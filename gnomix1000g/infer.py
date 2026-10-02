"""Run the pretrained Gnomix model, with and without Gnofix, on every sample.

For each chromosome and batch of samples:
  1. read phased haplotypes at the model SNPs (X),
  2. base model -> per-window ancestry probabilities (B),
  3. smoother -> posterior per window (the "raw" call, published 1000G phase),
  4. Gnofix (upstream gnofix() with the defaults of Gnomix.phase) per sample,
     then the model again on the corrected haplotypes (the "gnofix" call).

As in gnomix.py run_inference(phase=True), Gnofix window labels come from
Gnofix itself and posteriors from model.predict_proba on the corrected
haplotypes. The two disagree in about 1% of windows; tract posteriors are the
posterior of the Gnofix label. We also keep which windows Gnofix exchanged
("swap"), so the same correction can be applied to every variant (pfile.py).

Outputs: work/infer/chrN/batch_*.npz while running, then work/infer/chrN.npz.
"""

from __future__ import annotations

import gc
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from multiprocessing import get_context

import numpy as np

from .common import ADMIXED_POPULATIONS, WORK, atomic_path, log, single_thread_env, workers_for_memory
from .haplotypes import load_match, read_haplotypes
from .prepare import load_model, model_meta, variant_count, write_samples

# Peak memory of one worker: chromosome 1-2 models while unpickling.
WORKER_GB = 2.5
CALL_KEYS = ("sample_idx", "lab_raw", "prob_raw", "lab_fix", "prob_fix", "swap")

_STATE: dict = {}


def quantize(p: np.ndarray) -> np.ndarray:
    return np.rint(np.clip(p, 0, 1) * 255).astype(np.uint8)


def window_of_snp(C: int, M: int, W: int) -> np.ndarray:
    """Window index of each model SNP. The last window takes the remainder."""
    return np.minimum(np.arange(C) // M, W - 1)


def _clear_state() -> None:
    _STATE.clear()
    gc.collect()


def _init(variant_ct: int, raw_sample_ct: int) -> None:
    import atexit
    single_thread_env()
    atexit.register(_clear_state)
    _STATE.update(variant_ct=variant_ct, raw_sample_ct=raw_sample_ct)


def _use_chrom(chrom: int) -> None:
    """Load the model of `chrom` in this worker (one model at a time)."""
    if _STATE.get("chrom") == chrom:
        return
    _STATE.pop("model", None)
    gc.collect()
    model = load_model(chrom)
    _STATE.update(chrom=chrom, model=model, match=load_match(chrom),
                  snp_window=window_of_snp(model.C, model.M, model.W))


def gnomix_calls(model, X: np.ndarray, snp_window: np.ndarray) -> tuple[dict, int]:
    """Raw and Gnofix calls for phased haplotypes X (2n, C): labels, quantized
    posteriors and Gnofix swaps per window, plus the number of alleles upstream
    Gnofix would have placed differently (window remainder offset, see below)."""
    from src.Gnofix.gnofix import gnofix  # upstream Gnofix, unmodified

    n, W = X.shape[0] // 2, model.W
    B = model.base.predict_proba(X)  # (2n, W, A)
    P_raw = model.smooth.predict_proba(B)
    X_fix = X.copy()
    lab_fix = np.empty((2 * n, W), dtype=np.uint8)
    swap = np.zeros((n, W), dtype=bool)
    upstream_offset = 0
    for k, Bk in enumerate(B.reshape(n, 2, W, -1)):
        xm, _, ym, yp, _, tracker = gnofix(X[2 * k], X[2 * k + 1], B=Bk, smoother=model.smooth)
        lab_fix[2 * k], lab_fix[2 * k + 1] = ym, yp
        swap[k] = np.asarray(tracker[0]) == 1
        # Gnofix decides on window probabilities only, and its tracker records
        # which windows it exchanged. Upstream then applies the swaps at SNP
        # index w * (C // W). When the last window holds a remainder, C // W is
        # larger than the model window M (chr6, 7, 8, 19, 21) and the swaps
        # drift. We rebuild the haplotypes from the tracker with the true layout.
        s = swap[k][snp_window]
        X_fix[2 * k] = np.where(s, X[2 * k + 1], X[2 * k])
        X_fix[2 * k + 1] = np.where(s, X[2 * k], X[2 * k + 1])
        upstream_offset += int((X_fix[2 * k] != xm).sum())
    P_fix = model.predict_proba(X_fix)
    out = {"lab_raw": np.argmax(P_raw, axis=-1).astype(np.uint8), "prob_raw": quantize(P_raw),
           "lab_fix": lab_fix, "prob_fix": quantize(P_fix), "swap": swap}
    return out, upstream_offset


def run_batch(chrom: int, sample_idx: np.ndarray, out_path: str) -> dict:
    """Raw and Gnofix calls for one batch of samples on one chromosome."""
    _use_chrom(chrom)
    model = _STATE["model"]
    X = read_haplotypes(chrom, sample_idx, model.C, _STATE["match"],
                        _STATE["variant_ct"], _STATE["raw_sample_ct"])
    out, upstream_offset = gnomix_calls(model, X, _STATE["snp_window"])
    with atomic_path(out_path) as tmp:
        np.savez(tmp, sample_idx=sample_idx, **out)
    return {"chrom": chrom, "upstream_offset_snps": upstream_offset}


def batches_of(samples, batch_size: int) -> list[np.ndarray]:
    idx = samples.pgen_index.values
    return [idx[i:i + batch_size] for i in range(0, len(idx), batch_size)]


def _batch_dir(chrom: int, batch_size: int):
    """Batch files of a chromosome; refuses to mix runs with different batch sizes."""
    bdir = WORK / "infer" / f"chr{chrom}"
    bdir.mkdir(parents=True, exist_ok=True)
    marker = bdir / "batch_size"
    if marker.exists() and int(marker.read_text()) != batch_size:
        raise RuntimeError(f"{bdir} holds batches of {marker.read_text().strip()} samples, but BATCH_SIZE="
                           f"{batch_size}. Use the earlier BATCH_SIZE or move that folder away.")
    marker.write_text(f"{batch_size}\n")
    return bdir


def assemble(chrom: int, batch_size: int, idx: np.ndarray) -> None:
    """Join the batch files of a chromosome in panel order, then remove them."""
    bdir = WORK / "infer" / f"chr{chrom}"
    n_batches = -(-len(idx) // batch_size)
    parts = []
    for b in range(n_batches):
        with np.load(bdir / f"batch_{b:05d}.npz") as z:
            parts.append({k: z[k] for k in CALL_KEYS})
    out = {k: np.concatenate([p[k] for p in parts]) for k in CALL_KEYS}
    if not np.array_equal(out["sample_idx"], idx):
        raise RuntimeError(f"chr{chrom}: batch files do not cover the samples in panel order")
    with atomic_path(WORK / "infer" / f"chr{chrom}.npz") as tmp:
        np.savez(tmp, **out)
    shutil.rmtree(bdir)  # the assembled file holds the same data


def run(chroms: list[int], batch_size: int) -> None:
    """All chromosomes in one worker pool. Expensive tasks (Gnofix on admixed
    samples, long chromosomes) go first so that no core waits at the end."""
    samples = write_samples()
    idx = samples.pgen_index.values
    batches = batches_of(samples, batch_size)
    admixed = samples.Population.isin(ADMIXED_POPULATIONS).values
    has_admixed = [admixed[i:i + batch_size].any() for i in range(0, len(idx), batch_size)]
    tasks, todo = [], []
    for c in chroms:
        if (WORK / "infer" / f"chr{c}.npz").exists():
            log(f"  [skip] chr{c} inference")
            continue
        todo.append(c)
        bdir = _batch_dir(c, batch_size)
        W = int(model_meta(c)["W"])
        for b, s in enumerate(batches):
            out = bdir / f"batch_{b:05d}.npz"
            if not out.exists():
                # admixed batches (slow Gnofix) first, then longer chromosomes first
                tasks.append(((not has_admixed[b], -W, c, b), (c, s, str(out))))
    if not todo:
        return
    tasks.sort(key=lambda t: t[0])
    args = [a for _, a in tasks]
    log(f"  {len(args)} tasks ({len(batches)} batches of {batch_size} samples per chromosome)")
    remaining = {c: sum(1 for a in args if a[0] == c) for c in todo}
    offset = dict.fromkeys(todo, 0)
    n_workers = workers_for_memory(WORKER_GB)
    t0, done, attempts, total = time.time(), 0, 0, len(args)
    while args:
        # Fresh workers for every wave of tasks: bounds memory growth. A worker
        # that dies (for example killed for memory) breaks only its wave; the
        # unfinished tasks are retried in the next wave instead of hanging.
        wave, args = args[:n_workers * 8], args[n_workers * 8:]
        failed = []
        with ProcessPoolExecutor(n_workers, mp_context=get_context("spawn"), initializer=_init,
                                 initargs=(variant_count(), len(samples))) as ex:
            futures = {ex.submit(run_batch, *a): a for a in wave}
            for f in as_completed(futures):
                try:
                    r = f.result()
                except BrokenProcessPool:
                    failed.append(futures[f])
                    continue
                done += 1
                offset[r["chrom"]] += r["upstream_offset_snps"]
                remaining[r["chrom"]] -= 1
                if remaining[r["chrom"]] == 0:
                    log(f"    chr{r['chrom']} done ({time.time() - t0:.0f} s)")
                if done % max(1, total // 20) == 0:
                    log(f"    {done}/{total} tasks, {time.time() - t0:.0f} s")
        if failed:
            attempts += 1
            if attempts > 3:
                raise RuntimeError(f"{len(failed)} inference tasks failed repeatedly (worker processes died; "
                                   "out of memory?). Lower THREADS.")
            n_workers = max(1, n_workers // 2)
            log(f"    {len(failed)} tasks lost when a worker died; retrying them with {n_workers} workers")
            args = failed + args
    for c in todo:
        if offset[c]:
            log(f"    chr{c}: upstream Gnofix placed {offset[c]:,} haplotype alleles differently "
                f"(window remainder offset); corrected haplotypes rebuilt from the tracker")
        assemble(c, batch_size, idx)
        log(f"  chr{c}: inference assembled")
