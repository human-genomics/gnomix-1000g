"""Trio phasing benchmark: does Gnofix improve phase and haplotype ancestry?

Truth. The 1000 Genomes panel phased trio children with their parents
(SHAPEIT2-duohmm), so a child's haplotypes are close to the transmitted ones.
We treat that phase as the truth T.

Test. Children are split into TRIO_GROUPS groups (siblings together). For
each group we remove the children's phase and rephase their genotypes with
Beagle 5.5 against the rest of the panel; all families of the group are left
out of the reference. This gives the phase quality of an unrelated sample, R.
We then run Gnofix on R (F), and on T (TF: does Gnofix damage a correct phase?).

Metrics per child and chromosome, all at the model SNPs:
  - switch errors at the child's heterozygous sites, relative to T;
  - haplotype ancestry error: cM-weighted fraction of windows where the
    haplotype calls differ from the calls on T, after choosing the better of the
    two homolog orientations per chromosome;
  - the same error restricted to windows where the two truth haplotypes have
    different ancestry (the only windows where phase can matter);
  - diploid ancestry error (orientation-free).

Decision. Per population with at least MIN_CHILDREN trio children; other
populations follow their group (AMR, AFR-American, AFR, EUR, EAS, SAS). Gnofix
is used if it lowers the continental haplotype ancestry error (AFR + AHG → AFR,
EUR + WAS → EUR) versus R by at least MIN_GAIN, with paired one-sided
Wilcoxon p < 0.01. Trio children keep their published phase; their parents and
duo members follow the decision (their published phase is statistical quality,
see pedigree.py). The same samples get Gnofix phase in the corrected panel, so
tracts and phased genotypes agree; switch errors and phase accuracy by distance
are reported as the cost and gain of that choice. The 8-label error is reported
too, but EUR-WAS label noise, which phase correction cannot fix, dominates it.

Outputs (output/validation/): trio_benchmark_per_child.tsv,
trio_benchmark_per_child_genome.tsv, trio_benchmark_summary.tsv,
trio_benchmark.png, gnofix_policy.json
"""

from __future__ import annotations

import gzip
import os
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

from .common import (ADMIXED_POPULATIONS, ANCESTRIES, CONTINENTAL, OUTPUT, REFERENCE, TOOLS_DIR, WORK, atomic_path, group_of, log, read_json,
                     single_thread_env, threads, workers_for_memory, write_json)
from .haplotypes import load_match, open_reader, read_haplotypes
from .prepare import load_model, model_meta, variant_count, write_samples

# 4-byte genotype cells ("0|1" + tab) so rows are built with numpy, not Python strings.
PHASED_CELLS = np.frombuffer(b"0|0\t0|1\t1|0\t1|1\t.|.\t", dtype=np.uint8).reshape(5, 4)
UNPHASED_CELLS = np.frombuffer(b"0/0\t0/1\t0/1\t1/1\t./.\t", dtype=np.uint8).reshape(5, 4)


def use_published_policy() -> None:
    """Copy the committed policy and benchmark results of the published run into this run's outputs.

    RUN_TRIO_BENCHMARK=1 repeats the benchmark. Beagle's multithreaded phasing is not
    bit-identical between runs, so a rerun gives slightly different numbers; borderline
    decisions can change (PUR passed p < 0.01 in one of our two runs; Gnofix
    lowered its error in both, and the published policy uses Gnofix for PUR).
    """
    import shutil
    vdir = OUTPUT / "validation"
    current = vdir / "gnofix_policy.json"
    if current.exists() and read_json(current).get("source") == "benchmark":
        log("  keeping the policy of the benchmark already run in this data folder")
        return
    if not (REFERENCE / "gnofix_policy.json").exists():
        raise FileNotFoundError("reference/gnofix_policy.json not found; run with RUN_TRIO_BENCHMARK=1")
    for name in ("trio_benchmark_summary.tsv", "trio_benchmark.png"):
        with atomic_path(vdir / name) as tmp:
            shutil.copyfile(REFERENCE / name, tmp)
    policy = read_json(REFERENCE / "gnofix_policy.json")
    policy["source"] = "published"
    write_json(current, policy)


def trio_chroms() -> list[int]:
    from .common import chrom_list
    return chrom_list(os.environ.get("TRIO_CHROMS", "18-22"))


# ---- grouping ----------------------------------------------------------------

def child_groups(samples: pd.DataFrame, n_groups: int) -> pd.DataFrame:
    """Assign every trio child to a group. Children who share a parent stay together."""
    kids = samples[samples.trio_child].copy()
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for _, r in kids.iterrows():
        for p in (r.PAT, r.MAT):
            parent[find(r.IID)] = find(p)
    kids["family"] = kids.IID.map(find)
    fams = sorted(kids.family.unique())
    rng = np.random.default_rng(TRIO_SEED)
    rng.shuffle(fams)
    kids["group"] = kids.family.map({f: i % n_groups for i, f in enumerate(fams)})
    return kids


# ---- VCF writing -----------------------------------------------------------------

def _vcf_header(chrom: int, ids) -> bytes:
    return ("##fileformat=VCFv4.2\n"
            f"##contig=<ID={chrom}>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n").encode()


def _rows(fixed: list, codes: np.ndarray, cells: np.ndarray) -> bytes:
    """VCF body lines from genotype codes (0..3, 4 = missing)."""
    body = cells[codes].reshape(codes.shape[0], -1).copy()
    body[:, -1] = ord("\n")
    return b"".join(f + r.tobytes() for f, r in zip(fixed, body))


def write_vcfs(chrom: int, samples: pd.DataFrame, kids: pd.DataFrame, outdir) -> dict:
    """Reference VCF (all samples, phased) and one unphased query VCF per group.

    Markers are the model SNPs found in the panel, at their GRCh37 positions and
    with the model REF/ALT, so Beagle output lines up with the model directly.
    """
    meta = model_meta(chrom)
    match = load_match(chrom)
    pos, ref, alt = meta["snp_pos"], meta["snp_ref"], meta["snp_alt"]
    model_idx = match["model_idx"]
    # The model has a few duplicated positions; matching drops both copies.
    if not (np.all(np.diff(model_idx) > 0) and np.all(np.diff(pos[model_idx]) > 0)):
        raise RuntimeError("matched model SNPs must be sorted by GRCh37 position")
    order = np.argsort(match["pgen_idx"])          # match rows in panel order
    inv = np.empty(len(order), dtype=np.int64)
    inv[order] = np.arange(len(order))             # match row -> rank in panel order
    vids_sorted = match["pgen_idx"][order].astype(np.uint32)

    reader = open_reader(samples.pgen_index.values, variant_count(), len(samples))
    groups = sorted(kids.group.unique())
    kid_cols = {g: np.nonzero(samples.IID.isin(kids.IID[kids.group == g]).values)[0] for g in groups}
    paths = {"ref": outdir / "ref.vcf.gz", "map": outdir / "model.map"}
    fh_ref = gzip.open(str(paths["ref"]) + ".part", "wb", compresslevel=1)
    fh_ref.write(_vcf_header(chrom, samples.IID))
    fh_q = {}
    for g in groups:
        paths[g] = outdir / f"query_g{g}.vcf.gz"
        fh_q[g] = gzip.open(str(paths[g]) + ".part", "wb", compresslevel=1)
        fh_q[g].write(_vcf_header(chrom, ["Q_" + i for i in samples.IID.values[kid_cols[g]]]))

    n = len(samples)
    chunk = 20_000
    for s in range(0, len(model_idx), chunk):
        e = min(s + chunk, len(model_idx))
        ranks = np.sort(inv[s:e])
        buf = np.empty((e - s, 2 * n), dtype=np.int32)
        reader.read_alleles_list(vids_sorted[ranks], buf)
        a = buf[np.argsort(order[ranks])]          # back to model order s..e
        flip = match["flip"][s:e]
        a[flip] = np.where(a[flip] >= 0, 1 - a[flip], a[flip])
        codes = np.clip(a[:, 0::2], 0, 1) * 2 + np.clip(a[:, 1::2], 0, 1)
        codes[(a[:, 0::2] < 0) | (a[:, 1::2] < 0)] = 4
        mi = model_idx[s:e]
        fixed = [f"{chrom}\t{pos[m]}\t{chrom}:{pos[m]}\t{ref[m]}\t{alt[m]}\t.\tPASS\t.\tGT\t".encode() for m in mi]
        fh_ref.write(_rows(fixed, codes, PHASED_CELLS))
        for g in groups:
            fh_q[g].write(_rows(fixed, codes[:, kid_cols[g]], UNPHASED_CELLS))
    reader.close()
    fh_ref.close()
    os.replace(str(paths["ref"]) + ".part", paths["ref"])
    for g, f in fh_q.items():
        f.close()
        os.replace(str(paths[g]) + ".part", paths[g])

    cm = np.maximum.accumulate(np.interp(pos[model_idx], meta["gm_pos"], meta["gm_cm"]))
    # The map is written last: its existence marks the VCFs as complete.
    with atomic_path(paths["map"]) as tmp, open(tmp, "w") as fh:
        for p, c in zip(pos[model_idx], cm):
            fh.write(f"{chrom}\t.\t{c:.6f}\t{p}\n")
    return paths


def run_beagle(paths: dict, group: int, exclude_ids, outdir) -> tuple[str, str]:
    """Phase one group of children. Returns (output VCF, Beagle log text)."""
    out = outdir / f"phased_g{group}"
    if not (outdir / f"phased_g{group}.vcf.gz").exists():
        excl = outdir / f"exclude_g{group}.txt"
        excl.write_text("\n".join(exclude_ids) + "\n")
        mem = os.environ.get("BEAGLE_MEM_GB", "16")
        cmd = ["java", f"-Xmx{mem}g", "-jar", str(TOOLS_DIR / "beagle.jar"), f"gt={paths[group]}",
               f"ref={paths['ref']}", f"map={paths['map']}", f"out={out}.part", f"excludesamples={excl}",
               "impute=false", f"nthreads={threads()}", f"seed={TRIO_SEED}"]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
        os.replace(f"{out}.part.log", f"{out}.log")
        os.replace(f"{out}.part.vcf.gz", f"{out}.vcf.gz")
    return f"{out}.vcf.gz", open(f"{out}.log").read()


def beagle_reference_count(logtext: str) -> int:
    for line in logtext.splitlines():
        if line.startswith("Reference samples:"):
            return int(line.split(":")[1].replace(",", "").strip())
    raise RuntimeError("no 'Reference samples' line in Beagle log")


def read_phased(vcf: str, model_idx: np.ndarray, C: int, ids) -> np.ndarray:
    """Beagle output -> X (2k, C) in model order, absent model SNPs = REF. GT cells are 'a|b'."""
    from .haplotypes import FILL
    k = len(ids)
    X = np.full((2 * k, C), FILL, dtype=np.int8)
    rows = np.empty((len(model_idx), 4 * k - 1), dtype=np.uint8)
    r = 0
    with gzip.open(vcf, "rb") as fh:
        for line in fh:
            if line.startswith(b"##"):
                continue
            if line.startswith(b"#"):
                if line.rstrip(b"\n").split(b"\t")[9:] != [("Q_" + i).encode() for i in ids]:
                    raise RuntimeError("Beagle output samples differ from the query")
                continue
            tail = line.split(b"\t", 9)[9].rstrip(b"\n")
            if len(tail) != 4 * k - 1:
                raise RuntimeError("unexpected genotype field width in Beagle output")
            rows[r] = np.frombuffer(tail, dtype=np.uint8)
            r += 1
    if r != len(model_idx):
        raise RuntimeError(f"Beagle output has {r} markers, expected {len(model_idx)}")
    X[0::2][:, model_idx] = (rows[:, 0::4] - ord("0")).T
    X[1::2][:, model_idx] = (rows[:, 2::4] - ord("0")).T
    return X


# ---- ancestry and metrics -----------------------------------------------------------

_STATE: dict = {}


def _clear_state() -> None:
    import gc
    _STATE.clear()
    gc.collect()


def _init(chrom: int) -> None:
    import atexit
    single_thread_env()
    from .infer import window_of_snp
    m = load_model(chrom)
    atexit.register(_clear_state)
    _STATE.update(model=m, snp_window=window_of_snp(m.C, m.M, m.W))


def _gnofix(model, X, snp_window):
    from src.Gnofix.gnofix import gnofix
    B = model.base.predict_proba(X)
    _, _, ym, yp, _, tracker = gnofix(X[0], X[1], B=B.reshape(1, 2, model.W, -1)[0], smoother=model.smooth)
    swap = (np.asarray(tracker[0]) == 1)[snp_window]
    Xf = np.stack([np.where(swap, X[1], X[0]), np.where(swap, X[0], X[1])])
    return Xf, np.stack([ym, yp]), int((np.diff(np.asarray(tracker[0])) != 0).sum() + tracker[0][0])


LONG_RANGE_CM = (1, 5, 10, 20)
# Smallest mean reduction of the continental haplotype ancestry error that counts as a gain.
MIN_GAIN = 0.005
P_MAX = 0.01  # paired one-sided Wilcoxon test level
TRIO_SEED = 20261001  # Beagle seed and site-pair sampling
WORKER_GB = 1.5
LR_PAIRS = 2000  # heterozygous site pairs sampled per distance, child and chromosome
DISPLAY_ORDER = ["AFR-American", "ACB", "ASW", "AMR", "PEL", "MXL", "CLM", "PUR", "AFR", "EUR", "EAS", "SAS"]
# Populations with at least this many trio children get their own decision.
MIN_CHILDREN = 10
# Continental label index of each ancestry index (AFR + AHG → AFR, EUR + WAS → EUR).
CONTINENTAL_IDX = np.array([ANCESTRIES.index(CONTINENTAL[a]) for a in ANCESTRIES])


def long_range_accuracy(o: np.ndarray, cm: np.ndarray, d: float, rng) -> tuple[int, int]:
    """Pairs of heterozygous sites ~d cM apart: (correctly phased pairs, pairs).

    o[i] is True where the test haplotype 0 carries the truth haplotype 0 allele.
    A pair is correctly phased when both sites have the same orientation.
    """
    if len(o) < 2:
        return 0, 0
    i = rng.integers(0, len(o), size=min(LR_PAIRS, len(o)))
    j = np.searchsorted(cm, cm[i] + d)
    ok = j < len(o)
    i, j = i[ok], j[ok]
    return int((o[i] == o[j]).sum()), int(len(i))


def _child_task(XT: np.ndarray, XR: np.ndarray, wlen: np.ndarray, snp_cm: np.ndarray, seed: int) -> dict:
    model, snp_window = _STATE["model"], _STATE["snp_window"]
    LT = np.argmax(model.predict_proba(XT), axis=-1)
    LR = np.argmax(model.predict_proba(XR), axis=-1)
    XF, LF, nsw_f = _gnofix(model, XR, snp_window)
    XTF, LTF, nsw_tf = _gnofix(model, XT, snp_window)

    het = XT[0] != XT[1]  # absent model SNPs are REF on both haplotypes, so never heterozygous
    out = {"het_sites": int(het.sum()), "gnofix_swaps_on_rephased": nsw_f, "gnofix_swaps_on_truth": nsw_tf}
    het_cm = snp_cm[het]
    hetanc_snp = (LT[0] != LT[1])[snp_window][het]
    for name, X in (("R", XR), ("F", XF), ("TF", XTF)):
        o = X[0][het] == XT[0][het]
        out[f"switch_errors_{name}"] = int((o[1:] != o[:-1]).sum())
        for d in LONG_RANGE_CM:
            st = np.random.default_rng(seed + d)  # same pairs for R, F and TF
            good, n = long_range_accuracy(o, het_cm, d, st)
            out[f"lr{d}_good_{name}"], out[f"lr{d}_n_{name}"] = good, n
            oh, ch = o[hetanc_snp], het_cm[hetanc_snp]
            good, n = long_range_accuracy(oh, ch, d, np.random.default_rng(seed + 100 + d))
            out[f"lr{d}h_good_{name}"], out[f"lr{d}h_n_{name}"] = good, n
    hetanc = LT[0] != LT[1]
    tot = wlen.sum()
    for name, L in (("R", LR), ("F", LF), ("TF", LTF)):
        e_direct = ((L[0] != LT[0]).astype(float) + (L[1] != LT[1])) / 2
        e_swap = ((L[0] != LT[1]).astype(float) + (L[1] != LT[0])) / 2
        e = e_direct if (e_direct * wlen).sum() <= (e_swap * wlen).sum() else e_swap
        out[f"hap_err_{name}"] = float((e * wlen).sum() / tot)
        out[f"hap_err_hetanc_{name}"] = float((e * wlen)[hetanc].sum() / max(wlen[hetanc].sum(), 1e-9))
        dip = ~(((L[0] == LT[0]) & (L[1] == LT[1])) | ((L[0] == LT[1]) & (L[1] == LT[0])))
        out[f"dip_err_{name}"] = float((dip * wlen).sum() / tot)
        # continental labels (AFR + AHG → AFR, EUR + WAS → EUR)
        c, ct = CONTINENTAL_IDX[L], CONTINENTAL_IDX[LT]
        e1 = ((c[0] != ct[0]).astype(float) + (c[1] != ct[1])) / 2
        e2 = ((c[0] != ct[1]).astype(float) + (c[1] != ct[0])) / 2
        out[f"hap_err_cont_{name}"] = float(min((e1 * wlen).sum(), (e2 * wlen).sum()) / tot)
    out["hetanc_cM"] = float(wlen[hetanc].sum())
    out["chrom_cM"] = float(tot)
    for a, p in enumerate(ANCESTRIES):
        out[f"truth_{p}"] = float((((LT[0] == a).astype(float) + (LT[1] == a)) / 2 * wlen).sum() / tot)
    return out


# ---- main --------------------------------------------------------------------------

def run() -> None:
    from .tracts import window_lengths_cm, window_table

    vdir = OUTPUT / "validation"
    vdir.mkdir(parents=True, exist_ok=True)
    per_child_path = vdir / "trio_benchmark_per_child.tsv"
    samples = write_samples()
    n_groups = int(os.environ.get("TRIO_GROUPS", "6"))
    kids = child_groups(samples, n_groups)
    chroms = trio_chroms()
    # The benchmark chromosomes need models and matched SNPs even if CHROMS leaves them out.
    from .prepare import match_all, unpack_models, write_model_meta
    unpack_models(chroms)
    write_model_meta(chroms)
    match_all(chroms)
    log(f"  {len(kids)} trio children in {kids.group.nunique()} groups; chromosomes {chroms}")

    done = pd.read_csv(per_child_path, sep="\t") if per_child_path.exists() else pd.DataFrame()
    rows = [done] if len(done) else []
    for c in chroms:
        if len(done) and (done.chrom == c).any():
            log(f"  [skip] chr{c} benchmark")
            continue
        t0 = time.time()
        wdir = WORK / "trio" / f"chr{c}"
        wdir.mkdir(parents=True, exist_ok=True)
        marker = wdir / "groups"
        if marker.exists() and int(marker.read_text()) != n_groups:
            raise RuntimeError(f"{wdir} holds files for TRIO_GROUPS={marker.read_text().strip()}; "
                               f"use that value or move the folder away")
        marker.write_text(f"{n_groups}\n")
        if not (wdir / "model.map").exists():
            log(f"  chr{c}: writing reference and query VCFs")
            paths = write_vcfs(c, samples, kids, wdir)
        else:
            paths = {"ref": wdir / "ref.vcf.gz", "map": wdir / "model.map"}
            paths.update({g: wdir / f"query_g{g}.vcf.gz" for g in sorted(kids.group.unique())})
        meta = model_meta(c)
        match = load_match(c)
        wlen = window_lengths_cm(window_table(c))
        snp_cm = np.interp(meta["snp_pos"], meta["gm_pos"], meta["gm_cm"])
        XR_all = {}
        for g in sorted(kids.group.unique()):
            gk = kids[kids.group == g]
            gk_ids = samples.IID[samples.IID.isin(gk.IID)].tolist()  # panel order
            fam = set(gk.IID) | set(gk.PAT) | set(gk.MAT)
            vcf, logtext = run_beagle(paths, g, sorted(fam), wdir)
            n_ref = beagle_reference_count(logtext)
            if n_ref != len(samples) - len(fam):
                raise RuntimeError(f"chr{c} group {g}: Beagle used {n_ref} reference samples, expected "
                                   f"{len(samples) - len(fam)} (family excluded)")
            X = read_phased(vcf, match["model_idx"], int(meta["C"]), gk_ids)
            for k, iid in enumerate(gk_ids):
                XR_all[iid] = X[2 * k:2 * k + 2]
            log(f"    chr{c} group {g}: rephased {len(gk_ids)} children ({time.time() - t0:.0f} s)")
        ids = samples.IID[samples.IID.isin(kids.IID)].tolist()
        idx = samples.pgen_index[samples.IID.isin(kids.IID)].values
        XT = read_haplotypes(c, idx, int(meta["C"]), match, variant_count(), len(samples))
        # A dead worker raises BrokenProcessPool here instead of hanging; rerun to resume.
        with ProcessPoolExecutor(workers_for_memory(WORKER_GB), mp_context=get_context("spawn"),
                                 initializer=_init, initargs=(c,)) as ex:
            futs = [ex.submit(_child_task, XT[2 * k:2 * k + 2], XR_all[iid], wlen, snp_cm, 1000 * c + k)
                    for k, iid in enumerate(ids)]
            res_rows = []
            for iid, f in zip(ids, futs):
                r = f.result()
                r.update(sample=iid, chrom=c)
                res_rows.append(r)
        rows.append(pd.DataFrame(res_rows))
        cur = pd.concat(rows, ignore_index=True)
        with atomic_path(per_child_path) as tmp:
            cur.to_csv(tmp, sep="\t", index=False)
        log(f"  chr{c}: benchmark done ({time.time() - t0:.0f} s)")

    df = pd.concat(rows, ignore_index=True)
    info = samples.set_index("IID")
    df["Population"] = df["sample"].map(info.Population)
    df["SuperPop"] = df["sample"].map(info.SuperPop)
    df["group"] = [group_of(p, sp) for p, sp in zip(df.Population, df.SuperPop)]
    summarize(df)


def per_child(df: pd.DataFrame) -> pd.DataFrame:
    """Per child over all benchmark chromosomes: summed counts, cM-weighted error rates."""
    w = df.chrom_cM
    agg = df.assign(**{c + "_w": df[c] * w for c in df.columns if c.startswith(("hap_err_", "dip_err_"))},
                    **{c + "_wh": df[c] * df.hetanc_cM for c in df.columns if c.startswith("hap_err_hetanc_")})
    g = agg.groupby(["sample", "group", "Population"])
    child = g[[c for c in agg.columns if c.startswith(("switch_errors_", "het_sites", "gnofix_swaps", "lr"))]].sum()
    for name in ("R", "F", "TF"):
        for d in LONG_RANGE_CM:
            for h in ("", "h"):
                child[f"lr{d}{h}_acc_{name}"] = child[f"lr{d}{h}_good_{name}"] / child[f"lr{d}{h}_n_{name}"].clip(lower=1)
    for c in [c for c in df.columns if c.startswith(("hap_err_", "dip_err_")) and "hetanc" not in c]:
        child[c] = g[c + "_w"].sum() / g["chrom_cM"].sum()
    for c in [c for c in df.columns if c.startswith("hap_err_hetanc_")]:
        child[c] = g[c + "_wh"].sum() / g["hetanc_cM"].sum().clip(lower=1e-9)
    child["hetanc_fraction"] = g["hetanc_cM"].sum() / g["chrom_cM"].sum()
    return child.reset_index()


def improves(d: pd.DataFrame, worse: str, better: str) -> float:
    """Paired one-sided Wilcoxon p-value that column `better` is smaller than `worse`."""
    if ((d[worse] - d[better]) != 0).sum() < 5:
        return 1.0
    return float(wilcoxon(d[worse], d[better], alternative="greater", zero_method="wilcox").pvalue)


def summary_row(level: str, name: str, d: pd.DataFrame) -> dict:
    """Benchmark summary and Gnofix decision for one group or population."""
    p_cont = improves(d, "hap_err_cont_R", "hap_err_cont_F")
    use = p_cont < P_MAX and d.hap_err_cont_R.mean() - d.hap_err_cont_F.mean() >= MIN_GAIN
    row = {"level": level, "group": name, "children": len(d), "het_ancestry_fraction": d.hetanc_fraction.mean()}
    for col, key in (("hap_err", "hap_err_{}"), ("hap_err_hetanc", "hap_err_hetanc_{}"),
                     ("hap_err_continental", "hap_err_cont_{}"), ("dip_err", "dip_err_{}"),
                     ("switch_errors", "switch_errors_{}")):
        for tag, lab in (("R", "rephased"), ("F", "gnofix"), ("TF", "gnofix_on_truth")):
            if key.format(tag) in d:
                row[f"{col}_{lab}"] = d[key.format(tag)].mean()
    row["SER_rephased"] = d.switch_errors_R.sum() / d.het_sites.sum()
    row["SER_gnofix"] = d.switch_errors_F.sum() / d.het_sites.sum()
    for dd in LONG_RANGE_CM:
        for tag, lab in (("R", "rephased"), ("F", "gnofix"), ("TF", "gnofix_on_truth")):
            row[f"phase_acc_{dd}cM_{lab}"] = d[f"lr{dd}_good_{tag}"].sum() / max(d[f"lr{dd}_n_{tag}"].sum(), 1)
        for tag, lab in (("R", "rephased"), ("F", "gnofix")):
            row[f"phase_acc_{dd}cM_hetanc_{lab}"] = d[f"lr{dd}h_good_{tag}"].sum() / max(d[f"lr{dd}h_n_{tag}"].sum(), 1)
    row.update({"p_continental_err_improves": p_cont,
                "p_hap_err_improves": improves(d, "hap_err_R", "hap_err_F"),
                "p_10cM_phase_hetanc_improves": improves(d, "lr10h_acc_F", "lr10h_acc_R"),
                "p_switch_improves": improves(d, "switch_errors_R", "switch_errors_F"),
                "gnofix": bool(use)})
    return row


def summarize(df: pd.DataFrame) -> None:
    """Write the benchmark summary, the Gnofix policy and the figure."""
    vdir = OUTPUT / "validation"
    child = per_child(df)
    units = [("group", grp, d) for grp, d in child.groupby("group")]
    units += [("population", pop, d) for pop, d in child.groupby("Population") if len(d) >= MIN_CHILDREN]
    summ = pd.DataFrame([summary_row(*u) for u in units])
    summ["o"] = summ.group.map({k: i for i, k in enumerate(DISPLAY_ORDER)}).fillna(len(DISPLAY_ORDER))
    summ = summ.sort_values(["o", "level", "group"]).drop(columns="o")
    policy = {"rule": (f"Decided per population (populations with at least {MIN_CHILDREN} trio children; others "
                       "follow their group). Gnofix is used if it lowers the continental haplotype ancestry error "
                       f"(AFR + AHG → AFR, EUR + WAS → EUR) versus the rephased baseline by at least {MIN_GAIN:.3f} "
                       f"(mean over children; paired one-sided Wilcoxon p < {P_MAX}). The same samples get Gnofix "
                       "phase in the corrected panel. Trio children always keep their published phase (phased "
                       "with both parents); trio parents and duo members follow their population."),
              "apply_gnofix": {r.group: r.gnofix for r in summ.itertuples() if r.level == "group"},
              "apply_gnofix_population": {r.group: r.gnofix for r in summ.itertuples() if r.level == "population"},
              "chromosomes": sorted(int(c) for c in df.chrom.unique()),
              "source": "benchmark"}
    with atomic_path(vdir / "trio_benchmark_per_child_genome.tsv") as tmp:
        child.to_csv(tmp, sep="\t", index=False)
    with atomic_path(vdir / "trio_benchmark_summary.tsv") as tmp:
        summ.to_csv(tmp, sep="\t", index=False, float_format="%.5g")
    write_json(vdir / "gnofix_policy.json", policy)
    shown = summ[summ.group.isin(ADMIXED_POPULATIONS) |
                 ((summ.level == "group") & summ.group.isin(["AFR", "EUR", "EAS", "SAS"]))]
    plot_summary(shown, vdir / "trio_benchmark.png")
    log(summ[["level", "group", "children", "hap_err_continental_rephased", "hap_err_continental_gnofix",
              "p_continental_err_improves", "gnofix"]].to_string(index=False, float_format=lambda x: f"{x:.4g}"))


def plot_summary(summ: pd.DataFrame, path) -> None:
    """Two panels: continental haplotype ancestry error and 10 cM phase accuracy, rephased vs Gnofix."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator, MultipleLocator

    colors = {"rephased": "#A3A3A3", "gnofix": "#1F6F8B"}
    names = {"rephased": "Rephased without parents (statistical phase)", "gnofix": "Rephased, then Gnofix"}
    panels = [("hap_err_continental_{}", "Continental haplotype ancestry error", "%"),
              ("phase_acc_10cM_hetanc_{}", "Phase accuracy, het. sites 10 cM apart, mixed-ancestry regions", "%")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4), dpi=130)
    x = np.arange(len(summ))
    for ax, (col, title, unit) in zip(axes, panels):
        for i, cond in enumerate(("rephased", "gnofix")):
            ax.bar(x + (i - 0.5) * 0.38, summ[col.format(cond)].values * 100, width=0.36, color=colors[cond],
                   label=names[cond], linewidth=0)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{g}\n{n}" for g, n in zip(summ.group, summ.children)], fontsize=9)
        ax.set_title(title, fontsize=10.5, loc="left")
        ax.set_ylabel(unit, fontsize=10)
        ax.yaxis.set_major_locator(MultipleLocator(10) if "phase" in col else MaxNLocator(integer=True))
        ax.grid(axis="y", color="#e6e6e6", linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis="y", labelsize=9)
    axes[1].axhline(50, color="#666666", linewidth=0.8, linestyle="--")
    handles, labels = axes[0].get_legend_handles_labels()
    handles.append(Line2D([], [], color="#666666", linewidth=0.8, linestyle="--"))
    labels.append("50% = random phase")
    fig.text(0.005, 0.115, "children:", fontsize=9, color="#555555")
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, fontsize=9.5)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    with atomic_path(path) as tmp:
        fig.savefig(tmp, format="png")
    plt.close(fig)
