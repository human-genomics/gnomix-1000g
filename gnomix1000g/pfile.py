"""Write the 1000 Genomes panel (autosomes) with Gnofix phase corrections applied.

Gnofix decides, per sample, which model windows have their two haplotypes
exchanged (work/infer/chrN.npz, "swap"). Here we apply the same exchanges to
every variant of the panel, not only to the model SNPs:

  - Corrections are applied to exactly the samples whose final tracts use
    Gnofix (never trio or duo members), so haplotypes 1 and 2 of this panel are
    haplotypes A and B of the tracts and karyograms. Other samples are copied.
  - Each variant takes the swap state of the window of its nearest lifted model
    SNP (GRCh38). This does not assume that liftover keeps window order.
    Variants between two windows with different states get the state of the
    nearer one, so their phase relative to each other is uncertain.
  - After writing, every variant is read back and compared with the source:
    unordered genotypes must be identical, and phased genotypes must equal the
    planned swaps.

Output: output/pfile/chrN_gnofix.{pgen,pvar.zst,psam}, swaps.tsv.gz,
        verification.json, and the per-chromosome records they are merged from
        (output/pfile/parts/).
"""

from __future__ import annotations

import io
import shutil
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context

import numpy as np
import pandas as pd

from .common import (DOWNLOADS, OUTPUT, WORK, atomic_path, inputs_stamp, log, read_json, stamp_is_current,
                     workers_for_memory, write_json, write_stamp)
from .infer import window_of_snp
from .prepare import allele_offsets, model_meta, variant_count, write_samples
from .tracts import final_uses_gnofix

CHUNK = 8192  # variants per read/write block
WORKER_GB = 2.0
PFILE = OUTPUT / "pfile"


def nearest_window(snp_pos38: np.ndarray, snp_window: np.ndarray, pos: np.ndarray) -> np.ndarray:
    """Window of the nearest model SNP for each position (GRCh38; unlifted SNPs have pos <= 0).

    The SNPs need not be in GRCh38 order: liftover can reverse segments.
    """
    lifted = np.nonzero(snp_pos38 > 0)[0]
    order = lifted[np.argsort(snp_pos38[lifted], kind="stable")]
    sp, sw = snp_pos38[order], snp_window[order]
    right = np.clip(np.searchsorted(sp, pos), 0, len(sp) - 1)
    left = np.clip(right - 1, 0, None)
    use_left = np.abs(pos - sp[left]) <= np.abs(sp[right] - pos)
    return np.where(use_left, sw[left], sw[right])


def variant_windows(chrom: int, pos: np.ndarray) -> np.ndarray:
    meta = model_meta(chrom)
    with np.load(WORK / "match" / f"chr{chrom}.npz") as z:
        pos38 = z["pos38"]
    return nearest_window(pos38, window_of_snp(int(meta["C"]), int(meta["M"]), int(meta["W"])), pos)


def apply_swaps(alleles: np.ndarray, phased: np.ndarray, states: np.ndarray, samples: np.ndarray) -> None:
    """Exchange the two alleles of `samples` where `states` (variants x samples) is True (in place)."""
    for j, k in enumerate(samples):
        m = states[:, j] & (phased[:, k] == 1)
        a0 = alleles[m, 2 * k].copy()
        alleles[m, 2 * k] = alleles[m, 2 * k + 1]
        alleles[m, 2 * k + 1] = a0


def write_pvars(chroms: list[int]) -> dict:
    """Write chrN_gnofix.pvar.zst and return the record range [start, end) of each chromosome."""
    import zstandard

    ranges, writers, header = {}, {}, []
    i, cur, seen = 0, None, set()
    with open(DOWNLOADS / "all_hg38.pvar.zst", "rb") as fh:
        text = io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(fh), encoding="ascii")
        for line in text:
            if line.startswith("#"):
                header.append(line)
                continue
            c = line.split("\t", 1)[0]
            if c != cur:
                if c in seen:
                    raise RuntimeError(f"chromosome {c} is not contiguous in the pvar")
                seen.add(c)
                if cur in writers:
                    writers[cur][1].close()  # chromosomes are contiguous: this one is complete
                cur = c
                if c.isdigit() and int(c) in chroms:
                    path = PFILE / f".part.chr{c}_gnofix.pvar.zst"
                    w = zstandard.ZstdCompressor(level=6).stream_writer(open(path, "wb"))
                    w.write("".join(header).encode())
                    writers[c] = (path, w)
                    ranges[int(c)] = [i, i]
            if c in writers:
                writers[c][1].write(line.encode())
                ranges[int(c)][1] = i + 1
            i += 1
    if cur in writers:
        writers[cur][1].close()
    for c, (path, _) in writers.items():
        path.replace(PFILE / f"chr{c}_gnofix.pvar.zst")
    return ranges


def write_chrom(chrom: int, start: int, end: int, apply: np.ndarray) -> dict:
    import pgenlib
    import zstandard

    samples = write_samples()
    n = len(samples)
    with np.load(WORK / "infer" / f"chr{chrom}.npz") as z:
        swap = z["swap"]
    corrected = np.nonzero(apply & swap.any(axis=1))[0]
    with open(PFILE / f"chr{chrom}_gnofix.pvar.zst", "rb") as fh:
        text = io.TextIOWrapper(zstandard.ZstdDecompressor().stream_reader(fh), encoding="ascii")
        pos = pd.read_csv(text, sep="\t", comment="#", header=None, usecols=[1], dtype=np.int64)[1].values
    if len(pos) != end - start:
        raise RuntimeError(f"chr{chrom}: pvar has {len(pos)} records, expected {end - start}")
    vwin = variant_windows(chrom, pos)

    offsets = np.asarray(allele_offsets()[start:end + 1], dtype=np.int64)
    allele_cts = (offsets[1:] - offsets[:-1]).astype(np.uint32)
    source = pgenlib.PgenReader(str(WORK / "all_hg38.pgen").encode(), raw_sample_ct=n,
                                variant_ct=variant_count(), allele_idx_offsets=allele_offsets())
    out_pgen = PFILE / f"chr{chrom}_gnofix.pgen"
    tmp_pgen = PFILE / f".part.chr{chrom}_gnofix.pgen"
    writer = pgenlib.PgenWriter(str(tmp_pgen).encode(), n, end - start, False, int(allele_cts.max()), True)
    swapped = hets = 0
    for s in range(start, end, CHUNK):
        e = min(s + CHUNK, end)
        a = np.empty((e - s, 2 * n), dtype=np.int32)
        ph = np.empty((e - s, n), dtype=np.uint8)
        source.read_alleles_and_phasepresent_range(s, e, a, ph)
        states = swap[corrected][:, vwin[s - start:e - start]].T
        hets += int((a[:, 0::2] != a[:, 1::2]).sum())
        swapped += int((states & (ph[:, corrected] == 1) & (a[:, 2 * corrected] != a[:, 2 * corrected + 1])).sum())
        apply_swaps(a, ph, states, corrected)
        writer.append_partially_phased_batch(a, ph, allele_cts[s - start:e - start])
    writer.close()

    # Read back every variant and compare with the source and the plan.
    check = pgenlib.PgenReader(str(tmp_pgen).encode(), raw_sample_ct=n, variant_ct=end - start,
                               allele_idx_offsets=(offsets - offsets[0]).astype(np.uintp))
    for s in range(start, end, CHUNK):
        e = min(s + CHUNK, end)
        a = np.empty((e - s, 2 * n), dtype=np.int32)
        ph = np.empty((e - s, n), dtype=np.uint8)
        source.read_alleles_and_phasepresent_range(s, e, a, ph)
        b = np.empty_like(a)
        ph2 = np.empty_like(ph)
        check.read_alleles_and_phasepresent_range(s - start, e - start, b, ph2)
        if not np.array_equal(np.sort(np.stack([a[:, 0::2], a[:, 1::2]]), axis=0),
                              np.sort(np.stack([b[:, 0::2], b[:, 1::2]]), axis=0)) or not np.array_equal(ph, ph2):
            raise RuntimeError(f"chr{chrom}: genotypes differ from the source near record {s}")
        apply_swaps(a, ph, swap[corrected][:, vwin[s - start:e - start]].T, corrected)
        if not np.array_equal(a, b):
            raise RuntimeError(f"chr{chrom}: phase differs from the planned corrections near record {s}")
    check.close()
    source.close()

    # Records first, then the psam, and the pgen last: the pgen marks the chromosome as done.
    rows = []
    for k in corrected:
        st = swap[k][vwin]
        ch = np.nonzero(np.diff(st.astype(np.int8), prepend=0) != 0)[0]
        rows += [(samples.IID[k], chrom, int(pos[v]), bool(st[v])) for v in ch]
    parts = PFILE / "parts"
    with atomic_path(parts / f"chr{chrom}.swaps.tsv") as tmp:
        pd.DataFrame(rows, columns=["sample", "chrom", "pos_hg38", "swapped_from_here"]).to_csv(
            tmp, sep="\t", index=False)
    result = {"chrom": chrom, "variants": end - start, "samples_corrected": int(len(corrected)),
              "switch_points": len(rows), "exchanged_het_calls": swapped, "het_calls": hets,
              "read_back": "all variants: unordered genotypes equal the source; phase equals the plan"}
    write_json(parts / f"chr{chrom}.json", result)
    with atomic_path(PFILE / f"chr{chrom}_gnofix.psam") as tmp:
        shutil.copyfile(DOWNLOADS / "all_hg38.psam", tmp)
    tmp_pgen.replace(out_pgen)
    return result


def run(chroms: list[int]) -> None:
    PFILE.mkdir(parents=True, exist_ok=True)
    samples = write_samples()
    apply = final_uses_gnofix(samples)
    log(f"  phase corrections for {apply.sum():,} samples (the samples whose final tracts use Gnofix)")
    stamp = inputs_stamp(chroms, apply)
    if stamp_is_current(PFILE, stamp):
        todo = [c for c in chroms if not (PFILE / f"chr{c}_gnofix.pgen").exists()]
    else:
        todo = list(chroms)  # calls or policy changed: rewrite everything
        write_stamp(PFILE, "")
    if not todo:
        log("  [skip] pfile written")
    else:
        ranges = write_pvars(todo)
        with ProcessPoolExecutor(min(workers_for_memory(WORKER_GB), len(todo)), mp_context=get_context("spawn")) as ex:
            jobs = {c: ex.submit(write_chrom, c, *ranges[c], apply) for c in todo}
            for c, j in jobs.items():
                r = j.result()
                log(f"  chr{c}: {r['variants']:,} variants, {r['samples_corrected']:,} samples corrected, "
                    f"{r['exchanged_het_calls']:,} of {r['het_calls']:,} heterozygous calls exchanged; "
                    f"read-back OK")
        write_stamp(PFILE, stamp)
    parts = PFILE / "parts"
    with atomic_path(PFILE / "swaps.tsv.gz") as tmp:
        pd.concat([pd.read_csv(parts / f"chr{c}.swaps.tsv", sep="\t") for c in chroms]).to_csv(
            tmp, sep="\t", index=False, compression="gzip")
    write_json(PFILE / "verification.json", {str(c): read_json(parts / f"chr{c}.json") for c in chroms})
