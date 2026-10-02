"""Read phased haplotypes at the model SNPs from the PLINK 2 panel."""

from __future__ import annotations

import numpy as np

from .common import WORK

# Model SNPs absent from the panel are filled with the reference allele (0).
# Gnomix's own missing code (2) is out of range for the pretrained logistic
# base models: with 5% of model SNPs coded 2, calls on the Gnomix demo agree
# with the complete-data calls in only 79% of windows; with REF fill, 100%
# (validation/missing_snp_fill.tsv). Absent SNPs are mostly rare in 1000 Genomes.
FILL = 0


def load_match(chrom: int) -> dict:
    with np.load(WORK / "match" / f"chr{chrom}.npz") as z:
        return {k: z[k] for k in z.files}


def open_reader(sample_idx: np.ndarray, variant_ct: int, raw_sample_ct: int):
    import pgenlib

    from .prepare import allele_offsets

    return pgenlib.PgenReader(str(WORK / "all_hg38.pgen").encode(), raw_sample_ct=raw_sample_ct,
                              variant_ct=variant_ct, sample_subset=np.asarray(sample_idx, dtype=np.uint32),
                              allele_idx_offsets=allele_offsets())


def read_haplotypes(chrom: int, sample_idx: np.ndarray, n_model_snps: int, match: dict,
                    variant_ct: int, raw_sample_ct: int, chunk: int = 50_000) -> np.ndarray:
    """Return X with shape (2 * n_samples, n_model_snps), int8, in model SNP order.

    Rows are haplotypes: sample k has rows 2k and 2k+1 (first and second allele of
    the phased genotype). Alleles are coded 0/1 relative to the model REF/ALT.
    Model SNPs that are absent from the panel, and missing calls, get FILL (REF).
    """
    sample_idx = np.asarray(sample_idx, dtype=np.uint32)
    if np.any(np.diff(sample_idx.astype(np.int64)) <= 0):
        raise ValueError("sample_idx must be strictly increasing")
    n = len(sample_idx)
    reader = open_reader(sample_idx, variant_ct, raw_sample_ct)
    X = np.full((2 * n, n_model_snps), FILL, dtype=np.int8)
    order = np.argsort(match["pgen_idx"], kind="stable")
    vids = match["pgen_idx"][order].astype(np.uint32)
    cols = match["model_idx"][order]
    flips = match["flip"][order]
    buf = np.empty((chunk, 2 * n), dtype=np.int32)
    for s in range(0, len(vids), chunk):
        e = min(s + chunk, len(vids))
        out = buf[: e - s]
        reader.read_alleles_list(vids[s:e], out)
        a = out.T.astype(np.int8)  # (2n, k)
        miss = a < 0
        f = flips[s:e]
        a[:, f] = 1 - a[:, f]
        a[miss] = FILL
        X[:, cols[s:e]] = a
    reader.close()
    return X
