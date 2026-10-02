"""Checks of our model input against the official Gnomix path (chr22 demo data).

The Gnomix repository ships a GRCh37 query VCF (demo/data/small_query_chr22.vcf.gz,
nine 1000 Genomes Phase 3 samples, 100% of chr22 model SNPs). Official
inference on it (gnomix.py: vcf_to_npy + predict) gives the reference calls.

1. Missing-SNP handling. Mask the model SNPs that our panel lacks (and random
   5% / 20% subsets) and compare calls with Gnomix's missing code (2) versus
   filling with the reference allele (0).
2. Panel path. Build the same nine samples from the high-coverage GRCh38 panel
   with our liftover, matching and REF fill, and compare window calls with the
   official calls. Differences come from the different callsets and phasing of
   Phase 3 versus high coverage, so agreement is high but not exactly 100%.

Output: output/validation/missing_snp_fill.tsv, official_path_concordance.tsv
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .common import OUTPUT, TOOLS_DIR, atomic_path, import_gnomix, log
from .haplotypes import load_match, read_haplotypes
from .prepare import load_model, match_all, unpack_models, variant_count, write_model_meta, write_samples


def run() -> None:
    out = OUTPUT / "validation"
    out.mkdir(parents=True, exist_ok=True)
    if (out / "official_path_concordance.tsv").exists():
        log("  [skip] official-path checks")
        return
    import_gnomix()
    from src.utils import read_vcf, vcf_to_npy  # official Gnomix VCF encoding

    # The demo VCF is chr22: prepare chr22 even when CHROMS does not include it.
    unpack_models([22])
    write_model_meta([22])
    match_all([22])
    m = load_model(22)
    pops = np.asarray(m.population_order)
    vcf = read_vcf(str(TOOLS_DIR / "gnomix" / "demo" / "data" / "small_query_chr22.vcf.gz"), chm="22", fields="*")
    Xd = vcf_to_npy(vcf, m.snp_pos, m.snp_ref, verbose=False)  # official encoding
    Ld = np.argmax(m.predict_proba(Xd), -1)
    match = load_match(22)
    absent = np.ones(m.C, dtype=bool)
    absent[match["model_idx"]] = False

    rows = []
    rng = np.random.default_rng(1)
    masks = {"absent from 1000G high-coverage panel": absent,
             "random 1%": rng.random(m.C) < 0.01, "random 5%": rng.random(m.C) < 0.05,
             "random 20%": rng.random(m.C) < 0.20}
    for name, mask in masks.items():
        for fill, label in ((2, "Gnomix missing code 2"), (0, "reference allele 0")):
            X = Xd.copy()
            X[:, mask] = fill
            L = np.argmax(m.predict_proba(X), -1)
            rows.append({"masked_snps": name, "fraction_masked": round(mask.mean(), 4), "fill": label,
                         "window_agreement_with_complete_data": round(float((L == Ld).mean()), 4)})
    fill_df = pd.DataFrame(rows)
    with atomic_path(out / "missing_snp_fill.tsv") as tmp:
        fill_df.to_csv(tmp, sep="\t", index=False)
    log(fill_df.to_string(index=False))

    samples = write_samples().set_index("IID")
    ids = list(vcf["samples"])
    idx = samples.loc[ids, "pgen_index"].values
    order = np.argsort(idx)
    Xs = read_haplotypes(22, idx[order], m.C, match, variant_count(), len(samples))
    Xp = np.empty_like(Xs)
    for j, k in enumerate(order):
        Xp[2 * k:2 * k + 2] = Xs[2 * j:2 * j + 2]
    Lp = np.argmax(m.predict_proba(Xp), -1)
    # Unphased genotype agreement at model SNPs present in both (phase differs between callsets).
    both = np.zeros(m.C, dtype=bool)
    both[match["model_idx"]] = True
    gd, gp = Xd[0::2] + Xd[1::2], Xp[0::2] + Xp[1::2]
    ok = both & (Xd[0::2] < 2).all(axis=0) & (Xd[1::2] < 2).all(axis=0)
    rows = []
    for k, iid in enumerate(ids):
        dip = lambda L: np.sort(L[2 * k:2 * k + 2], axis=0)  # noqa: E731
        main_off = pops[np.bincount(Ld[2 * k:2 * k + 2].ravel(), minlength=len(pops)).argmax()]
        main_ours = pops[np.bincount(Lp[2 * k:2 * k + 2].ravel(), minlength=len(pops)).argmax()]
        rows.append({"sample": iid, "population": samples.loc[iid, "Population"],
                     "official_main_call": main_off, "our_main_call": main_ours,
                     "genotype_agreement": round(float((gd[k][ok] == gp[k][ok]).mean()), 4),
                     "diploid_window_agreement": round(float((dip(Ld) == dip(Lp)).all(axis=0).mean()), 4)})
    conc = pd.DataFrame(rows)
    log(conc.to_string(index=False))
    if (conc.official_main_call != conc.our_main_call).any() or conc.diploid_window_agreement.mean() < 0.9:
        raise RuntimeError("our model input disagrees with the official Gnomix path on the demo samples")
    with atomic_path(out / "official_path_concordance.tsv") as tmp:  # written only on success: skip marker
        conc.to_csv(tmp, sep="\t", index=False)
