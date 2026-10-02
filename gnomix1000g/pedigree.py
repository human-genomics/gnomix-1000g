"""Which samples have a pedigree-exact phase in the 1000 Genomes panel?

For every parent-child pair in the panel we read both phased genotypes at the
model SNPs (TRIO_CHROMS) and count phase switches that no meiosis explains:

  child   at sites where the child is heterozygous and the parent homozygous,
          which child haplotype carries the parent's allele. An exact phase
          never changes it.
  parent  at sites where the parent is heterozygous and the child homozygous,
          which parent haplotype carries the transmitted allele. An exact phase
          changes it only at the crossovers of that meiosis (about one per
          100 cM).

A switch counts only if the new state holds for at least MIN_RUN informative
sites, so isolated genotype errors are ignored.

Result in the published run: trio children have no such switches, while their
parents (and duo members) have about as many as statistically phased samples.
Only trio children keep their published phase; Gnofix decides for everyone
else (tracts.final_uses_gnofix).

Outputs (output/validation/): pedigree_phase_pairs.tsv, pedigree_phase_summary.tsv
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .common import OUTPUT, atomic_path, log
from .haplotypes import load_match, read_haplotypes
from .prepare import model_meta, variant_count, write_samples

MIN_RUN = 20  # informative sites a new phase state must hold to count as a switch


def switches(state: np.ndarray, min_run: int = MIN_RUN) -> int:
    """Changes of a 0/1 state sequence, after dropping runs shorter than min_run."""
    if len(state) < 2:
        return 0
    starts = np.r_[0, np.flatnonzero(np.diff(state)) + 1]
    lengths = np.diff(np.r_[starts, len(state)])
    kept = state[starts][lengths >= min_run]
    return int((np.diff(kept) != 0).sum())


def pair_switches(child: np.ndarray, parent: np.ndarray) -> tuple[int, int, int, int]:
    """(child switches, child informative sites, parent switches, parent informative sites).

    child and parent are (2, n_snps) haplotype pairs.
    """
    c_het, p_het = child[0] != child[1], parent[0] != parent[1]
    m = c_het & ~p_het
    child_state = (child[1, m] == parent[0, m]).astype(np.int8)
    m2 = p_het & ~c_het
    parent_state = (parent[1, m2] == child[0, m2]).astype(np.int8)
    return switches(child_state), int(m.sum()), switches(parent_state), int(m2.sum())


def family_pairs(samples: pd.DataFrame) -> pd.DataFrame:
    """One row per (child, parent) with both in the panel; kind = trio or duo (parents of the child in the panel)."""
    ids = set(samples.IID)
    rows = []
    for r in samples.itertuples():
        parents = [p for p in (r.PAT, r.MAT) if p in ids]
        for p in parents:
            rows.append({"child": r.IID, "parent": p, "Population": r.Population,
                         "kind": "trio" if len(parents) == 2 else "duo"})
    return pd.DataFrame(rows)


def run(chroms: list[int]) -> pd.DataFrame:
    vdir = OUTPUT / "validation"
    vdir.mkdir(parents=True, exist_ok=True)
    samples = write_samples()
    pairs = family_pairs(samples)
    members = samples[samples.IID.isin(set(pairs.child) | set(pairs.parent))]
    row_of = {iid: k for k, iid in enumerate(members.IID)}
    out = []
    expected = 0.0
    for c in chroms:
        meta = model_meta(c)
        cm = np.interp(meta["snp_pos"], meta["gm_pos"], meta["gm_cm"])
        expected += (cm[-1] - cm[0]) / 100
        X = read_haplotypes(c, members.pgen_index.values, int(meta["C"]), load_match(c), variant_count(),
                            len(samples))
        for r in pairs.itertuples():
            i, j = row_of[r.child], row_of[r.parent]
            cs, cn, ps, pn = pair_switches(X[2 * i:2 * i + 2], X[2 * j:2 * j + 2])
            out.append({"child": r.child, "parent": r.parent, "chrom": c, "child_switches": cs,
                        "child_sites": cn, "parent_switches": ps, "parent_sites": pn})
        del X
        log(f"    chr{c}: {len(pairs)} parent-child pairs checked")
    per = pd.DataFrame(out).groupby(["child", "parent"], as_index=False)[
        ["child_switches", "child_sites", "parent_switches", "parent_sites"]].sum()
    per = pairs.merge(per, on=["child", "parent"])
    per["expected_crossovers"] = round(expected, 2)
    summary = per.groupby("kind").agg(
        pairs=("child", "size"),
        child_switches_median=("child_switches", "median"),
        child_pairs_without_switches=("child_switches", lambda s: float((s == 0).mean())),
        parent_switches_median=("parent_switches", "median"),
        parent_switches_min=("parent_switches", "min"),
        expected_crossovers=("expected_crossovers", "first")).reset_index()
    summary["chromosomes"] = ",".join(map(str, chroms))
    with atomic_path(vdir / "pedigree_phase_pairs.tsv") as tmp:
        per.to_csv(tmp, sep="\t", index=False)
    with atomic_path(vdir / "pedigree_phase_summary.tsv") as tmp:
        summary.to_csv(tmp, sep="\t", index=False, float_format="%.4g")
    log(f"  pedigree phase, chromosomes {chroms} (switches not explained by a meiosis):")
    log(summary.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
    return summary
