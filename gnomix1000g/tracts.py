"""Turn window-level calls into msp files, ancestry tracts and global ancestry.

Three call sets exist for every haplotype window (see infer.py):
  raw     the pretrained model on the published 1000 Genomes phase
  gnofix  after Gnofix phase correction
  final   per sample, gnofix or raw according to output/validation/gnofix_policy.json.
          Trio children always keep the raw call: the panel phased them with
          both parents. Their parents and duo members get the population's
          decision like everyone else (validation/pedigree_phase_*.tsv).

Coordinates: the model works on GRCh37 SNP positions, so windows are defined in
GRCh37. Each window also gets GRCh38 coordinates from its lifted SNPs; windows
that liftover reverses, stretches or moves out of order have none
(hg38_discordant). A tract's GRCh38 span runs from its first to its last
window with GRCh38 coordinates (-1 if it has none). All positions are 1-based,
inclusive SNP positions. The msp "n snps" column counts the model SNPs found
in the panel (upstream Gnomix counts query VCF SNPs).

Outputs (output/):
  windows/chrN.windows.tsv                 window table: GRCh37 + GRCh38 + cM
  msp/chrN.{final,raw,gnofix}.msp.gz       Gnomix .msp format (GRCh37 spos/epos, as Gnomix writes)
  posteriors/chrN.final.npz                posterior per haplotype, window, ancestry (uint8, /255)
  tracts/tracts_{final,raw,gnofix}.tsv.gz  all haplotype tracts, both builds
  tracts/bed_hg38/IID_{A,B}.bed            Martin et al. (2017) format, UNK below posterior 0.9
  tracts/gnofix_switches.tsv.gz            where Gnofix exchanged the two haplotypes
  global_ancestry_{final,raw,gnofix}.tsv   per-sample ancestry fractions (cM-weighted)
  population_ancestry.tsv                  population means of the final calls
"""

from __future__ import annotations

import bisect
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from .common import OUTPUT, ANCESTRIES, UNK_THRESHOLD, WORK, atomic_path, group_of, log, read_json
from .infer import window_of_snp
from .prepare import model_meta, write_samples

CALL_SETS = ("final", "raw", "gnofix")
STRETCH_SLACK = 1_000_000  # bp: a lifted window may be at most 2x + this longer than in GRCh37


def sample_group(samples: pd.DataFrame) -> list[str]:
    return [group_of(p, sp) for p, sp in zip(samples.Population, samples.SuperPop)]


def pedigree_phased(samples: pd.DataFrame) -> np.ndarray:
    """Samples whose published phase is pedigree-exact: trio children (both parents in the panel).

    Parents and duo members are not: their phase switches as often as a statistical phase
    (pedigree.py), so Gnofix can help them as it helps unrelated samples.
    """
    return samples.trio_child.values.astype(bool)


def policy_applies(policy: dict, samples: pd.DataFrame) -> np.ndarray:
    """Per-population decision where the benchmark made one, otherwise the group decision."""
    per_pop = policy.get("apply_gnofix_population", {})
    return np.array([bool(per_pop[p]) if p in per_pop else bool(policy["apply_gnofix"].get(g, False))
                     for p, g in zip(samples.Population, sample_group(samples))])


def final_uses_gnofix(samples: pd.DataFrame) -> np.ndarray:
    """Per sample: do the final calls (and the corrected panel) use Gnofix?"""
    path = OUTPUT / "validation" / "gnofix_policy.json"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: run the 'trio' step first")
    return policy_applies(read_json(path), samples) & ~pedigree_phased(samples)


def calls(chrom: int, use_fix: np.ndarray) -> dict:
    """Labels (2N, W) and posteriors (2N, W, A) of the three call sets of a chromosome."""
    with np.load(WORK / "infer" / f"chr{chrom}.npz") as z:
        d = {k: z[k] for k in z.files}
    rows = np.repeat(use_fix, 2)
    return {"raw": (d["lab_raw"], d["prob_raw"]), "gnofix": (d["lab_fix"], d["prob_fix"]),
            "final": (np.where(rows[:, None], d["lab_fix"], d["lab_raw"]),
                      np.where(rows[:, None, None], d["prob_fix"], d["prob_raw"])),
            "swap": d["swap"]}


def called_posterior(labels: np.ndarray, probs: np.ndarray) -> np.ndarray:
    """Posterior (0-1) of the called ancestry in every haplotype window."""
    return np.take_along_axis(probs, labels[..., None].astype(np.int64), axis=2)[..., 0] / 255


# ---- windows -------------------------------------------------------------------

def longest_increasing(x: np.ndarray) -> np.ndarray:
    """Indices of a longest strictly increasing subsequence of x."""
    tails, tail_idx, prev = [], [], np.full(len(x), -1)
    for i, v in enumerate(x):
        k = bisect.bisect_left(tails, v)
        if k == len(tails):
            tails.append(v)
            tail_idx.append(i)
        else:
            tails[k], tail_idx[k] = v, i
        prev[i] = tail_idx[k - 1] if k else -1
    out, i = [], tail_idx[-1] if tail_idx else -1
    while i >= 0:
        out.append(i)
        i = prev[i]
    return np.array(out[::-1], dtype=np.int64)


def hg38_discordant(first: np.ndarray, last: np.ndarray, span19: np.ndarray) -> np.ndarray:
    """Windows without usable GRCh38 coordinates: no lifted SNP, reversed or stretched
    (> 2x + STRETCH_SLACK) by liftover, or out of order with the other windows (outside
    the longest increasing run of starts, or overlapping the previous window)."""
    s, e = np.minimum(first, last), np.maximum(first, last)
    bad = (first < 0) | (first > last) | (e - s > 2 * span19 + STRETCH_SLACK)
    ok = np.nonzero(~bad)[0]
    keep = np.zeros(len(first), dtype=bool)
    keep[ok[longest_increasing(s[ok])]] = True
    prev_end = -1
    for w in np.nonzero(keep)[0]:
        if s[w] <= prev_end:
            keep[w] = False
        else:
            prev_end = e[w]
    return ~keep


def window_table(chrom: int) -> pd.DataFrame:
    meta = model_meta(chrom)
    with np.load(WORK / "match" / f"chr{chrom}.npz") as z:
        pos38, model_idx = z["pos38"], z["model_idx"]
    C, M, W = int(meta["C"]), int(meta["M"]), int(meta["W"])
    pos19 = meta["snp_pos"]
    s_idx = np.arange(W) * M
    e_idx = np.append(s_idx[1:], C) - 1
    win = window_of_snp(C, M, W)
    in_panel = np.zeros(C, dtype=bool)
    in_panel[model_idx] = True

    df = pd.DataFrame({"chrom": chrom, "window": np.arange(W),
                       "spos_hg19": pos19[s_idx], "epos_hg19": pos19[e_idx]})
    # Gnomix: linear interpolation of the model's genetic map, clamped at the ends
    df["sgpos"] = np.round(np.interp(df.spos_hg19.values, meta["gm_pos"], meta["gm_cm"]), 5)
    df["egpos"] = np.round(np.interp(df.epos_hg19.values, meta["gm_pos"], meta["gm_cm"]), 5)
    df["n_model_snps"] = np.bincount(win, minlength=W)
    df["n_query_snps"] = np.bincount(win[in_panel], minlength=W)
    # GRCh38: first and last lifted SNP of each window, in model order.
    li = np.nonzero(pos38 > 0)[0]
    lw = win[li]
    starts = np.searchsorted(lw, np.arange(W), side="left")
    ends = np.searchsorted(lw, np.arange(W), side="right") - 1
    has = ends >= starts
    first = np.full(W, -1, dtype=np.int64)
    last = np.full(W, -1, dtype=np.int64)
    first[has] = pos38[li[starts[has]]]
    last[has] = pos38[li[ends[has]]]
    df["spos_hg38"] = np.minimum(first, last)
    df["epos_hg38"] = np.maximum(first, last)
    bad = hg38_discordant(first, last, (df.epos_hg19 - df.spos_hg19).values)
    df["hg38_discordant"] = bad
    return df


def window_lengths_cm(win: pd.DataFrame) -> np.ndarray:
    """Genetic length covered by each window: from its start to the next window's start."""
    s = win.sgpos.values
    return np.maximum(np.append(s[1:], win.egpos.values[-1]) - s, 0)


# ---- tracts ----------------------------------------------------------------------

def label_runs(labels: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Runs of equal labels along each row: (row, first window, last window)."""
    H, W = labels.shape
    change = np.ones((H, W), dtype=bool)
    change[:, 1:] = labels[:, 1:] != labels[:, :-1]
    rows, start = np.nonzero(change)  # row-major: sorted by row, then window
    last_in_row = np.append(rows[1:] != rows[:-1], True)
    end = np.where(last_in_row, W - 1, np.append(start[1:], 0) - 1)
    return rows, start, end


def hg38_span(win: pd.DataFrame, start_w: np.ndarray, end_w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """GRCh38 start/end of window runs: first and last concordant window inside each run (-1: none)."""
    good = ~win.hg38_discordant.values
    W = len(good)
    idx = np.arange(W)
    next_good = np.minimum.accumulate(np.where(good, idx, W)[::-1])[::-1]
    prev_good = np.maximum.accumulate(np.where(good, idx, -1))
    a, b = next_good[start_w], prev_good[end_w]
    ok = (a < W) & (b >= 0) & (a <= b)
    s38 = np.where(ok, win.spos_hg38.values[np.clip(a, 0, W - 1)], -1)
    e38 = np.where(ok, win.epos_hg38.values[np.clip(b, 0, W - 1)], -1)
    inverted = e38 < s38  # cannot happen after window_table's order check; keep the output valid
    return np.where(inverted, -1, s38), np.where(inverted, -1, e38)


def haplotype_tracts(labels: np.ndarray, probs: np.ndarray, win: pd.DataFrame, iids, chrom: int,
                     names=None) -> pd.DataFrame:
    """Collapse runs of equal labels into tracts, for every haplotype.

    `probs` gives each window's posterior of its label (0-1); `names` maps label
    codes to ancestry names (default: the model populations).
    """
    names = np.asarray(ANCESTRIES if names is None else names)
    rows, s_w, e_w = label_runs(labels)
    csum = np.concatenate([np.zeros((labels.shape[0], 1)), np.cumsum(probs, axis=1)], axis=1)
    mean_p = (csum[rows, e_w + 1] - csum[rows, s_w]) / (e_w - s_w + 1)
    start38, end38 = hg38_span(win, s_w, e_w)
    return pd.DataFrame({
        "sample": np.asarray(iids)[rows // 2],
        "haplotype": np.where(rows % 2 == 0, "A", "B"),
        "chrom": chrom,
        "start_hg38": start38, "end_hg38": end38,
        "start_hg19": win.spos_hg19.values[s_w], "end_hg19": win.epos_hg19.values[e_w],
        "start_cM": win.sgpos.values[s_w], "end_cM": win.egpos.values[e_w],
        "ancestry": names[labels[rows, s_w]],
        "n_windows": e_w - s_w + 1,
        "mean_posterior": np.round(mean_p, 4),
    })


def write_msp(path: Path, chrom: int, win: pd.DataFrame, labels: np.ndarray, iids) -> None:
    header = "#Subpopulation order/codes: " + "\t".join(f"{p}={i}" for i, p in enumerate(ANCESTRIES))
    cols = ["#chm", "spos", "epos", "sgpos", "egpos", "n snps"] + [f"{s}.{h}" for s in iids for h in (0, 1)]
    meta = [f"{chrom}\t{a}\t{b}\t{c:.5f}\t{d:.5f}\t{n}" for a, b, c, d, n in
            zip(win.spos_hg19, win.epos_hg19, win.sgpos, win.egpos, win.n_query_snps)]
    with atomic_path(path) as tmp, gzip.open(tmp, "wt", compresslevel=6) as fh:
        fh.write(header + "\n" + "\t".join(cols) + "\n")
        for w, m in enumerate(meta):
            fh.write(m + "\t" + "\t".join(map(str, labels[:, w])) + "\n")


def write_martin_beds(df: pd.DataFrame) -> None:
    """One BED per haplotype (GRCh38): chrom, spos, epos, ancestry, sgpos, egpos."""
    out = OUTPUT / "tracts" / "bed_hg38"
    out.mkdir(parents=True, exist_ok=True)
    df = df[df.start_hg38 > 0].rename(columns={"start_hg38": "spos", "end_hg38": "epos",
                                               "start_cM": "sgpos", "end_cM": "egpos"})
    for (sample, hap), g in df.groupby(["sample", "haplotype"], sort=False):
        with atomic_path(out / f"{sample}_{hap}.bed") as tmp:
            g.sort_values(["chrom", "spos"])[["chrom", "spos", "epos", "ancestry", "sgpos", "egpos"]].to_csv(
                tmp, sep="\t", index=False, header=False)


def write_table(df: pd.DataFrame, path: Path, **kw) -> None:
    with atomic_path(path) as tmp:
        df.to_csv(tmp, sep="\t", **kw)


# ---- main ----------------------------------------------------------------------

def run(chroms: list[int]) -> None:
    samples = write_samples()
    iids = samples.IID.values
    use_fix = final_uses_gnofix(samples)
    log(f"  final calls use Gnofix for {use_fix.sum():,} of {len(samples):,} samples")
    A = len(ANCESTRIES)
    totals = {k: np.zeros((len(samples) * 2, A)) for k in CALL_SETS}
    tracts = {k: [] for k in CALL_SETS}
    unk_tracts, switch_rows = [], []
    for c in chroms:
        win = window_table(c)
        write_table(win, OUTPUT / "windows" / f"chr{c}.windows.tsv", index=False)
        cs = calls(c, use_fix)
        wlen = window_lengths_cm(win)
        for k in CALL_SETS:
            lab, prob = cs[k]
            for a in range(A):
                totals[k][:, a] += ((lab == a) * wlen).sum(axis=1)
            tracts[k].append(haplotype_tracts(lab, called_posterior(lab, prob), win, iids, c))
            write_msp(OUTPUT / "msp" / f"chr{c}.{k}.msp.gz", c, win, lab, iids)
        lab, prob = cs["final"]
        with atomic_path(OUTPUT / "posteriors" / f"chr{c}.final.npz") as tmp:
            np.savez_compressed(tmp, prob=prob, samples=iids, populations=np.array(ANCESTRIES))
        # Martin et al. style: windows whose called-ancestry posterior is below 0.9 are UNK.
        post = called_posterior(lab, prob)
        lab_unk = np.where(post >= UNK_THRESHOLD, lab, A)
        unk_tracts.append(haplotype_tracts(lab_unk, post, win, iids, c, names=ANCESTRIES + ["UNK"]))
        # Gnofix switch points: window boundaries where a sample's swap state changes.
        swap = cs["swap"].astype(np.int8)
        k_idx, w_idx = np.nonzero(np.diff(swap, axis=1, prepend=0) != 0)
        switch_rows.append(pd.DataFrame({
            "sample": iids[k_idx], "chrom": c, "window": w_idx,
            "pos_hg19": win.spos_hg19.values[w_idx], "pos_hg38": win.spos_hg38.values[w_idx],
            "cM": win.sgpos.values[w_idx], "used_in_final": use_fix[k_idx]}))
        log(f"  chr{c}: tracts and msp written")

    for k in CALL_SETS:
        write_table(pd.concat(tracts[k], ignore_index=True), OUTPUT / "tracts" / f"tracts_{k}.tsv.gz",
                    index=False, compression="gzip")
    write_martin_beds(pd.concat(unk_tracts, ignore_index=True))
    write_table(pd.concat(switch_rows, ignore_index=True), OUTPUT / "tracts" / "gnofix_switches.tsv.gz",
                index=False, compression="gzip")

    meta_cols = ["IID", "SuperPop", "Population", "SEX", "trio_child", "has_parent_in_panel",
                 "has_child_in_panel", "phase3_unrelated_2504"]
    for k in CALL_SETS:
        t = totals[k]
        dip = t[0::2] + t[1::2]
        out = samples[meta_cols].copy()
        out["group"] = sample_group(samples)
        if k == "final":
            out["gnofix_in_final"] = use_fix
        for a, p in enumerate(ANCESTRIES):
            out[p] = np.round(dip[:, a] / dip.sum(axis=1), 5)
        for h, name in ((0, "A"), (1, "B")):
            th = t[h::2]
            for a, p in enumerate(ANCESTRIES):
                out[f"{p}_hap{name}"] = np.round(th[:, a] / th.sum(axis=1), 5)
        out["total_cM"] = np.round(dip.sum(axis=1) / 2, 3)
        write_table(out, OUTPUT / f"global_ancestry_{k}.tsv", index=False)
        if k == "final":
            pop = out.groupby(["SuperPop", "Population"])[ANCESTRIES].mean().round(4)
            pop.insert(0, "n", out.groupby(["SuperPop", "Population"]).size())
            write_table(pop, OUTPUT / "population_ancestry.tsv")
    log("  global ancestry written")
