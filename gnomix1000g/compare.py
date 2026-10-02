"""Compare our calls with Martin et al. (2017) and the 1000 Genomes Phase 1 LAI.

1. Martin et al. (2017), AJHG: population mean ancestry of the six admixed
   American populations (K=3), admixture timing, and individually named samples.
   The paper's per-sample tracts are not public, so we compare what it reports.
2. 1000 Genomes Phase 1 admixture working group (1000 Genomes Consortium
   2012): per-sample global proportions and diploid local-ancestry tracts
   (GRCh37) for 242 ASW, CLM, MXL and PUR samples. Consensus of LAMP-LD,
   HAPMIX, RFMix and MULTIMIX.

Both earlier analyses use three sources (AFR, EUR, NAT). Our model has eight.
For comparison we map AFR + AHG → AFR, EUR + WAS → EUR, NAT → NAT, and keep
EAS + SAS + OCE as OTHER (reported, and dropped when we renormalise to three ancestries).

Outputs: output/validation/compare_*.tsv, compare_*.png, COMPARISON.md
"""

from __future__ import annotations

import io
import zipfile

import numpy as np
import pandas as pd

from .common import ADMIXED_POPULATIONS as ADMIXED
from .common import ANCESTRIES, ANCESTRY_COLORS, DOWNLOADS, OUTPUT, REFERENCE, atomic_path, log
from .prepare import write_samples

THREE = {"AFR": ["AFR", "AHG"], "EUR": ["EUR", "WAS"], "NAT": ["NAT"], "OTHER": ["EAS", "SAS", "OCE"]}
P1_CODES = {1: ("EUR", "EUR"), 2: ("EUR", "AFR"), 3: ("AFR", "AFR"), 4: ("EUR", "NAT"),
            5: ("AFR", "NAT"), 6: ("NAT", "NAT")}


def three_way(df: pd.DataFrame, renormalise: bool) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    for k, cols in THREE.items():
        out[k] = df[cols].sum(axis=1)
    if renormalise:
        s = out[["AFR", "EUR", "NAT"]].sum(axis=1)
        out = out[["AFR", "EUR", "NAT"]].div(s, axis=0)
    return out


def read_ref(name: str) -> pd.DataFrame:
    return pd.read_csv(REFERENCE / name, sep="\t", comment="#")


# ---- Martin et al. 2017 ------------------------------------------------------------

def martin_population(glob: pd.DataFrame) -> pd.DataFrame:
    ref = read_ref("martin2017_population_ancestry.tsv")
    g = glob[glob.phase3_unrelated_2504 & glob.Population.isin(ADMIXED)]
    rows = []
    for _, r in ref.iterrows():
        d = g[g.Population == r.population]
        raw = three_way(d[ANCESTRIES], renormalise=False)
        ren = three_way(d[ANCESTRIES], renormalise=True)
        m = ren[r.ancestry]
        se = m.std(ddof=1) / np.sqrt(len(m))
        rows.append({"population": r.population, "ancestry": r.ancestry, "n_ours": len(d),
                     "martin_mean": r["mean"], "martin_ci": f"{r.ci_low:.2f}-{r.ci_high:.2f}",
                     "ours_mean_3way": round(m.mean(), 3),
                     "ours_ci_3way": f"{m.mean() - 1.96 * se:.3f}-{m.mean() + 1.96 * se:.3f}",
                     "ours_mean_8way_collapsed": round(raw[r.ancestry].mean(), 3),
                     "ours_other_EAS_SAS_OCE": round(raw["OTHER"].mean(), 3),
                     "difference": round(m.mean() - r["mean"], 3)})
    return pd.DataFrame(rows)


def martin_named(glob: pd.DataFrame) -> pd.DataFrame:
    ref = read_ref("martin2017_named_samples.tsv")
    g = glob.set_index("IID")
    rows = []
    for _, r in ref.iterrows():
        row = {"sample": r["sample"], "population": r.population, "martin_statement": r.statement,
               "expectation": r.expectation}
        if r["sample"] in g.index:
            for p in ANCESTRIES:
                row[p] = round(float(g.loc[r["sample"], p]), 3)
        rows.append(row)
    return pd.DataFrame(rows)


def admixture_timing(tracts: pd.DataFrame, glob: pd.DataFrame, min_cm: float = 2.0) -> pd.DataFrame:
    """Rough single-pulse admixture time from tract lengths.

    Tracts are first merged to continental labels (AFR + AHG → AFR, EUR + WAS → EUR),
    because EUR-WAS label changes would otherwise split them. Under a single
    pulse T generations ago, tracts of ancestry a are exponential with mean
    1 / (T (1 - m_a)) Morgans; tracts shorter than min_cm are dropped and the
    estimate uses the mean excess length over min_cm. Residual phase and label
    errors still split tracts, so this overestimates T. A sanity check only;
    Martin et al. fitted pulse models with Tracts.
    """
    ref = read_ref("martin2017_admixture_timing.tsv")
    g = glob[glob.phase3_unrelated_2504].set_index("IID")
    cont = {a: k for k, v in THREE.items() for a in v}
    t = tracts[tracts["sample"].isin(g.index[g.Population.isin(ADMIXED)])].copy()
    t["c"] = t.ancestry.map(cont)
    t = t.sort_values(["sample", "haplotype", "chrom", "start_cM"])
    key = t["sample"] + t.haplotype + t.chrom.astype(str)
    t["run"] = ((t.c != t.c.shift()) | (key != key.shift())).cumsum()
    runs = t.groupby("run").agg(sample=("sample", "first"), c=("c", "first"),
                                s=("start_cM", "min"), e=("end_cM", "max"))
    runs["L"] = (runs.e - runs.s) / 100
    runs["pop"] = runs["sample"].map(g.Population)
    rows = []
    for pop in ADMIXED:
        ids = g.index[g.Population == pop]
        props = three_way(g.loc[ids, ANCESTRIES], renormalise=False).mean()
        for a in ("AFR", "EUR", "NAT"):
            m = props[a]
            if m < 0.05 or m > 0.95:
                continue
            lens = runs.L[(runs["pop"] == pop) & (runs.c == a)]
            lens = lens[lens > min_cm / 100]
            if len(lens) < 50:
                continue
            T = 1 / ((lens.mean() - min_cm / 100) * (1 - m))
            rows.append({"population": pop, "ancestry": a, "proportion": round(m, 3), "tracts_over_2cM": len(lens),
                         "mean_tract_cM": round(100 * lens.mean(), 1), "single_pulse_T": round(T, 1)})
    out = pd.DataFrame(rows, columns=["population", "ancestry", "proportion", "tracts_over_2cM",
                                      "mean_tract_cM", "single_pulse_T"])
    out["martin_tracts_estimate"] = out.population.map(dict(zip(ref.population, ref.generations)))
    return out


# ---- 1000 Genomes Phase 1 -------------------------------------------------------------

def phase1_global(glob: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pop in ("ASW", "CLM", "MXL", "PUR"):
        p = pd.read_csv(DOWNLOADS / f"{pop}_phase1_globalproportions_and_unknown.txt", sep=r"\s+",
                        header=None, names=["sample", "EUR", "AFR", "NAT", "UNK"])
        p["population"] = pop
        rows.append(p)
    p1 = pd.concat(rows)
    s = p1[["EUR", "AFR", "NAT"]].sum(axis=1)
    for a in ("EUR", "AFR", "NAT"):
        p1[f"p1_{a}"] = p1[a] / s
    g = glob.set_index("IID")
    p1 = p1[p1["sample"].isin(g.index)].copy()
    ours = three_way(g.loc[p1["sample"], ANCESTRIES], renormalise=True)
    raw = three_way(g.loc[p1["sample"], ANCESTRIES], renormalise=False)
    for a in ("EUR", "AFR", "NAT"):
        p1[f"ours_{a}"] = ours[a].values
    p1["ours_OTHER"] = raw["OTHER"].values
    p1["p1_UNK"] = p1["UNK"]
    return p1[["sample", "population"] + [f"p1_{a}" for a in ("EUR", "AFR", "NAT", "UNK")]
              + [f"ours_{a}" for a in ("EUR", "AFR", "NAT", "OTHER")]]


def phase1_global_summary(p1: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for pop, group in list(p1.groupby("population")) + [("all", p1)]:
        for a in ("EUR", "AFR", "NAT"):
            # Phase 1 ran ASW with EUR and AFR only (NAT set to 0): leave ASW out of NAT.
            d = group[group.population != "ASW"] if a == "NAT" else group
            if not len(d):
                continue
            x, y = d[f"p1_{a}"], d[f"ours_{a}"]
            rows.append({"population": pop, "ancestry": a, "n": len(d),
                         "pearson_r": round(float(np.corrcoef(x, y)[0, 1]), 4) if x.std() > 0 else np.nan,
                         "mean_phase1": round(x.mean(), 4), "mean_ours": round(y.mean(), 4),
                         "mean_abs_diff": round(float((x - y).abs().mean()), 4)})
    return pd.DataFrame(rows)


def phase1_tracts() -> pd.DataFrame:
    rows = []
    for pop in ("ASW", "CLM", "MXL", "PUR"):
        with zipfile.ZipFile(DOWNLOADS / f"{pop}_phase1_ancestry_deconvolution.zip") as z:
            for name in z.namelist():
                if not name.endswith(".bed"):
                    continue
                d = pd.read_csv(io.BytesIO(z.read(name)), sep="\t", header=None,
                                names=["chrom", "start", "end", "code", "length"], dtype={"code": str})
                d["sample"] = name.split("/")[-1][:-4]
                d["population"] = pop
                rows.append(d)
    t = pd.concat(rows, ignore_index=True)
    t = t[t.code.isin([str(k) for k in P1_CODES])].copy()
    t["code"] = t.code.astype(int)
    t["chrom"] = t.chrom.astype(str).str.replace("chr", "").astype(int)
    return t


def phase1_tract_concordance(chroms: list[int]) -> pd.DataFrame:
    """Diploid concordance per window (GRCh37 midpoints), weighted by window cM."""
    from .tracts import calls, final_uses_gnofix, window_lengths_cm

    samples = write_samples()
    use_fix = final_uses_gnofix(samples)
    t = phase1_tracts()
    idx_of = {s: i for i, s in enumerate(samples.IID)}
    pops = np.array(ANCESTRIES)
    collapse = {a: k for k, v in THREE.items() for a in v}
    to3 = np.array([collapse[p] for p in pops])
    res = []
    for c in chroms:
        win = pd.read_csv(OUTPUT / "windows" / f"chr{c}.windows.tsv", sep="\t")
        mid = ((win.spos_hg19 + win.epos_hg19) // 2).values
        wlen = window_lengths_cm(win)
        lab, _ = calls(c, use_fix)["final"]
        tc = t[t.chrom == c]
        for (sample, pop), d in tc.groupby(["sample", "population"]):
            if sample not in idx_of:
                continue
            k = idx_of[sample]
            a, b = to3[lab[2 * k]], to3[lab[2 * k + 1]]
            d = d.sort_values("start")
            # Phase 1 BED: 0-based start, 1-based end; window midpoints are 1-based
            j = np.searchsorted(d.start.values, mid, side="left") - 1
            inside = (j >= 0) & (mid <= d.end.values[np.clip(j, 0, None)])
            codes = np.where(inside, d.code.values[np.clip(j, 0, None)], 0)
            called = codes > 0
            exp = np.array([P1_CODES.get(x, ("", "")) for x in codes])
            agree = ((a == exp[:, 0]) & (b == exp[:, 1])) | ((a == exp[:, 1]) & (b == exp[:, 0]))
            res.append({"sample": sample, "population": pop, "chrom": c,
                        "cM_compared": wlen[called].sum(), "cM_agree": wlen[called & agree].sum(),
                        "cM_ours_other": wlen[called & ((a == "OTHER") | (b == "OTHER"))].sum()})
    r = pd.DataFrame(res)
    per = r.groupby(["sample", "population"])[["cM_compared", "cM_agree", "cM_ours_other"]].sum().reset_index()
    per["concordance"] = per.cM_agree / per.cM_compared
    return per


# ---- figure and report ------------------------------------------------------------------

def scatter(p1: pd.DataFrame, path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"ASW": "#2a78d6", "CLM": "#eb6834", "MXL": "#1baf7a", "PUR": "#4a3aa7"}  # fixed per population
    fig, axes = plt.subplots(1, 3, figsize=(11, 4), dpi=130)
    for ax, a in zip(axes, ("AFR", "EUR", "NAT")):
        d = p1 if a != "NAT" else p1[p1.population != "ASW"]
        ax.plot([0, 1], [0, 1], color="#999999", lw=0.8, zorder=1)
        for pop, g in d.groupby("population"):
            ax.scatter(g[f"p1_{a}"], g[f"ours_{a}"], s=12, color=colors[pop], label=pop, zorder=2, linewidths=0)
        r = np.corrcoef(d[f"p1_{a}"], d[f"ours_{a}"])[0, 1]
        ax.set_title(f"{a}: r = {r:.4f}, n = {len(d)}", fontsize=11, loc="left")
        ax.set_xlabel("Phase 1 consensus LAI", fontsize=10)
        ax.set_ylabel("This pipeline", fontsize=10)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect("equal")
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.tick_params(labelsize=9, colors="#444444")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(path)
    plt.close(fig)


def bars(glob: pd.DataFrame, path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    g = glob[glob.Population.isin(ADMIXED) & glob.phase3_unrelated_2504]
    fig, axes = plt.subplots(1, len(ADMIXED), figsize=(14, 3.2), dpi=120, sharey=True)
    for ax, pop in zip(axes, ["ACB", "ASW", "PUR", "CLM", "MXL", "PEL"]):
        d = g[g.Population == pop].sort_values(["AFR", "NAT"])
        bottom = np.zeros(len(d))
        for p in ["AFR", "AHG", "EUR", "WAS", "NAT", "EAS", "SAS", "OCE"]:
            ax.bar(np.arange(len(d)), d[p], bottom=bottom, width=1.0, color=ANCESTRY_COLORS[p], label=p, linewidth=0)
            bottom += d[p].values
        ax.set_title(f"{pop} (n={len(d)})", fontsize=10)
        ax.set_xticks([])
        ax.set_xlim(-0.5, len(d) - 0.5)
    axes[0].set_ylabel("Global ancestry")
    axes[-1].legend(loc="center left", bbox_to_anchor=(1, 0.5), frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _write(df: pd.DataFrame, path) -> None:
    with atomic_path(path) as tmp:
        df.to_csv(tmp, sep="\t", index=False)


def run(chroms: list[int] | None = None) -> None:
    from .common import chrom_list

    chroms = chroms or chrom_list()
    vdir = OUTPUT / "validation"
    vdir.mkdir(parents=True, exist_ok=True)
    glob = pd.read_csv(OUTPUT / "global_ancestry_final.tsv", sep="\t")
    tracts = pd.read_csv(OUTPUT / "tracts" / "tracts_final.tsv.gz", sep="\t")

    mp = martin_population(glob)
    mn = martin_named(glob)
    tim = admixture_timing(tracts, glob)
    p1 = phase1_global(glob)
    p1s = phase1_global_summary(p1)
    conc = phase1_tract_concordance(chroms)
    _write(mp, vdir / "compare_martin2017_population.tsv")
    _write(mn, vdir / "compare_martin2017_named_samples.tsv")
    _write(tim, vdir / "compare_martin2017_timing.tsv")
    _write(p1, vdir / "compare_phase1_global_per_sample.tsv")
    _write(p1s, vdir / "compare_phase1_global_summary.tsv")
    _write(conc, vdir / "compare_phase1_tract_concordance.tsv")
    scatter(p1, vdir / "compare_phase1_global.png")
    bars(glob, vdir / "admixed_global_ancestry.png")

    cs = conc.groupby("population").apply(lambda d: pd.Series({
        "samples": len(d), "concordance_cM_weighted": d.cM_agree.sum() / d.cM_compared.sum(),
        "median_sample_concordance": d.concordance.median(),
        "cM_with_our_EAS_SAS_OCE_call": d.cM_ours_other.sum() / d.cM_compared.sum()})).reset_index()
    _write(cs, vdir / "compare_phase1_tract_concordance_summary.tsv")
    md = ["# Comparison with Martin et al. (2017) and 1000 Genomes Phase 1 LAI", "",
          "Mapping of our 8 ancestries to 3: AFR + AHG → AFR, EUR + WAS → EUR, NAT → NAT; EAS + SAS + OCE → other "
          "(dropped and renormalised for 3-way comparisons).", "",
          "### Martin et al. 2017, population means (1000 Genomes Phase 3 unrelated samples)", "",
          mp.to_markdown(index=False), "",
          "### Samples named by Martin et al. 2017", "", mn.drop(columns=["martin_statement"]).to_markdown(index=False), "",
          "### Admixture timing (rough single-pulse estimate vs Martin et al. Tracts estimates)", "",
          tim.to_markdown(index=False), "",
          "### 1000 Genomes Phase 1 consensus LAI: global proportions", "", p1s.to_markdown(index=False), "",
          "### 1000 Genomes Phase 1 consensus LAI: diploid local ancestry concordance", "",
          cs.to_markdown(index=False, floatfmt=".4f"), ""]
    with atomic_path(vdir / "COMPARISON.md") as tmp:
        tmp.write_text("\n".join(md))
    log(mp.to_string(index=False))
    log(p1s.to_string(index=False))
    log(cs.to_string(index=False))
