"""Summary figures for the report and the README (output/figures/)."""

from __future__ import annotations

import shutil

import numpy as np
import pandas as pd

from .common import ADMIXED_POPULATIONS, ANCESTRY_COLORS, ANCESTRY_NAMES, OUTPUT, atomic_path, log

SUPER_ORDER = ["AFR", "AMR", "EUR", "SAS", "EAS"]
POP_ORDER = {"AFR": ["YRI", "ESN", "GWD", "MSL", "LWK", "ACB", "ASW"],
             "AMR": ["PEL", "MXL", "CLM", "PUR"],
             "EUR": ["GBR", "CEU", "FIN", "IBS", "TSI"],
             "SAS": ["PJL", "GIH", "BEB", "ITU", "STU"],
             "EAS": ["CHB", "CHS", "CDX", "JPT", "KHV"]}
STACK = ["AFR", "AHG", "EUR", "WAS", "SAS", "EAS", "NAT", "OCE"]
# HG01893 is the example karyogram of Martin et al. (2017, Fig. 1). The other examples are the
# samples without relatives in the panel whose global ancestry is closest to their population mean.
KARYOGRAM_FIXED = [("PEL", "HG01893")]
KARYOGRAM_TYPICAL = ["MXL", "ASW", "PUR"]
GNOFIX_EXAMPLE_QUANTILE = 0.8  # README before/after example: more change than average, not the extreme
NO_DATA = "#E9E9E9"
ACROCENTRIC = [13, 14, 15, 21, 22]


def typical_samples(glob: pd.DataFrame) -> list[tuple[str, str]]:
    out = []
    free = glob[~glob.has_parent_in_panel & ~glob.has_child_in_panel]
    for pop in KARYOGRAM_TYPICAL:
        d = free[free.Population == pop]
        dist = ((d[STACK] - glob.loc[glob.Population == pop, STACK].mean()) ** 2).sum(axis=1)
        out.append((pop, d.loc[dist.idxmin(), "IID"]))
    return out


def _style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#AAAAAA")
    ax.spines["bottom"].set_color("#AAAAAA")
    ax.tick_params(colors="#444444", labelsize=8)


def structure_plot(glob: pd.DataFrame, path) -> None:
    """All samples as stacked bars of global ancestry, grouped by population."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    rows, labels, x = [], [], 0
    for sp in SUPER_ORDER:
        for pop in POP_ORDER[sp]:
            d = glob[glob.Population == pop]
            # sort by the population's main ancestry, then the second
            main = d[STACK].mean().sort_values(ascending=False).index[:2]
            d = d.sort_values(list(main), ascending=False)
            rows.append(d[STACK].values)
            labels.append((pop, sp, x, x + len(d)))
            x += len(d) + 6
    fig, ax = plt.subplots(figsize=(16, 4.2), dpi=130)
    for (pop, sp, a, b), vals in zip(labels, rows):
        bottom = np.zeros(len(vals))
        xs = np.arange(a, b)
        for i, anc in enumerate(STACK):
            ax.bar(xs, vals[:, i], bottom=bottom, width=1.0, color=ANCESTRY_COLORS[anc], linewidth=0)
            bottom += vals[:, i]
        ax.text((a + b) / 2, -0.04, pop, ha="center", va="top", fontsize=9, color="#333333",
                fontweight="bold" if pop in ADMIXED_POPULATIONS else "normal")
    for sp in SUPER_ORDER:
        xs = [(a, b) for pop, s, a, b in labels if s == sp]
        ax.text((xs[0][0] + xs[-1][1]) / 2, -0.17, sp, ha="center", va="top", fontsize=11, color="#111111")
        ax.plot([xs[0][0], xs[-1][1]], [-0.135, -0.135], color="#999999", linewidth=0.8, clip_on=False)
    ax.set_xlim(-3, x - 3)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0%", "50%", "100%"])
    _style(ax)
    ax.spines["bottom"].set_visible(False)
    ax.legend(handles=[Patch(color=ANCESTRY_COLORS[a], label=f"{ANCESTRY_NAMES[a]} ({a})") for a in STACK],
              loc="upper center", bbox_to_anchor=(0.5, 1.2), ncol=8, frameon=False, fontsize=10,
              handlelength=1.2, columnspacing=1.2)
    ax.set_title(f"Genome-wide ancestry of all {len(glob):,} samples (pretrained Gnomix; admixed populations "
                 "in bold)", fontsize=12, loc="left", pad=32, color="#222222")
    ax.tick_params(axis="y", labelsize=9)
    fig.tight_layout()
    with atomic_path(path) as tmp:
        fig.savefig(tmp, facecolor="white", format="png")
    plt.close(fig)


def _category_runs(lab: np.ndarray, to_cat: np.ndarray) -> int:
    """Number of runs of equal displayed category, summed over haplotypes (rows)."""
    c = to_cat[lab]
    return int((np.diff(c, axis=1) != 0).sum() + len(c))


def gnofix_example_choice(glob: pd.DataFrame, chroms: list[int]) -> tuple[str, str, int, pd.DataFrame]:
    """An admixed sample that Gnofix changed more than average, and its most changed chromosome.

    Change = tracts (runs of the displayed karyogram category) removed by Gnofix. Among the samples
    whose final calls use Gnofix we take the one at GNOFIX_EXAMPLE_QUANTILE of the genome-wide change.
    """
    from .common import WORK
    from .karyogram import categories

    use = glob[glob.gnofix_in_final].reset_index()
    to_cat = {}
    rows = []
    for c in chroms:
        with np.load(WORK / "infer" / f"chr{c}.npz") as z:
            lab_raw, lab_fix = z["lab_raw"], z["lab_fix"]
        for k, r in use.iterrows():
            if r.Population not in to_cat:
                t = np.zeros(8, dtype=np.int64)
                for i, (_, _, codes) in enumerate(categories(r.Population)):
                    t[codes] = i
                to_cat[r.Population] = t
            i = int(r["index"])
            rows.append((r.IID, r.Population, c, _category_runs(lab_raw[2 * i:2 * i + 2], to_cat[r.Population]),
                         _category_runs(lab_fix[2 * i:2 * i + 2], to_cat[r.Population])))
    df = pd.DataFrame(rows, columns=["IID", "Population", "chrom", "tracts_raw", "tracts_gnofix"])
    df["removed"] = df.tracts_raw - df.tracts_gnofix
    per = df.groupby(["IID", "Population"]).removed.sum().sort_values().reset_index()
    pick = per.iloc[int(round(GNOFIX_EXAMPLE_QUANTILE * (len(per) - 1)))]
    d = df[df.IID == pick.IID]
    if (~d.chrom.isin(ACROCENTRIC)).any():  # their short arms have no model SNPs: mostly grey
        d = d[~d.chrom.isin(ACROCENTRIC)]
    d = d.sort_values(["removed", "chrom"], ascending=[False, True])
    return pick.IID, pick.Population, int(d.chrom.iloc[0]), per


def gnofix_example(glob: pd.DataFrame, chroms: list[int], path) -> str:
    """One chromosome of one sample, published phase vs Gnofix (as in the Gnofix README figure)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, Patch, Rectangle

    from .common import WORK
    from .karyogram import categories, cytobands, drawn_spans, runs

    iid, pop, c, per = gnofix_example_choice(glob, chroms)
    k = int(np.nonzero(glob.IID.values == iid)[0][0])  # global_ancestry rows follow the panel order
    with np.load(WORK / "infer" / f"chr{c}.npz") as z:
        labs = {"raw": z["lab_raw"][2 * k:2 * k + 2], "fix": z["lab_fix"][2 * k:2 * k + 2]}
    cats = categories(pop)
    to_cat = np.zeros(8, dtype=np.int64)
    for i, (_, _, codes) in enumerate(cats):
        to_cat[codes] = i
    s, e, ok = drawn_spans(c)
    L = cytobands()[0][c] / 1e6
    H = 0.62  # bar height (y units); bars at y = 1 (A) and 0 (B)

    fig, axes = plt.subplots(1, 2, figsize=(12, 2.9), dpi=130)
    fig.subplots_adjust(left=0.04, right=0.98, top=0.72, bottom=0.25, wspace=0.16)
    used = set()
    for ax, key, title in zip(axes, ("raw", "fix"), ("Published phase", "After Gnofix")):
        ax.set_xlim(-0.02 * L, 1.02 * L)
        ax.set_ylim(-0.6, 1.6)
        bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        y_per_x = (2.2 / bbox.height) / (1.04 * L / bbox.width)
        rounded = dict(boxstyle=f"round,pad=0,rounding_size={H / 2 / y_per_x}", mutation_aspect=y_per_x)
        n = 0
        for h, y in ((0, 1), (1, 0)):
            shape = FancyBboxPatch((0, y - H / 2), L, H, facecolor=NO_DATA, edgecolor="none", **rounded)
            ax.add_patch(shape)
            cat = to_cat[labs[key][h]]
            n += _category_runs(labs[key][h:h + 1], to_cat)
            for xs, xe, ci, _ in runs(cat, np.ones(len(cat)), ok, s, e):
                r = Rectangle((xs / 1e6, y - H / 2), (xe - xs) / 1e6, H, facecolor=cats[ci][1], edgecolor="none",
                              antialiased=False)
                ax.add_patch(r)
                r.set_clip_path(shape)
                used.add(ci)
            ax.add_patch(FancyBboxPatch((0, y - H / 2), L, H, facecolor="none", edgecolor="#5A5A5A",
                                        linewidth=0.7, **rounded))
            ax.text(-0.015 * L, y, "AB"[h], ha="right", va="center", fontsize=10, color="#555555")
        ax.set_title(f"{title}: {n} tracts", fontsize=11.5, loc="left", color="#222222")
        ax.axis("off")
    fig.text(0.04, 0.93, f"{iid} ({pop}), chromosome {c}", fontsize=13, fontweight="bold", color="#111111",
             va="top")
    fig.text(0.505, 0.47, "→", fontsize=26, color="#444444", ha="center", va="center")
    fig.legend(handles=[Patch(color=cats[i][1], label=cats[i][0]) for i in sorted(used)] +
               [Patch(color=NO_DATA, label="no model SNPs")], loc="lower center", ncol=len(used) + 1,
               frameon=False, fontsize=9.5)
    with atomic_path(path) as tmp:
        fig.savefig(tmp, facecolor="white", format="png")
    plt.close(fig)
    with atomic_path(OUTPUT / "validation" / "gnofix_example_choice.tsv") as tmp:
        per.to_csv(tmp, sep="\t", index=False)
    return f"{iid} ({pop}) chr{c}"


def run(chroms: list[int]) -> None:
    fdir = OUTPUT / "figures"
    fdir.mkdir(parents=True, exist_ok=True)
    glob = pd.read_csv(OUTPUT / "global_ancestry_final.tsv", sep="\t")
    structure_plot(glob, fdir / "structure_all_samples.png")
    vdir = OUTPUT / "validation"
    for name in ("trio_benchmark.png", "compare_phase1_global.png"):
        if (vdir / name).exists():
            shutil.copyfile(vdir / name, fdir / name)
    if glob.gnofix_in_final.any() and (OUTPUT / "windows").exists():
        log(f"  Gnofix example: {gnofix_example(glob, chroms, fdir / 'gnofix_example.png')}")
    for pop, iid in KARYOGRAM_FIXED + typical_samples(glob):
        src = OUTPUT / "karyograms" / pop / f"{iid}.png"
        if src.exists():
            shutil.copyfile(src, fdir / f"karyogram_{pop}_{iid}.png")
