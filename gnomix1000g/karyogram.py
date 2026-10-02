"""Karyograms for every sample, in the style of Martin et al. (2017).

Each chromosome is drawn as two vertical bars, one per homolog (A, B). The
homolog order is arbitrary: A and B are not maternal and paternal. Colours are
the final ancestry calls (the calls of tracts_final); coordinates are GRCh38.
Opacity follows the posterior of the displayed category (for West Eurasian:
P(EUR) + P(WAS)): 1.0 at >= 0.9, 0.95 at 0.8, 0.72 at 0.6, 0.5 at <= 0.35.
Admixed populations (ACB, ASW, CLM, MXL, PEL, PUR) show EUR + WAS as one West
Eurasian category, because the model separates them poorly, and the rare EAS,
SAS, AHG and OCE calls as "Other". Karyograms are redrawn when the calls or the
Gnofix policy change, or with REDRAW=1.

Output: output/karyograms/<Population>/<IID>.png
"""

from __future__ import annotations

import gzip
import os
from multiprocessing import get_context

import numpy as np
import pandas as pd

from .common import (ADMIXED_POPULATIONS, ANCESTRIES, ANCESTRY_COLORS, ANCESTRY_NAMES, DOWNLOADS, OTHER,
                     OUTPUT, POPULATION_NAMES, WEST_EURASIAN, atomic_path, inputs_stamp, log, stamp_is_current,
                     workers_for_memory, write_stamp)
from .prepare import write_samples
from .tracts import calls, final_uses_gnofix

MIN_ALPHA, ALPHA_STEP = 0.5, 0.025  # opacity floor; quantization so runs of windows merge
MAX_DRAWN_SPAN = 5e6  # bp: longer window spans (gaps between lifted SNPs) stay grey
WORKER_GB = 0.5
NO_DATA = "#E9E9E9"
FOOTER = "github.com/human-genomics/gnomix-1000g"
BAR, GAP, STEP = 0.34, 0.07, 1.12  # bar width, gap between homologs, distance between chromosomes

_DATA: dict = {}


def cytobands() -> tuple[dict, dict]:
    """Chromosome lengths and centromere (start, end) in GRCh38 from UCSC cytoBand."""
    lengths, cen = {}, {}
    with gzip.open(DOWNLOADS / "cytoBand.hg38.txt.gz", "rt") as fh:
        for line in fh:
            c, s, e, _, stain = line.rstrip("\n").split("\t")
            if not c[3:].isdigit():
                continue
            k = int(c[3:])
            lengths[k] = max(lengths.get(k, 0), int(e))
            if stain == "acen":
                a, b = cen.get(k, (np.inf, 0))
                cen[k] = (min(a, int(s)), max(b, int(e)))
    return lengths, cen


def drawn_spans(chrom: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """GRCh38 start and end to draw for each window, and whether to draw it."""
    win = pd.read_csv(OUTPUT / "windows" / f"chr{chrom}.windows.tsv", sep="\t")
    good = ~win.hg38_discordant.values
    s = win.spos_hg38.values.astype(float)
    # Draw each window up to the next window's start (no gaps), unless that window
    # has no GRCh38 coordinates; then up to its own last SNP.
    nxt = np.append(s[1:], -1)
    nxt_ok = np.append(good[1:], False) & (nxt > s)
    e = np.where(nxt_ok, nxt, win.epos_hg38.values).astype(float)
    return s, e, good & (e > s) & (e - s < MAX_DRAWN_SPAN)


def load_chrom(chrom: int, use_fix: np.ndarray) -> dict:
    """Final calls of a chromosome, with drawable GRCh38 spans per window."""
    lab, prob = calls(chrom, use_fix)["final"]
    s, e, ok = drawn_spans(chrom)
    return {"lab": lab, "prob": prob, "s": s, "e": e, "ok": ok}


def opacity(post: np.ndarray) -> np.ndarray:
    """Opacity for posteriors in [0, 1]: smooth, nearly solid down to 0.8, MIN_ALPHA at <= 0.35."""
    x = np.clip((post - 0.35) / 0.55, 0, 1)
    alpha = MIN_ALPHA + (1 - MIN_ALPHA) * x * x * (3 - 2 * x)  # smoothstep
    return np.round(alpha / ALPHA_STEP) * ALPHA_STEP


def categories(population: str) -> list[tuple[str, str, list[int]]]:
    """Displayed categories: (label, colour, model ancestry codes), in legend order."""
    code = {a: i for i, a in enumerate(ANCESTRIES)}
    one = lambda a: [(f"{ANCESTRY_NAMES[a]} ({a})", ANCESTRY_COLORS[a], [code[a]])]  # noqa: E731
    if population in ADMIXED_POPULATIONS:
        # Three sources plus "other": four colours that stay distinct for colour-blind readers.
        return (one("AFR") + [("West Eurasian (EUR + WAS)", WEST_EURASIAN, [code["EUR"], code["WAS"]])]
                + one("NAT") + [("Other (EAS, SAS, AHG, OCE)", OTHER,
                                 [code["EAS"], code["SAS"], code["AHG"], code["OCE"]])])
    return sum((one(a) for a in ("AFR", "AHG", "EUR", "WAS", "SAS", "EAS", "NAT", "OCE")), [])


def runs(cat: np.ndarray, alpha: np.ndarray, ok: np.ndarray, s: np.ndarray, e: np.ndarray) -> list:
    """Merge consecutive windows with the same category and opacity into (start, end, cat, alpha)."""
    idx = np.nonzero(ok)[0]
    if not len(idx):
        return []
    key = cat[idx] * 1000 + np.rint(alpha[idx] / ALPHA_STEP)
    brk = np.nonzero((np.diff(key) != 0) | (np.diff(idx) != 1) | (s[idx[1:]] != e[idx[:-1]]))[0] + 1
    starts = np.concatenate([[0], brk])
    ends = np.concatenate([brk, [len(idx)]]) - 1
    return [(s[idx[a]], e[idx[b]], int(cat[idx[a]]), float(alpha[idx[a]])) for a, b in zip(starts, ends)]


def _draw_chromosomes(ax, fig, k: int, cats: list) -> None:
    from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle

    d = _DATA
    lengths, cen, chroms = d["lengths"], d["cen"], d["chroms"]
    to_cat = np.zeros(len(ANCESTRIES), dtype=np.int64)
    member = np.zeros((len(ANCESTRIES), len(cats)))  # ancestry -> displayed category
    for i, (_, _, codes) in enumerate(cats):
        to_cat[codes] = i
        member[codes, i] = 1
    xmax = (max(chroms) - 1) * STEP + 2 * BAR + GAP
    ymax = max(lengths[c] for c in chroms) / 1e6
    ax.set_xlim(-0.35, xmax + 0.35)
    ax.set_ylim(ymax + 6, -14)
    # radius of the rounded ends in data units: equal in x and y on screen
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    y_per_x = ((ymax + 20) / bbox.height) / ((xmax + 0.7) / bbox.width)
    rounded = dict(boxstyle=f"round,pad=0,rounding_size={BAR / 2}", mutation_aspect=y_per_x)
    for c in chroms:
        x0 = (c - 1) * STEP
        L = lengths[c] / 1e6
        ch = d["data"][c]
        for h in (0, 1):
            x = x0 + h * (BAR + GAP)
            shape = FancyBboxPatch((x, 0), BAR, L, facecolor=NO_DATA, edgecolor="none", zorder=1, **rounded)
            ax.add_patch(shape)
            cat = to_cat[ch["lab"][2 * k + h]]
            cat_post = (ch["prob"][2 * k + h] / 255) @ member  # posterior of each displayed category
            alpha = opacity(np.take_along_axis(cat_post, cat[:, None], axis=1)[:, 0])
            for ys, ye, ci, al in runs(cat, alpha, ch["ok"], ch["s"], ch["e"]):
                r = Rectangle((x, ys / 1e6), BAR, (ye - ys) / 1e6, facecolor=cats[ci][1], alpha=al,
                              edgecolor="none", zorder=2, antialiased=False)
                ax.add_patch(r)
                r.set_clip_path(shape)
            if c in cen:  # centromere: white notches on both sides
                cs, ce = cen[c][0] / 1e6, cen[c][1] / 1e6
                for side, sign in ((x - 0.004, 1), (x + BAR + 0.004, -1)):
                    ax.add_patch(Polygon([(side, cs - 1.5), (side + sign * BAR * 0.3, (cs + ce) / 2),
                                          (side, ce + 1.5)], closed=True, facecolor="white",
                                         edgecolor="none", zorder=3))
            ax.add_patch(FancyBboxPatch((x, 0), BAR, L, facecolor="none", edgecolor="#5A5A5A", linewidth=0.6,
                                        zorder=4, **rounded))
        ax.text(x0 + BAR + GAP / 2, -4, str(c), ha="center", va="bottom", fontsize=11, color="#333333")
    ax.set_xticks([])
    ax.set_yticks(np.arange(0, ymax + 1, 50))
    ax.tick_params(axis="y", labelsize=10, colors="#555555", length=3)
    ax.set_ylabel("Mb (GRCh38)", fontsize=10, color="#555555")
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#AAAAAA")


def _draw_titles(fig, row, gnofix: bool) -> None:
    title = fig.text(0.045, 0.955, row.IID, fontsize=18, fontweight="bold", color="#111111", va="top")
    fig.canvas.draw()  # measure the title so the subtitle follows it
    x = fig.transFigure.inverted().transform(title.get_window_extent())[1, 0] + 0.012
    pop_name = POPULATION_NAMES.get(row.Population, row.Population)
    fig.text(x, 0.951, f"{pop_name} ({row.Population}, {row.SuperPop})", fontsize=13, color="#333333", va="top")
    phase = "with Gnofix phase correction" if gnofix else "on the published 1000 Genomes phase"
    fig.text(0.045, 0.905, f"Pretrained Gnomix {phase}. Two bars per chromosome: homologs A and B "
             "(not maternal/paternal).", fontsize=10, color="#666666", va="top")


def _draw_legend(fig, cats: list, glob_row) -> None:
    from matplotlib.patches import Rectangle

    def swatch(x, y, color, alpha=1.0, width=0.016):
        fig.patches.append(Rectangle((x, y - 0.012), width, 0.03, transform=fig.transFigure, facecolor=color,
                                     alpha=alpha, edgecolor="none"))

    lx, ly = 0.775, 0.80
    fig.text(lx, ly + 0.05, "Genome-wide ancestry", fontsize=12, fontweight="bold", color="#222222")
    shares = [100 * sum(glob_row[ANCESTRIES[i]] for i in codes) for _, _, codes in cats]
    for (label, color, _), pct in zip(cats, shares):
        if pct < 0.05:
            continue
        swatch(lx, ly, color)
        fig.text(lx + 0.024, ly, f"{pct:5.1f}%", fontsize=11, color="#111111", va="center",
                 family="DejaVu Sans Mono")
        fig.text(lx + 0.078, ly, label, fontsize=10, color="#333333", va="center")
        ly -= 0.06
    ly -= 0.02
    fig.text(lx, ly, "Posterior of the shown ancestry", fontsize=10, color="#555555", va="center")
    ly -= 0.045
    main_color = cats[int(np.argmax(shares))][1]
    grid = np.linspace(0.3, 1.0, 36)
    for j, al in enumerate(opacity(grid)):
        swatch(lx + j * 0.0055, ly, main_color, al, width=0.0056)
    for p in (0.4, 0.6, 0.8, 1.0):
        fig.text(lx + (p - 0.3) / 0.7 * 35 * 0.0055, ly - 0.035, f"{p:.1f}", fontsize=8.5, color="#555555",
                 va="center", ha="center")
    ly -= 0.1
    swatch(lx, ly, NO_DATA)
    fig.text(lx + 0.024, ly, "no model SNPs / not lifted to GRCh38", fontsize=10, color="#555555", va="center")
    fig.text(lx, 0.06, FOOTER, fontsize=9, color="#999999")


def plot_sample(k: int) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = _DATA
    row = d["samples"].iloc[k]
    out = OUTPUT / "karyograms" / row.Population / f"{row.IID}.png"
    if out.exists() and os.environ.get("REDRAW") != "1":
        return str(out)
    cats = categories(row.Population)
    fig = plt.figure(figsize=(15, 6.4), dpi=120, facecolor="white")
    ax = fig.add_axes([0.045, 0.06, 0.705, 0.80])
    _draw_chromosomes(ax, fig, k, cats)
    _draw_titles(fig, row, bool(d["use_fix"][k]))
    _draw_legend(fig, cats, d["global_ancestry"].iloc[k])
    out.parent.mkdir(parents=True, exist_ok=True)
    with atomic_path(out) as tmp:
        fig.savefig(tmp, facecolor="white", format="png")
    plt.close(fig)
    return str(out)


def run(chroms: list[int], sample_ids=None) -> None:
    samples = write_samples()
    use_fix = final_uses_gnofix(samples)
    stamp = inputs_stamp(chroms, use_fix)
    kdir = OUTPUT / "karyograms"
    if sample_ids is None and not stamp_is_current(kdir, stamp) and any(kdir.glob("*/*.png")):
        log("  calls or Gnofix policy changed since the karyograms were drawn: redrawing all")
        os.environ["REDRAW"] = "1"
    lengths, cen = cytobands()
    glob = pd.read_csv(OUTPUT / "global_ancestry_final.tsv", sep="\t").set_index("IID").loc[samples.IID]
    _DATA.update(samples=samples, use_fix=use_fix, lengths=lengths, cen=cen, chroms=chroms,
                 global_ancestry=glob.reset_index(), data={c: load_chrom(c, use_fix) for c in chroms})
    todo = range(len(samples)) if sample_ids is None else \
        [int(i) for i in np.nonzero(samples.IID.isin(sample_ids).values)[0]]
    n = workers_for_memory(WORKER_GB)
    log(f"  drawing {len(todo)} karyograms with {n} workers")
    with get_context("fork").Pool(n) as pool:  # fork: workers share the loaded calls
        for i, _ in enumerate(pool.imap_unordered(plot_sample, todo, chunksize=8), 1):
            if i % 500 == 0:
                log(f"    {i}/{len(todo)} karyograms")
    if sample_ids is None:
        write_stamp(kdir, stamp)
