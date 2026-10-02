"""Build the browsable karyogram site (`site` step) for GitHub Pages.

Writes output/site/: index.html, data.js and karyograms/<POP>/<IID>.png (copied from output/karyograms).
The page browses one population at a time, with global ancestry, explicit pedigree relationships and
Gnofix status. Sample links open karyograms and their relatives without leaving the page.
"""

from __future__ import annotations

import json
import shutil

import pandas as pd

from .common import ANCESTRIES, ANCESTRY_COLORS, ANCESTRY_NAMES, OUTPUT, REFERENCE, WORK, atomic_path, log

SITE = OUTPUT / "site"
MARTIN_URL = "https://personal.broadinstitute.org/armartin/tgp_admixture/karyograms/{}.pdf"
REPO_URL = "https://github.com/human-genomics/gnomix-1000g"
RELEASE_URL = f"{REPO_URL}/releases/tag/v1.0.0"

SUPERPOP_NAMES = {"AFR": "African", "AMR": "Admixed American", "EAS": "East Asian",
                  "EUR": "European", "SAS": "South Asian"}
SUPERPOP_ORDER = ["AMR", "AFR", "EUR", "SAS", "EAS"]


def population_names() -> dict:
    """Population code -> full description (reference/kg_population_names.tsv)."""
    d = pd.read_csv(REFERENCE / "kg_population_names.tsv", sep="\t")
    d = d[d["Population Code"].notna()]
    return dict(zip(d["Population Code"], d["Population Description"]))


def sample_rows(glob: pd.DataFrame, samples: pd.DataFrame, martin: set) -> list[dict]:
    """One record per sample: ancestry fractions, family links and whether Martin et al. drew it."""
    fam = samples.set_index("IID")[["PAT", "MAT"]]
    ids = set(samples.IID)
    rows = []
    for r in glob.itertuples():
        pat, mat = fam.loc[r.IID, "PAT"], fam.loc[r.IID, "MAT"]
        parents = [p for p in (pat, mat) if p in ids]
        rows.append({
            "id": r.IID, "pop": r.Population, "sup": r.SuperPop,
            "a": [round(float(getattr(r, a)), 4) for a in ANCESTRIES],
            "par": parents, "fix": bool(r.gnofix_in_final), "m": r.IID in martin,
        })
    kids = {r["id"]: r["par"] for r in rows if r["par"]}
    children = {}
    for kid, par in kids.items():
        for p in par:
            children.setdefault(p, []).append(kid)
    for r in rows:
        r["kids"] = children.get(r["id"], [])
    by_id = {r["id"]: r for r in rows}
    for r in rows:
        roles = []
        if r["par"]:
            roles.append("Trio child" if len(r["par"]) == 2 else "Duo child")
        for size, label in ((2, "Trio parent"), (1, "Duo parent")):
            if any(len(by_id[k]["par"]) == size for k in r["kids"]):
                roles.append(label)
        r["family"] = roles or ["Unrelated"]
    return rows


def copy_karyograms() -> int:
    """Copy the karyograms into the site folder (the site is self-contained and publishable)."""
    src, dst = OUTPUT / "karyograms", SITE / "karyograms"
    n = 0
    for p in sorted(src.glob("*/*.png")):
        out = dst / p.parent.name / p.name
        if out.exists() and out.stat().st_size == p.stat().st_size:
            n += 1
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, out)
        n += 1
    return n


def run() -> None:
    SITE.mkdir(parents=True, exist_ok=True)
    glob = pd.read_csv(OUTPUT / "global_ancestry_final.tsv", sep="\t")
    samples = pd.read_csv(WORK / "samples.tsv", sep="\t", dtype={"PAT": str, "MAT": str})
    martin = set()
    mf = REFERENCE / "martin2017_karyogram_samples.txt"
    if mf.exists():
        martin = {x.strip() for x in mf.read_text().split() if x.strip()}
    rows = sample_rows(glob, samples, martin)
    pops = population_names()
    pop_meta = []
    for (sup, pop), d in glob.groupby(["SuperPop", "Population"]):
        pop_meta.append({"pop": pop, "sup": sup, "name": pops.get(pop, pop), "n": len(d),
                         "mean": [round(float(d[a].mean()), 4) for a in ANCESTRIES]})
    pop_meta.sort(key=lambda x: (SUPERPOP_ORDER.index(x["sup"]), x["pop"]))
    data = {
        "samples": rows, "pops": pop_meta,
        "anc": [{"code": a, "name": ANCESTRY_NAMES[a], "color": ANCESTRY_COLORS[a]} for a in ANCESTRIES],
        "sup": {k: SUPERPOP_NAMES.get(k, k) for k in SUPERPOP_ORDER},
        "martin_url": MARTIN_URL, "n_martin": sum(r["m"] for r in rows),
        "n_trio": sum(1 for r in rows if len(r["par"]) == 2),
    }
    with atomic_path(SITE / "data.js") as tmp:
        tmp.write_text("window.DATA = " + json.dumps(data, separators=(",", ":")) + ";\n")
    with atomic_path(SITE / "index.html") as tmp:
        tmp.write_text(page_html())
    n = copy_karyograms()
    log(f"  site: {len(rows)} samples, {n} karyograms, {data['n_martin']} with a Martin et al. karyogram")
    log(f"  wrote {SITE}/index.html")


def page_html() -> str:
    from .site_page import HTML
    return HTML
