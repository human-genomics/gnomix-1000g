"""Write output/REPORT.md and output/MANIFEST.sha256."""

from __future__ import annotations

import platform
from datetime import datetime, timezone

import pandas as pd

from . import __version__
from .common import ANCESTRIES, OUTPUT, TOOLS_DIR, WORK, atomic_path, log, read_json, sha256_file

EXAMPLES = ["HG01893", "NA20314", "HG01880", "NA19700", "HG00096", "NA18525", "HG03006"]


def _versions() -> dict:
    import numpy
    import sklearn
    import xgboost
    v = {"gnomix1000g": __version__, "python": platform.python_version(), "numpy": numpy.__version__,
         "scikit-learn": sklearn.__version__, "xgboost": xgboost.__version__}
    commit = TOOLS_DIR / "gnomix" / ".commit"
    if commit.exists():
        v["gnomix commit"] = commit.read_text().strip()
    return v


def run(chroms: list[int]) -> None:
    from . import figures
    figures.run(chroms)
    md = ["# gnomix-1000g results", "",
          f"Generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}. "
          f"Chromosomes: {', '.join(map(str, chroms))}.", ""]
    md += ["## Software", "", pd.DataFrame([_versions()]).T.rename(columns={0: "version"}).to_markdown(), ""]

    cov = pd.read_csv(WORK / "match" / "summary.tsv", sep="\t")
    cov = cov[cov.chrom.isin(chroms)]
    md += ["## Model SNPs found in the 1000 Genomes panel", "",
           f"Total: {cov.matched.sum():,} of {cov.model_snps.sum():,} model SNPs "
           f"({100 * cov.matched.sum() / cov.model_snps.sum():.2f}%). Gnomix asks for at least 80%.", "",
           cov.to_markdown(index=False), ""]

    vdir = OUTPUT / "validation"
    if (vdir / "trio_benchmark_summary.tsv").exists():
        s = pd.read_csv(vdir / "trio_benchmark_summary.tsv", sep="\t")
        pol = read_json(vdir / "gnofix_policy.json")
        cols = {"level": "level", "group": "group or population", "children": "children",
                "het_ancestry_fraction": "genome with 2 ancestries",
                "hap_err_rephased": "hap. ancestry error, rephased", "hap_err_gnofix": "hap. error, Gnofix",
                "hap_err_gnofix_on_truth": "hap. error, Gnofix on true phase",
                "hap_err_continental_rephased": "continental error, rephased",
                "hap_err_continental_gnofix": "continental error, Gnofix",
                "phase_acc_10cM_hetanc_rephased": "10 cM phase acc., rephased",
                "phase_acc_10cM_hetanc_gnofix": "10 cM phase acc., Gnofix",
                "SER_rephased": "switch error rate, rephased", "SER_gnofix": "switch error rate, Gnofix",
                "gnofix": "Gnofix used"}
        t = s[list(cols)].rename(columns=cols)
        md += ["## Trio phasing benchmark (Gnofix decision)", "",
               ("Benchmark run in this data folder. " if pol.get("source") == "benchmark" else
                "Benchmark not re-run (RUN_TRIO_BENCHMARK=0): policy and results of the published run "
                "(repository folder reference/). "),
               f"Chromosomes {pol['chromosomes']}. Rule: {pol['rule']}", "",
               "Errors are fractions of the genome (cM-weighted); 10 cM phase accuracy is for heterozygous "
               "site pairs in regions where the two homologs have different ancestry. Full table: "
               "`validation/trio_benchmark_summary.tsv`.", "",
               t.to_markdown(index=False, floatfmt=".3f"), "",
               "![trio benchmark](validation/trio_benchmark.png)", ""]
    if (vdir / "COMPARISON.md").exists():
        body = (vdir / "COMPARISON.md").read_text().split("\n", 1)[1]
        md += ["## Comparison with published results", body, "",
               "![global ancestry](validation/admixed_global_ancestry.png)", "",
               "![phase 1](validation/compare_phase1_global.png)", ""]
    md += ["![structure](figures/structure_all_samples.png)", ""]
    pop = pd.read_csv(OUTPUT / "population_ancestry.tsv", sep="\t")
    md += ["## Mean global ancestry per population (final calls)", "",
           pop.to_markdown(index=False, floatfmt=".3f"), ""]

    glob = pd.read_csv(OUTPUT / "global_ancestry_final.tsv", sep="\t").set_index("IID")
    md += ["## Example karyograms", ""]
    for iid in EXAMPLES:
        if iid in glob.index:
            p = glob.loc[iid, "Population"]
            md += [f"**{iid}** ({p}): " + ", ".join(f"{a} {100 * glob.loc[iid, a]:.1f}%" for a in ANCESTRIES
                                                  if glob.loc[iid, a] >= 0.01), "",
                   f"![{iid}](karyograms/{p}/{iid}.png)", ""]
    pf = OUTPUT / "pfile" / "verification.json"
    if pf.exists():
        md += ["## Phase-corrected panel", "", "```", pf.read_text(), "```", ""]
    with atomic_path(OUTPUT / "REPORT.md") as tmp:
        tmp.write_text("\n".join(md))

    # Karyograms (3,202 PNGs) and release copies are left out; so are hidden state files.
    files = sorted(p for p in OUTPUT.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256"
                   and not {"karyograms", "release"} & set(p.relative_to(OUTPUT).parts)
                   and not any(part.startswith(".") for part in p.relative_to(OUTPUT).parts))
    with atomic_path(OUTPUT / "MANIFEST.sha256") as tmp:
        tmp.write_text("".join(f"{sha256_file(p)}  {p.relative_to(OUTPUT)}\n" for p in files))
    log(f"  wrote {OUTPUT / 'REPORT.md'} and MANIFEST.sha256 ({len(files)} files)")
