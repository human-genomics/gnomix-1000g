"""Command line: python -m gnomix1000g <step> [--chroms 1-22]. main.sh runs "all"."""

from __future__ import annotations

import argparse
import os

from .common import chrom_list, log, single_thread_env


def _prepare(chroms, args):
    from . import checks, prepare
    prepare.unpack_models(chroms)
    prepare.write_model_meta(chroms)
    prepare.write_samples()
    prepare.split_pvar()
    prepare.decompress_pgen()
    prepare.match_all(chroms)
    log("  checks against the official Gnomix path (chr22 demo data)")
    checks.run()


def _infer(chroms, args):
    from . import infer
    infer.run(chroms, args.batch_size)


def _trio(chroms, args):
    from . import trio
    if os.environ.get("RUN_TRIO_BENCHMARK", "0") == "1":
        trio.run()
    else:
        log("  [skip] benchmark (RUN_TRIO_BENCHMARK=0): using the published policy, reference/gnofix_policy.json")
        trio.use_published_policy()


def _tracts(chroms, args):
    from . import tracts
    tracts.run(chroms)


def _karyograms(chroms, args):
    from . import karyogram
    karyogram.run(chroms)


def _compare(chroms, args):
    from . import compare
    compare.run(chroms)


def _pfile(chroms, args):
    from . import pfile
    if os.environ.get("MAKE_PFILE", "1") == "1":
        pfile.run(chroms)
    else:
        log("  [skip] phase-corrected panel (MAKE_PFILE=0)")


def _report(chroms, args):
    from . import report
    report.run(chroms)


def _site(chroms, args):
    from . import site
    site.run()


def _release(chroms, args):
    from . import release
    release.run(chroms)


# (name, description, function, part of "all")
STEPS = [
    ("prepare", "models, panel index, model SNP matching", _prepare, True),
    ("infer", "pretrained Gnomix with and without Gnofix, all samples", _infer, True),
    ("trio", "Gnofix policy (trio phasing benchmark)", _trio, True),
    ("tracts", "msp files, ancestry tracts, global ancestry", _tracts, True),
    ("karyograms", "one karyogram per sample", _karyograms, True),
    ("compare", "Martin et al. (2017) and 1000 Genomes Phase 1", _compare, True),
    ("pfile", "phase-corrected panel", _pfile, True),
    ("report", "REPORT.md, figures, checksums", _report, True),
    ("site", "browsable karyogram site (GitHub Pages)", _site, False),
    ("release", "release assets", _release, False),
]


def main(argv=None) -> None:
    single_thread_env()
    p = argparse.ArgumentParser(prog="gnomix1000g")
    p.add_argument("step", choices=[s[0] for s in STEPS] + ["all"])
    p.add_argument("--chroms", default=os.environ.get("CHROMS", "1-22"))
    p.add_argument("--batch-size", type=int, default=int(os.environ.get("BATCH_SIZE", "16")))
    args = p.parse_args(argv)
    chroms = chrom_list(args.chroms)
    for name, description, run, in_all in STEPS:
        if args.step == name or (args.step == "all" and in_all):
            log(f"==> {name}: {description}")
            run(chroms, args)


if __name__ == "__main__":
    main()
