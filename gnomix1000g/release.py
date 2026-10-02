"""Package the outputs as release assets (output/release/), with SHA256SUMS.

GitHub release assets must be smaller than 2 GB each, so the corrected panel is
shipped per chromosome, with the .pgen compressed by zstd (decompress with
`plink2 --zst-decompress` or `zstd -d`, as for the official PLINK 2 files).
"""

from __future__ import annotations

import shutil
import tarfile
import zipfile

from .common import OUTPUT, atomic_path, log, sha256_file, threads

RELEASE = OUTPUT / "release"
MAX_ASSET = 2_000_000_000


def _tar(name: str, members: list, gz: bool = False) -> None:
    with atomic_path(RELEASE / name) as tmp, tarfile.open(tmp, "w:gz" if gz else "w") as tar:
        for path in members:
            tar.add(path, arcname=str(path.relative_to(OUTPUT)))


def _copy(src, name: str) -> None:
    with atomic_path(RELEASE / name) as tmp:
        shutil.copyfile(src, tmp)


def run(chroms: list[int]) -> None:
    import zstandard

    RELEASE.mkdir(parents=True, exist_ok=True)
    for k in ("final", "raw", "gnofix"):
        _copy(OUTPUT / "tracts" / f"tracts_{k}.tsv.gz", f"tracts_{k}.tsv.gz")
    _copy(OUTPUT / "tracts" / "gnofix_switches.tsv.gz", "gnofix_switches.tsv.gz")
    _tar("global_ancestry.tar.gz", [OUTPUT / f"global_ancestry_{k}.tsv" for k in ("final", "raw", "gnofix")]
         + [OUTPUT / "population_ancestry.tsv"], gz=True)
    _tar("msp.tar", sorted((OUTPUT / "msp").glob("*.msp.gz")) + sorted((OUTPUT / "windows").glob("*.tsv")))
    _tar("posteriors.tar", sorted((OUTPUT / "posteriors").glob("*.npz")))
    _tar("bed_hg38.tar.gz", sorted((OUTPUT / "tracts" / "bed_hg38").glob("*.bed")), gz=True)
    _tar("validation.tar.gz", sorted(p for p in (OUTPUT / "validation").rglob("*") if p.is_file()), gz=True)
    _copy(OUTPUT / "REPORT.md", "REPORT.md")
    with atomic_path(RELEASE / "karyograms.zip") as tmp, zipfile.ZipFile(tmp, "w", zipfile.ZIP_STORED) as z:
        for p in sorted((OUTPUT / "karyograms").glob("*/*.png")):
            z.write(p, arcname=str(p.relative_to(OUTPUT)))
    log("  packaged tracts, ancestry, msp, posteriors, BEDs, karyograms, validation")

    pfile = OUTPUT / "pfile"
    if any(pfile.glob("chr*_gnofix.pgen")):
        cctx = zstandard.ZstdCompressor(level=6, threads=threads())
        for c in chroms:
            src_pgen, dst = pfile / f"chr{c}_gnofix.pgen", RELEASE / f"chr{c}_gnofix.pgen.zst"
            if not dst.exists() or dst.stat().st_mtime < src_pgen.stat().st_mtime:
                with atomic_path(dst) as tmp, open(src_pgen, "rb") as src, open(tmp, "wb") as out:
                    cctx.copy_stream(src, out)
            _copy(pfile / f"chr{c}_gnofix.pvar.zst", f"chr{c}_gnofix.pvar.zst")
            _copy(pfile / f"chr{c}_gnofix.psam", f"chr{c}_gnofix.psam")
        _copy(pfile / "verification.json", "pfile_verification.json")
        _copy(pfile / "swaps.tsv.gz", "pfile_swaps.tsv.gz")
        log("  packaged the phase-corrected panel")

    assets = sorted(p for p in RELEASE.iterdir()
                    if p.is_file() and p.name != "SHA256SUMS" and not p.name.startswith("."))
    too_big = [p.name for p in assets if p.stat().st_size >= MAX_ASSET]
    if too_big:
        raise RuntimeError(f"release assets of 2 GB or more (GitHub limit): {too_big}")
    with atomic_path(RELEASE / "SHA256SUMS") as tmp:
        tmp.write_text("".join(f"{sha256_file(p)}  {p.name}\n" for p in assets))
    total = sum(p.stat().st_size for p in assets)
    log(f"  {len(assets)} release assets, {total / 1e9:.1f} GB, in {RELEASE}")
