"""Prepare inputs: unpack models, index the PLINK 2 panel, match model SNPs.

Outputs in work/:
  models/chrN/model_chm_N.pkl   pretrained model (from the tarball)
  models/chrN/meta.npz          model SNPs (GRCh37), window layout, genetic map
  samples.tsv                   the 3,202 samples in .pgen order
  pvar/chrN.npz                 biallelic SNV records of the panel: pgen index, pos, ref, alt
  all_hg38.pgen                 decompressed genotype file
  match/chrN.npz                model SNP -> panel record, allele flip, hg38 position
  match/summary.tsv             SNP coverage per chromosome
"""

from __future__ import annotations

import io
import pickle
import tarfile

import numpy as np
import pandas as pd

from .common import DOWNLOADS, ANCESTRIES, WORK, atomic_path, import_gnomix, log, read_json, write_json
from .liftover import Chain, complement


# ---- models ------------------------------------------------------------------

def model_dir(chrom: int):
    return WORK / "models" / f"chr{chrom}"


def model_path(chrom: int):
    return model_dir(chrom) / f"model_chm_{chrom}.pkl"


def unpack_models(chroms: list[int]) -> None:
    todo = [c for c in chroms if not model_path(c).exists()]
    if not todo:
        log("  [skip] models unpacked")
        return
    wanted = {f"pretrained_gnomix_models/chr{c}/model_chm_{c}.pkl": c for c in todo}
    log(f"  unpacking models for chr{','.join(map(str, todo))}")
    with tarfile.open(DOWNLOADS / "pretrained_gnomix_models.tar.gz", "r:gz") as tar:
        for member in tar:
            c = wanted.get(member.name)
            if c is None:
                continue
            dest = model_path(c)
            dest.parent.mkdir(parents=True, exist_ok=True)
            part = dest.with_suffix(".pkl.part")
            with tar.extractfile(member) as src, open(part, "wb") as dst:
                while True:
                    buf = src.read(1 << 24)
                    if not buf:
                        break
                    dst.write(buf)
            part.replace(dest)
    missing = [c for c in todo if not model_path(c).exists()]
    if missing:
        raise RuntimeError(f"models missing from tarball: {missing}")


def _quiet_xgboost_teardown() -> None:
    """xgboost 1.1.1 boosters freed during interpreter teardown raise
    "'NoneType' object has no attribute 'XGBoosterFree'" (harmless, but it
    prints a traceback). Skip the free call once the library is gone."""
    import xgboost.core as xc
    if getattr(xc.Booster.__del__, "_gnomix1000g", False):
        return
    original = xc.Booster.__del__

    def __del__(self):
        if getattr(xc, "_LIB", None) is not None:
            original(self)

    __del__._gnomix1000g = True
    xc.Booster.__del__ = __del__


def load_model(chrom: int):
    """Load a pretrained Gnomix model and make it single-threaded."""
    _quiet_xgboost_teardown()
    import_gnomix()
    with open(model_path(chrom), "rb") as fh:
        model = pickle.load(fh)
    # Pickled with n_jobs=None and base_multithread=True: that starts a process
    # pool with one worker per core inside every call. We parallelise ourselves.
    model.base.base_multithread = False
    model.base.log_inference = False
    model.base.n_jobs = 1
    model.smooth.n_jobs = 1
    model.smooth.model.n_jobs = 1
    model.smooth.model.get_booster().set_param("nthread", 1)
    return model


def write_model_meta(chroms: list[int]) -> None:
    for c in chroms:
        out = model_dir(c) / "meta.npz"
        if out.exists():
            continue
        m = load_model(c)
        if list(m.population_order) != ANCESTRIES:
            raise RuntimeError(f"chr{c}: unexpected population order {list(m.population_order)}")
        gm = m.gen_map_df
        with atomic_path(out) as tmp:
            np.savez(
                tmp,
                snp_pos=np.asarray(m.snp_pos, dtype=np.int64),
                snp_ref=np.asarray(m.snp_ref).astype("U1"),
                snp_alt=np.asarray(m.snp_alt).astype("U1"),
                C=m.C, M=m.M, W=m.W, A=m.A, S=m.smooth.S, context=m.context,
                populations=np.asarray(m.population_order).astype("U3"),
                gm_pos=np.asarray(gm.pos, dtype=np.float64),
                gm_cm=np.asarray(gm.pos_cm, dtype=np.float64),
            )
        log(f"  chr{c}: model C={m.C} SNPs, window M={m.M} SNPs, W={m.W} windows, smoother S={m.smooth.S}")
        del m, gm


def model_meta(chrom: int) -> dict:
    with np.load(model_dir(chrom) / "meta.npz") as z:
        return {k: z[k] for k in z.files}


# ---- panel -------------------------------------------------------------------

def write_samples() -> pd.DataFrame:
    out = WORK / "samples.tsv"
    if not out.exists():
        psam = pd.read_csv(DOWNLOADS / "all_hg38.psam", sep="\t")
        psam = psam.rename(columns={"#IID": "IID"})
        ids = set(psam.IID)
        psam["pgen_index"] = np.arange(len(psam))
        has_pat, has_mat = psam.PAT.isin(ids), psam.MAT.isin(ids)
        parents = set(psam.PAT[has_pat]) | set(psam.MAT[has_mat])
        psam["trio_child"] = has_pat & has_mat
        psam["has_parent_in_panel"] = has_pat | has_mat
        psam["has_child_in_panel"] = psam.IID.isin(parents)
        panel = pd.read_csv(DOWNLOADS / "integrated_call_samples_v3.20130502.ALL.panel", sep="\t")
        psam["phase3_unrelated_2504"] = psam.IID.isin(set(panel["sample"]))
        with atomic_path(out) as tmp:
            psam.to_csv(tmp, sep="\t", index=False)
    return pd.read_csv(out, sep="\t")


def split_pvar() -> None:
    """Stream all_hg38.pvar.zst once and keep biallelic SNV records per autosome."""
    out = WORK / "pvar"
    if (out / "summary.json").exists() and (out / "allele_idx_offsets.npy").exists():
        log("  [skip] pvar indexed")
        return
    import zstandard

    out.mkdir(parents=True, exist_ok=True)
    keep = {str(c): {"idx": [], "pos": [], "ref": [], "alt": []} for c in range(1, 23)}
    allele_cts = []  # pgenlib needs allele offsets because some records are multiallelic
    offset = 0
    with open(DOWNLOADS / "all_hg38.pvar.zst", "rb") as fh:
        reader = zstandard.ZstdDecompressor().stream_reader(fh)
        text = io.TextIOWrapper(reader, encoding="ascii")
        while True:
            line = text.readline()
            if not line:
                raise RuntimeError("no #CHROM header in pvar")
            if line.startswith("#CHROM"):
                cols = line[1:].rstrip("\n").split("\t")
                break
        chunks = pd.read_csv(text, sep="\t", header=None, names=cols, usecols=["CHROM", "POS", "REF", "ALT"],
                             dtype={"CHROM": str, "POS": np.int64, "REF": str, "ALT": str},
                             chunksize=2_000_000, na_filter=False)
        for chunk in chunks:
            n = len(chunk)
            idx = np.arange(offset, offset + n, dtype=np.int64)
            offset += n
            allele_cts.append((2 + chunk.ALT.str.count(",").values).astype(np.uint8))
            snv = (chunk.REF.str.len().values == 1) & (chunk.ALT.str.len().values == 1)
            chrom = chunk.CHROM.values
            for c, store in keep.items():
                m = snv & (chrom == c)
                if m.any():
                    store["idx"].append(idx[m])
                    store["pos"].append(chunk.POS.values[m])
                    store["ref"].append(chunk.REF.values[m].astype("U1"))
                    store["alt"].append(chunk.ALT.values[m].astype("U1"))
            log(f"    pvar records read: {offset:,}")
    counts = {}
    for c, store in keep.items():
        arrays = {k: np.concatenate(v) for k, v in store.items()}
        np.savez(out / f"chr{c}.npz", **arrays)
        counts[c] = int(len(arrays["idx"]))
    cts = np.concatenate(allele_cts)
    offsets = np.zeros(len(cts) + 1, dtype=np.uintp)
    np.cumsum(cts, out=offsets[1:])
    np.save(out / "allele_idx_offsets.npy", offsets)
    write_json(out / "summary.json", {"variant_ct": offset, "multiallelic_records": int((cts > 2).sum()),
                                      "biallelic_snv_per_autosome": counts})


def variant_count() -> int:
    return int(read_json(WORK / "pvar" / "summary.json")["variant_ct"])


def allele_offsets() -> np.ndarray:
    """Allele index offsets for every panel record (memory-mapped, shared between workers)."""
    return np.load(WORK / "pvar" / "allele_idx_offsets.npy", mmap_mode="r")


def decompress_pgen() -> None:
    dest = WORK / "all_hg38.pgen"
    if dest.exists():
        log("  [skip] pgen decompressed")
        return
    import zstandard

    log("  decompressing all_hg38.pgen.zst (about 9.5 GB)")
    part = dest.with_suffix(".pgen.part")
    with open(DOWNLOADS / "all_hg38.pgen.zst", "rb") as src, open(part, "wb") as dst:
        zstandard.ZstdDecompressor().copy_stream(src, dst, read_size=1 << 24, write_size=1 << 24)
    part.replace(dest)


# ---- model SNP matching --------------------------------------------------------

def match_snps(model_idx: np.ndarray, pos38: np.ndarray, ref: np.ndarray, alt: np.ndarray,
               panel: pd.DataFrame) -> pd.DataFrame:
    """Match model SNPs (alleles on the GRCh38 + strand) to panel SNV records.

    A model SNP matches a record with the same position and the same alleles,
    possibly with REF and ALT exchanged (flip). Model SNPs that match several
    records, and records matched by several model SNPs, are dropped.
    Returns columns i (model index), pidx (panel record), flip, sorted by i.
    """
    mod = pd.DataFrame({"i": model_idx, "pos": pos38, "mref": ref, "malt": alt})
    j = mod.merge(panel, on="pos", how="inner")
    same = (j.mref == j.pref) & (j.malt == j.palt)
    swap = (j.mref == j.palt) & (j.malt == j.pref)
    j = j[same | swap].assign(flip=swap[same | swap].values)
    j = j[~j.i.duplicated(keep=False) & ~j.pidx.duplicated(keep=False)]
    return j.sort_values("i")[["i", "pidx", "flip"]]


def match_chrom(chrom: int, chain: Chain) -> dict:
    meta = model_meta(chrom)
    pos19 = meta["snp_pos"]
    ref, alt = meta["snp_ref"].astype(str), meta["snp_alt"].astype(str)
    dst_chrom, pos38, minus = chain.lift(f"chr{chrom}", pos19)
    same_chrom = dst_chrom == f"chr{chrom}"
    pos38 = np.where(same_chrom, pos38, -1)
    mref = np.where(minus, complement(ref), ref)
    malt = np.where(minus, complement(alt), alt)

    with np.load(WORK / "pvar" / f"chr{chrom}.npz") as z:
        panel = pd.DataFrame({"pidx": z["idx"], "pos": z["pos"], "pref": z["ref"], "palt": z["alt"]})
    lifted = np.nonzero(same_chrom)[0]
    j = match_snps(lifted, pos38[lifted], mref[lifted], malt[lifted], panel)
    out = WORK / "match"
    out.mkdir(parents=True, exist_ok=True)
    with atomic_path(out / f"chr{chrom}.npz") as tmp:
        np.savez(tmp, model_idx=j.i.values.astype(np.int64), pgen_idx=j.pidx.values.astype(np.int64),
                 flip=j.flip.values.astype(bool), pos38=pos38.astype(np.int64))
    stats = {
        "chrom": chrom,
        "model_snps": int(len(pos19)),
        "lifted_same_chrom": int(same_chrom.sum()),
        "lifted_minus_strand": int((minus & same_chrom).sum()),
        "matched": int(len(j)),
        "matched_allele_swapped": int(j.flip.sum()),
        "fraction_matched": round(len(j) / len(pos19), 5),
    }
    return stats


def match_all(chroms: list[int]) -> pd.DataFrame:
    out = WORK / "match" / "summary.tsv"
    done = pd.read_csv(out, sep="\t") if out.exists() else pd.DataFrame()
    have = set(done.chrom) if len(done) else set()
    todo = [c for c in chroms if c not in have or not (WORK / "match" / f"chr{c}.npz").exists()]
    if todo:
        chain = Chain(DOWNLOADS / "hg19ToHg38.over.chain.gz")
        rows = []
        for c in todo:
            s = match_chrom(c, chain)
            log(f"  chr{c}: {s['matched']:,}/{s['model_snps']:,} model SNPs in panel "
                f"({100 * s['fraction_matched']:.2f}%), {s['matched_allele_swapped']:,} with REF/ALT swapped")
            rows.append(s)
        done = pd.concat([done[~done.chrom.isin(todo)] if len(done) else done, pd.DataFrame(rows)])
        done = done.sort_values("chrom")
        with atomic_path(out) as tmp:
            done.to_csv(tmp, sep="\t", index=False)
    low = done[done.fraction_matched < 0.8]
    if len(low):
        raise RuntimeError(f"model SNP coverage below 80% (Gnomix minimum): {low.to_dict('records')}")
    return done
