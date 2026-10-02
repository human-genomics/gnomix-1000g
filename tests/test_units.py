"""Unit tests that need no downloaded data."""

from __future__ import annotations

import gzip
import importlib
import json

import numpy as np
import pandas as pd
import pytest

from gnomix1000g.common import ANCESTRIES, REFERENCE, chrom_list, group_of
from gnomix1000g.infer import window_of_snp
from gnomix1000g.liftover import Chain, complement
from gnomix1000g.pfile import apply_swaps, nearest_window
from gnomix1000g.prepare import match_snps
from gnomix1000g.tracts import (haplotype_tracts, hg38_discordant, hg38_span, label_runs, longest_increasing,
                                policy_applies, window_lengths_cm)

MODULES = ["common", "liftover", "prepare", "haplotypes", "checks", "infer", "trio", "tracts", "karyogram",
           "compare", "pfile", "figures", "report", "release", "__main__"]


@pytest.mark.parametrize("name", MODULES)
def test_modules_import(name):
    importlib.import_module(f"gnomix1000g.{name}")


def test_chrom_list():
    assert chrom_list("1-3 22") == [1, 2, 3, 22]
    assert chrom_list("22,21") == [21, 22]
    with pytest.raises(ValueError):
        chrom_list("23")


def test_groups():
    assert group_of("ASW", "AFR") == "AFR-American"
    assert group_of("PEL", "AMR") == "AMR"
    assert group_of("YRI", "AFR") == "AFR"


# ---- liftover and matching ------------------------------------------------------

def write_chain(path):
    # chr1 (size 1000) -> chrA (+): blocks [100,200)->[0,100), gap 50 / 10, [250,300)->[110,160)
    # chr2 (size 1000) -> chrB (-, size 500): block [0,100) -> minus-strand [0,100)
    text = ("chain 1000 chr1 1000 + 100 300 chrA 400 + 0 160 1\n100\t50\t10\n50\n\n"
            "chain 1000 chr2 1000 + 0 100 chrB 500 - 0 100 2\n100\n\n")
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def test_liftover(tmp_path):
    p = tmp_path / "t.chain.gz"
    write_chain(p)
    ch = Chain(p)
    c, pos, minus = ch.lift("chr1", [100, 101, 200, 201, 251, 300, 301])
    assert list(pos) == [-1, 1, 100, -1, 111, 160, -1]
    assert list(c[[1, 2]]) == ["chrA", "chrA"]
    c, pos, minus = ch.lift("chr2", [1, 100])
    assert list(pos) == [500, 401] and minus.all()  # minus strand: size - 1 - offset
    assert list(complement(np.array(["A", "C", "G", "T"]))) == ["T", "G", "C", "A"]


def test_match_snps():
    panel = pd.DataFrame({"pidx": [10, 11, 12, 13, 14], "pos": [100, 200, 300, 300, 400],
                          "pref": ["A", "C", "G", "G", "T"], "palt": ["G", "T", "A", "C", "C"]})
    j = match_snps(np.array([0, 1, 2, 3, 4, 5]), np.array([100, 200, 300, 400, 400, 999]),
                   np.array(["A", "T", "G", "T", "C", "A"]), np.array(["G", "C", "C", "C", "T", "G"]), panel)
    # 0: same alleles; 1: REF/ALT exchanged (flip); 2: the second record at 300;
    # 3 and 4 match the same record (both dropped); 5: no record at 999.
    assert list(j.i) == [0, 1, 2]
    assert list(j.pidx) == [10, 11, 13]
    assert list(j.flip) == [False, True, False]


def test_read_haplotypes_flip_and_fill(tmp_path, monkeypatch):
    pgenlib = pytest.importorskip("pgenlib")
    from gnomix1000g import haplotypes, prepare

    # 3 samples, 3 biallelic records; sample 1 is missing at record 2.
    alleles = np.array([[0, 1, 1, 1, 0, 0],
                        [1, 0, 0, 0, 1, 1],
                        [0, 0, -9, -9, 1, 0]], dtype=np.int32)
    w = pgenlib.PgenWriter(str(tmp_path / "all_hg38.pgen").encode(), 3, 3, False, 2, True)
    w.append_alleles_batch(alleles[:2], all_phased=True)
    w.append_partially_phased_batch(alleles[2:], np.array([[1, 0, 1]], dtype=np.uint8), np.array([2], np.uint32))
    w.close()
    (tmp_path / "pvar").mkdir()
    np.save(tmp_path / "pvar" / "allele_idx_offsets.npy", np.array([0, 2, 4, 6], dtype=np.uintp))
    monkeypatch.setattr(haplotypes, "WORK", tmp_path)
    monkeypatch.setattr(prepare, "WORK", tmp_path)
    # model has 4 SNPs; SNP 3 is absent from the panel; record 1 is allele-swapped
    match = {"model_idx": np.array([0, 1, 2]), "pgen_idx": np.array([0, 1, 2]),
             "flip": np.array([False, True, False])}
    X = haplotypes.read_haplotypes(1, np.array([0, 1, 2]), 4, match, 3, 3)
    assert X.shape == (6, 4)
    assert list(X[:, 0]) == [0, 1, 1, 1, 0, 0]
    assert list(X[:, 1]) == [0, 1, 1, 1, 0, 0]  # flipped
    assert list(X[:, 2]) == [0, 0, 0, 0, 1, 0]  # missing -> REF
    assert (X[:, 3] == 0).all()  # absent -> REF


# ---- windows and tracts -----------------------------------------------------------

def test_window_of_snp_remainder():
    assert list(window_of_snp(C=10, M=3, W=3)) == [0, 0, 0, 1, 1, 1, 2, 2, 2, 2]


def toy_windows(W=6, discordant=()):
    df = pd.DataFrame({
        "chrom": 1, "window": np.arange(W),
        "spos_hg19": np.arange(W) * 100 + 1, "epos_hg19": np.arange(W) * 100 + 90,
        "sgpos": np.arange(W) * 1.0, "egpos": np.arange(W) * 1.0 + 0.9,
        "spos_hg38": np.arange(W) * 100 + 1001, "epos_hg38": np.arange(W) * 100 + 1090,
        "hg38_discordant": False})
    df.loc[list(discordant), "hg38_discordant"] = True
    return df


def test_label_runs():
    rows, s, e = label_runs(np.array([[0, 0, 3, 3, 3, 0], [2, 2, 2, 2, 2, 2]]))
    assert list(rows) == [0, 0, 0, 1] and list(s) == [0, 2, 5, 0] and list(e) == [1, 4, 5, 5]


def test_haplotype_tracts():
    win = toy_windows()
    labels = np.array([[0, 0, 3, 3, 3, 0], [2, 2, 2, 2, 2, 2]], dtype=np.uint8)
    post = np.array([[1, 1, 1, 1, 1, 1], [.8, .8, .8, .8, .8, .8]])
    t = haplotype_tracts(labels, post, win, ["S1"], 1)
    assert list(t.ancestry) == ["EUR", "AFR", "EUR", "NAT"]
    assert list(t.haplotype) == ["A", "A", "A", "B"]
    assert list(t.n_windows) == [2, 3, 1, 6]
    assert list(t.start_hg19) == [1, 201, 501, 1] and list(t.end_hg19) == [190, 490, 590, 590]
    assert list(t.start_hg38) == [1001, 1201, 1501, 1001]
    assert np.allclose(t.mean_posterior, [1, 1, 1, 0.8])
    unk = haplotype_tracts(np.where(labels == 3, len(ANCESTRIES), labels), post, win, ["S1"], 1,
                           names=ANCESTRIES + ["UNK"])
    assert list(unk.ancestry) == ["EUR", "UNK", "EUR", "NAT"]


def test_hg38_span_skips_discordant_windows():
    win = toy_windows(discordant=(0, 3))
    s38, e38 = hg38_span(win, np.array([0, 3, 3]), np.array([2, 3, 5]))
    assert list(s38) == [1101, -1, 1401]  # run 0-2 starts at window 1; run 3-3 has none
    assert list(e38) == [1290, -1, 1590]


def test_longest_increasing():
    x = np.array([10, 20, 900, 30, 40, 5, 50])
    assert list(x[longest_increasing(x)]) == [10, 20, 30, 40, 50]
    assert len(longest_increasing(np.array([]))) == 0


def test_hg38_discordant():
    first = np.array([100, 300, 9000, 500, -1, 900, 1300, 1500])
    last = np.array([200, 400, 9100, 450, -1, 1000, 5_000_000, 1600])
    span19 = np.full(8, 100)
    # 2: moved far ahead (out of order); 3: reversed; 4: not lifted; 6: stretched
    assert list(np.nonzero(hg38_discordant(first, last, span19))[0]) == [2, 3, 4, 6]


def test_window_lengths():
    assert np.allclose(window_lengths_cm(toy_windows(3)), [1.0, 1.0, 0.9])


# ---- Gnofix policy and the corrected panel ------------------------------------------

def test_policy_applies_population_then_group():
    samples = pd.DataFrame({"Population": ["PUR", "MXL", "YRI", "CHB", "ACB"],
                            "SuperPop": ["AMR", "AMR", "AFR", "EAS", "AFR"]})
    policy = {"apply_gnofix": {"AMR": True, "AFR": False, "EAS": False, "AFR-American": True},
              "apply_gnofix_population": {"PUR": False, "MXL": True}}
    assert list(policy_applies(policy, samples)) == [False, True, False, False, True]


def test_reference_policy_and_summary_match_the_code():
    policy = json.loads((REFERENCE / "gnofix_policy.json").read_text())
    assert {"rule", "apply_gnofix", "apply_gnofix_population", "chromosomes"} <= set(policy)
    summary = pd.read_csv(REFERENCE / "trio_benchmark_summary.tsv", sep="\t")
    report_columns = ["level", "group", "children", "het_ancestry_fraction", "hap_err_rephased", "hap_err_gnofix",
                      "hap_err_gnofix_on_truth", "hap_err_continental_rephased", "hap_err_continental_gnofix",
                      "phase_acc_10cM_hetanc_rephased", "phase_acc_10cM_hetanc_gnofix", "SER_rephased",
                      "SER_gnofix", "gnofix"]
    assert set(report_columns) <= set(summary.columns)


def test_nearest_window_with_reversed_liftover():
    # windows 0,1,2 in model order; liftover reverses window 1 to the far end
    pos38 = np.array([100, 110, 900, 910, 300, 310, -1])
    win = np.array([0, 0, 1, 1, 2, 2, 2])
    got = nearest_window(pos38, win, np.array([90, 105, 200, 250, 305, 800, 905, 2000]))
    assert list(got) == [0, 0, 0, 2, 2, 1, 1, 1]


def test_apply_swaps():
    a = np.array([[0, 1, 1, 0], [0, 1, 1, 1]], dtype=np.int32)  # 2 variants, 2 samples
    ph = np.array([[1, 1], [1, 0]], dtype=np.uint8)
    states = np.array([[True], [True]])  # sample 1 swapped at both variants
    apply_swaps(a, ph, states, np.array([1]))
    assert a.tolist() == [[0, 1, 0, 1], [0, 1, 1, 1]]  # second variant unphased: unchanged


# ---- trio benchmark VCF I/O -----------------------------------------------------------

def test_vcf_rows_roundtrip(tmp_path):
    from gnomix1000g.trio import PHASED_CELLS, UNPHASED_CELLS, _rows, read_phased
    rng = np.random.default_rng(0)
    a = rng.integers(0, 2, size=(50, 6))  # 50 variants, 3 samples
    codes = a[:, 0::2] * 2 + a[:, 1::2]
    fixed = [f"1\t{i + 1}\t.\tA\tG\t.\tPASS\t.\tGT\t".encode() for i in range(50)]
    body = _rows(fixed, codes, PHASED_CELLS)
    assert body.count(b"\n") == 50
    first = body.split(b"\n")[0].split(b"\t")[9:]
    assert first == [f"{a[0, 2 * k]}|{a[0, 2 * k + 1]}".encode() for k in range(3)]
    unph = _rows(fixed, codes, UNPHASED_CELLS).split(b"\n")[0].split(b"\t")[9:]
    assert all(c in (b"0/0", b"0/1", b"1/1") for c in unph)
    vcf = tmp_path / "x.vcf.gz"
    with gzip.open(vcf, "wb") as fh:
        fh.write(b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQ_a\tQ_b\tQ_c\n")
        fh.write(body)
    model_idx = np.arange(0, 100, 2)  # 50 of 100 model SNPs present
    X = read_phased(str(vcf), model_idx, 100, ["a", "b", "c"])
    assert X.shape == (6, 100)
    assert np.array_equal(X[:, model_idx], a.T)
    assert (X[:, 1::2] == 0).all()  # absent SNPs get REF


# ---- inference batches and drawing helpers ------------------------------------------------

def test_assemble_joins_batches_in_panel_order(tmp_path, monkeypatch):
    from gnomix1000g import infer
    monkeypatch.setattr(infer, "WORK", tmp_path)
    bdir = tmp_path / "infer" / "chr1"
    bdir.mkdir(parents=True)
    idx = np.arange(5)
    for b, s in enumerate((idx[:2], idx[2:4], idx[4:])):
        n = len(s)
        np.savez(bdir / f"batch_{b:05d}.npz", sample_idx=s, lab_raw=np.full((2 * n, 3), b, np.uint8),
                 prob_raw=np.zeros((2 * n, 3, 8), np.uint8), lab_fix=np.zeros((2 * n, 3), np.uint8),
                 prob_fix=np.zeros((2 * n, 3, 8), np.uint8), swap=np.zeros((n, 3), bool))
    infer.assemble(1, 2, idx)
    with np.load(tmp_path / "infer" / "chr1.npz") as z:
        assert list(z["sample_idx"]) == list(idx)
        assert list(z["lab_raw"][:, 0]) == [0, 0, 0, 0, 1, 1, 1, 1, 2, 2]
    assert not bdir.exists()


def test_opacity_scale():
    from gnomix1000g.karyogram import opacity
    a = opacity(np.array([0.0, 0.35, 0.6, 0.8, 0.9, 1.0]))
    assert a[0] == a[1] == 0.5 and a[-1] == a[-2] == 1.0
    assert 0.7 <= a[2] <= 0.75 and a[3] >= 0.9
    assert (np.diff(opacity(np.linspace(0, 1, 101))) >= 0).all()


def test_karyogram_runs_merge_and_split():
    from gnomix1000g.karyogram import runs
    cat = np.array([0, 0, 1, 1, 1])
    alpha = np.array([1.0, 1.0, 1.0, 0.5, 0.5])
    s = np.array([0.0, 10, 20, 30, 40])
    e = np.array([10.0, 20, 30, 40, 50])
    ok = np.array([True, True, True, True, False])
    assert runs(cat, alpha, ok, s, e) == [(0.0, 20.0, 0, 1.0), (20.0, 30.0, 1, 1.0), (30.0, 40.0, 1, 0.5)]


def test_category_runs():
    from gnomix1000g.figures import _category_runs
    to_cat = np.array([0, 1, 1, 2, 0, 0, 0, 0])  # ancestries 1 and 2 shown as one category
    assert _category_runs(np.array([[1, 2, 1, 3], [0, 0, 0, 0]]), to_cat) == 2 + 1
