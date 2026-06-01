"""
Basit testler. Calistirmak icin:
    PYTHONPATH=src python -m pytest tests/      (pytest varsa)
veya:
    PYTHONPATH=src python tests/test_basic.py   (pytest yoksa)
"""

import fakhermsa as msa


def test_identical_sequences_have_max_similarity():
    s = "MKQLEDKVEEL"
    assert msa.fft_similarity(s, s) > 0.99


def test_alignment_columns_equal_length():
    seqs = ["MKQLED", "MKKLED", "MTQLEDR"]
    aligned = msa.align(seqs)
    lengths = {len(r) for r in aligned}
    assert len(lengths) == 1   # hepsi esit uzunlukta olmali


def test_alignment_preserves_residues():
    # Bosluklar cikarilinca orijinal dizi geri gelmeli
    seqs = ["MKQLED", "MKKLED", "MTQLEDR"]
    aligned = msa.align(seqs)
    for original, row in zip(seqs, aligned):
        assert row.replace("-", "") == original


def test_single_sequence():
    assert msa.align(["MKQL"]) == ["MKQL"]


def test_score_increases_for_similar():
    similar = msa.sum_of_pairs_score(msa.align(["MKQLED", "MKQLED", "MKQLED"]))
    assert similar > 0


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    for t in tests:
        t()
        print(f"[OK] {t.__name__}")
        passed += 1
    print(f"\n{passed}/{len(tests)} test gecti.")
