"""序列工具与方向约定测试。

这一组测试守的是全项目最容易出错的地方：**指头 N→C 与 DNA 5'→3' 反平行**。
如果这里错了，所有设计都会"看起来对、实际反了"。
"""

from __future__ import annotations

import pytest

from zfn_design.sequtil import (
    SequenceError,
    gc_content,
    hamming,
    finger_order,
    is_self_complementary,
    max_homopolymer,
    revcomp,
    tandem_repeat_count,
    triplets_by_finger,
)


def test_revcomp_basic():
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("AAAA") == "TTTT"
    assert revcomp("GCGTGGGCG") == "CGCCCACGC"


def test_revcomp_is_involution():
    s = "ACGTTGCAACCGGTTA"
    assert revcomp(revcomp(s)) == s


def test_triplets_by_finger_is_antiparallel():
    """半位点 ABC-DEF-GHI：F1(N 端) 读 GHI，F3(C 端) 读 ABC。"""
    m = triplets_by_finger("AAACCCGGG")
    assert m["F1"] == "GGG"
    assert m["F2"] == "CCC"
    assert m["F3"] == "AAA"


def test_triplets_by_finger_matches_zif268():
    """1AAY: F1→GCG, F2→TGG, F3→GCG（外侧两个同为 GCG，方向无关）。"""
    m = triplets_by_finger("GCGTGGGCG")
    assert m["F1"] == "GCG"
    assert m["F2"] == "TGG"
    assert m["F3"] == "GCG"


def test_triplets_rejects_wrong_length():
    with pytest.raises(SequenceError):
        triplets_by_finger("GCGTGGA")   # 7 bp，不是 3 的倍数


def test_triplets_by_finger_eighteen_bp():
    """六指半位点（18 bp）：F1–F3 读 3' 那 9 bp，F4–F6 读 5' 那 9 bp。"""
    m = triplets_by_finger("GCGTGGGCGAAATTTCCC")
    assert [m[f"F{i}"] for i in range(1, 7)] == [
        "CCC", "TTT", "AAA", "GCG", "TGG", "GCG",
    ]
    assert finger_order("GCGTGGGCGAAATTTCCC") == ["CCC", "TTT", "AAA", "GCG", "TGG", "GCG"]


def test_gc_content():
    assert gc_content("GCGC") == 1.0
    assert gc_content("ATAT") == 0.0
    assert gc_content("ACGT") == 0.5
    assert gc_content("") == 0.0


def test_hamming():
    assert hamming("ACGT", "ACGT") == 0
    assert hamming("ACGT", "ACTT") == 1
    with pytest.raises(ValueError):
        hamming("ACG", "AC")


def test_max_homopolymer():
    assert max_homopolymer("AAACGT") == 3
    assert max_homopolymer("AAAAAAA") == 7
    assert max_homopolymer("") == 0


def test_tandem_repeat_count():
    assert tandem_repeat_count("AAAAAAAA", 2) >= 1
    assert tandem_repeat_count("ACACAC", 2) >= 1
    assert tandem_repeat_count("ACGTACGT", 2) == 0   # 每个 2-mer 只出现 2 次，未达阈值


def test_is_self_complementary():
    # AAAAAA 与 TTTTTT 互为反向互补
    assert is_self_complementary("AAAAAATTTTTT", min_arm=6)
    assert not is_self_complementary("ACGATCGATC", min_arm=6)
