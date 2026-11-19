"""靶点枚举测试：数量、坐标、方向、可重建性。"""

from __future__ import annotations

import pytest

from conftest import DEMO_SITE
from zfn_design.models import SequenceRecord
from zfn_design.sequtil import revcomp
from zfn_design.targets import expected_site_count, scan_sequence


def test_site_count_matches_formula():
    seq = "ACGT" * 40          # 160 bp
    rec = SequenceRecord(name="t", seq=seq)
    sites = scan_sequence(rec, spacers=(5, 6, 7), both_strands=True)
    assert len(sites) == expected_site_count(len(seq), spacers=(5, 6, 7), both_strands=True)


def test_single_strand_count():
    seq = "ACGT" * 40
    rec = SequenceRecord(name="t", seq=seq)
    sites = scan_sequence(rec, spacers=(6,), both_strands=False)
    assert len(sites) == expected_site_count(len(seq), spacers=(6,), both_strands=False)


def test_demo_site_is_recovered(demo_record):
    sites = scan_sequence(demo_record, spacers=(6,), both_strands=False)
    target = [s for s in sites if s.left_half == "GCGTGGGCG" and s.right_half == "GCGTGGGCG"]
    assert target, "应当能找到左=右=GCGTGGGCG 的那个位点"
    s = target[0]
    assert s.spacer == "GATATC"
    assert s.spacer_len == 6


def test_site_sequence_reconstructs(demo_record):
    """site_seq 必须恒等于 left + spacer + revcomp(right)。"""
    for s in scan_sequence(demo_record, spacers=(5, 6, 7), both_strands=True):
        assert s.site_seq == s.left_half + s.spacer + revcomp(s.right_half)


def test_plus_strand_coordinates(demo_record):
    sites = scan_sequence(demo_record, spacers=(6,), both_strands=False)
    s = [x for x in sites if x.left_half == "GCGTGGGCG"][0]
    expected_start = demo_record.seq.index(DEMO_SITE) + 1
    assert s.start == expected_start
    assert s.strand == "+"
    assert demo_record.seq[s.start - 1 : s.end] == s.site_seq


def test_minus_strand_coordinates_are_on_input(demo_record):
    """负链位点的 start/end 必须仍然落在输入序列坐标上。"""
    for s in scan_sequence(demo_record, spacers=(6,), both_strands=True):
        if s.strand != "-":
            continue
        region = demo_record.seq[s.start - 1 : s.end]
        assert revcomp(region) == s.site_seq
        assert len(region) == len(s.site_seq)


def test_cut_site_inside_spacer(demo_record):
    for s in scan_sequence(demo_record, spacers=(5, 6, 7), both_strands=True):
        assert s.start <= s.cut_site_input <= s.end


def test_n_windows_are_skipped():
    seq = "ACGT" * 20 + "N" * 30 + "ACGT" * 20
    rec = SequenceRecord(name="n", seq=seq)
    sites = scan_sequence(rec, spacers=(6,), both_strands=True)
    for s in sites:
        region = rec.seq[s.start - 1 : s.end]
        assert "N" not in region
        assert "N" not in s.site_seq


def test_short_sequence_yields_nothing():
    rec = SequenceRecord(name="s", seq="ACGTACGT")
    assert scan_sequence(rec, spacers=(6,), both_strands=True) == []
