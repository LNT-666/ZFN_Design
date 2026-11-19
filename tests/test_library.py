"""指头库测试：证据等级、覆盖统计、查询过滤、导入格式。"""

from __future__ import annotations

import pytest

from conftest import ZIF268_PROTEIN
from zfn_design.library import FingerLibrary, LibraryError
from zfn_design.models import EvidenceLevel


def test_seed_library_loads(seed_library):
    assert len(seed_library) == 3
    assert seed_library.triplets == {"GCG", "TGG"}


def test_seed_library_is_pdb_verifiable(seed_library):
    """种子库里每条螺旋都必须真的出现在 1AAY 蛋白序列里。"""
    for rec in seed_library:
        assert rec.source_type == "pdb"
        assert rec.source_id == "1AAY"
        assert rec.helix7 in ZIF268_PROTEIN, f"{rec.helix_id} 的螺旋不在 1AAY 序列中"


def test_seed_coverage_is_intentionally_low(seed_library):
    stats = seed_library.coverage_stats()
    assert stats["triplets_covered"] == 2
    assert stats["triplets_total"] == 64
    assert stats["coverage_fraction"] < 0.05


def test_candidates_filter_by_split_position(seed_library):
    f1 = seed_library.candidates("GCG", split_position="F1")
    f3 = seed_library.candidates("GCG", split_position="F3")
    assert [r.helix_id for r in f1] == ["1AAY_F1"]
    assert [r.helix_id for r in f3] == ["1AAY_F3"]
    assert seed_library.candidates("GCG") == sorted(
        seed_library.candidates("GCG"),
        key=lambda r: (r.evidence.tier, r.affinity_nM if r.affinity_nM is not None else 1e9, r.helix_id),
    )


def test_missing_triplet_returns_empty(seed_library):
    assert seed_library.candidates("AAA") == []


def test_max_tier_filters_inferred(tiny_library):
    assert any(r.evidence is EvidenceLevel.CODE_INFERENCE for r in tiny_library.candidates("AAA"))
    assert tiny_library.candidates("AAA", max_tier=1) == []


def test_candidates_sorted_curated_first(tiny_library):
    cands = tiny_library.candidates("GCG")
    assert cands[0].evidence.tier == 1
    assert cands[-1].evidence.tier == 2


def test_evidence_level_parse():
    assert EvidenceLevel.parse("measured") is EvidenceLevel.MEASURED
    assert EvidenceLevel.parse("CODE-INFERENCE") is EvidenceLevel.CODE_INFERENCE
    assert EvidenceLevel.parse("literature") is EvidenceLevel.LITERATURE
    with pytest.raises(ValueError):
        EvidenceLevel.parse("i-made-this-up")


def test_evidence_tiers():
    assert EvidenceLevel.MEASURED.tier == 1
    assert EvidenceLevel.LITERATURE.tier == 1
    assert EvidenceLevel.VARIANT_UNVALIDATED.tier == 2
    assert EvidenceLevel.CODE_INFERENCE.tier == 2
    assert EvidenceLevel.NONE.tier == 3


def test_helix6_and_core4():
    rec = next(iter(FingerLibrary.load_csv(
        __import__("pathlib").Path(__file__).resolve().parents[1]
        / "zfn_design" / "data" / "finger_library_seed.csv")))
    # RSDELTR: -1=R 1=S 2=D 3=E 5=T 6=R
    assert rec.helix6 == "RSDETR"
    assert rec.core4 == "RDER"


def test_invalid_row_rejected(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text(
        "helix_id,triplet,helix7,evidence_level\n"
        "x,GCG,RSDELTR,not-a-level\n",
        encoding="utf-8",
    )
    with pytest.raises(LibraryError):
        FingerLibrary.load_csv(p)


def test_bad_triplet_rejected(tmp_path):
    p = tmp_path / "bad2.csv"
    p.write_text("helix_id,triplet,helix7,evidence_level\nx,GCX,RSDELTR,measured\n", encoding="utf-8")
    with pytest.raises(LibraryError):
        FingerLibrary.load_csv(p)


def test_csv_roundtrip(tmp_path, tiny_library):
    out = tiny_library.to_csv(tmp_path / "lib.csv")
    again = FingerLibrary.load_csv(out)
    assert len(again) == len(tiny_library)
    assert {r.helix_id for r in again} == {r.helix_id for r in tiny_library}


def test_merge_prefers_first(seed_library, tiny_library):
    merged = seed_library.merged(tiny_library)
    assert len(merged) == len(seed_library) + len(tiny_library)


def test_import_b1h_dir(tmp_path):
    """构造一个符合 B1H 格式的目录，验证导入器。"""
    d = tmp_path / "b1h"
    d.mkdir()
    (d / "GCG_F2.txt").write_text(
        "# helix(-1,1,2,3,5,6)  freq\n"
        "R S D E T R   0.91\n"
        "R S D E T K   0.40\n"
        "Q S S N L R   0.12\n",
        encoding="utf-8",
    )
    lib = FingerLibrary.import_b1h_dir(d, top_n=2)
    assert len(lib) == 2
    recs = lib.candidates("GCG")
    assert recs[0].helix7[0] == "R"
    assert all(r.evidence is EvidenceLevel.MEASURED for r in recs)
    assert all(r.split_position == "F2" for r in recs)
    # 位置 4 数据里没有（不接触 DNA），默认补 L，并在 notes 里留证
    assert recs[0].helix7[4] == "L"
    assert len(recs[0].helix7) == 7
    assert "第 4 位未测" in recs[0].notes
    # 名次写在 id 前面，便于跨上下文按名次交错排序
    assert recs[0].helix_id == "B1H_r01_F2_GCG"


def test_import_b1h_official_format(tmp_path):
    """Princeton 官方压缩包格式：6 残基连写 + 空格 + 频率，文件名带后缀。"""
    d = tmp_path / "F3_union"
    d.mkdir()
    (d / "TGG_amino_acid_union_cut.txt").write_text(
        "RRDTTN 0.0069884\nRKSTTN 0.0048866\n", encoding="utf-8"
    )
    lib = FingerLibrary.import_b1h_dir(d, top_n=5, split_position="F3", pos4="L")
    recs = lib.candidates("TGG", split_position="F3")
    assert len(recs) == 2
    assert recs[0].helix7 == "RRDT" + "L" + "TN"
    assert "源 6-mer=RRDTTN" in recs[0].notes
    assert lib.candidates("TGG", split_position="F2") == []
