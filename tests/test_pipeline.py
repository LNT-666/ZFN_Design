"""端到端测试：跑完整流水线，检查产物与不变量。"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from zfn_design import pipeline
from zfn_design.library import load_default_library
from zfn_design.models import SequenceRecord


@pytest.fixture
def payload(demo_record):
    return pipeline.run(
        [demo_record],
        library=load_default_library(),
        spacers=(6,),
        max_sites=5,
        n_designs=2,
        fusion_keys=["single_activator"],
        command="pytest",
    )


def test_pipeline_produces_designs(payload):
    assert payload["designs"], "demo 序列上应当至少能设计出 1 个 ZFN"
    assert payload["targets_summary"]["n_sites"] > 0


def test_design_schema(payload):
    d = payload["designs"][0]
    for key in ("design_id", "target", "monomer_left", "monomer_right", "score"):
        assert key in d
    assert len(d["monomer_left"]["fingers"]) == 3
    assert len(d["monomer_right"]["fingers"]) == 3
    assert d["monomer_left"]["protein"].startswith("M")


def test_triplet_orientation_matches_halfsite(payload):
    """左单体的三联体 (N→C) 必须是半位点 3'→5' 的三个三联体。"""
    from zfn_design.sequtil import triplets_by_finger

    for d in payload["designs"]:
        hs = d["monomer_left"]["halfsite"]
        expect = [triplets_by_finger(hs)[f] for f in ("F1", "F2", "F3")]
        assert d["monomer_left"]["triplets_N_to_C"] == expect
        assert d["monomer_left"]["fingers"][0]["triplet"] == expect[0]


def test_site_status_totals(payload):
    ts = payload["targets_summary"]
    assert sum(ts["by_status"].values()) == ts["n_sites"]


def test_gap_sites_list_missing_triplets(payload):
    gaps = [s for s in payload["sites"] if s["designability"] == "gap"]
    for s in gaps:
        assert s["missing_triplets"], "缺口位点必须报告缺哪些三联体"
        assert s["n_combos"] == 0


def test_offtarget_self_is_excluded(payload):
    for d in payload["designs"]:
        ot = d["offtarget"]
        assert ot is not None
        own = (d["target"]["seq_name"], d["target"]["start"], d["target"]["strand"])
        for h in ot["hits"]:
            assert (h["seq_name"], h["start"], h["strand"]) != own


def test_write_outputs(payload, tmp_path):
    paths = pipeline.write_outputs(payload, tmp_path)
    assert paths["json"].exists()
    assert paths["html"].exists()
    assert paths["sites_tsv"].exists()
    assert paths["designs_tsv"].exists()
    data = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert data["meta"]["version"]
    html = paths["html"].read_text(encoding="utf-8")
    assert "实验验证" in html
    assert "ZFN 设计报告" in html


def test_library_warning_surfaced():
    """种子库覆盖极低时，必须显式提示用户这是预期行为。"""
    from zfn_design.library import FingerLibrary

    seed_path = Path(__file__).resolve().parents[1] / "zfn_design" / "data" / "finger_library_seed.csv"
    rec = SequenceRecord(name="x", seq="ACGT" * 30)
    p = pipeline.run([rec], library=FingerLibrary.load_csv(seed_path, name="seed_1AAY"),
                     spacers=(6,), max_sites=3)
    assert any("64" in w for w in p["warnings"])


def test_context_fallback_is_reported():
    """B1H 数据只有 F2/F3 上下文时，F1 位退用别的上下文必须被记录。"""
    rec = SequenceRecord(name="x", seq="GCTTAGCTAGCTAGGCTTAGCTAGCTAAGCTTAG")
    p = pipeline.run([rec], library=load_default_library(), spacers=(6,), max_sites=2)
    tsum = p["targets_summary"]
    if tsum.get("n_context_fallback_sites"):
        assert any("上下文" in w for w in p["warnings"])
    else:  # 只有种子库时不应该有 fallback 标记
        assert "n_context_fallback_sites" in tsum


def test_fusion_marks_placeholders(payload):
    fusions = payload["fusions"]
    assert fusions
    f = fusions[0]
    assert "draft_protein" in f
    assert f["architecture"] == "single_activator"
