"""S1（DesignScore / all-target 报告）的验收测试。

对应宋成义 2026-11-18 09:37 定的口径：六项综合分叫 Design_score、按它排序、
报告能整包下载。
"""

from __future__ import annotations

import json
import zipfile

from zfn_design import gui, pipeline


def test_payload_is_stage1_and_ranked():
    payload = gui.run_request({"seq": gui.demo_payload()["seq"],
                               "n_designs": 2, "max_sites": 5})
    assert payload["ok"] is True

    total = payload["n_sites"]
    assert total > 0

    designs = payload["designs"]
    assert designs, "示例序列应该产出设计"
    totals = [d["total"] for d in designs]
    assert totals == sorted(totals, reverse=True), "设计必须按 Design_score 降序"
    assert [d["rank"] for d in designs] == list(range(1, len(designs) + 1))


def test_bundle_contains_every_artifact(tmp_path, monkeypatch):
    monkeypatch.setattr(gui, "OUT_ROOT", tmp_path / "gui")
    r = gui.run_request({"seq": gui.demo_payload()["seq"], "n_designs": 1, "max_sites": 1})

    files = r["files"]
    assert "bundle" in files and files["bundle"].endswith("zfn_bundle.zip")

    bundle = tmp_path / "gui" / r["out_dir"].split("\\")[-1] / "zfn_bundle.zip"
    assert bundle.is_file()
    with zipfile.ZipFile(bundle) as z:
        names = set(z.namelist())
    assert {"report.html", "report.json", "sites.tsv", "designs.tsv"} <= names
    assert "zfn_bundle.zip" not in names, "包不能把自己装进去"


def test_report_meta_declares_design_score(tmp_path):
    payload = pipeline.run(_records(), n_designs=1, max_sites=5,
                           command="pytest-s1")
    out = pipeline.write_outputs(payload, tmp_path / "out")

    data = json.loads(out["json"].read_text(encoding="utf-8"))
    assert data["meta"]["stage"] == "S1"
    assert data["meta"]["report_kind"] == "all_target"
    assert data["meta"]["design_score_name"] == "Design_score"
    assert data["scoring"]["design_score_keys"] == [
        "evidence", "subsite", "gc_balance", "finger_diversity",
        "complexity", "self_complementarity",
    ]

    html = out["html"].read_text(encoding="utf-8")
    assert "Design_score" in html
    assert "DesignScore" in html and "Stage 1" in html


def _records():
    from zfn_design.models import SequenceRecord
    return [SequenceRecord(name="demo", seq=gui.demo_payload()["seq"],
                           description="", source="pytest")]
