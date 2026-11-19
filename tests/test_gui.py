"""图形界面的后端测试。

界面本身是浏览器里的 HTML，没法在这里点按钮；但按钮点下去调用的那个函数
(:func:`zfn_design.gui.run_request`) 是纯 Python，可以像普通函数一样测。
这样"界面里跑出来的数 == 命令行跑出来的数"就不只是口头承诺。
"""

from __future__ import annotations

import json

import pytest

from zfn_design import gui, pipeline
from zfn_design.ingest import IngestError
from zfn_design.sequtil import SequenceError


def test_demo_payload_is_plain_sequence():
    d = gui.demo_payload()
    assert set(d) == {"name", "seq"}
    assert d["seq"] and ">" not in d["seq"]
    assert set(d["seq"]) <= set("ACGTN")


def test_run_request_matches_pipeline():
    """界面后端的结果必须和直接调 pipeline.run 完全一致。"""
    seq = gui.demo_payload()["seq"]
    via_gui = gui.run_request({"seq": seq, "spacers": [5, 6, 7], "n_designs": 3, "max_sites": 20})

    from zfn_design.ingest import parse_fasta
    records = parse_fasta(">input\n" + seq)
    direct = pipeline.run(records, spacers=(5, 6, 7), n_designs=3, max_sites=20)

    assert via_gui["ok"] is True
    assert via_gui["n_sites"] == len(direct["sites"])
    assert via_gui["by_status"] == direct["targets_summary"]["by_status"]
    assert [d["id"] for d in via_gui["designs"]] == [d["design_id"] for d in direct["designs"]]
    assert [d["total"] for d in via_gui["designs"]] == [d["score"]["total"] for d in direct["designs"]]


def test_run_request_shape():
    r = gui.run_request({"seq": gui.demo_payload()["seq"], "n_designs": 2, "max_sites": 5})
    for key in ("elapsed_ms", "run_id", "inputs", "library", "designs", "files",
                "disclaimers", "warnings", "out_dir"):
        assert key in r
    # n_designs 是"每个位点"的上限，总量还要乘位点数（max_sites=5）
    assert len(r["designs"]) <= 2 * 5
    assert set(r["files"]) >= {"json", "html", "sites_tsv", "designs_tsv"}
    # 每个设计带 7 个打分分项（六项设计分 + spacer 诊断项；脱靶是独立分数，不进这里）
    d = r["designs"][0]
    assert len(d["components"]) == 7
    keys = {c["key"] for c in d["components"]}
    assert keys == {"evidence", "subsite", "spacer", "gc_balance",
                    "finger_diversity", "complexity", "self_complementarity"}
    assert next(c for c in d["components"] if c["key"] == "spacer")["weight"] == 0.0
    assert abs(sum(c["weighted"] for c in d["components"]) - d["total"]) < 0.02
    assert len(r["sites"]) <= 60  # 只回传前 60 条给前端


def test_files_are_written():
    r = gui.run_request({"seq": gui.demo_payload()["seq"], "n_designs": 1, "max_sites": 1})
    from pathlib import Path
    out = Path(r["out_dir"])
    assert (out / "report.html").exists()
    assert (out / "sites.tsv").exists()
    payload = json.loads((out / "report.json").read_text(encoding="utf-8"))
    assert payload["meta"]["version"]


def test_run_dir_is_findable(tmp_path, monkeypatch):
    """输出目录必须落在固定根下、名字能看出时间与序列长度，而不是随机临时目录。"""
    from pathlib import Path

    monkeypatch.setattr(gui, "OUT_ROOT", tmp_path / "gui")
    r = gui.run_request({"seq": gui.demo_payload()["seq"], "n_designs": 1, "max_sites": 1})
    d = Path(r["out_dir"])
    assert d.is_dir()
    assert d.parent == tmp_path / "gui"
    assert "bp" in d.name and d.name[:8].isdigit()  # 20261114-151710_300bp_xxxxxx
    assert len(list(d.iterdir())) >= 5


def test_reveal_reports_failure_instead_of_raising(tmp_path):
    """不存在的目录：必须返回 (False, 原因)，不能把异常抛给请求线程。"""
    ok, msg = gui._reveal(tmp_path / "definitely-not-here-9f3a")
    assert ok is False
    assert msg


def test_empty_sequence_is_rejected():
    with pytest.raises(IngestError):
        gui.run_request({"seq": "   "})


def test_bad_characters_are_rejected():
    with pytest.raises(SequenceError):
        gui.run_request({"seq": "ACGTZZZACGT"})


def test_fasta_header_gives_name():
    r = gui.run_request({"seq": ">my_gene my description\nGCGTGGGCGGATATCCGCCCACGC",
                         "n_designs": 1, "max_sites": 3})
    assert r["inputs"][0]["name"] == "my_gene"


def test_params_are_clamped():
    """前端传什么数字都不能把服务打挂。"""
    r = gui.run_request({"seq": gui.demo_payload()["seq"], "n_designs": 999,
                         "max_sites": 99999, "max_mm_left": 99, "max_mm_right": -5})
    assert r["params"]["n_designs"] == 10
    assert r["params"]["max_sites"] == 500
    assert r["params"]["max_mm_left"] == 3
    assert r["params"]["max_mm_right"] == 0


def test_info_payload():
    info = gui.info_payload()
    assert info["version"]
    assert info["library"]["triplets_total"] == 64
    assert any(a["key"] == "paired_nuclease" for a in info["architectures"])


def test_html_page_is_self_contained():
    html = gui.INDEX_HTML
    assert "/api/run" in html and "/api/info" in html
    assert "http://" not in html.replace("http://127.0.0.1", "")  # 没有外链资源
    assert "<script src=" not in html and "<link rel=\"stylesheet\"" not in html
