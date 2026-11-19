"""S4：Fusion design + 人源化 CDS 的测试。"""
from __future__ import annotations

import json
from pathlib import Path

from zfn_design import codonopt, fusion_stage, pipeline
from zfn_design.models import SequenceRecord

from conftest import requires_full_library

SEQ = "TACGTAGCTAGCTTAGCTAGGCTA" "GCGTGGGCGGATATCCGCCCACGC" "TTACGGATCCTTAGGCAT"

# 实验室自备片段的替身：只要长度与组成合理，拼装逻辑与它们无关。
KRAB_STUB = "RSVTFKDVAVDFTQEEWQQLDPSQKRLYRDVMQENYKNMFQEVAGQIKQLTESLRQNLDQLKS"
FOKI_STUB = ("QLVKSELEEKKSELRHKLKYVPHEYIELIEIARNSTQDRILEMKVMEFFMKVYGYRGKHLGGSRK"
             "PDGAIYTVGSPIDYGVIVDTKAYSGGYNLPIGQADEMQRYVEENQTRNKHINPNEWWKVYPSSVT")


def _payload(tmp_path: Path) -> dict:
    recs = [SequenceRecord(name="s4", seq=SEQ, source="test")]
    payload = pipeline.run(recs, spacers=(6,), fingers=3, n_designs=2, max_sites=12,
                           fusion_keys=[], command="test")
    out = tmp_path / "run"
    pipeline.write_outputs(payload, out)
    return {"out": out, "payload": payload}


def test_codonopt_is_deterministic_and_constrained():
    protein = "MRPHHHHHHGSGGGGSLEPGEKPYKCPECGKSFSQSSNLVRHQRTHTGA*"
    a = codonopt.optimize(protein, seed=11)
    b = codonopt.optimize(protein, seed=11)
    assert a["cds"] == b["cds"]
    assert a["ok"] is True
    assert 0.40 <= a["gc"] <= 0.60
    assert a["forbidden_sites"] == []
    assert a["homopolymer_max"] <= codonopt.MAX_HOMOPOLYMER
    assert codonopt.translate(a["coding"]) == protein
    assert len(a["coding"]) == 3 * len(protein)
    assert len(a["cds"]) == len(a["coding"]) + len(a["kozak"])


def test_codonopt_refuses_placeholders():
    r = codonopt.optimize("MRP[FokI_CD:user_supplied]")
    assert r["ok"] is False
    assert r["cds"] == ""
    assert any("占位符" in n for n in r["notes"])


def test_codonopt_requires_start_and_stop():
    ok = codonopt.optimize("MPKKKRKVMERPYAC")
    assert ok["ok"] is True
    assert ok["cds"].startswith(ok["kozak"] + "ATG")
    assert ok["coding"].startswith("ATG")
    assert ok["coding"][-3:] in codonopt.STOP_CODONS
    assert ok["cds"][-3:] in codonopt.STOP_CODONS
    assert ok["start_codon"] == "ATG"
    assert ok["cds_coding_start"] == len(ok["kozak"])
    assert codonopt.translate(ok["coding"]) == "MPKKKRKVMERPYAC*"
    bad = codonopt.optimize("PKKKRKVMERPYAC")  # 没有起始 M
    assert bad["ok"] is False
    assert bad["cds"] == ""


@requires_full_library
def test_fusion_stage_defaults_all_ok(tmp_path):
    p = _payload(tmp_path)
    res = fusion_stage.run_fusion(p["out"], top=2)
    s = res["summary"]
    assert s["n_entries"] == 6  # 2 个设计 × 3 种架构
    # 成对核酸酶每条设计出左右两条链：2×2 条单链 + 2×2 条双链
    assert s["n_chains"] == 8
    assert s["n_ok"] == 8, [e.get("cds_notes") for e in res["entries"]]
    assert s["n_blocked"] == 0
    for e in res["entries"]:
        assert e["cds_status"] == "ok"
        assert e["protein"].startswith("M"), e["protein"][:12]
        assert "PKKKRKV" in e["protein"]
        assert e["has_nls"] is True
        assert e["cds"].startswith(e["cds_kozak"] + "ATG")
        assert e["cds"][-3:] in codonopt.STOP_CODONS
        assert e["cds_start_codon"] == "ATG"
        for ch in e["chains"]:
            assert ch["cds_status"] == "ok", ch.get("cds_notes")
            assert codonopt.translate(ch["cds_coding"]).rstrip("*") == ch["protein"].rstrip("*")
        if e["architecture"] == "paired_nuclease":
            assert e["n_chains"] == 2
            sides = {ch["side"] for ch in e["chains"]}
            assert sides == {"left", "right"}
            parts = {ch["side"]: ch["fok_part"] for ch in e["chains"]}
            # published obligate heterodimer pair is ELD:KKK (Doyon 2011)
            assert parts == {"left": "FokI_CD_ELD", "right": "FokI_CD_KKK"}
    saved = json.loads((p["out"] / "fusions.json").read_text(encoding="utf-8"))
    assert saved["meta"]["stage"] == "S4"


def test_fusion_stage_flags_placeholder_parts(tmp_path):
    p = _payload(tmp_path)
    res = fusion_stage.run_fusion(p["out"], top=1, arch_keys=["paired_nuclease"],
                                  parts={"FokI_CD": "[FokI_CD:user_supplied]"})
    s = res["summary"]
    assert s["n_entries"] == 1
    assert s["n_ok"] == 0
    assert s["n_blocked"] + s["n_warn"] == 2  # 左右两条链都缺 FokI
    entry = res["entries"][0]
    assert entry["cds_status"] != "ok"
    assert not entry["cds"]
    assert any("占位符" in str(e) for e in (entry.get("cds_notes") or []))


def test_fusion_stage_full_with_supplied_parts(tmp_path):
    p = _payload(tmp_path)
    res = fusion_stage.run_fusion(p["out"], top=1, arch_keys=list(fusion_stage.DEFAULT_ARCHS),
                                  parts={"KRAB": KRAB_STUB, "FokI_CD": FOKI_STUB})
    s = res["summary"]
    assert s["n_ok"] == 4, [e.get("cds_notes") for e in res["entries"]]  # 2 单链 + 1 双链两条
    for e in res["entries"]:
        assert e["cds_status"] == "ok"
        for ch in e["chains"]:
            assert 0.40 <= ch["cds_gc"] <= 0.60
            assert ch["cds_forbidden_sites"] == []
    assert (p["out"] / "cds.fasta").exists()
    assert (p["out"] / "order_form.csv").exists()
    import zipfile
    with zipfile.ZipFile(p["out"] / "zfn_bundle.zip") as z:
        names = set(z.namelist())
    assert {"fusions.json", "proteins.tsv", "cds.tsv", "cds.fasta", "order_form.csv"} <= names


def test_fusion_stage_rejects_bad_design_id(tmp_path):
    p = _payload(tmp_path)
    try:
        fusion_stage.run_fusion(p["out"], design_ids=["nope"])
    except ValueError as exc:
        assert "不在这次作业里" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("应当报错")
