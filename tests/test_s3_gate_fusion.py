"""S3 门禁（只允许通过过滤的候选）+ S3 → S4 串联。"""
from __future__ import annotations

import json
import random

import pytest

from zfn_design import fusion_stage, offtarget_stage, pipeline
from zfn_design.models import SequenceRecord

from conftest import requires_full_library

GENOME = 4000


def _genome(tmp_path):
    rnd = random.Random(11)
    genome = "".join(rnd.choice("ACGT") for _ in range(GENOME))
    p = tmp_path / "mini.fna"
    p.write_text(">chrT\n" + genome + "\n", encoding="utf-8")
    return p, genome


def _run_with_filter(tmp_path, seq):
    genome_path, _g = _genome(tmp_path)
    out = tmp_path / "out"
    payload = pipeline.run(
        [SequenceRecord(name="t1", seq=seq)],
        spacers=(6,), fingers=3, n_designs=2, max_sites=8,
        genome_fasta=str(genome_path), command="test",
    )
    pipeline.write_outputs(payload, out)
    return out, payload


def test_s3_requires_target_filter(tmp_path):
    # 模拟一份「带了基因组、但没有走 Target filter」的旧作业：把 genome_filter 摘掉
    _gpath, genome = _genome(tmp_path)
    out, payload = _run_with_filter(tmp_path, genome[1200:1800])
    rj = out / "report.json"
    data = json.loads(rj.read_text(encoding="utf-8"))
    data["targets_summary"].pop("genome_filter", None)
    rj.write_text(json.dumps(data, ensure_ascii=False), encoding="utf-8")
    with pytest.raises(ValueError) as ei:
        offtarget_stage.run_offtarget(out, top=1)
    assert "Target filter" in str(ei.value)


@requires_full_library
def test_s3_then_fusion(tmp_path):
    _gpath, genome = _genome(tmp_path)
    seq = genome[1200:1800]
    out, payload = _run_with_filter(tmp_path, seq)
    assert (payload.get("targets_summary") or {}).get("genome_filter")

    res = offtarget_stage.run_offtarget(out, top=1)
    assert res["ok"] is True
    assert res["n_designs"] == 1
    refreshed = offtarget_stage.load_payload(out)
    assert refreshed["meta"]["stage"] == "S3"
    assert refreshed.get("fusion_ready") is True
    with_ot = [d for d in refreshed["designs"] if d.get("offtarget_score")]
    assert with_ot and 0 <= with_ot[0]["offtarget_score"]["score"] <= 100

    fus = fusion_stage.run_fusion(out, top=1, arch_keys=["single_activator"])
    assert fus["summary"]["n_entries"] >= 1
    entry = fus["entries"][0]
    assert entry.get("offtarget_score") is not None
    assert entry["protein"]
    assert out.joinpath("order_form.csv").exists()
