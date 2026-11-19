"""S3：全基因组 off-target score 与报告接线。"""

from __future__ import annotations

import json
import random

from zfn_design import genome_offtarget, pipeline
from zfn_design.models import SequenceRecord

from conftest import requires_full_library


def _demo_seq(n: int = 180, seed: int = 7) -> str:
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def _write_fasta(path, chrom: str, seq: str) -> str:
    path.write_text(f">{chrom}\n{seq}\n", encoding="utf-8")
    return str(path)


def test_select_top_orders_by_design_score():
    designs = [
        type("D", (), {"design_id": "a", "score": type("S", (), {"total": 41.0})()})(),
        type("D", (), {"design_id": "b", "score": type("S", (), {"total": 88.5})()})(),
    ]
    picked = genome_offtarget.select_top(designs, 1)
    assert [d.design_id for d in picked] == ["b"]


@requires_full_library
def test_s3_hook_sets_stage_and_offtarget_score(tmp_path):
    seq = _demo_seq()
    # 基因组里放一份完全相同的序列：候选会命中自己，burden 很小但不为 0
    genome = _write_fasta(tmp_path / "mini.fa", "chrT", _demo_seq(2000, seed=3) + seq)
    payload = pipeline.run(
        [SequenceRecord(name="input", seq=seq)],
        spacers=(6,),
        fingers=3,
        n_designs=2,
        max_sites=8,
        genome_fasta=genome,
        offtarget_top=2,
        offtarget_max_mm=2,
        fusion_keys=[],
        command="test-s3",
    )
    assert payload["meta"]["stage"] == "S3"
    assert payload["meta"]["report_kind"] == "offtarget_score"

    ot = payload["offtarget_genome"]
    assert ot and ot["designs"], "S3 必须给出逐设计脱靶记录"
    for rec in ot["designs"]:
        assert 0.0 <= rec["offtarget_score"] <= 100.0
        assert rec["burden"] >= 0.0
        assert rec["scanned_bp"] > 0
    ranks = [r["burden_rank"] for r in ot["designs"]]
    assert ranks == sorted(ranks)
    assert ot["summary"]["n_designs"] == len(ot["designs"])


@requires_full_library
def test_s3_outputs_include_offtarget_tsv_and_report_section(tmp_path):
    seq = _demo_seq(seed=11)
    genome = _write_fasta(tmp_path / "mini2.fa", "chrT", _demo_seq(1500, seed=5) + seq)
    payload = pipeline.run(
        [SequenceRecord(name="input", seq=seq)],
        spacers=(6,),
        fingers=3,
        n_designs=1,
        max_sites=5,
        genome_fasta=genome,
        offtarget_top=1,
        fusion_keys=[],
        command="test-s3-files",
    )
    out = pipeline.write_outputs(payload, tmp_path / "run")
    assert (out["offtarget_tsv"]).is_file()
    html = (tmp_path / "run" / "report.html").read_text(encoding="utf-8")
    assert "Off-target score" in html
    assert "Stage 3" in html
    tsv = (tmp_path / "run" / "offtarget_genome.tsv").read_text(encoding="utf-8-sig")
    first = tsv.splitlines()[0].split("\t")
    assert "offtarget_score" in first and "burden" in first
    # zip 里也要带上
    import zipfile

    with zipfile.ZipFile(tmp_path / "run" / pipeline.BUNDLE_NAME) as z:
        assert "offtarget_genome.tsv" in z.namelist()
    meta = json.loads((tmp_path / "run" / "report.json").read_text(encoding="utf-8"))
    assert meta["meta"]["stage"] == "S3"
