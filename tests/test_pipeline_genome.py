"""全基因组模式的端到端：sites.tsv、重复过滤、multi_copy、truncated。"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from zfn_design import genome, pipeline
from zfn_design.library import load_default_library

from conftest import requires_full_library

CHR1 = ("TACGTAGCTAGCTTAGCTAGGCTAGCGTGGGCGGATATCCGCCCACGCTTACGGATCCTTAG"
        "TAGCTAGCGTGGGCGGATATCCGCCCACGCTTACGGATCCTTAG")
CHR2 = "GGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCGGGCCCAAATTT"


@pytest.fixture()
def mini_genome(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(f">chr1\n{CHR1}\n>chr2\n{CHR2}\n", encoding="utf-8")
    return str(fa)


@pytest.fixture()
def repeats_bed(tmp_path):
    p = tmp_path / "r.bed"
    p.write_text("chr1\t0\t20\n", encoding="utf-8")
    return str(p)


def _rows(path: Path) -> list[dict]:
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_run_genome_writes_sites_tsv(mini_genome, tmp_path):
    out = tmp_path / "o"
    payload = pipeline.run_genome(mini_genome, out, spacers=(6,), limit_sites=5000)

    assert payload["meta"]["genome_mode"] is True
    ts = payload["targets_summary"]
    assert ts["n_sites"] > 0
    # n_sites = 可设计行（full/partial/gap）；multi_copy 行单独计数、单独落盘
    assert ts["n_sites"] == sum(ts["by_status"].values())

    rows = _rows(out / "sites.tsv")
    assert len(rows) == ts["n_sites"] + ts["n_multi_copy"]
    assert list(rows[0].keys())[:4] == ["site_id", "seq_name", "start", "end"]
    assert all(r["designability"] for r in rows)


def test_repeats_are_removed(mini_genome, repeats_bed, tmp_path):
    lib = load_default_library()
    out = tmp_path / "o"
    payload = pipeline.run_genome(
        mini_genome, out, library=lib, spacers=(6,),
        repeats_path=repeats_bed, limit_sites=5000,
    )
    rows = _rows(out / "sites.tsv")
    assert rows
    assert all(r["seq_name"] != "chr1" or int(r["start"]) > 20 for r in rows)
    assert payload["targets_summary"]["scan_stats"]["skipped_repeat"] > 0
    assert any("重复" in w for w in payload["warnings"])


def test_without_repeats_keeps_them(mini_genome, tmp_path):
    out = tmp_path / "o"
    pipeline.run_genome(mini_genome, out, spacers=(6,), limit_sites=5000)
    rows = _rows(out / "sites.tsv")
    assert any(r["seq_name"] == "chr1" and int(r["start"]) <= 20 for r in rows)
    # 没有注释时必须显式提醒
    payload = pipeline.run_genome(mini_genome, tmp_path / "o2", spacers=(6,), limit_sites=1000)
    assert any("未提供重复序列注释" in w for w in payload["warnings"])


def test_multi_copy_is_flagged(mini_genome, tmp_path):
    """两条染色体上出现同一段序列 → 只能留一处，其余标 multi_copy。"""
    out = tmp_path / "o"
    payload = pipeline.run_genome(mini_genome, out, spacers=(6,), limit_sites=5000)
    rows = _rows(out / "sites.tsv")
    dup = [r for r in rows if r["designability"] == "multi_copy"]
    assert dup, "应该能识别出重复出现的同一序列"
    assert payload["targets_summary"]["n_multi_copy"] == len(dup)
    # 重复行只做记录：不进 by_status、不进 n_sites、也不进设计
    assert payload["targets_summary"]["n_sites"] == len(rows) - len(dup)


def test_limit_marks_truncated(mini_genome, tmp_path):
    payload = pipeline.run_genome(mini_genome, tmp_path / "o", spacers=(6,), limit_sites=10)
    assert payload["targets_summary"]["n_sites"] == 10
    assert payload["targets_summary"]["truncated"] is True
    assert any("上限" in w for w in payload["warnings"])


def test_regions_restrict_scan(mini_genome, tmp_path):
    bed = tmp_path / "t.bed"
    bed.write_text("chr2\t0\t40\n", encoding="utf-8")
    payload = pipeline.run_genome(
        mini_genome, tmp_path / "o", spacers=(6,), regions_path=str(bed), limit_sites=5000,
    )
    rows = _rows(tmp_path / "o" / "sites.tsv")
    assert rows and all(r["seq_name"] == "chr2" for r in rows)
    assert payload["targets_summary"]["regions"]["n_intervals_merged"] == 1


@requires_full_library
def test_genome_designs_have_no_offtarget(mini_genome, tmp_path):
    payload = pipeline.run_genome(
        mini_genome, tmp_path / "o", spacers=(6,), limit_sites=300,
        n_designs=1, max_design_sites=1,
    )
    assert payload["designs"], "应该有设计"
    d = payload["designs"][0]
    assert d["score"]["total"] > 0
    # 全基因组脱靶是独立分数：不进 Design_score 的分项
    assert "offtarget" not in {c["key"] for c in d["score"]["components"]}


def test_estimate_matches_full_scan(mini_genome, tmp_path):
    payload = pipeline.run_genome(mini_genome, tmp_path / "o",
                                  spacers=(6,), limit_sites=10_000_000)
    ts = payload["targets_summary"]
    total = sum(ts["by_status"].values()) + ts["n_multi_copy"]
    # estimate_sites 只是数量级预估（双链、不扣重叠），允许偏大，但不能偏小
    assert ts["estimated_total"] >= total
    assert ts["estimated_total"] <= total * 1.5
