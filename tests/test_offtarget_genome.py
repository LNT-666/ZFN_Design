"""全基因组脱靶扫描的测试。

原则：**不比较"分数好不好"**，只锁"找得到 / 数得对 / 两种引擎一致"。
"""

from __future__ import annotations

import random

import pytest

from zfn_design import offtarget_genome, pipeline
from zfn_design.sequtil import revcomp

LEFT = "GCGTGGGCG"
RIGHT = "GCGTGGGCG"
SPACER = "GATATC"


def _rnd(n: int, seed: int = 11) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


@pytest.fixture(scope="module")
def planted(tmp_path_factory) -> dict:
    """在一个随机基因组里植入：0 错配位点、左侧 1 错配位点。"""
    site = LEFT + SPACER + revcomp(RIGHT)
    mm1 = LEFT[:8] + ("A" if LEFT[8] != "A" else "C") + SPACER + revcomp(RIGHT)
    # 每段用不同 seed：同一 seed 会让三段互为前缀，整条序列自带大量重复，
    # 全基因组去重（site_seq 计数）会把候选清空——那是 fixture 的坑，不是引擎的。
    seq = _rnd(3000, 1) + site + _rnd(2000, 2) + mm1 + _rnd(2500, 3)
    p = tmp_path_factory.mktemp("ot") / "plant.fa"
    p.write_text(">chr1\n" + seq + "\n", encoding="utf-8")
    return {"path": str(p), "seq": seq,
            "exact_at": seq.index(site), "mm1_at": seq.index(mm1)}


def test_finds_planted_site(planted):
    r = offtarget_genome.genome_search(
        planted["path"], LEFT, RIGHT, spacers=(6,), max_mm_left=2, max_mm_right=2)
    assert r["n_total_hits"] >= 2
    exact = [h for h in r["hits"] if h["mm_total"] == 0]
    assert len(exact) == 1
    assert exact[0]["start0"] == planted["exact_at"]
    assert exact[0]["weight"] == 1.0


def test_one_mismatch_is_counted_on_the_right_side(planted):
    r = offtarget_genome.genome_search(
        planted["path"], LEFT, RIGHT, spacers=(6,), max_mm_left=2, max_mm_right=2)
    one = [h for h in r["hits"] if h["mm_total"] == 1]
    assert len(one) == 1
    assert one[0]["start0"] == planted["mm1_at"]
    assert one[0]["mm_left"] == 1 and one[0]["mm_right"] == 0
    assert one[0]["weight"] == 0.25


def test_identical_halves_are_not_counted_twice(planted):
    """左右半位点相同时，两种方向模式完全一样，不能把同一个位点数两遍。"""
    r = offtarget_genome.genome_search(
        planted["path"], LEFT, LEFT, spacers=(6,), max_mm_left=2, max_mm_right=2)
    starts = [(h["chrom"], h["start0"]) for h in r["hits"]]
    assert len(starts) == len(set(starts))


def test_burden_matches_weights(planted):
    r = offtarget_genome.genome_search(
        planted["path"], LEFT, RIGHT, spacers=(6,), max_mm_left=2, max_mm_right=2)
    assert r["burden"] == pytest.approx(sum(h["weight"] for h in r["hits"]), abs=1e-6)


def test_python_and_numpy_agree(planted):
    a = offtarget_genome.genome_search(
        planted["path"], LEFT, RIGHT, spacers=(6,), allow_numpy=True)
    b = offtarget_genome.genome_search(
        planted["path"], LEFT, RIGHT, spacers=(6,), allow_numpy=False)
    assert a["by_mm"] == b["by_mm"]
    assert a["burden"] == pytest.approx(b["burden"], abs=1e-9)
    assert sorted(h["start0"] for h in a["hits"]) == sorted(h["start0"] for h in b["hits"])


def test_max_hits_truncates_the_list_but_not_the_count(planted):
    r = offtarget_genome.genome_search(
        planted["path"], LEFT, RIGHT, spacers=(6,), max_mm_left=2,
        max_mm_right=2, max_hits=1)
    assert len(r["hits"]) == 1
    assert r["truncated"] is True
    assert r["n_total_hits"] >= 2  # 统计仍是全量


def test_spacer_长度不同_命中不同(planted):
    """spacer 差一个碱基就不是同一种位点，不该命中。"""
    r = offtarget_genome.genome_search(
        planted["path"], LEFT, RIGHT, spacers=(7,), max_mm_left=0, max_mm_right=0)
    assert r["by_mm"] == {}


def test_pipeline_reports_offtarget_separately(planted, tmp_path):
    """全基因组脱靶是独立一步：进 offtarget_score 与 offtargets.tsv，不进 Design_score。"""
    out = tmp_path / "g"
    # limit 要够大，才能扫到植入在 ~3000 bp 处的那个位点（每次窗口 ≈ 2 条链）
    p = pipeline.run_genome(
        planted["path"], out, spacers=(6,), fingers=3, limit_sites=9000,
        n_designs=1, ot_designs=1, ot_max_mm=2)
    assert p["genome_offtarget"], "应该至少有一个设计做了全基因组脱靶"
    did = next(iter(p["genome_offtarget"]))
    assert (out / "offtargets.tsv").exists()
    txt = (out / "offtargets.tsv").read_text(encoding="utf-8-sig")
    assert did in txt
    d = next(x for x in p["designs"] if x["design_id"] == did)
    ot = d["offtarget_score"]
    assert ot["source"] == "genome_wide"
    assert 0.0 <= ot["score"] <= 100.0
    assert ot["scanned_bp"] > 0
    # 两套分并列，不相加
    assert "offtarget" not in {c["key"] for c in d["score"]["components"]}
    assert any("全基因组" in w for w in p["warnings"])
    # 做了脱靶就不该再警告"未做脱靶"
    assert not any("未做脱靶" in w for w in p["warnings"])
