"""单体型（single）与组合型（paired）必须是两条不同的流水线，不能只换个标签。

这些测试锁住的是「导师说的结果一样」那个 bug：在同一段序列上，
两种构型的位点数、单体数、打分分项、脱靶判据、融合校验都必须不同。
"""

from __future__ import annotations

import pytest

from zfn_design import offtarget, pipeline, targets
from zfn_design.build import build_for_site
from zfn_design.designability import evaluate
from zfn_design.fusion import FusionError, build_fusion
from zfn_design.models import SequenceRecord

SEQ = "TACGTAGCTAGCTTAGCTAGGCTA" + "GCGTGGGCGGATATCCGCCCACGC" + "TTACGGATCCTTAGGCAT"


@pytest.fixture()
def record():
    return SequenceRecord(name="demo", seq=SEQ, description="", source="test")


def _designable_site(record, lib, mode: str, spacers=(6,)):
    """取第一个"真的能设计"的位点。

    不能用 `scan_sequence(...)[0]`：序列开头那些位点通常缺三联体，
    在只带 3 条种子的发布包里会直接变成空列表（早期版本就是这么挂的）。
    """
    for site in targets.scan_sequence(record, spacers=spacers, mode=mode):
        if evaluate(site, lib).status == "full":
            return site
    pytest.skip(f"{mode} 模式下没有完全可设计的位点（当前指头库覆盖不足）")


def test_paired_site_has_two_halves_and_single_has_one(record):
    paired = targets.scan_sequence(record, spacers=(6,), mode="paired")
    single = targets.scan_sequence(record, spacers=(6,), mode="single")

    assert paired, "组合型应当有候选位点"
    assert single, "单体型应当有候选位点"
    p, s = paired[0], single[0]

    assert p.kind == "paired" and p.is_paired
    assert p.spacer and p.right_half and len(p.half_sites) == 2

    assert s.kind == "single" and not s.is_paired
    assert s.spacer == "" and s.right_half == ""
    assert len(s.half_sites) == 1
    # 不变量在单体型下依然成立：site_seq == left + spacer + revcomp(right)
    assert s.site_seq == s.left_half + s.spacer


def test_single_site_count_matches_formula(record):
    hsl = 9
    got = len(targets.scan_sequence(record, spacers=(5, 6, 7), mode="single"))
    per_strand = max(0, len(SEQ) - hsl + 1)
    assert got == per_strand * 2
    # 单体型不依赖 spacer，所以三种 spacer 设置得到同样的位点数
    assert got == len(targets.scan_sequence(record, spacers=(6,), mode="single"))
    # 组合型更多：每个 spacer 一条
    paired = len(targets.scan_sequence(record, spacers=(5, 6, 7), mode="paired"))
    assert paired > got


def test_single_mode_rejects_bad_mode(record):
    with pytest.raises(ValueError):
        targets.scan_sequence(record, mode="dimer")


def test_single_design_has_no_partner(record):
    from zfn_design.library import load_default_library

    lib = load_default_library()
    site = _designable_site(record, lib, "single")
    dz = evaluate(site, lib)
    assert dz.right is None and len(dz.halves) == 1

    designs = build_for_site(site, lib, dz, n_designs=1)
    assert designs
    d = designs[0]
    assert d.is_paired is False
    assert d.partner_fingers == [] and d.partner_protein == ""
    assert len(d.fingers) == 3


def test_single_scoring_drops_spacer_and_renormalizes(record):
    from zfn_design.library import load_default_library
    from zfn_design.score import score_design

    lib = load_default_library()

    def scored(mode):
        site = _designable_site(record, lib, mode)
        dz = evaluate(site, lib)
        design = build_for_site(site, lib, dz, n_designs=1)[0]
        return score_design(design)

    paired = scored("paired")
    single = scored("single")

    keys_p = [c.key for c in paired.components]
    keys_s = [c.key for c in single.components]
    assert "spacer" in keys_p
    assert "spacer" not in keys_s
    assert len(keys_s) == len(keys_p) - 1

    # 权重重新归一化后总分仍是 0–100，且权重和为 1
    assert 0.0 <= single.total <= 100.0
    assert abs(sum(c.weight for c in single.components) - 1.0) < 1e-9
    assert any("重归一化" in n for n in single.notes)


def test_single_offtarget_counts_unpaired_single_sites(record):
    """同样一段序列、同样一个 9 bp 半位点：单侧命中的一定不少于成对。"""
    half = "GCGTGGGCG"
    single = offtarget.search_single([record], half_site=half, max_mm=2)
    paired = offtarget.search(
        [record], left_half=half, right_half=half, spacers=(6,), max_mm_left=2, max_mm_right=2
    )
    assert single.hits, "单侧搜索应当有命中"
    assert len(single.hits) >= len(paired.hits)
    assert all(h.spacer == "" and h.right_half == "" for h in single.hits)


def test_fusion_rejects_dimer_architecture_on_single_design(record):
    from zfn_design.library import load_default_library

    lib = load_default_library()
    site = _designable_site(record, lib, "single")
    design = build_for_site(site, lib, evaluate(site, lib), n_designs=1)[0]

    with pytest.raises(FusionError):
        build_fusion(design, arch_key="paired_nuclease")
    ok = build_fusion(design, arch_key="single_activator")
    assert ok["is_complete"] is True


def test_pipeline_modes_differ_end_to_end(record):
    paired = pipeline.run([record], spacers=(6,), architecture="paired_nuclease", max_sites=3, n_designs=1)
    single = pipeline.run([record], spacers=(6,), architecture="single_activator", max_sites=3, n_designs=1)

    assert paired["meta"]["params"]["mode"] == "paired"
    assert single["meta"]["params"]["mode"] == "single"

    assert paired["targets_summary"]["n_sites"] != single["targets_summary"]["n_sites"]

    pd = paired["designs"][0]
    sd = single["designs"][0]
    assert pd["monomer_right"] is not None
    assert sd["monomer_right"] is None
    assert "spacer" not in [c["key"] for c in sd["score"]["components"]]
    # 同样的处理位点数，分数不应恰好相同——这是那个 bug 的直接判据
    assert pd["score"]["total"] != sd["score"]["total"]
