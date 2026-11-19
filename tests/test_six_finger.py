"""六指阵列（18 bp 半位点）与 spacer 自定义解析。

六指阵列是本工具里最容易"看起来对、其实是编的"一块，所以这里把
每一项都被迫做出的假设都钉成断言：

* 半位点长度 = 3 × 指头数，18 bp 半位点在扫描时确实按 18 bp 切；
* 六指蛋白 = 两个三指模块 + 一个 TGEKP 连接肽，长度 184 aa；
* 阵列 N→C 顺序里 F1–F3 读 3' 那 9 bp，F4–F6 读 5' 那 9 bp；
* F1/F4 在 B1H 库里没有本位上下文，必须走 fallback 且被记进报告；
* spacer 文本既能写 ``6``，也能写 ``5,6`` 或 ``5-7``。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from zfn_design import framework, pipeline, targets
from zfn_design.cli import _parse_spacers
from zfn_design.designability import context_for_slot, plan_halfsite
from zfn_design.models import SequenceRecord

# 发布包默认不带 B1H 库（申请制数据），这几条测的是"64/64 全覆盖时仍然要走 fallback"，
# 种库下没有意义，直接 skip 而不是让发布包跑出一片红。
from zfn_design.library import has_full_library
needs_b1h = pytest.mark.skipif(not has_full_library(), reason="需要 B1H 指头库（发布包默认不含，见 data/finger_library/README.md）")


def test_half_site_length_follows_finger_count():
    rec = SequenceRecord(name="t", seq="ACGT" * 40, description="", source="test")
    three = targets.scan_many([rec], spacers=(6,), half_site_len=9)
    six = targets.scan_many([rec], spacers=(6,), half_site_len=18)
    assert all(len(s.left_half) == 9 and len(s.right_half) == 9 for s in three)
    assert all(len(s.left_half) == 18 and len(s.right_half) == 18 for s in six)
    assert len(three) == targets.expected_site_count(160, half_site_len=9, spacers=(6,))
    assert len(six) == targets.expected_site_count(160, half_site_len=18, spacers=(6,))
    assert len(six) < len(three)


def test_six_finger_protein_shape():
    h = ["RSDELTR", "RSDHLTT", "RSDERKR"] * 2
    prot = framework.build_array_protein(h)
    assert len(prot) == 184
    # 两个模块之间必须是 TGEKP（模块 A 的尾巴 + 连接肽 + 模块 B 去掉起始 Met）
    assert "HTKIHLRQKDTGEKPERPYACPVESCDRRFS" in prot
    assert prot.startswith("MERPYACPVESCDRRFSRSDELTR")
    assert prot.endswith("HTKIHLRQKD")
    # 六个螺旋按 N→C 顺序出现，且各自位置正确
    assert prot.count("RSDELTR") == 2
    assert prot.count("RSDHLTT") == 2
    assert prot.count("RSDERKR") == 2


def test_template_selection_by_finger_count():
    assert framework.template_for_fingers(3) == "zif268_1aay"
    assert framework.template_for_fingers(6) == "zif268_x2_junction"
    # 其它指数（2–12）走通用拼接骨架
    assert framework.template_for_fingers(4) == framework.GENERIC_KEY
    assert framework.template_for_fingers(12) == framework.GENERIC_KEY
    with pytest.raises(framework.FrameworkError):
        framework.template_for_fingers(1)
    with pytest.raises(framework.FrameworkError):
        framework.template_for_fingers(13)


def test_context_mapping_for_six_fingers():
    """六指阵列里 2/5 指用 F2 上下文、3/6 指用 F3，1/4 指没有本位数据。"""
    assert context_for_slot("F1", 3) == "F1"
    assert context_for_slot("F2", 6) == "F2"
    assert context_for_slot("F3", 6) == "F3"
    assert context_for_slot("F5", 6) == "F2"
    assert context_for_slot("F6", 6) == "F3"
    assert context_for_slot("F1", 6) == "F1"
    assert context_for_slot("F4", 6) == "F4"


@needs_b1h
def test_six_finger_plan_reports_context_fallback(default_library):
    """F1/F4 没有本位数据 → 必须走 fallback 并记录下来。"""
    plan = plan_halfsite("GCGTGGGCGAAATTTCCC", default_library)
    assert plan.n_fingers == 6
    assert set(plan.slots) == {f"F{i}" for i in range(1, 7)}
    assert "F1" in plan.context_fallback
    assert "F4" in plan.context_fallback
    assert not plan.missing            # 64 个三联体全覆盖


@needs_b1h
def test_pipeline_six_fingers_end_to_end(demo_record, default_library):
    payload = pipeline.run(
        [demo_record],
        library=default_library,
        spacers=(6,),
        fingers=6,
        n_designs=1,
        max_sites=2,
    )
    assert payload["meta"]["params"]["fingers"] == 6
    assert payload["meta"]["params"]["half_site_bp"] == 18
    assert payload["designs"], "六指阵列应当能产出设计"
    d = payload["designs"][0]
    assert len(d["monomer_left"]["protein"]) == 184
    assert len(d["monomer_left"]["triplets_N_to_C"]) == 6
    assert len(d["target"]["left_half"]) == 18
    # 「三指只读 9 bp」那条限制不该出现在六指运行里，取而代之的是六指的连接警告
    warnings = " ".join(payload["warnings"])
    assert "9 bp" not in warnings
    assert "两个三指模块" in warnings


def test_pipeline_three_fingers_carries_uniqueness_warning(demo_record, default_library):
    payload = pipeline.run([demo_record], library=default_library, spacers=(6,),
                           fingers=3, n_designs=1, max_sites=1)
    warnings = " ".join(payload["warnings"])
    assert "只读 9 bp" in warnings
    assert "并不唯一" in warnings


@pytest.mark.parametrize("text,expected", [
    ("6", (6,)),
    ("5,6", (5, 6)),
    ("5-7", (5, 6, 7)),
    ("5,7-8", (5, 7, 8)),
    (" 6 , 5 ", (5, 6)),
])
def test_parse_spacers(text, expected):
    assert _parse_spacers(text) == expected


@pytest.mark.parametrize("text", ["", "  ", "abc"])
def test_parse_spacers_rejects_junk(text):
    with pytest.raises(SystemExit):
        _parse_spacers(text)
