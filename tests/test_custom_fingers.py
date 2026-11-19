"""任意指数（2–12）阵列：骨架拼接、上下文映射、端到端。

这三件事以前是"三指或六指二选一"，现在开放成 Custom。
开放的前提是能说清"多出来的部分从哪来、哪里没有证据"：

* N 端、两种原生连接区、C 端都逐字符取自 1AAY；
* 三指一模块，模块之间用 ``TGEKP`` + 去掉 Met 的 N 端接起来
  —— 和 ``zif268_x2_junction`` 用的是同一个接法；
* 因此 **n 指蛋白长度恒为 28n + 6**（n=3 → 90 aa 正好还原 1AAY，
  n=6 → 184 aa 与专门写的六指模板逐字符相同）。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from zfn_design import designability, framework, pipeline, sequtil

H3 = ["RSDELTR", "RSDHLTT", "RSDERKR"]
H12 = (H3 * 4)[:12]

# 发布包默认不带 B1H 库（申请制数据）。四指阵列要 8 个三联体，种子库只覆盖 2 个，
# 端到端那条必须有 B1H 才有设计可出——在发布包里 skip，而不是留一片红。
from zfn_design.library import has_full_library
needs_b1h = pytest.mark.skipif(not has_full_library(), reason="需要 B1H 指头库（发布包默认不含，见 data/finger_library/README.md）")


def _fill(template_key: str, helices) -> str:
    s = framework.get_template(template_key)["sequence_template"]
    for i, h in enumerate(helices, start=1):
        s = s.replace(f"{{h{i}}}", h)
    return s


def test_generic_builder_reproduces_the_two_known_templates():
    """通用拼接在 n=3 / n=6 时必须与手写模板逐字符一致——这是它可信的唯一理由。"""
    assert framework.build_array_protein(H3) == _fill("zif268_1aay", H3)
    assert framework.build_array_protein(H12[:6]) == _fill("zif268_x2_junction", H12[:6])


@pytest.mark.parametrize("n", [2, 4, 5, 7, 8, 9, 10, 11, 12])
def test_generic_length_is_28n_plus_6(n):
    prot = framework.build_array_protein(H12[:n])
    assert len(prot) == 28 * n + 6
    assert prot.startswith("MERPYACPVESCDRRFS")
    assert prot.endswith("HTKIHLRQKD")
    # 每个螺旋按 N→C 顺序只出现一次
    for i, h in enumerate(H12[:n]):
        assert prot.count(h) >= 1


def test_generic_module_junction_is_tgekp():
    """第 3→4 指的接法必须是 TGEKP + 模块 N 端（去掉起始 Met）。"""
    prot = framework.build_array_protein(H12[:4])
    assert framework._G_MODULE_START in prot
    assert "HTKIHLRQKD" not in prot[:-10]  # C 端只出现一次，在末尾


def test_finger_count_bounds():
    assert framework.MIN_FINGERS == 2 and framework.MAX_FINGERS == 12
    with pytest.raises(framework.FrameworkError):
        framework.build_array_protein(H12[:1])
    with pytest.raises(framework.FrameworkError):
        framework.build_array_protein(H12 + H12[:1])


def test_context_mapping_generalises_by_module_position():
    """任意指数：模块内第 2 指→F2、第 3 指→F3、第 1 指→本位（库里无数据）。"""
    for n in (4, 5, 7, 9):
        slots = designability._slot_keys(n)
        got = [designability.context_for_slot(s, n) for s in slots]
        expect = []
        for pos in range(1, n + 1):
            within = ((pos - 1) % 3) + 1
            expect.append("F2" if within == 2 else "F3" if within == 3 else f"F{pos}")
        assert got == expect, f"n={n}"


def test_half_site_length_follows_finger_count():
    assert len(sequtil.triplets_by_finger("GCGTGGGCG")) == 3
    assert len(sequtil.triplets_by_finger("GCGTGGGCG" * 2)) == 6
    four = sequtil.triplets_by_finger("GCGTGGGCGTGG")
    assert list(four) == ["F1", "F2", "F3", "F4"]
    assert four["F1"] == "TGG"          # F1 读 3' 端那 3 bp


@needs_b1h
def test_pipeline_runs_with_four_fingers(demo_record, default_library):
    payload = pipeline.run(
        [demo_record], library=default_library,
        spacers=(6,), fingers=4, max_sites=2, n_designs=1,
    )
    assert payload["meta"]["params"]["fingers"] == 4
    assert payload["meta"]["params"]["half_site_bp"] == 12
    for s in payload["sites"]:
        assert len(s["left_half"]) == 12
        assert len(s["right_half"]) == 12
    designs = payload["designs"]
    assert designs, "4 指阵列在 demo 上应当有设计"
    d = designs[0]
    assert len(d["monomer_left"]["fingers"]) == 4
    assert len(d["monomer_left"]["protein"]) == 28 * 4 + 6
