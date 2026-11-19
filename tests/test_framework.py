"""骨架装配测试。

最重要的一条：把 1AAY 的三个螺旋填回模板，必须**一字不差**地还原
PDB 里那条链的序列。这是整个"输出蛋白序列"功能的真值校验。
"""

from __future__ import annotations

import pytest

from conftest import ZIF268_PROTEIN
from zfn_design import framework


def test_roundtrip_reproduces_pdb_1aay():
    seq = framework.build_array_protein(["RSDELTR", "RSDHLTT", "RSDERKR"])
    assert seq == ZIF268_PROTEIN
    assert len(seq) == 90


def test_helix_substitution_changes_only_helix():
    base = framework.build_array_protein(["RSDELTR", "RSDHLTT", "RSDERKR"])
    mutant = framework.build_array_protein(["RSDELTR", "RSDHATT", "RSDERKR"])
    assert len(base) == len(mutant)
    diffs = [i for i, (a, b) in enumerate(zip(base, mutant)) if a != b]
    # RSDHLTT -> RSDHATT 只改 1 个残基
    assert len(diffs) == 1


def test_template_metadata():
    d = framework.describe("zif268_1aay")
    assert d["n_fingers"] == 3
    assert d["source_id"] == "PDB:1AAY"
    assert d["template_length"] == 90


def test_wrong_helix_count_rejected():
    # 1 指低于下限（2 指起才拼得出像样的阵列）
    with pytest.raises(framework.FrameworkError):
        framework.build_array_protein(["RSDELTR"])
    # 指头数与指定模板不匹配
    with pytest.raises(framework.FrameworkError):
        framework.build_array_protein(["RSDELTR"] * 4, template_key="zif268_1aay")


def test_non_amino_acid_rejected():
    with pytest.raises(framework.FrameworkError):
        framework.build_array_protein(["RSDELTR", "RSDHLT1", "RSDERKR"])


def test_placeholder_x_allowed_but_warned():
    # 位置 4（索引 4）不接触 DNA，允许 X 占位
    seq = framework.build_array_protein(["RSDELTR", "RSDLXTT", "RSDERKR"])
    assert "X" in seq
    warns = framework.warnings_for(["RSDELTR", "RSDLXTT", "RSDERKR"])
    assert warns and "占位符" in warns[0]


def test_lowercase_input_normalised():
    a = framework.build_array_protein(["rsd eltr".replace(" ", ""), "RSDHLTT", "RSDERKR"])
    assert a == ZIF268_PROTEIN
