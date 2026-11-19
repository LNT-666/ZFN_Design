"""把三条识别螺旋装回蛋白骨架，得到可合成的氨基酸序列。

这是本方案里唯一"输出蛋白序列"的地方，因此必须极其保守：

* 骨架模板来自 **PDB 1AAY**（Zif268, Mus musculus）真实序列，逐字符抄录；
* 只做 backbone-swap：保留骨架与连接区，只替换识别螺旋的 7 个残基；
* 若某个螺旋的第 4 位是不接触 DNA 的占位符 ``X``，输出里会原样保留，
  并在 :func:`warnings_for` 里给出提醒——绝不悄悄猜一个残基填进去。

自检：把 ``RSDELTR/RSDHLTT/RSDERKR`` 填回模板，必须**一字不差**地还原
1AAY 的序列。这条断言在 ``tests/test_framework.py`` 里。
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Sequence

DATA_DIR = Path(__file__).resolve().parent / "data"
TEMPLATES_PATH = DATA_DIR / "framework_templates.json"

AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
#: 位置 4 不接触 DNA，允许用 X 占位
AA_ALLOW_X = AA | {"X"}

# ---- 任意指数阵列的通用骨架（构造体，不是晶体结构） ---------------------
# 全部片段逐字符取自 1AAY，只有"模块之间怎么接"是构造的：
#   N 端      MERPYACPVESCDRRFS
#   F1→F2     HIRIHTGQKPFQCRICMRNFS     （1AAY 原生连接区，13+8 aa）
#   F2→F3     HIRTHTGEKPFACDICGRKFA     （1AAY 原生连接区，13+8 aa）
#   F3→下一模块 TGEKP + 去掉 Met 的 N 端  （和 zif268_x2_junction 用的是同一个接法）
#   C 端      HTKIHLRQKD
GENERIC_KEY = "repeat_cassette"
_G_NTERM = "MERPYACPVESCDRRFS"
_G_TERM = "HTKIHLRQKD"
_G_LINK = ("HIRIHTGQKPFQCRICMRNFS", "HIRTHTGEKPFACDICGRKFA")
_G_MODULE_START = "TGEKP" + _G_NTERM[1:]      # 第二个模块起（无起始 Met）
MIN_FINGERS, MAX_FINGERS = 2, 12

GENERIC_META = {
    "key": GENERIC_KEY,
    "label": "1AAY 单指区段重复拼接（任意指数）",
    "n_fingers": None,          # 运行时决定
    "source_id": "PDB:1AAY ×k",
    "evidence_level": "constructed",
    "notes": (
        "由 1AAY 的 N 端、两种原生连接区与 C 端逐字符拼成，按三指一模块重复；"
        "**模块间界面、连接肽长度、相对取向都没有任何结构或活性证据**，"
        "只能当可合成的草稿，必须单独验证。"
    ),
}


def generic_sequence_template(n_fingers: int) -> str:
    """生成 n 指阵列的骨架模板字符串（``{h1}``…``{hn}`` 插槽）。"""
    if not (MIN_FINGERS <= n_fingers <= MAX_FINGERS):
        raise FrameworkError(
            f"指头数必须在 {MIN_FINGERS}–{MAX_FINGERS} 之间，收到 {n_fingers}")
    parts = [_G_NTERM]
    for i in range(1, n_fingers + 1):
        parts.append(f"{{h{i}}}")
        if i == n_fingers:
            break
        within = ((i - 1) % 3) + 1
        parts.append(_G_LINK[within - 1] if within < 3 else _G_MODULE_START)
    parts.append(_G_TERM)
    return "".join(parts)


class FrameworkError(ValueError):
    pass


@lru_cache(maxsize=1)
def load_templates() -> dict:
    return json.loads(TEMPLATES_PATH.read_text(encoding="utf-8"))


def template_keys() -> list[str]:
    return list(load_templates()["templates"].keys())


def get_template(key: str | None = None) -> dict:
    db = load_templates()
    key = key or db["default"]
    if key not in db["templates"]:
        raise FrameworkError(f"未知骨架模板 {key!r}；可用: {list(db['templates'])}")
    return db["templates"][key]


def template_for_fingers(n_fingers: int) -> str:
    """按指头数挑模板：3 指用 1AAY，6 指用两个模块串联的，其它指数用通用拼接。"""
    for key, tpl in load_templates()["templates"].items():
        if tpl["n_fingers"] == n_fingers:
            return key
    if MIN_FINGERS <= n_fingers <= MAX_FINGERS:
        return GENERIC_KEY
    raise FrameworkError(
        f"没有 {n_fingers} 指的骨架模板（支持 {MIN_FINGERS}–{MAX_FINGERS} 指）；现有固定模板: "
        + ", ".join(f"{k}({v['n_fingers']} 指)" for k, v in load_templates()["templates"].items())
    )


def _check_helix(helix7: str, index: int) -> str:
    h = (helix7 or "").strip().upper()
    if len(h) != 7:
        raise FrameworkError(f"第 {index} 个识别螺旋必须 7 残基，得到 {helix7!r}")
    bad = sorted(set(h) - AA_ALLOW_X)
    if bad:
        raise FrameworkError(f"第 {index} 个识别螺旋含非氨基酸字符 {bad!r}")
    return h


def build_array_protein(helices: Sequence[str], template_key: str | None = None) -> str:
    """``helices`` 必须是 N→C 顺序的 7 残基螺旋（任意指数，2–12）。

    ``template_key`` 省略时按螺旋条数自动挑模板：3 指→1AAY 真实骨架，
    6 指→双模块串联，其它→通用拼接（构造体）。
    """
    if not helices:
        raise FrameworkError("至少要给一条识别螺旋")
    if template_key is None:
        template_key = template_for_fingers(len(helices))
    if template_key == GENERIC_KEY:
        seq = generic_sequence_template(len(helices))
        for i in range(1, len(helices) + 1):
            seq = seq.replace(f"{{h{i}}}", _check_helix(helices[i - 1], i))
        return seq
    tpl = get_template(template_key)
    n = tpl["n_fingers"]
    if len(helices) != n:
        raise FrameworkError(f"模板 {template_key or 'default'} 需要 {n} 条螺旋，得到 {len(helices)}")
    seq = tpl["sequence_template"]
    for i, slot in enumerate(tpl["helix_slots"], start=1):
        seq = seq.replace(slot, _check_helix(helices[i - 1], i))
    return seq


def warnings_for(helices: Sequence[str]) -> list[str]:
    """输出蛋白序列时需要提醒用户的事项。"""
    out = []
    for i, h in enumerate(helices, start=1):
        if "X" in h.upper():
            out.append(f"指 {i} 的识别螺旋含占位符 X（第 4 位，不接触 DNA）：合成前请替换为真实残基")
    return out


def describe(template_key: str | None = None) -> dict:
    if template_key == GENERIC_KEY:
        return dict(GENERIC_META, template_length=None,
                    template_note=f"按指头数动态生成，长度 = 28n + 6（n 为指头数）")
    tpl = get_template(template_key)
    return {
        "key": template_key or load_templates()["default"],
        "label": tpl["label"],
        "n_fingers": tpl["n_fingers"],
        "source_id": tpl["source_id"],
        "evidence_level": tpl["evidence_level"],
        "template_length": len(tpl["sequence_template"]) - sum(len(s) for s in tpl["helix_slots"]) + 7 * tpl["n_fingers"],
        "notes": tpl["notes"],
    }


__all__ = ["FrameworkError", "GENERIC_KEY", "MIN_FINGERS", "MAX_FINGERS",
           "build_array_protein", "generic_sequence_template", "warnings_for", "describe",
           "template_keys", "template_for_fingers", "get_template", "load_templates"]
