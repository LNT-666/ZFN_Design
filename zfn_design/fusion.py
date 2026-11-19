"""阶段 7：fusion design——把"结合"变成"有功能"。

锌指阵列本身只是"手"。要变成能用的工具，C 端（或 N 端）必须挂上效应域：

===============  ==========================  ==========================
架构             组成                         何时用
===============  ==========================  ==========================
paired_nuclease  NLS + ZF + linker + FokI     要切割（基因敲除、敲入）
single_activator NLS + ZF + VP64              要激活（把基因打开）
single_repressor NLS + ZF + KRAB              要抑制（把基因关掉）
===============  ==========================  ==========================

前一种需要**成对位点**（两个单体夹住 5-7 bp spacer，FokI 才能二聚化）；
后两种只需要单个 9 bp 半位点，因此可用位点更多、也更容易挑到安全的。

诚实边界：本模块只负责"按架构把元件接起来 + 标出每个元件的证据状态"。
凡是没有可靠序列的元件（如 FokI 催化域），输出里保留 ``[占位符]``，
绝不凭记忆填一段序列——那是最容易造成实验事故的地方。
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any

from .models import ZfnDesign

DATA_DIR = Path(__file__).resolve().parent / "data"
PARTS_PATH = DATA_DIR / "fusion_parts.json"

PLACEHOLDER_FMT = "[{name}:{status}]"


class FusionError(ValueError):
    pass


@lru_cache(maxsize=1)
def load_parts_db() -> dict:
    return json.loads(PARTS_PATH.read_text(encoding="utf-8"))


def architectures() -> dict[str, dict]:
    return load_parts_db()["architectures"]


def parts() -> dict[str, dict]:
    return load_parts_db()["parts"]


def describe_architectures() -> list[dict[str, Any]]:
    out = []
    for key, a in architectures().items():
        out.append({
            "key": key,
            "label": a["label"],
            "description": a["description"],
            "order": a["order"],
            "requires_pair": a["dimer"],
            "requires": a["requires"],
        })
    return out


def build_fusion(design: ZfnDesign, arch_key: str = "single_activator",
                 *, zf_protein: str | None = None,
                 extra_parts: dict[str, str] | None = None,
                 part_overrides: dict[str, str] | None = None) -> dict[str, Any]:
    """按架构拼出一条（草案）蛋白序列，并列出所有红旗。"""
    archs = architectures()
    if arch_key not in archs:
        raise FusionError(f"未知架构 {arch_key!r}；可用 {list(archs)}")
    arch = archs[arch_key]
    pdb = parts()
    extra_parts = extra_parts or {}
    # token -> 实际元件键（成对核酸酶要把左右两条链换成不同的 FokI 变体）
    overrides = part_overrides or {}

    # 架构与设计的形态必须一致：单体型设计不能拼 FokI 二聚体架构（没有第二个
    # 单体可以二聚化）。反过来，成对设计的**单侧**单体拿去做 ZF-TF 是可行的
    # ——那正是单体型策略，但必须标明"唯一性只按一个半位点算"。
    if arch["dimer"] and not design.target.is_paired:
        raise FusionError(
            f"架构 {arch_key} 需要成对位点（FokI 二聚化），但这是单体型设计")
    single_from_pair = (not arch["dimer"]) and design.target.is_paired

    zf = zf_protein if zf_protein is not None else design.protein
    if not zf:
        raise FusionError("该设计没有可用的 ZF 阵列蛋白序列（骨架装配失败？）")

    seq_chunks: list[str] = []
    used: list[dict[str, Any]] = []
    flags: list[str] = []
    if single_from_pair:
        flags.append("这是成对设计的单侧单体：单独做成 ZF-TF 后只剩一个半位点的特异性，"
                     "不再有 FokI 二聚化带来的额外约束")

    for token in arch["order"]:
        if token == "ZF_ARRAY":
            seq_chunks.append(zf)
            used.append({
                "token": token,
                "label": "ZF 阵列（3 指）",
                "sequence": zf,
                "evidence_level": "derived",
                "source_id": "framework template (PDB:1AAY backbone)",
                "length": len(zf),
            })
            continue

        if token == "LINKER":
            lp = pdb.get("linker_flex") or pdb.get("linker_GGGGS") or {}
            s = (extra_parts.get("linker_flex") or extra_parts.get("linker_GGGGS")
                 or lp.get("sequence") or "")
            seq_chunks.append(s)
            used.append({"token": token, "label": lp.get("label", "linker"),
                         "sequence": s, "evidence_level": lp.get("evidence_level", ""),
                         "source_id": lp.get("source_id", ""), "length": len(s)})
            continue

        part_key = overrides.get(token, token)
        part = pdb.get(part_key)
        if part is None:
            flags.append(f"架构要求的元件 {part_key} 在 fusion_parts.json 里没有定义")
            seq_chunks.append(PLACEHOLDER_FMT.format(name=part_key, status="undefined"))
            used.append({"token": part_key, "label": part_key, "sequence": None,
                         "evidence_level": "missing", "source_id": "", "length": 0})
            continue

        supplied = (extra_parts.get(part_key) or extra_parts.get(token)
                    or (extra_parts.get("FokI_CD") if part_key.startswith("FokI_CD") else None))
        s = supplied or part.get("sequence")
        if supplied:
            ev = "user_supplied"
        else:
            ev = part.get("evidence_level", "")
        if not s:
            flags.append(f"{part_key} 没有内置序列——请通过 --parts 传入你实验室验证过的序列")
            seq_chunks.append(PLACEHOLDER_FMT.format(name=part_key, status=ev or "missing"))
            used.append({"token": part_key, "label": part.get("label", part_key), "sequence": None,
                         "evidence_level": ev, "source_id": part.get("source_id", ""), "length": 0})
            continue
        if ev == "needs_review":
            flags.append(f"{part_key} 序列需逐碱基核对（evidence_level=needs_review）")
        seq_chunks.append(s)
        used.append({"token": part_key, "label": part.get("label", part_key), "sequence": s,
                     "evidence_level": ev, "source_id": part.get("source_id", ""), "length": len(s)})

    complete = all(u["sequence"] is not None for u in used)
    if not arch["dimer"] and design.target.spacer_len not in (5, 6, 7):
        flags.append("单体型架构不使用 spacer 约束，但当前设计是按成对位点挑的")

    protein = "".join(seq_chunks)
    # 真核表达要有起始密码子：蛋白 N 端不是 M 时补一个 Met。
    # NLS 在最前面（PKKKRKV 起）时正是这种情况。
    if complete and protein and not protein.startswith("M"):
        protein = "M" + protein
        used.insert(0, {
            "token": "START_MET", "label": "起始甲硫氨酸（提供 ATG）",
            "sequence": "M", "evidence_level": "required",
            "source_id": "start codon requirement", "length": 1,
        })
        flags.append("N 端补了 Met（起始密码子 ATG）；若实验室构型不同请自行调整")

    has_nls = any(u["token"].startswith("NLS") for u in used)
    if not has_nls:
        flags.append("这条融合蛋白没有 NLS，进不了核")

    return {
        "design_id": design.design_id,
        "architecture": arch_key,
        "architecture_label": arch["label"],
        "requires_pair": arch["dimer"],
        "parts": used,
        "draft_protein": protein,
        "is_complete": complete,
        "length": len(protein),
        "has_nls": has_nls,
        "red_flags": flags,
        "notes": [
            "draft 里出现 [NAME:status] 的方括号占位符 = 该元件没有可靠序列，不能直接拿去合成",
            "单体型架构（activator/repressor）对每个半位点都适用，不要求成对",
            "蛋白以 M 起始（ATG），CDS 末尾自动带终止密码子",
        ],
    }


def build_fusions(designs: list[ZfnDesign], arch_keys: list[str] | None = None,
                  *, extra_parts: dict[str, str] | None = None) -> list[dict[str, Any]]:
    keys = arch_keys or list(architectures())
    out = []
    for d in designs:
        for k in keys:
            try:
                out.append(build_fusion(d, k, extra_parts=extra_parts))
            except FusionError as exc:
                out.append({"design_id": d.design_id, "architecture": k, "error": str(exc)})
    return out


__all__ = ["FusionError", "architectures", "parts", "describe_architectures",
           "build_fusion", "build_fusions", "load_parts_db"]
