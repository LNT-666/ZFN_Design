"""阶段 3：designability——先判断"设计得出来吗"，再决定要不要排队。

这一步的价值常被低估：如果指头库里根本没有某个三联体的可靠指头，
后面所有的打分都只是在给一个**造不出来的东西**打分。所以本阶段是
一道硬门槛，输出三态：

``full``
    两侧 6 个三联体全部有 Tier1 证据（实测/文献）候选。
``partial``
    能拼出阵列，但至少有一个三联体只有 Tier2（推断/未验证）候选。
``gap``
    至少有一个三联体一条候选都没有 —— 该位点直接淘汰，并列出缺哪些。
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Iterable, Mapping

from .library import FingerLibrary
from .models import EvidenceLevel, FingerRecord, TargetSite
from .sequtil import triplets_by_finger


def slots_for_halfsite(halfsite: str) -> dict[str, str]:
    """9 bp 半位点 → ``{'F1': 3bp, 'F2': 3bp, 'F3': 3bp}``（F1 = N 端）。

    18 bp 的六指半位点返回 F1–F6，见 :func:`zfn_design.sequtil.triplets_by_finger`。
    """
    return triplets_by_finger(halfsite)


def slots_for_site(site: TargetSite) -> dict[str, dict[str, str]]:
    return {
        "left": slots_for_halfsite(site.left_half),
        "right": slots_for_halfsite(site.right_half),
    }


def _slot_keys(n: int) -> list[str]:
    return [f"F{i}" for i in range(1, n + 1)]


def context_for_slot(slot: str, n_slots: int) -> str:
    """把"阵列里的第几个指头"映射到"库里的哪个上下文"。

    Princeton B1H 只测过三指阵列的**第 2 位与第 3 位**（F2 / F3）。任意指数的
    阵列都可以看成若干个三指模块首尾相接，所以按"模块内位置"取上下文：

    * 模块内第 2 指 → 用 F2 的上下文数据；
    * 模块内第 3 指 → 用 F3 的上下文数据；
    * 模块内第 1 指（N 端指）→ **库里没有本位数据**，只能靠 fallback。

    举例：三指 → F2/F3 命中第 2、3 位；六指 → 第 2、3、5、6 位命中，
    第 1、4 位走 fallback；四指 → 第 2、3 位命中，第 1、4 位走 fallback。

    这个映射是启发式的（模块内几何相似），不是实验结论。
    """
    if n_slots <= 3:
        return slot
    pos = int(slot[1:])
    within_module = ((pos - 1) % 3) + 1     # 1 / 2 / 3
    if within_module == 2:
        return "F2"
    if within_module == 3:
        return "F3"
    return slot


@dataclass
class HalfSitePlan:
    """一条半位点的可设计性。"""

    halfsite: str
    side: str                                  # left / right
    slots: dict[str, str] = field(default_factory=dict)          # F1..F3 -> triplet
    candidates: dict[str, list[FingerRecord]] = field(default_factory=dict)
    missing: list[str] = field(default_factory=list)             # 缺候选的 F 位
    best_tier: dict[str, int] = field(default_factory=dict)      # F 位 -> 最优 tier
    context_fallback: list[str] = field(default_factory=list)    # 用了非本位上下文的 F 位

    @property
    def n_fingers(self) -> int:
        return len(self.slots) if self.slots else 3

    @property
    def n_combos(self) -> int:
        if self.missing:
            return 0
        n = 1
        for f in _slot_keys(self.n_fingers):
            n *= max(1, len(self.candidates.get(f, [])))
        return n

    @property
    def status(self) -> str:
        if self.missing:
            return "gap"
        return "full" if all(t == 1 for t in self.best_tier.values()) else "partial"

    def to_dict(self) -> dict[str, Any]:
        return {
            "side": self.side,
            "halfsite": self.halfsite,
            "n_fingers": self.n_fingers,
            "slots": dict(self.slots),
            "n_candidates": {f: len(self.candidates.get(f, [])) for f in _slot_keys(self.n_fingers)},
            "best_tier": dict(self.best_tier),
            "missing": list(self.missing),
            "context_fallback": list(self.context_fallback),
            "n_combos": self.n_combos,
            "status": self.status,
        }


@dataclass
class Designability:
    site: TargetSite
    left: HalfSitePlan
    right: HalfSitePlan | None = None      # 单体型没有右半位点
    status: str = "gap"
    reasons: list[str] = field(default_factory=list)

    @property
    def halves(self) -> list[HalfSitePlan]:
        return [self.left] if self.right is None else [self.left, self.right]

    @property
    def n_fingers(self) -> int:
        return self.left.n_fingers + (self.right.n_fingers if self.right else 0)

    @property
    def n_combos(self) -> int:
        if self.right is None:
            return self.left.n_combos
        return min(self.left.n_combos, self.right.n_combos)

    @property
    def missing_triplets(self) -> list[str]:
        out = [self.left.slots[f] for f in self.left.missing]
        if self.right is not None:
            out += [self.right.slots[f] for f in self.right.missing]
        return out

    @property
    def context_fallback(self) -> list[str]:
        out = [f"L:{f}" for f in self.left.context_fallback]
        if self.right is not None:
            out += [f"R:{f}" for f in self.right.context_fallback]
        return out

    def to_dict(self) -> dict[str, Any]:
        d = self.site.to_dict()
        halves: dict[str, Any] = {"left": self.left.to_dict()}
        if self.right is not None:
            halves["right"] = self.right.to_dict()
        d.update({
            "designability": self.status,
            "n_fingers": self.n_fingers,
            "n_combos": self.n_combos,
            "missing_triplets": self.missing_triplets,
            "context_fallback": self.context_fallback,
            "reasons": list(self.reasons),
            "halves": halves,
        })
        return d


def plan_halfsite(halfsite: str, library: FingerLibrary, *, side: str = "left",
                  max_tier: int = 2) -> HalfSitePlan:
    slots = slots_for_halfsite(halfsite)
    plan = HalfSitePlan(halfsite=halfsite, side=side, slots=slots)
    n = len(slots)
    for f in _slot_keys(n):
        t = slots[f]
        ctx = context_for_slot(f, n)
        cands = library.candidates(t, split_position=ctx, max_tier=max_tier)
        if not cands:
            # 库里的实测数据只覆盖部分阵列位置（Princeton B1H 只有 F2/F3）。
            # 同一位点上"另一个上下文里的指头"仍然是该三联体的实测结合者，
            # 只是没在本位验证过——所以退一步用它，但必须记下来。
            cands = library.candidates(t, split_position=None, max_tier=max_tier)
            if cands:
                plan.context_fallback.append(f)
        plan.candidates[f] = cands
        if not cands:
            plan.missing.append(f)
        else:
            plan.best_tier[f] = min(c.evidence.tier for c in cands)
    return plan


def evaluate(site: TargetSite, library: FingerLibrary, *, max_tier: int = 2) -> Designability:
    left = plan_halfsite(site.left_half, library, side="left", max_tier=max_tier)
    right = (plan_halfsite(site.right_half, library, side="right", max_tier=max_tier)
             if site.is_paired else None)
    halves = [("左", left)] if right is None else [("左", left), ("右", right)]

    reasons: list[str] = []
    if any(p.missing for _, p in halves):
        status = "gap"
        for side_name, p in halves:
            for f in p.missing:
                reasons.append(f"{side_name}半位点 {f} 需要 {p.slots[f]}，指头库无候选")
    elif all(p.status == "full" for _, p in halves):
        status = "full"
        n = sum(len(p.slots) for _, p in halves)
        reasons.append(f"{n} 个三联体全部有 Tier1 候选")
    else:
        status = "partial"
        weak = [f"{s}:{f}:{p.slots[f]}" for s, p in halves
                for f, t in p.best_tier.items() if t > 1]
        reasons.append("以下三联体只有 Tier2（未验证）候选: " + ", ".join(sorted(weak)))

    fb = [f"{s}:{f}" for s, p in halves for f in p.context_fallback]
    if fb:
        reasons.append("以下指位库中没有本位上下文的数据，已用同一三联体在别的阵列位置实测的指头: "
                       + ", ".join(sorted(fb)))

    return Designability(site=site, left=left, right=right, status=status, reasons=reasons)


def evaluate_many(sites: Iterable[TargetSite], library: FingerLibrary,
                  *, max_tier: int = 2) -> list[Designability]:
    return [evaluate(s, library, max_tier=max_tier) for s in sites]


_STATUS_RANK = {"full": 0, "partial": 1, "gap": 2}


def triage(results: Iterable[Designability]) -> list[Designability]:
    """按可设计性排序：全 Tier1 优先，其次缺候选的排最后。"""
    return sorted(results, key=lambda d: (
        _STATUS_RANK[d.status],
        -d.n_combos,
        d.site.seq_name,
        d.site.start,
        d.site.strand,
    ))


def triplet_demand(results: Iterable[Designability], top: int = 24) -> list[dict[str, Any]]:
    """统计"输入里最需要哪些三联体"。

    这是本工具最实用的一张表：种子库覆盖不够时，用户不需要一次性导入
    全部 64 个三联体，只要按这张表的顺序补最有价值的那十几个即可。
    每行给出：三联体、出现次数、当前是否有候选、缺口次数。
    """
    total: dict[str, int] = {}
    missing_hits: dict[str, int] = {}
    for r in results:
        for plan in r.halves:
            for f in _slot_keys(plan.n_fingers):
                t = plan.slots[f]
                total[t] = total.get(t, 0) + 1
                if f in plan.missing:
                    missing_hits[t] = missing_hits.get(t, 0) + 1
    rows = [
        {
            "triplet": t,
            "occurrences": n,
            "n_candidate_slots_gap": missing_hits.get(t, 0),
            "covered": t not in missing_hits,
        }
        for t, n in total.items()
    ]
    rows.sort(key=lambda r: (r["covered"], -r["occurrences"]))
    return rows[:top]


def summary(results: Iterable[Designability]) -> dict[str, Any]:
    results = list(results)
    by_status: Mapping[str, int] = {
        s: sum(1 for r in results if r.status == s) for s in ("full", "partial", "gap")
    }
    missing: dict[str, int] = {}
    for r in results:
        for t in r.missing_triplets:
            missing[t] = missing.get(t, 0) + 1
    return {
        "n_sites": len(results),
        "by_status": dict(by_status),
        "most_often_missing_triplets": sorted(missing.items(), key=lambda kv: -kv[1])[:12],
        "n_context_fallback_sites": sum(1 for r in results if r.context_fallback),
        "triplet_demand": triplet_demand(results),
    }


__all__ = [
    "HalfSitePlan", "Designability", "plan_halfsite", "evaluate", "evaluate_many",
    "triage", "summary", "slots_for_halfsite", "slots_for_site", "EvidenceLevel",
]
