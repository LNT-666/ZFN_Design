"""阶段 5：score——给每个设计打一个**可审计**的分，而不是"预测活性"。

三条硬规矩
----------
1. **不叫亲和力预测**。每个分项都是把文献里的定性设计规则写成公式，
   分数只表示"按这些规则看，它更像一个靠谱的位点"。
2. **每个分项必须能打印出它的算法和扣分原因**。报告中给出 raw / weight /
   weighted 三列，任何人拿计算器都能复核总分。
3. **分数不能提升证据等级**。一个全 Tier2 的设计无论打多少分，
   报告里仍然标"未验证"。

分项与依据
----------
=================  ======  ============================================
分项               权重    依据
=================  ======  ============================================
evidence           0.30   指头是否有实测/文献支持（Tier1/Tier2）
offtarget          ——     已拆出：设计分不含脱靶项，全基因组脱靶由 score_offtarget 单独给分
subsite            0.18   三联体质量：GNN 偏好、CpG 甲基化风险
spacer             0.10   5-7 bp，6 最常用
gc_balance         0.07   全位点 GC 40-60% 更稳
finger_diversity   0.07   相邻指头螺旋完全相同 → 模块组装上下文风险
complexity         0.04   同碱基长串 / 串联重复 → 合成与特异性问题
self_complementarity 0.04 位点内部反向重复 → 自二聚/发夹风险
=================  ======  ============================================
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Sequence

from .models import ScoreBreakdown, ScoreComponent, TargetSite, ZfnDesign
from .offtarget import burden_to_score
from .sequtil import gc_content, is_self_complementary, max_homopolymer, tandem_repeat_count

DATA_DIR = Path(__file__).resolve().parent / "data"
WEIGHTS_PATH = DATA_DIR / "scoring_weights.json"

_EVIDENCE_POINTS = {1: 100.0, 2: 45.0, 3: 0.0}

COMPONENT_LABELS = {
    "evidence": "指头证据强度",
    "offtarget": "脱靶负担（in silico）",
    "subsite": "三联体质量",
    "spacer": "spacer 长度",
    "gc_balance": "GC 平衡",
    "finger_diversity": "指头多样性",
    "complexity": "序列复杂度",
    "self_complementarity": "自身互补风险",
}

# 阶段 1 的对外口径：这六项加权总分叫 Design_score（脱靶另算，不进这个分）。
DESIGN_SCORE_NAME = "Design_score"
DESIGN_SCORE_KEYS = ("evidence", "subsite", "gc_balance", "finger_diversity",
                     "complexity", "self_complementarity")


@lru_cache(maxsize=4)
def _load_cached(path_str: str) -> dict:
    return json.loads(Path(path_str).read_text(encoding="utf-8"))


def load_weights(path: str | Path | None = None) -> dict:
    p = Path(path).resolve() if path else WEIGHTS_PATH
    data = _load_cached(str(p))
    weights = data["weights"]
    total = sum(weights.values())
    if abs(total - 1.0) > 1e-6:
        raise ValueError(f"权重之和必须为 1.0，当前 {total:.4f}（{p}）")
    return data


# ---------------------------------------------------------------- 各分项
def _triplet_score(triplet: str, rules: dict) -> tuple[float, str]:
    t = triplet.upper()
    if t[0] == "G":
        base, why = rules["first_base_G"], "首位 G（GNN 偏好）"
    elif t[0] in "AT":
        base, why = rules["first_base_A_or_T"], f"首位 {t[0]}"
    else:
        base, why = rules["first_base_C"], "首位 C"
    notes = [why]
    if "CG" in t:
        base -= rules["penalty_contains_CG"]
        notes.append("含 CpG（甲基化风险，-%.0f）" % rules["penalty_contains_CG"])
    if len(set(t)) == 1:
        base -= rules["penalty_homopolymer_triplet"]
        notes.append("三联体同碱基（-%.0f）" % rules["penalty_homopolymer_triplet"])
    return max(0.0, base), "；".join(notes)


def _subsite_component(site: TargetSite, rules: dict) -> ScoreComponent:
    from .sequtil import finger_order  # 延迟导入，避免循环

    triplets = finger_order(site.left_half)
    if site.is_paired:
        triplets = triplets + finger_order(site.right_half)
    scores: list[float] = []
    cache: dict[str, tuple[float, str]] = {}
    counts: dict[str, int] = {}
    for t in triplets:
        if t not in cache:
            cache[t] = _triplet_score(t, rules)
        scores.append(cache[t][0])
        counts[t] = counts.get(t, 0) + 1
    raw = sum(scores) / len(scores)
    detail = " | ".join(
        f"{t}×{counts[t]}={cache[t][0]:.0f}({cache[t][1]})" for t in counts
    )
    return ScoreComponent("subsite", COMPONENT_LABELS["subsite"], raw, 0.0,
                          f"{len(triplets)} 个指位三联体（均值 {raw:.1f}）：{detail}")


def _gc_component(site: TargetSite, rules: dict) -> ScoreComponent:
    gc = gc_content(site.site_seq)
    if rules["ideal_low"] <= gc <= rules["ideal_high"]:
        raw = 100.0
    elif gc < rules["acceptable_low"] or gc > rules["acceptable_high"]:
        raw = 30.0
    else:
        raw = 70.0
    return ScoreComponent("gc_balance", COMPONENT_LABELS["gc_balance"], raw, 0.0,
                          f"GC={gc:.1%}（理想 {rules['ideal_low']:.0%}-{rules['ideal_high']:.0%}）")


def _complexity_component(site: TargetSite, rules: dict) -> ScoreComponent:
    seq = site.site_seq
    hp = max_homopolymer(seq)
    tand = tandem_repeat_count(seq, rules["tandem_unit"])
    hp_ok, hp_bad = rules["max_homopolymer_ok"], rules["max_homopolymer_bad"]
    if hp <= hp_ok and tand <= rules["tandem_count_ok"]:
        raw = 100.0
    elif hp >= hp_bad or tand >= rules["tandem_count_bad"]:
        raw = 25.0
    else:
        raw = 65.0
    return ScoreComponent("complexity", COMPONENT_LABELS["complexity"], raw, 0.0,
                          f"最长同碱基串={hp}；{rules['tandem_unit']}-mer 重复数={tand}")


def _selfcomp_component(site: TargetSite) -> ScoreComponent:
    risky = is_self_complementary(site.site_seq)
    return ScoreComponent("self_complementarity", COMPONENT_LABELS["self_complementarity"],
                          0.0 if risky else 100.0, 0.0,
                          "位点内部存在 >=6 bp 反向重复（自二聚/发夹风险）" if risky else "无内部反向重复")


def _diversity_component(design: ZfnDesign, rules: dict) -> ScoreComponent:
    helices = [f.helix7 for f in design.all_fingers]
    adj = sum(1 for a, b in zip(helices, helices[1:]) if a == b)
    raw = 100.0
    detail = "无相邻重复螺旋"
    if adj:
        raw = rules["adjacent_identical_helix"]
        detail = f"{adj} 处相邻指头螺旋完全相同（模块组装上下文风险）"
    return ScoreComponent("finger_diversity", COMPONENT_LABELS["finger_diversity"], raw, 0.0, detail)


def _evidence_component(design: ZfnDesign) -> ScoreComponent:
    fingers = design.all_fingers
    pts = [_EVIDENCE_POINTS[f.evidence.tier] for f in fingers]
    raw = sum(pts) / len(pts) if pts else 0.0
    n_t2 = sum(1 for f in fingers if f.evidence.tier == 2)
    detail = (f"{len(fingers)} 条指头：Tier1={len(fingers)-n_t2}, Tier2={n_t2}；"
              f"来源={sorted({f.source_id for f in fingers if f.source_id})}")
    return ScoreComponent("evidence", COMPONENT_LABELS["evidence"], raw, 0.0, detail)


def _offtarget_component(design: ZfnDesign, rules: dict) -> ScoreComponent:
    if design.offtarget is None:
        return ScoreComponent("offtarget", COMPONENT_LABELS["offtarget"], 50.0, 0.0,
                              "未做脱靶分析（按中间分 50 处理）")
    burden = design.offtarget.burden
    raw = burden_to_score(burden, scale=rules["burden_scale"])
    s = design.offtarget.summary()
    return ScoreComponent("offtarget", COMPONENT_LABELS["offtarget"], raw, 0.0,
                          f"burden={burden:.3f}；命中 {s['n_hits']} 个（exact {s['n_exact']}）；"
                          f"搜索 {s['searched_bp']} bp")


# ---------------------------------------------------------------- 主入口
def score_offtarget(design: ZfnDesign, weights_data: dict | None = None) -> dict[str, Any]:
    """脱靶负担单独给分（0-100），**不并进** Design score。

    只有做过脱靶搜索（局部 or 全基因组）时才有分；没做时返回 ``score=None``。
    """
    data = weights_data or load_weights()
    rules = data["rules"]["offtarget"]
    if design.offtarget is None:
        return {"score": None, "burden": None, "n_hits": None,
                "scale": rules["burden_scale"], "detail": "未做脱靶分析"}
    c = _offtarget_component(design, rules)
    return {"score": round(c.raw, 2), "burden": round(design.offtarget.burden, 4),
            "n_hits": design.offtarget.summary()["n_hits"], "scale": rules["burden_scale"],
            "detail": c.detail}


def grade_of(total: float, grades: dict) -> str:
    if total >= grades["A"]:
        return "A"
    if total >= grades["B"]:
        return "B"
    if total >= grades["C"]:
        return "C"
    return "D"


def score_design(design: ZfnDesign, weights_data: dict | None = None) -> ScoreBreakdown:
    data = weights_data or load_weights()
    w = data["weights"]
    r = data["rules"]

    components = [
        _evidence_component(design),
        _subsite_component(design.target, r["subsite"]),
        _gc_component(design.target, r["gc_balance"]),
        _diversity_component(design, r["finger_diversity"]),
        _complexity_component(design.target, r["complexity"]),
        _selfcomp_component(design.target),
    ]
    # spacer 只对组合型有意义：单体型没有 spacer，也没有二聚化，
    # 所以这一项不进分，其余权重按比例重归一化（总分仍是 0-100）。
    if design.is_paired:
        components.insert(3, ScoreComponent(
            "spacer", COMPONENT_LABELS["spacer"],
            r["spacer"].get(f"len_{design.target.spacer_len}", r["spacer"]["len_5_or_7"]),
            0.0, f"spacer={design.target.spacer}（{design.target.spacer_len} bp）"))

    for c in components:
        c.weight = w[c.key]

    weight_sum = sum(c.weight for c in components)
    renormalized = abs(weight_sum - 1.0) > 1e-9
    if renormalized:
        for c in components:
            c.weight = c.weight / weight_sum

    total = sum(c.contribution for c in components)

    warnings: list[str] = []
    if any(f.evidence.tier > 1 for f in design.all_fingers):
        warnings.append("该设计含 Tier2（未验证）指头——必须实验验证后才能使用")
    if design.offtarget is None:
        warnings.append("未做脱靶分析：设计分不含脱靶项，脱靶负担见单独的 Off-target score")
    elif design.offtarget.exact_hits > 0:
        warnings.append(f"在搜索空间内发现 {design.offtarget.exact_hits} 个 0 错配位点，需人工核查")
    if not design.is_paired:
        warnings.append("单体型（ZF-TF/单阵列）：无 FokI 二聚化，不产生双链断裂；"
                        "脱靶判据是「一个 9-18 bp 半位点像不像」，与组合型不可直接比较")
    warnings.append("所有计算输出都必须实验验证；本分数不是活性预测")

    notes = [
        "总分 = Σ(raw × weight)，权重见 scoring_weights.json",
        "分数只反映'按公开设计规则看更像靠谱位点'，与实测活性无定量关系",
    ]
    dropped = sorted(set(w) - {c.key for c in components})
    if dropped:
        # 单体型没有 spacer/二聚化：该项不进分母。若被剔除项权重本来就是 0，
        # 总分不变，但「按比例重归一化」这句话仍然成立（分母 = 其余权重之和）。
        notes.append(f"单体型不适用 {'/'.join(dropped)} 分项（该项权重为 0 时总分不变），"
                     f"其余权重按比例重归一化")
    return ScoreBreakdown(components=components, total=total,
                          grade=grade_of(total, data["grades"]),
                          warnings=warnings, notes=notes)


def rank(designs: Iterable[ZfnDesign], weights_data: dict | None = None) -> list[ZfnDesign]:
    data = weights_data or load_weights()
    out = []
    for d in designs:
        d.score = score_design(d, data)
        out.append(d)
    out.sort(key=lambda d: (-(d.score.total if d.score else 0.0), d.design_id))
    return out


def weights_documentation(weights_data: dict | None = None) -> dict[str, Any]:
    data = weights_data or load_weights()
    return {
        "weights": data["weights"],
        "grades": data["grades"],
        "rules": data["rules"],
        "labels": COMPONENT_LABELS,
        "design_score_name": DESIGN_SCORE_NAME,
        "design_score_keys": list(DESIGN_SCORE_KEYS),
    }


__all__ = ["load_weights", "score_design", "score_offtarget", "rank", "grade_of",
           "weights_documentation", "COMPONENT_LABELS",
           "DESIGN_SCORE_NAME", "DESIGN_SCORE_KEYS"]
