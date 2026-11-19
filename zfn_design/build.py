"""阶段 4：build——把一个"可设计"的位点变成若干个具体的 ZFN 设计。

输入是可设计性分析（已经确认每个三联体都有候选），输出是排好序的
:class:`~zfn_design.models.ZfnDesign` 列表。一个"设计"= 一对单体
（左阵列 + 右阵列），因为 ZFN 只有成对且方向正确时才能切 DNA。

排序原则（在设计阶段就定，不打分）：
1. 每个指位取可靠度最高的前 ``max_per_finger`` 个候选；
2. 组合按"指标等级之和"升序——Tier1 越多的组合越靠前；
3. 其余交给 score 阶段，避免这里重复实现一遍权重逻辑。
"""

from __future__ import annotations

from itertools import product
from typing import Iterator, Sequence

from . import framework
from .designability import Designability, HalfSitePlan
from .library import FingerLibrary
from .models import FingerRecord, TargetSite, ZfnDesign
from .sequtil import revcomp


class BuildError(ValueError):
    pass


def _enumerate_arrays(plan: HalfSitePlan, *, max_per_finger: int = 2,
                      max_tier: int = 2, library: FingerLibrary | None = None
                      ) -> Iterator[list[FingerRecord]]:
    if plan.missing:
        return
    pools: list[Sequence[FingerRecord]] = []
    for f in [f"F{i}" for i in range(1, plan.n_fingers + 1)]:
        cands = plan.candidates.get(f, [])
        pool = [c for c in cands if c.evidence.tier <= max_tier][:max_per_finger]
        if not pool:
            return
        pools.append(pool)
    for combo in product(*pools):
        yield list(combo)


def _tier_sum(fingers: Sequence[FingerRecord]) -> int:
    return sum(f.evidence.tier for f in fingers)


def _protein(fingers: Sequence[FingerRecord], template_key: str | None) -> tuple[str, list[str]]:
    try:
        seq = framework.build_array_protein([f.helix7 for f in fingers], template_key)
    except framework.FrameworkError:
        return "", ["骨架装配失败：识别螺旋含非法字符"]
    return seq, framework.warnings_for([f.helix7 for f in fingers])


def build_for_site(
    site: TargetSite,
    library: FingerLibrary,
    designability: Designability,
    *,
    n_designs: int = 5,
    max_per_finger: int = 2,
    max_tier: int = 2,
    template_key: str | None = None,
    architecture: str = "paired_nuclease",
) -> list[ZfnDesign]:
    """为单个位点生成最多 ``n_designs`` 个设计。

    组合型产出成对设计（左单体 + 右单体）；单体型只产出**一条**阵列的设计。
    """
    if designability.status == "gap":
        return []

    left_arrays = list(_enumerate_arrays(designability.left, max_per_finger=max_per_finger,
                                        max_tier=max_tier, library=library))
    if not left_arrays:
        return []

    if not site.is_paired:
        singles = sorted(left_arrays, key=_tier_sum)[:n_designs]
        designs: list[ZfnDesign] = []
        for idx, la in enumerate(singles, start=1):
            protein, warns = _protein(la, template_key)
            designs.append(
                ZfnDesign(
                    design_id=f"{site.site_id}__d{idx}",
                    target=site,
                    fingers=list(la),
                    protein=protein,
                    architecture=architecture,
                    notes=list(warns),
                )
            )
        return designs

    right_arrays = list(_enumerate_arrays(designability.right, max_per_finger=max_per_finger,
                                         max_tier=max_tier, library=library))
    if not right_arrays:
        return []

    pairs: list[tuple[int, list[FingerRecord], list[FingerRecord]]] = []
    for la in left_arrays:
        for ra in right_arrays:
            pairs.append((_tier_sum(la) + _tier_sum(ra), la, ra))
    pairs.sort(key=lambda t: t[0])

    designs = []
    for idx, (_, la, ra) in enumerate(pairs[:n_designs], start=1):
        protein, warns = _protein(la, template_key)
        partner_protein, pwarns = _protein(ra, template_key)
        notes = list(warns) + list(pwarns)
        notes.append(f"右单体按反向互补结合；其上链形态 = {revcomp(site.right_half)}")
        designs.append(
            ZfnDesign(
                design_id=f"{site.site_id}__d{idx}",
                target=site,
                fingers=list(la),
                partner_fingers=list(ra),
                protein=protein,
                partner_protein=partner_protein,
                architecture=architecture,
                notes=notes,
            )
        )
    return designs


def build_all(
    sites: Sequence[TargetSite],
    library: FingerLibrary,
    designabilities: Sequence[Designability],
    *,
    n_designs: int = 5,
    max_per_finger: int = 2,
    max_tier: int = 2,
    template_key: str | None = None,
    architecture: str = "paired_nuclease",
    max_sites: int | None = None,
) -> list[ZfnDesign]:
    by_id = {d.site.site_id: d for d in designabilities}
    out: list[ZfnDesign] = []
    used = 0
    for site in sites:
        dz = by_id.get(site.site_id)
        if dz is None or dz.status == "gap":
            continue
        out.extend(build_for_site(site, library, dz, n_designs=n_designs,
                                  max_per_finger=max_per_finger, max_tier=max_tier,
                                  template_key=template_key, architecture=architecture))
        used += 1
        if max_sites is not None and used >= max_sites:
            break
    return out


__all__ = ["BuildError", "build_for_site", "build_all"]
