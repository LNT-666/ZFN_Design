"""S3 · 全基因组脱靶负担（off-target score）。

口径（宋成义 2026-11-18）
------------------------
对 S1/S2 出来的设计，取 Design_score 前 N 个（默认 10），每个设计的
左/右半位点 + spacer 拿去扫整个基因组，每侧允许 0/1/2 错配，数出
近似位点数与加权 burden，再折成 0–100 的 off-target score。
**两套分并列展示，不相加**：off-target score 不进 Design_score。

为什么不直接用 :func:`zfn_design.offtarget_genome.genome_search`
---------------------------------------------------------------
可以直接用，本模块就是它的薄封装：把设计的半位点、spacer、错配上限拆成
参数，逐个跑、汇总成报告能直接吃的结构。真正的扫描逻辑仍在那一个函数里，
避免两份实现漂移。

成本
----
每扫一遍全基因组（3.1 Gb）约 60–120 s（numpy 路径，2 个朝向 × 1 个 spacer）；
10 个设计就是 10–20 分钟。所以默认 ``top=10``，且必须显式给 genome_fasta。
"""

from __future__ import annotations

import time
from typing import Callable, Iterable

from . import offtarget_genome
from .models import ZfnDesign

__all__ = ["select_top", "score_designs", "score_from_burden", "summary"]

score_from_burden = offtarget_genome.score_from_burden


def select_top(designs: Iterable[ZfnDesign], top: int) -> list[ZfnDesign]:
    """按 Design_score 总分降序取前 ``top`` 个（``top<=0`` 表示全取）。"""
    ranked = sorted(
        designs,
        key=lambda d: (d.score.total if d.score is not None else 0.0),
        reverse=True,
    )
    return ranked if top <= 0 else ranked[:top]


def score_designs(
    designs: Iterable[ZfnDesign],
    genome_fasta: str,
    *,
    top: int = 10,
    max_mm: int = 2,
    max_hits: int = 500,
    chunk: int = 4_000_000,
    progress: Callable[[str], None] | None = None,
    on_progress: Callable[[str, float, float, str], None] | None = None,
) -> list[dict]:
    """逐个设计做全基因组脱靶扫描，返回报告用的记录（按 burden 升序）。

    每条记录的 ``burden`` 越小越好，``offtarget_score`` 是
    ``100 - scale * burden``（scale 见 scoring_weights.json 的 offtarget 段）。

    ``on_progress(phase, done, total, message)`` 给前端画进度条用：``done``
    是"已完成设计数 + 当前设计扫过的碱基比例"，``total`` 是设计总数。
    """
    picked = select_top(designs, top)
    out: list[dict] = []
    n = len(picked)
    for i, design in enumerate(picked, 1):
        tgt = design.target
        left, right = tgt.left_half, tgt.right_half
        if progress:
            progress(f"[{i}/{n}] {design.design_id} {tgt.seq_name}:{tgt.start}{tgt.strand}")
        if on_progress:
            on_progress("全基因组脱靶扫描", i - 1, n, f"{design.design_id} 起扫")
        t0 = time.perf_counter()
        res = offtarget_genome.genome_search(
            genome_fasta,
            left,
            right,
            spacers=(tgt.spacer_len,),
            max_mm_left=max_mm,
            max_mm_right=max_mm,
            max_hits=max_hits,
            chunk=chunk,
            on_bp=(
                (lambda done_bp, total_bp: on_progress(
                    "全基因组脱靶扫描",
                    (i - 1) + (done_bp / total_bp if total_bp else 0.0),
                    n,
                    f"{design.design_id} {done_bp:,}/{total_bp:,} bp",
                ))
                if on_progress else None
            ),
        )
        score = score_from_burden(res["burden"])
        rec = {
            "design_id": design.design_id,
            "rank": i,
            "seq_name": tgt.seq_name,
            "start": tgt.start,
            "strand": tgt.strand,
            "spacer": tgt.spacer_len,
            "max_mm": max_mm,
            "left_half": left,
            "right_half": right,
            "design_score": round(design.score.total, 2) if design.score else None,
            "offtarget_score": round(score, 2),
            "burden": res["burden"],
            "by_mm": {str(k): v for k, v in sorted(res["by_mm"].items())},
            "n_hits_kept": res["n_hits"],
            "n_hits_total": res["n_total_hits"],
            "truncated": res["truncated"],
            "engine": res["engine"],
            "scanned_bp": res["scanned_bp"],
            "elapsed_s": round(time.perf_counter() - t0, 1),
            "hits": res["hits"],
        }
        design.offtarget_score = {
            "score": rec["offtarget_score"],
            "burden": rec["burden"],
            "by_mm": rec["by_mm"],
            "n_hits": rec["n_hits_total"],
            "engine": rec["engine"],
            "scanned_bp": rec["scanned_bp"],
        }
        out.append(rec)
    out.sort(key=lambda r: r["burden"])
    for i, rec in enumerate(out, 1):
        rec["burden_rank"] = i
    return out


def summary(records: list[dict]) -> dict:
    """报告顶部的汇总卡。"""
    if not records:
        return {"n_designs": 0}
    burdens = [r["burden"] for r in records]
    scores = [r["offtarget_score"] for r in records]
    return {
        "n_designs": len(records),
        "burden_min": min(burdens),
        "burden_max": max(burdens),
        "burden_median": sorted(burdens)[len(burdens) // 2],
        "score_min": min(scores),
        "score_max": max(scores),
        "engine": records[0]["engine"],
        "scanned_bp": records[0]["scanned_bp"],
        "by_mm_total": _sum_mm(records),
    }


def _sum_mm(records: list[dict]) -> dict[str, int]:
    total: dict[str, int] = {}
    for r in records:
        for k, v in r["by_mm"].items():
            total[k] = total.get(k, 0) + int(v)
    return {k: total[k] for k in sorted(total, key=lambda x: int(x))}
