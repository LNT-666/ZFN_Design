"""S3 · 全基因组脱靶评分，单独按键执行（宋济 2026-11-18）。

为什么单独成一步
----------------
Design Score（S1 + 可选 S2）只出设计与打分，几秒钟到几分钟；
全基因组脱靶每个设计要扫 31 亿碱基，2–3 分钟。混在同一次点击里，
按钮会长时间没有反馈，也没法只重算脱靶。所以拆成两个按键：

* 「Design Score」→ :func:`zfn_design.pipeline.run`（S1/S2）
* 「Off-target Score」→ 本模块，读已有 report.json，只补 S3

产物
----
``offtarget_genome.tsv``、更新后的 ``report.json`` / ``report.html``，
以及重打的 ``zfn_bundle.zip``。
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Callable

from . import genome_offtarget, pipeline, report

__all__ = ["STAGE", "STAGE_TITLE", "REPORT_KIND", "run_offtarget", "load_payload"]

STAGE = "S3"
STAGE_TITLE = "Stage 3 · Off-target score（全基因组脱靶评分）"
REPORT_KIND = "offtarget_score"

_BUNDLE_FILES = (
    "report.json",
    "report.html",
    "sites.tsv",
    "designs.tsv",
    "offtarget_genome.tsv",
    "fusions.json",
    "proteins.tsv",
    "cds.tsv",
    "cds.fasta",
    "order_form.csv",
)


def load_payload(out_dir: str | Path) -> dict[str, Any]:
    p = Path(out_dir) / "report.json"
    if not p.exists():
        raise ValueError(f"这个目录里没有 report.json：{out_dir}")
    with p.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _rebuild_bundle(out: Path) -> Path:
    """把目录里现成的产物重打一次 zip（只收存在且非空的文件）。"""
    files: dict[str, Path] = {}
    for name in _BUNDLE_FILES:
        p = out / name
        if p.exists() and p.stat().st_size > 0:
            files[name] = p
    return report.write_bundle(out / pipeline.BUNDLE_NAME, files)


def _shim(design: dict[str, Any]) -> SimpleNamespace:
    """把 report.json 里的设计字典包成 score_designs 认的形状。"""
    tgt = design.get("target") or {}
    sc = design.get("score") or {}
    return SimpleNamespace(
        design_id=design.get("design_id") or design.get("id"),
        target=SimpleNamespace(
            seq_name=tgt.get("seq_name"),
            start=tgt.get("start"),
            strand=tgt.get("strand"),
            left_half=tgt.get("left_half"),
            right_half=tgt.get("right_half"),
            spacer_len=int(tgt.get("spacer_len") or len(tgt.get("spacer") or "")),
        ),
        score=SimpleNamespace(total=float(sc.get("total") or 0.0)),
        offtarget_score=None,
    )


def run_offtarget(
    out_dir: str | Path,
    *,
    top: int = 10,
    max_mm: int = 2,
    max_hits: int = 500,
    chunk: int = 4_000_000,
    progress: Callable[[str, float, float, str], None] | None = None,
) -> dict[str, Any]:
    """对已有作业的 Design_score 前 ``top`` 个设计补算 off-target score。"""
    out = Path(out_dir)
    payload = load_payload(out)
    meta = payload.setdefault("meta", {})
    params = meta.get("params") or {}
    genome_fasta = params.get("genome_fasta")
    if not genome_fasta:
        raise ValueError(
            "这次作业没带基因组 FASTA，无法做全基因组脱靶；请先勾选基因组重跑 Design Score"
        )

    # 宋济 2026-11-18：只有通过 Target filter（全基因组重复过滤）的候选，才能做 S3。
    gf = (payload.get("targets_summary") or {}).get("genome_filter")
    if not gf:
        raise ValueError(
            "只有通过 Target filter（全基因组重复过滤）的候选才能做 Off-target Score："
            "请先在 Target filter 里选基因组 FASTA（可选再选 RepeatMasker 注释）重跑 Design Score"
        )

    designs = payload.get("designs") or []
    if not designs:
        raise ValueError("这次作业没有设计（候选可能被过滤器清光了），S3 无从评分")

    picked = designs if not top or top <= 0 else sorted(
        designs, key=lambda d: float((d.get("score") or {}).get("total") or 0.0), reverse=True
    )[:top]
    shims = [_shim(d) for d in picked]

    if progress:
        progress("全基因组脱靶扫描", 0, len(shims), f"{len(shims)} 个设计")

    records = genome_offtarget.score_designs(
        shims,
        genome_fasta,
        top=0,
        max_mm=max_mm,
        max_hits=max_hits,
        chunk=chunk,
        on_progress=progress,
    )
    summ = genome_offtarget.summary(records)

    by_id = {r["design_id"]: r for r in records}
    for d in designs:
        r = by_id.get(d.get("design_id") or d.get("id"))
        if r is None:
            continue
        d["offtarget_score"] = {
            "score": r["offtarget_score"],
            "burden": r["burden"],
            "by_mm": r["by_mm"],
            "n_hits": r["n_hits_total"],
            "engine": r["engine"],
            "scanned_bp": r["scanned_bp"],
        }

    meta = payload.setdefault("meta", {})
    meta["stage"] = STAGE
    meta["stage_title"] = STAGE_TITLE
    meta["report_kind"] = REPORT_KIND
    payload["offtarget_genome"] = {"summary": summ, "designs": records}
    # 宋济 2026-11-18：S3 的结果也要能接着做 Fusion design（S4）。
    payload["fusion_ready"] = True

    if progress:
        progress("写报告与打包", 0, 3, "")
    report.write_json(out / "report.json", payload)
    report.write_html(out / "report.html", payload)
    pipeline._write_offtarget_genome_tsv(out / "offtarget_genome.tsv", records)
    bundle = _rebuild_bundle(out)
    if progress:
        progress("完成", 1, 1, "")

    return {
        "ok": True,
        "out_dir": str(out),
        "stage": STAGE,
        "n_designs": len(records),
        "summary": summ,
        "designs": records,
        "files": {
            "report_html": "/files/report.html",
            "report_json": "/files/report.json",
            "offtarget_tsv": "/files/offtarget_genome.tsv",
            "bundle": f"/files/{pipeline.BUNDLE_NAME}",
        },
    }
