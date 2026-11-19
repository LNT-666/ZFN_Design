"""阶段 S4：Fusion design——从已有作业里挑设计，拼蛋白、出 CDS、出订单表。

与 S1–S3 分开按键：S4 不改动上一次的评分结果，只读 ``report.json``，
把选中的设计接上效应结构，再走 :mod:`codonopt` 生成 CDS 与合成订单。

输入是**已有作业目录**（GUI 的 out 目录），不是序列——因为"挑哪个设计"
这件事发生在设计完成之后，宋成义 2026-11-18 明确要求这一步单独按键。
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from . import codonopt, fusion, report

STAGE = "S4"
STAGE_TITLE = "Stage 4 · Fusion design（蛋白 + 人源化 CDS 订单）"
REPORT_KIND = "fusion_design"

DEFAULT_ARCHS = ("single_activator", "single_repressor", "paired_nuclease")


def _shim(design: dict[str, Any]):
    """把 report.json 里的设计字典包成 build_fusion 需要的最小对象。"""
    target = design.get("target") or {}
    left = design.get("monomer_left") or {}
    return SimpleNamespace(
        design_id=design.get("design_id", ""),
        protein=left.get("protein") or "",
        target=SimpleNamespace(
            is_paired=bool(design.get("is_paired", True)),
            spacer_len=target.get("spacer_len") or target.get("spacer") and len(str(target.get("spacer"))) or 6,
        ),
    )


def _chain_protein(design: dict[str, Any], side: str) -> str:
    """取某条链的 ZF 蛋白（左侧/右侧各 3 指）。"""
    key = "monomer_left" if side == "left" else "monomer_right"
    mono = design.get(key) or {}
    return mono.get("protein") or ""


def select_designs(designs: list[dict[str, Any]], *, design_ids: list[str] | None = None,
                   top: int = 3) -> list[dict[str, Any]]:
    """按 rank 取前 N 个；给了 id 就按 id 取（顺序按调用方给的）。"""
    if design_ids:
        by_id = {d.get("design_id"): d for d in designs}
        missing = [i for i in design_ids if i not in by_id]
        if missing:
            raise ValueError(f"这些设计 id 不在这次作业里：{', '.join(missing)}")
        return [by_id[i] for i in design_ids]
    return list(designs)[: max(0, int(top))]


def _chain_specs(arch_key: str) -> list[dict[str, Any]]:
    """架构要出几条链：成对核酸酶左右各一条（FokI ELD / KK）。"""
    arch = fusion.architectures().get(arch_key) or {}
    chains = arch.get("chains")
    if chains:
        return [dict(c, fok_part=c.get("fok_part")) for c in chains]
    return [{"side": "single", "label": arch.get("label") or arch_key,
             "fok_part": None}]


def _one_chain(design: dict[str, Any], arch_key: str, spec: dict[str, Any], *,
               parts: dict[str, str] | None, seed: int) -> dict[str, Any]:
    side = spec.get("side") or "single"
    zf = _chain_protein(design, side) if side in ("left", "right") else \
        (_chain_protein(design, "left") or _shim(design).protein)
    overrides: dict[str, str] = {}
    if spec.get("fok_part") and spec["fok_part"] != "FokI_CD_ELD":
        # 右侧链换成 KK 变体；左侧保持架构默认的 ELD
        overrides["FokI_CD_ELD"] = spec["fok_part"]
    shim = _shim(design)
    try:
        built = fusion.build_fusion(shim, arch_key, zf_protein=zf,
                                    extra_parts=parts, part_overrides=overrides)
    except fusion.FusionError as exc:
        return {"side": side, "label": spec.get("label") or side,
                "error": str(exc), "cds_status": "skipped"}
    if "error" in built:
        built.update({"side": side, "label": spec.get("label") or side,
                      "cds_status": "skipped"})
        return built

    protein = built.get("draft_protein") or ""
    out: dict[str, Any] = {
        "side": side,
        "label": spec.get("label") or side,
        "fok_part": spec.get("fok_part"),
        "protein": protein,
        "protein_len": len(protein),
        "has_nls": bool(built.get("has_nls")),
        "complete": bool(built.get("is_complete")),
        "red_flags": built.get("red_flags") or [],
        "parts": built.get("parts"),
    }
    if not out["complete"]:
        out["cds_status"] = "blocked"
        out["cds"] = ""
        out["cds_notes"] = ["拼装里有占位符：缺元件序列，CDS 不生成"]
        return out

    cd = codonopt.optimize(protein, seed=seed)
    cds = cd["cds"]
    coding = cd.get("coding") or ""
    starts_atg = coding.startswith("ATG")
    ends_stop = bool(coding) and coding[-3:] in codonopt.STOP_CODONS
    ok = bool(cd["ok"] and starts_atg and ends_stop)
    notes = list(cd["notes"])
    if not starts_atg:
        notes.append("编码区首密码子不是 ATG，不能直接表达")
    if not ends_stop:
        notes.append("编码区末密码子不是终止密码子")
    out.update({
        "cds": cds,
        "cds_coding": coding,
        "cds_len": len(cds),
        "cds_gc": cd["gc"],
        "cds_coding_start": cd.get("cds_coding_start", 0),
        "cds_kozak": cd.get("kozak", ""),
        "cds_start_codon": coding[:3] if coding else "",
        "cds_stop_codon": coding[-3:] if coding else "",
        "cds_has_nls": bool(built.get("has_nls")),
        "cds_homopolymer_max": cd["homopolymer_max"],
        "cds_forbidden_sites": cd["forbidden_sites"],
        "cds_ok": ok,
        "cds_notes": notes,
        "back_translation": cd["back_translation"],
        "cds_status": "ok" if ok else "warn",
    })
    return out


def _entry(design: dict[str, Any], arch_key: str, *, parts: dict[str, str] | None,
           seed: int) -> dict[str, Any]:
    shim = _shim(design)
    chains = [_one_chain(design, arch_key, spec, parts=parts, seed=seed)
              for spec in _chain_specs(arch_key)]
    first = chains[0] if chains else {}
    out: dict[str, Any] = {
        "design_id": shim.design_id,
        "rank": design.get("rank"),
        "architecture": arch_key,
        "architecture_label": (fusion.architectures().get(arch_key) or {}).get("label") or arch_key,
        "requires_pair": bool((fusion.architectures().get(arch_key) or {}).get("dimer")),
        "n_chains": len(chains),
        "chains": chains,
        "design_score": (design.get("score") or {}).get("total"),
        # S3 跑过之后 design 里有 offtarget_score；带上，订单里一并显示
        "offtarget_score": (design.get("offtarget_score") or {}).get("score"),
        "offtarget_burden": (design.get("offtarget_score") or {}).get("burden"),
    }
    # 顶层保留第一条链的字段，方便旧界面/旧调用方读取
    for k in ("protein", "protein_len", "has_nls", "complete", "red_flags", "parts",
              "cds", "cds_len", "cds_gc", "cds_coding_start", "cds_kozak",
              "cds_start_codon", "cds_stop_codon", "cds_has_nls",
              "cds_homopolymer_max", "cds_forbidden_sites", "cds_ok", "cds_notes",
              "back_translation", "cds_status", "error"):
        if k in first:
            out[k] = first[k]
    if "error" in first:
        return out
    out["cds_status"] = _worst_status([c.get("cds_status") for c in chains])
    return out


_STATUS_ORDER = {"ok": 0, "warn": 1, "blocked": 2, "skipped": 3}


def _worst_status(statuses: list[str | None]) -> str:
    known = [s for s in statuses if s]
    if not known:
        return "skipped"
    return max(known, key=lambda s: _STATUS_ORDER.get(s, 9))



def build_entries(designs: list[dict[str, Any]], *, design_ids: list[str] | None = None,
                  top: int = 3, arch_keys: list[str] | None = None,
                  parts: dict[str, str] | None = None, seed: int = 0) -> list[dict[str, Any]]:
    chosen = select_designs(designs, design_ids=design_ids, top=top)
    keys = list(arch_keys or DEFAULT_ARCHS)
    known = set(fusion.architectures())
    unknown = [k for k in keys if k not in known]
    if unknown:
        raise ValueError(f"未知架构：{', '.join(unknown)}；可用 {sorted(known)}")
    out: list[dict[str, Any]] = []
    for d in chosen:
        for k in keys:
            e = _entry(d, k, parts=parts, seed=seed)
            e["fusion_id"] = f"{e.get('design_id')}__{k}"
            out.append(e)
    return out


def _chain_rows(entries: list[dict[str, Any]]):
    """把每条融合拆成「一行一条链」，单链架构就是一行。"""
    for e in entries:
        chains = e.get("chains") or []
        if not chains:
            yield e, e
            continue
        for ch in chains:
            yield e, ch


def summary(entries: list[dict[str, Any]]) -> dict[str, Any]:
    rows = list(_chain_rows(entries))
    ok = [ch for _e, ch in rows if ch.get("cds_status") == "ok"]
    blocked = [ch for _e, ch in rows if ch.get("cds_status") == "blocked"]
    skipped = [ch for _e, ch in rows if ch.get("cds_status") == "skipped"]
    gcs = [ch["cds_gc"] for ch in ok if ch.get("cds_gc") is not None]
    return {
        "n_entries": len(entries),
        "n_chains": len(rows),
        "n_ok": len(ok),
        "n_chains_ok": len(ok),
        "n_warn": sum(1 for _e, ch in rows if ch.get("cds_status") == "warn"),
        "n_blocked": len(blocked),
        "n_skipped": len(skipped),
        "cds_total_bp": sum(int(ch.get("cds_len") or 0) for ch in ok),
        "gc_min": round(min(gcs), 4) if gcs else None,
        "gc_mean": round(sum(gcs) / len(gcs), 4) if gcs else None,
        "gc_max": round(max(gcs), 4) if gcs else None,
        "architectures": sorted({e.get("architecture") for e in entries if e.get("architecture")}),
    }


def write_outputs(out_dir: Path, meta: dict[str, Any], entries: list[dict[str, Any]],
                  summ: dict[str, Any]) -> dict[str, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    fj = out_dir / "fusions.json"
    fj.write_text(json.dumps({"meta": meta, "summary": summ, "entries": entries},
                             ensure_ascii=False, indent=2), encoding="utf-8")
    paths["fusions_json"] = fj

    pt = out_dir / "proteins.tsv"
    cols = ["fusion_id", "chain", "chain_label", "fok_part", "design_id", "rank",
            "architecture", "architecture_label", "requires_pair", "protein_len",
            "has_nls", "complete", "cds_status", "red_flags", "protein"]
    with pt.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for e, ch in _chain_rows(entries):
            row = {
                "fusion_id": e.get("fusion_id"), "chain": ch.get("side"),
                "chain_label": ch.get("label"), "fok_part": ch.get("fok_part"),
                "design_id": e.get("design_id"), "rank": e.get("rank"),
                "architecture": e.get("architecture"),
                "architecture_label": e.get("architecture_label"),
                "requires_pair": e.get("requires_pair"),
                "protein_len": ch.get("protein_len"), "has_nls": ch.get("has_nls"),
                "complete": ch.get("complete"), "cds_status": ch.get("cds_status"),
                "red_flags": " ; ".join(ch.get("red_flags") or []),
                "protein": ch.get("protein"),
            }
            w.writerow([row.get(c, "") for c in cols])
    paths["proteins_tsv"] = pt

    ct = out_dir / "cds.tsv"
    ccols = ["fusion_id", "chain", "chain_label", "design_id", "architecture",
             "cds_len", "cds_gc", "cds_kozak", "cds_start_codon", "cds_stop_codon",
             "cds_has_nls", "cds_homopolymer_max", "cds_sites", "cds_ok", "cds_status",
             "cds_notes"]
    with ct.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(ccols)
        for e, ch in _chain_rows(entries):
            sites = ch.get("cds_forbidden_sites") or []
            w.writerow([e.get("fusion_id"), ch.get("side"), ch.get("label"),
                        e.get("design_id"), e.get("architecture"),
                        ch.get("cds_len") or 0, ch.get("cds_gc") or "",
                        ch.get("cds_kozak") or "",
                        ch.get("cds_start_codon") or "", ch.get("cds_stop_codon") or "",
                        ch.get("cds_has_nls", ""),
                        ch.get("cds_homopolymer_max") or "",
                        " ; ".join(f"{s['enzyme']}@{s['pos']}" for s in sites),
                        ch.get("cds_ok", ""), ch.get("cds_status"),
                        " ; ".join(ch.get("cds_notes") or [])])
    paths["cds_tsv"] = ct

    fa = out_dir / "cds.fasta"
    with fa.open("w", encoding="utf-8", newline="\n") as fh:
        for e, ch in _chain_rows(entries):
            if ch.get("cds_status") not in ("ok", "warn") or not ch.get("cds"):
                continue
            head = (f">{e['fusion_id']} chain={ch.get('side')} design={e['design_id']} "
                    f"arch={e['architecture']} aa={ch.get('protein_len')} "
                    f"bp={ch.get('cds_len')} gc={ch.get('cds_gc')} "
                    f"kozak={ch.get('cds_kozak') or 'none'} status={ch.get('cds_status')}")
            fh.write(head + "\n")
            seq = ch["cds"]
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
    paths["cds_fasta"] = fa

    of = out_dir / "order_form.csv"
    with of.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["订单号", "融合编号", "链", "设计", "架构", "CDS(bp)", "GC",
                    "Kozak", "起始", "终止", "同聚物上限", "酶切位点", "状态", "备注"])
        for i, (e, ch) in enumerate(_chain_rows(entries), start=1):
            sites = ch.get("cds_forbidden_sites") or []
            w.writerow([f"S4-{i:03d}", e.get("fusion_id"), ch.get("label") or ch.get("side"),
                        e.get("design_id"),
                        e.get("architecture_label") or e.get("architecture"),
                        ch.get("cds_len") or 0,
                        f"{ch['cds_gc']:.1%}" if ch.get("cds_gc") is not None else "",
                        ch.get("cds_kozak") or "无",
                        ch.get("cds_start_codon") or "", ch.get("cds_stop_codon") or "",
                        ch.get("cds_homopolymer_max") or "",
                        " ; ".join(s["enzyme"] for s in sites) or "clean",
                        ch.get("cds_status"),
                        " ; ".join((ch.get("red_flags") or []) + (ch.get("cds_notes") or []))])
    paths["order_form"] = of
    return paths


def rebuild_bundle(out_dir: Path) -> Path:
    """把整包 zip 重打一遍，带上 S4 的四个文件。"""
    out_dir = Path(out_dir)
    files: dict[str, Path] = {}
    for name in ("report.html", "report.json", "sites.tsv", "designs.tsv",
                 "offtarget_genome.tsv", "fusions.json", "proteins.tsv", "cds.tsv",
                 "cds.fasta", "order_form.csv"):
        p = out_dir / name
        if p.exists():
            files[name] = p
    return report.write_bundle(out_dir / "zfn_bundle.zip", files)


def run_fusion(out_dir: str | Path, *, design_ids: list[str] | None = None, top: int = 3,
               arch_keys: list[str] | None = None, parts: dict[str, str] | None = None,
               seed: int = 0) -> dict[str, Any]:
    """S4 入口：读 report.json → 拼装 → 优化 → 落盘。"""
    out_dir = Path(out_dir)
    rj = out_dir / "report.json"
    if not rj.exists():
        raise ValueError(f"这个作业目录里没有 report.json：{out_dir}")
    payload = json.loads(rj.read_text(encoding="utf-8"))
    designs = payload.get("designs") or []
    if not designs:
        raise ValueError("这次作业没有设计（可能候选被你过滤光了），S4 无从挑选")

    entries = build_entries(designs, design_ids=design_ids, top=top,
                            arch_keys=arch_keys, parts=parts, seed=seed)
    summ = summary(entries)
    meta = {
        "stage": STAGE,
        "stage_title": STAGE_TITLE,
        "report_kind": REPORT_KIND,
        "source_run_id": (payload.get("meta") or {}).get("run_id"),
        "source_stage": (payload.get("meta") or {}).get("stage"),
        "design_ids": [e["design_id"] for e in entries[: len(entries)] if e.get("design_id")],
        "architectures": [k for k in (arch_keys or DEFAULT_ARCHS)],
        "top": top,
        "seed": seed,
    }
    paths = write_outputs(out_dir, meta, entries, summ)
    paths["bundle"] = rebuild_bundle(out_dir)
    return {"ok": True, "out_dir": str(out_dir), "meta": meta, "summary": summ,
            "entries": entries, "files": {k: f"/files/{v.name}" for k, v in paths.items()},
            "paths": {k: str(v) for k, v in paths.items()}}


def order_ready(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [e for e in entries if e.get("cds_status") in ("ok", "warn") and e.get("cds")]


__all__ = ["STAGE", "STAGE_TITLE", "REPORT_KIND", "DEFAULT_ARCHS", "select_designs",
           "build_entries", "summary", "write_outputs", "rebuild_bundle", "run_fusion",
           "order_ready"]
