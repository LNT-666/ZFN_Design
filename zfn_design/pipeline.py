"""把 7 个阶段串成一条流水线。

``cli.py`` 只负责解析参数，真正的编排都在这里，这样测试可以直接调用
:func:`run`，不必去 fork 一个进程。
"""

from __future__ import annotations

import uuid
from pathlib import Path
from typing import Any, Sequence

from . import (__version__, designability, fusion, genome, genome_landing,
               genome_offtarget,
               offtarget_genome, report, score, targets)
from .build import build_for_site
from .library import FingerLibrary, load_default_library
from .models import OffTargetHit, OffTargetReport, SequenceRecord
from .offtarget import search as offtarget_search
from .offtarget import search_single as offtarget_search_single

#: 组合型（成对 ZFN + FokI 二聚化）的架构名；其余都按单体型处理。
PAIRED_ARCHITECTURES = {"paired_nuclease"}

DISCLAIMERS = [
    "本工具输出的是候选设计，不是活性保证；所有结果必须实验验证",
    "打分不是亲和力预测，只反映'按公开设计规则看是否更像靠谱位点'",
    "脱靶分析是序列层面的启发式搜索，不等同于实际切割位点",
    "三指 ZFP 只读取 9 bp，单独作为单体在人类基因组中并不唯一；若需要基因组单拷贝结合，"
    "应将其作为 ZFN 配对的一半，或扩展为六指阵列并重新验证。本文件不保证活性或特异性，"
    "所有设计均需实验验证。",
]

#: 每个指头读 3 bp。三指阵列 = 9 bp 半位点，六指阵列 = 18 bp 半位点。
BASES_PER_FINGER = 3


def run(
    records: Sequence[SequenceRecord],
    *,
    library: FingerLibrary | None = None,
    spacers: Sequence[int] = (5, 6, 7),
    fingers: int = 3,
    both_strands: bool = True,
    max_tier: int = 2,
    n_designs: int = 5,
    max_per_finger: int = 2,
    max_sites: int | None = 20,
    repeats_path: str | None = None,
    genome_fasta: str | None = None,
    landing_chunk: int = 8_000_000,
    landing_max_hits: int = 2000,
    offtarget_top: int = 0,
    offtarget_max_mm: int = 2,
    offtarget_max_hits: int = 500,
    min_copies_remove: int = 3,
    template_key: str | None = None,
    architecture: str = "paired_nuclease",
    offtarget_records: Sequence[SequenceRecord] | None = None,
    max_mm_left: int = 2,
    max_mm_right: int = 2,
    weights_path: str | Path | None = None,
    fusion_keys: Sequence[str] | None = None,
    extra_parts: dict[str, str] | None = None,
    command: str = "",
    progress=None,
) -> dict[str, Any]:
    lib = library or load_default_library()
    weights_data = score.load_weights(weights_path)
    mode = "paired" if architecture in PAIRED_ARCHITECTURES else "single"

    def _p(phase: str, done: float = 0, total: float = 0, msg: str = "") -> None:
        if progress:
            progress(phase, done, total, msg)

    _p("扫描候选位点", 0, len(records), "逐条序列扫 27 bp 位点")
    sites = targets.scan_many(
        records, spacers=tuple(spacers), half_site_len=BASES_PER_FINGER * fingers,
        both_strands=both_strands, mode=mode,
    )

    # ---- Target filter（第 1 步）：**先过滤后评分**（宋成义 2026-11-18 口径）----
    # 给了 genome_fasta 就按全基因组口径判：
    #   ① 候选 motif 的任一基因组落点压在 RepeatMasker 注释区 → 剔除
    #   ② 同一 motif 的基因组落点数 >= min_copies_remove → 剔除（兜注释漏掉的低复杂度）
    # 只给 repeats_path 时退回旧口径（按输入序列自身坐标叠注释），仅作兼容保留。
    n_sites_raw = len(sites)
    n_repeat_filtered = 0
    landing_summary: dict[str, Any] | None = None
    if genome_fasta:
        total_bp = sum(n for _c, n in genome.fasta_index(genome_fasta))

        def _lp(info: dict) -> None:
            _p("全基因组重复过滤", info.get("scanned_bp", 0), total_bp,
               f"已扫 {info.get('chrom', '')}")

        sites, landing_summary = genome_landing.filter_sites_by_genome(
            sites,
            genome_fasta,
            repeats_path=repeats_path,
            min_copies=min_copies_remove,
            chunk=landing_chunk,
            max_hits=landing_max_hits,
            progress=_lp,
        )
        n_repeat_filtered = (
            landing_summary["n_removed_repeat"] + landing_summary["n_removed_multicopy"]
        )
    elif repeats_path:
        ridx = genome.parse_repeatmasker(repeats_path)
        kept = [s for s in sites if not ridx.overlaps(s.seq_name, s.start - 1, s.end)]
        n_repeat_filtered = n_sites_raw - len(kept)
        sites = kept

    # ---- 查重（第 2 步）：同一段 site_seq 出现 >= min_copies_remove 次 → 整段剔除；
    #      恰好 2 次的保留第一处，第二处不算可设计 ----
    counts: dict[str, int] = {}
    for s in sites:
        counts[s.site_seq] = counts.get(s.site_seq, 0) + 1
    removed = {k for k, v in counts.items() if v >= min_copies_remove}
    seen: set[str] = set()
    unique_sites = []
    n_multi_copy = 0
    for s in sites:
        if s.site_seq in removed:
            continue
        if s.site_seq in seen:
            n_multi_copy += 1
            continue
        seen.add(s.site_seq)
        unique_sites.append(s)
    n_removed_sites = sum(1 for s in sites if s.site_seq in removed)
    sites = unique_sites

    dz = designability.evaluate_many(sites, lib, max_tier=max_tier)
    _p("三联体可设计性评估", 1, 1, f"{len(dz)} 个候选")
    dz_sorted = designability.triage(dz)

    designable = [d for d in dz_sorted if d.status != "gap"]
    if max_sites is not None:
        designable = designable[:max_sites]

    designs = []
    for i, d in enumerate(designable, 1):
        designs.extend(
            build_for_site(d.site, lib, d, n_designs=n_designs,
                           max_per_finger=max_per_finger, max_tier=max_tier,
                           template_key=template_key, architecture=architecture)
        )
        if i == 1 or i == len(designable) or i % 10 == 0:
            _p("生成 ZFN 设计", i, len(designable), f"{len(designs)} 个设计")

    search_space = list(offtarget_records) if offtarget_records else list(records)
    for design in designs:
        tgt = design.target
        # on-target 的两个链是**同一段物理 DNA**，必须一起排除；
        # 只排除本链会让回文位点把自己的另一条链报成 0 错配脱靶。
        twin = "-" if tgt.strand == "+" else "+"
        exclude = {(tgt.seq_name, tgt.start, tgt.strand), (tgt.seq_name, tgt.start, twin)}
        if tgt.is_paired:
            design.offtarget = offtarget_search(
                search_space,
                left_half=tgt.left_half,
                right_half=tgt.right_half,
                spacers=tuple(spacers),
                max_mm_left=max_mm_left,
                max_mm_right=max_mm_right,
                exclude=exclude,
            )
        else:
            design.offtarget = offtarget_search_single(
                search_space,
                half_site=tgt.left_half,
                max_mm=max(max_mm_left, max_mm_right),
                exclude=exclude,
            )

    designs = score.rank(designs, weights_data)
    _p("打分排序", 1, 1, f"{len(designs)} 个设计")

    # ---- S3 · Off-target score：取 Design_score 前 N 个做全基因组脱靶负担 ----
    # 结果并列展示，**不加进** Design_score（宋成义 2026-11-18 口径）。
    offtarget_genome_records: list[dict] = []
    offtarget_genome_summary: dict[str, Any] | None = None
    if offtarget_top and genome_fasta and designs:
        offtarget_genome_records = genome_offtarget.score_designs(
            designs,
            genome_fasta,
            top=offtarget_top,
            max_mm=offtarget_max_mm,
            max_hits=offtarget_max_hits,
        )
        offtarget_genome_summary = genome_offtarget.summary(offtarget_genome_records)

    fusions: list[dict] = []
    if fusion_keys is None:
        fusion_keys = (["single_activator", "single_repressor"] if mode == "single"
                       else ["single_activator", "paired_nuclease"])
    top = designs[: max(1, min(len(designs), 5))]
    if top and fusion_keys:
        fusions = fusion.build_fusions(top, list(fusion_keys), extra_parts=extra_parts)

    payload = {
        "meta": {
            "run_id": uuid.uuid4().hex[:12],
            "generated_at": report.now_iso(),
            "version": __version__,
            "command": command,
            "stage": (
                "S3" if offtarget_genome_summary else ("S2" if genome_fasta else "S1")
            ),
            "stage_title": (
                "Stage 3 · Off-target score（全基因组脱靶评分）"
                if offtarget_genome_summary
                else (
                    "Stage 2 · Target filter（filtered-target 设计报告）"
                    if genome_fasta
                    else "Stage 1 · DesignScore（all-target 设计报告）"
                )
            ),
            "report_kind": (
                "offtarget_score"
                if offtarget_genome_summary
                else ("filtered_target" if genome_fasta else "all_target")
            ),
            "design_score_name": score.DESIGN_SCORE_NAME,
            "params": {
                "spacers": list(spacers),
                "fingers": fingers,
                "half_site_bp": BASES_PER_FINGER * fingers,
                "both_strands": both_strands,
                "max_tier": max_tier,
                "n_designs": n_designs,
                "max_per_finger": max_per_finger,
                "max_sites": max_sites,
                "template": template_key,
                "architecture": architecture,
                "mode": mode,
                "max_mm_left": max_mm_left,
                "max_mm_right": max_mm_right,
                "offtarget_space_bp": sum(r.length for r in search_space),
                "genome_fasta": genome_fasta,
                "repeats_path": repeats_path,
                "min_copies_remove": min_copies_remove,
                "offtarget_top": offtarget_top,
                "offtarget_max_mm": offtarget_max_mm,
            },
        },
        "inputs": [
            {"name": r.name, "length": r.length, "gc": round(
                sum(1 for c in r.seq if c in "GC") / max(1, r.length), 4),
             "n_bases": r.seq.count("N"), "source": r.source}
            for r in records
        ],
        "library": lib.coverage_stats(),
        "targets_summary": {
            **designability.summary(dz),
            "n_sites_raw": n_sites_raw,
            "n_repeat_filtered": n_repeat_filtered,
            "n_unique_sites": len(seen),
            "n_multi_copy": n_multi_copy,
            "n_removed_ge3": len(removed),
            "n_removed_ge3_sites": n_removed_sites,
            "min_copies_remove": min_copies_remove,
            "genome_filter": (
                {
                    **{k: v for k, v in landing_summary.items() if k != "per_site"},
                    "n_filtered_sites": n_repeat_filtered,
                    "removed_sites": [
                        {k: v for k, v in rec.items() if k != "landings"}
                        for rec in landing_summary.get("per_site", [])
                        if rec["verdict"] != "pass"
                    ][:200],
                }
                if landing_summary
                else None
            ),
        },
        "sites": [d.to_dict() for d in dz_sorted],
        "designs": [dict(d.to_dict(), rank=i + 1) for i, d in enumerate(designs)],
        "fusions": fusions,
        "offtarget_genome": (
            {"summary": offtarget_genome_summary, "designs": offtarget_genome_records}
            if offtarget_genome_summary
            else None
        ),
        "scoring": score.weights_documentation(weights_data),
        "disclaimers": list(DISCLAIMERS),
        "warnings": _warnings(lib, dz_sorted, designs, fingers=fingers, mode=mode),
    }
    return payload


# --------------------------------------------------------------------------
# 全基因组模式（流式）
# --------------------------------------------------------------------------

GENOME_SITE_COLUMNS = [
    "site_id", "seq_name", "start", "end", "strand", "kind", "site_seq",
    "left_half", "spacer", "spacer_len", "right_half",
    "designability", "n_combos", "missing_triplets",
]


def run_genome(
    fasta_path: str,
    outdir: str | Path,
    *,
    library: FingerLibrary | None = None,
    spacers: Sequence[int] = (5, 6, 7),
    fingers: int = 3,
    both_strands: bool = True,
    max_tier: int = 2,
    repeats_path: str | None = None,
    regions_path: str | None = None,
    limit_sites: int = 200_000,
    n_designs: int = 0,
    max_design_sites: int = 20,
    max_per_finger: int = 2,
    dedupe: bool = True,
    min_copies_remove: int = 3,
    ot_designs: int = 10,
    ot_max_mm: int = 2,
    ot_max_hits: int = 200,
    ot_chunk: int = 4_000_000,
    template_key: str | None = None,
    architecture: str = "paired_nuclease",
    weights_path: str | None = None,
    chunk: int = genome.DEFAULT_CHUNK,
    progress=None,
    command: str = "",
) -> dict:
    """在全基因组（或 BED 指定的区域）里流式扫靶点，边扫边写 TSV。

    与 :func:`run` 的区别：
    - 不物化 ``list[TargetSite]``，候选数可以到十亿级而内存不涨；
    - 重复区（RepeatMasker）与目标区间（BED）在构造 ``TargetSite`` 之前就过滤；
    - 同一段 ``site_seq`` 出现 ``>= min_copies_remove`` 次时，**整段剔除**
      （含第一处）；只出现 2 次的保留第一处，第二处标 ``multi_copy`` 不计入可设计统计；
    - 设计分（Design score）只由 6 个设计分项加权得到；全基因组脱靶是**另一步**，
      对 Design score 最高的前 ``ot_designs`` 个设计单独算 Off-target score，不再回填进设计分。
    """
    import time as _time

    t0 = _time.perf_counter()
    lib = library or load_default_library()
    mode = "paired" if architecture in PAIRED_ARCHITECTURES else "single"
    hsl = BASES_PER_FINGER * fingers
    spacers = tuple(sorted({int(s) for s in spacers}))

    repeats = genome.parse_repeatmasker(repeats_path) if repeats_path else None
    regions = genome.parse_bed(regions_path) if regions_path else None

    chroms = genome.fasta_index(fasta_path)
    wdata = score.load_weights(weights_path)
    total_bp = sum(n for _, n in chroms)
    repeat_bp = repeats.summary()["total_bp"] if repeats else 0
    estimate = genome.estimate_sites(
        chroms, spacers=spacers, half_site_len=hsl, both_strands=both_strands, mode=mode,
        repeat_fraction=(repeat_bp / total_bp) if (repeats and total_bp) else 0.0,
    )

    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    sites_path = out / "sites.tsv"

    by_status: dict[str, int] = {"full": 0, "partial": 0, "gap": 0}
    by_chrom: dict[str, int] = {}
    preview: list[dict] = []
    candidates: list[tuple] = []
    seq_counts: dict[str, int] = {}
    n_sites_raw = 0
    n_multi_copy = 0
    scan_stats: dict = {}

    fh = open(sites_path, "w", encoding="utf-8-sig", newline="")
    buf: list[str] = []

    def _flush() -> None:
        if buf:
            fh.write("".join(buf))
            buf.clear()

    try:
        fh.write("\t".join(GENOME_SITE_COLUMNS) + "\n")
        for site in genome.scan_genome(
            fasta_path, spacers=spacers, half_site_len=hsl, both_strands=both_strands,
            mode=mode, repeats=repeats, regions=regions, limit=limit_sites,
            chunk=chunk, progress=progress, stats=scan_stats,
        ):
            n_sites_raw += 1
            copies = seq_counts.get(site.site_seq, 0)
            seq_counts[site.site_seq] = copies + 1

            # 同一段序列第 2 次出现 = 不是单靶点：标 multi_copy，不参与设计、不计入可设计统计。
            # （>= min_copies_remove 次的，扫完后再整段剔除，见下方二次处理。）
            if dedupe and copies:
                buf.append("\t".join(str(x) for x in [
                    site.site_id, site.seq_name, site.start, site.end, site.strand,
                    site.kind, site.site_seq, site.left_half, site.spacer,
                    len(site.spacer), site.right_half, "multi_copy", 0, "",
                ]) + "\n")
                if len(buf) >= 2000:
                    _flush()
                continue

            dz = designability.evaluate(site, lib, max_tier=max_tier)
            buf.append("\t".join(str(x) for x in [
                site.site_id, site.seq_name, site.start, site.end, site.strand,
                site.kind, site.site_seq, site.left_half, site.spacer,
                len(site.spacer), site.right_half, dz.status, dz.n_combos,
                ",".join(dz.missing_triplets),
            ]) + "\n")
            if len(buf) >= 2000:
                _flush()

            if len(preview) < 200:
                preview.append(dz.to_dict())

            if n_designs > 0 and dz.status != "gap" and len(candidates) < max_design_sites:
                candidates.append((site, dz))
    finally:
        _flush()
        fh.close()

    # ---- 扫完之后才知道谁被复制了：site_seq 出现 >= min_copies_remove 次 → 整段剔除 ----
    removed_seqs = ({s for s, c in seq_counts.items() if c >= min_copies_remove}
                    if dedupe else set())
    n_removed_sites = 0
    text = sites_path.read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    header, body = lines[0], lines[1:]
    kept_rows: list[str] = []
    by_status = {"full": 0, "partial": 0, "gap": 0}
    by_chrom = {}
    n_multi_copy = 0
    for line in body:
        cols = line.split("\t")
        if len(cols) > 11 and cols[6] in removed_seqs:
            n_removed_sites += 1
            continue
        st = cols[11] if len(cols) > 11 else ""
        if st == "multi_copy":
            n_multi_copy += 1
            kept_rows.append(line)
            continue
        if st in by_status:
            by_status[st] += 1
        if len(cols) > 1:
            by_chrom[cols[1]] = by_chrom.get(cols[1], 0) + 1
        kept_rows.append(line)
    if removed_seqs:
        sites_path.write_text("\n".join([header] + kept_rows) + "\n", encoding="utf-8-sig")
    n_sites = sum(by_status.values())

    # ---- 设计：只对扫完仍存活的候选出设计 ----
    designs = []
    for site, dz in candidates:
        if site.site_seq in removed_seqs:
            continue
        try:
            built = build_for_site(
                site, lib, dz, n_designs=n_designs, max_per_finger=max_per_finger,
                max_tier=max_tier, template_key=template_key,
                architecture=architecture,
            )
            for d in built:
                d.score = score.score_design(d, weights_data=wdata)
                designs.append(d)
        except Exception as exc:  # 单个位点失败不该中断整条基因组扫描
            if preview:
                preview[0].setdefault("build_error", str(exc))

    designs.sort(key=lambda d: -((d.score.total if d.score else 0.0)))

    # ---- 单独一步：全基因组脱靶负担（只对 Design score 前 N 个设计） ----
    ot_rows: list[dict[str, Any]] = []
    genome_ot: dict[str, Any] = {}
    if ot_designs > 0 and designs:
        for d in designs[:ot_designs]:
            tgt = d.target
            spacers_ot = (len(tgt.spacer),) if tgt.spacer else (0,)
            res = offtarget_genome.genome_search(
                fasta_path, tgt.left_half, tgt.right_half,
                spacers=spacers_ot,
                max_mm_left=ot_max_mm, max_mm_right=ot_max_mm,
                max_hits=ot_max_hits, chunk=ot_chunk, progress=progress,
            )
            res["design_id"] = d.design_id
            res["score_100"] = round(offtarget_genome.score_from_burden(res["burden"]), 2)
            genome_ot[d.design_id] = res
            for h in res["hits"]:
                ot_rows.append({"design_id": d.design_id, **h})
            # 明细留在 offtarget 里，但**不回填**进设计分：脱靶单独给分。
            d.offtarget = OffTargetReport(
                hits=[
                    OffTargetHit(
                        seq_name=h["chrom"], start=h["start0"] + 1, strand=h["strand"],
                        left_half=tgt.left_half, spacer="N" * h["spacer"],
                        right_half=tgt.right_half,
                        mm_left=h["mm_left"], mm_right=h["mm_right"],
                    )
                    for h in res["hits"]
                ],
                searched_bp=res["scanned_bp"],
                max_mm_left=ot_max_mm, max_mm_right=ot_max_mm,
                spacers=tuple(spacers_ot), truncated=res["truncated"],
            )
            d.offtarget_score = {
                "score": res["score_100"],
                "burden": round(float(res["burden"]), 4),
                "n_hits": len(res["hits"]),
                "scanned_bp": res["scanned_bp"],
                "max_mm": ot_max_mm,
                "truncated": bool(res.get("truncated")),
                "engine": res.get("engine", "?"),
                "source": "genome_wide",
            }

    elapsed = _time.perf_counter() - t0
    n_designs_out = len(designs)
    scan_stats["removed_ge3_seqs"] = len(removed_seqs)
    scan_stats["removed_ge3_sites"] = n_removed_sites

    payload: dict[str, Any] = {
        "meta": {
            "tool": "zfn_design",
            "version": __version__,
            "run_id": uuid.uuid4().hex[:8],
            "generated_at": report.now_iso(),
            "command": command,
            "mode": mode,
            "genome_mode": True,
            "elapsed_s": round(elapsed, 3),
            "params": {
                "spacers": list(spacers), "fingers": fingers,
                "half_site_bp": hsl, "architecture": architecture,
                "max_tier": max_tier, "limit_sites": limit_sites,
                "chunk": chunk, "both_strands": both_strands,
                "min_copies_remove": min_copies_remove,
                "repeats": Path(repeats_path).name if repeats_path else None,
                "ot_designs": ot_designs,
            },
        },
        "inputs": [{"name": Path(fasta_path).name, "kind": "genome",
                    "n_chroms": len(chroms), "total_bp": total_bp}],
        "library": lib.coverage_stats(),
        "targets_summary": {
            "n_sites": n_sites,
            "n_sites_raw": n_sites_raw,
            "n_repeat_filtered": int(scan_stats.get("skipped_repeat", 0)),
            "n_unique_sites": len(seq_counts),
            "n_multi_copy": n_multi_copy,
            "n_removed_ge3": len(removed_seqs),
            "n_removed_ge3_sites": n_removed_sites,
            "min_copies_remove": min_copies_remove,
            "by_status": by_status,
            "by_chrom": by_chrom,
            "estimated_total": estimate,
            "scan_stats": scan_stats,
            "repeat_bp": repeat_bp,
            "repeat_intervals": repeats.summary() if repeats else None,
            "regions": regions.summary() if regions else None,
            "truncated": n_sites_raw >= limit_sites,
            "sites_tsv": str(sites_path),
        },
        "sites": preview,
        "designs": [d.to_dict() for d in designs],
        "genome_offtarget": genome_ot,
        "offtarget_score_info": {
            "source": "genome_wide",
            "top_n": ot_designs,
            "max_mm": ot_max_mm,
            "scale": wdata["rules"]["offtarget"]["burden_scale"],
            "note": "Off-target score 只对 Design score 前 N 个设计算，不并入 Design score",
        },
        "fusions": [],
        "scoring": score.weights_documentation(wdata),
        "disclaimers": list(DISCLAIMERS),
        "warnings": _genome_warnings(
            n_sites=n_sites, estimate=estimate, limit_sites=limit_sites,
            repeats=repeats, regions=regions, n_designs=n_designs_out,
            mode=mode, fingers=fingers, elapsed=elapsed,
            n_multi_copy=n_multi_copy, scan_stats=scan_stats,
            ot_note=(_ot_note(genome_ot, ot_max_mm, ot_max_hits) if genome_ot else None),
        ),
    }

    report.write_json(out / "report.json", payload)
    if ot_rows:
        _write_offtargets_tsv(out / "offtargets.tsv", ot_rows)
    if payload["designs"]:
        report.designs_tsv(out / "designs.tsv", payload["designs"])
        report.write_html(out / "report.html", payload)
    return payload


def _ot_note(genome_ot: dict[str, Any], ot_max_mm: int, ot_max_hits: int) -> str:
    n = len(genome_ot)
    worst = max(genome_ot.values(), key=lambda r: r["burden"])
    w0 = worst["by_mm"].get(0, 0)
    engine = worst.get("engine", "?")
    extra = "（numpy 加速）" if engine == "numpy" else "（纯 Python，慢）"
    msg = (
        f"全基因组脱靶：对分数最高的 {n} 个设计各扫了一遍基因组，"
        f"错配上限 {ot_max_mm}/侧。burden 最高的是 {worst['design_id']}"
        f"（{worst['n_total_hits']} 个近似位点，其中 0 错配 {w0} 个，"
        f"burden {worst['burden']}）。明细见 offtargets.tsv{extra}。"
    )
    if any(r["truncated"] for r in genome_ot.values()):
        msg += f" 有设计命中数超过 --ot-max-hits {ot_max_hits}，清单已截断，burden 仍是全量统计。"
    return msg


def _write_offtargets_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    cols = ["design_id", "chrom", "start0", "strand", "spacer", "mm_left", "mm_right",
            "mm_total", "weight"]
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def _genome_warnings(*, n_sites: int, estimate: int, limit_sites: int,
                     repeats, regions, n_designs: int, mode: str,
                     fingers: int, elapsed: float,
                     n_multi_copy: int = 0, scan_stats: dict | None = None,
                     ot_note: str | None = None) -> list[str]:
    out: list[str] = []
    out.append(
        f"基因组模式：扫到 {n_sites:,} 个候选窗口（预估上限 {estimate:,}），"
        f"用时 {elapsed:.1f} s。"
    )
    if scan_stats:
        sk = scan_stats.get("skipped_repeat", 0)
        win = scan_stats.get("windows", 0) or 1
        if sk:
            out.append(
                f"重复序列过滤剔除了 {sk:,} 个窗口，占扫过的 {sk / win:.1%}"
                "（RepeatMasker 注释覆盖的区间一律不产出候选）。"
            )
    if n_multi_copy:
        out.append(
            f"有 {n_multi_copy:,} 个候选的序列在基因组里**出现了不止一次**（已在 sites.tsv 标为 "
            "multi_copy，且不参与设计）。重复拷贝意味着不能当作单靶点——这正是要避开重复区的原因。"
        )
    if ot_note:
        out.append(ot_note)
    if n_sites >= limit_sites:
        out.append(
            f"已达 --limit {limit_sites:,} 上限而提前停止：结果**不是全量**，"
            "要全量请调大 --limit 或用 --regions 先缩小范围。"
        )
    if repeats is not None:
        st = repeats.summary()
        out.append(
            f"已按 RepeatMasker 注释剔除落在重复序列里的候选（{st['n_intervals_merged']:,} 个区间，"
            f"合计 {st['total_bp']:,} bp）。"
        )
    else:
        out.append(
            "**未提供重复序列注释**：结果里仍包含重复区（如 Alu/LINE）的候选，"
            "这些位置的特异性很差，正式使用前请加 --repeats。"
        )
    if regions is not None:
        out.append(f"只扫了 BED 指定的区域（{regions.summary()['n_intervals_merged']:,} 段）。")
    if n_designs and not ot_note:
        out.append("基因组模式下的设计**未做脱靶搜索**（需在基因组尺度再扫一遍），"
                   "脱靶请加 --ot-designs 或使用 --offtarget-fasta / 外部工具单独跑。")
    if mode == "single":
        out.append(
            f"当前是单体型（{fingers} 指，读 {BASES_PER_FINGER * fingers} bp）："
            "只扫一个半位点，不产生双链断裂，命中数天然多于组合型。"
        )
    return out


def _warnings(lib: FingerLibrary, dz, designs, *, fingers: int = 3,
              mode: str = "paired") -> list[str]:
    out: list[str] = []
    stats = lib.coverage_stats()
    if stats["triplets_covered"] < 16:
        out.append(
            f"指头库只覆盖 {stats['triplets_covered']}/64 个三联体，绝大多数位点会被判为缺口。"
            "这是随包种子库的预期表现——真实项目请导入 B1H 数据集（见文档）。"
        )
    n_fb = sum(1 for d in dz if d.context_fallback)
    if n_fb:
        out.append(
            f"{n_fb} 个位点用到了「非本位上下文」的指头：库里的实测数据只覆盖部分阵列位置"
            "（Princeton B1H 提供 F2/F3），缺口位置退而使用同一三联体在别的指位实测的指头。"
            "这些指头在该位点的行为未经验证。"
        )
    if not designs:
        out.append("没有生成任何设计——通常是位点全部落在缺口上。")
    t2 = sum(1 for d in designs if d.min_tier > 1)
    if t2:
        out.append(f"{t2} 个设计含 Tier2（未验证）指头，必须实验验证")
    if mode == "single":
        out.append(
            f"当前是**单体型**（{fingers} 指阵列，读 {BASES_PER_FINGER * fingers} bp）："
            "只有一个半位点、没有 spacer、没有 FokI 二聚化，因此不产生双链断裂；"
            "脱靶判据是「一个半位点像不像」而不是「一对半位点像不像」，"
            "命中数会明显多于组合型，两者的分数不可直接比较。"
        )
    elif fingers <= 3:
        out.append(
            "三指阵列只读 9 bp：单个三指 ZFP 在人类基因组中并不唯一。"
            "需要单拷贝级别的结合特异性时，应把它当作 ZFN 配对的一半（本工具默认成对输出），"
            "或改成六指阵列（18 bp）并重新验证。"
        )
    else:
        out.append(
            f"{fingers} 指阵列读 {BASES_PER_FINGER * fingers} bp，理论上比三指更可能唯一，"
            "但六指阵列由两个三指模块串联而成：模块各自的骨架有 1AAY 依据，"
            "模块之间的连接与界面没有任何结构或活性证据，必须单独验证。"
        )
    return out


def _write_offtarget_genome_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    """S3 的逐设计全基因组脱靶负担（一行一个设计）。"""
    cols = ["burden_rank", "design_id", "rank", "seq_name", "start", "strand", "spacer",
            "design_score", "offtarget_score", "burden", "mm0", "mm1", "mm2",
            "n_hits_kept", "n_hits_total", "truncated", "engine", "scanned_bp", "elapsed_s"]
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            by_mm = r.get("by_mm") or {}
            rec = dict(r)
            rec.update({"mm0": by_mm.get("0", 0), "mm1": by_mm.get("1", 0), "mm2": by_mm.get("2", 0)})
            fh.write("\t".join(str(rec.get(c, "")) for c in cols) + "\n")
    return path


BUNDLE_NAME = "zfn_bundle.zip"


def write_outputs(payload: dict, outdir: str | Path, *, bundle: bool = True) -> dict[str, Path]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {
        "json": report.write_json(out / "report.json", payload),
        "html": report.write_html(out / "report.html", payload),
        "sites_tsv": report.sites_tsv(out / "sites.tsv", payload["sites"]),
        "designs_tsv": report.designs_tsv(out / "designs.tsv", payload["designs"]),
    }
    if payload.get("fusions"):
        paths["fusions_json"] = report.write_json(out / "fusions.json", payload["fusions"])
    if payload.get("offtarget_genome"):
        paths["offtarget_tsv"] = _write_offtarget_genome_tsv(
            out / "offtarget_genome.tsv", payload["offtarget_genome"]["designs"]
        )
    if bundle:
        paths["bundle"] = report.write_bundle(out / BUNDLE_NAME, {k: v for k, v in paths.items()})
    return paths


__all__ = ["run", "write_outputs", "DISCLAIMERS"]
