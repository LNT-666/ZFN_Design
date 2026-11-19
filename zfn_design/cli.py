"""命令行入口。

::

    zfn-design run    --fasta in.fa --out out/
    zfn-design scan   --seq GCGTGGGCG... --out out/
    zfn-design library info
    zfn-design library validate --path mylib.csv
    zfn-design library import-b1h --dir b1h_dump --out mylib.csv
    zfn-design info

不写 ``if __name__`` 之外的花活：参数解析 -> 调 pipeline -> 打印一句人话总结。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import (__version__, designability, framework, fusion, genome, ingest,
               pipeline, report, score, targets)
from .library import FingerLibrary, LibraryError, load_default_library


def _add_input_args(p: argparse.ArgumentParser, *, genome_ok: bool = False) -> None:
    g = p.add_mutually_exclusive_group(required=not genome_ok)
    g.add_argument("--seq", help="直接给 DNA 序列")
    g.add_argument("--fasta", help="FASTA 文件（可多记录）")
    g.add_argument("--stdin", action="store_true", help="从标准输入读")
    if genome_ok:
        g.add_argument("--genome", metavar="FASTA",
                       help="全基因组 FASTA：流式扫描，不进内存")


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("--library", help="指头库 CSV；默认用随包种子库")
    p.add_argument("--spacers", default="5,6,7",
                   help="spacer 长度：逗号分隔或区间，如 6 / 5,6 / 5-7（默认 5,6,7）")
    p.add_argument("--fingers", type=int, default=3,
                   help="每个单体的指头数（2–12）。3 → 半位点 9 bp（1AAY 真实骨架）；"
                        "6 → 18 bp（两个三指模块串联）；其它指数走通用拼接骨架"
                        "（构造体，模块间界面未验证）")
    p.add_argument("--single-strand", action="store_true", help="只扫正链（默认双链）")
    p.add_argument("--max-tier", type=int, default=2, choices=[1, 2],
                   help="允许的最高指头证据 tier（1=只要实测/文献）")
    p.add_argument("--n-designs", type=int, default=5, help="每个位点输出几个设计")
    p.add_argument("--max-per-finger", type=int, default=2, help="每个指位最多取几个候选指头")
    p.add_argument("--max-sites", type=int, default=20, help="最多对几个位点做设计（0=不限）")
    p.add_argument("--template", default=None, help="骨架模板 key（默认 zif268_1aay）")
    p.add_argument("--architecture", default="paired_nuclease", metavar="KEY",
                   help="构型，同时决定单体型 / 组合型：paired_nuclease=组合型（左+spacer+右，"
                        "两单体 FokI 二聚化才切割）；single_activator / single_repressor=单体型"
                        "（只扫一个半位点，无 spacer、无搭档、不产生双链断裂）。默认 paired_nuclease")
    p.add_argument("--offtarget-fasta", help="额外的脱靶搜索序列（默认只搜输入序列）")
    p.add_argument("--max-mm-left", type=int, default=2)
    p.add_argument("--max-mm-right", type=int, default=2)
    p.add_argument("--weights", help="自定义 scoring_weights.json")


def _load_library(path: str | None) -> FingerLibrary:
    if path:
        return FingerLibrary.load_csv(path)
    return load_default_library()


def _parse_spacers(s: str) -> tuple[int, ...]:
    """``"6"`` / ``"5,6"`` / ``"5-7"`` / ``"5,6-7"`` 都能吃。"""
    vals: set[int] = set()
    try:
        for chunk in (s or "").replace(" ", "").split(","):
            if not chunk:
                continue
            if "-" in chunk:
                a, b = chunk.split("-", 1)
                lo, hi = int(a), int(b)
                if lo > hi:
                    raise SystemExit(f"--spacers 区间写反了: {chunk}")
                vals.update(range(lo, hi + 1))
            else:
                vals.add(int(chunk))
    except ValueError as exc:
        raise SystemExit(f"--spacers 写法不认识：{s!r}（试试 6 或 5,6 或 5-7）") from exc
    out = tuple(sorted(vals))
    if not out:
        raise SystemExit("--spacers 为空")
    if out[0] < 0:
        raise SystemExit("--spacers 不能为负数")
    return out


def _check_fingers(n: int) -> None:
    """指头数校验：2–12 指，超出就明确报错（骨架拼不出来）。"""
    lo, hi = framework.MIN_FINGERS, framework.MAX_FINGERS
    if not (lo <= n <= hi):
        raise SystemExit(
            f"--fingers 只能是 {lo}–{hi} 之间的整数，收到 {n}。\n"
            "  3 指 = 9 bp 半位点（1AAY 真实骨架，点这个最稳）\n"
            "  6 指 = 18 bp 半位点（两个三指模块串联）\n"
            "  其它指数走通用拼接骨架：可合成，但模块间界面完全未验证"
        )


def cmd_run(args: argparse.Namespace) -> int:
    _check_fingers(args.fingers)
    known = {a["key"] for a in fusion.describe_architectures()}
    if args.architecture not in known:
        sys.stderr.write(
            f"--architecture 取值不认识：{args.architecture!r}\n"
            f"可选：{', '.join(sorted(known))}\n"
            "  paired_nuclease                 → 组合型（左+spacer+右，FokI 二聚化）\n"
            "  single_activator / single_repressor → 单体型（只扫一个半位点）\n"
        )
        return 2

    records = ingest.read_input(sequence=args.seq, fasta=args.fasta, stdin=args.stdin)
    lib = _load_library(args.library)
    extra = None
    if args.offtarget_fasta:
        extra = ingest.read_input(fasta=args.offtarget_fasta)

    payload = pipeline.run(
        records,
        library=lib,
        spacers=_parse_spacers(args.spacers),
        fingers=args.fingers,
        both_strands=not args.single_strand,
        max_tier=args.max_tier,
        n_designs=args.n_designs,
        max_per_finger=args.max_per_finger,
        max_sites=None if args.max_sites == 0 else args.max_sites,
        template_key=args.template,
        architecture=args.architecture,
        offtarget_records=extra,
        repeats_path=getattr(args, "repeats", None),
        genome_fasta=getattr(args, "genome", None),
        offtarget_top=getattr(args, "offtarget_top", 0),
        offtarget_max_mm=getattr(args, "offtarget_max_mm", 2),
        max_mm_left=args.max_mm_left,
        max_mm_right=args.max_mm_right,
        weights_path=args.weights,
        command=" ".join(sys.argv[1:]),
    )
    paths = pipeline.write_outputs(payload, args.out)
    _print_summary(payload, paths)
    return 0


def cmd_scan(args: argparse.Namespace) -> int:
    _check_fingers(args.fingers)
    lib = _load_library(args.library)
    hsl = pipeline.BASES_PER_FINGER * args.fingers

    if getattr(args, "genome", None):
        return _scan_genome(args, lib, hsl)

    records = ingest.read_input(sequence=args.seq, fasta=args.fasta, stdin=args.stdin)
    sites = targets.scan_many(records, spacers=_parse_spacers(args.spacers),
                              half_site_len=hsl,
                              both_strands=not args.single_strand)

    n_before = len(sites)
    if getattr(args, "repeats", None):
        idx = genome.parse_repeatmasker(args.repeats)
        sites = [s for s in sites
                 if not idx.overlaps(s.seq_name, s.start - 1, s.end)]
        print(f"重复序列过滤：{n_before} → {len(sites)}（剔除 {n_before - len(sites)} 个）")

    dz = designability.evaluate_many(sites, lib, max_tier=args.max_tier)
    dz = designability.triage(dz)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = [d.to_dict() for d in dz]
    report.write_json(out / "sites.json", rows)
    report.sites_tsv(out / "sites.tsv", rows)
    s = designability.summary(dz)
    print(f"位点 {s['n_sites']} 个：可设计(全Tier1) {s['by_status']['full']}，"
          f"含未验证 {s['by_status']['partial']}，缺口 {s['by_status']['gap']}")
    print(f"输出: {out / 'sites.tsv'}")
    return 0


def _scan_genome(args: argparse.Namespace, lib, hsl: int) -> int:
    progress = None if args.quiet else genome.stderr_progress
    payload = pipeline.run_genome(
        args.genome,
        args.out,
        library=lib,
        spacers=_parse_spacers(args.spacers),
        fingers=args.fingers,
        both_strands=not args.single_strand,
        max_tier=args.max_tier,
        repeats_path=args.repeats,
        regions_path=args.regions,
        limit_sites=args.limit,
        n_designs=args.design,
        max_design_sites=args.max_design_sites,
        max_per_finger=args.max_per_finger,
        template_key=args.template,
        architecture=args.architecture,
        weights_path=args.weights,
        chunk=args.chunk,
        ot_designs=args.ot_designs,
        ot_max_mm=args.ot_max_mm,
        ot_max_hits=args.ot_max_hits,
        progress=progress,
        command=" ".join(sys.argv[1:]),
    )
    ts = payload["targets_summary"]
    print("=" * 68)
    print(f"基因组: {payload['inputs'][0]['n_chroms']} 条序列，"
          f"{payload['inputs'][0]['total_bp']:,} bp")
    print(f"候选位点: {ts['n_sites']:,}  "
          f"(可设计 {ts['by_status']['full']:,} / 含未验证 {ts['by_status']['partial']:,} / "
          f"缺口 {ts['by_status']['gap']:,})")
    if ts.get("repeat_intervals"):
        print(f"重复区剔除: 注释 {ts['repeat_intervals']['n_intervals_merged']:,} 段 / "
              f"{ts['repeat_intervals']['total_bp']:,} bp")
    print(f"生成设计: {len(payload['designs'])}")
    for w in payload["warnings"]:
        print(f"  ! {w}")
    print("-" * 68)
    print(f"  sites.tsv   {Path(args.out) / 'sites.tsv'}")
    print(f"  report.json {Path(args.out) / 'report.json'}")
    if payload["designs"]:
        print(f"  report.html {Path(args.out) / 'report.html'}")
    if payload.get("genome_offtarget"):
        print(f"  offtargets.tsv {Path(args.out) / 'offtargets.tsv'}")
    print("=" * 68)
    return 0


def cmd_library(args: argparse.Namespace) -> int:
    if args.action == "info":
        lib = _load_library(args.path)
        print(json.dumps(lib.coverage_stats(), ensure_ascii=False, indent=2))
        cov = lib.coverage()
        print("\n三联体覆盖:", {k: v for k, v in cov.items()})
        return 0
    if args.action == "validate":
        lib = _load_library(args.path)
        print(f"OK：{len(lib)} 条记录，{len(lib.triplets)} 个三联体")
        return 0
    if args.action == "import-b1h":
        if not args.dir:
            raise SystemExit("import-b1h 需要 --dir")
        lib = FingerLibrary.import_b1h_dir(
            args.dir, top_n=args.top_n, pos4=args.pos4,
            split_position=args.split_position, min_freq=args.min_freq,
        )
        dest = lib.to_csv(args.out or "finger_library_b1h.csv")
        print(f"导入 {len(lib)} 条记录 → {dest}")
        print(json.dumps(lib.coverage_stats(), ensure_ascii=False, indent=2))
        return 0
    raise SystemExit(f"未知 action {args.action}")


def cmd_info(args: argparse.Namespace) -> int:
    print(f"ZFN_Design {__version__}")
    try:
        cov = load_default_library().coverage_stats()
        print(f"\n指头库: {cov['library_name']}  {cov['n_records']} 条记录，"
              f"覆盖 {cov['triplets_covered']}/{cov['triplets_total']} 三联体")
        if cov["triplets_covered"] < cov["triplets_total"]:
            print("  仅种子库，缺完整 B1H 表时结果不可用于真实设计（见 data/finger_library/README.md）")
    except LibraryError as exc:
        print(f"\n指头库: 加载失败 — {exc}")
    print("\n骨架模板:")
    for k in framework.template_keys():
        d = framework.describe(k)
        print(f"  {k}: {d['label']} ({d['template_length']} aa, {d['source_id']})")
    print("\n融合架构:")
    for a in fusion.describe_architectures():
        print(f"  {a['key']}: {a['label']} — {'需要成对位点' if a['requires_pair'] else '单位点即可'}")
    w = score.weights_documentation()
    print("\n打分权重:")
    for k, v in w["weights"].items():
        print(f"  {k:22s} {v:.2f}  {w['labels'].get(k, '')}")
    print(f"  合计 {sum(w['weights'].values()):.2f}")
    print("\n元件证据状态:")
    for k, p in fusion.parts().items():
        print(f"  {k:14s} {p['evidence_level']:14s} {p['label']}")
    return 0


def cmd_fusion(args: argparse.Namespace) -> int:
    from . import fusion_stage  # 延迟导入：S4 只在真要用时才装进来

    parts: dict[str, str] = {}
    if args.parts:
        parts = json.loads(Path(args.parts).read_text(encoding="utf-8"))
    try:
        res = fusion_stage.run_fusion(
            Path(args.dir),
            design_ids=args.designs or None,
            top=args.top,
            arch_keys=args.arch or None,
            parts=parts,
            seed=args.seed,
        )
    except ValueError as exc:
        print(f"错误：{exc}")
        return 2

    s = res["summary"]
    print("=" * 68)
    print(f"融合条目 {s['n_entries']}：合格 {s['n_ok']} / 需复核 {s['n_warn']} / "
          f"缺元件 {s['n_blocked']} / 跳过 {s['n_skipped']}")
    print(f"CDS 总长 {s['cds_total_bp']} bp，GC {s['gc_min']:.1%}–{s['gc_max']:.1%}")
    for e in res["entries"]:
        print(f"  {e['fusion_id']:44s} {e.get('architecture', '-'):18s} "
              f"{e.get('cds_status', '-'):8s} "
              f"{len(e.get('cds') or ''):5d} bp  GC={e.get('cds_gc')}")
    print("-" * 68)
    for k, v in (res.get("files") or {}).items():
        print(f"  {k:12s} {v}")
    print("=" * 68)
    return 0


def cmd_gui(args: argparse.Namespace) -> int:
    from . import gui  # 延迟导入：只有真的要开界面才去碰 http.server

    return gui.main([
        "--host", args.host,
        "--port", str(args.port),
        *(["--no-browser"] if args.no_browser else []),
        *(["--out-root", args.out_root] if getattr(args, "out_root", None) else []),
    ])


def _print_summary(payload: dict, paths: dict[str, Path]) -> None:
    lib = payload["library"]
    ts = payload["targets_summary"]
    designs = payload["designs"]
    print("=" * 68)
    print(f"指头库: {lib['library_name']}  {lib['n_records']} 条记录，"
          f"覆盖 {lib['triplets_covered']}/{lib['triplets_total']} 三联体")
    print(f"候选位点: {ts['n_sites']}  "
          f"(可设计 {ts['by_status']['full']} / 含未验证 {ts['by_status']['partial']} / 缺口 {ts['by_status']['gap']})")
    print(f"生成设计: {len(designs)}")
    for i, d in enumerate(designs[:10], start=1):
        sc = d.get("score") or {}
        print(f"  {i:2d}. {d['design_id']:28s} score={sc.get('total', 0):6.1f} "
              f"grade={sc.get('grade', '-')}  tier={d['min_evidence_tier']} "
              f"triplets={'>'.join(d['monomer_left']['triplets_N_to_C'])}")
    for w in payload.get("warnings", []):
        print(f"  ! {w}")
    print("-" * 68)
    for k, v in paths.items():
        print(f"  {k:12s} {v}")
    print("=" * 68)


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="zfn-design",
        description="从 DNA 序列出发的 ZFN/ZFP 靶点筛选、设计、打分与报告工具",
    )
    ap.add_argument("--version", action="version", version=f"zfn-design {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_run = sub.add_parser("run", help="跑完整流水线并输出报告")
    _add_input_args(p_run)
    _add_common(p_run)
    p_run.add_argument("--out", default="out", help="输出目录（默认 out）")
    p_run.add_argument("--repeats", metavar="PATH",
                       help="RepeatMasker 注释（.out/.bed）：先剔除重复区候选，再评分")
    p_run.add_argument("--genome", metavar="FASTA",
                       help="基因组 FASTA：按全基因组落点过滤候选（S2）")
    p_run.add_argument("--offtarget-top", type=int, default=0, metavar="N",
                       help="对 Design_score 前 N 个设计做全基因组脱靶评估（S3，默认 0=关闭，需 --genome）")
    p_run.add_argument("--offtarget-max-mm", type=int, default=2,
                       help="脱靶扫描时每侧允许的最大错配数（默认 2）")
    p_run.set_defaults(func=cmd_run)

    p_scan = sub.add_parser("scan", help="只做靶点扫描 + 可设计性判定（支持全基因组）")
    _add_input_args(p_scan, genome_ok=True)
    _add_common(p_scan)
    p_scan.add_argument("--out", default="out_scan")
    p_scan.add_argument("--repeats", metavar="PATH",
                        help="RepeatMasker 注释（.out 或 .bed）：落在重复序列里的候选直接剔除")
    p_scan.add_argument("--regions", metavar="BED",
                        help="只在 BED 指定的区域里扫（启动子/外显子等），大幅提速")
    p_scan.add_argument("--limit", type=int, default=200_000,
                        help="最多记录多少个候选就停（默认 200000，防止内存/时间爆掉）")
    p_scan.add_argument("--design", type=int, default=0, metavar="N",
                        help="顺带为前若干个可设计位点各出 N 个设计（默认 0=只扫不设计）")
    p_scan.add_argument("--max-design-sites", type=int, default=20,
                        help="最多给几个位点做设计（默认 20）")
    p_scan.add_argument("--chunk", type=int, default=genome.DEFAULT_CHUNK,
                        help=f"流式扫描的块大小（默认 {genome.DEFAULT_CHUNK:,} bp）")
    p_scan.add_argument("--ot-designs", type=int, default=10, metavar="N",
                        help="对 Design score 最高的前 N 个设计做全基因组脱靶扫描（默认 10，0=关闭）")
    p_scan.add_argument("--ot-max-mm", type=int, default=2,
                        help="全基因组脱靶允许的每侧错配数（默认 2）")
    p_scan.add_argument("--ot-max-hits", type=int, default=200,
                        help="每个设计最多记录多少条脱靶命中（默认 200，burden 仍按全量算）")
    p_scan.add_argument("--quiet", action="store_true", help="不打印每条染色体的进度")
    p_scan.set_defaults(func=cmd_scan)

    p_lib = sub.add_parser("library", help="指头库查看/校验/导入")
    p_lib.add_argument("action", choices=["info", "validate", "import-b1h"])
    p_lib.add_argument("--path", help="库文件路径")
    p_lib.add_argument("--dir", help="import-b1h 的输入目录")
    p_lib.add_argument("--out", help="import-b1h 的输出 CSV")
    p_lib.add_argument("--top-n", type=int, default=10)
    p_lib.add_argument("--pos4", default="L",
                       help="B1H 数据不含识别螺旋第 4 位，用它填（默认 L）")
    p_lib.add_argument("--split-position", default="any",
                       help="数据目录对应的阵列位置：F1/F2/F3/any")
    p_lib.add_argument("--min-freq", type=float, default=0.0)
    p_lib.set_defaults(func=cmd_library)

    p_info = sub.add_parser("info", help="打印模板/架构/权重/元件状态")
    p_info.set_defaults(func=cmd_info)

    p_fus = sub.add_parser("fusion", help="对已有结果目录做 S4 融合设计（独立于开始设计）")
    p_fus.add_argument("dir", metavar="DIR", help="已有结果目录（含 report.json）")
    p_fus.add_argument("--top", type=int, default=3, metavar="N",
                       help="取 Design_score 前 N 个设计做融合（默认 3）")
    p_fus.add_argument("--designs", nargs="*", metavar="ID",
                       help="只做这些设计 id（与 --top 二选一，给了就按 id 来）")
    p_fus.add_argument("--arch", nargs="*", metavar="KEY",
                       choices=list(fusion.DEFAULT_ARCHS) if hasattr(fusion, "DEFAULT_ARCHS")
                       else ["single_activator", "single_repressor", "paired_nuclease"],
                       help="架构：single_activator / single_repressor / paired_nuclease（默认全开）")
    p_fus.add_argument("--parts", metavar="JSON",
                       help="元件序列 JSON，如 {\"KRAB\":\"...\",\"FokI_CD\":\"...\"}")
    p_fus.add_argument("--seed", type=int, default=0, help="密码子优化的随机种子（默认 0）")
    p_fus.set_defaults(func=cmd_fusion)

    p_gui = sub.add_parser("gui", help="启动本地图形界面（浏览器打开）")
    p_gui.add_argument("--host", default="127.0.0.1", help="默认只对本机开放")
    p_gui.add_argument("--port", type=int, default=8797,
                       help="默认 8797；0 = 自动挑空闲端口")
    p_gui.add_argument("--no-browser", action="store_true", help="不自动打开浏览器")
    p_gui.add_argument("--out-root", default=None,
                       help="结果写到哪儿，默认 <项目>/out/gui/")
    p_gui.set_defaults(func=cmd_gui)
    return ap


def _utf8_stdio() -> None:
    """Force UTF-8 on stdout/stderr.

    Labels below contain non-ASCII characters. On a Windows console that
    defaults to a legacy code page, a redirected stdout raises
    UnicodeEncodeError and the process exits 1 -- which is how the Windows
    CI runner failed the CLI self-check.
    """
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError):
            pass


def main(argv: list[str] | None = None) -> int:
    _utf8_stdio()
    ap = build_parser()
    args = ap.parse_args(argv)
    try:
        return args.func(args)
    except (LibraryError, ingest.IngestError, framework.FrameworkError, fusion.FusionError) as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 2
    except FileNotFoundError as exc:
        print(f"错误: 找不到文件 {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
