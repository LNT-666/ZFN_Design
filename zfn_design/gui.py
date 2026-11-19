"""图形界面（本地网页）——真正调用本程序，不是另写一套算法。

设计原则：

* **只用标准库**：``http.server`` + ``json``，和包本身一样零第三方依赖；
* **跑的是同一份代码**：按钮背后直接调 :func:`zfn_design.pipeline.run`，
  和命令行 ``python -m zfn_design run`` 是同一个函数、同一份权重、同一份指头库；
* **只监听 127.0.0.1**：序列不出本机，也不做公网服务。

用法::

    python -m zfn_design.gui              # 起服务并自动打开浏览器
    python -m zfn_design.gui --port 8765  # 指定端口
    python -m zfn_design.gui --no-browser # 只起服务
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
import traceback
import uuid
import webbrowser
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, unquote, urlparse

from . import __version__, fusion_stage, manual, offtarget_stage, pipeline
from .ingest import IngestError, parse_fasta
from .models import SequenceRecord
from .sequtil import SequenceError, clean, validate


# --------------------------------------------------------------------------
# 一次运行的产物放哪儿
# --------------------------------------------------------------------------
class _State:
    """服务进程级的可变状态（最近一次运行的输出目录）。"""

    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.out_dir: Path | None = None
        self.run_count = 0

    def set_dir(self, d: Path) -> None:
        with self.lock:
            self.out_dir = d
            self.run_count += 1

    def get_dir(self) -> Path | None:
        with self.lock:
            return self.out_dir


STATE = _State()

SERVED_FILES = {"report.html", "report.json", "sites.tsv", "designs.tsv", "fusions.json",
                "zfn_bundle.zip", "offtarget_genome.tsv", "proteins.tsv", "cds.tsv",
                "cds.fasta", "order_form.csv"}


def default_out_root() -> Path:
    """结果默认落在哪：<项目>/out/gui/，装成 wheel 且写不进去时退到本地数据目录。

    可用环境变量 ``ZFN_OUT_ROOT`` 覆盖（命令行也给了 ``--out-root``）。
    """
    env = os.environ.get("ZFN_OUT_ROOT")
    if env:
        p = Path(env).expanduser()
        p.mkdir(parents=True, exist_ok=True)
        return p
    project = Path(__file__).resolve().parent.parent
    for cand in (project / "out" / "gui",
                 Path(os.environ.get("LOCALAPPDATA")
                      or Path.home()) / "ZFN_Design" / "out" / "gui"):
        try:
            cand.mkdir(parents=True, exist_ok=True)
            return cand
        except OSError:
            continue
    return Path(tempfile.mkdtemp(prefix="zfn_gui_"))


OUT_ROOT = default_out_root()


# ---------------------------------------------------------------- 参考数据目录
GENOME_EXT = {".fna", ".fa", ".fasta", ".fas"}
REPEAT_EXT = {".out", ".rm", ".align", ".bed"}
REGION_EXT = {".bed"}


def default_ref_dirs() -> list:
    """参考数据目录：默认 ``<项目>/refs/``，可用 ``ZFN_REF_DIRS``（冒号/分号分隔）覆盖。"""
    env = os.environ.get("ZFN_REF_DIRS")
    if env:
        cands = [Path(x).expanduser() for x in env.replace(";", ":").split(":") if x.strip()]
    else:
        cands = [Path(__file__).resolve().parent.parent.parent / "refs"]
    return [p for p in cands if p.is_dir()]


def _human_size(n: int) -> str:
    if n < 1024:
        return f"{n} B"
    for unit in ("KB", "MB", "GB", "TB"):
        n /= 1024.0
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}"
    return f"{n:.1f} TB"


def refs_payload() -> dict:
    """列出 refs 目录下可选的基因组 / RepeatMasker 注释 / 区域 BED，供页面下拉框使用。"""
    genomes, repeats, regions = [], [], []
    seen = {"g": set(), "r": set(), "b": set()}
    dirs = default_ref_dirs()
    for d in dirs:
        try:
            entries = sorted(d.iterdir(), key=lambda p: p.name.lower())
        except OSError:
            continue
        for p in entries:
            try:
                if not p.is_file():
                    continue
                size = p.stat().st_size
            except OSError:
                continue
            ext = p.suffix.lower()
            item = {"path": str(p), "name": p.name, "size": size,
                    "size_h": _human_size(size)}
            if ext in GENOME_EXT and str(p) not in seen["g"]:
                seen["g"].add(str(p))
                genomes.append(item)
            if ext in REPEAT_EXT and str(p) not in seen["r"]:
                seen["r"].add(str(p))
                repeats.append(item)
            if ext in REGION_EXT and str(p) not in seen["b"]:
                seen["b"].add(str(p))
                regions.append(item)
    return {"dirs": [str(d) for d in dirs], "genomes": genomes,
            "repeats": repeats, "regions": regions}


def make_run_dir(n_bp: int, run_id: str) -> Path:
    """给一次运行建一个「看得懂、找得到」的目录：20261114-153012_300bp_5ea42c。"""
    stem = f"{_stamp()}_{n_bp}bp" + (f"_{(run_id or '')[:6]}")
    d = OUT_ROOT / stem
    n = 1
    while d.exists():
        n += 1
        d = OUT_ROOT / f"{stem}-{n}"
    d.mkdir(parents=True)
    return d


def _stamp() -> str:
    return datetime.now().strftime("%Y%m%d-%H%M%S")


def _reveal(path: Path) -> tuple[bool, str]:
    """在系统文件管理器里打开目录。失败返回 (False, 原因)，不抛异常。"""
    try:
        if sys.platform.startswith("win"):
            os.startfile(str(path))  # noqa: S606  (Windows 专用)
        elif sys.platform == "darwin":
            subprocess.Popen(["open", str(path)])  # noqa: S603, S607
        else:
            subprocess.Popen(["xdg-open", str(path)])  # noqa: S603, S607
        return True, ""
    except Exception as exc:  # noqa: BLE001
        return False, f"{type(exc).__name__}: {exc}"


# --------------------------------------------------------------------------
# 请求处理逻辑（与 HTTP 解耦，方便直接单测）
# --------------------------------------------------------------------------
def _clamp(value: Any, lo: int, hi: int, default: int) -> int:
    try:
        n = int(value)
    except (TypeError, ValueError):
        return default
    return max(lo, min(hi, n))


def _records_from_request(req: dict[str, Any]) -> list[SequenceRecord]:
    # 从 Windows 复制粘贴的文本常带 BOM，str.strip() 不会去掉它，
    # 会让下面的 startswith(">") 失效、把 FASTA 头当成碱基。先清掉。
    text = str(req.get("seq") or "").lstrip("\ufeff").strip()
    if not text:
        raise IngestError("序列是空的")
    if text.lstrip().startswith(">"):
        return parse_fasta(text, source="gui")
    seq = validate(clean(text), name="input", allow_n=True)
    name = str(req.get("name") or "input").strip() or "input"
    return [SequenceRecord(name=name, seq=seq, source="gui")]


def run_request(req: dict[str, Any], *, progress=None) -> dict[str, Any]:
    """执行一次「Design Score」（S1 + 可选 S2），返回给前端渲染的精简结果。

    S3 不在这一步跑：宋济 2026-11-18 要求脱靶分析独立成按键。
    """
    records = _records_from_request(req)

    spacers = tuple(sorted({
        s for s in (_clamp(x, 0, 30, 6) for x in (req.get("spacers") or [5, 6, 7]))
    })) or (5, 6, 7)
    fingers = _clamp(req.get("fingers"), 2, 12, 3)
    n_designs = _clamp(req.get("n_designs"), 1, 10, 3)
    max_sites = _clamp(req.get("max_sites"), 1, 500, 20)
    max_mm_left = _clamp(req.get("max_mm_left"), 0, 3, 2)
    max_mm_right = _clamp(req.get("max_mm_right"), 0, 3, 2)
    architecture = str(req.get("architecture") or "paired_nuclease")

    # S2：全基因组重复过滤（宋成义 2026-11-18 口径）。不给基因组就退回 S1 纯打分模式。
    genome_fasta = (req.get("genome") or "").strip().strip('"') or None
    if genome_fasta and not os.path.exists(genome_fasta):
        raise ValueError(f"这个基因组路径不存在：{genome_fasta}")
    repeats_path = (req.get("repeats") or "").strip().strip('"') or None
    if repeats_path and not os.path.exists(repeats_path):
        raise ValueError(f"这个 RepeatMasker 注释不存在：{repeats_path}")
    if repeats_path and not genome_fasta:
        raise ValueError("要给 RepeatMasker 注释，必须先选基因组 FASTA（过滤按全基因组落点判）")

    t0 = time.perf_counter()
    payload = pipeline.run(
        records,
        spacers=spacers,
        fingers=fingers,
        n_designs=n_designs,
        max_sites=max_sites,
        max_mm_left=max_mm_left,
        max_mm_right=max_mm_right,
        architecture=architecture,
        genome_fasta=genome_fasta,
        repeats_path=repeats_path,
        offtarget_top=0,  # S3 独立按键（宋济 2026-11-18），这里不跑
        fusion_keys=[],  # S4 单独按键，开始设计不自动拼融合蛋白（宋成义 2026-11-18）
        command="zfn_design.gui",
        progress=progress,
    )
    elapsed_ms = (time.perf_counter() - t0) * 1000

    out_dir = make_run_dir(len(records[0].seq) if records else 0,
                           (payload.get("meta") or {}).get("run_id", ""))
    files = pipeline.write_outputs(payload, out_dir)
    STATE.set_dir(out_dir)

    ts = payload["targets_summary"]
    designs = []
    for i, d in enumerate(payload["designs"], start=1):
        sc = d.get("score") or {}
        ot = (d.get("offtarget") or {}).get("summary") or {}
        left = d["monomer_left"]
        right = d.get("monomer_right")
        designs.append({
            "id": d["design_id"],
            "rank": d.get("rank", i),
            "architecture": d.get("architecture"),
            "is_paired": d.get("is_paired", True),
            "target": d.get("target"),
            "triplets_left": left["triplets_N_to_C"],
            "triplets_right": (right or {}).get("triplets_N_to_C") or [],
            "protein": left["protein"],
            "partner_protein": (right or {}).get("protein", ""),
            "min_evidence_tier": d.get("min_evidence_tier"),
            "total": sc.get("total"),
            "grade": sc.get("grade"),
            "components": sc.get("components") or [],
            "score_warnings": sc.get("warnings") or [],
            "burden": ot.get("burden"),
            "n_hits": ot.get("n_hits"),
            "n_exact": ot.get("n_exact"),
        })

    sites = payload["sites"]
    return {
        "ok": True,
        "elapsed_ms": round(elapsed_ms, 1),
        "run_id": payload["meta"]["run_id"],
        "generated_at": payload["meta"]["generated_at"],
        "version": payload["meta"]["version"],
        "params": payload["meta"]["params"],
        "inputs": payload["inputs"],
        "library": payload["library"],
        "n_sites": len(sites),
        "fingers": fingers,
        "half_site_bp": pipeline.BASES_PER_FINGER * fingers,
        "spacers": list(spacers),
        "by_status": ts.get("by_status", {}),
        "genome_filter": ts.get("genome_filter"),
        "offtarget_genome": payload.get("offtarget_genome"),
        "sites": sites[:60],
        "triplet_demand": (ts.get("triplet_demand") or [])[:12],
        "designs": designs,
        "warnings": payload.get("warnings", []),
        "disclaimers": payload.get("disclaimers", []),
        "files": {k: f"/files/{v.name}" for k, v in files.items()},
        "out_dir": str(out_dir),
    }


def fusion_request(req: dict[str, Any], *, progress=None) -> dict[str, Any]:
    """S4：只针对已有作业做融合设计，不重跑 S1–S3（单独按键，宋成义 2026-11-18）。"""
    if progress:
        progress("读取设计", 0, 1, "读 report.json")
    raw_dir = str(req.get("dir") or "").strip().strip('"')
    if not raw_dir:
        raise ValueError("缺少作业目录（先点「开始设计」跑出结果）")
    out_dir = Path(raw_dir)
    if not (out_dir / "report.json").exists():
        raise ValueError(f"这个目录里没有 report.json：{out_dir}")

    top = _clamp(req.get("top"), 1, 20, 3)
    ids = [str(x).strip() for x in (req.get("design_ids") or []) if str(x).strip()]
    archs = [str(x).strip() for x in (req.get("archs") or []) if str(x).strip()]
    if not archs:
        archs = ["single_activator", "single_repressor", "paired_nuclease"]
    parts = req.get("parts") if isinstance(req.get("parts"), dict) else None
    seed = _clamp(req.get("seed"), 0, 10_000, 0)

    t0 = time.perf_counter()
    res = fusion_stage.run_fusion(out_dir, design_ids=ids or None, top=top,
                                  arch_keys=archs, parts=parts, seed=seed)
    STATE.set_dir(out_dir)
    res["elapsed_ms"] = round((time.perf_counter() - t0) * 1000, 1)
    res["out_dir"] = str(out_dir)
    if progress:
        progress("完成", 1, 1, "")
    return res


def offtarget_request(req: dict[str, Any], *, progress=None) -> dict[str, Any]:
    """S3：对已有作业的设计补算全基因组脱靶（单独按键，宋济 2026-11-18）。"""
    raw_dir = str(req.get("dir") or "").strip().strip('"')
    if not raw_dir:
        raise ValueError("缺少作业目录（先点「Design Score」跑出设计）")
    out_dir = Path(raw_dir)
    if not (out_dir / "report.json").exists():
        raise ValueError(f"这个目录里没有 report.json：{out_dir}")

    top = _clamp(req.get("top"), 0, 50, 10)
    max_mm = _clamp(req.get("max_mm"), 0, 3, 2)
    t0 = time.perf_counter()
    res = offtarget_stage.run_offtarget(out_dir, top=top, max_mm=max_mm, progress=progress)
    STATE.set_dir(out_dir)
    res["elapsed_ms"] = round((time.perf_counter() - t0) * 1000, 1)
    return res


def genome_request(req: dict[str, Any], *, progress=None) -> dict[str, Any]:
    """全基因组扫描：吃本机路径（3 Gb 不可能走浏览器上传）。"""
    from . import genome as _genome

    fasta = (req.get("genome") or "").strip().strip('"')
    if not fasta:
        raise ValueError("请填基因组 FASTA 的完整路径")
    if not os.path.exists(fasta):
        raise ValueError(f"这个路径不存在：{fasta}")

    def _opt(key: str) -> str | None:
        v = (req.get(key) or "").strip().strip('"')
        if not v:
            return None
        if not os.path.exists(v):
            raise ValueError(f"这个路径不存在：{v}")
        return v

    repeats = _opt("repeats")
    skip_repeats = req.get("skip_repeats")
    skip_repeats = True if skip_repeats is None else bool(skip_repeats)
    repeats_used = repeats if (repeats and skip_repeats) else None
    regions = _opt("regions")
    limit = _clamp(req.get("limit"), 1, 2_000_000, 200_000)
    spacers = tuple(sorted({
        s for s in (_clamp(x, 0, 30, 6) for x in (req.get("spacers") or [5, 6, 7]))
    })) or (5, 6, 7)
    fingers = _clamp(req.get("fingers"), 2, 12, 3)
    architecture = str(req.get("architecture") or "paired_nuclease")
    n_designs = _clamp(req.get("n_designs"), 0, 10, 0)
    ot_designs = _clamp(req.get("ot_designs"), 0, 20, 10)
    ot_max_mm = _clamp(req.get("ot_max_mm"), 0, 3, 2)

    out_dir = OUT_ROOT / f"{_stamp()}_genome_{limit}"
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()
    payload = pipeline.run_genome(
        fasta, out_dir,
        spacers=spacers, fingers=fingers, architecture=architecture,
        repeats_path=repeats_used, regions_path=regions,
        limit_sites=limit, n_designs=n_designs, max_design_sites=10,
        ot_designs=ot_designs,
        ot_max_mm=ot_max_mm,
        command="zfn_design.gui",
    )
    elapsed_ms = (time.perf_counter() - t0) * 1000
    STATE.set_dir(out_dir)
    if repeats and not skip_repeats:
        payload.setdefault("warnings", []).append(
            "已按你的选择关闭重复序列过滤：RepeatMasker 注释未参与剔除（multi_copy 查重仍然生效）")

    ts = payload["targets_summary"]
    return {
        "ok": True,
        "genome_mode": True,
        "elapsed_ms": round(elapsed_ms, 1),
        "run_id": payload["meta"]["run_id"],
        "generated_at": payload["meta"]["generated_at"],
        "version": payload["meta"]["version"],
        "library": payload["library"],
        "inputs": payload["inputs"],
        "n_sites": ts["n_sites"],
        "n_multi_copy": ts.get("n_multi_copy", 0),
        "n_sites_raw": ts.get("n_sites_raw"),
        "n_repeat_filtered": ts.get("n_repeat_filtered"),
        "n_unique_sites": ts.get("n_unique_sites"),
        "n_removed_ge3": ts.get("n_removed_ge3"),
        "n_removed_ge3_sites": ts.get("n_removed_ge3_sites"),
        "skip_repeats": skip_repeats,
        "by_status": ts.get("by_status", {}),
        "scan_stats": ts.get("scan_stats", {}),
        "repeat_intervals": ts.get("repeat_intervals"),
        "estimated_total": ts.get("estimated_total"),
        "truncated": ts.get("truncated"),
        "sites": payload["sites"][:60],
        "designs": [{"id": d["design_id"],
                     "design_score": (d.get("score") or {}).get("total"),
                     "grade": (d.get("score") or {}).get("grade"),
                     "offtarget_score": (d.get("offtarget_score") or {}).get("score")}
                    for d in payload["designs"]],
        "genome_offtarget": payload.get("genome_offtarget") or {},
        "warnings": payload.get("warnings", []),
        "disclaimers": payload.get("disclaimers", []),
        "out_dir": str(out_dir),
        "sites_tsv": str(out_dir / "sites.tsv"),
        "offtargets_tsv": (str(out_dir / "offtargets.tsv")
                           if (out_dir / "offtargets.tsv").exists() else None),
    }


def info_payload() -> dict[str, Any]:
    from .framework import template_keys
    from .fusion import describe_architectures
    from .library import load_default_library
    from .score import weights_documentation

    lib = load_default_library()
    w = weights_documentation()
    return {
        "version": __version__,
        "library": lib.coverage_stats(),
        "templates": template_keys(),
        "architectures": [{"key": a["key"], "label": a.get("label", a["key"])}
                          for a in describe_architectures()],
        "n_weights": len(w["weights"]),
    }


def demo_payload() -> dict[str, str]:
    """示例序列：就是 examples/demo_target.fa 里那条。"""
    path = Path(__file__).resolve().parent.parent / "examples" / "demo_target.fa"
    if path.exists():
        recs = parse_fasta(path.read_text(encoding="utf-8"), source=str(path))
        r = recs[0]
        return {"name": r.name, "seq": r.seq}
    # 装成 wheel 后 examples 不一定在，给一条同样带 GCGTGGGCG-GATATC-GCGTGGGCG 的短序列
    seq = ("TACGTAGCTAGCTTAGCTAGGCTA" + "GCGTGGGCGGATATCCGCCCACGC" + "TTACGGATCCTTAGGCAT")
    return {"name": "demo_target", "seq": seq}



# --------------------------------------------------------------------------
# 后台作业：长任务不再占着一条 HTTP 连接（反代 90 s 会 504）
# --------------------------------------------------------------------------
_JOBS: dict[str, dict[str, Any]] = {}
_JOBS_LOCK = threading.Lock()
_JOB_KEEP = 32


def _job_start(kind: str, fn) -> str:
    """起一个后台线程跑 fn(progress)，立刻返回 job id。

    ``progress(phase, done, total, message="")`` 由作业内部回调，前端轮询
    ``/api/job`` 时读 ``progress`` 字段画进度条。
    """
    jid = uuid.uuid4().hex[:12]
    rec: dict[str, Any] = {
        "id": jid, "kind": kind, "state": "running",
        "t0": time.time(), "result": None, "error": None, "kind_detail": None,
        "progress": {"phase": "排队", "done": 0, "total": 0, "frac": 0.0, "message": ""},
    }

    def _progress(phase: str, done: float = 0, total: float = 0, message: str = "") -> None:
        frac = (float(done) / float(total)) if total else 0.0
        rec["progress"] = {
            "phase": phase,
            "done": round(float(done), 3),
            "total": round(float(total), 3),
            "frac": round(min(max(frac, 0.0), 1.0), 4),
            "message": message,
        }
    with _JOBS_LOCK:
        _JOBS[jid] = rec
        if len(_JOBS) > _JOB_KEEP:  # 只留最近的一批，免得常驻进程越攒越多
            for k in sorted(_JOBS, key=lambda x: _JOBS[x]["t0"])[: len(_JOBS) - _JOB_KEEP]:
                _JOBS.pop(k, None)

    def _work() -> None:
        try:
            rec["result"] = fn(_progress)
            rec["state"] = "done"
            _progress("完成", 1, 1)
        except ValueError as exc:
            rec["state"] = "error"
            rec["error"] = str(exc)
            rec["kind_detail"] = "input"
        except (IngestError, SequenceError) as exc:
            rec["state"] = "error"
            rec["error"] = str(exc)
            rec["kind_detail"] = "input"
        except Exception as exc:  # noqa: BLE001
            traceback.print_exc()
            rec["state"] = "error"
            rec["error"] = f"{type(exc).__name__}: {exc}"
            rec["kind_detail"] = "internal"
        finally:
            rec["elapsed_ms"] = int((time.time() - rec["t0"]) * 1000)

    threading.Thread(target=_work, daemon=True, name=f"job-{kind}").start()
    return jid


def _job_view(jid: str) -> dict[str, Any]:
    with _JOBS_LOCK:
        rec = _JOBS.get(jid)
        if rec is None:
            return {"ok": False, "error": "没有这个作业（可能已被清理，或 id 打错）"}
        out = {k: rec[k] for k in ("id", "kind", "state", "result", "error", "kind_detail")}
        out["ok"] = True
        out["progress"] = rec.get("progress") or {"phase": "", "frac": 0.0}
        out["elapsed_ms"] = rec.get("elapsed_ms") or int((time.time() - rec["t0"]) * 1000)
        return out


# --------------------------------------------------------------------------
# HTTP 层
# --------------------------------------------------------------------------
class _Handler(BaseHTTPRequestHandler):
    server_version = f"ZFN-Design-GUI/{__version__}"
    protocol_version = "HTTP/1.1"

    # 关掉默认的逐请求日志（控制台只留一行"已启动"）
    def log_message(self, fmt: str, *args: Any) -> None:  # noqa: A003
        return

    def _send(self, code: int, body: bytes, ctype: str,
              extra: dict[str, str] | None = None) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, code: int, obj: Any) -> None:
        self._send(code, json.dumps(obj, ensure_ascii=False).encode("utf-8"),
                   "application/json; charset=utf-8")

    def do_GET(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        if path in ("/", "/index.html"):
            self._send(200, INDEX_HTML.encode("utf-8"), "text/html; charset=utf-8")
        elif path == "/api/info":
            self._json(200, info_payload())
        elif path == "/api/demo":
            self._json(200, demo_payload())
        elif path == "/api/refs":
            self._json(200, refs_payload())
        elif path == "/api/job":
            qs = parse_qs(urlparse(self.path).query)
            self._json(200, _job_view((qs.get("id") or [""])[0]))
        elif path == "/manual":
            qs = parse_qs(urlparse(self.path).query)
            page = manual.render_page((qs.get("lang") or ["zh"])[0])
            self._send(200, page.encode("utf-8"), "text/html; charset=utf-8")
        elif path == "/manual/download":
            qs = parse_qs(urlparse(self.path).query)
            lang = (qs.get("lang") or ["zh"])[0]
            fmt = (qs.get("fmt") or ["docx"])[0]
            self._serve_manual(lang, fmt)
        elif path.startswith("/files/"):
            self._serve_file(unquote(path[len("/files/"):]))
        else:
            self._json(404, {"ok": False, "error": "没有这个地址"})

    def do_POST(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        if path == "/api/quit":
            self._json(200, {"ok": True})
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        if path == "/api/reveal":
            d = STATE.get_dir()
            if d is None:
                self._json(409, {"ok": False, "error": "还没有运行过"})
                return
            ok, msg = _reveal(d)
            self._json(200 if ok else 500, {"ok": ok, "dir": str(d), "error": msg})
            return
        if path == "/api/genome":
            length = _clamp(self.headers.get("Content-Length"), 0, 1_000_000, 0)
            raw = self.rfile.read(length) if length else b"{}"
            try:
                req = json.loads(raw.decode("utf-8") or "{}")
            except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                self._json(400, {"ok": False, "error": f"请求体不是合法 JSON：{exc}"})
                return
            jid = _job_start("genome", lambda progress: genome_request(req, progress=progress))
            self._json(200, {"ok": True, "job_id": jid, "async": True})
            return
        if path == "/api/fusion":
            length = _clamp(self.headers.get("Content-Length"), 0, 1_000_000, 0)
            raw = self.rfile.read(length) if length else b"{}"
            try:
                req = json.loads(raw.decode("utf-8") or "{}")
            except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                self._json(400, {"ok": False, "error": f"请求体不是合法 JSON：{exc}"})
                return
            jid = _job_start("fusion", lambda progress: fusion_request(req, progress=progress))
            self._json(200, {"ok": True, "job_id": jid, "async": True})
            return
        if path == "/api/offtarget":
            length = _clamp(self.headers.get("Content-Length"), 0, 1_000_000, 0)
            raw = self.rfile.read(length) if length else b"{}"
            try:
                req = json.loads(raw.decode("utf-8") or "{}")
            except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                self._json(400, {"ok": False, "error": f"请求体不是合法 JSON：{exc}"})
                return
            jid = _job_start("offtarget", lambda progress: offtarget_request(req, progress=progress))
            self._json(200, {"ok": True, "job_id": jid, "async": True})
            return
        if path != "/api/run":
            self._json(404, {"ok": False, "error": "没有这个地址"})
            return

        length = _clamp(self.headers.get("Content-Length"), 0, 8_000_000, 0)
        raw = self.rfile.read(length) if length else b"{}"
        try:
            req = json.loads(raw.decode("utf-8") or "{}")
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            self._json(400, {"ok": False, "error": f"请求体不是合法 JSON：{exc}"})
            return

        jid = _job_start("run", lambda progress: run_request(req, progress=progress))
        self._json(200, {"ok": True, "job_id": jid, "async": True})
        return

    def _serve_file(self, name: str) -> None:
        d = STATE.get_dir()
        if name not in SERVED_FILES or d is None:
            self._json(404, {"ok": False, "error": "文件不存在（先跑一次）"})
            return
        path = d / name
        if not path.exists():
            self._json(404, {"ok": False, "error": "文件不存在"})
            return
        ctype = {
            ".html": "text/html; charset=utf-8",
            ".json": "application/json; charset=utf-8",
            ".tsv": "text/tab-separated-values; charset=utf-8",
        }.get(path.suffix, "application/octet-stream")
        self._send(200, path.read_bytes(), ctype,
                   {"Content-Disposition": f'inline; filename="{name}"'})

    def _serve_manual(self, lang: str, fmt: str) -> None:
        lang = lang if lang in manual.LANGS else "zh"
        fmt = fmt if fmt in ("docx", "md") else "docx"
        path = manual.manual_path(lang).with_suffix("." + fmt)
        if not path.exists():
            self._json(404, {"ok": False, "error": "说明书文件缺失"})
            return
        ctype = ("application/vnd.openxmlformats-officedocument.wordprocessingml.document"
                 if fmt == "docx" else "text/markdown; charset=utf-8")
        name = f"ZFN-manual-{lang}.{fmt}"
        self._send(200, path.read_bytes(), ctype,
                   {"Content-Disposition": f'attachment; filename="{name}"'})


# --------------------------------------------------------------------------
# 页面
# --------------------------------------------------------------------------
INDEX_HTML = r"""<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ZFN 设计工作台</title>
<style>
  :root{
    --bg:#0d1016; --panel:#141922; --panel2:#10151d; --line:#232c3a;
    --txt:#dfe7f3; --dim:#8d9bb0; --dim2:#6b7889;
    --acc:#5aa2ff; --ok:#4fd08a; --warn:#e8c56a; --bad:#f08a7f;
    --mono:ui-monospace,SFMono-Regular,Consolas,"Liberation Mono",monospace;
  }
  *{box-sizing:border-box}
  body{margin:0;background:var(--bg);color:var(--txt);
    font:14px/1.65 system-ui,-apple-system,"Segoe UI","Microsoft YaHei",sans-serif}
  .wrap{max-width:1240px;margin:0 auto;padding:0 22px 70px}
  header{border-bottom:1px solid var(--line);background:linear-gradient(180deg,#131a25,#0d1016)}
  header .wrap{padding:20px 22px 16px}
  h1{margin:0;font-size:21px;letter-spacing:.2px}
  h1 small{color:var(--dim);font-weight:400;font-size:13px;margin-left:10px}
  .sub{color:var(--dim);font-size:12.5px;margin-top:6px}
  .grid{display:grid;grid-template-columns:minmax(320px,420px) 1fr;gap:18px;margin-top:20px}
  @media(max-width:900px){.grid{grid-template-columns:1fr}}
  .card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:16px}
  label.f{display:block;color:var(--dim);font-size:12.5px;margin:12px 0 5px}
  textarea{width:100%;height:190px;resize:vertical;background:#0a0e14;color:#cfe0f5;
    border:1px solid var(--line);border-radius:8px;padding:10px 12px;
    font-family:var(--mono);font-size:12px;line-height:1.7}
  input[type=text],input[type=number],select{
    background:#0a0e14;color:#cfe0f5;border:1px solid var(--line);border-radius:7px;
    padding:7px 9px;font:inherit;font-size:13px}
  input[type=number]{width:78px}
  .row{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
  .chk{display:inline-flex;align-items:center;gap:6px;background:#0f151e;border:1px solid var(--line);
    border-radius:7px;padding:5px 10px;font-size:12.5px;color:#c3d2e6;cursor:pointer}
  .chk input{accent-color:var(--acc)}
  button{font:inherit;font-size:13.5px;border-radius:8px;padding:9px 18px;cursor:pointer;
    border:1px solid #2a4b7c;background:linear-gradient(180deg,#20406b,#183049);
    color:#dce9ff;font-weight:600}
  button:hover{filter:brightness(1.12)}
  button:disabled{opacity:.5;cursor:default;filter:none}
  button.ghost{background:#111823;border-color:var(--line);color:#b9c8dc;font-weight:400}
  button.tab{background:#0f1620;border:1px solid var(--line);color:#93a4bb;font-weight:500;
    border-radius:8px 8px 0 0;padding:8px 14px}
  button.tab.on{background:#1b2838;border-color:#3d6ea8;color:#eaf2ff;font-weight:600}
  .stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(112px,1fr));gap:10px}
  .stat{background:var(--panel2);border:1px solid var(--line);border-radius:9px;padding:10px 12px}
  .stat .k{color:var(--dim2);font-size:11.5px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  .stat .v{font-size:20px;font-weight:600;margin-top:2px}
  .v.ok{color:var(--ok)} .v.warn{color:var(--warn)} .v.bad{color:var(--bad)}
  table{width:100%;border-collapse:collapse;font-size:12.5px}
  th{text-align:left;color:var(--dim);font-weight:600;font-size:12px;
    border-bottom:1px solid var(--line);padding:7px 8px;background:#111823}
  td{padding:6px 8px;border-bottom:1px solid #1a2130;vertical-align:top}
  td.mono,.mono{font-family:var(--mono);font-size:12px}
  .design table td:first-child{white-space:nowrap}
  tr:hover td{background:#141b26}
  .pill{display:inline-block;font-size:11px;padding:1px 7px;border-radius:5px;border:1px solid}
  .pill.full{color:#8ce7b6;border-color:#2a5b41;background:#12281c}
  .pill.gap{color:#f0a89f;border-color:#5b2a26;background:#1e1211}
  .pill.warnp{color:#f0d094;border-color:#5b4a24;background:#1e1a11}
  .design{background:var(--panel2);border:1px solid var(--line);border-radius:9px;padding:12px 14px;margin:10px 0}
  .design .h{display:flex;justify-content:space-between;gap:12px;align-items:baseline;flex-wrap:wrap}
  .design .id{font-family:var(--mono);font-size:13px}
  .design .sc{font-family:var(--mono);font-size:19px;font-weight:700}
  .grade{font-size:11px;padding:1px 7px;border-radius:5px;margin-left:7px;vertical-align:2px}
  .gA{background:#12281c;color:#5fe0a0;border:1px solid #2a5b41}
  .gB{background:#201c10;color:#e8c56a;border:1px solid #5b4a24}
  .gC{background:#1e1412;color:#e8a06a;border:1px solid #5b3a24}
  .gD{background:#1e1211;color:#f0a89f;border:1px solid #5b2a26}
  .seq{background:#0a0e14;border:1px solid var(--line);border-radius:7px;
    padding:8px 10px;font-family:var(--mono);font-size:11.5px;color:#a9bad4;
    overflow-wrap:anywhere;line-height:1.75;margin-top:6px}
  details summary{cursor:pointer;color:var(--dim);font-size:12.5px;margin:6px 0}
  .note{font-size:12px;color:var(--dim2)}
  .banner{border:1px solid #5b4a24;background:#1b1710;color:#e7cf9f;
    border-radius:9px;padding:11px 14px;font-size:12.5px;margin:16px 0}
  .err{border:1px solid #5b2a26;background:#1b1010;color:#f0a89f;
    border-radius:8px;padding:10px 12px;font-size:13px;margin:10px 0;white-space:pre-wrap}
  .files{display:flex;gap:10px;flex-wrap:wrap;margin-top:10px}
  .limit{margin-top:14px;padding:11px 14px;border-radius:8px;border:1px solid #7a5a1e;
    background:rgba(122,90,30,.14);color:#f0d9a6;font-size:12.5px;line-height:1.62}
  .limit b{color:#ffd98a}
  .files a{font-size:12.5px;color:#cfe0f5;text-decoration:none;border:1px solid var(--line);
    border-radius:7px;padding:6px 12px;background:#101823}
  .files a:hover{border-color:#2a4b7c;color:#fff}
  .files a.zip{border-color:#2f6b45;background:#12241a;color:#a9e5bd}
  .files a.zip:hover{border-color:#3f8c5c;color:#d6f6e2}
  .s4{margin-top:16px;border:1px solid #2b4a6f;border-radius:9px;padding:12px 14px;
    background:linear-gradient(180deg,rgba(33,66,102,.22),rgba(16,24,35,.2))}
  .s4h{font-weight:600;color:#bcd6f5;font-size:13.5px}
  .s4row{display:flex;gap:14px;align-items:center;flex-wrap:wrap;margin-top:9px}
  .s4row input[type=number]{background:#0e1620;border:1px solid var(--line);color:#e6eef8;
    border-radius:6px;padding:4px 7px}
  .s4row button{background:#1c4a74;border:1px solid #2f6b9f;color:#e8f2ff;border-radius:7px;
    padding:7px 14px;cursor:pointer;font-size:12.5px}
  .s4row button:hover{background:#245c8e}
  .s4row button:disabled{opacity:.55;cursor:default}
  .s4 textarea{width:100%;box-sizing:border-box;margin-top:8px;background:#0e1620;
    border:1px solid var(--line);color:#e6eef8;border-radius:6px;padding:7px 9px;
    font-family:ui-monospace,Consolas,monospace;font-size:11.5px;line-height:1.5}
  .fu{border:1px solid var(--line);border-radius:8px;padding:10px 12px;margin-top:10px;background:#0d151f}
  .fu .h{display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap;align-items:baseline}
  .fu .id{font-size:12.5px;color:#cfe0f5}
  .fu .pill2{font-size:11.5px;padding:1px 7px;border-radius:5px;border:1px solid var(--line)}
  .fu .ok2{color:#a9e5bd;border-color:#2f6b45;background:#12241a}
  .fu .wn2{color:#f0d9a6;border-color:#7a5a1e;background:#241f12}
  .fu .bd2{color:#ffb0b0;border-color:#7a2a2a;background:#241414}
  .fu .seq{max-height:120px;overflow:auto;font-family:ui-monospace,Consolas,monospace;
    font-size:11.5px;color:#b9c9dc;word-break:break-all;line-height:1.55;margin-top:6px}
  .outdir{display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin-top:10px;
    border-top:1px dashed var(--line);padding-top:10px}
  .outdir code{font-size:12px;color:#cfe0f5;background:#0d141d;border:1px solid var(--line);
    border-radius:6px;padding:5px 9px;flex:1 1 260px;min-width:0;max-width:60ch;
    white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  button.mini{font-size:12px;padding:5px 11px;font-weight:500;margin-left:auto;white-space:nowrap}
  .spin{display:inline-block;width:13px;height:13px;border:2px solid #2a4b7c;
    border-top-color:#7fb2ff;border-radius:50%;animation:sp .7s linear infinite;vertical-align:-2px}
  @keyframes sp{to{transform:rotate(360deg)}}
  .hidden{display:none}
  .muted{color:var(--dim2);font-size:12px}
  .pbar{height:6px;border-radius:4px;background:#16243c;overflow:hidden;margin-top:8px;border:1px solid var(--line)}
  .pbar .pfill{height:100%;width:0;background:#3f7fd0;transition:width .4s ease}
  .pbar.on .pfill{background:linear-gradient(90deg,#3f7fd0,#7fb2ff)}
  .ptxt{margin-top:4px;font-size:12px;color:var(--dim2);min-height:16px}
  .warnline{background:#2a2312;border-left:3px solid #c8a24a;padding:7px 9px;border-radius:5px;color:#e6d7b0}
  details.opt{margin-top:12px;border:1px solid var(--line);border-radius:8px;padding:8px 10px;background:#0e1725}
  details.opt>summary{cursor:pointer;font-size:13px;color:#c9d8ef}
  footer{margin-top:34px;border-top:1px solid var(--line);padding-top:14px;color:var(--dim2);font-size:12px}
  .langbtn{position:fixed;top:12px;left:12px;z-index:20;font-size:12px;font-weight:600;
    padding:6px 12px;border-radius:7px;border:1px solid #2a4b7c;background:#16233a;
    color:#cfe0f5;cursor:pointer;letter-spacing:.3px}
  .langbtn:hover{background:#1d2e49;color:#fff}
  .doclink{position:fixed;top:12px;right:12px;z-index:20;font-size:12px;font-weight:600;
    padding:6px 12px;border-radius:7px;border:1px solid #2a4b7c;background:#16233a;
    color:#cfe0f5;text-decoration:none;letter-spacing:.3px}
  .doclink:hover{background:#1d2e49;color:#fff}
</style>
</head>
<body>

<header><div class="wrap">
  <button class="langbtn" id="lang" type="button" title="Switch language / 切换语言">EN</button>
  <a class="doclink" id="doclink" href="/manual?lang=zh">使用说明书</a>
  <h1>ZFN 设计工作台 <small>本地图形界面 · 后端就是 <span class="mono">zfn_design</span> 本体</small></h1>
  <div class="sub" id="env">正在读取引擎信息…</div>
</div></header>

<div class="wrap">
  <div class="banner">
    按钮背后直接调用 <span class="mono">zfn_design.pipeline.run()</span>——和命令行
    <span class="mono">python -m zfn_design run</span> 是同一个函数、同一份权重、同一份指头库。
    分数只表示「按公开设计规则看更像靠谱位点」，不是活性预测；任何设计进实验室前必须 T7E1 或
    amplicon 深度测序验证。
  </div>

  <div class="grid">
    <!-- 左：输入 -->
    <div class="card">
      <div style="font-weight:600">输入</div>
      <div class="note" id="modehint" style="margin:8px 0 10px">
        粘贴一条序列，只设计这一条，输出按 Design_score 排序的 all-target 报告。
        重复序列过滤在下一步接进来。
      </div>
      <div id="pane-single">

      <label class="f">DNA 序列（纯文本或 FASTA 都可以）</label>
      <textarea id="seq" spellcheck="false" placeholder="粘贴 DNA 序列，或点下面的「载入示例」"></textarea>

      <div class="row" style="margin-top:8px">
        <button class="ghost" id="demo">载入示例</button>
        <button class="ghost" id="clear">清空</button>
        <span class="muted" id="seqinfo"></span>
      </div>

      <label class="f">构型（决定单体型 / 组合型）</label>
      <div class="row">
        <select id="arch" style="flex:1;min-width:0"></select>
      </div>
      <div class="note" id="archnote" style="margin-top:6px"></div>

      <label class="f" id="spacerlabel">spacer 长度（两半位点之间）</label>
      <div class="row" id="spacerrow">
        <label class="chk"><input type="radio" name="spmode" value="default" checked> 默认 5–7</label>
        <label class="chk"><input type="radio" name="spmode" value="custom"> 自定义</label>
        <input type="text" id="spcustom" value="6" placeholder="6 或 5,6,8 或 5-7"
               style="display:none;flex:1;min-width:0">
      </div>

      <label class="f">指头阵列</label>
      <div class="row">
        <select id="fingersel" style="flex:1;min-width:0">
          <option value="3" selected>三指阵列（读 9 bp）</option>
          <option value="6">六指阵列（读 18 bp，两个三指模块串联）</option>
          <option value="custom">自定义指数（Custom）</option>
        </select>
        <input type="number" id="fingcustom" value="4" min="2" max="12" step="1"
               title="2–12 指；非 3/6 用通用拼接骨架"
               style="display:none;width:74px">
      </div>
      <div class="note" id="archnote2" style="margin-top:4px">
        需要基因组单拷贝结合时再加长：三指单独不唯一，六指/更长阵列能提高唯一性，但模块间界面从未验证过。
      </div>

      <label class="f">参数</label>
      <div class="row">
        <span class="note">每靶点数</span>
        <input type="number" id="ndesign" value="3" min="1" max="10">
        <span class="note">最多处理位点</span>
        <input type="number" id="maxsites" value="20" min="1" max="500">
      </div>
      <div class="row" style="margin-top:8px">
        <span class="note">每侧最大错配</span>
        <input type="number" id="mm" value="2" min="0" max="3">
      </div>

      <details class="opt">
        <summary>Target filter · 全基因组重复过滤（可选，默认关闭）</summary>
        <div class="note" style="margin:6px 0 4px">
          选了基因组就按基因组口径过滤：候选 motif 的落点压在 RepeatMasker 注释区，
          或在基因组里出现 ≥3 次，整条剔除、<b>不参与打分</b>。只留通过过滤的候选做 DesignScore。
        </div>
        <div class="row">
          <select id="genome" style="flex:1;min-width:0">
            <option value="">（不过滤：只做 DesignScore）</option>
          </select>
        </div>
        <div class="row" style="margin-top:6px">
          <select id="repeats" style="flex:1;min-width:0">
            <option value="">（不叠 RepeatMasker 注释）</option>
          </select>
        </div>
      </details>

      <div class="row" style="margin-top:16px">
        <button id="run">Design Score（S1 + S2）</button>
        <span id="status" class="muted"></span>
      </div>
      <div class="pbar" id="run-bar"><div class="pfill" id="run-fill"></div></div>
      <div class="ptxt" id="run-txt"></div>

      <label class="f" style="margin-top:16px">Off-target Score（S3 · 单独执行）</label>
      <div class="note" style="margin-bottom:4px">
        对上一次 Design Score 出来的结果，取 Design_score 前 N 个设计做全基因组脱靶扫描
        （每侧 0–2 错配）。结果<b>并列展示，不加进 Design_score</b>。每扫一遍全基因组约 1–2 分钟。
      </div>
      <div class="note warnline" style="margin-bottom:6px">
        <b>前提：只有通过 Target filter（全基因组重复过滤）的候选，才能做 Off-target Score。</b>
        所以请先在「Target filter」里选基因组 FASTA（可选再加 RepeatMasker 注释），跑一次 Design Score，
        再点本按钮。没走过过滤的作业会被拒绝，并提示原因。
      </div>
      <div class="row">
        <span class="note">Top N</span>
        <input type="number" id="ot-top" value="10" min="1" max="50">
        <span class="note">每侧错配上限</span>
        <input type="number" id="ot-mm" value="2" min="0" max="3">
      </div>
      <div class="row" style="margin-top:8px">
        <button id="ot-run">Off-target Score（S3）</button>
      </div>
      <div class="pbar" id="ot-bar"><div class="pfill" id="ot-fill"></div></div>
      <div class="ptxt" id="ot-txt"></div>

      <div class="note" style="margin-top:10px">
        序列只在本机内存里处理，不联网、不上传。
      </div>

      </div><!-- /pane-single -->

    </div>

    <!-- 右：结果 -->
    <div class="card" id="result">
      <div class="muted">运行一次，这里会出现统计、设计清单和完整报告。</div>
    </div>
  </div>

  <footer>
    引擎 <span class="mono" id="ver"></span> ·
    <span class="mono">python -m zfn_design.gui</span> ·
    输出文件在临时目录，关掉服务不会自动删除。
  </footer>
</div>

<script>
const $ = s => document.querySelector(s);
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));

/* ---------- i18n：左上角按钮切换 中文 / English ---------- */
const I18N_EN = {
  "ZFN 设计工作台":"ZFN Design Workbench",
  "使用说明书":"User manual",
  "本地图形界面 · 后端就是":"Local web UI · powered by",
  "本体":"itself",
  "正在读取引擎信息…":"Reading engine info…",
  "按钮背后直接调用":"The buttons call",
  "——和命令行":"— identical to the command line",
  "是同一个函数、同一份权重、同一份指头库。 分数只表示「按公开设计规则看更像靠谱位点」，不是活性预测；任何设计进实验室前必须 T7E1 或 amplicon 深度测序验证。":"same function, same weights, same finger library. Scores only mean “looks plausible under the public design rules”, not predicted activity; every design must be validated by T7E1 or amplicon deep sequencing before it enters the lab.",
  "输入":"Input",
  "粘贴一条序列，只设计这一条，输出按 Design_score 排序的 all-target 报告。 重复序列过滤在下一步接进来。":"Paste one sequence; only this one gets designed. Output is an all-target report ranked by Design_score.",
  "DNA 序列（纯文本或 FASTA 都可以）":"DNA sequence (plain text or FASTA)",
  "粘贴 DNA 序列，或点下面的「载入示例」":"Paste a DNA sequence, or click “Load example” below",
  "载入示例":"Load example",
  "清空":"Clear",
  "构型（决定单体型 / 组合型）":"Architecture (single monomer / paired)",
  "spacer 长度（两半位点之间）":"Spacer length (between the two half-sites)",
  "默认 5–7":"Default 5–7",
  "自定义":"Custom",
  "6 或 5,6,8 或 5-7":"6, or 5,6,8, or 5-7",
  "指头阵列":"Finger array",
  "三指阵列（读 9 bp）":"3-finger array (reads 9 bp)",
  "六指阵列（读 18 bp，两个三指模块串联）":"6-finger array (reads 18 bp, two 3-finger modules)",
  "自定义指数（Custom）":"Custom finger count",
  "2–12 指；非 3/6 用通用拼接骨架":"2–12 fingers; counts other than 3/6 use the generic scaffold",
  "需要基因组单拷贝结合时再加长：三指单独不唯一，六指/更长阵列能提高唯一性，但模块间界面从未验证过。":"Lengthen the array when single-copy genomic binding is required: a 3-finger ZFP is not unique on its own; 6-finger or longer arrays improve uniqueness, but module interfaces have never been validated.",
  "参数":"Parameters",
  "每靶点数":"Designs per site",
  "最多处理位点":"Max sites",
  "每侧最大错配":"Max mismatches per half-site",
  "Target filter · 全基因组重复过滤（可选，默认关闭）":"Target filter · genome-wide repeat filtering (optional, off by default)",
  "选了基因组就按基因组口径过滤：候选 motif 的落点压在 RepeatMasker 注释区， 或在基因组里出现 ≥3 次，整条剔除、":"With a genome selected, filtering is genome-wide: a candidate motif is dropped entirely if any landing falls inside a RepeatMasker annotation, or if it occurs ≥3 times in the genome,",
  "不参与打分":"not scored",
  "。只留通过过滤的候选做 DesignScore。":". Only survivors go on to DesignScore.",
  "（不过滤：只做 DesignScore）":"(no filtering: DesignScore only)",
  "（不叠 RepeatMasker 注释）":"(no RepeatMasker annotation)",
  "Design Score（S1 + S2）":"Design Score (S1 + S2)",
  "Off-target Score（S3 · 单独执行）":"Off-target Score (S3 · separate run)",
  "对上一次 Design Score 出来的结果，取 Design_score 前 N 个设计做全基因组脱靶扫描 （每侧 0–2 错配）。结果":"Uses the last Design Score run: the top N designs by Design_score are scanned genome-wide (0–2 mismatches per half-site). Results",
  "并列展示，不加进 Design_score":"are shown side by side and never added to Design_score",
  "。每扫一遍全基因组约 1–2 分钟。":". One genome pass takes about 1–2 minutes.",
  "前提：只有通过 Target filter（全基因组重复过滤）的候选，才能做 Off-target Score。 所以请先在「Target filter」里选基因组 FASTA（可选再加 RepeatMasker 注释），跑一次 Design Score， 再点本按钮。没走过过滤的作业会被拒绝，并提示原因。":"Prerequisite: only candidates that passed Target filter (genome-wide repeat filtering) can be scored for off-target effects. Pick a genome FASTA under “Target filter” (RepeatMasker annotation optional), run Design Score, then click this button. Jobs without filtering are rejected with the reason.",
  "每侧错配上限":"Max mismatches per half-site",
  "Off-target Score（S3）":"Off-target Score (S3)",
  "序列只在本机内存里处理，不联网、不上传。":"Sequences are processed in local memory only — no network, no upload.",
  "运行一次，这里会出现统计、设计清单和完整报告。":"Run once; statistics, designs and the full report appear here.",
  "引擎":"Engine",
  "输出文件在临时目录，关掉服务不会自动删除。":"Output files live in a temporary directory; stopping the server does not delete them.",
  "单体型：只扫":"Single monomer: scans only",
  "一个半位点":"one half-site",
  "（无 spacer、无搭档、不产生双链断裂），":"(no spacer, no partner, no double-strand break),",
  "脱靶按「单点位」判定，结果不能与组合型直接比较。":"Off-target is judged per single half-site and cannot be compared with paired results.",
  "组合型：扫":"Paired: scans",
  "左半位点 + spacer + 右半位点":"left half-site + spacer + right half-site",
  "，两个单体靠 FokI 二聚化才切割。":"; the two monomers only cut after FokI dimerization.",
  "完成":"Done",
  "失败":"Failed",
  "运行失败":"Run failed",
  "提交失败":"Submit failed",
  "作业不存在":"Job not found",
  "先跑一次设计。":"Run a Design Score first.",
  "先跑一次 Design Score，再点 Off-target Score。":"Run Design Score first, then Off-target Score.",
  "先粘一条序列，或点「载入示例」。":"Paste a sequence first, or click “Load example”.",
  "至少选一种效应结构。":"Select at least one effector.",
  "生成融合蛋白与 CDS":"Build fusion proteins and CDS",
  "生成设计":"Generate designs",
  "整包下载 (zip)":"Download bundle (zip)",
  "打开文件夹":"Open folder",
  "输出目录":"Output directory",
  "输出链数":"Chains",
  "融合条目":"Fusion entries",
  "链":"Chain",
  "含 NLS":"NLS",
  "无 NLS":"No NLS",
  "无 Kozak":"No Kozak",
  "脱靶 burden":"Off-target burden",
  "burden 排名":"Burden rank",
  "burden 最低":"Min burden",
  "burden 最高":"Max burden",
  "off-target score 最高":"Max off-target score",
  "得分":"Score",
  "位点":"Site",
  "设计":"Design",
  "近似位点（错配:条数）":"Near sites (mismatches:count)",
  "参与设计":"Designs scored",
  "可设计":"Designable",
  "缺口":"Gap",
  "缺":"Missing",
  "有":"Yes",
  "未验证":"Unverified",
  "合格":"OK",
  "CDS 合格":"CDS OK",
  "需复核":"Review",
  "缺元件":"Missing part",
  "跳过":"Skipped",
  "正在做 Design Score…":"Running Design Score…",
  "接回后台作业…":"Reconnecting to the job…",
  "正在拼装效应结构、优化密码子…":"Assembling effectors, optimizing codons…",
  "结果":"Result",
  "refs 加载失败":"failed to load refs",
  "refs 目录：":"refs dir:",
  "（refs 目录里没有基因组 FASTA）":"(no genome FASTA in refs)",
  "（refs 目录里没有 RepeatMasker 注释）":"(no RepeatMasker annotation in refs)",
  "引擎信息读取失败":"failed to read engine info",
  "打不开文件管理器：":"cannot open the file manager:"
};
const I18N_PHASE = {
  "扫描候选位点":"Scanning candidate sites",
  "全基因组重复过滤":"Genome-wide repeat filtering",
  "三联体可设计性评估":"Triplet designability",
  "生成 ZFN 设计":"Building ZFN designs",
  "打分排序":"Scoring and ranking",
  "全基因组脱靶扫描":"Genome-wide off-target scan",
  "准备":"Preparing"
};
const CN_RE = /[\u4e00-\u9fff]/;
let LANG = localStorage.getItem("zfn_lang") === "en" ? "en" : "zh";
let __langBusy = false;

function _tr(s){
  const k = String(s).replace(/\s+/g, " ").trim();
  if (!k) return s;
  if (I18N_EN[k]) return I18N_EN[k];
  let out = k;
  for (const [zh, en] of Object.entries(I18N_PHASE)) out = out.split(zh).join(en);
  out = out.replace(/（已用\s*(\d+)\s*分\s*(\d+)\s*秒）/g, " (elapsed $1 min $2 s)");
  out = out.replace(/（已用\s*(\d+)\s*秒）/g, " (elapsed $1 s)");
  out = out.replace(/^完成，用时\s*([\d.]+)\s*秒$/, "Done in $1 s");
  return out;
}

function applyLang(lang){
  __langBusy = true;
  LANG = lang === "en" ? "en" : "zh";
  document.documentElement.lang = LANG === "en" ? "en" : "zh-CN";
  document.title = LANG === "en" ? "ZFN Design Workbench" : "ZFN 设计工作台";
  const btn = document.getElementById("lang");
  if (btn) btn.textContent = LANG === "en" ? "中文" : "EN";
  const dl = document.getElementById("doclink");
  if (dl) dl.href = "/manual?lang=" + LANG;
  const tw = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT, {
    acceptNode(n){
      const tag = n.parentNode && n.parentNode.nodeName;
      if (tag === "SCRIPT" || tag === "STYLE") return NodeFilter.FILTER_REJECT;
      if (n.parentNode && n.parentNode.id === "lang") return NodeFilter.FILTER_REJECT;
      if (LANG === "en") return CN_RE.test(n.nodeValue || "") ? NodeFilter.FILTER_ACCEPT : NodeFilter.FILTER_SKIP;
      return n.__zh != null ? NodeFilter.FILTER_ACCEPT : NodeFilter.FILTER_SKIP;
    }
  });
  const nodes = [];
  while (tw.nextNode()) nodes.push(tw.currentNode);
  for (const n of nodes){
    if (LANG === "en"){
      if (n.__zh == null) n.__zh = n.nodeValue;
      const raw = n.__zh;
      const m = raw.match(/^(\s*)([\s\S]*?)(\s*)$/);
      if (!m) continue;
      const next = m[1] + _tr(m[2]) + m[3];
      if (n.nodeValue !== next) n.nodeValue = next;
    } else if (n.__zh != null){
      if (n.nodeValue !== n.__zh) n.nodeValue = n.__zh;
    }
  }
  document.querySelectorAll("input,textarea,select,button,option").forEach(el => {
    ["placeholder","title"].forEach(a => {
      const cur = el.getAttribute(a);
      const key = "__zh_" + a;
      if (cur == null){ return; }
      if (el[key] == null) el[key] = cur;
      el.setAttribute(a, LANG === "en" ? _tr(el[key]) : el[key]);
    });
  });
  setTimeout(() => { __langBusy = false; }, 40);
}

(function(){
  const btn = document.getElementById("lang");
  if (btn) btn.onclick = () => {
    clearTimeout(window.__langT);
    const next = LANG === "en" ? "zh" : "en";
    localStorage.setItem("zfn_lang", next);
    applyLang(next);
  };
  const obs = new MutationObserver(() => {
    if (__langBusy || LANG !== "en") return;
    clearTimeout(window.__langT);
    window.__langT = setTimeout(() => { if (LANG === "en") applyLang("en"); }, 60);
  });
  obs.observe(document.body, {childList:true, subtree:true, characterData:true});
  setTimeout(() => applyLang(LANG), 0);
})();

function parseSpacerText(text){
  const out = new Set();
  String(text || "").replace(/\s+/g, "").split(",").forEach(chunk => {
    if (!chunk) return;
    const m = chunk.match(/^(\d+)-(\d+)$/);
    if (m){
      const a = +m[1], b = +m[2];
      for (let v = Math.min(a,b); v <= Math.max(a,b); v++) out.add(v);
    } else if (/^\d+$/.test(chunk)){
      out.add(+chunk);
    }
  });
  return [...out].sort((a,b) => a-b);
}
function spacerMode(){
  const el = document.querySelector('input[name="spmode"]:checked');
  return el ? el.value : "default";
}
function spacers(){
  if (spacerMode() === "default") return [5,6,7];
  const v = parseSpacerText($("#spcustom").value);
  return v.length ? v : [5,6,7];
}
function syncSpacerInput(){
  $("#spcustom").style.display = spacerMode() === "custom" ? "" : "none";
}
document.querySelectorAll('input[name="spmode"]').forEach(el => el.onchange = syncSpacerInput);

// 构型决定单体型还是组合型：单体型没有 spacer，也没有 FokI 二聚化。
const SINGLE_ARCHS = new Set(["single_activator", "single_repressor"]);
function syncArch(){
  const sel = $("#arch");
  if (!sel || !sel.value) return;
  const single = SINGLE_ARCHS.has(sel.value);
  $("#spacerlabel").style.opacity = single ? ".45" : "";
  $("#spacerrow").style.opacity = single ? ".45" : "";
  document.querySelectorAll('input[name="spmode"]').forEach(el => { el.disabled = single; });
  $("#spcustom").disabled = single;
  $("#archnote").innerHTML = single
    ? "单体型：只扫<b>一个半位点</b>（无 spacer、无搭档、不产生双链断裂），" +
      "脱靶按「单点位」判定，结果不能与组合型直接比较。"
    : "组合型：扫<b>左半位点 + spacer + 右半位点</b>，两个单体靠 FokI 二聚化才切割。";
}
function updateSeqInfo(){
  const raw = $("#seq").value;
  const n = raw.replace(/^>.*$/gm, "").replace(/\s+/g, "").length;
  $("#seqinfo").textContent = n ? n + " bp" : "";
}
$("#seq").addEventListener("input", updateSeqInfo);

// 指头阵列：3 / 6 / 自定义（2–12）。自定义时把数字框露出来。
function syncFingers(){
  const sel = $("#fingersel"), box = $("#fingcustom");
  if (!sel) return;
  const custom = sel.value === "custom";
  box.style.display = custom ? "" : "none";
  if (custom && !(box.value >= 2 && box.value <= 12)) box.value = "4";
}
function fingerCount(){
  const sel = $("#fingersel");
  if (sel && sel.value === "custom"){
    let n = parseInt($("#fingcustom").value, 10);
    if (!Number.isFinite(n)) n = 4;
    return Math.min(12, Math.max(2, n));
  }
  return sel ? (parseInt(sel.value, 10) || 3) : 3;
}
$("#fingersel").addEventListener("change", syncFingers);
$("#fingcustom").addEventListener("input", () => {
  const v = parseInt($("#fingcustom").value, 10);
  if (Number.isFinite(v) && v >= 2 && v <= 12) $("#archnote2").textContent =
    `自定义 ${v} 指：半位点 ${v * 3} bp，骨架长度 ${28 * v + 6} aa（通用拼接，模块间界面未验证）`;
});
syncFingers();

// 纯序列没有名字，用 FASTA 头或示例名时才有意义的名字；默认 "input"
let seqName = "input";
async function loadDemo(){
  const r = await fetch("/api/demo").then(r => r.json());
  $("#seq").value = r.seq;
  seqName = r.name || "input";
  updateSeqInfo();
  setStatus("");
}
$("#demo").onclick = loadDemo;
$("#clear").onclick = () => { $("#seq").value = ""; seqName = "input"; updateSeqInfo(); };
$("#run").onclick = run;
$("#ot-run").onclick = runOfftarget;

// refs 下拉：基因组 FASTA / RepeatMasker 注释（S2 全基因组重复过滤）
async function loadRefsOptions(){
  try{
    const r = await fetch("/api/refs").then(x => x.json());
    const fill = (sel, items, keep, label) => {
      const el = $(sel);
      if (!el) return;
      for (const it of items || []){
        const o = document.createElement("option");
        o.value = it.path;
        o.textContent = `${it.name} · ${it.size_h}`;
        el.appendChild(o);
      }
      if (!(items || []).length){
        const o = document.createElement("option");
        o.value = "";
        o.textContent = keep;
        o.disabled = true;
        el.appendChild(o);
      }
    };
    fill("#genome", r.genomes, "（refs 目录里没有基因组 FASTA）");
    fill("#repeats", r.repeats, "（refs 目录里没有 RepeatMasker 注释）");
    if (r.dirs && r.dirs.length) console.log("refs 目录：", r.dirs);
  }catch(e){ console.warn("refs 加载失败", e); }
}
loadRefsOptions();

function setStatus(html){ $("#status").innerHTML = html; }


// ---------------------------------------------------------------- 后台作业
async function readJson(res){
  const text = await res.text();
  try{ return JSON.parse(text); }
  catch(e){
    const head = text.replace(/\s+/g, " ").slice(0, 160);
    throw new Error(`服务端返回的不是 JSON（HTTP ${res.status}）：${head}`);
  }
}

async function postJob(url, body){
  const res = await fetch(url, {method:"POST", headers:{"Content-Type":"application/json"},
                                body: JSON.stringify(body)});
  const r = await readJson(res);
  if (!r.ok) throw new Error(r.error || "提交失败");
  if (!r.job_id) return r;           // 兼容旧的同步返回
  const kind = url.indexOf("genome") >= 0 ? "genome" : "single";
  rememberJob(r.job_id, kind);
  return {job_id: r.job_id};
}

// 轮询到结束；onTick(秒, 作业对象) 用来刷新界面上的计时与进度条
async function waitJob(jid, onTick){
  const t0 = Date.now();
  for (;;){
    await new Promise(r => setTimeout(r, 1500));
    const sec = Math.round((Date.now() - t0) / 1000);
    const res = await fetch(`/api/job?id=${encodeURIComponent(jid)}`);
    const r = await readJson(res);
    if (!r.ok) throw new Error(r.error || "作业不存在");
    if (onTick) onTick(sec, r);
    if (r.state === "running") continue;
    if (r.state === "error") throw new Error(r.error || "运行失败");
    return r.result || {};
  }
}

// 进度条：pct 是 0–100；没给 pct 就按 frac
function setBar(bar, fill, txt, prog, sec){
  prog = prog || {};
  const pct = Math.max(0, Math.min(100, Math.round((prog.frac || 0) * 100)));
  if (fill) fill.style.width = pct + "%";
  if (bar) bar.classList.toggle("on", pct > 0 && pct < 100);
  if (txt){
    const left = prog.phase ? `${esc(prog.phase)} · ${pct}%` : "";
    const detail = prog.message ? ` · ${esc(prog.message)}` : "";
    txt.textContent = `${left}${detail}${sec != null ? `（已用 ${fmtDur(sec)}）` : ""}`;
  }
}

function fmtDur(sec){
  if (sec < 60) return `${sec} 秒`;
  const m = Math.floor(sec / 60), s = sec % 60;
  return s ? `${m} 分 ${s} 秒` : `${m} 分钟`;
}

async function run(){
  const seq = $("#seq").value.trim();
  if (!seq){ $("#result").innerHTML = '<div class="err">先粘一条序列，或点「载入示例」。</div>'; return; }
  $("#run").disabled = true;
  $("#ot-run").disabled = true;
  const t0 = Date.now();
  setStatus('<span class="spin"></span> 正在做 Design Score…');
  setBar($("#run-bar"), $("#run-fill"), $("#run-txt"), {phase:"准备", frac:0}, 0);
  const body = {
    seq,
    name: seqName,
    spacers: spacers(),
    fingers: fingerCount(),
    n_designs: +$("#ndesign").value || 3,
    max_sites: +$("#maxsites").value || 20,
    max_mm_left: +$("#mm").value || 0,
    max_mm_right: +$("#mm").value || 0,
    architecture: $("#arch").value,
    genome: $("#genome").value || "",
    repeats: $("#repeats").value || ""
  };
  try{
    const j = await postJob("/api/run", body);
    const r = j.job_id ? await waitJob(j.job_id, (s, job) => {
      setStatus(`<span class="spin"></span> 正在做 Design Score…`);
      setBar($("#run-bar"), $("#run-fill"), $("#run-txt"), job.progress, s);
    }) : j;
    render(r);
    setBar($("#run-bar"), $("#run-fill"), $("#run-txt"), {phase:"完成", frac:1}, null);
    setStatus(`完成，用时 ${((r.elapsed_ms||0)/1000).toFixed(1)} 秒`);
  }catch(e){
    $("#result").innerHTML = `<div class="err">${esc(e.message)}</div>`;
    setStatus("");
    setBar($("#run-bar"), $("#run-fill"), $("#run-txt"), {phase:"失败", frac:0}, null);
  }finally{
    forgetJob();
    $("#run").disabled = false;
    $("#ot-run").disabled = false;
  }
}

// ---------------------------------------------------------------- Off-target Score（S3 单独按键）
function lastOutDir(){
  const r = window.__lastRun || {};
  return r.out_dir || "";
}

async function runOfftarget(){
  if (!$("#ot-run")) return;
  const dir = lastOutDir();
  if (!dir){
    $("#result").innerHTML = '<div class="err">先跑一次 Design Score，再点 Off-target Score。</div>';
    return;
  }
  $("#ot-run").disabled = true;
  const body = { dir, top: +$("#ot-top").value || 10, max_mm: +$("#ot-mm").value || 2 };
  setBar($("#ot-bar"), $("#ot-fill"), $("#ot-txt"), {phase:"准备", frac:0}, 0);
  try{
    const j = await postJob("/api/offtarget", body);
    const r = j.job_id ? await waitJob(j.job_id, (s, job) => {
      setBar($("#ot-bar"), $("#ot-fill"), $("#ot-txt"), job.progress, s);
    }) : j;
    renderOfftarget(r);
    setBar($("#ot-bar"), $("#ot-fill"), $("#ot-txt"), {phase:"完成", frac:1}, null);
  }catch(e){
    $("#result").innerHTML = `<div class="err">${esc(e.message)}</div>`;
    setBar($("#ot-bar"), $("#ot-fill"), $("#ot-txt"), {phase:"失败", frac:0}, null);
  }finally{
    $("#ot-run").disabled = false;
  }
}

function renderOfftarget(r){
  const s = r.summary || {};
  window.__lastRun = Object.assign({}, window.__lastRun || {}, {out_dir: r.out_dir});
  const rows = (r.designs || []).map(d => `<tr>
      <td class="mono">${esc(d.burden_rank ?? "")}</td>
      <td class="mono">${esc(d.design_id)}</td>
      <td class="mono">${esc(d.seq_name)}:${esc(d.start)}${esc(d.strand)}</td>
      <td class="mono">${d.design_score ?? ""}</td>
      <td class="mono"><b>${d.offtarget_score ?? ""}</b></td>
      <td class="mono">${d.burden ?? ""}${d.truncated ? "（下限）" : ""}</td>
      <td class="mono">${esc(Object.entries(d.by_mm || {}).map(([k,v]) => `${k}:${v}`).join(" "))}</td>
    </tr>`).join("");
  $("#result").innerHTML = `
    <div style="display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap;align-items:baseline">
      <div style="font-weight:600">结果</div>
      <div class="note">S3 · Off-target Score</div>
    </div>
    <div class="cards" style="margin-top:12px">
      <div class="k"><div class="kv">${s.n_designs ?? 0}</div><div class="kl">参与设计</div></div>
      <div class="k"><div class="kv">${s.burden_min ?? ""}</div><div class="kl">burden 最低</div></div>
      <div class="k"><div class="kv">${s.burden_max ?? ""}</div><div class="kl">burden 最高</div></div>
      <div class="k"><div class="kv">${s.score_max ?? ""}</div><div class="kl">off-target score 最高</div></div>
    </div>
    <table class="t" style="margin-top:12px"><thead><tr>
      <th>burden 排名</th><th>设计</th><th>位点</th><th>Design_score</th>
      <th>off-target score</th><th>burden</th><th>近似位点（错配:条数）</th>
    </tr></thead><tbody>${rows}</tbody></table>
    <div class="note" style="margin-top:8px">
      两个分并列，<b>不相加</b>。burden 越小越好；命中数触到上限时 burden 是下限，标「（下限）」。
    </div>
    <div class="files">${fileLinks(r.files)}</div>
    ${s4PanelHtml()}
    <div class="outdir">
      <span class="note">输出目录</span>
      <code title="${esc(r.out_dir || "")}">${esc(shortPath(r.out_dir))}</code>
      <button class="mini" onclick="revealFolder()">打开文件夹</button>
    </div>
  `;
  wireS4();
}

function shortPath(p){
  p = p || "";
  return p.length > 54 ? "…" + p.slice(-54) : p;
}

// ---------------------------------------------------------------- 刷新后接回
function rememberJob(jid, kind){
  try{ localStorage.setItem("zfn_job", JSON.stringify({id: jid, kind, t: Date.now()})); }catch(e){}
}
function forgetJob(){
  try{ localStorage.removeItem("zfn_job"); }catch(e){}
}

async function resumeJob(){
  let saved = null;
  try{ saved = JSON.parse(localStorage.getItem("zfn_job") || "null"); }catch(e){}
  if (!saved || !saved.id) return;
  try{
    const r = await readJson(await fetch(`/api/job?id=${encodeURIComponent(saved.id)}`));
    if (!r.ok){ forgetJob(); return; }
    if (r.state === "running"){
      $("#run").disabled = true;
      setStatus('<span class="spin"></span> 接回后台作业…');
      const res = await waitJob(saved.id, (s) =>
        setStatus(`<span class="spin"></span> 正在设计… 已用 ${fmtDur(s)}`));
      render(res);
      setStatus(`完成，用时 ${((res.elapsed_ms||0)/1000).toFixed(1)} 秒`);
      $("#run").disabled = false;
      forgetJob();
    } else if (r.state === "done" && r.result){
      render(r.result);
      setStatus("已接回上次完成的结果");
      forgetJob();
    } else {
      forgetJob();
    }
  }catch(e){ forgetJob(); }
}

setTimeout(resumeJob, 700);

async function revealFolder(){
  try{
    const r = await fetch("/api/reveal", {method:"POST"}).then(x => x.json());
    if (!r.ok){ alert("打不开文件管理器：" + (r.error || "") + "\n\n目录：" + r.dir); }
  }catch(e){ alert("打不开文件管理器：" + e.message); }
}

const FILE_LABELS = {json:"report.json", html:"report.html", sites_tsv:"sites.tsv",
                     designs_tsv:"designs.tsv", fusions_json:"fusions.json",
                     offtarget_tsv:"offtarget_genome.tsv", proteins_tsv:"proteins.tsv",
                     cds_tsv:"cds.tsv", cds_fasta:"cds.fasta", order_form:"order_form.csv",
                     report_html:"report.html", report_json:"report.json",
                     bundle:"整包下载 (zip)"};

function fileLinks(files){
  return Object.entries(files || {}).map(([k,v]) =>
    `<a href="${v}" target="_blank" rel="noopener"${k === "bundle" ? ' class="zip"' : ""}>` +
    `${esc(FILE_LABELS[k] || k)}</a>`).join("");
}

function s4PanelHtml(){
  return `
    <div class="s4">
      <div class="s4h">④ S4 · Fusion design（单独执行，不重跑设计）</div>
      <div class="note">从这次结果里挑设计，接上效应结构，生成蛋白序列与人源化 CDS，并出一张合成订单表。</div>
      <div class="note">S3 跑完也可以直接点这里：会带上该设计的 off-target 分（burden）一起进订单表。</div>
      <div class="s4row">
        <label class="note">取前 <input id="fu-top" type="number" value="3" min="1" max="20" style="width:64px"> 个设计</label>
        <label class="note"><input type="checkbox" id="fu-act" checked> 激活（VP64）</label>
        <label class="note"><input type="checkbox" id="fu-rep" checked> 抑制（KRAB）</label>
        <label class="note"><input type="checkbox" id="fu-nuc"> 核酸酶（FokI）</label>
        <button id="fu-run">生成融合蛋白与 CDS</button>
      </div>
      <details style="margin-top:8px"><summary class="note">需要换元件就填序列（JSON：{"KRAB":"…","FokI_CD":"…"}）</summary>
        <textarea id="fu-parts" rows="3" placeholder='{"KRAB":"你的KRAB片段","FokI_CD":"你的FokI催化域"}'></textarea>
      </details>
      <div id="fusion-out"></div>
    </div>`;
}

function wireS4(){
  const b = $("#fu-run");
  if (b) b.onclick = runFusion;
}

function render(r){
  window.__lastRun = r;
  const b = r.by_status || {};
  const lib = r.library || {};
  const cards = [
    ["候选位点", r.n_sites, ""],
    ["可设计·Tier1", b.full ?? 0, "ok"],
    ["含未验证", b.partial ?? 0, "warn"],
    ["缺口", b.gap ?? 0, "bad"],
    ["生成设计", r.designs.length, ""],
    ["库覆盖", `${lib.triplets_covered ?? 0}/${lib.triplets_total ?? 64}`, ""]
  ].map(([k,v,cls]) => `<div class="stat"><div class="k" title="${k}">${k}</div><div class="v ${cls}">${v}</div></div>`).join("");

  // 横幅已经讲了阵列长度的限制，别在警告里重复一遍
  const warn = (r.warnings||[])
    .filter(w => !/只读 9 bp|两个三指模块|当前是\*\*单体型\*\*/.test(w))
    .map(w => `<div class="note" style="color:#e7cf9f">⚠ ${esc(w)}</div>`).join("");

  const f = r.fingers || 3, hs = r.half_site_bp || (f * 3);
  const single = (r.params||{}).mode === "single";
  const limitBanner = single
    ? `<div class="limit"><b>单体型</b>（当前构型）：只扫<b>一个 ${hs} bp 半位点</b>，没有 spacer、没有搭档、` +
      `没有 FokI 二聚化，因此<b>不产生双链断裂</b>——它是 ZF-TF（激活/抑制），不是剪刀。` +
      `脱靶判据是「一个半位点像不像」，命中数天然多于组合型，<b>两类分数不可直接比较</b>。` +
      (f <= 3
        ? `另外：三指 ZFP 只读取 <b>9 bp</b>，单独作为单体在人类基因组中<b>并不唯一</b>；`
        : ``) +
      `本文件不保证活性或特异性，所有设计均需实验验证。</div>`
    : (f <= 3
    ? `<div class="limit"><b>最重要的限制</b>：三指 ZFP 只读取 <b>9 bp</b>，` +
      `单独作为单体在人类基因组中<b>并不唯一</b>。若需要基因组单拷贝结合，` +
      `应把它作为 ZFN 配对的一半（本工具默认成对输出），或加长阵列（六指读 18 bp）并重新验证。` +
      `本文件不保证活性或特异性，所有设计均需实验验证。</div>`
    : (() => {
        // 模块 N 端指：1, 4, 7…
        const heads = [];
        for (let i = 1; i <= f; i += 3) heads.push("F" + i);
        return `<div class="limit">当前是 <b>${f} 指阵列</b>，每个单体读 <b>${hs} bp</b>。` +
          `阵列按三指一模块首尾相接：每个模块的骨架逐字符取自 1AAY，` +
          `但<b>模块之间的 TGEKP 连接、相对取向与界面没有任何结构或活性证据</b>；` +
          `且模块 N 端指（${heads.join("/")}）在 B1H 库里没有本位上下文数据，只能借用别的位置的实测指头。` +
          `本文件不保证活性或特异性，所有设计均需实验验证。</div>`;
      })());

  const designCards = r.designs.map((d, i) => {
    const cls = "g" + (d.grade || "D");
    const comp = (d.components||[]).map(c =>
      `<tr><td>${esc(c.label_zh)}</td><td class="mono">${c.weight.toFixed(2)}</td>
       <td class="mono">${c.raw.toFixed(1)}</td><td class="mono">${c.weighted.toFixed(2)}</td>
       <td class="note">${esc(c.detail)}</td></tr>`).join("");
    const paired = d.is_paired !== false;
    const trip = paired
      ? `三联体 ${esc((d.triplets_left||[]).join("&gt;"))} / ${esc((d.triplets_right||[]).join("&gt;"))}`
      : `三联体 ${esc((d.triplets_left||[]).join("&gt;"))}（单阵列，无右单体）`;
    const rightSeq = paired
      ? `<details><summary>右单体蛋白序列（${(d.partner_protein||"").length} aa）</summary>
         <div class="seq">${esc(d.partner_protein)}</div></details>`
      : "";
    return `<div class="design">
      <div class="h">
        <div>
          <div class="id">#${d.rank ?? (i + 1)} ${esc(d.id)}</div>
          <div class="note">${esc(d.target?.site_seq || "")} ·
            ${paired ? "组合型（成对）" : "单体型（单阵列）"} · ${trip}</div>
        </div>
        <div><span class="sc">${(d.total ?? 0).toFixed(1)}</span><span class="grade ${cls}">${esc(d.grade)}</span></div>
      </div>
      <table style="margin-top:8px">
        <thead><tr><th>分项</th><th>权重</th><th>得分</th><th>加权</th><th>依据</th></tr></thead>
        <tbody>${comp}</tbody>
      </table>
      <div class="note" style="margin-top:6px">
        ${paired ? "脱靶 burden" : "单点位脱靶 burden"}=${d.burden ?? "—"}，命中 ${d.n_hits ?? 0} 个（exact ${d.n_exact ?? 0}）· 证据层级 Tier${d.min_evidence_tier}
      </div>
      <details><summary>${paired ? "左单体" : "阵列"}蛋白序列（${(d.protein||"").length} aa）</summary>
        <div class="seq">${esc(d.protein)}</div></details>
      ${rightSeq}
    </div>`;
  }).join("");

  const files = fileLinks(r.files);

  const demand = (r.triplet_demand||[]).map(t =>
    `<tr><td class="mono">${esc(t.triplet)}</td><td class="mono">${t.occurrences}</td>
     <td>${t.covered ? '<span class="pill full">有</span>' : '<span class="pill gap">缺</span>'}</td>
     <td class="mono">${t.n_candidate_slots_gap ?? 0}</td></tr>`).join("");

  const siteRows = (r.sites||[]).slice(0, 25).map(s =>
    `<tr><td class="mono">${esc(s.site?.site_id || s.site_id || "")}</td>
     <td class="mono">${esc(s.site?.left_half || "")}</td>
     <td class="mono">${esc(s.site?.spacer || "")}</td>
     <td class="mono">${esc(s.site?.right_half || "")}</td>
     <td>${s.status === "full" ? '<span class="pill full">可设计</span>'
            : s.status === "partial" ? '<span class="pill warnp">未验证</span>'
            : '<span class="pill gap">缺口</span>'}</td>
     <td class="mono">${s.n_combos ?? 0}</td></tr>`).join("");

  $("#result").innerHTML = `
    <div style="display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap;align-items:baseline">
      <div style="font-weight:600">结果</div>
      <div class="note">run ${esc(r.run_id)} · ${esc(r.generated_at)}</div>
    </div>
    <div class="stats" style="margin-top:12px">${cards}</div>
    ${limitBanner}
    ${warn}
    <div class="files">${files}</div>
    ${s4PanelHtml()}
    <div class="outdir">
      <span class="note">输出目录</span>
      <code title="${esc(r.out_dir || "")}">${esc(shortPath(r.out_dir))}</code>
      <button class="mini" onclick="revealFolder()">打开文件夹</button>
    </div>
    <div style="font-weight:600;margin-top:18px">设计（按分数排序）</div>
    ${designCards || '<div class="note">没有可输出的设计——位点全部落在指头库缺口上。</div>'}
    ${demand ? `<div style="font-weight:600;margin-top:18px">这段输入最缺的三联体</div>
      <table style="margin-top:6px"><thead><tr><th>三联体</th><th>出现次数</th><th>当前库</th><th>缺口次数</th></tr></thead>
      <tbody>${demand}</tbody></table>` : ""}
    ${siteRows ? `<div style="font-weight:600;margin-top:18px">候选位点（前 25 条）</div>
      <table style="margin-top:6px"><thead><tr><th>位点</th><th>左半位点</th><th>spacer</th><th>右半位点</th><th>可设计性</th><th>组合数</th></tr></thead>
      <tbody>${siteRows}</tbody></table>
      <div class="note" style="margin-top:6px">共 ${r.n_sites} 条，完整清单见 sites.tsv</div>` : ""}
  `;
  wireS4();
}


// ---------------------------------------------------------------- S4 融合设计
function fusionBody(){
  const archs = [];
  if ($("#fu-act") && $("#fu-act").checked) archs.push("single_activator");
  if ($("#fu-rep") && $("#fu-rep").checked) archs.push("single_repressor");
  if ($("#fu-nuc") && $("#fu-nuc").checked) archs.push("paired_nuclease");
  const r0 = window.__lastRun || {};
  let parts = null;
  const pt = $("#fu-parts");
  if (pt && pt.value.trim()){
    try{ parts = JSON.parse(pt.value); }catch(e){ throw new Error("元件 JSON 解析失败：" + e.message); }
  }
  return {dir: r0.out_dir, top: +($("#fu-top") ? $("#fu-top").value : 3) || 3, archs, parts};
}

async function runFusion(){
  const r0 = window.__lastRun || {};
  if (!r0.out_dir){ $("#fusion-out").innerHTML = '<div class="err">先跑一次设计。</div>'; return; }
  let body;
  try{ body = fusionBody(); }catch(e){ $("#fusion-out").innerHTML = `<div class="err">${esc(e.message)}</div>`; return; }
  if (!body.archs.length){ $("#fusion-out").innerHTML = '<div class="err">至少选一种效应结构。</div>'; return; }
  const btn = $("#fu-run");
  btn.disabled = true;
  $("#fusion-out").innerHTML = '<div class="note">正在拼装效应结构、优化密码子…</div>';
  try{
    const j = await postJob("/api/fusion", body);
    const res = j.job_id ? await waitJob(j.job_id) : j;
    renderFusion(res);
  }catch(e){
    $("#fusion-out").innerHTML = `<div class="err">${esc(e.message)}</div>`;
  }finally{
    btn.disabled = false;
  }
}

function renderFusion(res){
  if (!res || res.ok === false || res.state === "error"){
    $("#fusion-out").innerHTML = `<div class="err">${esc((res && (res.error || res.message)) || "失败")}</div>`;
    return;
  }
  const s = res.summary || {};
  const cards = [
    ["融合条目", s.n_entries ?? 0, ""],
    ["输出链数", s.n_chains ?? 0, ""],
    ["CDS 可用", s.n_ok ?? 0, "ok"],
    ["需复核", s.n_warn ?? 0, "warn"],
    ["缺元件", s.n_blocked ?? 0, "bad"],
    ["CDS 总长", (s.cds_total_bp ?? 0) + " bp", ""],
  ].map(([k,v,cls]) => `<div class="stat"><div class="k">${k}</div><div class="v ${cls}">${v}</div></div>`).join("");

  const FILES = {fusions_json:"fusions.json", proteins_tsv:"proteins.tsv", cds_tsv:"cds.tsv",
                 cds_fasta:"cds.fasta", order_form:"order_form.csv", bundle:"整包下载 (zip)"};
  const files = Object.entries(res.files || {}).map(([k,v]) =>
    `<a href="${v}" target="_blank" rel="noopener"${k === "bundle" ? ' class="zip"' : ""}>${esc(FILES[k] || k)}</a>`).join("");

  const entries = (res.entries || []).map(e => {
    if (e.error){
      return `<div class="fu"><div class="h"><div class="id">${esc(e.design_id)} · ${esc(e.architecture)}</div>
        <span class="pill2 bd2">跳过</span></div>
        <div class="note">${esc(e.error)}</div></div>`;
    }
    const st = e.cds_status;
    const pill = st === "ok" ? '<span class="pill2 ok2">CDS 合格</span>'
      : st === "warn" ? '<span class="pill2 wn2">需复核</span>'
      : '<span class="pill2 bd2">缺元件</span>';
    const chains = (e.chains && e.chains.length) ? e.chains : [e];
    const chainHtml = chains.map(ch => {
      const cst = ch.cds_status || "";
      const cpill = cst === "ok" ? '<span class="pill2 ok2">合格</span>'
        : cst === "warn" ? '<span class="pill2 wn2">需复核</span>'
        : cst === "skipped" ? '<span class="pill2 bd2">跳过</span>'
        : '<span class="pill2 bd2">缺元件</span>';
      const chNotes = [...(ch.red_flags || []), ...(ch.cds_notes || [])]
        .map(x => `<div class="note" style="color:#e7cf9f">⚠ ${esc(x)}</div>`).join("");
      const chSites = (ch.cds_forbidden_sites || []).length
        ? `<div class="note" style="color:#ffb0b0">酶切位点：${esc(ch.cds_forbidden_sites.map(x => x.enzyme).join(", "))}</div>` : "";
      const kozakTxt = ch.cds_kozak ? `Kozak ${esc(ch.cds_kozak)}` : "无 Kozak";
      const proteinBlock = ch.protein ? `<details><summary>融合蛋白（${ch.protein_len} aa）</summary><div class="seq">${esc(ch.protein)}</div></details>` : "";
      const cdsBlock = ch.cds ? `<details><summary>CDS（${ch.cds_len} bp，GC ${((ch.cds_gc ?? 0)*100).toFixed(1)}%，${kozakTxt}，起始 ${esc(ch.cds_start_codon || "—")} / 终止 ${esc(ch.cds_stop_codon || "—")}）</summary>
        <div class="seq">${esc(ch.cds)}</div></details>` : "";
      return `<div style="margin-top:8px;padding-top:8px;border-top:1px dashed #24344a">
        <div class="h"><div class="id">${esc(ch.label || ch.side || "链")}${
          ch.fok_part ? ` <span class="note">· ${esc(ch.fok_part)}</span>` : ""}
          <span class="note">· ${ch.protein_len ?? "—"} aa</span></div><div>${cpill}</div></div>
        <div class="note">起始 ${esc(ch.cds_start_codon || "—")} · 终止 ${esc(ch.cds_stop_codon || "—")} · ${
          ch.has_nls ? "含 NLS" : "无 NLS"}</div>
        ${chNotes}${chSites}${proteinBlock}${cdsBlock}
      </div>`;
    }).join("");
    const headLine = `${e.n_chains > 1 ? `成对输出 ${e.n_chains} 条链 · ` : ""}起始 ${esc(e.cds_start_codon || "—")} · 终止 ${esc(e.cds_stop_codon || "—")} · ${
      e.has_nls ? "含 NLS" : "无 NLS"}`;
    return `<div class="fu">
      <div class="h">
        <div class="id">${esc(e.fusion_id || (e.design_id + "__" + e.architecture))}
          <span class="note">· ${esc(e.architecture_label || e.architecture)} · ${e.protein_len} aa</span></div>
        <div>${pill}</div>
      </div>
      <div class="note">${headLine}</div>
      <div class="note">设计 #${esc(String(e.rank ?? "-"))} ${esc(e.design_id)} · Design_score ${e.design_score ?? "—"}${
        e.offtarget_score != null ? ` · off-target ${e.offtarget_score}（burden ${e.offtarget_burden ?? "—"}）` : " · off-target 未算"}${
        e.requires_pair ? " · 需成对位点" : ""}</div>
      ${chainHtml}
    </div>`;
  }).join("");

  $("#fusion-out").innerHTML = `
    <div class="stats" style="margin-top:12px">${cards}</div>
    <div class="files" style="margin-top:10px">${files}</div>
    <div class="note" style="margin-top:8px">S4 只读这次作业的 report.json，不重算 S1–S3；
      蛋白里的 [元件:状态] 方括号 = 缺该元件的可靠序列，那种条目不会出 CDS。</div>
    ${entries}`;
  const b = $("#fu-run");
  if (b) b.onclick = runFusion;
}

fetch("/api/info").then(r => r.json()).then(async info => {
  $("#ver").textContent = info.version;
  $("#env").innerHTML = `引擎 <span class="mono">zfn_design ${esc(info.version)}</span> ·
    指头库 <span class="mono">${info.library.n_records}</span> 条，
    覆盖 <span class="mono">${info.library.triplets_covered}/${info.library.triplets_total}</span> 个三联体 ·
    打分权重 <span class="mono">8 分项</span>`;
  const sel = $("#arch");
  sel.innerHTML = (info.architectures||[]).map(a => `<option value="${esc(a.key)}">${esc(a.label)}</option>`).join("");
  if ([...sel.options].some(o => o.value === "paired_nuclease")) sel.value = "paired_nuclease";
  sel.onchange = syncArch;
  // URL 参数：?fingers=6&spacers=6 或 5-7&arch=single_activator（也用于截图与自检）
  const qs = new URLSearchParams(location.search);
  const qf = qs.get("fingers");
  if (qf === "3" || qf === "6") $("#fingersel").value = qf;
  else if (qf && +qf >= 2 && +qf <= 12) { $("#fingersel").value = "custom"; $("#fingcustom").value = qf; }
  syncFingers();
  if (qs.get("arch") && [...sel.options].some(o => o.value === qs.get("arch"))){
    sel.value = qs.get("arch");
  }
  if (qs.get("spacers")){
    document.querySelector('input[name="spmode"][value="custom"]').checked = true;
    $("#spcustom").value = qs.get("spacers");
    syncSpacerInput();
  }
  syncArch();
  // ?autorun：打开就把示例跑一遍（截图、演示、自检时用）
  if (qs.has("autorun")){ await loadDemo(); run(); }
}).catch(() => { $("#env").textContent = "引擎信息读取失败"; });
</script>
</body>
</html>
"""


# --------------------------------------------------------------------------
# 启动
# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    global OUT_ROOT
    ap = argparse.ArgumentParser(description="ZFN 设计的本地图形界面")
    ap.add_argument("--host", default="127.0.0.1", help="只对本机开放，默认 127.0.0.1")
    ap.add_argument("--port", type=int, default=0, help="0 = 自动挑一个空闲端口")
    ap.add_argument("--no-browser", action="store_true", help="不自动打开浏览器")
    ap.add_argument("--out-root", default=None,
                    help="结果写到哪儿，默认 <项目>/out/gui/")
    args = ap.parse_args(argv)

    if args.out_root:
        OUT_ROOT = Path(args.out_root).expanduser()
        OUT_ROOT.mkdir(parents=True, exist_ok=True)

    httpd = ThreadingHTTPServer((args.host, args.port), _Handler)
    httpd.daemon_threads = True
    port = httpd.server_address[1]
    url = f"http://{args.host}:{port}/"

    print("=" * 62)
    print(f"ZFN 设计工作台  v{__version__}")
    print(f"地址: {url}")
    print(f"结果目录: {OUT_ROOT}")
    print("按 Ctrl+C 停止（或关掉这个窗口）")
    print("=" * 62)

    if not args.no_browser:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()

    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n已停止。")
    finally:
        httpd.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
