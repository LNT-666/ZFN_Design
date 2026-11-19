"""校验 docs/pipeline.html 与 Python 引擎的口径是否完全一致。

页面自己声称"和命令行同一套算法"，这个脚本负责把这句话变成可执行的证据：

1. 用 Python 跑一遍示例序列（直接调 ``zfn_design.pipeline.run``，不开进程）；
2. 用无头 Edge 打开 ``docs/pipeline.html?selftest=1``，把页面算出的关键数字抓回来；
3. 逐项比对：候选位点数、三态分布、设计数、每个设计的 ID 与总分、库覆盖。

任何一项不一致就非零退出。改规则（比如调打分权重）后忘了重建页面时，
这个脚本会先报错——这正是它存在的意义。

用法::

    python tools/check_pipeline_html.py
    python tools/check_pipeline_html.py --html docs/pipeline.html --verbose
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

EDGE_CANDIDATES = [
    Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"),
    Path(r"C:\Program Files\Microsoft\Edge\Application\msedge.exe"),
    Path("/usr/bin/microsoft-edge"),
    Path("/usr/bin/google-chrome"),
]


def find_browser() -> Path | None:
    for p in EDGE_CANDIDATES:
        if p.exists():
            return p
    for name in ("msedge", "microsoft-edge", "google-chrome", "chromium"):
        found = shutil.which(name)
        if found:
            return Path(found)
    return None


def python_side() -> dict:
    from zfn_design.pipeline import run
    from zfn_design.ingest import parse_fasta
    from zfn_design.models import SequenceRecord

    text = (ROOT / "examples" / "demo_target.fa").read_text(encoding="utf-8")
    records = parse_fasta(text, source="examples/demo_target.fa")
    payload = run(records, spacers=(5, 6, 7), n_designs=3, max_sites=20)
    return {
        "n_sites": len(payload["sites"]),
        "by_status": payload["targets_summary"]["by_status"],
        "n_designs": len(payload["designs"]),
        "designs": [
            {"id": d["design_id"], "total": round(d["score"]["total"], 6), "grade": d["score"]["grade"]}
            for d in payload["designs"]
        ],
        "triplets_covered": payload["library"]["triplets_covered"],
    }


def browser_side(html: Path, timeout_ms: int = 12000) -> dict:
    browser = find_browser()
    if browser is None:
        raise RuntimeError("找不到 Edge/Chrome，无法运行页面自检（可手动打开页面加 ?selftest=1）")
    url = html.resolve().as_uri() + "?selftest=1"
    cmd = [
        str(browser), "--headless=new", "--disable-gpu", "--no-sandbox",
        "--user-data-dir=" + str(Path.home() / ".cache" / "zfn-selftest-profile"),
        f"--virtual-time-budget={timeout_ms}", "--dump-dom", url,
    ]
    out = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                         errors="replace", timeout=90)
    dom = out.stdout or ""
    m = re.search(r"SELFTEST:(\{.*?\})</pre>", dom, re.S)
    if not m:
        raise RuntimeError(f"页面没有输出 selftest（DOM {len(dom)} 字符）。"
                           "可能 JS 报错——先手动打开页面看控制台。")
    return json.loads(m.group(1))


def compare(py: dict, js: dict) -> list[str]:
    diffs: list[str] = []
    for key in ("n_sites", "n_designs", "triplets_covered"):
        if py[key] != js[key]:
            diffs.append(f"{key}: python={py[key]}  html={js[key]}")
    for k in ("full", "partial", "gap"):
        a, b = py["by_status"][k], js["by_status"][k]
        if a != b:
            diffs.append(f"by_status.{k}: python={a}  html={b}")
    if [d["id"] for d in py["designs"]] != [d["id"] for d in js["designs"]]:
        diffs.append(f"design ids 不一致:\n  python={[d['id'] for d in py['designs']]}\n  html={[d['id'] for d in js['designs']]}")
    else:
        for a, b in zip(py["designs"], js["designs"]):
            if abs(a["total"] - b["total"]) > 1e-6:
                diffs.append(f"{a['id']} 总分: python={a['total']}  html={b['total']}")
            if a["grade"] != b["grade"]:
                diffs.append(f"{a['id']} 等级: python={a['grade']}  html={b['grade']}")
    # 排版体检：长序列块不该横向溢出（溢出说明 CSS 破了，文字会被切掉）
    overflow = js.get("seq_overflow") or []
    if overflow:
        diffs.append(f"排版: {len(overflow)} 个长序列块横向溢出 {overflow}")
    return diffs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="比对页面引擎与 Python 引擎")
    ap.add_argument("--html", default=str(ROOT / "docs" / "pipeline.html"))
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)

    html = Path(args.html)
    if not html.exists():
        print(f"× 找不到 {html}；先跑 python tools/build_pipeline_html.py")
        return 2

    py = python_side()
    if args.verbose:
        print("python:", json.dumps(py, ensure_ascii=False))
    try:
        js = browser_side(html)
    except Exception as exc:  # noqa: BLE001
        print(f"× 页面自检失败: {exc}")
        return 3
    if args.verbose:
        print("html  :", json.dumps(js, ensure_ascii=False))

    diffs = compare(py, js)
    if diffs:
        print("× 口径不一致：")
        for d in diffs:
            print("   -", d)
        return 1

    print(f"√ 口径一致：{py['n_sites']} 个位点（{py['by_status']['full']}/"
          f"{py['by_status']['partial']}/{py['by_status']['gap']}），"
          f"{py['n_designs']} 个设计，最高分 "
          f"{py['designs'][0]['total'] if py['designs'] else '—'} "
          f"({py['designs'][0]['grade'] if py['designs'] else '—'})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
