"""把整条流水线编译成一个**单文件、零依赖**的可交互 HTML。

为什么要单独做一页
------------------
`out/report.html` 是"跑完一次的结果快照"，只能看不能动；这一页是
"流水线本身的说明书 + 试算器"：

* 左边讲清 7 个阶段各自吃什么、吐什么、会在哪里失败；
* 右边可以贴一段自己的 DNA，在**浏览器里**跑出同一套结果（不用装 Python）。

口径必须和 Python 完全一致，否则这页就成了骗人的演示。所以：
1. 指头库、骨架模板、打分权重、融合元件全部从这里读 data/*.json，
   不在这里复制一份常量；
2. JS 引擎逐行对照 `targets.py` / `designability.py` / `build.py` /
   `score.py` / `offtarget.py` 实现；
3. 页面带 `?selftest=1` 自检模式，把关键数字写进 DOM，供外部脚本与
   `report.json` 比对（`tools/check_pipeline_html.py` 用的就是它）。

用法::

    python tools/build_pipeline_html.py            # 写 docs/pipeline.html
    python tools/build_pipeline_html.py --out x.html
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "zfn_design" / "data"


def load_library() -> list[dict]:
    """和 Python 侧同一个库：有 B1H 全库就用它，否则退回 3 条种子。"""
    path = DATA / "finger_library_b1h.csv"
    if not path.exists():
        path = DATA / "finger_library_seed.csv"
    rows = []
    with path.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            rows.append({
                "id": r["helix_id"],
                "triplet": r["triplet"],
                "helix7": r["helix7"],
                "tier": {"measured": 1, "literature": 1,
                         "variant_unvalidated": 2, "code_inference": 2}.get(
                             (r["evidence_level"] or "").strip().lower(), 3),
                "split": (r["split_position"] or "any").strip(),
                "source": r["source_id"],
                "assay": r["assay"],
                "affinity": float(r["affinity_nM"]) if r["affinity_nM"] else None,
            })
    return rows


def load_demo() -> tuple[str, str]:
    path = ROOT / "examples" / "demo_target.fa"
    name, chunks = "demo_target", []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            head = line[1:].strip()
            name = head.split()[0] if head else "demo_target"
        elif line.strip():
            chunks.append(line.strip().upper())
    return name, "".join(chunks)


def build_payload() -> dict:
    templates = json.loads((DATA / "framework_templates.json").read_text(encoding="utf-8"))
    weights = json.loads((DATA / "scoring_weights.json").read_text(encoding="utf-8"))
    parts = json.loads((DATA / "fusion_parts.json").read_text(encoding="utf-8"))
    demo_name, demo_seq = load_demo()
    return {
        "library": load_library(),
        "templates": templates,
        "weights": weights,
        "fusion": parts,
        "demo": {"name": demo_name, "seq": demo_seq},
    }


TEMPLATE_PATH = Path(__file__).resolve().parent / "pipeline_template.html"


def build(out: Path) -> Path:
    payload = build_payload()
    html = TEMPLATE_PATH.read_text(encoding="utf-8")
    html = html.replace("/*__ZFN_DATA__*/null", json.dumps(payload, ensure_ascii=False))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="生成可交互的流水线 HTML")
    ap.add_argument("--out", default=str(ROOT / "docs" / "pipeline.html"))
    args = ap.parse_args(argv)
    path = build(Path(args.out))
    print(f"已生成 {path}  ({path.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
