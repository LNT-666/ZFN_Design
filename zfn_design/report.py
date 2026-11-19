"""阶段 6：report——把结果写成机器能读、人能看的三份东西。

* ``report.json``  全量、结构化，给下游脚本用；
* ``*.tsv``        表格，直接拖进 Excel；
* ``report.html``  单文件、零外部依赖，浏览器打开即可，给实验同事看。

HTML 刻意做成"自解释"的：每个分数旁边都有分项拆解和算法说明，
任何人不看代码也能判断这个排序是怎么来的。
"""

from __future__ import annotations

import html
import json
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence

from .models import ScoreBreakdown  # noqa: F401  (类型提示用)
from .sequtil import wrap

STATUS_LABEL = {"full": "可设计（全 Tier1）", "partial": "可设计（含未验证）", "gap": "缺口"}
STATUS_CLASS = {"full": "ok", "partial": "warn", "gap": "bad"}
GRADE_CLASS = {"A": "ok", "B": "warn", "C": "warn", "D": "bad"}


# ---------------------------------------------------------------- JSON / TSV
def write_json(path: str | Path, payload: dict) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return p


def _tsv_escape(v: Any) -> str:
    if v is None:
        return ""
    s = str(v)
    return s.replace("\t", " ").replace("\n", " ")


def write_tsv(path: str | Path, header: Sequence[str], rows: Iterable[Sequence[Any]]) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(header)]
    lines += ["\t".join(_tsv_escape(c) for c in row) for row in rows]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8-sig")
    return p


def sites_tsv(path: str | Path, sites: Sequence[dict]) -> Path:
    header = ["site_id", "seq_name", "start", "end", "strand", "left_half", "spacer",
              "spacer_len", "right_half", "designability", "n_combos", "missing_triplets"]
    rows = [
        [s["site_id"], s["seq_name"], s["start"], s["end"], s["strand"], s["left_half"],
         s["spacer"], s["spacer_len"], s["right_half"], s["designability"], s["n_combos"],
         ",".join(s["missing_triplets"])]
        for s in sites
    ]
    return write_tsv(path, header, rows)


def designs_tsv(path: str | Path, designs: Sequence[dict]) -> Path:
    header = ["design_id", "seq_name", "start", "strand", "architecture",
              "left_half", "spacer", "right_half",
              "triplets_N_to_C", "partner_triplets_N_to_C", "min_tier",
              "design_score", "grade", "offtarget_score", "offtarget_burden",
              "protein", "partner_protein"]
    rows = []
    for d in designs:
        sc = d.get("score") or {}
        ot = (d.get("offtarget") or {}).get("summary") or {}
        right = d.get("monomer_right") or {}
        rows.append([
            d["design_id"], d["target"]["seq_name"], d["target"]["start"], d["target"]["strand"],
            d["architecture"], d["target"]["left_half"], d["target"]["spacer"],
            d["target"]["right_half"],
            ",".join(d["monomer_left"]["triplets_N_to_C"]),
            ",".join(right.get("triplets_N_to_C") or []),
            d["min_evidence_tier"], sc.get("total"), sc.get("grade"),
            (d.get("offtarget_score") or {}).get("score"), ot.get("burden"),
            d["monomer_left"]["protein"], right.get("protein", ""),
        ])
    return write_tsv(path, header, rows)


# ---------------------------------------------------------------- HTML
_CSS = """
:root{--bg:#0f1115;--panel:#171a21;--panel2:#1e222b;--fg:#e6e8ee;--dim:#9aa3b2;
--ok:#3fb950;--warn:#d29922;--bad:#f85149;--acc:#58a6ff;--line:#2a2f3a}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
font:14px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI","Microsoft YaHei",sans-serif}
header{padding:28px 32px;border-bottom:1px solid var(--line);background:linear-gradient(180deg,#1b2130,#0f1115)}
h1{margin:0 0 6px;font-size:22px;letter-spacing:.5px}
h2{margin:34px 0 12px;font-size:17px;border-left:3px solid var(--acc);padding-left:10px}
h3{margin:18px 0 8px;font-size:14px;color:var(--dim);text-transform:uppercase;letter-spacing:1px}
main{padding:0 32px 80px;max-width:1400px}
.sub{color:var(--dim);font-size:13px}
.banner{margin:18px 0;padding:14px 16px;border-radius:8px;background:#2a1f14;border:1px solid #5c3d1a;color:#f0c674}
.banner b{color:#ffd479}
.cards{display:flex;gap:14px;flex-wrap:wrap;margin:16px 0}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;padding:14px 18px;flex:1 1 120px;min-width:120px}
.card .k{color:var(--dim);font-size:12px}
.card .v{font-size:24px;font-weight:600;margin-top:4px}
table{border-collapse:collapse;width:100%;margin:10px 0 24px;background:var(--panel);font-size:13px}
th,td{padding:8px 10px;border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{background:var(--panel2);color:var(--dim);font-weight:600;position:sticky;top:0}
tr:hover td{background:#1b1f28}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.mono{font-family:ui-monospace,Consolas,"Courier New",monospace}
.ok{color:var(--ok)}.warn{color:var(--warn)}.bad{color:var(--bad)}.dim{color:var(--dim)}
.pill{display:inline-block;padding:1px 8px;border-radius:10px;font-size:12px;
border:1px solid currentColor}
.design{background:var(--panel);border:1px solid var(--line);border-radius:10px;
padding:16px 18px;margin:14px 0}
.design .top{display:flex;justify-content:space-between;align-items:baseline;gap:12px;flex-wrap:wrap}
.design .id{font-family:ui-monospace,monospace;font-size:13px;color:var(--acc)}
.score{font-size:26px;font-weight:700}
.bar{height:6px;background:#252a34;border-radius:3px;overflow:hidden;margin:8px 0 14px}
.bar>i{display:block;height:100%;background:linear-gradient(90deg,#58a6ff,#3fb950)}
pre{background:#0b0d11;border:1px solid var(--line);border-radius:8px;padding:12px;
overflow-x:auto;font-size:12px;line-height:1.5;margin:8px 0}
details summary{cursor:pointer;color:var(--acc);font-size:13px;margin:6px 0}
footer{color:var(--dim);font-size:12px;padding:24px 32px;border-top:1px solid var(--line)}
.tag{font-size:11px;padding:0 6px;border-radius:4px;background:#252a34;color:var(--dim);margin-left:6px}
"""


def _e(x: Any) -> str:
    return html.escape("" if x is None else str(x))


def _pill(text: str, cls: str) -> str:
    return f'<span class="pill {cls}">{_e(text)}</span>'


def _num(x: Any, nd: int = 2) -> str:
    """None 安全的小数格式化（脱靶分可能缺失）。"""
    if x is None or x == "":
        return "—"
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return _e(x)


def _score_table(score: dict) -> str:
    rows = []
    for c in score.get("components", []):
        rows.append(
            "<tr>"
            f'<td>{_e(c["label_zh"])}<span class="tag">{_e(c["key"])}</span></td>'
            f'<td class="num">{c["raw"]:.1f}</td>'
            f'<td class="num dim">{c["weight"]:.2f}</td>'
            f'<td class="num">{c["weighted"]:.1f}</td>'
            f'<td class="dim">{_e(c["detail"])}</td>'
            "</tr>"
        )
    return (
        "<table><thead><tr><th>分项</th><th>原始分</th><th>权重</th>"
        "<th>加权贡献</th><th>依据 / 扣分原因</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>"
    )


def _design_card(d: dict) -> str:
    sc = d.get("score") or {}
    total = sc.get("total", 0)
    grade = sc.get("grade", "-")
    ot = (d.get("offtarget") or {}).get("summary") or {}
    warns = sc.get("warnings") or []
    left = d["monomer_left"]
    right = d.get("monomer_right")
    warn_html = ""
    if warns:
        warn_html = '<div class="banner">' + "<br>".join("⚠ " + _e(w) for w in warns) + "</div>"
    hits = (d.get("offtarget") or {}).get("hits") or []
    hit_rows = "".join(
        f'<tr><td class="mono">{_e(h["seq_name"])}:{h["start"]}{h["strand"]}</td>'
        f'<td class="mono">{_e(h["left_half"])}</td><td class="mono">{_e(h["spacer"])}</td>'
        f'<td class="mono">{_e(h["right_half"])}</td>'
        f'<td class="num">{h["mm_left"]}/{h["mm_right"]}</td>'
        f'<td class="num">{h["weight"]}</td></tr>'
        for h in hits[:10]
    )
    hits_html = (
        "<h3>脱靶命中（前 10）</h3><table><thead><tr><th>位置</th><th>左半位点</th>"
        "<th>spacer</th><th>右半位点</th><th>错配 L/R</th><th>权重</th></tr></thead><tbody>"
        + hit_rows + "</tbody></table>"
    ) if hits else "<p class='dim'>搜索空间内未发现脱靶命中。</p>"

    right_row = ""
    if right:
        right_row = f"""
    <tr><td>右</td><td class="mono">{_e(right["halfsite"])}</td>
        <td class="mono">{_e(", ".join(right["triplets_N_to_C"]))}</td>
        <td class="mono">{_e(" / ".join(f["helix7"] for f in right["fingers"]))}</td>
        <td class="dim">{_e(" / ".join(f["helix_id"] for f in right["fingers"]))}</td></tr>"""

    kind_txt = ("组合型（成对 ZFN）" if d.get("is_paired", True) else "单体型（单阵列）")
    spacer_txt = (f"spacer {d['target']['spacer_len']} bp"
                  if d.get("is_paired", True) else "无 spacer")
    label_left = "左单体" if d.get("is_paired", True) else "阵列"
    right_pre = (f"\n\n右单体: {_e(wrap((right or {}).get('protein') or '（无）', 70))}"
                 if right else "")

    return f"""
<div class="design">
  <div class="top">
    <div>
      <div class="id">{_e(d["design_id"])}</div>
      <div class="sub">{_e(d["target"]["seq_name"])}:{d["target"]["start"]}{_e(d["target"]["strand"])}
        · {_e(d["architecture"])} · {kind_txt} · {spacer_txt}
        · 最低证据 Tier {d["min_evidence_tier"]}</div>
    </div>
    <div style="text-align:right">
      <div class="score {GRADE_CLASS.get(grade, "")}">{total:.1f} <span style="font-size:14px">Design_score ({_e(grade)})</span></div>
      <div class="sub">脱靶负担 {_e(ot.get("burden"))} · 命中 {_e(ot.get("n_hits"))} · Off-target score {_e((d.get("offtarget_score") or {}).get("score", "—"))}</div>
    </div>
  </div>
  <div class="bar"><i style="width:{min(100, max(0, total)):.1f}%"></i></div>
  {warn_html}
  <table><thead><tr><th>单体</th><th>半位点</th><th>三联体 (N→C)</th><th>识别螺旋</th><th>指头来源</th></tr></thead><tbody>
    <tr><td>左</td><td class="mono">{_e(left["halfsite"])}</td>
        <td class="mono">{_e(", ".join(left["triplets_N_to_C"]))}</td>
        <td class="mono">{_e(" / ".join(f["helix7"] for f in left["fingers"]))}</td>
        <td class="dim">{_e(" / ".join(f["helix_id"] for f in left["fingers"]))}</td></tr>{right_row}
  </tbody></table>
  <h3>打分拆解</h3>
  {_score_table(sc)}
  <details><summary>蛋白序列（草案）</summary>
    <pre>{_e(label_left)}: {_e(wrap(left["protein"] or "（无）", 70))}
{right_pre}</pre>
  </details>
  <details><summary>脱靶明细</summary>{hits_html}</details>
</div>"""


def render_html(payload: dict) -> str:
    meta = payload.get("meta", {})
    lib = payload.get("library", {})
    tsum = payload.get("targets_summary", {})
    sites = payload.get("sites", [])
    designs = payload.get("designs", [])
    fusions = payload.get("fusions", [])
    scoring = payload.get("scoring", {})
    inputs = payload.get("inputs", [])
    warnings = payload.get("warnings", [])

    def _input_row(i: dict) -> str:
        # 基因组模式没有 length/gc，用 n_chroms/total_bp 顶上
        if "length" in i:
            return (f'<tr><td class="mono">{_e(i["name"])}</td><td class="num">{i["length"]}</td>'
                    f'<td class="num">{i["gc"]:.1%}</td><td class="num">{i["n_bases"]}</td>'
                    f'<td class="dim">{_e(i["source"])}</td></tr>')
        return (f'<tr><td class="mono">{_e(i["name"])}</td>'
                f'<td class="num">{i.get("total_bp", 0):,} bp</td>'
                f'<td class="num">{i.get("n_chroms", 1)} 条</td><td class="num">—</td>'
                f'<td class="dim">全基因组</td></tr>')

    input_rows = "".join(_input_row(i) for i in inputs)

    site_rows = []
    for s in sites[:60]:
        cls = STATUS_CLASS.get(s["designability"], "")
        site_rows.append(
            f'<tr><td class="mono">{_e(s["site_id"])}</td>'
            f'<td class="mono">{_e(s["left_half"])}</td>'
            f'<td class="mono">{_e(s["spacer"])}</td>'
            f'<td class="mono">{_e(s["right_half"])}</td>'
            f'<td class="num">{s["spacer_len"]}</td>'
            f'<td>{_pill(STATUS_LABEL.get(s["designability"], s["designability"]), cls)}</td>'
            f'<td class="num">{s["n_combos"]}</td>'
            f'<td class="mono dim">{_e(",".join(s["missing_triplets"]))}</td></tr>'
        )
    if len(sites) > 60:
        site_rows.append(f'<tr><td colspan="8" class="dim">…共 {len(sites)} 个位点，此处只显示前 60 个（完整清单见 sites.tsv）</td></tr>')

    demand_rows = "".join(
        f'<tr><td class="mono">{_e(d["triplet"])}</td><td class="num">{d["occurrences"]}</td>'
        f'<td>{"有候选" if d["covered"] else _pill("缺", "bad")}</td>'
        f'<td class="num">{d["n_candidate_slots_gap"]}</td></tr>'
        for d in (tsum.get("triplet_demand") or [])
    )

    weight_rows = "".join(
        f'<tr><td  class="mono">{_e(k)}</td><td>{_e(scoring.get("labels", {}).get(k, ""))}</td>'
        f'<td class="num">{v:.2f}</td></tr>'
        for k, v in (scoring.get("weights") or {}).items()
    )

    fusion_rows = "".join(
        f'<tr><td class="mono">{_e(f.get("design_id"))}</td><td>{_e(f.get("architecture"))}</td>'
        f'<td class="num">{f.get("length", 0)}</td>'
        f'<td>{"完整" if f.get("is_complete") else "含占位符"}</td>'
        f'<td class="dim">{_e("；".join(f.get("red_flags") or []))}</td></tr>'
        for f in fusions
    )

    warn_html = ""
    if warnings:
        warn_html = '<div class="banner">' + "<br>".join("⚠ " + _e(w) for w in warnings) + "</div>"

    gf = tsum.get("genome_filter") or None
    filter_html = ""
    if gf:
        rows = []
        for rec in gf.get("removed_sites") or []:
            verdict = "重复区" if rec.get("verdict") == "repeat" else "多拷贝"
            rows.append(
                f'<tr><td class="mono">{_e(rec.get("site_id") or "")}</td>'
                f'<td class="mono">{_e(rec.get("motif") or "")}</td>'
                f'<td class="num">{rec.get("genome_copies", 0)}</td>'
                f'<td>{_pill(verdict, "bad")}</td></tr>'
            )
        filter_html = f"""
<h2>2 · Target filter（全基因组重复过滤）</h2>
<p class="sub">判据（宋成义 2026-11-18 口径）：候选 motif 的<b>全基因组落点</b>落在 RepeatMasker
注释的重复区，或同一 motif 在基因组里出现 ≥ {gf.get("min_copies", 3)} 次，即整条剔除。
<b>被剔除的位点不参与打分</b>——所以这份报告里只有通过过滤的候选才有 Design_score。</p>
<div class="cards">
  <div class="card"><div class="k">候选（过滤前）</div><div class="v">{gf.get("n_sites_in", 0)}</div></div>
  <div class="card"><div class="k">落在重复区剔除</div><div class="v bad">{gf.get("n_removed_repeat", 0)}</div></div>
  <div class="card"><div class="k">多拷贝剔除</div><div class="v bad">{gf.get("n_removed_multicopy", 0)}</div></div>
  <div class="card"><div class="k">通过</div><div class="v ok">{gf.get("n_kept", 0)}</div></div>
</div>
<p class="sub">基因组：<span class="mono">{_e(gf.get("genome_fasta"))}</span>
 · RepeatMasker：<span class="mono">{_e(gf.get("repeats_path") or "（未提供）")}</span>
 · motif 数 {gf.get("n_motifs", 0)} · 落点总数 {gf.get("n_landings_total", 0)}</p>
<table><thead><tr><th>位点</th><th>motif</th><th>基因组落点</th><th>判定</th></tr></thead>
<tbody>{"".join(rows) or '<tr><td colspan="4" class="dim">没有候选被剔除</td></tr>'}</tbody></table>
"""

    no = 3 if gf else 2

    ot = payload.get("offtarget_genome") or None
    ot_html = ""
    if ot:
        ots = ot.get("summary") or {}
        ot_rows = []
        for r in ot.get("designs") or []:
            by = r.get("by_mm") or {}
            mm_txt = " / ".join(
                f"{k}:{v}" for k, v in sorted(by.items(), key=lambda kv: int(kv[0]))
            ) or "0"
            ot_rows.append(
                f'<tr><td class="num">{r.get("burden_rank", "")}</td>'
                f'<td class="mono">{_e(r.get("design_id"))}</td>'
                f'<td class="mono">{_e(r.get("seq_name"))}:{r.get("start")}{_e(r.get("strand"))}</td>'
                f'<td class="num">{_num(r.get("design_score"))}</td>'
                f'<td class="num"><b>{_num(r.get("offtarget_score"))}</b></td>'
                f'<td class="num">{_num(r.get("burden"), 3)}</td>'
                f'<td class="num">{mm_txt}</td>'
                f'<td class="num">{r.get("n_hits_total", 0)}</td></tr>'
            )
        ot_html = f"""
<h2>{no + 4} · Off-target score（全基因组）</h2>
<p class="sub">本轮参与脱靶评估的，<b>只有通过 Target filter（全基因组重复过滤）后留下的候选</b>：
落在 RepeatMasker 注释重复区的、以及全基因组落点 ≥3 的，都已在上一步剔除，不会进这里。
对被选中的设计（默认 Design_score 前 {ots.get("n_designs", 0)} 个），把左/右半位点
放到整个基因组里找近似位点：每侧允许 0–2 个错配，按错配数加权（0 mm ×1、1 mm ×0.25、2 mm ×0.06）
累加成 <b>burden</b>，再折成 <b>off-target score = 100 − 25 × burden</b>（下限 0）。
<b>两套分并列，不相加</b>：off-target score 不进 Design_score。
burden 是 in silico 近似位点密度，<b>不是</b>切割概率。</p>
<div class="cards">
  <div class="card"><div class="k">参与评估</div><div class="v">{ots.get("n_designs", 0)}</div></div>
  <div class="card"><div class="k">burden 最低</div><div class="v ok">{_num(ots.get("burden_min"), 3)}</div></div>
  <div class="card"><div class="k">burden 中位</div><div class="v">{_num(ots.get("burden_median"), 3)}</div></div>
  <div class="card"><div class="k">burden 最高</div><div class="v bad">{_num(ots.get("burden_max"), 3)}</div></div>
  <div class="card"><div class="k">扫描规模</div><div class="v">{ots.get("scanned_bp", 0):,} bp</div></div>
  <div class="card"><div class="k">近似位点合计</div><div class="v">{sum((ots.get("by_mm_total") or {}).values())}</div></div>
</div>
<table><thead><tr><th>burden 排名</th><th>设计</th><th>位点</th><th>Design_score</th>
<th>off-target score</th><th>burden</th><th>近位点（总错配:条数）</th><th>命中</th></tr></thead>
<tbody>{"".join(ot_rows) or '<tr><td colspan="8" class="dim">没有设计参与脱靶评估</td></tr>'}</tbody></table>
<p class="sub">每侧允许 ≤{_e((ot.get("designs") or [{}])[0].get("max_mm") or "2")} 个错配，所以总错配是 0–4；
权重 0 mm ×1、1 mm ×0.25、2 mm ×0.06、3 mm ×0.015、4 mm ×0.004。
命中数超过上限时只累计前 N 条，此时 burden 是<b>下限</b>（本报告 {len(ot.get("designs") or [])} 个设计里
{"有" if any(r.get("truncated") for r in (ot.get("designs") or [])) else "没有"}命中被截断）。
burden 里含目标位点自身那 1 条精确命中。引擎 {_e(ots.get("engine") or "")} · 完整逐设计数值见
<span class="mono">offtarget_genome.tsv</span>，近似位点清单见 <span class="mono">report.json</span>。</p>
"""
    fus_no = no + 5 if ot else no + 4
    w_no = no + 6 if ot else no + 5

    return f"""<!DOCTYPE html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ZFN 设计报告 · {_e(meta.get("run_id", ""))}</title>
<style>{_CSS}</style></head>
<body>
<header>
  <h1>ZFN 设计报告 · {_e(meta.get("stage_title") or "S1 DesignScore（all-target）")}</h1>
  <div class="sub">生成时间 {_e(meta.get("generated_at"))} · 引擎版本 {_e(meta.get("version"))}
   · 运行 ID <span class="mono">{_e(meta.get("run_id"))}</span></div>
</header>
<main>
<div class="banner">
  <b>请先读这段。</b> 本报告全部结果是<b>计算预测</b>，不是活性保证。
  分数只表示"按公开设计规则看起来更像靠谱位点"，与实测亲和力/切割效率<b>没有定量关系</b>。
  任何设计进入实验前必须走：T7E1 或 amplicon 深度测序（on-target）+
  in silico 脱靶复核 + 至少两次独立重复。含 Tier2 指头的设计尤其如此。
</div>
{warn_html}

<div class="cards">
  <div class="card"><div class="k">输入序列</div><div class="v">{len(inputs)}</div></div>
  <div class="card"><div class="k">候选位点</div><div class="v">{tsum.get("n_sites", 0)}</div></div>
  <div class="card"><div class="k">可设计（全 Tier1）</div><div class="v ok">{(tsum.get("by_status") or {}).get("full", 0)}</div></div>
  <div class="card"><div class="k">可设计（含未验证）</div><div class="v warn">{(tsum.get("by_status") or {}).get("partial", 0)}</div></div>
  <div class="card"><div class="k">缺口</div><div class="v bad">{(tsum.get("by_status") or {}).get("gap", 0)}</div></div>
  <div class="card"><div class="k">生成设计</div><div class="v">{len(designs)}</div></div>
  <div class="card"><div class="k">指头库记录</div><div class="v">{lib.get("n_records", 0)}</div></div>
  <div class="card"><div class="k">三联体覆盖</div><div class="v">{lib.get("triplets_covered", 0)}/{lib.get("triplets_total", 64)}</div></div>
</div>

<h2>1 · 输入</h2>
<table><thead><tr><th>名称</th><th>长度</th><th>GC</th><th>N 数</th><th>来源</th></tr></thead>
<tbody>{input_rows}</tbody></table>
{filter_html}
<h2>{"3" if gf else "2"} · 指头库</h2>
<table><thead><tr><th>库</th><th>记录数</th><th>Tier1 记录</th><th>覆盖三联体</th><th>按证据分布</th></tr></thead>
<tbody><tr><td class="mono">{_e(lib.get("library_name"))}</td>
<td class="num">{lib.get("n_records")}</td><td class="num">{lib.get("n_curated_records")}</td>
<td class="num">{lib.get("triplets_covered")}/{lib.get("triplets_total")}</td>
<td class="dim">{_e(lib.get("by_evidence"))}</td></tr></tbody></table>
<p class="sub">覆盖不足会直接体现为下面的"缺口"。补齐覆盖的办法见 docs/数据来源与许可.md：
申请 Persikov 2015 的 B1H 数据集后，用 <span class="mono">zfn-design library import-b1h</span> 导入。</p>

<h2>{no + 1} · 候选位点与可设计性</h2>
<table><thead><tr><th>位点</th><th>左半位点</th><th>spacer</th><th>右半位点</th>
<th>spacer bp</th><th>可设计性</th><th>组合数</th><th>缺失三联体</th></tr></thead>
<tbody>{"".join(site_rows)}</tbody></table>

<h2>{no + 2} · 这份输入最需要的三联体</h2>
<p class="sub">按出现次数排序，未覆盖的排在前面。种子库覆盖不足时，<b>不必一次补全 64 个</b>——
照这张表从前往后补最有价值。</p>
<table><thead><tr><th>三联体</th><th>出现次数（指位）</th><th>当前库</th><th>缺口次数</th></tr></thead>
<tbody>{demand_rows}</tbody></table>

<h2>{no + 3} · 设计（按 Design_score 排序）</h2>{"".join(_design_card(d) for d in designs) if designs else '<p class="dim">没有可输出的设计——通常是位点全部落在缺口上。加大输入长度或补充指头库。</p>'}
{ot_html}
<h2>{fus_no} · 融合架构（fusion）</h2>
<table><thead><tr><th>设计</th><th>架构</th><th>长度 aa</th><th>完整性</th><th>红旗</th></tr></thead>
<tbody>{fusion_rows}</tbody></table>
<p class="sub">S3 跑完之后，这批设计同样可以接着做 Fusion design（S4）：挑设计、接 VP64 / KRAB / FokI，
生成蛋白与人源化 CDS，并在订单表里带上该设计的 off-target score 与 burden。
本次 {("已" if ot else "还没有")}算过全基因组脱靶。</p>

<h2>{w_no} · 打分权重（可审计）</h2>
<table><thead><tr><th>key</th><th>含义</th><th>权重</th></tr></thead><tbody>{weight_rows}</tbody></table>
<p class="sub">权重之和 = {sum((scoring.get("weights") or {}).values()):.2f}。
这些是工程默认值，<b>不是</b>从实验数据拟合的；改动它们会直接改变排序，
所以放在 <span class="mono">scoring_weights.json</span> 里而不是写死在代码里。</p>
</main>
<footer>
本报告由 ZFN_Design 生成 · 计算预测，需实验验证 ·
指头数据来源见 docs/数据来源与许可.md
</footer>
</body></html>"""


def write_html(path: str | Path, payload: dict) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(render_html(payload), encoding="utf-8")
    return p


def write_bundle(path: str | Path, files: dict[str, Path]) -> Path:
    """把一次运行的全部产物打成一个 zip，方便整包交付。"""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    written: set[str] = set()
    with zipfile.ZipFile(p, "w", zipfile.ZIP_DEFLATED) as z:
        for f in files.values():
            if f is None:
                continue
            f = Path(f)
            if not f.is_file() or f == p or f.name in written:
                continue
            z.write(f, arcname=f.name)
            written.add(f.name)
    return p


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d %H:%M:%S %z")


__all__ = ["write_json", "write_tsv", "sites_tsv", "designs_tsv",
           "render_html", "write_html", "write_bundle", "now_iso"]
