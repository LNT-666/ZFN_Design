"""使用说明书：把 data/manual_zh.md / manual_en.md 渲染成网页，并提供下载。

单一真源是 data/ 下的两个 markdown 文件；网页版本与 docx 下载都由它们生成。
"""
from __future__ import annotations

import html
import re
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"

LANGS = ("zh", "en")
TITLES = {
    "zh": "ZFN 设计工作台 · 程序使用说明书",
    "en": "ZFN Design Workbench · User Manual",
}
SWITCH = {"zh": "EN", "en": "中文"}
SWITCH_HREF = {"zh": "en", "en": "zh"}
BACK_LABEL = {"zh": "← 回到工作台", "en": "← Back to the workbench"}
DL_LABEL = {
    "zh": {"docx": "下载 Word 版", "md": "下载 Markdown 版"},
    "en": {"docx": "Download Word", "md": "Download Markdown"},
}

_CSS = """
:root{--bg:#0b1220;--card:#111a2b;--line:#1e2b44;--fg:#dbe6f5;--dim:#8fa3c0;--acc:#4a8fdc}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font:15px/1.75 -apple-system,BlinkMacSystemFont,"Segoe UI","Microsoft YaHei",sans-serif}
.wrap{max-width:980px;margin:0 auto;padding:28px 22px 80px}
.bar{position:sticky;top:0;background:rgba(11,18,32,.94);backdrop-filter:blur(6px);
  border-bottom:1px solid var(--line);padding:12px 22px;display:flex;gap:14px;
  align-items:center;flex-wrap:wrap;z-index:9}
.bar a{color:var(--acc);text-decoration:none;font-size:14px}
.bar .sp{flex:1}
.bar .tag{color:var(--dim);font-size:13px}
h1{font-size:26px;margin:26px 0 6px;letter-spacing:.3px}
h2{font-size:20px;margin:34px 0 10px;padding-bottom:8px;border-bottom:1px solid var(--line)}
h3{font-size:16px;margin:22px 0 6px;color:#a9c4e6}
p{margin:10px 0}
ul,ol{margin:10px 0 10px 22px;padding:0}
li{margin:4px 0}
code{background:#0e1726;border:1px solid var(--line);border-radius:4px;
  padding:1px 5px;font-size:13px;color:#9fd0ff}
pre{background:#0e1726;border:1px solid var(--line);border-radius:8px;padding:12px 14px;
  overflow:auto;font-size:13px}
pre code{border:0;background:none;padding:0}
table{border-collapse:collapse;width:100%;margin:14px 0;font-size:14px}
th,td{border:1px solid var(--line);padding:7px 10px;text-align:left;vertical-align:top}
th{background:#14203a;color:#cfe0f7;font-weight:600}
tr:nth-child(even) td{background:#0e1726}
blockquote{margin:12px 0;padding:8px 14px;border-left:3px solid var(--acc);
  background:#0e1726;color:var(--dim)}
.foot{margin-top:46px;color:var(--dim);font-size:13px;border-top:1px solid var(--line);padding-top:14px}
"""


def _inline(text: str) -> str:
    """行内标记：`code`、**bold**、[a](b)。先转义再替换。"""
    out = html.escape(text, quote=False)
    out = re.sub(r"`([^`]+)`", r"<code>\1</code>", out)
    out = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", out)
    out = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', out)
    return out


def md_to_html(md: str) -> str:
    """够用就行的 markdown → HTML：标题、表格、列表、段落、代码块、引用。"""
    lines = md.replace("\r\n", "\n").split("\n")
    out: list[str] = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        s = line.strip()

        if not s:
            i += 1
            continue

        # 代码块
        if s.startswith("```"):
            i += 1
            buf = []
            while i < n and not lines[i].strip().startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            out.append("<pre><code>" + html.escape("\n".join(buf)) + "</code></pre>")
            continue

        # 标题
        m = re.match(r"^(#{1,6})\s+(.*)$", s)
        if m:
            level = len(m.group(1))
            out.append(f"<h{level}>{_inline(m.group(2))}</h{level}>")
            i += 1
            continue

        # 表格
        if s.startswith("|") and i + 1 < n and re.match(r"^\|[\s:|-]+\|$", lines[i + 1].strip()):
            header = [c.strip() for c in s.strip("|").split("|")]
            i += 2
            rows = []
            while i < n and lines[i].strip().startswith("|"):
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            out.append("<table><thead><tr>")
            out.extend(f"<th>{_inline(c)}</th>" for c in header)
            out.append("</tr></thead><tbody>")
            for r in rows:
                out.append("<tr>")
                out.extend(f"<td>{_inline(c)}</td>" for c in r)
                out.append("</tr>")
            out.append("</tbody></table>")
            continue

        # 列表
        if re.match(r"^[-*]\s+", s) or re.match(r"^\d+\.\s+", s):
            ordered = bool(re.match(r"^\d+\.\s+", s))
            tag = "ol" if ordered else "ul"
            out.append(f"<{tag}>")
            while i < n:
                cur = lines[i].strip()
                mm = re.match(r"^[-*]\s+(.*)$", cur) or re.match(r"^\d+\.\s+(.*)$", cur)
                if not mm:
                    break
                out.append(f"<li>{_inline(mm.group(1))}</li>")
                i += 1
            out.append(f"</{tag}>")
            continue

        # 引用
        if s.startswith(">"):
            buf = []
            while i < n and lines[i].strip().startswith(">"):
                buf.append(lines[i].strip()[1:].strip())
                i += 1
            out.append("<blockquote>" + _inline(" ".join(buf)) + "</blockquote>")
            continue

        # 段落（连续非空行合并）
        buf = [s]
        i += 1
        while i < n:
            nxt = lines[i].strip()
            if (not nxt or nxt.startswith(("#", "|", ">", "```"))
                    or re.match(r"^[-*]\s+", nxt) or re.match(r"^\d+\.\s+", nxt)):
                break
            buf.append(nxt)
            i += 1
        out.append("<p>" + _inline(" ".join(buf)) + "</p>")

    return "\n".join(out)


def _normalise(lang: str | None) -> str:
    lang = (lang or "zh").lower()
    return lang if lang in LANGS else "zh"


def manual_path(lang: str | None = None) -> Path:
    return DATA_DIR / f"manual_{_normalise(lang)}.md"


def read_manual(lang: str | None = None) -> str:
    p = manual_path(lang)
    if not p.exists():
        return "# Manual not found\n\n`manual_{}_md` is missing from the package data directory.".format(lang)
    return p.read_text(encoding="utf-8")


def render_page(lang: str | None = None) -> str:
    """整页 HTML：顶部有返回、语言切换、下载链接。"""
    lang = _normalise(lang)
    other = SWITCH_HREF[lang]
    body = md_to_html(read_manual(lang))
    dl = DL_LABEL[lang]
    return f"""<!doctype html>
<html lang="{lang}">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(TITLES[lang])}</title>
<style>{_CSS}</style>
</head>
<body>
<div class="bar">
  <a href="/">{html.escape(BACK_LABEL[lang])}</a>
  <span class="sp"></span>
  <a href="/manual?lang={other}">{SWITCH[lang]}</a>
  <a href="/manual/download?lang={lang}&amp;fmt=docx">{html.escape(dl['docx'])}</a>
  <a href="/manual/download?lang={lang}&amp;fmt=md">{html.escape(dl['md'])}</a>
</div>
<div class="wrap">
{body}
<div class="foot">zfn_design 0.7.0 · {html.escape(TITLES[lang])}</div>
</div>
</body>
</html>
"""


__all__ = ["LANGS", "TITLES", "manual_path", "read_manual", "md_to_html", "render_page"]
