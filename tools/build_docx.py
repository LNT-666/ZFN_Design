#!/usr/bin/env python3
"""把 docs/*.md 编译成一份排版好的 Word 文档。

用法::

    python tools/build_docx.py                # 默认输出 docs/ZFN设计方案.docx
    python tools/build_docx.py --out x.docx

支持的 Markdown 子集（刻意保持小，够用即可）：

* 标题 ``#`` / ``##`` / ``###``
* 段落（含 ``**粗体**``、``*斜体*``、`` `代码` ``、``[文字](url)``、``<sup>`` 标签）
* 无序列表 ``-`` / ``*``（两层缩进）与有序列表 ``1.``
* 管道表格（带 ``|---|`` 对齐行）
* fenced code block（```）
* 引用 ``>``
* 分隔线 ``---``

刻意不做的：图片、脚注、LaTeX 数学（公式在 md 里已经写成等宽文本了）。
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

# ---------------------------------------------------------------- 版式常量
BODY_LATIN = "Calibri"
BODY_CJK = "宋体"
HEAD_LATIN = "Segoe UI Semibold"
HEAD_CJK = "微软雅黑"
MONO = "Consolas"
ACCENT = RGBColor(0x1F, 0x4E, 0x79)
GREY = RGBColor(0x59, 0x59, 0x59)
CODE_SHADE = "F2F4F7"
HEAD_SHADE = "E8EEF6"


# ---------------------------------------------------------------- 低层工具
def _set_run_font(run, latin: str, cjk: str, size: float | None = None,
                  bold: bool | None = None, color: RGBColor | None = None) -> None:
    run.font.name = latin
    rpr = run._element.get_or_add_rPr()
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts")
        rpr.append(rfonts)
    rfonts.set(qn("w:ascii"), latin)
    rfonts.set(qn("w:hAnsi"), latin)
    rfonts.set(qn("w:eastAsia"), cjk)
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.font.bold = bold
    if color is not None:
        run.font.color.rgb = color


def _shade(element, fill: str) -> None:
    """给段落或单元格铺底色。"""
    pr = element.get_or_add_pPr() if hasattr(element, "get_or_add_pPr") else element
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), fill)
    pr.append(shd)


def _para_border(paragraph, side: str = "bottom", color: str = "D0D7E2",
                 size: int = 6) -> None:
    ppr = paragraph._p.get_or_add_pPr()
    borders = ppr.find(qn("w:pBdr"))
    if borders is None:
        borders = OxmlElement("w:pBdr")
        ppr.append(borders)
    el = OxmlElement(f"w:{side}")
    el.set(qn("w:val"), "single")
    el.set(qn("w:sz"), str(size))
    el.set(qn("w:space"), "4")
    el.set(qn("w:color"), color)
    borders.append(el)


def _add_field(paragraph, instr: str, placeholder: str = "") -> None:
    run = paragraph.add_run()
    fld_begin = OxmlElement("w:fldChar")
    fld_begin.set(qn("w:fldCharType"), "begin")
    instr_el = OxmlElement("w:instrText")
    instr_el.set(qn("xml:space"), "preserve")
    instr_el.text = instr
    fld_sep = OxmlElement("w:fldChar")
    fld_sep.set(qn("w:fldCharType"), "separate")
    txt = OxmlElement("w:t")
    txt.text = placeholder
    fld_end = OxmlElement("w:fldChar")
    fld_end.set(qn("w:fldCharType"), "end")
    for el in (fld_begin, instr_el, fld_sep, txt, fld_end):
        run._r.append(el)


# ---------------------------------------------------------------- 行内解析
INLINE_RE = re.compile(
    r"(\*\*.+?\*\*"          # **bold**
    r"|`[^`]+`"              # `code`
    # [text](url)：URL 里允许一层成对括号（DOI 形如 10.1016/S0969-2126(96)00125-6）
    r"|\[[^\]]+\]\((?:https?://(?:[^()\s]|\([^()\s]*\))*|[^)]+)\)"
    r"|\*[^*\n]+?\*"         # *italic*
    r"|<sup>.*?</sup>"       # <sup>..</sup>
    r"|<https?://[^>]+>"     # <https://...> 自动链接
    r"|<[a-zA-Z/][^>]*>)"    # 其它裸标签，直接丢掉
)


def add_inline(paragraph, text: str, *, size: float = 10.5,
               color: RGBColor | None = None) -> None:
    """把一段带行内标记的文本写进 paragraph。"""
    pos = 0
    for m in INLINE_RE.finditer(text):
        if m.start() > pos:
            _plain(paragraph, text[pos:m.start()], size, color)
        token = m.group(0)
        if token.startswith("**"):
            inner = token[2:-2]
            # 粗体里常嵌行内代码（如 **`mode="paired"`（组合型）**），
            # 直接把反引号去掉，用等宽字体渲染整段
            if "`" in inner:
                r = paragraph.add_run(inner.replace("`", ""))
                _set_run_font(r, MONO, MONO, size - 0.5, bold=True,
                              color=RGBColor(0xA6, 0x2A, 0x2A))
            else:
                r = paragraph.add_run(inner)
                _set_run_font(r, BODY_LATIN, BODY_CJK, size, bold=True, color=color)
        elif token.startswith("`"):
            r = paragraph.add_run(token[1:-1])
            _set_run_font(r, MONO, MONO, size - 0.5, color=RGBColor(0xA6, 0x2A, 0x2A))
        elif token.startswith("["):
            mm = re.match(r"\[([^\]]+)\]\((.*)\)$", token, re.S)
            if not mm:
                _plain(paragraph, token, size, color)
                pos = m.end()
                continue
            label, url = mm.group(1), mm.group(2)
            r = paragraph.add_run(label)
            _set_run_font(r, BODY_LATIN, BODY_CJK, size, color=ACCENT)
            r.font.underline = True
            # 链接文字本身已经说明了去处（DOI / 完整 URL）时不再重复打印一遍
            shown = label.strip().rstrip("/").lower()
            target = url.strip().rstrip("/").lower()
            redundant = shown == target or target in shown or shown in target
            if url.startswith("http") and not redundant:
                r2 = paragraph.add_run(f"（{url}）")
                _set_run_font(r2, BODY_LATIN, BODY_CJK, size - 1.5, color=GREY)
        elif token.startswith("*"):
            r = paragraph.add_run(token[1:-1])
            _set_run_font(r, BODY_LATIN, BODY_CJK, size, color=color)
            r.font.italic = True
        elif token.startswith("<http"):
            r = paragraph.add_run(token[1:-1])
            _set_run_font(r, BODY_LATIN, BODY_CJK, size, color=ACCENT)
            r.font.underline = True
        elif token.startswith("<sup>"):
            r = paragraph.add_run(token[5:-6])
            _set_run_font(r, BODY_LATIN, BODY_CJK, size - 1.5, color=color)
            r.font.superscript = True
        pos = m.end()
    if pos < len(text):
        _plain(paragraph, text[pos:], size, color)


def _plain(paragraph, text: str, size: float, color: RGBColor | None) -> None:
    if not text:
        return
    r = paragraph.add_run(text)
    _set_run_font(r, BODY_LATIN, BODY_CJK, size, color=color)


# ---------------------------------------------------------------- 块级解析
def parse_markdown(text: str) -> list[tuple]:
    """极简 block parser。返回 ('kind', payload, ...) 元组列表。"""
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    blocks: list[tuple] = []
    i = 0
    n = len(lines)

    def is_table_sep(s: str) -> bool:
        return bool(re.fullmatch(r"\|[\s:\-|]+\|", s.strip())) and "-" in s

    while i < n:
        line = lines[i]
        stripped = line.strip()

        # 空行
        if not stripped:
            i += 1
            continue

        # 代码块
        if stripped.startswith("```"):
            lang = stripped[3:].strip()
            i += 1
            buf = []
            while i < n and not lines[i].strip().startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            blocks.append(("code", lang, buf))
            continue

        # 分隔线
        if re.fullmatch(r"-{3,}|\*{3,}|_{3,}", stripped):
            blocks.append(("hr",))
            i += 1
            continue

        # 标题
        m = re.match(r"^(#{1,4})\s+(.*)$", stripped)
        if m:
            level = len(m.group(1))
            blocks.append(("h", level, m.group(2).strip()))
            i += 1
            continue

        # 表格
        if stripped.startswith("|") and i + 1 < n and is_table_sep(lines[i + 1]):
            header = _split_row(stripped)
            i += 2
            rows = []
            while i < n and lines[i].strip().startswith("|"):
                rows.append(_split_row(lines[i].strip()))
                i += 1
            blocks.append(("table", header, rows))
            continue

        # 引用
        if stripped.startswith(">"):
            buf = []
            while i < n and lines[i].strip().startswith(">"):
                buf.append(lines[i].strip().lstrip(">").strip())
                i += 1
            blocks.append(("quote", " ".join(buf)))
            continue

        # 列表
        m = re.match(r"^(\s*)([-*+]|\d+\.)\s+(.*)$", line)
        if m:
            items = []
            ordered = m.group(2)[0].isdigit()
            while i < n:
                mm = re.match(r"^(\s*)([-*+]|\d+\.)\s+(.*)$", lines[i])
                if not mm:
                    # 列表项的续行
                    if lines[i].strip() and items and not lines[i].strip().startswith("|"):
                        prev = items[-1]
                        items[-1] = (prev[0], prev[1], prev[2] + " " + lines[i].strip())
                        i += 1
                        continue
                    break
                indent = len(mm.group(1).expandtabs(4)) // 2
                marker = mm.group(2) if ordered else None
                items.append((min(indent, 1), marker, mm.group(3).strip()))
                i += 1
            blocks.append(("list", ordered, items))
            continue

        # 普通段落：吃到空行或下一个块起始
        buf = [stripped]
        i += 1
        while i < n:
            nxt = lines[i]
            s = nxt.strip()
            if (not s or s.startswith("```") or s.startswith("|")
                    or s.startswith(">") or re.match(r"^#{1,4}\s", s)
                    or re.fullmatch(r"-{3,}", s)
                    or re.match(r"^(\s*)([-*+]|\d+\.)\s+", nxt)):
                break
            buf.append(s)
            i += 1
        blocks.append(("p", " ".join(buf)))

    return blocks


def _display_width(s: str) -> int:
    """粗略显示宽度：CJK 记 2，其余记 1（用于算列宽）。"""
    return sum(2 if ord(ch) > 0x2E80 else 1 for ch in s)


def _split_row(line: str) -> list[str]:
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [c.strip() for c in line.split("|")]


# ---------------------------------------------------------------- 渲染
class Builder:
    def __init__(self, doc: Document) -> None:
        self.doc = doc
        self.first_block_done = False

    def render(self, blocks: list[tuple], *, skip_h1: bool = False) -> None:
        for blk in blocks:
            kind = blk[0]
            if kind == "h":
                self._heading(blk[1], blk[2], skip_h1=skip_h1)
            elif kind == "p":
                self._paragraph(blk[1])
            elif kind == "code":
                self._code(blk[2])
            elif kind == "table":
                self._table(blk[1], blk[2])
            elif kind == "list":
                self._list(blk[1], blk[2])
            elif kind == "quote":
                self._quote(blk[1])
            elif kind == "hr":
                self._hr()
            self.first_block_done = True

    # -- 各类元素
    def _heading(self, level: int, text: str, *, skip_h1: bool = False) -> None:
        if level == 1 and skip_h1 and not self.first_block_done:
            return
        if level == 1 and self.first_block_done:
            self.doc.add_page_break()
        p = self.doc.add_paragraph(style=f"Heading {min(level, 3)}")
        sizes = {1: 18, 2: 14, 3: 12, 4: 11}
        add_inline(p, text, size=sizes.get(level, 11))
        for r in p.runs:
            _set_run_font(r, HEAD_LATIN, HEAD_CJK, sizes.get(level, 11), bold=True, color=ACCENT)
        p.paragraph_format.space_before = Pt(18 if level == 1 else 12)
        p.paragraph_format.space_after = Pt(6 if level == 1 else 4)
        p.paragraph_format.keep_with_next = True
        if level >= 2:
            _para_border(p, "bottom", "DCE4F0", 4)

    def _paragraph(self, text: str) -> None:
        p = self.doc.add_paragraph(style="Body Text")
        add_inline(p, text)
        p.paragraph_format.space_after = Pt(6)
        p.paragraph_format.line_spacing = 1.25
        # 左对齐而非两端对齐：参考文献里长 DOI / 作者串会把两端对齐的行拉成大洞
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT

    def _code(self, lines: list[str]) -> None:
        p = self.doc.add_paragraph()
        pf = p.paragraph_format
        pf.left_indent = Cm(0.3)
        pf.right_indent = Cm(0.3)
        pf.space_before = Pt(4)
        pf.space_after = Pt(8)
        pf.line_spacing = 1.0
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT
        _shade(p._p, CODE_SHADE)
        _para_border(p, "left", "9DB2D0", 18)
        for idx, line in enumerate(lines):
            r = p.add_run(line if line else " ")
            _set_run_font(r, MONO, MONO, 8.5)
            r.font.color.rgb = RGBColor(0x22, 0x2A, 0x35)
            if idx < len(lines) - 1:
                r.add_break(WD_BREAK.LINE)

    def _table(self, header: list[str], rows: list[list[str]]) -> None:
        ncol = len(header)
        table = self.doc.add_table(rows=1, cols=ncol)
        table.style = "Table Grid"
        table.alignment = WD_TABLE_ALIGNMENT.CENTER
        table.autofit = False

        # 按内容显示宽度分配列宽（总宽 16.2 cm）
        raw = []
        for j in range(ncol):
            cells_w = [_display_width(header[j])]
            cells_w += [_display_width(r[j]) for r in rows if j < len(r)]
            raw.append(max(4, max(cells_w)))
        total = sum(raw)
        avail = 16.2
        widths = [max(1.5, avail * w / total) for w in raw]
        scale = avail / sum(widths)
        widths = [w * scale for w in widths]

        tbl_pr = table._tbl.tblPr
        layout = OxmlElement("w:tblLayout")
        layout.set(qn("w:type"), "fixed")
        tbl_pr.append(layout)

        hdr = table.rows[0].cells
        for j, cell_text in enumerate(header):
            cell = hdr[j]
            cell.width = Cm(widths[j])
            cell.text = ""
            p = cell.paragraphs[0]
            p.paragraph_format.space_after = Pt(2)
            p.paragraph_format.space_before = Pt(2)
            add_inline(p, cell_text, size=9.5)
            for r in p.runs:
                _set_run_font(r, HEAD_LATIN, HEAD_CJK, 9.5, bold=True)
            _shade(cell._tc.get_or_add_tcPr(), HEAD_SHADE)

        for row in rows:
            cells = table.add_row().cells
            for j in range(ncol):
                cell = cells[j]
                cell.width = Cm(widths[j])
                cell.text = ""
                p = cell.paragraphs[0]
                p.paragraph_format.space_after = Pt(2)
                p.paragraph_format.space_before = Pt(2)
                add_inline(p, row[j] if j < len(row) else "", size=9)
                for r in p.runs:
                    _set_run_font(r, BODY_LATIN, BODY_CJK, 9)

        # 表头行跨页重复；所有行都不允许在行内被页边界劈开
        for row_obj in table.rows:
            tr_pr = row_obj._tr.get_or_add_trPr()
            cant = OxmlElement("w:cantSplit")
            cant.set(qn("w:val"), "true")
            tr_pr.append(cant)
        tr_pr = table.rows[0]._tr.get_or_add_trPr()
        tbl_header = OxmlElement("w:tblHeader")
        tbl_header.set(qn("w:val"), "true")
        tr_pr.append(tbl_header)

        self.doc.add_paragraph().paragraph_format.space_after = Pt(2)

    def _list(self, ordered: bool, items: list[tuple[int, str | None, str]]) -> None:
        """列表项自己画 marker，不用 List Number/List Bullet 样式。

        原因：python-docx 的 List Number 是全文档共享的同一个 numbering 实例，
        多个有序列表会连续编号（§2 于是从 6. 开始）；样式缩进也与正文栅格对不齐。
        这里直接写 "1." / "•" 文本并设悬挂缩进，编号忠实于 Markdown 源。
        """
        for indent, marker, text in items:
            p = self.doc.add_paragraph()
            pf = p.paragraph_format
            base = 0.85 if indent == 0 else 1.45
            pf.left_indent = Cm(base)
            pf.first_line_indent = Cm(-0.5)
            pf.space_after = Pt(2)
            pf.line_spacing = 1.15
            pf.tab_stops.add_tab_stop(Cm(base))
            label = marker if (ordered and marker) else ("•" if indent == 0 else "–")
            run = p.add_run(f"{label}\t")
            _set_run_font(run, BODY_LATIN, BODY_CJK, 10.5)
            add_inline(p, text)

    def _quote(self, text: str) -> None:
        p = self.doc.add_paragraph()
        pf = p.paragraph_format
        pf.left_indent = Cm(0.6)
        pf.space_before = Pt(4)
        pf.space_after = Pt(8)
        _shade(p._p, "F7F9FC")
        _para_border(p, "left", "9DB2D0", 18)
        add_inline(p, text, size=10, color=RGBColor(0x33, 0x3A, 0x44))
        for r in p.runs:
            r.font.italic = True

    def _hr(self) -> None:
        p = self.doc.add_paragraph()
        p.paragraph_format.space_before = Pt(2)
        p.paragraph_format.space_after = Pt(2)
        _para_border(p, "bottom", "DCE4F0", 6)


# ---------------------------------------------------------------- 文档骨架
def _setup_styles(doc: Document) -> None:
    normal = doc.styles["Normal"]
    normal.font.name = BODY_LATIN
    normal.font.size = Pt(10.5)
    normal.element.rPr.rFonts.set(qn("w:eastAsia"), BODY_CJK)
    normal.paragraph_format.line_spacing = 1.25

    for name in ("Body Text", "List Bullet", "List Bullet 2", "List Number", "List Number 2"):
        try:
            st = doc.styles[name]
        except KeyError:
            continue
        st.font.name = BODY_LATIN
        st.font.size = Pt(10.5)


def _title_page(doc: Document, title: str, subtitle: str, meta: list[str]) -> None:
    for _ in range(5):
        doc.add_paragraph()

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(title)
    _set_run_font(r, HEAD_LATIN, HEAD_CJK, 28, bold=True, color=ACCENT)

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_before = Pt(10)
    r = p.add_run(subtitle)
    _set_run_font(r, BODY_LATIN, BODY_CJK, 12.5, color=GREY)

    for _ in range(6):
        doc.add_paragraph()

    for line in meta:
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_after = Pt(3)
        r = p.add_run(line)
        _set_run_font(r, BODY_LATIN, BODY_CJK, 10.5, color=GREY)

    doc.add_page_break()


def _toc_page(doc: Document) -> None:
    p = doc.add_paragraph()
    r = p.add_run("目　录")
    _set_run_font(r, HEAD_LATIN, HEAD_CJK, 16, bold=True, color=ACCENT)
    p.paragraph_format.space_after = Pt(10)

    p = doc.add_paragraph()
    _add_field(p, r'TOC \o "1-2" \h \z \u', "［目录将在 Word 打开时自动生成；若未生成，请按 Ctrl+A 再按 F9］")

    p = doc.add_paragraph()
    r = p.add_run("提示：本文件开启了“打开时更新域”，Word 一般会自动填好目录与页码；"
                  "若你的 Word 禁用了自动更新，按 Ctrl+A 再按 F9 手动刷新即可。")
    _set_run_font(r, BODY_LATIN, BODY_CJK, 9, color=GREY)
    doc.add_page_break()


def _enable_update_fields(doc: Document) -> None:
    """让 Word 打开文档时自动更新域（目录、页码），省得用户手动按 F9。"""
    settings = doc.settings.element
    el = OxmlElement("w:updateFields")
    el.set(qn("w:val"), "true")
    settings.append(el)


def _footer_page_numbers(doc: Document) -> None:
    for section in doc.sections:
        footer = section.footer
        p = footer.paragraphs[0] if footer.paragraphs else footer.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        _add_field(p, "PAGE", "1")
        for r in p.runs:
            _set_run_font(r, BODY_LATIN, BODY_CJK, 9, color=GREY)


# ---------------------------------------------------------------- main
def build(md_paths: list[Path], out_path: Path) -> Path:
    doc = Document()
    _setup_styles(doc)

    sec = doc.sections[0]
    sec.page_width = Cm(21.0)   # A4
    sec.page_height = Cm(29.7)
    sec.top_margin = Cm(2.2)
    sec.bottom_margin = Cm(2.0)
    sec.left_margin = Cm(2.4)
    sec.right_margin = Cm(2.4)

    _title_page(
        doc,
        "ZFN 锌指核酸酶设计与评分系统",
        "从 DNA 序列到可审计设计方案的完整说明",
        [
            "版本 0.1.0",
            "2026-11-14",
            "",
            "配套代码：ZFN_Design（零第三方依赖，57 项测试）",
            "本文档由 docs/*.md 自动编译，请以 Markdown 源文件为准",
        ],
    )
    _toc_page(doc)

    builder = Builder(doc)
    for idx, path in enumerate(md_paths):
        blocks = parse_markdown(path.read_text(encoding="utf-8"))
        builder.render(blocks, skip_h1=(idx == 0))
        if idx < len(md_paths) - 1:
            doc.add_page_break()

    _footer_page_numbers(doc)
    _enable_update_fields(doc)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(out_path))
    return out_path


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser(description="把 docs/*.md 编译成 docx")
    ap.add_argument("--out", default=str(root / "docs" / "ZFN设计方案.docx"))
    ap.add_argument("--docs", nargs="*", default=None,
                    help="默认按顺序编译 设计方案/术语表/数据来源与许可")
    args = ap.parse_args(argv)

    docs = ([Path(p) for p in args.docs] if args.docs else [
        root / "docs" / "设计方案.md",
        root / "docs" / "术语表.md",
        root / "docs" / "数据来源与许可.md",
    ])
    missing = [p for p in docs if not p.exists()]
    if missing:
        print("找不到: " + ", ".join(str(p) for p in missing), file=sys.stderr)
        return 1

    out = build(docs, Path(args.out))
    size = out.stat().st_size
    print(f"已生成 {out}  ({size:,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
