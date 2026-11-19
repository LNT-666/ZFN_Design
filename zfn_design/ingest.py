"""阶段 1：ingest——把用户手里的东西变成程序能吃的序列。

支持三种输入：

* ``--seq ACGT...``           直接给序列
* ``--fasta path``            单个或多个 FASTA 记录（.fa/.fasta/.fna/.txt）
* ``--stdin``                 从标准输入读（管道用）

产出 :class:`~zfn_design.models.SequenceRecord` 列表。这里刻意**不**做任何
生物学解释，只做"清洗 + 校验 + 记录 N 的位置"，因为后面每一步都要能倒回
到原始坐标。
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable, TextIO

from .models import SequenceRecord
from .sequtil import SequenceError, clean, validate


class IngestError(ValueError):
    """输入无法解析。"""


def parse_fasta(text: str, *, source: str = "") -> list[SequenceRecord]:
    """解析 FASTA 文本。没有 ``>`` 头时整段当一条序列，名为 ``seq1``。"""
    records: list[SequenceRecord] = []
    name: str | None = None
    desc = ""
    buf: list[str] = []

    def flush() -> None:
        if name is None and not buf:
            return
        raw = "".join(buf)
        if not raw.strip():
            return
        seq = validate(clean(raw), name=name or "seq1", allow_n=True)
        records.append(SequenceRecord(name=name or "seq1", seq=seq, description=desc, source=source))

    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            parts = header.split(None, 1)
            name = parts[0] if parts else "seq1"
            desc = parts[1] if len(parts) > 1 else ""
            buf = []
        else:
            buf.append(line)
    flush()
    if not records:
        raise IngestError("FASTA 里没有可用序列")
    return records


def read_input(*, sequence: str | None = None, fasta: str | Path | None = None,
               stdin: bool = False, stdin_stream: TextIO | None = None) -> list[SequenceRecord]:
    """统一入口。三个来源必须且只能给一个。"""
    given = [x is not None and x is not False for x in (sequence, fasta, stdin or None)]
    if sum(1 for g in given if g) != 1:
        raise IngestError("必须且只能指定 --seq / --fasta / --stdin 之一")

    if sequence is not None:
        seq = validate(clean(sequence), name="<--seq>", allow_n=True)
        return [SequenceRecord(name="input", seq=seq, source="cli")]

    if fasta is not None:
        p = Path(fasta)
        if not p.exists():
            raise IngestError(f"文件不存在: {p}")
        text = p.read_text(encoding="utf-8-sig", errors="replace")
        if text.lstrip().startswith(">"):
            return parse_fasta(text, source=str(p))
        return [SequenceRecord(name=p.stem, seq=validate(clean(text), name=p.name), source=str(p))]

    stream = stdin_stream or sys.stdin
    text = stream.read()
    if text.lstrip().startswith(">"):
        return parse_fasta(text, source="stdin")
    return [SequenceRecord(name="stdin", seq=validate(clean(text), name="stdin"), source="stdin")]


def summarize(records: Iterable[SequenceRecord]) -> list[dict[str, object]]:
    """给报告用的输入概览。"""
    out = []
    for r in records:
        n_count = r.seq.count("N")
        out.append({
            "name": r.name,
            "length": r.length,
            "n_bases": n_count,
            "gc": round(sum(1 for c in r.seq if c in "GC") / max(1, r.length), 4),
            "source": r.source,
            "description": r.description,
        })
    return out


def check_length(records: Iterable[SequenceRecord], *, warn_below: int = 100) -> list[str]:
    """给出长度相关的提醒（不改数据）。"""
    msgs = []
    for r in records:
        if r.length < warn_below:
            msgs.append(f"{r.name}: 仅 {r.length} bp，做脱靶分析时统计意义很弱")
        n = r.seq.count("N")
        if n:
            msgs.append(f"{r.name}: 含 {n} 个 N，含 N 的窗口会被跳过")
    return msgs


__all__ = ["IngestError", "SequenceError", "parse_fasta", "read_input", "summarize", "check_length"]
