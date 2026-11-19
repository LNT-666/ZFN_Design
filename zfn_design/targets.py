"""阶段 2：targetscan——在 DNA 里找"能被 ZFN 夹住"的位置。

ZFN 的几何约束
--------------
一个 ZFN 二聚体位点由三段组成，从 5' 到 3'（上链）::

    左半位点 (9 bp) | spacer (5-7 bp) | 右半位点 (9 bp，以反向互补形式出现)

为什么右半位点是"反向互补"：右单体是**倒过来**结合同一条双链的，
它的锌指阵列 5'→3' 方向与上链相反。所以如果我们把右单体"读到"的 9 bp
记作 ``right_half``(5'→3')，那么上链对应位置写出来就是
``revcomp(right_half)``。

本模块只负责"枚举"，不做任何过滤或打分——过滤在 designability，
打分在 score。这样每个阶段都能单独测。
"""

from __future__ import annotations

from typing import Iterable, Sequence

from .models import SequenceRecord, TargetSite
from .sequtil import revcomp

#: 半位点默认长度（3 指）。4 指阵列可改成 12，见 ``half_site_len`` 参数。
DEFAULT_HALF_SITE_LEN = 9

#: spacer 允许长度。6 最常用；5/7 也有活性报道。
DEFAULT_SPACERS = (5, 6, 7)


def scan_sequence(
    record: SequenceRecord,
    *,
    spacers: Sequence[int] = DEFAULT_SPACERS,
    half_site_len: int = DEFAULT_HALF_SITE_LEN,
    both_strands: bool = True,
    skip_n: bool = True,
    mode: str = "paired",
) -> list[TargetSite]:
    """枚举一条序列上的所有候选位点。

    ``mode="paired"``
        组合型：左半位点 + spacer + 右半位点，每个 spacer 长度各扫一遍。
    ``mode="single"``
        单体型：只找一个半位点（ZF-TF / 单 ZF 阵列），没有 spacer 与搭档，
        ``spacers`` 参数被忽略。
    """
    if mode not in ("paired", "single"):
        raise ValueError(f"mode 只能是 'paired' 或 'single'，收到 {mode!r}")
    sites: list[TargetSite] = []
    strands = ("+", "-") if both_strands else ("+",)
    for strand in strands:
        seq = record.seq if strand == "+" else revcomp(record.seq)
        if mode == "single":
            sites.extend(_scan_single(seq, record.name, strand, half_site_len,
                                     len(record.seq), skip_n))
        else:
            sites.extend(
                _scan_one(seq, record.name, strand, spacers, half_site_len,
                          len(record.seq), skip_n)
            )
    return sites


def _scan_single(seq: str, name: str, strand: str, hsl: int,
                 original_len: int, skip_n: bool) -> list[TargetSite]:
    """单体型：一个半位点就是一个位点（没有 spacer、没有右半位点）。"""
    out: list[TargetSite] = []
    n = len(seq)
    if n < hsl:
        return out
    for i in range(0, n - hsl + 1):
        window = seq[i : i + hsl]
        if skip_n and "N" in window:
            continue
        start_input = i + 1 if strand == "+" else original_len - (i + hsl) + 1
        out.append(
            TargetSite(
                seq_name=name,
                start=start_input,
                strand=strand,
                left_half=window,
                spacer="",
                right_half="",
                site_seq=window,
                kind="single",
            )
        )
    return out


def _scan_one(seq: str, name: str, strand: str, spacers: Sequence[int],
              hsl: int, original_len: int, skip_n: bool) -> list[TargetSite]:
    out: list[TargetSite] = []
    n = len(seq)
    for s in spacers:
        span = hsl + s + hsl
        if n < span:
            continue
        for i in range(0, n - span + 1):
            window = seq[i : i + span]
            if skip_n and "N" in window:
                continue
            left = seq[i : i + hsl]
            spacer = seq[i + hsl : i + hsl + s]
            right_on_top = seq[i + hsl + s : i + span]
            right = revcomp(right_on_top)

            # 换算回输入序列坐标（负链要从右往左数）
            start_input = i + 1 if strand == "+" else original_len - (i + span) + 1
            out.append(
                TargetSite(
                    seq_name=name,
                    start=start_input,
                    strand=strand,
                    left_half=left,
                    spacer=spacer,
                    right_half=right,
                    site_seq=window,
                )
            )
    return out


def scan_many(records: Iterable[SequenceRecord], **kwargs) -> list[TargetSite]:
    out: list[TargetSite] = []
    for r in records:
        out.extend(scan_sequence(r, **kwargs))
    return out


def expected_site_count(length: int, *, half_site_len: int = DEFAULT_HALF_SITE_LEN,
                        spacers: Sequence[int] = DEFAULT_SPACERS,
                        both_strands: bool = True,
                        mode: str = "paired") -> int:
    """理论位点数，用于自检（单元测试会拿它和实际枚举比对）。"""
    if mode == "single":
        per_strand = max(0, length - half_site_len + 1)
    else:
        per_strand = sum(max(0, length - (2 * half_site_len + s) + 1) for s in spacers)
    return per_strand * (2 if both_strands else 1)


def group_by_position(sites: Iterable[TargetSite]) -> dict[str, list[TargetSite]]:
    """按 ``seq_name`` 分组，方便报告分节。"""
    groups: dict[str, list[TargetSite]] = {}
    for s in sites:
        groups.setdefault(s.seq_name, []).append(s)
    return groups
