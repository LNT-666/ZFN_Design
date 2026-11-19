"""脱靶分析（in silico）。

这是本方案里最容易被误读的一步，先把边界写清楚：

* 我们搜的是 **"序列上像 ZFN 位点的东西"**：左半位点 + 5-7 bp spacer +
  右半位点（右半位点以反向互补出现），每侧允许 ``max_mm`` 个错配；
* 这不等于"会在这里切割"。染色质状态、甲基化、DNA 可及性、局部结构
  都会改变实际切割；**结合脱靶 ≠ 切割脱靶**；
* 因此本模块的输出叫 ``burden``（负担/风险），不被解释为"脱靶切割概率"。

算法是暴力扫描 + 每侧错配上限，语义参照 ZFN-Site / Tagger
(Cradick et al. 2011, Nucleic Acids Res) 的"半位点 + spacer + 错配"模型，
但**没有**实现其完整的打分矩阵，属于简化启发式。
"""

from __future__ import annotations

from typing import Iterable, Sequence

from .models import OffTargetHit, OffTargetReport, SequenceRecord
from .sequtil import hamming, revcomp

DEFAULT_SPACERS = (5, 6, 7)


class OffTargetError(ValueError):
    pass


def search(
    records: Iterable[SequenceRecord],
    *,
    left_half: str,
    right_half: str,
    spacers: Sequence[int] = DEFAULT_SPACERS,
    max_mm_left: int = 2,
    max_mm_right: int = 2,
    exclude: Iterable[tuple[str, int, str]] = (),
    max_bp: int | None = 5_000_000,
) -> OffTargetReport:
    """在给定序列集合里找近似位点。

    参数
    ----
    left_half, right_half : 设计所指的两个 9 bp 半位点（各自 5'→3'）
    max_mm_left / right   : 每侧允许的错配数（0-3 之间才有意义）
    exclude               : 需要跳过的 on-target 坐标，元素为 ``(seq_name, start, strand)``
    """
    left_half = left_half.upper()
    right_half = right_half.upper()
    if len(left_half) != len(right_half):
        raise OffTargetError("左右半位点长度必须一致")
    hsl = len(left_half)
    right_on_top = revcomp(right_half)
    excluded = set(exclude)

    report = OffTargetReport(max_mm_left=max_mm_left, max_mm_right=max_mm_right,
                             spacers=tuple(spacers))
    total_bp = 0
    truncated = False

    for rec in records:
        total_bp += rec.length
        if max_bp is not None and total_bp > max_bp:
            truncated = True
            break
        for strand in ("+", "-"):
            seq = rec.seq if strand == "+" else revcomp(rec.seq)
            n = len(seq)
            for s in spacers:
                span = hsl + s + hsl
                if n < span:
                    continue
                for i in range(0, n - span + 1):
                    cand_left = seq[i : i + hsl]
                    if hamming(cand_left, left_half) > max_mm_left:
                        continue
                    cand_right = seq[i + hsl + s : i + span]
                    if hamming(cand_right, right_on_top) > max_mm_right:
                        continue
                    start_input = i + 1 if strand == "+" else rec.length - (i + span) + 1
                    key = (rec.name, start_input, strand)
                    if key in excluded:
                        continue
                    report.hits.append(
                        OffTargetHit(
                            seq_name=rec.name,
                            start=start_input,
                            strand=strand,
                            left_half=cand_left,
                            spacer=seq[i + hsl : i + hsl + s],
                            right_half=revcomp(cand_right),
                            mm_left=hamming(cand_left, left_half),
                            mm_right=hamming(cand_right, right_on_top),
                        )
                    )

    report.searched_bp = total_bp
    report.truncated = truncated
    report.hits.sort(key=lambda h: (h.total_mm, h.seq_name, h.start))
    return report


def search_single(
    records: Iterable[SequenceRecord],
    *,
    half_site: str,
    max_mm: int = 2,
    exclude: Iterable[tuple[str, int, str]] = (),
    max_bp: int | None = 5_000_000,
) -> OffTargetReport:
    """单体型脱靶：只找一个近似半位点，没有 spacer、没有搭档。

    与组合型的关键差别不是算法，而是**判据**：单体型只要一个 9 bp 半位点
    像就够了，所以同样的错配上限下命中会显著更多——这正是"三指单独不唯一"
    在数字上的体现。
    """
    half_site = half_site.upper()
    hsl = len(half_site)
    excluded = set(exclude)
    report = OffTargetReport(max_mm_left=max_mm, max_mm_right=0, spacers=())
    total_bp = 0
    truncated = False

    for rec in records:
        total_bp += rec.length
        if max_bp is not None and total_bp > max_bp:
            truncated = True
            break
        for strand in ("+", "-"):
            seq = rec.seq if strand == "+" else revcomp(rec.seq)
            n = len(seq)
            if n < hsl:
                continue
            for i in range(0, n - hsl + 1):
                cand = seq[i : i + hsl]
                mm = hamming(cand, half_site)
                if mm > max_mm:
                    continue
                start_input = i + 1 if strand == "+" else rec.length - (i + hsl) + 1
                key = (rec.name, start_input, strand)
                if key in excluded:
                    continue
                report.hits.append(
                    OffTargetHit(
                        seq_name=rec.name,
                        start=start_input,
                        strand=strand,
                        left_half=cand,
                        spacer="",
                        right_half="",
                        mm_left=mm,
                        mm_right=0,
                    )
                )

    report.searched_bp = total_bp
    report.truncated = truncated
    report.hits.sort(key=lambda h: (h.total_mm, h.seq_name, h.start))
    return report


def burden_to_score(burden: float, *, scale: float = 25.0) -> float:
    """把负担映射成 0-100 的"安全性"分项（越高越安全）。"""
    return max(0.0, 100.0 - scale * burden)


__all__ = ["OffTargetError", "search", "search_single", "burden_to_score", "DEFAULT_SPACERS"]
