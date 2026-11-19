"""全基因组输入与重复序列过滤——targetscan 的两块"入口"。

为什么单独一个模块
------------------
真核基因组不能像质粒那样整条读进内存再枚举：hg38 是 3.1 Gbp，
按 9-spacer-9、双链、三个 spacer 全扫一遍会产生 **约 186 亿** 个候选窗口。
所以这一层必须做到三件事：

1. **流式**读 FASTA（一次只驻留一条染色体，约 250 MB 峰值）；
2. **生成器**出结果（``scan_genome`` 是 generator，不物化列表）；
3. **先过滤后产出**（重复区、目标区间在被 ``yield`` 之前就问掉）。

过滤优先于枚举的收益是数量级的：RepeatMasker 通常覆盖人类基因组 50% 上下，
先去重复，后面的 designability 与打分就少做一半的活。

数据格式
--------
- **RepeatMasker ``.out``**：SW score / perc div / perc del / perc ins /
  query / begin / end / (left) / strand / ID / class-family / repeat…
  其中 ``begin``/``end`` 是 **1-based 闭区间**，这里统一转成 0-based 半开。
- **RepeatMasker ``.bed``**：chrom / start / end / name / score / strand，
  本来就是 0-based 半开。
- 两种都认，按扩展名判断，扩展名不认识就嗅探首行。
"""

from __future__ import annotations

import bisect
import os
import sys
from dataclasses import dataclass, field
from typing import Callable, Iterator, Sequence

from .models import SequenceRecord, TargetSite
from .sequtil import revcomp

#: 一次读多少 bp 做扫描（含跨块重叠）。4 Mb 时单块扫描约 1–2 秒。
DEFAULT_CHUNK = 4_000_000


class GenomeError(RuntimeError):
    """基因组 / 注释文件读不了、格式不对。"""


# --------------------------------------------------------------------------
# FASTA 流式读取
# --------------------------------------------------------------------------

def iter_fasta(path: str, *, upper: bool = True) -> Iterator[tuple[str, str]]:
    """逐条产出 ``(染色体名, 序列)``；一次只驻留一条序列。

    名字取 ``>`` 之后第一段（空格/制表符前），与 faidx 习惯一致。
    空行、行内空白都忽略。
    """
    if not os.path.exists(path):
        raise GenomeError(f"找不到 FASTA：{path}")
    name: str | None = None
    buf: list[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf).upper() if upper else "".join(buf)
                head = line[1:].strip()
                name = head.split()[0] if head else "unnamed"
                buf = []
            else:
                if name is None:
                    continue
                buf.append(line.strip())
    if name is not None:
        yield name, "".join(buf).upper() if upper else "".join(buf)


def fasta_index(path: str) -> list[tuple[str, int]]:
    """只读 ``>`` 行，快速拿到每条染色体的长度（不加载序列）。"""
    out: list[tuple[str, int]] = []
    name: str | None = None
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, n))
                head = line[1:].strip()
                name = head.split()[0] if head else "unnamed"
                n = 0
            elif name is not None:
                n += len(line.strip())
    if name is not None:
        out.append((name, n))
    return out


# --------------------------------------------------------------------------
# 区间索引（重复区 / 目标区）
# --------------------------------------------------------------------------

def _merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals.sort()
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    return [(s, e) for s, e in merged]


@dataclass
class IntervalIndex:
    """按染色体合并好的半开区间集合，``overlaps`` 走 bisect，O(log n)。"""

    spans: dict[str, list[tuple[int, int]]] = field(default_factory=dict)
    total_bp: dict[str, int] = field(default_factory=dict)
    n_raw: int = 0

    def add(self, chrom: str, start0: int, end0: int) -> None:
        if end0 <= start0:
            return
        self.spans.setdefault(chrom, []).append((start0, end0))
        self.n_raw += 1

    def finalize(self) -> "IntervalIndex":
        for chrom, iv in list(self.spans.items()):
            merged = _merge(iv)
            self.spans[chrom] = merged
            self.total_bp[chrom] = sum(e - s for s, e in merged)
        return self

    # -- 查询 --------------------------------------------------------------
    def overlaps(self, chrom: str, start0: int, end0: int) -> bool:
        iv = self.spans.get(chrom)
        if not iv:
            return False
        i = bisect.bisect_left(iv, (start0, -1))
        # 前一个区间可能罩住起点
        if i > 0 and iv[i - 1][1] > start0:
            return True
        return i < len(iv) and iv[i][0] < end0

    def contains(self, chrom: str, start0: int, end0: int) -> bool:
        """严格被某个区间包含（用于「只扫这些区域」）。"""
        iv = self.spans.get(chrom)
        if not iv:
            return False
        i = bisect.bisect_right(iv, (start0, 1 << 62))
        if i == 0:
            return False
        s, e = iv[i - 1]
        return s <= start0 and end0 <= e

    def spans_of(self, chrom: str) -> list[tuple[int, int]]:
        return self.spans.get(chrom, [])

    def chroms(self) -> list[str]:
        return sorted(self.spans)

    def summary(self) -> dict:
        return {
            "n_intervals_raw": self.n_raw,
            "n_intervals_merged": sum(len(v) for v in self.spans.values()),
            "n_chroms": len(self.spans),
            "total_bp": sum(self.total_bp.values()),
        }


def parse_repeatmasker(path: str) -> IntervalIndex:
    """读 RepeatMasker ``.out`` 或 ``.bed``，返回区间索引（自动判格式）。"""
    if not os.path.exists(path):
        raise GenomeError(f"找不到重复序列注释文件：{path}")
    idx = IntervalIndex()
    ext = os.path.splitext(path)[1].lower()
    fmt = "out" if ext in (".out", ".rm", ".align") else ("bed" if ext == ".bed" else None)
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if fmt == "out" and not (len(parts) >= 7 and parts[5].lstrip("+-").isdigit() and parts[6].isdigit()):
                continue
            if fmt is None:  # 嗅探：.out 的第 5、6 列是 query/begin
                is_out = len(parts) >= 15 and parts[8] in ("+", "C", "-")
                is_bed = len(parts) >= 3 and parts[1].lstrip("+-").isdigit()
                if not is_out and not is_bed:
                    continue
                fmt = "out" if is_out else "bed"
            try:
                if fmt == "out":
                    chrom = parts[4]
                    start0 = int(parts[5]) - 1
                    end0 = int(parts[6])
                else:
                    chrom, start0, end0 = parts[0], int(parts[1]), int(parts[2])
            except (IndexError, ValueError) as exc:
                raise GenomeError(f"{path}:{lineno} 看不懂这一行：{line.strip()[:80]}") from exc
            idx.add(chrom, start0, end0)
    if not idx.spans:
        raise GenomeError(f"{path} 里没解析出任何区间（格式对不对？）")
    return idx.finalize()


def parse_bed(path: str) -> IntervalIndex:
    """读普通 BED（0-based 半开），用于「只在启动子/外显子里扫」。"""
    if not os.path.exists(path):
        raise GenomeError(f"找不到 BED：{path}")
    idx = IntervalIndex()
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise GenomeError(f"{path}:{lineno} BED 至少要 3 列")
            try:
                idx.add(parts[0], int(parts[1]), int(parts[2]))
            except ValueError as exc:
                raise GenomeError(f"{path}:{lineno} 坐标不是整数") from exc
    if not idx.spans:
        raise GenomeError(f"{path} 里没解析出任何区间")
    return idx.finalize()


# --------------------------------------------------------------------------
# 流式扫描
# --------------------------------------------------------------------------

def _windows_for(spacers: Sequence[int], hsl: int, mode: str) -> list[int]:
    if mode == "single":
        return [0]
    return [s for s in spacers]


def scan_genome(
    fasta_path: str,
    *,
    spacers: Sequence[int] = (5, 6, 7),
    half_site_len: int = 9,
    both_strands: bool = True,
    mode: str = "paired",
    repeats: IntervalIndex | None = None,
    regions: IntervalIndex | None = None,
    limit: int | None = None,
    chunk: int = DEFAULT_CHUNK,
    progress: Callable[[str], None] | None = None,
    stats: dict | None = None,
) -> Iterator[TargetSite]:
    """在全基因组上流式枚举候选位点。

    生成器：调用方可以边拿边写文件，不会在内存里堆出十亿个对象。
    过滤顺序是「区间 → 窗口 → 重复」，重复区的窗口根本不会构造 ``TargetSite``。
    """
    if mode not in ("paired", "single"):
        raise ValueError(f"mode 只能是 'paired' 或 'single'，收到 {mode!r}")
    spacers = tuple(sorted({int(s) for s in spacers})) if mode == "paired" else (0,)
    spans = _windows_for(spacers, half_site_len, mode)
    longest = max(half_site_len * 2 + s for s in spans) if mode == "paired" else half_site_len

    emitted = 0
    strands = ("+", "-") if both_strands else ("+",)
    if stats is not None:
        for k in ("windows", "skipped_n", "skipped_repeat", "emitted"):
            stats.setdefault(k, 0)

    for chrom, seq in iter_fasta(fasta_path):
        L = len(seq)
        if L < longest:
            continue
        targets = regions.spans_of(chrom) if regions is not None else [(0, L)]
        if not targets:
            continue
        if progress:
            progress(f"{chrom}  {L:,} bp")

        for rs, re_ in targets:
            rs = max(0, rs)
            re_ = min(L, re_)
            block_start = rs
            while block_start < re_:
                core_end = min(re_, block_start + chunk)
                # 向右多读 longest-1，让起点落在块内的窗口都完整
                ext_end = min(L, core_end + longest - 1)
                block = seq[block_start:ext_end]
                nb = len(block)
                n_core = core_end - block_start
                for strand in strands:
                    sub = block if strand == "+" else revcomp(block)
                    for s in spans:
                        span = half_site_len if mode == "single" else half_site_len * 2 + s
                        if strand == "+":
                            # 起点在 core 内的窗口：[0, n_core-1]，且窗口要能装进 block
                            i_lo, i_hi = 0, min(n_core - 1, nb - span)
                        else:
                            # 负链上 sub 的第 i 个窗口对应正链起点 nb - i - span，
                            # 要让它落在 core 内：[nb - span - n_core + 1, nb - span]
                            i_lo = max(0, nb - span - n_core + 1)
                            i_hi = nb - span
                        if i_hi < i_lo:
                            continue
                        for i in range(i_lo, i_hi + 1):
                            if stats is not None:
                                stats["windows"] = stats.get("windows", 0) + 1
                            if "N" in sub[i:i + span]:
                                if stats is not None:
                                    stats["skipped_n"] = stats.get("skipped_n", 0) + 1
                                continue
                            if strand == "+":
                                p_start0 = block_start + i
                            else:
                                p_start0 = block_start + nb - i - span
                            if repeats is not None and repeats.overlaps(chrom, p_start0, p_start0 + span):
                                if stats is not None:
                                    stats["skipped_repeat"] = stats.get("skipped_repeat", 0) + 1
                                continue
                            start_input = p_start0 + 1
                            if mode == "single":
                                left = sub[i:i + half_site_len]
                                site = TargetSite(
                                    seq_name=chrom, start=start_input, strand=strand,
                                    left_half=left, spacer="", right_half="",
                                    site_seq=left, kind="single",
                                )
                            else:
                                left = sub[i:i + half_site_len]
                                sp = sub[i + half_site_len:i + half_site_len + s]
                                right = revcomp(sub[i + half_site_len + s:i + span])
                                site = TargetSite(
                                    seq_name=chrom, start=start_input, strand=strand,
                                    left_half=left, spacer=sp, right_half=right,
                                    site_seq=sub[i:i + span],
                                )
                            yield site
                            emitted += 1
                            if stats is not None:
                                stats["emitted"] = emitted
                            if limit is not None and emitted >= limit:
                                return
                block_start = core_end
        if progress:
            progress(f"{chrom} 完成，累计 {emitted:,} 个候选")


def estimate_sites(
    lengths: Sequence[int] | list[tuple[str, int]],
    *,
    spacers: Sequence[int] = (5, 6, 7),
    half_site_len: int = 9,
    both_strands: bool = True,
    mode: str = "paired",
    repeat_fraction: float = 0.0,
) -> int:
    """粗估全基因组候选数（用来先给用户一个数量级，避免跑到一半才发现太多）。"""
    total = 0
    for item in lengths:
        L = item[1] if isinstance(item, tuple) else int(item)
        if mode == "single":
            per = max(0, L - half_site_len + 1)
        else:
            per = sum(max(0, L - (2 * half_site_len + s) + 1) for s in spacers)
        total += per * (2 if both_strands else 1)
    if repeat_fraction:
        total = int(total * (1.0 - max(0.0, min(1.0, repeat_fraction))))
    return total


def describe_speed(n_bp: int, seconds: float) -> str:
    if seconds <= 0:
        return "—"
    rate = n_bp / seconds
    return f"{rate/1e6:.2f} Mb/s（{seconds:.1f} s / {n_bp/1e6:.1f} Mb）"


def stderr_progress(msg: str) -> None:
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()
