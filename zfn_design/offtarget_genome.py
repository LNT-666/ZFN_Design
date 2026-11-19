"""全基因组尺度的脱靶扫描。

和 :mod:`zfn_design.offtarget` 的区别：那边是在一小段输入序列上暴力枚举，
这边是**流式扫整个基因组 FASTA**，逐个位置算左右半位点的错配数。

为什么必须另写一个模块
----------------------
1. 内存：人类基因组 3 Gbp 不能一次读进内存，必须按块流转。
2. 速度：逐位置 Python 循环约 2 µs/位点，3 Gbp 要一个多小时；
   装了 numpy 时走向量化路径（每个 4 Mb 块几十毫秒），没有 numpy 就退回纯 Python
   并明确警告"慢，建议先用 --regions 缩小范围"。
3. 语义：这里的 burden 是**在用户提供的整个基因组上**数出来的近似位点，
   仍然是"序列相似度"，不是切割概率（见文档 §5）。

坐标约定
--------
``start0`` 是 **0-based**、正链上的窗口起点；一条 ZFN 位点在负链上出现时，
用 (revcomp(right), spacer, revcomp(left)) 这样的"交换模式"在同一正链序列上命中，
所以 ``strand`` 记录的是**该位点相对于输入序列的方向**，不是基因组的正负链。
"""

from __future__ import annotations

import time
from typing import Any, Callable, Iterator

from . import genome as _genome
from .sequtil import revcomp

# 错配 → 权重。总和最多 4（每侧 ≤2），所以表要覆盖到 4；
# 3/4 错配的权重只用于排序，不代表"可能被切"。
MM_WEIGHT = {0: 1.0, 1: 0.25, 2: 0.06, 3: 0.015, 4: 0.004}

try:  # 可选加速：没有 numpy 也能跑
    import numpy as _np
except Exception:  # pragma: no cover - 取决于安装环境
    _np = None


HAS_NUMPY = _np is not None


def _encode(s: str) -> bytes:
    return s.encode("ascii", "ignore")


def _scan_block_numpy(block: bytes | bytearray, pat_left: bytes, pat_right: bytes,
                      gap: int, max_mm_l: int, max_mm_r: int):
    """返回 (起始下标数组, mm_left 数组, mm_right 数组)。"""
    arr = _np.frombuffer(block, dtype=_np.uint8)
    n = arr.size
    span = len(pat_left) + gap + len(pat_right)
    if n < span:
        return _np.empty(0, dtype=_np.int64), _np.empty(0, dtype=_np.int8), _np.empty(0, dtype=_np.int8)
    m = n - span + 1

    mm_l = _np.zeros(m, dtype=_np.int8)
    for k, ch in enumerate(pat_left):
        mm_l += (arr[k:k + m] != ch)
    mm_r = _np.zeros(m, dtype=_np.int8)
    off = len(pat_left) + gap
    for k, ch in enumerate(pat_right):
        mm_r += (arr[off + k:off + k + m] != ch)

    keep = (mm_l <= max_mm_l) & (mm_r <= max_mm_r)
    idx = _np.flatnonzero(keep)
    return idx, mm_l[idx], mm_r[idx]


def _scan_block_py(block: bytes, pat_left: bytes, pat_right: bytes, gap: int,
                   max_mm_l: int, max_mm_r: int):
    span = len(pat_left) + gap + len(pat_right)
    out = []
    for i in range(0, len(block) - span + 1):
        w_l = block[i:i + len(pat_left)]
        if sum(a != b for a, b in zip(w_l, pat_left)) > max_mm_l:
            continue
        off = i + len(pat_left) + gap
        w_r = block[off:off + len(pat_right)]
        if sum(a != b for a, b in zip(w_r, pat_right)) > max_mm_r:
            continue
        out.append((i,
                    sum(a != b for a, b in zip(w_l, pat_left)),
                    sum(a != b for a, b in zip(w_r, pat_right))))
    return out


def genome_search(
    fasta_path: str,
    left_half: str,
    right_half: str,
    *,
    spacers: tuple[int, ...] = (6,),
    max_mm_left: int = 2,
    max_mm_right: int = 2,
    max_hits: int = 500,
    chunk: int = 4_000_000,
    progress: Callable[[str], None] | None = None,
    on_bp: Callable[[int, int], None] | None = None,
    allow_numpy: bool = True,
) -> dict[str, Any]:
    """在整个基因组上找某个 ZFN 的近似位点。

    返回 ``hits``（最多 ``max_hits`` 条）、``burden``（权重和）、``by_mm``
    （0/1/2 错配计数）、``scanned_bp``、``truncated``、``engine``（numpy / python）。
    """
    use_np = bool(allow_numpy and _np is not None)
    left = left_half.upper()
    right = right_half.upper()

    # 两条方向：位点在正链 / 位点在负链（等价于交换左右半位点再取反向互补）。
    # 左右半位点相同时两种模式完全一样，只扫一遍，否则同一个物理位点会被数两次。
    orientations = [
        ("+", _encode(left), _encode(revcomp(right))),
        ("-", _encode(right), _encode(revcomp(left))),
    ]
    if (orientations[0][1], orientations[0][2]) == (orientations[1][1], orientations[1][2]):
        orientations = orientations[:1]

    hits: list[dict[str, Any]] = []
    by_mm: dict[int, int] = {}
    scanned = 0
    truncated = False
    t0 = time.perf_counter()
    _total_bp = 0
    if on_bp is not None:
        try:
            _total_bp = sum(n for _c, n in _genome.fasta_index(fasta_path))
        except OSError:
            _total_bp = 0

    for chrom, seq in _genome.iter_fasta(fasta_path):
        data = seq.encode("ascii", "ignore")
        n = len(data)
        scanned += n
        for gap in spacers:
            for strand, pl, pr in orientations:
                span = len(pl) + gap + len(pr)
                if n < span:
                    continue
                for start in range(0, max(1, n - span + 1), chunk):
                    block = data[start:start + chunk + span - 1]
                    if use_np:
                        idx, ml, mr = _scan_block_numpy(block, pl, pr, gap,
                                                        max_mm_left, max_mm_right)
                        pairs = zip(idx.tolist(), ml.tolist(), mr.tolist())
                    else:
                        pairs = _scan_block_py(block, pl, pr, gap,
                                               max_mm_left, max_mm_right)
                    for i, a, b in pairs:
                        mm = int(a) + int(b)
                        by_mm[mm] = by_mm.get(mm, 0) + 1
                        if len(hits) < max_hits:
                            hits.append({
                                "chrom": chrom,
                                "start0": start + i,
                                "strand": strand,
                                "spacer": gap,
                                "mm_left": int(a),
                                "mm_right": int(b),
                                "mm_total": mm,
                                "weight": MM_WEIGHT.get(mm, 0.0),
                                "site": (left if strand == "+" else right),
                            })
                        else:
                            truncated = True
                    if progress and start and start % (chunk * 25) == 0:
                        progress(f"  {chrom} +{start:,} / {n:,} · {len(hits)} hits")
        if on_bp is not None:
            on_bp(scanned, _total_bp)

    burden = sum(h["weight"] for h in hits)
    elapsed = time.perf_counter() - t0
    return {
        "hits": hits,
        "n_hits": len(hits),
        "by_mm": by_mm,
        "n_total_hits": sum(by_mm.values()),
        "burden": round(burden, 4),
        "scanned_bp": scanned,
        "elapsed_s": round(elapsed, 2),
        "truncated": truncated,
        "engine": "numpy" if use_np else "python",
        "max_mm_left": max_mm_left,
        "max_mm_right": max_mm_right,
        "spacers": list(spacers),
    }


def score_from_burden(burden: float, scale: float = 25.0) -> float:
    """把 burden 折成 0–100 的分项分（和高通量那套一致）。"""
    return max(0.0, 100.0 - scale * burden)
