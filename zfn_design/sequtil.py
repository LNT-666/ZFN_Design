"""DNA 序列基础工具。

设计原则
--------
* 全部函数是纯函数、无全局状态，便于单元测试；
* 只做 DNA 层面的字符串运算，不引入生物信息学依赖（本机无 biopython）；
* 对非法输入抛 ``SequenceError``，绝不静默把 ``N`` 当作普通碱基丢进打分。
"""

from __future__ import annotations

from typing import Iterable, Iterator

#: 允许的碱基字符（大写）。IUPAC 简并码目前不支持，避免误判。
VALID_BASES = frozenset("ACGT")

#: 简并 / 模糊字符，出现即视为"不可设计"区域。
AMBIGUOUS = frozenset("NRYSWKMBDHV")

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


class SequenceError(ValueError):
    """序列非法（空、含非 DNA 字符等）。"""


def clean(raw: str) -> str:
    """去掉空白、数字、FASTA 注释痕迹，转大写。"""
    if raw is None:
        raise SequenceError("序列为 None")
    out = []
    for ch in raw.strip():
        if ch.isspace():
            continue
        out.append(ch.upper())
    seq = "".join(out)
    if not seq:
        raise SequenceError("空序列")
    return seq


def validate(seq: str, *, name: str = "<seq>", allow_n: bool = True) -> str:
    """校验并返回规范化后的大写序列。"""
    bad = sorted({c for c in seq if c not in VALID_BASES and not (allow_n and c == "N")})
    if bad:
        raise SequenceError(f"{name}: 含非法碱基 {bad!r}；只接受 A/C/G/T" + ("/N" if allow_n else ""))
    return seq


def revcomp(seq: str) -> str:
    """反向互补。"""
    return seq.translate(_COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    """GC 比例；序列为空时返回 0.0。"""
    if not seq:
        return 0.0
    gc = sum(1 for c in seq.upper() if c in "GC")
    return gc / len(seq)


def hamming(a: str, b: str) -> int:
    """等长序列的汉明距离；长度不等抛错。"""
    if len(a) != len(b):
        raise ValueError(f"hamming 需要等长序列: {len(a)} vs {len(b)}")
    return sum(1 for x, y in zip(a, b) if x != y)


def triplets(halfsite: str) -> tuple[str, str, str]:
    """把一个 9 bp 半位点切成 3 个三联体（按 5'→3' 顺序）。

    返回 ``(t_5p, t_mid, t_3p)``，注意这与指头编号 **不是** 同序，
    见 :func:`triplets_by_finger`。
    """
    if len(halfsite) != 9:
        raise SequenceError(f"半位点必须 9 bp，得到 {len(halfsite)}")
    return halfsite[0:3], halfsite[3:6], halfsite[6:9]


def triplets_by_finger(halfsite: str) -> dict[str, str]:
    """按 C2H2 阵列的反平行几何，把半位点分配给 F1..Fk（k = 长度/3）。

    锌指阵列在肽链上走 N→C，而 DNA 双链上走 5'→3'，两者**反平行**
    （Persikov 2015 NAR 原文：*"N-terminus-poolGGG-poolCCC-poolAAA-C-terminus
    (zinc fingers bind DNA anti-parallel to the 5'-3' sequence of DNA)"*）。

    因此对一个 5'-ABC-DEF-GHI-3' 的 9 bp 半位点：

    ==========  ==========  ==============================
    指头        在阵列中    识别的三联体
    ==========  ==========  ==============================
    F1          N 端        GHI（3' 端那个）
    F2          中间        DEF
    F3          C 端        ABC（5' 端那个）
    ==========  ==========  ==============================

    18 bp 的六指半位点同样处理：F1–F3 读 3' 那 9 bp，F4–F6 读 5' 那 9 bp。
    写成 DNA 蛋白序列时，阵列顺序应是 ``[F1, F2, ..., Fk]``。
    """
    n = len(halfsite)
    if n % 3 or n < 3:
        raise SequenceError(f"半位点长度必须是 3 的倍数，得到 {n}")
    blocks = [halfsite[i : i + 3] for i in range(0, n, 3)]  # 5'→3'
    k = len(blocks)
    return {f"F{i + 1}": blocks[k - 1 - i] for i in range(k)}


def finger_order(halfsite: str) -> list[str]:
    """返回 ``[F1, F2, ...]`` 各自的三联体，顺序即蛋白 N→C。"""
    m = triplets_by_finger(halfsite)
    return [m[k] for k in sorted(m, key=lambda s: int(s[1:]))]


def max_homopolymer(seq: str) -> int:
    """最长同碱基连续长度。"""
    if not seq:
        return 0
    best = run = 1
    for prev, cur in zip(seq, seq[1:]):
        run = run + 1 if cur == prev else 1
        best = max(best, run)
    return best


def tandem_repeat_count(seq: str, unit: int = 2) -> int:
    """统计长度为 ``unit`` 的串联重复出现次数（用于低复杂度判定）。"""
    if len(seq) < unit * 2:
        return 0
    counts: dict[str, int] = {}
    for i in range(len(seq) - unit + 1):
        k = seq[i : i + unit]
        counts[k] = counts.get(k, 0) + 1
    return sum(1 for v in counts.values() if v >= 3)


def is_self_complementary(seq: str, min_arm: int = 6) -> bool:
    """检测位点内部是否存在长度 >= ``min_arm`` 的反向重复（可形成发夹/自身二聚）。"""
    n = len(seq)
    for arm in range(min_arm, n // 2 + 1):
        for i in range(0, n - 2 * arm + 1):
            left = seq[i : i + arm]
            right = seq[i + arm : i + 2 * arm]
            if left == revcomp(right):
                return True
    return False


def iter_windows(seq: str, width: int) -> Iterator[tuple[int, str]]:
    """滑窗，返回 ``(0-based 起点, 窗口序列)``。"""
    for i in range(0, len(seq) - width + 1):
        yield i, seq[i : i + width]


def wrap(seq: str, width: int = 60) -> str:
    """按固定宽度折行，便于在报告里显示长序列。"""
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def dedupe(items: Iterable[str]) -> list[str]:
    """保序去重。"""
    seen: set[str] = set()
    out: list[str] = []
    for it in items:
        if it not in seen:
            seen.add(it)
            out.append(it)
    return out
