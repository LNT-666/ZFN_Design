"""全基因组落点检索（S2）。

候选位点 `site_seq = left_half + spacer + revcomp(right_half)` 在参考基因组里的
所有出现位置，正负链都算。用来做两件事：

* **multi_copy 计数**：按全基因组落点计，不按输入序列里的重复计；
* **重复区判定**：落点是否压在 RepeatMasker 注释的重复区上。

实现是分块流式扫，每条染色体按 ``chunk`` 切开，一次扫完全部 motif：

* 装了 ``pyahocorasick`` → 建 AC 自动机，单趟扫完（ms01 实测 33 MB/s，
  3.1 Gb ≈ 94 s，**耗时与 motif 数无关**，正负链各入一条）；
* 没装 → 退回逐 motif 的 ``bytes.find``（聚合 ~700 MB/s，但耗时随 motif 数线性涨，
  几百个 motif 扫 3.1 Gb 要几十分钟），只适合小基因组。两条路结果等价。
"""

from __future__ import annotations

import os
from typing import Callable, Iterable, Iterator

__all__ = [
    "LandingError",
    "revcomp",
    "iter_fasta_chunks",
    "find_landings",
    "filter_sites_by_genome",
    "landing_stats",
]

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


class LandingError(RuntimeError):
    """落点检索失败（文件缺失、motif 非法等）。"""


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def _clean_motif(motif: str) -> str:
    s = (motif or "").strip().upper()
    if len(s) < 8:
        raise LandingError(f"motif 太短，无法在全基因组可靠定位：{motif!r}")
    if set(s) - set("ACGTN"):
        raise LandingError(f"motif 含非法碱基：{motif!r}")
    return s


def iter_fasta_chunks(
    fasta_path: str, chunk: int = 8_000_000
) -> Iterator[tuple[str, int, bytes]]:
    """按染色体分块产出 ``(chrom, start1, bytes)``，序列统一大写。"""
    if not os.path.exists(fasta_path):
        raise LandingError(f"找不到基因组 FASTA：{fasta_path}")
    if chunk < 1_000_000:
        raise LandingError("chunk 太小，至少 1 Mb")
    name: str | None = None
    buf = bytearray()
    pos = 0
    with open(fasta_path, "rb") as fh:
        for raw in fh:
            if raw.startswith(b">"):
                if name is not None and buf:
                    yield name, pos + 1, bytes(buf).upper()
                head = raw[1:].strip().decode("utf-8", "replace")
                name = head.split()[0] if head else "unnamed"
                pos = 0
                buf = bytearray()
                continue
            if name is None:
                continue
            line = raw.strip()
            if not line:
                continue
            buf += line
            if len(buf) >= chunk:
                yield name, pos + 1, bytes(buf).upper()
                pos += len(buf)
                buf = bytearray()
    if name is not None and buf:
        yield name, pos + 1, bytes(buf).upper()


def _build_automaton(pats: dict[str, tuple[bytes, bytes]]):
    """有 pyahocorasick 就建 AC 自动机（一次扫全基因组，耗时与 motif 数无关）。

    没装就返回 ``None``，退回逐 motif 的 ``bytes.find``：扫描量随 motif 数线性涨，
    549 个 motif 在 3.1 Gb 上要 40 分钟，所以只适合小基因组/少量候选。
    """
    try:
        import ahocorasick  # type: ignore
    except ImportError:
        return None
    a = ahocorasick.Automaton()
    for motif in pats:
        a.add_word(motif, (motif, "+"))
        a.add_word(revcomp(motif), (motif, "-"))
    a.make_automaton()
    return a


def find_landings(
    fasta_path: str,
    motifs: Iterable[str],
    *,
    chunk: int = 8_000_000,
    max_hits: int = 2000,
    progress: Callable[[dict], None] | None = None,
) -> dict[str, list[dict]]:
    """返回 ``{motif: [{chrom, start, end, strand}, ...]}``（1-based 闭区间）。

    命中数超过 ``max_hits`` 的 motif 只保留前 ``max_hits`` 条并标 ``copy_count_capped``
    ——重复家族动辄上万份，全存没意义。
    """
    pats: dict[str, tuple[bytes, bytes]] = {}
    for m in motifs:
        s = _clean_motif(m)
        pats[s] = (s.encode(), revcomp(s).encode())
    if not pats:
        return {}
    hits: dict[str, list[dict]] = {m: [] for m in pats}
    truncated: dict[str, bool] = {m: False for m in pats}
    automaton = _build_automaton(pats)
    scanned = 0
    for chrom, start1, data in iter_fasta_chunks(fasta_path, chunk):
        scanned += len(data)
        if automaton is not None:
            text = data.decode("ascii", "replace")
            for end, (motif, strand) in automaton.iter(text):
                bucket = hits[motif]
                if len(bucket) >= max_hits:
                    truncated[motif] = True
                    continue
                s0 = end - len(motif) + 1
                bucket.append(
                    {
                        "chrom": chrom,
                        "start": start1 + s0,
                        "end": start1 + end,
                        "strand": strand,
                    }
                )
        else:
            for motif, (fwd, rev) in pats.items():
                bucket = hits[motif]
                if len(bucket) >= max_hits:
                    truncated[motif] = True
                    continue
                i = data.find(fwd)
                while i >= 0 and len(bucket) < max_hits:
                    bucket.append(
                        {
                            "chrom": chrom,
                            "start": start1 + i,
                            "end": start1 + i + len(motif) - 1,
                            "strand": "+",
                        }
                    )
                    i = data.find(fwd, i + 1)
                j = data.find(rev)
                while j >= 0 and len(bucket) < max_hits:
                    bucket.append(
                        {
                            "chrom": chrom,
                            "start": start1 + j,
                            "end": start1 + j + len(motif) - 1,
                            "strand": "-",
                        }
                    )
                    j = data.find(rev, j + 1)
                if len(bucket) >= max_hits:
                    truncated[motif] = True
        if progress is not None:
            progress({"scanned_bp": scanned, "chrom": chrom})
    for motif, bucket in hits.items():
        if truncated[motif]:
            for rec in bucket:
                rec["copy_count_capped"] = True
    return hits


def engine_name() -> str:
    """当前用的扫描后端：``ahocorasick``（一次扫完）或 ``bytes.find``（逐 motif）。"""
    try:
        import ahocorasick  # type: ignore  # noqa: F401
    except ImportError:
        return "bytes.find"
    return "ahocorasick"


def landing_stats(landings: dict[str, list[dict]]) -> dict:
    """落点检索的汇总（写进 report.json，不用原始落点表）。"""
    return {
        "engine": engine_name(),
        "n_motifs": len(landings),
        "n_landings": sum(len(v) for v in landings.values()),
        "n_truncated_motifs": sum(1 for v in landings.values() if any(r.get("copy_count_capped") for r in v)),
        "motifs_multicopy": sum(1 for v in landings.values() if len(v) >= 3),
    }


def filter_sites_by_genome(
    sites,
    genome_fasta: str,
    *,
    repeats_path: str | None = None,
    min_copies: int = 3,
    chunk: int = 8_000_000,
    max_hits: int = 2000,
    progress: Callable[[dict], None] | None = None,
    repeats_index=None,
) -> tuple[list, dict]:
    """按全基因组落点过滤候选位点（宋成义 2026-11-18 口径）。

    判据（任一命中即剔除）：
      1. 落点压在 RepeatMasker 注释的重复区上（注释只覆盖已知重复）；
      2. 同一 motif 的全基因组落点数 ``>= min_copies``（兜注释漏掉的低复杂度序列）。

    返回 ``(kept_sites, stats)``；``stats`` 里按 motif 记下每条候选的判定结果，
    便于报告里逐条解释为什么被去掉。
    """
    if min_copies < 2:
        raise LandingError("min_copies 至少 2")
    motif_of: dict[str, list] = {}
    for s in sites:
        motif_of.setdefault(s.site_seq, []).append(s)
    landings = find_landings(
        genome_fasta, motif_of.keys(), chunk=chunk, max_hits=max_hits, progress=progress
    )

    ridx = repeats_index
    if repeats_path and ridx is None:
        from . import genome as _genome

        ridx = _genome.parse_repeatmasker(repeats_path)

    kept: list = []
    per_site: list[dict] = []
    n_multi = n_rep = 0
    for motif, group in motif_of.items():
        recs = landings.get(motif, [])
        copies = len(recs)
        multi = copies >= min_copies
        in_repeat = False
        repeat_at = None
        if ridx is not None:
            for rec in recs:
                if ridx.overlaps(rec["chrom"], rec["start"] - 1, rec["end"]):
                    in_repeat = True
                    repeat_at = rec
                    break
        verdict = "multi_copy" if multi else ("repeat" if in_repeat else "pass")
        if multi:
            n_multi += len(group)
        elif in_repeat:
            n_rep += len(group)
        else:
            kept.extend(group)
        for s in group:
            per_site.append(
                {
                    "site_id": getattr(s, "site_id", None),
                    "seq_name": s.seq_name,
                    "start": s.start,
                    "strand": s.strand,
                    "motif": motif,
                    "genome_copies": copies,
                    "copies_capped": any(r.get("copy_count_capped") for r in recs),
                    "verdict": verdict,
                    "landings": recs[:10],
                }
            )
    stats = {
        "genome_fasta": genome_fasta,
        "repeats_path": repeats_path,
        "engine": engine_name(),
        "min_copies": min_copies,
        "n_sites_in": sum(len(v) for v in motif_of.values()),
        "n_motifs": len(motif_of),
        "n_removed_multicopy": n_multi,
        "n_removed_repeat": n_rep,
        "n_kept": len(kept),
        "n_landings_total": sum(len(v) for v in landings.values()),
        "per_site": per_site,
    }
    return kept, stats
