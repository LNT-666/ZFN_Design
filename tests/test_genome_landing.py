"""S2：全基因组落点过滤的测试（宋成义 2026-11-18 口径）。

覆盖三件事：
1. 同一 motif 在基因组里出现 >=3 次 → 剔除（multi_copy）；
2. 落点压在 RepeatMasker 注释区 → 剔除；
3. 只出现一次且不在重复区 → 通过，并且报告 meta 变成 S2 / filtered_target。
"""

from __future__ import annotations

import random

from zfn_design import genome_landing, pipeline
from zfn_design.library import load_default_library
from zfn_design.models import SequenceRecord, TargetSite

RNG = random.Random(7)


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _random_dna(n: int) -> str:
    return "".join(RNG.choice("ACGT") for _ in range(n))


def _write_fasta(path, chrom: str, seq: str) -> None:
    path.write_text(f">{chrom}\n" + "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60)) + "\n",
                    encoding="utf-8")


def _scaffold():
    """建一个 90 kb 的假基因组：motif A 出现 3 次，motif B 在重复区出现 1 次。"""
    filler = list(_random_dna(90_000))
    motif_a = "GGT" + "ACG" + "TGC" + "ATG" + "CGT"  # 15 bp
    motif_b = "CCT" + "AGA" + "TCA" + "GTC" + "AGC"
    placed = {}
    for name, motif, positions in (("a", motif_a, (1000, 30_000, 60_000)), ("b", motif_b, (80_000,))):
        for p in positions:
            filler[p : p + len(motif)] = list(motif)
        placed[name] = motif
    seq = "".join(filler)
    return seq, placed


def _site(motif: str) -> TargetSite:
    half = (len(motif) - 6) // 2
    return TargetSite(
        seq_name="input",
        start=1,
        strand="+",
        left_half=motif[:half],
        spacer=motif[half : half + 6],
        right_half=revcomp(motif[half + 6 :]),
        site_seq=motif,
    )


def test_find_landings_counts_both_strands(tmp_path):
    seq, placed = _scaffold()
    fasta = tmp_path / "genome.fa"
    _write_fasta(fasta, "chr1", seq)
    hits = genome_landing.find_landings(str(fasta), [placed["a"], placed["b"]], chunk=1_000_000)
    assert len(hits[placed["a"]]) == 3
    assert len(hits[placed["b"]]) == 1


def test_filter_removes_multicopy_and_repeat(tmp_path):
    seq, placed = _scaffold()
    fasta = tmp_path / "genome.fa"
    _write_fasta(fasta, "chr1", seq)
    bed = tmp_path / "rep.bed"
    bed.write_text("chr1\t79900\t80100\n", encoding="utf-8")

    sites = [_site(placed["a"]), _site(placed["b"])]
    kept, stats = genome_landing.filter_sites_by_genome(
        sites, str(fasta), repeats_path=str(bed), min_copies=3, chunk=1_000_000
    )
    assert kept == [], "两条都该被剔除：A 是多拷贝，B 落在重复区"
    assert stats["n_removed_multicopy"] == 1
    assert stats["n_removed_repeat"] == 1
    verdicts = {rec["motif"]: rec["verdict"] for rec in stats["per_site"]}
    assert verdicts == {placed["a"]: "multi_copy", placed["b"]: "repeat"}


def test_pipeline_run_with_genome_marks_stage_s2(tmp_path, demo_record):
    # 把 demo 序列整条塞进基因组，候选位点的落点都能找到
    seq = demo_record.seq
    fasta = tmp_path / "g.fa"
    _write_fasta(fasta, "chr1", "ACGT" * 2000 + seq + "ACGT" * 2000)
    payload = pipeline.run(
        [SequenceRecord(name="input", seq=seq)],
        library=load_default_library(),
        spacers=(6,),
        max_sites=3,
        n_designs=1,
        fusion_keys=[],
        genome_fasta=str(fasta),
        command="pytest",
    )
    assert payload["meta"]["stage"] == "S2"
    assert payload["meta"]["report_kind"] == "filtered_target"
    assert payload["targets_summary"]["genome_filter"] is not None
