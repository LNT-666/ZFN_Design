"""全基因组入口 + 重复序列过滤（targetscan 的过滤器）。"""

from __future__ import annotations

import pytest

from zfn_design import genome, targets
from zfn_design.models import SequenceRecord

# 两条"染色体"：一条 60 bp、一条 40 bp
CHR1 = "TACGTAGCTAGCTTAGCTAGGCTAGCGTGGGCGGATATCCGCCCACGCTTACGGATCCTTAG"
CHR2 = "GGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGC"


@pytest.fixture()
def fasta_file(tmp_path):
    p = tmp_path / "mini.fa"
    p.write_text(
        ">chr1 测试第一条\n"
        + "\n".join(CHR1[i:i + 50] for i in range(0, len(CHR1), 50))
        + "\n>chr2\n"
        + CHR2
        + "\n",
        encoding="utf-8",
    )
    return str(p)


def _sites_key(sites):
    return sorted((s.seq_name, s.start, s.strand, s.site_seq, s.spacer) for s in sites)


# ---------------------------------------------------------------- FASTA

def test_iter_fasta_names_and_lengths(fasta_file):
    got = list(genome.iter_fasta(fasta_file))
    assert [n for n, _ in got] == ["chr1", "chr2"]
    assert got[0][1] == CHR1 and got[1][1] == CHR2


def test_fasta_index_matches_sequences(fasta_file):
    assert genome.fasta_index(fasta_file) == [("chr1", len(CHR1)), ("chr2", len(CHR2))]


def test_iter_fasta_missing_file_raises():
    with pytest.raises(genome.GenomeError):
        list(genome.iter_fasta("nope.fa"))


# ------------------------------------------------- 区间索引 / RepeatMasker

def test_repeatmasker_out_and_bed_agree(tmp_path):
    out = tmp_path / "r.out"
    # .out 的 begin/end 是 1-based 闭区间
    out.write_text(
        "  100 1.2 0.3 0.1 chr1 10 20 (0) + 1 ALU SINE/Alu\n"
        "  200 2.0 0.5 0.2 chr1 30 35 (0) C 2 LINE/L1\n",
        encoding="utf-8",
    )
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t9\t20\tAlu\t0\t+\nchr1\t29\t35\tL1\t0\t-\n", encoding="utf-8")

    a = genome.parse_repeatmasker(str(out))
    b = genome.parse_repeatmasker(str(bed))
    assert a.spans_of("chr1") == b.spans_of("chr1") == [(9, 20), (29, 35)]
    assert a.summary()["total_bp"] == b.summary()["total_bp"] == 11 + 6


def test_interval_overlaps_and_contains(fasta_file, tmp_path):
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t10\t20\nchr1\t30\t40\n", encoding="utf-8")
    idx = genome.parse_repeatmasker(str(bed))

    assert idx.overlaps("chr1", 19, 25)      # 搭到 20 的边
    assert not idx.overlaps("chr1", 20, 30)  # 两段之间的空档
    assert idx.overlaps("chr1", 5, 11)       # 只搭到开头
    assert not idx.overlaps("chr2", 0, 99)   # 别的染色体
    assert idx.contains("chr1", 12, 18)
    assert not idx.contains("chr1", 18, 22)


def test_merge_overlapping_intervals(tmp_path):
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t0\t10\nchr1\t8\t20\nchr1\t15\t25\n", encoding="utf-8")
    idx = genome.parse_repeatmasker(str(bed))
    assert idx.spans_of("chr1") == [(0, 25)]


def test_parse_bed_rejects_garbage(tmp_path):
    p = tmp_path / "x.bed"
    p.write_text("chr1\tnot_a_number\t20\n", encoding="utf-8")
    with pytest.raises(genome.GenomeError):
        genome.parse_bed(str(p))


# ------------------------------------------------------------- 流式扫描

def test_scan_genome_equals_scan_sequence(fasta_file):
    """分块扫描必须与整条扫描逐位一致——包括跨块的窗口和负链坐标。"""
    got = _sites_key(genome.scan_genome(fasta_file, spacers=(6,), chunk=17))
    want = []
    for name, seq in genome.iter_fasta(fasta_file):
        rec = SequenceRecord(name=name, seq=seq, description="", source="test")
        want.extend(targets.scan_sequence(rec, spacers=(6,)))
    assert got == _sites_key(want)
    assert got, "应该有候选位点"


def test_scan_genome_single_mode(fasta_file):
    got = _sites_key(genome.scan_genome(fasta_file, mode="single", chunk=13))
    want = []
    for name, seq in genome.iter_fasta(fasta_file):
        rec = SequenceRecord(name=name, seq=seq, description="", source="test")
        want.extend(targets.scan_sequence(rec, mode="single"))
    assert got == _sites_key(want)


def test_scan_genome_drops_repeat_sites(fasta_file, tmp_path):
    bed = tmp_path / "r.bed"
    bed.write_text("chr1\t0\t30\n", encoding="utf-8")
    idx = genome.parse_repeatmasker(str(bed))
    filtered = list(genome.scan_genome(fasta_file, spacers=(6,), repeats=idx, chunk=11))
    assert filtered
    assert all(s.seq_name != "chr1" or s.start > 30 for s in filtered)


def test_scan_genome_regions_restrict(fasta_file, tmp_path):
    bed = tmp_path / "t.bed"
    bed.write_text("chr1\t0\t25\n", encoding="utf-8")
    reg = genome.parse_bed(str(bed))
    got = list(genome.scan_genome(fasta_file, spacers=(6,), regions=reg))
    assert got
    assert all(s.seq_name == "chr1" and s.start <= 25 for s in got)


def test_scan_genome_limit(fasta_file):
    got = list(genome.scan_genome(fasta_file, spacers=(5, 6, 7), limit=5))
    assert len(got) == 5


def test_estimate_sites_matches_actual(fasta_file):
    want = len(list(genome.scan_genome(fasta_file, spacers=(6,), mode="single")))
    got = genome.estimate_sites([(n, len(s)) for n, s in genome.iter_fasta(fasta_file)],
                                mode="single")
    assert got == want


def test_scan_genome_rejects_bad_mode(fasta_file):
    with pytest.raises(ValueError):
        list(genome.scan_genome(fasta_file, mode="triple"))
