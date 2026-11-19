"""阶段 S4：人源化密码子优化 + 合成订单约束。

把蛋白序列变成能直接下单的 CDS。三个硬约束：

* **人源密码子偏好**：按 Kazusa 人源密码子使用频率表（每千碱基频率）挑同义密码子；
* **GC 平衡**：整条 CDS 的 GC 落在 40%–60%；
* **可合成性**：不出现 ≥6 的同聚物、不引入常见 II 型酶切位点
  （EcoRI/BamHI/HindIII/XbaI/SalI/PstI/NdeI/XhoI/NotI/BsaI）。

诚实边界：这是*密码子优化*，不是表达量预测。序列避开了内部酶切位点，
但外源基因在真核细胞里的实际表达还取决于染色质、UTR、polyA 等——
本模块不承诺表达水平，只承诺"能按这条序列合成、能克隆、能翻译回原蛋白"。
"""

from __future__ import annotations

import random
from typing import Iterable

# Kazusa 人源密码子使用频率（freq per thousand，Codon Usage Database）。
HUMAN_CODONS: dict[str, list[tuple[str, float]]] = {
    "A": [("GCC", 27.7), ("GCT", 18.4), ("GCA", 15.8), ("GCG", 7.4)],
    "R": [("AGG", 12.3), ("AGA", 12.1), ("CGG", 11.4), ("CGC", 10.4),
          ("CGA", 6.2), ("CGT", 4.5)],
    "N": [("AAC", 19.1), ("AAT", 17.0)],
    "D": [("GAC", 25.1), ("GAT", 21.8)],
    "C": [("TGC", 12.6), ("TGT", 10.6)],
    "Q": [("CAG", 34.2), ("CAA", 12.3)],
    "E": [("GAG", 39.6), ("GAA", 29.0)],
    "G": [("GGC", 22.2), ("GGA", 16.5), ("GGG", 16.5), ("GGT", 10.8)],
    "H": [("CAC", 15.1), ("CAT", 10.9)],
    "I": [("ATC", 20.8), ("ATT", 16.0), ("ATA", 7.5)],
    "L": [("CTG", 39.6), ("CTC", 19.6), ("CTT", 13.2), ("TTG", 12.9),
          ("TTA", 7.7), ("CTA", 7.2)],
    "K": [("AAG", 31.9), ("AAA", 24.4)],
    "M": [("ATG", 22.0)],
    "F": [("TTC", 20.3), ("TTT", 17.6)],
    "P": [("CCC", 19.8), ("CCT", 17.5), ("CCA", 16.9), ("CCG", 6.9)],
    "S": [("AGC", 19.5), ("TCC", 17.7), ("TCT", 15.2), ("TCA", 12.2),
          ("AGT", 12.1), ("TCG", 4.4)],
    "T": [("ACC", 18.9), ("ACA", 15.1), ("ACT", 13.1), ("ACG", 6.1)],
    "W": [("TGG", 13.2)],
    "Y": [("TAC", 15.3), ("TAT", 12.2)],
    "V": [("GTG", 28.1), ("GTC", 14.5), ("GTT", 11.0), ("GTA", 7.1)],
    "*": [("TGA", 1.6), ("TAA", 1.0), ("TAG", 0.8)],
}

STOP_CODONS = ("TGA", "TAA", "TAG")

# 反查表：密码子 -> (氨基酸, 该密码子的人源频率‰)，用于回译与稀有罚分。
_BASES: dict[str, tuple[str, float]] = {
    codon: (aa, freq) for aa, cands in HUMAN_CODONS.items() for codon, freq in cands
}

# 常见 II 型限制酶位点：合成订单里最常被要求 6 个字母内不出现的那些。
FORBIDDEN_SITES: dict[str, str] = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "XbaI": "TCTAGA",
    "SalI": "GTCGAC",
    "PstI": "CTGCAG",
    "NdeI": "CATATG",
    "XhoI": "CTCGAG",
    "NotI": "GCGGCCGC",
    "BsaI": "GGTCTC",
}

MAX_HOMOPOLYMER = 5  # 连续同碱基上限（6 个及以上算失败）
GC_LOW, GC_HIGH = 0.40, 0.60
MAX_TRIES = 60


def protein_problems(protein: str) -> list[str]:
    """粗查蛋白序列：占位符、非标准残基。"""
    bad: list[str] = []
    if not protein:
        return ["蛋白序列为空"]
    if "[" in protein or "]" in protein:
        bad.append("蛋白里还有 [元件:状态] 占位符——该架构缺元件序列，不能生成 CDS")
    unknown = sorted({c for c in protein if c not in HUMAN_CODONS and c not in "[]:"})
    if unknown:
        bad.append(f"蛋白里有非标准残基：{''.join(unknown)}")
    return bad


def _gc(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def _max_homopolymer(seq: str) -> int:
    best = run = 1 if seq else 0
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best


def _longest_run(seq: str) -> tuple[int, int] | None:
    """最长同聚物 (起点, 长度)；用于定点修复。"""
    if not seq:
        return None
    best_start, best_len = 0, 1
    start, run = 0, 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
        else:
            start, run = i, 1
        if run > best_len:
            best_start, best_len = start, run
    return best_start, best_len


def _hp_repair(codons: list[str], protein: str, gc_low: float, gc_high: float,
               protect: tuple[int, ...] = (), rounds: int = 400) -> list[str]:
    """定点修同聚物：只动覆盖最长 run 的密码子的同义变体。

    随机爬山容易被卡住——7 连 G 只有极少数替换能一次降到 6，
    所以这里直接瞄准最长 run 所在的密码子逐个试。
    """
    best = _score("".join(codons), gc_low, gc_high)
    for _ in range(rounds):
        run = _longest_run("".join(codons))
        if run is None or run[1] <= MAX_HOMOPOLYMER:
            break
        start, length = run
        touched = sorted({(start + k) // 3 for k in range(length)})
        cand = [i for i in touched
                if i not in protect and i < len(codons)
                and len(HUMAN_CODONS.get(protein[i], [])) > 1]
        improved = False
        for i in cand:
            old = codons[i]
            for alt, _w in HUMAN_CODONS[protein[i]]:
                if alt == old:
                    continue
                codons[i] = alt
                sc = _score("".join(codons), gc_low, gc_high)
                if sc < best - 1.0:
                    best = sc
                    improved = True
                    break
                codons[i] = old
            if improved:
                break
        if not improved:
            break
    return codons


def _sites(seq: str) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    for name, site in FORBIDDEN_SITES.items():
        start = 0
        while True:
            i = seq.find(site, start)
            if i < 0:
                break
            hits.append({"enzyme": name, "site": site, "pos": i + 1})
            start = i + 1
    return hits


def _run_excess(seq: str) -> int:
    """所有超过 MAX_HOMOPOLYMER 的 run 的多余长度之和。

    用"和"而不是"最大 run"：后者在有多条 7 连时无从下降，
    爬山会卡住不动。
    """
    if not seq:
        return 0
    total = 0
    run = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            run += 1
        else:
            run = 1
        if run == MAX_HOMOPOLYMER + 1:
            total += 1
    return total


def _score(cds: str, gc_low: float, gc_high: float) -> float:
    """越小越好：位点/同聚物是硬伤，GC 越界次之，稀有密码子轻微罚分。"""
    pen = 0.0
    pen += 1_000_000.0 * len(_sites(cds))
    pen += 1_000_000.0 * _run_excess(cds)
    gc = _gc(cds)
    if gc < gc_low:
        pen += 1_000.0 * (gc_low - gc)
    elif gc > gc_high:
        pen += 1_000.0 * (gc - gc_high)
    # 稀有密码子轻微罚分：同样合格时优先人源高频密码子。
    rare = 0.0
    for i in range(0, len(cds) - 2, 3):
        _aa, freq = _BASES.get(cds[i:i + 3], ("X", 0.0))
        rare += max(0.0, 20.0 - freq) / 400.0
    return pen + rare


def _polish(codons: list[str], protein: str, rng: random.Random,
            gc_low: float, gc_high: float, steps: int = 2500,
            protect: tuple[int, ...] = ()) -> list[str]:
    """爬山：随机换同义密码子，只要总分下降就接受。"""
    cur = "".join(codons)
    best = _score(cur, gc_low, gc_high)
    idxs = [i for i in range(len(codons))
            if i not in protect and len(HUMAN_CODONS.get(protein[i], [])) > 1]
    if not idxs:
        return codons
    for _ in range(steps):
        if best <= 0.0:
            break
        i = idxs[rng.randrange(len(idxs))]
        alts = [c for c, _w in HUMAN_CODONS[protein[i]] if c != codons[i]]
        rng.shuffle(alts)
        old = codons[i]
        best_alt, best_alt_score = None, best
        for alt in alts:
            codons[i] = alt
            trial = "".join(codons)
            sc = _score(trial, gc_low, gc_high)
            if sc < best_alt_score:
                best_alt, best_alt_score = alt, sc
        if best_alt is None:
            codons[i] = old
            continue
        codons[i] = best_alt
        best = best_alt_score
    return codons


def optimize(protein: str, *, seed: int = 0, add_stop: bool = True,
             gc_range: tuple[float, float] = (GC_LOW, GC_HIGH),
             require_start: bool = True, kozak: str = "GCCACC") -> dict[str, object]:
    """蛋白 → CDS。返回 cds / gc / 位点 / 同聚物 / 是否合格 / 备注。

    require_start：蛋白必须以 M 起始，首密码子固定 ATG（真核翻译起点）。
    add_stop：末尾补一个终止密码子（不进入蛋白序列）。
    kozak：ATG 前的人源 Kozak 序列（默认 GCCACC）；传空串就不加。
    返回的 ``cds`` 是可直接下单的完整序列（含 Kozak），``coding`` 只有编码区。
    """
    kozak = (kozak or "").upper()
    problems = protein_problems(protein)
    if problems:
        return {
            "ok": False, "cds": "", "coding": "", "cds_coding_start": 0,
            "protein_len": len(protein), "gc": None,
            "forbidden_sites": [], "homopolymer_max": None,
            "back_translation": "", "start_codon": "", "stop_codon": "",
            "kozak": kozak, "notes": problems,
        }

    p = protein.rstrip("*").upper()
    if require_start and not p.startswith("M"):
        return {
            "ok": False, "cds": "", "coding": "", "cds_coding_start": 0,
            "protein_len": len(p), "gc": None,
            "forbidden_sites": [], "homopolymer_max": None,
            "back_translation": "", "start_codon": "", "stop_codon": "",
            "kozak": kozak,
            "notes": ["蛋白不是以 M 起始，CDS 没有起始密码子（N 端需补甲硫氨酸）"],
        }
    if add_stop:
        p = p + "*"

    rng = random.Random(seed)
    codons: list[str] = []
    for aa in p:
        cands = HUMAN_CODONS.get(aa)
        if not cands:
            raise ValueError(f"无法翻译残基 {aa!r}")
        head = cands[: max(1, min(3, len(cands)))]
        codons.append(head[rng.randrange(len(head))][0])

    # 起始密码子固定 ATG，且不允许爬山改动。
    protect: tuple[int, ...] = ()
    if p.startswith("M"):
        codons[0] = "ATG"
        protect = (0,)

    gc_low, gc_high = gc_range
    codons = _polish(codons, p, rng, gc_low, gc_high, protect=protect)
    codons = _hp_repair(codons, p, gc_low, gc_high, protect=protect)
    codons = _polish(codons, p, rng, gc_low, gc_high, protect=protect)
    codons = _hp_repair(codons, p, gc_low, gc_high, protect=protect)
    if protect:
        codons[0] = "ATG"
    coding = "".join(codons)
    cds = kozak + coding

    notes: list[str] = []
    hp = _max_homopolymer(cds)
    sites = _sites(cds)
    gc = _gc(coding)
    start_codon = coding[:3]
    stop_codon = coding[-3:]
    if not kozak:
        notes.append("没有加 Kozak 序列（ATG 前为空），人源细胞里翻译起始效率可能偏低")
    if start_codon != "ATG":
        notes.append(f"首密码子不是 ATG：{start_codon}")
    if add_stop and stop_codon not in STOP_CODONS:
        notes.append(f"末密码子不是终止密码子：{stop_codon}")
    if hp > MAX_HOMOPOLYMER:
        notes.append(f"存在 {hp} 个连续相同碱基")
    if not (gc_low <= gc <= gc_high):
        notes.append(f"GC {gc:.1%} 不在 {gc_low:.0%}–{gc_high:.0%}")
    if sites:
        notes.append("残留酶切位点：" + "、".join(sorted({s['enzyme'] for s in sites})))

    ok = (hp <= MAX_HOMOPOLYMER and not sites and gc_low <= gc <= gc_high
          and start_codon == "ATG" and (not add_stop or stop_codon in STOP_CODONS))
    return {
        "ok": ok,
        "cds": cds,
        "coding": coding,
        "cds_coding_start": len(kozak),
        "kozak": kozak,
        "protein_len": len(p),
        "gc": round(gc, 4),
        "forbidden_sites": sites,
        "homopolymer_max": hp,
        "start_codon": start_codon,
        "stop_codon": stop_codon,
        "back_translation": translate(coding),
        "notes": notes,
    }


def translate(cds: str) -> str:
    """翻译回蛋白，用来验证 CDS 没有走样。"""
    return "".join(_BASES.get(cds[i:i + 3], ("X", 0.0))[0]
                   for i in range(0, len(cds) - 2, 3))


def gc_histogram(seqs: Iterable[str]) -> dict[str, float]:
    vals = [_gc(s) for s in seqs if s]
    if not vals:
        return {"n": 0, "min": 0.0, "mean": 0.0, "max": 0.0}
    return {"n": len(vals), "min": round(min(vals), 4),
            "mean": round(sum(vals) / len(vals), 4), "max": round(max(vals), 4)}


__all__ = ["HUMAN_CODONS", "FORBIDDEN_SITES", "STOP_CODONS", "optimize", "translate",
           "protein_problems", "gc_histogram"]
