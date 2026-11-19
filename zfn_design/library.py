"""指头库（finger library）的读写、查询与导入。

本模块是整个方案里**唯一**允许"知识"进入系统的入口。规则是：

1. 任何一条记录都必须带 ``evidence_level``，没有等级的记录直接拒收；
2. 同一三联体可以有多条记录（不同指位、不同来源、不同亲和力）；
3. 查询接口按证据等级优先、亲和力次之排序，调用方不需要自己排序；
4. 缺失就是缺失——``candidates()`` 返回空列表，由上层报 Tier3 缺口，
   绝不在库里"造"一个凑数的指头。
"""

from __future__ import annotations

import csv
import io
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .models import EvidenceLevel, FingerRecord
from .sequtil import SequenceError

FIELDS = [
    "helix_id",
    "triplet",
    "helix7",
    "evidence_level",
    "source_type",
    "source_id",
    "assay",
    "affinity_nM",
    "split_position",
    "notes",
]

_EVIDENCE_RANK = {
    EvidenceLevel.MEASURED: 0,
    EvidenceLevel.LITERATURE: 1,
    EvidenceLevel.VARIANT_UNVALIDATED: 2,
    EvidenceLevel.CODE_INFERENCE: 3,
    EvidenceLevel.NONE: 4,
}


class LibraryError(ValueError):
    """指头库格式或内容错误。"""


def _to_float(raw: str | None) -> float | None:
    if raw is None:
        return None
    s = str(raw).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


class FingerLibrary:
    """不可变的指头集合（内部按三联体建索引）。"""

    def __init__(self, records: Sequence[FingerRecord], *, name: str = "library"):
        self.name = name
        self._records: list[FingerRecord] = list(records)
        self._by_triplet: dict[str, list[FingerRecord]] = defaultdict(list)
        for r in self._records:
            self._by_triplet[r.triplet].append(r)
        for lst in self._by_triplet.values():
            lst.sort(key=lambda r: (_EVIDENCE_RANK[r.evidence],
                                    r.affinity_nM if r.affinity_nM is not None else 1e9,
                                    r.helix_id))
        # 基因组模式每秒要问几十万次 candidates()，缓存是必须的
        self._cand_cache: dict[tuple[str, str | None, int], list[FingerRecord]] = {}

    # ---------------------------------------------------------------- 基本属性
    def __len__(self) -> int:
        return len(self._records)

    def __iter__(self):
        return iter(self._records)

    @property
    def records(self) -> list[FingerRecord]:
        return list(self._records)

    @property
    def triplets(self) -> set[str]:
        return set(self._by_triplet)

    def get(self, helix_id: str) -> FingerRecord | None:
        for r in self._records:
            if r.helix_id == helix_id:
                return r
        return None

    # ---------------------------------------------------------------- 查询
    def candidates(self, triplet: str, *, split_position: str | None = None,
                   max_tier: int = 2) -> list[FingerRecord]:
        """返回某三联体的候选指头，按可靠性排序。

        参数
        ----
        triplet: 3 bp，如 ``"GCG"``
        split_position: ``"F1"`` / ``"F2"`` / ``"F3"`` / ``None``。
            B1H 数据区分 F2/F3 位置；``any`` 的记录对任何位置都可用。
        max_tier: 允许的最高证据 tier（1 = 只要实测/文献；2 = 也接受推断）。
        """
        triplet = triplet.upper()
        key = (triplet, split_position, max_tier)
        cached = self._cand_cache.get(key)
        if cached is not None:
            return cached
        out = []
        for r in self._by_triplet.get(triplet, []):
            if r.evidence.tier > max_tier:
                continue
            if split_position and r.split_position not in ("any", split_position):
                continue
            out.append(r)
        self._cand_cache[key] = out
        return list(out)

    def coverage(self) -> dict[str, int]:
        """``{三联体: 可用记录数}``。"""
        return {t: len(v) for t, v in sorted(self._by_triplet.items())}

    def coverage_stats(self) -> dict[str, object]:
        total_triplets = 64
        covered = len(self._by_triplet)
        curated = sum(1 for r in self._records if r.evidence.tier == 1)
        return {
            "library_name": self.name,
            "n_records": len(self._records),
            "n_curated_records": curated,
            "triplets_covered": covered,
            "triplets_total": total_triplets,
            "coverage_fraction": round(covered / total_triplets, 4),
            "by_evidence": {
                lvl.value: sum(1 for r in self._records if r.evidence is lvl)
                for lvl in EvidenceLevel
                if any(r.evidence is lvl for r in self._records)
            },
        }

    def missing_triplets(self, triplets: Iterable[str]) -> list[str]:
        """给定一组三联体，返回库里没有候选的那些（保序去重）。"""
        out, seen = [], set()
        for t in triplets:
            t = t.upper()
            if t in seen:
                continue
            seen.add(t)
            if not self._by_triplet.get(t):
                out.append(t)
        return out

    # ---------------------------------------------------------------- 读写
    @classmethod
    def from_rows(cls, rows: Iterable[Mapping[str, str]], *, name: str = "library") -> "FingerLibrary":
        recs: list[FingerRecord] = []
        errors: list[str] = []
        for i, row in enumerate(rows, start=2):  # 2 = 表头之后的第一行
            try:
                triplet = (row.get("triplet") or "").strip().upper()
                helix7 = (row.get("helix7") or "").strip().upper()
                if not re.fullmatch(r"[ACGT]{3}", triplet):
                    raise LibraryError(f"triplet 非法: {triplet!r}")
                if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{7}", helix7):
                    raise LibraryError(f"helix7 必须是 7 个氨基酸: {helix7!r}")
                ev = EvidenceLevel.parse(row.get("evidence_level") or "")
                helix_id = (row.get("helix_id") or "").strip() or f"{triplet}_{helix7}"
                recs.append(
                    FingerRecord(
                        helix_id=helix_id,
                        triplet=triplet,
                        helix7=helix7,
                        evidence=ev,
                        source_type=(row.get("source_type") or "").strip(),
                        source_id=(row.get("source_id") or "").strip(),
                        assay=(row.get("assay") or "").strip(),
                        affinity_nM=_to_float(row.get("affinity_nM")),
                        split_position=(row.get("split_position") or "any").strip() or "any",
                        notes=(row.get("notes") or "").strip(),
                    )
                )
            except (LibraryError, ValueError, SequenceError) as exc:
                errors.append(f"第 {i} 行: {exc}")
        if errors:
            raise LibraryError("指头库有 %d 处错误:\n  %s" % (len(errors), "\n  ".join(errors)))
        return cls(recs, name=name)

    @classmethod
    def load_csv(cls, path: str | Path, *, name: str | None = None) -> "FingerLibrary":
        p = Path(path)
        if not p.exists():
            raise LibraryError(f"指头库文件不存在: {p}")
        text = p.read_text(encoding="utf-8-sig")
        reader = csv.DictReader(io.StringIO(text))
        missing = [f for f in ("triplet", "helix7", "evidence_level") if f not in (reader.fieldnames or [])]
        if missing:
            raise LibraryError(f"{p.name} 缺少必需列: {missing}")
        return cls.from_rows(reader, name=name or p.stem)

    def to_csv(self, path: str | Path) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        with p.open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader()
            for r in self._records:
                w.writerow({
                    "helix_id": r.helix_id,
                    "triplet": r.triplet,
                    "helix7": r.helix7,
                    "evidence_level": r.evidence.value,
                    "source_type": r.source_type,
                    "source_id": r.source_id,
                    "assay": r.assay,
                    "affinity_nM": "" if r.affinity_nM is None else r.affinity_nM,
                    "split_position": r.split_position,
                    "notes": r.notes,
                })
        return p

    # ---------------------------------------------------------------- 外部数据导入
    @classmethod
    def import_b1h_dir(cls, directory: str | Path, *, split_position: str = "any",
                       top_n: int = 10, min_freq: float = 0.0, pos4: str = "L",
                       source_id: str = "Persikov2015_B1H") -> "FingerLibrary":
        """导入 Persikov 2015 (NAR, PMC4330361) 的 B1H 数据集目录。

        该数据集从 https://zf.princeton.edu/b1h/download.html 免费申请获取
        （非商业用途，需填姓名/单位/邮箱）。每个三联体一个文本文件，每行形如::

            RRDTTN 0.0069884

        即识别螺旋 -1,1,2,3,5,6 六个位置的氨基酸 + 该螺旋在测序数据中的频率
        （频率已按三联体归一，和为 1）。第 4 位（位置 4）不接触 DNA、也不在
        该数据集中，本函数用 ``pos4``（默认 ``L``）填成完整 7 残基串，并在
        ``notes`` 里保留原始 6-mer 与这条假设——**打分不使用第 4 位**，但它会
        出现在蛋白序列里，合成前应人工确认。

        先导出的实测数据里，各列位置可用"位置 1 和 5 的熵显著更高"这条
        已知规律核对：本数据集第 2 列（索引 1）与第 5 列（索引 4）的平均熵
        最高（F2: 2.76 / 3.54 bit），与 (-1,1,2,3,5,6) 的约定一致。

        文件名可以是 ``GCG.txt``、``gcg.txt``、``GCG_F2.txt``、或 Princeton
        官方压缩包里的 ``GCG_amino_acid_union_cut.txt``（此时位置由
        ``split_position`` 参数给出，F2/F3 两个目录分别导入再合并）。
        """
        directory = Path(directory)
        if not directory.is_dir():
            raise LibraryError(f"不是目录: {directory}")
        files = [f for f in directory.rglob("*") if f.is_file() and f.suffix.lower() in (".txt", ".dat", ".tsv")]
        if not files:
            raise LibraryError(f"{directory} 下没有 .txt/.dat/.tsv 文件")

        recs: list[FingerRecord] = []
        for f in files:
            m = re.search(r"([ACGTacgt]{3})", f.stem)
            if not m:
                continue
            triplet = m.group(1).upper()
            pos_in_name = re.search(r"[_\-](F[123])", f.stem)
            rec_split = pos_in_name.group(1) if pos_in_name else split_position
            entries: list[tuple[float, str]] = []
            for line in f.read_text(encoding="utf-8", errors="replace").splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = re.split(r"[\s,\t]+", line)
                if len(parts) < 2:
                    continue
                # 两种写法都接受：
                #   "RRDTTN 0.0069884"        （Princeton 官方格式：6 残基连写）
                #   "R R D T T N 0.0069884"   （逐残基分列）
                if len(parts) >= 7 and all(len(p) == 1 for p in parts[:6]):
                    aa = "".join(parts[:6]).upper()
                    freq = _to_float(parts[6])
                else:
                    aa = parts[0].upper()
                    freq = _to_float(parts[-1])
                if not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]{6}", aa):
                    continue
                if freq is None or freq < min_freq:
                    continue
                entries.append((freq, aa))
            entries.sort(key=lambda x: -x[0])
            for rank, (freq, aa6) in enumerate(entries[:top_n], start=1):
                helix7 = aa6[:4] + pos4 + aa6[4:]
                recs.append(
                    FingerRecord(
                        helix_id=f"B1H_r{rank:02d}_{rec_split}_{triplet}",
                        triplet=triplet,
                        helix7=helix7,
                        evidence=EvidenceLevel.MEASURED,
                        source_type="b1h",
                        source_id=source_id,
                        assay="B1H",
                        affinity_nM=None,
                        split_position=rec_split,
                        notes=(f"freq={freq:.4g}（该三联体内归一）；"
                               f"源 6-mer={aa6}；第 4 位未测，按 {pos4} 填"),
                    )
                )
        if not recs:
            raise LibraryError("导入结果为空——请检查文件格式是否与 B1H 数据集一致")
        return cls(recs, name=f"b1h:{directory.name}")

    def merged(self, other: "FingerLibrary") -> "FingerLibrary":
        """合并两个库；同名 ``helix_id`` 保留先出现的。"""
        seen = {r.helix_id for r in self._records}
        recs = list(self._records) + [r for r in other if r.helix_id not in seen]
        return FingerLibrary(recs, name=f"{self.name}+{other.name}")


def default_library_path(data_dir: Path | None = None) -> Path:
    """指头库查找顺序：``data/`` 直放 → ``data/finger_library/`` 用户自备 → 3 条 PDB 种子。

    B1H 全库（Persikov et al. 2015）不再随仓库分发，需用户自行申请后放入
    ``data/finger_library/finger_library_b1h.csv``；缺失时退回到仅够跑通流程的种子库。
    """
    data = Path(data_dir) if data_dir is not None else Path(__file__).resolve().parent / "data"
    candidates = (
        data / "finger_library_b1h.csv",
        data / "finger_library" / "finger_library_b1h.csv",
    )
    for p in candidates:
        if p.exists():
            return p
    return data / "finger_library_seed.csv"


def has_full_library(data_dir: Path | None = None) -> bool:
    """是否已放入 B1H 全库（未放入时只能跑种子库）。"""
    return default_library_path(data_dir).name == "finger_library_b1h.csv"


def load_default_library() -> FingerLibrary:
    p = default_library_path()
    name = "b1h_persikov2015" if p.name == "finger_library_b1h.csv" else "seed_1AAY"
    return FingerLibrary.load_csv(p, name=name)
