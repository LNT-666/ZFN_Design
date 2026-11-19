"""核心数据模型。

所有跨模块传递的对象都在这里定义，避免各模块各造一套字典。
全部 ``dataclass`` 都可 JSON 序列化（:meth:`to_dict`）。
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from enum import Enum
from typing import Any


class EvidenceLevel(str, Enum):
    """指头证据等级——本方案最重要的一个字段。

    绝不允许把低等级冒充高等级；报告里会把等级原样打给用户看。
    """

    #: Tier 1：有实验实测（B1H / B2H / reporter / 晶体结构 + 亲和力），有 DOI。
    MEASURED = "measured"
    #: Tier 1：有文献表格/数据库出处，但本仓库未逐条复核原始数据。
    LITERATURE = "literature"
    #: Tier 2：由已知指头做定点突变衍生，**未验证**。
    VARIANT_UNVALIDATED = "variant_unvalidated"
    #: Tier 2：按识别密码规则枚举出来的候选，**纯计算、未验证**。
    CODE_INFERENCE = "code_inference"
    #: Tier 3：完全没有候选，报缺口。
    NONE = "none"

    @property
    def tier(self) -> int:
        return {
            EvidenceLevel.MEASURED: 1,
            EvidenceLevel.LITERATURE: 1,
            EvidenceLevel.VARIANT_UNVALIDATED: 2,
            EvidenceLevel.CODE_INFERENCE: 2,
            EvidenceLevel.NONE: 3,
        }[self]

    @property
    def label_zh(self) -> str:
        return {
            EvidenceLevel.MEASURED: "实测（Tier1）",
            EvidenceLevel.LITERATURE: "文献来源（Tier1）",
            EvidenceLevel.VARIANT_UNVALIDATED: "突变衍生·未验证（Tier2）",
            EvidenceLevel.CODE_INFERENCE: "规则推断·未验证（Tier2）",
            EvidenceLevel.NONE: "缺口（Tier3）",
        }[self]

    @classmethod
    def parse(cls, raw: str) -> "EvidenceLevel":
        key = (raw or "").strip().lower().replace("-", "_").replace(" ", "_")
        for m in cls:
            if m.value == key or m.name.lower() == key:
                return m
        raise ValueError(f"未知 evidence_level: {raw!r}；合法值 {[m.value for m in cls]}")


@dataclass(frozen=True)
class FingerRecord:
    """一条"指头—三联体"记录。"""

    helix_id: str
    triplet: str
    helix7: str
    evidence: EvidenceLevel
    source_type: str = ""          # pdb / b1h / zifdb / paper / variant
    source_id: str = ""            # 1AAY / DOI / 数据库主键
    assay: str = ""                # B1H / crystal / reporter / none
    affinity_nM: float | None = None
    split_position: str = "any"    # F2 / F3 / any，B1H 数据区分指位
    notes: str = ""

    @property
    def helix6(self) -> str:
        """识别螺旋的 -1,1,2,3,5,6 六个位置（Persikov 2015 的 core sequence 表示法）。"""
        if len(self.helix7) < 7:
            return ""
        return "".join(self.helix7[i] for i in (0, 1, 2, 3, 5, 6))

    @property
    def core4(self) -> str:
        """四个主要特异性决定位 -1,2,3,6。"""
        if len(self.helix7) < 7:
            return ""
        return "".join(self.helix7[i] for i in (0, 2, 3, 6))

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["evidence"] = self.evidence.value
        d["helix6"] = self.helix6
        d["core4"] = self.core4
        d["tier"] = self.evidence.tier
        return d


@dataclass
class SequenceRecord:
    """一条输入序列。"""

    name: str
    seq: str
    description: str = ""
    source: str = ""

    @property
    def length(self) -> int:
        return len(self.seq)

    def to_dict(self) -> dict[str, Any]:
        return {"name": self.name, "length": self.length,
                "description": self.description, "source": self.source}


@dataclass
class TargetSite:
    """一个 ZFN 位点。两种形态（``kind``）：

    ``"paired"``（组合型，经典 ZFN）
        左半位点 + spacer + 右半位点，两个单体各认一侧，FokI 在 spacer 上二聚化。
    ``"single"``（单体型，ZF-TF / 单 ZF 阵列）
        只有一个半位点，没有 spacer、没有搭档。``spacer`` 与 ``right_half`` 都是空串，
        因此下面的不变量依然成立。

    坐标约定
    --------
    * ``start`` 是 **输入序列** 上的 1-based 起点（位点占据的最左碱基）。
    * ``strand`` 说明位点朝哪边读。
    * ``site_seq`` 是位点按"读链方向"写出来的串，恒满足::

          site_seq == left_half + spacer + revcomp(right_half)

      ``strand == "+"`` 时 ``site_seq`` 等于输入序列对应片段；
      ``strand == "-"`` 时等于输入序列对应片段的反向互补。

    为什么要存两份：``right_half`` 是右单体自己的 5'→3' 序列，它在双链上
    以反向互补形式出现。只存一份最容易在"到底该用哪条链"上出错，
    所以两者都存下来，并有单元测试保证互相一致。
    """

    seq_name: str
    start: int
    strand: str
    left_half: str
    spacer: str
    right_half: str
    site_seq: str
    kind: str = "paired"

    @property
    def is_paired(self) -> bool:
        return self.kind == "paired"

    @property
    def spacer_len(self) -> int:
        return len(self.spacer)

    @property
    def half_sites(self) -> list[str]:
        if not self.is_paired:
            return [self.left_half]
        return [self.left_half, self.right_half]

    @property
    def end(self) -> int:
        return self.start + len(self.site_seq) - 1

    @property
    def input_region(self) -> tuple[int, int]:
        """位点在输入序列上占据的 1-based 闭区间。"""
        return self.start, self.end

    @property
    def cut_offset_in_site(self) -> int:
        """0-based 偏移（相对 ``site_seq``）。

        组合型：spacer 中点，FokI 在那一带切割（**近似**位置）。
        单体型：没有切割位点，返回结合位点中心，只用于报告展示。
        """
        if not self.is_paired:
            return len(self.left_half) // 2
        return len(self.left_half) + (self.spacer_len + 1) // 2 - 1

    @property
    def cut_site_input(self) -> int:
        """切割点换算回输入序列的 1-based 坐标（负链会自动翻转）。"""
        if self.strand == "+":
            return self.start + self.cut_offset_in_site
        return self.end - self.cut_offset_in_site

    @property
    def site_id(self) -> str:
        return f"{self.seq_name}:{self.start}:{self.strand}"

    def to_dict(self) -> dict[str, Any]:
        return {
            "site_id": self.site_id,
            "kind": self.kind,
            "is_paired": self.is_paired,
            "seq_name": self.seq_name,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "left_half": self.left_half,
            "spacer": self.spacer,
            "spacer_len": self.spacer_len,
            "right_half": self.right_half,
            "site_seq": self.site_seq,
            "cut_site_input_approx": self.cut_site_input,
        }


@dataclass
class OffTargetHit:
    """一个（潜在的）脱靶位点。"""

    seq_name: str
    start: int
    strand: str
    left_half: str
    spacer: str
    right_half: str
    mm_left: int
    mm_right: int

    @property
    def total_mm(self) -> int:
        return self.mm_left + self.mm_right

    @property
    def weight(self) -> float:
        """粗略的"像不像真位点"权重：0 错配 1.0，1 错配 0.25，2 错配 0.06，
        3/4 错配只用于排序（0.015 / 0.004），不代表"可能被切"。"""
        table = {0: 1.0, 1: 0.25, 2: 0.06, 3: 0.015, 4: 0.004}
        return table.get(self.total_mm, 0.0)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["total_mm"] = self.total_mm
        d["weight"] = self.weight
        return d


@dataclass
class OffTargetReport:
    hits: list[OffTargetHit] = field(default_factory=list)
    searched_bp: int = 0
    max_mm_left: int = 2
    max_mm_right: int = 2
    spacers: tuple[int, ...] = (5, 6, 7)
    truncated: bool = False

    @property
    def burden(self) -> float:
        """脱靶负担：所有命中权重之和（含 on-target 时请先剔除）。"""
        return sum(h.weight for h in self.hits)

    @property
    def exact_hits(self) -> int:
        return sum(1 for h in self.hits if h.total_mm == 0)

    def summary(self) -> dict[str, Any]:
        by_mm: dict[int, int] = {}
        for h in self.hits:
            by_mm[h.total_mm] = by_mm.get(h.total_mm, 0) + 1
        return {
            "searched_bp": self.searched_bp,
            "max_mm_left": self.max_mm_left,
            "max_mm_right": self.max_mm_right,
            "spacer_lengths": list(self.spacers),
            "n_hits": len(self.hits),
            "n_exact": self.exact_hits,
            "hits_by_mismatch": {str(k): by_mm[k] for k in sorted(by_mm)},
            "burden": round(self.burden, 4),
            "truncated": self.truncated,
        }

    def to_dict(self) -> dict[str, Any]:
        return {"summary": self.summary(), "hits": [h.to_dict() for h in self.hits]}


@dataclass
class ScoreComponent:
    """打分的一个可审计分项。"""

    key: str
    label_zh: str
    raw: float          # 0-100
    weight: float
    detail: str = ""

    @property
    def contribution(self) -> float:
        return self.raw * self.weight

    def to_dict(self) -> dict[str, Any]:
        return {"key": self.key, "label_zh": self.label_zh, "raw": round(self.raw, 2),
                "weight": self.weight, "weighted": round(self.contribution, 3),
                "detail": self.detail}


@dataclass
class ScoreBreakdown:
    components: list[ScoreComponent]
    total: float
    grade: str
    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "total": round(self.total, 2),
            "grade": self.grade,
            "components": [c.to_dict() for c in self.components],
            "warnings": list(self.warnings),
            "notes": list(self.notes),
        }


@dataclass
class ZfnDesign:
    """一个完整设计。

    * 组合型（``target.is_paired``）：一对 ZF 阵列（左单体 + 右单体）+ FokI 二聚体。
    * 单体型：只有一条 ZF 阵列（``partner_fingers`` / ``partner_protein`` 为空）。
    """

    design_id: str
    target: TargetSite
    fingers: list[FingerRecord]          # 左单体（单体型时就是唯一那条），N→C
    partner_fingers: list[FingerRecord] = field(default_factory=list)  # 右单体，N→C
    protein: str = ""
    partner_protein: str = ""
    score: ScoreBreakdown | None = None
    offtarget: OffTargetReport | None = None
    offtarget_score: dict[str, Any] | None = None  # 全基因组脱靶负担（单独评分，不进 Design score）
    architecture: str = "paired_nuclease"
    notes: list[str] = field(default_factory=list)

    @property
    def is_paired(self) -> bool:
        return self.target.is_paired

    @property
    def triplets(self) -> list[str]:
        return [f.triplet for f in self.fingers]

    @property
    def partner_triplets(self) -> list[str]:
        return [f.triplet for f in self.partner_fingers]

    @property
    def all_fingers(self) -> list[FingerRecord]:
        return list(self.fingers) + list(self.partner_fingers)

    @property
    def evidence_levels(self) -> list[str]:
        return [f.evidence.value for f in self.all_fingers]

    @property
    def min_tier(self) -> int:
        return max(f.evidence.tier for f in self.all_fingers)

    def to_dict(self) -> dict[str, Any]:
        return {
            "design_id": self.design_id,
            "architecture": self.architecture,
            "is_paired": self.is_paired,
            "target": self.target.to_dict(),
            "monomer_left": {
                "halfsite": self.target.left_half,
                "fingers": [f.to_dict() for f in self.fingers],
                "triplets_N_to_C": self.triplets,
                "protein": self.protein,
            },
            "monomer_right": {
                "halfsite": self.target.right_half,
                "fingers": [f.to_dict() for f in self.partner_fingers],
                "triplets_N_to_C": self.partner_triplets,
                "protein": self.partner_protein,
            } if self.is_paired else None,
            "min_evidence_tier": self.min_tier,
            "score": self.score.to_dict() if self.score else None,
            "offtarget": self.offtarget.to_dict() if self.offtarget else None,
            "offtarget_score": self.offtarget_score,
            "notes": list(self.notes),
        }
