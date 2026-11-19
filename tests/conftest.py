"""共享测试夹具。

刻意**不**下载任何网络数据：所有测试用内存里构造的指头库，
唯一例外是随包的种子库（它本身就是 PDB 1AAY 的可复核数据）。
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from zfn_design.library import FingerLibrary, load_default_library  # noqa: E402
from zfn_design.models import EvidenceLevel, FingerRecord, SequenceRecord  # noqa: E402

#: Zif268 (PDB 1AAY) 蛋白序列 —— 全项目的真值锚点
ZIF268_PROTEIN = (
    "MERPYACPVESCDRRFSRSDELTRHIRIHTGQKPFQCRICMRNFSRSDHLTTHIRTHTGEKPFACDICGRKFA"
    "RSDERKRHTKIHLRQKD"
)

#: 1AAY 上链 DNA；9 bp 结合位点 = AGCGTGGGCGT 中间那段
ZIF268_DNA_TOP = "AGCGTGGGCGT"
ZIF268_HALFSITE = "GCGTGGGCG"

#: 一个 9+6+9 的完整二聚体位点：左半位点 GCGTGGGCG，spacer GATATC，右半位点 GCGTGGGCG
DEMO_SITE = "GCGTGGGCGGATATCCGCCCACGC"


@pytest.fixture
def seed_library() -> FingerLibrary:
    """3 条 PDB 种子库（可再分发的那份），与默认库无关。"""
    from zfn_design.library import FingerLibrary as _FL

    path = Path(__file__).resolve().parents[1] / "zfn_design" / "data" / "finger_library_seed.csv"
    return _FL.load_csv(path, name="seed_1AAY")


@pytest.fixture
def default_library() -> FingerLibrary:
    """随包默认库：有 B1H 全库时就是它，否则退回种子库。"""
    return load_default_library()


@pytest.fixture
def tiny_library() -> FingerLibrary:
    """手工构造的小库：覆盖 GCG / TGG / AAA，含一个 Tier2 条目。"""
    recs = [
        FingerRecord("h_GCG_1", "GCG", "RSDELTR", EvidenceLevel.LITERATURE, "pdb", "1AAY", "crystal_structure", split_position="F1"),
        FingerRecord("h_GCG_2", "GCG", "RSDERKR", EvidenceLevel.LITERATURE, "pdb", "1AAY", "crystal_structure", split_position="F3"),
        FingerRecord("h_GCG_3", "GCG", "RSDKLTR", EvidenceLevel.VARIANT_UNVALIDATED, "variant", "synthetic", "none", split_position="any"),
        FingerRecord("h_TGG_1", "TGG", "RSDHLTT", EvidenceLevel.LITERATURE, "pdb", "1AAY", "crystal_structure", split_position="F2"),
        FingerRecord("h_AAA_1", "AAA", "QSGNLTE", EvidenceLevel.CODE_INFERENCE, "rule", "recognition_code", "none", split_position="any"),
    ]
    return FingerLibrary(recs, name="tiny")


@pytest.fixture
def demo_record() -> SequenceRecord:
    seq = "TACGTAGCTAGCTTAGCTAGGCTA" + DEMO_SITE + "TTACGGATCCTTAGGCAT"
    return SequenceRecord(name="demo_target", seq=seq, source="fixture")

#: 需要完整 B1H 指头库（发布包默认不带）才能跑出设计的集成用例。
from zfn_design.library import has_full_library  # noqa: E402

requires_full_library = pytest.mark.skipif(
    not has_full_library(),
    reason="需要完整 B1H 指头库才能生成设计；见 zfn_design/data/finger_library/README.md",
)
