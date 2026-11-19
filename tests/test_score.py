"""打分与排序测试：可审计、确定性、权重合法性。"""

from __future__ import annotations

import pytest

from conftest import DEMO_SITE
from zfn_design import score
from zfn_design.library import load_default_library
from zfn_design.models import EvidenceLevel, FingerRecord, SequenceRecord, ZfnDesign
from zfn_design.targets import scan_sequence


def test_weights_sum_to_one():
    data = score.load_weights()
    assert abs(sum(data["weights"].values()) - 1.0) < 1e-9


def test_bad_weights_rejected(tmp_path):
    p = tmp_path / "w.json"
    p.write_text('{"weights": {"evidence": 0.9}, "grades": {}, "rules": {}}', encoding="utf-8")
    with pytest.raises(ValueError):
        score.load_weights(p)


def _make_design(seed_library, demo_record):
    sites = [s for s in scan_sequence(demo_record, spacers=(6,), both_strands=False)
             if s.left_half == "GCGTGGGCG" and s.spacer == "GATATC"]
    site = sites[0]
    from zfn_design.designability import evaluate
    from zfn_design.build import build_for_site

    dz = evaluate(site, seed_library)
    assert dz.status == "full", dz.reasons
    return build_for_site(site, seed_library, dz, n_designs=3)[0]


def test_score_components_sum_to_total(seed_library, demo_record):
    d = _make_design(seed_library, demo_record)
    sb = score.score_design(d)
    recomputed = sum(c.raw * c.weight for c in sb.components)
    assert abs(recomputed - sb.total) < 1e-6
    assert abs(sum(c.weight for c in sb.components) - 1.0) < 1e-9
    assert 0.0 <= sb.total <= 100.0


def test_score_always_warns(seed_library, demo_record):
    d = _make_design(seed_library, demo_record)
    sb = score.score_design(d)
    assert any("实验验证" in w for w in sb.warnings)


def test_grade_thresholds():
    g = {"A": 85, "B": 70, "C": 55}
    assert score.grade_of(92, g) == "A"
    assert score.grade_of(85, g) == "A"
    assert score.grade_of(84.9, g) == "B"
    assert score.grade_of(70, g) == "B"
    assert score.grade_of(60, g) == "C"
    assert score.grade_of(54.9, g) == "D"


def test_tier2_design_flagged(seed_library, demo_record):
    """把某个指头换成 Tier2，必须出现'未验证'警告。"""
    d = _make_design(seed_library, demo_record)
    d.fingers[0] = FingerRecord("x", d.fingers[0].triplet, d.fingers[0].helix7,
                                EvidenceLevel.VARIANT_UNVALIDATED, "variant", "synthetic")
    sb = score.score_design(d)
    assert any("Tier2" in w for w in sb.warnings)
    # 证据分项应当低于全 Tier1 的情况
    d2 = _make_design(seed_library, demo_record)
    assert sb.total <= score.score_design(d2).total


def test_ranking_is_deterministic(seed_library, demo_record):
    from zfn_design.designability import evaluate, triage
    from zfn_design.build import build_all
    from zfn_design.targets import scan_many

    sites = scan_many([demo_record], spacers=(6,), both_strands=True)
    dz = triage([evaluate(s, seed_library) for s in sites])
    designs = build_all([d.site for d in dz], seed_library, dz, n_designs=2)
    a = [x.design_id for x in score.rank(designs)]
    b = [x.design_id for x in score.rank(designs)]
    assert a == b
