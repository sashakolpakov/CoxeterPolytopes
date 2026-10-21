from __future__ import annotations

import json

from scripts.certify_coble_vinberg_incidence_inverse_transfer import (
    OUTPUT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_live_replay_matches_archive() -> None:
    assert build_certificate() == _archive()


def test_denominator_is_surface_free_and_does_not_encode_classification() -> None:
    payload = _archive()
    raw = payload["surface_free_denominator"]
    assert raw["name"] == "D_CobVin"
    rendered = json.dumps(raw["definition"], sort_keys=True)
    assert "Kondo" not in rendered
    assert "Type VI" not in rendered
    assert "Type MII" not in rendered
    assert "a surface or an effective curve" in raw["not_a_hypothesis"]
    assert "the bound b in {1,2,5,8}" in raw["not_a_hypothesis"]
    assert raw["known_distinct_positive_inputs"] == 2


def test_finite_compiler_and_exact_return_are_stated_with_correct_quantifiers() -> None:
    payload = _archive()
    compiler = payload["finite_all_branch_compiler"]
    theorem = payload["promotion_and_return_theorem"]
    assert "P2 or a Hirzebruch surface" in compiler["branches"]
    assert "finite union of finite-type" in compiler["incidence_space"]
    assert "at most rank(N) contractions" in compiler["termination"]
    assert "if and only if" in theorem["intrinsic_iff"]
    assert "literally the oriented input chamber P" in theorem["same_P"]
    assert "Proposition Vinberg2" in theorem["root_exhaustion"]
    assert "primitive fiber multiplicity vector (2,1)" in theorem[
        "minus_one_exhaustion"
    ]


def test_high_dimensional_positive_calibrations() -> None:
    rows = {
        row["raw_id"]: row
        for row in _archive()["classification_pushback"][
            "positive_high_dimensional_calibrations"
        ]
    }
    assert set(rows) == {"H14/60", "H17/48"}
    assert rows["H14/60"]["boundary_count"] == 5
    assert rows["H14/60"]["CM_root_count"] == 20
    assert rows["H14/60"]["CM_cusp_count"] == 57
    assert rows["H14/60"]["single_boundary_lift_count"] == 35
    assert rows["H14/60"]["incidence_characteristics_not_two"] == [3]
    assert rows["H17/48"]["boundary_count"] == 8
    assert rows["H17/48"]["CM_root_count"] == 40
    assert rows["H17/48"]["CM_cusp_count"] == 202
    assert rows["H17/48"]["single_boundary_lift_count"] == 0
    assert rows["H17/48"]["incidence_characteristics_not_two"] == [3]
    assert all(_archive()["checks"].values())


def test_scope_is_l3_not_a_false_l4_claim() -> None:
    payload = _archive()
    score = payload["five_part_score"]
    assert score["C_surface_free_class"] == "pass"
    assert score["Car_finite_exhaustive_carrier_enumeration"] == "pass"
    assert score["Car_realized_branch_for_every_P_in_C"].startswith("fail:")
    assert score["Car_overall_under_strict_contract"] == "fail"
    assert "pass" in score["Return_exact_same_P"]
    assert score["Pushback_new_theorem_on_all_inputs"].startswith("fail:")
    assert score["rating"].startswith("carrier-relative L3")
    assert "strict programme L2-L3, not L4" in score["rating"]
    missing = payload["scope"]["exact_L4_missing_gate"]
    assert "force at least one incidence stack to be nonempty" in missing
    assert any(
        "new abstract H14 or H17" in claim
        for claim in payload["scope"]["not_claimed"]
    )
