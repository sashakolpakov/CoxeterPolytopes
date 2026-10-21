from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json

from scripts.certify_coble_nodot_maxdegree3_closure import OUTPUT


ARCHIVE_SHA256 = "56a0646809fff76514729a58887172eb916d28cbda699128ed2f96a23556e82d"


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_complete_subcubic_denominator_and_two_survivors() -> None:
    assert sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    payload = _archive()
    assert payload["all_checks_pass"] is True
    scope = payload["surface_free_denominator"]
    assert scope["surface_or_carrier_success_in_membership"] is False
    assert "between 11 and 18 edges" in scope["support_edge_range_is_exhaustive"]
    assert payload["combined_exact_census"] == {
        "weighted_trials": 10_456_309,
        "exact_Lorentzian_marked_presentations": 1_038,
        "weighted_root_isometry_classes": 837,
        "finite_volume_root_isometry_classes": 2,
        "finite_volume_root_types": ["Type I", "Type II"],
        "complete_cusp_lifted_polytope_isometry_classes": 2,
        "exact_Coble_same_P_returns": 2,
    }


def test_new_low_and_high_tails_are_exact() -> None:
    tails = _archive()["new_tail_censuses"]
    low = tails["low_11_13"]
    high = tails["high_16_18"]
    assert (
        low["weighted_trials"],
        low["exact_Lorentzian_marked_presentations"],
        low["weighted_isometry_classes"],
    ) == (1_394_476, 156, 130)
    assert (
        high["weighted_trials"],
        high["exact_Lorentzian_marked_presentations"],
        high["weighted_isometry_classes"],
    ) == (3_268_247, 408, 302)
    assert low["volume_ledger"] == {
        "finite_volume_noncompact": 1,
        "infinite_volume_hyperideal": 129,
    }
    assert high["volume_ledger"] == {"infinite_volume_hyperideal": 302}
    assert len(low["no_hyperideal_records"]) == 1
    assert low["no_hyperideal_records"][0]["known_root_diagram_matches"] == [
        "Type I"
    ]
    assert high["no_hyperideal_records"] == []
    assert Counter(record["classification"] for record in low["class_records"]) == {
        "infinite_volume_hyperideal": 129,
        "finite_volume_noncompact": 1,
    }
    assert Counter(record["classification"] for record in high["class_records"]) == {
        "infinite_volume_hyperideal": 302
    }
    assert all(
        record["positive_strict_interior"]["square"] > 0
        and record["positive_strict_interior"][
            "all_wall_slacks_strictly_positive"
        ]
        for tail in (low, high)
        for record in tail["class_records"]
    )


def test_p10_is_a_real_counterexample_to_maximal_automaticity() -> None:
    counterexample = _archive()["maximal_automaticity_counterexample"]
    assert counterexample["D_CobVin_arithmetic_gates_pass"] is True
    assert counterexample["complete_raw_package_regenerated"] is True
    assert counterexample["root_count"] == 10
    assert counterexample["complete_cusp_lift_count"] == 1
    assert counterexample["base_field_scope"] == (
        "complex numbers (characteristic zero)"
    )
    assert counterexample[
        "exact_complex_one_boundary_strong_Coble_full_nef_carrier_exists"
    ] is False
    assert counterexample["consequence"] == (
        "R_CobVin(C) is a proper subset of D_CobVin"
    )


def test_strict_l4_score_is_not_inflated() -> None:
    audit = _archive()["strict_L4_audit"]
    assert audit["class_wide_carrier_nonemptiness"].startswith("PASS")
    assert audit["exact_same_P_return"] == "PASS"
    assert audit["nonthin_new_AG_pushback"] == "FAIL"
    assert audit["rating"] == (
        "L3_PLUS_CLASSWIDE_SCOPED_ESSENTIAL_SURJECTIVITY_NOT_STRICT_L4"
    )
