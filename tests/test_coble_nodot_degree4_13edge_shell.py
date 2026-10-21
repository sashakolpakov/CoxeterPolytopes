from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json

from scripts.certify_coble_nodot_degree4_13edge_shell import OUTPUT


ARCHIVE_SHA256 = "815ea6a701888cb665e375d8fa42cebfd93221ddced48adf66c77fe3aa4287dd"


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_exact_degree_four_shell_census() -> None:
    assert sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert payload["exact_census"] == {
        "maximum_degree_at_most_four_support_count": 16_977,
        "maximum_degree_at_most_three_support_count": 2_888,
        "maximum_degree_exactly_four_support_count": 14_089,
        "weighted_trials": 5_325_642,
        "exact_Lorentzian_marked_presentations": 638,
        "weighted_diagram_isometry_classes": 500,
        "finite_volume_classes": 0,
        "infinite_volume_hyperideal_classes": 500,
    }


def test_pairing_shell_counts_are_exact() -> None:
    counts = _archive()["counts_by_pairing_two_count"]
    assert counts == {
        "0": {
            "weighted_trials": 14_089,
            "exact_Lorentzian_marked_presentations": 1,
            "weighted_diagram_isometry_classes": 1,
        },
        "1": {
            "weighted_trials": 183_157,
            "exact_Lorentzian_marked_presentations": 26,
            "weighted_diagram_isometry_classes": 24,
        },
        "2": {
            "weighted_trials": 1_098_942,
            "exact_Lorentzian_marked_presentations": 189,
            "weighted_diagram_isometry_classes": 145,
        },
        "3": {
            "weighted_trials": 4_029_454,
            "exact_Lorentzian_marked_presentations": 422,
            "weighted_diagram_isometry_classes": 330,
        },
    }


def test_every_weighted_class_has_an_exact_hyperideal_witness() -> None:
    payload = _archive()
    records = payload["class_records"]
    assert len(records) == 500
    assert Counter(record["pairing_two_count"] for record in records) == {
        0: 1,
        1: 24,
        2: 145,
        3: 330,
    }
    assert all(
        record["classification"] == "infinite_volume_hyperideal"
        and record["finite_ideal_hyperideal"][2] > 0
        and record["all_input_walls_irredundant"]
        and record["one_exact_hyperideal_ray_witness"]["type"]
        == "hyperideal"
        and record["one_exact_hyperideal_ray_witness"]["square"] < 0
        and record["positive_strict_interior"]["square"] > 0
        and record["positive_strict_interior"][
            "all_wall_slacks_strictly_positive"
        ]
        and record["literal_vinberg"]["literal_vinberg_passes"] is False
        and record["literal_vinberg"]["stopped_at_hyperideal_endpoint"]
        is True
        for record in records
    )
    assert payload["extreme_ray_type_totals_finite_ideal_hyperideal"] == [
        2_588,
        1_182,
        12_666,
    ]
    assert sum(payload["quotient_determinant_histogram"].values()) == 500


def test_no_algebraic_geometry_leverage_is_claimed() -> None:
    payload = _archive()
    scope = payload["intrinsic_surface_free_shell"]
    assert scope["membership_is_intrinsic_to_weighted_Coxeter_Gram_data"]
    assert scope["surface_or_carrier_success_in_membership"] is False
    assert payload["theorem"]["finite_volume_noncompact_classes"] == 0
    assert payload["theorem"]["compact_classes"] == 0
    assert payload["theorem"]["new_hyperbolic_finite_volume_examples"] == []
    audit = payload["polytope_to_surface_transition_audit"]
    assert audit["carrier_enumeration_reached"] is False
    assert audit["algebraic_geometry_invoked"] is False
    assert audit["independent_AG_pushback_reached"] is False
    assert audit["strict_L4_flagship_gate_crossed"] is False
    assert audit["rating"] == "SCOPED_P_SIDE_NEGATIVE_L2_L3_NOT_L4"


def test_adjacent_subcubic_archive_is_pinned() -> None:
    pinned = _archive()["pinned_sources"]
    assert pinned[
        "data/enumerations/coble_nodot_maxdegree3_closure.json"
    ]["sha256"] == (
        "56a0646809fff76514729a58887172eb916d28cbda699128ed2f96a23556e82d"
    )

