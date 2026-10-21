from __future__ import annotations

import json

from scripts.certify_coble_four_point_constructible_incidence import (
    OUTPUT,
    build_certificate,
)


def test_naive_free_incidence_predicate_fails_on_both_inputs() -> None:
    certificate = build_certificate()
    assert certificate["naive_free_incidence_verdict"]["decision"] == (
        "FAILS_ON_BOTH_CERTIFIED_HIGH_DIMENSIONAL_INPUTS"
    )
    assert certificate["checks"]["both_naive_free_predicates_fail"]


def test_complete_four_point_design_is_exactly_characteristic_three() -> None:
    theorem = build_certificate()["four_point_projective_line_theorem"]
    assert theorem["common_nonemptiness_iff"]["solution"] == (
        "char=3, lambda=2=-1=1/2"
    )
    assert theorem["grid_graph_design"][
        "all_24_evaluation_ranks_at_char3_lambda2"
    ] == {"3": 24}
    assert theorem["conic_polar_chord_design"][
        "all_six_evaluation_ranks_at_char3_lambda2"
    ] == {"2": 6}


def test_two_label_free_dependency_signatures_are_recovered() -> None:
    signatures = build_certificate()["P_derived_dependency_signatures"]
    rank15 = signatures["rank15_60wall_input"]
    rank18 = signatures["rank18_48wall_input"]
    assert rank15["raw_surface_labels_read"] is False
    assert rank18["raw_surface_labels_read"] is False
    assert rank15["boundary_count"] == 5
    assert rank15["total_blowup_centers"] == 14
    assert rank15["residual_root_curve_count"] == 6
    assert rank15["complete_pair_complement_set"]
    assert rank18["boundary_count"] == 8
    assert rank18["transverse_grid_centers"] == 16
    assert rank18["root_curve_count"] == 24
    assert rank18["complete_permutation_set"]


def test_core_class_is_two_member_L3_and_extension_gate_is_prospective() -> None:
    certificate = build_certificate()
    constructible = certificate["constructible_incidence_class"]
    assert constructible["certified_member_count"] == 2
    assert constructible["not_defined_as_union_of_two_archived_grams"]
    assert constructible["proved_scope"]["zero_extension_required"]
    prospective = constructible["prospective_extension_verifier"]
    assert prospective["status"] == "NOT_A_PROVED_NONEMPTINESS_THEOREM"
    assert prospective["certified_nonzero_extension_count"] == 0
    assert certificate["hostile_level"]["score"] == (
        "L3_TWO_MEMBER_CORE_CALIBRATION;_NOT_L3_PLUS_OR_L4"
    )
    assert certificate["all_checks_pass"]


def test_frozen_archive_matches_live_certificate() -> None:
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == build_certificate()
