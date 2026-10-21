from __future__ import annotations

import json

from scripts.certify_coble_remaining_seven_contraction_designs import (
    OUTPUT,
    build_certificate,
)


CERTIFICATE = build_certificate()


def test_all_seven_rows_are_excluded_from_the_two_ci4_cores() -> None:
    assert len(CERTIFICATE["seven_row_outcome"]) == 7
    assert all(
        not row["CI4_membership"]
        for row in CERTIFICATE["seven_row_outcome"].values()
    )
    status = CERTIFICATE["CI4_final_status"]
    assert not status["same_core_found_in_remaining_seven"]
    assert not status["nonzero_standard_free_extension_found"]
    assert "at L3" in status["rating"]


def test_four_small_rows_have_exhaustive_sequential_outcomes() -> None:
    rows = CERTIFICATE["seven_row_outcome"]
    assert (rows["one_boundary_22"]["exhaustive_sequential_P2_targets"], rows["one_boundary_22"]["endpoint_degree_signature_count"]) == (74, 16)
    assert (rows["one_boundary_30"]["exhaustive_sequential_P2_targets"], rows["one_boundary_30"]["endpoint_degree_signature_count"]) == (573, 28)
    assert (rows["two_boundary_26"]["exhaustive_sequential_P2_targets"], rows["two_boundary_26"]["endpoint_degree_signature_count"]) == (52, 9)
    assert (rows["two_boundary_36"]["exhaustive_sequential_P2_targets"], rows["two_boundary_36"]["endpoint_degree_signature_count"]) == (890, 28)


def test_raw_206_wall_input_recovers_the_mobius_plane() -> None:
    core = CERTIFICATE["two_boundary_206_mobius_core"]
    raw = core["raw_extraction"]
    assert raw["surface_labels_read"] is False
    assert raw["boundary_count"] == 2
    assert raw["crossing_minus_one_count"] == 10
    assert raw["contracted_endpoint_gram"] == [[0, 1], [1, 0]]
    assert raw["determinant_one_NS_basis"]
    assert raw["root_count"] == 30
    assert raw["raw_design_parameters"] == "3-(10,4,1)"


def test_f9_realization_and_strong_coble_numerics_are_exact() -> None:
    core = CERTIFICATE["two_boundary_206_mobius_core"]
    geometry = core["finite_geometry"]
    assert geometry["F9_graph_evaluation_rank_histogram"] == {"3": 30, "4": 180}
    assert geometry["all_thirty_unique_smooth_1_1_curves"]
    coble = core["strong_Coble_numerics"]
    assert coble["boundary_sum_equals_minus_2K"]
    assert coble["anti_canonical_O_2_2_evaluation_rank"] == 9
    assert coble["anti_canonical_kernel_dimension"] == 0
    assert core["exact_return_handoff"]["returned_wall_count"] == 206


def test_q_factory_is_coble_root_rigid_exactly_at_q_three() -> None:
    theorem = CERTIFICATE["mobius_plane_factory_rigidity_theorem"]
    assert theorem["formulas"] == {
        "B_square": "-(q-1)^2",
        "circle_square": "1-q",
        "K_dot_circle": "q-3",
        "B_plus_Bprime": "-2K+(q-3)(F+G)",
    }
    samples = theorem["sample_exact_arithmetic"]
    assert samples["3"]["Coble_root_package"]
    assert all(
        not record["Coble_root_package"]
        for q, record in samples.items()
        if q != "3"
    )


def test_vi_and_vii_are_recorded_only_as_unfinished_design_candidates() -> None:
    candidates = CERTIFICATE["different_finite_design_candidates"]
    assert "not_promoted" in candidates["one_boundary_78"]["status"]
    assert "not_promoted" in candidates["one_boundary_72"]["status"]


def test_frozen_archive_matches_live_certificate() -> None:
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == CERTIFICATE
