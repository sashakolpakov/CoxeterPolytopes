from __future__ import annotations

import json

from scripts.certify_polytope_first_coble_type_vi_h14 import (
    OUTPUT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_inverse_certificate_replays_from_only_the_raw_arithmetic_data() -> None:
    payload = build_certificate()
    assert payload == _archive()
    assert payload["input"]["discovery_fields_used"] == [
        "arithmetic.gram",
        "arithmetic.wall_vectors",
    ]
    assert "facets.labels" in payload["input"]["surface_fields_not_used"]
    assert "arithmetic.boundary_frame_fold" in payload["input"][
        "surface_fields_not_used"
    ]


def test_all_candidate_frames_are_enumerated_and_only_one_descales() -> None:
    payload = _archive()
    search = payload["exhaustive_boundary_frame_search"]
    assert search["prefix_survivor_counts"] == {
        "adjunction_alphabet": 1,
        "all_five_subsets": 3003,
        "all_physical_walls_descend": 1,
        "canonical_characteristic": 1,
        "canonical_square_minus_five": 1,
        "half_form_integral": 1,
        "odd_unimodular": 1,
        "pairwise_orthogonal": 338,
        "surviving_frames": 1,
    }
    assert search["unique_unordered_frame"] is True
    assert search["surface_labels_consulted"] is False
    ns = payload["reconstructed_Picard_lattice"]
    assert ns["identification"] == "odd unimodular I_(1,14)"
    assert ns["determinant"] == 1
    assert ns["canonical_square"] == -5
    assert ns["canonical_is_characteristic"] is True
    assert ns["minus_two_K_equals_boundary_sum"] is True


def test_root_fold_and_all_cusp_lifts_are_intrinsically_complete() -> None:
    payload = _archive()
    assert payload["wall_taxonomy"]["counts"] == {
        "boundary_minus4": 5,
        "crossing_minus1_curve": 10,
        "cusp_lift_minus1_curve": 35,
        "minus2_curve": 10,
    }
    cm = payload["Coble_Mukai_reconstruction"]
    assert cm["identification"] == "E10=U+E8(-1)"
    assert cm["roots_span_whole_CM"] is True
    assert cm["root_marking_counts"] == {"minus1_root": 10, "minus2": 10}
    assert cm["Type_VI_gate"]["labelled_gram_isomorphic_to_Kondo_Gamma_VI"]
    assert cm["Type_VI_gate"]["labelled_isometry_count"] == 120
    assert cm["exact_cone"]["finite_vertex_count"] == 600
    assert cm["exact_cone"]["ideal_vertex_count"] == 57
    assert cm["exact_cone"]["hyperideal_vertex_count"] == 0
    assert cm["exact_cone"]["all_causal_rays_in_one_time_sheet"] is True

    lifts = payload["cusp_lift_exhaustion"]
    assert lifts["lift_mask_histogram"] == {
        "00000": 22,
        "00001": 7,
        "00010": 7,
        "00100": 7,
        "01000": 7,
        "10000": 7,
    }
    assert lifts["integral_cusp_boundary_pairs"] == 35
    assert lifts["all_and_only_integral_lifts_are_raw_minus_one_walls"] is True


def test_exact_same_polytope_return_and_scope_are_not_overstated() -> None:
    payload = _archive()
    returned = payload["exact_same_P_return"]
    assert returned["oriented_mirror_sets_equal"] is True
    assert returned["reconstructed_mirror_count"] == 60
    assert returned["full_H14_cone"]["finite_vertex_count"] == 27_845
    assert returned["full_H14_cone"]["ideal_vertex_count"] == 1_027
    assert returned["full_H14_cone"]["hyperideal_vertex_count"] == 0
    assert returned["full_H14_cone"]["all_facets_irredundant"] is True
    assert returned["full_H14_cone"]["all_causal_rays_in_one_time_sheet"] is True

    gate = payload["algebro_geometric_gate"]
    assert gate["characteristic"] == 3
    assert gate["surface_isomorphism_classes"] == 1
    assert gate["R_invariant"] == "(E6 + D5, Z/2Z)"
    assert gate["full_nef_carrier_exists"] is True
    assert "automatic Coble realization for every numerical package passing C0-C2" in payload[
        "scope"
    ]["not_proved"]
    assert "a theorem about all arithmetic Coxeter polytopes" in payload["scope"][
        "not_proved"
    ]
