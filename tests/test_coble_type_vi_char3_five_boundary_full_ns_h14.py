from __future__ import annotations

from collections import Counter
import json

from scripts.certify_coble_type_vi_char3_five_boundary_full_ns_h14 import (
    OUTPUT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_frozen_certificate_replays_exactly() -> None:
    assert build_certificate() == _archive()


def test_five_boundary_surface_has_complete_finite_volume_full_nef_cone() -> None:
    payload = _archive()
    assert payload["all_checks_pass"]
    assert payload["full_NS_lattice"]["identification"] == "odd unimodular I_(1,14)"
    assert payload["full_NS_lattice"]["canonical_square"] == -5
    assert payload["curve_exhaustion"]["negative_curve_types"] == {
        "boundary_minus4": 5,
        "minus1": 45,
        "minus2": 10,
    }
    assert payload["curve_exhaustion"]["minus1_decomposition"] == {
        "ten_meeting_two_boundaries_once": 10,
        "thirty_five_meeting_one_boundary_twice": 35,
    }
    masks = payload["curve_exhaustion"]["boundary_lift_membership_by_sorted_CM_cusp"]
    assert Counter(sum(mask) for mask in masks) == {0: 22, 1: 35}
    assert Counter(
        index
        for mask in masks
        for index, accepted in enumerate(mask)
        if accepted
    ) == {0: 7, 1: 7, 2: 7, 3: 7, 4: 7}

    literal = payload["literal_vinberg"]
    assert literal["decision"] == "finite_volume_noncompact"
    assert (
        literal["facet_count"],
        literal["finite_vertex_count"],
        literal["ideal_vertex_count"],
        literal["elliptic_edge_count"],
    ) == (60, 27_845, 1_027, 244_330)
    assert literal["bad_edge_count"] == 0
    assert literal["bad_parabolic_neighborhood_count"] == 0

    cone = payload["exact_double_description"]
    assert (
        cone["extreme_ray_count"],
        cone["finite_vertex_count"],
        cone["ideal_vertex_count"],
        cone["hyperideal_vertex_count"],
    ) == (28_872, 27_845, 1_027, 0)
    assert cone["all_facets_irredundant"]
    assert set(cone["all_sixty_facet_face_ranks"]) == {14}


def test_arithmetic_lineage_and_marked_fold_are_exact_and_scoped() -> None:
    payload = _archive()
    arithmetic = payload["arithmetic"]
    assert arithmetic["arithmetic_lattice"] == "2*Pic~(S), Nikulin (15,7,1)"
    assert arithmetic["identification"] == "U(2)+E8(-1)+A1(-1)^5"
    assert arithmetic["primitive_square_divisibility_histogram"] == {
        "-2,1": 15,
        "-2,2": 35,
        "-4,2": 10,
    }
    assert arithmetic["walls_span_whole_arithmetic_lattice"]
    assert arithmetic["all_60_reflections_integral"]

    fold = arithmetic["boundary_frame_fold"]
    assert fold["boundary_indices_zero_based"] == [0, 1, 2, 3, 4]
    assert fold["boundary_gram"] == [
        [-4, 0, 0, 0, 0],
        [0, -4, 0, 0, 0],
        [0, 0, -4, 0, 0],
        [0, 0, 0, -4, 0],
        [0, 0, 0, 0, -4],
    ]
    assert fold["negative_primitive_restriction_count"] == 20
    assert fold["isotropic_redundant_restriction_count"] == 35
    assert fold["recovered_CM_chamber"] == "Kondo Type VI, H9/20"
    assert "not by itself" in fold["inverse_boundary"]

    lineage = payload["lineage"]
    assert lineage["exact_identification"] == "Alexeev (15,7,1) = Shimada D5 H14/60"
    assert sorted(lineage["type_VI_to_Alexeev_one_based"]) == list(range(1, 61))
    assert lineage["decision"] == {
        "new_abstract_polytope": False,
        "positive_gain": (
            "new project-level complete full-nef realization by Kondo's "
            "five-boundary characteristic-three Type-VI Coble surface"
        ),
        "same_abstract_polytope": True,
        "same_arithmetic_reflection_chamber": True,
        "worldwide_priority": (
            "closed known for the abstract chamber via Alexeev; priority for "
            "this full-NS Coble realization has not been asserted"
        ),
    }
    assert payload["checks"][
        "source_surface_is_char3_five_boundary_not_Martin_char5_one_boundary"
    ]
    assert "an unmarked polytope-to-surface theorem" in payload["claim"]["not_claimed"]
