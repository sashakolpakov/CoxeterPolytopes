from __future__ import annotations

import json

from scripts.enumerate_jacobian_kummer_h10_class0007_components import (
    OUTPUT as COMPONENT_CERTIFICATE,
)
from scripts.independent_certify_jacobian_kummer_h10_class0007 import (
    SUMMARY_OUTPUT as CONE_CERTIFICATE,
    load_ray_archive,
)
from scripts.screen_jacobian_kummer_h10_class0007 import OUTPUT as SCREEN_CERTIFICATE
from scripts.verify_jacobian_kummer_h10_class0007_vinberg import (
    OUTPUT as VINBERG_CERTIFICATE,
)


def _load(path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_class0007_exact_fixed_slice_facets_and_coxeter_gate() -> None:
    screen = _load(SCREEN_CERTIFICATE)

    assert screen["group_class"] == screen["group_class"] | {
        "gap_subgroup_class_index": 7,
        "subgroup_structure": "C2",
        "conjugacy_class_size": 60,
        "normalizer_order": 384,
        "representative_trace": 5,
        "fixed_rank": 11,
        "hyperbolic_dimension": 10,
    }
    assert screen["group_class"]["surface_automorphism_required"] is False
    assert screen["wall_orbit_count"] == 174
    assert screen["wall_orbit_size_histogram"] == {"1": 32, "2": 142}
    assert screen["orbit_gram_inertia"] == [1, 10, 163]
    assert screen["facet_count"] == 84
    assert screen["redundant_orbit_count"] == 90
    assert screen["facet_normal_square_histogram"] == {
        "-1": 64,
        "-2": 8,
        "-4": 12,
    }
    assert screen["facet_pair_count"] == 3_486
    assert screen["coxeter_edge_counts"] == {
        "2": 1_166,
        "3": 72,
        "4": 112,
        "infinity": 860,
        "ultraparallel": 1_276,
    }
    assert screen["non_coxeter_pair_count"] == 0
    assert screen["status"]["all_facet_pairs_are_coxeter"] is True


def test_class0007_complete_components_and_literal_vinberg_closure() -> None:
    components = _load(COMPONENT_CERTIFICATE)
    vinberg = _load(VINBERG_CERTIFICATE)

    assert components["elliptic_order_histogram"] == {
        "1": 84,
        "2": 184,
        "3": 456,
        "4": 528,
        "5": 480,
        "6": 384,
    }
    assert components["affine_order_histogram"] == {
        "2": 860,
        "3": 712,
        "4": 606,
        "5": 348,
        "6": 280,
        "7": 144,
    }
    assert components["distinct_full_zero_neighborhood_count"] == 862
    assert components["invalid_full_zero_neighborhood_count"] == 0

    assert vinberg["logically_independent_of_exact_double_description"] is True
    assert vinberg["finite_vertex_count"] == 15_698
    assert vinberg["ideal_vertex_count"] == 862
    assert vinberg["elliptic_rank_9_edge_count"] == 107_233
    assert vinberg["endpoint_count_histogram"] == {"2": 107_233}
    assert vinberg["endpoint_type_histogram"] == {
        "FF": 56_560,
        "FI": 43_860,
        "II": 6_813,
    }
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["unrecognized_endpoint_edge_count"] == 0
    assert vinberg["invalid_full_zero_neighborhood_count"] == 0
    assert vinberg["vinberg_decision"] == "finite_volume_noncompact"
    assert 2 * vinberg["elliptic_rank_9_edge_count"] == (
        vinberg["finite_endpoint_incidence_count"]
        + vinberg["ideal_endpoint_incidence_count"]
    )


def test_class0007_independent_exact_cone_agrees() -> None:
    cone = _load(CONE_CERTIFICATE)
    vinberg = _load(VINBERG_CERTIFICATE)
    exact = cone["exact_double_description"]
    edges = cone["vinberg_edge_endpoint_audit"]

    assert cone["cone"]["facet_count"] == 84
    assert cone["cone"]["inequality_rank"] == 11
    assert cone["cone"]["lineality_dimension"] == 0
    assert exact["extreme_ray_count"] == 16_560
    assert exact["finite_vertex_count"] == vinberg["finite_vertex_count"]
    assert exact["ideal_vertex_count"] == vinberg["ideal_vertex_count"]
    assert exact["hyperideal_vertex_count"] == 0
    assert exact["finite_vertex_active_facet_count_histogram"] == {"10": 15_698}
    assert cone["affine_full_neighborhood_cross_check"][
        "equals_exact_lightlike_ray_incidence_sets"
    ] is True
    assert edges["elliptic_rank_9_edge_count"] == vinberg["elliptic_rank_9_edge_count"]
    assert edges["endpoint_count_histogram"] == {"2": 107_233}
    assert edges["finite_finite_edge_count"] == 56_560
    assert edges["finite_ideal_edge_count"] == 43_860
    assert edges["ideal_ideal_edge_count"] == 6_813
    assert edges["bad_edge_count"] == 0
    assert edges["equals_complete_elliptic_rank_9_edge_set"] is True
    assert cone["decision"] == cone["decision"] | {
        "finite_volume": True,
        "compact": False,
        "classification": "finite_volume_noncompact",
    }

    header, rays = load_ray_archive()
    assert header["schema"] == "jacobian_kummer_h10_class0007_exact_cone_rays_v1"
    assert header["ray_count"] == len(rays) == 16_560
    assert len({active for _, active in rays}) == 16_560
