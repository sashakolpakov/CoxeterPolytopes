from __future__ import annotations

import json

from scripts.build_jacobian_kummer_h12_interim import OUTPUT as METRIC_CERTIFICATE
from scripts.enumerate_jacobian_kummer_h12_affine_components import (
    OUTPUT as AFFINE_CERTIFICATE,
)
from scripts.independent_certify_jacobian_kummer_h12 import (
    SUMMARY_OUTPUT as CONE_CERTIFICATE,
)
from scripts.verify_jacobian_kummer_h12_vinberg import (
    OUTPUT as VINBERG_CERTIFICATE,
)


def _load(path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_exact_h12_metric_and_complete_affine_census() -> None:
    metric = _load(METRIC_CERTIFICATE)
    affine = _load(AFFINE_CERTIFICATE)

    assert metric["group_class"]["fixed_rank"] == 13
    assert metric["group_class"]["hyperbolic_dimension"] == 12
    assert metric["facet_count"] == 172
    assert metric["redundant_orbit_count"] == 8
    assert metric["wall_orbit_size_histogram"] == {"1": 44, "2": 136}
    assert metric["facet_normal_square_histogram"] == {
        "-1": 124,
        "-2": 40,
        "-4": 8,
    }
    assert metric["coxeter_edge_counts"] == {
        "2": 2_820,
        "3": 144,
        "4": 800,
        "infinity": 1_750,
        "ultraparallel": 9_192,
    }

    assert len(affine["connected_elliptic_components"]) == 16_908
    assert len(affine["connected_affine_components"]) == 17_074
    assert affine["affine_order_histogram"] == {
        "2": 1_750,
        "3": 2_784,
        "4": 3_420,
        "5": 3_768,
        "6": 1_760,
        "7": 544,
        "8": 1_128,
        "9": 1_536,
        "10": 384,
    }
    assert affine["distinct_full_zero_neighborhood_count"] == 6_214
    assert affine["invalid_full_zero_neighborhood_count"] == 0


def test_pure_diagram_vinberg_certificate_is_independent_and_closes() -> None:
    vinberg = _load(VINBERG_CERTIFICATE)
    cone = _load(CONE_CERTIFICATE)

    assert vinberg["logically_independent_of_cdd_extreme_ray_certificate"] is True
    assert vinberg["finite_vertex_count"] == 113_274
    assert vinberg["ideal_vertex_count"] == 6_214
    assert vinberg["elliptic_rank_11_edge_count"] == 910_316
    assert vinberg["endpoint_count_histogram"] == {"2": 910_316}
    assert vinberg["endpoint_type_histogram"] == {
        "FF": 505_568,
        "FI": 348_152,
        "II": 56_596,
    }
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["unrecognized_endpoint_edge_count"] == 0
    assert vinberg["invalid_full_zero_neighborhood_count"] == 0
    assert vinberg["vinberg_decision"] == "finite_volume_noncompact"

    assert 2 * vinberg["elliptic_rank_11_edge_count"] == (
        vinberg["finite_endpoint_incidence_count"]
        + vinberg["ideal_endpoint_incidence_count"]
    )
    assert vinberg["finite_endpoint_incidence_count"] == (
        2 * vinberg["endpoint_type_histogram"]["FF"]
        + vinberg["endpoint_type_histogram"]["FI"]
    )
    assert vinberg["ideal_endpoint_incidence_count"] == (
        vinberg["endpoint_type_histogram"]["FI"]
        + 2 * vinberg["endpoint_type_histogram"]["II"]
    )

    independent = cone["exact_double_description"]
    assert independent["finite_vertex_count"] == vinberg["finite_vertex_count"]
    assert independent["ideal_vertex_count"] == vinberg["ideal_vertex_count"]
    assert independent["hyperideal_vertex_count"] == 0
    assert cone["vinberg_edge_endpoint_audit"]["elliptic_rank_11_edge_count"] == (
        vinberg["elliptic_rank_11_edge_count"]
    )
