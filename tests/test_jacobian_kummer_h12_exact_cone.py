from __future__ import annotations

import json

from scripts.independent_certify_jacobian_kummer_h12 import (
    INEQUALITY_OUTPUT,
    RAY_OUTPUT,
    SUMMARY_OUTPUT,
    load_ray_archive,
)


def _certificate() -> dict:
    return json.loads(SUMMARY_OUTPUT.read_text(encoding="utf-8"))


def test_exact_cone_certificate_has_no_hyperideal_vertex() -> None:
    payload = _certificate()
    cone = payload["cone"]
    vertices = payload["exact_double_description"]

    assert cone["hyperbolic_dimension"] == 12
    assert cone["facet_count"] == 172
    assert cone["inequality_rank"] == 13
    assert cone["full_dimensional"] is True
    assert cone["pointed"] is True
    assert vertices["extreme_ray_count"] == 119_488
    assert vertices["finite_vertex_count"] == 113_274
    assert vertices["ideal_vertex_count"] == 6_214
    assert vertices["hyperideal_vertex_count"] == 0
    assert vertices["finite_vertex_active_facet_count_histogram"] == {
        "12": 113_274
    }
    assert vertices["ideal_vertex_active_facet_count_histogram"] == {
        "13": 2_304,
        "14": 3_256,
        "15": 576,
        "16": 76,
        "18": 2,
    }


def test_exact_vinberg_endpoint_certificate() -> None:
    payload = _certificate()
    affine = payload["affine_full_neighborhood_cross_check"]
    edges = payload["vinberg_edge_endpoint_audit"]

    assert affine["connected_affine_component_count"] == 17_074
    assert affine["distinct_full_neighborhood_count"] == 6_214
    assert affine["invalid_full_neighborhood_count"] == 0
    assert affine["equals_exact_lightlike_ray_incidence_sets"] is True

    assert edges["elliptic_rank_11_edge_count"] == 910_316
    assert edges["endpoint_count_histogram"] == {"2": 910_316}
    assert edges["finite_finite_edge_count"] == 505_568
    assert edges["finite_ideal_edge_count"] == 348_152
    assert edges["ideal_ideal_edge_count"] == 56_596
    assert edges["bad_edge_count"] == 0
    assert payload["decision"]["classification"] == "finite_volume_noncompact"
    assert all(payload["checks"].values())


def test_replayable_gmp_input_and_compressed_ray_archive() -> None:
    inequality_lines = INEQUALITY_OUTPUT.read_text(encoding="utf-8").splitlines()
    assert inequality_lines[:3] == [
        "H-representation",
        "begin",
        "172 14 rational",
    ]
    assert inequality_lines[-1] == "end"

    header, records = load_ray_archive(RAY_OUTPUT)
    assert header["schema"] == "jacobian_kummer_h12_exact_cone_rays_v1"
    assert header["ray_count"] == 119_488
    assert header["coordinate_basis_facet_indices"] == list(range(1, 14))
    assert len(records) == 119_488
    assert all(len(ray) == 13 and len(active) >= 12 for ray, active in records)
