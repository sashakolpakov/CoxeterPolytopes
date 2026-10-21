from __future__ import annotations

from collections import Counter
import gzip
import hashlib
import json

import sympy as sp

from scripts.certify_coble_one_boundary_full_ns_h10 import (
    _calibration_walls,
    _load_templates,
)
from scripts.certify_coble_type_vii_char5_one_boundary_full_ns_h10 import (
    EDGE_OUTPUT,
    EXPECTED_P10_INDEX,
    EXPECTED_RAY_HASH,
    RAY_OUTPUT,
    SUMMARY_OUTPUT,
    build_certificate,
)


SUMMARY_SHA256 = "a528ca8c70f2af325d93404634afd0f9fd9800bfc9f979716051c32abee21242"
RAYS_SHA256 = "10662253f402c11fa3826e31ad2d05ead38463667d1cfd373f874fd6f7a84df9"
EDGES_SHA256 = "bc57acaaedf1600aa8792ae6074aebd8ef50c77c986d2483b8ec93dc696b9765"


def _sha256(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_archived_type_vii_full_nef_theorem() -> None:
    assert _sha256(SUMMARY_OUTPUT) == SUMMARY_SHA256
    assert _sha256(RAY_OUTPUT) == RAYS_SHA256
    assert _sha256(EDGE_OUTPUT) == EDGES_SHA256
    document = json.loads(SUMMARY_OUTPUT.read_text(encoding="utf-8"))
    assert document["schema"] == "coxeter.coble_type_vii_char5_one_boundary_full_ns_h10.v1"
    assert document["direction"] == "SURFACE_TO_POLYTOPE"
    assert document["inverse_carrier_score"]["level"] == "N/A/L0"
    assert document["surface"] == {
        "automorphism_group": "S5",
        "base_field": "algebraically closed, characteristic 5",
        "boundary_components": 1,
        "classification": "Kondo Type VII strong Coble surface",
        "construction": "Kondo Example 4.15 and Proposition 4.16",
        "effective_irreducible_CM_roots": 20,
        "uniqueness": "unique isomorphism class in Kondo Theorem 1.1/Table 2",
    }
    bridge = document["algebro_geometric_completeness"][
        "intrinsic_unique_even_E10_overlattice_bridge"
    ]
    assert bridge["root_lattice_determinant"] == -4
    assert bridge["root_lattice_discriminant_smith_diagonal"] == [1] * 8 + [2, 2]
    assert set(bridge["nonzero_q_value_set"]) == {"-5/2", "-2", "-7/2"}
    assert bridge["unique_nonzero_even_isotropic_class_count"] == 1
    assert bridge["unique_even_overlattice_is_even_unimodular_E10"] is True
    assert bridge["marking_uniqueness_not_claimed"] is True
    polytope = document["polytope"]
    cone = polytope["exact_double_description"]
    vinberg = polytope["literal_vinberg"]
    assert (polytope["dimension"], polytope["facet_count"]) == (10, 72)
    assert polytope["classification"] == "finite_volume_noncompact"
    assert (cone["finite_vertex_count"], cone["ideal_vertex_count"], cone["hyperideal_vertex_count"]) == (18_761, 1_731, 0)
    assert cone["all_facets_irredundant"] is True
    assert set(cone["all_facet_face_ranks"]) == {10}
    assert cone["primitive_ray_list_sha256"] == EXPECTED_RAY_HASH
    assert cone["all_extreme_rays_on_one_lorentz_time_sheet"] is True
    assert cone["minimum_strict_interior_to_extreme_ray_pairing"] > 0
    assert vinberg["elliptic_rank_9_edge_diagram_count"] == 143_565
    assert vinberg["elliptic_rank_10_finite_endpoint_diagram_count"] == 18_761
    assert vinberg["edge_endpoint_type_histogram"] == {
        "FF": 60_250,
        "FI": 67_110,
        "II": 16_205,
    }
    assert vinberg["connected_affine_component_count"] == 4_943
    assert vinberg["distinct_full_parabolic_neighbourhood_count"] == 1_731
    assert vinberg["bad_edge_count"] == vinberg["bad_parabolic_neighbourhood_count"] == 0
    assert polytope["solomon"]["index_in_W_P10"] == EXPECTED_P10_INDEX
    assert document["checks"]["finite_volume"] is True
    assert document["checks"]["noncompact"] is True
    assert document["all_checks_pass"] is True


def test_archived_complete_ray_and_edge_ledgers() -> None:
    with gzip.open(RAY_OUTPUT, "rt", encoding="utf-8") as handle:
        ray_header = json.loads(next(handle))
        rays = [json.loads(line) for line in handle]
    assert ray_header["ray_count"] == len(rays) == 20_492
    assert Counter(record["kind"] for record in rays) == {
        "finite": 18_761,
        "ideal": 1_731,
    }
    assert Counter(
        (record["kind"], len(record["active_facets"])) for record in rays
    ) == {
        ("finite", 10): 18_761,
        ("ideal", 10): 120,
        ("ideal", 11): 500,
        ("ideal", 12): 721,
        ("ideal", 13): 295,
        ("ideal", 14): 90,
        ("ideal", 15): 5,
    }
    canonical = sorted(tuple(record["primitive_coordinate_basis_vector"]) for record in rays)
    assert hashlib.sha256(
        json.dumps(canonical, separators=(",", ":")).encode()
    ).hexdigest() == EXPECTED_RAY_HASH

    with gzip.open(EDGE_OUTPUT, "rt", encoding="utf-8") as handle:
        edge_header = json.loads(next(handle))
        endpoint_histogram = Counter()
        edge_count = 0
        seen_supports = set()
        for line in handle:
            support, endpoints, kinds, signature = json.loads(line)
            assert len(support) == 9 and len(set(support)) == 9
            assert len(endpoints) == 2 and endpoints[0] != endpoints[1]
            assert all(kind in {"finite", "ideal"} for kind in kinds)
            assert sum(int(kind[1:]) for kind in signature) == 9
            seen_supports.add(tuple(support))
            endpoint_histogram["".join(sorted(kind[0].upper() for kind in kinds))] += 1
            edge_count += 1
    assert edge_header["edge_count"] == edge_count == len(seen_supports) == 143_565
    assert endpoint_histogram == {"FF": 60_250, "FI": 67_110, "II": 16_205}


def test_type_vii_wall_construction_rebuilds_entrywise() -> None:
    chamber_data, ray_data = _load_templates()
    lattice, walls, construction = _calibration_walls("VII", chamber_data, ray_data)
    document = json.loads(SUMMARY_OUTPUT.read_text(encoding="utf-8"))
    polytope = document["polytope"]
    assert construction == {
        "cm_cusp_count": 51,
        "cm_root_count": 20,
        "wall_order": "b=-K; CM (-2)-roots; cusp lifts E=f-b",
    }
    assert sp.Matrix(polytope["lattice_gram"]) == lattice
    assert polytope["wall_vectors"] == [[int(value) for value in wall] for wall in walls]
    wall_matrix = sp.Matrix.hstack(*walls).T
    assert sp.Matrix(polytope["negative_wall_gram"]) == wall_matrix * lattice * wall_matrix.T
    assert [int((wall.T * lattice * wall)[0]) for wall in walls] == [-1] + [-2] * 20 + [-1] * 51
    boundary = walls[0]
    for exceptional in walls[21:]:
        cusp = exceptional + boundary
        assert (cusp.T * lattice * cusp)[0] == 0
        assert (cusp.T * lattice * boundary)[0] == 0


def test_full_exact_cdd_diagram_and_vinberg_replay() -> None:
    rebuilt, ray_archive, edge_archive = build_certificate()
    archived = json.loads(SUMMARY_OUTPUT.read_text(encoding="utf-8"))
    assert rebuilt == archived
    assert ray_archive == RAY_OUTPUT.read_bytes()
    assert edge_archive == EDGE_OUTPUT.read_bytes()
