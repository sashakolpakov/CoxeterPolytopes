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
from scripts.certify_coble_type_vi_char5_one_boundary_full_ns_h10 import (
    EDGE_OUTPUT,
    EXPECTED_P10_INDEX,
    EXPECTED_RAY_HASH,
    RAY_OUTPUT,
    SUMMARY_OUTPUT,
    build_certificate,
)


SUMMARY_SHA256 = "2f58973a6e7305032588681cdec0f2186d9ce640207a72eeecadb6397a0bfffc"
RAYS_SHA256 = "72ec876fcc319d5366f90303fed2ce7c6d1b1d6f37d3f2a601898ab91de1c438"
EDGES_SHA256 = "d841d5d747041824be1d1b3c839cb64f5665130e75bd7cc94f07e5d7b342d6db"


def _sha256(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_archived_type_vi_full_nef_theorem() -> None:
    assert _sha256(SUMMARY_OUTPUT) == SUMMARY_SHA256
    assert _sha256(RAY_OUTPUT) == RAYS_SHA256
    assert _sha256(EDGE_OUTPUT) == EDGES_SHA256
    document = json.loads(SUMMARY_OUTPUT.read_text(encoding="utf-8"))
    assert document["schema"] == (
        "coxeter.coble_type_vi_char5_one_boundary_full_ns_h10.v1"
    )
    assert document["direction"] == "SURFACE_TO_POLYTOPE"
    assert document["inverse_carrier_score"]["level"] == "N/A/L0"
    assert document["surface"] == {
        "R_invariant": "(E6+A4,{0})",
        "automorphism_group": "S5",
        "base_field": "algebraically closed, characteristic 5",
        "boundary_components": 1,
        "classification": "Kondo Type VI strong Coble surface",
        "construction": "Kondo Example 4.12 and Proposition 4.13",
        "effective_irreducible_CM_roots": 20,
        "uniqueness": "unique isomorphism class in Kondo Theorem 1.1/Table 2",
    }
    bridge = document["algebro_geometric_completeness"][
        "intrinsic_full_E10_marking_bridge"
    ]
    assert bridge["root_coordinate_span_smith_diagonal"] == [1] * 10
    assert bridge["root_coordinate_span_is_all_E10"] is True
    assert abs(bridge["HNF_basis_determinant"]) == 1
    assert bridge["root_lattice_determinant"] == -1
    assert bridge["root_lattice_gram_smith_diagonal"] == [1] * 10
    assert bridge["graph_only_match_used"] is False

    polytope = document["polytope"]
    cone = polytope["exact_double_description"]
    vinberg = polytope["literal_vinberg"]
    assert (polytope["dimension"], polytope["facet_count"]) == (10, 78)
    assert polytope["classification"] == "finite_volume_noncompact"
    assert (
        cone["finite_vertex_count"],
        cone["ideal_vertex_count"],
        cone["hyperideal_vertex_count"],
    ) == (21_661, 1_797, 0)
    assert cone["all_facets_irredundant"] is True
    assert set(cone["all_facet_face_ranks"]) == {10}
    assert cone["primitive_ray_list_sha256"] == EXPECTED_RAY_HASH
    assert cone["all_extreme_rays_on_one_lorentz_time_sheet"] is True
    assert cone["minimum_strict_interior_to_extreme_ray_pairing"] == 362_134
    assert vinberg["all_pair_compatible_rank_9_cliques"] == 26_769_085
    assert vinberg["all_pair_compatible_rank_10_cliques"] == 40_619_993
    assert vinberg["elliptic_rank_9_edge_diagram_count"] == 166_320
    assert vinberg["elliptic_rank_10_finite_endpoint_diagram_count"] == 21_661
    assert vinberg["edge_endpoint_type_histogram"] == {
        "FF": 69_210,
        "FI": 78_190,
        "II": 18_920,
    }
    assert vinberg["connected_affine_component_count"] == 5_551
    assert vinberg["connected_affine_order_histogram"] == {
        "2": 622,
        "3": 1_890,
        "4": 1_410,
        "5": 804,
        "6": 625,
        "7": 200,
    }
    assert vinberg["cusp_affine_component_count_histogram"] == {
        "2": 480,
        "3": 772,
        "4": 460,
        "5": 75,
        "6": 10,
    }
    assert vinberg["distinct_full_parabolic_neighbourhood_count"] == 1_797
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighbourhood_count"] == 0
    assert polytope["solomon"]["rational_coxeter_euler_characteristic"] == (
        "-527/6480"
    )
    assert polytope["solomon"]["index_in_W_P10"] == EXPECTED_P10_INDEX
    assert polytope["solomon"]["H10_Gauss_Bonnet_volume"] == (
        "1054*pi^5/382725"
    )
    bounded = document["lineage"]["bounded_coxiter_corpus"]
    if bounded["available"]:
        assert bounded["ordinary_H10_78_header_matches"] == []
        assert bounded["doubling_H10_78_header_matches"] == []
    assert document["checks"]["finite_volume"] is True
    assert document["checks"]["noncompact"] is True
    assert document["all_checks_pass"] is True


def test_archived_complete_type_vi_ray_and_edge_ledgers() -> None:
    with gzip.open(RAY_OUTPUT, "rt", encoding="utf-8") as handle:
        ray_header = json.loads(next(handle))
        rays = [json.loads(line) for line in handle]
    assert ray_header["ray_count"] == len(rays) == 23_458
    assert Counter(record["kind"] for record in rays) == {
        "finite": 21_661,
        "ideal": 1_797,
    }
    assert Counter(
        (record["kind"], len(record["active_facets"])) for record in rays
    ) == {
        ("finite", 10): 21_661,
        ("ideal", 11): 480,
        ("ideal", 12): 772,
        ("ideal", 13): 460,
        ("ideal", 14): 75,
        ("ideal", 15): 10,
    }
    canonical = sorted(
        tuple(record["primitive_coordinate_basis_vector"]) for record in rays
    )
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
            endpoint_histogram[
                "".join(sorted(kind[0].upper() for kind in kinds))
            ] += 1
            edge_count += 1
    assert edge_header["edge_count"] == edge_count == len(seen_supports) == 166_320
    assert endpoint_histogram == {"FF": 69_210, "FI": 78_190, "II": 18_920}


def test_type_vi_wall_construction_rebuilds_entrywise() -> None:
    chamber_data, ray_data = _load_templates()
    lattice, walls, construction = _calibration_walls("VI", chamber_data, ray_data)
    document = json.loads(SUMMARY_OUTPUT.read_text(encoding="utf-8"))
    polytope = document["polytope"]
    assert construction == {
        "cm_cusp_count": 57,
        "cm_root_count": 20,
        "wall_order": "b=-K; CM (-2)-roots; cusp lifts E=f-b",
    }
    assert sp.Matrix(polytope["lattice_gram"]) == lattice
    assert polytope["wall_vectors"] == [
        [int(value) for value in wall] for wall in walls
    ]
    wall_matrix = sp.Matrix.hstack(*walls).T
    assert sp.Matrix(polytope["negative_wall_gram"]) == (
        wall_matrix * lattice * wall_matrix.T
    )
    assert [int((wall.T * lattice * wall)[0]) for wall in walls] == (
        [-1] + [-2] * 20 + [-1] * 57
    )
    boundary = walls[0]
    for exceptional in walls[21:]:
        cusp = exceptional + boundary
        assert (cusp.T * lattice * cusp)[0] == 0
        assert (cusp.T * lattice * boundary)[0] == 0


def test_full_exact_type_vi_cdd_diagram_and_vinberg_replay() -> None:
    rebuilt, ray_archive, edge_archive = build_certificate()
    archived = json.loads(SUMMARY_OUTPUT.read_text(encoding="utf-8"))
    assert rebuilt == archived
    assert ray_archive == RAY_OUTPUT.read_bytes()
    assert edge_archive == EDGE_OUTPUT.read_bytes()
