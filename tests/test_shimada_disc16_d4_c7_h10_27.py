from __future__ import annotations

from collections import Counter
import gzip
import json

import sympy as sp

from scripts.certify_shimada_disc16_d4_c7_h10_27_exact_cone import (
    CASE_INPUT,
    SUMMARY_OUTPUT as CONE_OUTPUT,
    _rebuild_exact_cone,
    load_ray_archive,
)
from scripts.certify_shimada_disc16_d4_c7_h10_27_vinberg_faces import (
    EDGE_OUTPUT,
    OUTPUT as VINBERG_OUTPUT,
)


def _load(path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_archived_d4_c7_exact_section_and_coxeter_gate() -> None:
    case = _load(CASE_INPUT)
    assert case["id"] == "D0004_C007"
    assert case["case"]["archived_stabilizer_order"] == 48
    assert case["case"]["subgroup_order"] == 4
    assert case["case"]["subgroup_conjugacy_class_size"] == 3
    assert case["saturated_fixed_lattice"]["inertia"] == [1, 10, 0]
    assert case["saturated_fixed_lattice"]["determinant"] == 512
    assert case["folded_walls"]["wall_orbit_count"] == 62
    reduction = case["exact_reduction"]
    assert reduction["hyperbolic_facet_count"] == 27
    assert reduction["conic_redundant_count"] == 35
    assert reduction["unresolved_conic_facet_count"] == 0
    assert reduction["causally_inactive_conic_facet_count"] == 0
    gate = case["coxeter_gate"]
    assert gate["all_pairs_coxeter"] is True
    assert gate["non_coxeter_pair_count"] == 0
    assert gate["coxeter_edge_counts"] == {
        "2": 235,
        "3": 15,
        "4": 14,
        "infinity": 48,
        "ultraparallel": 39,
    }


def test_gmp_cone_rays_replay_exact_incidence_and_lorentz_types() -> None:
    summary = _load(CONE_OUTPUT)
    exact = summary["exact_double_description"]
    assert exact["extreme_ray_count"] == 962
    assert exact["finite_vertex_count"] == 895
    assert exact["ideal_vertex_count"] == 67
    assert exact["hyperideal_vertex_count"] == 0
    assert exact["active_facet_count_histogram"] == {
        "finite_10": 895,
        "ideal_11": 6,
        "ideal_12": 58,
        "ideal_13": 3,
    }
    assert summary["decision"] == {
        "classification": "finite_volume_noncompact",
        "compact": False,
        "criterion": "complete pointed projective cone has no spacelike/hyperideal extreme ray",
        "finite_volume": True,
    }

    cone = _rebuild_exact_cone()
    inequality: sp.Matrix = cone["inequality_matrix"]
    fixed_gram: sp.Matrix = cone["fixed_gram"]
    header, rays = load_ray_archive()
    assert header["ray_count"] == len(rays) == 962
    replay_types = Counter()
    replay_active = Counter()
    for ray, archived_active, archived_square, archived_kind in rays:
        vector = sp.Matrix(ray)
        slacks = inequality * vector
        assert all(value >= 0 for value in slacks)
        active = tuple(index + 1 for index, value in enumerate(slacks) if value == 0)
        square = sp.Rational((vector.T * fixed_gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        assert active == tuple(archived_active)
        assert str(square) == archived_square
        assert kind == archived_kind
        replay_types[kind] += 1
        replay_active[(kind, len(active))] += 1
    assert replay_types == {"finite": 895, "ideal": 67}
    assert replay_active == {
        ("finite", 10): 895,
        ("ideal", 11): 6,
        ("ideal", 12): 58,
        ("ideal", 13): 3,
    }


def test_literal_vinberg_edges_endpoints_and_affine_neighborhoods() -> None:
    vinberg = _load(VINBERG_OUTPUT)
    diagram = vinberg["diagram_enumeration"]
    endpoint = vinberg["literal_endpoint_audit"]
    assert diagram["pair_compatible_clique_count_by_order"] == {
        "1": 27,
        "2": 264,
        "3": 1446,
        "4": 5238,
        "5": 13683,
        "6": 27177,
        "7": 42516,
        "8": 53658,
        "9": 55451,
        "10": 47223,
        "11": 33078,
        "12": 18870,
    }
    assert diagram["elliptic_rank_9_edge_diagram_count"] == 6221
    assert diagram["elliptic_rank_10_finite_endpoint_diagram_count"] == 895
    assert diagram["connected_affine_component_count"] == 198
    assert endpoint["edge_endpoint_count_histogram"] == {"2": 6221}
    assert endpoint["edge_endpoint_type_histogram"] == {
        "ff": 3161,
        "fi": 2628,
        "ii": 432,
    }
    assert endpoint["finite_vertex_count"] == 895
    assert endpoint["ideal_vertex_count"] == 67
    assert endpoint["bad_edge_count"] == 0
    assert endpoint["bad_parabolic_neighborhood_count"] == 0
    assert vinberg["decision"]["classification"] == "finite_volume_noncompact"

    _, rays = load_ray_archive()
    with gzip.open(EDGE_OUTPUT, "rt", encoding="utf-8") as stream:
        edge_header = json.loads(next(stream))
        edges = tuple(json.loads(line) for line in stream)
    assert edge_header["edge_count"] == len(edges) == 6221
    seen = set()
    endpoint_types = Counter()
    for edge_facets, endpoint_indices, archived_types, finite_components in edges:
        assert finite_components
        left, right = (index - 1 for index in endpoint_indices)
        intersection = sorted(set(rays[left][1]).intersection(rays[right][1]))
        assert intersection == edge_facets
        assert archived_types == [rays[left][3], rays[right][3]]
        assert len(edge_facets) == 9
        assert tuple(edge_facets) not in seen
        seen.add(tuple(edge_facets))
        endpoint_types["".join(sorted(kind[0] for kind in archived_types))] += 1
    assert endpoint_types == {"ff": 3161, "fi": 2628, "ii": 432}
