from __future__ import annotations

from collections import Counter
import gzip
import json

import sympy as sp

from scripts.certify_shimada_disc16_d5_full_subgroup_fishing import (
    OUTPUT as CENSUS_OUTPUT,
)
from scripts.certify_shimada_disc16_d5_h10_30_exact_cone import (
    SUMMARY_OUTPUT as CONE_OUTPUT,
    _rebuild_exact_cone,
    load_ray_archive,
)
from scripts.certify_shimada_disc16_d5_h10_30_vinberg_faces import (
    EDGE_OUTPUT,
    OUTPUT as VINBERG_OUTPUT,
)


def _load(path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_complete_d20_subgroup_denominator_and_above_h9_classes() -> None:
    census = _load(CENSUS_OUTPUT)
    stabilizer = census["weighted_wall_stabilizer"]
    subgroups = census["subgroup_census"]

    assert stabilizer["order"] == 20
    assert stabilizer["element_order_histogram"] == {
        "1": 1,
        "2": 11,
        "5": 4,
        "10": 4,
    }
    assert stabilizer["recognition"]["recognized_structure"] == (
        "dihedral group of order 20 (D_10)"
    )
    assert stabilizer[
        "archived_lattice_stabilizer_equals_full_exact_weighted_wall_graph_automorphism_group"
    ] is True
    assert subgroups["subgroup_count"] == 22
    assert subgroups["conjugacy_class_count"] == 10
    assert subgroups["order_dimension_histogram"] == {
        "order_1_H19": 1,
        "order_2_H10": 5,
        "order_2_H11": 5,
        "order_2_H14": 1,
        "order_4_H8": 5,
        "order_5_H3": 1,
        "order_10_H2": 2,
        "order_10_H3": 1,
        "order_20_H2": 1,
    }
    assert [
        (record["order"], record["hyperbolic_dimension"], record["class_size"])
        for record in census["above_H9"]["conjugacy_classes"]
    ] == [(1, 19, 1), (2, 14, 1), (2, 10, 5), (2, 11, 5)]


def test_all_above_h9_section_gates() -> None:
    screens = _load(CENSUS_OUTPUT)["above_H9"]["complete_section_screens"]
    assert {
        key: (
            screen.get("facet_count", screen.get("source_certified_irredundant_facet_count")),
            screen["pair_gate"]["non_coxeter_pair_count"],
        )
        for key, screen in screens.items()
    } == {
        "class_01_H19": (165, 2625),
        "class_02_H14": (60, 0),
        "class_03_H10": (30, 0),
        "class_04_H11": (31, 109),
    }

    h10 = screens["class_03_H10"]
    assert h10["wall_orbit_count"] == 86
    assert h10["redundant_orbit_count"] == 56
    assert h10["pair_gate"]["pair_count"] == 435
    assert h10["pair_gate"]["coxeter_edge_counts"] == {
        "2": 276,
        "3": 15,
        "4": 30,
        "infinity": 62,
        "ultraparallel": 52,
    }
    retained = set(h10["facet_orbit_indices"])
    assert len(retained) == 30
    assert all(
        all(term["facet_orbit_index"] in retained for term in record["terms"])
        for record in h10["redundancy_certificates"]
    )
    assert all(
        record["method"] == "orthogonal_projection_of_published_weyl_vector"
        for record in h10["facet_witnesses"]
    )

    h11 = screens["class_04_H11"]
    assert h11["wall_orbit_count"] == 91
    assert h11["redundant_orbit_count"] == 60
    assert sum(
        record["method"].endswith("exact_tangent_perturbation")
        for record in h11["facet_witnesses"]
    ) == 2


def test_gmp_rays_replay_all_exact_incidence_and_lorentz_types() -> None:
    cone_summary = _load(CONE_OUTPUT)
    exact = cone_summary["exact_double_description"]
    assert exact["extreme_ray_count"] == 1540
    assert exact["finite_vertex_count"] == 1333
    assert exact["ideal_vertex_count"] == 207
    assert exact["hyperideal_vertex_count"] == 0
    assert exact["active_facet_count_histogram"] == {
        "finite_10": 1333,
        "ideal_10": 48,
        "ideal_11": 98,
        "ideal_12": 52,
        "ideal_13": 9,
    }
    assert cone_summary["decision"] == {
        "classification": "finite_volume_noncompact",
        "compact": False,
        "criterion": "complete pointed projective cone has no spacelike/hyperideal extreme ray",
        "finite_volume": True,
    }

    cone = _rebuild_exact_cone()
    inequality: sp.Matrix = cone["inequality_matrix"]
    basis_gram: sp.Matrix = cone["basis_gram"]
    header, rays = load_ray_archive()
    assert header["ray_count"] == len(rays) == 1540
    replay_types = Counter()
    replay_active = Counter()
    for ray, archived_active, archived_square, archived_kind in rays:
        vector = sp.Matrix(ray)
        slacks = inequality * vector
        assert all(value >= 0 for value in slacks)
        active = tuple(index + 1 for index, value in enumerate(slacks) if value == 0)
        square = sp.Rational((vector.T * basis_gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        assert active == tuple(archived_active)
        assert str(square) == archived_square
        assert kind == archived_kind
        replay_types[kind] += 1
        replay_active[(kind, len(active))] += 1
    assert replay_types == {"finite": 1333, "ideal": 207}
    assert replay_active == {
        ("finite", 10): 1333,
        ("ideal", 10): 48,
        ("ideal", 11): 98,
        ("ideal", 12): 52,
        ("ideal", 13): 9,
    }


def test_literal_vinberg_edges_endpoints_and_full_affine_neighborhoods() -> None:
    vinberg = _load(VINBERG_OUTPUT)
    diagram = vinberg["diagram_enumeration"]
    endpoint = vinberg["literal_endpoint_audit"]
    assert diagram["pair_compatible_clique_count_by_order"] == {
        "1": 30,
        "2": 321,
        "3": 1929,
        "4": 7706,
        "5": 22293,
        "6": 49110,
        "7": 85038,
        "8": 118152,
        "9": 133463,
        "10": 123508,
        "11": 93930,
        "12": 58638,
    }
    assert diagram["elliptic_rank_9_edge_diagram_count"] == 10030
    assert diagram["elliptic_rank_10_finite_endpoint_diagram_count"] == 1333
    assert diagram["connected_affine_component_count"] == 436
    assert endpoint["edge_endpoint_count_histogram"] == {"2": 10030}
    assert endpoint["edge_endpoint_type_histogram"] == {
        "ff": 4423,
        "fi": 4484,
        "ii": 1123,
    }
    assert endpoint["finite_vertex_count"] == 1333
    assert endpoint["ideal_vertex_count"] == 207
    assert endpoint["bad_edge_count"] == 0
    assert endpoint["bad_parabolic_neighborhood_count"] == 0
    assert vinberg["decision"]["classification"] == "finite_volume_noncompact"

    ray_header, rays = load_ray_archive()
    del ray_header
    with gzip.open(EDGE_OUTPUT, "rt", encoding="utf-8") as stream:
        edge_header = json.loads(next(stream))
        edges = tuple(json.loads(line) for line in stream)
    assert edge_header["edge_count"] == len(edges) == 10030
    seen = set()
    endpoint_types = Counter()
    for edge_facets, endpoint_indices, archived_types, finite_components in edges:
        del finite_components
        left, right = (index - 1 for index in endpoint_indices)
        intersection = sorted(set(rays[left][1]).intersection(rays[right][1]))
        assert intersection == edge_facets
        assert archived_types == [rays[left][3], rays[right][3]]
        assert len(edge_facets) == 9
        assert tuple(edge_facets) not in seen
        seen.add(tuple(edge_facets))
        endpoint_types["".join(sorted(kind[0] for kind in archived_types))] += 1
    assert endpoint_types == {"ff": 4423, "fi": 4484, "ii": 1123}
