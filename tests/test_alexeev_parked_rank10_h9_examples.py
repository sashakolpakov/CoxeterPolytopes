from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_alexeev_parked_rank10_h9_examples import OUTPUT, build


ROOT = Path(__file__).resolve().parents[1]


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_five_parked_rank10_h9_examples_rebuild_exactly() -> None:
    archived = _archive()
    rebuilt, _ = build()
    assert rebuilt == archived
    assert archived["summary"] == {
        "compact_count": 0,
        "facet_count_histogram": {"12": 1, "13": 1, "15": 2, "16": 1},
        "finite_volume_noncompact_count": 5,
        "new_project_ledger_diagram_count": 5,
        "total_bad_edge_count": 0,
        "total_bad_parabolic_neighborhood_count": 0,
    }


def test_exact_volume_vinberg_and_lattice_gates() -> None:
    expected = {
        "A14_8_0_C01": (12, 44, 5, 260, 10, []),
        "A14_8_1_C02": (16, 165, 23, 1112, 55, [7, 15, 16]),
        "A14_10_1_C01": (13, 70, 7, 393, 14, [3, 5, 11]),
        "A14_10_1_C02": (15, 123, 15, 778, 35, [1, 3, 6, 15]),
        "A16_14_1_C01": (15, 100, 9, 588, 23, [1, 8, 14, 15]),
    }
    for case in _archive()["cases"]:
        facets, finite, ideal, edges, affine, doubling = expected[case["id"]]
        cone = case["exact_cone"]
        exact = cone["exact_double_description"]
        assert cone["decision"]["classification"] == "finite_volume_noncompact"
        assert case["section"]["exact_reduction"][
            "active_hyperbolic_facet_count"
        ] == facets
        assert exact["finite_vertex_count"] == finite
        assert exact["ideal_vertex_count"] == ideal
        assert exact["hyperideal_vertex_count"] == 0
        assert case["time_orientation"]["all_extreme_rays_same_time_cone"]
        assert case["literal_vinberg"]["elliptic_edge_count"] == edges
        assert case["literal_vinberg"][
            "connected_affine_component_count"
        ] == affine
        assert case["literal_vinberg"]["bad_edge_count"] == 0
        assert case["literal_vinberg"]["bad_parabolic_neighborhood_count"] == 0
        assert case["doubling_walls_one_based"] == doubling
        assert case["section"]["restricted_mirrors"][
            "all_active_reflections_integral_on_saturated_fixed_lattice"
        ]


def test_denominator_and_nonduplication_scope() -> None:
    archived = _archive()
    assert len(archived["source"]["complete_parked_rank10_subgroup_rows"]) == 6
    assert archived["source"]["previously_frozen_row"] == [[14, 8, 1], 1]
    assert not any(archived["pairwise_relation_diagram_isomorphisms"].values())
    assert not any(
        archived["matches_with_prior_project_H9_factory_or_cube_diagrams"].values()
    )
    assert not archived["scope"]["surface_automorphism_or_nef_claimed"]
    assert not archived["scope"]["worldwide_novelty_claimed"]
