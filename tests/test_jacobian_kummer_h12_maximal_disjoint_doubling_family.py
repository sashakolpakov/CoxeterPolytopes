from __future__ import annotations

from functools import lru_cache
import json

from scripts.certify_jacobian_kummer_h12_maximal_disjoint_doubling_family import (
    OUTPUT,
    PINNED_SHA256,
    build,
)


@lru_cache(maxsize=1)
def _payload() -> dict:
    return build()


def test_archive_is_an_exact_independent_raw_source_rebuild() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["pinned_inputs"] == PINNED_SHA256
    independence = payload["independence"]
    assert independence["published_source_reparsed"] is True
    assert independence["fixed_lattice_and_172_roots_reconstructed"] is True
    assert independence["other_JK_H12_scout_or_certifier_imported"] is False
    assert independence["graph_library_used"] is False


def test_exact_124_wall_doubling_stock() -> None:
    stock = _payload()["doubling_stock"]
    assert stock["doubling_wall_count"] == 124
    assert stock["doubling_walls_one_based"] == [
        31,
        32,
        33,
        36,
        39,
        42,
        48,
        51,
        56,
        *range(58, 173),
    ]
    assert stock["whole_seed_relation_histogram"] == {
        "m2": 2820,
        "m3": 144,
        "m4": 800,
        "parallel": 1750,
        "ultraparallel": 9192,
    }
    assert stock["doubling_stock_relation_histogram"] == {
        "m2": 264,
        "parallel": 1138,
        "ultraparallel": 6224,
    }
    assert all(stock["checks"].values())


def test_unique_maximal_112_wall_disjoint_stock() -> None:
    graph = _payload()["maximal_disjoint_graph_theorem"]
    assert graph["relation"] == "parallel_or_ultraparallel"
    assert graph["vertex_count"] == 124
    assert graph["edge_count"] == 7362
    assert graph["degree_histogram"] == {"99": 12, "120": 16, "121": 96}
    assert graph["maximal_clique_count"] == 100
    assert graph["maximal_clique_size_histogram"] == {
        "36": 9,
        "55": 36,
        "74": 42,
        "93": 12,
        "112": 1,
    }
    assert graph["maximum_size"] == 112
    assert graph["unique_maximum_one_based"] == list(range(61, 173))
    assert graph["unique_maximum_relation_histogram"] == {
        "parallel": 664,
        "ultraparallel": 5552,
    }
    assert graph["conflict_matching_size"] == 12
    assert graph["matching_upper_bound"] == 112
    assert set(
        graph["maximum_clique_size_containing_each_excluded_vertex"].values()
    ) == {93}
    assert len(graph["all_maximal_cliques_one_based"]) == 100
    assert all(graph["checks"].values())


def test_unique_maximal_96_wall_strict_ultraparallel_stock() -> None:
    graph = _payload()["maximal_ultraparallel_graph_theorem"]
    assert graph["relation"] == "strictly_ultraparallel"
    assert graph["vertex_count"] == 124
    assert graph["edge_count"] == 6224
    assert graph["degree_histogram"] == {"56": 12, "64": 16, "112": 96}
    assert graph["maximal_clique_count"] == 61
    assert graph["maximal_clique_size_histogram"] == {
        "32": 32,
        "57": 12,
        "61": 16,
        "96": 1,
    }
    assert graph["maximum_size"] == 96
    assert graph["unique_maximum_one_based"] == list(range(77, 173))
    assert graph["unique_maximum_relation_histogram"] == {
        "ultraparallel": 4560
    }
    assert graph["conflict_matching_size"] == 28
    assert graph["matching_upper_bound"] == 96
    assert set(
        graph["maximum_clique_size_containing_each_excluded_vertex"].values()
    ) == {57, 61}
    assert len(graph["all_maximal_cliques_one_based"]) == 61
    assert all(graph["checks"].values())


def test_allcock_finite_subtrees_and_full_Picard_scope() -> None:
    payload = _payload()
    theorem = payload["theorem"]
    assert theorem["hyperbolic_dimension"] == 12
    assert theorem["maximal_disjoint_doubling_stock_size"] == 112
    assert theorem["maximal_strictly_ultraparallel_stock_size"] == 96
    assert theorem["allcock_4_1_applicability"]["applies"] is True
    assert theorem["allcock_4_3_applicability"]["applies"] is True
    assert "|V(T)|" in theorem["exact_index_and_volume"]

    carrier = theorem["carrier"]
    assert carrier["seed_classification"] == "finite_volume_noncompact"
    assert carrier["seed_cusp_count"] == 6214
    assert carrier["carrier_lattice"] == "U(2)+D5(-1)^2+<-4>"
    assert carrier["carrier_rank"] == 13
    assert carrier["carrier_moduli_dimension"] == 7
    assert carrier["disjoint_112_boundary_core"]["walls_one_based"] == list(
        range(1, 61)
    )
    assert carrier["ultraparallel_96_boundary_core"]["walls_one_based"] == list(
        range(1, 77)
    )
    assert carrier["disjoint_112_boundary_core"]["Smith_diagonal"] == [1] * 13
    assert carrier["ultraparallel_96_boundary_core"]["Smith_diagonal"] == [1] * 13
    assert all(carrier["checks"].values())

    scope = payload["lineage_and_scope"]
    assert scope["worldwide_priority"] == "open"
    assert scope["worldwide_novelty_claimed"] is False
    assert scope["new_dimension_claimed"] is False
    assert scope["new_arithmetic_commensurability_class_claimed"] is False
    assert scope["full_nef_claimed"] is False
    assert scope["surface_automorphism_claimed"] is False
