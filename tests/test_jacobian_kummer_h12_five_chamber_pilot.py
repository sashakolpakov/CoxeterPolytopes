from __future__ import annotations

from functools import lru_cache
import json

from scripts.certify_jacobian_kummer_h12_five_chamber_pilot import (
    OUTPUT,
    PILOTS,
    TARGET_FACETS,
    build,
)


@lru_cache(maxsize=1)
def _payload() -> dict:
    return build()


def test_exact_replay_matches_frozen_archive() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["result"] == {
        "decision": "PASS_FREEZE",
        "hyperbolic_dimension": 12,
        "explicit_five_chamber_polytope_count": 5,
        "facet_counts": [752, 761, 770, 779, 788],
        "pairwise_nonisomorphic_ordinary_Coxeter_matrix_count": 5,
        "reflection_index_and_volume_multiplier": 5,
        "finite_volume": True,
        "compact": False,
        "worldwide_priority": "open",
    }


def test_five_explicit_properly_colored_paths_and_formula() -> None:
    rows = _payload()["pilots"]
    assert len(rows) == len(PILOTS) == 5
    assert [row["facet_count"] for row in rows] == list(TARGET_FACETS)
    assert [row["B_edge_occurrence_count"] for row in rows] == list(range(5))
    assert [row["edge_colors_one_based"] for row in rows] == [
        [color + 1 for color in colors] for _b_count, colors in PILOTS
    ]
    for row in rows:
        b_count = row["B_edge_occurrence_count"]
        assert row["tree_shape"] == "P5"
        assert row["proper_edge_coloring"] is True
        assert row["all_four_edge_colors_distinct"] is True
        assert row["chamber_count"] == 5
        assert row["internal_mirror_count"] == 4
        assert row["orientation_cancellation_internal_mirror_count"] == 4
        assert row["facet_count"] == 752 + 9 * b_count
        assert row[
            "expected_facet_count_172_plus_145_times_4_plus_9B"
        ] == 752 + 9 * b_count
        assert row["marked_chamber_adjacency"] == [
            [index + 1, index + 2, color]
            for index, color in enumerate(row["edge_colors_one_based"])
        ]


def test_ordinary_coxeter_integrality_interior_and_smith_gates() -> None:
    rows = _payload()["pilots"]
    boundary_hashes = set()
    coxeter_hashes = set()
    for row in rows:
        facets = row["facet_count"]
        assert row["primitive_integral_reflection_root_count"] == facets
        assert row["strict_interior_square"] > 0
        assert row["minimum_strict_interior_facet_pairing"] == 5
        assert row["fixed_boundary_core_Smith_diagonal"] == [1] * 13
        assert row["full_linear_rank"] == 13
        histogram = row["ordinary_Coxeter_relation_histogram"]
        assert set(histogram) == {"m2", "m3", "m4", "m6", "infinity"}
        assert sum(histogram.values()) == facets * (facets - 1) // 2
        assert histogram["m6"] == 0
        boundary_hashes.add(row["marked_boundary_root_sha256"])
        coxeter_hashes.add(row["marked_ordinary_Coxeter_matrix_sha256"])
    assert len(boundary_hashes) == len(coxeter_hashes) == 5
    # Distinct matrix orders alone are an unambiguous nonisomorphism proof.
    assert len({row["facet_count"] for row in rows}) == 5


def test_Allcock_volume_lineage_and_carrier_scope() -> None:
    payload = _payload()
    stock = payload["seed_and_Allcock_stock"]["stock"]
    assert stock["wall_count"] == stock["doubling_wall_count"] == 112
    assert stock["A_wall_count_W61_through_W76"] == 16
    assert stock["B_wall_count_W77_through_W172"] == 96
    assert stock["pairwise_disjoint_pair_count"] == 6216
    assert stock["all_doubling"] is True
    assert stock["all_pairwise_disjoint_in_Allcock_sense"] is True

    lineage = payload["lineage"]
    targets = set(TARGET_FACETS)
    assert targets.isdisjoint(lineage["pinned_project_predecessor_facet_counts"])
    assert targets.isdisjoint(lineage["pinned_project_ledger_H12_facet_counts"])
    assert targets.isdisjoint(lineage["pinned_CoxIter_H12_facet_counts"])
    assert lineage["target_counts_absent_from_both_project_sets"] is True
    assert lineage["target_headers_absent_from_pinned_CoxIter_lanes"] is True
    assert lineage["worldwide_priority"] == "open"

    scope = payload["carrier_scope"]
    assert scope["same_full_Picard_K3_family"] is True
    assert scope["fixed_index_one_divisor_lattice_core"] is True
    assert scope["weak_divisor_reflection_carrier_only"] is True
    assert scope["full_nef_realization_claimed"] is False
    assert scope["surface_automorphism_realization_claimed"] is False
    boundary = payload["claim_boundary"]
    assert boundary["full_five_vertex_coloring_or_shape_census_claimed"] is False
    assert boundary["general_facet_formula_for_every_five_vertex_tree_claimed"] is False
    assert boundary["worldwide_novelty_claimed"] is False
