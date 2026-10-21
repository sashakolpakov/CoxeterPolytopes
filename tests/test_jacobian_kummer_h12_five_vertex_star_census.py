from __future__ import annotations

from collections import Counter
import json
import math

import scripts.certify_jacobian_kummer_h12_small_subtree_census as base
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star


def _archive() -> dict:
    return json.loads(star.OUTPUT.read_text(encoding="utf-8"))


def test_frozen_archive_and_complete_K14_denominator() -> None:
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    denominator = payload["coloring_denominator"]
    assert denominator["shape"] == "K1,4"
    assert denominator["raw_unrooted_proper_colorings"] == math.comb(112, 4) == 6_210_820
    assert denominator["seed_action_order"] == 384
    assert denominator["Burnside_fixed_subset_sum"] == 6_478_848
    assert denominator["Burnside_orbit_count"] == 16_872
    assert denominator["explicit_orbit_cover_count"] == 16_872
    assert denominator["explicit_orbit_size_sum"] == 6_210_820
    assert sum(
        int(size) * count
        for size, count in denominator["seed_orbit_size_histogram"].items()
    ) == 6_210_820


def test_complete_type_count_facet_formula_and_raw_histogram() -> None:
    payload = _archive()
    theorem = payload["theorem"]
    rows = payload["orbit_representatives"]
    assert theorem["scope"] == "complete K1,4 five-vertex shape only"
    assert theorem["ordinary_Coxeter_type_count"] == len(rows) == 16_872
    assert theorem["pairwise_nonisomorphic"] is True
    assert theorem["reflection_index_and_volume_multiple"] == 5
    assert theorem["finite_volume"] is True
    assert theorem["compact"] is False
    assert Counter(row["facet_count"] for row in rows) == Counter(
        {int(key): value for key, value in theorem["type_facet_histogram"].items()}
    )
    expected_raw = {
        752 + 9 * b_count: math.comb(16, 4 - b_count) * math.comb(96, b_count)
        for b_count in range(5)
    }
    assert {int(key): value for key, value in theorem[
        "raw_unrooted_coloring_facet_histogram"
    ].items()} == expected_raw
    for row in rows:
        assert row["tree_shape"] == "K1,4"
        assert row["proper_edge_coloring"] is True
        assert row["all_four_incident_edge_colors_distinct"] is True
        assert row["facet_count"] == 752 + 9 * row["B_edge_occurrence_count"]
        assert row["expected_facet_count_752_plus_9B"] == row["facet_count"]


def test_every_archived_root_interior_smith_and_ordinary_gate() -> None:
    rows = _archive()["orbit_representatives"]
    type_hashes = set()
    for row in rows:
        facets = row["facet_count"]
        assert row["chamber_count"] == 5
        assert row["internal_mirror_count"] == 4
        assert row["strict_interior_square"] > 0
        assert row["minimum_strict_interior_facet_pairing"] > 0
        assert row["primitive_integral_reflection_root_count"] == facets
        assert row["fixed_boundary_core_Smith_diagonal"] == [1] * 13
        assert row["full_linear_rank"] == 13
        histogram = row["ordinary_Coxeter_relation_histogram"]
        assert set(histogram) == {"m2", "m3", "m4", "m6", "infinity"}
        assert sum(histogram.values()) == facets * (facets - 1) // 2
        type_hashes.add(row["ordinary_Coxeter_type_invariant_sha256"])
    assert len(type_hashes) == len(rows) == 16_872


def test_one_exact_boundary_replay_for_each_B_count() -> None:
    rows = _archive()["orbit_representatives"]
    selected = {}
    for row in rows:
        selected.setdefault(row["B_edge_occurrence_count"], row)
    assert set(selected) == set(range(5))
    _carrier, gram, roots, interior, _bounded = base._load_inputs()
    engine = star.ExactStarBoundaryEngine(gram, roots, interior)
    for b_count, row in selected.items():
        stock_colors = tuple(
            color - 1 - star.STOCK_START for color in row["edge_colors_one_based"]
        )
        boundary, boundary_data = engine.boundary(stock_colors)
        _labels, degree_key, ordinary_data = engine.ordinary_data(boundary)
        assert boundary_data["facet_count"] == 752 + 9 * b_count
        assert boundary_data["marked_boundary_root_sha256"] == row[
            "marked_boundary_root_sha256"
        ]
        assert ordinary_data["ordinary_degree_profile_multiset_sha256"] == row[
            "ordinary_degree_profile_multiset_sha256"
        ]
        assert degree_key


def test_claim_boundary_remains_shape_local_and_weak_carrier() -> None:
    payload = _archive()
    scope = payload["lineage_and_claim_boundary"]
    assert scope["P5_complete_census_claimed"] is False
    assert scope["degree_three_fork_complete_census_claimed"] is False
    assert scope["all_five_vertex_shapes_complete_census_claimed"] is False
    assert scope["worldwide_priority"] == "open"
    carrier = payload["Allcock_and_carrier"]
    assert carrier["same_full_Picard_K3_family"] is True
    assert carrier["fixed_index_one_divisor_lattice_core"] is True
    assert carrier["weak_divisor_reflection_carrier_only"] is True
    assert carrier["full_nef_effectivity_or_surface_action_claim"] is False
