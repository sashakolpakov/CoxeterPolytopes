from __future__ import annotations

from collections import Counter
import json

import scripts.certify_jacobian_kummer_h12_small_subtree_census as base
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star


def _archive() -> dict:
    return json.loads(p5.OUTPUT.read_text(encoding="utf-8"))


def test_frozen_scope_and_exact_binary_shard_denominators() -> None:
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    denominator = payload["coloring_denominator"]
    assert denominator["shape"] == "P5"
    assert tuple(denominator["selected_binary_class_word_shards"]) == p5.PRODUCTION_PATTERNS
    assert denominator["selected_shard_count"] == 7
    assert denominator["selected_raw_unrooted_proper_colorings"] == 7_562_616
    assert denominator["selected_seed_and_reversal_orbit_count"] == 20_693
    assert denominator["all_ten_pattern_orbit_count"] == 201_747
    assert denominator["remaining_disjoint_patterns"] == ["ABBB", "BABB", "BBBB"]
    assert denominator["remaining_disjoint_orbit_count"] == 181_054
    for pattern in p5.PRODUCTION_PATTERNS:
        row = denominator["per_pattern"][pattern]
        assert row["raw_unrooted_proper_colorings"] == p5._raw_pattern_count(pattern)
        assert row["orbit_count"] == p5.EXPECTED_PATTERN_ORBITS[pattern]
        assert row["explicit_orbit_cover_count"] == row["orbit_count"]
        assert row["explicit_orbit_size_sum"] == row["raw_unrooted_proper_colorings"]


def test_complete_type_count_formula_and_pairwise_invariants() -> None:
    payload = _archive()
    theorem = payload["theorem"]
    rows = payload["orbit_representatives"]
    assert theorem["ordinary_Coxeter_type_count"] == len(rows) == 20_693
    assert theorem["pairwise_nonisomorphic_within_selected_layer"] is True
    assert theorem["reflection_index_and_volume_multiple"] == 5
    assert theorem["finite_volume"] is True
    assert theorem["compact"] is False
    assert theorem["nonadjacent_color_repetition_type_histogram"] == {
        "0": 19_677,
        "1": 1_005,
        "2": 11,
    }
    assert Counter(row["edge_class_word"] for row in rows) == Counter(
        {
            pattern: p5.EXPECTED_PATTERN_ORBITS[pattern]
            for pattern in p5.PRODUCTION_PATTERNS
        }
    )
    type_hashes = set()
    for row in rows:
        assert row["tree_shape"] == "P5"
        assert row["proper_edge_coloring"] is True
        assert row["edge_class_word"] in p5.PRODUCTION_PATTERNS
        assert row["B_edge_occurrence_count"] <= 2
        assert row["facet_count"] == 752 + 9 * row["B_edge_occurrence_count"]
        assert row["expected_facet_count_752_plus_9B"] == row["facet_count"]
        type_hashes.add(row["ordinary_Coxeter_type_invariant_sha256"])
    assert len(type_hashes) == len(rows)


def test_every_archived_root_interior_smith_and_ordinary_gate() -> None:
    for row in _archive()["orbit_representatives"]:
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


def test_exact_burnside_counts_for_all_ten_resumable_shards() -> None:
    actions = star._restricted_actions()
    assert all(
        all((color < 16) == (action[color] < 16) for color in range(112))
        for action in actions
    )
    counts = {
        pattern: p5._burnside_pattern_orbits(pattern, actions)["orbit_count"]
        for pattern in p5.ALL_PATTERNS
    }
    assert counts == p5.EXPECTED_PATTERN_ORBITS
    assert sum(counts.values()) == 201_747


def test_one_boundary_replay_per_completed_binary_shard() -> None:
    rows = _archive()["orbit_representatives"]
    selected = {}
    for row in rows:
        selected.setdefault(row["edge_class_word"], row)
    assert set(selected) == set(p5.PRODUCTION_PATTERNS)
    _carrier, gram, roots, interior, _bounded = base._load_inputs()
    engine = p5.ExactP5BoundaryEngine(gram, roots, interior)
    for pattern, row in selected.items():
        stock_colors = tuple(
            color - p5.STOCK_START - 1 for color in row["edge_colors_one_based"]
        )
        boundary, boundary_data = engine.boundary(stock_colors)
        _labels, degree_key, ordinary_data = engine.ordinary_data(boundary)
        assert boundary_data["facet_count"] == 752 + 9 * pattern.count("B")
        assert boundary_data["marked_boundary_root_sha256"] == row[
            "marked_boundary_root_sha256"
        ]
        assert ordinary_data["ordinary_degree_profile_multiset_sha256"] == row[
            "ordinary_degree_profile_multiset_sha256"
        ]
        assert degree_key


def test_pilot_recovery_and_claim_boundary() -> None:
    payload = _archive()
    recovery = payload["pilot_recovery"]
    assert recovery["eligible_B_0_1_2_pilot_count"] == 3
    assert recovery[
        "project_positive_increment_over_the_five_path_pilot_within_this_layer"
    ] == 20_690
    assert {row["B_edge_occurrence_count"] for row in recovery["rows"]} == {0, 1, 2}
    assert all(row["same_seed_action_and_path_reversal_orbit"] for row in recovery["rows"])
    scope = payload["lineage_and_claim_boundary"]
    assert scope["complete_B_at_most_2_layer_claimed"] is True
    assert scope["complete_P5_census_claimed"] is False
    assert scope["cross_star_pure_exact_nonisomorphism_claimed"] is False
    assert scope["worldwide_priority"] == "open"
    carrier = payload["Allcock_and_carrier"]
    assert carrier["same_full_Picard_K3_family"] is True
    assert carrier["weak_divisor_reflection_carrier_only"] is True
    assert carrier["full_nef_effectivity_or_surface_action_claim"] is False
