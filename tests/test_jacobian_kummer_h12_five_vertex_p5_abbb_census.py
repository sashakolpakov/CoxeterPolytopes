from __future__ import annotations

from collections import Counter
import json

import scripts.certify_jacobian_kummer_h12_small_subtree_census as seed
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_abbb_census as abbb
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as frozen
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star


def _archive() -> dict:
    return json.loads(abbb.OUTPUT.read_text(encoding="utf-8"))


def test_frozen_predecessor_and_complete_ABBB_denominator() -> None:
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["pinned_inputs"]["P5_B_le_2_builder"] == abbb.PINNED_FROZEN_SHA256[
        "P5_B_le_2_builder"
    ]
    assert payload["pinned_inputs"]["P5_B_le_2_archive"] == abbb.PINNED_FROZEN_SHA256[
        "P5_B_le_2_archive"
    ]
    denominator = payload["coloring_denominator"]
    assert denominator["selected_binary_class_word_shards"] == ["ABBB"]
    assert denominator["selected_raw_unrooted_proper_colorings"] == 13_862_400
    assert denominator["selected_seed_and_reversal_orbit_count"] == 36_247
    row = denominator["per_pattern"]["ABBB"]
    assert row["raw_unrooted_proper_colorings"] == 13_862_400
    assert row["orbit_count"] == row["explicit_orbit_cover_count"] == 36_247
    assert row["explicit_orbit_size_sum"] == 13_862_400
    assert "remaining_disjoint_patterns" not in denominator
    assert "remaining_disjoint_orbit_count" not in denominator
    assert denominator["complement_of_selected_shard_patterns"] == [
        "AAAA",
        "AAAB",
        "AABA",
        "AABB",
        "ABAB",
        "ABBA",
        "BAAB",
        "BABB",
        "BBBB",
    ]
    assert denominator["complement_of_selected_shard_orbit_count"] == 165_500
    assert denominator["globally_remaining_binary_class_word_shards"] == [
        "BABB",
        "BBBB",
    ]
    assert denominator["globally_remaining_seed_and_reversal_orbits"] == 144_807
    checks = payload["checks"]
    assert "production_scope_is_exactly_B_at_most_2" not in checks
    assert "full_P5_census_not_claimed_for_production_scope" not in checks
    assert checks["predecessor_B_le_2_scope_pinned_and_green"] is True
    assert checks["production_scope_is_exactly_ABBB"] is True
    assert checks["complete_P5_census_not_claimed"] is True


def test_complete_ABBB_type_count_and_pairwise_invariants() -> None:
    payload = _archive()
    theorem = payload["theorem"]
    rows = payload["orbit_representatives"]
    assert theorem["scope"] == "complete ABBB five-vertex P5 class-word shard"
    assert theorem["ordinary_Coxeter_type_count"] == len(rows) == 36_247
    assert theorem["pairwise_nonisomorphic_within_ABBB_shard"] is True
    assert theorem["type_facet_histogram"] == {"779": 36_247}
    assert theorem["raw_unrooted_coloring_facet_histogram"] == {"779": 13_862_400}
    assert theorem["reflection_index_and_volume_multiple"] == 5
    assert theorem["finite_volume"] is True
    assert theorem["compact"] is False
    assert Counter(row["edge_class_word"] for row in rows) == {"ABBB": 36_247}
    type_hashes = set()
    for row in rows:
        assert row["tree_shape"] == "P5"
        assert row["proper_edge_coloring"] is True
        assert row["edge_class_word"] == "ABBB"
        assert row["B_edge_occurrence_count"] == 3
        assert row["facet_count"] == row["expected_facet_count_752_plus_9B"] == 779
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


def test_independent_ABBB_Burnside_value() -> None:
    actions = star._restricted_actions()
    burnside = frozen._burnside_pattern_orbits("ABBB", actions)
    assert burnside["orbit_count"] == 36_247
    assert burnside["reversal_fixed_sum"] == 0
    assert burnside["Burnside_denominator"] == 384


def test_selected_exact_boundary_replays() -> None:
    rows = _archive()["orbit_representatives"]
    selected = [rows[0]]
    repeated = next(
        (row for row in rows if row["nonadjacent_color_repetition_count"] > 0), None
    )
    if repeated is not None and repeated["orbit_number"] != selected[0]["orbit_number"]:
        selected.append(repeated)
    _carrier, gram, roots, interior, _bounded = seed._load_inputs()
    engine = frozen.ExactP5BoundaryEngine(gram, roots, interior)
    for row in selected:
        stock_colors = tuple(
            color - frozen.STOCK_START - 1 for color in row["edge_colors_one_based"]
        )
        boundary, boundary_data = engine.boundary(stock_colors)
        _labels, degree_key, ordinary_data = engine.ordinary_data(boundary)
        assert boundary_data["facet_count"] == 779
        assert boundary_data["marked_boundary_root_sha256"] == row[
            "marked_boundary_root_sha256"
        ]
        assert ordinary_data["ordinary_degree_profile_multiset_sha256"] == row[
            "ordinary_degree_profile_multiset_sha256"
        ]
        assert degree_key


def test_combined_exact_progress_and_claim_boundary() -> None:
    payload = _archive()
    combined = payload["combined_exact_P5_progress"]
    assert combined["frozen_B_le_2_type_count"] == 20_693
    assert combined["new_ABBB_type_count"] == 36_247
    assert combined["combined_completed_type_count"] == 56_940
    assert combined["combined_completed_type_facet_histogram"] == {
        "752": 142,
        "761": 2_091,
        "770": 18_460,
        "779": 36_247,
    }
    assert combined["combined_positive_increment_beyond_all_five_original_P5_pilots"] == 56_936
    recovery = payload["pilot_recovery"]
    assert recovery["eligible_ABBB_B3_pilot_count"] == 1
    assert recovery["ABBB_positive_increment_beyond_original_B3_pilot"] == 36_246
    scope = payload["lineage_and_claim_boundary"]
    assert scope["complete_P5_census_claimed"] is False
    assert scope["remaining_BABB_BBBB_orbit_denominator"] == 144_807
    assert scope["cross_star_pure_exact_nonisomorphism_claimed"] is False
    assert scope["worldwide_priority"] == "open"
    carrier = payload["Allcock_and_carrier"]
    assert carrier["weak_divisor_reflection_carrier_only"] is True
    assert carrier["full_nef_effectivity_or_surface_action_claim"] is False
