from __future__ import annotations

import gzip
import hashlib
import json

import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_complete_streaming_census as complete
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard as census
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star


SUMMARY, ROWS = census._default_paths("BBBB")


def _summary() -> dict:
    return json.loads(SUMMARY.read_text(encoding="utf-8"))


def test_BBBB_compact_archive_and_denominator() -> None:
    payload = _summary()
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert hashlib.sha256(census.Path(census.__file__).read_bytes()).hexdigest() == payload[
        "pinned_inputs"
    ]["builder"]
    denominator = payload["coloring_denominator"]
    assert denominator["raw_unrooted_proper_colorings"] == 41_154_000
    assert denominator["explicit_orbit_size_sum"] == 41_154_000
    assert denominator["orbit_count"] == denominator["explicit_orbit_cover_count"] == 108_159
    assert denominator["Burnside_denominator"] == 768
    assert denominator["reversal_fixed_sum"] == 659_328


def test_BBBB_rows_and_stream_digests() -> None:
    payload = _summary()
    digest = hashlib.sha256()
    count = orbit_sum = refined = 0
    with gzip.open(ROWS, "rb") as stream:
        for line in stream:
            digest.update(line)
            row = json.loads(line)
            count += 1
            orbit_sum += row["unrooted_coloring_orbit_size_under_seed_action"]
            refined += int(row["color_refinement_applied"])
            assert row["orbit_number"] == count
            assert row["edge_class_word"] == "BBBB"
            assert row["facet_count"] == 788
    transcript = payload["compact_exact_transcript"]
    assert count == transcript["row_count"] == 108_159
    assert orbit_sum == 41_154_000
    assert refined == 13_591
    assert digest.hexdigest() == transcript["uncompressed_canonical_JSONL_sha256"]
    assert hashlib.sha256(ROWS.read_bytes()).hexdigest() == transcript[
        "deterministic_gzip_sha256"
    ]


def test_BBBB_nonisomorphism_and_sole_star_pair() -> None:
    payload = _summary()
    theorem = payload["theorem"]
    assert theorem["ordinary_Coxeter_type_count"] == 108_159
    assert theorem["pairwise_nonisomorphic_within_shard"] is True
    assert theorem["type_facet_histogram"] == {"788": 108_159}
    exact = payload["ordinary_Coxeter_nonisomorphism_certificate"]
    assert exact["degree_digest_class_count"] == 100_935
    assert exact["degree_digest_collision_group_count"] == 6_366
    assert exact["degree_digest_collision_member_count"] == 13_590
    assert exact["literal_degree_collision_group_count"] == 6_366
    assert exact["literal_degree_collision_member_count"] == 13_590
    assert exact["color_refined_current_member_count"] == 13_591
    assert exact["color_refined_discrete_current_count"] == 13_232
    cross = payload["cross_prior_exact_separation"]
    assert cross["prior_same_facet_row_count"] == 8_985
    assert cross["archived_degree_digest_candidate_count"] == 1
    assert cross["prior_literal_degree_replay_count"] == 1
    assert cross["literal_cross_prior_degree_match_count"] == 1
    assert cross["cross_prior_exact_refinement_pair_count"] == 1
    assert cross["unresolved_cross_prior_invariant_pair_count"] == 0


def test_BBBB_independent_Burnside_and_selected_exact_replay() -> None:
    actions = star._restricted_actions()
    assert p5._burnside_pattern_orbits("BBBB", actions)["orbit_count"] == 108_159
    selected = []
    with gzip.open(ROWS, "rt", encoding="utf-8") as stream:
        for line in stream:
            row = json.loads(line)
            if not selected or row["color_refinement_applied"]:
                selected.append(row)
            if len(selected) == 2:
                break
    _carrier, gram, roots, interior, _bounded = p5.base._load_inputs()
    engine = p5.ExactP5BoundaryEngine(gram, roots, interior)
    for row in selected:
        colors = tuple(color - p5.STOCK_START - 1 for color in row["edge_colors_one_based"])
        boundary, boundary_data = engine.boundary(colors)
        _labels, degree_key, ordinary_data = engine.ordinary_data(boundary)
        assert boundary_data["facet_count"] == 788
        assert boundary_data["marked_boundary_root_sha256"] == row[
            "marked_boundary_root_sha256"
        ]
        assert hashlib.sha256(degree_key).hexdigest() == row[
            "ordinary_degree_profile_multiset_sha256"
        ]
        assert census._json_sha256({**boundary_data, **ordinary_data}) == row[
            "exact_gate_record_sha256"
        ]


def test_complete_P5_merge_and_claim_boundaries() -> None:
    payload = json.loads(complete.OUTPUT.read_text(encoding="utf-8"))
    assert payload == complete.build()
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    theorem = payload["theorem"]
    assert theorem["ordinary_Coxeter_type_count"] == 201_747
    assert theorem["pairwise_nonisomorphic"] is True
    assert theorem["type_facet_histogram"] == {
        "752": 142,
        "761": 2_091,
        "770": 18_460,
        "779": 72_895,
        "788": 108_159,
    }
    accounting = payload["project_accounting"]
    assert accounting["new_BABB_plus_BBBB_census_type_count"] == 144_807
    assert accounting["complete_P5_increment_beyond_original_pilots"] == 201_742
    assert accounting["new_BABB_plus_BBBB_increment_beyond_original_pilots"] == 144_806
    union = payload["exact_high_B_cross_shape_union"]
    assert union["pairwise_nonisomorphic_union_type_count"] == 197_926
    assert union["sole_star_BBBB_equal_degree_profile_pair_count"] == 1
    assert union["sole_pair_literal_degree_bytes_equal"] is True
    assert union["sole_pair_exact_refinement_invariants_equal"] is False
    assert union["star_plus_all_P5_claimed"] is False
    assert payload["claim_boundary"]["worldwide_priority"] == "open"
    assert payload["claim_boundary"]["all_five_vertex_shapes_complete_census_claimed"] is False
