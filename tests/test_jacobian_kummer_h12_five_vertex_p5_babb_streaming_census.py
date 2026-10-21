from __future__ import annotations

import gzip
import hashlib
import json

import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_streaming_shard as census
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star


SUMMARY, ROWS = census._default_paths("BABB")


def _summary() -> dict:
    return json.loads(SUMMARY.read_text(encoding="utf-8"))


def _rows() -> list[dict]:
    with gzip.open(ROWS, "rt", encoding="utf-8") as stream:
        return [json.loads(line) for line in stream]


def test_compact_archive_hashes_and_denominator() -> None:
    payload = _summary()
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["pattern"] == "BABB"
    assert hashlib.sha256(census.Path(census.__file__).read_bytes()).hexdigest() == payload[
        "pinned_inputs"
    ]["builder"]
    denominator = payload["coloring_denominator"]
    assert denominator["raw_unrooted_proper_colorings"] == 14_008_320
    assert denominator["orbit_count"] == 36_648
    assert denominator["explicit_orbit_cover_count"] == 36_648
    assert denominator["explicit_orbit_size_sum"] == 14_008_320
    assert denominator["Burnside_denominator"] == 384


def test_compact_rows_and_stream_digests() -> None:
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
            assert row["edge_class_word"] == "BABB"
            assert row["facet_count"] == 779
            assert len(row["ordinary_degree_profile_multiset_sha256"]) == 64
            assert len(row["exact_gate_record_sha256"]) == 64
    transcript = payload["compact_exact_transcript"]
    assert count == transcript["row_count"] == 36_648
    assert orbit_sum == 14_008_320
    assert refined == 4_015
    assert digest.hexdigest() == transcript["uncompressed_canonical_JSONL_sha256"]
    assert hashlib.sha256(ROWS.read_bytes()).hexdigest() == transcript[
        "deterministic_gzip_sha256"
    ]


def test_nonisomorphism_and_cross_prior_certificate() -> None:
    payload = _summary()
    theorem = payload["theorem"]
    assert theorem["ordinary_Coxeter_type_count"] == 36_648
    assert theorem["pairwise_nonisomorphic_within_shard"] is True
    assert theorem["type_facet_histogram"] == {"779": 36_648}
    assert theorem["separated_from_previously_completed_P5_rows"] is True
    assert theorem["separated_from_complete_same_facet_star_rows"] is True
    exact = payload["ordinary_Coxeter_nonisomorphism_certificate"]
    assert exact["degree_digest_class_count"] == 34_583
    assert exact["degree_digest_collision_group_count"] == 1_950
    assert exact["degree_digest_collision_member_count"] == 4_015
    assert exact["literal_degree_collision_group_count"] == 1_950
    assert exact["literal_degree_collision_member_count"] == 4_015
    assert exact["color_refined_discrete_current_count"] == 3_990
    cross = payload["cross_prior_exact_separation"]
    assert cross["prior_same_facet_row_count"] == 36_247 + 6_117
    assert cross["archived_degree_digest_candidate_count"] == 0
    assert cross["literal_cross_prior_degree_match_count"] == 0
    assert cross["unresolved_cross_prior_invariant_pair_count"] == 0


def test_independent_burnside_and_selected_exact_replay() -> None:
    actions = star._restricted_actions()
    assert p5._burnside_pattern_orbits("BABB", actions)["orbit_count"] == 36_648
    rows = _rows()
    selected = [rows[0], next(row for row in rows if row["color_refinement_applied"])]
    _carrier, gram, roots, interior, _bounded = p5.base._load_inputs()
    engine = p5.ExactP5BoundaryEngine(gram, roots, interior)
    for row in selected:
        colors = tuple(color - p5.STOCK_START - 1 for color in row["edge_colors_one_based"])
        boundary, boundary_data = engine.boundary(colors)
        _labels, degree_key, ordinary_data = engine.ordinary_data(boundary)
        gate_data = {**boundary_data, **ordinary_data}
        assert boundary_data["facet_count"] == 779
        assert boundary_data["marked_boundary_root_sha256"] == row[
            "marked_boundary_root_sha256"
        ]
        assert hashlib.sha256(degree_key).hexdigest() == row[
            "ordinary_degree_profile_multiset_sha256"
        ]
        assert census._json_sha256(gate_data) == row["exact_gate_record_sha256"]


def test_arithmetic_Allcock_resource_and_claim_boundaries() -> None:
    payload = _summary()
    gates = payload["aggregate_arithmetic_gates"]
    assert gates["minimum_strict_interior_facet_pairing"] == 5
    assert gates["ordinary_relation_histogram_inventory"] == [
        {
            "histogram": {
                "infinity": 285_024,
                "m2": 13_442,
                "m3": 656,
                "m4": 3_909,
                "m6": 0,
            },
            "type_count": 36_648,
        }
    ]
    allcock = payload["Allcock_and_carrier"]
    assert allcock["finite_volume"] is True
    assert allcock["compact"] is False
    assert allcock["reflection_index_and_volume_multiplier"] == 5
    assert allcock["weak_divisor_reflection_carrier_only"] is True
    assert allcock["full_nef_effectivity_or_surface_action_claim"] is False
    assert payload["resource_certificate"]["RSS_ceiling_bytes"] == 1 << 30
    scope = payload["lineage_and_claim_boundary"]
    assert scope["direction"] == "SURFACE_TO_POLYTOPE"
    assert scope["worldwide_priority"] == "open"
    assert scope["worldwide_novelty_claim"] is False
