from __future__ import annotations

import hashlib
import json

import scripts.certify_jacobian_kummer_h12_five_vertex_fork_b3_terminal_acceptance as terminal
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star
import scripts.certify_jacobian_kummer_h12_star_p5_b_le_2_exact_separation as exact


def test_frozen_exact_low_B_union_and_provenance() -> None:
    payload = json.loads(exact.OUTPUT.read_text())
    assert payload["all_checks_pass"] and all(payload["checks"].values())
    assert payload["provenance"]["this_builder_sha256"] == hashlib.sha256(
        exact.Path(exact.__file__).read_bytes()
    ).hexdigest()
    assert payload["provenance"]["pinned_frozen_inputs_sha256"] == exact.PINNED_SHA256
    assert payload["denominator"]["star_type_count_by_B"] == {
        "0": 19,
        "1": 169,
        "2": 1_582,
    }
    assert payload["denominator"]["P5_type_count_by_B"] == {
        "0": 142,
        "1": 2_091,
        "2": 18_460,
    }
    assert payload["exact_replay"]["literal_cross_shape_degree_match_pair_count"] == 0
    assert payload["exact_replay"]["unresolved_pair_count"] == 0
    theorem = payload["theorem"]
    assert theorem["complete_star_plus_complete_P5_pairwise_nonisomorphic_type_count"] == 218_619
    assert theorem[
        "complete_P5_star_plus_fork_B_le_3_pairwise_nonisomorphic_type_count"
    ] == 311_596
    assert theorem["finite_volume"] and not theorem["compact"]


def test_one_literal_exact_replay_per_low_B_header() -> None:
    star_rows: dict[int, tuple[tuple[int, ...], str]] = {}
    for _index, colors, digest in terminal._pretty_json_orbit_rows(exact.STAR, 16_872):
        b = exact._b_count(colors)
        if b <= 2:
            star_rows.setdefault(b, (colors, digest))
        if len(star_rows) == 3:
            break
    p5_rows: dict[int, tuple[tuple[int, ...], str]] = {}
    for _index, colors, digest in terminal._pretty_json_orbit_rows(exact.P5_LOW, 20_693):
        b = exact._b_count(colors)
        p5_rows.setdefault(b, (colors, digest))
        if len(p5_rows) == 3:
            break
    assert set(star_rows) == set(p5_rows) == {0, 1, 2}
    _carrier, gram, roots, interior, _bounded = p5.base._load_inputs()
    star_engine = star.ExactStarBoundaryEngine(gram, roots, interior)
    p5_engine = p5.ExactP5BoundaryEngine(gram, roots, interior)
    for b in range(3):
        star_key, _labels = exact._rebuild(
            star_engine, *star_rows[b], exact.FACETS[b]
        )
        p5_key, _labels = exact._rebuild(p5_engine, *p5_rows[b], exact.FACETS[b])
        assert star_key != p5_key


def test_only_B4_remains_in_the_five_vertex_denominator() -> None:
    payload = json.loads(exact.OUTPUT.read_text())
    assert payload["claim_boundary"]["only_unclassified_fork_B4_denominator"] == 106_385
    assert not payload["claim_boundary"]["complete_five_vertex_three_shape_census_claimed"]
    assert payload["claim_boundary"]["worldwide_priority"] == "open"
    assert not payload["claim_boundary"]["worldwide_novelty_claim"]
