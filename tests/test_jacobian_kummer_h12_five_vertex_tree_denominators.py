from __future__ import annotations

import json

from scripts.certify_jacobian_kummer_h12_five_vertex_tree_denominators import (
    OUTPUT,
    build,
)


def test_exact_replay_matches_frozen_denominator_archive() -> None:
    payload = build()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())


def test_all_three_shape_denominators_and_claim_boundary() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    shapes = payload["tree_shapes"]
    assert shapes["P5"]["raw_unrooted_colorings"] == 76_587_336
    assert shapes["P5"]["Burnside_numerator"] == 154_941_696
    assert shapes["P5"]["seed_and_tree_automorphism_orbits"] == 201_747
    assert shapes["K1,4"]["raw_unrooted_colorings"] == 6_210_820
    assert shapes["K1,4"]["Burnside_numerator"] == 6_478_848
    assert shapes["K1,4"]["seed_and_tree_automorphism_orbits"] == 16_872
    assert shapes["degree_3_fork"]["raw_unrooted_colorings"] == 75_897_360
    assert shapes["degree_3_fork"]["Burnside_numerator"] == 153_110_016
    assert shapes["degree_3_fork"]["seed_and_tree_automorphism_orbits"] == 199_362
    assert payload["total_raw_unrooted_colored_trees"] == 158_695_516
    assert payload["total_seed_and_tree_automorphism_orbits"] == 417_981
    boundary = payload["claim_boundary"]
    assert boundary["this_is_a_colored_tree_denominator_not_a_polytope_census"] is True
    assert boundary[
        "different_color_orbits_may_still_give_isomorphic_ordinary_Coxeter_matrices"
    ] is True
