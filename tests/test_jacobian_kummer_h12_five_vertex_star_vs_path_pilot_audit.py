from __future__ import annotations

import json

from scripts.audit_jacobian_kummer_h12_five_vertex_star_vs_path_pilot import (
    OUTPUT,
    build,
)


def test_exact_cross_shape_replay_matches_archive() -> None:
    payload = build()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())


def test_all_five_path_pilots_are_separated_from_complete_star_census() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    result = payload["result"]
    assert result["complete_K14_type_count"] == 16_872
    assert result["explicit_P5_pilot_type_count"] == 5
    assert result["pairwise_nonisomorphic_combined_five_chamber_type_count"] == 16_877
    assert result["finite_volume"] is True
    assert result["compact"] is False
    assert result["reflection_index_and_volume_multiplier"] == 5
    rows = payload["P5_rows"]
    assert [row["B_edge_occurrence_count"] for row in rows] == list(range(5))
    assert [row["facet_count"] for row in rows] == [752, 761, 770, 779, 788]
    assert all(
        row["equal_degree_profile_multiset_found_in_complete_K14_census"] is False
        for row in rows
    )
    boundary = payload["claim_boundary"]
    assert boundary["complete_P5_census_claimed"] is False
    assert boundary["complete_degree_3_fork_census_claimed"] is False
    assert boundary["worldwide_novelty_claimed"] is False
