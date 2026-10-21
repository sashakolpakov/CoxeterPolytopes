from __future__ import annotations

from scripts.verify_h9_14_record31_a1_orbit0_independent_transcript import verify


def test_source_orbit0_independent_incidence_replay() -> None:
    result = verify()
    assert result["production_search_node_count"] == 231_773
    assert result["independent_search_node_count"] == 242_124
    assert result["ridge_valid_target_orbit_count"] == 1
    assert result["terminal_pure_target_censuses_identical"] is True
    assert result["independent_literal_blocker_and_ridge_replay"] is True
