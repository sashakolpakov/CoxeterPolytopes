from __future__ import annotations

from scripts.verify_h9_14_record31_a1_orbit3_transcript import (
    EXPECTED_ARCHIVE_SHA256,
    EXPECTED_FAMILIES,
    EXPECTED_GENERATED_SOURCE_SHA256,
    verify,
)


def test_terminal_source_orbit3_transcript_and_literal_replay() -> None:
    result = verify()
    assert result["archive_sha256"] == EXPECTED_ARCHIVE_SHA256
    assert result["generated_source_sha256"] == EXPECTED_GENERATED_SOURCE_SHA256
    assert result["search_node_count"] == 469_990_577
    assert result["pure_120_facet_target_count"] == 4_445
    assert result["ridge_valid_literal_target_count"] == 11
    assert result["ridge_valid_target_orbit_count"] == len(EXPECTED_FAMILIES) == 11
    assert result["target_literal_multiplicities"] == [1] * 11
    assert result["independent_blocker_and_ridge_replay"] is True
