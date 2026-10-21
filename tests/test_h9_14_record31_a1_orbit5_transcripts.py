from __future__ import annotations

from scripts.verify_h9_14_record31_a1_orbit5_transcripts import (
    EXPECTED_INDEPENDENT_SHA256,
    EXPECTED_PRODUCTION_SHA256,
    verify,
)


def test_terminal_source_orbit5_two_quotient_no_go() -> None:
    result = verify()
    assert result["production_archive_sha256"] == EXPECTED_PRODUCTION_SHA256
    assert result["independent_archive_sha256"] == EXPECTED_INDEPENDENT_SHA256
    assert result["production_search_node_count"] == 41_488_368
    assert result["independent_search_node_count"] == 6_042_612
    assert result["pure_120_facet_target_count"] == 735
    assert result["pure_ridge_signature_count"] == 231
    assert result["ridge_valid_target_orbit_count"] == 0
    assert result["terminal_pure_target_censuses_identical"] is True
