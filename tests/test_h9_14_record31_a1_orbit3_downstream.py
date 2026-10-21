from __future__ import annotations

from scripts.verify_h9_14_record31_a1_orbit3_downstream import (
    EXPECTED_INDEPENDENT_SHA256,
    EXPECTED_PRODUCTION_SHA256,
    verify,
)


def test_source_orbit3_topology_and_independent_compact_label_no_go() -> None:
    result = verify(replay=True)
    assert result["production_downstream_sha256"] == EXPECTED_PRODUCTION_SHA256
    assert result["independent_downstream_sha256"] == EXPECTED_INDEPENDENT_SHA256
    assert result["target_orbit_count"] == 11
    assert result["topology_survivor_count"] == 11
    assert result["compact_label_survivor_count"] == 0
    assert result["gram_candidate_count"] == 0
    assert result["production_and_independent_replayed"] is True
