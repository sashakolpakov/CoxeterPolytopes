import json
from pathlib import Path

from scripts.audit_h9_14_record0_a1_source1_checkpoint52_increment import (
    ARCHIVE,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(Path(ARCHIVE).read_text(encoding="utf-8"))


def test_fresh_incremental_audit_matches_archive() -> None:
    assert build_certificate() == _archive()


def test_both_new_core_orbits_are_independently_replayed() -> None:
    rows = _archive()["new_cores"]
    assert [row["human_core_number"] for row in rows] == [51, 52]
    assert [row["full_source_orbit_size"] for row in rows] == [10080, 10080]
    assert [row["core_stabilizer_order"] for row in rows] == [24, 24]
    assert rows[0]["independent_local_passing_state_counts"] == [
        108,
        110,
        48,
        0,
    ]
    assert rows[1]["independent_local_passing_state_counts"] == [
        108,
        110,
        62,
        103,
        10,
        0,
    ]
    assert all(all(row["checks"].values()) for row in rows)


def test_checkpoint_is_still_bounded_not_terminal() -> None:
    certificate = _archive()
    assert all(certificate["checks"].values())
    status = certificate["bounded_checkpoint_status"]
    assert status["excluded_local_core_orbits"] == 52
    assert status["terminal_solver_status"] == "OPTIMAL"
    assert status["stopped_reason"] == "core_cutoff"
    assert status["fixed_source_compact_coxeter_no_go"] is False
    assert status["theorem_certificate"] is False
    assert status["finite_volume_noncompact_checked"] is False


def test_cumulative_orbit_size_sum_is_scoped_exactly() -> None:
    status = _archive()["bounded_checkpoint_status"]
    assert status["previous_cumulative_orbit_size_sum"] == 367836
    assert status["incremental_orbit_size_sum"] == 20160
    assert status["cumulative_orbit_size_sum"] == 387996
    assert status["orbit_size_sum_is_not_a_cross_core_distinct_clause_count"]
