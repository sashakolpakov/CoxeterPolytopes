import json
from pathlib import Path

from scripts.audit_h9_14_record0_a1_source1_checkpoint56_increment import (
    ARCHIVE,
    EXPECTED_DISTINCT_CLAUSE_COUNT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(Path(ARCHIVE).read_text(encoding="utf-8"))


def test_fresh_incremental_and_global_audit_matches_archive() -> None:
    assert build_certificate() == _archive()


def test_all_four_new_core_orbits_are_independently_replayed() -> None:
    rows = _archive()["new_cores"]
    assert [row["human_core_number"] for row in rows] == [53, 54, 55, 56]
    assert [row["full_source_orbit_size"] for row in rows] == [
        504,
        504,
        2016,
        10080,
    ]
    assert [row["core_stabilizer_order"] for row in rows] == [
        480,
        480,
        120,
        24,
    ]
    assert all(all(row["checks"].values()) for row in rows)


def test_all_56_orbits_are_distinct_with_exact_global_clause_count() -> None:
    census = _archive()["global_orbit_census"]
    assert census["orbit_family_count"] == 56
    assert (
        census["independently_rebuilt_orbit_size_sum"]
        == census["distinct_clause_signature_count"]
        == EXPECTED_DISTINCT_CLAUSE_COUNT
    )
    assert all(census["checks"].values())


def test_checkpoint_is_still_bounded_not_terminal() -> None:
    certificate = _archive()
    assert all(certificate["checks"].values())
    status = certificate["bounded_checkpoint_status"]
    assert status["excluded_local_core_orbits"] == 56
    assert status["terminal_solver_status"] == "OPTIMAL"
    assert status["stopped_reason"] == "core_cutoff"
    assert status["fixed_source_compact_coxeter_no_go"] is False
    assert status["theorem_certificate"] is False
    assert status["finite_volume_noncompact_checked"] is False


def test_increment_is_13104_distinct_clause_signatures() -> None:
    status = _archive()["bounded_checkpoint_status"]
    assert status["previous_distinct_clause_count"] == 387996
    assert status["incremental_distinct_clause_count"] == 13104
    assert status["distinct_clause_signature_count"] == 401100
    assert status["cross_core_distinctness_independently_checked"]
