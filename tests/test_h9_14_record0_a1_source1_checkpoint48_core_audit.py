import json
from pathlib import Path

from scripts.audit_h9_14_record0_a1_source1_checkpoint48_core import (
    ARCHIVE,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(Path(ARCHIVE).read_text(encoding="utf-8"))


def test_fresh_certificate_matches_frozen_archive() -> None:
    assert build_certificate() == _archive()


def test_new_core_is_an_independent_lossless_label_contradiction() -> None:
    certificate = _archive()
    assert all(certificate["checks"].values())
    assert certificate["new_core"]["independent_passing_state_counts"] == [
        108,
        110,
        16,
        0,
    ]
    assert certificate["new_core"]["full_source_orbit_size"] == 168
    assert certificate["new_core"]["core_stabilizer_order"] == 1440


def test_scope_remains_bounded_and_nonterminal() -> None:
    status = _archive()["bounded_checkpoint_status"]
    assert status["terminal_solver_status"] == "OPTIMAL"
    assert status["stopped_reason"] == "core_cutoff"
    assert status["excluded_local_core_orbits"] == 48
    assert status["fixed_source_compact_coxeter_no_go"] is False
    assert status["theorem_certificate"] is False
