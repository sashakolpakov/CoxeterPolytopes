from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_h9_14_record0_a1_source0_joint_branches import (
    EXPECTED_CERTIFICATE_SHA256,
    EXPECTED_MODEL_COUNTS,
    EXPECTED_ORBIT_SIZES,
    OUTPUT,
    independent_model_counts,
    partition_record,
    verify_record,
)


def archive() -> dict[str, object]:
    return json.loads(Path(OUTPUT).read_text(encoding="utf-8"))


def test_archive_replays_and_is_frozen() -> None:
    summary = verify_record(archive())
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["closed_branch_count"] == 9
    assert summary["covered_classical_A4_placement_count"] == 108
    assert summary["fixed_source_compact_coxeter_no_go"] is True


def test_partition_is_exact_and_complete() -> None:
    partition = partition_record()
    assert partition["orbit_count"] == 9
    assert partition["orbit_sizes"] == list(EXPECTED_ORBIT_SIZES)
    assert sum(partition["orbit_sizes"]) == 108
    assert partition["coverage_is_disjoint_and_complete"] is True


def test_independent_joint_model_counts() -> None:
    assert independent_model_counts() == EXPECTED_MODEL_COUNTS


def test_every_branch_is_terminal_without_a_survivor() -> None:
    document = archive()
    evidence = document["branch_evidence"]
    assert [row["branch_index"] for row in evidence] == list(range(9))
    assert [row["placement_orbit_size"] for row in evidence] == list(
        EXPECTED_ORBIT_SIZES
    )
    assert all(row["terminal_solver_status"] == "INFEASIBLE" for row in evidence)
    assert all(row["compact_label_survivor"] is None for row in evidence)
    assert all(row["branch_compact_coxeter_no_go"] is True for row in evidence)


def test_schedule_corollary_and_broader_scope_are_literal() -> None:
    document = archive()
    assert document["record0_a1_corollary"]["closed_source_orbit_indices"] == [
        0,
        1,
        2,
        3,
    ]
    assert document["theorem"]["whole_record0_a1_compact_schedule_no_go"] is True
    assert document["theorem"]["whole_record0_compact_shadow_no_go"] is False
    assert document["theorem"]["finite_volume_noncompact_checked"] is False
    assert document["theorem"]["compact_label_or_gram_survivor"] is False
    assert document["broader_frontier"]["global_wholly_excluded_shadow_count"] == 134
    assert document["broader_frontier"]["global_open_shadow_records"] == [
        0,
        1,
        3,
        4,
        5,
        7,
        9,
        31,
    ]
