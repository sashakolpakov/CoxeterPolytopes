from __future__ import annotations

import json

from scripts.certify_h9_14_record14_a2_schedule_no_go import (
    EXPECTED_CERTIFICATE_DIGEST,
    EXPECTED_TARGET_GROUPS,
    OUTPUT,
    TARGET_SOURCE_ORBITS,
    verify,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    audit_compact_labels,
    build_problem,
    enumerate_models,
    mask,
)
from scripts.independent_h9_14_record14_a2_audit import replay_orbit


def test_record14_a2_frozen_schedule_theorem_and_scope() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 14,
        "schedule_a": 2,
        "source_orbit_count": 36,
        "literal_target_count": 24,
        "literal_target_orbit_count": 6,
        "terminal_compact_label_state_count": 0,
        "closed_schedules": [2, 3, 4, 5, 6],
        "remaining_schedules": [0, 1],
        "whole_shadow_compact_coxeter_no_go": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_digest": EXPECTED_CERTIFICATE_DIGEST,
    }
    evidence = document["evidence"]
    assert evidence["source_target_count_histogram"] == {"0": 32, "6": 4}
    assert evidence["target_bearing_source_orbits"] == list(TARGET_SOURCE_ORBITS)
    assert evidence["all_target_orbits_unrestricted_compact_label_excluded"] is True


def test_record14_a2_literal_equality_and_target_orbits_are_complete() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rows = document["evidence"]["rows"]
    assert [row["source_positive_A4_orbit_index"] for row in rows] == list(range(36))
    assert all(row["production_independent_literal_equality"] for row in rows)
    assert all(
        row["production_z3_model_count"]
        == row["independent_cp_sat_model_count"]
        == row["literal_pure_target_count"]
        for row in rows
    )
    target_rows = [row for row in rows if row["literal_pure_target_count"]]
    assert {
        row["source_positive_A4_orbit_index"]: tuple(
            tuple(group) for group in row["literal_target_orbit_groups"]
        )
        for row in target_rows
    } == EXPECTED_TARGET_GROUPS
    assert sum(row["literal_target_orbit_count"] for row in target_rows) == 6


def test_record14_a2_two_solver_and_label_spot_replay() -> None:
    production_problem = build_problem(2, 5)
    production = enumerate_models(production_problem)
    independent = replay_orbit(5)
    assert production["literal_pure_targets"] == independent["literal_pure_targets"]
    assert production["literal_pure_target_count"] == 6
    representative = tuple(mask(face) for face in production["literal_pure_targets"][0])
    label = audit_compact_labels(production_problem["fours"], representative)
    assert [stage["passing_state_count"] for stage in label["stages"]] == [108, 84, 0]
    assert label["unrestricted_compact_label_no_go"] is True
