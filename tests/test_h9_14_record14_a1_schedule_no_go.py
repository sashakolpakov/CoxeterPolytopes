from __future__ import annotations

import json

from scripts.certify_h9_14_record14_a1_schedule_no_go import (
    EXPECTED_CERTIFICATE_DIGEST,
    EXPECTED_TARGET_GROUPS,
    OUTPUT,
    TARGET_SOURCE_ORBITS,
    _fixed_source_actions,
    verify,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import build_problem


def test_record14_a1_frozen_schedule_theorem_and_scope() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 14,
        "schedule_a": 1,
        "source_orbit_count": 5,
        "generator_lex_target_count": 19,
        "full_target_orbit_count": 10,
        "terminal_compact_label_state_count": 0,
        "closed_schedules": [1, 2, 3, 4, 5, 6],
        "remaining_schedules": [0],
        "whole_shadow_compact_coxeter_no_go": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_digest": EXPECTED_CERTIFICATE_DIGEST,
    }
    evidence = document["evidence"]
    assert evidence["source_target_count_histogram"] == {"0": 3, "8": 1, "11": 1}
    assert evidence["target_bearing_source_orbits"] == list(TARGET_SOURCE_ORBITS)
    assert evidence["all_target_orbits_unrestricted_compact_label_excluded"] is True


def test_record14_a1_two_solver_literal_equality_and_full_orbits() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rows = document["evidence"]["rows"]
    assert [row["source_positive_A4_orbit_index"] for row in rows] == list(range(5))
    assert all(row["production_independent_literal_equality"] for row in rows)
    assert all(
        row["production_z3_generator_lex_model_count"]
        == row["independent_cp_sat_generator_lex_model_count"]
        == row["generator_lex_slice_pure_target_count"]
        for row in rows
    )
    target_rows = [row for row in rows if row["generator_lex_slice_pure_target_count"]]
    observed = {
        row["source_positive_A4_orbit_index"]: tuple(
            sorted(tuple(group) for group in row["generator_lex_target_orbit_groups"])
        )
        for row in target_rows
    }
    assert observed == {
        source: tuple(sorted(groups)) for source, groups in EXPECTED_TARGET_GROUPS.items()
    }
    assert sum(row["target_orbit_count"] for row in target_rows) == 10


def test_record14_a1_full_stabilizer_orders_and_terminal_labels() -> None:
    problem0 = build_problem(1, 0)
    problem3 = build_problem(1, 3)
    assert len(_fixed_source_actions(problem0["fours"][0])) == 4320
    assert len(_fixed_source_actions(problem3["fours"][0])) == 10080
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    audits = document["evidence"]["unrestricted_compact_label_audits"]
    assert len(audits) == 10
    assert all(row["rank_nine_actual_vertex_count"] == 108 for row in audits)
    assert all(row["terminal_sparse_label_state_count"] == 0 for row in audits)
    assert all(row["unrestricted_compact_label_no_go"] for row in audits)
