from __future__ import annotations

import json

from scripts.certify_h9_14_record14_upper_schedules import (
    EXPECTED_CERTIFICATE_DIGEST,
    EXPECTED_SCHEDULE_ROWS,
    OUTPUT,
    verify,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    build_problem,
    enumerate_models,
)


def test_record14_upper_schedule_certificate_is_literal_and_scoped() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 14,
        "closed_schedules": [3, 4, 5, 6],
        "enumerated_source_orbit_count": 72,
        "production_target_count": 0,
        "independent_target_count": 0,
        "remaining_schedules": [0, 1, 2],
        "whole_shadow_compact_coxeter_no_go": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_digest": EXPECTED_CERTIFICATE_DIGEST,
    }
    evidence = document["evidence"]
    assert evidence["schedule_rows"] == [list(row) for row in EXPECTED_SCHEDULE_ROWS]
    assert evidence["source_orbit_count_by_schedule"] == {"3": 59, "4": 13}
    assert evidence["production_search_node_total"] == 6_566_617
    assert evidence["production_no_small_transversal_family_total"] == 30_399
    assert evidence["independent_z3_model_total"] == 0


def test_record14_upper_source_indices_cover_every_orbit_once() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rows = document["evidence"]["rows"]
    assert [(row["a"], row["source_positive_A4_orbit_index"]) for row in rows] == [
        *((3, index) for index in range(59)),
        *((4, index) for index in range(13)),
    ]
    assert all(row["production_literal_pure_target_count"] == 0 for row in rows)
    assert all(row["independent_literal_pure_target_count"] == 0 for row in rows)
    assert all(row["independent_enumeration_complete"] is True for row in rows)


def test_record14_upper_independent_replay_hits_hard_branches() -> None:
    for a, orbit_index in ((3, 6), (4, 0)):
        replay = enumerate_models(build_problem(a, orbit_index))
        assert replay["enumeration_complete"] is True
        assert replay["z3_model_count"] == 0
        assert replay["literal_pure_target_count"] == 0
