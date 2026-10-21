from __future__ import annotations

import json

from scripts.certify_h9_14_record31_positive_a5 import (
    ARCHIVE,
    fresh_record,
    verify,
)


def test_record31_a5_compact_schedule_is_exactly_closed() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 31,
        "schedule_count": 8,
        "strict_schedule_split": {
            "surviving": [0, 1, 2, 3, 4, 5],
            "excluded": [6, 7],
        },
        "newly_closed_schedule": 5,
        "remaining_open_schedules": [0, 1, 2, 3, 4],
        "production_a5_surviving_orbits": 0,
        "independent_a5_surviving_paths": 0,
        "compact_only": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }

    production = document["production"]
    assert production["strict_compact_excluded_a"] == [6, 7]
    assert production["local_link_constraints"]["upper_endpoint_histogram"] == {
        "0": 30,
        "1": 156,
        "2": 84,
        "3": 103,
        "4": 2,
        "5": 36,
        "6": 9,
    }
    a5 = production["branches"][-1]
    assert a5["source_a4_full_orbits"] == 62_234
    assert a5["interval_extension_paths"] == 34_615_149
    assert a5["shadow_feasible_extension_paths"] == 158_736
    assert a5["orbit_count"] == 0
    assert document["independent"]["viable_paths"] == 0


def test_record31_a5_fresh_production_and_independent_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    rebuilt = fresh_record()
    assert rebuilt == archived
    assert verify(rebuilt)["newly_closed_schedule"] == 5
