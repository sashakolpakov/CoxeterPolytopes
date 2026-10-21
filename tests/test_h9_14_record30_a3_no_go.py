from __future__ import annotations

import json

from scripts.certify_h9_14_record30_a3_no_go import (
    ARCHIVE,
    fresh_record,
    verify,
)


def test_record30_entire_compact_a3_schedule_is_excluded() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    summary = verify(document)
    assert summary == {
        "record_index": 30,
        "newly_closed_schedule": 3,
        "a3_orbit_count": 32,
        "no_small_transversal_family_count": 1_896_496,
        "scheduled_126_facet_target_count": 0,
        "remaining_open_schedules": [0, 1, 2, 4],
        "compact_only": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "projection_sha256": (
            "7afaebc75008dd16c6bf499df24de767c5f0a9cd9bc672b563e175c7de38c5a5"
        ),
    }
    for route in (document["production"], document["independent"]):
        assert all(
            "126" not in row["order_five_transversal_count_histogram"]
            for row in route["orbits"]
        )


def test_record30_a3_fresh_production_and_independent_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    rebuilt = fresh_record()
    assert rebuilt == archived
    assert verify(rebuilt)["remaining_open_schedules"] == [0, 1, 2, 4]
