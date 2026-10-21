from __future__ import annotations

import json

from scripts.certify_h9_14_record30_a4_no_go import (
    ARCHIVE,
    fresh_record,
    verify,
)


def test_record30_entire_compact_a4_schedule_is_excluded() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 30,
        "newly_closed_schedule": 4,
        "a4_orbit_count": 7,
        "no_small_transversal_family_count": 840,
        "scheduled_124_facet_target_count": 0,
        "remaining_open_schedules": [0, 1, 2],
        "compact_only": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "projection_sha256": (
            "e11696f63b8cf3105292295e14370a598594bf96dd702eb1a415c31cfc2630ac"
        ),
    }
    assert document["production"]["orbits"][0][
        "order_five_transversal_count_histogram"
    ] == {"92": 90, "97": 120, "102": 480, "107": 150}


def test_record30_a4_fresh_literal_and_pivoted_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    rebuilt = fresh_record()
    assert rebuilt == archived
    assert verify(rebuilt)["remaining_open_schedules"] == [0, 1, 2]
