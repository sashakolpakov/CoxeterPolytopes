from __future__ import annotations

import json

from scripts.certify_h9_14_record31_a4_no_go import (
    ARCHIVE,
    fresh_record,
    verify,
)


def test_record31_entire_compact_a4_schedule_is_excluded() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 31,
        "newly_closed_schedule": 4,
        "a4_orbit_count": 20,
        "scheduled_target_count": 0,
        "remaining_open_schedules": [0, 1, 2, 3],
        "compact_only": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }
    assert document["production_low_additions"][
        "literal_low_addition_target_count"
    ] == 0
    high = {
        row["orbit_representative_index"]: row
        for row in document["production_three_additions"][
            "representative_audits"
        ]
    }
    assert high[2]["no_small_transversal_addition_count"] == 0
    assert high[18]["facet_count_histogram"] == {"105": 150, "109": 120}
    assert high[19]["facet_count_histogram"]["115"] == 1
    assert all(
        row["scheduled_114_facet_target_count"] == 0
        for row in document["independent"]["representative_audits"]
    )


def test_record31_a4_fresh_production_and_independent_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    rebuilt = fresh_record()
    assert rebuilt == archived
    assert verify(rebuilt)["remaining_open_schedules"] == [0, 1, 2, 3]
