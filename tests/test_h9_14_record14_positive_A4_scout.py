from __future__ import annotations

import json

from scripts.certify_h9_14_record14_positive_A4_scout import (
    ARCHIVE,
    fresh_record,
    verify,
)


def test_record14_positive_scout_frozen_graph_derivation_and_scope() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 14,
        "schedule_count": 7,
        "local_constraint_count": 410,
        "strict_schedule_split": {
            "surviving": [0, 1, 2, 3, 4],
            "excluded": [5, 6],
        },
        "positive_branch_orbit_count": 113,
        "positive_branch_closed_count": 0,
        "collision_free_a4_ground_recovery": True,
        "compact_only": True,
        "finite_volume_noncompact_chamber_excluded": False,
    }

    vertices = {
        row["vertex"]: row["missing_four_incidence_interval"]
        for row in document["local_link_constraints"]["vertices"]
    }
    assert vertices == {
        0: [0, 3],
        1: [0, 3],
        2: [0, 4],
        **{vertex: [0, 5] for vertex in range(3, 12)},
        12: [0, 1],
        13: [0, 0],
    }
    assert document["local_link_constraints"]["upper_endpoint_histogram"] == {
        "0": 66,
        "1": 146,
        "2": 57,
        "3": 95,
        "4": 37,
        "5": 9,
    }


def test_record14_positive_scout_fresh_source_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    rebuilt = fresh_record()
    assert rebuilt == archived
    assert verify(rebuilt)["collision_free_a4_ground_recovery"] is True
