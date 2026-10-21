from __future__ import annotations

import json

from scripts.certify_h9_14_record30_positive_scout import ARCHIVE, verify


def test_record30_positive_scout_archive_counts_and_scope() -> None:
    assert verify(json.loads(ARCHIVE.read_text()), allow_unfrozen=True) == {
        "positive_orbit_counts": [3, 21, 32, 7, 1],
        "a5_compact_schedule_excluded": True,
        "strict_sum_excluded_schedules": [6, 7],
        "whole_shadow_compact_coxeter_excluded": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }
