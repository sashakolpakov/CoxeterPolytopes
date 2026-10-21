from __future__ import annotations

import json

from scripts.certify_independent_h9_14_record18_a3_audit import (
    OUTPUT,
    build,
    verify,
)


def test_record18_a3_independent_archive_closes_only_a3() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "input_positive_orbit_count": 40,
        "complete_addition_transcript_replayed": True,
        "no_small_transversal_addition_count": 6_011,
        "scheduled_pure_104_facet_target_count": 0,
        "a3_schedule_incidence_empty": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def test_record18_a3_no_project_sources_freshly_replay_all_orbits() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rebuilt = build()
    assert rebuilt == archived
    assert verify(rebuilt)["a3_schedule_incidence_empty"] is True
