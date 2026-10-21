from __future__ import annotations

import json

from scripts.certify_independent_h9_14_record18_positive_audit import (
    OUTPUT,
    build,
    verify,
)


def test_record18_independent_positive_archive_has_exact_scoped_result() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "positive_branch_counts_reproduced": True,
        "complete_orbit_key_sets_reproduced": True,
        "positive_orbit_count": 70,
        "a4_orbit_count": 13,
        "a4_schedule_incidence_empty": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def test_record18_independent_positive_sources_freshly_replay() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rebuilt = build()
    assert rebuilt == archived
    assert verify(rebuilt)["a4_schedule_incidence_empty"] is True
