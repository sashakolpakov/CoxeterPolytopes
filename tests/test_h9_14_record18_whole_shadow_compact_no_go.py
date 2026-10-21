from __future__ import annotations

import json

from scripts.certify_h9_14_record18_whole_shadow_compact_no_go import (
    OUTPUT,
    fresh_record,
    verify_record,
)


def _load():
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_record18_whole_shadow_compact_certificate() -> None:
    document = _load()
    summary = verify_record(document)
    assert summary["compact_record18_excluded"] is True
    assert summary["closed_schedules"] == list(range(7))
    assert summary["finite_volume_noncompact_checked"] is False
    assert summary["new_hyperbolic_example"] is False

    branches = document["branch_summary"]
    assert branches["a_ge_7_global_g4_excluded"] is True
    assert branches["a5_a6_summed_link_g4_excluded"] is True
    assert branches["a4_positive_orbit_count"] == 13
    assert branches["a4_incidence_target_count"] == 0
    assert branches["a3_positive_orbit_count"] == 40
    assert branches["a3_incidence_target_count"] == 0
    assert branches["a2_positive_orbit_count"] == 14
    assert branches["a2_add4_literal_target_count"] == 8
    assert branches["a2_label_survivor_count"] == 0
    assert branches["a1_positive_orbit_count"] == 3
    assert branches["a1_literal_target_count"] == 1080
    assert branches["a1_label_survivor_count"] == 0
    assert branches["a0_literal_target_count"] == 196
    assert branches["a0_label_survivor_count"] == 0
    assert branches["all_nonnegative_a_schedules_closed"] is True


def test_record18_aggregate_fresh_composition_matches_archive() -> None:
    assert fresh_record() == _load()
