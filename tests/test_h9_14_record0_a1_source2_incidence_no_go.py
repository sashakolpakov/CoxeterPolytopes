from __future__ import annotations

import json

from scripts.certify_h9_14_record0_a1_source2_incidence_no_go import (
    EXPECTED_CERTIFICATE_DIGEST,
    EXPECTED_COUNTS,
    OUTPUT,
    verify,
)


def test_frozen_record0_a1_source2_theorem_and_scope() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 0,
        "schedule_a": 1,
        "source_orbit_index": 2,
        "representative_A4": [0, 1, 2, 3],
        "production_terminal_status": "INFEASIBLE",
        "independent_terminal_status": "unsat",
        "fixed_source_compact_incidence_no_go": True,
        "remaining_a1_source_orbits": [0, 1, 3],
        "whole_record0_a1_schedule_excluded": False,
        "whole_shadow_compact_coxeter_no_go": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_digest": EXPECTED_CERTIFICATE_DIGEST,
    }


def test_two_independent_models_match_every_frozen_count() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    evidence = document["evidence"]
    assert evidence["production_terminal_status"] == "INFEASIBLE"
    assert evidence["independent_terminal_status"] == "unsat"
    assert evidence["production_counts"] == EXPECTED_COUNTS
    assert evidence["independent_counts"] == EXPECTED_COUNTS
    assert evidence["independent_reconstruction_exact_count_agreement"] is True


def test_certificate_does_not_overclaim_schedule_shadow_or_volume() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    fixed = document["fixed_source"]
    theorem = document["theorem"]
    assert fixed["a1_source_orbit_count"] == 4
    assert fixed["other_source_orbits_not_closed_by_this_certificate"] == [0, 1, 3]
    assert theorem["excluded_a1_source_orbits"] == [2]
    assert theorem["remaining_a1_source_orbits"] == [0, 1, 3]
    assert theorem["whole_record0_a1_schedule_excluded"] is False
    assert theorem["whole_shadow_compact_coxeter_no_go"] is False
    assert theorem["finite_volume_noncompact_chamber_excluded"] is False
