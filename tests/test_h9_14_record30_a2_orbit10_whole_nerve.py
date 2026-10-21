from __future__ import annotations

import json

from scripts.explore_h9_14_record30_a2_orbit10_cross_supersets import (
    ARCHIVE,
    EXPECTED_CERTIFICATE_SHA256,
    verify,
)
from scripts.independent_h9_14_record30_a2_orbit10_whole_nerve import replay


def test_record30_a2_orbit10_whole_nerve_certificate() -> None:
    frozen = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert verify(frozen) == {
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
        "first_empty_connected_shell": 5,
        "literal_connected_shell_counts": {
            "2": 36_288,
            "3": 82_944,
            "4": 24_624,
        },
        "residual_internal_passing_count": 0,
        "sparse_state_count": 648,
        "whole_orbit10_nerve_compact_coxeter_impossible": True,
    }


def test_record30_a2_orbit10_hostile_whole_nerve_replay() -> None:
    hostile = replay()
    assert hostile["agrees_with_frozen_production"] is True
    assert hostile["implementation"]["imports_production_enumerator"] is False
    assert hostile["implementation"]["imports_coxeter_label_module"] is False
    assert hostile["residual_internal_passing_count"] == 0
    assert hostile["connected_shell_counts_per_state"] == {
        "2": 56,
        "3": 128,
        "4": 38,
    }
    assert hostile["first_empty_connected_shell"] == 5
    assert hostile["zero_second_prime_fixed_minor_count"] == 0
    assert hostile["whole_orbit10_nerve_compact_coxeter_impossible"] is True


def test_record30_a2_orbit10_whole_nerve_scope() -> None:
    frozen = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    scope = frozen["scope"]
    assert scope["whole_orbit10_nerve_compact_excluded"] is True
    assert scope["whole_record30_a2_schedule_excluded"] is False
    assert scope["whole_record30_shadow_excluded"] is False
    assert scope["finite_volume_noncompact_statement"] is False
    assert frozen["conclusion"]["new_hyperbolic_polytope"] is False
