from __future__ import annotations

import json

from scripts.certify_h9_14_record14_a0_whole_shadow_no_go import (
    EXPECTED_CERTIFICATE_DIGEST,
    EXPECTED_LABEL_STATE_SEQUENCES,
    EXPECTED_TARGET_GROUPS,
    OUTPUT,
    verify,
)
from scripts.explore_h9_14_record14_a0_downstream import full_source_orbit_key
from scripts.explore_h9_14_record14_fixed_orbit_z3 import mask
from scripts.independent_h9_14_record14_a0_downstream_audit import (
    EXPECTED_GROUPS as INDEPENDENT_EXPECTED_GROUPS,
    EXPECTED_SUMMARY_SHA256 as INDEPENDENT_EXPECTED_SUMMARY_SHA256,
    OUTPUT as INDEPENDENT_OUTPUT,
    replay as independent_replay,
)


def _move(value: int, action: tuple[int, ...]) -> int:
    return sum(1 << action[vertex] for vertex in range(14) if value & (1 << vertex))


def test_record14_frozen_whole_shadow_theorem_and_scope() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 14,
        "schedule_a": 0,
        "generator_lex_target_count": 52,
        "full_target_orbit_count": 3,
        "terminal_compact_label_state_count": 0,
        "closed_schedules": [0, 1, 2, 3, 4, 5, 6],
        "remaining_schedules": [],
        "whole_shadow_compact_coxeter_no_go": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_digest": EXPECTED_CERTIFICATE_DIGEST,
    }


def test_record14_a0_literal_equality_and_full_orbit_partition() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    evidence = document["evidence"]
    assert evidence["production_z3_generator_lex_model_count"] == 52
    assert evidence["independent_cp_sat_generator_lex_model_count"] == 52
    assert evidence["generator_lex_slice_pure_target_count"] == 52
    assert evidence["production_independent_literal_equality"] is True
    assert tuple(
        tuple(group) for group in evidence["full_target_orbit_groups"]
    ) == EXPECTED_TARGET_GROUPS
    assert sorted(index for group in EXPECTED_TARGET_GROUPS for index in group) == list(range(52))
    assert evidence["target_orbit_size_histogram"] == {"14": 1, "18": 1, "20": 1}


def test_record14_a0_unrestricted_label_extinction() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    audits = document["evidence"]["unrestricted_compact_label_audits"]
    assert tuple(
        tuple(stage["passing_state_count"] for stage in audit["stages"])
        for audit in audits
    ) == EXPECTED_LABEL_STATE_SEQUENCES
    assert all(audit["rank_nine_actual_vertex_count"] == 110 for audit in audits)
    assert all(audit["terminal_sparse_label_state_count"] == 0 for audit in audits)
    assert all(audit["unrestricted_compact_label_no_go"] for audit in audits)


def test_record14_a0_canonical_key_is_invariant_under_source_generators() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    raw = document["evidence"]["generator_lex_slice_pure_targets"][0]
    family = tuple(sorted(mask(face) for face in raw))
    active = list(range(14))
    active[0], active[1] = active[1], active[0]
    isolate = list(range(14))
    isolate[3], isolate[4] = isolate[4], isolate[3]
    key = full_source_orbit_key(family)
    assert full_source_orbit_key(tuple(reversed(family))) == key
    assert full_source_orbit_key(tuple(_move(face, tuple(active)) for face in family)) == key
    assert full_source_orbit_key(tuple(_move(face, tuple(isolate)) for face in family)) == key


def test_record14_a0_hostile_downstream_replay() -> None:
    frozen = json.loads(INDEPENDENT_OUTPUT.read_text(encoding="utf-8"))
    assert independent_replay() == frozen
    assert frozen["summary_sha256"] == INDEPENDENT_EXPECTED_SUMMARY_SHA256
    assert tuple(
        tuple(group) for group in frozen["full_target_orbit_groups"]
    ) == INDEPENDENT_EXPECTED_GROUPS
    assert frozen["full_target_orbit_count"] == 3
    assert frozen["generator_lex_target_count"] == 52
    assert frozen["all_target_orbits_unrestricted_compact_label_excluded"] is True
    assert [
        [stage["passing_state_count"] for stage in audit["stages"]]
        for audit in frozen["label_audits"]
    ] == [[420, 96, 144, 16, 0], [420, 96, 8, 0], [420, 96, 8, 0]]
