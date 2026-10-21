from __future__ import annotations

from hashlib import sha256
import json

from scripts.certify_h9_14_record30_a1_schedule_no_go import ARCHIVE, verify


EXPECTED_ARCHIVE_SHA256 = (
    "65152e8f4b687f487852932a5dae3896d189a56d1ed99c6764f7389aaea8d341"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "932d69e5646106785171562174e5423329ac09137e0224ac65e6ce436187067d"
)
EXPECTED_LITERAL_SET_SHA256 = (
    "628d18323217fba80cea165f79e2ff047f9bc9ed06c06aa62265d89c5a43fa07"
)


def _certificate() -> dict[str, object]:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def test_frozen_a1_certificate_and_scope() -> None:
    certificate = _certificate()
    assert sha256(ARCHIVE.read_bytes()).hexdigest() == EXPECTED_ARCHIVE_SHA256
    assert verify(certificate) == {
        "record_index": 30,
        "fixed_positive_A4_count": 1,
        "necessary_incidence_target_orbit_count": 3,
        "literal_target_count": 168,
        "terminal_finite_label_state_count": 0,
        "record30_a1_compact_schedule_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }


def test_nauty_and_no_nauty_exhaustions_are_literal_set_identical() -> None:
    certificate = _certificate()
    comparison = certificate["incidence_cross_audit"]
    assert comparison["production_and_independent_common_counters_identical"]
    assert comparison["production_and_independent_literal_target_sets_identical"]
    assert comparison["literal_target_count"] == 168
    assert comparison["literal_target_set_sha256_all_orbits"] == (
        EXPECTED_LITERAL_SET_SHA256
    )
    production = certificate["production_nauty_incidence_replay"]
    independent = certificate["independent_explicit_incidence_replay"]
    assert [row["canonical_state_counts"] for row in production] == [
        [1, 11, 256, 9_229, 309_175],
        [1, 9, 153, 4_123, 159_275, 1_399_930],
        [0, 0, 0, 0, 0, 0, 0],
    ]
    assert [row["canonical_state_counts"] for row in independent] == [
        row["canonical_state_counts"] for row in production
    ]
    assert [row["target_orbit_sizes"] for row in independent] == [
        [48],
        [80, 40],
        [],
    ]


def test_two_active_branch_and_all_three_label_targets_are_killed() -> None:
    certificate = _certificate()
    witness = certificate["independent_explicit_incidence_replay"][2][
        "strict_link_witness"
    ]
    assert witness == {
        "subset": [0, 1],
        "polytope_dimension": 7,
        "f_prefix": [1, 10, 43, 104],
        "h": [1, 3, 4, 4],
        "g_without_g0": [2, 1, 0],
        "strict_entry": 3,
        "observed": 0,
        "equality_forced_missing_face_minimum_order": 6,
        "inherited_maximum_missing_face_order": 5,
    }
    labels = certificate["independent_complete_label_replay"]
    assert labels["terminal_finite_label_state_count_across_targets"] == 0
    assert labels["passing_state_count_transcripts"] == [
        [108, 84, 0],
        [108, 84, 36, 108, 0],
        [108, 84, 0],
    ]
    assert certificate["theorem_scope"] == {
        "record30_a1_compact_schedule_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "whole_record30_shadow_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
    }
