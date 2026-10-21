import json

from coxeter_polytopes.combinatorics.h9_14_record18_coxeter_labels import (
    audit_record18_compact_coxeter_labels,
)
from scripts.certify_h9_14_record18_coxeter_labels import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)
from scripts.independent_h9_14_record18_coxeter_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay,
)


def test_record18_selected_nerve_compact_gate_is_exhaustive() -> None:
    audit = audit_record18_compact_coxeter_labels()
    assert audit.actual_vertex_count == 110
    assert audit.symmetry_group_order == 8
    assert audit.full_left_state_count == audit.full_right_state_count == 420
    assert audit.left_orbit_representative_count == 57
    assert audit.left_state_orbit_size_histogram == ((4, 9), (8, 48))
    assert audit.quotient_single_edge_trial_count == 598_500
    assert audit.full_single_edge_trial_count == 4_410_000
    assert audit.quotient_single_edge_passing_count == 3_780
    assert audit.passing_triple_orbit_count == 3_660
    assert audit.full_single_edge_passing_count == 28_800
    assert audit.stronger_label_survivors == ((4, 0), (5, 0))
    assert audit.quotient_residual_double_edge_trial_count == 0
    assert audit.unrestricted_compact_coxeter_no_go is True


def test_record18_residual_blocks_have_exact_negative_index_two() -> None:
    audit = audit_record18_compact_coxeter_labels()
    assert audit.exceptional_full_class_histogram == (
        ("A5_A5_type_1_1", 14_400),
        ("A5_A5_type_4_1", 14_400),
    )
    assert all(
        witness.signs == "++++-----+" and witness.inertia == (8, 2, 0)
        for witness in audit.exact_index_two_witnesses
    )


def test_record18_compact_gate_hostile_replay_and_frozen_certificate() -> None:
    independent = replay()
    assert independent["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    assert independent["unrestricted_compact_coxeter_no_go"] is True
    assert independent["finite_volume_noncompact_checked"] is False
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 18,
        "selected_nerve_only": True,
        "unrestricted_compact_coxeter_no_go": True,
        "finite_volume_noncompact_checked": False,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
