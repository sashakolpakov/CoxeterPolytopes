import json

from coxeter_polytopes.combinatorics.h9_14_record31_coxeter_labels import (
    audit_record31_compact_coxeter_labels,
)
from scripts.certify_h9_14_record31_coxeter_labels import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    verify,
)
from scripts.independent_h9_14_record31_coxeter_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay,
)


def test_record31_selected_nerve_compact_gate_is_exhaustive() -> None:
    audit = audit_record31_compact_coxeter_labels()
    assert audit.actual_vertex_count == 122
    assert audit.symmetry_group_order == 16
    assert audit.full_left_state_count == audit.full_right_state_count == 420
    assert audit.left_orbit_representative_count == 110
    assert audit.right_orbit_representative_count == 110
    assert audit.left_state_orbit_size_histogram == ((2, 10), (4, 100))
    assert audit.right_state_orbit_size_histogram == ((2, 10), (4, 100))
    assert audit.quotient_gate.label_three_test_count == 302_500
    assert len(audit.quotient_gate.passing_trials) == 1_932
    assert audit.stronger_label_survivors == ((4, 0), (5, 0))
    assert audit.quotient_residual_double_edge_trial_count == 36
    assert audit.quotient_residual_double_edge_passing_count == 0
    assert audit.unrestricted_compact_coxeter_no_go is True


def test_record31_inertia_and_residual_lanner_branches_close() -> None:
    audit = audit_record31_compact_coxeter_labels()
    assert audit.exceptional_class_histogram == (
        ("A5_A5_type_0_1", 36),
        ("A5_A5_type_1_1", 900),
        ("A5_A5_type_1_1_opposite_roots", 36),
        ("A5_A5_type_4_1", 960),
    )
    assert audit.index_two_single_edge_trial_count == 1_860
    assert all(
        witness.signs == "++++-----+" and witness.inertia == (8, 2, 0)
        for witness in audit.exact_index_two_witnesses
    )
    assert audit.residual_lanner_csp_trial_count == 72
    assert audit.residual_lanner_csp_completing_count == 0
    assert audit.residual_lanner_csp_node_histogram == (
        (1, 50),
        (4, 12),
        (23, 2),
        (38, 2),
        (63, 2),
        (68, 2),
        (113, 2),
    )


def test_record31_compact_gate_hostile_replay_and_frozen_certificate() -> None:
    independent = replay()
    assert independent["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    assert independent["unrestricted_compact_coxeter_no_go"] is True
    assert independent["finite_volume_noncompact_checked"] is False
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert verify(document) == {
        "record_index": 31,
        "selected_nerve_only": True,
        "unrestricted_compact_coxeter_no_go": True,
        "finite_volume_noncompact_checked": False,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }

