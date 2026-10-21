import json

import coxeter_polytopes.combinatorics.h9_14_selected_nerve_coxeter_no_go as gate
from scripts.certify_h9_14_selected_nerve_coxeter_no_go import OUTPUT, verify
from scripts.independent_h9_14_selected_nerve_coxeter_no_go_audit import (
    compare_archive,
    replay,
)


def test_exact_radical_signs_derive_the_exceptional_block_inertias() -> None:
    record42 = gate._record42_witness((2, 3, 4, 5), (6, 7, 8, 9, 10))
    a5_classes = gate._a5_witnesses((2, 3, 4, 5, 6), (1, 7, 8, 9, 10))
    assert record42.signs == "+++++++-+"
    assert record42.inertia == (7, 2, 0)
    assert [(item.signs, item.inertia) for item in a5_classes] == [
        ("++++-----+", (8, 2, 0)),
        ("++++-----+", (8, 2, 0)),
    ]


def test_selected_support_audits_are_unrestricted_after_scope_reaudit() -> None:
    audits = gate.audit_selected_nerve_coxeter_no_gos()
    assert [audit.record_index for audit in audits] == [30, 32, 42, 43, 44, 47]
    assert [audit.vertex_count for audit in audits] == [128, 110, 100, 90, 90, 90]
    assert all(audit.unrestricted_compact_coxeter_no_go for audit in audits)

    by_index = {audit.record_index: audit for audit in audits}
    assert by_index[30].support_gates[0].label_three_test_count == 1_296
    assert by_index[30].support_gates[0].passing_trials == ()
    assert by_index[30].raw_double_edge_trial_count == 9_720
    assert by_index[30].monotonicity_residual_double_trial_count == 0
    assert by_index[30].normalizing_support_symmetry_orders == (24, 24)

    record42 = by_index[42]
    assert record42.support_gates[0].label_three_test_count == 900
    assert len(record42.support_gates[0].passing_trials) == 2
    assert record42.stronger_single_edge_survivors == ((4, 0), (5, 0))
    assert record42.raw_double_edge_trial_count == 8_550
    assert record42.monotonicity_residual_double_trial_count == 1
    assert record42.double_edge_passing_count == 0
    assert record42.exact_index_two_witnesses[0].inertia == (7, 2, 0)
    assert record42.normalizing_support_symmetry_orders == (24, 120)

    for index in (43, 44, 47):
        audit = by_index[index]
        assert len(audit.support_gates) == 3
        assert all(item.label_three_test_count == 52_500 for item in audit.support_gates)
        assert all(len(item.passing_trials) == 300 for item in audit.support_gates)
        assert audit.stronger_single_edge_survivors == ((4, 0), (5, 0))
        assert audit.raw_double_edge_trial_count == 1_890_000
        assert audit.monotonicity_residual_double_trial_count == 0
        assert audit.double_edge_passing_count == 0
        assert audit.exceptional_class_histogram == (
            ("A5_A5_type_1_1", 360),
            ("A5_A5_type_4_1", 540),
        )
        assert [item.inertia for item in audit.exact_index_two_witnesses] == [
            (8, 2, 0),
            (8, 2, 0),
        ]
        assert audit.normalizing_support_symmetry_orders == (120,)


def test_frozen_seven_selected_nerve_certificate_replays_after_scope_reaudit() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary["selected_nerve_indices"] == [30, 32, 35, 42, 43, 44, 47]
    assert summary["selected_nerve_count"] == 7
    assert summary["compact_coxeter_realization_exists"] is False
    assert "seven frozen/displayed nerve completions" in summary["scope"]
    assert summary["reaudited_indices"] == [30, 42, 43, 44, 47]
    nonclaims = manifest["theorem"]["nonclaims"]
    assert any("graph shadow" in claim for claim in nonclaims)
    assert any("finite-volume noncompact" in claim for claim in nonclaims)


def test_hostile_replay_agrees_with_the_frozen_aggregate() -> None:
    result = replay()
    compare_archive(result)
    assert result["selected_nerve_indices"] == [30, 32, 35, 42, 43, 44, 47]
    assert result["reaudited_support_indices"] == [30, 42, 43, 44, 47]
    assert result["record32"]["A4_passing_count"] == 48
    assert result["record32"]["distinct_FIRST_A4_projection_count"] == 24
    assert result["record32"]["FIRST_A4_projection_multiplicity_histogram"] == {
        "2": 24
    }
    assert result["all_selected_nerves_compact_coxeter_excluded"] is True
    assert result["principal_block_scope_reaudit_passed"] is True
