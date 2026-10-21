from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.h9_14_record43 import TARGET_MISSING_FIVES
from coxeter_polytopes.combinatorics.h9_14_record43_all_completions import (
    EXPECTED_SCHEDULES,
    SECOND_TARGET_MISSING_FIVES,
    audit_record43_all_completions,
)
from coxeter_polytopes.combinatorics.h9_14_record43_second_nerve_coxeter_no_go import (
    audit_record43_second_nerve_labels,
)
from scripts.certify_h9_14_record43_all_completions import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    SUPPORT,
    verify,
)
from scripts.independent_h9_14_record43_all_completions_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay,
)


ROOT = Path(__file__).resolve().parents[1]
HOSTILE_SOURCE = (
    ROOT / "scripts/independent_h9_14_record43_all_completions_audit.py"
)
EXPECTED_ARCHIVE_FILE_SHA256 = (
    "68ae21860f4eea6a6a1eca978f4a07df30415b698de6d4ab0783d89cb3d7a85e"
)


@pytest.fixture(scope="session")
def production_audit():
    source = json.loads(SUPPORT.read_text(encoding="utf-8"))
    record = next(
        row for row in source["records"] if row["graph6"] == "M?????????O?o?w??"
    )
    schedules_json = json.dumps(
        record["dehn_sommerville_schedule"],
        sort_keys=True,
        separators=(",", ":"),
    )
    return audit_record43_all_completions(schedules_json)


@pytest.fixture(scope="session")
def second_label_audit():
    return audit_record43_second_nerve_labels()


@pytest.fixture(scope="session")
def hostile_replay():
    return replay()


def test_all_schedules_reduce_to_two_full_graph_orbits(production_audit) -> None:
    assert production_audit.verified
    assert production_audit.exact_residual_orbit_count == 2
    assert production_audit.schedule.rows == EXPECTED_SCHEDULES
    assert production_audit.schedule.surviving_missing_four_counts == (0, 1)
    assert production_audit.schedule.excluded_missing_four_counts == (2, 3)
    assert production_audit.a1.independent_a4_count == 690
    assert production_audit.a1.interval_feasible_a4_count == 378
    assert production_audit.a1.aggregate_a5_feasible_count == 28
    assert production_audit.a1.representative_vertex_feasible_pair_count == 540
    assert production_audit.a1.representative_small_transversal_pair_count == 540
    assert production_audit.a1.representative_transversal_size_histogram == ((4, 540),)
    assert production_audit.a0.enumeration_path_count == 1073
    assert production_audit.a0.exact_mask_state_count == 887
    assert production_audit.a0.canonical_quotient_orbit_count == 394
    assert production_audit.a0.target_path_count == 2
    assert production_audit.a0.target_exact_mask_state_count == 2
    assert production_audit.a0.target_quotient_orbit_count == 2
    assert production_audit.a0.active_graph_automorphism_count == 2
    assert production_audit.a0.target_full_graph_orbit_count == 2
    decoded = tuple(
        tuple(
            tuple(vertex for vertex in range(14) if mask >> vertex & 1)
            for mask in state
        )
        for state in production_audit.a0.target_missing_five_states
    )
    assert decoded == (TARGET_MISSING_FIVES, SECOND_TARGET_MISSING_FIVES)


def test_second_nerve_uses_correct_distinguished_support_orbits(
    second_label_audit,
) -> None:
    assert second_label_audit.unrestricted_compact_coxeter_no_go
    assert second_label_audit.verified_nerve_automorphism_count == 576
    assert second_label_audit.full_left_placement_count == 420
    assert second_label_audit.full_right_placement_count == 420
    assert second_label_audit.left_orbit_size_histogram == ((12, 5), (24, 15))
    assert second_label_audit.right_orbit_size_histogram == ((12, 5), (24, 15))
    assert second_label_audit.cross_gate.left_state_count == 20
    assert second_label_audit.cross_gate.right_state_count == 20
    assert second_label_audit.cross_gate.label_three_test_count == 10_000
    assert len(second_label_audit.cross_gate.passing_trials) == 55
    assert second_label_audit.stronger_single_edge_survivors == ((4, 0), (5, 0))
    assert second_label_audit.raw_double_edge_trial_count == 120_000
    assert second_label_audit.monotonicity_residual_double_trial_count == 0
    assert second_label_audit.double_edge_passing_count == 0
    assert second_label_audit.exceptional_class_histogram == (
        ("A5_A5_type_1_1", 25),
        ("A5_A5_type_4_1", 30),
    )
    assert tuple(w.inertia for w in second_label_audit.exact_index_two_witnesses) == (
        (8, 2, 0),
        (8, 2, 0),
    )


def test_no_project_import_hostile_replay_closes_the_same_theorem(
    hostile_replay,
) -> None:
    source = HOSTILE_SOURCE.read_text(encoding="utf-8")
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source
    assert hostile_replay["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    assert hostile_replay["a0_census"]["scheduled_target_full_graph_orbit_count"] == 2
    assert hostile_replay["a1_census"]["surviving_addition_pair_count"] == 0
    gate = hostile_replay["second_target_independent_label_gate"]
    assert gate["verified_nerve_automorphism_count"] == 576
    assert gate["representative_count_per_support"] == 20
    assert gate["label_three_test_count"] == 10_000
    assert gate["label_three_passing_orbit_count"] == 55
    assert gate["monotonicity_residual_two_edge_trial_count"] == 0
    assert gate["unrestricted_compact_coxeter_no_go"] is True
    first_gate = hostile_replay["first_target_compact_coxeter_archive"]
    assert first_gate["normalizing_support_symmetry_orders"] == [120]
    assert first_gate["unrestricted_compact_coxeter_no_go"] is True
    assert "principal submatrix" in first_gate["principal_block_lemma"]
    assert hostile_replay["scope"]["finite_volume_noncompact_checked"] is False


def test_frozen_whole_shadow_certificate_rebuilds_compact_only() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert sha256(OUTPUT.read_bytes()).hexdigest() == EXPECTED_ARCHIVE_FILE_SHA256
    summary = verify(document)
    assert summary == {
        "record_index": 43,
        "exact_residual_incidence_orbits": 2,
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert document["scope"]["compact_only"] is True
    assert document["scope"]["finite_volume_noncompact_chamber_excluded"] is False
    chain = document["first_unrestricted_compact_label_chain"]
    assert chain["aggregate_kind"] == (
        "h9_14_seven_selected_nerve_compact_coxeter_no_go_reaudited"
    )
    assert chain["normalizing_support_symmetry_orders"] == [120]
    assert chain["unrestricted_compact_coxeter_no_go"] is True
    assert "principal submatrix" in chain["principal_block_lemma"]
    assert "inertia" in chain["principal_block_lemma"]
