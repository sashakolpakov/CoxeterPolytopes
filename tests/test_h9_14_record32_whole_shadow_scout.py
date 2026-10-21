from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import pytest

from scripts.certify_h9_14_record32_all_completions_candidate import (
    EXPECTED_CANDIDATE_FILE_SHA256,
    OUTPUT as CANDIDATE_OUTPUT,
    verify as verify_candidate,
)
from scripts.explore_h9_14_record32_a0_target_gate import audit_family
from scripts.explore_h9_14_record32_a1_completions import (
    A1_REPRESENTATIVES,
    audit_representative as audit_a1_representative,
)
from scripts.explore_h9_14_record32_a1_cross_gate import (
    TARGETS as A1_CROSS_TARGETS,
    audit_target as audit_a1_cross_target,
)
from scripts.explore_h9_14_record32_a2_additions import (
    A2_REPRESENTATIVES,
    audit_representative as audit_a2_representative,
)
from scripts.explore_h9_14_record32_a3_additions import audit_a3_representative
from scripts.explore_h9_14_record32_all_completions import (
    EXPECTED_SCHEDULES,
    schedule_audit,
    state_space_audit,
)


ROOT = Path(__file__).resolve().parents[1]
TARGETS = (
    ROOT / "data/enumerations/h9_14_record32_a0_all_targets_scout.json"
)
EXPECTED_TARGET_FILE_SHA256 = (
    "ae84f9c8dde0f66cd4d8fc8c0b9dee2afab985c7d484f2e052709ad3c3734938"
)


@pytest.fixture(scope="session")
def target_document() -> dict[str, object]:
    return json.loads(TARGETS.read_text(encoding="utf-8"))


def test_exact_schedules_and_lossless_a0_state_space() -> None:
    schedules = schedule_audit()
    assert tuple(
        (
            row["missing_four_count"],
            row["total_missing_independent_five_sets"],
            row["rank_nine_facet_count"],
            row["summed_vertex_link_g4"],
        )
        for row in schedules["rows"]
    ) == EXPECTED_SCHEDULES
    assert schedules["strict_compact_surviving_a"] == [0, 1, 2, 3]
    assert schedules["strict_compact_excluded_a"] == [4, 5]

    state = state_space_audit()
    assert state["active_graph_description"] == "P4 disjoint union K2"
    assert state["active_graph_automorphism_count"] == 4
    assert state["independent_active_mask_count"] == 24
    assert state["independent_active_mask_size_histogram"] == {
        0: 1,
        1: 6,
        2: 11,
        3: 6,
    }
    assert state["active_row_multisets_after_caps_before_cover_filter"] == 22_513
    assert state["active_row_multisets_after_exact_active_cover_filter"] == 14_731
    assert "No singleton-active-row normalization" in state["lossless_statement"]


def test_a0_scout_target_archive_is_frozen_but_not_a_theorem(
    target_document: dict[str, object],
) -> None:
    assert sha256(TARGETS.read_bytes()).hexdigest() == EXPECTED_TARGET_FILE_SHA256
    assert target_document["scheduled_facet_count"] == 114
    assert target_document["scheduled_target_path_count"] == 1_263
    assert target_document["exact_unordered_target_family_count"] == 1_025
    families = target_document["target_family_masks"]
    assert len(families) == 1_025
    assert len({tuple(row) for row in families}) == 1_025


def test_one_a0_ridge_candidate_has_an_uncapped_unrestricted_label_no_go(
    target_document: dict[str, object],
) -> None:
    family = target_document["target_family_masks"][6]
    audit = audit_family(6, family, state_cap=250_000, incidence_only=False)
    assert audit["minimal_transversal_size_histogram"] == {5: 114}
    assert audit["ridge_degree_histogram"] == {2: 513}
    assert audit["all_ridge_degrees_two"] is True
    gate = audit["sparse_label_gate"]
    assert gate["state_cap_exceeded"] is False
    assert gate["unrestricted_compact_label_no_go"] is True
    assert [
        (
            stage["support_index"],
            stage["compatible_extension_count_before_cap"],
            stage["passing_state_count"],
        )
        for stage in gate["stages"]
    ] == [(0, 420, 420), (1, 1_332, 96), (2, 384, 0)]


@pytest.fixture(scope="session")
def a2_audits() -> tuple[dict[str, object], ...]:
    return tuple(
        audit_a2_representative(index, representative)
        for index, representative in enumerate(A2_REPRESENTATIVES)
    )


def test_a2_additional_a5_closure_reaches_only_the_selected_nerve(
    a2_audits: tuple[dict[str, object], ...],
) -> None:
    assert [row["scheduled_target_count"] for row in a2_audits] == [0, 1, 0, 0, 0, 0]
    selected = a2_audits[1]
    assert selected["forced_five_shadow_count"] == 17
    assert selected["additional_minimal_A5_count"] == 1
    assert selected["scheduled_target_additions"] == [[[6, 7, 8, 9, 10]]]
    large = a2_audits[2]
    assert large["raw_unordered_addition_count"] == 16_328_620
    assert large["vertex_cap_feasible_count"] == 1_794_605
    assert large["no_transversal_of_order_at_most_four_count"] == 840
    assert large["facet_count_histogram"] == {78: 90, 81: 120, 84: 480, 87: 150}


def test_a3_unique_shadow_orbit_is_empty_before_labels() -> None:
    audit = audit_a3_representative()
    assert audit["forced_five_shadow_count"] == 22
    assert audit["additional_minimal_A5_count"] == 2
    assert audit["individual_candidate_count"] == 56
    assert audit["raw_unordered_addition_count"] == 1_540
    assert audit["vertex_cap_feasible_count"] == 540
    assert audit["no_transversal_of_order_at_most_four_count"] == 0
    assert audit["scheduled_target_count"] == 0


def test_a1_first_orbit_lossless_quotient_finds_the_unique_ridge_target() -> None:
    audit = audit_a1_representative(
        0, A1_REPRESENTATIVES[0], incidence_only=True
    )
    assert audit["individual_candidate_count"] == 664
    assert audit["raw_unordered_addition_count"] == 220_116
    assert audit["isolate_column_multiset_pair_path_count"] == 1_161
    assert audit["distinct_candidate_row_path_count"] == 1_102
    assert audit["no_transversal_of_order_at_most_four_path_count"] == 9
    assert audit["fixed_A4_stabilizer_target_orbit_count"] == 1
    target = audit["target_audits"][0]
    assert target["minimal_transversal_size_histogram"] == {5: 112}
    assert target["ridge_degree_histogram"] == {2: 504}


def test_second_a1_ridge_target_has_zero_single_cross_edge_survivors() -> None:
    audit = audit_a1_cross_target(A1_CROSS_TARGETS[1])
    assert audit["sparse_locally_elliptic_state_count"] == 864
    assert audit["minimum_nonright_cross_edge_count_from_inertia"] == 1
    assert audit["cross_pair_count"] == 20
    assert audit["single_m3_test_count"] == 17_280
    assert audit["single_m3_passing_count"] == 0
    assert audit["unrestricted_compact_label_no_go"] is True


def test_candidate_archive_refuses_theorem_promotion_before_hostile_replay() -> None:
    document = json.loads(CANDIDATE_OUTPUT.read_text(encoding="utf-8"))
    assert sha256(CANDIDATE_OUTPUT.read_bytes()).hexdigest() == (
        EXPECTED_CANDIDATE_FILE_SHA256
    )
    assert verify_candidate(document) == {
        "record_index": 32,
        "production_candidate_complete": True,
        "independent_replay_complete": False,
        "whole_graph_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "candidate_file_sha256": EXPECTED_CANDIDATE_FILE_SHA256,
    }
