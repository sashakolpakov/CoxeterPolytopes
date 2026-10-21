from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.h9_14_record42_all_completions import (
    EXPECTED_LABEL_STAGE_TRANSCRIPT_SHA256,
    EXPECTED_SCHEDULES,
    EXPECTED_TARGET_FAMILY_SHA256,
    EXPECTED_TOPOLOGY_ORBIT_KEY_SHA256,
    audit_record42_all_completions,
)
from scripts.certify_h9_14_record42_all_completions import (
    EXPECTED_CERTIFICATE_SHA256,
    EXPECTED_RECORD42_SUPPORT_SYMMETRY_ORDERS,
    EXPECTED_SCOPE_REAUDIT_STATUS,
    EXPECTED_SELECTED_LABEL_KIND,
    OUTPUT,
    SUPPORT,
    verify,
)
from scripts.independent_h9_14_record42_all_completions_audit import (
    EXPECTED_SUMMARY_SHA256,
    replay,
)


ROOT = Path(__file__).resolve().parents[1]
HOSTILE_SOURCE = (
    ROOT / "scripts/independent_h9_14_record42_all_completions_audit.py"
)
EXPECTED_ARCHIVE_FILE_SHA256 = (
    "c9beee7a93965774bff273b89c95430a2b9504b9d094b4730f78033aa659e74c"
)
EXPECTED_HOSTILE_LABEL_TRANSCRIPT_SHA256 = (
    "0e726439a3f85977f4456b3184e74106050804c0f2e106456aad115df47b298e"
)


@pytest.fixture(scope="session")
def production_audit():
    source = json.loads(SUPPORT.read_text(encoding="utf-8"))
    record = [row for row in source["records"] if row["survives_second_layer"]][42]
    serialized = json.dumps(record, sort_keys=True, separators=(",", ":"))
    return audit_record42_all_completions(serialized)


@pytest.fixture(scope="session")
def hostile_replay():
    return replay()


def test_all_schedules_and_positive_branches_are_exhausted(production_audit) -> None:
    assert production_audit.verified
    assert production_audit.schedule.rows == EXPECTED_SCHEDULES
    assert production_audit.schedule.surviving_a == (0, 1, 2)
    assert production_audit.schedule.excluded_a == (3, 4)
    positive = production_audit.positive
    assert positive.individually_feasible_a4_count == 546
    assert positive.a2_unordered_pair_count == 148_785
    assert positive.a2_family_interval_count == 131_355
    assert positive.a2_forced_shadow_overrun_count == 58_275
    assert positive.a2_strict_count == 0
    assert positive.a1_strict_count == 126
    assert positive.a1_forced_product_count == 126
    assert positive.a1_incidence_orbit_count == 1
    assert positive.a1_representative_four == (2, 3, 4, 5)
    assert positive.a1_representative_added_five == (6, 7, 8, 9, 10)


def test_corrected_a0_census_covers_all_82_active_row_types(production_audit) -> None:
    a0 = production_audit.a0
    assert a0.active_row_configuration_count == 82
    assert a0.graph_cover_prefilter_path_count == 1_098_455
    assert a0.prefilter_repeated_row_path_count == 28_805
    assert a0.small_transversal_failure_count == 0
    assert a0.blocker_valid_path_count == 1_069_650
    assert a0.target_path_count == 35
    assert a0.exact_target_family_count == 15
    assert a0.target_family_sha256 == EXPECTED_TARGET_FAMILY_SHA256
    assert a0.topology_valid_family_indices == (0, 1, 2, 3, 6)
    assert a0.topology_orbit_count == 1
    assert a0.topology_orbit_key_sha256 == EXPECTED_TOPOLOGY_ORBIT_KEY_SHA256
    assert all(
        audit.minimal_transversal_histogram == ((5, 102),)
        for audit in a0.family_audits
    )


def test_complete_unrestricted_label_gate_kills_all_15_targets(
    production_audit,
) -> None:
    a0 = production_audit.a0
    assert len(a0.label_audits) == 15
    assert all(audit.compact_no_go for audit in a0.label_audits)
    assert a0.label_stage_transcript_sha256 == (
        EXPECTED_LABEL_STAGE_TRANSCRIPT_SHA256
    )
    representative = a0.label_audits[0]
    assert tuple(
        (
            stage.raw_extension_count,
            stage.compatible_extension_count,
            stage.passing_state_count,
        )
        for stage in representative.stages
    ) == ((420, 420, 420), (176_400, 1_332, 96), (40_320, 384, 0))


def test_hostile_no_project_import_replay_closes_the_same_theorem(
    hostile_replay,
) -> None:
    source = HOSTILE_SOURCE.read_text(encoding="utf-8")
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source
    assert hostile_replay["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    a0 = hostile_replay["a0_census"]
    assert a0["active_row_configuration_count"] == 82
    assert a0["cover_prefilter_paths"] == 1_098_455
    assert a0["small_transversal_failure_count"] == 0
    assert a0["exact_target_family_count"] == 15
    assert a0["ridge_valid_family_indices"] == [0, 1, 2, 3, 6]
    assert a0["ridge_valid_orbit_count"] == 1
    assert hostile_replay["a0_label_stage_transcript_sha256"] == (
        EXPECTED_HOSTILE_LABEL_TRANSCRIPT_SHA256
    )
    assert all(
        audit["unrestricted_compact_coxeter_no_go"]
        for audit in hostile_replay["a0_all_target_label_audits"]
    )
    product = hostile_replay["a1_product_label_chain"]
    assert product["aggregate_kind"] == EXPECTED_SELECTED_LABEL_KIND
    assert product["scope_reaudit_status"] == EXPECTED_SCOPE_REAUDIT_STATUS
    assert product["normalizing_support_symmetry_orders"] == (
        EXPECTED_RECORD42_SUPPORT_SYMMETRY_ORDERS
    )
    assert "principal submatrix" in product["principal_block_lemma"]
    assert "inertia" in product["principal_block_lemma"]
    assert product["single_edge_test_count"] == 900
    assert product["raw_two_edge_trial_count"] == 8_550
    assert product["two_edge_passing_count"] == 0
    assert product["exceptional_inertia"] == [7, 2, 0]
    assert product["unrestricted_compact_coxeter_no_go"] is True
    assert hostile_replay["theorem"] == {
        "whole_record42_graph_shadow_compact_coxeter_impossible": True,
        "a0_exact_target_families_excluded": 15,
        "a1_product_orbit_excluded": True,
        "finite_volume_noncompact_excluded": False,
    }


def test_frozen_whole_shadow_certificate_is_compact_only() -> None:
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert EXPECTED_ARCHIVE_FILE_SHA256 != "TO_BE_FROZEN"
    assert sha256(OUTPUT.read_bytes()).hexdigest() == EXPECTED_ARCHIVE_FILE_SHA256
    assert verify(document) == {
        "record_index": 42,
        "a0_exact_target_families_excluded": 15,
        "a1_incidence_orbits_excluded": 1,
        "whole_graph_shadow_compact_coxeter_excluded": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
    assert document["scope"]["compact_only"] is True
    assert document["scope"][
        "finite_volume_noncompact_chamber_excluded"
    ] is False
    product = document["a1_unrestricted_product_label_chain"]
    assert product["aggregate_kind"] == EXPECTED_SELECTED_LABEL_KIND
    assert product["scope_reaudit_status"] == EXPECTED_SCOPE_REAUDIT_STATUS
    assert product["normalizing_support_symmetry_orders"] == (
        EXPECTED_RECORD42_SUPPORT_SYMMETRY_ORDERS
    )
    assert "principal submatrix" in product["principal_block_lemma"]
    assert product["unrestricted_compact_coxeter_no_go"] is True
