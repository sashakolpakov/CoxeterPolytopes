from __future__ import annotations

import json

from scripts.certify_h9_14_record30_whole_shadow_compact_no_go import (
    COMPONENTS,
    OPEN_COMPACT_SHADOWS_AFTER_RECORD30,
    OUTPUT,
    fresh_record,
    schedule_table,
    verify_record,
)


def _load() -> dict[str, object]:
    value = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert isinstance(value, dict)
    return value


def test_record30_whole_shadow_compact_certificate() -> None:
    document = _load()
    summary = verify_record(document, require_archive_file=True)
    assert summary["compact_record30_excluded"] is True
    assert summary["closed_schedules"] == list(range(8))
    assert summary["all_nonnegative_a_closed"] is True
    assert summary["finite_volume_noncompact_checked"] is False
    assert summary["new_hyperbolic_example"] is False
    assert summary["compact_frontier"] == "133 + 9"
    assert summary["open_compact_shadows"] == (
        OPEN_COMPACT_SHADOWS_AFTER_RECORD30
    )
    assert summary["ordinary_incidence_frontier"] == "123 + 19"


def test_fresh_semantic_composition_matches_frozen_archive() -> None:
    assert fresh_record() == _load()


def test_schedule_table_and_cover_are_exact() -> None:
    rows = schedule_table()
    assert [row["a"] for row in rows] == list(range(8))
    assert [row["b"] for row in rows] == [8, 14, 20, 26, 32, 38, 44, 50]
    assert [row["facets"] for row in rows] == [132, 130, 128, 126, 124, 122, 120, 118]
    assert [row["global_g4"] for row in rows] == [8, 7, 6, 5, 4, 3, 2, 1]
    assert [row["summed_vertex_link_g4"] for row in rows] == [48, 42, 36, 30, 24, 18, 12, 6]
    covered = sorted(
        schedule
        for component in COMPONENTS
        for schedule in component["schedules"]
    )
    assert covered == list(range(8))

    exhaustion = _load()["schedule_exhaustion"]
    assert exhaustion["strict_global_g4_excludes_a_at_least_8"] is True
    assert exhaustion["strict_summed_vertex_link_g4_excludes"] == [6, 7]
    assert exhaustion["component_certificates_close"] == list(range(8))
    assert exhaustion["all_nonnegative_a_closed"] is True


def test_component_archives_are_file_and_internal_digest_pinned() -> None:
    archived = _load()["component_certificates"]
    assert len(archived) == len(COMPONENTS) == 6
    for expected, observed in zip(COMPONENTS, archived, strict=True):
        assert observed["id"] == expected["id"]
        assert observed["schedules"] == expected["schedules"]
        assert observed["archive"] == expected["archive"]
        assert observed["archive_sha256"] == expected["archive_sha256"]
        assert observed["certificate_sha256"] == expected["certificate_sha256"]
        assert observed["verifier"] == expected["verifier"]
        assert observed["verifier_sha256"] == expected["verifier_sha256"]


def test_branch_semantics_and_scope_are_not_overpromoted() -> None:
    document = _load()
    components = {
        row["id"]: row for row in document["component_certificates"]
    }
    assert components["a0_dual_exhaustion"]["semantic_summary"]["drat_verified"] is True
    assert components["a0_dual_exhaustion"]["semantic_summary"]["compact_label_survivor_count"] == 0
    assert components["a1_schedule"]["semantic_summary"]["literal_target_count"] == 168
    assert components["a1_schedule"]["semantic_summary"]["terminal_finite_label_state_count"] == 0
    assert components["a2_schedule"]["semantic_summary"]["pure_target_orbit_count"] == 77
    assert components["a2_schedule"]["semantic_summary"]["exceptional_target_count"] == 3
    assert components["a3_schedule"]["semantic_summary"]["scheduled_126_facet_target_count"] == 0
    assert components["a4_schedule"]["semantic_summary"]["scheduled_124_facet_target_count"] == 0
    assert components["a5_a7_positive_scout"]["semantic_summary"]["a5_compact_schedule_excluded"] is True
    assert components["a5_a7_positive_scout"]["semantic_summary"]["strict_sum_excluded_schedules"] == [6, 7]

    assert document["compact_record30_excluded"] is True
    assert document["whole_shadow_compact_coxeter_theorem"] is True
    assert document["compact_only"] is True
    assert document["finite_volume_noncompact_checked"] is False
    assert document["finite_volume_noncompact_chamber_excluded"] is False
    assert document["new_hyperbolic_coxeter_polytope"] is False


def test_a0_independent_orbit_and_deep_replay_provenance_are_present() -> None:
    document = _load()
    hostile = document["a0_no_import_orbit_audit"]
    orbit = hostile["orbit_and_clause_replay"]
    assert orbit["target_orbit_count"] == 8
    assert orbit["blocked_seed_form_count"] == 63_960
    assert orbit["all_orbit_stabilizer_double_counts_exact"] is True
    assert orbit["cnf_tail_exactly_reconstructed"] is True

    final_hostile = hostile["final_hostile_audit"]
    assert final_hostile["audit_status"] == "PASS"
    assert final_hostile["remaining_fatal_or_material_gap"] is False
    assert final_hostile["independent_drat_replay_status"] == "VERIFIED"
    assert final_hostile["independent_drat_replay_seconds"] == 488.306
    assert final_hostile["checker_resolution_step_count"] == 34_109_825
    assert final_hostile["compact_only"] is True
    assert final_hostile["finite_volume_noncompact_statement"] is False

    a0 = document["component_certificates"][0]
    provenance = a0["deep_replay_provenance"]
    assert provenance["solver"] == "Kissat 4.0.4"
    assert provenance["solver_status"] == "UNSATISFIABLE"
    assert provenance["checker"] == "drat-trim"
    assert provenance["checker_status"] == "VERIFIED"
    assert provenance["checker_resolution_step_count"] == 34_109_825
    assert provenance["deep_mode_enforces_exact_solver_and_checker_binary_hashes"] is True
