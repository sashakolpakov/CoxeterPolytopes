"""Hostile regressions for the independent unconditional synthesis audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

from scripts.audit_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness import (
    ALL_INPUT_PINS,
    AUDIT_REPORT,
    BLOWUP_MAX,
    CLASSICAL_AUDIT_PINS,
    CLASSICAL_OWNER_PINS,
    FIELD_DEGREE_MAX,
    OUTPUT,
    PRIMARY_SOURCE_PINS,
    RANK_MAX,
    ROOT,
    SYNTHESIS_OWNER_PINS,
    TRACE_RANK_MAX,
    WEAK_AUDIT_PINS,
    WEAK_OWNER_PINS,
    arithmetic_bound_and_box_replay,
    branch_archive_cross_checks,
    build_payload,
    chow_composition_replay,
    digest,
    encoded_audit,
    finite_gallery_calibration,
    local_link_audit,
    provenance_replay,
    quadratic_eigenslice_replay,
    sha256_file,
    synthesis_scope_audit,
)


AUDITOR = (
    ROOT
    / "scripts/audit_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.py"
)
AUDIT_ARTIFACT_PINS = {
    "report": "75e93b1938eaf0d89ba4d382f61ca089eeb86ea08102918845d30d6d582d5a4f",
    "script": "5b581f6e39ce085ca2faebeda6453bd08d21854dab78be21bd7fcfcb71cae623",
    "archive": "004a48b347a3b21163bdd008be6ba00e6dd2e03cded900f7fe366333642f8b37",
}


def test_all_five_terminal_four_artifact_packages_are_pinned() -> None:
    assert len(SYNTHESIS_OWNER_PINS) == 4
    assert len(CLASSICAL_OWNER_PINS) == 4
    assert len(CLASSICAL_AUDIT_PINS) == 4
    assert len(WEAK_OWNER_PINS) == 4
    assert len(WEAK_AUDIT_PINS) == 4
    assert len(ALL_INPUT_PINS) == 21  # five packages plus the scope no-go
    assert len(set(ALL_INPUT_PINS)) == 21


def test_literal_input_scope_and_repaired_gram_wording() -> None:
    replay = synthesis_scope_audit()
    assert replay["all_checks_pass"]
    assert "admissible Vinberg lattice-form" in replay["input"]
    assert all(replay["marker_checks"].values())
    assert all(replay["machine_checks"].values())
    assert replay["link_audit"]["all_local_links_resolve"]


def test_noncommuting_chow_replay_fixes_composition_order() -> None:
    replay = chow_composition_replay()
    checks = replay["checks"]
    assert replay["all_checks_pass"]
    assert replay["composition_rule"] == "C_A J C_B = C_(BA)"
    assert replay["standard_action_rule"] == "Gamma_g o Gamma_h = Gamma_(gh)"
    assert checks["noncommuting_pair"]
    assert checks["standard_composition_order"]
    assert checks["reversed_order_is_detectably_wrong"]
    assert checks["reflection_cycles_self_transpose"]
    assert checks["reflection_cycles_square_to_projector_cycle"]
    assert checks["m2_relation"] and checks["m3_relation"]
    assert checks["parabolic_nilpotency_index_three"]


def test_quadratic_field_slice_separates_raw_and_projected_grams() -> None:
    replay = quadratic_eigenslice_replay()
    checks = replay["checks"]
    assert replay["all_checks_pass"]
    assert replay["field"] == "Q(sqrt(3))"
    assert replay["native_rank"] == 3
    assert replay["trace_rank"] == 6
    assert replay["trace_signature"] == [1, 5, 0]
    assert replay["raw_wall_gram"] != replay["projected_wall_gram"]
    assert replay["raw_normalized_pair_square"] != (
        replay["projected_normalized_pair_square"]
    )
    assert checks["trace_pairings_equal_distinguished_pairings"]
    assert checks["projected_wall_gram_is_distinguished_gram"]
    assert checks["raw_trace_gram_differs_from_projected_gram"]


def test_29_25_750_2997_and_finite_vieta_box_are_recomputed() -> None:
    replay = arithmetic_bound_and_box_replay()
    assert replay["all_checks_pass"]
    assert RANK_MAX == 30
    assert FIELD_DEGREE_MAX == 25
    assert TRACE_RANK_MAX == 750
    assert BLOWUP_MAX == 2997
    assert replay["simple_cycle_slots"] == 24893912605687593731774060679329
    assert all(replay["checks"].values())


def test_group_to_chamber_gallery_contract_is_calibrated() -> None:
    replay = finite_gallery_calibration()
    assert replay["all_checks_pass"]
    assert replay["group_order"] == 6
    assert len(replay["chamber_sign_patterns"]) == 6
    assert replay["checks"]["one_orbit_meets_all_six_chambers"]
    assert "reflection maximality forces W(P)=R(L)" in replay["general_bridge"]
    assert "complete chambers form one reflection-group orbit" in replay["general_bridge"]


def test_provenance_graph_has_no_weak_carrier_to_finiteness_path() -> None:
    replay = provenance_replay()
    checks = replay["checks"]
    assert replay["all_checks_pass"]
    assert replay["forbidden_edge"] == [
        "weak_chow_carrier",
        "finite_maximal_chambers",
    ]
    assert replay["finiteness_derived_from_weak_carrier"] is False
    assert replay["uniform_2997_imports_classical_bounds"] is True
    assert checks["weak_carrier_does_not_reach_native_pairs"]
    assert checks["weak_carrier_does_not_reach_similarity_hosts"]
    assert checks["weak_carrier_does_not_reach_chamber_finiteness"]
    assert checks["factorized_finish_requires_maximality"]


def test_branch_archives_and_full_payload_pass() -> None:
    assert all(branch_archive_cross_checks().values())
    payload = build_payload()
    assert payload["verdict"] == "PASS_AFTER_SCOPE_REPAIRS"
    assert payload["all_checks_pass"]
    assert all(payload["checks"].values())
    assert len(payload["repairs_verified"]) == 3
    assert payload["replay_boundary"] == {
        "imports_synthesis_owner_producer": False,
        "imports_classical_owner_or_audit_producer": False,
        "imports_weak_owner_or_audit_producer": False,
        "published_literature_machine_proved": False,
        "temporary_primary_sources_required_for_clean_replay": False,
    }


def test_primary_source_hashes_match_when_cached() -> None:
    for record in PRIMARY_SOURCE_PINS.values():
        path = Path(record["path"])
        if path.exists():
            assert sha256_file(path) == record["sha256"]


def test_auditor_does_not_import_any_owner_or_branch_producer() -> None:
    source = AUDITOR.read_text(encoding="utf-8")
    forbidden = (
        "from scripts.certify_unconditional_ns_chow",
        "from scripts.certify_unconditional_maximal",
        "from scripts.certify_universal_weak",
        "from scripts.audit_unconditional_maximal",
        "from scripts.audit_universal_weak",
        "rational_surface_carrier import",
    )
    assert not any(needle in source for needle in forbidden)


def test_report_links_math_delimiters_and_bytes_are_clean() -> None:
    links = local_link_audit(AUDIT_REPORT)
    assert links["all_local_links_resolve"]
    assert links["broken_local_links"] == []
    assert links["links_escaping_repository"] == []
    text = AUDIT_REPORT.read_text(encoding="utf-8")
    lines = [line.strip() for line in text.splitlines()]
    assert lines.count("\\[") == lines.count("\\]")
    raw = AUDIT_REPORT.read_bytes()
    forbidden = set(range(0, 9)) | {11, 12} | set(range(14, 32)) | {127}
    assert not any(byte in forbidden for byte in raw)
    observed = {
        "report": hashlib.sha256(AUDIT_REPORT.read_bytes()).hexdigest(),
        "script": hashlib.sha256(AUDITOR.read_bytes()).hexdigest(),
        "archive": hashlib.sha256(OUTPUT.read_bytes()).hexdigest(),
    }
    assert observed == AUDIT_ARTIFACT_PINS


def test_archive_is_exact_replay_and_does_not_self_pin() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert archived["payload"] == build_payload()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert OUTPUT.read_text(encoding="utf-8") == encoded_audit()
    assert str(OUTPUT.relative_to(ROOT)) not in archived["payload"][
        "terminal_input_sha256"
    ]


def test_cli_check_output() -> None:
    subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.audit_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness",
            "--check-output",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
