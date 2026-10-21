"""Independent hostile tests for the universal weak Chow-action package."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys

from scripts.audit_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile import (
    DEPENDENCY_PINS,
    OUTPUT,
    OWNER_PINS,
    ROOT,
    build_audit,
    chow_orientation_replay,
    digest,
    eigenslice_replay,
    encoded_audit,
    four_square_bound_replay,
    leverage_replay,
    pin_audit,
    stable_lattice_replay,
)


AUDITOR = (
    ROOT
    / "scripts/audit_universal_weak_chow_action_carrier_and_leverage_no_go_"
    "independent_hostile.py"
)
AUDIT_REPORT = (
    ROOT
    / "results/universal_weak_chow_action_carrier_and_leverage_no_go_"
    "independent_hostile_audit.md"
)
AUDIT_PINS = {
    "results/universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.md": (
        "91345590b506ef7c7d4287be5a4d94c7563ba7613167ffbe34272e6b7fc422c1"
    ),
    "scripts/audit_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile.py": (
        "98765eaae927d9265cf69591f02c31d1d0adc71653936d2199f53eb75b16a347"
    ),
    "data/enumerations/universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.json": (
        "a9c19a9f2dce89c03182c9dc7ad8248a0ab75e32e5e58cd1bca6fac7683764f4"
    ),
}


def test_terminal_owner_and_dependency_bytes_are_pinned() -> None:
    assert pin_audit(OWNER_PINS) == OWNER_PINS
    assert pin_audit(DEPENDENCY_PINS) == DEPENDENCY_PINS
    source = AUDITOR.read_text(encoding="utf-8")
    assert "from scripts.certify_universal_weak" not in source
    assert "rational_surface_carrier import" not in source


def test_audit_report_script_and_archive_are_frozen() -> None:
    observed = {
        relative: hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        for relative in AUDIT_PINS
    }
    assert observed == AUDIT_PINS


def test_stable_lattice_sum_of_translates_is_exact() -> None:
    replay = stable_lattice_replay()
    assert replay["all_checks_pass"]
    assert replay["stable_action"] == [["0", "1"], ["1", "0"]]
    assert replay["cleared_stable_gram"] == [["0", "1"], ["1", "0"]]
    assert all(replay["checks"].values())


def test_four_square_embedding_reconstructs_4R_minus_3() -> None:
    replay = four_square_bound_replay()
    assert replay["all_checks_pass"]
    assert replay["source_diagonal"] == ["7/3", "-5/6", "-11/10"]
    assert replay["common_scale"] == 30
    assert replay["blowup_count"] == 9
    assert replay["general_count"] == "1+4(R-1)=4R-3"
    assert all(replay["checks"].values())


def test_noncommuting_chow_pair_fixes_composition_orientation() -> None:
    replay = chow_orientation_replay()
    checks = replay["checks"]
    assert replay["all_checks_pass"]
    assert replay["composition_rule"] == "C_A J C_B = C_(BA)"
    assert replay["action_rule"] == "C_A^t J = A"
    assert replay["homomorphism_rule"] == "Gamma_g o Gamma_h = Gamma_(gh)"
    assert checks["noncommuting_pair_really_noncommutes"]
    assert checks["middle_contraction_is_first_left_then_right"]
    assert checks["middle_contraction_is_not_wrong_order"]
    assert checks["reflection_cycles_self_transpose"]
    assert checks["reflection_cycles_square_to_projector"]
    assert checks["m2_relation"] and checks["m3_relation"]
    assert checks["parabolic_nilpotency_index_three"]


def test_quadratic_field_eigenslice_returns_exact_oriented_walls() -> None:
    replay = eigenslice_replay()
    checks = replay["checks"]
    assert replay["all_checks_pass"]
    assert replay["field"] == "Q(sqrt(2))"
    assert replay["trace_rank"] == 4
    assert replay["trace_signature"] == [1, 3, 0]
    assert replay["distinguished_signature"] == [1, 1]
    assert replay["conjugate_signature"] == [0, 2]
    assert replay["wall_pairings"] == ["1", "2 + 2*sqrt(2)"]
    assert checks["field_action_is_trace_self_adjoint"]
    assert checks["field_action_splits_by_real_embeddings"]
    assert checks["native_reflection_is_trace_isometry"]
    assert checks["native_reflection_is_involution"]
    assert checks["native_reflection_commutes_with_field_action"]
    assert checks["trace_pairing_equals_distinguished_pairing"]
    assert checks["oriented_wall_pairings_positive"]


def test_imported_2997_bound_and_leverage_no_go_are_separate() -> None:
    replay = leverage_replay()
    assert replay["all_checks_pass"]
    assert replay["only_internal_rank_statement"] == (
        "d(n+1)=R <= rho(Bl_b(P2))=b+1"
    )
    assert replay["circularity_verdict"] == (
        "2997 follows only after importing n+1<=30 and [k:Q]<=25"
    )
    assert replay["checks"]["quadratic_field_relations_hold_on_fixed_rank_two_tate_part"]
    assert replay["checks"]["imported_rank_bound"]
    assert replay["checks"]["imported_blowup_bound"]


def test_complete_audit_scope_is_pass_without_strong_category_overclaim() -> None:
    payload = build_audit()
    assert payload["verdict"] == "PASS"
    assert payload["all_checks_pass"]
    assert all(payload["checks"].values())
    assert payload["independence"] == {
        "imports_owner_certifier": False,
        "imports_owner_embedding_implementation": False,
        "noncommuting_pair_makes_composition_order_observable": True,
        "quadratic_field_eigenslice_is_separate_from_owner_triangle": True,
    }
    assert "effective/nef/Mori" in payload["scope_verdict"]["not_proved"]
    assert "AG-derived ABSW" in payload["scope_verdict"]["not_proved"]


def test_report_records_proof_pressure_points_and_has_no_controls() -> None:
    text = AUDIT_REPORT.read_text(encoding="utf-8")
    for needle in (
        "PASS, with the owner's marked weak-category and imported-bound scope",
        "The stable trace lattice is legitimate",
        "The \\(4R-3\\) embedding count is exact",
        "Literal Chow composition, transpose, and action orientation",
        "C_AJC_B=C_{BA}",
        "There is no orientation defect",
        "The distinguished eigenslice returns the chamber exactly",
        "The number \\(2997\\) is imported, not derived",
        "The leverage no-go survives hostile pressure",
        "No owner byte was edited",
    ):
        assert needle in text
    assert text.count("\\[") == text.count("\\]")
    raw = AUDIT_REPORT.read_bytes()
    forbidden = set(range(0, 9)) | {11, 12} | set(range(14, 32)) | {127}
    assert not any(byte in forbidden for byte in raw)
    for relative in (
        "scripts/audit_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile.py",
        "data/enumerations/universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.json",
        "tests/test_universal_weak_chow_action_carrier_and_leverage_no_go_independent_hostile_audit.py",
    ):
        assert (ROOT / relative).exists()


def test_archive_is_exact_replay() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert archived["payload"] == build_audit()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert OUTPUT.read_text(encoding="utf-8") == encoded_audit()


def test_cli_check() -> None:
    subprocess.run(
        [sys.executable, str(AUDITOR), "--check"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
