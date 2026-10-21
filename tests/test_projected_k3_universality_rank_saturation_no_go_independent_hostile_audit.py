from __future__ import annotations

import hashlib
import json
from pathlib import Path

from scripts.audit_projected_k3_universality_rank_saturation_no_go_independent_hostile import (
    OUTPUT,
    OWNER_PINS,
    ROOT,
    build_audit,
)


REPORT = (
    ROOT
    / "results/projected_k3_universality_rank_saturation_no_go_independent_hostile_audit.md"
)
SCRIPT = (
    ROOT
    / "scripts/audit_projected_k3_universality_rank_saturation_no_go_independent_hostile.py"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "0235d2d46b8d7badc55b87e7f74532a80fb227f48085a1538152729fc87a9c80"
)
EXPECTED_AUDIT_FILE_SHA256 = {
    "report": "89296f04bac6b0b4ed68d9b86936a08ede515ba5bf9c847abeaadbc628926071",
    "script": "2e9b68a1049a1e796453972af9b4c679e298007e507defdde0731aac2d545595",
    "archive": "0ce5e3675079c3188cbd3575860b9f8130e51b7a362783f8d2df90a7f8a6337a",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_independent_audit_passes_and_is_frozen() -> None:
    payload = build_audit()
    assert payload["all_checks_pass"]
    assert payload["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256


def test_archive_is_exact_replay() -> None:
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == build_audit()


def test_terminal_owner_pins_are_checked() -> None:
    for relative, expected in OWNER_PINS.items():
        assert _sha256(ROOT / relative) == expected


def test_audit_artifacts_are_pinned() -> None:
    assert _sha256(REPORT) == EXPECTED_AUDIT_FILE_SHA256["report"]
    assert _sha256(SCRIPT) == EXPECTED_AUDIT_FILE_SHA256["script"]
    assert _sha256(OUTPUT) == EXPECTED_AUDIT_FILE_SHA256["archive"]


def test_x3_rank_saturation_and_square_class_are_independent() -> None:
    replay = build_audit()["x3_independent_replay"]
    assert replay["hyperbolic_dimension"] == 19
    assert replay["raw_rank"] == 20
    assert replay["inertia"] == [1, 19, 16]
    assert replay["smith_nonzero"] == [1] * 19 + [3]
    assert replay["facet_norm_counts"] == {"-2": 24, "-6": 12}
    assert replay["complex_k3_picard_ceiling"] == 20
    assert replay["forced_projector_complement_rank"] == 0
    assert replay["square_class_ratio"] == "3"
    assert not replay["square_class_ratio_is_rational_square"]
    assert replay["three_adic_valuation"] == 1


def test_separate_block_checks_projector_shell_and_glue_calibration() -> None:
    block = build_audit()["independent_block_calibration"]
    assert block["all_checks_pass"]
    assert block["ambient_signature"] == [1, 4, 0]
    assert block["schur_signature"] == [0, 3, 0]
    assert block["scale"] > block["row_shell_bound"]
    assert block["facet_rank_mod_prime"] == 2
    assert block["complement_rank"] >= block["facet_rank_mod_prime"]
    assert all(block["checks"].values())


def test_all_five_hostile_logic_pressure_points_are_recorded() -> None:
    logic = build_audit()["logic_audit"]
    assert set(logic) == {
        "generic_facet_root_necessity",
        "exact_module_similarity",
        "compact_nonrational_shell",
        "rational_cusp_completeness",
        "compression_discriminant_bound",
    }
    assert "projection is the identity" in logic["generic_facet_root_necessity"]
    assert "one rational similarity" in logic["exact_module_similarity"]
    assert "Minkowski lattice" in logic["compact_nonrational_shell"]
    assert "integral reflection vector" in logic["rational_cusp_completeness"]
    assert "injects into A_K" in logic["compression_discriminant_bound"]


def test_scope_is_all_finite_volume_not_maximal_only() -> None:
    scope = build_audit()["scope_verdict"]
    assert scope["status"] == "PASS"
    assert "all finite-volume arithmetic chambers" in scope["proved"]
    assert "absolutely maximal" in scope["not_proved"]
    assert "mixed-divisor" in scope["not_proved"]


def test_report_preserves_owner_scope_and_no_repair_verdict() -> None:
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for needle in (
        "PASS, with the owner's stated exact-module and non-maximal-only scope",
        "generic positive point",
        "common-\\(\\mathbf Q\\)-similarity",
        "Compact nonrational residue shells",
        "rational cusped host-completeness argument",
        "Facet compression and primitive glue",
        "does not address a theorem whose inputs are required to be absolutely",
        "No scope repair to the owner bytes is required",
    ):
        assert needle in text
