"""Tests for the unconditional NS/Chow and maximal-chamber synthesis."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

from scripts.certify_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness import (
    ARCHIVE,
    BLOWUP_MAX,
    REPORT,
    TRACE_RANK_MAX,
    build_certificate,
    digest,
    encoded_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFIER = ROOT / "scripts/certify_unconditional_ns_chow_encoding_and_maximal_chamber_finiteness.py"


def test_unconditional_weak_carrier_scope() -> None:
    carrier = build_certificate()["unconditional_weak_carrier"]
    assert "admissible Vinberg lattice-form" in carrier["input"]
    assert carrier["requires_reflection_maximality"] is False
    assert carrier["full_trace_NS_injection"] is True
    assert carrier["faithful_rational_Chow_reflection_action"] is True
    assert carrier["exact_semialgebraic_eigenslice_return"] is True
    assert carrier["effective_or_nef_return"] is False
    assert TRACE_RANK_MAX == 750
    assert BLOWUP_MAX == 2997


def test_unconditional_maximal_chamber_scope() -> None:
    theorem = build_certificate()["unconditional_maximal_chambers"]
    assert theorem["finite_native_pair_set"] is True
    assert theorem["finite_native_similarity_hosts"] is True
    assert theorem["finite_unmarked_normalized_chambers"] is True
    assert theorem["ABSW_free_factorization"] is True
    assert theorem["AG_carrier_hypothesis"] is False


def test_provenance_does_not_turn_weak_carrier_into_an_ag_proof() -> None:
    scope = build_certificate()["provenance_boundary"]
    assert scope["finiteness_derived_from_weak_carrier"] is False
    assert scope["classical_dimension_and_field_inputs_imported"] is True
    assert scope["global_AG_only_finiteness_claimed"] is False
    assert scope["universal_effective_nef_or_mixed_action_compiler_proved"] is False
    assert scope["current_complex_high_leverage_stack_universal"] is False


def test_report_answers_the_attached_question_literally() -> None:
    text = REPORT.read_text(encoding="utf-8")
    assert "## 4. Exact answer to the attached conditional statement" in text
    assert r"most \(2997\)" in text
    assert "finite set of unmarked normalized classes" in text
    assert "unconditional **AG-only**" in text
    assert "standard admissible Vinberg lattice-form model" in text
    assert "raw ambient trace" in text


def test_archive_is_exact_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_cli_check_output() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check-output"],
        cwd=ROOT,
        check=True,
    )
