"""Tests for the semiampleness-to-ABSW weak-log composition theorem."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from scripts.certify_pfirst_weak_log_complete_carrier_transfer import (
    build_certificate,
    verify_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/pfirst_weak_log_complete_carrier_transfer.md"
CERTIFIER = ROOT / "scripts/certify_pfirst_weak_log_complete_carrier_transfer.py"
ARCHIVE = ROOT / "data/enumerations/pfirst_weak_log_complete_carrier_transfer.json"
PINS = {
    REPORT: "bff338b5e2f4dde2c918742a853656661c2ac074684938031ebbc95df8fd14b2",
    CERTIFIER: "1be1ceaf48fcd4cd9d719a552f2bb9f0c37b80dd3531c3a94767cb05d5786802",
    ARCHIVE: "b8426d30d2e5cdfa8224b76567c004926b77d58e53a0e4dbdbbaadfdb4b67e74",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def archived() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def redigest(payload: dict) -> None:
    payload.pop("content_sha256_without_this_field", None)
    payload["content_sha256_without_this_field"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def test_frozen_bytes_and_exact_archive_replay() -> None:
    assert {path: sha256(path) for path in PINS} == PINS
    payload = archived()
    assert payload == build_certificate()
    assert verify_certificate(payload)


def test_rank_26_elimination_is_a_global_linear_combination() -> None:
    row = archived()["rank_elimination"]
    assert row["nonnegative_linear_combination"] == (
        "7*(10-rho+b)+(4*rho+8-7*b)=78-3*rho>=0"
    )
    assert row["maximum_picard_rank"] == 26
    assert row["maximum_hyperbolic_dimension"] == 25
    assert row["empty_boundary_maximum_picard_rank"] == 10
    for rho in range(0, 200):
        feasible_b = [
            b
            for b in range(0, 200)
            if 10 - rho + b >= 0 and 4 * rho + 8 - 7 * b >= 0
        ]
        assert bool(feasible_b) is (rho <= 26)


def test_composed_f4_packet_is_nonzero_and_uses_the_complete_provider() -> None:
    row = archived()["f4_composed_calibration"]
    assert row["weak_log_class"] == [[3, 2], [6, 1], [-1, 1]]
    assert row["twice_weak_log_class"] == [3, 12, -2]
    assert row["weak_log_square"] == [8, 1]
    assert row["affine_chart_count"] == 5
    assert row["global_section_count"] == 4
    assert row["coxeter_orders"] == {
        "B_half_E": 2,
        "B_half_G": 3,
        "E_G": "infinity",
    }
    assert row["area_over_pi"] == [1, 6]
    assert row["exact_full_nef_return"] is True


def test_scope_rejects_totality_or_circular_nef_promotions() -> None:
    for key in (
        "composition_producer_accepts_arbitrary_carrier_packet",
        "arbitrary_PB_Atlas_adapter_implemented",
        "raw_surface_free_essential_surjectivity_proved",
        "nonrational_eigen_nef_branch_exhausted",
        "global_ABSW_replacement_proved",
    ):
        mutated = deepcopy(archived())
        mutated["scope"][key] = True
        redigest(mutated)
        with pytest.raises(ValueError, match="promoted"):
            verify_certificate(mutated)

    for key in (
        "complete_Mori_generator_list_assumed",
        "nefness_inferred_from_proposed_walls",
    ):
        mutated = deepcopy(archived())
        mutated["scope"][key] = True
        redigest(mutated)
        with pytest.raises(ValueError):
            verify_certificate(mutated)


def test_false_geometric_gate_and_stale_source_are_rejected() -> None:
    mutated = deepcopy(archived())
    mutated["certificate_gates"]["exact_effective_irreducible_oriented_facets"] = False
    redigest(mutated)
    with pytest.raises(ValueError, match="gate"):
        verify_certificate(mutated)

    mutated = deepcopy(archived())
    mutated["sources"]["semiampleness_theorem"]["sha256"] = "0" * 64
    redigest(mutated)
    with pytest.raises(ArithmeticError, match="stale source"):
        verify_certificate(mutated)


def test_report_states_exact_payoff_and_universal_debt() -> None:
    report = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "does not assert that every raw maximal arithmetic chamber",
        "does not infer nefness by testing",
        r"\boxed{P=\mathbf P(\operatorname {Nef}S)}",
        r"\boxed{n\leq25}",
        r"\boxed{k(P)=\mathbf Q}",
        "does not assume the conclusion it is meant to return",
        "arbitrary PB-Atlas B1/B2a-to-cover adapter",
        "not an implementation",
        "does **not** ingest and validate an arbitrary complete carrier packet",
        "does not require them to be listed or to generate the Mori cone",
    ):
        assert phrase in report

    scope = archived()["scope"]
    assert scope["composition_producer_is_dependency_and_calibration_replay"] is True
    assert scope["composition_producer_accepts_arbitrary_carrier_packet"] is False


def test_cli_checks_frozen_archive() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
    )
