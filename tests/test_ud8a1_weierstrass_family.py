from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_ud8a1_weierstrass_family import OUTPUT, build_certificate


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_explicit_open_family_and_witness(certificate: dict) -> None:
    family = certificate["weierstrass_family"]
    assert family["parameters"] == ["r", "a", "b", "c", "d", "e", "f", "g", "h"]
    assert family["fiber_configuration"] == "I4*+I2+12I1"
    assert family["nonempty_witness"]["parameters"] == [1, 0, 0, 0, 1, 0, 0, 0, 0]
    assert family["nonempty_witness"]["residual_polynomial_discriminant"] != 0
    assert family["nonempty_witness"]["resultant_with_P"] != 0
    assert all(family["checks"].values())


def test_very_general_picard_and_mordell_weil_gate(certificate: dict) -> None:
    lattice = certificate["very_general_picard_lattice"]
    assert lattice["isometry_class"] == "U+D8(-1)+A1(-1)"
    assert lattice["rank_signature_determinant"] == [11, [1, 10, 0], 8]
    assert lattice["smith"] == [1] * 8 + [2] * 3
    assert lattice["torsion_saturation"]["maximum_total_correction"] == "5/2"
    assert lattice["torsion_saturation"]["nonzero_section_height_lower_bound"] == "3/2"
    assert all(lattice["checks"].values())


def test_unimodular_bridge_to_published_curve_basis(certificate: dict) -> None:
    bridge = certificate["explicit_fiber_basis_bridge"]
    assert bridge["determinant_T"] == -1
    assert bridge["standard_basis"] == [
        "F=A9+A10",
        "F+A11",
        "A1",
        "A2",
        "A3",
        "A4",
        "A5",
        "A6",
        "A7",
        "A13",
        "A9",
    ]
    assert all(bridge["checks"].values())


def test_scope_keeps_classical_family_and_new_carrier_distinct(
    certificate: dict,
) -> None:
    scope = certificate["scope"]
    handoff = certificate["future_A18_C04_H10_28_handoff"]
    assert "Roulleau already published" in scope["classical_part"]
    assert "equation-level K3 family" in scope["positive_gain"]
    assert scope["no_automorphism_claim"] is True
    assert scope["no_H10_28_claim_before_composition"] is True
    assert handoff["status"] == "equation_target_ready_carrier_composition_pending"
    assert certificate["all_checks_pass"] is True


def test_report_freezes_equation_and_lineage() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "ud8a1_weierstrass_family.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "nine-parameter" in normalized
    assert "i4* + i2 + 12i1" in normalized
    assert "u+d8(-1)+a1(-1)" in normalized
    assert "height" in normalized
    assert "roulleau" in normalized
    assert "not yet" in normalized
