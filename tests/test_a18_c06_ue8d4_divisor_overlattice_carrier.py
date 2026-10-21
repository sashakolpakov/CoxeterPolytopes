from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_a18_c06_ue8d4_divisor_overlattice_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_unique_even_index_four_overclosure(certificate: dict) -> None:
    lattice = certificate["lattice_and_embedding"]
    census = lattice["stable_even_overclosure_census"]
    overclosure = lattice["overclosure_N"]
    assert census["stable_integral_orders"] == [1, 4]
    assert census["stable_even_orders"] == [1, 4]
    assert census["unique_nontrivial_even_index_four"] is True
    assert overclosure["isometry_class"] == "U+E8(-1)+D4(-1)"
    assert overclosure["determinant"] == -4
    assert overclosure["smith"] == [1] * 12 + [2, 2]
    assert all(lattice["checks"].values())


def test_explicit_embedding_and_all_reflections(certificate: dict) -> None:
    lattice = certificate["lattice_and_embedding"]
    inclusion = lattice["explicit_M_into_family_NS"]
    mirrors = certificate["primitive_mirrors"]
    assert inclusion["index"] == 4
    assert inclusion["identity"] == "I^T Gram(NS) I = Gram(M)"
    assert mirrors["count"] == 31
    assert mirrors["M_presentation_square_histogram"] == {
        "-2": 14,
        "-4": 5,
        "-8": 12,
    }
    assert mirrors["NS_coordinate_content_histogram"] == {"1": 19, "2": 12}
    assert mirrors["primitive_NS_square_histogram"] == {"-2": 26, "-4": 5}
    assert mirrors["primitive_NS_square_divisibility_histogram"] == {
        "-2,1": 26,
        "-4,2": 5,
    }
    assert mirrors["all_31_reflections_extend_integrally"] is True
    assert all(mirrors["normalization_checks"].values())
    assert all(all(row["checks"].values()) for row in mirrors["records"])


def test_actual_weierstrass_family_and_effective_roots(certificate: dict) -> None:
    family = certificate["explicit_K3_family"]
    effectivity = certificate["effectivity_and_scope"]
    assert family["very_general_NS"] == "U+E8(-1)+D4(-1)"
    assert "t^3*(t-1)*(t-a)*(t-b)*(t-c)" in family["equation"]
    assert all(family["checks"].values())
    assert effectivity["minimum_oriented_pairing"] == "1"
    assert all(effectivity["checks"].values())
    assert "26 primitive square -2" in effectivity["minus_two_effectivity"]
    assert "five primitive square -4" in effectivity["mixed_warning"]


def test_literal_vinberg_gate_and_scope(certificate: dict) -> None:
    vinberg = certificate["vinberg_finite_volume"]
    scope = certificate["scope"]
    assert vinberg["classification"] == "finite_volume_noncompact_Coxeter_polytope"
    assert (vinberg["dimension"], vinberg["facets"]) == (13, 31)
    assert (vinberg["finite_vertices"], vinberg["cusps"]) == (3433, 238)
    assert vinberg["elliptic_edges"] == 28766
    assert vinberg["connected_affine_components"] == 529
    assert vinberg["bad_edges"] == vinberg["bad_parabolic_neighborhoods"] == 0
    assert all(vinberg["checks"].values())
    assert "refine" in scope["not_full_nef"]
    assert "no claim" in scope["not_automorphism_fixed"]
    assert "no new abstract-polytope" in scope["novelty"]
    assert certificate["all_checks_pass"] is True


def test_report_freezes_positive_gain_and_boundary() -> None:
    report_path = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "a18_c06_ue8d4_divisor_overlattice_carrier.md"
    )
    report = " ".join(report_path.read_text(encoding="utf-8").lower().split())
    assert "equation-level k3 divisor carrier" in report
    assert "26 of square `-2`" in report
    assert "5 of square `-4`" in report
    assert "not the full nef cone" in report
    assert "not a new abstract polytope" in report
