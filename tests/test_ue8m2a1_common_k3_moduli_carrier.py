from __future__ import annotations

import json
from pathlib import Path

import pytest
import sympy as sp

from scripts.certify_ue8m2a1_common_k3_moduli_carrier import (
    F_GRAM,
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_F_and_exhaustive_reflection_alphabet(certificate: dict) -> None:
    lattice = certificate["lattice"]
    assert lattice["rank"] == 11
    assert lattice["signature"] == [1, 10, 0]
    assert lattice["determinant"] == 512
    assert lattice["smith"] == [1, 1] + [2] * 9
    assert lattice["discriminant_exponent"] == 2
    assert lattice["primitive_integral_reflection_alphabet"] == [
        {"square": -2, "divisibility": 1},
        {"square": -2, "divisibility": 2},
        {"square": -4, "divisibility": 2},
    ]
    witnesses = lattice["explicit_negative_type_witnesses"]
    assert {
        (row["square"], row["divisibility"])
        for row in witnesses.values()
    } == {(-2, 1), (-2, 2), (-4, 2)}
    assert all(
        row["primitive"]
        and row["dual_class_exact_order"] == row["divisibility"]
        and row["reflection_integral"]
        and row["reflection_preserves_F"]
        and row["hyperbolic_mirror_normal"]
        for row in witnesses.values()
    )


def test_positive_norm_integral_reflection_is_not_a_wall_type(
    certificate: dict,
) -> None:
    """The U-vector e+f is the adversarial scope counterexample."""

    boundary = certificate["lattice"][
        "positive_norm_integral_reflection_boundary_witness"
    ]
    assert boundary["coordinates"] == [1, 1] + [0] * 9
    assert boundary["square"] == 2
    assert boundary["divisibility"] == 1
    assert boundary["primitive"] is True
    assert boundary["reflection_integral"] is True
    assert boundary["reflection_preserves_F"] is True
    assert boundary["hyperbolic_mirror_normal"] is False

    vector = sp.Matrix(boundary["coordinates"])
    reflection = sp.eye(11) - vector * (vector.T * F_GRAM)
    assert reflection.T * F_GRAM * reflection == F_GRAM


def test_equation_family_is_full_nine_dimensional_carrier(
    certificate: dict,
) -> None:
    family = certificate["equation_family"]
    assert family["parameter_count"] == family["moduli_dimension"] == 9
    assert family["very_general_full_Picard"] == (
        "NS(X)=U+E8(-2)+<-2> for a very general parameter"
    )
    assert family["nonempty_open_witness"]
    assert "A(s^2)" in family["equation"]
    assert "ambient K3 surface" in family["weak_carrier_scope"]
    assert "does not provide equations" in family["weak_carrier_scope"]


def test_three_calibrations_share_F_and_pass_volume_gates(
    certificate: dict,
) -> None:
    rows = certificate["calibrations"]
    assert list(rows) == [
        "A15_11-C2-H10-21",
        "K3-S3-C2-H10-13",
        "X2-RAT-C2-H10-37",
    ] or set(rows) == {
        "K3-S3-C2-H10-13",
        "A15_11-C2-H10-21",
        "X2-RAT-C2-H10-37",
    }
    by_facets = {row["facets"]: row for row in rows.values()}
    assert set(by_facets) == {13, 21, 37}
    assert [by_facets[count]["wall_span_index"] for count in [13, 21, 37]] == [
        1,
        1,
        1,
    ]
    assert [by_facets[count]["cusps"] for count in [13, 21, 37]] == [
        4,
        50,
        261,
    ]
    assert [
        by_facets[count]["elliptic_rank9_edges"] for count in [13, 21, 37]
    ] == [208, 4159, 24350]
    assert [
        by_facets[count]["index_in_full_reflection_group"]
        for count in [13, 21, 37]
    ] == [1, 6720, 241920]
    assert all(row["all_checks_pass"] for row in rows.values())


def test_full_group_criterion_and_compact_no_go(certificate: dict) -> None:
    group = certificate["full_reflection_group"]
    assert group["fundamental_chamber"] == "K3-S3-C2-H10-13"
    assert group["fundamental_facets"] == 13
    assert group["cusps"] == 4
    assert group["nonuniform"] is True
    proof = group["proof_of_subgroup_criterion"]
    assert "Coxeter-Poincare" in proof
    assert "disjoint union" in proof
    assert "finite iff the index is finite" in proof
    assert "powers of parabolic" in proof
    theorem = certificate["theorem"]
    assert "primitive negative" in theorem["common_equation_family_criterion"]
    assert "Positive-norm" in theorem["converse_for_this_marked_family"]
    assert "fundamental domain" in theorem["finite_volume_subgroup_criterion"]
    assert "primitive negative" in theorem["compact_no_go_inside_F"]


def test_promotion_and_scope_are_not_overclaimed(certificate: dict) -> None:
    boundary = certificate["promotion_boundary"]
    assert "separate complete Vinberg" in boundary["volume_gate"]
    assert "finite index" in boundary["volume_gate"]
    assert "may be reducible" in boundary["reducibility"]
    assert boundary["full_nef_criterion"].startswith(
        "After the chamber-specific Weyl re-marking"
    )
    assert "cannot be a full K3 nef wall" in boundary["minus_four_obstruction"]
    scope = certificate["scope"]
    assert scope["new_abstract_polytopes"] == "none"
    assert scope["surface_automorphisms"] == "not asserted"
    assert "separate Weyl markings" in scope["simultaneous_ample_marking"]
    assert scope["equation_level_wall_representatives"].startswith(
        "not asserted"
    )


def test_all_source_files_and_payloads_are_pinned(certificate: dict) -> None:
    sources = certificate["sources"]
    assert set(sources) == {
        "family",
        "h10_13",
        "h10_21",
        "h10_37",
        "h10_37_lineage",
    }
    for key in ["family", "h10_13", "h10_21", "h10_37"]:
        assert len(sources[key]["file_sha256"]) == 64
        assert len(sources[key]["payload_sha256"]) == 64
    assert len(sources["h10_37_lineage"]["file_sha256"]) == 64
    assert "payload_sha256" not in sources["h10_37_lineage"]


def test_report_states_the_bidirectional_boundary(certificate: dict) -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results/ue8m2a1_common_k3_moduli_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "carrier factory for **every genuine finite-sided coxeter" in normalized
    assert "primitive-negative" in normalized
    assert "exhaustive arithmetic alphabet" in normalized
    assert "e+f" in normalized
    assert "finite volume, finite index" in normalized
    assert "no compact `h10` coxeter polytope" in normalized
    assert "complete vinberg" in normalized
    assert "not a no-go for the abstract polytope" in normalized
    assert certificate["all_checks_pass"] is True
    assert all(certificate["checks"].values())
