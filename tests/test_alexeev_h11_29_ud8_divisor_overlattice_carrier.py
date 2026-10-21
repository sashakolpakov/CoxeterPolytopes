from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_alexeev_h11_29_ud8_divisor_overlattice_carrier import (
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
    census = certificate["lattice_and_embedding"][
        "reflection_stable_even_overclosure_census"
    ]
    overclosure = certificate["lattice_and_embedding"]["overclosure_N"]
    assert census["stable_even_orders"] == [1, 4]
    assert census["nontrivial_even_index_four_count"] == 1
    assert census["unique"] is True
    assert overclosure["isometry_class"] == "U+D8(-1)+A1(-1)^2"
    assert overclosure["determinant"] == -16
    assert all(overclosure["checks"].values())


def test_explicit_embedding_and_all_reflections(certificate: dict) -> None:
    lattice = certificate["lattice_and_embedding"]
    inclusion = lattice["explicit_M_into_family_NS"]
    walls = lattice["facet_reflection_extensions"]
    assert inclusion["index"] == 4
    assert inclusion["identity"] == "I^T Gram(NS) I = Gram(M)"
    assert walls["count"] == 29
    assert walls["M_presentation_square_histogram"] == {
        "-2": 11,
        "-4": 12,
        "-8": 6,
    }
    assert walls["NS_coordinate_content_histogram"] == {"1": 23, "2": 6}
    assert walls["primitive_NS_square_histogram"] == {"-2": 17, "-4": 12}
    assert walls["primitive_NS_square_divisibility_histogram"] == {
        "-2,1": 16,
        "-2,2": 1,
        "-4,2": 12,
    }
    assert all(walls["normalization_checks"].values())
    assert walls["all_29_extend_integrally"] is True
    assert all(all(row["checks"].values()) for row in walls["records"])
    assert all(
        row["checks"]["presentation_and_primitive_reflections_equal"]
        and row["checks"]["same_hyperplane"]
        for row in walls["records"]
    )


def test_actual_weierstrass_and_scoped_action(certificate: dict) -> None:
    family = certificate["explicit_K3_family"]
    action = certificate["family_action_on_divisor_sublattice"]
    assert family["BH_id"] == "0.4.4.2"
    assert family["very_general_NS"] == "U+D8(-1)+A1(-1)^2"
    assert all(family["checks"].values())
    assert action["checks"]["M_is_stable_under_family_Picard_involution"] is True
    assert action["checks"]["selected_oriented_simple_wall_set_not_preserved"] is True
    assert action["selected_oriented_simple_walls_preserved"] == 14
    assert "chamber is not preserved" in action["conclusion"]


def test_bounded_primary_source_lineage(certificate: dict) -> None:
    audit = certificate["bounded_primary_source_lineage_audit"]
    assert "118-root" in audit["alexeev"]["records"]
    assert audit["brandhorst_hofmann"]["records"] == "classification row 0.4.4.2"
    assert "19 curve walls" in audit["roulleau"]["records"]
    assert "appears unrecorded" in audit["conclusion"]
    assert audit["scope"] == "bounded primary-source evidence, not worldwide novelty"


def test_literal_vinberg_gate_and_scope(certificate: dict) -> None:
    vinberg = certificate["vinberg_finite_volume"]
    scope = certificate["scope"]
    effectivity = certificate["effectivity_and_divisor_interpretation"]
    assert vinberg["classification"] == "finite_volume_noncompact_Coxeter_polytope"
    assert (vinberg["dimension"], vinberg["facets"]) == (11, 29)
    assert (vinberg["finite_vertices"], vinberg["cusps"]) == (3876, 253)
    assert vinberg["elliptic_edges"] == 28798
    assert vinberg["bad_edges"] == vinberg["bad_parabolic_neighborhoods"] == 0
    assert all(vinberg["checks"].values())
    assert "refine" in scope["not_full_nef"]
    assert "not the selected 29-simple-wall set" in scope[
        "not_automorphism_fixed_section"
    ]
    assert "seventeen square -2" in effectivity["minus_two_wall_effectivity"]
    assert "twelve square -4" in effectivity["mixed_mirror_warning"]
    assert certificate["all_checks_pass"] is True


def test_report_freezes_primitive_normalization_and_scope() -> None:
    report_path = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "alexeev_h11_29_ud8_divisor_overlattice_carrier.md"
    )
    report = " ".join(report_path.read_text(encoding="utf-8").lower().split())
    assert "for a very general member" in report
    assert "17 of square -2" in report
    assert "0 of square -8" in report
    assert "not an automorphism-fixed-section claim" in report
    assert "not the full nef cone" in report
