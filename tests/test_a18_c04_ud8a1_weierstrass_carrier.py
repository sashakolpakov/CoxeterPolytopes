from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_a18_c04_ud8a1_weierstrass_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_pinned_equation_family_and_exact_composition(certificate: dict) -> None:
    family = certificate["equation_family"]
    composition = certificate["exact_composition"]
    assert family["fiber_configuration"] == "I4*+I2+12I1"
    assert family["very_general_NS"] == "U+D8(-1)+A1(-1)"
    assert family["very_general_MW"] == "trivial"
    assert len(composition["nef_reduction_wall_sequence_one_based"]) == 226
    assert composition["determinant_N_to_equation_fiber_basis"] == -1
    assert composition["maximum_absolute_transport_entry"] == 58
    assert all(composition["checks"].values())


def test_all_28_mirrors_and_reflections_are_literal(certificate: dict) -> None:
    mirrors = certificate["divisor_mirrors"]
    assert mirrors["count"] == 28
    assert mirrors["type_histogram_square_divisibility"] == {
        "-4,2": 4,
        "-2,1": 16,
        "-2,2": 8,
    }
    assert mirrors["maximum_absolute_fiber_basis_coordinate"] == 6
    assert mirrors["generation_index"] == 1
    assert len(mirrors["oriented_28_by_28_Gram"]) == 28
    assert len(mirrors["records"]) == 28
    assert all(
        len(row["reflection_matrix_in_equation_fiber_basis"]) == 11
        and all(row["checks"].values())
        for row in mirrors["records"]
    )


def test_explicit_ample_divisor_lies_inside_both_chambers(
    certificate: dict,
) -> None:
    interior = certificate["algebro_geometric_interior"]
    assert interior["transported_chamber_interior_square"] == 508
    assert interior["transported_chamber_minimum_pairing"] == 3
    assert interior["status"].startswith("big_and_nef")
    assert interior["ample_inside_both_chambers_square"] == 1408
    assert interior["minimum_complete_nef_wall_pairing"] == 1
    assert interior["minimum_imported_wall_pairing"] == 7
    assert "nonempty intersection" in interior["conclusion"]
    assert all(interior["checks"].values())


def test_literal_vinberg_gate_and_bounded_scope(certificate: dict) -> None:
    vinberg = certificate["literal_Vinberg_gate"]
    scope = certificate["lineage_and_scope"]
    assert vinberg["classification"] == "finite_volume_noncompact_Coxeter_polytope"
    assert (vinberg["dimension"], vinberg["facet_count"]) == (10, 28)
    assert (vinberg["finite_vertex_count"], vinberg["ideal_vertex_count"]) == (
        1087,
        74,
    )
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighborhood_count"] == 0
    assert all(vinberg["checks"].values())
    assert "equation-to-28-wall" in scope["new_positive_gain"]
    assert scope["not_the_complete_nef_chamber"] is True
    assert scope["no_surface_automorphism_claim"] is True
    assert scope["not_a_new_abstract_polytope"] is True
    assert certificate["all_checks_pass"] is True


def test_report_freezes_positive_theorem_and_boundary() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "a18_c04_ud8a1_weierstrass_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "toot" in normalized
    assert "nine-parameter" in normalized
    assert "28 explicit primitive divisor mirrors" in normalized
    assert "square `1408`" in normalized
    assert "ample cone" in normalized
    assert "not the complete nef chamber" in normalized
    assert "no automorphism" in normalized

