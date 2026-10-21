from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_ghess_h10_27_full_picard_k3_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_fixed_lattice_and_discriminant_form(certificate: dict) -> None:
    lattice = certificate["full_Picard_lattice"]
    discriminant = lattice["discriminant_form"]
    split = lattice["standard_identification"]
    assert lattice["rank"] == 11
    assert lattice["signature"] == [1, 10, 0]
    assert lattice["determinant"] == 256
    assert lattice["smith_invariants"] == [1] * 4 + [2] * 6 + [4]
    assert discriminant["group"] == "A_S=(Z/4)+(Z/2)^6"
    assert discriminant["quadratic_form"] == (
        "q_S=<3/4> on Z/4 + u(2)^2 + v(2)"
    )
    assert discriminant["generator_orders"] == [4, 2, 2, 2, 2, 2, 2]
    assert all(discriminant["checks"].values())
    assert split["identification"].startswith("S_H=U+K_H")
    assert all(split["checks"].values())


def test_explicit_nonstable_primitive_k3_embedding(certificate: dict) -> None:
    embedding = certificate["full_Picard_lattice"]["primitive_K3_embedding"]
    period = certificate["full_Picard_lattice"]["period_consequence"]
    assert embedding["criterion"].startswith(
        "explicit full-discriminant graph gluing"
    )
    assert embedding["orthogonal_transcendental_lattice"] == (
        "T_H=U+U(2)+A2(-2)"
    )
    assert embedding["generator_orders"] == [6, 2, 2, 2]
    assert embedding["ambient_identification"].endswith("Lambda_K3")
    assert embedding["fixed_orthogonal_complement_signature"] == [2, 9, 0]
    assert all(embedding["checks"].values())
    assert period["K3_moduli_dimension"] == 9
    assert period["very_general_full_Picard"] == "NS(X)=S_H"


def test_all_27_walls_generate_full_picard(certificate: dict) -> None:
    chamber = certificate["transported_chamber"]
    assert chamber["id"] == "GHESS-C2-H10-27"
    assert chamber["known_abstract_id"] == "D16-D4-C7-H10-27"
    assert chamber["hyperbolic_dimension"] == 10
    assert chamber["facet_count"] == len(chamber["wall_records"]) == 27
    assert chamber["wall_span_rank"] == 11
    assert chamber["wall_span_smith"] == [1] * 11
    assert chamber["wall_span_index_in_full_Picard"] == 1
    assert chamber["primitive_mirror_type_histogram"] == {
        "-2,1": 8,
        "-4,2": 10,
        "-4,4": 9,
    }
    assert all(
        all(record["checks"].values())
        for record in chamber["wall_records"]
    )


def test_strict_orientation(certificate: dict) -> None:
    chamber = certificate["transported_chamber"]
    assert chamber["strict_interior_square"] == 20
    assert chamber["minimum_wall_pairing"] == 1
    assert chamber["strict_wall_pairings"] == (
        [1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1]
        + [4, 8, 8, 8, 4, 4, 4, 16, 16, 16, 16, 16, 16]
    )


def test_literal_vinberg_and_exact_cone(certificate: dict) -> None:
    finite = certificate["finite_volume"]
    vinberg = finite["literal_Vinberg"]
    cone = finite["exact_GMP_cone"]
    assert finite["classification"] == (
        "finite_volume_noncompact_Coxeter_polytope"
    )
    assert (
        vinberg["finite_vertex_count"],
        vinberg["ideal_vertex_count"],
        vinberg["elliptic_edge_count"],
    ) == (895, 67, 6221)
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighborhood_count"] == 0
    assert (
        cone["extreme_ray_count"],
        cone["finite_vertex_count"],
        cone["ideal_vertex_count"],
        cone["hyperideal_extreme_ray_count"],
    ) == (962, 895, 67, 0)
    assert cone["active_facet_count_histogram"] == {
        "finite_10": 895,
        "ideal_11": 6,
        "ideal_12": 58,
        "ideal_13": 3,
    }
    assert all(finite["checks"].values())


def test_lineage_and_scope(certificate: dict) -> None:
    lineage = certificate["lineage"]
    scope = certificate["scope"]
    assert lineage["normalized_diagram"] == "D16-D4-C7-H10-27"
    assert lineage["known_reflection_subgroup_index"] == 2016
    assert lineage["new_abstract_polytope"] is False
    assert lineage["new_arithmetic_class"] is False
    assert scope["not_equation_level"] is True
    assert "not K3 negative-curve walls" in scope["not_full_nef"]
    assert "not claimed" in scope["not_surface_automorphism"]
    assert scope["worldwide_priority_claim"] == "none"
    assert all(certificate["checks"].values())
    assert certificate["all_checks_pass"] is True


def test_report_freezes_positive_and_caveats() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results/ghess_h10_27_full_picard_k3_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "toot" in normalized
    assert "nine-dimensional" in normalized
    assert "895" in normalized and "67" in normalized and "6,221" in normalized
    assert "not the full nef cone" in normalized
    assert "not a surface-automorphism" in normalized
    assert "not a new abstract polytope" in normalized
