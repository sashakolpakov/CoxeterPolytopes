from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_k3_d4_h12_full_picard_k3_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_saturated_fixed_lattice_and_standard_decomposition(
    certificate: dict,
) -> None:
    handoff = certificate["saturated_fixed_lattice_handoff"]
    lattice = certificate["full_Picard_lattice"]
    assert handoff["saturation_check"] is True
    assert handoff["source_ambient_Picard"] == (
        "U+E8(-1)+D4(-1)+A1(-1)"
    )
    assert lattice["name"] == "S=U+E8(-1)+A2(-1)+A1(-1)"
    assert lattice["ordered_standard_decomposition"] == (
        "U+A2(-1)+A1(-1)+E8(-1)"
    )
    assert lattice["rank"] == 13
    assert lattice["signature"] == [1, 12, 0]
    assert lattice["determinant"] == 6
    assert lattice["smith_invariants"] == [1] * 12 + [6]
    assert all(lattice["checks"].values())


def test_explicit_primitive_k3_embedding_and_generic_period(
    certificate: dict,
) -> None:
    lattice = certificate["full_Picard_lattice"]
    embedding = lattice["primitive_K3_embedding"]
    period = lattice["period_consequence"]
    assert embedding["ambient_lattice"] == "Lambda_K3=U^3+E8(-1)^2"
    assert embedding["orthogonal_complement_signature"] == [2, 7, 0]
    assert len(embedding["visible_unit_minor_rows_zero_based"]) == 13
    assert all(embedding["checks"].values())
    assert period["K3_moduli_dimension"] == 7
    assert period["very_general_full_Picard"] == "NS(X)=S"


def test_all_17_walls_are_integral_and_span_full_picard(
    certificate: dict,
) -> None:
    chamber = certificate["transported_chamber"]
    assert chamber["id"] == "K3-D4-H12"
    assert chamber["hyperbolic_dimension"] == 12
    assert chamber["facet_count"] == len(chamber["wall_records"]) == 17
    assert chamber["wall_span_rank"] == 13
    assert chamber["wall_span_smith"] == [1] * 13
    assert chamber["wall_span_index_in_full_Picard"] == 1
    assert chamber["primitive_mirror_type_histogram"] == {
        "-2,1": 13,
        "-2,2": 2,
        "-6,3": 2,
    }
    assert all(
        all(record["checks"].values())
        for record in chamber["wall_records"]
    )


def test_strict_interior_orients_complete_chamber(certificate: dict) -> None:
    chamber = certificate["transported_chamber"]
    assert chamber["strict_interior_square"] == 1110
    assert chamber["strict_wall_pairings"] == [2] * 13 + [54, 54, 58, 58]
    assert chamber["minimum_wall_pairing"] == 2


def test_literal_and_optimized_volume_gates(certificate: dict) -> None:
    volume = certificate["finite_volume"]
    assert volume["classification"] == (
        "finite_volume_noncompact_Coxeter_polytope"
    )
    assert (
        volume["finite_vertices"],
        volume["cusps"],
        volume["hyperideal_extreme_rays"],
    ) == (93, 3, 0)
    assert volume["elliptic_rank11_edges"] == 623
    assert volume["endpoint_type_histogram"] == {
        "0,2": 3,
        "1,1": 124,
        "2,0": 496,
    }
    assert all(volume["checks"].values())


def test_full_group_lineage_and_scope(certificate: dict) -> None:
    lineage = certificate["reflection_group_lineage"]
    scope = certificate["scope"]
    assert lineage["full_chamber"] is True
    assert lineage["index"] == 1
    assert lineage["known_source"] == (
        "Scharlau--Kirschmer List12/Gram12 entry 6"
    )
    assert all(lineage["checks"].values())
    assert scope["not_new_abstract_polytope"] is True
    assert "square-minus-six" in scope["not_full_nef"]
    assert "not claimed" in scope["triality_not_surface_automorphism"]
    assert "not asserted" in scope["reflections_not_surface_automorphisms"]
    assert scope["equation_level"] is False
    assert scope["worldwide_priority_claim"] == "none"
    assert all(certificate["checks"].values())
    assert certificate["all_checks_pass"] is True


def test_report_freezes_positive_and_caveats() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results/k3_d4_h12_full_picard_k3_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "toot" in normalized
    assert "seven-dimensional" in normalized
    assert "93" in normalized and "three cusps" in normalized
    assert "not a surface automorphism" in normalized
    assert "not the full nef cone" in normalized

