from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_three_h11_full_picard_k3_carriers import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def _chamber(certificate: dict, chamber_id: str) -> dict:
    return next(row for row in certificate["chambers"] if row["id"] == chamber_id)


def _walls(chamber: dict) -> dict:
    return chamber.get(
        "oriented_wall_system",
        chamber.get("oriented_wall_system_on_closure"),
    )


def test_c2_explicit_full_picard_k3_carrier(certificate: dict) -> None:
    chamber = _chamber(certificate, "F3-C2-H11-206")
    lattice = chamber["original_wall_lattice"]
    embedding = lattice["primitive_K3_embedding"]
    assert lattice["name"] == "U+E8(-2)+<-2>^2"
    assert lattice["signature"] == [1, 11, 0]
    assert lattice["determinant_absolute"] == 1024
    assert lattice["smith_invariants"] == [1, 1] + [2] * 10
    assert lattice["discriminant"]["length"] == 10
    assert embedding["criterion"] == "explicit full primitive embedding"
    assert embedding["orthogonal_complement"] == "E8(-2)+<2>^2"
    assert embedding["combined_index"] == 1024
    assert all(embedding["checks"].values())
    assert lattice["very_general_period_consequence"]["moduli_dimension"] == 8
    assert lattice["minimal_even_scaled_dual"]["minimal_scale"] == 4
    assert lattice["minimal_even_scaled_dual"]["determinant_absolute"] == 16384


def test_c2_all_206_walls_and_literal_vinberg(certificate: dict) -> None:
    chamber = _chamber(certificate, "F3-C2-H11-206")
    walls = _walls(chamber)
    assert chamber["facet_count"] == len(walls["wall_records"]) == 206
    assert walls["strict_interior_square"] == 4
    assert walls["minimum_wall_pairing"] == 1
    assert walls["orientation_sign_histogram"] == {"-1": 103, "1": 103}
    assert walls["primitive_mirror_type_histogram"] == {
        "-2,1": 12,
        "-2,2": 164,
        "-4,2": 30,
    }
    assert all(all(row["checks"].values()) for row in walls["wall_records"])
    assert walls["wall_span"]["rank"] == 12
    assert walls["wall_span"]["index_in_carrier_lattice"] == 1
    volume = chamber["volume"]
    assert (volume["finite_vertices"], volume["cusps"]) == (327566, 18354)
    assert volume["elliptic_rank10_edges"] == 2688957
    assert volume["bad_edges"] == volume["bad_parabolic_neighborhoods"] == 0


def test_v4_mandatory_closure_is_the_full_picard_carrier(
    certificate: dict,
) -> None:
    chamber = _chamber(certificate, "F3-V4-H11-328")
    original = chamber["original_wall_lattice"]
    closure = chamber["certified_reflection_stable_overclosures"][1]
    lattice = closure["carrier_lattice"]
    embedding = lattice["primitive_K3_embedding"]
    assert original["determinant_absolute"] == 768
    assert original["smith_invariants"] == [1] * 6 + [2] * 4 + [4, 12]
    assert original["reflection_failure"] == "36 of 328 reflections are nonintegral"
    assert closure["index"] == 4
    assert closure["unique_minimal_for_all_328_reflections"] is True
    assert lattice["determinant_absolute"] == 48
    assert lattice["smith_invariants"] == [1] * 8 + [2] * 3 + [6]
    assert lattice["discriminant"]["length"] == 4
    assert embedding["inequality"] == "10 >= 4+2"
    assert embedding["primitive_embedding_exists"] is True
    assert all(embedding["checks"].values())
    assert lattice["very_general_period_consequence"]["moduli_dimension"] == 8
    assert lattice["minimal_even_scaled_dual"]["minimal_scale"] == 6
    assert original["minimal_even_scaled_dual"]["minimal_scale"] == 12


def test_v4_all_328_walls_and_exact_cone(certificate: dict) -> None:
    chamber = _chamber(certificate, "F3-V4-H11-328")
    walls = _walls(chamber)
    assert chamber["facet_count"] == len(walls["wall_records"]) == 328
    assert walls["strict_interior_square"] == 4
    assert walls["minimum_wall_pairing"] == 1
    assert walls["orientation_sign_histogram"] == {"-1": 158, "1": 170}
    assert walls["primitive_mirror_type_histogram"] == {
        "-2,1": 28,
        "-4,2": 24,
        "-6,3": 276,
    }
    assert all(all(row["checks"].values()) for row in walls["wall_records"])
    assert walls["wall_span"]["rank"] == 12
    assert walls["wall_span"]["index_in_carrier_lattice"] == 1
    volume = chamber["volume"]
    assert (volume["finite_vertices"], volume["cusps"]) == (225408, 11554)
    assert volume["hyperideal_vertices"] == 0


def test_ss2_3b_stable_range_full_picard_carrier_and_vinberg(
    certificate: dict,
) -> None:
    chamber = _chamber(certificate, "SS2-3B-H11-22")
    lattice = chamber["original_wall_lattice"]
    embedding = lattice["primitive_K3_embedding"]
    walls = _walls(chamber)
    assert lattice["name"] == "U+A2(-1)^5"
    assert lattice["signature"] == [1, 11, 0]
    assert lattice["determinant_absolute"] == 243
    assert lattice["smith_invariants"] == [1] * 7 + [3] * 5
    assert lattice["discriminant"]["length"] == 5
    assert embedding["inequality"] == "10 >= 5+2"
    assert embedding["primitive_embedding_exists"] is True
    assert all(embedding["checks"].values())
    assert lattice["very_general_period_consequence"]["moduli_dimension"] == 8
    assert lattice["minimal_even_scaled_dual"]["minimal_scale"] == 3
    assert chamber["facet_count"] == len(walls["wall_records"]) == 22
    assert walls["primitive_mirror_type_histogram"] == {"-2,1": 12, "-6,3": 10}
    assert walls["wall_span"]["rank"] == 12
    assert walls["wall_span"]["index_in_carrier_lattice"] == 1
    assert all(all(row["checks"].values()) for row in walls["wall_records"])
    volume = chamber["volume"]
    assert (volume["finite_vertices"], volume["cusps"]) == (2304, 72)
    assert volume["elliptic_rank10_edges"] == 14580
    assert volume["bad_edges"] == volume["bad_parabolic_neighborhoods"] == 0


def test_joint_positive_theorem_has_exact_boundary(certificate: dict) -> None:
    assert certificate["all_checks_pass"] is True
    assert all(certificate["checks"].values())
    scope = certificate["scope"]
    assert "three known cofinite H11" in scope["positive_gain"]
    assert scope["not_new_abstract_polytopes"] is True
    assert "mixed square-minus-four" in scope["not_full_nef"]
    assert "not asserted" in scope["not_surface_automorphisms"]
    assert scope["all_wall_sets_generate_full_Picard"] is True
    assert scope["worldwide_priority_claim"] == "none"


def test_report_freezes_the_three_h11_positive_theorem() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results/three_h11_full_picard_k3_carriers.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "positive theorem" in normalized
    assert "f3-c2-h11-206" in normalized
    assert "f3-v4-h11-328" in normalized
    assert "ss2-3b-h11-22" in normalized
    assert "eight-dimensional" in normalized
    assert "index-four" in normalized
    assert "not full-nef" in normalized
    assert "not surface automorphisms" in normalized
