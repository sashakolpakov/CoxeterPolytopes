from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_f3_v4_jk_h10_shared_k3_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_shared_lattice_and_explicit_primitive_k3_embedding(
    certificate: dict,
) -> None:
    lattice = certificate["shared_Picard_lattice"]
    embedding = lattice["primitive_K3_embedding"]
    assert lattice["name"] == "S=U+E8(-2)+<-4>"
    assert lattice["signature"] == [1, 10, 0]
    assert lattice["determinant"] == 1024
    assert lattice["smith_invariants"] == [1, 1] + [2] * 8 + [4]
    assert lattice["discriminant_length"] == 9
    assert lattice["discriminant_form"] == "q_plus(F2^8) + <-1/4>"
    assert embedding["orthogonal_complement"] == "U+E8(-2)+<4>"
    assert embedding["combined_embedding_complement_index"] == 1024
    assert all(embedding["checks"].values())
    assert lattice["period_consequence"]["K3_moduli_dimension"] == 9
    assert lattice["period_consequence"]["very_general_full_Picard_lattice"] == (
        "NS(X)=S"
    )
    assert all(lattice["checks"].values())


def test_f3_v4_all_123_oriented_walls_and_volume(certificate: dict) -> None:
    chamber = certificate["chambers"][0]
    assert chamber["id"] == "F3-V4-H10-123"
    assert chamber["facet_count"] == len(chamber["wall_records"]) == 123
    assert chamber["strict_interior_square"] == 4
    assert chamber["minimum_wall_pairing"] == 1
    assert chamber["orientation_sign_histogram"] == {"-1": 64, "1": 59}
    assert chamber["primitive_mirror_type_histogram"] == {
        "-2,1": 10,
        "-4,2": 31,
        "-4,4": 82,
    }
    assert all(all(row["checks"].values()) for row in chamber["wall_records"])
    assert chamber["wall_span"]["index_in_full_Picard_lattice"] == 2
    assert chamber["wall_span"]["rank"] == 11
    volume = chamber["literal_Vinberg"]
    assert (volume["finite_vertices"], volume["cusps"], volume["hyperideal_vertices"]) == (
        37451,
        1732,
        0,
    )
    assert volume["elliptic_rank9_edges"] == 246610
    assert volume["bad_edges"] == volume["bad_parabolic_neighborhoods"] == 0
    assert all(chamber["checks"].values())


def test_jacobian_kummer_all_84_oriented_walls_and_volume(
    certificate: dict,
) -> None:
    chamber = certificate["chambers"][1]
    assert chamber["id"] == "JK-H10-0007"
    assert chamber["facet_count"] == len(chamber["wall_records"]) == 84
    assert chamber["strict_interior_square"] == 8
    assert chamber["minimum_wall_pairing"] == 1
    assert chamber["primitive_mirror_type_histogram"] == {
        "-2,1": 8,
        "-4,2": 33,
        "-4,4": 43,
    }
    assert all(all(row["checks"].values()) for row in chamber["wall_records"])
    assert chamber["wall_span"]["index_in_full_Picard_lattice"] == 2
    assert chamber["wall_span"]["rank"] == 11
    volume = chamber["literal_Vinberg"]
    assert (volume["finite_vertices"], volume["cusps"], volume["hyperideal_vertices"]) == (
        15698,
        862,
        0,
    )
    assert volume["elliptic_rank9_edges"] == 107233
    assert volume["bad_edges"] == volume["bad_parabolic_neighborhoods"] == 0
    assert all(chamber["checks"].values())


def test_minimal_even_scaled_dual_is_recorded_but_not_used(
    certificate: dict,
) -> None:
    dual = certificate["shared_Picard_lattice"]["minimal_even_scaled_dual"]
    assert dual["minimal_positive_even_integral_scale"] == 8
    assert dual["no_smaller_scale_works"] is True
    assert dual["determinant"] == 2**23
    assert dual["discriminant_length"] == 11
    assert "positive K3 carrier uses S itself" in dual["note"]


def test_scope_is_positive_but_exactly_bounded(certificate: dict) -> None:
    scope = certificate["scope"]
    assert "two known abstract" in scope["positive_gain"]
    assert scope["not_new_abstract_polytopes"] is True
    assert "mixed divisor-reflection" in scope["not_full_nef"]
    assert "not asserted" in scope["not_surface_automorphisms"]
    assert "index-two" in scope["wall_span_index"]
    assert scope["worldwide_priority_claim"] == "none"
    assert certificate["all_checks_pass"] is True


def test_report_freezes_the_joint_positive_theorem() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results/f3_v4_jk_h10_shared_k3_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "positive theorem" in normalized
    assert "u+e8(-2)+<-4>" in normalized
    assert "123" in normalized and "84" in normalized
    assert "nine-dimensional" in normalized
    assert "not full-nef" in normalized
    assert "not surface automorphisms" in normalized

