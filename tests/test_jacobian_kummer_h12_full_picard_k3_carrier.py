from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_jacobian_kummer_h12_full_picard_k3_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_rank13_carrier_and_nikulin_stable_range(certificate: dict) -> None:
    carrier = certificate["chamber"]["carrier"]
    embedding = carrier["primitive_K3_embedding"]
    assert carrier["name"] == "U(2)+D5(-1)^2+<-4>"
    assert carrier["signature"] == [1, 12, 0]
    assert carrier["determinant_absolute"] == 256
    assert carrier["smith_invariants"] == [1] * 8 + [2] * 2 + [4] * 3
    assert carrier["discriminant"]["generator_orders"] == [2, 2, 4, 4, 4]
    assert carrier["discriminant"]["length"] == 5
    assert carrier["canonical_discriminant_form"] == (
        "u(2) + <3/4> + <3/4> + <7/4>"
    )
    assert embedding["criterion"] == "Nikulin primitive-embedding stable range"
    assert embedding["complement_rank"] == 9
    assert embedding["complement_signature"] == [2, 7]
    assert embedding["inequality"] == "9 >= 5+2"
    assert embedding["primitive_embedding_exists"] is True
    assert all(embedding["checks"].values())
    period = carrier["very_general_period_consequence"]
    assert period["moduli_dimension"] == 7
    assert period["full_Picard_lattice"] == "NS(X)=U(2)+D5(-1)^2+<-4>"


def test_all_172_oriented_walls_act_on_and_generate_picard(
    certificate: dict,
) -> None:
    chamber = certificate["chamber"]
    walls = chamber["oriented_wall_system"]
    assert chamber["facet_count"] == len(walls["wall_records"]) == 172
    assert walls["strict_interior_coordinates"] == [
        -6,
        -6,
        -6,
        -6,
        -6,
        -6,
        11,
        13,
        11,
        14,
        11,
        13,
        17,
    ]
    assert walls["strict_interior_square"] == 8
    assert walls["minimum_wall_pairing"] == 1
    assert walls["orientation_sign_histogram"] == {"1": 172}
    assert walls["primitive_mirror_type_histogram"] == {
        "-2,1": 16,
        "-4,2": 36,
        "-4,4": 96,
        "-8,4": 24,
    }
    assert all(all(row["checks"].values()) for row in walls["wall_records"])
    assert walls["wall_span"]["rank"] == 13
    assert walls["wall_span"]["index_in_carrier_lattice"] == 1
    assert walls["wall_span"]["smith_nonzero_invariants"] == [1] * 13


def test_literal_vinberg_volume_gate(certificate: dict) -> None:
    volume = certificate["chamber"]["finite_volume"]["literal_Vinberg"]
    assert volume["classification"] == "finite_volume_noncompact"
    assert volume["elliptic_rank11_edges"] == 910316
    assert volume["endpoint_count_histogram"] == {"2": 910316}
    assert (volume["finite_vertices"], volume["cusps"]) == (113274, 6214)
    assert volume["bad_edges"] == 0
    assert volume["unrecognized_endpoint_edges"] == 0
    assert volume["invalid_full_zero_neighborhoods"] == 0


def test_independent_exact_cone_volume_gate(certificate: dict) -> None:
    volume = certificate["chamber"]["finite_volume"][
        "independent_exact_GMP_cone"
    ]
    assert volume["classification"] == "finite_volume_noncompact"
    assert volume["extreme_rays"] == 119488
    assert (volume["finite_vertices"], volume["cusps"]) == (113274, 6214)
    assert volume["hyperideal_rays"] == 0
    assert all(volume["all_exact_replay_checks"].values())


def test_scaled_dual_is_censused_but_not_needed(certificate: dict) -> None:
    dual = certificate["chamber"]["carrier"]["minimal_even_scaled_dual"]
    assert dual["minimal_scale"] == 8
    assert dual["determinant_absolute"] == 2**31
    assert dual["smith_invariants"] == [2] * 3 + [4] * 2 + [8] * 8
    assert dual["discriminant_length"] == 13


def test_positive_gain_and_boundary_are_exact(certificate: dict) -> None:
    assert certificate["all_checks_pass"] is True
    assert all(certificate["checks"].values())
    scope = certificate["scope"]
    assert "full Picard lattice" in scope["positive_gain"]
    assert scope["not_new_abstract_polytope"] is True
    assert scope["not_first_algebraic_occurrence"] is True
    assert "square-minus-four" in scope["not_full_nef"]
    assert "not asserted" in scope["not_surface_automorphisms"]
    assert scope["worldwide_priority_claim"] == "none"


def test_report_freezes_the_jk_h12_carrier_theorem() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results/jacobian_kummer_h12_full_picard_k3_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "positive theorem" in normalized
    assert "jk-h12-172" in normalized
    assert "u(2)+d5(-1)^2+<-4>" in normalized
    assert "seven-dimensional" in normalized
    assert "910,316" in normalized
    assert "119,488" in normalized
    assert "not full-nef" in normalized
    assert "not surface automorphisms" in normalized
