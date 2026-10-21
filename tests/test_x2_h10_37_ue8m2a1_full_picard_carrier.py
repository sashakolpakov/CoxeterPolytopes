from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_x2_h10_37_ue8m2a1_full_picard_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_explicit_nine_parameter_family_has_full_generic_picard(
    certificate: dict,
) -> None:
    family = certificate["explicit_family"]
    assert family["parameter_count"] == 9
    assert family["very_general_full_Picard"] == (
        "NS(X)=U+E8(-2)+<-2> for a very general parameter"
    )
    assert family["fibers"] == {
        "I1": 22,
        "I2": 1,
        "smooth_at_infinity": 1,
    }
    assert family["nonempty_open_witness"]
    assert all(family["checks"].values())


def test_lattice_and_explicit_primitive_k3_embedding(
    certificate: dict,
) -> None:
    lattice = certificate["full_Picard_lattice"]
    embedding = lattice["primitive_K3_embedding"]
    assert lattice["name"] == "F=U+E8(-2)+<-2>"
    assert lattice["rank"] == 11
    assert lattice["signature"] == [1, 10, 0]
    assert lattice["determinant"] == 512
    assert lattice["smith_invariants"] == [1, 1] + [2] * 9
    assert lattice["nikulin_invariants"] == {"r": 11, "a": 9, "delta": 1}
    assert lattice["K3_moduli_dimension"] == 9
    assert embedding["orthogonal_complement"] == "U+E8(-2)+<2>"
    assert embedding["combined_embedding_complement_index"] == 512
    assert all(embedding["checks"].values())


def test_all_37_walls_transport_and_generate_full_picard(
    certificate: dict,
) -> None:
    chamber = certificate["transported_H10_37_chamber"]
    assert chamber["id"] == "X2-RAT-C2-H10-37"
    assert chamber["hyperbolic_dimension"] == 10
    assert chamber["facet_count"] == len(chamber["wall_records"]) == 37
    assert chamber["walls_span_full_Picard"] is True
    assert chamber["wall_span_rank"] == 11
    assert chamber["wall_span_smith"] == [1] * 11
    assert chamber["wall_span_index"] == 1
    assert chamber["primitive_mirror_type_histogram"] == {
        "-2,1": 6,
        "-2,2": 21,
        "-4,2": 10,
    }
    assert all(
        all(record["checks"].values())
        for record in chamber["wall_records"]
    )
    assert all(chamber["checks"].values())


def test_square12_class_orients_every_halfspace(certificate: dict) -> None:
    chamber = certificate["transported_H10_37_chamber"]
    assert chamber["square12_strict_interior_coordinates"] == [
        226,
        11,
        -44,
        -48,
        -69,
        -97,
        -75,
        -43,
        -12,
        -5,
        36,
    ]
    assert chamber["minimum_wall_pairing"] == 1
    assert len(chamber["square12_strict_interior_pairings"]) == 37
    assert min(chamber["square12_strict_interior_pairings"]) == 1


def test_literal_vinberg_finite_volume_gate(certificate: dict) -> None:
    volume = certificate["finite_volume"]
    assert volume["classification"] == (
        "finite_volume_noncompact_Coxeter_polytope"
    )
    assert (volume["finite_vertices"], volume["cusps"]) == (3236, 261)
    assert volume["elliptic_rank9_edges"] == 24350
    assert volume["bad_edges"] == 0
    assert volume["bad_parabolic_neighborhoods"] == 0
    assert all(volume["checks"].values())


def test_positive_gain_and_scope_are_exact(certificate: dict) -> None:
    gain = certificate["shared_carrier_gain"]
    scope = certificate["scope"]
    assert gain["same_explicit_family_already_carries"] == (
        "A15_11-C2-H10-21"
    )
    assert gain["newly_attached_chamber"] == "X2-RAT-C2-H10-37"
    assert "two distinct certified H10" in gain["meaning"]
    assert "equation" in scope["positive_gain"]
    assert scope["not_new_abstract_polytope"] is True
    assert "no claim" in scope["not_full_nef"]
    assert "not asserted" in scope["not_surface_automorphisms"]
    assert scope["worldwide_priority_claim"] == "none"
    assert all(certificate["checks"].values())
    assert certificate["all_checks_pass"] is True


def test_report_freezes_the_positive_without_overclaim() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results/x2_h10_37_ue8m2a1_full_picard_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "toot" in normalized
    assert "nine-parameter" in normalized
    assert "37" in normalized and "21" in normalized
    assert "3,236" in normalized and "261" in normalized
    assert "not the full nef cone" in normalized
    assert "not surface automorphisms" in normalized

