from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_explicit_family_has_very_general_full_picard(certificate: dict) -> None:
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


def test_explicit_GL11Z_bridge_and_primitive_embedding(
    certificate: dict,
) -> None:
    lattice = certificate["full_Picard_lattice"]
    assert lattice["name"] == "F=U+E8(-2)+<-2>"
    assert lattice["signature"] == [1, 10, 0]
    assert lattice["determinant"] == 512
    assert lattice["smith_invariants"] == [1, 1] + [2] * 9
    assert lattice["isometry_determinant"] == 1
    assert lattice["standard_to_archived_fixed_GL11Z"][0] == [
        0,
        0,
        1,
        0,
        2,
        0,
        -2,
        1,
        0,
        0,
        1,
    ]
    assert all(lattice["checks"].values())
    embedding = lattice["primitive_K3_embedding"]
    assert embedding["orthogonal_complement"] == "U+E8(-2)+<2>"
    assert embedding["combined_embedding_complement_index"] == 512
    assert all(embedding["checks"].values())


def test_all_13_walls_transport_and_span_full_picard(certificate: dict) -> None:
    chamber = certificate["transported_H10_13_chamber"]
    assert chamber["id"] == "K3-S3-C2-H10-13"
    assert chamber["hyperbolic_dimension"] == 10
    assert chamber["facet_count"] == len(chamber["wall_records"]) == 13
    assert chamber["known_full_reflection_chamber"] is True
    assert chamber["wall_span_rank"] == 11
    assert chamber["wall_span_smith"] == [1] * 11
    assert chamber["wall_span_index"] == 1
    assert chamber["primitive_mirror_type_histogram"] == {
        "-2,1": 2,
        "-2,2": 2,
        "-4,2": 9,
    }
    assert all(
        all(record["checks"].values()) for record in chamber["wall_records"]
    )
    assert all(chamber["checks"].values())


def test_square358_class_strictly_orients_the_chamber(certificate: dict) -> None:
    chamber = certificate["transported_H10_13_chamber"]
    assert chamber["square358_strict_interior_standard_F_coordinates"] == [
        -42,
        -40,
        2,
        -6,
        -7,
        -11,
        -8,
        -14,
        3,
        21,
        19,
    ]
    assert chamber["square358_strict_interior_pairings"] == [
        2,
        2,
        2,
        2,
        2,
        6,
        2,
        2,
        2,
        2,
        2,
        2,
        38,
    ]
    assert chamber["minimum_wall_pairing"] == 2


def test_fresh_literal_and_optimized_finite_volume_replays(
    certificate: dict,
) -> None:
    volume = certificate["finite_volume"]
    assert volume["classification"] == "finite_volume_noncompact"
    assert (
        volume["finite_vertices"],
        volume["cusps"],
        volume["hyperideal_vertices"],
    ) == (34, 4, 0)
    assert volume["elliptic_rank9_edges"] == 208
    assert volume["endpoint_type_histogram"] == {
        "0,2": 5,
        "1,1": 66,
        "2,0": 137,
    }
    assert volume["literal_bad_edges"] == 0
    assert volume["optimized_bad_edges"] == 0
    assert volume["parabolic_neighborhood_audits"] == 7
    assert all(volume["checks"].values())


def test_three_chambers_share_one_equation_family(certificate: dict) -> None:
    gain = certificate["shared_carrier_gain"]
    assert gain["one_explicit_family_carries"] == [
        "K3-S3-C2-H10-13",
        "A15_11-C2-H10-21",
        "X2-RAT-C2-H10-37",
    ]
    assert gain["facet_counts"] == [13, 21, 37]
    assert "three distinct certified H10" in gain["meaning"]


def test_scope_and_report_do_not_overclaim(certificate: dict) -> None:
    scope = certificate["scope"]
    assert scope["not_new_abstract_polytope"] is True
    assert "not asserted" in scope["special_surface_geometry_kept_separate"]
    assert "no claim" in scope["not_full_nef"]
    assert "not asserted" in scope["not_surface_automorphisms"]
    assert scope["worldwide_priority_claim"] == "none"
    assert all(certificate["checks"].values())
    assert certificate["all_checks_pass"] is True

    report = (
        Path(__file__).resolve().parents[1]
        / "results/k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "toot" in normalized
    assert "determinant-one" in normalized
    assert "13, 21, and 37" in normalized
    assert "34 finite" in normalized and "four cusps" in normalized
    assert "not the full nef cone" in normalized
    assert "not surface automorphisms" in normalized

