from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_a18_c03_ud4a1x6_weierstrass_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_equation_open_and_fibers(certificate: dict) -> None:
    family = certificate["weierstrass_family"]
    assert len(family["parameters"]) == 8
    assert family["fibers"] == (
        "I0* at t=0; I2 at the six roots of S; I1 at the six roots of H"
    )
    assert family["nonempty_witness"]["H_at_zero"] == 940
    assert family["nonempty_witness"]["two_torsion_root_ideal_groebner_basis"] == [
        "1"
    ]
    assert all(family["checks"].values())


def test_exact_eight_moduli_and_full_picard(certificate: dict) -> None:
    tangent = certificate["weierstrass_family"]["moduli_tangent_certificate"]
    picard = certificate["very_general_picard_lattice"]
    assert tangent == {
        "family_rank": 8,
        "coordinate_orbit_rank": 4,
        "combined_rank": 12,
        "minor_rows_zero_based": [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15],
        "minor_determinant": -37705049539909728,
    }
    assert picard["lattice"] == "U+D4(-1)+A1(-1)^6"
    assert picard["rank_signature_determinant"] == [12, [1, 11, 0], 256]
    assert picard["smith"] == [1] * 4 + [2] * 8
    assert "unit ideal" in picard["torsion_saturation"]
    assert all(picard["checks"].values())


def test_exact_marking_and_all_28_walls(certificate: dict) -> None:
    bridge = certificate["equation_to_H11_28_bridge"]
    assert bridge["determinant_Q"] == -1
    assert len(bridge["standard_to_archived_N_columns_Q"]) == 12
    assert len(bridge["mirror_records"]) == 28
    assert bridge["generation_index"] == 1
    assert bridge["root_matrix_HNF"] == [
        [1 if i == j else 0 for j in range(12)] for i in range(12)
    ]
    assert bridge["strict_interior_square"] == 608
    assert bridge["minimum_wall_pairing"] == 2
    assert all(bridge["checks"].values())
    assert all(
        all(record["checks"].values()) for record in bridge["mirror_records"]
    )


def test_literal_vinberg_and_scope(certificate: dict) -> None:
    vinberg = certificate["literal_Vinberg_audit"]
    scope = certificate["scope"]
    assert vinberg["classification"] == "finite_volume_noncompact_Coxeter_polytope"
    assert (vinberg["dimension"], vinberg["facet_count"]) == (11, 28)
    assert (vinberg["finite_vertex_count"], vinberg["ideal_vertex_count"]) == (
        1252,
        64,
    )
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighborhood_count"] == 0
    assert scope["avoids_global_two_torsion_overclosure"] is True
    assert scope["not_the_complete_nef_chamber"] is True
    assert scope["strict_interior_not_asserted_ample"] is True
    assert scope["no_surface_automorphism_claim"] is True
    assert certificate["all_checks_pass"] is True


def test_report(certificate: dict) -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "a18_c03_ud4a1x6_weierstrass_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "positive theorem" in normalized
    assert "eight-moduli" in normalized
    assert "i0*+6i2+6i1" in normalized
    assert "unit ideal" in normalized
    assert "all 28" in normalized
    assert "finite-volume noncompact" in normalized
    assert "not the complete nef" in normalized

