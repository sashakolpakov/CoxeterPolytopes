from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_ss2_2c_h14_98_hilb2_scaled_dual_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_unique_odd_closure_and_actual_hilbert_square_target(
    certificate: dict,
) -> None:
    closure = certificate["stable_closure"]
    bridge = certificate["scaled_dual_bridge"]
    target = bridge["target"]
    assert closure["operation"].startswith("unique reflection-stable odd")
    assert closure["determinant"] == 64
    assert closure["smith_invariants"] == [1] * 9 + [2] * 6
    assert closure["parity"] == "odd"
    assert bridge["operation"] == "N=M^vee(2), with Gram 2*Gram(M)^(-1)"
    assert bridge["N_determinant"] == 512
    assert bridge["N_smith_invariants"] == [1] * 6 + [2] * 9
    assert bridge["N_parity"] == "even"
    assert target["divisor_lattice"] == "Pic(Hilb^2(X))=(15,9,1)"
    assert bridge["isometry_determinant"] == 1
    assert all(bridge["checks"].values())


def test_all_98_positive_halfspaces_and_reflections_transport(
    certificate: dict,
) -> None:
    records = certificate["facet_records"]
    assert len(records) == 98
    assert {row["positive_halfspace_factor_from_M"] for row in records} == {
        "1",
        "2",
    }
    assert {row["positive_halfspace_factor_from_original_L"] for row in records} == {
        "1",
        "2",
    }
    assert all(all(row["checks"].values()) for row in records)
    assert sorted(row["alexeev_wall_zero_based"] for row in records) == list(range(98))


def test_primitive_type_exchange_and_full_lattice_generation(
    certificate: dict,
) -> None:
    closure = certificate["stable_closure"]
    bridge = certificate["scaled_dual_bridge"]
    generation = certificate["full_lattice_generation"]
    identity = [[1 if i == j else 0 for j in range(15)] for i in range(15)]
    assert closure["primitive_mirror_types_square_divisibility"] == {
        "-4,2": 14,
        "-2,1": 14,
        "-1,1": 70,
    }
    assert bridge["primitive_target_mirror_types_square_divisibility"] == {
        "-4,2": 14,
        "-2,1": 14,
        "-2,2": 70,
    }
    assert generation["generation_index"] == 1
    assert generation["N_root_matrix_HNF"] == identity
    assert generation["Hilb2_root_matrix_HNF"] == identity
    assert generation["strict_interior_square"] == 40432
    assert generation["minimum_wall_pairing"] == 76
    assert all(generation["checks"].values())


def test_literal_vinberg_gate_and_scope(certificate: dict) -> None:
    vinberg = certificate["inherited_literal_Vinberg_audit"]
    scope = certificate["scope"]
    assert vinberg["classification"] == "finite_volume_noncompact_Coxeter_polytope"
    assert (vinberg["dimension"], vinberg["facet_count"]) == (14, 98)
    assert (vinberg["finite_vertex_count"], vinberg["ideal_vertex_count"]) == (
        75898,
        2114,
    )
    assert vinberg["elliptic_edge_count"] == 646786
    assert vinberg["connected_affine_component_count"] == 5838
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighborhood_count"] == 0
    assert all(vinberg["checks"].values())
    assert scope["known_abstract_polytope"] is True
    assert scope["not_a_new_geometric_chamber_slot"] is True
    assert scope["not_a_Kahler_or_movable_cone_claim"] is True
    assert scope["not_an_automorphism_fixed_chamber_claim"] is True
    assert scope["worldwide_priority_claim"] == "none"
    assert certificate["all_checks_pass"] is True


def test_report_freezes_positive_gain_and_boundaries() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "ss2_2c_h14_98_hilb2_scaled_dual_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "positive theorem" in normalized
    assert "unique integral index-two overclosure" in normalized
    assert "pic(hilb^2(x))" in normalized
    assert "70 * (-2,div2)" in normalized
    assert "75,898" in normalized
    assert "finite-volume and noncompact" in normalized
    assert "does **not** add a new abstract" in normalized
    assert "no claim is made" in normalized
    assert "kähler or movable cone" in normalized
