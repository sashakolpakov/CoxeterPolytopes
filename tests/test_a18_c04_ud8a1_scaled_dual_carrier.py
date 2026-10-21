from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_a18_c04_ud8a1_scaled_dual_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_odd_closure_dual_is_actual_roulleau_k3_picard_lattice(
    certificate: dict,
) -> None:
    closure = certificate["stable_closure"]
    bridge = certificate["scaled_dual_bridge"]
    target = bridge["target"]
    assert closure["operation"].startswith("unique reflection-stable odd")
    assert closure["determinant"] == 256
    assert closure["smith_invariants"] == [1] * 3 + [2] * 8
    assert closure["parity"] == "odd"
    assert bridge["operation"] == "N=M^vee(2), with Gram 2*Gram(M)^(-1)"
    assert bridge["N_determinant"] == 8
    assert bridge["N_smith_invariants"] == [1] * 8 + [2] * 3
    assert bridge["N_parity"] == "even"
    assert target["divisor_lattice"] == "NS(X)=U+D8(-1)+A1(-1)"
    assert target["published_complete_nef_wall_count"] == 14
    assert bridge["isometry_determinant"] == 1
    assert all(bridge["checks"].values())


def test_all_28_positive_halfspaces_and_reflections_transport(
    certificate: dict,
) -> None:
    records = certificate["facet_records"]
    assert len(records) == 28
    assert sum(row["source_presentation_orientation"] == -1 for row in records) == 19
    assert sum(row["source_presentation_orientation"] == 1 for row in records) == 9
    assert {row["positive_halfspace_factor_from_M"] for row in records} == {
        "1",
        "2",
    }
    assert {row["positive_halfspace_factor_from_original_L"] for row in records} == {
        "1",
        "2",
    }
    assert all(all(row["checks"].values()) for row in records)


def test_mirror_exchange_and_full_picard_generation(certificate: dict) -> None:
    closure = certificate["stable_closure"]
    bridge = certificate["scaled_dual_bridge"]
    generation = certificate["full_lattice_generation"]
    identity = [[1 if i == j else 0 for j in range(11)] for i in range(11)]
    assert closure["primitive_mirror_types_square_divisibility"] == {
        "-4,2": 16,
        "-2,1": 4,
        "-1,1": 8,
    }
    assert bridge["primitive_target_mirror_types_square_divisibility"] == {
        "-4,2": 4,
        "-2,1": 16,
        "-2,2": 8,
    }
    assert generation["generation_index"] == 1
    assert generation["N_root_matrix_HNF"] == identity
    assert generation["K3_root_matrix_HNF"] == identity
    assert generation["strict_interior_square"] == 508
    assert generation["minimum_wall_pairing"] == 3
    assert all(generation["checks"].values())


def test_literal_vinberg_gate_and_exact_scope(certificate: dict) -> None:
    vinberg = certificate["inherited_literal_Vinberg_audit"]
    scope = certificate["lineage_and_scope"]
    assert vinberg["classification"] == "finite_volume_noncompact_Coxeter_polytope"
    assert (vinberg["dimension"], vinberg["facet_count"]) == (10, 28)
    assert (vinberg["finite_vertex_count"], vinberg["ideal_vertex_count"]) == (
        1087,
        74,
    )
    assert vinberg["elliptic_edge_count"] == 6956
    assert vinberg["connected_affine_component_count"] == 204
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighborhood_count"] == 0
    assert all(vinberg["checks"].values())
    assert "determinant 1024" in scope["existing_geometric_carrier"]
    assert "determinant-eight" in scope["new_positive_gain"]
    assert scope["not_a_new_abstract_polytope"] is True
    assert "four primitive square-minus-four" in scope[
        "not_the_complete_K3_nef_chamber"
    ]
    assert scope["not_an_automorphism_fixed_chamber_claim"] is True
    assert scope["worldwide_priority_claim"] == "none"
    assert certificate["all_checks_pass"] is True


def test_report_freezes_positive_gain_and_boundaries() -> None:
    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "a18_c04_ud8a1_scaled_dual_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "positive theorem" in normalized
    assert "u+d8(-1)+a1(-1)" in normalized
    assert "16 | `(-2, div 1)`" in normalized
    assert "1,087" in normalized
    assert "second, arithmetically different k3" in normalized
    assert "not the k3's complete nef chamber" in normalized
    assert "no surface automorphism" in normalized
    assert "no new abstract-polytope" in normalized

