from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_a16_h11_29_fermat_char3_scaled_dual_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_odd_closure_dual_is_actual_fermat_fixed_divisor_lattice(
    certificate: dict,
) -> None:
    closure = certificate["stable_closure"]
    bridge = certificate["scaled_dual_bridge"]
    target = bridge["target"]
    assert closure["operation"].endswith("odd index-eight closure")
    assert closure["determinant"] == -4
    assert closure["smith_invariants"] == [1] * 10 + [2] * 2
    assert closure["parity"] == "odd"
    assert bridge["operation"] == "N=M^vee(2), with Gram 2*Gram(M)^(-1)"
    assert bridge["N_determinant"] == -1024
    assert bridge["N_smith_invariants"] == [1] * 2 + [2] * 10
    assert bridge["N_parity"] == "even"
    assert target["divisor_lattice"] == "NS(X)^C2=U+E8(-2)+<-2>^2"
    assert target["existing_Fermat_section_wall_count"] == 206
    assert target["known_full_reflection_chamber_wall_count"] == 14
    assert bridge["isometry_determinant"] == 1
    assert all(bridge["checks"].values())


def test_exceptional_U_split_and_explicit_complement_bridge(
    certificate: dict,
) -> None:
    chain = certificate["scaled_dual_bridge"]["constructive_isometry_chain"]
    assert chain["exceptional_isotropic_vector_u"] == [
        3,
        4,
        2,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
    ]
    assert chain["isotropic_partner_v"] == [
        4,
        3,
        2,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        0,
    ]
    assert len(chain["positive_complement_gram"]) == 10
    assert len(chain["positive_complement_to_E8_2_A1_2_squared"]) == 10
    assert len(chain["standard_basis_columns_in_N"]) == 12


def test_all_29_halfspaces_reflections_and_full_generation(
    certificate: dict,
) -> None:
    records = certificate["facet_records"]
    generation = certificate["full_lattice_generation"]
    identity = [[1 if i == j else 0 for j in range(12)] for i in range(12)]
    assert len(records) == 29
    assert {row["source_presentation_orientation"] for row in records} == {1}
    assert {row["positive_halfspace_factor_from_M"] for row in records} == {
        "1",
        "2",
    }
    assert {row["positive_halfspace_factor_from_original_L"] for row in records} == {
        "1",
        "2",
    }
    assert all(all(row["checks"].values()) for row in records)
    assert generation["generation_index"] == 1
    assert generation["N_root_matrix_HNF"] == identity
    assert generation["Fermat_root_matrix_HNF"] == identity
    assert generation["strict_interior_square"] == 436
    assert generation["minimum_wall_pairing"] == 4
    assert all(generation["checks"].values())


def test_mirror_exchange_and_literal_vinberg_gate(certificate: dict) -> None:
    closure = certificate["stable_closure"]
    bridge = certificate["scaled_dual_bridge"]
    vinberg = certificate["inherited_literal_Vinberg_audit"]
    assert closure["primitive_mirror_types_square_divisibility"] == {
        "-4,2": 5,
        "-2,1": 17,
        "-1,1": 7,
    }
    assert bridge["primitive_target_mirror_types_square_divisibility"] == {
        "-4,2": 17,
        "-2,1": 5,
        "-2,2": 7,
    }
    assert vinberg["classification"] == "finite_volume_noncompact_Coxeter_polytope"
    assert (vinberg["dimension"], vinberg["facet_count"]) == (11, 29)
    assert (vinberg["finite_vertex_count"], vinberg["ideal_vertex_count"]) == (
        3876,
        253,
    )
    assert vinberg["elliptic_edge_count"] == 28798
    assert vinberg["connected_affine_component_count"] == 668
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighborhood_count"] == 0
    assert all(vinberg["checks"].values())


def test_exact_scope_and_report(certificate: dict) -> None:
    scope = certificate["lineage_and_scope"]
    assert "determinant-16" in scope["prior_full_rank_complex_K3_carrier"]
    assert "determinant-1024" in scope["new_positive_gain"]
    assert "206-wall" in scope["target_already_has_other_chambers"]
    assert scope["not_a_new_abstract_polytope"] is True
    assert scope["not_a_new_lattice_or_maximal_reflection_group"] is True
    assert "seventeen primitive square-minus-four" in scope[
        "not_the_complete_Fermat_nef_chamber"
    ]
    assert scope["not_an_automorphism_fixed_chamber_claim"] is True
    assert scope["worldwide_priority_claim"] == "none"
    assert certificate["all_checks_pass"] is True

    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "a16_h11_29_fermat_char3_scaled_dual_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "positive theorem" in normalized
    assert "supersingular fermat quartic" in normalized
    assert "exceptional primitive isotropic" in normalized
    assert "3,876" in normalized
    assert "additional algebraic carrier" in normalized
    assert "not the complete nef chamber" in normalized
    assert "no surface automorphism" in normalized
    assert "no new abstract-polytope" in normalized
