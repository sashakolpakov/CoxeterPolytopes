from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_a18_c03_ud4a1x6_scaled_dual_carrier import (
    OUTPUT,
    build_certificate,
)


@pytest.fixture(scope="module")
def certificate() -> dict:
    rebuilt = build_certificate()
    archived = json.loads(Path(OUTPUT).read_text(encoding="utf-8"))
    assert rebuilt == archived
    return rebuilt


def test_scaled_dual_is_full_picard_k3_lattice(certificate: dict) -> None:
    closure = certificate["stable_closure"]
    bridge = certificate["scaled_dual_K3_bridge"]
    assert closure["determinant"] == -16
    assert closure["smith_invariants"] == [1] * 8 + [2] * 4
    assert closure["parity"] == "odd"
    assert bridge["operation"] == "N=M^vee(2), with Gram 2*Gram(M)^(-1)"
    assert bridge["N_determinant"] == -256
    assert bridge["N_smith_invariants"] == [1] * 4 + [2] * 8
    assert bridge["N_parity"] == "even"
    assert bridge["Nikulin_invariants"] == {"r": 12, "a": 8, "delta": 1}
    assert bridge["discriminant_group"] == "(Z/2)^8"
    assert bridge["standard_representative"]["lattice"] == (
        "U+D4(-1)+A1(-1)^6"
    )
    assert bridge["period_consequence"]["N_polarized_moduli_dimension"] == 8
    assert bridge["period_consequence"]["not_merely_a_sublattice"] is True
    assert all(bridge["checks"].values())


def test_explicit_primitive_k3_graph_glue(certificate: dict) -> None:
    glue = certificate["scaled_dual_K3_bridge"]["primitive_K3_embedding"]
    assert glue["complement"]["lattice"] == "T=U+U(2)+A1(-1)^6"
    assert glue["complement"]["signature"] == [2, 8, 0]
    assert glue["complement"]["determinant"] == 256
    assert len(glue["glue_vectors_in_N_plus_T_coordinates"]) == 8
    assert glue["discriminant_anti_isometry"]["all_elements_checked"] == 256
    assert "unique even unimodular lattice" in glue["classification"]
    assert all(glue["checks"].values())


def test_all_28_oriented_mirrors_transport_and_generate(certificate: dict) -> None:
    records = certificate["facet_records"]
    closure = certificate["stable_closure"]
    bridge = certificate["scaled_dual_K3_bridge"]
    generation = certificate["full_lattice_generation"]
    identity = [[1 if i == j else 0 for j in range(12)] for i in range(12)]
    assert len(records) == 28
    assert sum(row["source_presentation_orientation"] == -1 for row in records) == 17
    assert sum(row["source_presentation_orientation"] == 1 for row in records) == 11
    assert {row["positive_halfspace_factor_from_M"] for row in records} == {
        "1",
        "2",
    }
    assert all(all(row["checks"].values()) for row in records)
    assert closure["primitive_mirror_types_square_divisibility"] == {
        "-4,2": 8,
        "-2,1": 8,
        "-1,1": 12,
    }
    assert bridge["primitive_target_mirror_types_square_divisibility"] == {
        "-4,2": 8,
        "-2,1": 8,
        "-2,2": 12,
    }
    assert generation["generation_index"] == 1
    assert generation["N_root_matrix_HNF"] == identity
    assert generation["strict_interior_square"] == 608
    assert generation["minimum_wall_pairing"] == 2
    assert all(generation["checks"].values())


def test_literal_vinberg_finite_volume_gate(certificate: dict) -> None:
    vinberg = certificate["inherited_literal_Vinberg_audit"]
    assert vinberg["classification"] == "finite_volume_noncompact_Coxeter_polytope"
    assert (vinberg["dimension"], vinberg["facet_count"]) == (11, 28)
    assert (vinberg["finite_vertex_count"], vinberg["ideal_vertex_count"]) == (
        1252,
        64,
    )
    assert vinberg["elliptic_edge_count"] == 8851
    assert vinberg["connected_affine_component_count"] == 186
    assert vinberg["bad_edge_count"] == 0
    assert vinberg["bad_parabolic_neighborhood_count"] == 0
    assert all(vinberg["checks"].values())


def test_scope_and_report(certificate: dict) -> None:
    scope = certificate["lineage_and_scope"]
    assert "Smith invariants" in scope["arithmetically_different_from_original_carrier"]
    assert scope["not_a_new_abstract_polytope"] is True
    assert "eight primitive square-minus-four" in scope[
        "not_the_complete_K3_nef_chamber"
    ]
    assert scope["not_an_automorphism_fixed_chamber_claim"] is True
    assert scope["no_named_or_equation_family_claim"] is True
    assert scope["worldwide_priority_claim"] == "none"
    assert certificate["all_checks_pass"] is True

    report = (
        Path(__file__).resolve().parents[1]
        / "results"
        / "a18_c03_ud4a1x6_scaled_dual_carrier.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "positive theorem" in normalized
    assert "u + d4(-1) + a1(-1)^6" in normalized
    assert "explicit primitive embedding" in normalized
    assert "1,252" in normalized
    assert "second full-picard" in normalized
    assert "not the complete nef chamber" in normalized
    assert "no surface automorphism" in normalized
    assert "no new abstract-polytope" in normalized
