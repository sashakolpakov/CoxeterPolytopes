from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_r12_h10_16_order4_divisor_overlattice_carrier import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "data/enumerations/r12_h10_16_order4_divisor_overlattice_carrier.json"
)
REPORT = ROOT / "results/r12_h10_16_order4_divisor_overlattice_carrier.md"


def archived() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_certificate_replays_exactly() -> None:
    payload = archived()
    assert build_certificate() == payload
    assert payload["all_checks_pass"]
    assert payload["payload_sha256"] == (
        "53a4d077584ea71528e496eace2e24ac0f699f520118dd2c14f6900e455db829"
    )


def test_inverse_overextension_is_complete_and_unique() -> None:
    census = archived()["inverse_overextension_census"]
    assert census["all_subspace_count"] == 67
    assert census["dimension_histogram"] == {
        "0": 1,
        "1": 15,
        "2": 35,
        "3": 15,
        "4": 1,
    }
    assert census["totally_q_isotropic_subspace_count"] == 4
    assert census["action_stable_subspace_count"] == 27
    assert census["action_stable_totally_q_isotropic_subspace_count"] == 2
    assert census["nontrivial_action_stable_even_overgroup_count"] == 1
    assert census["unique_nontrivial_subgroup_elements"] == [0, 1]
    assert census["unique_glue_half_vector_bits"] == 344


def test_unique_saturation_is_the_actual_row39_action_pair() -> None:
    saturation = archived()["unique_saturation"]
    assert saturation["index_over_M"] == 2
    assert saturation["isometry_class"] == "U+E8(-1)+A1(-1)^2"
    assert saturation["rank_signature_determinant"] == [12, [1, 11, 0], -4]
    assert saturation["smith"] == [*[1] * 10, 2, 2]
    conjugacy = saturation["exact_action_pair_conjugacy_to_row39"]
    assert conjugacy["nikulin_stable_lifting_applies"]
    assert conjugacy["all_glue_elements_carried"]
    assert conjugacy["full_integral_lattice_action_pairs_conjugate"]
    assert conjugacy["coinvariant_gram"] == [[-4]]


def test_direct_divisor_sublattice_has_the_c2b_action() -> None:
    divisor = archived()["direct_row39_divisor_sublattice"]
    assert divisor["index_in_NS"] == 2
    assert divisor["rank_signature_determinant"] == [12, [1, 11, 0], -16]
    assert divisor["smith"] == [*[1] * 8, 2, 2, 2, 2]
    assert divisor["parity_functional"] == [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
    conjugacy = divisor["exact_action_pair_conjugacy_to_archived_C2b"]
    assert conjugacy["fixed_discriminant_isometry"][
        "matrix_on_nontrivial_smith_coordinates"
    ] == [[1, 1, 1], [0, 1, 1], [0, 0, 1]]
    assert conjugacy["all_glue_elements_carried"]
    assert conjugacy["full_integral_lattice_action_pairs_conjugate"]


def test_actual_surface_and_chamber_scope_are_explicit() -> None:
    payload = archived()
    carrier = payload["explicit_algebraic_carrier"]
    assert carrier["automorphism_global_order"] == 4
    assert carrier["automorphism_picard_order"] == 2
    assert carrier["brandhorst_hofmann_id"] == "0.4.4.1"
    assert carrier["very_general_NS"] == "U+E8(-1)+A1(-1)^2"
    chamber = payload["divisor_chamber_handoff"]
    assert chamber["all_archived_roots_remain_algebraic_minus_two_classes"]
    assert chamber["archived_complete_M_reflection_chamber_retained"]
    assert chamber["full_NS_nef_chamber_may_be_a_strict_refinement"]
    assert "not claimed to be the full nef cone" in chamber["scope"]


def test_literal_vinberg_gate_is_inherited_without_downgrade() -> None:
    polytope = archived()["hyperbolic_polytope"]
    assert polytope == {
        "bad_edge_count": 0,
        "canonical_full_audit_sha256": (
            "900df6d000ca91669cea117474cda559e4188600fb2744bafa0f69fd1541c226"
        ),
        "classification": "finite_volume_noncompact",
        "cusps": 6,
        "dimension": 10,
        "edge_extension_histogram": {"2": 426},
        "elliptic_edges": 426,
        "endpoint_type_histogram": {"0,2": 10, "1,1": 132, "2,0": 284},
        "facets": 16,
        "finite_vertices": 70,
        "hyperideal_vertices": 0,
        "id": "R12-D8A1x2-C2b-H10-16",
        "parabolic_neighborhood_audit_count": 14,
        "simple_roots_complete": True,
    }


def test_report_states_positive_result_and_precise_limit() -> None:
    text = " ".join(REPORT.read_text(encoding="utf-8").lower().split())
    assert "toot" in text
    assert "explicit algebro-geometric carrier" in text
    assert "unique nonzero survivor" in text
    assert "effective divisor classes" in text
    assert "not claimed to be the full nef cone" in text
    assert "finite-volume noncompact" in text
    assert "no contradiction" in text
