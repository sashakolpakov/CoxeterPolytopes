from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_a15_7_c01_special_k3_symplectic_lift import (
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT / "data" / "enumerations" / "a15_7_c01_special_k3_symplectic_lift.json"
)
REPORT = ROOT / "results" / "a15_7_c01_special_k3_symplectic_lift.md"
LEDGER = ROOT / "results" / "higher_dimensional_polytope_lineage_ledger.md"
SOURCE_REPORT = ROOT / "results" / "alexeev_remaining_fixed_section_positive_examples.md"


def _archived() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_special_k3_lift_replays_exactly() -> None:
    archived = _archived()
    assert build_certificate() == archived
    assert archived["all_checks_pass"]
    assert archived["case_id"] == "A15_7_1_C01"
    assert archived["certificate_payload_sha256"] == (
        "98304fee498980acae78915903155f5310778dc29cb718d842ac8e41cb2213e9"
    )


def test_discriminant_graph_and_nikulin_lattice_are_exact() -> None:
    archived = _archived()
    glue = archived["transcendental_glue"]
    assert glue["lattice"] == "U(2)+U(2)+<-2>^3"
    assert glue["anti_invariant_lattice"] == "D4(-2)"
    assert glue["actions_equivariant"]
    assert glue["all_128_discriminant_classes_checked"]
    assert len(glue["quadratic_form_table"]) == 128
    assert glue["source_q2_histogram_mod_4"] == {
        "0": 28,
        "1": 36,
        "2": 36,
        "3": 28,
    }

    action = archived["k3_lattice_action"]
    assert action["even_unimodular_signature"] == [3, 19]
    assert action["action_order"] == 2
    assert action["action_characteristic_polynomial"] == "(x-1)^14*(x+1)^8"
    assert action["fixed_lattice_rank_signature_determinant"] == [
        14,
        [3, 11],
        -256,
    ]
    assert action["anti_invariant_lattice"] == "E8(-2)"
    assert action["anti_invariant_rank"] == 8
    assert action["anti_invariant_determinant"] == 256
    assert action["anti_invariant_has_no_minus_two_vectors"]


def test_period_picard_and_fixed_divisors_are_exact() -> None:
    archived = _archived()
    period = archived["period_and_picard_lattice"]
    assert period["period_eigenvalue"] == "1"
    assert period["period_square"] == "0"
    assert period["period_positive_norm"] == "4"
    assert period["rational_T_orthogonal"] == "exactly the anti-invariant D4(-2)"
    assert period["picard_rank"] == 19
    assert period["picard_signature"] == [1, 18]
    assert period["picard_determinant"] == 32
    assert period["picard_smith_diagonal_absolute"] == [*[1] * 16, 2, 4, 4]
    assert period["ambient_N_is_contained"]
    assert period["added_algebraic_rank"] == 4
    assert period["all_added_algebraic_directions_are_anti_invariant"]
    assert period["fixed_NS_is_exactly_archived_N_fixed_lattice"]
    assert period["fixed_NS_rank"] == 11
    assert "genuine symplectic Nikulin involution" in archived[
        "torelli_and_chamber"
    ]["global_torelli_conclusion"]


def test_literal_vinberg_and_duplicate_lineage_are_not_overclaimed() -> None:
    archived = _archived()
    assert archived["hyperbolic_fixed_section"] == {
        "bad_edge_count": 0,
        "bad_parabolic_neighborhood_count": 0,
        "compact": False,
        "criterion": (
            "Vinberg edge-endpoint criterion plus every irreducible affine "
            "component's full zero-neighborhood"
        ),
        "dimension": 10,
        "elliptic_edge_count": 6221,
        "facet_count": 27,
        "finite_vertex_count": 895,
        "finite_volume": True,
        "ideal_vertex_count": 67,
    }
    lineage = archived["lineage"]
    assert lineage["same_signed_normalized_gram_as"] == "D16-D4-C7-H10-27"
    assert len(lineage["facet_permutation_one_based"]) == 27
    assert not lineage["abstract_polytope_is_new"]
    assert not lineage["dimension_is_new"]
    assert not lineage["stable_lattice_or_maximal_reflection_group_is_new"]


def test_reports_and_ledger_state_the_geometric_and_nef_scope() -> None:
    report = REPORT.read_text(encoding="utf-8")
    ledger = LEDGER.read_text(encoding="utf-8")
    source_report = SOURCE_REPORT.read_text(encoding="utf-8")
    normalized_report = " ".join(report.lower().split())
    normalized_source_report = " ".join(source_report.split())
    assert "Picard rank 19" in report
    assert "symplectic Nikulin involution" in report
    assert "not a new abstract" in normalized_report
    assert "not asserted to be the literal nef cone" in normalized_report
    assert "A15_7-C2-H10-27" in ledger
    assert "special Picard-rank-19 K3" in ledger
    assert "A15_7_1_C01" in source_report
    assert "genuine symplectic Nikulin involution" in normalized_source_report
