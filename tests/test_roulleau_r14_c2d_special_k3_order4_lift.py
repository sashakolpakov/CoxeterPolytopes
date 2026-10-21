from __future__ import annotations

import json
from pathlib import Path

import pytest

from scripts.certify_roulleau_r14_c2d_special_k3_order4_lift import (
    BH_ORDER4_ROW_COUNT,
    BH_ORDER4_SHA256,
    OUTPUT,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results" / "roulleau_r14_c2d_special_k3_order4_lift.md"
BH_LOCAL = Path(
    "/private/tmp/K3Groups_extract/K3Groups/results/purely_ns/order04.txt"
)


def _archived() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_positive_order_four_k3_lift_replays_exactly() -> None:
    archived = _archived()
    assert build_certificate() == archived
    assert archived["all_checks_pass"]
    assert archived["certificate_payload_sha256"] == (
        "053a15e0d9c8200f30cdecb2308f5d60d682de4e141647dc41e978ec45082b3f"
    )
    assert "There exists a complex projective K3 surface" in archived["theorem"]
    assert "automorphism sigma of order four" in archived["theorem"]


def test_picard_pair_and_rootless_kahler_gate_are_exact() -> None:
    archived = _archived()
    pair = archived["archived_picard_pair"]
    assert pair["rank_signature_determinant"] == [14, [1, 13], -16]
    assert pair["canonical_action"].endswith(
        "(12)(34) on the four A1(-1) summands"
    )
    assert pair["fixed_rank_signature_determinant"] == [12, [1, 11], -16]
    assert pair["fixed_smith"] == [*[1] * 10, 4, 4]
    assert pair["anti_invariant_lattice"] == "<-4>+<-4>"
    assert pair["anti_invariant_rootless"]

    kahler = archived["kahler_and_torelli"]
    assert kahler["all_27_walls_permuted"]
    assert kahler["invariant_ample_square"] == 528
    assert kahler["minimum_wall_pairing"] == "2"
    assert kahler["strictly_inside_source_chamber"]
    assert not kahler["anti_invariant_NS_has_minus_two_vectors"]
    assert "genuine K3 automorphism" in kahler["strong_torelli_conclusion"]


def test_equivariant_unimodular_glue_and_exact_period_close() -> None:
    archived = _archived()
    glue = archived["equivariant_k3_glue"]
    assert glue["transcendental_lattice"] == "U+U+A1(-1)^4"
    assert glue["transcendental_action_square"] == "-identity"
    assert glue["glued_lattice"] == "even unimodular II_(3,19)"
    assert glue["glued_action_order"] == 4
    assert glue["square_fixed_lattice_is_exactly_N"]
    assert glue["global_fixed_lattice_is_exactly_archived_fixed_lattice"]
    assert len(glue["quadratic_form_table"]) == 16
    assert all(
        (row["4q_N_mod_8"] + row["4q_T_mod_8"]) % 8 == 0
        for row in glue["quadratic_form_table"]
    )

    period = archived["explicit_period_and_picard"]
    assert period["minimal_polynomial_irreducible_over_Q_i"]
    assert period["i_eigenbasis_hermitian_gram"] == [
        ["4", "0", "0", "0"],
        ["0", "-4", "0", "0"],
        ["0", "0", "-4", "0"],
        ["0", "0", "0", "-4"],
    ]
    assert period["period_square"] == "0"
    assert period["period_eigenvalue"] == "i"
    assert period["period_positive_norm_strict_lower_bound"] == "244949/62500"
    assert period["no_nonzero_rational_T_orthogonal"]
    assert period["picard_lattice_is_exactly_N"]
    assert period["picard_rank"] == 14


def test_literal_vinberg_finite_volume_h11_19_handoff_is_frozen() -> None:
    chamber = _archived()["hyperbolic_fixed_section"]
    assert (chamber["dimension"], chamber["facet_count"]) == (11, 19)
    assert chamber["finite_volume"] and not chamber["compact"]
    assert chamber["elliptic_edge_count"] == 795
    assert chamber["finite_vertex_count"] == 114
    assert chamber["ideal_vertex_count"] == 13
    assert chamber["bad_edge_count"] == 0
    assert chamber["bad_parabolic_neighborhood_count"] == 0
    assert chamber["edge_endpoint_type_histogram"] == {
        "0,2": 31,
        "1,1": 274,
        "2,0": 490,
    }


def test_classification_lineage_and_elliptic_chart_are_stated_precisely() -> None:
    archived = _archived()
    bh = archived["brandhorst_hofmann_cross_identification"]
    assert bh["member_sha256"] == BH_ORDER4_SHA256
    assert bh["complete_row_count"] == BH_ORDER4_ROW_COUNT == 79
    assert bh["cohomological_character_r_l_m"] == [12, 2, 4]
    assert [row["id"] for row in bh["all_rows_with_that_character"]] == [
        "0.4.3.4",
        "0.4.3.12",
    ]
    assert bh["unique_row_with_fixed_determinant_minus_16"] == "0.4.3.4"
    assert bh["target_NS_pair_exactly_U_E8_A1x4_with_double_transposition"]

    elliptic = archived["elliptic_algebro_geometric_chart"]
    assert elliptic["equation_pullback"] == "sigma^*(F)=-F"
    assert elliptic["square"].endswith("elliptic inversion")
    assert elliptic["generic_fiber_configuration"] == [
        "I4* at infinity",
        "4 I2",
        "6 I1",
    ]
    assert "not by point counting" in elliptic["scope_caution"]

    report = REPORT.read_text(encoding="utf-8")
    assert "TOOT" in report
    assert "genuine order-four K3" in report and "automorphism" in report
    assert "finite-volume and noncompact" in report
    assert "not a new" in report and "abstract Coxeter polytope" in report


@pytest.mark.skipif(not BH_LOCAL.exists(), reason="local published BH table absent")
def test_optional_full_brandhorst_hofmann_replay_is_exact() -> None:
    replayed = build_certificate(BH_LOCAL)
    bh = replayed["brandhorst_hofmann_cross_identification"]
    assert bh["external_full_replay_performed"]
    replay = bh["external_replay"]
    assert replay["member_sha256"] == BH_ORDER4_SHA256
    assert replay["row_count"] == 79
    assert replay["target_pair_exactly_integrally_conjugate_to_R14_C2D"]
    assert [row["id"] for row in replay["matching_character_rows"]] == [
        "0.4.3.4",
        "0.4.3.12",
    ]
