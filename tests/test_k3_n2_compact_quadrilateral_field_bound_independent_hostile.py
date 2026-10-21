"""Independent hostile tests for the all-degree compact quadrilateral theorem."""

from __future__ import annotations

from scripts.audit_k3_n2_compact_quadrilateral_field_bound_independent_hostile import (
    PINS,
    ROOT,
    build_audit,
    determinant_certificate,
    diagonal_certificate,
    field_packet_certificate,
    norm_gap_certificate,
    quadratic_certificate,
    sha256,
)


def test_frozen_owner_bytes() -> None:
    for relative, expected in PINS.values():
        assert sha256(ROOT / relative) == expected


def test_gram_determinant_is_a_formal_polynomial_identity() -> None:
    record = determinant_certificate()
    assert record["formal_polynomial_identity"] is True
    assert record["determinant_monomial_count"] == 17
    assert record["growth_inequality"] == "(U^2-1)(V^2-1)<4U+4V+4UV+8"


def test_full_trace_degree_and_norm_gap_are_all_degree() -> None:
    record = norm_gap_certificate()
    assert record["degree_upper_bound"] == 6
    assert record["degree_rows"][-1] == {
        "degree": 6,
        "full_trace_rank": 18,
        "strict_shift_lower_bound": "1/4096",
    }
    assert record["norm_gap"] == "A-1>4^(-d)"


def test_degree_six_diagonal_constant_rederives_exactly() -> None:
    record = diagonal_certificate(6, 32_772)
    assert record["V_lower_bound"] == "sqrt(4097)/64"
    assert record["linear_coefficient"] == "16384+256sqrt(4097)"
    assert record["constant_coefficient"] == "32769+256sqrt(4097)"
    assert record["remainder"] == "537034767-8389888sqrt(4097)"
    assert record["remainder_square_difference"] == 17_606_953_271_521
    assert record["both_opposite_squared_ratios_strictly_below"] == 1_074_003_984


def test_complete_packet_not_one_coordinate_generates_general_field() -> None:
    record = field_packet_certificate()
    assert record["simple_cycle_lengths"] == [2, 3, 4]
    assert record["complete_cycle_packet_generates_k_is_hypothesis"] is True
    assert record["single_coordinate_generation_used_for_degree_gt_2"] is False
    assert record["minimal_generator_subpacket_size_upper_bound"] == 2
    assert record["primitive_grid"] == "theta1+c*theta2 for c=0,...,15"
    assert record["Vandermonde_difference_height"] == 590_584_093_515_710_595_072
    assert record["Vandermonde_exponent"] == 30
    assert record["conclusion_is_finite_field_set"] is True


def test_quadratic_specialization_uses_disjoint_conjugate_intervals() -> None:
    diagonal = diagonal_certificate(2, 132)
    assert diagonal["shift_gap"] == "1/16"
    assert diagonal["remainder_square_difference"] == 1_286_881
    record = quadratic_certificate()
    assert record["disjoint_intervals_force_nonrational"] is True
    assert record["one_nonrational_coordinate_generates_only_because_degree_is_two"] is True
    assert record["trace_upper"] == 69_699
    assert record["field_discriminant_upper"] == 4_857_950_597


def test_combined_audit_preserves_exact_scope() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_CONDITIONAL_EXACT_SCOPE"
    assert all(audit["scope"].values())
