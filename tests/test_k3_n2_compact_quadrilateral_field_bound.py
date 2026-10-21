"""Regression tests for the compact K3 quadrilateral field theorem."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import pytest

from scripts.certify_k3_n2_compact_quadrilateral_field_bound import (
    build_certificate,
    determinant,
    determinant_formula,
    determinant_record,
    diagonal_bound_record,
    digest,
    encoded_certificate,
    general_field_bound_record,
    gram_matrix,
    permutation_sign,
    quadratic_field_bound_record,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_compact_quadrilateral_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_compact_quadrilateral_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_compact_quadrilateral_field_bound.json"
AUDIT = ROOT / "results/k3_n2_compact_quadrilateral_field_bound_hostile_scope_audit.md"


def test_determinant_expansion_matches_direct_leibniz_on_exact_grid():
    tested = 0
    adjacent_values = (Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4))
    opposite_values = (Fraction(5, 4), Fraction(3, 2), Fraction(2), Fraction(7, 3))
    for a in adjacent_values:
        for b in adjacent_values:
            for c in adjacent_values:
                for d in adjacent_values:
                    for U in opposite_values:
                        for V in opposite_values:
                            direct = determinant(gram_matrix(a, b, c, d, U, V))
                            assert direct == determinant_formula(a, b, c, d, U, V)
                            tested += 1
    assert tested == 4096


def test_coefficient_bounds_hold_on_dense_rational_angle_grid():
    for denominator in range(2, 40):
        for a_numerator in range(denominator):
            for b_numerator in range(denominator):
                a = Fraction(a_numerator, denominator)
                b = Fraction(b_numerator, denominator)
                c = Fraction((3 * a_numerator + 1) % denominator, denominator)
                d = Fraction((5 * b_numerator + 2) % denominator, denominator)
                P = a * b + c * d
                Q = a * d + b * c
                R = a * c + b * d
                K = (
                    -a * a
                    - b * b
                    - c * c
                    - d * d
                    + a * a * c * c
                    + b * b * d * d
                    - 2 * a * b * c * d
                )
                assert 0 <= P < 2
                assert 0 <= Q < 2
                assert 0 <= R < 2
                assert abs(K) < 8


def test_quarter_integral_diagonal_bound_constants_are_exact():
    record = diagonal_bound_record(2)
    assert record["shift_lower_bound"] == "1/16"
    assert record["opposite_magnitude_lower_bound"] == "sqrt(17)/4"
    assert record["linear_coefficient_upper_bound"] == "64+16sqrt(17)"
    assert record["constant_coefficient_upper_bound"] == "129+16sqrt(17)"
    assert record["trial_opposite_magnitude"] == 132
    assert record["trial_quadratic_remainder"] == "8847-2128sqrt(17)"
    assert record["remainder_square_difference"] == 1_286_881
    assert record["opposite_squared_ratio_strict_upper_bound"] == 17_424
    general = diagonal_bound_record(6)
    assert general["shift_lower_bound"] == "1/4096"
    assert general["opposite_magnitude_lower_bound"] == "sqrt(4097)/64"
    assert general["linear_coefficient_upper_bound"] == "16384+256sqrt(4097)"
    assert general["constant_coefficient_upper_bound"] == "32769+256sqrt(4097)"
    assert general["trial_opposite_magnitude"] == 32_772
    assert general["trial_quadratic_remainder"] == "537034767-8389888sqrt(4097)"
    assert general["remainder_square_difference"] == 17_606_953_271_521


def test_trial_132_is_beyond_the_positive_quadratic_root_exactly():
    assert 8847 > 0
    assert 8847**2 - 17 * 2128**2 == 1_286_881 > 0
    # The vertex lies below 132 because (64+16*sqrt(17))/2 < 72.
    assert (64 + 16 * 5) // 2 < 132


def test_field_trace_norm_and_discriminant_constants_recompute():
    record = quadratic_field_bound_record()
    assert record["pair_cycle_distinguished_strict_upper_bound"] == 69_696
    assert record["pair_cycle_conjugate_strict_upper_bound"] == 4
    assert record["integral_trace_upper_bound"] == 69_699
    assert record["integral_norm_upper_bound"] == 278_783
    assert record["field_discriminant_upper_bound"] == 4_857_950_597
    assert record["opposite_pair_cycle_generates_quadratic_field"] is True
    assert 69_699**2 - 4 == 4_857_950_597
    general = general_field_bound_record()
    assert general["field_degree_upper_bound"] == 6
    assert general["opposite_magnitude_strict_upper_bound"] == 32_772
    assert general["opposite_squared_ratio_strict_upper_bound"] == 1_074_003_984
    assert general["simple_cyclic_coordinate_conjugate_bound"] == 18_455_752_922_365_956_096
    assert general["selected_generator_count_upper_bound"] == 2
    assert general["bad_primitive_combination_coefficient_count_upper_bound"] == 15
    assert general["Vandermonde_difference_bound"] == 590_584_093_515_710_595_072
    assert general["field_discriminant_bound"] == "590584093515710595072^30"


def test_archive_replays_and_freezes_exact_scope():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    conclusions = payload["conclusions"]
    assert conclusions["K3_full_trace_field_degree_at_most_6"] is True
    assert conclusions["all_degree_both_opposite_magnitudes_below_32772"] is True
    assert conclusions["bounded_simple_cycles_generate_ground_field"] is True
    assert conclusions["quadratic_field_discriminant_at_most_4857950597"] is True
    assert conclusions["bounded_height_premise_required"] is False
    assert conclusions["angle_order_classification_required"] is False
    scope = payload["scope"]
    assert scope["primitive_exact_full_trace_Picard_required"] is True
    assert scope["returned_field_is_cyclic_Cartan_ground_field_required"] is True
    assert scope["exact_same_compact_Coxeter_quadrilateral_eigen_nef_chamber_required"] is True
    assert scope["complete_saturated_reflection_packet_retained"] is True
    assert scope["reflection_maximality_in_declared_order_host_retained"] is True
    assert scope["reflection_maximality_used_in_bound"] is False
    assert scope["noncompact_quadrilateral_closed"] is False
    assert scope["arbitrary_n2_polygon_closed"] is False


def test_report_states_normalization_proof_and_nonclaims():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    required = [
        "proved for every full-trace complex-K3 field degree",
        "minimal complex projective K3 surface",
        "primitive full trace",
        "exactly a compact Coxeter quadrilateral",
        "complete and saturated",
        "reflection-maximal",
        "Only the squared ratios",
        "no square root is silently adjoined",
        r"r,s>4^{-d}\ge4^{-6}={1\over4096}",
        r"0=rs+K-2PU-2QV-2RUV",
        r"|K|<8",
        r"rs<4U+4V+4UV+8",
        r"U^2 <(16384+256\sqrt{4097})U +(32769+256\sqrt{4097})",
        r"537034767^2-4097\cdot8389888^2 =17,606,953,271,521>0",
        "Theorem 6.1 (all-degree compact-quadrilateral field theorem)",
        r"\bigl(590,584,093,515,710,595,072\bigr)^{30}",
        r"|\operatorname {Disc}(k)| \le69,699^2-4 =4,857,950,597",
        "The determinant alone allows the geometric escape",
        r"The theorem does not: - prove that every exact \(n=2\) chamber is a quadrilateral",
    ]
    for phrase in required:
        assert phrase in text


def test_no_angle_or_reflective_lattice_classification_imported():
    combined = " ".join(
        (REPORT.read_text(encoding="utf-8") + CERTIFIER.read_text(encoding="utf-8")).split()
    )
    assert "No imposed height, power relation, angle classification" in combined
    assert "uses_reflective_lattice_classification\": True" not in combined


def test_hostile_audit_keeps_determinant_and_scope_literal():
    text = " ".join(AUDIT.read_text(encoding="utf-8").split())
    required = [
        "PASS, all-degree compact-quadrilateral scope",
        "single determinant expansion",
        r"(U^2-1)(V^2-1)+K-2PU-2QV-2RUV",
        r"r,s>4^{-d}",
        r"U,V<32772",
        "bounded simple cyclic packet",
        r"590,584,093,515,710,595,072",
        "Maximality is retained but unused",
        "does not cover noncompact quadrilaterals",
    ]
    for phrase in required:
        assert phrase in text


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        determinant(((Fraction(1), Fraction(0)),))
    assert permutation_sign((0, 1, 2, 3)) == 1
    assert permutation_sign((1, 0, 2, 3)) == -1
    record = determinant_record(0, 0, 0, 0, 2, 3)
    assert record["expansion_identity"] is True
