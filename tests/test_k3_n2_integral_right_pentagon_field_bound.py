"""Regression tests for the general K3 right-pentagon field theorem."""

from __future__ import annotations

from fractions import Fraction
import json
from math import sqrt
from pathlib import Path

import pytest

from scripts.certify_k3_n2_integral_right_pentagon_field_bound import (
    build_certificate,
    closure_record,
    digest,
    encoded_certificate,
    general_field_bound_record,
    integral_coordinate_candidate_packet,
    phi_add,
    phi_multiply,
    phi_scale,
    poison_packet_record,
    quarter_integral_bound_record,
    recurrence_residuals,
    shifted_pentagon,
    squarefree_part,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_integral_right_pentagon_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_integral_right_pentagon_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_integral_right_pentagon_field_bound.json"
AUDIT = ROOT / "results/k3_n2_integral_right_pentagon_field_bound_hostile_scope_audit.md"


def test_all_five_gram_minors_and_A2_relations_vanish_exactly():
    tested = 0
    for left_denominator in range(1, 11):
        for right_denominator in range(1, 11):
            for left_numerator in range(1, 5 * left_denominator):
                for right_numerator in range(1, 5 * right_denominator):
                    left = Fraction(left_numerator, left_denominator)
                    right = Fraction(right_numerator, right_denominator)
                    if left * right <= 1:
                        continue
                    record = closure_record(left, right)
                    assert record["A2_recurrence_residuals"] == ["0"] * 5
                    assert record["principal_Gram_minor_residuals"] == ["0"] * 5
                    assert Fraction(record["moving_denominator"]) == left * right - 1
                    tested += 1
    assert tested > 40_000


def test_quarter_integral_norm_and_recurrence_constants_are_exact():
    record = quarter_integral_bound_record()
    assert record["integral_scale"] == 4
    assert Fraction(record["strict_shifted_coordinate_lower_bound"]) == Fraction(1, 16)
    assert Fraction(record["strict_shifted_coordinate_upper_bound"]) == 288
    assert Fraction(record["moving_denominator_upper_bound"]) == 288
    assert record["proof_identity"] == "x[i+3]=(z+1)/(x[i]*z-1)"


def test_dense_positive_A2_grid_obeys_universal_upper_bound_when_lower_gate_holds():
    passed = 0
    epsilon = Fraction(1, 16)
    for left_denominator in range(1, 16):
        for right_denominator in range(1, 16):
            for left_numerator in range(1, 20):
                for right_numerator in range(1, 20):
                    left = Fraction(left_numerator, left_denominator)
                    right = Fraction(right_numerator, right_denominator)
                    if left * right <= 1:
                        continue
                    coordinates = shifted_pentagon(left, right)
                    if min(coordinates) <= epsilon:
                        continue
                    assert max(coordinates) < 288
                    assert recurrence_residuals(coordinates) == (0, 0, 0, 0, 0)
                    passed += 1
    assert passed > 10_000


def test_general_pair_generator_and_discriminant_constants_recompute():
    record = general_field_bound_record()
    assert record["pair_cycle_distinguished_strict_upper_bound"] == 1156
    assert record["pair_cycle_conjugate_strict_upper_bound"] == 4
    assert record["integral_trace_upper_bound"] == 1159
    assert record["integral_norm_upper_bound"] == 4623
    assert record["polynomial_and_field_discriminant_upper_bound"] == 1_343_277
    assert record["each_pair_cycle_generates_quadratic_field"] is True
    assert 1159**2 - 4 == 1_343_277


def test_integral_coordinate_corollary_has_528_rows_and_336_candidate_fields():
    packet = integral_coordinate_candidate_packet()
    assert packet["trace_range"] == [3, 34]
    assert packet["norm_range"] == [1, 33]
    assert packet["candidate_trace_norm_row_count"] == 528
    assert packet["candidate_squarefree_radicand_count"] == 336
    assert len(packet["candidate_squarefree_radicands"]) == 336
    assert max(row["polynomial_discriminant"] for row in packet["candidate_rows"]) == 1152
    assert packet["candidate_not_realization_count"] is True


def test_squarefree_reduction_for_candidate_fields():
    assert squarefree_part(8) == 2
    assert squarefree_part(12) == 3
    assert squarefree_part(45) == 5
    assert squarefree_part(1152) == 2


def test_poison_packet_is_exactly_quarter_integral_but_not_integral():
    packet = poison_packet_record()
    assert packet["recurrence_residuals"] == [["0", "0"]] * 5
    assert packet["four_times_each_coordinate_integral"] == [True] * 5
    assert packet["each_coordinate_integral"] == [False, False, False, True, True]
    assert packet["first_three_coordinates_not_integral"] is True
    assert packet["directed_acute_five_cycle_in_1_phi_basis"] == ["-408", "-664"]
    assert packet["asserts_exact_K3_or_maximal_carrier"] is False


def test_poison_packet_real_embeddings_have_the_claimed_strict_signs():
    phi0 = (1 + sqrt(5)) / 2
    phi1 = (1 - sqrt(5)) / 2
    coordinates = [
        (Fraction(-1, 2), Fraction(1, 2)),
        (Fraction(-1, 4), Fraction(1, 4)),
        (Fraction(1, 2), Fraction(2)),
        (Fraction(8), Fraction(14)),
        (Fraction(2), Fraction(4)),
    ]
    distinguished = [float(a) + float(b) * phi0 for a, b in coordinates]
    conjugate = [float(a) + float(b) * phi1 for a, b in coordinates]
    assert min(distinguished) > 0
    assert max(distinguished) < 288
    assert all(-1 < value < 0 for value in conjugate)
    one = (Fraction(1), Fraction(0))
    for index in range(5):
        product = phi_multiply(coordinates[(index - 1) % 5], coordinates[(index + 1) % 5])
        assert product == phi_add(coordinates[index], one)
        assert all(coefficient.denominator == 1 for coefficient in phi_scale(4, coordinates[index]))


def test_archive_replays_and_freezes_general_and_sharp_scopes():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    conclusions = payload["conclusions"]
    assert conclusions["moving_denominator_AC_minus_A_minus_C_below_288"] is True
    assert conclusions["general_field_discriminant_at_most_1343277"] is True
    assert conclusions["power_coupling_required"] is False
    assert conclusions["bounded_mixed_height_premise_required"] is False
    assert conclusions["integral_A_C_corollary_candidate_squareclasses"] == 336
    scope = payload["scope"]
    assert scope["exact_same_compact_right_pentagon_eigen_nef_chamber_required"] is True
    assert scope["complete_saturated_reflection_packet_retained"] is True
    assert scope["reflection_maximality_in_declared_order_host_retained"] is True
    assert scope["reflection_maximality_used_in_bound"] is False
    assert scope["A_and_C_algebraic_integers_required_for_general_theorem"] is False
    assert scope["A_and_C_algebraic_integers_required_for_1152_corollary"] is True
    assert scope["arbitrary_n2_chamber_closed"] is False


def test_report_states_every_hypothesis_constant_and_nonclaim():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    required = [
        "proved general compact-right-pentagon field theorem",
        "minimal complex projective K3 surface",
        "primitive full trace",
        "exactly a compact right-angled pentagon",
        "complete and saturated",
        "reflection-maximal",
        r"0<A_1,B_1,C_1,D_1,E_1<1",
        r"All five Gram identities are the \(A_2\) pentagon",
        r"x_{i-1}x_{i+1}=x_i+1",
        r"y_i:=4x_i\in\mathcal O",
        r"(x_i)_0>{1\over16}",
        r"x_i<288",
        r"0<T^2-4N\le1159^2-4=1,343,277",
        "Theorem 6.1 (general integral right-pentagon field bound)",
        r"|\operatorname {Disc}(k)|\le34^2-4=1152",
        "528 possible trace--norm rows and 336 distinct squarefree radicands",
        "These are candidates, not asserted realized K3 fields",
        "Poison packet against division by four",
        "not a claimed K3 surface",
        r"does not: - prove that every exact \(n=2\) chamber is a right pentagon",
    ]
    for phrase in required:
        assert phrase in text


def test_no_pair_generation_gap_or_classification_import():
    combined = " ".join(
        (
            REPORT.read_text(encoding="utf-8")
            + CERTIFIER.read_text(encoding="utf-8")
        ).split()
    )
    assert "a rational element would have" in combined
    assert "No classification of reflective lattices is imported" in combined
    assert "uses_reflective_lattice_classification\": True" not in combined


def test_hostile_audit_keeps_the_quarter_integral_scope_exact():
    text = " ".join(AUDIT.read_text(encoding="utf-8").split())
    required = [
        "PASS, general right-pentagon scope",
        r"The automatic statement is \(4x_i\in\mathcal O\)",
        r"x_{i-1}x_{i+1}=x_i+1",
        r"1/16<x_i<288",
        "every pair cycle is nonrational",
        r"field discriminant bound \(1,343,277\)",
        "336 candidate squareclasses, not 336 realized fields",
        "poison packet",
        "Maximality is retained but unused",
        r"does not settle arbitrary \(n=2\) incidence",
    ]
    for phrase in required:
        assert phrase in text


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        shifted_pentagon(0, 2)
    with pytest.raises(ValueError):
        shifted_pentagon(Fraction(1, 2), 2)
    with pytest.raises(ValueError):
        recurrence_residuals((Fraction(1),) * 4)
    with pytest.raises(ValueError):
        quarter_integral_bound_record(0)
    with pytest.raises(ValueError):
        squarefree_part(0)
