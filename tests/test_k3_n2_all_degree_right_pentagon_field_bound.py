"""Regression tests for the all-degree compact K3 right-pentagon theorem."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import pytest

from scripts.certify_k3_n2_all_degree_right_pentagon_field_bound import (
    all_degree_field_record,
    build_certificate,
    closure_record,
    degree_bound_record,
    digest,
    encoded_certificate,
    quadratic_record,
    recurrence_residuals,
    shifted_pentagon,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_all_degree_right_pentagon_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_all_degree_right_pentagon_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_all_degree_right_pentagon_field_bound.json"


def test_all_five_A2_relations_replay_on_exact_grid():
    tested = 0
    for left_denominator in range(1, 11):
        for right_denominator in range(1, 11):
            for left_numerator in range(1, 6 * left_denominator):
                for right_numerator in range(1, 6 * right_denominator):
                    left = Fraction(left_numerator, left_denominator)
                    right = Fraction(right_numerator, right_denominator)
                    if left * right <= 1:
                        continue
                    coordinates = shifted_pentagon(left, right)
                    assert recurrence_residuals(coordinates) == (0, 0, 0, 0, 0)
                    record = closure_record(left, right)
                    assert record["recurrence_residuals"] == ["0"] * 5
                    assert record["principal_Gram_minor_residuals"] == ["0"] * 5
                    tested += 1
    assert tested > 60_000


def test_degree_dependent_norm_and_recurrence_bounds_are_exact():
    expected_upper = {
        2: 288,
        3: 4_224,
        4: 66_048,
        5: 1_050_624,
        6: 16_785_408,
    }
    for degree in range(2, 7):
        record = degree_bound_record(degree)
        reciprocal_floor = 4**degree
        assert record["strict_shifted_coordinate_lower_bound"] == f"1/{reciprocal_floor}"
        assert record["strict_shifted_coordinate_upper_bound"] == expected_upper[degree]
        assert record["strict_squared_pair_coordinate_upper_bound"] == (reciprocal_floor + 1) ** 2
        assert record["strict_pair_cycle_upper_bound"] == 4 * (reciprocal_floor + 1) ** 2
        assert record["strict_length_five_cycle_absolute_upper_bound"] == 32 * (reciprocal_floor + 1) ** 5


def test_corrected_bound_does_not_retain_quadratic_floor_in_degree_six():
    epsilon = Fraction(1, 4**6)
    correct = Fraction(1, epsilon**2) + Fraction(2, epsilon)
    invalid_mixed_floor = 16 + 17 / epsilon
    assert correct == 16_785_408
    assert invalid_mixed_floor == 69_648
    assert correct > invalid_mixed_floor


def test_dense_positive_packet_obeys_degree_six_bound_under_degree_six_floor():
    epsilon = Fraction(1, 4096)
    upper = 16_785_408
    tested = 0
    for left_denominator in range(1, 15):
        for right_denominator in range(1, 15):
            for left_numerator in range(1, 25):
                for right_numerator in range(1, 25):
                    left = Fraction(left_numerator, left_denominator)
                    right = Fraction(right_numerator, right_denominator)
                    if left * right <= 1:
                        continue
                    coordinates = shifted_pentagon(left, right)
                    if min(coordinates) <= epsilon:
                        continue
                    assert max(coordinates) < upper
                    tested += 1
    assert tested > 20_000


def test_cyclic_height_and_vandermonde_constants_recompute():
    record = all_degree_field_record()
    assert record["field_degree_upper_bound"] == 6
    assert record["uniform_shifted_coordinate_strict_upper_bound"] == 16_785_408
    assert record["uniform_squared_pair_coordinate_strict_upper_bound"] == 16_785_409
    assert record["simple_cyclic_generator_conjugate_bound"] == 36_938_546_139_294_728_224
    assert record["selected_generator_count_upper_bound"] == 2
    assert record["bad_primitive_combination_coefficient_count_upper_bound"] == 15
    assert record["Vandermonde_difference_bound"] == 1_182_033_476_457_431_303_168
    assert record["field_discriminant_bound"] == "1182033476457431303168^30"


def test_quadratic_compatibility_retains_frozen_sharp_constant():
    record = quadratic_record()
    assert record["strict_shifted_coordinate_lower_bound"] == "1/16"
    assert record["strict_shifted_coordinate_upper_bound"] == 288
    assert record["pair_cycle_generates_quadratic_field"] is True
    assert record["quadratic_field_discriminant_upper_bound"] == 1_343_277


def test_archive_replays_and_freezes_scope():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    conclusions = payload["conclusions"]
    assert conclusions["K3_full_trace_field_degree_at_most_6"] is True
    assert conclusions["automatic_shift_floor_is_4_to_minus_d"] is True
    assert conclusions["all_degree_every_shift_below_16785408"] is True
    assert conclusions["suggested_69648_bound_is_not_claimed"] is True
    assert conclusions["bounded_mixed_height_premise_required"] is False
    scope = payload["scope"]
    assert scope["returned_field_is_cyclic_Cartan_ground_field_required"] is True
    assert scope["exact_same_compact_right_pentagon_eigen_nef_chamber_required"] is True
    assert scope["complete_saturated_reflection_packet_retained"] is True
    assert scope["reflection_maximality_in_declared_order_host_retained"] is True
    assert scope["reflection_maximality_used_in_bound"] is False
    assert scope["arbitrary_n2_polygon_closed"] is False


def test_report_states_every_hypothesis_constant_and_nonclaim():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    required = [
        "proved for all full-trace complex-K3 degrees",
        "minimal complex projective K3 surface",
        "primitive full trace",
        "exactly a compact right-angled Coxeter pentagon",
        "complete and saturated",
        "reflection-maximal",
        "ground field generated by the cyclic Cartan products",
        r"4^{-d}<x_i<(4^d+1)^2-1",
        r"x_i<16,785,408",
        r"x_i+1<4097^2=16,785,409",
        r"One cannot retain \(1/16\)",
        r"16+17\cdot4^6=69,648",
        r"c_5^2=4^5\prod_{i=0}^4 A_i",
        "Theorem 6.1 (all-degree compact-right-pentagon field theorem)",
        "Hermite--Minkowski proves the finite field set",
        r"does not: - prove that every exact \(n=2\) chamber is a right pentagon",
    ]
    for phrase in required:
        assert phrase in text


def test_no_classification_or_hidden_height_is_imported():
    combined = " ".join(
        (REPORT.read_text(encoding="utf-8") + CERTIFIER.read_text(encoding="utf-8")).split()
    )
    assert "No imposed correspondence height" in combined
    assert "reflective-lattice classification is used" in combined
    assert '"bounded_mixed_height_premise_required": True' not in combined


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        shifted_pentagon(Fraction(1, 2), 2)
    with pytest.raises(ValueError):
        recurrence_residuals((Fraction(1),) * 4)
    with pytest.raises(ValueError):
        degree_bound_record(1)
    with pytest.raises(ValueError):
        degree_bound_record(7)
