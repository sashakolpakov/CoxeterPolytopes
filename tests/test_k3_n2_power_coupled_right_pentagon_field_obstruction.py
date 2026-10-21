"""Regression tests for the power-coupled rank-six K3 pentagon obstruction."""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
import json
from pathlib import Path

import pytest

from scripts.certify_k3_n2_power_coupled_right_pentagon_field_obstruction import (
    build_certificate,
    definite_place_record,
    digest,
    encoded_certificate,
    exact_height_inequality,
    norm_divisor_records,
    pentagon_closure_record,
    power_coupling_record,
    proper_divisors,
    quadratic_candidate_records,
    squarefree_part,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_power_coupled_right_pentagon_field_obstruction.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_power_coupled_right_pentagon_field_obstruction.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_power_coupled_right_pentagon_field_obstruction.json"
AUDIT = ROOT / "results/k3_n2_power_coupled_right_pentagon_field_obstruction_hostile_scope_audit.md"


def test_all_five_rank_three_principal_minors_vanish_on_exact_grid():
    for numerator_a in range(3, 15):
        for numerator_c in range(3, 15):
            A = Fraction(numerator_a, 2)
            C = Fraction(numerator_c, 2)
            if A * C - A - C <= 0:
                continue
            record = pentagon_closure_record(A, C)
            assert record["all_five_principal_minors_vanish"] is True
            assert record["principal_four_by_four_minors"] == ["0"] * 5
            assert Fraction(record["B"]) > 1
            assert Fraction(record["D"]) > 1
            assert Fraction(record["E"]) > 1


def test_power_coupling_gives_the_fixed_divisor_identity_exactly():
    for A in (Fraction(3, 2), Fraction(7, 4), Fraction(2), Fraction(5, 2), Fraction(3)):
        for exponent in range(2, 18):
            u = A**exponent - A ** (exponent - 1) - 1
            if u <= 0:
                continue
            record = power_coupling_record(A, exponent)
            assert record["Delta_equals_A_times_u"] is True
            assert Fraction(record["E"]) == 1 + 1 / u
            assert Fraction(record["D"]) == 1 + A ** (exponent - 1) / u
            assert Fraction(record["u_times_v"]) == 4


def test_definite_place_intervals_are_strict_on_dense_rational_grid():
    for denominator in range(2, 30):
        for numerator in range(1, denominator):
            conjugate = Fraction(numerator, denominator)
            for exponent in (2, 3, 4, 9, 20):
                record = definite_place_record(conjugate, exponent)
                assert -2 < Fraction(record["u_conjugate"]) < -1
                assert -4 < Fraction(record["v_conjugate"]) < -2


def test_norm_gate_has_exactly_the_divisors_of_sixteen():
    records = norm_divisor_records()
    assert proper_divisors(16) == (1, 2, 4, 8, 16)
    assert [record["negative_norm_magnitude"] for record in records] == [1, 2, 4, 8, 16]
    for record in records:
        n = record["negative_norm_magnitude"]
        assert Fraction(record["distinguished_u_strict_upper_bound"]) == Fraction(16, n)
        assert record["uniform_u_strict_upper_bound"] == 16


def test_quadratic_trace_norm_enumeration_is_the_claimed_five_rows():
    records = quadratic_candidate_records()
    assert [
        (
            record["trace"],
            record["norm"],
            record["discriminant"],
            record["field_squarefree_radicand"],
        )
        for record in records
    ] == [
        (3, 1, 5, 5),
        (4, 1, 12, 3),
        (4, 2, 8, 2),
        (5, 2, 17, 17),
        (5, 3, 13, 13),
    ]


def test_exact_radical_comparison_agrees_with_high_precision_control():
    getcontext().prec = 80
    envelope = (Decimal(1) + Decimal(69).sqrt()) / 2
    for trace in range(2, 10):
        for norm in range(1, 15):
            discriminant = trace * trace - 4 * norm
            if discriminant <= 0:
                continue
            root = (Decimal(trace) + Decimal(discriminant).sqrt()) / 2
            assert exact_height_inequality(trace, norm) is (root < envelope)


def test_squarefree_field_labels_are_reduced():
    assert squarefree_part(5) == 5
    assert squarefree_part(8) == 2
    assert squarefree_part(12) == 3
    assert squarefree_part(13) == 13
    assert squarefree_part(17) == 17


def test_archive_replays_and_freezes_every_scope_gate():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    conclusions = payload["conclusions"]
    assert conclusions["power_coupling_fixed_divisor_identity"] == "u*v=4"
    assert conclusions["possible_field_radicands"] == [2, 3, 5, 13, 17]
    assert conclusions["finite_field_set_inside_power_coupled_ansatz"] is True
    scope = payload["scope"]
    assert scope["minimal_complex_projective_K3_required"] is True
    assert scope["primitive_exact_full_trace_Picard_rank_six_required"] is True
    assert scope["exact_same_compact_right_pentagon_eigen_nef_chamber_required"] is True
    assert scope["complete_saturated_integral_reflection_packet_required"] is True
    assert scope["reflection_maximality_in_declared_order_host_retained"] is True
    assert scope["additional_power_coupling_C_equals_A_to_q_required"] is True
    assert scope["additional_A_algebraic_integer_field_generator_required"] is True
    assert scope["reflection_maximality_used_to_derive_height_bound"] is False
    assert scope["arbitrary_independent_A_C_coupling_closed"] is False
    assert scope["existence_of_a_carrier_in_the_ansatz_asserted"] is False


def test_report_visually_separates_carrier_hypotheses_from_extra_ansatz():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    required = [
        "## 1. The exact carrier package retained",
        "## 2. The additional power-coupled ansatz",
        "minimal complex projective K3 surface",
        "primitive full trace lattice",
        "exactly a compact right-angled pentagon",
        "complete and saturated",
        "reflection-maximal among the integral",
        r"A\in\mathcal O,\qquad k=\mathbf Q(A),\qquad C=A^q",
        r"The algebraic integrality of \(A\) is an actual extra hypothesis",
        r"\Delta=AC-A-C",
        r"v:=4E-4={4\over u}\in\mathcal O",
        r"N(u)N(v)=16,\qquad n\mid16",
        r"A_0<{1+\sqrt{69}\over2}",
        "Theorem 6.1 (power-coupled right-pentagon obstruction)",
        r"\mathbf Q(\sqrt{17})",
        "Reflection maximality is not used in deriving it",
        r"without the relation \(C=A^q\)",
        "The arbitrary nonorthogonal coupling problem remains open",
    ]
    for phrase in required:
        assert phrase in text


def test_no_forbidden_classification_or_general_n2_claim_appears():
    combined = REPORT.read_text(encoding="utf-8") + CERTIFIER.read_text(encoding="utf-8")
    assert "No reflective-lattice classification is used" in combined
    assert "arbitrary_independent_A_C_coupling_closed\": True" not in combined
    assert "existence_of_a_carrier_in_the_ansatz_asserted\": True" not in combined


def test_hostile_audit_keeps_the_general_gap_open():
    text = " ".join(AUDIT.read_text(encoding="utf-8").split())
    required = [
        "PASS, scoped",
        "The exact K3/maximal package is retained",
        "The power relation and algebraic integrality are extra",
        r"B=(A-1)(C-1)",
        r"u\,v=4",
        r"v_1\in(-4,-2)",
        r"n\mid16",
        r"2\le T\le5",
        r"1\le N\le4",
        r"\mathbf Q(\sqrt2)",
        r"independent \(A,C\)",
        r"does not solve the general \(n=2\) automatic-height problem",
    ]
    for phrase in required:
        assert phrase in text


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        pentagon_closure_record(1, 3)
    with pytest.raises(ValueError):
        pentagon_closure_record(Fraction(3, 2), Fraction(3, 2))
    with pytest.raises(ValueError):
        power_coupling_record(2, 1)
    with pytest.raises(ValueError):
        power_coupling_record(Fraction(6, 5), 2)
    with pytest.raises(ValueError):
        definite_place_record(0, 2)
    with pytest.raises(ValueError):
        definite_place_record(1, 2)
    with pytest.raises(ValueError):
        proper_divisors(0)
    with pytest.raises(ValueError):
        squarefree_part(0)
