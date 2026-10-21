"""Regression tests for the bounded chordless-cycle K3 field gate."""

from __future__ import annotations

from fractions import Fraction
import json
from math import comb
from pathlib import Path

import pytest

from scripts.certify_k3_chordless_cycle_height_field_bound import (
    affine_cycle_record,
    build_certificate,
    chord_decomposition_record,
    cyclic_tridiagonal,
    digest,
    encoded_certificate,
    k3_cycle_bound_record,
    matrix_rank,
    trace_induction_record,
    tridiagonal_rank_record,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_chordless_cycle_height_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_chordless_cycle_height_field_bound.py"
ARCHIVE = (
    ROOT / "data/enumerations/k3_chordless_cycle_height_field_bound.json"
)
AUDIT = (
    ROOT / "results/k3_chordless_cycle_height_field_bound_hostile_scope_audit.md"
)


def test_cyclic_tridiagonal_samples_have_nullity_at_most_two():
    for length in range(3, 41):
        record = tridiagonal_rank_record(length)
        assert record["sample_rank"] >= length - 2
        assert record["universal_nullity_bound"] == 2
        assert record["universal_kernel_injection"] == "ker(A)->(x_1,x_2)"


def test_reflection_trace_expansion_isolates_one_full_cycle():
    for length in range(1, 21):
        record = trace_induction_record(length)
        assert record["nonempty_subsequence_term_count"] == 2**length - 1
        assert record["full_cycle_coefficient"] == (-1) ** length
        assert record["full_cycle_occurs_once"] is True
        assert record["inductive_integrality_isolation"] is True
        assert record["terms_by_cycle_length"] == {
            str(size): comb(length, size)
            for size in range(1, length + 1)
        }


def test_nullity_bound_survives_many_exact_coefficient_variants():
    for length in range(3, 15):
        for shift in range(-3, 8):
            matrix = [list(row) for row in cyclic_tridiagonal(length, shift)]
            for index in range(length):
                matrix[index][index] += (index * index + 3 * shift) % 7
                successor = (index + 1) % length
                value = 1 + ((index + shift) % 4)
                matrix[index][successor] = value
                matrix[successor][index] = value
            assert matrix_rank(matrix) >= length - 2


def test_every_chord_decomposition_has_the_exact_edge_exponents():
    records = []
    for length in range(4, 25):
        for left in range(length):
            for step in range(2, length - 1):
                right = (left + step) % length
                # Rotate endpoints so the certifier's canonical ordering applies.
                if left < right:
                    record = chord_decomposition_record(length, left, right)
                else:
                    record = chord_decomposition_record(length, right, left)
                assert record["edge_exponent_identity"] is True
                assert record["identity"] == (
                    "c(C1)c(C2)=c(C)c(chord-pair)"
                )
                records.append(record)
    assert len(records) > 1000


def test_affine_signed_cycle_has_the_small_or_dominant_dichotomy():
    for length in range(2, 11):
        unit = 2**length
        others = (-unit + 1, -1, 0, 1, unit)
        for beta0 in (0, unit - 1, unit, unit + 1, 3 * unit):
            record = affine_cycle_record(length, beta0, others)
            conjugates = [
                Fraction(value)
                for value in record["affine_coordinate_conjugates"]
            ]
            assert min(conjugates) > 0
            if beta0 > unit:
                assert record["large_at_distinguished_place"] is True
                assert conjugates[0] >= max(conjugates[1:])
            else:
                assert record["large_at_distinguished_place"] is False
                assert max(conjugates) <= 2 * unit


def test_n2_constant_and_discriminant_exponent_are_exact():
    record = k3_cycle_bound_record(
        hyperbolic_dimension=2,
        correspondence_height_bound=256,
        square_lower_bound=2,
    )
    assert record["field_module_rank"] == 3
    assert record["field_degree_bound"] == 6
    assert record["chordless_cycle_length_bound"] == 5
    assert record["cycle_Cauchy_absolute_bound"] == 32
    assert record["affine_small_coordinate_bound"] == 64
    assert record["selected_generator_count_bound"] == 2
    assert record["primitive_combination_coefficient_bound"] == 16
    assert record["individual_affine_coordinate_bound"] == "128"
    assert record["primitive_generator_conjugate_bound"] == "4096"
    assert int(record["field_discriminant_bound"]) == 4096**30


def test_k3_dimension_ladder_closes_exactly_at_nine_for_nonrational_fields():
    for dimension in range(2, 10):
        record = k3_cycle_bound_record(dimension, 256)
        assert record["field_degree_bound"] == 20 // (dimension + 1)
        assert record["chordless_cycle_length_bound"] == dimension + 3
        assert record["selected_generator_count_bound"] <= 2
        assert record["finite_polynomial_box_size"] > 0
    with pytest.raises(ValueError):
        k3_cycle_bound_record(10, 256)


def test_archive_replays_and_freezes_the_scope():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    conclusions = payload["conclusions"]
    assert (
        conclusions[
            "ground_field_cycles_reduce_to_pairs_and_chordless_cycles"
        ]
        is True
    )
    assert conclusions["cyclic_integrality_from_reflection_trace_induction"]
    assert conclusions["individual_Cartan_entries_in_order_required"] is False
    assert conclusions["chordless_cycle_length_at_most_n_plus_3"] is True
    assert conclusions["pair_generation_required"] is False
    assert (
        conclusions["bounded_large_cycle_mixed_heights_give_finite_field_set"]
        is True
    )
    scope = payload["scope"]
    assert scope["complex_projective_K3_required"] is True
    assert scope["exact_full_Picard_trace_K3_required"] is True
    assert (
        scope[
            "finite_projective_order_module_or_explicit_trace_certificate_required"
        ]
        is True
    )
    assert scope["returned_reflections_preserve_order_module_required"] is True
    assert scope["bounded_mixed_degree_derived_from_maximality"] is False
    assert scope["uses_reflective_lattice_classification"] is False


def test_report_states_every_gate_and_nonclaim():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    required = [
        "Theorem 2.1 (the \\(m+2\\) cycle theorem)",
        "Lemma 1.1 (trace induction gives cyclic integrality)",
        "finite projective rank-\\(m\\) module",
        r"q(x,y)=\operatorname {Tr}_{k/\mathbf Q}",
        "repeated and non-induced proper subsequences are covered",
        r"\operatorname {tr}_k(s_{i_1}\cdots s_{i_\ell})",
        r"\ell\le m+2",
        r"\ker A\longrightarrow\mathbf R^2",
        "Lemma 3.1 (small-or-dominant dichotomy)",
        r"\theta_C=\beta_C+U_C",
        "Theorem 4.1 (bounded cyclic correspondence field gate)",
        "complex projective",
        "uniform constants \\(s_0,C>0\\)",
        r"D&=\left\lfloor{20\over n+1}\right\rfloor",
        r"\left(32\max\{64,C/s_0\}\right)^{30}",
        "finite field set directly",
        r"H\cdot\tau_t(H)=32a-10",
        "pair-generation restriction",
        "Reflection maximality has not yet been proved",
    ]
    for phrase in required:
        assert phrase in text


def test_no_reflection_finiteness_classification_is_imported():
    combined = REPORT.read_text(encoding="utf-8") + CERTIFIER.read_text(
        encoding="utf-8"
    )
    assert "classification of reflective lattices" in combined
    assert "uses_reflective_lattice_classification\": True" not in combined


def test_hostile_audit_keeps_the_height_premise_literal():
    text = " ".join(AUDIT.read_text(encoding="utf-8").split())
    required = [
        "PASS, conditionally",
        "Individual Cartan entries need not lie in the returned order",
        "finite projective returned order-module",
        r"c(C_1)c(C_2)=c(C)\,a_{ij}a_{ji}",
        r"\ell-2\le m",
        r"|\sigma(c(C))|\le2^\ell",
        r"\theta_C=(-1)^\ell c(C)+2^\ell",
        r"d(n+1)=\rho(X)\le20",
        r"(32R_0)^{30}",
        "finite field set, not only a bounded field degree",
        r"H\cdot\tau_{\theta_C}(H)\le C",
        "does not rule out one whose canonical mixed degree escapes",
    ]
    for phrase in required:
        assert phrase in text


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        cyclic_tridiagonal(2)
    with pytest.raises(ValueError):
        trace_induction_record(0)
    with pytest.raises(ValueError):
        chord_decomposition_record(3, 0, 2)
    with pytest.raises(ValueError):
        chord_decomposition_record(5, 0, 1)
    with pytest.raises(ValueError):
        affine_cycle_record(1, 2, (1,))
    with pytest.raises(ValueError):
        affine_cycle_record(3, -1, (1,))
    with pytest.raises(ValueError):
        affine_cycle_record(3, 1, (9,))
    with pytest.raises(ValueError):
        affine_cycle_record(3, 1, (-8,))
    with pytest.raises(ValueError):
        k3_cycle_bound_record(1, 10)
