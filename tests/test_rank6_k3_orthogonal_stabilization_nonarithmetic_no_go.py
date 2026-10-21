"""Regression tests for the rank-six orthogonal-stabilization no-go."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import pytest

from scripts.certify_rank6_k3_orthogonal_stabilization_nonarithmetic_no_go import (
    build_certificate,
    determinant,
    digest,
    embedding_record,
    five_roots,
    mixed_midpoint_coefficients,
    mixed_midpoint_sign,
    multiply,
    obstruction_record,
    q,
    sector_roots,
    sign_b_minus_ta,
    trace_lattice_record,
    transpose,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.md"
CERTIFIER = (
    ROOT / "scripts/certify_rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.py"
)
ARCHIVE = (
    ROOT
    / "data/enumerations/rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.json"
)


def test_trace_lattice_is_unimodularly_diagonal_and_self_adjoint():
    for a in range(2, 301):
        row = trace_lattice_record(a)
        assert abs(row["unimodular_change_determinant"]) == 1
        assert row["tau_self_adjoint"] is True
        assert row["tau_minimal_polynomial"] == "x^2-2*a*x+1"
        assert row["field_module_rank"] == 3
        assert row["real_signatures"] == [[1, 2], [3, 0]]
        d = a - 1
        assert row["diagonal_gram"] == [
            [2, 0, 0, 0, 0, 0],
            [0, -2, 0, 0, 0, 0],
            [0, 0, -2, 0, 0, 0],
            [0, 0, 0, -4 * d, 0, 0],
            [0, 0, 0, 0, -2, 0],
            [0, 0, 0, 0, 0, -4 * d],
        ]


def test_primitive_k3_embedding_reconstructs_the_gram():
    for a in range(2, 201):
        row = embedding_record(a)
        assert abs(row["primitive_minor_determinant"]) == 1
        assert row["primitive_embedding"] is True
        assert row["orthogonal_complement_signature"] == [2, 14]
        vectors = tuple(tuple(vector) for vector in row["embedding_vectors"])
        target = (
            (0, 1, 0, 0, 0, 0, 0, 0, 0),
            (1, 0, 0, 0, 0, 0, 0, 0, 0),
            (0, 0, 0, 1, 0, 0, 0, 0, 0),
            (0, 0, 1, 0, 0, 0, 0, 0, 0),
            (0, 0, 0, 0, 0, 1, 0, 0, 0),
            (0, 0, 0, 0, 1, 0, 0, 0, 0),
            (0, 0, 0, 0, 0, 0, -2, 0, 0),
            (0, 0, 0, 0, 0, 0, 0, -2, 0),
            (0, 0, 0, 0, 0, 0, 0, 0, -2),
        )
        columns = transpose(vectors)
        assert multiply(multiply(transpose(columns), target), columns) == tuple(
            tuple(entry for entry in gram_row) for gram_row in row["diagonal_gram"]
        )


def test_five_visible_vectors_are_roots_for_every_parameter():
    roots = five_roots()
    assert roots[-1] == (1, 0, 1, 0, 1, 0)
    for d in range(1, 1001):
        assert [q(root, root, d) for root in roots] == [1] * 5


def test_visible_pentagon_vertices_are_strictly_inside_the_klein_disk():
    # Put x=R^2-1.  The only nontrivial vertex test is
    # 1+x^2 < 1+x, which follows exactly from 0<x<1.
    for a in range(2, 2001):
        from math import sqrt

        t = a + sqrt(a * a - 1)
        radius_squared = t / (t - 1)
        x = radius_squared - 1
        assert 0 < x < 1
        assert 1 + x * x < radius_squared


def test_mixed_midpoint_has_one_equality_on_a_large_exact_root_shell():
    gamma = (1, 0, 1, 0, 1, 0)
    for a in range(2, 23):
        d = a - 1
        equalities = []
        for root in sector_roots(d, 36):
            coefficient_a, coefficient_b = mixed_midpoint_coefficients(root)
            assert coefficient_a >= 0
            sign = mixed_midpoint_sign(root, a)
            assert sign <= 0
            if sign == 0:
                equalities.append(root)
        assert equalities == [gamma]


def test_height_one_exposure_cases_are_exhaustive_without_floats():
    gamma = (1, 0, 1, 0, 1, 0)
    for d in range(1, 501):
        roots = list(sector_roots(d, 1))
        equalities = []
        for root in roots:
            coefficient_a, coefficient_b = mixed_midpoint_coefficients(root)
            assert coefficient_a >= 0
            sign = sign_b_minus_ta(coefficient_a, coefficient_b, d + 1)
            assert sign <= 0
            if sign == 0:
                equalities.append(root)
        assert equalities == [gamma]


def test_global_cauchy_gap_used_for_heights_at_least_two():
    for d in range(1, 1001):
        coefficient = 8 * d * d + 12 * d + 6
        for height in (2, 3, 7, 31, 127):
            assert coefficient * (height * height + 1) < (
                4 * (2 * d + 1) ** 2 * height * height
            )


def test_cyclic_product_norm_is_never_integral():
    for d in range(1, 10001):
        row = obstruction_record(d + 1)
        denominator = 4 * d - 1
        assert row["unavoidable_denominator"] == denominator
        norm = Fraction(row["cyclic_product_field_norm"])
        assert norm == Fraction(32 * d, denominator)
        assert norm == 8 + Fraction(8, denominator)
        assert norm.denominator != 1
        assert row["cyclic_product_norm_is_integer"] is False
        assert row["any_order_lattice_reflection_packet_possible"] is False


def test_cyclic_product_norm_recomputes_in_the_quadratic_field():
    for d in range(1, 1001):
        a = d + 1
        denominator = 4 * d - 1
        c0 = Fraction(4 * (2 * d - 1), denominator)
        c1 = Fraction(4, denominator)
        # If x=c0+c1*t and t*t'=1, t+t'=2a, then
        # N(x)=c0^2+2a*c0*c1+c1^2.
        norm = c0 * c0 + 2 * a * c0 * c1 + c1 * c1
        assert norm == Fraction(32 * d, denominator)


def test_every_higher_orthogonal_packet_contains_the_same_pair_obstruction():
    # In q=-Y^2+X^2+sum(U_i^2+2*d*Z_i^2), the vector with
    # Y=U_i=U_j=1 has norm one independently of r and of the chosen pair.
    for rank_over_field in range(3, 11):
        slope_count = rank_over_field - 1
        for first in range(slope_count):
            for second in range(first + 1, slope_count):
                coordinates = [1, 0] + [0, 0] * slope_count
                coordinates[2 + 2 * first] = 1
                coordinates[2 + 2 * second] = 1
                for d in (1, 2, 7, 31):
                    norm = -coordinates[0] ** 2 + coordinates[1] ** 2
                    norm += sum(
                        coordinates[2 + 2 * index] ** 2
                        + 2 * d * coordinates[3 + 2 * index] ** 2
                        for index in range(slope_count)
                    )
                    assert norm == 1


def test_archive_replays_and_keeps_every_scope_gate_literal():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate(max_a=64)
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    conclusions = payload["conclusions"]
    assert conclusions["primitive_exact_Picard_K3_carrier_available"] is True
    assert conclusions["complete_restricted_component_compact"] is True
    assert conclusions["mixed_closing_root_is_exposed"] is True
    assert conclusions["mixed_cyclic_product_algebraic_integral"] is False
    assert conclusions["integral_reflection_saturation_possible"] is False
    assert conclusions["overorder_or_overlattice_repairs_obstruction"] is False
    assert conclusions["all_orthogonal_stabilizations_r_at_least_2_obstructed"] is True
    scope = payload["scope"]
    assert scope["orthogonal_trace_stabilizations_only"] is True
    assert scope["all_nonorthogonal_rank_six_couplings_excluded"] is False
    assert scope["primitive_K3_embedding_beyond_rank_six_claimed"] is False
    assert scope["compactness_beyond_rank_six_claimed"] is False
    assert scope["n_at_least_two_global_field_theorem_resolved"] is False


def test_report_states_the_positive_packet_failure_and_surviving_problem():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "primitive rank-six K3 trace carrier",
        "compact right-angled pentagon",
        "mixed K3 root required to close the new direction",
        "nonintegral cyclic Cartan product",
        "any commensurable full restriction-of-scalars lattice",
        "applies even when the commensurable \\(\\mathbf Z\\)-lattice is not",
        "cannot be a product or any orthogonal stabilization",
        "every orthogonal stabilization fails",
        "no uniform K3 embedding or full-chamber compactness is asserted",
        "No impossibility theorem for all such nonorthogonal couplings is claimed",
        "exact surviving \\(n\\ge2\\) field-height problem",
    ):
        assert phrase in text


def test_no_classification_or_global_n_two_claim_is_smuggled_in():
    combined = REPORT.read_text(encoding="utf-8") + CERTIFIER.read_text(
        encoding="utf-8"
    )
    for forbidden in (
        "Nikulin classification",
        "Vinberg classification",
        "ABSW theorem",
        "all nonorthogonal couplings are impossible",
    ):
        assert forbidden not in combined


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        trace_lattice_record(1)
    with pytest.raises(ValueError):
        embedding_record(1)
    with pytest.raises(ValueError):
        obstruction_record(1)
    with pytest.raises(ValueError):
        build_certificate(max_a=1)
    with pytest.raises(ValueError):
        sign_b_minus_ta(-1, 0, 2)
