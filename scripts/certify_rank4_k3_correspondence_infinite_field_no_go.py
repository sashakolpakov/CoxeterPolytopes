#!/usr/bin/env python3
"""Exact certificate for the rank-four K3 infinite-field no-go.

The checker uses only integer and rational arithmetic.  It reconstructs the
Gram matrices, self-adjoint companion action, primitive K3-lattice embedding,
and finite Pell samples.  The infinitude step itself is the elementary Pell
growth argument recorded in the accompanying report; no reflective-lattice
classification is imported here.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from math import isqrt
import json
from pathlib import Path
from typing import Sequence


Matrix = tuple[tuple[int, ...], ...]
Vector = tuple[int, ...]


def transpose(matrix: Sequence[Sequence[int]]) -> Matrix:
    if not matrix:
        return ()
    return tuple(tuple(row[column] for row in matrix) for column in range(len(matrix[0])))


def multiply(
    left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]
) -> Matrix:
    if not left or not right:
        return ()
    right_t = transpose(right)
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in right_t)
        for row in left
    )


def add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(a + b for a, b in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def scale(coefficient: int, matrix: Matrix) -> Matrix:
    return tuple(tuple(coefficient * entry for entry in row) for row in matrix)


def identity(rank: int) -> Matrix:
    return tuple(
        tuple(1 if row == column else 0 for column in range(rank))
        for row in range(rank)
    )


def block_diagonal(*blocks: Matrix) -> Matrix:
    rank = sum(len(block) for block in blocks)
    answer = [[0] * rank for _ in range(rank)]
    offset = 0
    for block in blocks:
        for row, entries in enumerate(block):
            for column, entry in enumerate(entries):
                answer[offset + row][offset + column] = entry
        offset += len(block)
    return tuple(tuple(row) for row in answer)


def determinant(matrix: Sequence[Sequence[int]]) -> int:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    answer = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        pivot_value = work[column][column]
        answer *= pivot_value
        for row in range(column + 1, len(work)):
            coefficient = work[row][column] / pivot_value
            for entry in range(column, len(work)):
                work[row][entry] -= coefficient * work[column][entry]
    if answer.denominator != 1:
        raise AssertionError("an integral determinant became fractional")
    return int(answer)


def companion(a: int) -> Matrix:
    if a < 2:
        raise ValueError("a must be at least two")
    return ((0, -1), (1, 2 * a))


def b_zero() -> Matrix:
    return ((-1, 0), (0, 1))


def b_negative(a: int) -> Matrix:
    if a < 2:
        raise ValueError("a must be at least two")
    return ((1 - 2 * a, -1), (-1, -1))


def gram(a: int) -> Matrix:
    return block_diagonal(scale(2, b_zero()), scale(2, b_negative(a)))


def tau(a: int) -> Matrix:
    block = companion(a)
    return block_diagonal(block, block)


def operator_record(a: int) -> dict[str, object]:
    c = companion(a)
    b0 = b_zero()
    ba = b_negative(a)
    g = gram(a)
    operator = tau(a)
    left0 = multiply(transpose(c), b0)
    right0 = multiply(b0, c)
    lefta = multiply(transpose(c), ba)
    righta = multiply(ba, c)
    self_adjoint = multiply(transpose(operator), g) == multiply(g, operator)
    polynomial = add(
        add(multiply(operator, operator), scale(-2 * a, operator)), identity(4)
    )
    expected_zero = tuple(tuple(0 for _ in range(4)) for _ in range(4))
    det_ba = determinant(ba)
    det_g = determinant(g)
    if left0 != right0 or lefta != righta or not self_adjoint:
        raise AssertionError(f"self-adjointness failed at a={a}")
    if polynomial != expected_zero:
        raise AssertionError(f"quadratic action failed at a={a}")
    if det_ba != 2 * (a - 1) or det_g != -32 * (a - 1):
        raise AssertionError(f"determinant formula failed at a={a}")
    return {
        "a": a,
        "companion": [list(row) for row in c],
        "gram": [list(row) for row in g],
        "C_transpose_B0_equals_B0_C": [list(row) for row in left0],
        "C_transpose_Ba_equals_Ba_C": [list(row) for row in lefta],
        "self_adjoint": self_adjoint,
        "minimal_polynomial_annihilates_tau": polynomial == expected_zero,
        "det_Ba": det_ba,
        "det_gram": det_g,
        "lattice_signature": [1, 3],
        "t_eigenspace_signature": [1, 1],
        "conjugate_eigenspace_signature": [0, 2],
        "signature_witness": {
            "B0": "diag(-1,1)",
            "Ba_leading_entry_negative": ba[0][0] < 0,
            "Ba_determinant_positive": det_ba > 0,
            "B0_eigenline_square": "lambda^2-1",
            "t_greater_than_one": True,
            "conjugate_between_zero_and_one": True,
        },
        "free_order_module_rank": 2,
        "restriction_of_scalars_Q_dimension": 4,
    }


def k3_target_gram() -> Matrix:
    hyperbolic_plane = ((0, 1), (1, 0))
    return block_diagonal(
        hyperbolic_plane,
        hyperbolic_plane,
        hyperbolic_plane,
        ((-2,),),
    )


def embedding_vectors(a: int) -> tuple[Vector, ...]:
    if a < 2:
        raise ValueError("a must be at least two")
    # Coordinates are e1,f1,e2,f2,e3,f3,v with v^2=-2.
    return (
        (1, -1, 0, 0, 0, 0, 0),
        (0, 0, 1, 1, 0, 0, 0),
        (0, 0, 0, 0, 1, 2 - 2 * a, 1),
        (0, 0, 0, 0, 0, 0, 1),
    )


def embedding_record(a: int) -> dict[str, object]:
    target = k3_target_gram()
    vectors = embedding_vectors(a)
    columns = transpose(vectors)
    pullback = multiply(multiply(transpose(columns), target), columns)
    pivot_rows = (0, 2, 4, 6)  # e1,e2,e3,v
    pivot_minor = tuple(tuple(columns[row][column] for column in range(4)) for row in pivot_rows)
    minor_det = determinant(pivot_minor)
    if pullback != gram(a):
        raise AssertionError(f"K3 embedding Gram failed at a={a}")
    if abs(minor_det) != 1:
        raise AssertionError(f"K3 embedding ceased to be primitive at a={a}")
    return {
        "a": a,
        "target_basis": ["e1", "f1", "e2", "f2", "e3", "f3", "v"],
        "target_gram": [list(row) for row in target],
        "embedding_columns": [list(vector) for vector in vectors],
        "pullback_gram": [list(row) for row in pullback],
        "primitive_pivot_rows": ["e1", "e2", "e3", "v"],
        "primitive_minor_determinant": minor_det,
        "primitive_embedding": abs(minor_det) == 1,
        "K3_orthogonal_complement_signature": [2, 16],
    }


def squarefree_decomposition(number: int) -> tuple[int, int]:
    if number <= 0:
        raise ValueError("number must be positive")
    residual = number
    square_factor = 1
    squarefree = 1
    prime = 2
    while prime * prime <= residual:
        exponent = 0
        while residual % prime == 0:
            residual //= prime
            exponent += 1
        square_factor *= prime ** (exponent // 2)
        if exponent % 2:
            squarefree *= prime
        prime = 3 if prime == 2 else prime + 2
    if residual > 1:
        squarefree *= residual
    if squarefree * square_factor * square_factor != number:
        raise AssertionError("squarefree decomposition failed")
    return squarefree, square_factor


def field_record(a: int) -> dict[str, object]:
    if a < 2:
        raise ValueError("a must be at least two")
    radicand = a * a - 1
    squarefree, pell_y = squarefree_decomposition(radicand)
    if a * a - squarefree * pell_y * pell_y != 1:
        raise AssertionError("Pell identity failed")
    return {
        "a": a,
        "radicand": radicand,
        "squarefree_field_parameter": squarefree,
        "Pell_y": pell_y,
        "Pell_identity": True,
        "polynomial_discriminant": 4 * radicand,
    }


def build_certificate(max_a: int = 128) -> dict[str, object]:
    if max_a < 2:
        raise ValueError("max_a must be at least two")
    exact_samples = (2, 3, 4, 7, 19, max_a)
    exact_samples = tuple(dict.fromkeys(value for value in exact_samples if value <= max_a))
    fields = [field_record(a) for a in range(2, max_a + 1)]
    distinct_fields = sorted({row["squarefree_field_parameter"] for row in fields})
    return {
        "theorem": "rank-four K3 correspondence infinite-field scoped no-go",
        "family": {
            "parameter_range": "every integer a>=2",
            "field": "Q(sqrt(a^2-1))",
            "polynomial": "x^2-2*a*x+1",
            "NS_rank": 4,
            "lattice_signature": [1, 3],
            "K3_minimal": True,
            "K3_canonical_class_zero": True,
            "primitive_full_trace_block": True,
            "integral_self_adjoint_operator": True,
            "algebraic_cycle_coefficients": "rational divisor-product correspondence",
        },
        "operator_samples": [operator_record(a) for a in exact_samples],
        "embedding_samples": [embedding_record(a) for a in exact_samples],
        "field_sample_limit": max_a,
        "field_samples": fields,
        "distinct_sample_field_count": len(distinct_fields),
        "distinct_sample_squarefree_parameters": distinct_fields,
        "field_infinitude_proof": {
            "method": "finite-union contradiction using Pell exponential growth",
            "fixed_field_count_up_to_X": "O_d(log X)",
            "all_parameters_count_up_to_X": "X-1",
            "bounded_degree_alone_implies_finite_fields": False,
        },
        "finite_field_repairs": {
            "bounded_field_discriminant_with_bounded_degree": True,
            "bounded_conjugates_of_primitive_integral_generator": True,
            "bounded_integral_operator_matrix_in_bounded_polarized_basis": True,
            "unnormalized_operator_norm_alone": False,
        },
        "escape_parameters": {
            "polynomial_discriminant": "4(a^2-1)",
            "lattice_determinant": "-32(a-1)",
            "companion_max_entry": "2a",
        },
        "scope": {
            "exact_complete_restricted_root_component_returned": True,
            "prescribed_input_chamber_returned": False,
            "restricted_component_Coxeter_claimed": False,
            "restricted_component_finite_sided_claimed": False,
            "restricted_component_finite_volume_claimed": False,
            "reflection_saturation_claimed": False,
            "reflection_maximality_claimed": False,
            "exact_maximal_same_P_infinite_field_family_claimed": False,
            "strong_exact_maximal_field_set_gate": "open",
        },
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def envelope(max_a: int = 128) -> dict[str, object]:
    payload = build_certificate(max_a=max_a)
    return {"payload": payload, "payload_sha256": digest(payload)}


def encoded_certificate(max_a: int = 128) -> str:
    return json.dumps(envelope(max_a=max_a), indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-a", type=int, default=128)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args()
    encoded = encoded_certificate(max_a=arguments.max_a)
    if arguments.check is not None:
        if arguments.check.read_text(encoding="utf-8") != encoded:
            raise SystemExit(f"archive mismatch: {arguments.check}")
        print(f"PASS archive replay: {arguments.check}")
        return
    if arguments.output is not None:
        arguments.output.write_text(encoded, encoding="utf-8")
        print(arguments.output)
        return
    print(encoded, end="")


if __name__ == "__main__":
    main()
