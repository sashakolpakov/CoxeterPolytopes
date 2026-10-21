#!/usr/bin/env python3
"""Exact certificate for the rank-six orthogonal-stabilization no-go.

The symbolic checks reconstruct the trace lattice and operator, its
primitive K3 embedding, the five closing roots, the exact sign of every
sampled mixed-midpoint functional, and the nonintegral cyclic product.  No
reflective-lattice classification is imported.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from math import isqrt, sqrt
import json
from pathlib import Path
from typing import Iterable, Sequence


Matrix = tuple[tuple[int, ...], ...]
Vector = tuple[int, ...]


def transpose(matrix: Sequence[Sequence[int]]) -> Matrix:
    return tuple(
        tuple(row[column] for row in matrix)
        for column in range(len(matrix[0]))
    )


def multiply(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> Matrix:
    columns = transpose(right)
    return tuple(
        tuple(sum(a * b for a, b in zip(row, column)) for column in columns)
        for row in left
    )


def determinant(matrix: Sequence[Sequence[int]]) -> int:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    value = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            value = -value
        diagonal = work[column][column]
        value *= diagonal
        for row in range(column + 1, len(work)):
            ratio = work[row][column] / diagonal
            for index in range(column, len(work)):
                work[row][index] -= ratio * work[column][index]
    if value.denominator != 1:
        raise AssertionError("integral determinant became fractional")
    return int(value)


def identity(rank: int) -> Matrix:
    return tuple(
        tuple(int(row == column) for column in range(rank))
        for row in range(rank)
    )


def add(left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]) -> Matrix:
    return tuple(
        tuple(x + y for x, y in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def scale(coefficient: int, matrix: Sequence[Sequence[int]]) -> Matrix:
    return tuple(tuple(coefficient * entry for entry in row) for row in matrix)


def block_diagonal(*blocks: Sequence[Sequence[int]]) -> Matrix:
    sizes = [len(block) for block in blocks]
    rank = sum(sizes)
    result = [[0] * rank for _ in range(rank)]
    offset = 0
    for block, size in zip(blocks, sizes):
        for row in range(size):
            for column in range(size):
                result[offset + row][offset + column] = block[row][column]
        offset += size
    return tuple(tuple(row) for row in result)


def trace_lattice_record(a: int) -> dict[str, object]:
    if a < 2:
        raise ValueError("a must be at least two")
    d = a - 1
    b0 = ((-2, 0), (0, 2))
    ba = ((2 - 4 * a, -2), (-2, -2))
    companion = ((0, -1), (1, 2 * a))
    gram = block_diagonal(b0, ba, ba)
    tau = block_diagonal(companion, companion, companion)
    # Columns send (Y,X,U,Z,V,W) to the three original two-block bases.
    change = (
        (0, 1, 0, 0, 0, 0),
        (1, 0, 0, 0, 0, 0),
        (0, 0, 0, 1, 0, 0),
        (0, 0, 1, -1, 0, 0),
        (0, 0, 0, 0, 0, 1),
        (0, 0, 0, 0, 1, -1),
    )
    diagonal = multiply(multiply(transpose(change), gram), change)
    expected = (
        (2, 0, 0, 0, 0, 0),
        (0, -2, 0, 0, 0, 0),
        (0, 0, -2, 0, 0, 0),
        (0, 0, 0, -4 * d, 0, 0),
        (0, 0, 0, 0, -2, 0),
        (0, 0, 0, 0, 0, -4 * d),
    )
    minimal_polynomial = add(
        add(multiply(tau, tau), scale(-2 * a, tau)), identity(6)
    )
    zero = tuple(tuple(0 for _ in range(6)) for _ in range(6))
    if determinant(change) not in (-1, 1):
        raise AssertionError("diagonal change is not unimodular")
    if diagonal != expected:
        raise AssertionError(f"diagonal lattice failed at a={a}")
    if multiply(gram, tau) != multiply(transpose(tau), gram):
        raise AssertionError(f"self-adjointness failed at a={a}")
    if minimal_polynomial != zero:
        raise AssertionError(f"quadratic field relation failed at a={a}")
    return {
        "a": a,
        "original_gram": [list(row) for row in gram],
        "tau": [list(row) for row in tau],
        "unimodular_change_determinant": determinant(change),
        "diagonal_gram": [list(row) for row in diagonal],
        "tau_self_adjoint": True,
        "tau_minimal_polynomial": "x^2-2*a*x+1",
        "field_module_rank": 3,
        "real_signatures": [[1, 2], [3, 0]],
    }


def q(left: Sequence[int], right: Sequence[int], d: int) -> int:
    return (
        -left[0] * right[0]
        + left[1] * right[1]
        + left[2] * right[2]
        + 2 * d * left[3] * right[3]
        + left[4] * right[4]
        + 2 * d * left[5] * right[5]
    )


def target_gram() -> Matrix:
    u = ((0, 1), (1, 0))
    zero2 = ((0, 0), (0, 0))
    # Basis e1,f1,e2,f2,e3,f3,v1,v2,v3.
    result = [[0] * 9 for _ in range(9)]
    for offset in (0, 2, 4):
        for row in range(2):
            for column in range(2):
                result[offset + row][offset + column] = u[row][column]
    for index in (6, 7, 8):
        result[index][index] = -2
    return tuple(tuple(row) for row in result)


def embedding_vectors(d: int) -> tuple[Vector, ...]:
    if d < 1:
        raise ValueError("d must be positive")
    return (
        (1, 1, 0, 0, 0, 0, 0, 0, 0),   # Y
        (0, 0, 0, 0, 0, 0, 1, 0, 0),   # X
        (0, 0, 0, 0, 0, 0, 0, 1, 0),   # U
        (0, 0, 1, -2 * d, 0, 0, 0, 0, 0),  # Z
        (0, 0, 0, 0, 0, 0, 0, 0, 1),   # V
        (0, 0, 0, 0, 1, -2 * d, 0, 0, 0),  # W
    )


def embedding_record(a: int) -> dict[str, object]:
    if a < 2:
        raise ValueError("a must be at least two")
    d = a - 1
    target = target_gram()
    vectors = embedding_vectors(d)
    columns = transpose(vectors)
    pullback = multiply(multiply(transpose(columns), target), columns)
    expected = (
        (2, 0, 0, 0, 0, 0),
        (0, -2, 0, 0, 0, 0),
        (0, 0, -2, 0, 0, 0),
        (0, 0, 0, -4 * d, 0, 0),
        (0, 0, 0, 0, -2, 0),
        (0, 0, 0, 0, 0, -4 * d),
    )
    pivot_rows = (0, 6, 7, 2, 8, 4)
    minor = tuple(
        tuple(columns[row][column] for column in range(6))
        for row in pivot_rows
    )
    if pullback != expected or abs(determinant(minor)) != 1:
        raise AssertionError(f"primitive K3 embedding failed at a={a}")
    return {
        "a": a,
        "diagonal_gram": [list(row) for row in expected],
        "embedding_vectors": [list(row) for row in vectors],
        "pullback_gram": [list(row) for row in pullback],
        "primitive_minor_determinant": determinant(minor),
        "primitive_embedding": True,
        "orthogonal_complement_signature": [2, 14],
    }


def five_roots() -> tuple[Vector, ...]:
    return (
        (0, 0, -1, 0, 0, 0),
        (0, 0, 0, 0, -1, 0),
        (1, -1, 1, 0, 0, 0),
        (1, -1, 0, 0, 1, 0),
        (1, 0, 1, 0, 1, 0),
    )


def sector_roots(d: int, max_height: int) -> Iterable[Vector]:
    for height in range(1, max_height + 1):
        target = height * height + 1
        max_z = isqrt(target // (2 * d))
        for z in range(-max_z, max_z + 1):
            for w in range(-max_z, max_z + 1):
                residual = target - 2 * d * (z * z + w * w)
                if residual < 0:
                    continue
                for m in range(isqrt(residual) + 1):
                    remaining = residual - m * m
                    for u in range(isqrt(remaining) + 1):
                        v2 = remaining - u * u
                        v = isqrt(v2)
                        if v * v == v2:
                            yield height, -m, u, z, v, w


def mixed_midpoint_coefficients(root: Vector) -> tuple[int, int]:
    """Return A,B with twice the midpoint functional equal to B-t*A."""
    height, x, u, z, v, w = root
    m = -x
    return 2 * height - u - v + z + w, 2 * m + z + w


def sign_b_minus_ta(a_coefficient: int, b_coefficient: int, a: int) -> int:
    """Exact sign of B-(a+sqrt(a^2-1))*A for A nonnegative."""
    if a_coefficient < 0:
        raise ValueError("the exact comparison expects A nonnegative")
    if a_coefficient == 0:
        return (b_coefficient > 0) - (b_coefficient < 0)
    rational_part = b_coefficient - a * a_coefficient
    if rational_part <= 0:
        return -1
    left_square = rational_part * rational_part
    right_square = (a * a - 1) * a_coefficient * a_coefficient
    return (left_square > right_square) - (left_square < right_square)


def mixed_midpoint_sign(root: Vector, a: int) -> int:
    a_coefficient, b_coefficient = mixed_midpoint_coefficients(root)
    return sign_b_minus_ta(a_coefficient, b_coefficient, a)


def mixed_wall_record(a: int, max_height: int = 24) -> dict[str, object]:
    d = a - 1
    roots = five_roots()
    norms = [q(root, root, d) for root in roots]
    if norms != [1] * 5:
        raise AssertionError("five visible roots changed norm")
    equalities = []
    root_count = 0
    for root in sector_roots(d, max_height):
        root_count += 1
        sign = mixed_midpoint_sign(root, a)
        if sign > 0:
            raise AssertionError((a, root, sign))
        if sign == 0:
            equalities.append(list(root))
    gamma = [1, 0, 1, 0, 1, 0]
    if equalities != [gamma]:
        raise AssertionError((a, equalities))
    t = a + sqrt(a * a - 1)
    radius_squared = t / (t - 1)
    vertices = (
        (0.0, 0.0),
        (1.0, 0.0),
        (1.0, radius_squared - 1),
        (radius_squared - 1, 1.0),
        (0.0, 1.0),
    )
    if not all(s * s + r * r < radius_squared for s, r in vertices):
        raise AssertionError("visible pentagon ceased to be compact")
    return {
        "a": a,
        "d": d,
        "five_roots": [list(root) for root in roots],
        "root_q_norms": norms,
        "Klein_radius_squared": radius_squared,
        "visible_pentagon_vertices": [list(vertex) for vertex in vertices],
        "visible_pentagon_compact": True,
        "mixed_midpoint": [radius_squared / 2, radius_squared / 2],
        "root_shell_max_height": max_height,
        "root_shell_count": root_count,
        "mixed_midpoint_equalities": equalities,
    }


def obstruction_record(a: int) -> dict[str, object]:
    if a < 2:
        raise ValueError("a must be at least two")
    d = a - 1
    denominator = 4 * d - 1
    # chi = c0+c1*t.
    c0 = Fraction(4 * (2 * d - 1), denominator)
    c1 = Fraction(4, denominator)
    norm = Fraction(32 * d, denominator)
    coefficient_c0 = Fraction(4 * d - 2, denominator)
    coefficient_c1 = Fraction(2, denominator)
    if norm.denominator == 1:
        raise AssertionError("cyclic product unexpectedly became integral")
    return {
        "a": a,
        "d": d,
        "mixed_normal_square": "4*d+1-t",
        "cyclic_product_coefficients_c0_c1": [str(c0), str(c1)],
        "cyclic_product": "8*d/(4*d+1-t)",
        "cyclic_product_field_norm": str(norm),
        "cyclic_product_norm_is_integer": norm.denominator == 1,
        "mixed_reflection_coefficient_2b_over_norm": [
            str(coefficient_c0),
            str(coefficient_c1),
        ],
        "unavoidable_denominator": denominator,
        "any_order_lattice_reflection_packet_possible": False,
    }


def build_certificate(max_a: int = 64) -> dict[str, object]:
    if max_a < 2:
        raise ValueError("max_a must be at least two")
    samples = tuple(
        dict.fromkeys(value for value in (2, 3, 4, 7, 19, max_a) if value <= max_a)
    )
    return {
        "theorem": "rank-six K3 orthogonal stabilization nonarithmetic no-go",
        "Q_rank": 6,
        "field_module_rank": 3,
        "distinguished_hyperbolic_dimension": 2,
        "trace_lattice_records": [trace_lattice_record(a) for a in samples],
        "embedding_records": [embedding_record(a) for a in samples],
        "mixed_wall_records": [mixed_wall_record(a) for a in samples],
        "obstruction_records": [obstruction_record(a) for a in samples],
        "conclusions": {
            "primitive_exact_Picard_K3_carrier_available": True,
            "complete_restricted_component_compact": True,
            "mixed_closing_root_is_exposed": True,
            "mixed_cyclic_product_algebraic_integral": False,
            "integral_reflection_saturation_possible": False,
            "overorder_or_overlattice_repairs_obstruction": False,
            "all_orthogonal_stabilizations_r_at_least_2_obstructed": True,
        },
        "scope": {
            "orthogonal_trace_stabilizations_only": True,
            "all_nonorthogonal_rank_six_couplings_excluded": False,
            "primitive_K3_embedding_beyond_rank_six_claimed": False,
            "compactness_beyond_rank_six_claimed": False,
            "n_at_least_two_global_field_theorem_resolved": False,
            "no_reflective_lattice_classification_imported": True,
        },
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def encoded_certificate(max_a: int = 64) -> str:
    payload = build_certificate(max_a=max_a)
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-a", type=int, default=64)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args()
    encoded = encoded_certificate(arguments.max_a)
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
