#!/usr/bin/env python3
"""Exact replay for the noncompact K3 cusp-divisibility host gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from math import gcd
from pathlib import Path
from typing import Sequence


Matrix = list[list[int]]
Vector = list[int]


def determinant(matrix: Matrix) -> int:
    """Bareiss determinant over the integers."""

    if not matrix:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    size = len(work)
    for column in range(size - 1):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        value = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = work[row][target] * value
                numerator -= work[row][column] * work[column][target]
                work[row][target] = numerator // previous
        previous = value
    return sign * work[-1][-1]


def transpose(matrix: Matrix) -> Matrix:
    return [list(row) for row in zip(*matrix)]


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return [
        [sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0]))]
        for i in range(len(left))
    ]


def matvec(matrix: Matrix, vector: Vector) -> Vector:
    return [sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix]


def bilinear(left: Vector, form: Matrix, right: Vector) -> int:
    return sum(
        left[i] * form[i][j] * right[j]
        for i in range(len(left))
        for j in range(len(right))
    )


def identity(size: int) -> Matrix:
    return [[int(row == column) for column in range(size)] for row in range(size)]


def reflection_matrix(form: Matrix, root: Vector) -> Matrix:
    if bilinear(root, form, root) != -2:
        raise ValueError("this certifier uses square-minus-two roots")
    covector = matvec(form, root)
    # s_r(x)=x+(x,r)r; vectors are columns.
    return [
        [int(row == column) + root[row] * covector[column] for column in range(len(root))]
        for row in range(len(root))
    ]


def gcd_of_maximal_minors(matrix: Matrix) -> int:
    """Return the gcd of full-column-rank minors."""

    rows = len(matrix)
    columns = len(matrix[0])
    result = 0
    for selected in itertools.combinations(range(rows), columns):
        minor = [[matrix[row][column] for column in range(columns)] for row in selected]
        result = gcd(result, abs(determinant(minor)))
    return result


def family_record(index: int) -> dict[str, object]:
    if index < 2:
        raise ValueError("the acute three-root display is used only for d>=2")

    d = index
    host = [[0, d, 0], [d, 0, 0], [0, 0, -2]]
    cusp = [1, 0, 0]
    alpha = [0, 0, 1]
    beta = [1, 0, -1]
    transverse_root = [0, 1, -1]
    roots = [alpha, beta, transverse_root]
    root_columns = transpose(roots)
    root_gram = matmul(transpose(root_columns), matmul(host, root_columns))

    expected_root_gram = [
        [-2, 2, 2],
        [2, -2, d - 2],
        [2, d - 2, -2],
    ]
    if root_gram != expected_root_gram:
        raise AssertionError("root Gram formula changed")
    if determinant(root_columns) != 1:
        raise AssertionError("the three displayed roots stopped being a host basis")
    if determinant(host) != 2 * d * d:
        raise AssertionError("host discriminant formula changed")
    if determinant(root_gram) != determinant(host):
        raise AssertionError("root and host discriminants differ")
    if [alpha[i] + beta[i] for i in range(3)] != cusp:
        raise AssertionError("the affine null relation changed")
    if bilinear(cusp, host, cusp) != 0:
        raise AssertionError("the cusp is not isotropic")
    if bilinear(alpha, host, beta) != 2:
        raise AssertionError("the fixed affine A1 packet changed")
    if any(bilinear(cusp, host, root) for root in (alpha, beta)):
        raise AssertionError("an affine root misses the cusp")
    cusp_pairings = matvec(host, cusp)
    cusp_divisibility = gcd(*(abs(value) for value in cusp_pairings))
    if cusp_divisibility != d:
        raise AssertionError("cusp divisibility formula changed")
    if bilinear(cusp, host, transverse_root) != d:
        raise AssertionError("transverse root degree formula changed")

    unit = identity(3)
    for root in roots:
        if bilinear(root, host, root) != -2:
            raise AssertionError("a displayed physical root changed norm")
        reflection = reflection_matrix(host, root)
        if matmul(reflection, reflection) != unit:
            raise AssertionError("a root reflection is not involutive")
        if matmul(transpose(reflection), matmul(host, reflection)) != host:
            raise AssertionError("a root reflection does not preserve the host")

    # Primitive embedding into U+U+<-2> inside the K3 lattice:
    # u=e1, v=d*f1+e2, a=a.  A unit maximal minor proves saturation.
    ambient = [
        [0, 1, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 0, -2],
    ]
    embedding = [
        [1, 0, 0],
        [0, d, 0],
        [0, 1, 0],
        [0, 0, 0],
        [0, 0, 1],
    ]
    pullback = matmul(transpose(embedding), matmul(ambient, embedding))
    if pullback != host:
        raise AssertionError("the explicit K3-lattice pullback changed")
    minor_gcd = gcd_of_maximal_minors(embedding)
    if minor_gcd != 1:
        raise AssertionError("the displayed embedding is not primitive")

    return {
        "d": d,
        "host_gram": host,
        "signature": [1, 2],
        "absolute_discriminant": 2 * d * d,
        "primitive_cusp": cusp,
        "cusp_divisibility": cusp_divisibility,
        "fixed_affine_A1_gram": [[-2, 2], [2, -2]],
        "affine_relation": "alpha+beta=u",
        "transverse_root_square": -2,
        "transverse_root_cusp_degree": d,
        "root_basis_gram": root_gram,
        "root_basis_change_determinant": 1,
        "all_three_root_reflections_integral": True,
        "primitive_K3_embedding_minor_gcd": minor_gcd,
    }


def build_certificate() -> dict[str, object]:
    sample_indices = [2, 3, 4, 5, 7, 11, 19, 31]
    return {
        "verdict": "cusp-local data do not bound the noncompact K3 integral host",
        "family": [family_record(index) for index in sample_indices],
        "symbolic_family": {
            "host": "M_d=U(d)+<-2>, d>=2",
            "signature": "(1,2)",
            "absolute_discriminant": "2*d^2",
            "primitive_cusp_divisibility": "div(u)=d",
            "fixed_incident_packet": "alpha^2=beta^2=-2, alpha.beta=2, alpha+beta=u",
            "transverse_square_minus_two_root": "c=v-alpha, c^2=-2, u.c=d",
            "primitive_K3_embedding": "u=e1, v=d*f1+e2, alpha=a in U+U+<-2>",
            "pairwise_nonisometric": "absolute discriminants 2*d^2 are distinct",
        },
        "conditional_global_host_lemma": {
            "hypothesis": "for each rho<=20 choose rho physical (-2)-facets spanning NS_Q with pairings bounded by B(rho)",
            "finite_gram_reason": "diagonal -2 and integral off-diagonal entries in a finite box",
            "finite_saturation_reason": "R is nonsingular and R is contained in NS which is contained in R-dual",
            "finite_overlattice_reason": "R-dual/R is finite",
            "maximality_role": "downstream one-complete-mirror-chamber orbit per host",
        },
        "scope": {
            "exact_Picard_K3_periods_exist_for_each_M_d": True,
            "primitive_nef_isotropic_class_gives_genus_one_fibration": True,
            "cusp_divisibility_is_the_numerical_multisection_index": True,
            "fixed_ADE_or_affine_cusp_packet_bounds_divisibility": False,
            "K3_rank_bound_bounds_divisibility": False,
            "K3_primitive_embedding_bounds_divisibility": False,
            "finite_volume_claimed_for_M_d": False,
            "reflection_maximality_claimed_for_M_d": False,
            "infinite_family_under_all_target_hypotheses_claimed": False,
            "global_bounded_spanning_root_lemma_proved_here": False,
            "finite_host_theorem_proved_without_that_lemma": False,
            "Property_4_proved": False,
        },
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args(argv)
    encoded = encoded_certificate()
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
