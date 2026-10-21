#!/usr/bin/env python3
"""Exact lattice replay for the fixed-Jacobian K3 cusp-index no-go."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
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
            (row for row in range(column, size) if work[row][column]),
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


def rank(matrix: Matrix) -> int:
    """Exact rational row rank."""

    if not matrix:
        return 0
    work = [[Fraction(value) for value in row] for row in matrix]
    rows = len(work)
    columns = len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][column]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            value = work[row][column]
            work[row] = [
                work[row][target] - value * work[pivot_row][target]
                for target in range(columns)
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def transpose(matrix: Matrix) -> Matrix:
    return [list(row) for row in zip(*matrix)]


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum(left[row][middle] * right[middle][column] for middle in range(len(right)))
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def matvec(matrix: Matrix, vector: Vector) -> Vector:
    return [sum(row[column] * vector[column] for column in range(len(vector))) for row in matrix]


def bilinear(left: Vector, form: Matrix, right: Vector) -> int:
    return sum(
        left[row] * form[row][column] * right[column]
        for row in range(len(left))
        for column in range(len(right))
    )


def block_diagonal(*blocks: Matrix) -> Matrix:
    size = sum(len(block) for block in blocks)
    result = [[0 for _ in range(size)] for _ in range(size)]
    offset = 0
    for block in blocks:
        for row in range(len(block)):
            for column in range(len(block)):
                result[offset + row][offset + column] = block[row][column]
        offset += len(block)
    return result


def e8_cartan() -> Matrix:
    """Bourbaki E8 Cartan matrix with branch at the fourth displayed node."""

    result = [[2 * int(row == column) for column in range(8)] for row in range(8)]
    edges = [(0, 2), (2, 3), (1, 3), (3, 4), (4, 5), (5, 6), (6, 7)]
    for left, right in edges:
        result[left][right] = -1
        result[right][left] = -1
    return result


def scaled_u(scale: int) -> Matrix:
    return [[0, scale], [scale, 0]]


def unit_vector(size: int, position: int) -> Vector:
    return [int(index == position) for index in range(size)]


def family_record(index: int) -> dict[str, object]:
    if index < 2:
        raise ValueError("the no-go family is recorded for d>=2")
    d = index

    e8_positive = e8_cartan()
    if determinant(e8_positive) != 1:
        raise AssertionError("the E8 Cartan matrix stopped being unimodular")
    e8_negative = [[-value for value in row] for row in e8_positive]
    host = block_diagonal(scaled_u(d), e8_negative, e8_negative)
    if determinant(host) != -(d * d):
        raise AssertionError("the host discriminant formula changed")

    # Lambda_K3 = U^3 + E8(-1)^2.  Embed u=e1 and v=d*f1+e2.
    ambient = block_diagonal(scaled_u(1), scaled_u(1), scaled_u(1), e8_negative, e8_negative)
    embedding = [[0 for _ in range(18)] for _ in range(22)]
    embedding[0][0] = 1
    embedding[1][1] = d
    embedding[2][1] = 1
    for coordinate in range(16):
        embedding[6 + coordinate][2 + coordinate] = 1
    pullback = matmul(transpose(embedding), matmul(ambient, embedding))
    if pullback != host:
        raise AssertionError("the explicit primitive K3 embedding changed")
    # Rows e1,e2 and the sixteen E8 rows form the identity maximal minor.
    selected_rows = [0, 2, *range(6, 22)]
    unit_minor = [[embedding[row][column] for column in range(18)] for row in selected_rows]
    if determinant(unit_minor) != 1:
        raise AssertionError("the primitive unit minor changed")

    # Its exact orthogonal complement is U(d)+U.
    complement = [[0 for _ in range(4)] for _ in range(22)]
    complement[0][0] = -1       # -e1+d*f2
    complement[3][0] = d
    complement[2][1] = 1        # e2
    complement[4][2] = 1        # e3
    complement[5][3] = 1        # f3
    complement_gram = matmul(transpose(complement), matmul(ambient, complement))
    expected_complement = block_diagonal(scaled_u(d), scaled_u(1))
    if complement_gram != expected_complement:
        raise AssertionError("the orthogonal complement formula changed")
    # The displayed T=U(d)+U embedding is independently primitive: the
    # e1,e2,e3,f3 rows give a determinant-minus-one maximal minor.
    complement_rows = [0, 2, 4, 5]
    complement_minor = [
        [complement[row][column] for column in range(4)]
        for row in complement_rows
    ]
    if determinant(complement_minor) != -1:
        raise AssertionError("the transcendental embedding stopped being primitive")
    if matmul(transpose(embedding), matmul(ambient, complement)) != [[0] * 4 for _ in range(18)]:
        raise AssertionError("the displayed host and complement are not orthogonal")

    fibre = unit_vector(18, 0)
    if bilinear(fibre, host, fibre):
        raise AssertionError("the fibre is not isotropic")
    fibre_divisibility = gcd(*(abs(value) for value in matvec(host, fibre)))
    if fibre_divisibility != d:
        raise AssertionError("the fibre divisibility formula changed")

    # Highest-root coefficients for the chosen E8 numbering.
    marks = [2, 3, 4, 6, 5, 4, 3, 2]
    if matvec(e8_positive, marks) != [0, 0, 0, 0, 0, 0, 0, 1]:
        raise AssertionError("the E8 highest-root marks changed")

    components: list[Vector] = []
    for block in (0, 1):
        offset = 2 + 8 * block
        simples = [unit_vector(18, offset + position) for position in range(8)]
        theta = [sum(marks[position] * simples[position][coordinate] for position in range(8)) for coordinate in range(18)]
        affine = [fibre[coordinate] - theta[coordinate] for coordinate in range(18)]
        if bilinear(theta, host, theta) != -2 or bilinear(affine, host, affine) != -2:
            raise AssertionError("an affine E8 component stopped being a root")
        relation = [affine[coordinate] + theta[coordinate] for coordinate in range(18)]
        if relation != fibre:
            raise AssertionError("the affine fibre relation changed")
        components.extend([affine, *simples])

    component_columns = transpose(components)
    component_gram = matmul(transpose(component_columns), matmul(host, component_columns))
    if rank(component_columns) != 17 or rank(component_gram) != 16:
        raise AssertionError("the full local parabolic ranks changed")
    common_fixed_equations = matmul(transpose(component_columns), host)
    if rank(common_fixed_equations) != 17:
        raise AssertionError("the incident components no longer have a one-dimensional common orthogonal")
    if any(bilinear(fibre, host, component) for component in components):
        raise AssertionError("an incident fibre component misses the cusp")

    theta_first = [sum(marks[position] * components[1 + position][coordinate] for position in range(8)) for coordinate in range(18)]
    transverse_root = [unit_vector(18, 1)[coordinate] - theta_first[coordinate] for coordinate in range(18)]
    if bilinear(transverse_root, host, transverse_root) != -2:
        raise AssertionError("the displayed transverse lattice root changed norm")
    if bilinear(fibre, host, transverse_root) != d:
        raise AssertionError("the transverse root degree formula changed")
    incident_pairings = [bilinear(transverse_root, host, component) for component in components]
    if min(incident_pairings) < 0:
        raise AssertionError("the transverse root stopped being acute-compatible with the cusp packet")

    return {
        "d": d,
        "picard_rank": 18,
        "host": "U(d)+E8(-1)^2",
        "signature": [1, 17],
        "absolute_discriminant": d * d,
        "discriminant_group_invariant_factors": [d, d],
        "discriminant_length": 2,
        "primitive_embedding_unit_minor": 1,
        "transcendental_embedding_unit_minor": -1,
        "orthogonal_complement": "U(d)+U",
        "orthogonal_complement_gram": complement_gram,
        "actual_fibre_packet_marking": "F=u and vertical roots=E8(-1)^2",
        "primitive_fibre_divisibility": fibre_divisibility,
        "multisection_index": d,
        "fixed_fibres": ["II*", "II*", "I1", "I1", "I1", "I1"],
        "euler_sum": 10 + 10 + 4,
        "jacobian_mordell_weil_rank": 0,
        "affine_packet": "affine_E8+affine_E8",
        "incident_curve_count": len(components),
        "incident_vector_span_rank": rank(component_columns),
        "incident_gram_rank": rank(component_gram),
        "finite_root_rank": 16,
        "expected_finite_root_rank": 18 - 2,
        "common_orthogonal_dimension": 18 - rank(common_fixed_equations),
        "transverse_lattice_root_square": -2,
        "transverse_lattice_root_fibre_degree": d,
        "transverse_root_incident_pairings": sorted(incident_pairings),
        "all_positive_divisor_fibre_degrees_divisible_by_d": True,
    }


def build_certificate() -> dict[str, object]:
    sample_indices = [2, 3, 4, 5, 7, 11, 19, 31]
    return {
        "verdict": "fixed-Jacobian elliptic and K3-embedding data do not bound cusp divisibility",
        "family": [family_record(index) for index in sample_indices],
        "symbolic_family": {
            "fixed_jacobian": "very general S with NS(S)=U+E8(-1)^2 and fibres 2*II*+4*I1",
            "brauer_class": "alpha_d:T(S)=U^2 -> Q/Z, alpha_d(a*e1+...)=a/d mod Z",
            "brauer_torsion": "Br(S)[d]=(Z/d)^4 and alpha_d has exact order d",
            "algebraic_torsor": "X_d is the projective Ogg-Shafarevich torsor associated with alpha_d",
            "period_index": "div(F_d)=multisection_index(X_d)=ord(alpha_d)=d",
            "transcendental_kernel": "T(X_d)=ker(alpha_d)=U(d)+U",
            "neron_severi_complement": "NS(X_d)=U(d)+E8(-1)^2",
            "local_cusp_packet": "2*affine_E8, finite root rank 16=rho-2, common orthogonal Q*F_d",
            "escaping_height": "every transverse divisor degree is a positive multiple of d",
        },
        "source_dependencies": {
            "Ogg_Shafarevich": "Br(S)=algebraic Sha(S) parametrizes projective K3 torsors with relative Jacobian S",
            "period_equals_index": "Meinsma-Shinder 2024, Lemma 4.4(iii)",
            "transcendental_kernel": "Meinsma-Shinder 2024, Lemma 4.4(iv)",
            "K3_complement": "Huybrechts Theorem 14.1.12/Nikulin Theorem 1.14.4 plus explicit primitive T embedding and exact orthogonal complement",
            "actual_packet_alignment": "unimodular vertical E8(-1)^2 splitting plus the discriminant-form parity exclusion of the even d/2 binary lattice",
            "physical_affine_packet": "Huybrechts Chapter 11 Section 4.1 equation (4.1): every closed fibre agrees with every relative Jacobian fibre",
        },
        "scope": {
            "same_fixed_relative_jacobian_for_all_d": True,
            "projective_K3_torsor_for_all_d": True,
            "same_Kodaira_and_Euler_packet_for_all_d": True,
            "full_local_parabolic_rank_for_selected_cusp": True,
            "selected_jacobian_MW_rank_zero": True,
            "selected_fibration_stabilizer_finite": True,
            "K3_primitive_embedding_and_discriminant_length_bound_divisibility": False,
            "Euler_24_and_Shioda_Tate_bound_divisibility": False,
            "full_finite_volume_chamber_claimed_for_X_d": False,
            "full_automorphism_group_finite_claimed_for_X_d": False,
            "displayed_transverse_lattice_root_claimed_irreducible_physical": False,
            "finite_volume_counterfamily_claimed": False,
            "global_bounded_spanning_root_theorem_proved": False,
            "published_118_host_classification_reproved": False,
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
