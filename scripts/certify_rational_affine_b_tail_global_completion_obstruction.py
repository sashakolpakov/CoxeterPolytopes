#!/usr/bin/env python3
"""Exact certificate for the ruled affine-B-tail completion obstruction.

The computation is deliberately elementary.  It verifies the integral
wall basis, all dual rays of the visible simplex, the m<=4 finite-volume
threshold, the ruling-swap reflection, and the first compatible truncating
root family.  The ample-cone crossing and reflection-overgroup arguments
are symbolic proofs in the accompanying report.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "rational_affine_b_tail_global_completion_obstruction.json"
)


Vector = tuple[int, ...]
Matrix = list[list[Fraction]]


def pair(left: Sequence[int | Fraction], right: Sequence[int | Fraction]) -> Fraction:
    """Intersection form U + <-1>^m in c,d,e_1,... coordinates."""

    return (
        Fraction(left[0]) * Fraction(right[1])
        + Fraction(left[1]) * Fraction(right[0])
        - sum(
            Fraction(a) * Fraction(b)
            for a, b in zip(left[2:], right[2:])
        )
    )


def identity(size: int) -> Matrix:
    return [
        [Fraction(int(row == column)) for column in range(size)]
        for row in range(size)
    ]


def transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix)]


def multiply(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum(a * b for a, b in zip(row, column))
            for column in transpose(right)
        ]
        for row in left
    ]


def inverse(matrix: Matrix) -> Matrix:
    size = len(matrix)
    work = [row[:] + unit for row, unit in zip(matrix, identity(size))]
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        value = work[column][column]
        work[column] = [entry / value for entry in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            coefficient = work[row][column]
            work[row] = [
                entry - coefficient * pivot_entry
                for entry, pivot_entry in zip(work[row], work[column])
            ]
    return [row[size:] for row in work]


def determinant(matrix: Matrix) -> Fraction:
    work = [row[:] for row in matrix]
    value = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            value = -value
        pivot_value = work[column][column]
        value *= pivot_value
        for row in range(column + 1, len(work)):
            if not work[row][column]:
                continue
            coefficient = work[row][column] / pivot_value
            for index in range(column, len(work)):
                work[row][index] -= coefficient * work[column][index]
    return value


def visible_walls(length: int) -> tuple[list[str], list[Vector]]:
    if length < 2:
        raise ValueError("the affine B_m packet starts at m=2")
    size = length + 2
    labels = ["D", "delta_0"]
    walls: list[Vector] = []

    section = [0] * size
    section[1] = 1
    section[2] = -1
    walls.append(tuple(section))

    delta_zero = [0] * size
    delta_zero[0] = 1
    delta_zero[2] = delta_zero[3] = -1
    walls.append(tuple(delta_zero))

    for index in range(1, length):
        delta = [0] * size
        delta[1 + index] = 1
        delta[2 + index] = -1
        walls.append(tuple(delta))
        labels.append(f"delta_{index}")

    endpoint = [0] * size
    endpoint[-1] = 1
    walls.append(tuple(endpoint))
    labels.append(f"delta_{length}")
    return labels, walls


def gram(vectors: Sequence[Vector]) -> Matrix:
    return [
        [pair(left, right) for right in vectors]
        for left in vectors
    ]


def linear_combination(coefficients: Sequence[Fraction], vectors: Sequence[Vector]) -> tuple[Fraction, ...]:
    return tuple(
        sum(coefficient * vector[index] for coefficient, vector in zip(coefficients, vectors))
        for index in range(len(vectors[0]))
    )


def expected_dual(label: str, length: int) -> Vector:
    size = length + 2
    if label == "D":
        return (1, 0, *((0,) * length))
    if label == "delta_0":
        return (0, 1, *((0,) * length))
    index = int(label.split("_")[1])
    if index == 1:
        return (1, 1, -1, *((0,) * (length - 1)))
    return (1, 2, *((-1,) * index), *((0,) * (length - index)))


def visible_certificate(length: int) -> dict[str, object]:
    labels, walls = visible_walls(length)
    matrix = gram(walls)
    matrix_inverse = inverse(matrix)
    duals = [
        linear_combination([matrix_inverse[row][column] for row in range(len(walls))], walls)
        for column in range(len(walls))
    ]

    expected = [expected_dual(label, length) for label in labels]
    if duals != expected:
        raise AssertionError((length, labels, duals, expected))
    if any(
        pair(dual, wall) != int(i == j)
        for i, dual in enumerate(duals)
        for j, wall in enumerate(walls)
    ):
        raise AssertionError("dual-wall Kronecker pairing failed")

    squares = [int(pair(vector, vector)) for vector in duals]
    expected_squares = [0, 0, 1] + [4 - index for index in range(2, length + 1)]
    if squares != expected_squares:
        raise AssertionError((length, squares, expected_squares))
    if abs(determinant(matrix)) != 1:
        raise AssertionError("the visible walls are not a unimodular basis")

    return {
        "m": length,
        "wall_count": len(walls),
        "gram_determinant": int(determinant(matrix)),
        "dual_squares": squares,
        "finite_volume_visible_simplex": all(square >= 0 for square in squares),
        "hyperideal_wall_labels": [
            label for label, square in zip(labels, squares) if square < 0
        ],
    }


def form_matrix(length: int) -> Matrix:
    size = length + 2
    result = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    result[0][1] = result[1][0] = Fraction(1)
    for index in range(2, size):
        result[index][index] = Fraction(-1)
    return result


def swap_reflection_certificate(length: int) -> dict[str, object]:
    size = length + 2
    root = (1, -1, *((0,) * length))
    if pair(root, root) != -2:
        raise AssertionError("the ruling-swap vector is not a root")
    q_root = [sum(form_matrix(length)[row][column] * root[column] for column in range(size)) for row in range(size)]
    reflection = identity(size)
    for row in range(size):
        for column in range(size):
            reflection[row][column] += Fraction(root[row]) * q_root[column]

    if any(entry.denominator != 1 for row in reflection for entry in row):
        raise AssertionError("the ruling-swap reflection is not integral")
    if multiply(reflection, reflection) != identity(size):
        raise AssertionError("the ruling-swap reflection is not an involution")
    if multiply(transpose(reflection), multiply(form_matrix(length), reflection)) != form_matrix(length):
        raise AssertionError("the ruling-swap reflection is not an isometry")

    c = [[Fraction(int(row == 0))] for row in range(size)]
    d = [[Fraction(int(row == 1))] for row in range(size)]
    if multiply(reflection, c) != d or multiply(reflection, d) != c:
        raise AssertionError("the root reflection does not swap the rulings")

    return {
        "root": list(root),
        "surface_square": -2,
        "doubled_host_square": -4,
        "reflection_integral": True,
        "reflection_swaps_c_and_d": True,
        "pullback_c_plus_d_pairing": int(pair((1, 1, *((0,) * length)), root)),
    }


def truncator_certificate(length: int) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    if length >= 5:
        vector = (1, 2, *((-1,) * 5), *((0,) * (length - 5)))
        records.append(
            {
                "label": "T_5",
                "vector": list(vector),
                "support_length": 5,
                "square": int(pair(vector, vector)),
                "canonical_pairing": -1,
                "arithmetic_genus": 0,
                "pairing_with_first_escape_v5": int(pair(vector, vector)),
                "status": (
                    "effective divisor class; irreducibility, smoothness, "
                    "and facet status are not automatic"
                ),
            }
        )
    for a in range(1, (length - 2) // 4 + 1):
        support = 4 * a + 2
        vector = (a, 2, *((-1,) * support), *((0,) * (length - support)))
        square = int(pair(vector, vector))
        canonical_pairing = -2 * a - 4 + support
        genus = 1 + (square + canonical_pairing) // 2
        if (square, canonical_pairing, genus) != (-2, 2 * a - 2, a - 1):
            raise AssertionError("truncator adjunction formula failed")

        labels, walls = visible_walls(length)
        intersections = {
            label: int(pair(vector, wall))
            for label, wall in zip(labels, walls)
            if pair(vector, wall)
        }
        expected_tail_label = f"delta_{support}"
        if intersections.get(expected_tail_label) != 1:
            raise AssertionError((length, a, intersections))
        if intersections.get("D", 0) != a - 1:
            raise AssertionError((length, a, intersections))
        if any(
            label not in {"D", expected_tail_label}
            for label in intersections
        ):
            raise AssertionError((length, a, intersections))

        v5 = (1, 2, *((-1,) * 5), *((0,) * (length - 5)))
        records.append(
            {
                "label": f"R_{a}",
                "vector": list(vector),
                "support_length": support,
                "square": square,
                "canonical_pairing": canonical_pairing,
                "arithmetic_genus": genus,
                "pairing_with_first_escape_v5": int(pair(vector, v5)),
                "nonzero_visible_intersections": intersections,
                "linear_system_h0": 3 * a + 3,
                "jet_condition_count": support,
                "minimum_jet_row_dependencies_for_effectivity": a,
                "irreducible_exact_multiplicity_locus_codimension_lower_bound": a,
            }
        )
    return records


def affine_h3_certificate() -> dict[str, object]:
    """Check the affine-H3 star and the 2-elementary saturation countermodel."""

    physical_star: Matrix = [
        [Fraction(-1), Fraction(1), Fraction(1), Fraction(1), Fraction(1)],
        [Fraction(1), Fraction(-4), Fraction(0), Fraction(0), Fraction(0)],
        [Fraction(1), Fraction(0), Fraction(-4), Fraction(0), Fraction(0)],
        [Fraction(1), Fraction(0), Fraction(0), Fraction(-4), Fraction(0)],
        [Fraction(1), Fraction(0), Fraction(0), Fraction(0), Fraction(-4)],
    ]
    fibre = [[Fraction(value)] for value in (4, 1, 1, 1, 1)]
    if multiply(physical_star, fibre) != [[Fraction(0)] for _ in range(5)]:
        raise AssertionError("the affine-H3 divisor is not null")
    canonical_pairings = (-1, 2, 2, 2, 2)
    canonical_fibre = sum(
        coefficient * pairing
        for coefficient, pairing in zip((4, 1, 1, 1, 1), canonical_pairings)
    )
    genus = 1 + canonical_fibre // 2
    if (canonical_fibre, genus) != (4, 3):
        raise AssertionError("affine-H3 adjunction failed")

    # Negative D4 Cartan matrix: the first finite simple root is central.
    d4: Matrix = [
        [Fraction(-2), Fraction(1), Fraction(1), Fraction(1)],
        [Fraction(1), Fraction(-2), Fraction(0), Fraction(0)],
        [Fraction(1), Fraction(0), Fraction(-2), Fraction(0)],
        [Fraction(1), Fraction(0), Fraction(0), Fraction(-2)],
    ]
    if abs(determinant(d4)) != 4:
        raise AssertionError("D4 determinant is not four")
    if any(
        (2 * entry).denominator != 1
        for row in inverse(d4)
        for entry in row
    ):
        raise AssertionError("the D4 discriminant does not have exponent two")

    samples = []
    for saturation_index in (1, 2, 3, 5, 11, 37):
        # Ambient basis u,v,alpha_central,alpha_leaf1,...,alpha_leaf3.
        ambient: Matrix = [[Fraction(0) for _ in range(6)] for _ in range(6)]
        ambient[0][1] = ambient[1][0] = Fraction(1)
        for row in range(4):
            for column in range(4):
                ambient[2 + row][2 + column] = d4[row][column]

        highest_root = (0, 0, 2, 1, 1, 1)
        affine_root = (
            saturation_index,
            0,
            -2,
            -1,
            -1,
            -1,
        )
        finite_roots = []
        for index in range(4):
            vector = [0] * 6
            vector[2 + index] = 1
            finite_roots.append(tuple(vector))
        roots = (affine_root, *finite_roots)
        affine_gram = [
            [
                sum(
                    Fraction(left[i]) * ambient[i][j] * Fraction(right[j])
                    for i in range(6)
                    for j in range(6)
                )
                for right in roots
            ]
            for left in roots
        ]
        expected_gram: Matrix = [
            [Fraction(-2), Fraction(1), Fraction(0), Fraction(0), Fraction(0)],
            [Fraction(1), Fraction(-2), Fraction(1), Fraction(1), Fraction(1)],
            [Fraction(0), Fraction(1), Fraction(-2), Fraction(0), Fraction(0)],
            [Fraction(0), Fraction(1), Fraction(0), Fraction(-2), Fraction(0)],
            [Fraction(0), Fraction(1), Fraction(0), Fraction(0), Fraction(-2)],
        ]
        if affine_gram != expected_gram:
            raise AssertionError((saturation_index, affine_gram))

        mark_sum = tuple(
            roots[0][coordinate]
            + 2 * roots[1][coordinate]
            + roots[2][coordinate]
            + roots[3][coordinate]
            + roots[4][coordinate]
            for coordinate in range(6)
        )
        if mark_sum != (saturation_index, 0, 0, 0, 0, 0):
            raise AssertionError((saturation_index, mark_sum, highest_root))
        samples.append(
            {
                "N": saturation_index,
                "affine_null_relation": list(mark_sum),
                "ambient_discriminant_exponent": 2,
                "affine_gram_independent_of_N": True,
            }
        )

    return {
        "physical_null_coefficients": [4, 1, 1, 1, 1],
        "canonical_pairing": canonical_fibre,
        "arithmetic_genus": genus,
        "semiample_index_over_C": 1,
        "index_one_mechanism": (
            "Pic0 of the nonreduced rational tree is unipotent and torsion-free; "
            "O_F(F) torsion therefore becomes trivial"
        ),
        "relatively_minimal_core_picard_bound": 26,
        "genus_two_branch": False,
        "lattice_saturation_countermodel": samples,
    }


def build_certificate() -> dict[str, object]:
    visible = [visible_certificate(length) for length in range(2, 41)]
    finite_lengths = [
        record["m"]
        for record in visible
        if record["finite_volume_visible_simplex"]
    ]
    if finite_lengths != [2, 3, 4]:
        raise AssertionError(finite_lengths)

    swap_samples = [swap_reflection_certificate(length) for length in (2, 5, 10, 25, 40)]
    truncators = {
        str(length): truncator_certificate(length)
        for length in range(5, 41)
    }
    return {
        "theorem": "P1xP1 affine-B-tail global completion obstruction",
        "visible_simplex_finite_volume_iff": "2 <= m <= 4",
        "absolute_reflection_maximality_on_full_NS": False,
        "obstructing_root": "c-d",
        "visible_samples": visible,
        "swap_reflection_samples": swap_samples,
        "first_truncating_root_samples": truncators,
        "affine_h3_index_gate": affine_h3_certificate(),
        "claim_boundary": (
            "This excludes absolute full-NS reflection maximality for the explicit "
            "P1xP1 blowup family.  It does not exclude finite-volume nonmaximal nef "
            "cones or maximality relative to a proper host omitting c-d."
        ),
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    payload = build_certificate()
    archive = {
        "payload": payload,
        "payload_sha256": digest(payload),
    }
    rendered = json.dumps(archive, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.exists():
            raise FileNotFoundError(arguments.output)
        if arguments.output.read_text() != rendered:
            raise AssertionError("archive differs from exact recomputation")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print(archive["payload_sha256"])


if __name__ == "__main__":
    main()
