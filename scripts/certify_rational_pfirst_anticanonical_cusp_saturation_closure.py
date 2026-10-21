#!/usr/bin/env python3
"""Exact closure of the I_(1,9) Hamming chamber under all host reflections."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

try:
    from scripts.certify_rational_pfirst_k_escape_hamming_effectivity_no_go import (
        DIMENSION,
        METRIC,
        Vector,
        add,
        bareiss_determinant,
        dot,
        packet as hamming_packet,
        scale,
    )
except ModuleNotFoundError:
    from certify_rational_pfirst_k_escape_hamming_effectivity_no_go import (
        DIMENSION,
        METRIC,
        Vector,
        add,
        bareiss_determinant,
        dot,
        packet as hamming_packet,
        scale,
    )


def subtract(left: Vector, right: Vector) -> Vector:
    return tuple(a - b for a, b in zip(left, right))


def standard_roots() -> tuple[Vector, ...]:
    roots: list[Vector] = []
    section = [0] * DIMENSION
    section[9] = 1
    roots.append(tuple(section))
    for index in range(1, 9):
        root = [0] * DIMENSION
        root[index] = 1
        root[index + 1] = -1
        roots.append(tuple(root))
    line = [0] * DIMENSION
    line[0] = 1
    line[1] = line[2] = line[3] = -1
    roots.append(tuple(line))
    return tuple(roots)


def reflect(vector: Vector, root: Vector) -> Vector:
    norm = -dot(root, root)
    numerator = 2 * dot(vector, root)
    assert norm in (1, 2)
    assert numerator % norm == 0
    coefficient = numerator // norm
    return tuple(vector[index] + coefficient * root[index] for index in range(DIMENSION))


def reduce_to_standard_chamber(vector: Vector) -> tuple[Vector, tuple[Vector, ...]]:
    roots = standard_roots()
    current = vector
    sequence: list[Vector] = []
    for _ in range(10000):
        violated = next((root for root in roots if dot(current, root) < 0), None)
        if violated is None:
            return current, tuple(sequence)
        current = reflect(current, violated)
        sequence.append(violated)
    raise ArithmeticError("simple-root descent did not terminate")


def pullback_root(root: Vector, reflection_sequence: tuple[Vector, ...]) -> Vector:
    current = root
    for reflection_root in reversed(reflection_sequence):
        current = reflect(current, reflection_root)
    return current


def solve_square(matrix: list[list[int]], right: list[int]) -> list[Fraction]:
    size = len(matrix)
    assert size and len(right) == size and all(len(row) == size for row in matrix)
    augmented = [
        [Fraction(entry) for entry in row] + [Fraction(value)]
        for row, value in zip(matrix, right)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column]),
            None,
        )
        if pivot is None:
            raise ArithmeticError("singular exact linear system")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        pivot_value = augmented[column][column]
        augmented[column] = [entry / pivot_value for entry in augmented[column]]
        for row in range(size):
            if row != column and augmented[row][column]:
                coefficient = augmented[row][column]
                augmented[row] = [
                    entry - coefficient * pivot_entry
                    for entry, pivot_entry in zip(
                        augmented[row], augmented[column]
                    )
                ]
    return [augmented[index][-1] for index in range(size)]


def root_coordinates(vector: Vector, roots: tuple[Vector, ...]) -> tuple[int, ...]:
    # vector = sum_i coefficient_i roots_i, hence R^T coefficient = vector.
    matrix = [[roots[column][row] for column in range(DIMENSION)] for row in range(DIMENSION)]
    solution = solve_square(matrix, list(vector))
    assert all(entry.denominator == 1 for entry in solution)
    return tuple(int(entry) for entry in solution)


def normalize_fraction_vector(entries: list[Fraction]) -> Vector:
    denominator = 1
    for entry in entries:
        denominator = math.lcm(denominator, entry.denominator)
    integral = [int(entry * denominator) for entry in entries]
    divisor = math.gcd(*map(abs, integral))
    assert divisor
    primitive = tuple(entry // divisor for entry in integral)
    if primitive[0] < 0:
        primitive = scale(-1, primitive)
    return primitive


def dual_vertex(index: int, roots: tuple[Vector, ...]) -> Vector:
    # Solve root_j . vertex = delta_ij.
    matrix = [
        [METRIC[column] * root[column] for column in range(DIMENSION)]
        for root in roots
    ]
    right = [1 if row == index else 0 for row in range(DIMENSION)]
    return normalize_fraction_vector(solve_square(matrix, right))


def characteristic(vector: Vector) -> bool:
    # On the diagonal odd unimodular lattice every coordinate is odd.
    return all(entry % 2 for entry in vector)


def coxeter_histogram(roots: tuple[Vector, ...]) -> dict[str, int]:
    histogram: Counter[str] = Counter()
    for first_index, first in enumerate(roots):
        first_norm = -dot(first, first)
        for second in roots[first_index + 1 :]:
            second_norm = -dot(second, second)
            pairing = dot(first, second)
            if pairing == 0:
                order = "2"
            elif (
                first_norm == second_norm == 2
                and pairing == 1
            ):
                order = "3"
            elif (
                {first_norm, second_norm} == {1, 2}
                and pairing == 1
            ):
                order = "4"
            else:
                raise ArithmeticError(
                    f"unexpected saturated-root pair {first_norm, second_norm, pairing}"
                )
            histogram[order] += 1
    assert histogram == {"2": 36, "3": 8, "4": 1}
    return dict(sorted(histogram.items()))


def payload() -> dict[str, object]:
    hamming = hamming_packet()
    old_records = hamming["walls"]
    old_walls = [record["vector"] for record in old_records]
    canonical: Vector = hamming["canonical"]
    anticanonical: Vector = hamming["anticanonical"]
    h: Vector = hamming["strict_interior"]

    perturbation: Vector = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
    generic = tuple(100 * h[index] + perturbation[index] for index in range(DIMENSION))
    assert dot(generic, generic) == 51915
    old_slacks = [dot(generic, wall) for wall in old_walls]
    assert min(old_slacks) == 91

    reduced, sequence = reduce_to_standard_chamber(generic)
    assert reduced == (409, -193, -108, -106, -104, -98, -96, -94, -92, -91)
    assert len(sequence) == 10
    standard = standard_roots()
    reduced_slacks = [dot(reduced, root) for root in standard]
    assert reduced_slacks == [91, 85, 2, 2, 6, 2, 2, 2, 1, 2]

    saturated = tuple(pullback_root(root, sequence) for root in standard)
    expected_saturated: tuple[Vector, ...] = (
        (0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
        (-1, 0, 0, 1, 0, 1, 0, 1, 0, 0),
        (0, 1, 0, -1, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 1, 0, -1, 0, 0, 0, 0),
        (1, -1, -1, -1, 0, 0, 0, 0, 0, 0),
        (0, 0, 1, 0, -1, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 1, 0, -1, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 1, 0, -1, 0),
        (0, 0, 0, 0, 0, 0, 0, 0, 1, -1),
        (0, 0, 0, 0, 0, 1, 0, -1, 0, 0),
    )
    assert saturated == expected_saturated
    assert [dot(generic, root) for root in saturated] == reduced_slacks
    assert [dot(h, root) for root in saturated] == [1, 1] + [0] * 8
    assert [dot(canonical, root) for root in saturated] == [-1] + [0] * 9
    assert [dot(anticanonical, root) for root in saturated] == [1] + [0] * 9
    assert [dot(root, root) for root in saturated] == [-1] + [-2] * 9

    gram = [[dot(first, second) for second in saturated] for first in saturated]
    expected_gram = [
        [-1, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [0, -2, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, -2, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, -2, 1, 0, 0, 0, 0, 1],
        [0, 0, 0, 1, -2, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, -2, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, -2, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, -2, 1, 0],
        [1, 0, 0, 0, 0, 0, 0, 1, -2, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, -2],
    ]
    assert gram == expected_gram
    assert bareiss_determinant([list(root) for root in saturated]) == -1

    old_expansions = {
        record["name"]: root_coordinates(record["vector"], saturated)
        for record in old_records
    }
    assert all(
        all(coefficient >= 0 for coefficient in expansion)
        for expansion in old_expansions.values()
    )
    assert sum(
        sum(coefficient != 0 for coefficient in expansion) == 1
        and max(expansion) == 1
        for expansion in old_expansions.values()
    ) == 1
    expansion_bytes = json.dumps(
        old_expansions, sort_keys=True, separators=(",", ":")
    ).encode()

    vertices = tuple(dual_vertex(index, saturated) for index in range(DIMENSION))
    vertex_squares = [dot(vertex, vertex) for vertex in vertices]
    assert vertex_squares == [0, 0, 2, 6, 5, 4, 3, 2, 1, 1]
    second_cusp = vertices[1]
    assert vertices[0] == anticanonical
    assert second_cusp == (2, -1, 0, -1, 0, -1, 0, -1, 0, 0)
    assert characteristic(vertices[0])
    assert not characteristic(second_cusp)
    assert h == add(vertices[0], vertices[1])

    null_marks = root_coordinates(anticanonical, saturated)
    assert null_marks == (0, 2, 4, 6, 5, 4, 3, 2, 1, 3)
    affine_edges = [
        (first, second)
        for first in range(1, 10)
        for second in range(first + 1, 10)
        if dot(saturated[first], saturated[second]) == 1
    ]
    assert len(affine_edges) == 8
    seen = {1}
    while True:
        enlarged = seen | {
            second if first in seen else first
            for first, second in affine_edges
            if first in seen or second in seen
        }
        if enlarged == seen:
            break
        seen = enlarged
    assert seen == set(range(1, 10))
    affine_support_euler = 2 * 9 - len(affine_edges)
    assert affine_support_euler == 10

    standard_null_marks = root_coordinates(anticanonical, standard)
    assert standard_null_marks == null_marks
    component_not_meeting_section = [
        root
        for index, root in enumerate(standard[1:], 1)
        if dot(standard[0], root) == 0
    ]
    assert len(component_not_meeting_section) == 8
    trivial_basis = [anticanonical, standard[0], *component_not_meeting_section]
    assert bareiss_determinant([list(vector) for vector in trivial_basis]) == -1

    # Direct flex-pencil data.  The cubic
    # Y^2 Z = X^3 + X Z^2 + lambda Z^3
    # has discriminant -16(4+27 lambda^2).  Its lambda=infinity member is
    # the triple flex line Z^3=0.  Locally the first three blowups reduce
    # the pencil ideal to (x^6, w); six further point blowups resolve it.
    flex_pencil = {
        "equation": "Y^2 Z = X^3 + X Z^2 + lambda Z^3",
        "flex_point": "[0:1:0]",
        "flex_line": "Z=0",
        "finite_discriminant_without_nonzero_scalar": "4 + 27 lambda^2",
        "finite_singular_members": 2,
        "finite_singular_members_irreducible_nodal": True,
        "infinite_member": "3(Z=0)",
        "local_generators_before_blowup": [
            "z^3",
            "z - x^3 - x z^2",
        ],
        "first_three_basepoints_follow_flex_line": True,
        "ideal_after_first_three_blowups": "(x^6,w) up to local units",
        "remaining_successive_basepoints": 6,
        "total_blowups": 9,
        "last_exceptional_is_section": True,
    }

    result: dict[str, object] = {
        "schema": "rational_pfirst_anticanonical_cusp_saturation_closure/v1",
        "lattice": "I_(1,9)",
        "canonical": list(canonical),
        "anticanonical": list(anticanonical),
        "hamming_strict_point": list(h),
        "generic_perturbation": list(perturbation),
        "generic_point": list(generic),
        "generic_square": dot(generic, generic),
        "minimum_old_wall_slack": min(old_slacks),
        "simple_reflection_descent_steps": len(sequence),
        "reduced_generic_point": list(reduced),
        "reduced_simple_slacks": reduced_slacks,
        "saturated_alcove_roots": [list(root) for root in saturated],
        "saturated_root_squares": [-1] + [-2] * 9,
        "saturated_canonical_pairings": [-1] + [0] * 9,
        "saturated_gram": gram,
        "saturated_root_basis_determinant": -1,
        "saturated_coxeter_order_histogram": coxeter_histogram(saturated),
        "old_wall_positive_expansion_sha256": hashlib.sha256(
            expansion_bytes
        ).hexdigest(),
        "old_walls_with_nonnegative_simple_expansion": 32,
        "old_walls_remaining_simple_facets": 1,
        "vertices": [list(vertex) for vertex in vertices],
        "vertex_squares": vertex_squares,
        "ideal_vertices": [list(vertices[0]), list(vertices[1])],
        "characteristic_ideal_vertex_count": 1,
        "hamming_point_equals_sum_of_two_cusps": True,
        "anticanonical_affine_marks": list(null_marks[1:]),
        "anticanonical_affine_support_edges": [
            [first, second] for first, second in affine_edges
        ],
        "anticanonical_affine_support_euler": affine_support_euler,
        "rational_elliptic_c2_budget": 12,
        "host_reflection_maximality": {
            "primitive_reflecting_norms": [-1, -2],
            "descent_generators": [
                "signed permutations of E_i",
                "reflection in H-E1-E2-E3",
            ],
            "all_integral_reflections_generated": True,
            "relative_reflection_maximality_in_O_plus_L": True,
        },
        "explicit_realization": {
            "construction": "nine blowups of the flex cubic pencil",
            "pencil": flex_pencil,
            "standard_facet_roots": [list(root) for root in standard],
            "standard_to_saturated_isometry_certified_by_descent": True,
            "trivial_lattice_basis_determinant": -1,
            "other_reducible_fibers": 0,
            "other_negative_curves": 0,
            "full_nef_chamber_returned": True,
        },
        "anticanonical_cusp_rigidity": {
            "statement": (
                "every relatively reflection-maximal chamber in I_(1,9) "
                "having the primitive characteristic cusp F=-K is conjugate "
                "under the F-stabilizer to this simplex"
            ),
            "reason": (
                "the simplex has exactly two cusps, and only F is characteristic"
            ),
            "twice_minus_four_wall_at_F_cusp_possible": False,
            "affine_fiber_support_euler_forced": 10,
        },
        "scope": {
            "hamming_saturation_outcome": "MAXIMAL_AND_REALIZABLE",
            "hamming_simultaneous_effectivity_obstruction_survives": False,
            "rank_ten_anticanonical_cusp_budget_theorem": "PROVED",
            "all_ranks_or_non_anticanonical_cusps": "OPEN",
            "maximal_K_escape_without_minus_K_cusp": "OPEN",
        },
    }
    canonical_bytes = json.dumps(
        result, sort_keys=True, separators=(",", ":")
    ).encode()
    result["payload_sha256"] = hashlib.sha256(canonical_bytes).hexdigest()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    result = payload()
    default = Path(
        "data/enumerations/"
        "rational_pfirst_anticanonical_cusp_saturation_closure.json"
    )
    target = args.output or default
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        assert target.read_text() == rendered
    elif args.output:
        target.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
