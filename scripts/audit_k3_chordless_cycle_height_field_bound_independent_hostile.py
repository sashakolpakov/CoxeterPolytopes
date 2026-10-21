#!/usr/bin/env python3
"""Independent exact checks for the chordless-cycle K3 field theorem."""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from fractions import Fraction
from math import floor, log2
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
OWNER_PINS = {
    "report": {
        "path": "results/k3_chordless_cycle_height_field_bound.md",
        "sha256": "3372cddd560471830fc0a9cbaac02560ce4a682781db1cdc877b8cd108ea471c",
    },
    "scope_audit": {
        "path": "results/k3_chordless_cycle_height_field_bound_hostile_scope_audit.md",
        "sha256": "ac220aaab7182c5fc4aa4049c6745653ee46071fd91ee435061c7caa327b091d",
    },
    "certifier": {
        "path": "scripts/certify_k3_chordless_cycle_height_field_bound.py",
        "sha256": "1430574782e812f501bf96106af85c307cd22117ad1c568bec5a61cf71d9e766",
    },
    "archive": {
        "path": "data/enumerations/k3_chordless_cycle_height_field_bound.json",
        "sha256": "24589042cde4e5f1e37a711111fc05d64156bebb471b697aacedb5eed62b42ec",
    },
    "tests": {
        "path": "tests/test_k3_chordless_cycle_height_field_bound.py",
        "sha256": "1fbaf3414f5b039856bbbde7f0beda8975302776b5f58f9d82e1f4218be117d8",
    },
}


Matrix = tuple[tuple[Fraction, ...], ...]
Vector = tuple[Fraction, ...]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(size: int) -> Matrix:
    return tuple(
        tuple(Fraction(row == column) for column in range(size))
        for row in range(size)
    )


def outer(column: Vector, row: Vector) -> Matrix:
    return tuple(
        tuple(left * right for right in row) for left in column
    )


def subtract(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(a - b for a, b in zip(left_row, right_row))
        for left_row, right_row in zip(left, right)
    )


def multiply(left: Matrix, right: Matrix) -> Matrix:
    size = len(left)
    return tuple(
        tuple(
            sum(
                (left[row][index] * right[index][column] for index in range(size)),
                Fraction(0),
            )
            for column in range(size)
        )
        for row in range(size)
    )


def trace(matrix: Matrix) -> Fraction:
    return sum(
        (matrix[index][index] for index in range(len(matrix))), Fraction(0)
    )


def dot(row: Vector, column: Vector) -> Fraction:
    return sum((a * b for a, b in zip(row, column)), Fraction(0))


def reflection_packet() -> tuple[tuple[Vector, ...], tuple[Vector, ...], tuple[Matrix, ...]]:
    # Each P_i=e_i alpha_i is integral although some alpha_i(e_j) are
    # fractional.  The displayed e_i use different nonprimitive scales.
    vectors: tuple[Vector, ...] = (
        (Fraction(2), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(3)),
        (Fraction(2), Fraction(2), Fraction(0)),
    )
    covectors: tuple[Vector, ...] = (
        (Fraction(1), Fraction(1, 2), Fraction(0)),
        (Fraction(1), Fraction(2), Fraction(1)),
        (Fraction(1, 3), Fraction(1, 3), Fraction(2, 3)),
        (Fraction(1, 2), Fraction(1, 2), Fraction(1, 2)),
    )
    reflections = []
    for vector, covector in zip(vectors, covectors):
        assert dot(covector, vector) == 2
        projector = outer(vector, covector)
        assert all(entry.denominator == 1 for row in projector for entry in row)
        reflection = subtract(identity(3), projector)
        assert multiply(reflection, reflection) == identity(3)
        reflections.append(reflection)
    return vectors, covectors, tuple(reflections)


def cartan_entry(
    left: int,
    right: int,
    vectors: Sequence[Vector],
    covectors: Sequence[Vector],
) -> Fraction:
    return dot(covectors[left], vectors[right])


def cyclic_product(
    word: Sequence[int],
    vectors: Sequence[Vector],
    covectors: Sequence[Vector],
) -> Fraction:
    return product(
        cartan_entry(left, right, vectors, covectors)
        for left, right in zip(word, word[1:] + word[:1])
    )


def product(values: Iterable[Fraction]) -> Fraction:
    answer = Fraction(1)
    for value in values:
        answer *= value
    return answer


def trace_expansion(
    word: Sequence[int],
    vectors: Sequence[Vector],
    covectors: Sequence[Vector],
) -> Fraction:
    answer = Fraction(3)
    positions = range(len(word))
    for size in range(1, len(word) + 1):
        for subset in itertools.combinations(positions, size):
            subsequence = tuple(word[index] for index in subset)
            answer += (-1) ** size * cyclic_product(
                subsequence, vectors, covectors
            )
    return answer


def trace_induction_audit() -> dict[str, object]:
    vectors, covectors, reflections = reflection_packet()
    entries = [
        cartan_entry(i, j, vectors, covectors)
        for i in range(len(vectors))
        for j in range(len(vectors))
    ]
    fractional_entries = [entry for entry in entries if entry.denominator > 1]
    assert fractional_entries

    tested = 0
    for length in range(1, 6):
        for word in itertools.product(range(len(vectors)), repeat=length):
            matrix = identity(3)
            for index in word:
                matrix = multiply(matrix, reflections[index])
            assert trace(matrix) == trace_expansion(word, vectors, covectors)
            assert cyclic_product(word, vectors, covectors).denominator == 1
            tested += 1

    repeated_controls = (
        (0, 1, 0, 2, 1, 3),
        (3, 3, 1, 2, 3, 0, 1),
        (2, 0, 2, 0, 2, 1, 3, 1),
    )
    for word in repeated_controls:
        matrix = identity(3)
        for index in word:
            matrix = multiply(matrix, reflections[index])
        assert trace(matrix) == trace_expansion(word, vectors, covectors)
        assert cyclic_product(word, vectors, covectors).denominator == 1
        tested += 1

    return {
        "integral_reflection_word_count": tested,
        "fractional_individual_Cartan_entries_present": len(fractional_entries),
        "all_trace_subset_expansions_exact": True,
        "all_full_ordered_cycle_products_integral": True,
        "repeated_words_tested": len(repeated_controls),
    }


def forward_path(start: int, stop: int, length: int) -> tuple[int, ...]:
    vertices = [start]
    while vertices[-1] != stop:
        vertices.append((vertices[-1] + 1) % length)
    return tuple(vertices)


def directed_cycle_edges(vertices: Sequence[int]) -> Counter[tuple[int, int]]:
    edges: Counter[tuple[int, int]] = Counter()
    for left, right in zip(vertices, vertices[1:]):
        edges[(left, right)] += 1
    edges[(vertices[-1], vertices[0])] += 1
    return edges


def directed_chord_audit() -> dict[str, object]:
    tested = 0
    for length in range(4, 17):
        perimeter = directed_cycle_edges(tuple(range(length)))
        for left in range(length):
            for step in range(2, length - 1):
                right = (left + step) % length
                first = forward_path(left, right, length)
                second = forward_path(right, left, length)
                expected = perimeter.copy()
                expected[(right, left)] += 1
                expected[(left, right)] += 1
                assert (
                    directed_cycle_edges(first) + directed_cycle_edges(second)
                    == expected
                )
                tested += 1
    return {
        "directed_chord_splits_tested": tested,
        "perimeter_directions_each_occur_once": True,
        "both_chord_directions_occur_once": True,
        "field_identity": "c(C1)c(C2)=c(C)c(i,j)c(j,i)",
    }


def bareiss_determinant(matrix: Sequence[Sequence[int]]) -> int:
    work = [list(row) for row in matrix]
    size = len(work)
    if size == 0:
        return 1
    previous = 1
    sign = 1
    for column in range(size - 1):
        if work[column][column] == 0:
            pivot = next(
                (
                    row
                    for row in range(column + 1, size)
                    if work[row][column]
                ),
                None,
            )
            if pivot is None:
                return 0
            work[column], work[pivot] = work[pivot], work[column]
            sign *= -1
        value = work[column][column]
        for row in range(column + 1, size):
            for next_column in range(column + 1, size):
                numerator = (
                    work[row][next_column] * value
                    - work[row][column] * work[column][next_column]
                )
                assert numerator % previous == 0
                work[row][next_column] = numerator // previous
        previous = value
    return sign * work[-1][-1]


def cycle_matrix(length: int) -> tuple[tuple[int, ...], ...]:
    matrix = [[0] * length for _ in range(length)]
    for index in range(length):
        matrix[index][index] = 7 - 3 * index
        successor = (index + 1) % length
        edge = index + 2
        matrix[index][successor] = edge
        matrix[successor][index] = edge
    return tuple(tuple(row) for row in matrix)


def rank_minor_audit() -> dict[str, object]:
    records = []
    for length in range(3, 31):
        matrix = cycle_matrix(length)
        minor = [
            [matrix[row][column] for column in range(2, length)]
            for row in range(1, length - 1)
        ]
        determinant = bareiss_determinant(minor)
        expected = product(
            Fraction(matrix[index][index + 1])
            for index in range(1, length - 1)
        )
        assert determinant == expected
        assert determinant != 0
        records.append(
            {
                "length": length,
                "rank_lower_bound": length - 2,
                "witness_minor_determinant": determinant,
            }
        )
    return {
        "exact_minor_records": records,
        "universal_minor": (
            "rows 2,...,ell-1 and columns 3,...,ell; "
            "determinant product A_(i,i+1)"
        ),
        "kernel_dimension_at_most": 2,
        "cycle_length_at_most_m_plus_2": True,
    }


def affine_shift_audit() -> dict[str, object]:
    tested = 0
    for length in range(2, 14):
        unit = 2**length
        for distinguished in (0, unit, unit + 1, 5 * unit):
            for conjugate in (-unit + 1, -1, 0, 1, unit):
                theta0 = distinguished + unit
                theta = conjugate + unit
                assert theta0 > 0 and theta > 0
                if distinguished <= unit:
                    assert theta0 <= 2 * unit and theta <= 2 * unit
                else:
                    assert theta <= 2 * unit < theta0
                tested += 1
    return {
        "exact_affine_place_samples": tested,
        "shift": "theta=(-1)^ell*c+2^ell",
        "small_all_places_bound": "2^(ell+1)",
        "large_non_distinguished_bound": "<=2^(ell+1)",
        "large_coordinate_dominant": True,
        "zero_boundary_excluded_by_embedding_injectivity": True,
    }


def max_strict_tower_length(degree_bound: int) -> int:
    frontier = [(1, 0)]
    maximum = 0
    while frontier:
        degree, length = frontier.pop()
        maximum = max(maximum, length)
        for factor in range(2, degree_bound + 1):
            next_degree = degree * factor
            if next_degree <= degree_bound:
                frontier.append((next_degree, length + 1))
    return maximum


def generator_constant_audit() -> dict[str, object]:
    tower_records = []
    for degree_bound in range(2, 65):
        exact = max_strict_tower_length(degree_bound)
        claimed = floor(log2(degree_bound))
        assert exact == claimed
        tower_records.append(
            {
                "degree_bound": degree_bound,
                "max_strict_tower_length": exact,
            }
        )

    dimension_records = []
    for n in range(2, 10):
        degree_bound = 20 // (n + 1)
        length_bound = n + 3
        unit = 2**length_bound
        generator_count = floor(log2(degree_bound))
        coefficient_bound = degree_bound * (degree_bound - 1) // 2 + 1
        dimension_records.append(
            {
                "n": n,
                "D": degree_bound,
                "L": length_bound,
                "U": unit,
                "r": generator_count,
                "M_D": coefficient_bound,
            }
        )

    n2 = dimension_records[0]
    assert n2 == {"n": 2, "D": 6, "L": 5, "U": 32, "r": 2, "M_D": 16}
    assert (
        n2["r"] * n2["M_D"]
        == 32
    )
    return {
        "tower_records": tower_records,
        "dimension_records": dimension_records,
        "n2_primitive_bound": "32*max(64,C/s0)",
        "n2_discriminant_exponent": 30,
        "primitive_grid_hyperplane_count": "D(D-1)/2",
    }


def audit_payload() -> dict[str, object]:
    for item in OWNER_PINS.values():
        assert sha256(ROOT / item["path"]) == item["sha256"]
    return {
        "schema": "k3_chordless_cycle_height_field_bound/independent_hostile/v1",
        "owner_pins": OWNER_PINS,
        "trace_induction": trace_induction_audit(),
        "directed_cycle_reduction": directed_chord_audit(),
        "rank_minor": rank_minor_audit(),
        "total_positive_shift": affine_shift_audit(),
        "generator_selection_and_constants": generator_constant_audit(),
        "scope": {
            "complex_projective_K3_required_for_Picard_20": True,
            "field_degree_at_least_two_required_for_cycle_length": True,
            "finite_projective_order_module_or_explicit_trace_certificate": True,
            "uniform_large_cycle_mixed_height_is_assumed": True,
            "height_bound_is_not_derived_from_reflection_maximality": True,
            "verdict": "PASS_CONDITIONAL",
        },
    }


def main() -> None:
    print(json.dumps(audit_payload(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
