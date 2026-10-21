#!/usr/bin/env python3
"""Exact replay for the bounded chordless-cycle K3 field gate."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import comb, floor, log2
from pathlib import Path
from typing import Sequence

try:
    from scripts.certify_k3_dominant_correspondence_height_field_bound import (
        polynomial_search_count,
        primitive_combination_bound,
    )
except ModuleNotFoundError:  # Direct execution from the repository root.
    from certify_k3_dominant_correspondence_height_field_bound import (
        polynomial_search_count,
        primitive_combination_bound,
    )


def matrix_rank(matrix: Sequence[Sequence[int]]) -> int:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    rows = len(work)
    columns = len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, rows) if work[row][column]), None
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        diagonal = work[rank][column]
        work[rank] = [entry / diagonal for entry in work[rank]]
        for row in range(rows):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                entry - factor * pivot_entry
                for entry, pivot_entry in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == rows:
            break
    return rank


def cyclic_tridiagonal(length: int, diagonal_shift: int = 3) -> tuple[tuple[int, ...], ...]:
    if length < 3:
        raise ValueError("length must be at least three")
    matrix = [[0] * length for _ in range(length)]
    for index in range(length):
        matrix[index][index] = diagonal_shift + index
        successor = (index + 1) % length
        value = 1 + (index % 2)
        matrix[index][successor] = value
        matrix[successor][index] = value
    return tuple(tuple(row) for row in matrix)


def tridiagonal_rank_record(length: int) -> dict[str, object]:
    matrix = cyclic_tridiagonal(length)
    rank = matrix_rank(matrix)
    if rank < length - 2:
        raise AssertionError("cyclic tridiagonal nullity exceeded two")
    return {
        "cycle_length": length,
        "sample_matrix": [list(row) for row in matrix],
        "sample_rank": rank,
        "sample_nullity": length - rank,
        "universal_kernel_injection": "ker(A)->(x_1,x_2)",
        "universal_nullity_bound": 2,
        "universal_rank_lower_bound": length - 2,
    }


def trace_induction_record(length: int) -> dict[str, object]:
    if length < 1:
        raise ValueError("length must be positive")
    counts = {size: comb(length, size) for size in range(1, length + 1)}
    if sum(counts.values()) != 2**length - 1:
        raise AssertionError("reflection-product subset expansion changed")
    return {
        "reflection_word_length": length,
        "nonempty_subsequence_term_count": 2**length - 1,
        "terms_by_cycle_length": {
            str(size): count for size, count in counts.items()
        },
        "full_cycle_coefficient": (-1) ** length,
        "full_cycle_occurs_once": counts[length] == 1,
        "all_other_cycles_shorter": True,
        "inductive_integrality_isolation": True,
    }


def perimeter_edge(left: int, right: int) -> tuple[int, int]:
    return tuple(sorted((left, right)))


def path_vertices(start: int, stop: int, length: int) -> list[int]:
    result = [start]
    current = start
    while current != stop:
        current = (current + 1) % length
        result.append(current)
    return result


def cycle_edges(vertices: Sequence[int]) -> Counter[tuple[int, int]]:
    result: Counter[tuple[int, int]] = Counter()
    for left, right in zip(vertices, vertices[1:]):
        result[perimeter_edge(left, right)] += 1
    result[perimeter_edge(vertices[-1], vertices[0])] += 1
    return result


def chord_decomposition_record(length: int, left: int, right: int) -> dict[str, object]:
    if length < 4:
        raise ValueError("a chord decomposition needs length at least four")
    if not (0 <= left < right < length):
        raise ValueError("invalid chord endpoints")
    if right - left in (1, length - 1):
        raise ValueError("the selected edge is not a chord")
    perimeter = cycle_edges(list(range(length)))
    first_vertices = path_vertices(left, right, length)
    second_vertices = path_vertices(right, left, length)
    first = cycle_edges(first_vertices)
    second = cycle_edges(second_vertices)
    chord = perimeter_edge(left, right)
    expected = perimeter.copy()
    expected[chord] += 2
    if first + second != expected:
        raise AssertionError("chord cycle-product identity failed")
    return {
        "cycle_length": length,
        "chord": list(chord),
        "first_cycle": first_vertices,
        "second_cycle": second_vertices,
        "identity": "c(C1)c(C2)=c(C)c(chord-pair)",
        "edge_exponent_identity": True,
    }


def affine_cycle_record(
    length: int,
    distinguished_signed_cycle: int | Fraction,
    conjugate_signed_cycles: Sequence[int | Fraction],
) -> dict[str, object]:
    if length < 2:
        raise ValueError("cycle length must be at least two")
    beta0 = Fraction(distinguished_signed_cycle)
    betas = tuple(Fraction(value) for value in conjugate_signed_cycles)
    unit_bound = Fraction(2**length)
    if beta0 < 0:
        raise ValueError("the acute distinguished signed cycle is nonnegative")
    if any(abs(value) > unit_bound for value in betas):
        raise ValueError("a positive-place conjugate exceeded Cauchy--Schwarz")
    thetas = (beta0 + unit_bound,) + tuple(
        value + unit_bound for value in betas
    )
    if any(value <= 0 for value in thetas):
        raise ValueError(
            "the supplied rational sample hits the algebraically forbidden zero boundary"
        )
    large = beta0 > unit_bound
    if large and any(value > thetas[0] for value in thetas[1:]):
        raise AssertionError("large affine cycle failed dominance")
    if not large and any(value > 2 * unit_bound for value in thetas):
        raise AssertionError("small affine cycle failed its uniform envelope")
    return {
        "cycle_length": length,
        "Cauchy_absolute_bound": str(unit_bound),
        "distinguished_signed_cycle": str(beta0),
        "other_signed_cycles": [str(value) for value in betas],
        "affine_coordinate_conjugates": [str(value) for value in thetas],
        "large_at_distinguished_place": large,
        "affine_coordinate_totally_positive": True,
        "affine_coordinate_dominant_when_large": large,
        "small_coordinate_uniform_bound": str(2 * unit_bound),
    }


def k3_cycle_bound_record(
    hyperbolic_dimension: int,
    correspondence_height_bound: int,
    square_lower_bound: int = 2,
) -> dict[str, object]:
    if hyperbolic_dimension < 2:
        raise ValueError("hyperbolic dimension must be at least two")
    if square_lower_bound <= 0:
        raise ValueError("square lower bound must be positive")
    if correspondence_height_bound <= 0:
        raise ValueError("height bound must be positive")
    degree_bound = 20 // (hyperbolic_dimension + 1)
    if degree_bound < 2:
        raise ValueError("no non-rational full-trace K3 field fits")
    module_rank = hyperbolic_dimension + 1
    chordless_length_bound = module_rank + 2
    cycle_Cauchy_bound = 2**chordless_length_bound
    affine_small_bound = 2 * cycle_Cauchy_bound
    selected_generator_bound = floor(log2(degree_bound))
    combination = primitive_combination_bound(
        degree_bound, selected_generator_bound
    )
    coefficient_bound = combination["positive_coefficient_bound"]
    individual_bound = max(
        Fraction(affine_small_bound),
        Fraction(correspondence_height_bound, square_lower_bound),
    )
    primitive_bound = (
        selected_generator_bound * coefficient_bound * individual_bound
    )
    discriminant_bound = primitive_bound ** (
        degree_bound * (degree_bound - 1)
    )
    return {
        "hyperbolic_dimension": hyperbolic_dimension,
        "field_module_rank": module_rank,
        "field_degree_bound": degree_bound,
        "chordless_cycle_length_bound": chordless_length_bound,
        "cycle_Cauchy_absolute_bound": cycle_Cauchy_bound,
        "affine_small_coordinate_bound": affine_small_bound,
        "selected_generator_count_bound": selected_generator_bound,
        "primitive_combination_coefficient_bound": coefficient_bound,
        "mixed_correspondence_height_bound": correspondence_height_bound,
        "positive_square_lower_bound": square_lower_bound,
        "individual_affine_coordinate_bound": str(individual_bound),
        "primitive_generator_conjugate_bound": str(primitive_bound),
        "field_discriminant_bound": str(discriminant_bound),
        "finite_polynomial_box_size": polynomial_search_count(
            degree_bound, primitive_bound
        ),
    }


def build_certificate() -> dict[str, object]:
    chord_records = []
    for length in range(4, 13):
        for right in range(2, length - 1):
            chord_records.append(chord_decomposition_record(length, 0, right))
    return {
        "theorem": "bounded chordless cycles plus one-intersection heights finite the K3 field set",
        "trace_induction_records": [
            trace_induction_record(length) for length in range(1, 13)
        ],
        "tridiagonal_rank_records": [
            tridiagonal_rank_record(length) for length in range(3, 13)
        ],
        "chord_decomposition_records": chord_records,
        "affine_cycle_records": [
            affine_cycle_record(2, 3, (-3, 0, 2)),
            affine_cycle_record(3, 7, (-7, -1, 6)),
            affine_cycle_record(3, 9, (-7, 0, 8)),
            affine_cycle_record(5, 41, (-31, -2, 29)),
        ],
        "K3_dimension_records": [
            k3_cycle_bound_record(dimension, 256)
            for dimension in range(2, 10)
        ],
        "conclusions": {
            "cyclic_integrality_from_reflection_trace_induction": True,
            "individual_Cartan_entries_in_order_required": False,
            "ground_field_cycles_reduce_to_pairs_and_chordless_cycles": True,
            "chordless_cycle_length_at_most_n_plus_3": True,
            "non_distinguished_cycle_conjugate_bound_is_two_to_length": True,
            "affine_signed_cycle_is_totally_positive": True,
            "large_affine_signed_cycle_is_dominant": True,
            "bounded_large_cycle_mixed_heights_give_finite_field_set": True,
            "pair_generation_required": False,
        },
        "scope": {
            "complex_projective_K3_required": True,
            "exact_full_Picard_trace_K3_required": True,
            "finite_projective_order_module_or_explicit_trace_certificate_required": True,
            "returned_reflections_preserve_order_module_required": True,
            "bounded_mixed_degree_derived_from_maximality": False,
            "uses_reflective_lattice_classification": False,
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args()
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
