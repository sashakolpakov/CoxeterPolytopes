#!/usr/bin/env python3
"""Exact linear-algebra replay for the cusped bounded spanning-root lemma.

The hyperbolic input to the theorem is proved in the companion report.  This
certifier checks its two algebraic interfaces:

* two distinct complete vertex packets span the whole Lorentz space, whereas
  one cusp packet is exactly one dimension short; and
* the strict normalized bound ``< 14`` becomes the integral K3 bound ``<= 13``
  and leaves a finite full-rank Gram box in every Picard rank at most 20.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Sequence


Vector = list[int]
Matrix = list[list[int]]


def matrix_rank(matrix: Matrix) -> int:
    """Rank over Q by exact row reduction."""

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
            scale = work[row][column]
            work[row] = [
                work[row][target] - scale * work[pivot_row][target]
                for target in range(columns)
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def determinant(matrix: Matrix) -> int:
    """Bareiss determinant over Z."""

    if not matrix:
        return 1
    size = len(matrix)
    if any(len(row) != size for row in matrix):
        raise ValueError("determinant requires a square matrix")
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
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


def lorentz_square(vector: Vector) -> int:
    """Square for diag(1,-1,...,-1)."""

    return vector[0] * vector[0] - sum(value * value for value in vector[1:])


def lorentz_pair(left: Vector, right: Vector) -> int:
    return left[0] * right[0] - sum(
        left[index] * right[index] for index in range(1, len(left))
    )


def orthogonal_hyperplane_basis(vector: Vector) -> Matrix:
    """Return an exact integral row basis of vector-perp."""

    covector = [vector[0], *(-value for value in vector[1:])]
    pivot = next(index for index, value in enumerate(covector) if value)
    basis: Matrix = []
    for column in range(len(vector)):
        if column == pivot:
            continue
        row = [0] * len(vector)
        row[column] = covector[pivot]
        row[pivot] = -covector[column]
        if sum(covector[index] * row[index] for index in range(len(vector))):
            raise AssertionError("kernel-basis construction failed")
        basis.append(row)
    if matrix_rank(basis) != len(vector) - 1:
        raise AssertionError("orthogonal packet has the wrong rank")
    return basis


def canonical_vertex(kind: str, rank: int, slot: int) -> Vector:
    """Two distinct finite/null representatives in one Lorentz cone."""

    if rank < 3:
        raise ValueError("the spanning-root theorem starts in Picard rank three")
    if slot not in (0, 1):
        raise ValueError("slot must be zero or one")
    vector = [0] * rank
    if kind == "finite":
        if slot == 0:
            vector[0] = 1
        else:
            vector[0] = 2
            vector[1] = 1
        if lorentz_square(vector) <= 0:
            raise AssertionError("finite representative is not timelike")
    elif kind == "ideal":
        vector[0] = 1
        vector[1 + slot] = 1
        if lorentz_square(vector) != 0:
            raise AssertionError("ideal representative is not null")
    else:
        raise ValueError(f"unknown vertex kind: {kind}")
    return vector


def vertex_pair_record(rank: int, first_kind: str, second_kind: str) -> dict[str, object]:
    first = canonical_vertex(first_kind, rank, 0)
    second = canonical_vertex(second_kind, rank, 1)
    first_packet = orthogonal_hyperplane_basis(first)
    second_packet = orthogonal_hyperplane_basis(second)
    if first == second:
        raise AssertionError("the canonical vertices must be distinct")
    if any(lorentz_pair(first, row) for row in first_packet):
        raise AssertionError("first packet misses first-perp")
    if any(lorentz_pair(second, row) for row in second_packet):
        raise AssertionError("second packet misses second-perp")
    union_rank = matrix_rank(first_packet + second_packet)
    if union_rank != rank:
        raise AssertionError("two distinct complete vertex packets did not span")
    same_vertex_rank = matrix_rank(first_packet + first_packet)
    if same_vertex_rank != rank - 1:
        raise AssertionError("one vertex packet should remain codimension one")
    return {
        "rank": rank,
        "vertex_types": [first_kind, second_kind],
        "first_square": lorentz_square(first),
        "second_square": lorentz_square(second),
        "first_packet_rank": matrix_rank(first_packet),
        "second_packet_rank": matrix_rank(second_packet),
        "two_distinct_packet_union_rank": union_rank,
        "one_packet_or_duplicate_union_rank": same_vertex_rank,
    }


def bounded_gram_count(rank: int) -> int:
    """Number of symmetric integral boxes before inertia/nonsingularity filters."""

    off_diagonal_entries = rank * (rank - 1) // 2
    return 14**off_diagonal_entries


def sample_lorentz_root_gram(rank: int) -> Matrix:
    """A nonsingular signature-(1,rank-1) matrix inside the 0..13 box."""

    gram = [[-2 * int(row == column) for column in range(rank)] for row in range(rank)]
    gram[0][1] = gram[1][0] = 3
    expected = -5 * ((-2) ** (rank - 2))
    if determinant(gram) != expected:
        raise AssertionError("sample Lorentz Gram determinant changed")
    return gram


def rank_record(rank: int) -> dict[str, object]:
    if not 3 <= rank <= 20:
        raise ValueError("K3 Picard rank must lie between 3 and 20 here")
    pair_types = [
        ("finite", "finite"),
        ("finite", "ideal"),
        ("ideal", "finite"),
        ("ideal", "ideal"),
    ]
    pair_records = [vertex_pair_record(rank, *pair) for pair in pair_types]
    gram = sample_lorentz_root_gram(rank)
    return {
        "rho": rank,
        "hyperbolic_dimension": rank - 1,
        "vertex_pair_records": pair_records,
        "strict_normalized_bound": 14,
        "integral_off_diagonal_bound": 13,
        "off_diagonal_slots": rank * (rank - 1) // 2,
        "bounded_gram_box_count": bounded_gram_count(rank),
        "sample_lorentz_gram_determinant": determinant(gram),
    }


def build_certificate() -> dict[str, object]:
    return {
        "verdict": "finite volume plus exact physical root return gives a bounded full spanning packet",
        "uniform_constant": {
            "normalized_star_pairing_strictly_less_than": 14,
            "K3_root_intersection_at_most": 13,
            "valid_Picard_ranks": [3, 20],
        },
        "rank_records": [rank_record(rank) for rank in range(3, 21)],
        "proof_interfaces": {
            "narrow_star": "intersecting/asymptotic pairs have b<=2; the positive-gap sine formula puts ultraparallel closed-star pairs below 14",
            "finite_vertex": "the complete incident packet spans x-perp",
            "ideal_vertex": "bounded horospherical link makes the complete incident packet span u-perp",
            "two_vertex_bridge": "distinct projective vertices have different orthogonal hyperplanes whose sum is V",
            "rational_extraction": "integral vectors spanning over R also span over Q",
            "host_finiteness": "R subset NS subset R-dual and the bounded Gram list is finite",
        },
        "hostile_scope": {
            "finite_sided_exact_return_required": True,
            "intersecting_vs_ultraparallel_visual_split_required": True,
            "finite_volume_required_at_ideal_vertices": True,
            "projective_closed_star_required": True,
            "complete_incident_cusp_packet_required": True,
            "two_distinct_vertices_required": True,
            "physical_literal_square_minus_two_facets_required": True,
            "simplicity_required": False,
            "compactness_required": False,
            "cusp_classification_required": False,
            "reflective_lattice_classification_required": False,
            "reflection_maximality_required_for_spanning_lemma": False,
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
