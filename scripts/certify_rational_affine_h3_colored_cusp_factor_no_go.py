#!/usr/bin/env python3
"""Exact certificate for the affine-H3 colored cusp-factor local no-go."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_affine_h3_colored_cusp_factor_no_go.md"


def dot(left: list[int], right: list[int]) -> int:
    return sum(x * y for x, y in zip(left, right))


def exact_rank(matrix: list[list[int]]) -> int:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if work else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            coefficient = work[row][column]
            work[row] = [
                entry - coefficient * pivot_entry
                for entry, pivot_entry in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def mat_vec(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [
        sum(entry * coefficient for entry, coefficient in zip(row, vector))
        for row in matrix
    ]


def affine_d_normals(rank: int) -> list[list[int]]:
    """Euclidean normals alpha_0,...,alpha_n for affine D_n."""
    assert rank >= 4
    zero = [0] * rank
    zero[0] = zero[1] = -1
    roots = [zero]
    for index in range(rank - 1):
        root = [0] * rank
        root[index] = 1
        root[index + 1] = -1
        roots.append(root)
    last = [0] * rank
    last[-2] = last[-1] = 1
    roots.append(last)
    return roots


def affine_marks(rank: int) -> list[int]:
    return [1, 1] + [2] * (rank - 3) + [1, 1]


def affine_gram(rank: int) -> list[list[int]]:
    normals = affine_d_normals(rank)
    return [[-dot(left, right) for right in normals] for left in normals]


def physical_types(rank: int) -> list[str]:
    assert rank % 2 == 0
    types = []
    for vertex in range(rank + 1):
        is_boundary = (
            vertex in {0, 1, rank - 1, rank}
            or (3 <= vertex <= rank - 3 and vertex % 2 == 1)
        )
        types.append("B" if is_boundary else "E")
    return types


def physical_multiplicities(rank: int) -> list[int]:
    return [
        mark if kind == "B" else 2 * mark
        for mark, kind in zip(affine_marks(rank), physical_types(rank))
    ]


def physical_intersection_gram(rank: int) -> list[list[int]]:
    root_gram = affine_gram(rank)
    kinds = physical_types(rank)
    matrix = [[0] * (rank + 1) for _ in range(rank + 1)]
    for i in range(rank + 1):
        matrix[i][i] = -4 if kinds[i] == "B" else -1
        for j in range(i + 1, rank + 1):
            if root_gram[i][j] == 1:
                matrix[i][j] = matrix[j][i] = 1
    return matrix


def proximity_clusters(rank: int) -> list[tuple[int, ...]]:
    """Downstairs branch multiplicities at the contracted E vertices."""
    kinds = physical_types(rank)
    marks = affine_marks(rank)
    physical = physical_intersection_gram(rank)
    return [
        tuple(
            marks[other]
            for other, kind in enumerate(kinds)
            if kind == "B" and physical[vertex][other] == 1
        )
        for vertex, kind in enumerate(kinds)
        if kind == "E"
    ]


def canonical_colors(rank: int) -> list[int]:
    return [2 if kind == "B" else -2 for kind in physical_types(rank)]


def coordinate_refinement(rank: int) -> list[int]:
    vector = [0] * rank
    vector[-1] = 2
    return vector


def finite_simple_roots(rank: int) -> list[list[int]]:
    return affine_d_normals(rank)[1:]


def coordinate_coefficients(rank: int) -> list[int]:
    # 2e_n = alpha_n - alpha_{n-1}.
    return [0] * (rank - 2) + [-1, 1]


def certify_rank(rank: int) -> dict[str, object]:
    assert rank >= 4 and rank % 2 == 0
    roots = affine_d_normals(rank)
    gram = affine_gram(rank)
    marks = affine_marks(rank)
    kinds = physical_types(rank)
    multiplicities = physical_multiplicities(rank)
    physical = physical_intersection_gram(rank)
    colors = canonical_colors(rank)

    # The displayed Euclidean vectors prove negative semidefiniteness.
    assert all(dot(root, root) == 2 for root in roots)
    assert exact_rank(gram) == rank
    assert mat_vec(gram, marks) == [0] * (rank + 1)
    null_normal = [
        sum(mark * root[column] for mark, root in zip(marks, roots))
        for column in range(rank)
    ]
    assert null_normal == [0] * rank

    # Exact physical divisor equations and adjunction.
    assert mat_vec(physical, multiplicities) == [0] * (rank + 1)
    assert sum(1 for kind in kinds if kind == "B") == rank // 2 + 2
    assert sum(1 for kind in kinds if kind == "E") == rank // 2 - 1
    positive_mass = sum(
        mark for mark, kind in zip(marks, kinds) if kind == "B"
    )
    negative_mass = sum(
        mark for mark, kind in zip(marks, kinds) if kind == "E"
    )
    assert positive_mass == rank
    assert negative_mass == rank - 2
    assert sum(mark * color for mark, color in zip(marks, colors)) == 4

    canonical_intersections = [
        2 if kind == "B" else -1 for kind in kinds
    ]
    assert (
        sum(
            coefficient * canonical
            for coefficient, canonical in zip(
                multiplicities, canonical_intersections
            )
        )
        == 4
    )
    physical_gcd = gcd(*multiplicities)
    assert physical_gcd == 1

    # The coordinate mirror is primitive, norm four, zero-colored, and
    # changes one sign.  It pairs evenly with the D_n root sublattice.
    coordinate = coordinate_refinement(rank)
    finite_roots = finite_simple_roots(rank)
    coefficients = coordinate_coefficients(rank)
    reconstructed = [
        sum(
            coefficient * root[column]
            for coefficient, root in zip(coefficients, finite_roots)
        )
        for column in range(rank)
    ]
    assert reconstructed == coordinate
    assert dot(coordinate, coordinate) == 4
    assert all(dot(coordinate, root) % 2 == 0 for root in finite_roots)
    finite_colors = colors[1:]
    assert sum(
        coefficient * color
        for coefficient, color in zip(coefficients, finite_colors)
    ) == 0

    sample = list(range(rank))
    reflected = sample.copy()
    reflected[-1] *= -1
    assert sum(sample) % 2 == sum(reflected) % 2

    # Proximity degrees: leaf B has one E neighbor, internal B has two;
    # the blowdown squares are consequently -3 and -2.
    degrees = [
        sum(
            physical[index][other]
            for other in range(rank + 1)
            if other != index
        )
        for index in range(rank + 1)
    ]
    for kind, mark, degree in zip(kinds, marks, degrees):
        if kind != "B":
            continue
        assert degree == (1 if mark == 1 else 2)
        assert -4 + degree == (-3 if mark == 1 else -2)

    clusters = proximity_clusters(rank)
    assert all(sum(cluster) == 4 for cluster in clusters)
    if rank == 4:
        assert clusters == [(1, 1, 1, 1)]
    else:
        assert sorted(len(cluster) for cluster in clusters).count(3) == 2
        assert all(len(cluster) in {2, 3} for cluster in clusters)

    return {
        "affine_rank": rank,
        "vertices": rank + 1,
        "boundary_vertices": kinds.count("B"),
        "exceptional_vertices": kinds.count("E"),
        "positive_mark_mass": positive_mass,
        "negative_mark_mass": negative_mass,
        "canonical_mark_sum": 4,
        "physical_gcd_witness": physical_gcd,
        "proximity_cluster_arities": [len(cluster) for cluster in clusters],
        "coordinate_refinement_norm": 4,
        "coordinate_refinement_color": 0,
    }


def certificate() -> dict[str, object]:
    rows = [certify_rank(rank) for rank in range(4, 42, 2)]
    report_text = REPORT.read_text()
    required = [
        "sum_{\\mathcal B}a_i=n",
        "sum_{\\mathcal E}a_i=n-2",
        "q_M(K_S,v)=0",
        "q_M(M,v)\\subset2\\mathbf Z",
        "B_n-B_{n-1}\\in2\\operatorname {NS}(S)",
        "v^\\perp=\\{x_n=0\\}",
        "1+1+1+1=4",
        "must still",
        "projective rational surface",
    ]
    assert all(token in report_text for token in required)
    return {
        "status": "pass",
        "theorem": "colored_affine_factor_local_no_go",
        "tested_even_ranks": [row["affine_rank"] for row in rows],
        "rows": rows,
        "scope": {
            "proved": "marked affine packets and adapted-host parity gate",
            "not_proved": [
                "arbitrary-rank projective rational genus-three realization",
                "reflection-maximal complete-nef counterfamily",
                "unbounded zero-color chains under H1/H2 exclusion",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compact", action="store_true")
    arguments = parser.parse_args()
    payload = certificate()
    print(
        json.dumps(
            payload,
            sort_keys=True,
            indent=None if arguments.compact else 2,
        )
    )


if __name__ == "__main__":
    main()
