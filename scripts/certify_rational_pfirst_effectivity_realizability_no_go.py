#!/usr/bin/env python3
"""Exact replay of the rank-three P-first realizability obstruction."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path


J = (1, -1, -1)


def dot(x: tuple[int, ...], y: tuple[int, ...]) -> int:
    return sum(s * a * b for s, a, b in zip(J, x, y))


def det3(columns: tuple[tuple[int, ...], ...]) -> int:
    a, b, c = columns
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - b[0] * (a[1] * c[2] - a[2] * c[1])
        + c[0] * (a[1] * b[2] - a[2] * b[1])
    )


def reflection_matrix(root: tuple[int, ...]) -> list[list[int]]:
    square = dot(root, root)
    out: list[list[int]] = []
    for j in range(3):
        basis = tuple(1 if i == j else 0 for i in range(3))
        numerator = 2 * dot(basis, root)
        assert numerator % square == 0
        image = tuple(
            basis[i] - (numerator // square) * root[i] for i in range(3)
        )
        out.append(list(image))
    # Store ordinary row-major matrices rather than image columns.
    return [list(row) for row in zip(*out)]


def matmul(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    return [
        [sum(a[i][k] * b[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def transpose(a: list[list[int]]) -> list[list[int]]:
    return [list(row) for row in zip(*a)]


def payload() -> dict[str, object]:
    K = (-3, 1, 1)
    physical = {
        "R": (0, 1, -1),
        "B1": (0, -2, 0),
        "B2": (-2, 2, 2),
    }
    primitive = {
        "R": physical["R"],
        "B1": (0, -1, 0),
        "B2": (-1, 1, 1),
    }
    order = ("R", "B1", "B2")
    h = (-4, 1, 2)
    vertices = {
        "R_B1": (-1, 0, 0),
        "R_B2": (-2, 1, 1),
        "B1_B2": (-1, 0, 1),
    }

    assert all(k % 2 for k in K)
    assert dot(K, K) == 7 == 10 - 3
    assert det3(tuple(primitive[name] for name in order)) == 1

    physical_gram = [
        [dot(physical[a], physical[b]) for b in order] for a in order
    ]
    assert physical_gram == [[-2, 2, 0], [2, -4, 4], [0, 4, -4]]
    assert [dot(K, physical[name]) for name in order] == [0, 2, 2]
    assert [dot(h, physical[name]) for name in order] == [1, 2, 2]
    assert dot(h, h) == 11

    vertex_records: dict[str, object] = {}
    incident = {
        "R_B1": ("R", "B1"),
        "R_B2": ("R", "B2"),
        "B1_B2": ("B1", "B2"),
    }
    for label, v in vertices.items():
        zeros = incident[label]
        assert all(dot(v, physical[name]) == 0 for name in zeros)
        assert all(dot(v, physical[name]) >= 0 for name in order)
        vertex_records[label] = {
            "vector": list(v),
            "square": dot(v, v),
            "incident": list(zeros),
        }
    assert sorted(record["square"] for record in vertex_records.values()) == [0, 1, 2]

    metric = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    reflection_records: dict[str, object] = {}
    for name in order:
        matrix = reflection_matrix(primitive[name])
        assert matmul(matrix, matrix) == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        assert matmul(matmul(transpose(matrix), metric), matrix) == metric
        reflection_records[name] = matrix

    rr = {str(q): 2 - q for q in (1, 2, 4)}
    assert rr == {"1": 1, "2": 0, "4": -2}

    result: dict[str, object] = {
        "schema": "rational_pfirst_effectivity_realizability_no_go/v1",
        "lattice": "I_(1,2)",
        "canonical": list(K),
        "canonical_square": dot(K, K),
        "physical_walls": {name: list(physical[name]) for name in order},
        "primitive_roots": {name: list(primitive[name]) for name in order},
        "primitive_basis_determinant": 1,
        "physical_gram": physical_gram,
        "canonical_pairings": [0, 2, 2],
        "strict_interior": list(h),
        "strict_interior_square": 11,
        "strict_slacks": [1, 2, 2],
        "vertices": vertex_records,
        "coxeter_orders": {"R_B1": 4, "R_B2": 2, "B1_B2": "infinity"},
        "area_over_pi": [1, 4],
        "proper_overgroup_area_upper_over_pi": [1, 8],
        "minimum_noncompact_coxeter_area_over_pi": [1, 6],
        "reflection_matrices": reflection_records,
        "riemann_roch_chi_by_negative_square": rr,
        "conclusion": (
            "the characteristic canonical vector lies in the proposed chamber, "
            "so no marked smooth projective rational surface can have this full nef cone"
        ),
    }
    canonical = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    result["payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    result = payload()
    default = Path("data/enumerations/rational_pfirst_effectivity_realizability_no_go.json")
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

