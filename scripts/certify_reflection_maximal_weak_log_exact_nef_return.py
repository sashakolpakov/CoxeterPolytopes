#!/usr/bin/env python3
"""Exact replay of the weak-log exact-nef-return F_4 calibration."""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import TypeAlias


Q: TypeAlias = Fraction
Vector: TypeAlias = tuple[Q, Q, Q]
Matrix: TypeAlias = list[list[Q]]

L_GRAM: Matrix = [
    [Q(-4), Q(1), Q(0)],
    [Q(1), Q(0), Q(0)],
    [Q(0), Q(0), Q(-1)],
]
M_GRAM: Matrix = [
    [Q(-2), Q(1), Q(0)],
    [Q(1), Q(0), Q(0)],
    [Q(0), Q(0), Q(-2)],
]
IDENTITY: Matrix = [
    [Q(1), Q(0), Q(0)],
    [Q(0), Q(1), Q(0)],
    [Q(0), Q(0), Q(1)],
]


def dot(x: Vector, y: Vector, gram: Matrix) -> Q:
    return sum(
        (x[i] * gram[i][j] * y[j] for i in range(3) for j in range(3)),
        Q(0),
    )


def matmul(a: Matrix, b: Matrix) -> Matrix:
    return [
        [sum((a[i][k] * b[k][j] for k in range(3)), Q(0)) for j in range(3)]
        for i in range(3)
    ]


def transpose(a: Matrix) -> Matrix:
    return [list(row) for row in zip(*a)]


def reflection_matrix(root: Vector, gram: Matrix) -> Matrix:
    """Return the row-major matrix of reflection in a square -2 root."""

    assert dot(root, root, gram) == -2
    columns: list[Vector] = []
    for j in range(3):
        basis = tuple(Q(1 if i == j else 0) for i in range(3))
        pairing = dot(basis, root, gram)
        image = tuple(basis[i] + pairing * root[i] for i in range(3))
        columns.append(image)
    return [list(row) for row in zip(*columns)]


def encode(value: object) -> object:
    if isinstance(value, Fraction):
        return [value.numerator, value.denominator]
    if isinstance(value, tuple):
        return [encode(item) for item in value]
    if isinstance(value, list):
        return [encode(item) for item in value]
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    return value


def payload() -> dict[str, object]:
    B: Vector = (Q(1), Q(0), Q(0))
    F: Vector = (Q(0), Q(1), Q(0))
    E: Vector = (Q(0), Q(0), Q(1))
    G: Vector = (Q(0), Q(1), Q(-1))
    K: Vector = (Q(-2), Q(-6), Q(1))
    D: Vector = (Q(3, 2), Q(6), Q(-1))

    assert dot(K, K, L_GRAM) == 7
    assert tuple(-K[i] - B[i] / 2 for i in range(3)) == D
    curves = {"B": B, "E": E, "G": G}
    assert [dot(curves[name], curves[name], L_GRAM) for name in curves] == [
        -4,
        -1,
        -1,
    ]
    assert [dot(K, curves[name], L_GRAM) for name in curves] == [2, -1, -1]
    assert [
        dot(curves[name], curves[name], L_GRAM)
        + dot(K, curves[name], L_GRAM)
        for name in curves
    ] == [-2, -2, -2]
    d_pairings = [dot(D, curves[name], L_GRAM) for name in curves]
    assert d_pairings == [0, 1, Q(1, 2)]
    assert dot(D, D, L_GRAM) == 8

    # The complete toric Mori cone is Cone(B,E,G).  These three exact
    # nonnegative pairings are therefore an independent nef certificate.
    positive_section_strict = tuple(B[i] + 4 * G[i] + 3 * E[i] for i in range(3))
    assert positive_section_strict == (Q(1), Q(4), Q(-1))
    assert F == tuple(G[i] + E[i] for i in range(3))

    roots: dict[str, Vector] = {
        "B_half": (Q(1), Q(0), Q(0)),
        "E": (Q(0), Q(0), Q(1)),
        "G": (Q(0), Q(1), Q(-1)),
    }
    order = ("B_half", "E", "G")
    root_gram = [[dot(roots[a], roots[b], M_GRAM) for b in order] for a in order]
    assert root_gram == [
        [Q(-2), Q(0), Q(1)],
        [Q(0), Q(-2), Q(2)],
        [Q(1), Q(2), Q(-2)],
    ]

    reflections: dict[str, Matrix] = {}
    for name in order:
        matrix = reflection_matrix(roots[name], M_GRAM)
        assert all(entry.denominator == 1 for row in matrix for entry in row)
        assert matmul(matrix, matrix) == IDENTITY
        assert matmul(matmul(transpose(matrix), M_GRAM), matrix) == M_GRAM
        reflections[name] = matrix

    assert matmul(
        matmul(reflections["B_half"], reflections["E"]),
        matmul(reflections["B_half"], reflections["E"]),
    ) == IDENTITY
    bg = matmul(reflections["B_half"], reflections["G"])
    assert matmul(matmul(bg, bg), bg) == IDENTITY
    eg = matmul(reflections["E"], reflections["G"])
    assert eg != IDENTITY and matmul(eg, eg) != IDENTITY

    raw: dict[str, object] = {
        "schema": "reflection_maximal_weak_log_exact_nef_return/v1",
        "surface": "Bl_p(F_4), p torus-fixed on the positive section and p not in B",
        "lattice_basis": ["B", "F", "E"],
        "lattice_gram": L_GRAM,
        "canonical": K,
        "canonical_square": dot(K, K, L_GRAM),
        "curves": curves,
        "curve_squares": {name: dot(curve, curve, L_GRAM) for name, curve in curves.items()},
        "canonical_pairings": {
            name: dot(K, curve, L_GRAM) for name, curve in curves.items()
        },
        "weak_log_class": D,
        "weak_log_square": dot(D, D, L_GRAM),
        "weak_log_pairings_on_complete_toric_mori_generators": {
            name: dot(D, curve, L_GRAM) for name, curve in curves.items()
        },
        "complete_toric_mori_generators": ["B", "E", "G"],
        "adapted_basis": ["B_half", "F", "E"],
        "adapted_doubled_gram": M_GRAM,
        "facet_roots": roots,
        "facet_root_gram": root_gram,
        "reflection_matrices": reflections,
        "coxeter_orders": {"B_half_E": 2, "B_half_G": 3, "E_G": "infinity"},
        "area_over_pi": Q(1, 6),
        "proper_overgroup_area_upper_over_pi": Q(1, 12),
        "minimum_noncompact_coxeter_area_over_pi": Q(1, 6),
        "conclusion": (
            "the independently certified nef weak-log class and adapted-host "
            "reflection maximality return the complete nef chamber"
        ),
    }
    result = encode(raw)
    assert isinstance(result, dict)
    canonical = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    result["payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    result = payload()
    default = Path(
        "data/enumerations/reflection_maximal_weak_log_exact_nef_return.json"
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
