#!/usr/bin/env python3
"""Certify rational-similarity lineage of the live H10/H12 lattices.

Arithmetic hyperbolic reflection groups defined by rational quadratic forms
are commensurable when the forms are similar over Q.  This script compares
determinant square classes and every local Hasse invariant, using an exact
implementation of the rational Hilbert symbol.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import prod
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "higher_dimensional_rational_lineage_certificate.json"
)


def _e8_cartan() -> sp.Matrix:
    return sp.Matrix(
        [
            [2, 0, -1, 0, 0, 0, 0, 0],
            [0, 2, 0, -1, 0, 0, 0, 0],
            [-1, 0, 2, -1, 0, 0, 0, 0],
            [0, -1, -1, 2, -1, 0, 0, 0],
            [0, 0, 0, -1, 2, -1, 0, 0],
            [0, 0, 0, 0, -1, 2, -1, 0],
            [0, 0, 0, 0, 0, -1, 2, -1],
            [0, 0, 0, 0, 0, 0, -1, 2],
        ]
    )


def _d5_negative() -> sp.Matrix:
    gram = sp.diag(*([-2] * 5))
    for left, right in ((0, 1), (1, 2), (2, 3), (2, 4)):
        gram[left, right] = gram[right, left] = 1
    return gram


def _lattices() -> dict[str, sp.Matrix]:
    u = sp.Matrix([[0, 1], [1, 0]])
    u2 = 2 * u
    e8 = _e8_cartan()
    a2_negative = sp.Matrix([[-2, 1], [1, -2]])
    return {
        "JK-H10-0007": sp.diag(u, -2 * e8, sp.Matrix([[-4]])),
        "K3-D4-H12": sp.diag(
            u, -e8, a2_negative, sp.Matrix([[-2]])
        ),
        "JK-H12": sp.diag(
            u2, _d5_negative(), _d5_negative(), sp.Matrix([[-4]])
        ),
    }


def _diagonalize(matrix: sp.Matrix) -> list[Fraction]:
    """Return an exact rational diagonal form under congruence."""

    work = sp.Matrix(matrix)
    diagonal: list[Fraction] = []
    while work.rows:
        pivot = next((i for i in range(work.rows) if work[i, i] != 0), None)
        if pivot is None:
            pair = next(
                (
                    (i, j)
                    for i in range(work.rows)
                    for j in range(i + 1, work.rows)
                    if work[i, j] != 0
                ),
                None,
            )
            if pair is None:
                raise ArithmeticError("quadratic form is degenerate")
            change = sp.eye(work.rows)
            change[pair[1], pair[0]] = 1
            work = change.T * work * change
            pivot = pair[0]
        order = [pivot] + [i for i in range(work.rows) if i != pivot]
        work = work.extract(order, order)
        value = sp.Rational(work[0, 0])
        diagonal.append(Fraction(int(value.p), int(value.q)))
        if work.rows == 1:
            break
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / value
    return diagonal


def _prime_factors(number: int) -> list[int]:
    return [int(prime) for prime in sp.factorint(abs(number))]


def _square_class(value: Fraction) -> int:
    """Signed squarefree integer representing a nonzero Q*/Q*2 class."""

    if value == 0:
        raise ArithmeticError("zero has no square class")
    result = -1 if value < 0 else 1
    numerator = abs(value.numerator)
    denominator = value.denominator
    valuations: dict[int, int] = {}
    for prime, exponent in sp.factorint(numerator).items():
        valuations[int(prime)] = valuations.get(int(prime), 0) + int(exponent)
    for prime, exponent in sp.factorint(denominator).items():
        valuations[int(prime)] = valuations.get(int(prime), 0) - int(exponent)
    for prime, exponent in valuations.items():
        if exponent % 2:
            result *= prime
    return result


def _legendre(unit: int, prime: int) -> int:
    residue = pow(unit % prime, (prime - 1) // 2, prime)
    if residue == 1:
        return 1
    if residue == prime - 1:
        return -1
    raise ArithmeticError("Legendre symbol requested for a nonunit")


def _hilbert(left: Fraction, right: Fraction, place: int | str) -> int:
    a = _square_class(left)
    b = _square_class(right)
    if place == "infinity":
        return -1 if a < 0 and b < 0 else 1
    prime = int(place)
    if prime == 2:
        alpha = int(a % 2 == 0)
        beta = int(b % 2 == 0)
        unit_a = a // (2 if alpha else 1)
        unit_b = b // (2 if beta else 1)
        exponent = (
            ((unit_a - 1) // 2) * ((unit_b - 1) // 2)
            + alpha * ((unit_b * unit_b - 1) // 8)
            + beta * ((unit_a * unit_a - 1) // 8)
        )
        return -1 if exponent % 2 else 1
    alpha = int(a % prime == 0)
    beta = int(b % prime == 0)
    unit_a = a // (prime if alpha else 1)
    unit_b = b // (prime if beta else 1)
    sign = -1 if alpha * beta * ((prime - 1) // 2) % 2 else 1
    if beta:
        sign *= _legendre(unit_a, prime)
    if alpha:
        sign *= _legendre(unit_b, prime)
    return sign


def _hasse(diagonal: list[Fraction], place: int | str) -> int:
    return prod(
        _hilbert(diagonal[left], diagonal[right], place)
        for left in range(len(diagonal))
        for right in range(left + 1, len(diagonal))
    )


def _relevant_primes(*diagonals: list[Fraction]) -> list[int]:
    primes = {2}
    for diagonal in diagonals:
        for value in diagonal:
            square_class = _square_class(value)
            primes.update(_prime_factors(square_class))
    return sorted(primes)


def _compare(
    lattice: sp.Matrix,
    reference: sp.Matrix,
    similarity_scalar: Fraction,
) -> dict[str, Any]:
    lattice_diagonal = _diagonalize(lattice)
    reference_diagonal = [
        similarity_scalar * value for value in _diagonalize(reference)
    ]
    places: list[int | str] = _relevant_primes(
        lattice_diagonal, reference_diagonal
    ) + ["infinity"]
    local = {
        str(place): {
            "lattice_hasse": _hasse(lattice_diagonal, place),
            "reference_hasse": _hasse(reference_diagonal, place),
        }
        for place in places
    }
    determinant_ratio = Fraction(
        int(lattice.det()), int(reference.det())
    )
    scaled_determinant_ratio = determinant_ratio / (
        similarity_scalar ** lattice.rows
    )
    determinant_square_class_matches = _square_class(scaled_determinant_ratio) == 1
    all_hasse_match = all(
        values["lattice_hasse"] == values["reference_hasse"]
        for values in local.values()
    )
    return {
        "rank": lattice.rows,
        "similarity_scalar": str(similarity_scalar),
        "lattice_determinant": int(lattice.det()),
        "reference_determinant": int(reference.det()),
        "scaled_determinant_ratio": str(scaled_determinant_ratio),
        "scaled_determinant_ratio_is_square": determinant_square_class_matches,
        "lattice_diagonal": [str(value) for value in lattice_diagonal],
        "scaled_reference_diagonal": [str(value) for value in reference_diagonal],
        "local_hasse_invariants": local,
        "all_local_hasse_invariants_match": all_hasse_match,
        "rationally_similar": determinant_square_class_matches and all_hasse_match,
    }


def build_certificate() -> dict[str, Any]:
    lattices = _lattices()
    odd_unimodular_h10 = sp.diag(1, *([-1] * 10))
    odd_unimodular_h12 = sp.diag(1, *([-1] * 12))
    mcleod_h10 = sp.diag(3, *([-1] * 10))
    mcleod_h12 = sp.diag(3, *([-1] * 12))

    comparisons = {
        "JK-H10-0007_vs_odd_unimodular_H10": _compare(
            lattices["JK-H10-0007"], odd_unimodular_h10, Fraction(1024)
        ),
        "JK-H10-0007_vs_McLeod_f3_H10": _compare(
            lattices["JK-H10-0007"], mcleod_h10, Fraction(1024, 3)
        ),
        "JK-H12_vs_odd_unimodular_H12": _compare(
            lattices["JK-H12"], odd_unimodular_h12, Fraction(256)
        ),
        "JK-H12_vs_McLeod_f3_H12": _compare(
            lattices["JK-H12"], mcleod_h12, Fraction(256, 3)
        ),
        "K3-D4-H12_vs_odd_unimodular_H12": _compare(
            lattices["K3-D4-H12"], odd_unimodular_h12, Fraction(6)
        ),
        "K3-D4-H12_vs_McLeod_f3_H12": _compare(
            lattices["K3-D4-H12"], mcleod_h12, Fraction(2)
        ),
    }
    expected = {
        "JK-H10-0007_vs_odd_unimodular_H10": True,
        "JK-H10-0007_vs_McLeod_f3_H10": False,
        "JK-H12_vs_odd_unimodular_H12": True,
        "JK-H12_vs_McLeod_f3_H12": True,
        "K3-D4-H12_vs_odd_unimodular_H12": False,
        "K3-D4-H12_vs_McLeod_f3_H12": False,
    }
    checks = {
        key: comparisons[key]["rationally_similar"] == value
        for key, value in expected.items()
    }
    return {
        "schema": "higher_dimensional_rational_lineage_certificate_v1",
        "criterion": {
            "statement": (
                "Nondegenerate rational forms are isometric iff dimension, "
                "determinant square class, real signature, and all local Hasse "
                "invariants agree; arithmetic orthogonal groups of similar "
                "forms are commensurable."
            ),
            "reference": (
                "O'Meara, Introduction to Quadratic Forms; "
                "Maclachlan--Reid, Arithmetic of Hyperbolic 3-Manifolds, Ch. 6"
            ),
        },
        "sign_convention": (
            "signature (1,n); McLeod's -3*x0^2+sum(xi^2) is multiplied by -1"
        ),
        "comparisons": comparisons,
        "conclusions": {
            "JK-H10-0007": (
                "same rational quadratic space and arithmetic commensurability "
                "class as Vinberg's odd-unimodular H10 chamber"
            ),
            "JK-H12": (
                "same arithmetic commensurability class as both the odd-unimodular "
                "H12 chamber and McLeod's f3^12 chamber"
            ),
            "K3-D4-H12": (
                "not rationally similar to either of those two displayed classical "
                "forms; its prior lineage is instead the determinant-6 "
                "Scharlau--Kirschmer census entry"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archived rational-lineage certificate is stale")
        print("archived higher-dimensional rational-lineage certificate verified")
    else:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)


if __name__ == "__main__":
    main()
