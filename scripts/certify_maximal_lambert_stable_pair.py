#!/usr/bin/env python3
"""Exact algebraic checks for the maximal Lambert stable-pair obstruction."""

from __future__ import annotations

import json
from fractions import Fraction

import sympy as sp


def triangle_candidates(index: int) -> list[tuple[int, int, int]]:
    target = Fraction(1) - Fraction(1, 6 * index)
    candidates: list[tuple[int, int, int]] = []

    # From 1/a+1/b+1/c=target and a<=b<=c, one has
    # a<=3/target.  After a is fixed, b<=2/(target-1/a).
    a_max = (Fraction(3, 1) / target).numerator // (
        Fraction(3, 1) / target
    ).denominator
    for a in range(2, a_max + 1):
        remainder_after_a = target - Fraction(1, a)
        if remainder_after_a <= 0:
            continue
        b_bound = Fraction(2, 1) / remainder_after_a
        b_max = b_bound.numerator // b_bound.denominator
        for b in range(a, b_max + 1):
            remainder = target - Fraction(1, a) - Fraction(1, b)
            if remainder <= 0 or remainder.numerator != 1:
                continue
            c = remainder.denominator
            if c >= b:
                candidates.append((a, b, c))
    return candidates


def main() -> None:
    x, y = sp.symbols("x y", positive=True)
    gram = sp.Matrix(
        [
            [1, 0, -x, -sp.Rational(1, 2)],
            [0, 1, 0, -y],
            [-x, 0, 1, 0],
            [-sp.Rational(1, 2), -y, 0, 1],
        ]
    )
    determinant = sp.factor(gram.det())
    expected_determinant = (
        4 * x**2 * y**2 - 4 * x**2 - 4 * y**2 + 3
    ) / 4
    assert sp.simplify(determinant - expected_determinant) == 0

    principal_minors = [
        sp.factor(gram.minor_submatrix(i, i).det()) for i in range(4)
    ]
    expected_minors = [
        1 - y**2,
        sp.Rational(3, 4) - x**2,
        sp.Rational(3, 4) - y**2,
        1 - x**2,
    ]
    assert all(
        sp.simplify(actual - expected) == 0
        for actual, expected in zip(principal_minors, expected_minors)
    )

    y_squared = (4 * x**2 - 3) / (4 * (x**2 - 1))
    assert sp.simplify(expected_determinant.subs(y**2, y_squared)) == 0

    z = sp.symbols("z")
    transform = lambda value: 1 + 1 / (4 * (value - 1))
    assert sp.simplify(transform(transform(z)) - z) == 0

    triangles = {str(d): triangle_candidates(d) for d in range(2, 8)}
    assert triangles == {
        "2": [(2, 3, 12), (2, 4, 6), (3, 3, 4)],
        "3": [(2, 3, 9)],
        "4": [(2, 3, 8)],
        "5": [],
        "6": [],
        "7": [(2, 3, 7)],
    }

    # Fujita's strict SNC inequalities for n=2 and four weights 1/2.
    total_weight = Fraction(2)
    fujita_margins = {
        str(k): k * total_weight - 3 * k * Fraction(1, 2)
        for k in (1, 2)
    }
    assert all(margin > 0 for margin in fujita_margins.values())

    # The Q-arithmetic integral locus of (A-4)(B-4)=4.
    arithmetic_pairs = [
        (4 + divisor, 4 + 4 // divisor) for divisor in (1, 2, 4)
    ]
    assert arithmetic_pairs == [(5, 8), (6, 6), (8, 5)]

    certificate = {
        "determinant": str(determinant),
        "principal_minors": [str(value) for value in principal_minors],
        "area": "pi/6",
        "proper_overgroup_index_interval": [2, 7],
        "triangle_candidates": triangles,
        "stable_pair_log_class": "H",
        "stable_pair_log_volume": 1,
        "fujita_strict_margins": {
            key: str(value) for key, value in fujita_margins.items()
        },
        "rational_arithmetic_cycle_pairs": arithmetic_pairs,
    }
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
