#!/usr/bin/env python3
"""Independent superset replay of the global zero-E8 no-go.

This file imports no project code and enumerates the covectors directly in
``(A,C)`` coordinates.  The primary proof gives ``|X_a|,|X_b|<=1`` for
``X=-A/phi``, hence conservatively ``|A_a|,|A_b|<=2``.  Together with its
``|C_c|<=3, |C_d|<=1`` bound, this direct 525-pair superset must contain every
separator.  It contains no separator.
"""

from __future__ import annotations

from collections import Counter


Pair = tuple[int, int]
INV_PHI: Pair = (-1, 1)
INTERIOR_XU: tuple[Pair, Pair] = ((601, 199), (2, 0))
RAY_XU: tuple[Pair, Pair] = ((13, 0), (0, 1))


def add(x: Pair, y: Pair) -> Pair:
    return x[0] + y[0], x[1] + y[1]


def neg(x: Pair) -> Pair:
    return -x[0], -x[1]


def mul(x: Pair, y: Pair) -> Pair:
    a, b = x
    c, d = y
    return a * c + b * d, a * d + b * c + b * d


def conjugate(x: Pair) -> Pair:
    return x[0] + x[1], -x[1]


def sign(x: Pair) -> int:
    a, b = x
    t = 2 * a + b
    if t == 0:
        return (b > 0) - (b < 0)
    if b == 0 or (t > 0) == (b > 0):
        return (t > 0) - (t < 0)
    if t * t > 5 * b * b:
        return (t > 0) - (t < 0)
    return (b > 0) - (b < 0)


def run() -> dict[str, object]:
    admissible_count = 0
    sign_histogram: Counter[str] = Counter()
    cutting_records = []
    for aa in range(-2, 3):
        for ab in range(-2, 3):
            A = (aa, ab)
            for ca in range(-3, 4):
                for cb in range(-1, 2):
                    C = (ca, cb)
                    q = add(neg(mul(mul(A, A), INV_PHI)), mul(C, C))
                    remainder = add((4, 0), neg(q))
                    if not (
                        sign(q) > 0
                        and sign(remainder) > 0
                        and sign(conjugate(q)) > 0
                        and sign(conjugate(remainder)) > 0
                    ):
                        continue
                    admissible_count += 1
                    interior = add(mul(A, INTERIOR_XU[0]), mul(C, INTERIOR_XU[1]))
                    ray = add(mul(A, RAY_XU[0]), mul(C, RAY_XU[1]))
                    signs = sign(interior), sign(ray)
                    sign_histogram[f"{signs[0]},{signs[1]}"] += 1
                    if signs[0] * signs[1] < 0:
                        cutting_records.append((A, C, q, interior, ray))
    return {
        "literal_superset_pair_count": 5 * 5 * 7 * 3,
        "signature_admissible_count": admissible_count,
        "sign_histogram": dict(sorted(sign_histogram.items())),
        "cutting_records": cutting_records,
    }


if __name__ == "__main__":
    print(run())
