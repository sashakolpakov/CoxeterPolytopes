#!/usr/bin/env python3
"""Exact chronological Schur probe for the candidate H2786 Lambda22 wall.

The fixed vectors are ``d=e1`` and ``e=-g738`` with Gram
``[[4,2],[2,10]]``.  This script deliberately prints an exact necessary
superset; positive extensions still require arithmetic resolution before the
candidate can be called accepted.
"""

from __future__ import annotations

from functools import lru_cache
from math import isqrt, lcm, sqrt

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    LAMBDA22_COEFFICIENT,
    enumerate_lambda22_controls,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import (
    has_short_binary_combination,
)


HEIGHT = Silver(2786, 1970)
CONTROL = Silver(82, 58)
H14_CONTROL = Silver(6, 4)
NORM_GATES = {
    Silver(1): (1, 1),
    Silver(2): (1, 1),
    Silver(4, 2): (2, 1),
    Silver(4): (2, 2),
    Silver(3): (3, 3),
    Silver(6): (3, 3),
    Silver(12, 6): (6, 3),
    Silver(12): (6, 6),
}
MINIMUM = {1: 4, 2: 12, 3: 24, 6: 96}


def norm_gate(square: int, modulus: int) -> bool:
    return (
        (not square or square >= MINIMUM[modulus])
        and square % (2 * modulus) == 0
        and (modulus % 2 or square % 8 in (0, 4))
    )


@lru_cache(maxsize=None)
def one_rows(
    square: int, modulus: int, kind: int
) -> tuple[tuple[float, Silver, int, int, int], ...]:
    rows = []
    for d_dot in range(-isqrt(4 * square), isqrt(4 * square) + 1):
        if d_dot % modulus:
            continue
        for e_dot in range(-isqrt(10 * square), isqrt(10 * square) + 1):
            if e_dot % modulus:
                continue
            schur = (
                36 * square
                - 10 * d_dot * d_dot
                + 4 * d_dot * e_dot
                - 4 * e_dot * e_dot
            )
            if schur < 0:
                continue
            value = (
                Silver(26 * d_dot + e_dot, 18 * d_dot + e_dot)
                if kind == 0
                else Silver(36 * d_dot + 2 * e_dot, 26 * d_dot + e_dot)
            )
            numeric = float(value.a) + sqrt(2) * float(value.b)
            rows.append((numeric, value, d_dot, e_dot, schur))
    return tuple(sorted(rows, reverse=True))


def cross_extensions(
    u_square: int,
    v_square: int,
    u_dot_v: int,
    u_modulus: int,
    v_modulus: int,
    other_control: Silver,
) -> tuple[tuple[Silver, int, int, int, int, int, int, int, Silver], ...]:
    threshold = LAMBDA22_COEFFICIENT * CONTROL * other_control
    threshold_numeric = float(threshold.a) + sqrt(2) * float(threshold.b)
    u_rows = one_rows(u_square, u_modulus, 0)
    v_rows = one_rows(v_square, v_modulus, 1)
    output = []
    if not u_rows or not v_rows:
        return ()
    for u_numeric, u_value, d_u, e_u, schur_a in u_rows:
        if u_numeric + v_rows[0][0] < threshold_numeric - 1e-8:
            break
        for v_numeric, v_value, d_v, e_v, schur_b in v_rows:
            if u_numeric + v_numeric < threshold_numeric - 1e-8:
                break
            pairing = u_value + v_value - threshold
            if pairing.sign(0) < 0:
                break
            schur_n = 36 * u_dot_v - (
                10 * d_u * d_v
                - 2 * (d_u * e_v + e_u * d_v)
                + 4 * e_u * e_v
            )
            if schur_a * schur_b < schur_n * schur_n:
                continue
            h14_e_pairing = (
                Silver(e_u + 2 * e_v, e_u + e_v)
                - LAMBDA22_COEFFICIENT * H14_CONTROL * other_control
            )
            if h14_e_pairing.sign(0) > 0:
                continue
            output.append(
                (
                    pairing,
                    d_u,
                    d_v,
                    e_u,
                    e_v,
                    schur_a,
                    schur_b,
                    schur_n,
                    h14_e_pairing,
                )
            )
    return tuple(output)


def main() -> None:
    type_count = 0
    row_count = 0
    extension_count = 0
    positive_count = 0
    for level in enumerate_lambda22_controls(HEIGHT):
        if level.height <= Silver(14, 10):
            continue
        u_modulus, v_modulus = NORM_GATES[level.norm]
        for item in candidate_leech_gram_triples(level.target):
            if has_short_binary_combination(item):
                continue
            gram = (item.u_square, item.v_square, item.u_dot_v)
            if (
                not norm_gate(gram[0], u_modulus)
                or not norm_gate(gram[1], v_modulus)
                or gram[2] % lcm(u_modulus, v_modulus)
            ):
                continue
            type_count += 1
            extensions = cross_extensions(
                *gram, u_modulus, v_modulus, level.control
            )
            if not extensions:
                continue
            row_count += 1
            extension_count += len(extensions)
            positives = [row for row in extensions if row[0].sign(0) > 0]
            positive_count += len(positives)
            print(
                "ROW",
                str(level.height),
                str(level.norm),
                str(level.control),
                gram,
                "extensions",
                len(extensions),
                "positives",
                len(positives),
                flush=True,
            )
            for extension in positives:
                print(
                    "POSITIVE",
                    str(level.height),
                    gram,
                    tuple(str(value) if isinstance(value, Silver) else value for value in extension),
                    flush=True,
                )
    print(
        "SUMMARY",
        f"types={type_count}",
        f"rows={row_count}",
        f"extensions={extension_count}",
        f"positives={positive_count}",
    )


if __name__ == "__main__":
    main()
