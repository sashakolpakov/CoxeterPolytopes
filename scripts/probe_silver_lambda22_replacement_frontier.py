#!/usr/bin/env python3
"""Probe the one-direction Lambda22 replacement ray past H1598.

This is an exact necessary-superset scan, not a Vinberg termination claim.
It reports every control whose relaxed, gate-compatible Gram data can pair
nonnegatively with ``(1,(10/31)d)``.  Arithmetic resolutions are deliberately
left to a subsequent certificate once the first new threat is isolated.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import lcm

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    LAMBDA22_COEFFICIENT,
    enumerate_lambda22_controls,
    gated_relaxed_witness_bound,
)
from coxeter_polytopes.search.silver_leech import candidate_leech_gram_triples
from coxeter_polytopes.search.silver_leech_manifest import (
    has_short_binary_combination,
)


START = Silver(1598, 1130)
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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--real-height",
        type=int,
        default=5000,
        help="identity-embedding height bound (a rational field element)",
    )
    parser.add_argument("--height-a", type=int)
    parser.add_argument("--height-b", type=int)
    parser.add_argument("--scale", default="10/31")
    arguments = parser.parse_args()
    if (arguments.height_a is None) != (arguments.height_b is None):
        parser.error("--height-a and --height-b must be supplied together")
    bound = (
        Silver(arguments.height_a, arguments.height_b)
        if arguments.height_a is not None
        else Silver(arguments.real_height)
    )
    scale = Silver(Fraction(arguments.scale))
    retained = 0
    controls = 0
    rows: list[tuple[object, ...]] = []
    for level in enumerate_lambda22_controls(bound):
        if level.height <= START:
            continue
        controls += 1
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
            try:
                d_u, d_v, maximum = gated_relaxed_witness_bound(
                    *gram, u_modulus, v_modulus
                )
            except AssertionError:
                continue
            retained += 1
            margin = scale * maximum - LAMBDA22_COEFFICIENT * level.control
            if margin.sign(0) >= 0:
                rows.append(
                    (
                        level.height,
                        level.norm,
                        level.control,
                        gram,
                        (d_u, d_v),
                        maximum,
                        margin,
                    )
                )
                print(
                    "THREAT",
                    *(str(value) for value in rows[-1][:3]),
                    gram,
                    (d_u, d_v),
                    str(maximum),
                    str(margin),
                    flush=True,
                )
    print(
        "SUMMARY",
        f"bound={bound}",
        f"scale={scale}",
        f"post_start_controls={controls}",
        f"retained_types={retained}",
        f"threats={len(rows)}",
    )


if __name__ == "__main__":
    main()
