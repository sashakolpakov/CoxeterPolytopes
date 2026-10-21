#!/usr/bin/env python3
"""Lightweight exact preflight for five-black BB critical equations.

This is not a certificate producer.  It measures the bidegrees and pairwise
resultants for selected records before a bounded census is designed.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time

from flint import nmod_mpoly_ctx

from certify_gram_belyi_m5_degree29_reduction import lambda_mod
from certify_gram_belyi_m5_degree29_four_black_ww_simple import (
    mpoly_to_nmod_poly,
    outer_add,
    outer_mul,
    outer_power,
    trim_outer,
)


def bb_remainder(
    prime: int, lam: int, weights: tuple[int, int, int, int, int]
):
    """Return R0,...,R3 for roots lambda,0,1,u,v over F_p."""

    context = nmod_mpoly_ctx.get(("u", "v"), prime, "degrevlex")
    u, v = context.gens()
    roots = (context.constant(lam), context.constant(0), context.constant(1), u, v)
    factors = [[-root, context.constant(1)] for root in roots]
    product_value = [context.constant(1)]
    for factor, exponent in zip(factors, weights):
        product_value = outer_mul(product_value, outer_power(factor, exponent))

    derivative_residual = [context.constant(0)]
    for omitted, weight in enumerate(weights):
        term = [context.constant(1)]
        for index, factor in enumerate(factors):
            if index != omitted:
                term = outer_mul(term, factor)
        derivative_residual = outer_add(
            derivative_residual, [weight * coefficient for coefficient in term]
        )

    inverse29 = pow(29, -1, prime)
    remainder = list(product_value)
    while len(remainder) >= 5:
        quotient_coefficient = inverse29 * remainder[-1]
        shift = len(remainder) - 5
        for index, coefficient in enumerate(derivative_residual):
            remainder[shift + index] -= quotient_coefficient * coefficient
        trim_outer(remainder)
    remainder += [context.constant(0)] * (4 - len(remainder))
    return context, remainder


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("weights", nargs=5, type=int)
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument("--singular", action="store_true")
    arguments = parser.parse_args()
    weights = tuple(arguments.weights)
    if sum(weights) != 29 or min(weights) <= 0:
        raise SystemExit("weights must be positive and sum to 29")
    _sqrt5, lam = lambda_mod(arguments.prime)
    started = time.monotonic()
    _context, remainder = bb_remainder(arguments.prime, lam, weights)
    equations = remainder[1:]
    rows = []
    resultant_polynomials = []
    for left, right, label in (
        (equations[0], equations[1], "R1_R2"),
        (equations[0], equations[2], "R1_R3"),
        (equations[1], equations[2], "R2_R3"),
    ):
        before = time.monotonic()
        resultant_mpoly = left.resultant(right, "v")
        resultant = mpoly_to_nmod_poly(resultant_mpoly, arguments.prime)
        resultant_polynomials.append(resultant)
        rows.append(
            {
                "pair": label,
                "degree_u": resultant.degree(),
                "elapsed_seconds": time.monotonic() - before,
            }
        )
    gcd12_13 = resultant_polynomials[0].gcd(resultant_polynomials[1])
    gcd_all = gcd12_13.gcd(resultant_polynomials[2])
    singular = None
    if arguments.singular:
        singular_input = "\n".join(
            [
                f"ring r={arguments.prime},(u,v),dp;",
                *(f"poly r{index}={value};" for index, value in enumerate(equations, 1)),
                "ideal I=r1,r2,r3;",
                "ideal J=std(I);",
                'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print(J); }',
                "quit;",
            ]
        )
        before = time.monotonic()
        completed = subprocess.run(
            ["Singular", "-q"],
            input=singular_input,
            text=True,
            capture_output=True,
            check=True,
        )
        singular = {
            "stdout": completed.stdout.strip(),
            "stderr": completed.stderr.strip(),
            "elapsed_seconds": time.monotonic() - before,
        }
    print(
        json.dumps(
            {
                "prime": arguments.prime,
                "weights": weights,
                "equation_total_degrees": [int(value.total_degree()) for value in equations],
                "equation_bidegrees": [
                    [int(degree) for degree in value.degrees()]
                    for value in equations
                ],
                "resultants": rows,
                "gcd_R12_R13_degree": gcd12_13.degree(),
                "gcd_all_degree": gcd_all.degree(),
                "singular": singular,
                "elapsed_seconds": time.monotonic() - started,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
