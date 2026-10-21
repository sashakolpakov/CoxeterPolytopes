#!/usr/bin/env python3
"""Bounded exact preflight for the five-black WW master equations.

This script is exploratory only.  It forms the residual-quartic common-value
equations and the two fixed endpoint-ratio equations for selected black-degree
records, then asks Singular whether the affine ideal is the unit ideal.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time

from flint import nmod_mpoly_ctx

from certify_gram_belyi_m5_degree29_reduction import lambda_mod
from certify_gram_belyi_m5_degree29_four_black_ww_simple import (
    outer_add,
    outer_mul,
    outer_power,
    trim_outer,
)


def canonical_remainder(
    prime: int, weights: tuple[int, int, int, int, int]
):
    """Return R0,...,R3 for roots 0,1,u,v,w over F_p."""

    context = nmod_mpoly_ctx.get(("u", "v", "w"), prime, "degrevlex")
    u, v, w = context.gens()
    roots = (context.constant(0), context.constant(1), u, v, w)
    factors = [[-root, context.constant(1)] for root in roots]
    product_value = [context.constant(1)]
    for factor, exponent in zip(factors, weights):
        product_value = outer_mul(
            product_value, outer_power(factor, exponent)
        )
    derivative_residual = [context.constant(0)]
    for omitted, weight in enumerate(weights):
        term = [context.constant(1)]
        for index, factor in enumerate(factors):
            if index != omitted:
                term = outer_mul(term, factor)
        derivative_residual = outer_add(
            derivative_residual, [weight * coefficient for coefficient in term]
        )
    remainder = list(product_value)
    inverse29 = pow(29, -1, prime)
    while len(remainder) >= 5:
        quotient_coefficient = inverse29 * remainder[-1]
        shift = len(remainder) - 5
        for index, coefficient in enumerate(derivative_residual):
            remainder[shift + index] -= quotient_coefficient * coefficient
        trim_outer(remainder)
    remainder += [context.constant(0)] * (4 - len(remainder))
    return context, remainder, derivative_residual


def singular_affine_unit(
    prime: int,
    q: int,
    weights: tuple[int, int, int, int, int],
    remainder,
    derivative_residual,
    algorithm: str,
    status: str,
    variable_order: str,
) -> tuple[str, float, str]:
    r0, r1, r2, r3 = remainder
    a, b, c, d, e = weights
    g_y = f"y^{a}*(y-1)^{b}*(y-u)^{c}*(y-v)^{d}*(y-w)^{e}"
    g_qy = (
        f"({q}*y)^{a}*({q}*y-1)^{b}*({q}*y-u)^{c}*"
        f"({q}*y-v)^{d}*({q}*y-w)^{e}"
    )
    d_y = "+".join(
        f"({coefficient})*y^{power}"
        for power, coefficient in enumerate(derivative_residual)
    )
    d_qy = "+".join(
        f"({coefficient})*({q}*y)^{power}"
        for power, coefficient in enumerate(derivative_residual)
    )
    endpoint_equations = {
        "SS": [f"{g_y}-r0", f"{g_qy}-r0"],
        "SC": [f"{g_y}-r0", d_qy],
        "CS": [d_y, f"{g_qy}-r0"],
        "CC": [d_y, d_qy],
    }[status]
    program = "\n".join(
        [
            'LIB "standard.lib";',
            f"ring r={prime},({variable_order}),dp;",
            "option(redSB);",
            f"poly r0={r0};",
            f"poly r1={r1};",
            f"poly r2={r2};",
            f"poly r3={r3};",
            f"poly h1={endpoint_equations[0]};",
            f"poly h2={endpoint_equations[1]};",
            "ideal I=r1,r2,r3,h1,h2;",
            f"ideal J={algorithm}(I);",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(dim(J)); print(vdim(J)); }',
            "quit;",
        ]
    )
    before = time.monotonic()
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    return completed.stdout.strip(), time.monotonic() - before, completed.stderr.strip()


def singular_critical_summary(prime: int, remainder, algorithm: str):
    _r0, r1, r2, r3 = remainder
    program = "\n".join(
        [
            'LIB "standard.lib";',
            f"ring r={prime},(u,v,w),dp;",
            "option(redSB);",
            f"poly r1={r1};",
            f"poly r2={r2};",
            f"poly r3={r3};",
            "ideal I=r1,r2,r3;",
            f"ideal J={algorithm}(I);",
            'print("SIZE",size(J),"DIM",dim(J),"VDIM",vdim(J));',
            "quit;",
        ]
    )
    before = time.monotonic()
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    return completed.stdout.strip(), time.monotonic() - before, completed.stderr.strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("weights", nargs=5, type=int)
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument(
        "--algorithm",
        choices=("std", "slimgb", "groebner", "stdhilb"),
        default="std",
    )
    parser.add_argument("--critical-only", action="store_true")
    parser.add_argument("--status", choices=("SS", "SC", "CS", "CC"), default="SS")
    parser.add_argument(
        "--variable-order",
        choices=("u,v,w,y", "y,w,v,u", "y,u,v,w", "w,y,u,v"),
        default="u,v,w,y",
    )
    arguments = parser.parse_args()
    weights = tuple(arguments.weights)
    if min(weights) <= 0 or sum(weights) != 29:
        raise SystemExit("positive weights must sum to 29")
    _sqrt5, lam = lambda_mod(arguments.prime)
    q = (29 * lam - 13) % arguments.prime
    started = time.monotonic()
    _context, remainder, derivative_residual = canonical_remainder(
        arguments.prime, weights
    )
    if arguments.critical_only:
        result, singular_seconds, stderr = singular_critical_summary(
            arguments.prime, remainder, arguments.algorithm
        )
    else:
        result, singular_seconds, stderr = singular_affine_unit(
            arguments.prime,
            q,
            weights,
            remainder,
            derivative_residual,
            arguments.algorithm,
            arguments.status,
            arguments.variable_order,
        )
    print(
        json.dumps(
            {
                "weights": weights,
                "prime": arguments.prime,
                "q": q,
                "q_polynomial": (q * q + 12 * q + 16) % arguments.prime,
                "critical_total_degrees": [
                    int(value.total_degree()) for value in remainder[1:]
                ],
                "critical_multidegrees": [
                    [int(degree) for degree in value.degrees()]
                    for value in remainder[1:]
                ],
                "algorithm": arguments.algorithm,
                "status": arguments.status,
                "variable_order": arguments.variable_order,
                "affine_result": result,
                "singular_stderr": stderr,
                "singular_seconds": singular_seconds,
                "total_seconds": time.monotonic() - started,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
