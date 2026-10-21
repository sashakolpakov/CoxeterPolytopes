#!/usr/bin/env python3
"""Exact preflight for five-black WW in the original endpoint coordinates."""

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


def original_remainder(
    prime: int, lam: int, weights: tuple[int, int, int, int, int]
):
    context = nmod_mpoly_ctx.get(("u", "v", "w", "z"), prime, "degrevlex")
    u, v, w, z = context.gens()
    roots = (context.constant(lam), u, v, w, z)
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


def singular_unit(prime: int, lam: int, weights, remainder, algorithm: str):
    r0, r1, r2, r3 = remainder
    a, b, c, d, e = weights
    g0 = f"(-L)^{a}*(-u)^{b}*(-v)^{c}*(-w)^{d}*(-z)^{e}"
    g1 = (
        f"(1-L)^{a}*(1-u)^{b}*(1-v)^{c}*(1-w)^{d}*(1-z)^{e}"
    )
    program = "\n".join(
        [
            f"ring r={prime},(u,v,w,z),dp;",
            "option(redSB);",
            f"number L={lam};",
            f"poly r0={r0};",
            f"poly r1={r1};",
            f"poly r2={r2};",
            f"poly r3={r3};",
            f"poly h0={g0}-r0;",
            f"poly h1={g1}-r0;",
            "ideal I=r1,r2,r3,h0,h1;",
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


def singular_partition_five_unit(
    prime: int,
    lam: int,
    weights,
    derivative_residual,
    algorithm: str,
    projective: bool = False,
):
    a, b, c_weight, d, e = weights
    g0 = f"(-L)^{a}*(-u)^{b}*(-v)^{c_weight}*(-w)^{d}*(-z)^{e}"
    g1 = (
        f"(1-L)^{a}*(1-u)^{b}*(1-v)^{c_weight}*"
        f"(1-w)^{d}*(1-z)^{e}"
    )
    gc = (
        f"(c-L)^{a}*(c-u)^{b}*(c-v)^{c_weight}*"
        f"(c-w)^{d}*(c-z)^{e}"
    )
    desired = (
        "29*c^4",
        "-116*c^3",
        "174*c^2",
        "-116*c",
    )
    variables = "(u,v,w,z,c,s)" if projective else "(u,v,w,z,c)"
    commands = [
            f"ring r={prime},{variables},dp;",
            "option(redSB);",
            f"number L={lam};",
            *(
                f"poly d{index}=({derivative_residual[index]})-({desired[index]});"
                for index in range(4)
            ),
            f"poly h0={g0}-({gc});",
            f"poly h1={g1}-({gc});",
            "ideal I=d0,d1,d2,d3,h0,h1;",
            f"ideal J={algorithm}(I);",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(dim(J)); print(vdim(J)); }',
    ]
    if projective:
        commands += [
            "poly D0=homog(d0,s); poly D1=homog(d1,s);",
            "poly D2=homog(d2,s); poly D3=homog(d3,s);",
            "poly H0=homog(h0,s); poly H1=homog(h1,s);",
        ]
        for chart in ("u", "v", "w", "z", "c"):
            commands += [
                f"ideal K{subst_index(chart)}="
                + ",".join(
                    f"subst(subst({name},s,0),{chart},1)"
                    for name in ("D0", "D1", "D2", "D3", "H0", "H1")
                )
                + ";",
                f"ideal L{subst_index(chart)}={algorithm}(K{subst_index(chart)});",
                f'if (size(L{subst_index(chart)})==1 && string(L{subst_index(chart)}[1])=="1") {{ print("INF_{chart}_UNIT"); }} else {{ print("INF_{chart}_NONUNIT"); print(dim(L{subst_index(chart)})); }}',
            ]
    commands.append("quit;")
    program = "\n".join(commands)
    before = time.monotonic()
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    return completed.stdout.strip(), time.monotonic() - before, completed.stderr.strip()


def subst_index(variable: str) -> int:
    return {"u": 0, "v": 1, "w": 2, "z": 3, "c": 4}[variable]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("weights", nargs=5, type=int)
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="std")
    parser.add_argument("--partition-five", action="store_true")
    parser.add_argument("--projective", action="store_true")
    arguments = parser.parse_args()
    weights = tuple(arguments.weights)
    if min(weights) <= 0 or sum(weights) != 29:
        raise SystemExit("positive weights must sum to 29")
    _sqrt5, lam = lambda_mod(arguments.prime)
    started = time.monotonic()
    _context, remainder, derivative_residual = original_remainder(
        arguments.prime, lam, weights
    )
    if arguments.partition_five:
        result, elapsed, stderr = singular_partition_five_unit(
            arguments.prime,
            lam,
            weights,
            derivative_residual,
            arguments.algorithm,
            arguments.projective,
        )
    else:
        result, elapsed, stderr = singular_unit(
            arguments.prime, lam, weights, remainder, arguments.algorithm
        )
    print(
        json.dumps(
            {
                "weights": weights,
                "prime": arguments.prime,
                "critical_total_degrees": [
                    int(value.total_degree()) for value in remainder[1:]
                ],
                "critical_multidegrees": [
                    [int(degree) for degree in value.degrees()]
                    for value in remainder[1:]
                ],
                "result": result,
                "partition_five": arguments.partition_five,
                "projective": arguments.projective,
                "stderr": stderr,
                "singular_seconds": elapsed,
                "total_seconds": time.monotonic() - started,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
