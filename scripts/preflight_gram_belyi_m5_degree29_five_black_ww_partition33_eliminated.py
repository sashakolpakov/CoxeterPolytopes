#!/usr/bin/env python3
"""Quadratic-factor elimination preflight for partition (3^2,1^23)."""

from __future__ import annotations

import argparse
import json
import subprocess
import time

from certify_gram_belyi_m5_degree29_reduction import lambda_mod
from certify_gram_belyi_m5_degree29_five_black_ww_partition5 import (
    residual_derivative,
)
from certify_gram_belyi_m5_degree29_four_black_ww_simple import (
    outer_mul,
    outer_power,
    trim_outer,
)


def eliminated_equations(prime: int, lam: int, weights):
    context, derivative = residual_derivative(prime, lam, weights)
    u, v, w, z = context.gens()
    inv2 = pow(2, -1, prime)
    inv29 = pow(29, -1, prime)

    # If D=29*(x^2-S*x+P)^2, its x^3 and x^2 coefficients determine S,P.
    S = -pow(58, -1, prime) * derivative[3]
    P = inv2 * (inv29 * derivative[2] - S * S)
    square_equations = [
        derivative[1] + 58 * S * P,
        derivative[0] - 29 * P * P,
    ]

    roots = (context.constant(lam), u, v, w, z)
    product = [context.constant(1)]
    for root, exponent in zip(roots, weights):
        product = outer_mul(product, outer_power([-root, context.constant(1)], exponent))

    # Divide G by the monic quadratic C=x^2-S*x+P.
    remainder = list(product)
    quotient = [context.constant(0) for _ in range(len(product) - 2)]
    divisor = [P, -S, context.constant(1)]
    while len(remainder) >= 3:
        quotient_coefficient = remainder[-1]
        shift = len(remainder) - 3
        quotient[shift] = quotient_coefficient
        for index, coefficient in enumerate(divisor):
            remainder[shift + index] -= quotient_coefficient * coefficient
        trim_outer(remainder)
    remainder += [context.constant(0)] * (2 - len(remainder))
    B, A = remainder
    Q0 = quotient[0]
    Q1 = sum(quotient, context.constant(0))
    return square_equations + [A, Q0, Q1], S, P


def singular_unit(prime: int, equations, algorithm: str, incremental: bool):
    commands = [
            f"ring R={prime},(u,v,w,z),dp;",
            "option(redSB);",
            *(f"poly f{i}={value};" for i, value in enumerate(equations)),
    ]
    if incremental:
        commands += [
            "ideal I=f0,f1,f2,f3;",
            f"ideal B={algorithm}(I);",
            'print("BASE"); print(dim(B)); print(vdim(B));',
            f"ideal J={algorithm}(B+ideal(f4));",
        ]
    else:
        commands += [
            "ideal I=" + ",".join(f"f{i}" for i in range(len(equations))) + ";",
            f"ideal J={algorithm}(I);",
        ]
    commands += [
        'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(dim(J)); print(vdim(J)); }',
        "quit;",
    ]
    program = "\n".join(commands)
    before = time.monotonic()
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    return completed.stdout.strip(), completed.stderr.strip(), time.monotonic() - before


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("weights", nargs=5, type=int)
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument("--base", action="store_true", help="omit the final Q(1) test")
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="std")
    parser.add_argument("--incremental", action="store_true")
    arguments = parser.parse_args()
    weights = tuple(arguments.weights)
    if min(weights) <= 0 or sum(weights) != 29:
        raise SystemExit("positive weights must sum to 29")
    sqrt5, lam = lambda_mod(arguments.prime)
    before = time.monotonic()
    equations, S, P = eliminated_equations(arguments.prime, lam, weights)
    if arguments.base:
        equations = equations[:-1]
    stdout, stderr, singular_seconds = singular_unit(
        arguments.prime, equations, arguments.algorithm, arguments.incremental
    )
    print(
        json.dumps(
            {
                "weights": weights,
                "prime": arguments.prime,
                "sqrt5": sqrt5,
                "lambda": lam,
                "equation_total_degrees": [int(value.total_degree()) for value in equations],
                "S_total_degree": int(S.total_degree()),
                "P_total_degree": int(P.total_degree()),
                "result": stdout,
                "stderr": stderr,
                "singular_seconds": singular_seconds,
                "total_seconds": time.monotonic() - before,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
