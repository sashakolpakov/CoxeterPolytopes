#!/usr/bin/env python3
"""Lightweight exact preflight for five-black BW/WB critical endpoints."""

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


def critical_system(
    prime: int, alpha: int, weights: tuple[int, int, int, int, int]
):
    context = nmod_mpoly_ctx.get(("u", "v", "w"), prime, "degrevlex")
    u, v, w = context.gens()
    roots = (context.constant(alpha), context.constant(0), u, v, w)
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
    endpoint_critical = sum(derivative_residual, context.constant(0))
    return context, remainder[1:] + [endpoint_critical]


def singular_unit(prime: int, equations, algorithm: str, ordering: str = "dp") -> tuple[str, float]:
    variables = "(u,v,w)" if ordering == "dp" else "(w,u,v)"
    program = "\n".join(
        [
            f"ring r={prime},{variables},{ordering};",
            "option(redSB);",
            *(f"poly r{i}={value};" for i, value in enumerate(equations, 1)),
            "ideal I=" + ",".join(f"r{i}" for i in range(1, len(equations) + 1)) + ";",
            f"ideal J={algorithm}(I);",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(J); }',
            "quit;",
        ]
    )
    before = time.monotonic()
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    if completed.stderr.strip():
        raise ArithmeticError(completed.stderr)
    return completed.stdout.strip(), time.monotonic() - before


def convert_w_linear(value, target_context):
    u, v = target_context.gens()
    constant = target_context.constant(0)
    linear = target_context.constant(0)
    for (u_power, v_power, w_power), coefficient in value.terms():
        term = target_context.constant(int(coefficient)) * (u**u_power) * (v**v_power)
        if w_power == 0:
            constant += term
        elif w_power == 1:
            linear += term
        else:
            raise ArithmeticError(("not linear in w", value.degrees()))
    return linear, constant


def clear_linear_substitution(equations, endpoint_equation, prime: int):
    """Return C^26 R_i(u,v,-K/C), plus the exceptional C,K."""

    target = nmod_mpoly_ctx.get(("u", "v"), prime, "degrevlex")
    u, v = target.gens()
    c_value, k_value = convert_w_linear(endpoint_equation, target)
    c_powers = [target.constant(1)]
    k_powers = [target.constant(1)]
    for _ in range(26):
        c_powers.append(c_powers[-1] * c_value)
        k_powers.append(k_powers[-1] * (-k_value))
    cleared = []
    for equation in equations:
        if int(equation.degrees()[2]) > 26:
            raise ArithmeticError(("w degree", equation.degrees()))
        answer = target.constant(0)
        for (u_power, v_power, w_power), coefficient in equation.terms():
            answer += (
                target.constant(int(coefficient))
                * (u**u_power)
                * (v**v_power)
                * k_powers[w_power]
                * c_powers[26 - w_power]
            )
        cleared.append(answer)
    return cleared, c_value, k_value


def singular_unit_two_variables(prime: int, equations, algorithm: str) -> tuple[str, float]:
    program = "\n".join(
        [
            f"ring r={prime},(u,v),dp;",
            "option(redSB);",
            *(f"poly r{i}={value};" for i, value in enumerate(equations, 1)),
            "ideal I=" + ",".join(f"r{i}" for i in range(1, len(equations) + 1)) + ";",
            f"ideal J={algorithm}(I);",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); }',
            "quit;",
        ]
    )
    before = time.monotonic()
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    if completed.stderr.strip():
        raise ArithmeticError(completed.stderr)
    return completed.stdout.strip(), time.monotonic() - before


def singular_saturated_unit(prime: int, equations, denominator) -> tuple[str, float]:
    program = "\n".join(
        [
            'LIB "elim.lib";',
            f"ring r={prime},(u,v),dp;",
            "option(redSB);",
            *(f"poly r{i}={value};" for i, value in enumerate(equations, 1)),
            f"poly c={denominator};",
            "ideal I=" + ",".join(f"r{i}" for i in range(1, len(equations) + 1)) + ";",
            "ideal J=sat(I,ideal(c))[1];",
            "J=std(J);",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); }',
            "quit;",
        ]
    )
    before = time.monotonic()
    completed = subprocess.run(
        ["Singular", "-q"], input=program, text=True, capture_output=True, check=True
    )
    if completed.stderr.strip():
        raise ArithmeticError(completed.stderr)
    return completed.stdout.strip(), time.monotonic() - before


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("weights", nargs=5, type=int)
    parser.add_argument("--prime", type=int, default=101)
    parser.add_argument("--orientation", choices=("BW", "WB"), default="BW")
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="slimgb")
    parser.add_argument("--substitute", action="store_true")
    parser.add_argument("--saturate", action="store_true")
    parser.add_argument("--ordering", choices=("dp", "lp"), default="dp")
    arguments = parser.parse_args()
    weights = tuple(arguments.weights)
    if min(weights) <= 0 or sum(weights) != 29:
        raise SystemExit("positive weights must sum to 29")
    _sqrt5, lam = lambda_mod(arguments.prime)
    alpha = lam if arguments.orientation == "BW" else (1 - lam) % arguments.prime
    started = time.monotonic()
    _context, equations = critical_system(arguments.prime, alpha, weights)
    if arguments.substitute:
        cleared, c_value, k_value = clear_linear_substitution(
            equations[:3], equations[3], arguments.prime
        )
        if arguments.saturate:
            result, elapsed = singular_saturated_unit(
                arguments.prime, cleared, c_value
            )
        else:
            result, elapsed = singular_unit_two_variables(
                arguments.prime, cleared, arguments.algorithm
            )
        exceptional, exceptional_elapsed = singular_unit(
            arguments.prime,
            equations[:3] + [
                # Lift C,K back to the three-variable context as polynomials
                # independent of w; their string syntax is unchanged.
                c_value,
                k_value,
            ],
            arguments.algorithm,
            arguments.ordering,
        )
        substituted = {
            "cleared_total_degrees": [int(value.total_degree()) for value in cleared],
            "cleared_multidegrees": [
                [int(degree) for degree in value.degrees()] for value in cleared
            ],
            "generic_standard_basis": result,
            "generic_seconds": elapsed,
            "exceptional_standard_basis": exceptional,
            "exceptional_seconds": exceptional_elapsed,
            "C_degree": int(c_value.total_degree()),
            "K_degree": int(k_value.total_degree()),
        }
    else:
        result, elapsed = singular_unit(
            arguments.prime, equations, arguments.algorithm, arguments.ordering
        )
        substituted = None
    print(
        json.dumps(
            {
                "orientation": arguments.orientation,
                "weights": weights,
                "prime": arguments.prime,
                "total_degrees": [int(value.total_degree()) for value in equations],
                "multidegrees": [
                    [int(degree) for degree in value.degrees()] for value in equations
                ],
                "affine_standard_basis": result,
                "algorithm": arguments.algorithm,
                "singular_seconds": elapsed,
                "substitution": substituted,
                "total_seconds": time.monotonic() - started,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
