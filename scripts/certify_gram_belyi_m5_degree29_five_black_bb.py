#!/usr/bin/env python3
"""Exact no-go certificate for the degree-29 five-black BB block.

With lambda, 0, and 1 all black, write

    G(x)=(x-lambda)^a*x^b*(x-1)^c*(x-u)^d*(x-v)^e.

After removing the forced black-root powers from G', the residual derivative
is a quartic D.  If G mod D = R3*x^3+R2*x^2+R1*x+R0, the Shabat condition
requires R1=R2=R3=0.  For each of the 11,011 canonical degree records this
script proves that the three bivariate equations have no common projective
point after reduction at a good split prime:

* Singular computes a standard basis [1] on the affine chart;
* the three leading homogeneous forms have gcd 1 on u != 0 at infinity,
  while at least one is nonzero at [u:v]=[0:1].

Projective properness then prevents a characteristic-zero solution from
escaping to infinity.  This is a bounded coefficient-system census, not a
plane-tree search and not a complete five-black no-go.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

from flint import nmod_mpoly_ctx, nmod_poly

from certify_gram_belyi_m5_degree29_reduction import lambda_mod
from certify_gram_belyi_m5_degree29_five_black_quotient import (
    canonical_black_records,
)
from certify_gram_belyi_m5_degree29_four_black_ww_simple import (
    outer_add,
    outer_mul,
    outer_power,
    trim_outer,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/gram_belyi_m5_degree29_five_black_bb.json"
SCHEMA = "coxeter.gram_belyi.m5_degree29_five_black_bb.v1"
SPLIT_PRIMES = (101, 109, 131, 139, 149, 151, 179, 181, 191, 199)
HOMOGENIZATION_DEGREES = (27, 27, 26)


def bb_remainder(
    prime: int, lam: int, weights: tuple[int, int, int, int, int]
):
    """Return R0,...,R3 for roots lambda,0,1,u,v over F_p."""

    context = nmod_mpoly_ctx.get(("u", "v"), prime, "degrevlex")
    u, v = context.gens()
    roots = (
        context.constant(lam),
        context.constant(0),
        context.constant(1),
        u,
        v,
    )
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
            derivative_residual,
            [weight * coefficient for coefficient in term],
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


def leading_form_chart_u_one(value, degree: int, prime: int) -> nmod_poly:
    """Dehomogenize the degree-``degree`` leading form by u=1."""

    coefficients = [0] * (degree + 1)
    for (u_power, v_power), coefficient in value.terms():
        if u_power + v_power == degree:
            coefficients[v_power] = int(coefficient)
    return nmod_poly(coefficients, prime)


def point_u_zero_v_one(value, degree: int) -> int:
    """Evaluate the fixed-degree leading form at [u:v]=[0:1]."""

    return int(value.to_dict().get((0, degree), 0))


def affine_standard_basis_is_unit(
    prime: int, equations
) -> tuple[bool, str, str]:
    """Use an exact Singular standard basis on F_p[u,v]."""

    program = "\n".join(
        [
            f"ring r={prime},(u,v),dp;",
            "option(redSB);",
            *(f"poly r{index}={value};" for index, value in enumerate(equations, 1)),
            "ideal I=r1,r2,r3;",
            "ideal J=std(I);",
            'if (size(J)==1 && string(J[1])=="1") { print("UNIT"); } else { print("NONUNIT"); print(J); }',
            "quit;",
        ]
    )
    completed = subprocess.run(
        ["Singular", "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    stdout = completed.stdout.strip()
    return stdout == "UNIT", stdout, completed.stderr.strip()


def certify_one(
    prime: int, weights: tuple[int, int, int, int, int]
) -> dict[str, Any] | None:
    sqrt5, lam = lambda_mod(prime)
    _context, remainder = bb_remainder(prime, lam, weights)
    equations = remainder[1:]
    total_degrees = tuple(int(value.total_degree()) for value in equations)
    bidegrees = tuple(
        tuple(int(degree) for degree in value.degrees()) for value in equations
    )
    if any(
        actual > bound
        for actual, bound in zip(total_degrees, HOMOGENIZATION_DEGREES)
    ):
        raise ArithmeticError(("structural degree bound failed", weights, total_degrees))

    leading = [
        leading_form_chart_u_one(value, degree, prime)
        for value, degree in zip(equations, HOMOGENIZATION_DEGREES)
    ]
    nonzero_leading = [value for value in leading if value]
    if not nonzero_leading:
        return None
    infinity_gcd = nonzero_leading[0]
    for value in nonzero_leading[1:]:
        infinity_gcd = infinity_gcd.gcd(value)
    point_values = [
        point_u_zero_v_one(value, degree) % prime
        for value, degree in zip(equations, HOMOGENIZATION_DEGREES)
    ]
    if infinity_gcd.degree() != 0 or not any(point_values):
        return None

    affine_unit, stdout, stderr = affine_standard_basis_is_unit(prime, equations)
    if stderr:
        raise ArithmeticError(("Singular stderr", weights, prime, stderr))
    if not affine_unit:
        return None

    equation_digest = hashlib.sha256(
        "\n".join(str(value) for value in equations).encode()
    ).hexdigest()
    return {
        "a_lambda": weights[0],
        "b_x0": weights[1],
        "c_x1": weights[2],
        "d": weights[3],
        "e": weights[4],
        "prime": prime,
        "sqrt5_mod_prime": sqrt5,
        "lambda_mod_prime": lam,
        "critical_equation_total_degrees": list(total_degrees),
        "critical_equation_bidegrees": [list(value) for value in bidegrees],
        "homogenization_degrees": list(HOMOGENIZATION_DEGREES),
        "infinity_chart_gcd": [1],
        "infinity_point_u0_v1_values": point_values,
        "affine_standard_basis": [1],
        "singular_stdout": stdout,
        "critical_equation_sha256": equation_digest,
    }


def histogram(rows: list[dict[str, Any]]) -> dict[str, int]:
    return {
        str(prime): count
        for prime, count in sorted(Counter(row["prime"] for row in rows).items())
    }


def build_archive(progress_every: int = 500) -> dict[str, Any]:
    records = [tuple(value) for value in canonical_black_records(2)]
    if len(records) != 11011:
        raise ArithmeticError(("five-black BB record count", len(records)))
    rows = []
    started = time.monotonic()
    for index, weights in enumerate(records, 1):
        found = None
        for prime in SPLIT_PRIMES:
            found = certify_one(prime, weights)
            if found is not None:
                break
        if found is None:
            raise ArithmeticError(("unresolved five-black BB system", weights))
        rows.append(found)
        if progress_every and index % progress_every == 0:
            print(
                json.dumps(
                    {
                        "progress": index,
                        "total": len(records),
                        "elapsed_seconds": round(time.monotonic() - started, 3),
                    },
                    sort_keys=True,
                ),
                file=sys.stderr,
                flush=True,
            )

    return {
        "schema": SCHEMA,
        "status": "EXACT_FIVE_BLACK_BB_NO_GO",
        "global_degree_29_status": "OPEN",
        "five_black_marked_passports_before": 190047,
        "five_black_binary_status_master_systems_before": 33079,
        "eliminated_block": "BB",
        "eliminated_five_black_master_systems": 11011,
        "eliminated_five_black_marked_passports": 55055,
        "remaining_five_black_master_systems": 22068,
        "remaining_five_black_marked_passports": 134992,
        "global_marked_passports_eliminated": 103398,
        "global_marked_passports_remaining": 42507280,
        "prime_histogram": histogram(rows),
        "rows": rows,
        "remaining_blocks": {
            "BW_simple": {"master_systems": 4225, "marked_passports": 21125},
            "BW_critical": {"master_systems": 4225, "marked_passports": 29575},
            "WB_simple": {"master_systems": 4225, "marked_passports": 21125},
            "WB_critical": {"master_systems": 4225, "marked_passports": 29575},
            "WW_all_binary_statuses": {
                "master_systems": 5168,
                "marked_passports": 33592,
            },
        },
        "degree_29_witness_found": False,
        "complete_five_black_no_go_claimed": False,
        "complete_degree_29_no_go_claimed": False,
        "broad_plane_tree_search_used": False,
    }


def verify_archive(archive: dict[str, Any]) -> None:
    if archive["eliminated_five_black_master_systems"] != 11011:
        raise ArithmeticError("five-black BB master count changed")
    if archive["eliminated_five_black_marked_passports"] != 55055:
        raise ArithmeticError("five-black BB passport count changed")
    if archive["remaining_five_black_master_systems"] != 22068:
        raise ArithmeticError("five-black system remainder changed")
    if archive["remaining_five_black_marked_passports"] != 134992:
        raise ArithmeticError("five-black passport remainder changed")
    if archive["global_marked_passports_remaining"] != 42507280:
        raise ArithmeticError("global passport remainder changed")
    if archive["complete_five_black_no_go_claimed"]:
        raise ArithmeticError("five-black stratum remains open")
    if archive["complete_degree_29_no_go_claimed"]:
        raise ArithmeticError("degree 29 remains open")


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--progress-every", type=int, default=500)
    arguments = parser.parse_args()
    archive = build_archive(arguments.progress_every)
    verify_archive(archive)
    encoded = canonical_bytes(archive)
    digest = hashlib.sha256(encoded).hexdigest()
    if arguments.check:
        if not arguments.output.exists() or arguments.output.read_bytes() != encoded:
            raise SystemExit("five-black BB archive is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": archive["status"],
                "eliminated_master_systems": archive[
                    "eliminated_five_black_master_systems"
                ],
                "eliminated_marked_passports": archive[
                    "eliminated_five_black_marked_passports"
                ],
                "remaining_master_systems": archive[
                    "remaining_five_black_master_systems"
                ],
                "remaining_marked_passports": archive[
                    "remaining_five_black_marked_passports"
                ],
                "prime_histogram": archive["prime_histogram"],
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
