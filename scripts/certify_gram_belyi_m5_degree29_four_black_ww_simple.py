#!/usr/bin/env python3
"""Exact no-go certificate for every four-black WW endpoint status.

For each of the 640 marked black-degree records, translate the marked black
root to 0 and scale one other root to 1:

    G(y) = y^a (y-1)^b (y-u)^c (y-v)^d.

If D is the residual derivative cubic and

    G mod D = R2*y^2 + R1*y + R0,

then the Shabat condition is R1=R2=0.  The two normalized white endpoints
have fixed ratio q=29*lambda-13, so they require a common root of

    H(y)=G(y)-R0,  H(q*y).

At a good split prime this script:

1. forms A(u)=Res_v(R1,R2), of the sharp Bezout degree 28*27=756;
2. factors A and, over each residue field F_p[u]/(f), computes
   gcd_v(R1,R2);
3. factors that gcd and, in each resulting tower field, computes the full
   gcd_y(H(y),H(q*y)).

Every endpoint gcd is constant.  Degree 756 makes u integral locally; at
least one critical equation retains degree 27 in v on every relevant
component, making v integral; H is monic of degree 29, making y integral.
Thus reduction cannot lose a characteristic-zero solution at infinity.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

from flint import (
    fmpz_mod_poly_ctx,
    fq_default_ctx,
    fq_default_poly_ctx,
    nmod_mpoly_ctx,
    nmod_poly,
)

from certify_gram_belyi_m5_degree29_reduction import FpPoly, lambda_mod
from certify_gram_belyi_m5_degree29_four_black_reduction import (
    four_black_degree_records,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/gram_belyi_m5_degree29_four_black_ww_simple.json"
SCHEMA = "coxeter.gram_belyi.m5_degree29_four_black_ww_simple.v1"
SPLIT_PRIMES = (101, 109, 131, 139, 149, 151, 179, 181, 191, 199)


def trim_outer(value):
    while len(value) > 1 and not value[-1]:
        value.pop()
    return value


def outer_add(left, right):
    context = (left or right)[0].context()
    answer = [
        (left[index] if index < len(left) else context.constant(0))
        + (right[index] if index < len(right) else context.constant(0))
        for index in range(max(len(left), len(right)))
    ]
    return trim_outer(answer)


def outer_mul(left, right):
    context = left[0].context()
    answer = [context.constant(0) for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] += a * b
    return trim_outer(answer)


def outer_power(value, exponent: int):
    context = value[0].context()
    answer = [context.constant(1)]
    base = value
    while exponent:
        if exponent & 1:
            answer = outer_mul(answer, base)
        base = outer_mul(base, base)
        exponent //= 2
    return answer


def canonical_remainder(
    prime: int, lam: int, weights: tuple[int, int, int, int]
):
    """Return the coefficient context and R0,R1,R2 in F_p[u,v]."""

    context = nmod_mpoly_ctx.get(("u", "v"), prime, "degrevlex")
    u, v = context.gens()
    roots = (context.constant(0), context.constant(1), u, v)
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

    inverse29 = pow(29, -1, prime)
    remainder = list(product_value)
    while len(remainder) >= 4:
        quotient_coefficient = inverse29 * remainder[-1]
        shift = len(remainder) - 4
        for index, coefficient in enumerate(derivative_residual):
            remainder[shift + index] -= quotient_coefficient * coefficient
        trim_outer(remainder)
    if len(remainder) != 3:
        remainder += [context.constant(0)] * (3 - len(remainder))
    return context, remainder


def mpoly_to_nmod_poly(value, prime: int, variable: int = 0) -> nmod_poly:
    terms = value.to_dict()
    degree = max((powers[variable] for powers in terms), default=0)
    coefficients = [0] * (degree + 1)
    for powers, coefficient in terms.items():
        if any(
            powers[index]
            for index in range(len(powers))
            if index != variable
        ):
            raise ValueError((powers, value))
        coefficients[powers[variable]] = int(coefficient)
    return nmod_poly(coefficients, prime)


def specialize_u_as_v_poly(value, theta, field, polynomial_context):
    """Specialize an F_p[u,v] polynomial at u=theta, retaining v."""

    coefficients = [field(0)] * (value.degrees()[1] + 1)
    for (u_power, v_power), coefficient in value.terms():
        coefficients[v_power] += field(int(coefficient)) * (theta**u_power)
    return polynomial_context(coefficients)


def tower_trim(value):
    while len(value) > 1 and not value[-1]:
        value.pop()
    return value


def tower_add(left, right, zero):
    answer = [
        (left[index] if index < len(left) else zero)
        + (right[index] if index < len(right) else zero)
        for index in range(max(len(left), len(right)))
    ]
    return tower_trim(answer)


def tower_neg(value):
    return [-coefficient for coefficient in value]


def tower_sub(left, right, zero):
    return tower_add(left, tower_neg(right), zero)


def tower_mul(left, right, modulus, zero):
    answer = [zero for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] = (answer[i + j] + a.mul_mod(b, modulus)) % modulus
    return tower_trim(answer)


def tower_power(value, exponent: int, modulus, one, zero):
    answer = [one]
    base = value
    while exponent:
        if exponent & 1:
            answer = tower_mul(answer, base, modulus, zero)
        base = tower_mul(base, base, modulus, zero)
        exponent //= 2
    return answer


def tower_remainder(dividend, divisor, modulus, zero):
    left = tower_trim(list(dividend))
    right = tower_trim(list(divisor))
    inverse_leading = right[-1].inverse_mod(modulus)
    while not (len(left) == 1 and not left[0]) and len(left) >= len(right):
        coefficient = left[-1].mul_mod(inverse_leading, modulus)
        shift = len(left) - len(right)
        for index, value in enumerate(right):
            left[shift + index] = (
                left[shift + index] - coefficient.mul_mod(value, modulus)
            ) % modulus
        tower_trim(left)
    return left


def tower_gcd(left, right, modulus, one, zero):
    a, b = tower_trim(list(left)), tower_trim(list(right))
    while not (len(b) == 1 and not b[0]):
        a, b = b, tower_remainder(a, b, modulus, zero)
    inverse = a[-1].inverse_mod(modulus)
    return [coefficient.mul_mod(inverse, modulus) for coefficient in a]


def endpoint_gcd_degree(
    r0,
    weights: tuple[int, int, int, int],
    q: int,
    theta,
    field,
    v_modulus,
    polynomial_context,
) -> int:
    """Compute gcd_y(H(y),H(qy)) in the tower K[v]/(v_modulus)."""

    zero = polynomial_context(0)
    one = polynomial_context(1)
    u_value = theta
    v_value = polynomial_context.gen() % v_modulus
    u_tower = polynomial_context([u_value])

    factors = (
        [zero, one],
        [-one, one],
        [-u_tower, one],
        [-v_value, one],
    )
    g_value = [one]
    for factor, exponent in zip(factors, weights):
        g_value = tower_mul(
            g_value,
            tower_power(factor, exponent, v_modulus, one, zero),
            v_modulus,
            zero,
        )

    critical_value = polynomial_context(0)
    for (u_power, v_power), coefficient in r0.terms():
        base = polynomial_context(
            [field(int(coefficient)) * (theta**u_power)]
        )
        critical_value = (
            critical_value
            + base.mul_mod(v_value.pow_mod(v_power, v_modulus), v_modulus)
        ) % v_modulus
    h_value = list(g_value)
    h_value[0] = (h_value[0] - critical_value) % v_modulus

    q_tower = polynomial_context([field(q)])
    h_scaled = []
    q_power = one
    for coefficient in h_value:
        h_scaled.append(coefficient.mul_mod(q_power, v_modulus))
        q_power = q_power.mul_mod(q_tower, v_modulus)
    gcd = tower_gcd(h_value, h_scaled, v_modulus, one, zero)
    return len(gcd) - 1


def certify_one(
    prime: int, weights: tuple[int, int, int, int]
) -> dict[str, Any] | None:
    sqrt5, lam = lambda_mod(prime)
    q = (29 * lam - 13) % prime
    _context, remainder = canonical_remainder(prime, lam, weights)
    r0, r1, r2 = remainder

    if r1.total_degree() != 28 or r2.total_degree() != 27:
        return None
    if r1.degrees() != (27, 27) or r2.degrees() != (27, 27):
        return None

    resultant_mpoly = r1.resultant(r2, "v")
    if resultant_mpoly.degrees() != (756, 0):
        return None
    resultant = mpoly_to_nmod_poly(resultant_mpoly, prime)
    if resultant.degree() != 756:
        return None
    resultant_monic = resultant * pow(
        int(resultant[resultant.degree()]), -1, prime
    )
    _unit, u_factors = resultant.factor()

    integer_poly_context = fmpz_mod_poly_ctx(prime)
    component_rows = []
    endpoint_orbits = 0
    for u_factor, u_exponent in u_factors:
        modulus = integer_poly_context([int(value) for value in u_factor])
        field = fq_default_ctx(modulus=modulus, var="u", fq_type="FQ_NMOD")
        polynomial_context = fq_default_poly_ctx(field)
        theta = field.gen()
        s1 = specialize_u_as_v_poly(r1, theta, field, polynomial_context)
        s2 = specialize_u_as_v_poly(r2, theta, field, polynomial_context)
        if max(s1.degree(), s2.degree()) != 27:
            return None
        common_v = s1.gcd(s2)
        if common_v.degree() == 0:
            component_rows.append(
                {
                    "u_factor_degree": u_factor.degree(),
                    "u_factor_exponent": u_exponent,
                    "specialized_R1_R2_degrees": [s1.degree(), s2.degree()],
                    "common_v_degree": 0,
                    "v_factor_degrees": [],
                    "endpoint_gcd_degrees": [],
                }
            )
            continue

        _v_unit, v_factors = common_v.factor()
        v_factor_degrees = []
        endpoint_degrees = []
        for v_factor, v_exponent in v_factors:
            degree = endpoint_gcd_degree(
                r0,
                weights,
                q,
                theta,
                field,
                v_factor,
                polynomial_context,
            )
            v_factor_degrees.append([v_factor.degree(), v_exponent])
            endpoint_degrees.append(degree)
            endpoint_orbits += 1
            if degree != 0:
                return None
        component_rows.append(
            {
                "u_factor_degree": u_factor.degree(),
                "u_factor_exponent": u_exponent,
                "specialized_R1_R2_degrees": [s1.degree(), s2.degree()],
                "common_v_degree": common_v.degree(),
                "v_factor_degrees": v_factor_degrees,
                "endpoint_gcd_degrees": endpoint_degrees,
            }
        )

    return {
        "a": weights[0],
        "b": weights[1],
        "c": weights[2],
        "d": weights[3],
        "prime": prime,
        "sqrt5_mod_prime": sqrt5,
        "lambda_mod_prime": lam,
        "endpoint_ratio_q_mod_prime": q,
        "critical_equation_total_degrees": [28, 27],
        "critical_equation_bidegrees": [[27, 27], [27, 27]],
        "critical_resultant_degree": 756,
        "critical_resultant_monic_coefficients": [
            int(value) for value in resultant_monic
        ],
        "u_factor_count": len(u_factors),
        "endpoint_orbit_count": endpoint_orbits,
        "components": component_rows,
        "all_endpoint_gcds_are_one": True,
    }


def build_archive() -> dict[str, Any]:
    rows = []
    for weights in four_black_degree_records():
        found = None
        for prime in SPLIT_PRIMES:
            found = certify_one(prime, weights)
            if found is not None:
                break
        if found is None:
            raise ArithmeticError(("unresolved four-black WW simple", weights))
        rows.append(found)

    return {
        "schema": SCHEMA,
        "status": "EXACT_FOUR_BLACK_WW_NO_GO",
        "global_degree_29_status": "OPEN",
        "black_degree_master_system_count": 640,
        "represented_endpoint_statuses": ["SS", "SC", "CS", "CC"],
        "represented_endpoint_status_system_count": 2560,
        "represented_marked_passport_count": 8960,
        "eliminated_black_degree_master_system_count": len(rows),
        "eliminated_endpoint_status_system_count": 4 * len(rows),
        "eliminated_marked_passport_count": 14 * len(rows),
        "remaining_black_degree_master_system_count": 640 - len(rows),
        "remaining_endpoint_status_system_count": 2560 - 4 * len(rows),
        "remaining_marked_passport_count": 8960 - 14 * len(rows),
        "four_black_coefficient_systems_remaining": 3458,
        "four_black_marked_passports_remaining": 10374,
        "global_marked_passports_eliminated": 37969,
        "global_marked_passports_remaining": 42572709,
        "prime_histogram": dict(
            (str(prime), count)
            for prime, count in sorted(Counter(row["prime"] for row in rows).items())
        ),
        "total_endpoint_orbits_checked": sum(
            row["endpoint_orbit_count"] for row in rows
        ),
        "rows": rows,
        "degree_29_witness_found": False,
        "complete_four_black_no_go_claimed": False,
        "complete_degree_29_no_go_claimed": False,
        "broad_plane_tree_search_used": False,
    }


def verify_archive(archive: dict[str, Any]) -> None:
    if archive["black_degree_master_system_count"] != 640:
        raise ArithmeticError("system count changed")
    if archive["eliminated_black_degree_master_system_count"] != 640:
        raise ArithmeticError("WW block did not close")
    if archive["eliminated_endpoint_status_system_count"] != 2560:
        raise ArithmeticError("endpoint-status conversion changed")
    if archive["eliminated_marked_passport_count"] != 8960:
        raise ArithmeticError("passport conversion changed")
    if archive["remaining_black_degree_master_system_count"] != 0:
        raise ArithmeticError("WW systems remain")
    if archive["four_black_marked_passports_remaining"] != 10374:
        raise ArithmeticError("four-black boundary changed")
    if archive["global_marked_passports_remaining"] != 42572709:
        raise ArithmeticError("global boundary changed")
    if archive["complete_four_black_no_go_claimed"]:
        raise ArithmeticError("other four-black blocks remain")


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    archive = build_archive()
    verify_archive(archive)
    encoded = canonical_bytes(archive)
    digest = hashlib.sha256(encoded).hexdigest()
    if arguments.check:
        if not arguments.output.exists() or arguments.output.read_bytes() != encoded:
            raise SystemExit("four-black WW simple archive is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": archive["status"],
                "eliminated_black_degree_master_systems": archive[
                    "eliminated_black_degree_master_system_count"
                ],
                "eliminated_endpoint_status_systems": archive[
                    "eliminated_endpoint_status_system_count"
                ],
                "eliminated_marked_passports": archive[
                    "eliminated_marked_passport_count"
                ],
                "total_endpoint_orbits_checked": archive[
                    "total_endpoint_orbits_checked"
                ],
                "prime_histogram": archive["prime_histogram"],
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
