#!/usr/bin/env python3
"""Exact attack on the final four-black BW/WB simple-white block.

Orient the black endpoint to x=0.  If the marked black root is alpha, scale
y=x/alpha.  The black roots become 0,1,u,v and the white endpoint becomes

    z_BW = 1/lambda,       z_WB = 1/(1-lambda).

The critical configurations depend only on the four multiplicities, not on
the orientation.  They are enumerated exactly through the degree-756
resultant and finite-field tower factorization used by the four-black WW
certificate.  On each irreducible (u,v) orbit this script tests whether

    G(z_orientation) - R0(u,v)

vanishes in the tower field.  A nonzero value on every orbit excludes the
orientation.  Endpoint simplicity is not assumed, so a complete no-go here
closes both simple and critical refinements (the critical refinements were
already eliminated independently).
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx

from certify_gram_belyi_m5_degree29_reduction import lambda_mod
from certify_gram_belyi_m5_degree29_four_black_reduction import (
    one_black_endpoint_records,
)
from certify_gram_belyi_m5_degree29_four_black_ww_simple import (
    canonical_remainder,
    mpoly_to_nmod_poly,
    specialize_u_as_v_poly,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/gram_belyi_m5_degree29_four_black_one_black_simple.json"
SCHEMA = "coxeter.gram_belyi.m5_degree29_four_black_one_black_simple.v1"
SPLIT_PRIMES = (101, 109, 131, 139, 149, 151, 179, 181, 191, 199)


def specialize_r0_in_tower(
    r0, theta, field, v_value, v_modulus, polynomial_context
):
    answer = polynomial_context(0)
    for (u_power, v_power), coefficient in r0.terms():
        base = polynomial_context(
            [field(int(coefficient)) * (theta**u_power)]
        )
        answer = (
            answer
            + base.mul_mod(v_value.pow_mod(v_power, v_modulus), v_modulus)
        ) % v_modulus
    return answer


def endpoint_difference(
    r0,
    canonical_weights: tuple[int, int, int, int],
    endpoint: int,
    theta,
    field,
    v_modulus,
    polynomial_context,
):
    """Return G(endpoint)-R0 as an element of K[v]/(v_modulus)."""

    v_value = polynomial_context.gen() % v_modulus
    one = polynomial_context(1)
    z_value = polynomial_context([field(endpoint)])
    u_value = polynomial_context([theta])
    factors = (
        z_value,
        z_value - one,
        z_value - u_value,
        z_value - v_value,
    )
    g_value = one
    for factor, exponent in zip(factors, canonical_weights):
        g_value = g_value.mul_mod(
            factor.pow_mod(exponent, v_modulus), v_modulus
        )
    critical_value = specialize_r0_in_tower(
        r0, theta, field, v_value, v_modulus, polynomial_context
    )
    return (g_value - critical_value) % v_modulus


def certify_one(
    prime: int, degrees: tuple[int, int, int, int]
) -> dict[str, Any] | None:
    """Try one split-prime certificate for both endpoint orientations."""

    a, b, c, d = degrees
    sqrt5, lam = lambda_mod(prime)
    endpoints = {
        "BW": pow(lam, -1, prime),
        "WB": pow((1 - lam) % prime, -1, prime),
    }
    # Canonical roots are 0 (black endpoint), 1 (marked root), u, v.
    canonical_weights = (b, a, c, d)
    _context, remainder = canonical_remainder(
        prime, lam, canonical_weights
    )
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
    critical_orbits = 0
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
                    "endpoint_nonzero_BW_WB": [],
                }
            )
            continue

        _v_unit, v_factors = common_v.factor()
        v_factor_degrees = []
        endpoint_checks = []
        for v_factor, v_exponent in v_factors:
            differences = [
                endpoint_difference(
                    r0,
                    canonical_weights,
                    endpoints[orientation],
                    theta,
                    field,
                    v_factor,
                    polynomial_context,
                )
                for orientation in ("BW", "WB")
            ]
            nonzero = [bool(value) for value in differences]
            v_factor_degrees.append([v_factor.degree(), v_exponent])
            endpoint_checks.append(nonzero)
            critical_orbits += 1
            if not all(nonzero):
                return None
        component_rows.append(
            {
                "u_factor_degree": u_factor.degree(),
                "u_factor_exponent": u_exponent,
                "specialized_R1_R2_degrees": [s1.degree(), s2.degree()],
                "common_v_degree": common_v.degree(),
                "v_factor_degrees": v_factor_degrees,
                "endpoint_nonzero_BW_WB": endpoint_checks,
            }
        )

    return {
        "a_lambda": a,
        "b_black_endpoint": b,
        "c": c,
        "d": d,
        "canonical_weights_0_1_u_v": list(canonical_weights),
        "prime": prime,
        "sqrt5_mod_prime": sqrt5,
        "lambda_mod_prime": lam,
        "endpoint_coordinates_mod_prime": endpoints,
        "critical_equation_total_degrees": [28, 27],
        "critical_equation_bidegrees": [[27, 27], [27, 27]],
        "critical_resultant_degree": 756,
        "critical_resultant_monic_coefficients": [
            int(value) for value in resultant_monic
        ],
        "u_factor_count": len(u_factors),
        "critical_orbit_count": critical_orbits,
        "components": component_rows,
        "all_BW_WB_endpoint_differences_nonzero": True,
    }


def build_archive() -> dict[str, Any]:
    rows = []
    for degrees in one_black_endpoint_records():
        found = None
        for prime in SPLIT_PRIMES:
            found = certify_one(prime, degrees)
            if found is not None:
                break
        if found is None:
            raise ArithmeticError(("unresolved one-black system", degrees))
        rows.append(found)

    return {
        "schema": SCHEMA,
        "status": "EXACT_FOUR_BLACK_COMPLETE_NO_GO",
        "global_degree_29_status": "OPEN",
        "black_endpoint_degree_record_count": 1729,
        "endpoint_orientations": ["BW", "WB"],
        "broad_orientation_system_count": 3458,
        "remaining_simple_endpoint_status_system_count_before": 3458,
        "remaining_simple_marked_passport_count_before": 10374,
        "eliminated_simple_endpoint_status_system_count": 3458,
        "eliminated_simple_marked_passport_count": 10374,
        "four_black_coefficient_systems_remaining": 0,
        "four_black_marked_passports_remaining": 0,
        "global_marked_passports_eliminated": 48343,
        "global_marked_passports_remaining": 42562335,
        "prime_histogram": dict(
            (str(prime), count)
            for prime, count in sorted(Counter(row["prime"] for row in rows).items())
        ),
        "total_critical_orbits_checked": sum(
            row["critical_orbit_count"] for row in rows
        ),
        "rows": rows,
        "degree_29_witness_found": False,
        "complete_four_black_no_go_claimed": True,
        "complete_degree_29_no_go_claimed": False,
        "broad_plane_tree_search_used": False,
    }


def verify_archive(archive: dict[str, Any]) -> None:
    if archive["black_endpoint_degree_record_count"] != 1729:
        raise ArithmeticError("degree-record count changed")
    if archive["broad_orientation_system_count"] != 3458:
        raise ArithmeticError("orientation count changed")
    if archive["eliminated_simple_endpoint_status_system_count"] != 3458:
        raise ArithmeticError("final four-black block did not close")
    if archive["eliminated_simple_marked_passport_count"] != 10374:
        raise ArithmeticError("passport conversion changed")
    if archive["four_black_marked_passports_remaining"] != 0:
        raise ArithmeticError("four-black passports remain")
    if archive["global_marked_passports_remaining"] != 42562335:
        raise ArithmeticError("global boundary changed")
    if not archive["complete_four_black_no_go_claimed"]:
        raise ArithmeticError("four-black completion flag missing")
    if archive["complete_degree_29_no_go_claimed"]:
        raise ArithmeticError("higher black-vertex strata remain")


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
            raise SystemExit("four-black one-black/simple archive is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": archive["status"],
                "eliminated_systems": archive[
                    "eliminated_simple_endpoint_status_system_count"
                ],
                "eliminated_marked_passports": archive[
                    "eliminated_simple_marked_passport_count"
                ],
                "total_critical_orbits_checked": archive[
                    "total_critical_orbits_checked"
                ],
                "prime_histogram": archive["prime_histogram"],
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
