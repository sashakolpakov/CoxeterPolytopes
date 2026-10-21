#!/usr/bin/env python3
"""Bounded exact elimination in the three-black/WW degree-29 stratum.

This file does not enumerate plane trees.  The 1,372 marked passport records
collapse to 784 distinct coefficient systems.  First, the 980 records in
which at least one normalised endpoint is a non-root critical point collapse
to 588 systems:

* 196 systems with both endpoints critical;
* 196 systems with only x=0 required critical;
* 196 systems with only x=1 required critical.

For each system, one root is eliminated linearly from the residual quadratic
factor of the logarithmic derivative.  The remaining necessary equations are
univariate.  Degree-preserving reduction at a split prime followed by a full
polynomial gcd is an exact non-vanishing-resultant certificate, not a search
for F_p-rational points.

The other 392 records (the two possible white partitions with simple/simple
endpoints) are the same 196 coefficient systems.  Translation and scaling
put the black roots at 0,1,t.  The endpoint ratio becomes

    q = -(1-lambda)/lambda = 29*lambda-13,  q^2+12*q+16=0.

The common-critical-value condition is one degree-28 equation S(t)=0.  At a
good split prime we factor S completely.  In each residue field factor we
compute the full gcd of H(y)=F(y)-R_0(t) and H(q*y).  Every gcd is one.

The modular conclusion is exact for a particularly simple reason: S has its
full characteristic-zero degree 28 at the displayed prime, so t is integral
at that prime after scaling S to monic; H is monic of degree 29 in y, so y is
integral too.  Any characteristic-zero solution would therefore reduce to a
solution over the algebraic closure of the residue field.  None exists.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from flint import fmpz_mod_poly_ctx, fq_default_ctx, fq_default_poly_ctx

from certify_gram_belyi_m5_degree29_reduction import (
    FpPoly,
    lambda_mod,
    xr_poly_product,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/gram_belyi_m5_degree29_three_black_ww.json"
SCHEMA = "coxeter.gram_belyi.m5_degree29_three_black_ww.v1"

# All are split for Q(sqrt(5)); primes dividing 29 or a possible constant
# denominator are rejected case by case.
SPLIT_PRIMES = (
    101,
    109,
    131,
    139,
    149,
    151,
    179,
    181,
    191,
    199,
    211,
    229,
    239,
    241,
    251,
    269,
    271,
    281,
    311,
    349,
)


def black_degree_records() -> list[tuple[int, int, int]]:
    """Return (a,b,c), with lambda of degree a and b<=c."""

    answer = []
    for a in range(1, 28):
        for b in range(1, (29 - a) // 2 + 1):
            c = 29 - a - b
            if b <= c:
                answer.append((a, b, c))
    if len(answer) != 196:
        raise ArithmeticError(len(answer))
    return answer


def product(ring: FpPoly, factors: Iterable[tuple[list[int], int]]) -> list[int]:
    answer = [1]
    for factor, exponent in factors:
        answer = ring.mul(answer, ring.power(factor, exponent))
    return answer


def both_endpoints_critical_equations(
    prime: int, lam: int, a: int, b: int, c: int
) -> tuple[list[int], list[int], tuple[int, int]]:
    """Equations after D=29*x*(x-1) and linear elimination of s.

    If F=(x-lambda)^a(x-r)^b(x-s)^c, write

      D=29*x^2-A*x+B,
      A=(b+c)lambda+(a+c)r+(a+b)s,
      B=a*r*s+b*lambda*s+c*lambda*r.

    Both endpoints critical forces A=29 and B=0.  The first equation solves
    s=S/(a+b); substitution in the second gives a quadratic q(r).  Equality
    F(0)=F(1) gives e(r).
    """

    ring = FpPoly(prime)
    denominator = (a + b) % prime
    if denominator == 0:
        raise ZeroDivisionError((prime, a, b))
    r = [0, 1]
    # S = 29-(b+c)lambda-(a+c)r, so s=S/(a+b).
    s_num = [(29 - (b + c) * lam) % prime, -(a + c)]

    q = ring.add(
        ring.add(
            ring.scale(ring.mul(r, s_num), a),
            ring.scale(s_num, b * lam),
        ),
        ring.scale(r, c * lam * denominator),
    )

    f0_num = product(
        ring,
        (
            ([-lam], a),
            ([0, -1], b),
            (ring.neg(s_num), c),
        ),
    )
    f1_num = product(
        ring,
        (
            ([(1 - lam) % prime], a),
            ([1, -1], b),
            (ring.sub([denominator], s_num), c),
        ),
    )
    endpoint = ring.sub(f1_num, f0_num)
    # The endpoint polynomial has degree b+c=29-a: its leading coefficient
    # is a nonzero multiple of (1-lambda)^a-(-lambda)^a.  This cannot vanish
    # because 0<lambda<1 and a>0 (for odd a the signs differ; for even a
    # equality would force lambda=1/2).
    return q, endpoint, (2, 29 - a)


def endpoint_critical_equations(
    prime: int, lam: int, a: int, b: int, c: int, endpoint: int
) -> tuple[list[int], list[int], list[int], tuple[int, int], int]:
    """Eliminate s from D(endpoint)=0 and return two univariate equations.

    Write D(e)=C(r)s+K(r).  On C!=0, s=-K/C.  The returned equations are the
    cleared numerators of F(1-e)-F(e) and F(z)-F(e), where z is the other
    root of D.  The scalar is Res_r(C,K); its non-vanishing excludes the only
    branch lost by division by C.
    """

    if endpoint not in (0, 1):
        raise ValueError(endpoint)
    ring = FpPoly(prime)
    r = [0, 1]
    a0 = [(b + c) * lam % prime, a + c]

    if endpoint == 0:
        # D(0)=B=(a*r+b*lambda)s+c*lambda*r.
        c_poly = [b * lam, a]
        k_poly = [0, c * lam]
    else:
        # D(1)=29-A+B.
        c_poly = [(b * lam - a - b) % prime, a]
        k_poly = [
            (29 - (b + c) * lam) % prime,
            (c * lam - a - c) % prime,
        ]

    # Resultant of C=c0+c1*r and K=k0+k1*r.  If it is nonzero, C and K
    # cannot vanish together over the algebraic closure.
    exceptional_resultant = (
        c_poly[1] * k_poly[0] - c_poly[0] * k_poly[1]
    ) % prime

    other = 1 - endpoint
    # x-s=(x*C+K)/C after s=-K/C.
    def fixed_value_numerator(x_value: int) -> list[int]:
        return product(
            ring,
            (
                ([(x_value - lam) % prime], a),
                ([x_value, -1], b),
                (ring.add(ring.scale(c_poly, x_value), k_poly), c),
            ),
        )

    f_endpoint = fixed_value_numerator(endpoint)
    f_other = fixed_value_numerator(other)
    endpoint_equation = ring.sub(f_other, f_endpoint)

    # z=A/29-e and A=A0+(a+b)s.  Thus
    # z=Z/(29*C), Z=A0*C-(a+b)K-29e*C.
    z_num = ring.sub(
        ring.sub(ring.mul(a0, c_poly), ring.scale(k_poly, a + b)),
        ring.scale(c_poly, 29 * endpoint),
    )
    z_minus_lambda = ring.sub(z_num, ring.scale(c_poly, 29 * lam))
    z_minus_r = ring.sub(z_num, ring.scale(ring.mul(r, c_poly), 29))
    z_minus_s = ring.add(z_num, ring.scale(k_poly, 29))
    f_z_num = product(
        ring,
        ((z_minus_lambda, a), (z_minus_r, b), (z_minus_s, c)),
    )
    scaled_f_endpoint = ring.scale(
        ring.mul(ring.power(c_poly, a + b), f_endpoint),
        pow(29, 29, prime),
    )
    critical_equation = ring.sub(f_z_num, scaled_f_endpoint)

    return (
        endpoint_equation,
        critical_equation,
        c_poly,
        (29 - a, 58),
        exceptional_resultant,
    )


def separating_prime_both(a: int, b: int, c: int) -> dict[str, Any] | None:
    for prime in SPLIT_PRIMES:
        if (a + b) % prime == 0:
            continue
        sqrt5, lam = lambda_mod(prime)
        ring = FpPoly(prime)
        first, second, degree_bounds = both_endpoints_critical_equations(
            prime, lam, a, b, c
        )
        degrees = (len(first) - 1, len(second) - 1)
        if degrees != degree_bounds:
            continue
        if ring.gcd(first, second) == [1]:
            return {
                "prime": prime,
                "sqrt5_mod_prime": sqrt5,
                "lambda_mod_prime": lam,
                "equation_degrees": list(degrees),
                "gcd": [1],
            }
    return None


def separating_prime_one(
    a: int, b: int, c: int, endpoint: int
) -> dict[str, Any] | None:
    for prime in SPLIT_PRIMES:
        sqrt5, lam = lambda_mod(prime)
        ring = FpPoly(prime)
        first, second, _c_poly, degree_bounds, exceptional = (
            endpoint_critical_equations(prime, lam, a, b, c, endpoint)
        )
        degrees = (len(first) - 1, len(second) - 1)
        if exceptional == 0 or degrees != degree_bounds:
            continue
        if ring.gcd(first, second) == [1]:
            return {
                "prime": prime,
                "sqrt5_mod_prime": sqrt5,
                "lambda_mod_prime": lam,
                "equation_degrees": list(degrees),
                "exceptional_branch_resultant_mod_prime": exceptional,
                "gcd": [1],
            }
    return None


def canonical_critical_remainder(
    prime: int, a: int, b: int, c: int
) -> tuple[list[int], list[int]]:
    """Return R_0(t),R_1(t) for y^a(y-1)^b(y-t)^c modulo D.

    Here

      D=29*y^2-(a+c+(a+b)t)y+a*t.

    Thus, modulo the monic form of D,

      y^2=((a+c+(a+b)t)/29)y-a*t/29.
    """

    ring = FpPoly(prime)
    inverse29 = pow(29, -1, prime)
    t = [0, 1]
    u = [(a + c) * inverse29 % prime, (a + b) * inverse29 % prime]
    v = [0, -a * inverse29 % prime]
    return xr_poly_product(
        (
            (([0], [1]), a),
            (([-1], [1]), b),
            ((ring.neg(t), [1]), c),
        ),
        u,
        v,
        ring,
    )


def evaluate_in_extension(coefficients, theta, field):
    answer = field(0)
    for coefficient in reversed(coefficients):
        answer = answer * theta + field(int(coefficient))
    return answer


def separating_prime_simple_simple(
    a: int, b: int, c: int
) -> dict[str, Any] | None:
    """Exclude one canonical simple/simple system over a split residue field."""

    for prime in SPLIT_PRIMES:
        sqrt5, lam = lambda_mod(prime)
        q = (29 * lam - 13) % prime
        if q == 0:
            continue
        r0, r1 = canonical_critical_remainder(prime, a, b, c)
        # Homogeneous quadratic reduction gives deg S<=28.  Full degree here
        # proves that its characteristic-zero leading coefficient is a unit
        # at the selected prime.
        if len(r1) - 1 != 28:
            continue

        residue_poly_context = fmpz_mod_poly_ctx(prime)
        shabat = residue_poly_context(r1)
        _unit, factors = shabat.factor()
        factor_rows = []
        prime_separates = True

        for factor, exponent in factors:
            field = fq_default_ctx(
                modulus=factor, var="t", fq_type="FQ_NMOD"
            )
            polynomial_context = fq_default_poly_ctx(field)
            theta = field.gen()
            y = polynomial_context.gen()
            f_value = (y**a) * ((y - field(1)) ** b) * ((y - theta) ** c)
            critical_value = evaluate_in_extension(r0, theta, field)
            h_value = f_value - critical_value
            q_field = field(q)
            h_scaled = polynomial_context(
                [coefficient * (q_field**index) for index, coefficient in enumerate(h_value)]
            )
            gcd = h_value.gcd(h_scaled)
            factor_rows.append(
                {
                    "factor_coefficients": [int(value) for value in factor],
                    "factor_degree": factor.degree(),
                    "factor_exponent": exponent,
                    "H_degrees": [h_value.degree(), h_scaled.degree()],
                    "gcd_degree": gcd.degree(),
                }
            )
            if h_value.degree() != 29 or h_scaled.degree() != 29 or gcd.degree() != 0:
                prime_separates = False
                break

        if prime_separates and sum(
            row["factor_degree"] * row["factor_exponent"] for row in factor_rows
        ) == 28:
            return {
                "prime": prime,
                "sqrt5_mod_prime": sqrt5,
                "lambda_mod_prime": lam,
                "endpoint_ratio_q_mod_prime": q,
                "shabat_equation_degree": 28,
                "critical_remainder_coefficients_mod_prime": r0,
                "shabat_equation_coefficients_mod_prime": r1,
                "irreducible_factor_checks": factor_rows,
                "all_extension_gcds_are_one": True,
            }
    return None


def build_archive() -> dict[str, Any]:
    records = black_degree_records()
    rows_both = []
    rows_one = []
    rows_simple = []
    unresolved = []

    for a, b, c in records:
        certificate = separating_prime_both(a, b, c)
        row = {"a": a, "b": b, "c": c, "endpoint_case": "both critical"}
        if certificate is None:
            unresolved.append(row)
        else:
            rows_both.append({**row, **certificate})

        for endpoint in (0, 1):
            certificate = separating_prime_one(a, b, c, endpoint)
            row = {
                "a": a,
                "b": b,
                "c": c,
                "endpoint_case": f"x={endpoint} critical",
            }
            if certificate is None:
                unresolved.append(row)
            else:
                rows_one.append({**row, **certificate})

        certificate = separating_prime_simple_simple(a, b, c)
        row = {"a": a, "b": b, "c": c, "endpoint_case": "both simple"}
        if certificate is None:
            unresolved.append(row)
        else:
            rows_simple.append({**row, **certificate})

    eliminated_systems = len(rows_both) + len(rows_one)
    eliminated_records = 0
    if eliminated_systems == 588:
        # A one-critical coefficient system covers two marked passports:
        # the exceptional endpoint may have white degree 2 or 3.  A
        # both-critical system covers the single (2,2) endpoint record.
        eliminated_records = 2 * len(rows_one) + len(rows_both)

    eliminated_simple_records = 2 * len(rows_simple)

    all_rows = rows_both + rows_one
    archive = {
        "schema": SCHEMA,
        "status": (
            "EXACT_THREE_BLACK_WW_NO_GO"
            if not unresolved
            else "PARTIAL_THREE_BLACK_WW_REDUCTION"
        ),
        "global_degree_29_status": "OPEN",
        "marked_three_black_WW_records_before": 1372,
        "marked_black_degree_records": len(records),
        "distinct_three_black_WW_coefficient_systems": 784,
        "distinct_endpoint_critical_systems": 588,
        "endpoint_critical_marked_records": 980,
        "eliminated_endpoint_critical_systems": eliminated_systems,
        "eliminated_endpoint_critical_marked_records": eliminated_records,
        "simple_simple_marked_records": 392,
        "distinct_simple_simple_systems": 196,
        "eliminated_simple_simple_systems": len(rows_simple),
        "eliminated_simple_simple_marked_records": eliminated_simple_records,
        "eliminated_three_black_WW_coefficient_systems": (
            eliminated_systems + len(rows_simple)
        ),
        "eliminated_three_black_WW_marked_records": (
            eliminated_records + eliminated_simple_records
        ),
        "remaining_three_black_WW_marked_records": (
            1372 - eliminated_records - eliminated_simple_records
        ),
        "remaining_degree_29_marked_records_with_at_least_four_black_vertices": 42605329,
        "degree_29_witness_found": False,
        "complete_degree_29_no_go_claimed": False,
        "method": (
            "linear root elimination; degree-preserving split-prime "
            "univariate gcd/resultant certificates"
        ),
        "endpoint_critical_prime_histogram": dict(
            (str(prime), count)
            for prime, count in sorted(
                Counter(row["prime"] for row in all_rows).items()
            )
        ),
        "simple_simple_prime_histogram": dict(
            (str(prime), count)
            for prime, count in sorted(
                Counter(row["prime"] for row in rows_simple).items()
            )
        ),
        "both_endpoints_critical": rows_both,
        "exactly_one_endpoint_required_critical": rows_one,
        "simple_simple": rows_simple,
        "unresolved_three_black_WW_systems": unresolved,
    }
    return archive


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def write_archive(archive: dict[str, Any], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(canonical_bytes(archive))


def verify_archive(archive: dict[str, Any]) -> None:
    if archive["marked_black_degree_records"] != 196:
        raise ArithmeticError("black-degree census changed")
    if archive["distinct_endpoint_critical_systems"] != 588:
        raise ArithmeticError("system count changed")
    if archive["endpoint_critical_marked_records"] != 980:
        raise ArithmeticError("passport count changed")
    if archive["eliminated_endpoint_critical_systems"] != 588:
        raise ArithmeticError(archive["unresolved_endpoint_critical_systems"])
    if archive["eliminated_endpoint_critical_marked_records"] != 980:
        raise ArithmeticError("record conversion changed")
    if archive["eliminated_simple_simple_systems"] != 196:
        raise ArithmeticError(archive["unresolved_three_black_WW_systems"])
    if archive["eliminated_simple_simple_marked_records"] != 392:
        raise ArithmeticError("simple/simple record conversion changed")
    if archive["eliminated_three_black_WW_marked_records"] != 1372:
        raise ArithmeticError("three-black/WW census did not close")
    if archive["eliminated_three_black_WW_coefficient_systems"] != 784:
        raise ArithmeticError("three-black/WW system census did not close")
    if archive["remaining_three_black_WW_marked_records"] != 0:
        raise ArithmeticError("three-black/WW remainder is nonempty")
    if archive["complete_degree_29_no_go_claimed"]:
        raise ArithmeticError("this bounded certificate is not a global no-go")


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
        if not arguments.output.exists():
            raise SystemExit(f"missing archive: {arguments.output}")
        if arguments.output.read_bytes() != encoded:
            raise SystemExit("archive is stale")
    else:
        write_archive(archive, arguments.output)
    print(
        json.dumps(
            {
                "status": archive["status"],
                "eliminated_systems": archive[
                    "eliminated_three_black_WW_coefficient_systems"
                ],
                "eliminated_marked_records": archive[
                    "eliminated_three_black_WW_marked_records"
                ],
                "remaining_three_black_WW_marked_records": archive[
                    "remaining_three_black_WW_marked_records"
                ],
                "endpoint_critical_prime_histogram": archive[
                    "endpoint_critical_prime_histogram"
                ],
                "simple_simple_prime_histogram": archive[
                    "simple_simple_prime_histogram"
                ],
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
