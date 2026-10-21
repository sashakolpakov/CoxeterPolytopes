#!/usr/bin/env python3
"""Bounded exact reduction of the degree-29 four-black stratum.

The full four-black marked-passport block has 42,994 records but only 12,752
distinct coefficient systems.  This certificate deliberately attacks two
uniform subloci, without enumerating plane trees:

* both normalised endpoints black (3,276 systems / 9,828 records);
* exactly one endpoint black and the other endpoint a non-root critical point
  (3,458 systems / 13,832 records).

Both blocks reduce to univariate resultants over Q(lambda).  Modular gcd one
at a degree-preserving split prime proves the characteristic-zero resultant
nonzero.  The global degree-29 problem remains open.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from certify_gram_belyi_m5_degree29_reduction import FpPoly, lambda_mod


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/gram_belyi_m5_degree29_four_black_reduction.json"
SCHEMA = "coxeter.gram_belyi.m5_degree29_four_black_reduction.v1"
SPLIT_PRIMES = (101, 109, 131, 139, 149, 151, 179, 181, 191, 199)


def four_black_degree_records() -> list[tuple[int, int, int, int]]:
    """lambda degree a, followed by an unordered positive triple b<=c<=d."""

    answer = []
    for a in range(1, 27):
        for b in range(1, 29 - a):
            for c in range(b, 29 - a - b):
                d = 29 - a - b - c
                if c <= d:
                    answer.append((a, b, c, d))
    if len(answer) != 640:
        raise ArithmeticError(len(answer))
    return answer


def one_black_endpoint_records() -> list[tuple[int, int, int, int]]:
    """lambda degree a, endpoint degree b, remaining roots c<=d."""

    answer = []
    for a in range(1, 27):
        for b in range(1, 28 - a):
            for c in range(1, (29 - a - b) // 2 + 1):
                d = 29 - a - b - c
                if c <= d:
                    answer.append((a, b, c, d))
    if len(answer) != 1729:
        raise ArithmeticError(len(answer))
    return answer


def both_black_endpoint_records() -> list[tuple[int, int, int, int]]:
    """lambda degree a, ordered endpoint degrees b,c, final degree d."""

    answer = []
    for a in range(1, 27):
        for b in range(1, 28 - a):
            for c in range(1, 29 - a - b):
                d = 29 - a - b - c
                if d > 0:
                    answer.append((a, b, c, d))
    if len(answer) != 3276:
        raise ArithmeticError(len(answer))
    return answer


def trim_bivariate(value: list[list[int]]) -> list[list[int]]:
    while len(value) > 1 and value[-1] == [0]:
        value.pop()
    return value


def bivariate_add(
    ring: FpPoly, left: list[list[int]], right: list[list[int]]
) -> list[list[int]]:
    answer = [
        ring.add(
            left[index] if index < len(left) else [0],
            right[index] if index < len(right) else [0],
        )
        for index in range(max(len(left), len(right)))
    ]
    return trim_bivariate(answer)


def bivariate_scale(
    ring: FpPoly, value: list[list[int]], scalar: int
) -> list[list[int]]:
    return trim_bivariate([ring.scale(coefficient, scalar) for coefficient in value])


def bivariate_mul(
    ring: FpPoly, left: list[list[int]], right: list[list[int]]
) -> list[list[int]]:
    answer = [[0] for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[i + j] = ring.add(answer[i + j], ring.mul(a, b))
    return trim_bivariate(answer)


def bivariate_power(
    ring: FpPoly, value: list[list[int]], exponent: int
) -> list[list[int]]:
    answer = [[1]]
    base = value
    while exponent:
        if exponent & 1:
            answer = bivariate_mul(ring, answer, base)
        base = bivariate_mul(ring, base, base)
        exponent //= 2
    return answer


def parameter_remainder(
    prime: int,
    roots: tuple[list[int], ...],
    weights: tuple[int, ...],
) -> list[list[int]]:
    """Remainder in x, with coefficients in F_p[t], modulo residual D."""

    ring = FpPoly(prime)
    factors = [(ring.neg(root), [1]) for root in roots]
    product_value = [[1]]
    for factor, exponent in zip(factors, weights):
        product_value = bivariate_mul(
            ring, product_value, bivariate_power(ring, list(factor), exponent)
        )

    derivative_residual = [[0]]
    for omitted, weight in enumerate(weights):
        term = [[1]]
        for index, factor in enumerate(factors):
            if index != omitted:
                term = bivariate_mul(ring, term, list(factor))
        derivative_residual = bivariate_add(
            ring, derivative_residual, bivariate_scale(ring, term, weight)
        )

    inverse_leading = pow(derivative_residual[-1][0], -1, prime)
    remainder = [coefficient[:] for coefficient in product_value]
    while len(remainder) >= len(derivative_residual):
        quotient_coefficient = ring.scale(remainder[-1], inverse_leading)
        shift = len(remainder) - len(derivative_residual)
        for index, coefficient in enumerate(derivative_residual):
            remainder[shift + index] = ring.sub(
                remainder[shift + index],
                ring.mul(quotient_coefficient, coefficient),
            )
        trim_bivariate(remainder)
    return remainder


def certify_both_black_endpoints() -> list[dict[str, Any]]:
    rows = []
    for a, b, c, d in both_black_endpoint_records():
        found = None
        for prime in SPLIT_PRIMES:
            sqrt5, lam = lambda_mod(prime)
            ring = FpPoly(prime)
            remainder = parameter_remainder(
                prime,
                ([lam], [0], [1], [0, 1]),
                (a, b, c, d),
            )
            if len(remainder) != 3:
                continue
            r1, r2 = remainder[1], remainder[2]
            # R_2 is homogeneous of degree 27.  R_1 is homogeneous of
            # degree 28, but its pure t^28 term vanishes when the other
            # three roots coalesce at zero, so both t-degrees are <=27.
            if len(r1) - 1 != 27 or len(r2) - 1 != 27:
                continue
            if ring.gcd(r1, r2) == [1]:
                found = {
                    "a": a,
                    "b_x0": b,
                    "c_x1": c,
                    "d_moving": d,
                    "prime": prime,
                    "sqrt5_mod_prime": sqrt5,
                    "lambda_mod_prime": lam,
                    "remainder_degrees_R1_R2": [27, 27],
                    "gcd": [1],
                }
                break
        if found is None:
            raise ArithmeticError(("both endpoints black", a, b, c, d))
        rows.append(found)
    return rows


def scalar_product(values: Iterable[int], prime: int) -> int:
    answer = 1
    for value in values:
        answer = answer * value % prime
    return answer


def residual_value_at(
    endpoint: int,
    roots: tuple[int, int, int, int],
    weights: tuple[int, int, int, int],
    prime: int,
) -> int:
    answer = 0
    for omitted, weight in enumerate(weights):
        answer += weight * scalar_product(
            (
                (endpoint - roots[index]) % prime
                for index in range(4)
                if index != omitted
            ),
            prime,
        )
    return answer % prime


def scalar_remainder(
    prime: int,
    roots: tuple[int, int, int, int],
    weights: tuple[int, int, int, int],
) -> list[int]:
    ring = FpPoly(prime)
    factors = [[-root % prime, 1] for root in roots]
    product_value = [1]
    for factor, exponent in zip(factors, weights):
        product_value = ring.mul(product_value, ring.power(factor, exponent))
    derivative_residual = [0]
    for omitted, weight in enumerate(weights):
        term = [1]
        for index, factor in enumerate(factors):
            if index != omitted:
                term = ring.mul(term, factor)
        derivative_residual = ring.add(
            derivative_residual, ring.scale(term, weight)
        )
    _quotient, remainder = ring.divmod(product_value, derivative_residual)
    return remainder + [0] * (3 - len(remainder))


def interpolate(prime: int, points: list[tuple[int, int]]) -> list[int]:
    """Newton interpolation over F_p, returning ascending coefficients."""

    ring = FpPoly(prime)
    xs = [point[0] for point in points]
    divided = [point[1] for point in points]
    for width in range(1, len(points)):
        for index in range(len(points) - 1, width - 1, -1):
            divided[index] = (
                (divided[index] - divided[index - 1])
                * pow((xs[index] - xs[index - width]) % prime, -1, prime)
            ) % prime
    answer = [0]
    basis = [1]
    for index, coefficient in enumerate(divided):
        answer = ring.add(answer, ring.scale(basis, coefficient))
        basis = ring.mul(basis, [-xs[index] % prime, 1])
    return answer


def critical_endpoint_equations_by_interpolation(
    prime: int,
    alpha: int,
    weights: tuple[int, int, int, int],
) -> tuple[list[int], list[int], int]:
    """Eliminate s from D(1)=0 for roots alpha,0,r,s.

    If D(1)=C(r)s+K(r), the returned polynomials are

      P_1=C^28 R_1(r,-K/C),  P_2=C^27 R_2(r,-K/C).

    Homogeneity gives deg P_1<=55 and deg P_2<=54.  Exact interpolation uses
    56 ordinary values and is verified at every other r with C(r)!=0.
    """

    ring = FpPoly(prime)

    def c_and_k(r_value: int) -> tuple[int, int]:
        roots0 = (alpha, 0, r_value, 0)
        roots1 = (alpha, 0, r_value, 1)
        k_value = residual_value_at(1, roots0, weights, prime)
        c_value = (
            residual_value_at(1, roots1, weights, prime) - k_value
        ) % prime
        return c_value, k_value

    c0, k0 = c_and_k(0)
    c_at_1, k_at_1 = c_and_k(1)
    c1 = (c_at_1 - c0) % prime
    k1 = (k_at_1 - k0) % prime
    exceptional_resultant = (c1 * k0 - c0 * k1) % prime

    values = []
    for r_value in range(prime):
        c_value, k_value = c_and_k(r_value)
        if c_value == 0:
            continue
        s_value = -k_value * pow(c_value, -1, prime) % prime
        remainder = scalar_remainder(
            prime, (alpha, 0, r_value, s_value), weights
        )
        values.append(
            (
                r_value,
                pow(c_value, 28, prime) * remainder[1] % prime,
                pow(c_value, 27, prime) * remainder[2] % prime,
            )
        )

    if len(values) < 56:
        raise ArithmeticError((prime, weights, len(values)))
    p1 = interpolate(prime, [(r, first) for r, first, _second in values[:56]])
    p2 = interpolate(prime, [(r, second) for r, _first, second in values[:56]])
    if any(
        ring.evaluate(p1, r) != first or ring.evaluate(p2, r) != second
        for r, first, second in values
    ):
        raise ArithmeticError(("interpolation bound failed", prime, weights))
    return p1, p2, exceptional_resultant


def certify_one_black_one_critical() -> list[dict[str, Any]]:
    rows = []
    for orientation in ("BW", "WB"):
        for a, b, c, d in one_black_endpoint_records():
            found = None
            for prime in SPLIT_PRIMES:
                sqrt5, lam = lambda_mod(prime)
                alpha = lam if orientation == "BW" else (1 - lam) % prime
                ring = FpPoly(prime)
                p1, p2, exceptional = critical_endpoint_equations_by_interpolation(
                    prime, alpha, (a, b, c, d)
                )
                if exceptional == 0:
                    continue
                if len(p1) - 1 != 55 or len(p2) - 1 != 54:
                    continue
                if ring.gcd(p1, p2) == [1]:
                    found = {
                        "orientation": orientation,
                        "a_lambda": a,
                        "b_black_endpoint": b,
                        "c": c,
                        "d": d,
                        "prime": prime,
                        "sqrt5_mod_prime": sqrt5,
                        "marked_root_mod_prime": alpha,
                        "cleared_remainder_degrees_R1_R2": [55, 54],
                        "exceptional_branch_resultant_mod_prime": exceptional,
                        "gcd": [1],
                    }
                    break
            if found is None:
                raise ArithmeticError((orientation, a, b, c, d))
            rows.append(found)
    return rows


def histogram(rows: list[dict[str, Any]]) -> dict[str, int]:
    return dict(
        (str(prime), count)
        for prime, count in sorted(Counter(row["prime"] for row in rows).items())
    )


def build_archive() -> dict[str, Any]:
    black_records = four_black_degree_records()
    one_black = one_black_endpoint_records()
    both_black = both_black_endpoint_records()
    rows_bb = certify_both_black_endpoints()
    rows_bc = certify_one_black_one_critical()

    return {
        "schema": SCHEMA,
        "status": "EXACT_QUANTIFIED_FOUR_BLACK_REDUCTION",
        "global_degree_29_status": "OPEN",
        "four_black_marked_passports_before": 42994,
        "four_black_coefficient_systems_before": 12752,
        "passport_system_quotient": {
            "marked_black_degree_records": len(black_records),
            "one_black_endpoint_degree_records": len(one_black),
            "both_black_endpoint_degree_records": len(both_black),
            "BB": {"marked_passports": 9828, "coefficient_systems": 3276},
            "BW": {"marked_passports": 12103, "coefficient_systems": 3458},
            "WB": {"marked_passports": 12103, "coefficient_systems": 3458},
            "WW": {"marked_passports": 8960, "coefficient_systems": 2560},
            "white_partitions": ["(4,1^25)", "(3,2,1^24)", "(2,2,2,1^23)"],
            "white_endpoint_status_multiplicities": {
                "one_endpoint": {"simple": 3, "critical": 4},
                "two_ordered_endpoints": {
                    "SS": 3,
                    "SC": 4,
                    "CS": 4,
                    "CC": 3,
                },
            },
            "at_least_one_white_endpoint_critical": {
                "coefficient_systems": 5378,
                "marked_passports": 20872,
            },
        },
        "eliminated_blocks": {
            "both_endpoints_black": {
                "coefficient_system_count": len(rows_bb),
                "marked_passport_count": 9828,
                "prime_histogram": histogram(rows_bb),
                "rows": rows_bb,
            },
            "one_black_other_white_critical": {
                "coefficient_system_count": len(rows_bc),
                "marked_passport_count": 13832,
                "prime_histogram": histogram(rows_bc),
                "rows": rows_bc,
            },
        },
        "eliminated_four_black_coefficient_systems": len(rows_bb) + len(rows_bc),
        "eliminated_four_black_marked_passports": 9828 + 13832,
        "remaining_four_black_coefficient_systems": 12752 - len(rows_bb) - len(rows_bc),
        "remaining_four_black_marked_passports": 42994 - 9828 - 13832,
        "global_marked_passports_eliminated_through_four_black": 29009,
        "global_marked_passports_remaining": 42581669,
        "remaining_blocks": {
            "one_black_other_white_simple": {
                "coefficient_systems": 3458,
                "marked_passports": 10374,
            },
            "both_white_at_least_one_critical": {
                "coefficient_systems": 1920,
                "marked_passports": 7040,
            },
            "both_white_simple": {
                "coefficient_systems": 640,
                "marked_passports": 1920,
            },
        },
        "degree_29_witness_found": False,
        "complete_degree_29_no_go_claimed": False,
        "complete_four_black_no_go_claimed": False,
        "broad_plane_tree_search_used": False,
    }


def verify_archive(archive: dict[str, Any]) -> None:
    if archive["four_black_marked_passports_before"] != 42994:
        raise ArithmeticError("four-black passport count changed")
    if archive["four_black_coefficient_systems_before"] != 12752:
        raise ArithmeticError("four-black system count changed")
    if archive["eliminated_four_black_coefficient_systems"] != 6734:
        raise ArithmeticError("eliminated system count changed")
    if archive["eliminated_four_black_marked_passports"] != 23660:
        raise ArithmeticError("eliminated passport count changed")
    if archive["remaining_four_black_coefficient_systems"] != 6018:
        raise ArithmeticError("remaining system count changed")
    if archive["remaining_four_black_marked_passports"] != 19334:
        raise ArithmeticError("remaining passport count changed")
    if archive["global_marked_passports_eliminated_through_four_black"] != 29009:
        raise ArithmeticError("global eliminated count changed")
    if archive["global_marked_passports_remaining"] != 42581669:
        raise ArithmeticError("global remainder changed")
    if archive["complete_four_black_no_go_claimed"]:
        raise ArithmeticError("four-black stratum is not closed")


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
            raise SystemExit("four-black archive is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_bytes(encoded)
    print(
        json.dumps(
            {
                "status": archive["status"],
                "eliminated_systems": archive[
                    "eliminated_four_black_coefficient_systems"
                ],
                "eliminated_marked_passports": archive[
                    "eliminated_four_black_marked_passports"
                ],
                "remaining_systems": archive[
                    "remaining_four_black_coefficient_systems"
                ],
                "remaining_marked_passports": archive[
                    "remaining_four_black_marked_passports"
                ],
                "sha256": digest,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
