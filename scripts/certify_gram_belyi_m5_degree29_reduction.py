#!/usr/bin/env python3
"""Exact finite reduction for sharpness of the m=5 Belyi lower bound.

Let lambda=(7+2*sqrt(5))/29.  Rodriguez gives h_Bel(lambda)>=29.
This certificate asks whether equality can hold.  After replacing B by 1-B,
we may assume B(lambda)=0 and call that fibre black.

The script exhausts all ramification passports in which the black fibre has
at most two vertices, and all three-black-vertex passports in which at least
one of the normalized source points 0,1 is black.  It also exhausts the
Galois-symmetric factor ansatz

    B = scalar*(29*x^2-14*x+1)^a*x^b*(x-1)^c,
    2*a+b+c=29.

Nonexistence is certified by exact polynomial identities reduced at good
split primes.  For one-parameter families, both equations retain their
proved characteristic-zero degree and have gcd one modulo the displayed
prime; equivalently their resultant is nonzero.  A true characteristic-zero
common factor would make every such good reduction resultant zero.

The remaining passports are explicitly counted.  A clean run is therefore a
rigorous finite reduction, NOT a complete degree-29 no-go.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/gram_belyi_m5_degree29_reduction.json"
SCHEMA = "coxeter.gram_belyi.m5_degree29_reduction.v1"
DEGREE = 29
SPLIT_PRIMES = (101, 109, 131, 139, 149, 151)


class FpPoly:
    """Small univariate polynomial engine over F_p, ascending coefficients."""

    def __init__(self, prime: int):
        self.p = prime

    def trim(self, values: Iterable[int]) -> list[int]:
        answer = [value % self.p for value in values]
        if not answer:
            answer = [0]
        while len(answer) > 1 and answer[-1] == 0:
            answer.pop()
        return answer

    def add(self, left: list[int], right: list[int]) -> list[int]:
        return self.trim(
            (left[i] if i < len(left) else 0)
            + (right[i] if i < len(right) else 0)
            for i in range(max(len(left), len(right)))
        )

    def neg(self, value: list[int]) -> list[int]:
        return self.trim(-coefficient for coefficient in value)

    def sub(self, left: list[int], right: list[int]) -> list[int]:
        return self.add(left, self.neg(right))

    def scale(self, value: list[int], scalar: int) -> list[int]:
        return self.trim(scalar * coefficient for coefficient in value)

    def mul(self, left: list[int], right: list[int]) -> list[int]:
        answer = [0] * (len(left) + len(right) - 1)
        for i, a in enumerate(left):
            for j, b in enumerate(right):
                answer[i + j] = (answer[i + j] + a * b) % self.p
        return self.trim(answer)

    def power(self, value: list[int], exponent: int) -> list[int]:
        answer = [1]
        base = value
        while exponent:
            if exponent & 1:
                answer = self.mul(answer, base)
            base = self.mul(base, base)
            exponent //= 2
        return answer

    def divmod(self, dividend: list[int], divisor: list[int]) -> tuple[list[int], list[int]]:
        left = self.trim(dividend)
        right = self.trim(divisor)
        if right == [0]:
            raise ZeroDivisionError
        if len(left) < len(right):
            return [0], left
        quotient = [0] * (len(left) - len(right) + 1)
        inverse = pow(right[-1], -1, self.p)
        while left != [0] and len(left) >= len(right):
            shift = len(left) - len(right)
            coefficient = left[-1] * inverse % self.p
            quotient[shift] = coefficient
            for index, value in enumerate(right):
                left[shift + index] = (
                    left[shift + index] - coefficient * value
                ) % self.p
            left = self.trim(left)
        return self.trim(quotient), left

    def exact_div(self, dividend: list[int], divisor: list[int]) -> list[int]:
        quotient, remainder = self.divmod(dividend, divisor)
        if remainder != [0]:
            raise ArithmeticError((dividend, divisor, remainder))
        return quotient

    def gcd(self, left: list[int], right: list[int]) -> list[int]:
        a, b = self.trim(left), self.trim(right)
        while b != [0]:
            _quotient, remainder = self.divmod(a, b)
            a, b = b, remainder
        if a == [0]:
            return a
        return self.scale(a, pow(a[-1], -1, self.p))

    def derivative(self, value: list[int]) -> list[int]:
        return self.trim(index * value[index] for index in range(1, len(value)))

    def evaluate(self, value: list[int], point: int) -> int:
        answer = 0
        for coefficient in reversed(value):
            answer = (answer * point + coefficient) % self.p
        return answer


def sqrt5_mod(prime: int) -> int:
    roots = [value for value in range(prime) if value * value % prime == 5]
    if len(roots) != 2:
        raise ArithmeticError((prime, roots))
    return roots[0]


def lambda_mod(prime: int) -> tuple[int, int]:
    sqrt5 = sqrt5_mod(prime)
    lam = (7 + 2 * sqrt5) * pow(29, -1, prime) % prime
    if (29 * lam * lam - 14 * lam + 1) % prime:
        raise ArithmeticError((prime, sqrt5, lam))
    return sqrt5, lam


def xr_scalar_mul(
    left: tuple[int, int],
    right: tuple[int, int],
    u: int,
    v: int,
    prime: int,
) -> tuple[int, int]:
    ac = left[0] * right[0]
    bd = left[1] * right[1]
    return (ac + bd * v) % prime, (
        left[0] * right[1] + left[1] * right[0] + bd * u
    ) % prime


def xr_scalar_power(
    base: tuple[int, int], exponent: int, u: int, v: int, prime: int
) -> tuple[int, int]:
    answer = (1, 0)
    while exponent:
        if exponent & 1:
            answer = xr_scalar_mul(answer, base, u, v, prime)
        base = xr_scalar_mul(base, base, u, v, prime)
        exponent //= 2
    return answer


def xr_scalar_product(
    factors: Iterable[tuple[tuple[int, int], int]],
    u: int,
    v: int,
    prime: int,
) -> tuple[int, int]:
    answer = (1, 0)
    for base, exponent in factors:
        answer = xr_scalar_mul(
            answer, xr_scalar_power(base, exponent, u, v, prime), u, v, prime
        )
    return answer


def xr_poly_mul(
    left: tuple[list[int], list[int]],
    right: tuple[list[int], list[int]],
    u: list[int],
    v: list[int],
    ring: FpPoly,
) -> tuple[list[int], list[int]]:
    ac = ring.mul(left[0], right[0])
    bd = ring.mul(left[1], right[1])
    return ring.add(ac, ring.mul(bd, v)), ring.add(
        ring.add(ring.mul(left[0], right[1]), ring.mul(left[1], right[0])),
        ring.mul(bd, u),
    )


def xr_poly_power(
    base: tuple[list[int], list[int]],
    exponent: int,
    u: list[int],
    v: list[int],
    ring: FpPoly,
) -> tuple[list[int], list[int]]:
    answer = ([1], [0])
    while exponent:
        if exponent & 1:
            answer = xr_poly_mul(answer, base, u, v, ring)
        base = xr_poly_mul(base, base, u, v, ring)
        exponent //= 2
    return answer


def xr_poly_product(
    factors: Iterable[tuple[tuple[list[int], list[int]], int]],
    u: list[int],
    v: list[int],
    ring: FpPoly,
) -> tuple[list[int], list[int]]:
    answer = ([1], [0])
    for base, exponent in factors:
        answer = xr_poly_mul(
            answer, xr_poly_power(base, exponent, u, v, ring), u, v, ring
        )
    return answer


def find_scalar_obstruction(builder) -> dict[str, int]:
    for prime in SPLIT_PRIMES:
        sqrt5, lam = lambda_mod(prime)
        value = builder(prime, lam)
        if value % prime:
            return {
                "prime": prime,
                "sqrt5_mod_prime": sqrt5,
                "lambda_mod_prime": lam,
                "nonzero_obstruction": value % prime,
            }
    raise ArithmeticError("the prime bank did not separate a scalar obstruction")


def certify_two_black_endpoint_black() -> list[dict[str, Any]]:
    rows = []
    for side in ("x=0 black", "x=1 black"):
        for a in range(1, 29):
            b = 29 - a

            def obstruction(prime: int, lam: int) -> int:
                marked = lam if side == "x=0 black" else (1 - lam) % prime
                z = b * marked * pow(29, -1, prime) % prime
                critical = pow(z - marked, a, prime) * pow(z, b, prime)
                endpoint = pow(1 - marked, a, prime)
                return endpoint - critical

            rows.append(
                {
                    "side": side,
                    "marked_root_multiplicity": a,
                    "other_root_multiplicity": b,
                    **find_scalar_obstruction(obstruction),
                }
            )
    return rows


def two_black_white_equations(
    prime: int, lam: int, a: int
) -> tuple[list[int], list[int]]:
    ring = FpPoly(prime)
    b = 29 - a
    t = [0, 1]
    f0 = ring.scale(ring.power(t, b), -pow(lam, a, prime))
    f1 = ring.scale(ring.power([1, -1], b), pow(1 - lam, a, prime))
    first = ring.sub(f0, f1)
    # z-lambda=a(t-lambda)/29 and z-t=-b(t-lambda)/29.
    coefficient = (
        pow(a, a, prime)
        * pow(-b, b, prime)
        * pow(pow(29, 29, prime), -1, prime)
    ) % prime
    fz = ring.scale(ring.power([-lam, 1], 29), coefficient)
    second = ring.sub(f0, fz)
    return first, second


def certify_two_black_both_white() -> list[dict[str, Any]]:
    rows = []
    for a in range(1, 29):
        b = 29 - a
        found = None
        for prime in SPLIT_PRIMES:
            sqrt5, lam = lambda_mod(prime)
            ring = FpPoly(prime)
            first, second = two_black_white_equations(prime, lam, a)
            # In characteristic zero the degrees are b and 29.  The first
            # leading coefficient cannot vanish: it would make the positive
            # real number lambda/(1-lambda) a root of unity.
            if len(first) - 1 != b or len(second) - 1 != 29:
                continue
            gcd = ring.gcd(first, second)
            if gcd == [1]:
                found = {
                    "marked_root_multiplicity": a,
                    "other_root_multiplicity": b,
                    "prime": prime,
                    "sqrt5_mod_prime": sqrt5,
                    "lambda_mod_prime": lam,
                    "equation_degrees": [b, 29],
                    "gcd": [1],
                }
                break
        if found is None:
            raise ArithmeticError(("two black, both endpoints white", a, b))
        rows.append(found)
    return rows


def three_black_one_endpoint_equations(
    prime: int, marked: int, a: int, b: int, c: int
) -> tuple[list[int], list[int]]:
    ring = FpPoly(prime)
    t = [0, 1]
    inverse29 = pow(29, -1, prime)
    # Roots marked^a, 0^b, t^c.  The non-root derivative factor is
    # 29*x^2-((a+b)t+(b+c)marked)x+b*marked*t.
    u = ring.scale(
        ring.add(ring.scale(t, a + b), [(b + c) * marked]), inverse29
    )
    v = ring.scale(t, -b * marked * inverse29)
    remainder = xr_poly_product(
        (
            (([-marked], [1]), a),
            (([0], [1]), b),
            ((ring.neg(t), [1]), c),
        ),
        u,
        v,
        ring,
    )
    endpoint = ring.scale(ring.power([1, -1], c), pow(1 - marked, a, prime))
    return remainder[1], ring.sub(endpoint, remainder[0])


def certify_three_black_one_endpoint() -> list[dict[str, Any]]:
    rows = []
    for side in ("x=0 black", "x=1 black"):
        for a in range(1, 28):
            for b in range(1, 29 - a):
                c = 29 - a - b
                found = None
                for prime in SPLIT_PRIMES:
                    sqrt5, lam = lambda_mod(prime)
                    marked = lam if side == "x=0 black" else (1 - lam) % prime
                    ring = FpPoly(prime)
                    first, second = three_black_one_endpoint_equations(
                        prime, marked, a, b, c
                    )
                    # Homogeneous division by the quadratic derivative factor
                    # gives degree at most 28 in t for both equations; the
                    # displayed reduction proves both leading coefficients
                    # are nonzero in characteristic zero.
                    if len(first) - 1 != 28 or len(second) - 1 != 28:
                        continue
                    if ring.gcd(first, second) == [1]:
                        found = {
                            "side": side,
                            "multiplicities_marked_endpoint_free": [a, b, c],
                            "prime": prime,
                            "sqrt5_mod_prime": sqrt5,
                            "lambda_mod_prime": lam,
                            "equation_degrees": [28, 28],
                            "gcd": [1],
                        }
                        break
                if found is None:
                    raise ArithmeticError(("three black, one endpoint", side, a, b, c))
                rows.append(found)
    return rows


def certify_three_black_both_endpoints() -> list[dict[str, Any]]:
    rows = []
    for a in range(1, 28):
        for b in range(1, 29 - a):
            c = 29 - a - b

            def obstruction(prime: int, lam: int) -> int:
                inverse29 = pow(29, -1, prime)
                u = (a + b + (b + c) * lam) * inverse29 % prime
                v = -b * lam * inverse29 % prime
                remainder = xr_scalar_product(
                    (
                        (((-lam) % prime, 1), a),
                        (((0, 1)), b),
                        (((-1) % prime, 1), c),
                    ),
                    u,
                    v,
                    prime,
                )
                return remainder[1]

            rows.append(
                {
                    "multiplicities_lambda_0_1": [a, b, c],
                    **find_scalar_obstruction(obstruction),
                }
            )
    return rows


def certify_galois_symmetric_ansatz() -> list[dict[str, Any]]:
    rows = []
    for a in range(1, 15):
        for b in range(0, 30 - 2 * a):
            c = 29 - 2 * a - b
            found = None
            for prime in SPLIT_PRIMES:
                ring = FpPoly(prime)
                q = [1, -14, 29]
                x = [0, 1]
                xm1 = [-1, 1]
                raw = ring.mul(
                    ring.power(q, a),
                    ring.mul(ring.power(x, b), ring.power(xm1, c)),
                )
                forced = ring.power(q, a - 1)
                if b:
                    forced = ring.mul(forced, ring.power(x, b - 1))
                if c:
                    forced = ring.mul(forced, ring.power(xm1, c - 1))
                derivative = ring.derivative(raw)
                residual = ring.exact_div(derivative, forced)
                _quotient, remainder = ring.divmod(raw, residual)
                obstruction = None
                for index in range(1, len(remainder)):
                    if remainder[index]:
                        obstruction = {
                            "kind": "critical_values_not_equal",
                            "remainder_coefficient_index": index,
                            "nonzero_obstruction": remainder[index],
                        }
                        break
                if obstruction is None and b == 0:
                    value = (ring.evaluate(raw, 0) - remainder[0]) % prime
                    if value:
                        obstruction = {
                            "kind": "x=0_not_in_common_white_fibre",
                            "nonzero_obstruction": value,
                        }
                if obstruction is None and c == 0:
                    value = (ring.evaluate(raw, 1) - remainder[0]) % prime
                    if value:
                        obstruction = {
                            "kind": "x=1_not_in_common_white_fibre",
                            "nonzero_obstruction": value,
                        }
                if obstruction is not None:
                    found = {
                        "exponents_q_x_xminus1": [a, b, c],
                        "residual_derivative_degree": len(residual) - 1,
                        "prime": prime,
                        **obstruction,
                    }
                    break
            if found is None:
                raise ArithmeticError(("Galois-symmetric survivor", a, b, c))
            rows.append(found)
    return rows


def partitions(total: int, maximum: int | None = None):
    if total == 0:
        yield ()
        return
    largest = min(total, maximum or total)
    for first in range(largest, 0, -1):
        for tail in partitions(total - first, first):
            yield (first,) + tail


def choice_counts(multiplicities: Counter[int]) -> tuple[int, int]:
    distinct = sum(value > 0 for value in multiplicities.values())
    repeatable = sum(value >= 2 for value in multiplicities.values())
    return distinct, distinct * (distinct - 1) + repeatable


def passport_census() -> dict[str, Any]:
    by_length: dict[int, list[tuple[int, ...]]] = defaultdict(list)
    for partition in partitions(29):
        by_length[len(partition)].append(partition)

    ordinary_by_black_vertices = {}
    marked_by_black_vertices = {}
    marked_total = 0
    ordinary_total = 0
    for black_vertices in range(1, 30):
        black = by_length[black_vertices]
        white = by_length[30 - black_vertices]
        ordinary = len(black) * len(white)
        ordinary_by_black_vertices[str(black_vertices)] = ordinary
        ordinary_total += ordinary

        marked_black_count = marked_black_one = marked_black_two = 0
        for partition in black:
            multiplicities = Counter(partition)
            for lambda_degree in multiplicities:
                available = multiplicities.copy()
                available[lambda_degree] -= 1
                one, two = choice_counts(available)
                marked_black_count += 1
                marked_black_one += one
                marked_black_two += two

        white_one = white_two = 0
        for partition in white:
            one, two = choice_counts(Counter(partition))
            white_one += one
            white_two += two

        color_rows = {
            "BB": marked_black_two * len(white),
            "BW": marked_black_one * white_one,
            "WB": marked_black_one * white_one,
            "WW": marked_black_count * white_two,
        }
        subtotal = sum(color_rows.values())
        marked_by_black_vertices[str(black_vertices)] = {
            **color_rows,
            "total": subtotal,
        }
        marked_total += subtotal

    if ordinary_total != 235_465 or marked_total != 42_610_678:
        raise ArithmeticError((ordinary_total, marked_total))
    return {
        "ordinary_passport_pair_count": ordinary_total,
        "ordinary_by_lambda_color_vertex_count": ordinary_by_black_vertices,
        "marked_ramification_passport_count": marked_total,
        "marked_by_lambda_color_vertex_count_and_endpoint_colors": marked_by_black_vertices,
        "marking_convention": (
            "B(lambda)=0; B/W record the colours of source points x=0,x=1; "
            "cycle identities with equal valency are quotient out, while endpoint order is retained"
        ),
    }


def row_digest(rows: list[dict[str, Any]]) -> str:
    encoded = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def summarized(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "trial_count": len(rows),
        "witness_prime_histogram": dict(
            (str(prime), count)
            for prime, count in sorted(Counter(row["prime"] for row in rows).items())
        ),
        "canonical_rows_sha256": row_digest(rows),
        "rows": rows,
    }


def build_certificate() -> dict[str, Any]:
    census = passport_census()
    two_endpoint = certify_two_black_endpoint_black()
    two_white = certify_two_black_both_white()
    three_one = certify_three_black_one_endpoint()
    three_both = certify_three_black_both_endpoints()
    symmetric = certify_galois_symmetric_ansatz()

    low = census["marked_by_lambda_color_vertex_count_and_endpoint_colors"]
    if low["1"] != {"BB": 0, "BW": 0, "WB": 0, "WW": 1, "total": 1}:
        raise ArithmeticError(low["1"])
    if low["2"] != {"BB": 0, "BW": 56, "WB": 56, "WW": 84, "total": 196}:
        raise ArithmeticError(low["2"])
    if low["3"] != {
        "BB": 756,
        "BW": 1512,
        "WB": 1512,
        "WW": 1372,
        "total": 5152,
    }:
        raise ArithmeticError(low["3"])

    eliminated_marked = 1 + 196 + 756 + 1512 + 1512
    if eliminated_marked != 3_977:
        raise ArithmeticError(eliminated_marked)
    remaining = census["marked_ramification_passport_count"] - eliminated_marked
    if remaining != 42_606_701:
        raise ArithmeticError(remaining)

    return {
        "schema": SCHEMA,
        "all_checks_pass": True,
        "problem": {
            "coordinate": "lambda=(7+2*sqrt(5))/29",
            "minimal_polynomial": "29*x^2-14*x+1",
            "Rodriguez_lower_bound": 29,
            "question": "does a normalized degree-29 Belyi polynomial map lambda to 0 or 1?",
        },
        "passport_census": census,
        "exact_reductions": {
            "one_black_vertex": {
                "marked_passports_eliminated": 1,
                "proof": (
                    "B=C*(x-lambda)^29 would force the negative real ratio "
                    "(1-lambda)/(-lambda) to have 29th power one"
                ),
            },
            "two_black_one_endpoint_black": summarized(two_endpoint),
            "two_black_both_endpoints_white": summarized(two_white),
            "three_black_both_endpoints_black": summarized(three_both),
            "three_black_exactly_one_endpoint_black": summarized(three_one),
            "Galois_symmetric_q_power_ansatz": summarized(symmetric),
        },
        "theorem": {
            "degree_29_global_status": "OPEN",
            "proved": (
                "every degree-29 witness has at least three vertices in the colour "
                "containing lambda; if it has exactly three, both x=0 and x=1 "
                "have the opposite colour"
            ),
            "marked_ramification_passports_eliminated": eliminated_marked,
            "marked_ramification_passports_remaining": remaining,
            "three_black_remaining_endpoint_colour": "WW",
            "three_black_remaining_marked_passport_count": 1372,
            "four_or_more_black_remaining_marked_passport_count": 42_605_329,
            "Galois_symmetric_ansatz_trials_eliminated": len(symmetric),
        },
        "claim_boundary": {
            "degree_29_witness_found": False,
            "complete_degree_29_no_go_claimed": False,
            "three_black_WW_coefficient_systems_solved": False,
            "four_or_more_black_strata_solved": False,
            "ordinary_passports_are_not_dessins": True,
            "modular_timeout_used_as_evidence": False,
        },
    }


def canonical(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    rendered = canonical(payload)
    if args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("frozen m=5 degree-29 reduction differs from exact replay")
    print(rendered, end="")


if __name__ == "__main__":
    main()
