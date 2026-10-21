#!/usr/bin/env python3
"""Exact finite checks for the numerical constants in the AGCat theorem.

The accompanying proof is
``results/coxeter_finiteness_constant_audit_and_h3_cutoff.md``.  This
program certifies only the finite integer sieves in that note.  In
particular, a surviving label pair is a survivor of a necessary norm
inequality, not a realizable Coxeter packet.
"""

from __future__ import annotations

import json
from math import gcd


OLD_FORK_CUTOFF = 1_000_000
FORK_INDIVIDUAL_SCAN_START = 841
FORK_INDIVIDUAL_BOUND = 840
FORK_PAIR_BOUND = 420
TAIL_OLD_CUTOFF = 20_000
TAIL_BOUND = 126


def arithmetic_sieves(limit: int) -> tuple[list[int], list[int]]:
    """Return Euler phi and the least prime factor through ``limit``."""
    phi = list(range(limit + 1))
    least = list(range(limit + 1))
    for prime in range(2, limit + 1):
        if phi[prime] != prime:
            continue
        for value in range(prime, limit + 1, prime):
            phi[value] -= phi[value] // prime
        if prime * prime <= limit:
            for value in range(prime * prime, limit + 1, prime):
                if least[value] == value:
                    least[value] = prime
    return phi, least


def gamma_from_least(n: int, least: list[int]) -> int:
    """Cyclotomic norm gamma(n): p for a prime power p^a, and 1 otherwise."""
    prime = least[n]
    value = n
    while value % prime == 0:
        value //= prime
    return prime if value == 1 else 1


def greater_than_power_of_two(value: int, exponent: int) -> bool:
    """Decide ``value > 2**exponent`` exactly, usually without making 2**e."""
    top = value.bit_length() - 1
    if top != exponent:
        return top > exponent
    return value & (value - 1) != 0


def fork_individual_cutoff_certificate() -> dict[str, object]:
    """Exclude 841 <= M < 10^6 by the one-label norm envelope.

    Put d=phi(M)/2.  The necessary inequality is

        4**d < 4*M**4*gamma(M)*3**d.

    Most values are certified without large powers by
    (4/3)**10 > 2**4.  The few misses are compared as exact integers.
    The already frozen H3 proof excludes M >= 10^6.
    """
    limit = OLD_FORK_CUTOFF - 1
    phi, least = arithmetic_sieves(limit)
    exact_fallbacks: list[dict[str, object]] = []
    failures: list[int] = []
    contraction_failures: list[int] = []
    for label in range(FORK_INDIVIDUAL_SCAN_START, limit + 1):
        half_phi = phi[label] // 2
        gamma = gamma_from_least(label, least)
        exceptional = 4 * label**4 * gamma

        # Check 3*h(label)/4 < 1 exactly.  After raising to phi(label),
        # this is gamma(label)^2*3^phi < 4^phi.  For gamma=1 it is
        # immediate; for a prime-power label the same exact binary lower
        # bound avoids constructing enormous powers.
        if gamma != 1:
            contraction_exponent = 4 * (phi[label] // 10)
            contraction_coarse = (gamma**2).bit_length() <= contraction_exponent
            if not contraction_coarse:
                contraction_exact = 4 ** phi[label] > gamma**2 * 3 ** phi[label]
                if not contraction_exact:
                    contraction_failures.append(label)

        # (4/3)^half_phi > 2^(4*floor(half_phi/10)).  If this
        # lower bound already dominates the exceptional factor, the desired
        # strict integer comparison follows.
        binary_exponent = 4 * (half_phi // 10)
        coarse_pass = exceptional.bit_length() <= binary_exponent
        if coarse_pass:
            continue

        exact_pass = 4**half_phi > exceptional * 3**half_phi
        exact_fallbacks.append(
            {
                "label": label,
                "phi": phi[label],
                "gamma": gamma,
                "passed": exact_pass,
            }
        )
        if not exact_pass:
            failures.append(label)

    return {
        "scan": [FORK_INDIVIDUAL_SCAN_START, limit],
        "necessary_inequality": "4^d < 4*M^4*gamma(M)*3^d, d=phi(M)/2",
        "excluded_count": limit - FORK_INDIVIDUAL_SCAN_START + 1,
        "exact_fallbacks": exact_fallbacks,
        "failures": failures,
        "contraction_failures": contraction_failures,
        "passed": not failures and not contraction_failures,
        "conclusion_with_old_analytic_tail": "M<=840",
    }


def real_cyclotomic_compositum_degree(m: int, n: int, phi: list[int]) -> int:
    """Degree of Q(s_m,s_n), s_r=4*sin(pi/r)^2, for m,n >= 3."""
    common = gcd(m, n)
    multiple = m // common * n
    # The simultaneous sign kernel in (Z/lcm)^* has size two when
    # gcd(m,n)>2 and size four when gcd(m,n)<=2.
    divisor = 2 if common > 2 else 4
    assert phi[multiple] % divisor == 0
    return phi[multiple] // divisor


def fork_pair_sieve_certificate() -> dict[str, object]:
    """Enumerate the exact necessary norm filter after M <= 840."""
    label_limit = FORK_INDIVIDUAL_BOUND
    lcm_limit = label_limit * (label_limit - 1)
    phi, least = arithmetic_sieves(lcm_limit)
    survivors: list[dict[str, int]] = []
    divisibility_failures: list[tuple[int, int]] = []

    for m in range(3, label_limit + 1):
        degree_m = phi[m] // 2
        gamma_m = gamma_from_least(m, least)
        for n in range(m, label_limit + 1):
            degree_n = phi[n] // 2
            gamma_n = gamma_from_least(n, least)
            degree = real_cyclotomic_compositum_degree(m, n, phi)
            if degree % degree_m or degree % degree_n:
                divisibility_failures.append((m, n))
                continue

            numerator = (
                4
                * m**2
                * n**2
                * gamma_m ** (degree // degree_m)
                * gamma_n ** (degree // degree_n)
            )
            # The necessary strict inequality is numerator > 4**degree.
            if greater_than_power_of_two(numerator, 2 * degree):
                survivors.append({"m": m, "n": n, "degree": degree})

    maximum = max(row["n"] for row in survivors)
    endpoint = [row for row in survivors if row["n"] == maximum]
    return {
        "input_label_box": [3, label_limit],
        "degree_formula": (
            "phi(lcm(m,n))/2 if gcd(m,n)>2, otherwise phi(lcm(m,n))/4"
        ),
        "necessary_inequality": (
            "4*m^2*n^2*gamma(m)^(D/d_m)*gamma(n)^(D/d_n) > 4^D"
        ),
        "divisibility_failures": divisibility_failures,
        "survivor_count": len(survivors),
        "maximum_surviving_label": maximum,
        "endpoint_survivors": endpoint,
        "survivors": survivors,
        "passed": (
            not divisibility_failures
            and len(survivors) == 525
            and maximum == FORK_PAIR_BOUND
            and endpoint == [{"m": 3, "n": 420, "degree": 48}]
        ),
    }


def tail_sieve_certificate() -> dict[str, object]:
    """Sharpen the tail-chart cutoff below the old analytic threshold."""
    phi, least = arithmetic_sieves(TAIL_OLD_CUTOFF - 1)
    survivors: list[int] = []
    for label in range(3, TAIL_OLD_CUTOFF):
        half_phi = phi[label] // 2
        gamma = gamma_from_least(label, least)
        numerator = 601 * label**4 * gamma**2
        # A label can survive only if numerator > 4**half_phi.
        if greater_than_power_of_two(numerator, 2 * half_phi):
            survivors.append(label)
    return {
        "scan": [3, TAIL_OLD_CUTOFF - 1],
        "necessary_inequality": "601*d^4*gamma(d)^2 > 4^(phi(d)/2)",
        "survivor_count": len(survivors),
        "maximum_surviving_label": max(survivors),
        "survivors": survivors,
        "passed": len(survivors) == 66 and max(survivors) == TAIL_BOUND,
        "conclusion_with_old_analytic_tail": "d<=126",
    }


def build_certificate() -> dict[str, object]:
    individual = fork_individual_cutoff_certificate()
    pairs = fork_pair_sieve_certificate()
    tail = tail_sieve_certificate()
    return {
        "schema": "coxeter-finiteness-constant-improvements-v1",
        "all_checks_pass": all(
            bool(row["passed"]) for row in (individual, pairs, tail)
        ),
        "fork_individual_cutoff": individual,
        "fork_pair_sieve": pairs,
        "tail_sieve": tail,
        "scope_warning": (
            "Every survivor passes a necessary norm inequality only; this "
            "certificate does not assert field or Coxeter-packet existence."
        ),
    }


def main() -> None:
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))
    if not certificate["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
