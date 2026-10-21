#!/usr/bin/env python3
"""Exact checks for the compact-H3 general-fork field funnel.

The mathematical proof is in
``results/compact_h3_general_fork_arrangement_field_finiteness.md``.
This certifier deliberately checks only its finite algebraic/arithmetic
claims; it is not a substitute for the geometric existence argument.
"""

from __future__ import annotations

import json
from math import gcd

import sympy as sp


LABEL_CUTOFF = 1_000_000


def euler_phi(n: int) -> int:
    result = n
    value = n
    p = 2
    while p * p <= value:
        if value % p == 0:
            while value % p == 0:
                value //= p
            result -= result // p
        p += 1
    if value > 1:
        result -= result // value
    return result


def three_node_determinant_identity() -> dict[str, object]:
    u, a, b, x = sp.symbols("u a b x")
    gram = sp.Matrix([[-2, u, a], [u, -2, b], [a, b, -2]])
    determinant = sp.expand(gram.det())
    beta = sp.expand(determinant / 2)
    discriminant = sp.factor(sp.discriminant(x**2 + a * b * x + a**2 + b**2 - 4, x))
    expected_beta = u**2 + a * b * u + a**2 + b**2 - 4
    expected_discriminant = (4 - a**2) * (4 - b**2)
    return {
        "determinant": str(determinant),
        "beta": str(beta),
        "determinant_passed": sp.expand(beta - expected_beta) == 0,
        "discriminant": str(discriminant),
        "discriminant_passed": sp.expand(discriminant - expected_discriminant) == 0,
        "cartan_cycle_formula": "beta=t12+t13+t23-c123-4",
    }


def totient_lower_bound_regression(limit: int = 100_000) -> bool:
    # phi(n)^2 >= n/2, checked here and proved multiplicatively in the report.
    return all(2 * euler_phi(n) ** 2 >= n for n in range(3, limit + 1))


def compositum_degree_lower_bound(m: int, n: int) -> int:
    """Elementary lower bound [Q(s_m,s_n):Q] >= max(phi(m),phi(n))/2."""
    return max(1, max(euler_phi(m), euler_phi(n)) // 2)


def capacity_ratio_bound_at_cutoff() -> dict[str, object]:
    # For M=max(m,n)>=10^6, D=max(phi(m),phi(n))/2 >=353 by
    # phi(M)^2>=M/2.  The report's beta norm gives
    #   1 < 4 M^4 (13/16)^D,
    # while the following exact integer inequality gives the contradiction.
    M = LABEL_CUTOFF
    D = 353
    lhs = 16**D
    rhs = 4 * M**4 * 13**D
    return {
        "M": M,
        "D_lower_bound": D,
        "comparison": "16^D > 4*M^4*13^D",
        "lhs": str(lhs),
        "rhs": str(rhs),
        "passed": lhs > rhs,
    }


def monotonicity_certificate() -> dict[str, object]:
    # f(x)=sqrt(x/2)/2*log(16/13)-4log x-log 4.
    # Since log(16/13)>1/5, f'(x)>0 once sqrt(x/2)>80.
    return {
        "log_bound": "log(16/13)>1/5",
        "derivative_positive_if_sqrt_x_over_2_gt": 80,
        "applied_from": LABEL_CUTOFF,
        "passed": LABEL_CUTOFF > 2 * 80**2,
    }


def field_degree_sanity(limit: int = 100) -> dict[str, object]:
    # The base E=Q(4 sin^2(pi/m),4 sin^2(pi/n)) contains each factor,
    # so its degree is at least max(phi(m),phi(n))/2.  This finite loop only
    # guards integer rounding.
    rows = []
    passed = True
    for m in range(3, limit + 1):
        for n in range(3, limit + 1):
            lower = compositum_degree_lower_bound(m, n)
            target = max(euler_phi(m), euler_phi(n)) / 2
            ok = lower >= target and lower >= 1
            passed = passed and ok
            if (m, n) in {(3, 3), (5, 7), (99, 100)}:
                rows.append({"m": m, "n": n, "lower": lower, "passed": ok})
    return {"checked_through": limit, "samples": rows, "passed": passed}


def build_certificate() -> dict[str, object]:
    determinant = three_node_determinant_identity()
    cutoff = capacity_ratio_bound_at_cutoff()
    monotonicity = monotonicity_certificate()
    degree_sanity = field_degree_sanity()
    totient = totient_lower_bound_regression()
    passed = all(
        [
            determinant["determinant_passed"],
            determinant["discriminant_passed"],
            cutoff["passed"],
            monotonicity["passed"],
            degree_sanity["passed"],
            totient,
        ]
    )
    return {
        "schema": "compact-h3-general-fork-field-finiteness-v1",
        "all_checks_pass": passed,
        "three_node_packet": determinant,
        "label_cutoff": cutoff,
        "monotonicity": monotonicity,
        "totient_bound": {
            "formula": "2*phi(d)^2 >= d",
            "checked_through": 100_000,
            "passed": totient,
        },
        "base_field_degree": degree_sanity,
        "field_primitivity_gate": (
            "beta has beta_0>0 and beta_sigma<0 at every nonselected place; "
            "any embedding fixing beta would send the positive selected value "
            "to itself and is therefore selected, hence Q(beta)=k"
        ),
    }


def main() -> None:
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))
    if not certificate["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
