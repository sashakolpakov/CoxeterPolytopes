#!/usr/bin/env python3
"""Independent exact audit for the compact-H3 general-fork field funnel.

This checks the finite algebra used in the hostile audit.  The accompanying
report, not this script, proves the geometric packet-existence and Fekete
steps.
"""

from __future__ import annotations

import json
from fractions import Fraction

import sympy as sp


LABEL_CUTOFF = 1_000_000
ENDPOINT_EXPONENT = 353


def determinant_and_cartan_sign() -> dict[str, object]:
    u, a, b = sp.symbols("u a b", real=True)
    gram = sp.Matrix([[-2, u, a], [u, -2, b], [a, b, -2]])
    beta = sp.expand(gram.det() / 2)

    # With c_ij=2(r_i,r_j)/(r_i,r_i) and square -2 roots,
    # c_ij=-G_ij.  Hence z=c12*c23*c31=-uab.
    t12, t13, t23, z = u**2, a**2, b**2, -u * a * b
    invariant_beta = sp.expand(t12 + t13 + t23 - z - 4)
    wrong_plus_sign = sp.expand(t12 + t13 + t23 + z - 4)

    polynomial = sp.Poly(beta, u)
    discriminant = sp.factor(sp.discriminant(polynomial.as_expr(), u))
    minimum = sp.factor(polynomial.eval(-a * b / 2))
    equality_example = gram.subs({a: 1, b: 1, u: sp.Rational(-1, 2)})
    equality_leading_minors = [
        sp.factor(equality_example[:size, :size].det()) for size in (1, 2, 3)
    ]

    return {
        "determinant": str(sp.expand(gram.det())),
        "beta": str(beta),
        "cycle_formula": "t12+t13+t23-z-4",
        "cycle_formula_passed": sp.expand(beta - invariant_beta) == 0,
        "plus_z_formula_is_wrong": sp.expand(beta - wrong_plus_sign) != 0,
        "discriminant": str(discriminant),
        "discriminant_passed": sp.expand(
            discriminant - (4 - a**2) * (4 - b**2)
        )
        == 0,
        "minimum": str(minimum),
        "minimum_passed": sp.expand(
            minimum + (4 - a**2) * (4 - b**2) / 4
        )
        == 0,
        "lower_endpoint_negative_definite_example": {
            "parameters": "a=b=1,u=-1/2",
            "leading_principal_minors": [str(v) for v in equality_leading_minors],
            "passed": equality_leading_minors
            == [sp.Integer(-2), sp.Rational(15, 4), sp.Rational(-9, 2)],
        },
    }


def corrected_norm_identity() -> dict[str, object]:
    """Audit the exceptional-place factor in the beta norm.

    If L_rho=s_m^rho*s_n^rho/4, there are e extensions above every
    base embedding rho.  Exactly one extension above rho_0 is selected.
    Thus the product bound is

      256 * L_0^(e-1) * prod_{rho != rho_0} L_rho^e
      = 256/L_0 * (prod_rho L_rho)^e.
    """

    e = sp.symbols("e", integer=True, positive=True)
    L0, product_L = sp.symbols("L0 product_L", positive=True)
    lhs = 256 * L0 ** (e - 1) * (product_L / L0) ** e
    rhs = 256 / L0 * product_L**e
    return {
        "identity": "256*L0^(e-1)*(product_L/L0)^e=(256/L0)*product_L^e",
        "passed": sp.simplify(lhs - rhs) == 0,
        "substitution": "L0=s_m0*s_n0/4",
        "selected_factor": "1024/(s_m0*s_n0)",
        "polynomial_envelope": "<=4*M^4",
    }


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


def gamma(d: int) -> int:
    """Phi_d(1): p for a prime power p^r and 1 otherwise."""
    value = d
    prime = None
    p = 2
    while p * p <= value:
        if value % p == 0:
            if prime is not None:
                return 1
            prime = p
            while value % p == 0:
                value //= p
        p += 1
    if value > 1:
        if prime is not None:
            return 1
        prime = value
    return int(prime) if prime is not None else 1


def h_power_bound_regression(limit: int = 20_000) -> dict[str, object]:
    # h(d)=gamma(d)^(2/phi(d)).  To avoid floating point, h(d)<=3 is
    # equivalent to gamma(d)^2 <= 3^phi(d).
    failures = []
    equality = []
    for d in range(3, limit + 1):
        left = gamma(d) ** 2
        right = 3 ** euler_phi(d)
        if left > right:
            failures.append(d)
        if left == right:
            equality.append(d)
    return {
        "checked_through": limit,
        "failures": failures,
        "equality_cases": equality,
        "passed": not failures and equality == [3],
    }


def endpoint_cutoff() -> dict[str, object]:
    # For M>=10^6, h(M)<13/12 and h(other)<=3, so the cyclotomic
    # ratio is <(13/16)^D.  We use the weaker D>=353.
    lhs = 16**ENDPOINT_EXPONENT
    rhs = 4 * LABEL_CUTOFF**4 * 13**ENDPOINT_EXPONENT

    # A rational derivative certificate for
    # F(x)=4*x^4*(13/16)^(sqrt(x/2)/2):
    # log(16/13)>3/16, and sqrt(x)>=1000 at the cutoff.  The derivative
    # is negative once sqrt(x)>16*sqrt(2)/log(16/13), whose right side
    # is <512/3<1000.
    derivative_threshold_upper = Fraction(512, 3)
    return {
        "M": LABEL_CUTOFF,
        "D_weakened_lower_bound": ENDPOINT_EXPONENT,
        "comparison": "16^353 > 4*M^4*13^353",
        "lhs": str(lhs),
        "rhs": str(rhs),
        "endpoint_passed": lhs > rhs,
        "derivative_threshold_upper": str(derivative_threshold_upper),
        "monotonicity_passed": LABEL_CUTOFF**2 > derivative_threshold_upper**2,
    }


def large_label_h_bound() -> dict[str, object]:
    # phi(M)>=sqrt(M/2).  The function 2 log(x)/sqrt(x/2) decreases for
    # x>e^2. At 10^6, log(10)<7/3 and sqrt(500000)>700 give an upper
    # bound 1/25, while log(13/12)>1/13.
    return {
        "log_h_upper_at_cutoff": "<1/25",
        "log_13_over_12_lower": ">1/13",
        "endpoint_passed": Fraction(1, 25) < Fraction(1, 13),
        "decreasing_from_cutoff": True,
        "conclusion": "h(M)<13/12 for M>=10^6",
    }


def build_certificate() -> dict[str, object]:
    determinant = determinant_and_cartan_sign()
    norm = corrected_norm_identity()
    h_global = h_power_bound_regression()
    h_large = large_label_h_bound()
    endpoint = endpoint_cutoff()
    all_checks = [
        determinant["cycle_formula_passed"],
        determinant["plus_z_formula_is_wrong"],
        determinant["discriminant_passed"],
        determinant["minimum_passed"],
        determinant["lower_endpoint_negative_definite_example"]["passed"],
        norm["passed"],
        h_global["passed"],
        h_large["endpoint_passed"],
        endpoint["endpoint_passed"],
        endpoint["monotonicity_passed"],
    ]
    return {
        "schema": "compact-h3-general-fork-independent-hostile-v1",
        "all_checks_pass": all(all_checks),
        "determinant_and_cartan_sign": determinant,
        "corrected_norm": norm,
        "global_cyclotomic_factor": h_global,
        "large_label_factor": h_large,
        "label_cutoff": endpoint,
        "fekete_scope": {
            "base": "E=Q(s_m,s_n)",
            "contracting_intervals": "[-s_m^rho*s_n^rho/4,0]",
            "exceptional_intervals": 1,
            "exceptional_bound": "[0,256]",
        },
    }


def main() -> None:
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))
    if not certificate["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
