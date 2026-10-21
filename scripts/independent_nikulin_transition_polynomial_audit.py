#!/usr/bin/env python3
"""Exact hostile audit of equation (24) in arXiv:0910.5217v4.

This script does not use Nikulin's degree theorem as an input.  It checks the
single polynomial certificate printed on pages 15--16, exhibits an exact
rational counterexample to the printed sup norm, computes its true sup norm,
and verifies the minimal one-exponent repair suggested by the phrase
"similar to (13)" in the paper.

The mathematical conclusion is deliberately narrow: the printed certificate
does not prove the stated degree-23 subcase.  Replacing the final exponent 2
by 3 repairs that subcase.  This does not independently prove the global
N(14) theorem, whose remaining geometric and graph reductions require their
own certificates.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/nikulin_transition_polynomial_audit.json"
SCHEMA = "coxeter.nikulin_transition_polynomial_hostile_audit.v1"


def _rational_string(value: sp.Expr) -> str:
    value = sp.Rational(value)
    return f"{value.p}/{value.q}"


def _assert_positive(value: sp.Expr, label: str) -> None:
    if sp.simplify(value).is_positive is not True:
        raise ArithmeticError(f"could not certify {label} > 0: {value}")


def _assert_negative(value: sp.Expr, label: str) -> None:
    if sp.simplify(value).is_negative is not True:
        raise ArithmeticError(f"could not certify {label} < 0: {value}")


def build_certificate() -> dict[str, Any]:
    x, y, z = sp.symbols("x y z")

    printed = (
        (x + 2) ** 3
        * (x + 1) ** 4
        * x**4
        * (x - 1) ** 3
        * (x**2 + x - 1) ** 2
    )
    claimed_delta = sp.Rational(884736, 9765625)  # 0.0905969664 exactly.

    # A one-point exact refutation already suffices: x=3/4 lies in [-2,1].
    witness_x = sp.Rational(3, 4)
    witness_value = sp.Abs(printed.subs(x, witness_x))
    _assert_positive(witness_value - claimed_delta, "rational witness excess")

    # The symmetry x -> -1-x reduces the full optimization to a degree-nine
    # polynomial in y=x(x+1), whose range on [-2,1] is [-1/4,2].
    absolute_reduction = y**4 * (y - 1) ** 2 * (2 - y) ** 3
    if sp.expand(printed.subs(x, -1 - x) - printed) != 0:
        raise ArithmeticError("the printed polynomial lost its reflection symmetry")
    if sp.expand(printed + absolute_reduction.subs(y, x * (x + 1))) != 0:
        raise ArithmeticError("the y=x(x+1) reduction failed")

    derivative = sp.factor(sp.diff(absolute_reduction, y))
    expected_derivative = -y**3 * (y - 2) ** 2 * (y - 1) * (9 * y**2 - 19 * y + 8)
    if sp.expand(derivative - expected_derivative) != 0:
        raise ArithmeticError("unexpected derivative factorization")

    sqrt73 = sp.sqrt(73)
    y_minus = (19 - sqrt73) / 18
    y_plus = (19 + sqrt73) / 18
    true_delta = sp.collect(
        sp.radsimp(sp.expand(absolute_reduction.subs(y, y_plus))), sqrt73
    )
    expected_true_delta = (sp.Integer(42029518) + 2327386 * sqrt73) / 387420489
    if sp.simplify(true_delta - expected_true_delta) != 0:
        raise ArithmeticError("true maximum expression changed")

    candidate_values = {
        "y=-1/4": absolute_reduction.subs(y, -sp.Rational(1, 4)),
        "y=0": sp.Integer(0),
        "y=(19-sqrt(73))/18": absolute_reduction.subs(y, y_minus),
        "y=1": sp.Integer(0),
        "y=(19+sqrt(73))/18": true_delta,
        "y=2": sp.Integer(0),
    }
    for name, value in candidate_values.items():
        if name != "y=(19+sqrt(73))/18":
            _assert_positive(true_delta - value, f"maximum minus {name}")

    minpoly = sp.Poly(sp.minimal_polynomial(true_delta, z), z)
    if minpoly.all_coeffs() != [387420489, -84059036, 3538944]:
        raise ArithmeticError(f"unexpected maximum minpoly: {minpoly.as_expr()}")

    printed_at_14 = sp.Integer(printed.subs(x, 14))
    # B=1+log(P(14))/(-log(delta)).  B>27 iff P(14)*delta^26>1;
    # B<28 iff P(14)*delta^27<1.  These are algebraic sign tests.
    printed_lower_test = sp.expand(printed_at_14 * true_delta**26 - 1)
    printed_upper_test = sp.expand(printed_at_14 * true_delta**27 - 1)
    _assert_positive(printed_lower_test, "printed bound minus 27")
    _assert_negative(printed_upper_test, "printed bound minus 28")

    # Minimal repair: exponent 2 -> 3.  It is exactly equation (13)
    # translated from [0,3] to [-2,1].
    repaired = printed * (x**2 + x - 1)
    equation_13 = (
        z**3
        * (z - 1) ** 4
        * (z - 2) ** 4
        * (z - 3) ** 3
        * (z**2 - 3 * z + 1) ** 3
    )
    if sp.expand(equation_13.subs(z, x + 2) - repaired) != 0:
        raise ArithmeticError("the repaired polynomial is not translated equation (13)")

    repaired_reduction = y**4 * (y - 1) ** 3 * (y - 2) ** 3
    if sp.expand(repaired - repaired_reduction.subs(y, x * (x + 1))) != 0:
        raise ArithmeticError("repaired y-reduction failed")

    repaired_critical_y = sp.Rational(8, 5)
    repaired_delta = sp.Abs(repaired_reduction.subs(y, repaired_critical_y))
    if repaired_delta != claimed_delta:
        raise ArithmeticError((repaired_delta, claimed_delta))

    repaired_derivative = sp.factor(sp.diff(-repaired_reduction, y))
    expected_repaired_derivative = (
        -y**3
        * (y - 2) ** 2
        * (y - 1) ** 2
        * (2 * y - 1)
        * (5 * y - 8)
    )
    if sp.expand(repaired_derivative - expected_repaired_derivative) != 0:
        raise ArithmeticError(
            f"unexpected repaired derivative factorization: {repaired_derivative}"
        )

    repaired_at_14 = sp.Integer(repaired.subs(x, 14))
    repaired_candidates = {
        "y=-1/4": sp.Abs(repaired_reduction.subs(y, -sp.Rational(1, 4))),
        "y=0": sp.Integer(0),
        "y=1/2": sp.Abs(repaired_reduction.subs(y, sp.Rational(1, 2))),
        "y=1": sp.Integer(0),
        "y=8/5": repaired_delta,
        "y=2": sp.Integer(0),
    }
    for name, value in repaired_candidates.items():
        if name != "y=8/5":
            _assert_positive(repaired_delta - value, f"repaired maximum minus {name}")
    # The repaired real bound lies strictly between 23 and 24.
    repaired_lower_test = repaired_at_14 * repaired_delta**22 - 1
    repaired_upper_test = repaired_at_14 * repaired_delta**23 - 1
    _assert_positive(repaired_lower_test, "repaired bound minus 23")
    _assert_negative(repaired_upper_test, "repaired bound minus 24")

    checks = {
        "printed_degree_is_18": bool(sp.degree(printed, x) == 18),
        "claimed_delta_is_exact_decimal": str(sp.N(claimed_delta, 10)) == "0.09059696640",
        "rational_point_refutes_claimed_supremum": bool(witness_value > claimed_delta),
        "printed_true_supremum_is_certified": bool(all(
            true_delta >= value for value in candidate_values.values()
        )),
        "printed_lemma7_bound_is_between_27_and_28": (
            printed_lower_test.is_positive is True
            and printed_upper_test.is_negative is True
        ),
        "printed_certificate_does_not_prove_degree_23": True,
        "repair_changes_only_last_exponent_2_to_3": (
            sp.simplify(repaired / printed) == x**2 + x - 1
        ),
        "repair_is_translated_equation_13": (
            sp.expand(equation_13.subs(z, x + 2) - repaired) == 0
        ),
        "repaired_degree_is_20": bool(sp.degree(repaired, x) == 20),
        "repaired_supremum_is_claimed_delta": bool(repaired_delta == claimed_delta),
        "repaired_lemma7_bound_is_between_23_and_24": (
            bool(repaired_lower_test > 0 and repaired_upper_test < 0)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(checks)

    return {
        "schema": SCHEMA,
        "status": "arxiv_v4_certificate_fails_journal_repair_passes",
        "source": {
            "paper": "V. V. Nikulin, The transition constant for arithmetic hyperbolic reflection groups",
            "arxiv": "0910.5217v4",
            "arxiv_date": "2010-01-23",
            "source_pdf_sha256": "e5391ab018e6e3f92ba9e2af8c35aff79332f973f3021f9c085148f52336f54b",
            "printed_location": "pages 15--16, equation (24), case s=k=3",
            "journal_version": {
                "citation": "Izvestiya: Mathematics 75:5 (2011), 971--1005",
                "doi": "10.1070/IM2011v075n05ABEH002561",
                "source_pdf_sha256": "716e6e7333efbfaa22060b60dc97e8303363a4914409ccdd56ef7a90c923c7af",
                "corrected_location": "page 986, equation (4.19), case s=k=3",
                "corrected_exponent": 3,
            },
        },
        "printed_certificate": {
            "polynomial": "(x+2)^3*(x+1)^4*x^4*(x-1)^3*(x^2+x-1)^2",
            "degree": 18,
            "interval": ["-2", "1"],
            "claimed_supremum": _rational_string(claimed_delta),
            "rational_counterexample": {
                "x": _rational_string(witness_x),
                "absolute_value": _rational_string(witness_value),
                "strictly_exceeds_claimed_supremum": True,
            },
            "exact_reduction": {
                "substitution": "y=x*(x+1)",
                "y_interval": ["-1/4", "2"],
                "absolute_polynomial": "y^4*(y-1)^2*(2-y)^3",
                "derivative_factorization": "-y^3*(y-2)^2*(y-1)*(9*y^2-19*y+8)",
            },
            "true_supremum": {
                "critical_y": "(19+sqrt(73))/18",
                "critical_x_pair": [
                    "-1/2-sqrt(47+2*sqrt(73))/6",
                    "-1/2+sqrt(47+2*sqrt(73))/6",
                ],
                "exact_value": "(42029518+2327386*sqrt(73))/387420489",
                "minimal_polynomial_coefficients": [387420489, -84059036, 3538944],
                "decimal_50": str(sp.N(true_delta, 50)),
            },
            "lemma_7_at_t_14": {
                "P_14": int(printed_at_14),
                "exact_sign_test_lower": "P(14)*delta^26>1",
                "exact_sign_test_upper": "P(14)*delta^27<1",
                "real_bound_decimal_30": str(
                    sp.N(1 + sp.log(printed_at_14) / (-sp.log(true_delta)), 30)
                ),
                "integer_degree_bound_from_this_certificate": 27,
                "printed_claimed_degree_bound": 23,
                "printed_claim_is_established_by_printed_certificate": False,
            },
        },
        "minimal_repair": {
            "change": "replace the final exponent 2 in equation (24) by 3",
            "polynomial": "(x+2)^3*(x+1)^4*x^4*(x-1)^3*(x^2+x-1)^3",
            "degree": 20,
            "identity": "P_repaired(x)=P_equation_13(x+2)",
            "exact_supremum": _rational_string(repaired_delta),
            "supremum_decimal": "0.0905969664",
            "critical_y": _rational_string(repaired_critical_y),
            "lemma_7_at_t_14": {
                "P_14": int(repaired_at_14),
                "exact_sign_test_lower": "P(14)*delta^22>1",
                "exact_sign_test_upper": "P(14)*delta^23<1",
                "real_bound_decimal_30": str(
                    sp.N(1 + sp.log(repaired_at_14) / (-sp.log(repaired_delta)), 30)
                ),
                "integer_degree_bound": 23,
            },
            "verdict": "independently_verified_repair_of_this_subcase",
        },
        "propagation_scope": {
            "printed_dependency_chain": [
                "equation (24), pages 15--16",
                "Theorem 11, page 19, s=k=3 branch",
                "Propositions 12 and 13, page 20",
                "Theorem 4 and its printed N(14)<=25 derivation",
                "Applications Theorems 18 and 19",
            ],
            "literal_v4_proof_status": "invalid_as_printed",
            "peer_reviewed_journal_status": (
                "equation (4.19) silently makes the exponent-3 repair; this "
                "specific gap does not affect the published journal theorem"
            ),
            "journal_dependency_chain": [
                "equation (4.19), page 986",
                "Theorem 4.6, page 989",
                "Propositions 4.1 and 4.2, pages 990--991",
                "Theorem 3.2 proof closure, page 998",
            ],
            "global_N14_theorem_disproved": False,
            "remaining_global_dependencies_independently_certified_here": False,
        },
        "methodological_scope": {
            "role_of_Nikulin": "hostile regression target, not a dependency of the project dimension theorem",
            "proved_here": (
                "exact failure of the arXiv-v4 polynomial and an independent "
                "verification of the repair later printed in the journal"
            ),
            "not_proved_here": "the global N(14) field theorem or an ABSW-type dimension theorem",
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
