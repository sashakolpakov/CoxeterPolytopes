#!/usr/bin/env python3
"""Exact arithmetic/interface checks for the minimality-14 field funnel."""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TRANSITION_ARCHIVE = (
    ROOT / "data/enumerations/nikulin_four_node_transition_families.json"
)


def euler_phi(n: int) -> int:
    result = n
    p = 2
    value = n
    while p * p <= value:
        if value % p == 0:
            while value % p == 0:
                value //= p
            result -= result // p
        p += 1
    if value > 1:
        result -= result // value
    return result


def prime_power_base(n: int) -> int | None:
    """Return p for n=p^a, and None otherwise."""
    value = n
    found: int | None = None
    p = 2
    while p * p <= value:
        if value % p == 0:
            if found is not None:
                return None
            found = p
            while value % p == 0:
                value //= p
        p += 1
    if value > 1:
        if found is not None:
            return None
        found = value
    return found


def gamma_value(d: int) -> int:
    base = prime_power_base(d)
    return base if base is not None else 1


def lower_bound_totient_exact(n: int) -> bool:
    # Avoid square roots: phi(n)^2 >= n/2 iff 2 phi(n)^2 >= n.
    return 2 * euler_phi(n) ** 2 >= n


def tail_capacity_strict(d: int) -> bool:
    degree = euler_phi(d) // 2
    return gamma_value(d) < 4**degree


def cutoff_integer_comparison() -> bool:
    return 2**100 > 601 * 20_000**6


def monotonicity_certificate() -> dict[str, object]:
    # For x >= 20000, log(2)>2/3 gives
    # d/dx(log(2)*sqrt(x/2)-6log x)>0 once x>648.
    return {
        "log2_lower_bound": "log(2)>2/3",
        "derivative_positive_if_x_gt": 648,
        "applied_from": 20_000,
        "passed": 20_000 > 648,
    }


def transition_routes() -> dict[str, dict[str, object]]:
    archive = json.loads(TRANSITION_ARCHIVE.read_text())
    matrices = archive["symbolic_gram_matrices"]
    # Node numbers are one-based.  A tail meets exactly one of endpoints 1,2;
    # a fork meets both.  Gamma1/2 fork labels are bounded by the independently
    # certified spherical rank-three path grammar.
    declared = {
        "Gamma1": {"node": 3, "kind": "fork", "bounded_labels": True},
        "Gamma2": {"node": 3, "kind": "fork", "bounded_labels": True},
        "Gamma3": {"node": 3, "kind": "tail", "bounded_labels": False},
        "Gamma4": {"node": 4, "kind": "tail", "bounded_labels": False},
        "Gamma5": {"node": 3, "kind": "tail", "bounded_labels": False},
    }
    result: dict[str, dict[str, object]] = {}
    for family, route in declared.items():
        matrix = matrices[family]["gram_matrix"]
        j = int(route["node"]) - 1
        endpoint_entries = [matrix[0][j], matrix[1][j]]
        incidence = sum(entry != "0" for entry in endpoint_entries)
        expected = 1 if route["kind"] == "tail" else 2
        passed = incidence == expected
        if route["kind"] == "fork":
            passed = passed and archive["parameter_grammar_census"][
                "rank_three_path_lemma"
            ]["solutions_for_p_q_at_least_3_up_to_order"] == [
                [3, 3],
                [3, 4],
                [3, 5],
            ]
        result[family] = {
            **route,
            "endpoint_entries": endpoint_entries,
            "endpoint_incidence": incidence,
            "passed": passed,
        }
    return result


def h_vector_forcing_checks(max_n: int = 5000) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for n in range(10, max_n + 1):
        m = n - 1
        choose = m * (m - 1) // 2
        lower = Fraction(5 * choose, choose + m // 2)
        if n % 2 == 0:
            upper = Fraction(4, 1) + Fraction(4, n - 2)
        else:
            upper = Fraction(4, 1) + Fraction(4, n - 3)
        rows.append(
            {
                "n": n,
                "lower": str(lower),
                "strict_upper": str(upper),
                "contradiction": lower >= upper,
            }
        )
    return rows


def small_face_forcing_checks(max_n: int = 5000) -> list[dict[str, object]]:
    """Strict h-vector upper bound is at most 5 from n=6 onward."""
    rows: list[dict[str, object]] = []
    for n in range(6, max_n + 1):
        if n % 2 == 0:
            upper = Fraction(4, 1) + Fraction(4, n - 2)
        else:
            upper = Fraction(4, 1) + Fraction(4, n - 3)
        rows.append(
            {
                "n": n,
                "strict_upper": str(upper),
                "forces_triangle_or_quadrilateral": upper <= 5,
            }
        )
    return rows


def build_certificate() -> dict[str, object]:
    totient_checks = all(lower_bound_totient_exact(n) for n in range(3, 100_001))
    capacity_checks = all(tail_capacity_strict(d) for d in range(3, 20_000))
    routes = transition_routes()
    forcing = h_vector_forcing_checks()
    small_faces = small_face_forcing_checks()
    exact_at_cutoff = cutoff_integer_comparison()
    monotonicity = monotonicity_certificate()
    all_checks_pass = all(
        [
            totient_checks,
            capacity_checks,
            exact_at_cutoff,
            bool(monotonicity["passed"]),
            all(row["passed"] for row in routes.values()),
            all(row["contradiction"] for row in forcing),
            all(row["forces_triangle_or_quadrilateral"] for row in small_faces),
        ]
    )
    return {
        "schema": "global-minimality14-field-funnel-v1",
        "all_checks_pass": all_checks_pass,
        "tail_label_cutoff": {
            "excluded_from": 20_000,
            "largest_unexcluded_integer": 19_999,
            "exact_integer_comparison": exact_at_cutoff,
            "lhs": str(2**100),
            "rhs": str(601 * 20_000**6),
            "monotonicity": monotonicity,
        },
        "totient_bound": {
            "formula": "2*phi(d)^2 >= d",
            "exact_regression_through": 100_000,
            "passed": totient_checks,
            "proof": (
                "phi(d)^2/d=product_{p^a||d} p^(a-2)(p-1)^2; "
                "the sole possible factor below one is 1/2 from 2^1"
            ),
        },
        "tail_capacity": {
            "formula": "gamma(d)/4^(phi(d)/2) < 1",
            "checked_labels": [3, 19_999],
            "passed": capacity_checks,
            "samples": {
                str(d): {
                    "phi": euler_phi(d),
                    "gamma": gamma_value(d),
                    "capacity_numerator": gamma_value(d),
                    "capacity_denominator": 4 ** (euler_phi(d) // 2),
                }
                for d in [3, 4, 5, 6, 7, 8, 12, 30, 125]
            },
        },
        "transition_routes": routes,
        "dimension_forcing": {
            "range_checked": [10, 5000],
            "all_contradictions": all(row["contradiction"] for row in forcing),
            "endpoint_rows": [forcing[0], forcing[1], forcing[-1]],
        },
        "dimension_six_small_face_forcing": {
            "range_checked": [6, 5000],
            "strict_average_bound": "a<5",
            "all_force_triangle_or_quadrilateral": all(
                row["forces_triangle_or_quadrilateral"] for row in small_faces
            ),
            "endpoint_rows": [small_faces[0], small_faces[1], small_faces[-1]],
        },
    }


def main() -> None:
    certificate = build_certificate()
    print(json.dumps(certificate, indent=2, sort_keys=True))
    if not certificate["all_checks_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
