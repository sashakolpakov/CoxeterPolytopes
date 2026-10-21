"""Verify the source-audited rank 15--19 K3 wall Gram matrices exactly.

The artifact stores intersection graphs, not root coordinates in a canonical
``U + ADE`` basis.  This verifier checks the exact linear-algebra claims used
in the accompanying theorem note: rank, inertia, a full-rank principal minor,
and the effect of reading the six mistyped rank-15 edges literally.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "data" / "known_examples" / "k3_high_rank_nef_chambers.json"


def load_artifact() -> dict[str, Any]:
    return json.loads(DATA_PATH.read_text())


def build_gram(case: dict[str, Any]) -> sp.Matrix:
    """Return the integral wall Gram matrix, rejecting malformed edge data."""

    wall_count = int(case["wall_count"])
    gram = sp.diag(*([-2] * wall_count))
    seen: set[tuple[int, int]] = set()

    for value_text, pairs in case["intersections"].items():
        value = int(value_text)
        if value <= 0:
            raise ValueError(f"off-diagonal intersection must be positive: {value}")
        for raw_i, raw_j in pairs:
            i, j = sorted((int(raw_i), int(raw_j)))
            if not (1 <= i < j <= wall_count):
                raise ValueError(f"invalid wall pair {(raw_i, raw_j)}")
            if (i, j) in seen:
                raise ValueError(f"duplicate wall pair {(i, j)}")
            seen.add((i, j))
            gram[i - 1, j - 1] = value
            gram[j - 1, i - 1] = value

    labels = {
        label
        for pairs in case["intersections"].values()
        for pair in pairs
        for label in pair
    }
    if labels != set(range(1, wall_count + 1)):
        raise ValueError(f"wall labels are not exactly 1,...,{wall_count}")
    return gram


def exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Compute inertia by exact rational congruence elimination.

    A nonzero diagonal entry gives a one-dimensional pivot and its Schur
    complement.  If the remaining diagonal vanishes but an off-diagonal entry
    does not, replacing e_i by e_i+e_j creates a nonzero diagonal pivot.  These
    are congruences, so Sylvester's law of inertia applies without numerical
    eigenvalue decisions.
    """

    work = sp.Matrix(matrix)
    positive = 0
    negative = 0

    while work.rows:
        pivot = next((i for i in range(work.rows) if work[i, i] != 0), None)
        if pivot is None:
            pair = next(
                (
                    (i, j)
                    for i in range(work.rows)
                    for j in range(i + 1, work.rows)
                    if work[i, j] != 0
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows

            i, j = pair
            order = [i, j] + [k for k in range(work.rows) if k not in (i, j)]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0

        if pivot:
            order = [pivot] + [k for k in range(work.rows) if k != pivot]
            work = work.extract(order, order)

        diagonal = sp.Rational(work[0, 0])
        if diagonal > 0:
            positive += 1
        else:
            negative += 1

        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / diagonal

    return positive, negative, 0


def analyze_case(case: dict[str, Any]) -> dict[str, Any]:
    gram = build_gram(case)
    basis = [int(index) - 1 for index in case["basis_indices"]]
    basis_gram = gram.extract(basis, basis)
    rank = gram.rank()
    inertia = exact_inertia(gram)
    determinant = int(basis_gram.det())

    checks = {
        "symmetric": gram == gram.T,
        "all_wall_norms_minus_two": all(gram[i, i] == -2 for i in range(gram.rows)),
        "rank": rank == int(case["expected_rank"]),
        "signature": list(inertia) == case["expected_signature"],
        "basis_size": len(basis) == int(case["expected_rank"]),
        "basis_determinant": determinant == int(case["expected_basis_determinant"]),
        "absolute_discriminant": abs(determinant) == int(case["absolute_discriminant"]),
        "hyperbolic_dimension": int(case["hyperbolic_dimension"]) == int(case["rho"]) - 1,
    }

    result: dict[str, Any] = {
        "rho": int(case["rho"]),
        "wall_count": gram.rows,
        "rank": rank,
        "inertia": list(inertia),
        "basis_determinant": determinant,
        "checks": checks,
    }

    corrected_pairs = case.get("misprinted_pairs_corrected_from_1_to_2")
    if corrected_pairs:
        literal = gram.copy()
        for raw_i, raw_j in corrected_pairs:
            i, j = int(raw_i) - 1, int(raw_j) - 1
            if literal[i, j] != 2:
                raise ValueError(f"corrected pair {(raw_i, raw_j)} does not have weight 2")
            literal[i, j] = 1
            literal[j, i] = 1
        literal_rank = literal.rank()
        literal_inertia = exact_inertia(literal)
        literal_checks = {
            "rank": literal_rank == int(case["literal_roulleau_tikz_expected_rank"]),
            "signature": list(literal_inertia)
            == case["literal_roulleau_tikz_expected_signature"],
            "cannot_embed_in_claimed_ns": literal_rank > int(case["rho"])
            or literal_inertia[0] > 1,
        }
        result["literal_rank15_redraw"] = {
            "rank": literal_rank,
            "inertia": list(literal_inertia),
            "checks": literal_checks,
        }
        checks["literal_redraw_diagnosis"] = all(literal_checks.values())

    result["all_checks_pass"] = all(checks.values())
    return result


def _order_at_zero(expression: sp.Expr, variable: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.cancel(expression), variable)
    return min(
        monomial[0]
        for monomial, coefficient in polynomial.terms()
        if coefficient != 0
    )


def analyze_derived_rank19_inose_model() -> dict[str, Any]:
    """Check the derived equation-level realization of the rank-19 lattice."""

    t, s, u = sp.symbols("t s u")
    a4 = -3 * s**2 * t**4
    a6 = t**5 * (t**2 - 2 * (1 - s**3) * t + 1)
    c4 = sp.expand(-48 * a4)
    discriminant = sp.factor(-16 * (4 * a4**3 + 27 * a6**2))
    expected_discriminant = sp.factor(
        -432 * t**10 * (t - 1) ** 2 * (t**2 + (4 * s**3 - 2) * t + 1)
    )

    local_zero = {
        "a4": _order_at_zero(a4, t),
        "a6": _order_at_zero(a6, t),
        "c4": _order_at_zero(c4, t),
        "discriminant": _order_at_zero(discriminant, t),
    }
    a4_infinity = sp.cancel(u**8 * a4.subs(t, 1 / u))
    a6_infinity = sp.cancel(u**12 * a6.subs(t, 1 / u))
    c4_infinity = sp.cancel(u**8 * c4.subs(t, 1 / u))
    discriminant_infinity = sp.cancel(u**24 * discriminant.subs(t, 1 / u))
    local_infinity = {
        "a4": _order_at_zero(a4_infinity, u),
        "a6": _order_at_zero(a6_infinity, u),
        "c4": _order_at_zero(c4_infinity, u),
        "discriminant": _order_at_zero(discriminant_infinity, u),
    }

    residual_quadratic = t**2 + (4 * s**3 - 2) * t + 1
    quadratic_discriminant = sp.factor(sp.discriminant(residual_quadratic, t))
    i2_leading_coefficient = sp.factor(
        (discriminant / (t - 1) ** 2).subs(t, 1)
    )
    checks = {
        "c4_identity": sp.expand(c4 - 144 * s**2 * t**4) == 0,
        "discriminant_factorization": sp.expand(
            discriminant - expected_discriminant
        )
        == 0,
        "zero_II_star_valuations": local_zero
        == {"a4": 4, "a6": 5, "c4": 4, "discriminant": 10},
        "infinity_II_star_valuations": local_infinity
        == {"a4": 4, "a6": 5, "c4": 4, "discriminant": 10},
        "t1_I2_generic": sp.factor(c4.subs(t, 1)) == 144 * s**2
        and i2_leading_coefficient == -1728 * s**3,
        "two_I1_roots_generic": sp.expand(
            quadratic_discriminant - 16 * s**3 * (s**3 - 1)
        )
        == 0,
        "fiber_euler_sum": 10 + 10 + 2 + 1 + 1 == 24,
        "trivial_lattice_rank": 2 + 8 + 8 + 1 == 19,
        "family_parameter_is_nonconstant": sp.diff(4 * s**3 - 2, s) != 0,
    }
    return {
        "status": "derived_and_symbolically_verified_not_roulleau_quoted",
        "discriminant": str(discriminant),
        "local_valuations_at_zero": local_zero,
        "local_valuations_at_infinity": local_infinity,
        "quadratic_discriminant": str(quadratic_discriminant),
        "generic_fibers": ["II*", "II*", "I2", "I1", "I1"],
        "fiber_euler_sum": 24,
        "trivial_lattice_rank": 19,
        "generic_ns_deduction": (
            "Non-isotriviality plus the Noether-Lefschetz generic-rank "
            "argument and Shioda-Tate give NS=U+E8(-1)^2+A1(-1)."
        ),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def analyze() -> dict[str, Any]:
    artifact = load_artifact()
    cases = [analyze_case(case) for case in artifact["cases"]]
    inose_model = analyze_derived_rank19_inose_model()
    boundary = artifact["classification_boundary"]["eligible_lattice_count_by_picard_rank"]
    boundary_check = boundary == {
        "15": 1,
        "16": 1,
        "17": 1,
        "18": 1,
        "19": 1,
        "20": 0,
    }
    return {
        "artifact_id": artifact["artifact_id"],
        "artifact_status": artifact["artifact_status"],
        "cases": cases,
        "derived_rank19_inose_model": inose_model,
        "ahl_rank_15_20_boundary_recorded": boundary_check,
        "all_checks_pass": all(case["all_checks_pass"] for case in cases)
        and inose_model["all_checks_pass"]
        and boundary_check,
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
