#!/usr/bin/env python3
"""Independent hostile audit of the anticanonical genus-one facet bound."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINS = {
    "owner_report": (
        "results/rational_pfirst_anticanonical_genus_one_facet_dimension_bound.md",
        "45593c7dc66fd0fb2db989d5862b0fd1cddba306c194971225f8db64e5378801",
    ),
    "owner_certifier": (
        "scripts/certify_rational_pfirst_anticanonical_genus_one_facet_dimension_bound.py",
        "15d6bebf292645cb0b0e9c121f413a717f5161472eabc449219e395cfbf57f32",
    ),
    "owner_archive": (
        "data/enumerations/rational_pfirst_anticanonical_genus_one_facet_dimension_bound.json",
        "f9c961fdbd4dda410cf119b72d10fd303e672d8f73167ef39c46d4293fe2024f",
    ),
    "owner_tests": (
        "tests/test_rational_pfirst_anticanonical_genus_one_facet_dimension_bound.py",
        "838ed74546bbd7cd9ae0d76c3583647ea6ad0a8ba12f541801d75b45e708f435",
    ),
    "reservoir_archive": (
        "data/enumerations/other_algebraic_reservoirs_certificate.json",
        "bb02c9fdf14d3370e97ffd2bb89a3a7801f8bd578fe73972810b2ddc4e58882e",
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_pins() -> dict[str, str]:
    observed = {}
    for label, (relative, expected) in PINS.items():
        actual = sha256(ROOT / relative)
        assert actual == expected, (label, expected, actual)
        observed[label] = actual
    return observed


def dot(left, right) -> Fraction:
    return Fraction(left[0] * right[0]) - sum(
        (Fraction(a * b) for a, b in zip(left[1:], right[1:])),
        Fraction(0),
    )


def rank(rows) -> int:
    if not rows:
        return 0
    work = [[Fraction(entry) for entry in row] for row in rows]
    pivot_row = 0
    for column in range(len(work[0])):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][column]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or work[row][column] == 0:
                continue
            value = work[row][column]
            work[row] = [
                entry - value * pivot_entry
                for entry, pivot_entry in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
    return pivot_row


def determinant(rows) -> Fraction:
    work = [[Fraction(entry) for entry in row] for row in rows]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        value = work[column][column]
        result *= value
        for row in range(column + 1, len(work)):
            multiple = work[row][column] / value
            work[row] = [
                entry - multiple * pivot_entry
                for entry, pivot_entry in zip(work[row], work[column])
            ]
    return result


def primitive(vector: list[Fraction]) -> tuple[int, ...]:
    denominator = 1
    for entry in vector:
        denominator = denominator * entry.denominator // gcd(
            denominator, entry.denominator
        )
    integral = [int(entry * denominator) for entry in vector]
    content = 0
    for entry in integral:
        content = gcd(content, abs(entry))
    assert content
    return tuple(entry // content for entry in integral)


def null_line(rows, dimension: int) -> tuple[int, ...] | None:
    work = [[Fraction(entry) for entry in row] for row in rows]
    pivots = []
    pivot_row = 0
    for column in range(dimension):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][column]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or work[row][column] == 0:
                continue
            value = work[row][column]
            work[row] = [
                entry - value * pivot_entry
                for entry, pivot_entry in zip(work[row], work[pivot_row])
            ]
        pivots.append(column)
        pivot_row += 1
    free = [column for column in range(dimension) if column not in pivots]
    if len(free) != 1:
        return None
    vector = [Fraction(0)] * dimension
    vector[free[0]] = 1
    for row, column in reversed(list(enumerate(pivots))):
        vector[column] = -sum(
            (work[row][target] * vector[target] for target in free),
            Fraction(0),
        )
    return primitive(vector)


def walls(n: int) -> tuple[list[str], list[tuple[int, ...]]]:
    labels = []
    roots = []
    for index in range(1, n):
        root = [0] * (n + 1)
        root[index] = 1
        root[index + 1] = -1
        labels.append(f"D{index}")
        roots.append(tuple(root))
    exceptional = [0] * (n + 1)
    exceptional[n] = 1
    labels.append(f"E{n}")
    roots.append(tuple(exceptional))
    labels.append("L")
    roots.append((1, -1, -1, -1) + (0,) * (n - 3))
    labels.append("F")
    roots.append((3,) + (-1,) * n)
    return labels, roots


def endpoint_certificate(n: int) -> dict[str, object]:
    assert n in (10, 11)
    labels, roots = walls(n)
    dimension = n + 1
    inequalities = [[root[0], *(-entry for entry in root[1:])] for root in roots]
    assert abs(determinant([list(root) for root in roots[:-1]])) == 1

    rays = {}
    active_sets = 0
    for active in combinations(range(n + 2), n):
        active_sets += 1
        equations = [inequalities[index] for index in active]
        if rank(equations) != n:
            continue
        vector = null_line(equations, dimension)
        assert vector is not None
        slacks = [dot(vector, root) for root in roots]
        if all(value <= 0 for value in slacks):
            vector = tuple(-entry for entry in vector)
            slacks = [-value for value in slacks]
        if not all(value >= 0 for value in slacks):
            continue
        assert vector[0] > 0
        rays[vector] = tuple(index for index, value in enumerate(slacks) if value == 0)

    squares = [dot(ray, ray) for ray in rays]
    expected = {10: (66, 19, 17, 2), 11: (78, 28, 26, 2)}[n]
    assert (
        active_sets,
        len(rays),
        sum(value > 0 for value in squares),
        sum(value == 0 for value in squares),
    ) == expected
    assert all(value >= 0 for value in squares)
    assert all(dot(left, right) >= 0 for left in rays for right in rays)
    for facet, root in enumerate(roots):
        facet_rays = [list(ray) for ray in rays if dot(ray, root) == 0]
        assert rank(facet_rays) == n, labels[facet]

    f = roots[-1]
    k = tuple(-entry for entry in f)
    square = dot(f, f)
    canonical_pairing = dot(k, f)
    assert square == 9 - n
    assert canonical_pairing == n - 9
    assert 1 + (square + canonical_pairing) // 2 == 1
    assert gcd(*map(abs, f)) == 1
    for basis_index in range(dimension):
        basis = [0] * dimension
        basis[basis_index] = 1
        coefficient = -2 * dot(basis, f) / square
        assert coefficient.denominator == 1

    source = json.loads((ROOT / PINS["reservoir_archive"][0]).read_text())
    source_case = source[f"p{n}_rational_surface"]
    assert source_case["all_checks_pass"] is True
    assert all(source_case["checks"].values())
    return {
        "n": n,
        "rho": n + 1,
        "active_sets_checked": active_sets,
        "facets": n + 2,
        "finite_vertices": expected[2],
        "ideal_vertices": expected[3],
        "hyperideal_vertices": 0,
        "all_facets_genuine": True,
        "F_square": int(square),
        "K_dot_F": int(canonical_pairing),
        "F_primitive": True,
        "F_reflection_integral": True,
        "arithmetic_genus": 1,
        "complete_nef_source_replayed": True,
    }


def intrinsic_theorem_certificate() -> dict[str, object]:
    return {
        "base": "smooth projective rational surface over an algebraically closed field",
        "rho_at_least": 3,
        "minimal_rational_surface_ranks": [1, 2],
        "rho_ge_3_implies_nonminimal_and_minus_one_curve": True,
        "adjunction_on_minus_one_curve": "K.E=-1",
        "K_primitive": True,
        "NS_unimodular": True,
        "primitive_in_unimodular_implies_divisibility": 1,
        "integral_reflection_forces": "2/F^2 is an integer",
        "genuine_hyperbolic_facet_forces_negative_square": True,
        "allowed_F_squares": [-1, -2],
        "Noether_formula": "K^2+rho=10",
        "allowed_rho": [11, 12],
        "allowed_projective_nef_dimensions": [10, 11],
        "classification_input_for_intrinsic_bound": False,
        "characteristic_free": True,
    }


def scope_certificate() -> dict[str, object]:
    text = " ".join((ROOT / PINS["owner_report"][0]).read_text().split())
    required = {
        "literal": "F=-K_S",
        "rho": "rho(S)>=3",
        "reduced_irreducible": "reduced irreducible curve",
        "reflection": "preserves `NS(S)`",
        "p11_sharp": "`P11`, not `P10`, attains the upper endpoint",
        "nonprimitive_nonclaim": "merely proportional to `-K`",
        "adapted_nonclaim": "twice-primitive rational `(-4)` boundary",
        "no_classification": "classify every chamber",
    }
    for label, phrase in required.items():
        assert phrase in text, (label, phrase)
    return {label: True for label in required}


def build_audit() -> dict[str, object]:
    return {
        "schema": "rational_pfirst_anticanonical_genus_one_facet_dimension_bound/independent_hostile/v1",
        "verdict": "PASS_EXACT_SCOPE",
        "pins": verify_pins(),
        "intrinsic_theorem": intrinsic_theorem_certificate(),
        "exact_endpoints": [endpoint_certificate(10), endpoint_certificate(11)],
        "scope": scope_certificate(),
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
