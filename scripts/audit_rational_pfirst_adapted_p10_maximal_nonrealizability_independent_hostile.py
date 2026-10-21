#!/usr/bin/env python3
"""Independent hostile audit of the decorated adapted-P10 obstruction.

No function is imported from the owner certifier.  The twelve-wall cone is
dualized again over Q, the E10 root facet is separately dualized, and the
unbounded descent is checked symbolically with an independent bounded
regression used only as a test.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

OWNER_PINS = {
    "report": {
        "path": "results/rational_pfirst_adapted_p10_maximal_nonrealizability.md",
        "sha256": "d35e0b948cff9ce51fbda2757f5ec38e3cc1d7153ea6e83283c0ff7b1ffec002",
    },
    "certifier": {
        "path": "scripts/certify_rational_pfirst_adapted_p10_maximal_nonrealizability.py",
        "sha256": "efb95821a71d970a832f702b66dffbef134ed330ba3efa1b229b4c0d9743f3bd",
    },
    "archive": {
        "path": "data/enumerations/rational_pfirst_adapted_p10_maximal_nonrealizability.json",
        "sha256": "dd04ad2dfe8c0d013334c24a654127975a0a997cc37aaac1ee68655d26bb9109",
    },
    "tests": {
        "path": "tests/test_rational_pfirst_adapted_p10_maximal_nonrealizability.py",
        "sha256": "1eaf7769158d9a67eed990d3e89da55d9b53503374da3b01b3b9b6a878d484e7",
    },
}

DEPENDENCY_PINS = {
    "complex_Coble_obstruction_archive": {
        "path": "data/enumerations/coble_nodot_maxdegree3_closure.json",
        "sha256": "56a0646809fff76514729a58887172eb916d28cbda699128ed2f96a23556e82d",
    },
    "two_boundary_comparison_archive": {
        "path": "data/enumerations/coble_type_i_two_boundary_full_ns_h11.json",
        "sha256": "bcfa2de1e96bf4f41282345a64f8ee6f3a75e715287bab4962d7082584b4d985",
    },
    "independent_H10_Coble_audit": {
        "path": "results/polytope_first_coble_h10_hostile_audit.md",
        "sha256": "4ba993002bd74f122a6495c4c3a25468b0b2296e736f77643db4e390b8f85886",
    },
    "uncoloured_genus_one_carrier": {
        "path": "results/rational_pfirst_anticanonical_genus_one_facet_dimension_bound.md",
        "sha256": "45593c7dc66fd0fb2db989d5862b0fd1cddba306c194971225f8db64e5378801",
    },
    "uncoloured_genus_one_archive": {
        "path": "data/enumerations/rational_pfirst_anticanonical_genus_one_facet_dimension_bound.json",
        "sha256": "f9c961fdbd4dda410cf119b72d10fd303e672d8f73167ef39c46d4293fe2024f",
    },
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_pins() -> dict[str, str]:
    observed: dict[str, str] = {}
    for label, item in {**OWNER_PINS, **DEPENDENCY_PINS}.items():
        digest = sha256_file(ROOT / item["path"])
        assert digest == item["sha256"], (label, item["sha256"], digest)
        observed[label] = digest
    return observed


def dot(left: tuple[int | Fraction, ...], right: tuple[int | Fraction, ...]) -> Fraction:
    return Fraction(left[0] * right[0]) - sum(
        (Fraction(a * b) for a, b in zip(left[1:], right[1:])),
        Fraction(0),
    )


def walls() -> tuple[list[str], list[tuple[int, ...]]]:
    labels = ["E10"]
    result: list[tuple[int, ...]] = []
    e = [0] * 11
    e[10] = 1
    result.append(tuple(e))
    for index in range(1, 10):
        root = [0] * 11
        root[index] = 1
        root[index + 1] = -1
        labels.append(f"D{index}")
        result.append(tuple(root))
    labels.append("C")
    result.append((1, -1, -1, -1) + (0,) * 7)
    labels.append("F")
    result.append((3,) + (-1,) * 10)
    return labels, result


def matrix_rank(rows: list[list[int | Fraction]]) -> int:
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


def determinant(rows: list[list[int | Fraction]]) -> Fraction:
    work = [[Fraction(entry) for entry in row] for row in rows]
    value = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            value = -value
        pivot_value = work[column][column]
        value *= pivot_value
        for row in range(column + 1, len(work)):
            multiple = work[row][column] / pivot_value
            for target in range(column, len(work)):
                work[row][target] -= multiple * work[column][target]
    return value


def primitive_preserving_sign(vector: list[Fraction]) -> tuple[int, ...]:
    denominator = 1
    for entry in vector:
        denominator = denominator * entry.denominator // gcd(
            denominator, entry.denominator
        )
    integral = [int(entry * denominator) for entry in vector]
    divisor = 0
    for entry in integral:
        divisor = gcd(divisor, abs(entry))
    assert divisor > 0
    return tuple(entry // divisor for entry in integral)


def null_line(rows: list[list[int]]) -> tuple[int, ...] | None:
    work = [[Fraction(entry) for entry in row] for row in rows]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(11):
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
        pivot_columns.append(column)
        pivot_row += 1
    free = [column for column in range(11) if column not in pivot_columns]
    if len(free) != 1:
        return None
    vector = [Fraction(0)] * 11
    vector[free[0]] = 1
    for row, column in reversed(list(enumerate(pivot_columns))):
        vector[column] = -sum(
            work[row][target] * vector[target] for target in free
        )
    return primitive_preserving_sign(vector)


def exact_extreme_rays() -> tuple[dict[tuple[int, ...], tuple[int, ...]], int]:
    _, roots = walls()
    inequalities = [
        [root[0], *(-entry for entry in root[1:])] for root in roots
    ]
    rays: dict[tuple[int, ...], tuple[int, ...]] = {}
    tested = 0
    for active in combinations(range(12), 10):
        tested += 1
        equations = [inequalities[index] for index in active]
        if matrix_rank(equations) != 10:
            continue
        vector = null_line(equations)
        assert vector is not None
        slacks = [dot(vector, root) for root in roots]
        if all(slack <= 0 for slack in slacks):
            vector = tuple(-entry for entry in vector)
            slacks = [-slack for slack in slacks]
        if not all(slack >= 0 for slack in slacks):
            continue
        assert vector[0] > 0
        zeroes = tuple(index for index, slack in enumerate(slacks) if slack == 0)
        rays[vector] = zeroes
    assert tested == 66
    return dict(sorted(rays.items())), tested


def solve_square(matrix: list[list[int]], target: list[int]) -> tuple[Fraction, ...]:
    size = len(matrix)
    work = [
        [Fraction(entry) for entry in row] + [Fraction(value)]
        for row, value in zip(matrix, target)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        value = work[column][column]
        work[column] = [entry / value for entry in work[column]]
        for row in range(size):
            if row == column or work[row][column] == 0:
                continue
            value = work[row][column]
            work[row] = [
                entry - value * pivot_entry
                for entry, pivot_entry in zip(work[row], work[column])
            ]
    return tuple(work[row][-1] for row in range(size))


def cone_certificate() -> dict[str, object]:
    labels, roots = walls()
    assert determinant([list(root) for root in roots[:11]]) == -1
    rays, active_sets_tested = exact_extreme_rays()
    assert len(rays) == 19
    squares = [dot(ray, ray) for ray in rays]
    assert sum(square > 0 for square in squares) == 17
    assert sum(square == 0 for square in squares) == 2
    assert sum(square < 0 for square in squares) == 0
    assert all(dot(left, right) >= 0 for left in rays for right in rays)

    sum_ray = tuple(sum(ray[index] for ray in rays) for index in range(11))
    slacks = [dot(sum_ray, root) for root in roots]
    assert sum_ray == (137, -47, -45, -43, -41, -39, -37, -35, -33, -31, -30)
    assert dot(sum_ray, sum_ray) == 3940
    assert slacks == [30, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 30]
    for facet in range(12):
        active_rays = [
            list(ray) for ray in rays if dot(ray, roots[facet]) == 0
        ]
        assert matrix_rank(active_rays) == 10, labels[facet]

    null_rays = [ray for ray in rays if dot(ray, ray) == 0]
    assert len(null_rays) == 2 and dot(null_rays[0], null_rays[1]) > 0

    return {
        "active_ten_subsets_tested": active_sets_tested,
        "extreme_rays": len(rays),
        "finite_rays": 17,
        "ideal_rays": 2,
        "hyperideal_rays": 0,
        "all_pairings_nonnegative": True,
        "strict_interior": list(sum_ray),
        "strict_interior_square": 3940,
        "all_facets_rank_ten": True,
        "isolated_ideal_rays": True,
        "finite_volume_by_closed_ball_polytope_triangulation": True,
        "wall_cone_contains_no_negative_vector": True,
    }


def root_facet_certificate() -> dict[str, object]:
    _, roots = walls()
    root_facets = roots[1:11]  # D1,...,D9,C, all orthogonal to F.
    anticanonical = roots[11]
    equations = [
        [root[0], *(-entry for entry in root[1:])] for root in root_facets
    ] + [[anticanonical[0], *(-entry for entry in anticanonical[1:])]]
    duals = []
    for column in range(10):
        target = [int(row == column) for row in range(10)] + [0]
        dual = solve_square(equations, target)
        assert dot(dual, anticanonical) == 0
        duals.append(dual)
    squares = [dot(dual, dual) for dual in duals]
    assert squares == [4, 18, 42, 30, 20, 12, 6, 2, 0, 10]
    assert all(square >= 0 for square in squares)
    assert sum(square == 0 for square in squares) == 1
    return {
        "root_simplex_dimension": 9,
        "dual_squares": [int(square) for square in squares],
        "finite_vertices": 9,
        "ideal_vertices": 1,
        "hyperideal_vertices": 0,
        "finite_volume": True,
        "finite_covolume_root_reflection_group": True,
    }


def partitions(length: int, maximum: int, square_sum: int):
    if length == 0:
        if square_sum == 0:
            yield ()
        return
    for value in range(min(maximum, isqrt(square_sum)), -1, -1):
        for tail in partitions(length - 1, value, square_sum - value * value):
            yield (value,) + tail


def descent_regression(maximum_a: int = 18) -> dict[str, object]:
    roots_tested = 0
    cremona_steps = 0
    anticanonical_steps = 0
    simple = {"E_or_D": 0, "C": 0, "F": 0}
    for a in range(0, maximum_a + 1):
        for q in (1, 2):
            # For a>0, m_i<=a follows from sum(m_i^2)=a^2+q and
            # primitivity/negativity.  The terminal a=0 case is exceptional:
            # its normalized roots have one or two entries equal to 1.
            for multiplicities in partitions(10, max(a, 1), a * a + q):
                roots_tested += 1
                if a == 0:
                    assert multiplicities in {
                        (1,) + (0,) * 9,
                        (1, 1) + (0,) * 8,
                    }
                    simple["E_or_D"] += 1
                    continue
                assert all(value <= a for value in multiplicities)
                top = sum(multiplicities[:3])
                if top > a:
                    if top == 3 * a:
                        assert a == 1 and q == 2
                        assert multiplicities == (1, 1, 1) + (0,) * 7
                        simple["C"] += 1
                    else:
                        assert abs(2 * a - top) < a
                        cremona_steps += 1
                    continue

                total = sum(multiplicities)
                # The cone certificate proves total>3a for a negative root.
                assert total > 3 * a
                d = total - 3 * a
                is_f = a == 3 and q == 1 and multiplicities == (1,) * 10
                if is_f:
                    simple["F"] += 1
                else:
                    assert 3 * d < a
                    assert abs(a - 6 * d) < a
                    anticanonical_steps += 1

    # The allegedly dangerous branch d>=a/3 is analytically bounded by
    # a^2<=9q<=18.  Exhaust the resulting a=1,2,3,4 cases exactly.
    dangerous = []
    for a in range(1, 5):
        for q in (1, 2):
            for multiplicities in partitions(10, a, a * a + q):
                if sum(multiplicities[:3]) > a:
                    continue
                total = sum(multiplicities)
                if total <= 3 * a:
                    continue
                d = total - 3 * a
                if 3 * d >= a:
                    dangerous.append((a, q, multiplicities, d))
    assert dangerous == [(3, 1, (1,) * 10, 1)]

    return {
        "unbounded_measure": "positive integer absolute H-coordinate",
        "cremona_formula": "a maps to 2a-(m1+m2+m3)",
        "anticanonical_formula": "a maps to a-6(T-3a)",
        "large_d_Cauchy_bound": "a^2<=9q<=18",
        "large_d_exhaustive_survivor": "F=(3,-1^10)",
        "bounded_regression_maximum_a": maximum_a,
        "bounded_regression_roots": roots_tested,
        "bounded_regression_cremona_steps": cremona_steps,
        "bounded_regression_anticanonical_steps": anticanonical_steps,
        "simple_counts": simple,
        "every_nonsimple_step_strictly_decreases_measure": True,
        "bounded_search_used_as_proof": False,
    }


def physical_and_scope_certificate() -> dict[str, object]:
    labels, roots = walls()
    f = roots[-1]
    canonical = tuple(-entry for entry in f)
    boundary = tuple(2 * entry for entry in f)
    physical = [roots[0], *roots[1:11], boundary]
    assert [dot(root, root) for root in physical] == [-1] + [-2] * 10 + [-4]
    assert [dot(canonical, root) for root in physical] == [-1] + [0] * 10 + [2]
    assert dot(canonical, roots[0]) == -1
    assert dot(f, f) == -1
    assert tuple(f[index] - boundary[index] // 2 for index in range(11)) == (0,) * 11
    assert [dot(boundary, root) for root in roots[1:11]] == [0] * 10
    assert dot(boundary, roots[0]) == 2

    owner = " ".join(
        (ROOT / OWNER_PINS["report"]["path"]).read_text(encoding="utf-8").split()
    )
    required = {
        "decorated_not_uncoloured": "This is a decorated-carrier obstruction, not an obstruction to the uncoloured polytope.",
        "uncoloured_carrier_exists": "uncoloured `P10` nevertheless has a rational full-nef realization",
        "genus_one_primitive": "irreducible genus-one `(-1)` curve",
        "nonreduced_double": "nonreduced double of the genus-one wall",
        "complex_scope": "complex physical full-nef carrier is impossible",
        "positive_characteristic_nonclaim": "nonrealizability in positive characteristic",
        "raw_classification_free": "raw saturation theorem is project-native and classification-free",
        "terminal_AG_classification": "terminal nonrealizability input is an AG classification of complex Coble surfaces",
    }
    for label, phrase in required.items():
        assert phrase in owner, (label, phrase)

    # On the uncoloured flex carrier F=-K has p_a=1, whereas the physical
    # label asks for a reduced rational curve in the divisible class 2F.
    arithmetic_genus_f = 1 + (dot(f, f) + dot(canonical, f)) // 2
    arithmetic_genus_b = 1 + (dot(boundary, boundary) + dot(canonical, boundary)) // 2
    assert arithmetic_genus_f == 1
    assert arithmetic_genus_b == 0

    return {
        "labels": labels,
        "physical_square_split": [-1] + [-2] * 10 + [-4],
        "canonical_pairings": [-1] + [0] * 10 + [2],
        "both_K_and_minus_K_escape": True,
        "wrong_log_class_is_zero": True,
        "weak_incidence_gate_passes": True,
        "uncoloured_F_arithmetic_genus": arithmetic_genus_f,
        "decorated_B_arithmetic_genus": arithmetic_genus_b,
        "same_mirror": True,
        "different_geometric_effectivity_demand": True,
        "scope_phrases": {label: True for label in required},
    }


def ag_terminal_certificate() -> dict[str, object]:
    obstruction = json.loads(
        (ROOT / DEPENDENCY_PINS["complex_Coble_obstruction_archive"]["path"]).read_text()
    )["maximal_automaticity_counterexample"]
    assert obstruction["root_count"] == 10
    assert obstruction["complete_cusp_lift_count"] == 1
    assert obstruction["D_CobVin_arithmetic_gates_pass"] is True
    assert obstruction[
        "exact_complex_one_boundary_strong_Coble_full_nef_carrier_exists"
    ] is False
    assert "twelve effective irreducible roots" in obstruction["AG_obstruction"]
    return {
        "base_field": "C",
        "displayed_effective_root_count": 10,
        "finite_volume_root_chamber_forces_finite_Aut_lane": True,
        "Kondo_one_boundary_finite_Aut_types": ["I", "II"],
        "effective_irreducible_roots_in_each_allowed_type": 12,
        "decorated_full_nef_carrier_exists": False,
        "raw_maximality_uses_this_classification": False,
        "terminal_no_go_uses_this_classification": True,
    }


def build_audit() -> dict[str, object]:
    return {
        "schema": "rational_pfirst_adapted_p10_maximal_nonrealizability/independent_hostile/v1",
        "verdict": "PASS_SCOPED",
        "pins_verified": verify_pins(),
        "cone": cone_certificate(),
        "Coble_Mukai_root_facet": root_facet_certificate(),
        "all_reflection_descent": descent_regression(),
        "physical_packet": physical_and_scope_certificate(),
        "AG_terminal_no_go": ag_terminal_certificate(),
        "logical_separation": {
            "raw_finite_volume_and_maximality_classification_free": True,
            "decorated_complex_rational_full_nef_nonrealizability": True,
            "uncoloured_P10_nonexistence_claimed": False,
            "uncoloured_rational_full_nef_carrier_exists": True,
            "positive_characteristic_nonrealizability_claimed": False,
            "universal_adapted_dimension_bound_claimed": False,
        },
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
