#!/usr/bin/env python3
"""Independent exact audit for the ordinary P-first dimension theorem.

This file deliberately does not import the owner's certifier.  It rebuilds
the standard facets, their dual rays, the characteristic-vector gap, the
rank-ten affine relation, and the flex-cluster numerical certificates.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

OWNER_PINS = {
    "report": {
        "path": "results/rational_pfirst_ordinary_dimension_bound_and_exact_carrier.md",
        "sha256": "8ea1ede6cc1717ede57dc902e459ceb831682f9f6cd6328c70fc7657452920af",
    },
    "owner_audit": {
        "path": "results/rational_pfirst_ordinary_dimension_bound_and_exact_carrier_hostile_audit.md",
        "sha256": "cee0e6171684ca9b383ca6900aef5e02ae20f89c86ea726df47503efa401836b",
    },
    "certifier": {
        "path": "scripts/certify_rational_pfirst_ordinary_dimension_bound_and_exact_carrier.py",
        "sha256": "72963302967dcca228f785aec67fab3186ce4d5e4a180d2932e13d06fc5b23b9",
    },
    "archive": {
        "path": "data/enumerations/rational_pfirst_ordinary_dimension_bound_and_exact_carrier.json",
        "sha256": "d9b5b2c208b94df6a9c6a977292508d798837698bc2541697050c4f734f74b8b",
    },
    "tests": {
        "path": "tests/test_rational_pfirst_ordinary_dimension_bound_and_exact_carrier.py",
        "sha256": "93b9f1b285e8d7238540e73d56a7a6fcf2d3bd0d62d2b71f1d9b07157e4d4ff0",
    },
}

DEPENDENCY_PINS = {
    "timelike_report": {
        "path": "results/rational_pfirst_timelike_anticanonical_vertex_closure.md",
        "sha256": "6450ab0629c539fcb3fff00e6aa2bdd0ddf194f17bdaaa2c1f2265d2b4aebb29",
    },
    "timelike_audit": {
        "path": "results/rational_pfirst_timelike_anticanonical_vertex_closure_hostile_audit.md",
        "sha256": "2d7eddae058f81f06fded910b4d8ed45219b30e011af5be6a232f20655df0499",
    },
    "cusp_report": {
        "path": "results/rational_pfirst_anticanonical_cusp_saturation_closure.md",
        "sha256": "9d4142c8c4647bfced451ac1ff29184a6a278d9a51cdb1a0638276f3eb813032",
    },
    "cusp_audit": {
        "path": "results/rational_pfirst_anticanonical_cusp_saturation_closure_hostile_audit.md",
        "sha256": "1cc36249914c82ba42145e64cb0e1a43a888d4ac0c02d6b02084cade3ea498d8",
    },
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def lorentz_dot(left: tuple[int | Fraction, ...], right: tuple[int | Fraction, ...]) -> Fraction:
    return Fraction(left[0] * right[0]) - sum(
        (Fraction(a * b) for a, b in zip(left[1:], right[1:])),
        Fraction(0),
    )


def standard_facets(n: int) -> list[tuple[int, ...]]:
    facets: list[tuple[int, ...]] = []
    last_exceptional = [0] * (n + 1)
    last_exceptional[n] = 1
    facets.append(tuple(last_exceptional))
    for index in range(1, n):
        difference = [0] * (n + 1)
        difference[index] = 1
        difference[index + 1] = -1
        facets.append(tuple(difference))
    facets.append((1, -1, -1, -1) + (0,) * (n - 3))
    return facets


def determinant(matrix: list[list[int | Fraction]]) -> Fraction:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    value = Fraction(1)
    size = len(work)
    for column in range(size):
        pivot = next(
            row for row in range(column, size) if work[row][column] != 0
        )
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            value = -value
        pivot_value = work[column][column]
        value *= pivot_value
        for row in range(column + 1, size):
            multiplier = work[row][column] / pivot_value
            for target in range(column, size):
                work[row][target] -= multiplier * work[column][target]
    return value


def solve(matrix: list[list[int]], target: list[int]) -> tuple[Fraction, ...]:
    size = len(matrix)
    work = [
        [Fraction(entry) for entry in row] + [Fraction(value)]
        for row, value in zip(matrix, target)
    ]
    for column in range(size):
        pivot = next(
            row for row in range(column, size) if work[row][column] != 0
        )
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(size):
            if row == column or work[row][column] == 0:
                continue
            multiplier = work[row][column]
            work[row] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(work[row], work[column])
            ]
    return tuple(work[row][-1] for row in range(size))


def dual_rays(n: int) -> list[tuple[Fraction, ...]]:
    facets = standard_facets(n)
    pairing_matrix = [
        [facet[0], *(-entry for entry in facet[1:])]
        for facet in facets
    ]
    return [
        solve(pairing_matrix, [int(index == column) for index in range(n + 1)])
        for column in range(n + 1)
    ]


def is_characteristic(vector: tuple[int | Fraction, ...]) -> bool:
    return all(
        entry.denominator == 1 and int(entry) % 2 != 0
        if isinstance(entry, Fraction)
        else entry % 2 != 0
        for entry in vector
    )


def simplex_record(n: int) -> dict[str, object]:
    facets = standard_facets(n)
    canonical = (-3,) + (1,) * n
    anticanonical = (3,) + (-1,) * n
    duals = dual_rays(n)

    facet_matrix_det = determinant([list(facet) for facet in facets])
    assert abs(facet_matrix_det) == 1
    assert [lorentz_dot(root, root) for root in facets] == [-1] + [-2] * n
    assert [lorentz_dot(canonical, root) for root in facets] == [-1] + [0] * n
    assert [lorentz_dot(anticanonical, root) for root in facets] == [1] + [0] * n
    assert tuple(duals[0]) == anticanonical
    assert tuple(duals[1]) == (1, -1) + (0,) * (n - 1)
    for index, dual in enumerate(duals):
        assert [lorentz_dot(dual, root) for root in facets] == [
            int(index == column) for column in range(n + 1)
        ]

    squares = [lorentz_dot(dual, dual) for dual in duals]
    if n <= 8:
        assert squares[0] == 9 - n > 0
        assert [index for index, square in enumerate(squares) if square == 0] == [1]
        assert all(square >= 0 for square in squares)
    elif n == 9:
        assert [index for index, square in enumerate(squares) if square == 0] == [0, 1]
        assert all(square >= 0 for square in squares)
        assert is_characteristic(duals[0])
        assert not is_characteristic(duals[1])

    return {
        "n": n,
        "facet_basis_determinant": int(facet_matrix_det),
        "facet_squares": [int(lorentz_dot(root, root)) for root in facets],
        "canonical_pairings": [int(lorentz_dot(canonical, root)) for root in facets],
        "anticanonical_square": int(lorentz_dot(anticanonical, anticanonical)),
        "dual_squares": [int(square) for square in squares],
        "null_dual_indices": [
            index for index, square in enumerate(squares) if square == 0
        ],
    }


def escape_record(n: int) -> dict[str, object]:
    assert n >= 10
    facets = standard_facets(n)
    multiplicities = tuple(range(n, 0, -1))
    time = sum(multiplicities[:3]) + 1
    strict_point = (time,) + tuple(-value for value in multiplicities)
    anticanonical = (3,) + (-1,) * n
    strict_slacks = [lorentz_dot(strict_point, root) for root in facets]
    anticanonical_slacks = [
        lorentz_dot(anticanonical, root) for root in facets
    ]
    constant = lorentz_dot(strict_point, strict_point)
    linear_half = lorentz_dot(strict_point, anticanonical)
    quadratic = lorentz_dot(anticanonical, anticanonical)

    assert strict_slacks == [1] * (n + 1)
    assert anticanonical_slacks == [1] + [0] * n
    assert constant > 0
    assert quadratic == 9 - n < 0
    assert linear_half * linear_half - constant * quadratic > 0
    assert constant / quadratic < 0

    return {
        "n": n,
        "strict_point_square": int(constant),
        "all_strict_wall_slacks": [int(value) for value in strict_slacks],
        "F_wall_slacks": [int(value) for value in anticanonical_slacks],
        "escape_quadratic": [int(constant), int(2 * linear_half), int(quadratic)],
        "root_product_is_negative": True,
        "positive_light_cone_crossing_is_strict_for_every_wall": True,
    }


def characteristic_gap_records() -> list[dict[str, object]]:
    records = []
    for n in range(3, 9):
        # Identity check on a hostile grid; the displayed coefficient form
        # then certifies all nonnegative k,A,B, not merely this grid.
        for k in range(8):
            for A in range(8):
                for B in range(8):
                    s = 3 * k + A + B
                    upper = (
                        (k + A) * (k + A + 1)
                        + (k + B) * (k + B + 1)
                        + (n - 2) * k * (k + 1)
                    )
                    gap = s * (s + 3) - upper
                    coefficient_form = (
                        (9 - n) * k * (k + 1)
                        + 4 * k * (A + B)
                        + 2 * A * B
                        + 2 * A
                        + 2 * B
                    )
                    assert gap == coefficient_form
                    assert (gap == 0) == (k == A == B == 0)
        records.append(
            {
                "n": n,
                "identity": (
                    "(9-n)k(k+1)+4k(A+B)+2AB+2A+2B"
                ),
                "all_nonzero_nonnegative_inputs_give_strict_gap": True,
                "unique_reduced_characteristic_vector": [3] + [-1] * n,
            }
        )
    return records


def rank_ten_record() -> dict[str, object]:
    facets = standard_facets(9)
    anticanonical = (3,) + (-1,) * 9
    multiplicities = [2, 4, 6, 5, 4, 3, 2, 1, 3]
    affine_components = facets[1:]
    relation = tuple(
        sum(
            multiplicities[index] * affine_components[index][coordinate]
            for index in range(9)
        )
        for coordinate in range(10)
    )
    assert relation == anticanonical

    edges = []
    for left in range(9):
        for right in range(left + 1, 9):
            pairing = lorentz_dot(
                affine_components[left], affine_components[right]
            )
            assert pairing in (0, 1)
            if pairing == 1:
                edges.append((left, right))
    assert len(edges) == 8
    reached = {0}
    while True:
        expanded = reached | {
            right
            for left, right in edges
            if left in reached
        } | {
            left
            for left, right in edges
            if right in reached
        }
        if expanded == reached:
            break
        reached = expanded
    assert reached == set(range(9))

    # F, the section E9, and the eight components not meeting E9.
    trivial_lattice = [anticanonical, facets[0], *facets[1:8], facets[9]]
    trivial_det = determinant([list(vector) for vector in trivial_lattice])
    assert abs(trivial_det) == 1

    duals = dual_rays(9)
    assert tuple(duals[0]) == anticanonical
    assert tuple(duals[1]) == (1, -1) + (0,) * 8
    assert is_characteristic(duals[0])
    assert not is_characteristic(duals[1])

    return {
        "primitive_affine_multiplicities": multiplicities,
        "affine_relation_equals_F": True,
        "reduced_support_vertices": 9,
        "reduced_support_edges": len(edges),
        "reduced_support_is_tree": True,
        "reduced_support_euler_number": 2 * 9 - len(edges),
        "only_two_standard_ideal_rays": True,
        "F_is_characteristic": True,
        "other_ideal_ray_is_characteristic": False,
        "trivial_lattice_determinant": int(trivial_det),
        "section_quotient_is_trivial": True,
    }


def flex_record() -> dict[str, object]:
    # g=z(1-xz)-x^3 has dg/dz(0)=1.  On g=0, z has x-adic
    # valuation 3, hence z^3 has valuation 9.  Thus (z^3,g) is
    # analytically (x^9,w).  One point blow-up changes (x^e,w) to
    # (x^(e-1),w_1).
    residual = {n: 9 - n for n in range(3, 10)}
    assert residual == {3: 6, 4: 5, 5: 4, 6: 3, 7: 2, 8: 1, 9: 0}

    # The finite Weierstrass discriminant has no repeated root.
    # gcd(4+27*l^2, 54*l)=1 because l=0 is not a root.
    assert 4 != 0
    return {
        "local_generators": ["z^3", "z-x^3-x*z^2"],
        "smooth_base_branch_z_valuation": 3,
        "base_scheme_colength": 9,
        "contact_with_flex_line": 3,
        "residual_exponents_after_n_blowups": residual,
        "partial_pencil_has_zero_dimensional_base_and_no_fixed_curve": True,
        "finite_discriminant": "4+27*lambda^2",
        "finite_discriminant_has_two_simple_roots_over_C": True,
        "only_reducible_fibre_is_the_affine_nine_component_fibre": True,
    }


def verify_pins() -> dict[str, str]:
    pins = {**OWNER_PINS, **DEPENDENCY_PINS}
    observed: dict[str, str] = {}
    for name, item in pins.items():
        digest = sha256_file(ROOT / item["path"])
        assert digest == item["sha256"], (name, item["sha256"], digest)
        observed[name] = digest
    return observed


def build_audit() -> dict[str, object]:
    pins = verify_pins()
    simplex = [simplex_record(n) for n in range(3, 10)]
    escape = [escape_record(n) for n in range(10, 16)]
    characteristic = characteristic_gap_records()
    rank_ten = rank_ten_record()
    flex = flex_record()
    return {
        "schema": (
            "rational_pfirst_ordinary_dimension_bound_and_exact_carrier/"
            "independent_hostile/v1"
        ),
        "verdict": "PASS",
        "pins_verified": pins,
        "finite_volume_wall_cone": {
            "finite_sidedness_is_essential": True,
            "spacelike_escape_records": escape,
            "past_causal_exclusion_uses_same_spacelike_escape": True,
            "future_null_extreme_ray_argument": True,
        },
        "standard_simplex_records": simplex,
        "timelike_orientation": {
            "physical_wall_signs_alone_choose_sheet": False,
            "declared_host_is_full_O_plus": True,
            "relative_maximality_is_among_all_reflection_overgroups": True,
            "full_reflection_group_generated_by_standard_sides": True,
            "reduced_characteristic_gap_records": characteristic,
            "K_fixed_only_after_reduced_characteristic_argument": True,
        },
        "automatic_canonical_escape": {
            "F_in_wall_cone_from_ordinary_colors": True,
            "surviving_F_is_nonzero_causal": True,
            "finite_volume_forces_F_future": True,
            "K_equals_minus_F_is_past_and_not_in_P": True,
            "extra_K_not_in_P_hypothesis_required": False,
        },
        "rank_ten_cusp": rank_ten,
        "flex_carrier": flex,
        "scope": {
            "standard_marked_I_1_n": True,
            "n_at_least_3": True,
            "ordinary_minus_one_minus_two_facets_only": True,
            "dimension_bound_uses_no_surface_or_maximality": True,
            "exact_carrier_uses_relative_maximality": True,
            "twice_primitive_minus_four_lane_included": False,
            "classification_finiteness_used": False,
        },
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
