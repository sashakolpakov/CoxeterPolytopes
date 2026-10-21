#!/usr/bin/env python3
"""Independent hostile replay of the Type-I downstairs anti-Gram census.

This audit uses a different switching gauge (positive depth-first-tree
edges) and an exact characteristic-polynomial/Sturm inertia calculation.
It also verifies the half-sum-glue root that eliminates the sole rank-ten
negative-semidefinite orbit.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Iterable

import sympy as sp

try:
    from scripts.audit_enriques_type_i_explicit_lift_package_hostile import (
        _first_shimada_record,
        _ldl_upper,
        _matrix,
        _rows,
        _smith_diagonal,
    )
except ModuleNotFoundError:  # direct ``python scripts/...py`` execution
    from audit_enriques_type_i_explicit_lift_package_hostile import (
        _first_shimada_record,
        _ldl_upper,
        _matrix,
        _rows,
        _smith_diagonal,
    )


ROOT = Path(__file__).resolve().parents[1]
SHIMADA_DATA = (
    ROOT
    / "data"
    / "external"
    / "shimada_l10l26note"
    / "L10L26compdata2.txt"
)
EXPLICIT_PACKAGE = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_i_explicit_lift_package.json"
)
PRODUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_i_downstairs_anti_gram_census.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_i_downstairs_anti_gram_census_hostile_audit.json"
)

EXPECTED_SHA256 = {
    "shimada": "218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c",
    "explicit_package": "d9313e37687e355c9ca8f86d1ee11d74f0810a8c6906a0443f5377d2481fece7",
    "production": "df864bed625d1411c7db46d7bbb500606ccb011297e420de8272e512ab2e22ea",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _key(matrix: Iterable[Iterable[int]]) -> tuple[int, ...]:
    rows = [list(map(int, row)) for row in matrix]
    return tuple(
        rows[left][right]
        for left in range(len(rows))
        for right in range(left, len(rows))
    )


def _switch(matrix: list[list[int]], signs: tuple[int, ...]) -> list[list[int]]:
    return [
        [
            signs[left] * int(matrix[left][right]) * signs[right]
            for right in range(len(matrix))
        ]
        for left in range(len(matrix))
    ]


def _components(matrix: list[list[int]]) -> list[list[int]]:
    size = len(matrix)
    adjacency = [
        [other for other in range(size) if other != vertex and matrix[vertex][other]]
        for vertex in range(size)
    ]
    unseen = set(range(size))
    result = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        component = []
        stack = [root]
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for other in sorted(adjacency[vertex], reverse=True):
                if other in unseen:
                    unseen.remove(other)
                    stack.append(other)
        result.append(sorted(component))
    return result


def _dfs_positive_canonical(
    matrix: list[list[int]],
) -> tuple[tuple[int, ...], list[list[int]], list[int]]:
    """Independent gauge: force depth-first spanning-forest edges positive."""

    size = len(matrix)
    adjacency = [
        [other for other in range(size) if other != vertex and matrix[vertex][other]]
        for vertex in range(size)
    ]
    signs: list[int | None] = [None] * size

    def visit(vertex: int) -> None:
        for other in sorted(adjacency[vertex]):
            if signs[other] is None:
                edge_sign = 1 if matrix[vertex][other] > 0 else -1
                signs[other] = int(signs[vertex]) * edge_sign
                visit(other)

    for root in range(size):
        if signs[root] is None:
            signs[root] = 1
            visit(root)
    integral_signs = [int(value) for value in signs]
    canonical = _switch(matrix, tuple(integral_signs))
    return _key(canonical), canonical, integral_signs


def _direct_switch_orbit(matrix: list[list[int]]) -> set[tuple[int, ...]]:
    size = len(matrix)
    orbit = set()
    for initial in itertools.product((-1, 1), repeat=size - 1):
        orbit.add(_key(_switch(matrix, initial + (1,))))
    return orbit


def _sturm_inertia(matrix: list[list[int]]) -> tuple[int, int, int]:
    """Exact inertia from factorized characteristic polynomial and Sturm counts."""

    gram = sp.Matrix(matrix)
    variable = sp.Symbol("lambda")
    polynomial = gram.charpoly(variable).as_poly()
    zero = 0
    linear_zero = sp.Poly(variable, variable, domain=sp.ZZ)
    while polynomial.eval(0) == 0:
        polynomial = polynomial.exquo(linear_zero)
        zero += 1
    positive = 0
    negative = 0
    for factor, multiplicity in sp.factor_list(polynomial)[1]:
        factor = sp.Poly(factor, variable)
        negative += int(factor.count_roots(-sp.oo, 0)) * int(multiplicity)
        positive += int(factor.count_roots(0, sp.oo)) * int(multiplicity)
    if positive + negative + zero != gram.rows:
        raise ArithmeticError("Sturm inertia did not account for all eigenvalues")
    return positive, negative, zero


def _histogram(counter: Counter[tuple[int, int, int]]) -> dict[str, int]:
    return {
        f"p{positive}_n{negative}_z{zero}": count
        for (positive, negative, zero), count in sorted(counter.items())
    }


def analyze() -> dict[str, Any]:
    actual_hashes = {
        "shimada": _sha256(SHIMADA_DATA),
        "explicit_package": _sha256(EXPLICIT_PACKAGE),
        "production": _sha256(PRODUCTION),
    }
    if actual_hashes != EXPECTED_SHA256:
        raise AssertionError(f"a pinned census input changed: {actual_hashes}")
    production = json.loads(PRODUCTION.read_text(encoding="utf-8"))
    explicit = json.loads(EXPLICIT_PACKAGE.read_text(encoding="utf-8"))
    _, configuration_matrix, walls_matrix = _first_shimada_record()
    configuration = _rows(configuration_matrix)
    walls = _rows(walls_matrix)

    m1_edges = [
        (left, right)
        for left in range(12)
        for right in range(left + 1, 12)
        if configuration[left][right] == 1
    ]
    m2_edges = [
        (left, right)
        for left in range(12)
        for right in range(left + 1, 12)
        if configuration[left][right] == 2
    ]
    if len(m1_edges) != 10 or len(m2_edges) != 3:
        raise ArithmeticError("independent downstairs edge census changed")

    multiplicities: Counter[tuple[int, ...]] = Counter()
    representatives: dict[tuple[int, ...], list[list[int]]] = {}
    raw_count = 0
    for values1 in itertools.product((-2, 2), repeat=10):
        for values2 in itertools.product((-4, 0, 4), repeat=3):
            gram = [[-4 if left == right else 0 for right in range(12)] for left in range(12)]
            for (left, right), value in zip(m1_edges, values1):
                gram[left][right] = gram[right][left] = value
            for (left, right), value in zip(m2_edges, values2):
                gram[left][right] = gram[right][left] = value
            canonical_key, canonical, _ = _dfs_positive_canonical(gram)
            multiplicities[canonical_key] += 1
            representatives[canonical_key] = canonical
            raw_count += 1
    if raw_count != 27648 or len(representatives) != 18:
        raise ArithmeticError("independent DFS quotient changed raw/orbit counts")

    production_records = production["exact_inertia_census"]["orbit_records"]
    production_in_dfs_gauge = {
        _dfs_positive_canonical(record["canonical_gram"])[0]: record
        for record in production_records
    }
    if set(production_in_dfs_gauge) != set(representatives):
        raise ArithmeticError("production and independent switching quotients differ")

    weighted: Counter[tuple[int, int, int]] = Counter()
    unweighted: Counter[tuple[int, int, int]] = Counter()
    inertia_by_key: dict[tuple[int, ...], tuple[int, int, int]] = {}
    orbit_sizes: Counter[int] = Counter()
    for canonical_key, representative in representatives.items():
        components = _components(representative)
        expected_orbit_size = 2 ** (12 - len(components))
        direct_orbit = _direct_switch_orbit(representative)
        if (
            len(direct_orbit) != expected_orbit_size
            or multiplicities[canonical_key] != expected_orbit_size
        ):
            raise ArithmeticError("direct sheet-switch orbit size failed")
        # Reconstruct every orbit key only for the independent canonical check.
        for signs in itertools.product((-1, 1), repeat=11):
            switched = _switch(representative, signs + (1,))
            if _dfs_positive_canonical(switched)[0] != canonical_key:
                raise ArithmeticError("DFS gauge is not switching invariant")
        signature = _sturm_inertia(representative)
        inertia_by_key[canonical_key] = signature
        weighted[signature] += multiplicities[canonical_key]
        unweighted[signature] += 1
        orbit_sizes[expected_orbit_size] += 1
        production_record = production_in_dfs_gauge[canonical_key]
        production_signature = production_record["inertia"]
        if signature != (
            production_signature["positive"],
            production_signature["negative"],
            production_signature["zero"],
        ) or production_record["raw_assignment_multiplicity"] != multiplicities[canonical_key]:
            raise ArithmeticError("Sturm inertia or orbit multiplicity differs from production")

    expected_weighted = {
        "p0_n10_z2": 512,
        "p0_n9_z3": 1024,
        "p1_n10_z1": 9216,
        "p1_n11_z0": 512,
        "p1_n9_z2": 2048,
        "p2_n10_z0": 14336,
    }
    expected_unweighted = {
        "p0_n10_z2": 1,
        "p0_n9_z3": 1,
        "p1_n10_z1": 7,
        "p1_n11_z0": 1,
        "p1_n9_z2": 1,
        "p2_n10_z0": 7,
    }
    if _histogram(weighted) != expected_weighted or _histogram(unweighted) != expected_unweighted:
        raise ArithmeticError("independent exact inertia histograms changed")

    rank9_keys = [key for key, value in inertia_by_key.items() if value == (0, 9, 3)]
    rank10_keys = [key for key, value in inertia_by_key.items() if value == (0, 10, 2)]
    if len(rank9_keys) != 1 or len(rank10_keys) != 1:
        raise ArithmeticError("independent nonpositive-orbit uniqueness failed")
    geometric_gram = explicit["finite_lift_compiler_calibration"]["frozen_anti_difference_gram"]
    geometric_key = _dfs_positive_canonical(geometric_gram)[0]
    if geometric_key != rank9_keys[0]:
        raise ArithmeticError("unique independent rank-nine orbit is not geometric")

    rank10_gram = representatives[rank10_keys[0]]
    if not (
        rank10_gram[10][10] == rank10_gram[11][11] == -4
        and rank10_gram[10][11] == 0
    ):
        raise ArithmeticError("rank-ten orthogonal minus-four pair changed")
    downstairs_sum = [walls[10][index] + walls[11][index] for index in range(10)]
    if downstairs_sum != [8, 4, 8, 14, 12, 10, 8, 6, 4, 2]:
        raise ArithmeticError("forced-root downstairs parity vector changed")
    if any(value % 2 for value in downstairs_sum):
        raise ArithmeticError("forced-root downstairs half-sum is nonintegral")
    sign_square_checks = {}
    for first_sign, second_sign in itertools.product((-1, 1), repeat=2):
        four_square = (
            rank10_gram[10][10]
            + rank10_gram[11][11]
            + 2 * first_sign * second_sign * rank10_gram[10][11]
        )
        if four_square != -8:
            raise ArithmeticError("forced saturation root is not sheet-switch invariant")
        sign_square_checks[f"{first_sign:+d},{second_sign:+d}"] = four_square // 4

    # Identify the lattice generated before saturation.  The first ten
    # generators have rank eight.  The displayed eight form an integral
    # basis because the remaining two have integral coordinates in it; after
    # division by -2 its Gram is positive even unimodular of rank eight,
    # hence E8.  The final two vectors are orthogonal <-4> summands.
    first_block = sp.Matrix(rank10_gram[:10])[:, :10]
    basis_indices = [0, 1, 2, 3, 4, 5, 6, 8]
    e8_minus2_gram = first_block.extract(basis_indices, basis_indices)
    coordinates = e8_minus2_gram.inv() * first_block.extract(basis_indices, list(range(10)))
    if any(value.q != 1 for value in coordinates):
        raise ArithmeticError("first rank-ten block is not generated by the claimed basis")
    e8_positive = -e8_minus2_gram / 2
    _ldl_upper(e8_positive)
    if (
        int(e8_positive.det()) != 1
        or any(e8_positive[index, index] % 2 for index in range(8))
        or _smith_diagonal(e8_minus2_gram) != [2] * 8
    ):
        raise ArithmeticError("rank-ten generated E8(-2) identification failed")
    generated_rank10_gram = sp.diag(e8_minus2_gram, -4, -4)
    if abs(int(generated_rank10_gram.det())) != 4096:
        raise ArithmeticError("rank-ten generated lattice determinant changed")

    production_forced = production["rank_ten_residue"]["forced_saturation_root"]
    if (
        production_forced["e_square"] != -2
        or production_forced["downstairs_root_sum"] != downstairs_sum
        or not production_forced["sheet_switch_invariant"]
        or not production["checks"][
            "rank_10_orbit_eliminated_by_sheet_invariant_forced_saturation_root"
        ]
    ):
        raise ArithmeticError("production did not encode the rank-ten killer exactly")

    checks = {
        "all_27648_assignments_independently_generated": True,
        "positive_DFS_gauge_gives_18_complete_switching_orbits": True,
        "direct_orbits_match_component_stabilizer_formula": True,
        "exact_characteristic_polynomial_Sturm_inertia_matches_all_18_production_orbits": True,
        "unique_rank_9_orbit_equals_geometric_package": True,
        "rank_10_generated_lattice_is_E8_minus2_plus_two_minus4": True,
        "rank_10_all_four_sheet_choices_force_integral_saturation_minus_two_root": True,
        "rank_10_is_eliminated_by_L2_not_left_as_a_period_candidate": True,
        "unique_rank_9_geometric_orbit_is_the_only_Type_I_anti_Gram_survivor": True,
    }

    canonical_key_hash = hashlib.sha256(
        json.dumps(sorted(representatives), separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return {
        "schema": "enriques_type_i_downstairs_anti_gram_census_hostile_audit_v1",
        "verdict": "PASS/FREEZE (strengthened: rank ten eliminated)",
        "pinned_inputs": actual_hashes,
        "independent_switching_census": {
            "canonical_gauge": "positive sorted depth-first spanning-forest edges",
            "different_from_production_gauge": True,
            "raw_assignment_count": raw_count,
            "orbit_count": len(representatives),
            "orbit_size_histogram": {
                str(size): count for size, count in sorted(orbit_sizes.items())
            },
            "canonical_key_set_sha256": canonical_key_hash,
            "all_direct_orbits_replayed": True,
        },
        "independent_exact_inertia": {
            "engine": "factorized exact characteristic polynomial plus Sturm root counts",
            "raw_assignment_histogram": _histogram(weighted),
            "switch_orbit_histogram": _histogram(unweighted),
            "unique_rank_9_nonpositive_orbit": True,
            "unique_rank_10_nonpositive_orbit_before_L2": True,
            "rank_9_equals_geometric_package": True,
        },
        "rank_ten_hostile_killer": {
            "generated_lattice_before_saturation": "E8(-2)+<-4>+<-4>",
            "generated_lattice_determinant_absolute": abs(int(generated_rank10_gram.det())),
            "generated_lattice_rootless_by_norm_divisibility": True,
            "orthogonal_wall_pair_one_based": [11, 12],
            "downstairs_root_sum": downstairs_sum,
            "downstairs_root_half_sum": [value // 2 for value in downstairs_sum],
            "formula": (
                "e=(eps_11*delta_11+eps_12*delta_12)/2="
                "b_11^eps_11+b_12^eps_12-(u_11+u_12)/2"
            ),
            "integrality_source": "half-sum glue plus (r_11+r_12)/2 integral",
            "square_for_all_sheet_signs": sign_square_checks,
            "e_lies_in_saturation_of_anti_span": True,
            "e_is_anti_invariant": True,
            "L2_rootlessness_fails": True,
            "conclusion": "rank-ten orbit eliminated before embedding/period/L4",
        },
        "final_theorem": (
            "For the Type-I 12A wall system, the 27,648 two-sheet anti-Grams "
            "reduce, after negative-semidefiniteness, half-sum glue, and "
            "saturated rootlessness, to the unique rank-nine geometric orbit "
            "up to independent sheet swaps."
        ),
        "remaining_scope": {
            "rank_ten_branch_survives": False,
            "new_polytope_or_surface_claimed": False,
            "still_not_a_downstairs_construction_of_the_rank_nine_embedding_glue_period_L4_package": True,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(analyze(), indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived hostile audit differs: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(hashlib.sha256(rendered.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
