#!/usr/bin/env python3
"""Exhaust the Type-I two-sheet anti-Grams from downstairs data alone.

For two lifted Enriques curves over roots with downstairs intersection ``m``,
nonnegative same/cross intersections force

``delta_i.delta_j in {-2m,-2m+4,...,2m}``.

The Type-I graph has ten ``m=1`` pairs and three ``m=2`` pairs, hence 27,648
raw assignments.  This script constructs all of them, quotients exactly by
independent sheet swaps ``delta_i -> -delta_i``, computes exact inertia only
on canonical orbit representatives, and compares the unique rank-nine
negative-semidefinite orbit with the independently reconstructed geometric
lift package.  The sole rank-ten negative-semidefinite orbit is then
eliminated by the half-sum glue: it forces a root in the saturated
anti-invariant span.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Iterable

import sympy as sp

try:
    from scripts.certify_enriques_variable_lift_package import parse_gap_records
except ModuleNotFoundError:  # direct execution
    from certify_enriques_variable_lift_package import parse_gap_records
from coxeter_polytopes.exact.matrix import inertia


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
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_i_downstairs_anti_gram_census.json"
)

SHIMADA_SHA256 = (
    "218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c"
)
EXPLICIT_PACKAGE_SHA256 = (
    "d9313e37687e355c9ca8f86d1ee11d74f0810a8c6906a0443f5377d2481fece7"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _matrix_key(matrix: Iterable[Iterable[int]]) -> tuple[int, ...]:
    rows = [list(map(int, row)) for row in matrix]
    return tuple(
        rows[left][right]
        for left in range(len(rows))
        for right in range(left, len(rows))
    )


def _matrix_from_key(key: tuple[int, ...], size: int) -> list[list[int]]:
    matrix = [[0] * size for _ in range(size)]
    cursor = 0
    for left in range(size):
        for right in range(left, size):
            matrix[left][right] = matrix[right][left] = int(key[cursor])
            cursor += 1
    if cursor != len(key):
        raise ArithmeticError("bad upper-triangle key length")
    return matrix


def _support_components(matrix: list[list[int]]) -> list[list[int]]:
    size = len(matrix)
    adjacency = [
        [other for other in range(size) if other != vertex and matrix[vertex][other] != 0]
        for vertex in range(size)
    ]
    components = []
    unseen = set(range(size))
    while unseen:
        root = min(unseen)
        component = []
        queue = [root]
        unseen.remove(root)
        while queue:
            vertex = queue.pop(0)
            component.append(vertex)
            for other in sorted(adjacency[vertex]):
                if other in unseen:
                    unseen.remove(other)
                    queue.append(other)
        components.append(sorted(component))
    return components


def _switch_canonical(
    matrix: list[list[int]],
) -> tuple[tuple[int, ...], list[list[int]], list[int]]:
    """Canonicalize a signed support graph under diagonal sign switches.

    In each support component, sorted breadth-first tree edges are forced
    negative.  Cycle-edge signs are switching invariants, so the resulting
    matrix is a complete canonical representative for this fixed labelled
    support graph.
    """

    size = len(matrix)
    adjacency = [
        [other for other in range(size) if other != vertex and matrix[vertex][other] != 0]
        for vertex in range(size)
    ]
    signs: list[int | None] = [None] * size
    for root in range(size):
        if signs[root] is not None:
            continue
        signs[root] = 1
        queue = [root]
        while queue:
            vertex = queue.pop(0)
            for other in sorted(adjacency[vertex]):
                if signs[other] is None:
                    edge_sign = 1 if matrix[vertex][other] > 0 else -1
                    signs[other] = -int(signs[vertex]) * edge_sign
                    queue.append(other)
    integral_signs = [int(sign) for sign in signs]
    canonical = [
        [
            integral_signs[left] * matrix[left][right] * integral_signs[right]
            for right in range(size)
        ]
        for left in range(size)
    ]
    return _matrix_key(canonical), canonical, integral_signs


def _switch_orbit(matrix: list[list[int]]) -> set[tuple[int, ...]]:
    size = len(matrix)
    # A global sign is invisible, so fix the final vertex sign to +1.
    orbit = set()
    for initial in itertools.product((-1, 1), repeat=size - 1):
        signs = initial + (1,)
        switched = [
            [signs[left] * matrix[left][right] * signs[right] for right in range(size)]
            for left in range(size)
        ]
        orbit.add(_matrix_key(switched))
    return orbit


def _key_hash(keys: Iterable[tuple[int, ...]]) -> str:
    rendered = json.dumps(sorted(keys), separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(rendered).hexdigest()


def analyze() -> dict[str, Any]:
    if _sha256(SHIMADA_DATA) != SHIMADA_SHA256:
        raise AssertionError("Shimada input changed")
    if _sha256(EXPLICIT_PACKAGE) != EXPLICIT_PACKAGE_SHA256:
        raise AssertionError("explicit Type-I package changed")

    record = parse_gap_records(SHIMADA_DATA)["12A"]
    configuration = sp.Matrix(record["configuration"])
    if configuration.shape != (12, 12):
        raise ArithmeticError("Type-I configuration is not 12 by 12")
    m1_edges = [
        (left, right)
        for left in range(12)
        for right in range(left + 1, 12)
        if configuration[left, right] == 1
    ]
    m2_edges = [
        (left, right)
        for left in range(12)
        for right in range(left + 1, 12)
        if configuration[left, right] == 2
    ]
    if len(m1_edges) != 10 or len(m2_edges) != 3:
        raise ArithmeticError("Type-I nonorthogonal edge census changed")
    if any(
        configuration[left, right] not in (0, 1, 2)
        for left in range(12)
        for right in range(left + 1, 12)
    ):
        raise ArithmeticError("unexpected downstairs intersection")

    orbit_multiplicity: Counter[tuple[int, ...]] = Counter()
    orbit_representatives: dict[tuple[int, ...], list[list[int]]] = {}
    raw_count = 0
    for m1_values in itertools.product((-2, 2), repeat=len(m1_edges)):
        for m2_values in itertools.product((-4, 0, 4), repeat=len(m2_edges)):
            matrix = [[-4 if left == right else 0 for right in range(12)] for left in range(12)]
            for (left, right), value in zip(m1_edges, m1_values):
                matrix[left][right] = matrix[right][left] = value
            for (left, right), value in zip(m2_edges, m2_values):
                matrix[left][right] = matrix[right][left] = value
            key, representative, _ = _switch_canonical(matrix)
            orbit_multiplicity[key] += 1
            orbit_representatives[key] = representative
            raw_count += 1
    if raw_count != 2**10 * 3**3 or sum(orbit_multiplicity.values()) != raw_count:
        raise ArithmeticError("raw anti-Gram denominator changed")
    if len(orbit_multiplicity) != 18:
        raise ArithmeticError("sheet-switch orbit count changed")

    weighted_inertia: Counter[tuple[int, int, int]] = Counter()
    orbit_inertia: Counter[tuple[int, int, int]] = Counter()
    orbit_records = []
    for key in sorted(orbit_representatives):
        representative = orbit_representatives[key]
        exact_inertia = inertia(representative)
        signature = (
            exact_inertia.positive,
            exact_inertia.negative,
            exact_inertia.zero,
        )
        multiplicity = orbit_multiplicity[key]
        components = _support_components(representative)
        predicted_orbit_size = 2 ** (12 - len(components))
        direct_orbit = _switch_orbit(representative)
        if len(direct_orbit) != predicted_orbit_size or multiplicity != predicted_orbit_size:
            raise ArithmeticError("canonical switching multiplicity is not a full orbit")
        if {_switch_canonical(_matrix_from_key(item, 12))[0] for item in direct_orbit} != {key}:
            raise ArithmeticError("switching orbit does not return to one canonical key")
        weighted_inertia[signature] += multiplicity
        orbit_inertia[signature] += 1
        orbit_records.append(
            {
                "canonical_gram": representative,
                "canonical_key_sha256": _key_hash([key]),
                "raw_assignment_multiplicity": multiplicity,
                "support_components_one_based": [
                    [vertex + 1 for vertex in component] for component in components
                ],
                "support_component_count": len(components),
                "switch_orbit_size_formula": f"2^(12-{len(components)})",
                "inertia": asdict(exact_inertia),
                "rank": exact_inertia.positive + exact_inertia.negative,
            }
        )

    expected_weighted = Counter(
        {
            (2, 10, 0): 14336,
            (1, 9, 2): 2048,
            (1, 10, 1): 9216,
            (1, 11, 0): 512,
            (0, 9, 3): 1024,
            (0, 10, 2): 512,
        }
    )
    expected_orbits = Counter(
        {
            (2, 10, 0): 7,
            (1, 9, 2): 1,
            (1, 10, 1): 7,
            (1, 11, 0): 1,
            (0, 9, 3): 1,
            (0, 10, 2): 1,
        }
    )
    if weighted_inertia != expected_weighted or orbit_inertia != expected_orbits:
        raise ArithmeticError("exact anti-Gram inertia census changed")

    nonpositive = [
        record
        for record in orbit_records
        if record["inertia"]["positive"] == 0
    ]
    rank_nine = [record for record in nonpositive if record["rank"] == 9]
    rank_ten = [record for record in nonpositive if record["rank"] == 10]
    if len(nonpositive) != 2 or len(rank_nine) != 1 or len(rank_ten) != 1:
        raise ArithmeticError("nonpositive orbit uniqueness changed")

    # The rank-ten residue is not a good-package candidate.  In its fixed
    # wall labelling delta_11 and delta_12 are orthogonal square-minus-four
    # vectors, while the corresponding downstairs roots have even sum.  In
    # any half-sum-glued realization and for any sheet signs eps_i, put
    #
    #   b_i^eps=(u_i+eps_i delta_i)/2,
    #   e=(eps_11 delta_11+eps_12 delta_12)/2
    #    =b_11^eps_11+b_12^eps_12-(u_11+u_12)/2.
    #
    # The last expression proves integrality, and orthogonality gives
    # e^2=(-4-4)/4=-2.  Thus e belongs to the saturation of the anti span
    # and violates L2.  This argument is invariant under every sheet swap.
    rank_ten_gram = rank_ten[0]["canonical_gram"]
    forced_pair = (10, 11)
    if not (
        rank_ten_gram[10][10] == -4
        and rank_ten_gram[11][11] == -4
        and rank_ten_gram[10][11] == 0
    ):
        raise ArithmeticError("rank-ten forced-root pair changed")
    walls = sp.Matrix(record["walls"])
    forced_downstairs_sum = [
        int(walls[10, column] + walls[11, column]) for column in range(10)
    ]
    if any(value % 2 for value in forced_downstairs_sum):
        raise ArithmeticError("rank-ten forced-root downstairs sum is not even")
    forced_root_square_numerator = (
        rank_ten_gram[10][10]
        + rank_ten_gram[11][11]
        + 2 * rank_ten_gram[10][11]
    )
    if forced_root_square_numerator != -8:
        raise ArithmeticError("rank-ten forced-root norm numerator changed")

    explicit = json.loads(EXPLICIT_PACKAGE.read_text(encoding="utf-8"))
    geometric_gram = explicit["finite_lift_compiler_calibration"]["frozen_anti_difference_gram"]
    geometric_key, geometric_canonical, geometric_signs = _switch_canonical(geometric_gram)
    if geometric_key != _matrix_key(rank_nine[0]["canonical_gram"]):
        raise ArithmeticError("unique rank-nine orbit is not the geometric package")
    geometric_orbit = _switch_orbit(geometric_gram)
    if len(geometric_orbit) != 1024 or geometric_orbit != {
        item
        for item in geometric_orbit
        if _switch_canonical(_matrix_from_key(item, 12))[0] == geometric_key
    }:
        raise ArithmeticError("geometric sheet-switch orbit is inconsistent")

    checks = {
        "all_27648_anti_grams_generated_from_downstairs_alphabet": True,
        "all_raw_assignments_partitioned_into_18_complete_sheet_switch_orbits": True,
        "exact_inertia_computed_on_every_orbit_representative": True,
        "unique_negative_semidefinite_rank_9_orbit": True,
        "rank_9_orbit_equals_explicit_geometric_package": True,
        "unique_negative_semidefinite_rank_10_orbit": True,
        "rank_10_orbit_eliminated_by_sheet_invariant_forced_saturation_root": True,
    }

    def histogram(counter: Counter[tuple[int, int, int]]) -> dict[str, int]:
        return {
            f"p{positive}_n{negative}_z{zero}": count
            for (positive, negative, zero), count in sorted(counter.items())
        }

    return {
        "schema": "enriques_type_i_downstairs_anti_gram_census_v1",
        "scope": {
            "input": "only Shimada Type-I downstairs intersections and curve nonnegativity alphabet",
            "result": (
                "the rank-nine geometric anti-Gram is the only Type-I orbit "
                "surviving negative-semidefiniteness and saturated rootlessness"
            ),
            "rank_ten_status": (
                "eliminated before period/embedding/L4: half-sum glue forces "
                "a minus-two vector in the saturated anti span"
            ),
            "remaining_rank_nine_compiler_work": [
                "derive the compatible primitive embedding and half-sum glue orbit downstairs",
                "recover the exact period and ample orientation downstairs",
                "replay L4 and cone completeness without importing the OSCAR carrier",
            ],
            "new_abstract_polytope_claimed": False,
        },
        "pinned_inputs": {
            "shimada_path": str(SHIMADA_DATA.relative_to(ROOT)),
            "shimada_sha256": SHIMADA_SHA256,
            "explicit_package_path": str(EXPLICIT_PACKAGE.relative_to(ROOT)),
            "explicit_package_sha256": EXPLICIT_PACKAGE_SHA256,
        },
        "downstairs_alphabet": {
            "formula": "delta_i.delta_j=2(a-c)=4a-2m with a,c>=0 and a+c=m",
            "m1_edges_zero_based": [list(edge) for edge in m1_edges],
            "m2_edges_zero_based": [list(edge) for edge in m2_edges],
            "m1_values": [-2, 2],
            "m2_values": [-4, 0, 4],
            "raw_assignment_count": raw_count,
        },
        "sheet_switch_quotient": {
            "action": "delta_i -> epsilon_i delta_i, epsilon_i in {+1,-1}",
            "canonicalization": (
                "force sorted breadth-first spanning-tree edges negative in every support component"
            ),
            "orbit_count": len(orbit_records),
            "orbit_size_histogram": {
                str(size): count
                for size, count in sorted(Counter(orbit_multiplicity.values()).items())
            },
            "all_canonical_keys_sha256": _key_hash(orbit_representatives),
        },
        "exact_inertia_census": {
            "raw_assignment_histogram": histogram(weighted_inertia),
            "sheet_switch_orbit_histogram": histogram(orbit_inertia),
            "orbit_records": orbit_records,
        },
        "generic_rank_nine_theorem": {
            "unique_orbit": True,
            "raw_assignment_count": rank_nine[0]["raw_assignment_multiplicity"],
            "canonical_gram": rank_nine[0]["canonical_gram"],
            "canonical_gram_inertia": rank_nine[0]["inertia"],
            "equals_explicit_Type_I_package_up_to_sheet_swaps": True,
            "geometric_switching_signs_to_canonical": geometric_signs,
            "geometric_canonical_gram": geometric_canonical,
            "interpretation": (
                "for Picard rank 19 (anti rank 9), the downstairs curve-intersection "
                "alphabet already forces the geometric anti-Gram up to relabelling the two sheets"
            ),
        },
        "rank_ten_residue": {
            "unique_signature_candidate_orbit": True,
            "raw_assignment_count": rank_ten[0]["raw_assignment_multiplicity"],
            "canonical_gram": rank_ten[0]["canonical_gram"],
            "canonical_gram_inertia": rank_ten[0]["inertia"],
            "forced_saturation_root": {
                "wall_pair_one_based": [forced_pair[0] + 1, forced_pair[1] + 1],
                "downstairs_root_sum": forced_downstairs_sum,
                "downstairs_half_sum_integral": True,
                "formula_for_arbitrary_sheet_signs": (
                    "e=(eps_11*delta_11+eps_12*delta_12)/2="
                    "b_11^eps_11+b_12^eps_12-(u_11+u_12)/2"
                ),
                "four_times_e_square": forced_root_square_numerator,
                "e_square": forced_root_square_numerator // 4,
                "e_is_integral_anti_invariant_and_in_saturation": True,
                "sheet_switch_invariant": True,
            },
            "status": "ELIMINATED: forced saturated anti-invariant minus-two root violates L2",
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
            raise SystemExit(f"archived certificate differs: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(hashlib.sha256(rendered.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
