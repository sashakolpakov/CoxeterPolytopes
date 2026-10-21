#!/usr/bin/env python3
"""Independent hostile audit of the source-29 triple-123 cut.

This deliberately imports no project mathematics.  It reconstructs the
relevant source ground from the two missing edges and two missing four-sets,
derives the link g-vector from face counts, uses Burnside's lemma for the
configuration-orbit counts, and checks the frozen near seed directly from
its integer masks.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from math import comb
from typing import Iterable, Iterator, Mapping


VERTICES = tuple(range(14))
T = frozenset((1, 2, 3))
Q = tuple(range(4, 13))
FREE_FIBRE = frozenset((5, 6, 7))
REMAINDER = tuple(range(8, 13))
MISSING_EDGES = (frozenset((0, 13)), frozenset((1, 13)))
MISSING_FOURS = (
    frozenset((0, 1, 2, 3)),
    frozenset((0, 1, 2, 4)),
)
ACTIVE_PATTERNS = (
    (227,),
    (229,),
    (227, 229),
    (229, 230),
    (227, 229, 230),
)
NEAR_SEED = (
    3296, 7328, 6337, 1578, 2094, 4944, 1872, 6277, 6157, 302,
    15488, 595, 6188, 7312, 789, 3114, 12332, 9040, 856, 342, 3176,
)
EXPECTED_AUDIT_DIGEST = "97d9990e9ae3f4837aca1f5fd9c4e4becb689ad78abda8fff8b53bf4771af5ab"


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _decode(value: int) -> frozenset[int]:
    return frozenset(vertex for vertex in VERTICES if value >> vertex & 1)


def _allowed(face: frozenset[int]) -> bool:
    return not any(row <= face for row in (*MISSING_EDGES, *MISSING_FOURS))


def _macaulay_raise(value: int, degree: int) -> int:
    """Compute the Macaulay successor without importing the project helper."""

    if value < 0 or degree < 1:
        raise ValueError("Macaulay arguments must satisfy value >= 0, degree >= 1")
    remainder = value
    previous_top: int | None = None
    successor = 0
    for lower in range(degree, 0, -1):
        top = lower
        ceiling = previous_top - 1 if previous_top is not None else None
        while (ceiling is None or top + 1 <= ceiling) and comb(top + 1, lower) <= remainder:
            top += 1
        if comb(top, lower) <= remainder:
            remainder -= comb(top, lower)
            successor += comb(top + 1, lower + 1)
            previous_top = top
        else:
            previous_top = top
    if remainder:
        raise AssertionError("failed to construct the Macaulay representation")
    return successor


def _base_link() -> dict[str, object]:
    # T has a six-dimensional simplicial link.  Its j-faces are the allowed
    # j-vertex extensions of T before the additional missing A5s are chosen.
    extension_counts = [1]
    outside = tuple(vertex for vertex in VERTICES if vertex not in T)
    for order in (1, 2, 3):
        extension_counts.append(
            sum(
                _allowed(T | frozenset(extension))
                for extension in combinations(outside, order)
            )
        )
    link_dimension = 6
    h = []
    for index in range(4):
        h.append(
            sum(
                (-1) ** (index - order)
                * comb(link_dimension - order, link_dimension - index)
                * extension_counts[order]
                for order in range(index + 1)
            )
        )
    g = tuple(h[index] - h[index - 1] for index in range(1, 4))
    return {
        "extension_counts_f_minus_one_through_f2": extension_counts,
        "h0_through_h3": h,
        "g1_through_g3": list(g),
    }


def _ground_and_fibres() -> dict[str, object]:
    candidates = tuple(
        frozenset(row)
        for row in combinations(VERTICES, 5)
        if _allowed(frozenset(row))
    )
    containing = tuple(row for row in candidates if T <= row)
    pairs = tuple(row - T for row in containing)
    fibres = tuple(T | frozenset(row) for row in combinations(Q, 3))
    killers = {
        fibre: tuple(
            row for row in candidates if not T <= row and row < fibre
        )
        for fibre in fibres
    }
    killer_union = set().union(*map(set, killers.values()))

    coverage = Counter()
    intersection = Counter()
    for degree in (0, 1, 2):
        for selected in combinations(pairs, degree):
            covered = {
                fibre
                for fibre in fibres
                if any(pair <= fibre - T for pair in selected)
            }
            coverage[(degree, len(covered))] += 1
            if degree == 2:
                intersection[
                    "intersecting" if selected[0] & selected[1] else "disjoint"
                ] += 1

    return {
        "candidate_A5_count": len(candidates),
        "allowed_six_fibre_count": len(fibres),
        "T_containing_candidate_A5_count": len(containing),
        "T_containing_pair_bijection": set(pairs)
        == {frozenset(row) for row in combinations(Q, 2)},
        "killers_per_fibre_histogram": dict(sorted(Counter(map(len, killers.values())).items())),
        "noncontaining_killer_union_count": len(killer_union),
        "coverage_census": [
            {"d": key[0], "base_coverage": key[1], "configuration_count": value}
            for key, value in sorted(coverage.items())
        ],
        "pair_pair_census": dict(sorted(intersection.items())),
    }


def _strict_budget() -> list[dict[str, int | str]]:
    base_g = (2, 3, 4)
    coverages = ((0, "none", 0), (1, "one pair", 7),
                 (2, "intersecting pairs", 13), (2, "disjoint pairs", 14))
    rows: list[dict[str, int | str]] = []
    for degree, kind, base_coverage in coverages:
        final_g2 = base_g[1] - degree
        constant = base_g[2] + 5 * degree
        minimum = constant - _macaulay_raise(final_g2, 2)
        maximum = constant
        rows.append(
            {
                "d": degree,
                "configuration": kind,
                "final_g2": final_g2,
                "allowed_missing_six_min": minimum,
                "allowed_missing_six_max": maximum,
                "base_coverage": base_coverage,
                "outside_fibre_budget": maximum - base_coverage,
            }
        )
    return rows


def _block_group(*blocks: tuple[int, ...]) -> Iterator[dict[int, int]]:
    permutation_rows: Iterable[tuple[tuple[int, ...], ...]] = product(
        *(permutations(block) for block in blocks)
    )
    for rows in permutation_rows:
        action: dict[int, int] = {}
        for block, row in zip(blocks, rows, strict=True):
            action.update(zip(block, row, strict=True))
        yield action


EDGES = tuple(combinations(Q, 2))
EDGE_INDEX = {frozenset(edge): index for index, edge in enumerate(EDGES)}


def _edge_action(action: Mapping[int, int]) -> tuple[int, ...]:
    return tuple(
        EDGE_INDEX[frozenset((action.get(left, left), action.get(right, right)))]
        for left, right in EDGES
    )


def _burnside_edge_subset_orbits(*blocks: tuple[int, ...]) -> tuple[int, int]:
    fixed_one_sum = fixed_two_sum = group_order = 0
    for action in _block_group(*blocks):
        group_order += 1
        moved = _edge_action(action)
        fixed_edges = sum(index == target for index, target in enumerate(moved))
        edge_two_cycles = sum(
            index < target and moved[target] == index
            for index, target in enumerate(moved)
        )
        fixed_one_sum += fixed_edges
        fixed_two_sum += comb(fixed_edges, 2) + edge_two_cycles
    if fixed_one_sum % group_order or fixed_two_sum % group_order:
        raise AssertionError("Burnside average was not integral")
    return fixed_one_sum // group_order, fixed_two_sum // group_order


def _active_two_orbits() -> tuple[int, Counter[tuple[int, bool]]]:
    group = tuple(_block_group(tuple(FREE_FIBRE), REMAINDER))

    def moved(configuration: tuple[int, int], action: Mapping[int, int]) -> tuple[int, int]:
        edge_action = _edge_action(action)
        return tuple(sorted(edge_action[index] for index in configuration))  # type: ignore[return-value]

    configurations = tuple(combinations(range(len(EDGES)), 2))
    representatives = {
        min(moved(configuration, action) for action in group)
        for configuration in configurations
    }
    profile: Counter[tuple[int, bool]] = Counter()
    for left_index, right_index in representatives:
        left = frozenset(EDGES[left_index])
        right = frozenset(EDGES[right_index])
        profile[(len(left & right), left <= FREE_FIBRE or right <= FREE_FIBRE)] += 1
    return len(representatives), profile


def _symmetry_and_branches() -> dict[str, object]:
    zero_d1, zero_d2 = _burnside_edge_subset_orbits(tuple(range(5, 13)))
    active_d1, active_d2_burnside = _burnside_edge_subset_orbits(
        tuple(FREE_FIBRE), REMAINDER
    )
    active_d2_direct, profile = _active_two_orbits()
    if active_d2_burnside != active_d2_direct:
        raise AssertionError("independent Burnside and canonical-image counts disagree")

    decoded_patterns = tuple(tuple(_decode(value) for value in row) for row in ACTIVE_PATTERNS)
    active_group = tuple(_block_group(tuple(FREE_FIBRE), REMAINDER))
    if any(
        {
            frozenset(action.get(vertex, vertex) for vertex in block)
            for block in pattern
        }
        != set(pattern)
        for pattern in decoded_patterns
        for action in active_group
    ):
        raise AssertionError("the claimed active T-stabilizer moves a fixed mask pattern")

    fixed_fibre = T | FREE_FIBRE
    fixed_mask_kills = {
        value: (_decode(value) < fixed_fibre and not T <= _decode(value))
        for value in (227, 229, 230)
    }
    excluded = profile[(0, False)]
    branch_rows = [
        {
            "triple012_branch": "zero_support",
            "T_stabilizer_on_Q": "S8 on 5..12; vertex 4 fixed",
            "orbits_by_d": {"0": 1, "1": zero_d1, "2": zero_d2},
            "configuration_branch_count": 1 + zero_d1 + zero_d2,
        }
    ]
    for pattern in ACTIVE_PATTERNS:
        has_fixed_killer = 230 in pattern
        d2 = active_d2_direct - excluded if has_fixed_killer else active_d2_direct
        branch_rows.append(
            {
                "triple012_branch": list(pattern),
                "T_stabilizer_on_Q": "S3 on 567 times S5 on 8..12; vertex 4 fixed",
                "fixed_123_killer_on_567": has_fixed_killer,
                "orbits_by_d": {"0": 1, "1": active_d1, "2": d2},
                "excluded_disjoint_noncovering_d2_orbits": excluded if has_fixed_killer else 0,
                "configuration_branch_count": 1 + active_d1 + d2,
            }
        )

    return {
        "zero_support_orbits_by_d": {"d0": 1, "d1": zero_d1, "d2": zero_d2},
        "active_support_orbits_by_d": {"d0": 1, "d1": active_d1, "d2": active_d2_direct},
        "active_d2_profile": [
            {
                "pair_intersection_size": key[0],
                "covers_567": key[1],
                "orbit_count": value,
            }
            for key, value in sorted(profile.items())
        ],
        "fixed_mask_kills_123567": {str(key): value for key, value in fixed_mask_kills.items()},
        "mask230_excluded_disjoint_orbits_per_affected_branch": excluded,
        "branch_rows": branch_rows,
        "combined_configuration_branch_count": sum(
            int(row["configuration_branch_count"]) for row in branch_rows
        ),
        "residual_stabilizer_check": (
            "After fixing T=123, the S_{0,1} and S_{3,4} factors cannot act. "
            "Thus zero support leaves S8 on 5..12, while normalized active "
            "support leaves S3(567)xS5(8..12); each active mask pattern is "
            "invariant under the latter action."
        ),
    }


def _near_seed() -> dict[str, object]:
    decoded = tuple((value, _decode(value)) for value in NEAR_SEED)
    containing = tuple((value, row - T) for value, row in decoded if T <= row)
    all_outside_fibres = tuple(frozenset(row) for row in combinations(Q, 3))
    base_covered = {
        outside
        for outside in all_outside_fibres
        if any(pair <= outside for _, pair in containing)
    }
    new_killers: list[tuple[int, frozenset[int]]] = []
    for value, row in decoded:
        if T <= row:
            continue
        for outside in all_outside_fibres:
            if row < T | outside and outside not in base_covered:
                new_killers.append((value, outside))
    new_supports = {outside for _, outside in new_killers}
    missing = base_covered | new_supports
    final_g2 = 3 - len(containing)
    final_g3 = 4 + 5 * len(containing) - len(missing)
    return {
        "T_containing_masks": [value for value, _ in containing],
        "outside_pairs": [sorted(pair) for _, pair in containing],
        "base_covered_fibre_count": len(base_covered),
        "new_outside_killers": [
            {"mask": value, "outside_fibre": sorted(outside)}
            for value, outside in sorted(new_killers)
        ],
        "new_outside_fibre_count": len(new_supports),
        "actual_missing_six_fibre_count": len(missing),
        "final_g2": final_g2,
        "final_g3": final_g3,
        "allowed_g3_interval": [0, _macaulay_raise(final_g2, 2)],
        "negative_g3_defect": max(0, -final_g3),
    }


def audit() -> dict[str, object]:
    result = {
        "base_link": _base_link(),
        "ground_and_fibres": _ground_and_fibres(),
        "strict_budget": _strict_budget(),
        "symmetry_and_branches": _symmetry_and_branches(),
        "near_seed": _near_seed(),
        "scope": {
            "source29_status": "OPEN",
            "enumerated_objects": "pair-configuration orbits, not completed strict-link survivors",
            "configuration_rows_are_survivors": False,
            "outside_killer_supports_still_unresolved": True,
            "compact_incidence_or_hyperbolic_witness_produced": False,
        },
    }
    digest = _digest(result)
    if digest != EXPECTED_AUDIT_DIGEST:
        raise AssertionError(f"hostile-audit transcript changed: {digest}")
    return {**result, "audit_digest": digest}


def main() -> None:
    print(json.dumps(audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
