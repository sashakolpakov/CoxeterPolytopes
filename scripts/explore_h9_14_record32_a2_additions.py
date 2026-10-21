#!/usr/bin/env python3
"""Exact additional-A5 census for the six viable record-32 ``a=2`` orbits.

The six representatives are produced by the lossless local/forced-shadow
scout in ``explore_h9_14_record32_all_completions.py``.  This script exhausts
every additional minimal A5 choice allowed by the exact vertex caps.  It then
applies the exact order-at-most-four transversal gate and counts rank-nine
facet complements by intersecting precomputed five-cover hit words.

This is an incidence scout, not a topology or Coxeter theorem.  Scheduled
targets still require literal blocker dualization, ridge/full-link audits,
and unrestricted Lanner/actual-vertex gates.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json
from math import comb

from coxeter_polytopes.combinatorics.h9_14_record32 import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
)
from coxeter_polytopes.combinatorics.h9_14_record35_all_completions import (
    _forced_five_shadow,
    _shadow_gate,
)
from coxeter_polytopes.combinatorics.h9_14_record46 import _base_link_g3_g4
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
)


VERTEX_COUNT = 14
EDGE_MASKS = tuple(
    sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES
)
TOTAL_FIVES = 18
SCHEDULED_FACETS = 110
A2_REPRESENTATIVES = (
    ((2, 3, 4, 11), (5, 6, 7, 8)),
    ((2, 3, 4, 5), (3, 4, 5, 11)),
    ((1, 2, 3, 4), (2, 3, 4, 11)),
    ((1, 2, 3, 4), (3, 4, 11, 12)),
    ((1, 2, 3, 4), (3, 5, 11, 12)),
    ((1, 2, 3, 4), (5, 6, 11, 12)),
)


def mask(values) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1)


def small_base_transversals(hyperedges: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        candidate
        for order in range(1, 5)
        for chosen in combinations(range(VERTEX_COUNT), order)
        for candidate in (mask(chosen),)
        if all(candidate & edge for edge in hyperedges)
    )


def hit_word(face: int, candidates: tuple[int, ...]) -> int:
    return sum(bool(face & candidate) << index for index, candidate in enumerate(candidates))


def at_least_masks(additions: tuple[int, ...]) -> tuple[int, ...]:
    if not additions:
        return ()
    if len(additions) == 1:
        return additions
    if len(additions) == 2:
        first, second = additions
        return first | second, first & second
    if len(additions) == 3:
        first, second, third = additions
        return (
            first | second | third,
            (first & second) | (first & third) | (second & third),
            first & second & third,
        )
    raise ValueError("record-32 a=2 additional count escaped 0..3")


def audit_representative(index: int, representative) -> dict[str, object]:
    family = tuple(mask(face) for face in representative)
    fives = independent_masks(TARGET_GRAPH6, 5)
    base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)
    shadow = _forced_five_shadow(family, fives)
    strict, aggregate, additional, lower, upper = _shadow_gate(
        family, shadow, base_g3, base_g4, TOTAL_FIVES
    )
    if not strict or not aggregate or additional not in range(4):
        raise AssertionError("record-32 a=2 representative gate changed")
    candidates = tuple(
        five
        for five in fives
        if five not in shadow
        and not any(five & four == four for four in family)
        and all(
            not (five & (1 << vertex)) or upper[vertex] >= 1
            for vertex in range(VERTEX_COUNT)
        )
    )
    lower_at_least = tuple(
        mask(
            vertex
            for vertex in range(VERTEX_COUNT)
            if lower[vertex] >= threshold
        )
        for threshold in range(1, additional + 1)
    )
    upper_below = tuple(
        mask(
            vertex
            for vertex in range(VERTEX_COUNT)
            if upper[vertex] < threshold
        )
        for threshold in range(1, additional + 1)
    )
    small = small_base_transversals((*EDGE_MASKS, *family, *shadow))
    small_words = {five: hit_word(five, small) for five in candidates}
    five_covers = graph_vertex_cover_masks(TARGET_GRAPH6, 5)
    five_words = {face: hit_word(face, five_covers) for face in (*family, *candidates)}
    fixed_five_word = (1 << len(five_covers)) - 1
    for four in family:
        fixed_five_word &= five_words[four]

    cap_feasible = 0
    no_small = 0
    facet_histogram: Counter[int] = Counter()
    targets = []
    choices = ((),) if additional == 0 else combinations(candidates, additional)
    for additions in choices:
        additions = tuple(additions)
        present = at_least_masks(additions)
        if any(
            required & actual != required or forbidden & actual
            for required, forbidden, actual in zip(
                lower_at_least, upper_below, present, strict=True
            )
        ):
            continue
        cap_feasible += 1
        if additions:
            small_word = (1 << len(small)) - 1
            for face in additions:
                small_word &= small_words[face]
            if small_word:
                continue
        elif small:
            continue
        no_small += 1
        word = fixed_five_word
        for face in additions:
            word &= five_words[face]
        facet_count = word.bit_count()
        facet_histogram[facet_count] += 1
        if facet_count == SCHEDULED_FACETS:
            targets.append(additions)
    return {
        "orbit_representative_index": index,
        "missing_fours": [list(face) for face in representative],
        "forced_five_shadow_count": len(shadow),
        "additional_minimal_A5_count": additional,
        "individual_candidate_count": len(candidates),
        "raw_unordered_addition_count": comb(len(candidates), additional),
        "vertex_cap_feasible_count": cap_feasible,
        "no_transversal_of_order_at_most_four_count": no_small,
        "facet_count_histogram": dict(sorted(facet_histogram.items())),
        "scheduled_target_count": len(targets),
        "scheduled_target_addition_masks": [list(row) for row in targets],
        "scheduled_target_additions": [
            [list(vertices(face)) for face in row] for row in targets
        ],
    }


def main() -> None:
    audits = tuple(
        audit_representative(index, representative)
        for index, representative in enumerate(A2_REPRESENTATIVES)
    )
    print(
        json.dumps(
            {
                "kind": "h9_14_record32_a2_additional_A5_scout",
                "theorem_certificate": False,
                "orbit_representative_count": len(audits),
                "scheduled_target_count_across_representatives": sum(
                    audit["scheduled_target_count"] for audit in audits
                ),
                "representative_audits": audits,
                "scope": (
                    "Exact for additional A5 choices on the six strict-shadow "
                    "orbit representatives. Scheduled targets still need blocker, "
                    "topology, orbit, and unrestricted Coxeter-label audits."
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
