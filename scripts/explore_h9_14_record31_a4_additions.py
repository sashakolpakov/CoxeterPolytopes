#!/usr/bin/env python3
"""Exact low-addition census for Record 31's twenty compact ``a=4`` orbits.

The three representatives requiring three additions are emitted for the
separate pivoted C++ search.  The remaining seventeen representatives are
exhausted literally here.  This is a compact-only production scout; target
topology and unrestricted labels are downstream gates.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json
from math import comb
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record35_all_completions import (
    _forced_five_shadow,
    _shadow_gate,
)
from coxeter_polytopes.combinatorics.h9_14_record46 import _base_link_g3_g4
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_record31_positive_a5.json"
GRAPH6 = "M?????????O?O?o??"
VERTEX_COUNT = 14
TOTAL_FIVES = 33
SCHEDULED_FACETS = 114
MISSING_EDGES = ((0, 11), (1, 12), (0, 13), (1, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in MISSING_EDGES)


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
    return sum(
        bool(face & candidate) << index
        for index, candidate in enumerate(candidates)
    )


def at_least_masks(additions: tuple[int, ...]) -> tuple[int, ...]:
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
    raise ValueError("Record 31 a=4 additions escaped 1..3")


def representatives() -> tuple[tuple[int, ...], ...]:
    document = json.loads(SOURCE.read_text(encoding="utf-8"))
    branch = document["production"]["branches"][3]
    if branch["a"] != 4 or branch["orbit_count"] != 20:
        raise AssertionError("Record 31 a=4 orbit archive changed")
    return tuple(
        tuple(mask(face) for face in orbit["representative"])
        for orbit in branch["orbits"]
    )


def audit_representative(index: int, family: tuple[int, ...]) -> dict[str, object]:
    fives = independent_masks(GRAPH6, 5)
    base_g3, base_g4 = _base_link_g3_g4(GRAPH6)
    shadow = _forced_five_shadow(family, fives)
    strict, aggregate, additional, lower, upper = _shadow_gate(
        family, shadow, base_g3, base_g4, TOTAL_FIVES
    )
    if not strict or not aggregate or additional not in (1, 2, 3):
        raise AssertionError("Record 31 a=4 strict shadow gate changed")
    candidates = tuple(
        five
        for five in fives
        if five not in shadow
        and all(
            not (five & (1 << vertex)) or upper[vertex] >= 1
            for vertex in range(VERTEX_COUNT)
        )
    )
    base = (*EDGE_MASKS, *family, *shadow)
    small = small_base_transversals(base)
    result: dict[str, object] = {
        "orbit_representative_index": index,
        "missing_four_masks": list(family),
        "missing_fours": [list(vertices(face)) for face in family],
        "forced_five_shadow_count": len(shadow),
        "additional_minimal_A5_count": additional,
        "individual_candidate_count": len(candidates),
        "raw_unordered_addition_count": comb(len(candidates), additional),
        "addition_lower": list(lower),
        "addition_upper": list(upper),
        "base_small_transversal_count": len(small),
    }
    if additional == 3:
        result["deferred_to_pivoted_cpp"] = True
        return result

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
    small_words = {five: hit_word(five, small) for five in candidates}
    five_covers = graph_vertex_cover_masks(GRAPH6, 5)
    five_words = {
        face: hit_word(face, five_covers) for face in (*family, *candidates)
    }
    fixed_five_word = (1 << len(five_covers)) - 1
    for four in family:
        fixed_five_word &= five_words[four]

    cap_feasible = 0
    no_small = 0
    facet_histogram: Counter[int] = Counter()
    targets: list[tuple[int, ...]] = []
    for additions in combinations(candidates, additional):
        present = at_least_masks(additions)
        if any(
            required & actual != required or forbidden & actual
            for required, forbidden, actual in zip(
                lower_at_least, upper_below, present, strict=True
            )
        ):
            continue
        cap_feasible += 1
        small_word = (1 << len(small)) - 1
        for face in additions:
            small_word &= small_words[face]
        if small_word:
            continue
        no_small += 1
        word = fixed_five_word
        for face in additions:
            word &= five_words[face]
        facet_count = word.bit_count()
        facet_histogram[facet_count] += 1
        if facet_count == SCHEDULED_FACETS:
            targets.append(tuple(additions))
    result.update(
        {
            "deferred_to_pivoted_cpp": False,
            "vertex_cap_feasible_count": cap_feasible,
            "no_transversal_of_order_at_most_four_count": no_small,
            "facet_count_histogram": dict(sorted(facet_histogram.items())),
            "scheduled_target_count": len(targets),
            "scheduled_target_addition_masks": [list(row) for row in targets],
            "scheduled_target_additions": [
                [list(vertices(face)) for face in row] for row in targets
            ],
        }
    )
    return result


def main() -> None:
    audits = tuple(
        audit_representative(index, family)
        for index, family in enumerate(representatives())
    )
    print(
        json.dumps(
            {
                "kind": "h9_14_record31_a4_low_addition_scout",
                "compact_only": True,
                "orbit_representative_count": len(audits),
                "literal_low_addition_target_count": sum(
                    audit.get("scheduled_target_count", 0) for audit in audits
                ),
                "deferred_three_addition_indices": [
                    audit["orbit_representative_index"]
                    for audit in audits
                    if audit["deferred_to_pivoted_cpp"]
                ],
                "representative_audits": audits,
                "scope": (
                    "Exact literal addition coverage for every a=4 representative "
                    "requiring one or two A5s. Three-addition representatives are "
                    "metadata-only here and require the pivoted C++ search. "
                    "Compact-only; no finite-volume inference."
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
