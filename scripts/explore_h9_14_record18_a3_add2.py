#!/usr/bin/env python3
"""Exact two-A5 addition gate for all record-18 a=3 shadow-24 orbits.

The input representatives and their full-action orbit keys come from the
production positive-A4 scout.  Every candidate pair is checked against the
exact vertex incidence intervals and every base transversal of order at most
four.  No topology, Coxeter, or finite-volume conclusion is made here.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable, Sequence

from coxeter_polytopes.combinatorics.h9_14_record18 import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
)
from coxeter_polytopes.combinatorics.h9_14_record46 import _base_link_g3_g4
from coxeter_polytopes.combinatorics.h9_14_record92 import independent_masks
from coxeter_polytopes.combinatorics.missing_edges import macaulay_raise


ROOT = Path(__file__).resolve().parents[1]
POSITIVE = ROOT / "data/enumerations/h9_14_record18_positive_A4_scout.json"
VERTEX_COUNT = 14
FULL_MASK = (1 << VERTEX_COUNT) - 1
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES)


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1)


def minimal_transversals(clutter: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        candidate
        for candidate in range(1, FULL_MASK + 1)
        if all(candidate & edge for edge in clutter)
        and all(
            not all((candidate ^ (1 << vertex)) & edge for edge in clutter)
            for vertex in vertices(candidate)
        )
    )


def audit_representative(
    orbit_index: int, orbit_key: Sequence[int], family: Sequence[Sequence[int]]
) -> dict[str, object]:
    fours = tuple(mask(face) for face in family)
    fives = independent_masks(TARGET_GRAPH6, 5)
    shadow = tuple(
        five for five in fives if any(five & four == four for four in fours)
    )
    if len(shadow) != 24:
        raise ValueError("record-18 a3 add2 audit received a non-shadow-24 orbit")
    base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)
    lower = []
    upper = []
    for vertex in range(VERTEX_COUNT):
        a4_incidence = sum(bool(four & (1 << vertex)) for four in fours)
        previous = base_g3[vertex] - a4_incidence
        forced = sum(bool(five & (1 << vertex)) for five in shadow)
        before = base_g4[vertex] + 6 * a4_incidence - forced
        lower.append(max(0, before - macaulay_raise(previous, 3)))
        upper.append(before - 1)
    candidates = tuple(
        five
        for five in fives
        if five not in shadow
        and all(
            not (five & (1 << vertex)) or upper[vertex] >= 1
            for vertex in range(VERTEX_COUNT)
        )
    )
    base = (*EDGE_MASKS, *fours, *shadow)
    small = tuple(
        candidate
        for order in range(1, 5)
        for chosen in combinations(range(VERTEX_COUNT), order)
        for candidate in (mask(chosen),)
        if all(candidate & edge for edge in base)
    )
    cap_count = 0
    no_small_count = 0
    facet_histogram: Counter[int] = Counter()
    target_count = 0
    for additions in combinations(candidates, 2):
        incidences = tuple(
            sum(bool(five & (1 << vertex)) for five in additions)
            for vertex in range(VERTEX_COUNT)
        )
        if any(
            not lower[vertex] <= incidences[vertex] <= upper[vertex]
            for vertex in range(VERTEX_COUNT)
        ):
            continue
        cap_count += 1
        if any(all(transversal & five for five in additions) for transversal in small):
            continue
        no_small_count += 1
        blocker = minimal_transversals((*base, *additions))
        facet_histogram[len(blocker)] += 1
        target_count += len(blocker) == 104
    return {
        "orbit_index": orbit_index,
        "positive_orbit_key": list(map(int, orbit_key)),
        "representative_A4s": [list(map(int, face)) for face in family],
        "forced_A5_shadow_count": len(shadow),
        "additional_A5_count": 2,
        "candidate_count": len(candidates),
        "raw_unordered_pair_count": len(candidates) * (len(candidates) - 1) // 2,
        "base_small_transversal_count": len(small),
        "vertex_cap_feasible_pair_count": cap_count,
        "no_transversal_of_order_at_most_four_pair_count": no_small_count,
        "facet_count_histogram": dict(sorted(facet_histogram.items())),
        "scheduled_104_facet_target_count": target_count,
    }


def audit_all() -> dict[str, object]:
    document = json.loads(POSITIVE.read_text(encoding="utf-8"))
    branch = document["branches"][2]
    rows = []
    for index, orbit in enumerate(branch["orbits"]):
        family = orbit["representative"]
        fours = tuple(mask(face) for face in family)
        fives = independent_masks(TARGET_GRAPH6, 5)
        shadow_count = sum(
            any(five & four == four for four in fours) for five in fives
        )
        if shadow_count == 24:
            rows.append(audit_representative(index, orbit["key"], family))
    if len(rows) != 29:
        raise AssertionError("record-18 a3 shadow-24 orbit count changed")
    if any(row["no_transversal_of_order_at_most_four_pair_count"] for row in rows):
        raise AssertionError("record-18 a3 add2 branch unexpectedly survived")
    return {
        "kind": "h9_14_record18_a3_add2_scout",
        "theorem_certificate": False,
        "compact_only": True,
        "input_positive_orbit_count": 40,
        "shadow_24_add2_orbit_count": len(rows),
        "orbit_audits": rows,
        "aggregate_cap_feasible_pair_count": sum(
            int(row["vertex_cap_feasible_pair_count"]) for row in rows
        ),
        "aggregate_no_small_transversal_pair_count": 0,
        "scheduled_target_count": 0,
        "scope": (
            "Necessary compact incidence addition gate only; the other a3 "
            "orbits, whole shadow, and finite-volume noncompact case remain open."
        ),
    }


def main() -> None:
    print(json.dumps(audit_all(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
