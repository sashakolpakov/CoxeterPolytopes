#!/usr/bin/env python3
"""Exact exploratory compact-label gate for the seven open selected nerves.

This script deliberately makes no whole-shadow claim.  It applies the complete
order-four/order-five Lanner placement tables to each already-certified
ordinary-polytopal selection and checks every actual rank-nine vertex after
each deterministic join.  Undeclared finite edges are set to label two; by
Perron--Frobenius monotonicity, failure at that sparse completion is an
unrestricted compact-label obstruction.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from importlib import import_module
from itertools import combinations
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)


TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
RECORDS = (0, 1, 3, 4, 5, 7, 9)
LabelState = tuple[tuple[Pair, int], ...]


@lru_cache(maxsize=None)
def placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(support))
    )


@lru_cache(maxsize=262_144)
def finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return is_finite_coxeter_labeling(range(order), labels)


def restriction_is_finite(
    restriction: tuple[int, ...], labels: Mapping[Pair, int]
) -> bool:
    code = tuple(
        labels[pair(left, right)]
        for left, right in combinations(restriction, 2)
    )
    return finite_code(len(restriction), code)


def next_support_index(
    supports: Sequence[tuple[int, ...]],
    remaining: Iterable[int],
    cluster_union: Iterable[int],
) -> int:
    choices = tuple(sorted(remaining))
    union = set(cluster_union)
    if not union:
        return choices[0]
    return max(
        choices,
        key=lambda index: (
            len(set(supports[index]) & union),
            -len(set(supports[index]) - union),
            -index,
        ),
    )


def audit_record(record: int) -> bool:
    module = import_module(
        f"coxeter_polytopes.combinatorics.h9_14_record{record}"
    )
    missing_edges = tuple(tuple(face) for face in module.TARGET_MISSING_EDGES)
    missing_fours = tuple(
        tuple(face) for face in getattr(module, "TARGET_MISSING_FOURS", ())
    )
    missing_fives = tuple(tuple(face) for face in module.TARGET_MISSING_FIVES)
    supports = missing_fours + missing_fives
    complements = tuple(sorted(module.target_facet_complements()))
    vertices = tuple(
        sorted(
            tuple(v for v in range(TARGET_VERTEX_COUNT) if (FULL_MASK ^ value) >> v & 1)
            for value in complements
        )
    )
    if any(len(vertex) != TARGET_DIMENSION for vertex in vertices):
        raise AssertionError(f"record {record}: non-rank-nine selected facet")

    remaining = set(range(len(supports)))
    states: tuple[LabelState, ...] = ((),)
    cluster_union: tuple[int, ...] = ()
    print(
        "RECORD",
        record,
        "VERTICES",
        len(vertices),
        "SUPPORTS",
        Counter(map(len, supports)),
        flush=True,
    )
    stage = 0
    while states and remaining:
        support_index = next_support_index(supports, remaining, cluster_union)
        support = supports[support_index]
        cluster_union = tuple(sorted(set(cluster_union) | set(support)))
        cluster_set = set(cluster_union)
        restrictions = tuple(
            sorted(
                {
                    tuple(v for v in vertex if v in cluster_set)
                    for vertex in vertices
                },
                key=lambda value: (len(value), value),
            )
        )
        if any(
            set(dotted) <= set(restriction)
            for dotted in missing_edges
            for restriction in restrictions
        ):
            raise AssertionError(f"record {record}: actual vertex contains dotted pair")

        options = placement_options(tuple(support))
        compatible = 0
        passing: set[LabelState] = set()
        for state in states:
            old = dict(state)
            for option in options:
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    pair(left, right): 2
                    for left, right in combinations(cluster_union, 2)
                }
                labels.update(merged)
                if all(
                    restriction_is_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        stage += 1
        print(
            "STAGE",
            stage,
            "INDEX",
            support_index,
            "SUPPORT",
            support,
            "UNION",
            len(cluster_union),
            "RESTRICTIONS",
            Counter(map(len, restrictions)),
            "OPTIONS",
            len(options),
            "RAW",
            len(states) * len(options),
            "COMPATIBLE",
            compatible,
            "PASS",
            len(passing),
            flush=True,
        )
        states = tuple(sorted(passing))
        remaining.remove(support_index)
    print(
        "RESULT",
        record,
        "NO_GO" if not states else "SURVIVORS",
        len(states),
        "UNJOINED",
        len(remaining),
        flush=True,
    )
    return not states


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("records", nargs="*", type=int, default=list(RECORDS))
    args = parser.parse_args()
    if any(record not in RECORDS for record in args.records):
        parser.error(f"records must lie in {RECORDS}")
    for record in args.records:
        audit_record(record)


if __name__ == "__main__":
    main()
