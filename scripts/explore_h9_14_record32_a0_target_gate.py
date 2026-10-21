#!/usr/bin/env python3
"""Stream record-32 ``a=0`` targets into exact incidence and label gates.

Input is the compact mask-only JSON emitted by
``explore_h9_14_record32_a0.cpp --target-output``.  Every target is first
dualized directly.  The script records the complete minimal-transversal size
histogram and the ridge-degree histogram, then joins complete order-five
Lanner placement tables while imposing every actual rank-nine vertex
restriction after each support.

The sparse label gate is unrestricted only when its terminal state count is
zero: changing any unassigned right angle to a non-right finite label can
only increase the nonnegative cosine adjacency, so it cannot repair a failed
elliptic principal restriction.  A positive or capped terminal count is not a
survivor theorem; it is handed to the cross-edge/inertia stage.

This is an exact scout, not a frozen whole-shadow certificate.  It does not
infer ordinary polytopality, PL topology, compact realizability, or finite
volume.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_record32 import TARGET_MISSING_EDGES


VERTEX_COUNT = 14
DIMENSION = 9
FULL_MASK = (1 << VERTEX_COUNT) - 1
EDGE_MASKS = tuple(
    sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES
)

LabelState = tuple[tuple[Pair, int], ...]


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1)


def minimal_transversals(hyperedges: Sequence[int]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, FULL_MASK + 1):
        if not all(candidate & edge for edge in hyperedges):
            continue
        if any(
            all((candidate ^ (1 << vertex)) & edge for edge in hyperedges)
            for vertex in vertices(candidate)
        ):
            continue
        answers.append(candidate)
    return tuple(answers)


def ridge_degree_histogram(facet_complements: Sequence[int]) -> tuple[tuple[int, int], ...]:
    ridge_degrees: Counter[int] = Counter()
    for complement in facet_complements:
        facet = FULL_MASK ^ complement
        for removed in vertices(facet):
            ridge_degrees[facet ^ (1 << removed)] += 1
    return tuple(sorted(Counter(ridge_degrees.values()).items()))


def rank_nine_vertices(missing_fives: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(frozenset(edge) for edge in TARGET_MISSING_EDGES) + tuple(
        frozenset(vertices(face)) for face in missing_fives
    )
    return tuple(
        chosen
        for chosen in combinations(range(VERTEX_COUNT), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


@lru_cache(maxsize=None)
def placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(support))
    )


@lru_cache(maxsize=524_288)
def finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return is_finite_coxeter_labeling(range(order), labels)


def restriction_finite(
    restriction: Sequence[int], labels: Mapping[Pair, int]
) -> bool:
    return finite_code(
        len(restriction),
        tuple(
            labels[pair(left, right)]
            for left, right in combinations(restriction, 2)
        ),
    )


def next_support(
    supports: Sequence[tuple[int, ...]], remaining: set[int], union: set[int]
) -> int:
    if not union:
        return min(remaining)
    return max(
        remaining,
        key=lambda index: (
            len(set(supports[index]) & union),
            -len(set(supports[index]) - union),
            -index,
        ),
    )


def sparse_label_gate(
    missing_fives: Sequence[int], state_cap: int
) -> dict[str, object]:
    supports = tuple(vertices(face) for face in missing_fives)
    actual_vertices = rank_nine_vertices(missing_fives)
    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[LabelState, ...] = ((),)
    stages = []
    capped = False
    while states and remaining:
        support_index = next_support(supports, remaining, union)
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {
                    tuple(facet for facet in vertex if facet in union)
                    for vertex in actual_vertices
                },
                key=lambda row: (len(row), row),
            )
        )
        options = placement_options(support)
        passing: set[LabelState] = set()
        compatible = 0
        raw = len(states) * len(options)
        for state in states:
            old = dict(state)
            for option in options:
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(
                    restriction_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
                    if len(passing) > state_cap:
                        capped = True
                        break
            if capped:
                break
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "raw_extension_count": raw,
                "compatible_extension_count_before_cap": compatible,
                "passing_state_count": len(states),
                "restriction_size_histogram": dict(
                    sorted(Counter(map(len, restrictions)).items())
                ),
                "state_cap_exceeded": capped,
            }
        )
        remaining.remove(support_index)
        if capped:
            break
    return {
        "rank_nine_vertex_count": len(actual_vertices),
        "stages": stages,
        "terminal_sparse_state_count": len(states),
        "state_cap_exceeded": capped,
        "unrestricted_compact_label_no_go": not capped and not states,
    }


def audit_family(
    index: int,
    missing_fives: Sequence[int],
    state_cap: int,
    incidence_only: bool,
) -> dict[str, object]:
    family = tuple(sorted(map(int, missing_fives)))
    transversals = minimal_transversals((*EDGE_MASKS, *family))
    histogram = tuple(sorted(Counter(value.bit_count() for value in transversals).items()))
    five_transversals = tuple(
        value for value in transversals if value.bit_count() == 5
    )
    ridge_histogram = ridge_degree_histogram(five_transversals)
    result: dict[str, object] = {
        "family_index": index,
        "missing_five_masks": list(family),
        "minimal_transversal_size_histogram": dict(histogram),
        "pure_order_five_blocker": histogram == ((5, 114),),
        "facet_complement_sha256": sha256(
            b"".join(value.to_bytes(2, "little") for value in five_transversals)
        ).hexdigest(),
        "ridge_degree_histogram": dict(ridge_histogram),
        # A ridge histogram computed only from order-five transversals is not a
        # pseudomanifold certificate when larger minimal transversals also
        # exist.  Require the blocker itself to be pure before exposing the
        # ridge flag or entering the label gate.
        "all_ridge_degrees_two": (
            histogram == ((5, 114),) and ridge_histogram == ((2, 513),)
        ),
    }
    if not incidence_only and result["all_ridge_degrees_two"]:
        result["sparse_label_gate"] = sparse_label_gate(family, state_cap)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("targets", type=Path)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--state-cap", type=int, default=250_000)
    parser.add_argument(
        "--incidence-only",
        action="store_true",
        help="run only blocker purity and the exact ridge-degree gate",
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="omit per-family audit rows but retain exact passing indices",
    )
    args = parser.parse_args()
    document = json.loads(args.targets.read_text(encoding="utf-8"))
    families = tuple(tuple(map(int, row)) for row in document["target_family_masks"])
    if len(families) != int(document["exact_unordered_target_family_count"]):
        raise AssertionError("record-32 a=0 target file count changed")
    stop = len(families) if args.limit is None else min(
        len(families), args.start + args.limit
    )
    audits = [
        audit_family(index, families[index], args.state_cap, args.incidence_only)
        for index in range(args.start, stop)
    ]
    summary: dict[str, object] = {
                "kind": "h9_14_record32_a0_target_incidence_and_label_scout",
                "theorem_certificate": False,
                "target_file_family_count": len(families),
                "audited_range": [args.start, stop],
                "pure_blocker_count": sum(
                    audit["pure_order_five_blocker"] for audit in audits
                ),
                "ridge_pseudomanifold_count": sum(
                    audit["all_ridge_degrees_two"] for audit in audits
                ),
                "ridge_pseudomanifold_indices": [
                    audit["family_index"]
                    for audit in audits
                    if audit["all_ridge_degrees_two"]
                ],
                "unrestricted_sparse_label_no_go_count": sum(
                    audit.get("sparse_label_gate", {}).get(
                        "unrestricted_compact_label_no_go", False
                    )
                    for audit in audits
                ),
                "capped_label_count": sum(
                    audit.get("sparse_label_gate", {}).get(
                        "state_cap_exceeded", False
                    )
                    for audit in audits
                ),
                "ridge_candidate_label_transcripts": [
                    {
                        "family_index": audit["family_index"],
                        "stage_support_and_passing_counts": [
                            [
                                stage["support_index"],
                                stage["compatible_extension_count_before_cap"],
                                stage["passing_state_count"],
                            ]
                            for stage in audit.get("sparse_label_gate", {}).get(
                                "stages", ()
                            )
                        ],
                        "state_cap_exceeded": audit.get(
                            "sparse_label_gate", {}
                        ).get("state_cap_exceeded", False),
                        "unrestricted_compact_label_no_go": audit.get(
                            "sparse_label_gate", {}
                        ).get("unrestricted_compact_label_no_go", False),
                    }
                    for audit in audits
                    if audit["all_ridge_degrees_two"]
                    and "sparse_label_gate" in audit
                ],
                "scope": (
                    "Exact blocker/ridge and sparse Lanner/actual-vertex scout. "
                    "Only zero uncapped sparse terminals are unrestricted compact "
                    "label no-gos. Nonzero/capped terminals require cross-edge and "
                    "principal-inertia completion. No finite-volume claim."
                ),
            }
    if not args.summary_only:
        summary["family_audits"] = audits
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
