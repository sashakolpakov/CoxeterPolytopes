#!/usr/bin/env python3
"""Standalone compact-label gate for the three Record-30 a=1 targets.

The incidence representatives are written literally below.  This file uses
only the Python standard library and rebuilds the complete 108 order-four and
420 order-five labelled Lanner placement tables, every one of the 130 actual
rank-nine vertices of each target, and the finite Coxeter classifier.

The conclusion concerns the displayed incidence targets only until a second
implementation has certified the exhaustive a=1 incidence census.  It makes
no PL, convex-polytopal, ideal-vertex, or finite-volume assertion.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb
from typing import Iterable, Mapping, Sequence


N = 14
DIMENSION = 9
FULL = (1 << N) - 1
EDGES = ((0, 11), (1, 12), (2, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in EDGES)

# One isolate-A4 target, followed by the two one-active-A4 targets.
TARGETS = (
    (
        (3, 4, 5, 6),
        (
            (1, 7, 8, 9, 10),
            (4, 7, 8, 9, 10),
            (7, 8, 9, 10, 11),
            (0, 3, 5, 6, 12),
        ),
    ),
    (
        (0, 3, 4, 5),
        (
            (3, 4, 5, 7, 8),
            (6, 7, 8, 9, 10),
            (1, 6, 9, 10, 11),
            (6, 8, 9, 10, 11),
            (3, 4, 5, 7, 12),
        ),
    ),
    (
        (0, 3, 4, 5),
        (
            (3, 4, 5, 7, 8),
            (6, 7, 8, 9, 10),
            (2, 6, 9, 10, 11),
            (6, 9, 10, 11, 12),
            (1, 3, 4, 5, 13),
        ),
    ),
)

Pair = tuple[int, int]
State = tuple[tuple[Pair, int], ...]

LANNER_REPRESENTATIVES: dict[
    int, tuple[tuple[tuple[Pair, int], ...], ...]
] = {
    4: (
        (((0, 1), 3), ((0, 2), 3), ((0, 3), 5)),
        (((0, 1), 5), ((0, 2), 3), ((1, 3), 3)),
        (((0, 1), 3), ((0, 2), 4), ((1, 3), 5)),
        (((0, 1), 3), ((0, 2), 5), ((1, 3), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 3), 4)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 3), 5)),
        (((0, 1), 3), ((0, 2), 4), ((1, 3), 4), ((2, 3), 3)),
        (((0, 1), 3), ((0, 2), 4), ((1, 3), 5), ((2, 3), 3)),
        (((0, 1), 3), ((0, 2), 5), ((1, 3), 5), ((2, 3), 3)),
    ),
    5: (
        (((0, 1), 3), ((0, 2), 3), ((0, 3), 3), ((1, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 4), ((2, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 5), ((2, 4), 5)),
        (
            ((0, 1), 3),
            ((0, 2), 3),
            ((1, 3), 3),
            ((2, 4), 3),
            ((3, 4), 4),
        ),
    ),
}


def pair(left: int, right: int) -> Pair:
    if left == right:
        raise ValueError("a Coxeter pair needs distinct vertices")
    return (left, right) if left < right else (right, left)


def mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value >> vertex & 1)


def canonical_sha(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def permuted_code(
    order: int, code: tuple[int, ...], permutation: Sequence[int]
) -> tuple[int, ...]:
    pairs = tuple(combinations(range(order), 2))
    source = dict(zip(pairs, code, strict=True))
    moved = {
        pair(permutation[left], permutation[right]): label
        for (left, right), label in source.items()
    }
    return tuple(moved[edge] for edge in pairs)


@lru_cache(maxsize=None)
def lanner_codes(order: int) -> tuple[tuple[int, ...], ...]:
    pairs = tuple(combinations(range(order), 2))
    result = set()
    for representative in LANNER_REPRESENTATIVES[order]:
        labels = {edge: 2 for edge in pairs}
        labels.update(representative)
        code = tuple(labels[edge] for edge in pairs)
        for permutation in permutations(range(order)):
            result.add(permuted_code(order, code, permutation))
    expected = {4: 108, 5: 420}[order]
    if len(result) != expected:
        raise AssertionError(f"order-{order} Lanner placement count changed")
    return tuple(sorted(result))


def path_labels(
    component: set[int], adjacency: Mapping[int, tuple[tuple[int, int], ...]]
) -> tuple[int, ...]:
    current = min(vertex for vertex in component if len(adjacency[vertex]) == 1)
    previous = None
    result = []
    while True:
        following = [row for row in adjacency[current] if row[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("finite path branched")
        other, label = following[0]
        result.append(label)
        previous, current = current, other
    forward = tuple(result)
    return min(forward, tuple(reversed(forward)))


def finite_component(component: set[int], labels: Mapping[Pair, int]) -> bool:
    if len(component) == 1:
        return True
    adjacency = {
        vertex: tuple(
            (other, labels[pair(vertex, other)])
            for other in component
            if other != vertex and labels[pair(vertex, other)] > 2
        )
        for vertex in component
    }
    if sum(map(len, adjacency.values())) // 2 != len(component) - 1:
        return False
    degree_counts = Counter(map(len, adjacency.values()))
    if max(degree_counts) <= 2:
        sequence = path_labels(component, adjacency)
        order = len(component)
        return (
            (order == 2 and sequence[0] in (3, 4, 5))
            or all(label == 3 for label in sequence)
            or sequence == (3,) * (order - 2) + (4,)
            or (order == 3 and sequence == (3, 5))
            or (order == 4 and sequence in ((3, 3, 5), (3, 4, 3)))
        )
    if degree_counts[3] != 1 or any(degree > 3 for degree in degree_counts):
        return False
    if any(
        labels[pair(left, right)] != 3
        for left, right in combinations(component, 2)
        if labels[pair(left, right)] > 2
    ):
        return False
    branch = next(vertex for vertex in component if len(adjacency[vertex]) == 3)
    arms = []
    for neighbor, _label in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(
                other for other, _ in adjacency[current] if other != previous
            )
            previous, current = current, following
            length += 1
        arms.append(length)
    shape = tuple(sorted(arms))
    return shape == (1, 1, len(component) - 3) or shape in {
        (1, 2, 2),
        (1, 2, 3),
        (1, 2, 4),
    }


@lru_cache(maxsize=1_048_576)
def finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    unseen = set(range(order))
    components = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        component = {root}
        queue = [root]
        while queue:
            current = queue.pop()
            neighbors = {
                other
                for other in unseen
                if labels[pair(current, other)] > 2
            }
            unseen.difference_update(neighbors)
            component.update(neighbors)
            queue.extend(neighbors)
        components.append(component)
    return all(finite_component(component, labels) for component in components)


def minimal_transversals(clutter: Sequence[int]) -> tuple[int, ...]:
    result = []
    for candidate in range(1, FULL + 1):
        if not all(candidate & face for face in clutter):
            continue
        if any(
            all((candidate ^ (1 << vertex)) & face for face in clutter)
            for vertex in vertices(candidate)
        ):
            continue
        result.append(candidate)
    return tuple(result)


def actual_vertices(supports: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(map(frozenset, (*EDGES, *supports)))
    return tuple(
        chosen
        for chosen in combinations(range(N), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


def placement_options(support: Sequence[int]) -> tuple[State, ...]:
    ordered = tuple(sorted(support))
    edges = tuple(combinations(ordered, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_codes(len(ordered))
    )


def state_digest(states: Sequence[State]) -> str:
    return canonical_sha(
        [
            [[left, right, label] for (left, right), label in state]
            for state in sorted(states)
        ]
    )


def audit_target(
    target_index: int,
    four: Sequence[int],
    additions: Sequence[Sequence[int]],
) -> dict[str, object]:
    supports = (tuple(sorted(four)),) + tuple(
        tuple(sorted(face)) for face in additions
    )
    addition_masks = tuple(mask(face) for face in additions)
    clutter = (*EDGE_MASKS, mask(four), *addition_masks)
    blocker = minimal_transversals(clutter)
    if len(blocker) != 130 or {face.bit_count() for face in blocker} != {5}:
        raise AssertionError(f"target {target_index} lost its pure 130-blocker")
    rank_nine_vertices = actual_vertices(supports)
    blocker_vertices = {
        tuple(vertex for vertex in range(N) if not complement >> vertex & 1)
        for complement in blocker
    }
    if blocker_vertices != set(rank_nine_vertices):
        raise AssertionError(f"target {target_index} blocker duality failed")

    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[State, ...] = ((),)
    stages = []
    while states and remaining:
        support_index = (
            0
            if not union
            else max(
                remaining,
                key=lambda index: (
                    len(set(supports[index]) & union),
                    len(supports[index]),
                    -len(set(supports[index]) - union),
                    -index,
                ),
            )
        )
        support = supports[support_index]
        old_state_count = len(states)
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {
                    tuple(vertex for vertex in actual if vertex in union)
                    for actual in rank_nine_vertices
                },
                key=lambda row: (len(row), row),
            )
        )
        options = placement_options(support)
        compatible = 0
        passing: set[State] = set()
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
                    finite_code(
                        len(restriction),
                        tuple(
                            labels[pair(left, right)]
                            for left, right in combinations(restriction, 2)
                        ),
                    )
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "option_count": len(options),
                "raw_extension_count": old_state_count * len(options),
                "compatible_extension_count": compatible,
                "restriction_count": len(restrictions),
                "restriction_size_histogram": {
                    str(size): count
                    for size, count in sorted(Counter(map(len, restrictions)).items())
                },
                "passing_state_count": len(states),
                "passing_state_sha256": state_digest(states),
            }
        )
        remaining.remove(support_index)
    if states:
        raise AssertionError(f"target {target_index} survived the label gate")
    return {
        "target_index": target_index,
        "fixed_A4": list(four),
        "additional_A5s": [list(face) for face in additions],
        "actual_rank_nine_vertex_count": len(rank_nine_vertices),
        "actual_vertices_sha256": canonical_sha(rank_nine_vertices),
        "minimal_transversal_histogram": {"5": len(blocker)},
        "stages": stages,
        "terminal_finite_label_state_count": 0,
        "displayed_target_compact_coxeter_impossible": True,
    }


def replay() -> dict[str, object]:
    rows = [
        audit_target(index, four, additions)
        for index, (four, additions) in enumerate(TARGETS)
    ]
    stage_counts = [
        [int(stage["passing_state_count"]) for stage in row["stages"]]
        for row in rows
    ]
    if stage_counts != [[108, 84, 0], [108, 84, 36, 108, 0], [108, 84, 0]]:
        raise AssertionError("Record-30 a=1 label-stage transcript changed")
    result = {
        "kind": "h9_14_record30_a1_displayed_targets_label_replay",
        "compact_only": True,
        "imports_project_module": False,
        "complete_order_four_lanner_placement_count": len(lanner_codes(4)),
        "complete_order_five_lanner_placement_count": len(lanner_codes(5)),
        "displayed_target_count": len(rows),
        "label_audits": rows,
        "passing_state_count_transcripts": stage_counts,
        "terminal_finite_label_state_count_across_targets": 0,
        "displayed_target_compact_no_go_count": len(rows),
        "monotonicity_extension": (
            "Every unspecified meeting pair is assigned the minimal label 2. "
            "Once an actual-vertex principal diagram is nonfinite, increasing "
            "an unspecified Coxeter label only enlarges the nonnegative cosine "
            "adjacency matrix and cannot restore positive definiteness."
        ),
        "theorem_scope": {
            "all_three_displayed_targets_compact_coxeter_impossible": True,
            "whole_record30_a1_schedule_excluded": False,
            "whole_record30_shadow_excluded": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
        "scope": (
            "Standalone complete finite Lanner/actual-rank-nine-vertex gate for "
            "the three exploratory Record-30 a=1 incidence representatives. "
            "Promotion to a whole-schedule theorem requires an independent "
            "exhaustion replay of the incidence census. No finite-volume claim."
        ),
    }
    result["summary_sha256"] = canonical_sha(result)
    return result


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
