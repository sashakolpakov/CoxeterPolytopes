#!/usr/bin/env python3
"""No-project-import Coxeter-label replay for Record 30 a=2 orbit 14.

The input is a JSON document containing ``pure_target_representatives`` as
five integer masks.  This script uses only the standard library.  It rebuilds
the complete 108 order-four and 420 order-five labelled Lanner placements,
the finite Coxeter classifier, all actual rank-nine vertices, and the sparse
monotonicity join.  It also expands every target through the explicit
``C2 x C2 x S6`` stabilizer to audit quotient disjointness and literal count.

The conclusion is compact-only.  Ideal vertices and finite-volume noncompact
chambers are not tested.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import argparse
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


VERTEX_COUNT = 14
DIMENSION = 9
FULL_MASK = (1 << VERTEX_COUNT) - 1
EDGES = ((0, 11), (1, 12), (2, 13))
FOURS = ((0, 1, 3, 4), (0, 2, 3, 4))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in EDGES)
FOUR_MASKS = tuple(sum(1 << vertex for vertex in face) for face in FOURS)

Pair = tuple[int, int]
LabelState = tuple[tuple[Pair, int], ...]


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


def _pair(left: int, right: int) -> Pair:
    if left == right:
        raise ValueError("a Coxeter pair needs two distinct vertices")
    return (left, right) if left < right else (right, left)


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _vertices(value: int) -> tuple[int, ...]:
    return tuple(
        vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1
    )


def _canonical_sha(value: object) -> str:
    return sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode()
    ).hexdigest()


def _permuted_code(
    order: int, code: tuple[int, ...], permutation: Sequence[int]
) -> tuple[int, ...]:
    pairs = tuple(combinations(range(order), 2))
    source = dict(zip(pairs, code, strict=True))
    target = {
        _pair(permutation[left], permutation[right]): label
        for (left, right), label in source.items()
    }
    return tuple(target[edge] for edge in pairs)


@lru_cache(maxsize=None)
def _lanner_codes(order: int) -> tuple[tuple[int, ...], ...]:
    if order not in LANNER_REPRESENTATIVES:
        raise ValueError("only order-four and order-five Lanner tables are used")
    pairs = tuple(combinations(range(order), 2))
    placements = set()
    for representative in LANNER_REPRESENTATIVES[order]:
        labels = {edge: 2 for edge in pairs}
        labels.update(representative)
        code = tuple(labels[edge] for edge in pairs)
        for permutation in permutations(range(order)):
            placements.add(_permuted_code(order, code, permutation))
    expected = {4: 108, 5: 420}[order]
    if len(placements) != expected:
        raise AssertionError(f"order-{order} Lanner placement count changed")
    return tuple(sorted(placements))


def _path_labels(
    component: set[int],
    adjacency: Mapping[int, tuple[tuple[int, int], ...]],
) -> tuple[int, ...]:
    endpoint = min(
        vertex for vertex in component if len(adjacency[vertex]) == 1
    )
    labels = []
    previous = None
    current = endpoint
    while True:
        following = [row for row in adjacency[current] if row[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("finite path classifier branched")
        next_vertex, label = following[0]
        labels.append(label)
        previous, current = current, next_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def _arm_lengths(
    component: set[int],
    adjacency: Mapping[int, tuple[tuple[int, int], ...]],
) -> tuple[int, ...]:
    branch = next(
        vertex for vertex in component if len(adjacency[vertex]) == 3
    )
    lengths = []
    for neighbor, _label in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(
                vertex
                for vertex, _ in adjacency[current]
                if vertex != previous
            )
            previous, current = current, following
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def _finite_component(
    component: set[int], labels: Mapping[Pair, int]
) -> bool:
    if len(component) == 1:
        return True
    adjacency = {
        vertex: tuple(
            (other, labels[_pair(vertex, other)])
            for other in component
            if other != vertex and labels[_pair(vertex, other)] > 2
        )
        for vertex in component
    }
    if sum(map(len, adjacency.values())) // 2 != len(component) - 1:
        return False
    degrees = Counter(map(len, adjacency.values()))
    if max(degrees) <= 2:
        sequence = _path_labels(component, adjacency)
        order = len(component)
        if order == 2:
            return sequence[0] in (3, 4, 5)
        return (
            all(label == 3 for label in sequence)
            or sequence == (3,) * (order - 2) + (4,)
            or (order == 3 and sequence == (3, 5))
            or (order == 4 and sequence in ((3, 3, 5), (3, 4, 3)))
        )
    if degrees[3] != 1 or any(degree > 3 for degree in degrees):
        return False
    if any(
        labels[_pair(left, right)] != 3
        for left in component
        for right in component
        if left < right and labels[_pair(left, right)] > 2
    ):
        return False
    arms = _arm_lengths(component, adjacency)
    return arms == (1, 1, len(component) - 3) or arms in {
        (1, 2, 2),
        (1, 2, 3),
        (1, 2, 4),
    }


@lru_cache(maxsize=1_048_576)
def _finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    unseen = set(range(order))
    components = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        frontier = [start]
        while frontier:
            vertex = frontier.pop()
            neighbors = {
                other
                for other in unseen
                if labels[_pair(vertex, other)] > 2
            }
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        components.append(component)
    return all(_finite_component(component, labels) for component in components)


def _minimal_transversals(clutter: Sequence[int]) -> tuple[int, ...]:
    answer = []
    for candidate in range(1, FULL_MASK + 1):
        if not all(candidate & face for face in clutter):
            continue
        if all(
            not all((candidate ^ (1 << vertex)) & face for face in clutter)
            for vertex in _vertices(candidate)
        ):
            answer.append(candidate)
    return tuple(answer)


def _actual_vertices(supports: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(map(frozenset, (*EDGES, *supports)))
    return tuple(
        chosen
        for chosen in combinations(range(VERTEX_COUNT), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


def _placement_options(support: Sequence[int]) -> tuple[LabelState, ...]:
    ordered = tuple(sorted(support))
    edges = tuple(combinations(ordered, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in _lanner_codes(len(ordered))
    )


def _state_digest(states: Sequence[LabelState]) -> str:
    return _canonical_sha(
        [
            [[left, right, label] for (left, right), label in state]
            for state in sorted(states)
        ]
    )


def _label_family(
    index: int,
    additions: Sequence[int],
    *,
    fours: Sequence[Sequence[int]] = FOURS,
    require_no_go: bool = True,
) -> dict[str, object]:
    normalized_fours = tuple(tuple(sorted(face)) for face in fours)
    supports = normalized_fours + tuple(_vertices(face) for face in additions)
    actual_vertices = _actual_vertices(supports)
    clutter = (
        EDGE_MASKS
        + tuple(_mask(face) for face in normalized_fours)
        + tuple(additions)
    )
    blocker = _minimal_transversals(clutter)
    if len(blocker) != 128 or {face.bit_count() for face in blocker} != {5}:
        raise AssertionError(f"target {index} lost pure 128-facet blocker")
    complements = {
        tuple(vertex for vertex in range(VERTEX_COUNT) if not face >> vertex & 1)
        for face in blocker
    }
    if complements != set(actual_vertices):
        raise AssertionError(f"target {index} blocker/vertex duality failed")

    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[LabelState, ...] = ((),)
    stages = []
    while states and remaining:
        support_index = (
            0
            if not union
            else max(
                remaining,
                key=lambda item: (
                    len(set(supports[item]) & union),
                    len(supports[item]),
                    -len(set(supports[item]) - union),
                    -item,
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
                    for actual in actual_vertices
                },
                key=lambda row: (len(row), row),
            )
        )
        options = _placement_options(support)
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
                    _pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(
                    _finite_code(
                        len(restriction),
                        tuple(
                            labels[_pair(left, right)]
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
                "restriction_size_histogram": dict(
                    sorted(Counter(map(len, restrictions)).items())
                ),
                "passing_state_count": len(states),
                "passing_state_sha256": _state_digest(states),
            }
        )
        remaining.remove(support_index)
    if states and require_no_go:
        raise AssertionError(f"target {index} survived the sparse label gate")
    return {
        "target_index": index,
        "fixed_A4_masks": [_mask(face) for face in normalized_fours],
        "additional_A5_masks": list(additions),
        "actual_rank_nine_vertex_count": len(actual_vertices),
        "actual_vertices_sha256": _canonical_sha(actual_vertices),
        "minimal_transversal_histogram": {"5": 128},
        "stages": stages,
        "terminal_sparse_state_count": len(states),
        "terminal_sparse_states_sha256": _state_digest(states),
        "unrestricted_compact_label_no_go": not states,
    }


def _stabilizer() -> tuple[tuple[int, ...], ...]:
    answer = []
    for swap_12 in (False, True):
        for swap_34 in (False, True):
            for tail in permutations(range(5, 11)):
                permutation = list(range(VERTEX_COUNT))
                if swap_12:
                    permutation[1], permutation[2] = 2, 1
                    permutation[12], permutation[13] = 13, 12
                if swap_34:
                    permutation[3], permutation[4] = 4, 3
                for source, target in zip(range(5, 11), tail, strict=True):
                    permutation[source] = target
                answer.append(tuple(permutation))
    if len(answer) != 2880 or len(set(answer)) != 2880:
        raise AssertionError("explicit stabilizer changed")
    fixed_edges = {frozenset(edge) for edge in EDGES}
    fixed_fours = {frozenset(face) for face in FOURS}
    for permutation in answer:
        if {
            frozenset(permutation[vertex] for vertex in edge) for edge in EDGES
        } != fixed_edges:
            raise AssertionError("claimed stabilizer does not preserve edges")
        if {
            frozenset(permutation[vertex] for vertex in face) for face in FOURS
        } != fixed_fours:
            raise AssertionError("claimed stabilizer does not preserve A4s")
    return tuple(answer)


def _permuted_mask(value: int, permutation: Sequence[int]) -> int:
    return _mask(permutation[vertex] for vertex in _vertices(value))


def _family_orbit(
    family: Sequence[int], stabilizer: Sequence[Sequence[int]]
) -> frozenset[tuple[int, ...]]:
    return frozenset(
        tuple(sorted(_permuted_mask(face, permutation) for face in family))
        for permutation in stabilizer
    )


def _target_families(document: Mapping[str, object]) -> tuple[tuple[int, ...], ...]:
    raw = document.get("pure_target_representatives")
    if not isinstance(raw, list) or len(raw) != 37:
        raise ValueError("expected exactly 37 target representatives")
    families = []
    for index, family in enumerate(raw):
        if not isinstance(family, list) or len(family) != 5:
            raise ValueError(f"target family {index} is malformed")
        normalized = tuple(sorted(int(face) for face in family))
        if len(set(normalized)) != 5 or any(face.bit_count() != 5 for face in normalized):
            raise ValueError(f"target family {index} is not five distinct A5s")
        families.append(normalized)
    return tuple(families)


def replay(document: Mapping[str, object]) -> dict[str, object]:
    families = _target_families(document)
    stabilizer = _stabilizer()
    literal_union: set[tuple[int, ...]] = set()
    orbit_sizes = []
    for index, family in enumerate(families):
        orbit = _family_orbit(family, stabilizer)
        if literal_union.intersection(orbit):
            raise AssertionError(f"target orbit {index} overlaps an earlier orbit")
        literal_union.update(orbit)
        orbit_sizes.append(len(orbit))
    histogram = dict(sorted(Counter(orbit_sizes).items()))
    expected_histogram = {60: 17, 90: 1, 120: 14, 180: 1, 240: 1, 360: 3}
    if histogram != expected_histogram or len(literal_union) != 4290:
        raise AssertionError("literal target orbit expansion changed")

    label_rows = [_label_family(index, family) for index, family in enumerate(families)]
    transcript = [row["stages"] for row in label_rows]
    result = {
        "kind": "h9_14_record30_a2_orbit14_no_project_label_replay",
        "compact_only": True,
        "imports_project_module": False,
        "complete_order_four_lanner_placement_count": len(_lanner_codes(4)),
        "complete_order_five_lanner_placement_count": len(_lanner_codes(5)),
        "stabilizer_order": len(stabilizer),
        "target_orbit_count": len(families),
        "target_family_sha256": _canonical_sha(families),
        "target_orbit_sizes": orbit_sizes,
        "orbit_size_histogram": {str(key): value for key, value in histogram.items()},
        "literal_target_count": len(literal_union),
        "literal_target_union_sha256": _canonical_sha(sorted(literal_union)),
        "label_audits": label_rows,
        "label_stage_transcript_sha256": _canonical_sha(transcript),
        "terminal_sparse_state_count_across_orbits": sum(
            int(row["terminal_sparse_state_count"]) for row in label_rows
        ),
        "unrestricted_compact_label_no_go_orbit_count": sum(
            bool(row["unrestricted_compact_label_no_go"]) for row in label_rows
        ),
        "theorem": {
            "whole_record30_a2_orbit14_nerve_compact_coxeter_impossible": True,
            "whole_record30_a2_schedule_excluded": False,
            "whole_record30_shadow_excluded": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
        "scope": (
            "Complete sparse Lanner and actual-rank-nine-vertex replay for all "
            "37 exact incidence target orbits in Record-30 compact a=2 "
            "positive-A4 orbit 14. The zero-state conclusion is unrestricted "
            "within finite labels 2,3,4,5 by Perron--Frobenius monotonicity "
            "and principal-submatrix inheritance. No finite-volume inference."
        ),
    }
    result["summary_sha256"] = _canonical_sha(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    args = parser.parse_args()
    document = json.loads(args.input.read_text(encoding="utf-8"))
    print(json.dumps(replay(document), sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
