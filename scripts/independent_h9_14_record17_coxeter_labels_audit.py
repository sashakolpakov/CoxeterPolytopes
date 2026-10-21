#!/usr/bin/env python3
"""Independent replay of the selected record-17 compact-Coxeter no-go.

This verifier does not import the production record-17 gate, the reusable
common-four-core gate, or the project's Coxeter classifier/Lanner tables.  It
rebuilds the 420 placements from five classical representatives, joins them
on the common core, independently classifies every finite vertex restriction,
and compares the resulting census with the frozen archive.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record17_coxeter_labels.json"
POLYTOPALITY = ROOT / "data/enumerations/h9_14_record17_polytopality.json"
FACET_COUNT = 14
DIMENSION = 9
OUTER = (0, 1, 2)
CORE = (3, 4, 5, 6)
SUPPORTS = tuple(tuple(sorted((*CORE, outer))) for outer in OUTER)
MISSING_EDGES = ((0, 12), (1, 12), (0, 13), (2, 13), (12, 13))
MISSING_FIVES = (
    (0, 3, 4, 5, 6),
    (1, 3, 4, 5, 6),
    (2, 3, 4, 5, 6),
    (0, 5, 6, 7, 8),
    (0, 7, 8, 9, 10),
    (1, 2, 3, 4, 11),
    (1, 2, 9, 10, 11),
    (7, 8, 9, 10, 11),
)

# Five order-five Lanner isomorphism representatives.  Unlisted pairs have
# label two.  Permutation closure below is independent of the production table.
LANNER_REPRESENTATIVES = (
    (((0, 1), 3), ((0, 2), 3), ((0, 3), 3), ((1, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 4), ((2, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 5), ((2, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 3), ((3, 4), 4)),
)


def _pair(left: int, right: int) -> tuple[int, int]:
    if left == right:
        raise ValueError("a pair needs distinct entries")
    return (left, right) if left < right else (right, left)


def _path_labels(
    component: frozenset[int],
    adjacency: Mapping[int, tuple[tuple[int, int], ...]],
) -> tuple[int, ...]:
    endpoint = min(vertex for vertex in component if len(adjacency[vertex]) == 1)
    previous: int | None = None
    current = endpoint
    labels: list[int] = []
    while True:
        following = [item for item in adjacency[current] if item[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("independent path traversal branched")
        next_vertex, label = following[0]
        labels.append(label)
        previous, current = current, next_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def _arm_lengths(
    component: frozenset[int],
    adjacency: Mapping[int, tuple[tuple[int, int], ...]],
) -> tuple[int, ...]:
    branch = next(vertex for vertex in component if len(adjacency[vertex]) == 3)
    lengths = []
    for neighbor, _label in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(
                vertex for vertex, _ in adjacency[current] if vertex != previous
            )
            previous, current = current, following
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def _finite_component(
    component: frozenset[int], labels: Mapping[tuple[int, int], int]
) -> bool:
    order = len(component)
    if order == 1:
        return True
    adjacency = {
        vertex: tuple(
            (other, labels[_pair(vertex, other)])
            for other in component
            if other != vertex and labels[_pair(vertex, other)] > 2
        )
        for vertex in component
    }
    if sum(map(len, adjacency.values())) // 2 != order - 1:
        return False
    degrees = Counter(len(adjacency[vertex]) for vertex in component)
    if max(degrees) <= 2:
        sequence = _path_labels(component, adjacency)
        if order == 2:
            return sequence[0] in (3, 4, 5)
        if all(label == 3 for label in sequence):
            return True
        if sequence == (3,) * (order - 2) + (4,):
            return True
        if order == 3 and sequence == (3, 5):
            return True
        if order == 4 and sequence in ((3, 3, 5), (3, 4, 3)):
            return True
        return False
    if degrees[3] != 1 or any(degree > 3 for degree in degrees):
        return False
    if any(
        labels[_pair(left, right)] != 3
        for left, right in combinations(component, 2)
        if labels[_pair(left, right)] > 2
    ):
        return False
    arms = _arm_lengths(component, adjacency)
    return arms == (1, 1, order - 3) or arms in {
        (1, 2, 2),
        (1, 2, 3),
        (1, 2, 4),
    }


def _finite(
    facets: Iterable[int], labels: Mapping[tuple[int, int], int]
) -> bool:
    vertices = tuple(sorted(set(facets)))
    unseen = set(vertices)
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
        components.append(frozenset(component))
    return all(_finite_component(component, labels) for component in components)


def _representative_code(
    representative: Sequence[tuple[tuple[int, int], int]],
) -> tuple[int, ...]:
    pairs = tuple(combinations(range(5), 2))
    labels = {edge: 2 for edge in pairs}
    labels.update(representative)
    return tuple(labels[edge] for edge in pairs)


def _permuted_code(
    code: Sequence[int], permutation: Sequence[int]
) -> tuple[int, ...]:
    pairs = tuple(combinations(range(5), 2))
    source = dict(zip(pairs, code, strict=True))
    moved = {
        _pair(permutation[left], permutation[right]): label
        for (left, right), label in source.items()
    }
    return tuple(moved[edge] for edge in pairs)


def _lanner_codes() -> tuple[tuple[int, ...], ...]:
    codes = tuple(
        sorted(
            {
                _permuted_code(_representative_code(representative), image)
                for representative in LANNER_REPRESENTATIVES
                for image in permutations(range(5))
            }
        )
    )
    if len(codes) != 420:
        raise AssertionError("independent Lanner placement count changed")
    canonical = {
        min(_permuted_code(code, image) for image in permutations(range(5)))
        for code in codes
    }
    if len(canonical) != 5:
        raise AssertionError("independent Lanner quotient count changed")
    local_pairs = tuple(combinations(range(5), 2))
    for code in codes:
        labels = dict(zip(local_pairs, code, strict=True))
        if _finite(range(5), labels):
            raise AssertionError("a declared Lanner placement is finite")
        for deleted in range(5):
            if not _finite((v for v in range(5) if v != deleted), labels):
                raise AssertionError("a declared Lanner placement is not minimal")
    return codes


def _support_states(
    support: Sequence[int], codes: Sequence[Sequence[int]]
) -> tuple[tuple[tuple[tuple[int, int], int], ...], ...]:
    pairs = tuple(combinations(tuple(sorted(support)), 2))
    return tuple(tuple(zip(pairs, code, strict=True)) for code in codes)


def _merge(
    left: Sequence[tuple[tuple[int, int], int]],
    right: Sequence[tuple[tuple[int, int], int]],
) -> tuple[tuple[tuple[int, int], int], ...] | None:
    labels = dict(left)
    for edge, label in right:
        previous = labels.get(edge)
        if previous is not None and previous != label:
            return None
        labels[edge] = label
    return tuple(sorted(labels.items()))


def _join_tables(
    codes: Sequence[Sequence[int]],
) -> tuple[tuple[tuple[tuple[tuple[int, int], int], ...], ...], tuple[int, ...]]:
    states: tuple[tuple[tuple[tuple[int, int], int], ...], ...] = ((),)
    counts = []
    for support in SUPPORTS:
        options = _support_states(support, codes)
        states = tuple(
            sorted(
                {
                    merged
                    for state in states
                    for option in options
                    if (merged := _merge(state, option)) is not None
                }
            )
        )
        counts.append(len(states))
    return states, tuple(counts)


def _vertices() -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(frozenset(face) for face in MISSING_EDGES + MISSING_FIVES)
    return tuple(
        facets
        for facets in combinations(range(FACET_COUNT), DIMENSION)
        if not any(face <= set(facets) for face in forbidden)
    )


def replay() -> dict[str, object]:
    vertices = _vertices()
    if len(vertices) != 110:
        raise AssertionError("independent record-17 vertex count changed")
    cluster = tuple(sorted((*OUTER, *CORE)))
    cluster_set = set(cluster)
    restrictions = tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in cluster_set)
                for vertex in vertices
            },
            key=lambda value: (len(value), value),
        )
    )
    codes = _lanner_codes()
    states, join_counts = _join_tables(codes)
    labels = []
    for state in states:
        full = {_pair(left, right): 2 for left, right in combinations(cluster, 2)}
        full.update(state)
        labels.append(full)

    stage_counts = []
    survivors = tuple(labels)
    rank_six_failure_histogram: tuple[tuple[int, int], ...] = ()
    for order in sorted(set(map(len, restrictions))):
        tests = tuple(value for value in restrictions if len(value) == order)
        if order == 6:
            rank_six_failure_histogram = tuple(
                sorted(
                    Counter(
                        sum(not _finite(test, state) for test in tests)
                        for state in survivors
                    ).items()
                )
            )
        survivors = tuple(
            state
            for state in survivors
            if all(_finite(test, state) for test in tests)
        )
        stage_counts.append((order, len(tests), len(survivors)))

    polytope = json.loads(POLYTOPALITY.read_text(encoding="utf-8"))
    target = polytope["target"]
    archived_edges = tuple(tuple(face) for face in target["missing_edges"])
    archived_fives = tuple(tuple(face) for face in target["missing_fives"])
    if archived_edges != MISSING_EDGES or archived_fives != MISSING_FIVES:
        raise AssertionError("independent polytopality target mismatch")
    dotted_in_cluster = tuple(
        edge
        for edge in MISSING_EDGES
        if edge[0] in cluster_set and edge[1] in cluster_set
    )
    result: dict[str, object] = {
        "record_index": 17,
        "vertex_count": len(vertices),
        "lanner_placement_count": len(codes),
        "lanner_isomorphism_type_count": 5,
        "supports": [list(support) for support in SUPPORTS],
        "join_counts": list(join_counts),
        "restriction_size_histogram": {
            str(order): count
            for order, count in sorted(Counter(map(len, restrictions)).items())
        },
        "stage_counts": [
            {
                "order": order,
                "restriction_count": tests,
                "survivors": count,
            }
            for order, tests, count in stage_counts
        ],
        "rank_six_failure_multiplicity_histogram": {
            str(multiplicity): count
            for multiplicity, count in rank_six_failure_histogram
        },
        "dotted_pairs_in_cluster_union": [list(edge) for edge in dotted_in_cluster],
        "unrestricted_compact_coxeter_no_go": (
            not dotted_in_cluster and stage_counts[-1][-1] == 0
        ),
    }
    expected = {
        "vertex_count": 110,
        "lanner_placement_count": 420,
        "join_counts": [420, 1_332, 4_620],
        "restriction_size_histogram": {"4": 9, "5": 16, "6": 4},
        "stage_counts": [
            {"order": 4, "restriction_count": 9, "survivors": 4_620},
            {"order": 5, "restriction_count": 16, "survivors": 80},
            {"order": 6, "restriction_count": 4, "survivors": 0},
        ],
        "dotted_pairs_in_cluster_union": [],
        "unrestricted_compact_coxeter_no_go": True,
    }
    for key, value in expected.items():
        if result[key] != value:
            raise AssertionError(f"independent record-17 field {key} changed")
    return result


def compare_archive(result: Mapping[str, object]) -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    theorem = document["theorem"]
    gate = theorem["decisive_common_four_core_gate"]
    archived_stages = [
        {
            "order": stage["order"],
            "restriction_count": stage["restriction_count"],
            "survivors": stage["cumulative_passing_state_count"],
        }
        for stage in gate["cumulative_restriction_stages"]
    ]
    comparisons = {
        "vertex_count": theorem["rank_nine_vertex_count"],
        "lanner_placement_count": theorem["lanner_placements_per_support"],
        "supports": gate["supports"],
        "join_counts": gate["join_counts"],
        "restriction_size_histogram": {
            str(stage["order"]): stage["restriction_count"]
            for stage in gate["cumulative_restriction_stages"]
        },
        "stage_counts": archived_stages,
        "dotted_pairs_in_cluster_union": theorem[
            "dotted_pairs_in_decisive_cluster_union"
        ],
        "unrestricted_compact_coxeter_no_go": theorem[
            "unrestricted_compact_coxeter_no_go"
        ]["proved"],
    }
    for key, value in comparisons.items():
        if result[key] != value:
            raise AssertionError(f"archive disagrees on independent field {key}")
    scope = document["scope"]
    if not scope["selected_nerve_only"] or not scope["whole_graph_shadow_excluded"] is False:
        raise AssertionError("record-17 archive scope widened")
    if scope["finite_volume_noncompact_chamber_excluded"]:
        raise AssertionError("record-17 archive makes a noncompact volume claim")


def main() -> None:
    result = replay()
    compare_archive(result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
