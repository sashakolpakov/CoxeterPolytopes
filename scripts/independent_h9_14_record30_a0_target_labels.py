#!/usr/bin/env python3
"""Standalone compact Coxeter-label consumer for Record-30 ``a=0`` targets.

The input may be the connected-ridge census JSON (with
``target_orbit_representatives``) or a single audited candidate JSON (with
``missing_A5s``).  Every target must consist of the fixed three missing edges
and eight missing five-sets.  This script uses only the Python standard
library.  It reconstructs the literal 132 rank-nine vertices, regenerates all
420 labelled order-five Lanner placements, and performs a lossless sparse
join with every actual-vertex ellipticity restriction.

A zero terminal state count is an unrestricted compact-label no-go in the
finite alphabet {2,3,4,5}: omitted meeting labels were tested at their minimal
value 2, and Perron--Frobenius monotonicity prevents larger labels from
repairing a nonelliptic principal diagram.  A nonzero terminal count is only
a necessary-condition survivor.  Ideal faces and finite-volume noncompact
chambers are not tested.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys
from typing import Iterable, Mapping, Sequence


VERTEX_COUNT = 14
DIMENSION = 9
FULL_MASK = (1 << VERTEX_COUNT) - 1
EDGES = ((0, 11), (1, 12), (2, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in EDGES)
EXPECTED_A5_COUNT = 8
EXPECTED_ACTUAL_VERTEX_COUNT = 132

Pair = tuple[int, int]
LabelState = tuple[tuple[Pair, int], ...]

# The five classical order-five Lanner diagrams.  Missing pairs have label 2.
LANNER_REPRESENTATIVES = (
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
)


def _pair(left: int, right: int) -> Pair:
    if left == right:
        raise ValueError("a Coxeter pair needs two distinct facets")
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
    code: tuple[int, ...], permutation: Sequence[int]
) -> tuple[int, ...]:
    pairs = tuple(combinations(range(5), 2))
    source = dict(zip(pairs, code, strict=True))
    moved = {
        _pair(permutation[left], permutation[right]): label
        for (left, right), label in source.items()
    }
    return tuple(moved[edge] for edge in pairs)


@lru_cache(maxsize=1)
def _lanner_codes() -> tuple[tuple[int, ...], ...]:
    pairs = tuple(combinations(range(5), 2))
    placements = set()
    for representative in LANNER_REPRESENTATIVES:
        labels = {edge: 2 for edge in pairs}
        labels.update(representative)
        code = tuple(labels[edge] for edge in pairs)
        for permutation in permutations(range(5)):
            placements.add(_permuted_code(code, permutation))
    canonical_types = {
        min(_permuted_code(code, permutation) for permutation in permutations(range(5)))
        for code in placements
    }
    if (len(placements), len(canonical_types)) != (420, 5):
        raise AssertionError("the complete order-five Lanner table changed")
    return tuple(sorted(placements))


def _path_labels(
    component: set[int],
    adjacency: Mapping[int, tuple[tuple[int, int], ...]],
) -> tuple[int, ...]:
    current = min(
        vertex for vertex in component if len(adjacency[vertex]) == 1
    )
    previous: int | None = None
    labels = []
    while True:
        following = [row for row in adjacency[current] if row[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("finite Coxeter path classifier branched")
        other, label = following[0]
        labels.append(label)
        previous, current = current, other
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
        for left, right in combinations(component, 2)
        if labels[_pair(left, right)] > 2
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
            current = frontier.pop()
            neighbors = {
                other
                for other in unseen
                if labels[_pair(current, other)] > 2
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
        minimal = True
        for vertex in _vertices(candidate):
            smaller = candidate ^ (1 << vertex)
            if all(smaller & face for face in clutter):
                minimal = False
                break
        if minimal:
            answer.append(candidate)
    return tuple(answer)


def _literal_actual_vertices(
    supports: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(map(frozenset, (*EDGES, *supports)))
    answer = []
    for chosen in combinations(range(VERTEX_COUNT), DIMENSION):
        chosen_set = set(chosen)
        if not any(face <= chosen_set for face in forbidden):
            answer.append(chosen)
    return tuple(answer)


@lru_cache(maxsize=1024)
def _placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    edges = tuple(combinations(sorted(support), 2))
    return tuple(
        tuple(zip(edges, code, strict=True)) for code in _lanner_codes()
    )


def _state_digest(states: Sequence[LabelState]) -> str:
    return _canonical_sha(
        [
            [[left, right, label] for (left, right), label in state]
            for state in sorted(states)
        ]
    )


def _normalize_face(raw: object, *, target_index: int) -> int:
    if isinstance(raw, int):
        value = raw
    elif isinstance(raw, list) and all(isinstance(vertex, int) for vertex in raw):
        if len(set(raw)) != len(raw):
            raise ValueError(f"target {target_index} has a repeated facet")
        if any(vertex < 0 or vertex >= VERTEX_COUNT for vertex in raw):
            raise ValueError(f"target {target_index} has an out-of-range facet")
        value = _mask(raw)
    else:
        raise ValueError(f"target {target_index} has a malformed A5 support")
    if value <= 0 or value > FULL_MASK or value.bit_count() != 5:
        raise ValueError(f"target {target_index} has a non-five-element support")
    if any(value & edge == edge for edge in EDGE_MASKS):
        raise ValueError(f"target {target_index} A5 contains a fixed missing edge")
    return value


def _target_families(document: object) -> tuple[tuple[int, ...], ...]:
    declared_count: int | None = None
    if isinstance(document, Mapping):
        if "target_orbit_representatives" in document:
            raw_families = document["target_orbit_representatives"]
            if "target_orbit_count" in document:
                declared_count = int(document["target_orbit_count"])
        elif "missing_A5s" in document:
            raw_families = [document["missing_A5s"]]
        else:
            raise ValueError(
                "input needs target_orbit_representatives or missing_A5s"
            )
    elif isinstance(document, list):
        raw_families = document
    else:
        raise ValueError("input JSON must be an object or a target-family list")
    if not isinstance(raw_families, list):
        raise ValueError("target family collection is malformed")
    families = []
    for target_index, raw_family in enumerate(raw_families):
        if not isinstance(raw_family, list) or len(raw_family) != EXPECTED_A5_COUNT:
            raise ValueError(
                f"target {target_index} must have exactly eight A5 supports"
            )
        family = tuple(
            sorted(
                _normalize_face(face, target_index=target_index)
                for face in raw_family
            )
        )
        if len(set(family)) != EXPECTED_A5_COUNT:
            raise ValueError(f"target {target_index} repeats an A5 support")
        families.append(family)
    if len(set(families)) != len(families):
        raise ValueError("the input repeats a target representative literally")
    if declared_count is not None and declared_count != len(families):
        raise ValueError(
            "declared target_orbit_count does not match the representative list"
        )
    return tuple(families)


def _extend_states(
    states: Sequence[LabelState],
    support: Sequence[int],
    union: set[int],
    actual_vertices: Sequence[Sequence[int]],
) -> tuple[tuple[LabelState, ...], int, tuple[tuple[int, ...], ...]]:
    cluster = tuple(sorted(union | set(support)))
    restrictions = tuple(
        sorted(
            {
                tuple(vertex for vertex in actual if vertex in cluster)
                for actual in actual_vertices
            },
            key=lambda row: (len(row), row),
        )
    )
    passing: set[LabelState] = set()
    compatible = 0
    options = _placement_options(tuple(support))
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
    return tuple(sorted(passing)), compatible, restrictions


def _audit_target(index: int, family: Sequence[int]) -> dict[str, object]:
    supports = tuple(_vertices(value) for value in family)
    clutter = EDGE_MASKS + tuple(family)
    blocker = _minimal_transversals(clutter)
    histogram = Counter(value.bit_count() for value in blocker)
    if histogram != {5: EXPECTED_ACTUAL_VERTEX_COUNT}:
        raise ValueError(
            f"target {index} is not a pure 5^{EXPECTED_ACTUAL_VERTEX_COUNT} schedule: "
            f"{dict(sorted(histogram.items()))}"
        )
    actual_vertices = _literal_actual_vertices(supports)
    blocker_vertices = tuple(
        sorted(_vertices(FULL_MASK ^ complement) for complement in blocker)
    )
    if blocker_vertices != actual_vertices:
        raise AssertionError(f"target {index} blocker/vertex duality failed")
    if len(actual_vertices) != EXPECTED_ACTUAL_VERTEX_COUNT:
        raise AssertionError(f"target {index} lost a literal actual vertex")

    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[LabelState, ...] = ((),)
    stages = []
    while states and remaining:
        if not union and len(remaining) > 1:
            first_pair = max(
                combinations(sorted(remaining), 2),
                key=lambda items: (
                    len(set(supports[items[0]]) & set(supports[items[1]])),
                    -items[0],
                    -items[1],
                ),
            )
            support_index = first_pair[0]
        elif not union:
            support_index = min(remaining)
        else:
            support_index = max(
                remaining,
                key=lambda item: (
                    len(set(supports[item]) & union),
                    -len(set(supports[item]) - union),
                    -item,
                ),
            )
        support = supports[support_index]
        old_count = len(states)
        states, compatible, restrictions = _extend_states(
            states, support, union, actual_vertices
        )
        union.update(support)
        remaining.remove(support_index)
        stages.append(
            {
                "support_index": support_index,
                "support_mask": family[support_index],
                "support": list(support),
                "option_count": len(_lanner_codes()),
                "raw_extension_count": old_count * len(_lanner_codes()),
                "compatible_extension_count": compatible,
                "restriction_count": len(restrictions),
                "restriction_size_histogram": {
                    str(size): count
                    for size, count in sorted(
                        Counter(map(len, restrictions)).items()
                    )
                },
                "passing_state_count": len(states),
                "passing_state_sha256": _state_digest(states),
            }
        )

    result = {
        "target_index": index,
        "A5_masks": list(family),
        "A5_supports": [list(support) for support in supports],
        "minimal_transversal_histogram": {
            str(size): count for size, count in sorted(histogram.items())
        },
        "facet_complement_sha256": sha256(
            b"".join(value.to_bytes(2, "little") for value in blocker)
        ).hexdigest(),
        "actual_rank_nine_vertex_count": len(actual_vertices),
        "actual_vertices_sha256": _canonical_sha(actual_vertices),
        "A5_pair_intersection_size_histogram": {
            str(size): count
            for size, count in sorted(
                Counter(
                    len(set(left) & set(right))
                    for left, right in combinations(supports, 2)
                ).items()
            )
        },
        "stages": stages,
        "processed_support_indices": [row["support_index"] for row in stages],
        "unprocessed_support_indices": sorted(remaining),
        "terminal_sparse_state_count": len(states),
        "terminal_sparse_states_sha256": _state_digest(states),
        "unrestricted_compact_label_no_go": not states,
        "decisive_processed_support_masks": (
            [family[int(row["support_index"])] for row in stages]
            if not states
            else []
        ),
    }
    return result


def replay(document: object, *, progress: bool = False) -> dict[str, object]:
    families = _target_families(document)
    audits = []
    for index, family in enumerate(families):
        if progress:
            print(
                f"record30 a=0 label target {index + 1}/{len(families)}",
                file=sys.stderr,
                flush=True,
            )
        audits.append(_audit_target(index, family))
    histogram = Counter(
        int(row["terminal_sparse_state_count"]) for row in audits
    )
    survivors = [
        {
            "target_index": int(row["target_index"]),
            "terminal_sparse_state_count": int(
                row["terminal_sparse_state_count"]
            ),
        }
        for row in audits
        if int(row["terminal_sparse_state_count"]) > 0
    ]
    result = {
        "kind": "h9_14_record30_a0_all_target_sparse_label_replay",
        "schema_version": "1.0.0",
        "compact_only": True,
        "imports_coxeter_polytopes": False,
        "fixed_missing_edges": [list(edge) for edge in EDGES],
        "A5_support_count_per_target": EXPECTED_A5_COUNT,
        "actual_rank_nine_vertex_count_per_target": (
            EXPECTED_ACTUAL_VERTEX_COUNT
        ),
        "complete_order_five_lanner_placement_count": len(_lanner_codes()),
        "target_orbit_count": len(families),
        "target_family_sha256": _canonical_sha(families),
        "target_audits": audits,
        "terminal_sparse_state_count_histogram": {
            str(count): multiplicity
            for count, multiplicity in sorted(histogram.items())
        },
        "unrestricted_compact_label_no_go_target_count": sum(
            bool(row["unrestricted_compact_label_no_go"]) for row in audits
        ),
        "sparse_label_survivor_count": len(survivors),
        "sparse_label_survivors": survivors,
        "all_nonempty_input_targets_have_unrestricted_compact_label_no_go": (
            bool(audits) and not survivors
        ),
        "finite_volume_noncompact_chamber_excluded": False,
        "new_hyperbolic_coxeter_polytope": False,
        "scope": (
            "Standalone complete order-five Lanner-placement and literal "
            "actual-rank-nine-vertex sparse join for the supplied compact "
            "Record-30 a=0 target representatives. A zero state count is "
            "unrestricted within finite labels 2,3,4,5 by Perron--Frobenius "
            "monotonicity. Surviving sparse states require further whole-nerve "
            "metric work. No ideal-face or finite-volume inference."
        ),
    }
    result["summary_sha256"] = _canonical_sha(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--pretty", action="store_true")
    parser.add_argument("--progress", action="store_true")
    args = parser.parse_args()
    document = json.loads(args.input.read_text(encoding="utf-8"))
    print(
        json.dumps(
            replay(document, progress=args.progress),
            indent=2 if args.pretty else None,
            sort_keys=True,
            separators=None if args.pretty else (",", ":"),
        )
    )


if __name__ == "__main__":
    main()
