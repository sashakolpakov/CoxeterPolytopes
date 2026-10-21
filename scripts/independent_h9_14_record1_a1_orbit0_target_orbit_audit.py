#!/usr/bin/env python3
"""Hostile replay of the first Record-1 ``a=1`` target orbit.

This file imports neither the production open-shadow completion model nor its
orbit certifier.  It reconstructs the three-edge missing graph, the fixed
four-nonface, the representative nineteen-row completion, the full
fixed-source automorphism group, and the literal blocker/ridge conditions.
The compact-label replay is the separately written classical-Lanner audit.

The result is deliberately bounded: it excludes one 60,480-element *labelled*
target orbit.  It does not claim that this is the only incidence orbit for the
source.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from typing import Iterable, Iterator, Sequence

from scripts.independent_h9_14_record1_a1_target0_label_audit import (
    ADDITIONS,
    replay as replay_labels,
)


N = 14
FULL = (1 << N) - 1
EDGE_ROWS = ((0, 13), (1, 13), (2, 13))
SOURCE_A4 = (0, 3, 4, 5)
LEAF_BLOCK = (1, 2)
SOURCE_ISOLATE_BLOCK = (3, 4, 5)
OUTSIDE_ISOLATE_BLOCK = (6, 7, 8, 9, 10, 11, 12)


def _mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def _vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value >> vertex & 1)


def _move(value: int, action: Sequence[int]) -> int:
    return _mask(action[vertex] for vertex in _vertices(value))


def _actions() -> Iterator[tuple[int, ...]]:
    """Enumerate ``S2 x S3 x S7`` without importing production symmetry code."""

    for leaves in permutations(LEAF_BLOCK):
        for source_isolates in permutations(SOURCE_ISOLATE_BLOCK):
            for outside_isolates in permutations(OUTSIDE_ISOLATE_BLOCK):
                action = list(range(N))
                action[1], action[2] = leaves
                action[3:6] = source_isolates
                action[6:13] = outside_isolates
                yield tuple(action)


def _hits_all(value: int, rows: Sequence[int]) -> bool:
    return all(value & row for row in rows)


def _minimal_transversals(rows: Sequence[int]) -> tuple[int, ...]:
    answer = []
    for value in range(1, FULL + 1):
        if not _hits_all(value, rows):
            continue
        if all(
            not _hits_all(value ^ (1 << vertex), rows)
            for vertex in _vertices(value)
        ):
            answer.append(value)
    return tuple(answer)


def _family_bytes(family: Sequence[int]) -> bytes:
    return (",".join(map(str, sorted(family))) + "\n").encode("ascii")


def replay() -> dict[str, object]:
    edge_masks = tuple(_mask(edge) for edge in EDGE_ROWS)
    source_mask = _mask(SOURCE_A4)
    target = tuple(sorted(_mask(face) for face in ADDITIONS))

    # The graph/source incidence signatures force precisely these five blocks:
    # the degree-three centre, the distinguished leaf in A4, the other two
    # leaves, the three A4 isolates, and the seven remaining isolates.
    graph_degrees = tuple(
        sum(vertex in edge for edge in EDGE_ROWS) for vertex in range(N)
    )
    source_partition = (
        (13,),
        (0,),
        LEAF_BLOCK,
        SOURCE_ISOLATE_BLOCK,
        OUTSIDE_ISOLATE_BLOCK,
    )
    signatures = tuple(
        (graph_degrees[vertex], int(vertex in SOURCE_A4))
        for vertex in range(N)
    )
    signature_classes = tuple(
        tuple(
            vertex for vertex, signature in enumerate(signatures)
            if signature == value
        )
        for value in sorted(set(signatures), reverse=True)
    )
    if set(signature_classes) != set(source_partition):
        raise AssertionError("independent graph/source partition changed")

    action_count = 0
    fixed_target_actions = 0
    stream = sha256()
    edge_set = frozenset(edge_masks)
    for action in _actions():
        action_count += 1
        if frozenset(_move(edge, action) for edge in edge_masks) != edge_set:
            raise AssertionError("derived action does not preserve the missing graph")
        if _move(source_mask, action) != source_mask:
            raise AssertionError("derived action does not preserve the A4 source")
        moved = tuple(sorted(_move(face, action) for face in target))
        fixed_target_actions += int(moved == target)
        stream.update(_family_bytes(moved))
    if action_count != 60_480 or fixed_target_actions != 1:
        raise AssertionError("independent target orbit-stabilizer count changed")

    minimal_nonfaces = (*edge_masks, source_mask, *target)
    blockers = _minimal_transversals(minimal_nonfaces)
    blocker_histogram = Counter(value.bit_count() for value in blockers)
    if blocker_histogram != {5: 158}:
        raise AssertionError("representative is not a literal pure 158-facet target")
    ridge_histogram: Counter[int] = Counter()
    for six in combinations(range(N), 6):
        six_mask = _mask(six)
        degree = sum(facet & six_mask == facet for facet in blockers)
        ridge_histogram[degree] += 1
    if any(degree not in (0, 2) for degree in ridge_histogram):
        raise AssertionError("representative violates the literal ridge equation")

    labels = replay_labels()
    if not labels["unrestricted_compact_label_no_go"]:
        raise AssertionError("independent compact-label obstruction disappeared")

    return {
        "kind": "independent_h9_14_record1_a1_orbit0_target_orbit_audit",
        "record_index": 1,
        "schedule_a": 1,
        "source_orbit_index": 0,
        "missing_edges": [list(edge) for edge in EDGE_ROWS],
        "representative_A4": list(SOURCE_A4),
        "source_stabilizer_partition": [list(block) for block in source_partition],
        "source_stabilizer_structure": "S2 x S3 x S7",
        "source_stabilizer_order": action_count,
        "representative_target_stabilizer_order": fixed_target_actions,
        "labelled_target_orbit_size": action_count // fixed_target_actions,
        "action_image_stream_sha256": stream.hexdigest(),
        "representative_additional_A5s": [list(face) for face in ADDITIONS],
        "representative_additional_A5_sha256": sha256(
            _family_bytes(target)
        ).hexdigest(),
        "representative_minimal_transversal_size_histogram": {
            str(size): count for size, count in sorted(blocker_histogram.items())
        },
        "representative_ridge_degree_histogram": {
            str(degree): count for degree, count in sorted(ridge_histogram.items())
        },
        "representative_compact_label_stage_counts": [
            stage["passing_state_count"] for stage in labels["stages"]
        ],
        "representative_compact_label_state_hashes": [
            stage["passing_state_sha256"] for stage in labels["stages"]
        ],
        "whole_labelled_target_orbit_compact_label_no_go": True,
        "quotient_incidence_census_complete": False,
        "whole_source_orbit_compact_no_go": False,
        "finite_volume_noncompact_excluded": False,
        "scope": (
            "Orbit-stabilizer and label-equivariance exclude one 60,480-member "
            "labelled incidence orbit; other incidence orbits for this source "
            "remain open."
        ),
    }


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
