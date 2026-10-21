#!/usr/bin/env python3
"""Exact blocker census for the ``a=0`` branch of H9/14 record 37.

The five missing A5s are rows of a 5-by-14 incidence matrix.  Local link
intervals force vertex 13 out, allow vertices 0, 1, and 11 in at most one
row, allow vertex 12 in at most two rows, and forbid any row from containing
two of those four active vertices.  The nine isolated graph vertices are
quotiented by storing their five-bit row columns as a multiset.  This script
enumerates that quotient exactly and uses integer bitsets to audit small and
five-element transversals.  It imports no SAT/SMT backend.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path
from typing import Iterator

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6
from scripts.enumerate_h9_14_four_a5_blocker import _local_interval_table


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
TARGET_INDEX = 37
TARGET_GRAPH6 = "M?????????O?O?o@?"
ROW_COUNT = 5
ACTIVE_VERTICES = (0, 1, 11, 12)
ISOLATED_VERTICES = tuple(range(2, 11))
FORBIDDEN_VERTEX = 13


def _mask(vertices) -> int:
    return sum(1 << vertex for vertex in vertices)


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if mask & (1 << vertex))


def _isolate_count_vectors(row_sums: tuple[int, ...]) -> Iterator[tuple[int, ...]]:
    """Yield all multisets of nine non-full five-bit columns with row sums."""

    patterns = tuple(sorted(range((1 << ROW_COUNT) - 1), key=lambda p: (-p.bit_count(), -p)))
    bits = tuple(
        tuple(int(bool(pattern & (1 << row))) for row in range(ROW_COUNT))
        for pattern in patterns
    )
    future_has_bit = []
    for index in range(len(patterns) + 1):
        future = patterns[index:]
        future_has_bit.append(
            tuple(any(pattern & (1 << row) for pattern in future) for row in range(ROW_COUNT))
        )
    counts = [0] * (1 << ROW_COUNT)

    def visit(index: int, columns_left: int, capacities: tuple[int, ...]):
        if any(value < 0 or value > columns_left for value in capacities):
            return
        if sum(capacities) > 4 * columns_left:
            return
        if any(
            capacities[row] and not future_has_bit[index][row]
            for row in range(ROW_COUNT)
        ):
            return
        if index == len(patterns):
            if columns_left == 0 and not any(capacities):
                yield tuple(counts)
            return
        pattern = patterns[index]
        pattern_bits = bits[index]
        maximum = columns_left
        for row, bit in enumerate(pattern_bits):
            if bit:
                maximum = min(maximum, capacities[row])
        for multiplicity in range(maximum + 1):
            counts[pattern] = multiplicity
            reduced = tuple(
                capacities[row] - multiplicity * pattern_bits[row]
                for row in range(ROW_COUNT)
            )
            yield from visit(index + 1, columns_left - multiplicity, reduced)
        counts[pattern] = 0

    yield from visit(0, len(ISOLATED_VERTICES), row_sums)


def _isolate_faces_from_counts(isolate_counts: tuple[int, ...]) -> tuple[int, ...]:
    isolate_patterns = tuple(
        pattern
        for pattern, multiplicity in enumerate(isolate_counts)
        for _ in range(multiplicity)
    )
    if len(isolate_patterns) != len(ISOLATED_VERTICES):
        raise AssertionError("isolate multiset lost a column")
    faces = []
    for row in range(ROW_COUNT):
        vertices = tuple(
            vertex
            for vertex, pattern in zip(
                ISOLATED_VERTICES, isolate_patterns, strict=True
            )
            if pattern & (1 << row)
        )
        faces.append(_mask(vertices))
    return tuple(faces)


def _cover_masks_of_sizes(
    edge_masks: tuple[int, ...], sizes: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(
        _mask(vertices)
        for size in sizes
        for vertices in combinations(range(14), size)
        if all(_mask(vertices) & edge for edge in edge_masks)
    )


def _hit_word(face: int, covers: tuple[int, ...]) -> int:
    return sum(bool(face & cover) << index for index, cover in enumerate(covers))


def enumerate_a0() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    record = survivors[TARGET_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("record-37 graph6 string changed")
    schedules = {
        schedule["missing_four_count"]: schedule
        for schedule in record["dehn_sommerville_schedule"]
    }
    schedule = schedules[0]
    if schedule["missing_five_count_offset"] != ROW_COUNT:
        raise ValueError("record-37 a=0 schedule no longer has five A5s")
    target = schedule["f_vector"][8]
    _, edges = decode_graph6(TARGET_GRAPH6)
    edge_masks = tuple(_mask(edge) for edge in edges)
    intervals, _ = _local_interval_table(TARGET_GRAPH6)
    vertex_intervals = tuple(intervals[1 << vertex] for vertex in range(14))
    if vertex_intervals != (
        (0, 1),
        (0, 1),
        *((0, 4),) * 9,
        (0, 1),
        (0, 2),
        (0, 0),
    ):
        raise ValueError("record-37 vertex A5 intervals changed")
    active_pair_masks = tuple(
        _mask(pair) for pair in combinations(ACTIVE_VERTICES, 2)
    )
    if any(
        pair not in intervals and not any(pair & edge == edge for edge in edge_masks)
        for pair in active_pair_masks
    ):
        raise AssertionError("an active pair escaped both graph and link bounds")
    for pair in active_pair_masks:
        if pair in intervals and intervals[pair] != (0, 0):
            raise AssertionError("an independent active pair can share an A5")

    small_covers = _cover_masks_of_sizes(edge_masks, (1, 2, 3, 4))
    five_covers = graph_vertex_cover_masks(TARGET_GRAPH6, 5)
    all_small = (1 << len(small_covers)) - 1
    all_five = (1 << len(five_covers)) - 1
    small_hits: dict[int, int] = {}
    five_hits: dict[int, int] = {}
    histogram: Counter[int] = Counter()
    active_type_counts: Counter[str] = Counter()
    active_type_blocker_counts: Counter[str] = Counter()
    quotient_states = 0
    distinct_row_states = 0
    blocker_states = 0
    target_state_count = 0
    target_states = []

    active_types_by_size: dict[int, list[tuple[int | None, ...]]] = {
        size: [] for size in range(ROW_COUNT + 1)
    }
    for used_0, used_1, used_11, used_12 in product(
        range(2), range(2), range(2), range(3)
    ):
        active = (
            [0] * used_0
            + [1] * used_1
            + [11] * used_11
            + [12] * used_12
        )
        if len(active) > ROW_COUNT:
            continue
        row_labels = tuple([None] * (ROW_COUNT - len(active)) + active)
        active_types_by_size[len(active)].append(row_labels)

    for active_size, active_types in active_types_by_size.items():
        if not active_types:
            continue
        row_sums = (5,) * (ROW_COUNT - active_size) + (4,) * active_size
        for isolate_counts in _isolate_count_vectors(row_sums):
            isolate_faces = _isolate_faces_from_counts(isolate_counts)
            for row_labels in active_types:
                quotient_states += 1
                type_key = ",".join(
                    "-" if item is None else str(item) for item in row_labels
                )
                active_type_counts[type_key] += 1
                faces = tuple(
                    isolate_face
                    | (0 if active is None else 1 << active)
                    for isolate_face, active in zip(
                        isolate_faces, row_labels, strict=True
                    )
                )
                if any(face.bit_count() != 5 for face in faces):
                    raise AssertionError("row compiler lost A5 cardinality")
                if any(face & (1 << FORBIDDEN_VERTEX) for face in faces):
                    raise AssertionError("forbidden vertex entered an A5")
                if any(
                    any(face & edge == edge for edge in edge_masks)
                    for face in faces
                ):
                    raise AssertionError(
                        "row compiler produced a graph-dependent A5"
                    )
                if len(set(faces)) != ROW_COUNT:
                    continue
                distinct_row_states += 1
                common_small = all_small
                for face in faces:
                    if face not in small_hits:
                        small_hits[face] = _hit_word(face, small_covers)
                    common_small &= small_hits[face]
                if common_small:
                    continue
                blocker_states += 1
                active_type_blocker_counts[type_key] += 1
                common_five = all_five
                for face in faces:
                    if face not in five_hits:
                        five_hits[face] = _hit_word(face, five_covers)
                    common_five &= five_hits[face]
                facet_count = common_five.bit_count()
                histogram[facet_count] += 1
                if facet_count == target:
                    target_state_count += 1
                if facet_count == target and len(target_states) < 100:
                    transversals = tuple(
                        cover
                        for index, cover in enumerate(five_covers)
                        if common_five & (1 << index)
                    )
                    target_states.append(
                        {
                            "row_labels": [
                                -1 if label is None else label
                                for label in row_labels
                            ],
                            "isolate_pattern_counts": list(isolate_counts),
                            "missing_five_masks": list(faces),
                            "missing_five_vertices": [
                                list(_vertices(face)) for face in faces
                            ],
                            "facet_complement_masks": list(transversals),
                        }
                    )

    return {
        "record_index": TARGET_INDEX,
        "graph6": TARGET_GRAPH6,
        "missing_edges": [list(edge) for edge in edges],
        "row_count": ROW_COUNT,
        "active_vertices": list(ACTIVE_VERTICES),
        "isolated_vertices": list(ISOLATED_VERTICES),
        "forbidden_vertex": FORBIDDEN_VERTEX,
        "vertex_a5_intervals": [list(item) for item in vertex_intervals],
        "small_graph_cover_count": len(small_covers),
        "five_graph_cover_count": len(five_covers),
        "row-normalized-isolate-quotient-states": quotient_states,
        "distinct-row-states": distinct_row_states,
        "blocker-valid-states": blocker_states,
        "active-type-state-counts": dict(sorted(active_type_counts.items())),
        "active-type-blocker-counts": dict(
            sorted(active_type_blocker_counts.items())
        ),
        "required-facet-count": target,
        "facet-count-histogram": dict(sorted(histogram.items())),
        "target-state-count": target_state_count,
        "target-state-count-stored": len(target_states),
        "target-states": target_states,
    }


def main() -> None:
    print(json.dumps(enumerate_a0(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
