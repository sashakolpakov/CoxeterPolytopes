"""Exact compressed blocker census for H9/14 ``a=0``, four-A5 branches.

The four missing five-faces are represented by four incidence rows.  Graph
vertices on an edge are kept as distinguished columns; the isolated graph
vertices are compressed to multiplicities of their four-bit membership
columns.  Row permutations are quotiented exactly.  Every vertex, edge, and
triple strict-link/Macaulay interval is checked before blocker duality.

This is a self-contained finite enumerator.  It imports no SAT/SMT backend.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from math import comb
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
)
from coxeter_polytopes.combinatorics.missing_edges import (
    decode_graph6,
    macaulay_raise,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
ROW_COUNT = 4
ROW_PERMUTATIONS = tuple(permutations(range(ROW_COUNT)))


def _mask(vertices: Sequence[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if mask & (1 << vertex))


def _permute_pattern(pattern: int, permutation: Sequence[int]) -> int:
    return sum(
        bool(pattern & (1 << permutation[new_row])) << new_row
        for new_row in range(ROW_COUNT)
    )


ROW_PATTERN_IMAGES = tuple(
    tuple(_permute_pattern(pattern, permutation) for pattern in range(1 << ROW_COUNT))
    for permutation in ROW_PERMUTATIONS
)
ROW_PATTERN_PREIMAGES = tuple(
    tuple(images.index(pattern) for pattern in range(1 << ROW_COUNT))
    for images in ROW_PATTERN_IMAGES
)


def _canonical_state_key(
    distinguished_columns: tuple[int, ...],
    isolate_counts: tuple[int, ...],
) -> tuple[int, ...]:
    return min(
        tuple(images[pattern] for pattern in distinguished_columns)
        + tuple(isolate_counts[preimages[pattern]] for pattern in range(1 << ROW_COUNT))
        for images, preimages in zip(
            ROW_PATTERN_IMAGES, ROW_PATTERN_PREIMAGES, strict=True
        )
    )


def _is_independent(vertices: Sequence[int], edge_set: set[tuple[int, int]]) -> bool:
    return all(tuple(sorted(pair)) not in edge_set for pair in combinations(vertices, 2))


def _local_a5_interval(
    subset: tuple[int, ...],
    edge_set: set[tuple[int, int]],
) -> tuple[int, int, int, int]:
    """Return ``previous, baseline high, lower, upper`` A5 incidence."""

    size = len(subset)
    link_dimension = 9 - size
    high_order = 5 - size
    rest = tuple(vertex for vertex in range(14) if vertex not in subset)
    link_f = {-1: 1}
    for order in range(1, high_order + 1):
        link_f[order - 1] = sum(
            _is_independent(subset + extension, edge_set)
            for extension in combinations(rest, order)
        )
    h = tuple(
        sum(
            (-1) ** (index - order)
            * comb(link_dimension - order, link_dimension - index)
            * link_f[order - 1]
            for order in range(index + 1)
        )
        for index in range(high_order + 1)
    )
    g = tuple(h[index] - h[index - 1] for index in range(1, len(h)))
    previous = g[-2]
    baseline_high = g[-1]
    lower = max(
        0,
        baseline_high - macaulay_raise(previous, high_order - 1),
    )
    upper = baseline_high - 1
    return previous, baseline_high, lower, upper


def _local_interval_table(
    graph6: str,
) -> tuple[dict[int, tuple[int, int]], set[tuple[int, int]]]:
    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != 14:
        raise ValueError("four-A5 census is specific to fourteen vertices")
    edge_set = {tuple(sorted(edge)) for edge in edges}
    intervals: dict[int, tuple[int, int]] = {}
    for size in (1, 2, 3):
        for subset in combinations(range(vertex_count), size):
            if not _is_independent(subset, edge_set):
                continue
            previous, baseline, lower, upper = _local_a5_interval(
                subset, edge_set
            )
            if previous < 1 or baseline < 1 or lower > upper:
                raise AssertionError("input graph failed its strict local precheck")
            intervals[_mask(subset)] = (lower, upper)
    return intervals, edge_set


def _active_assignments(
    distinguished: tuple[int, ...],
    intervals: Mapping[int, tuple[int, int]],
    edge_set: set[tuple[int, int]],
) -> dict[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    choices = []
    for vertex in distinguished:
        lower, upper = intervals[1 << vertex]
        choices.append(
            tuple(
                pattern
                for pattern in range(1 << ROW_COUNT)
                if lower <= pattern.bit_count() <= upper
            )
        )
    by_row_counts: dict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
    for columns in product(*choices):
        valid = True
        for row in range(ROW_COUNT):
            row_vertices = tuple(
                vertex
                for vertex, column in zip(distinguished, columns, strict=True)
                if column & (1 << row)
            )
            if not _is_independent(row_vertices, edge_set):
                valid = False
                break
        if not valid:
            continue
        for size in (2, 3):
            for indices in combinations(range(len(distinguished)), size):
                subset_vertices = tuple(distinguished[index] for index in indices)
                subset_mask = _mask(subset_vertices)
                if subset_mask not in intervals:
                    continue
                incidence = (1 << ROW_COUNT) - 1
                for index in indices:
                    incidence &= columns[index]
                lower, upper = intervals[subset_mask]
                if not lower <= incidence.bit_count() <= upper:
                    valid = False
                    break
            if not valid:
                break
        if not valid:
            continue
        row_counts = tuple(
            sum(bool(column & (1 << row)) for column in columns)
            for row in range(ROW_COUNT)
        )
        by_row_counts[row_counts].append(tuple(columns))
    return {
        key: tuple(values) for key, values in sorted(by_row_counts.items())
    }


def _full_local_intervals_hold(
    vertices: tuple[int, ...],
    columns: tuple[int, ...],
    intervals: Mapping[int, tuple[int, int]],
) -> bool:
    for size in (1, 2, 3):
        for indices in combinations(range(len(vertices)), size):
            subset_mask = _mask(tuple(vertices[index] for index in indices))
            if subset_mask not in intervals:
                # A graph edge in the subset means no independent A5 can
                # contain it; the row construction is checked separately.
                incidence = (1 << ROW_COUNT) - 1
                for index in indices:
                    incidence &= columns[index]
                if incidence:
                    return False
                continue
            incidence = (1 << ROW_COUNT) - 1
            for index in indices:
                incidence &= columns[index]
            lower, upper = intervals[subset_mask]
            if not lower <= incidence.bit_count() <= upper:
                return False
    return True


def _faces_from_columns(
    vertices: tuple[int, ...], columns: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(
        _mask(
            tuple(
                vertex
                for vertex, column in zip(vertices, columns, strict=True)
                if column & (1 << row)
            )
        )
        for row in range(ROW_COUNT)
    )


def _has_small_transversal(edge_masks: tuple[int, ...], faces: tuple[int, ...]) -> bool:
    clutter = edge_masks + faces
    for size in range(1, 5):
        for subset in combinations(range(14), size):
            candidate = _mask(subset)
            if all(candidate & item for item in clutter):
                return True
    return False


def _five_transversals(
    graph_cover_fives: tuple[int, ...], faces: tuple[int, ...]
) -> tuple[int, ...]:
    return tuple(
        cover
        for cover in graph_cover_fives
        if all(cover & face for face in faces)
    )


def enumerate_record(index: int) -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    record = survivors[index]
    graph6 = record["graph6"]
    schedules = {
        schedule["missing_four_count"]: schedule
        for schedule in record["dehn_sommerville_schedule"]
    }
    if schedules[0]["missing_five_count_offset"] != ROW_COUNT:
        raise ValueError("target does not have an a=0 four-A5 branch")
    target = schedules[0]["f_vector"][8]
    intervals, edge_set = _local_interval_table(graph6)
    if any(lower != 0 for lower, _ in intervals.values()):
        raise ValueError(
            "this compiler currently requires zero lower A5 incidences"
        )
    degrees = Counter(vertex for edge in edge_set for vertex in edge)
    isolates = tuple(vertex for vertex in range(14) if degrees[vertex] == 0)
    distinguished = tuple(
        vertex
        for vertex in range(14)
        if degrees[vertex] > 0 and intervals[1 << vertex][1] >= 1
    )
    forbidden = tuple(
        vertex
        for vertex in range(14)
        if degrees[vertex] > 0 and intervals[1 << vertex][1] == 0
    )
    if not isolates:
        raise ValueError("compressed census requires a nonempty isolate block")
    isolate_interval = {intervals[1 << vertex] for vertex in isolates}
    if len(isolate_interval) != 1:
        raise AssertionError("isolated vertices do not have a common interval")
    isolate_lower, isolate_upper = isolate_interval.pop()
    if isolate_lower != 0:
        raise ValueError("positive isolate lower degrees need a different compiler")

    active_by_counts = _active_assignments(distinguished, intervals, edge_set)
    allowed_isolate_patterns = tuple(
        pattern
        for pattern in range(1 << ROW_COUNT)
        if pattern.bit_count() <= isolate_upper
    )
    edge_masks = tuple(_mask(edge) for edge in edge_set)
    graph_cover_fives = graph_vertex_cover_masks(graph6, 5)
    canonical_states: set[tuple[int, ...]] = set()
    facet_histogram: Counter[int] = Counter()
    target_states = []
    raw_join_count = 0
    local_state_count = 0
    distinct_row_state_count = 0
    blocker_state_count = 0

    for isolate_columns in combinations_with_replacement(
        allowed_isolate_patterns, len(isolates)
    ):
        isolate_row_counts = tuple(
            sum(bool(column & (1 << row)) for column in isolate_columns)
            for row in range(ROW_COUNT)
        )
        required_active = tuple(5 - value for value in isolate_row_counts)
        if any(value < 0 for value in required_active):
            continue
        for distinguished_columns in active_by_counts.get(required_active, ()):
            raw_join_count += 1
            vertices = distinguished + isolates
            columns = distinguished_columns + isolate_columns
            if not _full_local_intervals_hold(vertices, columns, intervals):
                continue
            local_state_count += 1
            faces = _faces_from_columns(vertices, columns)
            if any(face.bit_count() != 5 for face in faces):
                raise AssertionError("column compiler lost an A5 row sum")
            if len(set(faces)) != ROW_COUNT:
                continue
            distinct_row_state_count += 1
            isolate_counts = tuple(
                isolate_columns.count(pattern) for pattern in range(1 << ROW_COUNT)
            )
            key = _canonical_state_key(distinguished_columns, isolate_counts)
            if key in canonical_states:
                continue
            canonical_states.add(key)
            if _has_small_transversal(edge_masks, faces):
                continue
            blocker_state_count += 1
            transversals = _five_transversals(graph_cover_fives, faces)
            facet_histogram[len(transversals)] += 1
            if len(transversals) == target:
                target_states.append(
                    {
                        "missing_five_masks": list(faces),
                        "missing_five_vertices": [
                            list(_vertices(face)) for face in faces
                        ],
                        "facet_complement_masks": list(transversals),
                    }
                )

    return {
        "record_index": index,
        "graph6": graph6,
        "missing_edges": [list(edge) for edge in sorted(edge_set)],
        "distinguished_allowed_vertices": list(distinguished),
        "forbidden_nonisolated_vertices": list(forbidden),
        "isolated_vertices": list(isolates),
        "local_a5_intervals": {
            str(mask): list(interval) for mask, interval in sorted(intervals.items())
        },
        "active_assignment_row_count_cells": len(active_by_counts),
        "raw_active_isolate_joins": raw_join_count,
        "full_local_states_before_row_quotient": local_state_count,
        "distinct_row_states_before_row_quotient": distinct_row_state_count,
        "canonical_local_states": len(canonical_states),
        "blocker_valid_states": blocker_state_count,
        "required_facet_count": target,
        "facet_count_histogram": dict(sorted(facet_histogram.items())),
        "target_state_count": len(target_states),
        "target_states": target_states,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("index", type=int)
    args = parser.parse_args(argv)
    print(json.dumps(enumerate_record(args.index), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
