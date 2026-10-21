#!/usr/bin/env python3
"""Exact compressed blocker exploration for H9/14 survivor record 43.

The missing-edge graph is the three-step Ferrers graph

    0--11,12,13;  1--12,13;  2--13

plus eight isolates.  This exploratory checker contains no SAT/SMT code.  Its
``a=0`` census quotients only the permutation symmetry of the isolates and
the unorderedness of the four missing five-faces; duplicated representatives
are harmless because the theorem-facing question is whether a facet count is
present in the exact support.
"""

from __future__ import annotations

import json
from collections import Counter
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from pathlib import Path

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_record54 import _base_link_g3_g4
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
    minimal_graph_vertex_covers,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


GRAPH6 = "M?????????O?o?w??"
ACTIVE = (0, 1, 2, 11, 12, 13)
ISOLATES = tuple(range(3, 11))
TARGET_A0_FACETS = 90


def mask(vertices) -> int:
    return sum(1 << vertex for vertex in vertices)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if value & (1 << vertex))


def is_independent_active(value: int, edge_masks: tuple[int, ...]) -> bool:
    return all(value & edge != edge for edge in edge_masks)


def has_small_transversal(
    faces: tuple[int, ...], covers: tuple[int, ...]
) -> bool:
    """Whether a size-three graph cover extends within the size-four budget."""

    for cover in covers:
        missed = tuple(face for face in faces if not face & cover)
        if not missed:
            return True
        common = missed[0]
        for face in missed[1:]:
            common &= face
        if common:
            return True
    return False


@lru_cache(maxsize=None)
def column_distributions(
    allowed: tuple[int, ...], row_sums: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    """Multiplicity vectors for eight indistinguishable isolate columns."""

    counts = [0] * len(allowed)
    answers: list[tuple[int, ...]] = []

    def recurse(position: int, left: int, sums: tuple[int, ...]) -> None:
        if position == len(allowed):
            if left == 0 and not any(sums):
                answers.append(tuple(counts))
            return
        pattern = allowed[position]
        maximum = left
        for row in range(4):
            if pattern & (1 << row):
                maximum = min(maximum, sums[row])
        for multiplicity in range(maximum + 1):
            next_left = left - multiplicity
            next_sums = tuple(
                sums[row] - multiplicity * bool(pattern & (1 << row))
                for row in range(4)
            )
            if all(0 <= value <= next_left for value in next_sums):
                counts[position] = multiplicity
                recurse(position + 1, next_left, next_sums)
        counts[position] = 0

    recurse(0, len(ISOLATES), row_sums)
    return tuple(answers)


def a0_census() -> None:
    order, edges = decode_graph6(GRAPH6)
    if order != 14:
        raise AssertionError("wrong graph order")
    edge_masks = tuple(mask(edge) for edge in edges)
    active_mask = mask(ACTIVE)
    active_edges = tuple(edge for edge in edge_masks if edge & active_mask == edge)
    covers = minimal_graph_vertex_covers(GRAPH6)
    if Counter(cover.bit_count() for cover in covers) != {3: 4}:
        raise AssertionError("record 43 should have four size-three graph covers")
    base_g3, base_g4 = _base_link_g3_g4(GRAPH6)
    if base_g3 != base_g4:
        raise AssertionError("unexpected unequal base link arrays")
    caps = tuple(value - 1 for value in base_g4)

    independent_active = tuple(
        value
        for size in range(5)
        for chosen in combinations(ACTIVE, size)
        for value in (mask(chosen),)
        if is_independent_active(value, active_edges)
        and all(bool(value & (1 << vertex)) <= caps[vertex] for vertex in ACTIVE)
    )

    active_configs = []
    for config in combinations_with_replacement(independent_active, 4):
        if any(
            sum(bool(row & (1 << vertex)) for row in config) > caps[vertex]
            for vertex in ACTIVE
        ):
            continue
        missed_rows = tuple(
            sum(
                1 << row_index
                for row_index, row in enumerate(config)
                if not row & cover
            )
            for cover in covers
        )
        # A graph cover plus one active vertex must not hit every row.
        if any(
            all(
                not (pattern & (1 << row_index))
                or config[row_index] & (1 << vertex)
                for row_index in range(4)
            )
            for pattern in missed_rows
            for vertex in ACTIVE
        ):
            continue
        active_configs.append((config, missed_rows))

    cover_fives = graph_vertex_cover_masks(GRAPH6, 5)
    histogram: Counter[int] = Counter()
    target_witnesses: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    states = 0
    bad_replays = 0
    for config, missed_rows in active_configs:
        allowed = tuple(
            column
            for column in range(1 << 4)
            if column.bit_count() <= 3
            and all(column & rows != rows for rows in missed_rows)
        )
        row_sums = tuple(5 - row.bit_count() for row in config)
        for multiplicities in column_distributions(allowed, row_sums):
            columns = tuple(
                column
                for column, multiplicity in zip(allowed, multiplicities)
                for _ in range(multiplicity)
            )
            if len(columns) != len(ISOLATES):
                raise AssertionError("wrong isolate-column total")
            faces = list(config)
            for isolate, column in zip(ISOLATES, columns):
                for row in range(4):
                    if column & (1 << row):
                        faces[row] |= 1 << isolate
            face_tuple = tuple(faces)
            if any(face.bit_count() != 5 for face in face_tuple):
                raise AssertionError("row-size error")
            if len(set(face_tuple)) != 4:
                continue
            states += 1
            if has_small_transversal(face_tuple, covers):
                bad_replays += 1
                continue
            facet_count = sum(
                all(transversal & face for face in face_tuple)
                for transversal in cover_fives
            )
            histogram[facet_count] += 1
            if facet_count == TARGET_A0_FACETS:
                selected = tuple(
                    transversal
                    for transversal in cover_fives
                    if all(transversal & face for face in face_tuple)
                )
                target_witnesses.append((face_tuple, selected))

    print("graph_edges", edges)
    print("base_g3", base_g3)
    print("base_g4", base_g4)
    print("minimal_covers", tuple(vertices(cover) for cover in covers))
    print("independent_active_masks", len(independent_active))
    print("active_configs", len(active_configs))
    print("column_signatures", column_distributions.cache_info().currsize)
    print("replayed_states", states)
    print("bad_replays", bad_replays)
    print("histogram", dict(sorted(histogram.items())))
    print("target_present", TARGET_A0_FACETS in histogram)
    print(
        "target_missing_fives",
        tuple(tuple(vertices(face)) for face in target_witnesses[0][0])
        if target_witnesses
        else None,
    )
    if target_witnesses:
        root = Path(__file__).resolve().parents[1]
        support = json.loads(
            (root / "data/enumerations/h9_14_nerve_support.json").read_text()
        )
        support_record = next(
            item for item in support["records"] if item["graph6"] == GRAPH6
        )
        for witness_index, (_, selected) in enumerate(target_witnesses):
            audit = audit_global_facet_selection(
                GRAPH6,
                support_record["ridge_support_core_masks"],
                support_record["dehn_sommerville_schedule"],
                selected,
            )
            print(
                "target_audit",
                witness_index,
                {
                    "support": audit.selected_inside_support_core,
                    "ridges": audit.all_ridge_degrees_two,
                    "dual_components": audit.dual_component_count,
                    "triples": audit.covers_all_independent_triples,
                    "face_vector": audit.face_vector,
                    "schedule": audit.matches_dehn_sommerville_schedule,
                    "strict_g": audit.first_strict_local_g_failure,
                    "macaulay": audit.first_local_macaulay_failure,
                    "euler": audit.first_local_euler_failure,
                    "homology": audit.first_mod2_homology_failure,
                    "passes_model": audit.passes_selection_model,
                    "passes_homology": audit.passes_connected_homology_audit,
                },
            )


if __name__ == "__main__":
    a0_census()
