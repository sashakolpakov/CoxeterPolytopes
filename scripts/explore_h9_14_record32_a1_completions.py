#!/usr/bin/env python3
"""Losslessly close the four strict-shadow record-32 ``a=1`` orbits.

Raw A5 combinations are as large as ``C(1252,4)``.  The exact enumeration
instead uses only genuine symmetries of a fixed A4 representative.  Its eight
graph isolates split into those inside and outside the A4; each class is
encoded by the complete multiset of its A5-row incidence columns.  Every
independent active subset is allowed in every row.  Thus the encoding is a
lossless quotient by the stabilizer's two isolate symmetric groups, not a
singleton-active-row ansatz or a sampling.

Each canonical completion is checked literally for distinct minimal rows,
vertex caps, transversals of order at most four, and the scheduled 112 facets.
Scheduled targets are dualized; pure ridge pseudomanifolds are immediately
sent through the sparse Lanner/actual-vertex gate.  This remains a scout until
the census and label transcripts have an independent replay.  It is strictly
compact-only and makes no finite-volume inference.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations_with_replacement, permutations
import argparse
import json
from math import comb
from typing import Iterable, Sequence

from coxeter_polytopes.combinatorics.h9_14_record32 import TARGET_GRAPH6
from coxeter_polytopes.combinatorics.h9_14_record35_all_completions import (
    _forced_five_shadow,
    _shadow_gate,
)
from coxeter_polytopes.combinatorics.h9_14_record46 import _base_link_g3_g4
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
)
from scripts.explore_h9_14_record32_a0_target_gate import (
    EDGE_MASKS,
    minimal_transversals,
    ridge_degree_histogram,
    sparse_label_gate,
)
from scripts.explore_h9_14_record32_a2_additions import (
    hit_word,
    mask,
    small_base_transversals,
    vertices,
)
from scripts.explore_h9_14_record32_all_completions import (
    ACTIVE_MASK,
    ACTIVE_VERTICES,
    ISOLATED_VERTICES,
    active_graph_automorphisms,
    independent_active_masks,
)


VERTEX_COUNT = 14
TOTAL_FIVES = 12
SCHEDULED_FACETS = 112
A1_REPRESENTATIVES = (
    (3, 4, 5, 6),
    (2, 3, 4, 5),
    (1, 2, 3, 4),
    (2, 3, 4, 11),
)
EXPECTED_CANDIDATE_COUNTS = (664, 672, 672, 1_252)


def _rows_from_columns(
    active_rows: Sequence[int],
    groups: Sequence[tuple[Sequence[int], Sequence[int]]],
) -> tuple[int, ...]:
    rows = list(active_rows)
    for group_vertices, columns in groups:
        if len(group_vertices) != len(columns):
            raise AssertionError("isolate column count changed")
        for vertex, column in zip(group_vertices, columns, strict=True):
            for row in range(len(rows)):
                if column & (1 << row):
                    rows[row] |= 1 << vertex
    return tuple(rows)


@lru_cache(maxsize=None)
def column_multisets_by_row_sums(
    row_count: int, column_count: int, column_cap: int
) -> dict[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    """All exact ``S_column_count`` column multisets, grouped by row sums."""

    allowed = tuple(
        column
        for column in range(1 << row_count)
        if column.bit_count() <= column_cap
    )
    grouped: defaultdict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
    for columns in combinations_with_replacement(allowed, column_count):
        sums = tuple(
            sum(bool(column & (1 << row)) for column in columns)
            for row in range(row_count)
        )
        grouped[sums].append(columns)
    return {key: tuple(value) for key, value in grouped.items()}


def _active_stabilizer(four: int) -> tuple[dict[int, int], ...]:
    answers = []
    for image_order in active_graph_automorphisms():
        action = dict(zip(ACTIVE_VERTICES, image_order, strict=True))
        image = mask(action.get(vertex, vertex) for vertex in vertices(four))
        if image == four:
            answers.append(action)
    return tuple(answers)


def family_orbit_key(four: int, additions: Sequence[int]) -> object:
    """Canonical key under the full graph stabilizer of the fixed A4."""

    inside = tuple(vertex for vertex in ISOLATED_VERTICES if four & (1 << vertex))
    outside = tuple(vertex for vertex in ISOLATED_VERTICES if not four & (1 << vertex))
    keys = []
    for action in _active_stabilizer(four):
        transformed = tuple(
            mask(action.get(vertex, vertex) for vertex in vertices(face))
            for face in additions
        )
        # Sorting rows removes their genuine permutation action.  Isolate
        # permutations are represented exactly by sorted incidence columns
        # within the two fixed-A4 classes.
        for ordered in permutations(transformed):
            active_parts = tuple(face & ACTIVE_MASK for face in ordered)
            inside_columns = tuple(
                sorted(
                    sum(bool(ordered[row] & (1 << vertex)) << row for row in range(len(ordered)))
                    for vertex in inside
                )
            )
            outside_columns = tuple(
                sorted(
                    sum(bool(ordered[row] & (1 << vertex)) << row for row in range(len(ordered)))
                    for vertex in outside
                )
            )
            keys.append((active_parts, inside_columns, outside_columns))
    return min(keys)


def audit_representative(
    index: int,
    representative: Iterable[int],
    *,
    state_cap: int = 250_000,
    incidence_only: bool = False,
) -> dict[str, object]:
    four = mask(representative)
    family = (four,)
    fives = independent_masks(TARGET_GRAPH6, 5)
    base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)
    shadow = _forced_five_shadow(family, fives)
    strict, aggregate, additional, lower, upper = _shadow_gate(
        family, shadow, base_g3, base_g4, TOTAL_FIVES
    )
    if not strict or not aggregate or any(lower):
        raise AssertionError("record-32 a=1 representative gate changed")
    candidate_set = {
        five
        for five in fives
        if five not in shadow
        and five & four != four
        and all(
            not (five & (1 << vertex)) or upper[vertex] >= 1
            for vertex in range(VERTEX_COUNT)
        )
    }
    if len(candidate_set) != EXPECTED_CANDIDATE_COUNTS[index]:
        raise AssertionError("record-32 a=1 candidate count changed")

    inside = tuple(vertex for vertex in ISOLATED_VERTICES if four & (1 << vertex))
    outside = tuple(vertex for vertex in ISOLATED_VERTICES if not four & (1 << vertex))
    inside_cap = min((upper[vertex] for vertex in inside), default=additional)
    outside_cap = min((upper[vertex] for vertex in outside), default=additional)
    if any(upper[vertex] != inside_cap for vertex in inside) or any(
        upper[vertex] != outside_cap for vertex in outside
    ):
        raise AssertionError("fixed-A4 isolate classes do not have uniform caps")
    inside_patterns = column_multisets_by_row_sums(
        additional, len(inside), inside_cap
    )
    outside_patterns = column_multisets_by_row_sums(
        additional, len(outside), outside_cap
    )

    allowed_active_rows = tuple(
        row
        for row in independent_active_masks()
        if all(not (row & (1 << vertex)) or upper[vertex] >= 1 for vertex in ACTIVE_VERTICES)
    )
    active_configuration_count = 0
    column_pair_path_count = 0
    distinct_candidate_row_path_count = 0
    no_small_path_count = 0
    facet_histogram: Counter[int] = Counter()
    target_families: set[tuple[int, ...]] = set()

    small = small_base_transversals((*EDGE_MASKS, four, *shadow))
    full_small_word = (1 << len(small)) - 1
    small_words = {five: hit_word(five, small) for five in candidate_set}
    five_covers = graph_vertex_cover_masks(TARGET_GRAPH6, 5)
    fixed_five_word = hit_word(four, five_covers)
    five_words = {five: hit_word(five, five_covers) for five in candidate_set}

    for active_rows in combinations_with_replacement(
        allowed_active_rows, additional
    ):
        if any(
            sum(bool(row & (1 << vertex)) for row in active_rows) > upper[vertex]
            for vertex in ACTIVE_VERTICES
        ):
            continue
        active_configuration_count += 1
        needed = tuple(5 - row.bit_count() for row in active_rows)
        for inside_sums, inside_columns_rows in inside_patterns.items():
            outside_sums = tuple(
                needed[row] - inside_sums[row] for row in range(additional)
            )
            if min(outside_sums, default=0) < 0:
                continue
            outside_columns_rows = outside_patterns.get(outside_sums, ())
            for inside_columns in inside_columns_rows:
                for outside_columns in outside_columns_rows:
                    column_pair_path_count += 1
                    rows = _rows_from_columns(
                        active_rows,
                        ((inside, inside_columns), (outside, outside_columns)),
                    )
                    if len(set(rows)) != additional or any(
                        row not in candidate_set for row in rows
                    ):
                        continue
                    distinct_candidate_row_path_count += 1
                    small_word = full_small_word
                    for row in rows:
                        small_word &= small_words[row]
                    if small_word:
                        continue
                    no_small_path_count += 1
                    word = fixed_five_word
                    for row in rows:
                        word &= five_words[row]
                    facets = word.bit_count()
                    facet_histogram[facets] += 1
                    if facets == SCHEDULED_FACETS:
                        target_families.add(tuple(sorted(rows)))

    target_orbits: defaultdict[object, list[tuple[int, ...]]] = defaultdict(list)
    for additions in sorted(target_families):
        target_orbits[family_orbit_key(four, additions)].append(additions)
    representatives = tuple(min(rows) for rows in target_orbits.values())
    target_audits = []
    for additions in representatives:
        minimal_faces = (four, *additions)
        transversals = minimal_transversals((*EDGE_MASKS, *minimal_faces))
        histogram = tuple(
            sorted(Counter(value.bit_count() for value in transversals).items())
        )
        five_transversals = tuple(
            value for value in transversals if value.bit_count() == 5
        )
        ridge_histogram = ridge_degree_histogram(five_transversals)
        pure_ridge = (
            histogram == ((5, SCHEDULED_FACETS),)
            and bool(ridge_histogram)
            and set(dict(ridge_histogram)) == {2}
        )
        target_audits.append(
            {
                "additional_A5_masks": list(additions),
                "additional_A5s": [list(vertices(face)) for face in additions],
                "minimal_transversal_size_histogram": dict(histogram),
                "ridge_degree_histogram": dict(ridge_histogram),
                "pure_ridge_pseudomanifold": pure_ridge,
                "sparse_label_gate": (
                    sparse_label_gate(minimal_faces, state_cap)
                    if pure_ridge and not incidence_only
                    else None
                ),
            }
        )

    return {
        "orbit_representative_index": index,
        "missing_four": list(representative),
        "forced_five_shadow_count": len(shadow),
        "additional_minimal_A5_count": additional,
        "individual_candidate_count": len(candidate_set),
        "raw_unordered_addition_count": comb(len(candidate_set), additional),
        "inside_A4_isolate_count": len(inside),
        "outside_A4_isolate_count": len(outside),
        "active_row_configuration_count": active_configuration_count,
        "isolate_column_multiset_pair_path_count": column_pair_path_count,
        "distinct_candidate_row_path_count": distinct_candidate_row_path_count,
        "no_transversal_of_order_at_most_four_path_count": no_small_path_count,
        "facet_count_histogram": dict(sorted(facet_histogram.items())),
        "literal_canonical_target_family_count": len(target_families),
        "fixed_A4_stabilizer_target_orbit_count": len(target_orbits),
        "target_audits": target_audits,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit", type=int, choices=range(4), action="append")
    parser.add_argument("--incidence-only", action="store_true")
    parser.add_argument("--state-cap", type=int, default=250_000)
    args = parser.parse_args()
    indices = tuple(args.orbit) if args.orbit else tuple(range(4))
    audits = tuple(
        audit_representative(
            index,
            A1_REPRESENTATIVES[index],
            state_cap=args.state_cap,
            incidence_only=args.incidence_only,
        )
        for index in indices
    )
    print(
        json.dumps(
            {
                "kind": "h9_14_record32_a1_lossless_completion_scout",
                "theorem_certificate": False,
                "audited_orbit_indices": list(indices),
                "representative_audits": audits,
                "scope": (
                    "Lossless fixed-A4 isolate-column stabilizer quotient with "
                    "all independent active row parts. Census and label stages "
                    "need independent replay. Compact-only; no finite-volume claim."
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
