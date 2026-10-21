#!/usr/bin/env python3
"""Lossless exact scout for all H9/14 record-32 completions.

This is deliberately a scout, not a theorem certificate.  It independently
reconstructs the Dehn--Sommerville schedule table, the strict compact schedule
cut, the full graph symmetry, and the proof-safe state spaces for the positive
and ``a=0`` branches.  In particular, the ``a=0`` encoding admits every
independent subset of the six active vertices as the active part of each of
the six A5 rows.  It does *not* impose a singleton-active-row normalization.

The default run is quick and prints the schedule/state-space audit.  Optional
positive-branch and a=0 distribution counts are exact but more expensive::

    PYTHONPATH=src:. python3 scripts/explore_h9_14_record32_all_completions.py
    PYTHONPATH=src:. python3 scripts/explore_h9_14_record32_all_completions.py --positive
    PYTHONPATH=src:. python3 scripts/explore_h9_14_record32_all_completions.py --a0-count

No polytopality, compact-Coxeter, or finite-volume conclusion is made by this
file.  Exact target families must still be dualized, topology-audited, and
sent through the unrestricted Lanner/actual-vertex gate.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations
import json
from math import comb
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record32 import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_RECORD_INDEX,
)
from coxeter_polytopes.combinatorics.h9_14_record46 import (
    _a4_interval_feasible,
    _base_link_g3_g4,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import independent_masks
from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6, macaulay_raise


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"

VERTEX_COUNT = 14
DIMENSION = 9
ACTIVE_VERTICES = (0, 1, 2, 11, 12, 13)
ISOLATED_VERTICES = tuple(range(3, 11))
ACTIVE_MASK = sum(1 << vertex for vertex in ACTIVE_VERTICES)
EDGE_MASKS = tuple(
    sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES
)
A0_ROW_COUNT = 6

EXPECTED_SCHEDULES = (
    (0, 6, 114, 36),
    (1, 12, 112, 30),
    (2, 18, 110, 24),
    (3, 24, 108, 18),
    (4, 30, 106, 12),
    (5, 36, 104, 6),
)
EXPECTED_MINIMAL_GRAPH_COVERS = (
    (0, 1, 2),
    (0, 2, 12),
    (0, 1, 13),
    (1, 11, 13),
    (0, 12, 13),
    (11, 12, 13),
)


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1)


def _support_record() -> Mapping[str, object]:
    document = json.loads(SUPPORT.read_text(encoding="utf-8"))
    survivors = [
        row for row in document["records"] if row["survives_second_layer"]
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise AssertionError("record-32 support index or graph6 changed")
    return record


def _summed_vertex_link_g4(f_vector: Sequence[int]) -> int:
    """Double-count the sum of g4 over all fourteen vertex links."""

    summed_link_f = {0: VERTEX_COUNT}
    for order in range(1, 5):
        summed_link_f[order] = (order + 1) * int(f_vector[order])

    def h_entry(index: int) -> int:
        return sum(
            (-1) ** (index - order)
            * comb(8 - order, 8 - index)
            * summed_link_f[order]
            for order in range(index + 1)
        )

    return h_entry(4) - h_entry(3)


def schedule_audit() -> dict[str, object]:
    record = _support_record()
    rows = tuple(
        (
            int(row["missing_four_count"]),
            int(row["missing_five_count_offset"]),
            int(row["f_vector"][-1]),
            _summed_vertex_link_g4(tuple(map(int, row["f_vector"]))),
        )
        for row in record["dehn_sommerville_schedule"]
    )
    if rows != EXPECTED_SCHEDULES:
        raise AssertionError("record-32 Dehn--Sommerville schedules changed")
    return {
        "rows": [
            {
                "missing_four_count": a,
                "total_missing_independent_five_sets": b,
                "rank_nine_facet_count": facets,
                "summed_vertex_link_g4": total,
            }
            for a, b, facets, total in rows
        ],
        "strict_compact_link_sum_minimum": VERTEX_COUNT,
        "strict_compact_surviving_a": [
            a for a, _b, _facets, total in rows if total >= VERTEX_COUNT
        ],
        "strict_compact_excluded_a": [
            a for a, _b, _facets, total in rows if total < VERTEX_COUNT
        ],
    }


def minimal_graph_covers() -> tuple[int, ...]:
    answers = []
    for candidate in range(1, 1 << VERTEX_COUNT):
        if not all(candidate & edge for edge in EDGE_MASKS):
            continue
        if any(
            all((candidate ^ (1 << vertex)) & edge for edge in EDGE_MASKS)
            for vertex in vertices(candidate)
        ):
            continue
        answers.append(candidate)
    return tuple(answers)


@lru_cache(maxsize=1)
def active_graph_automorphisms() -> tuple[tuple[int, ...], ...]:
    edge_set = {frozenset(edge) for edge in TARGET_MISSING_EDGES}
    answers = []
    for image_order in permutations(ACTIVE_VERTICES):
        action = dict(zip(ACTIVE_VERTICES, image_order, strict=True))
        if {
            frozenset((action[left], action[right]))
            for left, right in TARGET_MISSING_EDGES
        } == edge_set:
            answers.append(image_order)
    return tuple(answers)


def independent_active_masks() -> tuple[int, ...]:
    return tuple(
        value
        for order in range(len(ACTIVE_VERTICES) + 1)
        for chosen in combinations(ACTIVE_VERTICES, order)
        for value in (mask(chosen),)
        if all(value & edge != edge for edge in EDGE_MASKS)
    )


def a0_active_configurations() -> tuple[
    tuple[tuple[int, ...], tuple[int, ...]], ...
]:
    """Return the exact proof-safe active-row quotient for ``a=0``.

    Rows are first quotiented only by their genuine ``S6`` permutation action.
    The graph-cover filter is exact for whether an *active* fourth vertex can
    extend a minimal order-three graph cover to a forbidden transversal.
    Isolate extensions are handled later by the allowed-column condition.
    """

    covers = minimal_graph_covers()
    expected_covers = tuple(mask(values) for values in EXPECTED_MINIMAL_GRAPH_COVERS)
    if covers != expected_covers:
        raise AssertionError("record-32 minimal graph covers changed")
    _base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)
    caps = tuple(value - 1 for value in base_g4)
    active_rows = independent_active_masks()
    answers = []
    for configuration in combinations_with_replacement(active_rows, A0_ROW_COUNT):
        if any(
            sum(bool(row & (1 << vertex)) for row in configuration) > caps[vertex]
            for vertex in ACTIVE_VERTICES
        ):
            continue
        missed_rows = tuple(
            sum(
                1 << row_index
                for row_index, row in enumerate(configuration)
                if not row & cover
            )
            for cover in covers
        )
        if any(
            all(
                not (pattern & (1 << row_index))
                or bool(configuration[row_index] & (1 << vertex))
                for row_index in range(A0_ROW_COUNT)
            )
            for pattern in missed_rows
            for vertex in ACTIVE_VERTICES
        ):
            continue
        answers.append((configuration, missed_rows))
    return tuple(answers)


def state_space_audit() -> dict[str, object]:
    order, edges = decode_graph6(TARGET_GRAPH6)
    if order != VERTEX_COUNT or edges != TARGET_MISSING_EDGES:
        raise AssertionError("record-32 graph decode changed")
    active_rows = independent_active_masks()
    configurations = a0_active_configurations()
    _base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)
    return {
        "graph6": TARGET_GRAPH6,
        "decoded_missing_edges": [list(edge) for edge in edges],
        "active_graph_description": "P4 disjoint union K2",
        "active_vertices": list(ACTIVE_VERTICES),
        "isolated_vertices": list(ISOLATED_VERTICES),
        "active_graph_automorphism_count": len(active_graph_automorphisms()),
        "full_symmetry_description": "Aut(P4 disjoint union K2) x S8",
        "minimal_graph_covers": [list(vertices(value)) for value in minimal_graph_covers()],
        "base_vertex_link_g4": list(base_g4),
        "a0_A5_incidence_caps": [value - 1 for value in base_g4],
        "a0_row_count": A0_ROW_COUNT,
        "independent_active_masks": [list(vertices(value)) for value in active_rows],
        "independent_active_mask_size_histogram": dict(
            sorted(Counter(value.bit_count() for value in active_rows).items())
        ),
        "independent_active_mask_count": len(active_rows),
        "active_row_multisets_after_caps_before_cover_filter": (
            _active_configuration_count_before_cover_filter(active_rows, base_g4)
        ),
        "active_row_multisets_after_exact_active_cover_filter": len(configurations),
        "lossless_statement": (
            "Every independent active subset admitted by the strict singleton-link "
            "caps is available in every one of the six unordered A5 rows; the "
            "eight isolate columns are subsequently enumerated as exact six-bit "
            "multiplicity vectors. No singleton-active-row normalization is used."
        ),
    }


def _active_configuration_count_before_cover_filter(
    active_rows: Sequence[int], base_g4: Sequence[int]
) -> int:
    caps = tuple(value - 1 for value in base_g4)
    return sum(
        not any(
            sum(bool(row & (1 << vertex)) for row in configuration) > caps[vertex]
            for vertex in ACTIVE_VERTICES
        )
        for configuration in combinations_with_replacement(active_rows, A0_ROW_COUNT)
    )


def _family_interval_feasible(family: Sequence[int], local: object) -> bool:
    for figure in local.vertex_figures:
        incidence = sum(bool(face & (1 << figure.vertex)) for face in family)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.edge_figures:
        subset = mask(figure.edge)
        incidence = sum(face & subset == subset for face in family)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.triple_figures:
        subset = mask(figure.triple)
        incidence = sum(face & subset == subset for face in family)
        if incidence > figure.missing_four_incidence_max:
            return False
    return True


def _shadow_gate_from_word(
    family: Sequence[int],
    shadow_word: int,
    five_vertex_words: Sequence[int],
    base_g3: Sequence[int],
    base_g4: Sequence[int],
    total_fives: int,
) -> tuple[bool, bool, int, int]:
    shadow_count = shadow_word.bit_count()
    additional = total_fives - shadow_count
    if additional < 0:
        return False, False, shadow_count, additional
    lower = []
    upper = []
    strict = True
    for vertex in range(VERTEX_COUNT):
        a4_incidence = sum(bool(face & (1 << vertex)) for face in family)
        previous = base_g3[vertex] - a4_incidence
        forced_incidence = (shadow_word & five_vertex_words[vertex]).bit_count()
        before = base_g4[vertex] + 6 * a4_incidence - forced_incidence
        strict &= previous >= 1 and before >= 1
        lower.append(
            max(0, before - macaulay_raise(previous, 3))
            if previous >= 0
            else 0
        )
        upper.append(before - 1)
    aggregate = (
        strict
        and all(value >= 0 for value in upper)
        and sum(lower) <= 5 * additional <= sum(upper)
    )
    return strict, aggregate, shadow_count, additional


def _family_orbit_key(family: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    """Canonical A4-family key under the full graph action.

    Isolate permutations are handled by sorting their incidence columns; row
    permutations are exhausted literally.  This is used only on the small set
    surviving the exact shadow gate.
    """

    keys = []
    for image_order in active_graph_automorphisms():
        action = dict(zip(ACTIVE_VERTICES, image_order, strict=True))
        transformed = tuple(
            mask(action.get(vertex, vertex) for vertex in vertices(face))
            for face in family
        )
        for row_order in permutations(range(len(family))):
            active_parts = tuple(
                transformed[index] & ACTIVE_MASK for index in row_order
            )
            isolate_columns = tuple(
                sorted(
                    sum(
                        bool(transformed[row_order[row]] & (1 << isolate)) << row
                        for row in range(len(family))
                    )
                    for isolate in ISOLATED_VERTICES
                )
            )
            keys.append((active_parts, isolate_columns))
    return min(keys)


def positive_branch_scout(max_a: int = 3) -> dict[str, object]:
    """Run the exact local/forced-shadow gate through ``a=max_a``.

    This does not enumerate the additional minimal A5s, so its output is a
    rigorous branch reduction rather than an all-completion result.
    """

    local = analyze_local_link_graph(TARGET_GRAPH6)
    local_figures = (*local.vertex_figures, *local.edge_figures)
    if any(
        figure.missing_four_incidence_interval is None
        or figure.missing_four_incidence_interval[0] != 0
        for figure in local_figures
    ):
        raise AssertionError(
            "record-32 singleton A4 ground set is not proven lossless"
        )
    fours = tuple(
        face
        for face in independent_masks(TARGET_GRAPH6, 4)
        if _a4_interval_feasible(face, local)
    )
    fives = independent_masks(TARGET_GRAPH6, 5)
    shadow_words = tuple(
        sum((five & four == four) << index for index, five in enumerate(fives))
        for four in fours
    )
    five_vertex_words = tuple(
        sum(bool(five & (1 << vertex)) << index for index, five in enumerate(fives))
        for vertex in range(VERTEX_COUNT)
    )
    base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)
    schedules = {a: b for a, b, _facets, _total in EXPECTED_SCHEDULES}
    branches = []
    pair_neighbors = [0] * len(fours)
    for a in range(1, max_a + 1):
        raw_count = comb(len(fours), a)
        interval_count = 0
        gate_histogram: Counter[tuple[bool, bool, int, int]] = Counter()
        viable: list[tuple[int, ...]] = []
        if a == 3:
            cap_two_subsets = []
            for figure in local.vertex_figures:
                interval = figure.missing_four_incidence_interval
                if interval is not None and interval[1] == 2:
                    cap_two_subsets.append(1 << figure.vertex)
            for figure in local.edge_figures:
                interval = figure.missing_four_incidence_interval
                if interval is not None and interval[1] == 2:
                    cap_two_subsets.append(mask(figure.edge))
            for figure in local.triple_figures:
                if figure.missing_four_incidence_max == 2:
                    cap_two_subsets.append(mask(figure.triple))
            cap_two_signatures = tuple(
                sum(
                    (four & subset == subset) << constraint
                    for constraint, subset in enumerate(cap_two_subsets)
                )
                for four in fours
            )

            def triple_indices():
                full = (1 << len(fours)) - 1
                for first in range(len(fours)):
                    seconds = pair_neighbors[first] & (full ^ ((1 << (first + 1)) - 1))
                    while seconds:
                        second_bit = seconds & -seconds
                        second = second_bit.bit_length() - 1
                        seconds ^= second_bit
                        thirds = (
                            pair_neighbors[first]
                            & pair_neighbors[second]
                            & (full ^ ((1 << (second + 1)) - 1))
                        )
                        common_cap_two = (
                            cap_two_signatures[first] & cap_two_signatures[second]
                        )
                        while thirds:
                            third_bit = thirds & -thirds
                            third = third_bit.bit_length() - 1
                            thirds ^= third_bit
                            if common_cap_two & cap_two_signatures[third]:
                                continue
                            yield first, second, third

            index_families = triple_indices()
        else:
            index_families = combinations(range(len(fours)), a)
        for indices in index_families:
            family = tuple(fours[index] for index in indices)
            # For a=3 the generator is already the exact interval test:
            # individual cap-zero constraints were applied when ``fours`` was
            # built, pair cap-one constraints are encoded by pair_neighbors,
            # and cap-two constraints are rejected by cap_two_signatures.
            # Re-running the generic test on every admitted triple makes the
            # finite census needlessly quadratic in its constraint list.
            if a > 1 and a != 3 and not _family_interval_feasible(family, local):
                continue
            interval_count += 1
            if a == 2:
                first, second = indices
                pair_neighbors[first] |= 1 << second
                pair_neighbors[second] |= 1 << first
            shadow_word = 0
            for index in indices:
                shadow_word |= shadow_words[index]
            gate = _shadow_gate_from_word(
                family,
                shadow_word,
                five_vertex_words,
                base_g3,
                base_g4,
                schedules[a],
            )
            gate_histogram[gate] += 1
            if gate[0] and gate[1]:
                viable.append(family)
        orbits: defaultdict[object, list[tuple[int, ...]]] = defaultdict(list)
        for family in viable:
            orbits[_family_orbit_key(family)].append(family)
        representatives = tuple(
            min(orbit) for _key, orbit in sorted(orbits.items())
        )
        branches.append(
            {
                "missing_four_count": a,
                "raw_unordered_A4_family_count": raw_count,
                "family_interval_feasible_count": interval_count,
                "shadow_gate_histogram": [
                    {
                        "strict": gate[0],
                        "aggregate": gate[1],
                        "forced_shadow_count": gate[2],
                        "additional_minimal_A5_count": gate[3],
                        "multiplicity": count,
                    }
                    for gate, count in sorted(gate_histogram.items())
                ],
                "strict_aggregate_viable_family_count": len(viable),
                "full_graph_action_orbit_count": len(orbits),
                "orbit_size_histogram": dict(
                    sorted(Counter(map(len, orbits.values())).items())
                ),
                "orbit_representatives": [
                    [list(vertices(face)) for face in family]
                    for family in representatives
                ],
            }
        )
    return {
        "independent_A4_count": len(independent_masks(TARGET_GRAPH6, 4)),
        "individually_interval_feasible_A4_count": len(fours),
        "singleton_ground_set_lossless_for_families": True,
        "branches": branches,
        "scope": (
            "Exact through local intervals and forced-shadow strict/Macaulay "
            "caps. Additional minimal A5 families, blockers, topology, and "
            "Coxeter labels are not enumerated here."
        ),
    }


@lru_cache(maxsize=None)
def _distribution_count(
    allowed: tuple[int, ...], row_sums: tuple[int, ...]
) -> int:
    @lru_cache(maxsize=None)
    def visit(
        position: int, columns_left: int, sums: tuple[int, ...]
    ) -> int:
        if position == len(allowed):
            return int(columns_left == 0 and not any(sums))
        column = allowed[position]
        maximum = columns_left
        for row, value in enumerate(sums):
            if column & (1 << row):
                maximum = min(maximum, value)
        answer = 0
        for multiplicity in range(maximum + 1):
            next_left = columns_left - multiplicity
            next_sums = tuple(
                value - multiplicity * bool(column & (1 << row))
                for row, value in enumerate(sums)
            )
            if all(0 <= value <= next_left for value in next_sums):
                answer += visit(position + 1, next_left, next_sums)
        return answer

    return visit(0, len(ISOLATED_VERTICES), row_sums)


def a0_distribution_count() -> dict[str, object]:
    """Count the lossless ``a=0`` isolate-column quotient exactly.

    The count includes repeated six-row A5 families; those must be removed in
    the subsequent streaming enumeration.  It is nevertheless an exact size
    for the proof-safe column-multiset search before that postcheck.
    """

    _distribution_count.cache_clear()
    configurations = a0_active_configurations()
    path_count = 0
    positive_problem_count = 0
    allowed_histogram: Counter[int] = Counter()
    path_count_histogram: Counter[int] = Counter()
    for configuration, missed_rows in configurations:
        allowed = tuple(
            column
            for column in range(1 << A0_ROW_COUNT)
            if column.bit_count() <= A0_ROW_COUNT - 1
            and all(column & rows != rows for rows in missed_rows)
        )
        row_sums = tuple(5 - row.bit_count() for row in configuration)
        count = _distribution_count(allowed, row_sums)
        allowed_histogram[len(allowed)] += 1
        path_count_histogram[count] += 1
        path_count += count
        positive_problem_count += count > 0
    return {
        "active_row_configuration_count": len(configurations),
        "distinct_column_distribution_problem_count": (
            _distribution_count.cache_info().currsize
        ),
        "positive_column_distribution_problem_count": positive_problem_count,
        "allowed_column_count_histogram": dict(sorted(allowed_histogram.items())),
        "distribution_count_histogram": dict(sorted(path_count_histogram.items())),
        "total_column_multiset_path_count_before_repeated_row_postcheck": path_count,
        "scope": (
            "Exact lossless S8 column-multiset count before repeated-row, "
            "literal small-transversal, scheduled-facet, topology, and label gates."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--positive", action="store_true")
    parser.add_argument("--positive-max-a", type=int, choices=(1, 2, 3), default=3)
    parser.add_argument("--a0-count", action="store_true")
    args = parser.parse_args()
    result: dict[str, object] = {
        "kind": "h9_14_record32_lossless_all_completion_scout",
        "theorem_certificate": False,
        "schedule_audit": schedule_audit(),
        "state_space_audit": state_space_audit(),
    }
    if args.positive:
        result["positive_branch_scout"] = positive_branch_scout(args.positive_max_a)
    if args.a0_count:
        result["a0_distribution_count"] = a0_distribution_count()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
