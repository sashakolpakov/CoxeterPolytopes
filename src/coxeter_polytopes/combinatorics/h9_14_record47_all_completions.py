"""Exact compact-candidate completion census for the H9/14 record-47 shadow.

The record-47 missing-edge graph has six active vertices and eight isolates.
The four Dehn--Sommerville schedules are first reduced by strict positivity of
the fourteen vertex-link ``g_4`` values.  The remaining ``a=1`` schedule is
excluded by enumerating every possible missing four-face and its forced
five-shadow.  In the ``a=0`` schedule, the four unordered missing five-faces
are enumerated up to permutation of the eight graph isolates by their four-bit
column multiplicities.

This module is exhaustive under the compact-Coxeter missing-face necessities;
it does not classify arbitrary ordinary completions whose minimal nonfaces
could have order greater than five.  The compact Coxeter conclusion is
obtained by chaining its unique survivor to the separately frozen,
unrestricted selected-nerve label obstruction.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations
import json
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record47 import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_RECORD_INDEX,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.local_link import (
    analyze_local_link_graph,
)
from coxeter_polytopes.combinatorics.missing_edges import (
    decode_graph6,
    macaulay_raise,
)


VERTEX_COUNT = 14
FULL_MASK = (1 << VERTEX_COUNT) - 1
ACTIVE_VERTICES = (0, 1, 2, 11, 12, 13)
ISOLATED_VERTICES = tuple(range(3, 11))
EDGE_MASKS = tuple(
    sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES
)

# (missing-four count, independent five-nonface count, facet count, summed link g4)
EXPECTED_SCHEDULES = (
    (0, 4, 90, 24),
    (1, 10, 88, 18),
    (2, 16, 86, 12),
    (3, 22, 84, 6),
)
EXPECTED_A1_GATE_HISTOGRAM = (
    ((False, False, 8, 2), 84),
    ((False, False, 9, 1), 168),
    ((False, False, 10, 0), 70),
)
EXPECTED_A0_PATH_FACET_COUNT_HISTOGRAM = (
    (58, 3),
    (60, 30),
    (61, 45),
    (62, 27),
    (63, 12),
    (64, 195),
    (65, 30),
    (66, 27),
    (67, 87),
    (68, 144),
    (69, 1),
    (70, 48),
    (71, 18),
    (72, 31),
    (73, 3),
    (74, 6),
    (75, 6),
    (76, 12),
    (79, 6),
    (82, 3),
    (83, 9),
    (90, 1),
)
EXPECTED_A0_PATH_TRANSCRIPT_SHA256 = (
    "33d1820fec2a111c75049ad08757d48ef0e222ca19ec57aa6f99c167a258b1b2"
)
EXPECTED_A0_EXACT_MASK_FACET_COUNT_HISTOGRAM = (
    (58, 3),
    (60, 30),
    (61, 45),
    (62, 27),
    (63, 9),
    (64, 161),
    (65, 30),
    (66, 27),
    (67, 57),
    (68, 106),
    (69, 1),
    (70, 39),
    (71, 18),
    (72, 24),
    (73, 3),
    (74, 6),
    (75, 6),
    (76, 9),
    (79, 6),
    (82, 3),
    (83, 7),
    (90, 1),
)
EXPECTED_A0_EXACT_MASK_TRANSCRIPT_SHA256 = (
    "ab285a0f421ec7e86b75088b869b74e839baf0d787bd612ea38a274aea746b8f"
)
EXPECTED_A0_QUOTIENT_FACET_COUNT_HISTOGRAM = (
    (58, 3),
    (60, 9),
    (61, 15),
    (62, 12),
    (63, 6),
    (64, 49),
    (65, 16),
    (66, 11),
    (67, 29),
    (68, 26),
    (69, 1),
    (70, 21),
    (71, 9),
    (72, 9),
    (73, 3),
    (74, 3),
    (75, 3),
    (76, 9),
    (79, 3),
    (82, 3),
    (83, 6),
    (90, 1),
)
EXPECTED_A0_QUOTIENT_TRANSCRIPT_SHA256 = (
    "6f175b31b17bc8411ed841fd6f88537dc9676b3886b67366eda8755b64ef7d6d"
)
EXPECTED_A0_PATH_MULTIPLICITY_HISTOGRAM = (
    (1, 47),
    (2, 123),
    (3, 29),
    (4, 1),
    (6, 38),
    (12, 7),
    (24, 2),
)
EXPECTED_MINIMAL_GRAPH_COVERS = (4101, 8195, 12289, 14336)


def mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if value & (1 << vertex))


def _is_independent(value: int) -> bool:
    return all(value & edge != edge for edge in EDGE_MASKS)


@lru_cache(maxsize=None)
def _independent_masks(order: int) -> tuple[int, ...]:
    return tuple(
        mask(chosen)
        for chosen in combinations(range(VERTEX_COUNT), order)
        if _is_independent(mask(chosen))
    )


def _minimal_graph_covers() -> tuple[int, ...]:
    result = []
    for value in range(1, 1 << VERTEX_COUNT):
        if not all(value & edge for edge in EDGE_MASKS):
            continue
        if all(
            not all((value ^ (1 << vertex)) & edge for edge in EDGE_MASKS)
            for vertex in vertices(value)
        ):
            result.append(value)
    return tuple(result)


def _graph_cover_fives() -> tuple[int, ...]:
    return tuple(
        value
        for chosen in combinations(range(VERTEX_COUNT), 5)
        for value in (mask(chosen),)
        if all(value & edge for edge in EDGE_MASKS)
    )


def _base_link_g3_g4() -> tuple[tuple[int, ...], tuple[int, ...]]:
    independent = {order: _independent_masks(order) for order in range(2, 6)}
    g3_values = []
    g4_values = []
    for vertex in range(VERTEX_COUNT):
        counts = {
            order: sum(bool(value & (1 << vertex)) for value in independent[order])
            for order in independent
        }
        f0, f1, f2, f3 = (counts[order] for order in range(2, 6))
        g3_values.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4_values.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    return tuple(g3_values), tuple(g4_values)


def _a4_interval_feasible(value: int, local: object) -> bool:
    for figure in local.vertex_figures:
        incidence = int(bool(value & (1 << figure.vertex)))
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.edge_figures:
        subset = mask(figure.edge)
        incidence = int(value & subset == subset)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.triple_figures:
        subset = mask(figure.triple)
        incidence = int(value & subset == subset)
        if incidence > figure.missing_four_incidence_max:
            return False
    return True


def _a1_shadow_gate(
    four: int,
    independent_fives: Sequence[int],
    base_g3: Sequence[int],
    base_g4: Sequence[int],
) -> tuple[bool, bool, int, int]:
    shadow = tuple(value for value in independent_fives if value & four == four)
    additional = 10 - len(shadow)
    if additional < 0:
        return False, False, len(shadow), additional
    strict = True
    lower_sum = 0
    upper_sum = 0
    for vertex in range(VERTEX_COUNT):
        incidence = int(bool(four & (1 << vertex)))
        previous = base_g3[vertex] - incidence
        shadow_incidence = sum(bool(value & (1 << vertex)) for value in shadow)
        before_additional = base_g4[vertex] + 6 * incidence - shadow_incidence
        if previous < 1 or before_additional < 1:
            strict = False
        if previous >= 0:
            lower_sum += max(
                0,
                before_additional - macaulay_raise(previous, 3),
            )
        upper_sum += before_additional - 1
    aggregate = strict and lower_sum <= 5 * additional <= upper_sum
    return strict, aggregate, len(shadow), additional


@lru_cache(maxsize=None)
def _column_distributions(
    allowed: tuple[int, ...],
    row_sums: tuple[int, ...],
) -> tuple[tuple[int, ...], ...]:
    """Enumerate multisets of the eight isolate-incidence columns."""

    counts = [0] * len(allowed)
    answers: list[tuple[int, ...]] = []

    def recurse(position: int, left: int, sums: tuple[int, ...]) -> None:
        if position == len(allowed):
            if left == 0 and not any(sums):
                answers.append(tuple(counts))
            return
        column = allowed[position]
        maximum = left
        for row in range(4):
            if column & (1 << row):
                maximum = min(maximum, sums[row])
        for multiplicity in range(maximum + 1):
            next_left = left - multiplicity
            next_sums = tuple(
                sums[row] - multiplicity * bool(column & (1 << row))
                for row in range(4)
            )
            if all(0 <= value <= next_left for value in next_sums):
                counts[position] = multiplicity
                recurse(position + 1, next_left, next_sums)
        counts[position] = 0

    recurse(0, len(ISOLATED_VERTICES), row_sums)
    return tuple(answers)


def _has_small_transversal(
    missing_fives: Sequence[int], graph_covers: Sequence[int]
) -> bool:
    """Test whether a minimal graph cover extends inside the size-four budget."""

    for cover in graph_covers:
        missed = tuple(value for value in missing_fives if not value & cover)
        if not missed:
            return True
        common = missed[0]
        for value in missed[1:]:
            common &= value
        if common:
            return True
    return False


def _state_digest(rows: Sequence[tuple[tuple[int, ...], int]]) -> str:
    payload = [
        {"missing_five_masks": list(fives), "facet_count": facet_count}
        for fives, facet_count in rows
    ]
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _quotient_key(
    missing_fives: Sequence[int],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Canonicalize an unordered row family under the graph-isolate ``S8``."""

    active_mask = mask(ACTIVE_VERTICES)
    keys = []
    for row_order in permutations(range(4)):
        active_rows = tuple(
            missing_fives[index] & active_mask for index in row_order
        )
        isolate_columns = tuple(
            sorted(
                sum(
                    bool(missing_fives[row_order[position]] & (1 << isolate))
                    << position
                    for position in range(4)
                )
                for isolate in ISOLATED_VERTICES
            )
        )
        keys.append((active_rows, isolate_columns))
    return min(keys)


def _quotient_digest(
    rows: Sequence[
        tuple[tuple[tuple[int, ...], tuple[int, ...]], int]
    ],
) -> str:
    payload = [
        {
            "active_row_masks": list(key[0]),
            "isolate_columns": list(key[1]),
            "facet_count": facet_count,
        }
        for key, facet_count in rows
    ]
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


@dataclass(frozen=True)
class ScheduleAudit:
    rows: tuple[tuple[int, int, int, int], ...]
    strict_sum_required: int
    surviving_missing_four_counts: tuple[int, ...]
    excluded_missing_four_counts: tuple[int, ...]

    @property
    def verified(self) -> bool:
        return (
            self.rows == EXPECTED_SCHEDULES
            and self.strict_sum_required == VERTEX_COUNT
            and self.surviving_missing_four_counts == (0, 1)
            and self.excluded_missing_four_counts == (2, 3)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "rows": [
                {
                    "missing_four_count": a,
                    "independent_five_nonface_count": b,
                    "facet_count": facets,
                    "summed_vertex_link_g4": total,
                }
                for a, b, facets, total in self.rows
            ],
            "strict_sum_required": self.strict_sum_required,
            "surviving_missing_four_counts": list(
                self.surviving_missing_four_counts
            ),
            "excluded_missing_four_counts": list(
                self.excluded_missing_four_counts
            ),
            "verified": self.verified,
        }


@dataclass(frozen=True)
class A1Audit:
    independent_a4_count: int
    interval_feasible_a4_count: int
    gate_histogram: tuple[tuple[tuple[bool, bool, int, int], int], ...]
    strict_shadow_feasible_count: int
    aggregate_a5_feasible_count: int

    @property
    def verified(self) -> bool:
        return (
            self.independent_a4_count == 690
            and self.interval_feasible_a4_count == 322
            and self.gate_histogram == EXPECTED_A1_GATE_HISTOGRAM
            and self.strict_shadow_feasible_count == 0
            and self.aggregate_a5_feasible_count == 0
        )

    def to_record(self) -> dict[str, object]:
        return {
            "independent_a4_count": self.independent_a4_count,
            "interval_feasible_a4_count": self.interval_feasible_a4_count,
            "gate_histogram": [
                {
                    "strict_vertex_link_g_positive": gate[0],
                    "aggregate_additional_a5_feasible": gate[1],
                    "forced_a5_shadow_size": gate[2],
                    "additional_a5_count": gate[3],
                    "multiplicity": count,
                }
                for gate, count in self.gate_histogram
            ],
            "strict_shadow_feasible_count": self.strict_shadow_feasible_count,
            "aggregate_a5_feasible_count": self.aggregate_a5_feasible_count,
            "branch_excluded": self.aggregate_a5_feasible_count == 0,
            "verified": self.verified,
        }


@dataclass(frozen=True)
class A0Audit:
    active_vertices: tuple[int, ...]
    isolated_vertices: tuple[int, ...]
    minimal_graph_covers: tuple[int, ...]
    independent_active_mask_count: int
    active_row_configuration_count: int
    column_distribution_problem_count: int
    repeated_row_path_count: int
    enumeration_path_count: int
    exact_mask_state_count: int
    canonical_quotient_orbit_count: int
    small_transversal_path_count: int
    path_facet_count_histogram: tuple[tuple[int, int], ...]
    exact_mask_facet_count_histogram: tuple[tuple[int, int], ...]
    quotient_facet_count_histogram: tuple[tuple[int, int], ...]
    path_multiplicity_per_quotient_orbit_histogram: tuple[tuple[int, int], ...]
    path_transcript_sha256: str
    exact_mask_transcript_sha256: str
    quotient_transcript_sha256: str
    target_path_count: int
    target_exact_mask_state_count: int
    target_quotient_orbit_count: int
    unique_target_missing_fives: tuple[int, ...]

    @property
    def verified(self) -> bool:
        return (
            self.active_vertices == ACTIVE_VERTICES
            and self.isolated_vertices == ISOLATED_VERTICES
            and self.minimal_graph_covers == EXPECTED_MINIMAL_GRAPH_COVERS
            and self.independent_active_mask_count == 8
            and self.active_row_configuration_count == 78
            and self.repeated_row_path_count == 18
            and self.enumeration_path_count == 744
            and self.exact_mask_state_count == 618
            and self.canonical_quotient_orbit_count == 247
            and self.small_transversal_path_count == 0
            and self.path_facet_count_histogram
            == EXPECTED_A0_PATH_FACET_COUNT_HISTOGRAM
            and self.exact_mask_facet_count_histogram
            == EXPECTED_A0_EXACT_MASK_FACET_COUNT_HISTOGRAM
            and self.quotient_facet_count_histogram
            == EXPECTED_A0_QUOTIENT_FACET_COUNT_HISTOGRAM
            and self.path_multiplicity_per_quotient_orbit_histogram
            == EXPECTED_A0_PATH_MULTIPLICITY_HISTOGRAM
            and self.path_transcript_sha256 == EXPECTED_A0_PATH_TRANSCRIPT_SHA256
            and self.exact_mask_transcript_sha256
            == EXPECTED_A0_EXACT_MASK_TRANSCRIPT_SHA256
            and self.quotient_transcript_sha256
            == EXPECTED_A0_QUOTIENT_TRANSCRIPT_SHA256
            and self.target_path_count == 1
            and self.target_exact_mask_state_count == 1
            and self.target_quotient_orbit_count == 1
            and tuple(vertices(value) for value in self.unique_target_missing_fives)
            == TARGET_MISSING_FIVES
        )

    def to_record(self) -> dict[str, object]:
        return {
            "enumeration_and_quotient": (
                "78 partially row-sorted active configurations and S8 "
                "isolate-column multiplicities give 744 lossless paths. "
                "Exact-mask deduplication gives 618 states. Canonicalization "
                "over all 24 row permutations and the isolate-column "
                "multiset gives 247 S4-row x S8-isolate orbits. No "
                "active-graph automorphism quotient is used."
            ),
            "active_vertices": list(self.active_vertices),
            "isolated_vertices": list(self.isolated_vertices),
            "minimal_graph_covers": [
                list(vertices(value)) for value in self.minimal_graph_covers
            ],
            "independent_active_mask_count": self.independent_active_mask_count,
            "active_row_configuration_count": self.active_row_configuration_count,
            "column_distribution_problem_count": (
                self.column_distribution_problem_count
            ),
            "column_distribution_problem_count_scope": (
                "diagnostic for this deterministic recursion order, not a "
                "mathematical quotient invariant"
            ),
            "repeated_row_path_count": self.repeated_row_path_count,
            "enumeration_path_count": self.enumeration_path_count,
            "exact_mask_state_count": self.exact_mask_state_count,
            "canonical_quotient_orbit_count": (
                self.canonical_quotient_orbit_count
            ),
            "small_transversal_path_count": (
                self.small_transversal_path_count
            ),
            "path_facet_count_histogram": {
                str(facet_count): count
                for facet_count, count in self.path_facet_count_histogram
            },
            "exact_mask_facet_count_histogram": {
                str(facet_count): count
                for facet_count, count in self.exact_mask_facet_count_histogram
            },
            "quotient_facet_count_histogram": {
                str(facet_count): count
                for facet_count, count in self.quotient_facet_count_histogram
            },
            "path_multiplicity_per_quotient_orbit_histogram": {
                str(multiplicity): count
                for multiplicity, count in (
                    self.path_multiplicity_per_quotient_orbit_histogram
                )
            },
            "path_transcript_sha256": self.path_transcript_sha256,
            "exact_mask_transcript_sha256": self.exact_mask_transcript_sha256,
            "quotient_transcript_sha256": self.quotient_transcript_sha256,
            "target_facet_count": 90,
            "target_path_count": self.target_path_count,
            "target_exact_mask_state_count": self.target_exact_mask_state_count,
            "target_quotient_orbit_count": self.target_quotient_orbit_count,
            "unique_target_missing_fives": [
                list(vertices(value)) for value in self.unique_target_missing_fives
            ],
            "unique_target_equals_selected_record47_nerve": (
                tuple(vertices(value) for value in self.unique_target_missing_fives)
                == TARGET_MISSING_FIVES
            ),
            "verified": self.verified,
        }


@dataclass(frozen=True)
class Record47AllCompletionsAudit:
    decoded_edges: tuple[tuple[int, int], ...]
    schedule: ScheduleAudit
    a1: A1Audit
    a0: A0Audit

    @property
    def verified(self) -> bool:
        return (
            self.decoded_edges == TARGET_MISSING_EDGES
            and self.schedule.verified
            and self.a1.verified
            and self.a0.verified
        )

    @property
    def exact_residual_orbit_count(self) -> int:
        return int(self.verified and self.a0.target_quotient_orbit_count == 1)

    def to_record(self) -> dict[str, object]:
        return {
            "record_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
            "decoded_edges": [list(edge) for edge in self.decoded_edges],
            "schedule_audit": self.schedule.to_record(),
            "a1_audit": self.a1.to_record(),
            "a0_audit": self.a0.to_record(),
            "exact_residual_incidence_orbit_count_under_compact_necessities": (
                self.exact_residual_orbit_count
            ),
            "whole_shadow_incidence_census_complete": self.verified,
            "verified": self.verified,
        }


def _audit_schedules(
    schedules: Sequence[Mapping[str, object]],
) -> ScheduleAudit:
    rows = tuple(
        (
            int(row["missing_four_count"]),
            int(row["missing_five_count_offset"]),
            int(row["f_vector"][-1]),  # type: ignore[index]
            summed_vertex_link_g4(row["f_vector"]),  # type: ignore[arg-type]
        )
        for row in schedules
    )
    surviving = tuple(a for a, _b, _facets, total in rows if total >= VERTEX_COUNT)
    excluded = tuple(a for a, _b, _facets, total in rows if total < VERTEX_COUNT)
    audit = ScheduleAudit(
        rows=rows,
        strict_sum_required=VERTEX_COUNT,
        surviving_missing_four_counts=surviving,
        excluded_missing_four_counts=excluded,
    )
    if not audit.verified:
        raise AssertionError("record-47 Dehn--Sommerville schedule audit changed")
    return audit


def _audit_a1() -> A1Audit:
    independent_fours = _independent_masks(4)
    independent_fives = _independent_masks(5)
    local = analyze_local_link_graph(TARGET_GRAPH6)
    interval_fours = tuple(
        value for value in independent_fours if _a4_interval_feasible(value, local)
    )
    base_g3, base_g4 = _base_link_g3_g4()
    gates = tuple(
        _a1_shadow_gate(value, independent_fives, base_g3, base_g4)
        for value in interval_fours
    )
    audit = A1Audit(
        independent_a4_count=len(independent_fours),
        interval_feasible_a4_count=len(interval_fours),
        gate_histogram=tuple(sorted(Counter(gates).items())),
        strict_shadow_feasible_count=sum(gate[0] for gate in gates),
        aggregate_a5_feasible_count=sum(gate[1] for gate in gates),
    )
    if not audit.verified:
        raise AssertionError("record-47 a=1 exclusion changed")
    return audit


def _audit_a0() -> A0Audit:
    _column_distributions.cache_clear()
    graph_covers = _minimal_graph_covers()
    graph_cover_fives = _graph_cover_fives()
    active_edge_masks = tuple(
        edge for edge in EDGE_MASKS if edge & mask(ACTIVE_VERTICES) == edge
    )
    _base_g3, base_g4 = _base_link_g3_g4()
    incidence_caps = tuple(value - 1 for value in base_g4)
    independent_active = tuple(
        value
        for order in range(6)
        for chosen in combinations(ACTIVE_VERTICES, order)
        for value in (mask(chosen),)
        if all(value & edge != edge for edge in active_edge_masks)
        and all(
            int(bool(value & (1 << vertex))) <= incidence_caps[vertex]
            for vertex in ACTIVE_VERTICES
        )
    )

    active_configurations = []
    for configuration in combinations_with_replacement(independent_active, 4):
        if any(
            sum(bool(row & (1 << vertex)) for row in configuration)
            > incidence_caps[vertex]
            for vertex in ACTIVE_VERTICES
        ):
            continue
        missed_rows = tuple(
            sum(
                1 << row_index
                for row_index, row in enumerate(configuration)
                if not row & cover
            )
            for cover in graph_covers
        )
        if any(
            all(
                not (pattern & (1 << row_index))
                or configuration[row_index] & (1 << vertex)
                for row_index in range(4)
            )
            for pattern in missed_rows
            for vertex in ACTIVE_VERTICES
        ):
            continue
        active_configurations.append((configuration, missed_rows))

    states: list[tuple[tuple[int, ...], int]] = []
    repeated = 0
    small_transversal_failures = 0
    target_states: list[tuple[int, ...]] = []
    for configuration, missed_rows in active_configurations:
        allowed_columns = tuple(
            column
            for column in range(1 << 4)
            if column.bit_count() <= 3
            and all(column & rows != rows for rows in missed_rows)
        )
        row_sums = tuple(5 - row.bit_count() for row in configuration)
        for multiplicities in _column_distributions(allowed_columns, row_sums):
            columns = tuple(
                column
                for column, multiplicity in zip(
                    allowed_columns, multiplicities, strict=True
                )
                for _ in range(multiplicity)
            )
            if len(columns) != len(ISOLATED_VERTICES):
                raise AssertionError("record-47 isolate-column count changed")
            missing_fives = list(configuration)
            for isolate, column in zip(ISOLATED_VERTICES, columns, strict=True):
                for row in range(4):
                    if column & (1 << row):
                        missing_fives[row] |= 1 << isolate
            canonical_fives = tuple(sorted(missing_fives))
            if any(value.bit_count() != 5 for value in canonical_fives):
                raise AssertionError("record-47 a=0 row size changed")
            if len(set(canonical_fives)) != 4:
                repeated += 1
                continue
            if _has_small_transversal(canonical_fives, graph_covers):
                small_transversal_failures += 1
                continue
            facet_count = sum(
                all(transversal & value for value in canonical_fives)
                for transversal in graph_cover_fives
            )
            states.append((canonical_fives, facet_count))
            if facet_count == 90:
                target_states.append(canonical_fives)

    path_histogram = tuple(
        sorted(Counter(count for _state, count in states).items())
    )
    exact_mask_states = tuple(sorted(set(states)))
    exact_mask_histogram = tuple(
        sorted(Counter(count for _state, count in exact_mask_states).items())
    )
    quotient_path_multiplicities: Counter[
        tuple[tuple[int, ...], tuple[int, ...]]
    ] = Counter()
    quotient_facet_counts: dict[
        tuple[tuple[int, ...], tuple[int, ...]], int
    ] = {}
    for missing_fives, facet_count in states:
        key = _quotient_key(missing_fives)
        quotient_path_multiplicities[key] += 1
        previous = quotient_facet_counts.setdefault(key, facet_count)
        if previous != facet_count:
            raise AssertionError("facet count changed inside an a=0 quotient orbit")
    quotient_states = tuple(sorted(quotient_facet_counts.items()))
    quotient_histogram = tuple(
        sorted(Counter(count for _key, count in quotient_states).items())
    )
    path_multiplicity_histogram = tuple(
        sorted(Counter(quotient_path_multiplicities.values()).items())
    )
    target_exact_states = tuple(
        row for row in exact_mask_states if row[1] == 90
    )
    target_quotient_states = tuple(
        row for row in quotient_states if row[1] == 90
    )
    unique_target = target_states[0] if len(target_states) == 1 else ()
    audit = A0Audit(
        active_vertices=ACTIVE_VERTICES,
        isolated_vertices=ISOLATED_VERTICES,
        minimal_graph_covers=graph_covers,
        independent_active_mask_count=len(independent_active),
        active_row_configuration_count=len(active_configurations),
        column_distribution_problem_count=_column_distributions.cache_info().currsize,
        repeated_row_path_count=repeated,
        enumeration_path_count=len(states),
        exact_mask_state_count=len(exact_mask_states),
        canonical_quotient_orbit_count=len(quotient_states),
        small_transversal_path_count=small_transversal_failures,
        path_facet_count_histogram=path_histogram,
        exact_mask_facet_count_histogram=exact_mask_histogram,
        quotient_facet_count_histogram=quotient_histogram,
        path_multiplicity_per_quotient_orbit_histogram=(
            path_multiplicity_histogram
        ),
        path_transcript_sha256=_state_digest(states),
        exact_mask_transcript_sha256=_state_digest(exact_mask_states),
        quotient_transcript_sha256=_quotient_digest(quotient_states),
        target_path_count=len(target_states),
        target_exact_mask_state_count=len(target_exact_states),
        target_quotient_orbit_count=len(target_quotient_states),
        unique_target_missing_fives=unique_target,
    )
    if not audit.verified:
        raise AssertionError("record-47 a=0 exhaustive census changed")
    return audit


@lru_cache(maxsize=1)
def audit_record47_all_completions(
    schedules_json: str,
) -> Record47AllCompletionsAudit:
    """Audit all record-47 schedules from canonical serialized schedule rows."""

    schedules = json.loads(schedules_json)
    if not isinstance(schedules, list):
        raise ValueError("record-47 schedules must be a JSON list")
    order, edges = decode_graph6(TARGET_GRAPH6)
    if order != VERTEX_COUNT:
        raise AssertionError("record-47 graph order changed")
    audit = Record47AllCompletionsAudit(
        decoded_edges=edges,
        schedule=_audit_schedules(schedules),
        a1=_audit_a1(),
        a0=_audit_a0(),
    )
    if not audit.verified:
        raise AssertionError("record-47 all-completion audit changed")
    return audit


__all__ = [
    "A0Audit",
    "A1Audit",
    "EXPECTED_A0_PATH_FACET_COUNT_HISTOGRAM",
    "EXPECTED_A0_QUOTIENT_FACET_COUNT_HISTOGRAM",
    "EXPECTED_A1_GATE_HISTOGRAM",
    "EXPECTED_SCHEDULES",
    "Record47AllCompletionsAudit",
    "ScheduleAudit",
    "audit_record47_all_completions",
]
