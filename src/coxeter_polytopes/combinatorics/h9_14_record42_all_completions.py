"""Exact whole-shadow compact-Coxeter census for H9/14 record 42.

Record 42 has a five-cycle of missing edges and nine graph isolates.  The
five Dehn--Sommerville schedules are exhausted as follows.

* Strict positivity of the fourteen vertex-link ``g_4`` entries removes the
  ``a=3,4`` schedules.
* The ``a=2`` schedule is removed by the complete unordered pair census of
  locally admissible missing four-faces and their forced five-shadows.
* In the ``a=1`` schedule, all survivors form one ``S_9`` orbit and the last
  minimal five-face is forced.  This is exactly the frozen
  ``C_5 x Delta_3 x Delta_4`` product nerve.
* In the ``a=0`` schedule, a lossless row/column enumeration gives fifteen
  exact facet-count families.  Five pass the full incidence/topology audit
  and they form one exact ``Aut(C_5) x S_9 x S_5`` incidence-matrix orbit.

For extra safety the unrestricted compact Coxeter-label gate is run on all
fifteen ``a=0`` target families, not merely on the five topology survivors.
Every order-five Lanner placement is retained.  Unassigned meeting pairs are
temporarily right angled; Perron--Frobenius monotonicity makes rejection at
that value proof-safe for every later label in ``{3,4,5}``.

The remaining ``a=1`` product nerve is chained, by the certificate builder,
to the re-audited selected-support principal-block no-go.  That dependency
uses the full ``S4 x S5`` normalization and interlacing, so labels incident to
nodes outside the decisive block cannot repair its negative index two.
Nothing here addresses ideal vertices or finite-volume noncompact chambers.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations
import json
from math import comb
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_record46 import (
    _a4_interval_feasible,
    _base_link_g3_g4,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.local_link import (
    analyze_local_link_graph,
)
from coxeter_polytopes.combinatorics.missing_edges import (
    decode_graph6,
    macaulay_raise,
)


TARGET_RECORD_INDEX = 42
TARGET_GRAPH6 = "M?????????O?o?O@?"
TARGET_MISSING_EDGES = (
    (0, 11),
    (0, 12),
    (1, 12),
    (1, 13),
    (11, 13),
)
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
ACTIVE_VERTICES = (0, 1, 11, 12, 13)
ISOLATED_VERTICES = tuple(range(2, 11))
ROW_COUNT = 5
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES)

# The active vertices in cyclic order in the missing-edge graph.
ACTIVE_CYCLE = (0, 11, 13, 1, 12)

# The selected a=1 product completion.
PRODUCT_MISSING_FOUR = (2, 3, 4, 5)
PRODUCT_MISSING_FIVE = (6, 7, 8, 9, 10)

# Canonical representative among the five topology-valid a=0 families.
A0_REPRESENTATIVE_FIVES = (
    (0, 3, 4, 5, 6),
    (1, 3, 4, 5, 6),
    (2, 3, 4, 5, 6),
    (2, 7, 8, 9, 10),
    (7, 8, 9, 10, 12),
)

# (a, total missing five-sets, facet count, summed vertex-link g4).
EXPECTED_SCHEDULES = (
    (0, 5, 102, 30),
    (1, 11, 100, 24),
    (2, 17, 98, 18),
    (3, 23, 96, 12),
    (4, 29, 94, 6),
)
EXPECTED_A1_SHADOW_HISTOGRAM = ((8, 420), (10, 126))
EXPECTED_A2_SHADOW_HISTOGRAM = (
    (15, 420),
    (16, 70_140),
    (17, 2_520),
    (18, 50_400),
    (19, 1_260),
    (20, 6_615),
)
EXPECTED_TARGET_FAMILY_SHA256 = (
    "27412b42ae9cbe2c09e915c70b2986f3d5ca4b90a11b758d0df674727534d9d5"
)
EXPECTED_TOPOLOGY_ORBIT_KEY_SHA256 = (
    "05c1d904e9aab3257c140300f94b7a3d3d13029736f70ab730db0b286cd0b9e2"
)
EXPECTED_TOPOLOGY_VALID_FAMILY_INDICES = (0, 1, 2, 3, 6)
EXPECTED_A0_FACET_HISTOGRAM = (
    (40, 10), (45, 240), (50, 2_085), (51, 30), (53, 180),
    (54, 180), (55, 7_100), (56, 980), (57, 900), (58, 3_690),
    (59, 2_160), (60, 20_960), (61, 5_600), (62, 10_010),
    (63, 14_250), (64, 15_590), (65, 40_640), (66, 20_855),
    (67, 32_460), (68, 33_100), (69, 38_670), (70, 63_375),
    (71, 52_490), (72, 53_040), (73, 68_875), (74, 55_505),
    (75, 94_675), (76, 47_270), (77, 85_790), (78, 41_275),
    (79, 67_915), (80, 34_715), (81, 50_905), (82, 19_205),
    (83, 29_590), (84, 15_825), (85, 16_245), (86, 4_225),
    (87, 5_720), (88, 5_100), (89, 3_740), (90, 960), (91, 1_050),
    (92, 725), (93, 555), (94, 265), (95, 610), (96, 80), (97, 45),
    (98, 90), (100, 25), (101, 30), (102, 35), (105, 10),
)
EXPECTED_LABEL_STAGE_TRANSCRIPT_SHA256 = (
    "5c9bab8a166c58bdbaa1e3179d900a42059ba53ff2f3aa15dbdb349e882d3554"
)


def mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex))


def _canonical_digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _rank_nine_vertices(
    minimal_nonfaces: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(frozenset(face) for face in minimal_nonfaces)
    return tuple(
        chosen
        for chosen in combinations(range(TARGET_VERTEX_COUNT), TARGET_DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


def _minimal_graph_covers() -> tuple[int, ...]:
    answers = []
    for value in range(1, 1 << TARGET_VERTEX_COUNT):
        if not all(value & edge for edge in EDGE_MASKS):
            continue
        if all(
            not all((value ^ (1 << vertex)) & edge for edge in EDGE_MASKS)
            for vertex in vertices(value)
        ):
            answers.append(value)
    return tuple(answers)


EXPECTED_MINIMAL_GRAPH_COVERS = (2051, 6146, 8195, 12289, 14336)


def _local_a5_interval(subset: tuple[int, ...]) -> tuple[int, int]:
    """Strict local/Macaulay interval for A5 incidence through ``subset``."""

    size = len(subset)
    link_dimension = TARGET_DIMENSION - size
    high_order = 5 - size
    rest = tuple(vertex for vertex in range(TARGET_VERTEX_COUNT) if vertex not in subset)

    def independent(row: Sequence[int]) -> bool:
        value = mask(row)
        return all(value & edge != edge for edge in EDGE_MASKS)

    link_f = {-1: 1}
    for order in range(1, high_order + 1):
        link_f[order - 1] = sum(
            independent(subset + extension)
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
    previous, baseline = g[-2], g[-1]
    lower = max(0, baseline - macaulay_raise(previous, high_order - 1))
    upper = baseline - 1
    if previous < 1 or baseline < 1 or lower > upper:
        raise AssertionError("record-42 local A5 interval became infeasible")
    return lower, upper


def _active_row_configurations() -> tuple[tuple[int, ...], ...]:
    """All unordered disjoint independent active-row patterns.

    Strict singleton link intervals put each active vertex in at most one A5.
    A row must be independent in the missing-edge C5, hence has zero, one, or
    two active vertices.  Taking all multisets of five such rows and enforcing
    disjointness is a lossless row-permutation quotient.  Pair/triple local
    intervals are deliberately *not* used, so this is a proof-safe superset.
    """

    active_masks = tuple(
        value
        for order in range(3)
        for chosen in combinations(ACTIVE_VERTICES, order)
        for value in (mask(chosen),)
        if all(value & edge != edge for edge in EDGE_MASKS)
    )
    if len(active_masks) != 11:
        raise AssertionError("record-42 active C5 independent-mask count changed")
    configurations = tuple(
        rows
        for rows in combinations_with_replacement(active_masks, ROW_COUNT)
        if all(not (left & right) for left, right in combinations(rows, 2))
    )
    if len(configurations) != 82:
        raise AssertionError("record-42 active-row quotient count changed")
    return configurations


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


def _forced_shadow(family: Sequence[int], fives: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        five for five in fives if any(five & four == four for four in family)
    )


def _shadow_gate(
    family: Sequence[int],
    shadow: Sequence[int],
    base_g3: Sequence[int],
    base_g4: Sequence[int],
    total_fives: int,
) -> tuple[bool, bool, int, tuple[int, ...]]:
    additional = total_fives - len(shadow)
    if additional < 0:
        return False, False, additional, ()
    lower: list[int] = []
    upper: list[int] = []
    strict = True
    for vertex in range(TARGET_VERTEX_COUNT):
        a4_incidence = sum(bool(face & (1 << vertex)) for face in family)
        previous = base_g3[vertex] - a4_incidence
        forced_incidence = sum(bool(face & (1 << vertex)) for face in shadow)
        before = base_g4[vertex] + 6 * a4_incidence - forced_incidence
        if previous < 1 or before < 1:
            strict = False
        lower.append(
            max(0, before - macaulay_raise(previous, 3)) if previous >= 0 else 0
        )
        upper.append(before - 1)
    aggregate = (
        strict
        and all(value >= 0 for value in upper)
        and sum(lower) <= 5 * additional <= sum(upper)
    )
    return strict, aggregate, additional, tuple(upper)


@dataclass(frozen=True, slots=True)
class ScheduleAudit:
    rows: tuple[tuple[int, int, int, int], ...]
    surviving_a: tuple[int, ...]
    excluded_a: tuple[int, ...]

    @property
    def verified(self) -> bool:
        return (
            self.rows == EXPECTED_SCHEDULES
            and self.surviving_a == (0, 1, 2)
            and self.excluded_a == (3, 4)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "rows": [
                {
                    "missing_four_count": a,
                    "total_missing_five_sets": b,
                    "facet_count": facets,
                    "summed_vertex_link_g4": link_sum,
                }
                for a, b, facets, link_sum in self.rows
            ],
            "strict_link_sum_required": TARGET_VERTEX_COUNT,
            "surviving_missing_four_counts": list(self.surviving_a),
            "excluded_missing_four_counts": list(self.excluded_a),
            "verified": self.verified,
        }


@dataclass(frozen=True, slots=True)
class PositiveBranchAudit:
    individually_feasible_a4_count: int
    a1_shadow_histogram: tuple[tuple[int, int], ...]
    a1_strict_count: int
    a1_forced_product_count: int
    a1_incidence_orbit_count: int
    a1_representative_four: tuple[int, ...]
    a1_representative_added_five: tuple[int, ...]
    a2_unordered_pair_count: int
    a2_family_interval_count: int
    a2_shadow_histogram: tuple[tuple[int, int], ...]
    a2_forced_shadow_overrun_count: int
    a2_strict_count: int

    @property
    def verified(self) -> bool:
        return (
            self.individually_feasible_a4_count == 546
            and self.a1_shadow_histogram == EXPECTED_A1_SHADOW_HISTOGRAM
            and self.a1_strict_count == 126
            and self.a1_forced_product_count == 126
            and self.a1_incidence_orbit_count == 1
            and self.a1_representative_four == PRODUCT_MISSING_FOUR
            and self.a1_representative_added_five == PRODUCT_MISSING_FIVE
            and self.a2_unordered_pair_count == 148_785
            and self.a2_family_interval_count == 131_355
            and self.a2_shadow_histogram == EXPECTED_A2_SHADOW_HISTOGRAM
            and self.a2_forced_shadow_overrun_count == 58_275
            and self.a2_strict_count == 0
        )

    def to_record(self) -> dict[str, object]:
        return {
            "individually_interval_feasible_missing_fours": (
                self.individually_feasible_a4_count
            ),
            "a1": {
                "forced_shadow_size_histogram": dict(self.a1_shadow_histogram),
                "strict_shadow_survivors": self.a1_strict_count,
                "forced_product_completions": self.a1_forced_product_count,
                "exact_incidence_orbits_under_S9": self.a1_incidence_orbit_count,
                "representative_missing_four": list(self.a1_representative_four),
                "forced_additional_minimal_five": list(
                    self.a1_representative_added_five
                ),
                "equals_frozen_C5_x_Delta3_x_Delta4_nerve": True,
            },
            "a2": {
                "all_unordered_distinct_pairs": self.a2_unordered_pair_count,
                "family_interval_feasible_pairs": self.a2_family_interval_count,
                "forced_shadow_size_histogram": dict(self.a2_shadow_histogram),
                "forced_shadow_count_overrun_pairs": (
                    self.a2_forced_shadow_overrun_count
                ),
                "strict_shadow_survivors": self.a2_strict_count,
                "branch_excluded": self.a2_strict_count == 0,
            },
            "verified": self.verified,
        }


def _audit_positive_branches() -> PositiveBranchAudit:
    local = analyze_local_link_graph(TARGET_GRAPH6)
    fours = tuple(
        face
        for face in independent_masks(TARGET_GRAPH6, 4)
        if _a4_interval_feasible(face, local)
    )
    fives = independent_masks(TARGET_GRAPH6, 5)
    base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)

    a1_shadow_histogram: Counter[int] = Counter()
    strict_fours: list[int] = []
    product_completions: list[tuple[int, int]] = []
    for four in fours:
        shadow = _forced_shadow((four,), fives)
        strict, aggregate, additional, upper = _shadow_gate(
            (four,), shadow, base_g3, base_g4, 11
        )
        a1_shadow_histogram[len(shadow)] += 1
        if not (strict and aggregate):
            continue
        strict_fours.append(four)
        candidates = tuple(
            five
            for five in fives
            if five not in shadow
            and all(
                not (five & (1 << vertex)) or upper[vertex] >= 1
                for vertex in range(TARGET_VERTEX_COUNT)
            )
        )
        if additional != 1 or len(candidates) != 1:
            raise AssertionError("record-42 a=1 final A5 was not forced")
        product_completions.append((four, candidates[0]))

    expected_product_pairs = {
        (mask(four), mask(set(ISOLATED_VERTICES) - set(four)))
        for four in combinations(ISOLATED_VERTICES, 4)
    }
    if set(product_completions) != expected_product_pairs:
        raise AssertionError("record-42 a=1 survivors are not the product orbit")

    a2_local = 0
    a2_overrun = 0
    a2_strict = 0
    a2_shadow_histogram: Counter[int] = Counter()
    for family in combinations(fours, 2):
        if not _family_interval_feasible(family, local):
            continue
        a2_local += 1
        shadow = _forced_shadow(family, fives)
        strict, _aggregate, additional, _upper = _shadow_gate(
            family, shadow, base_g3, base_g4, 17
        )
        a2_shadow_histogram[len(shadow)] += 1
        a2_overrun += additional < 0
        a2_strict += strict

    audit = PositiveBranchAudit(
        individually_feasible_a4_count=len(fours),
        a1_shadow_histogram=tuple(sorted(a1_shadow_histogram.items())),
        a1_strict_count=len(strict_fours),
        a1_forced_product_count=len(product_completions),
        a1_incidence_orbit_count=int(bool(product_completions)),
        a1_representative_four=PRODUCT_MISSING_FOUR,
        a1_representative_added_five=PRODUCT_MISSING_FIVE,
        a2_unordered_pair_count=comb(len(fours), 2),
        a2_family_interval_count=a2_local,
        a2_shadow_histogram=tuple(sorted(a2_shadow_histogram.items())),
        a2_forced_shadow_overrun_count=a2_overrun,
        a2_strict_count=a2_strict,
    )
    if not audit.verified:
        raise AssertionError("record-42 positive-schedule census changed")
    return audit


def _distribution_iterator(
    allowed: tuple[int, ...], row_sums: tuple[int, ...]
):
    counts = [0] * len(allowed)

    def visit(position: int, columns_left: int, sums: tuple[int, ...]):
        if position == len(allowed):
            if columns_left == 0 and not any(sums):
                yield tuple(counts)
            return
        column = allowed[position]
        maximum = columns_left
        for row in range(ROW_COUNT):
            if column & (1 << row):
                maximum = min(maximum, sums[row])
        for multiplicity in range(maximum + 1):
            next_left = columns_left - multiplicity
            next_sums = tuple(
                sums[row] - multiplicity * bool(column & (1 << row))
                for row in range(ROW_COUNT)
            )
            if all(0 <= value <= next_left for value in next_sums):
                counts[position] = multiplicity
                yield from visit(position + 1, next_left, next_sums)
        counts[position] = 0

    yield from visit(0, len(ISOLATED_VERTICES), row_sums)


@lru_cache(maxsize=None)
def _column_distributions(
    allowed: tuple[int, ...], row_sums: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    return tuple(_distribution_iterator(allowed, row_sums))


def _faces_from_distribution(
    active_rows: Sequence[int],
    allowed: Sequence[int],
    multiplicities: Sequence[int],
) -> tuple[int, ...]:
    columns = tuple(
        column
        for column, multiplicity in zip(allowed, multiplicities, strict=True)
        for _ in range(multiplicity)
    )
    if len(columns) != len(ISOLATED_VERTICES):
        raise AssertionError("record-42 isolate-column count changed")
    faces = list(active_rows)
    for isolate, column in zip(ISOLATED_VERTICES, columns, strict=True):
        for row in range(ROW_COUNT):
            if column & (1 << row):
                faces[row] |= 1 << isolate
    return tuple(sorted(faces))


def _hit_word(face: int, covers: Sequence[int]) -> int:
    return sum(bool(face & cover) << index for index, cover in enumerate(covers))


def _has_small_transversal(fives: Sequence[int], graph_covers: Sequence[int]) -> bool:
    for cover in graph_covers:
        missed = tuple(five for five in fives if not five & cover)
        if not missed:
            return True
        common = missed[0]
        for five in missed[1:]:
            common &= five
        if common:
            return True
    return False


def _minimal_transversals(hyperedges: Sequence[int]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, FULL_MASK + 1):
        if not all(candidate & edge for edge in hyperedges):
            continue
        if any(
            all((candidate ^ (1 << vertex)) & edge for edge in hyperedges)
            for vertex in vertices(candidate)
        ):
            continue
        answers.append(candidate)
    return tuple(answers)


def _cycle_automorphisms() -> tuple[tuple[int, ...], ...]:
    answers = []
    for direction in (1, -1):
        for shift in range(5):
            image = list(range(TARGET_VERTEX_COUNT))
            for index, source in enumerate(ACTIVE_CYCLE):
                image[source] = ACTIVE_CYCLE[(shift + direction * index) % 5]
            answers.append(tuple(image))
    return tuple(answers)


def _orbit_key(family: Sequence[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Exact canonical key under D5 active, S9 isolates, and S5 rows.

    Once the row order and active automorphism are fixed, two incidence
    matrices are in the same isolate orbit iff their nine column multisets
    agree.  Exhausting all 10*120 outer choices therefore proves both
    directions of the quotient.
    """

    active_mask = mask(ACTIVE_VERTICES)
    keys = []
    for automorphism in _cycle_automorphisms():
        transformed = tuple(
            sum(1 << automorphism[vertex] for vertex in vertices(face))
            for face in family
        )
        for order in permutations(range(ROW_COUNT)):
            active_rows = tuple(transformed[index] & active_mask for index in order)
            columns = tuple(
                sorted(
                    sum(
                        bool(transformed[order[row]] & (1 << isolate)) << row
                        for row in range(ROW_COUNT)
                    )
                    for isolate in ISOLATED_VERTICES
                )
            )
            keys.append((active_rows, columns))
    return min(keys)


LabelState = tuple[tuple[Pair, int], ...]


@lru_cache(maxsize=None)
def _placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(support))
    )


@lru_cache(maxsize=262_144)
def _finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return is_finite_coxeter_labeling(range(order), labels)


def _restriction_finite(restriction: Sequence[int], labels: Mapping[Pair, int]) -> bool:
    code = tuple(
        labels[pair(left, right)] for left, right in combinations(restriction, 2)
    )
    return _finite_code(len(restriction), code)


def _label_state_digest(states: Sequence[LabelState]) -> str:
    return _canonical_digest(
        [
            [[left, right, label] for (left, right), label in state]
            for state in sorted(states)
        ]
    )


@dataclass(frozen=True, slots=True)
class LabelStage:
    support_index: int
    support: tuple[int, ...]
    cluster_union: tuple[int, ...]
    restriction_size_histogram: tuple[tuple[int, int], ...]
    raw_extension_count: int
    compatible_extension_count: int
    passing_state_count: int
    passing_state_sha256: str

    def to_record(self) -> dict[str, object]:
        return {
            "support_index": self.support_index,
            "support": list(self.support),
            "complete_lanner_placement_count": len(_placement_options(self.support)),
            "cluster_union": list(self.cluster_union),
            "actual_vertex_restriction_size_histogram": dict(
                self.restriction_size_histogram
            ),
            "raw_extension_count": self.raw_extension_count,
            "compatible_extension_count": self.compatible_extension_count,
            "passing_state_count": self.passing_state_count,
            "passing_state_sha256": self.passing_state_sha256,
        }


@dataclass(frozen=True, slots=True)
class FamilyLabelAudit:
    family_index: int
    missing_fives: tuple[tuple[int, ...], ...]
    rank_nine_vertex_count: int
    stages: tuple[LabelStage, ...]

    @property
    def compact_no_go(self) -> bool:
        return bool(self.stages) and self.stages[-1].passing_state_count == 0

    def to_record(self) -> dict[str, object]:
        return {
            "family_index": self.family_index,
            "missing_fives": [list(face) for face in self.missing_fives],
            "rank_nine_vertex_count": self.rank_nine_vertex_count,
            "decisive_support_indices": [stage.support_index for stage in self.stages],
            "stages": [stage.to_record() for stage in self.stages],
            "unrestricted_compact_coxeter_no_go": self.compact_no_go,
        }


def _next_support(
    supports: Sequence[tuple[int, ...]], remaining: set[int], union: set[int]
) -> int:
    if not union:
        return min(remaining)
    return max(
        remaining,
        key=lambda index: (
            len(set(supports[index]) & union),
            -len(set(supports[index]) - union),
            -index,
        ),
    )


def _audit_family_labels(
    family_index: int, missing_five_masks: Sequence[int]
) -> FamilyLabelAudit:
    missing_fives = tuple(vertices(value) for value in missing_five_masks)
    nonfaces = TARGET_MISSING_EDGES + missing_fives
    rank_nine = _rank_nine_vertices(nonfaces)
    supports = missing_fives
    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[LabelState, ...] = ((),)
    stages = []
    while states and remaining:
        support_index = _next_support(supports, remaining, union)
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {
                    tuple(facet for facet in vertex if facet in union)
                    for vertex in rank_nine
                },
                key=lambda row: (len(row), row),
            )
        )
        if any(
            set(edge) <= set(restriction)
            for edge in TARGET_MISSING_EDGES
            for restriction in restrictions
        ):
            raise AssertionError("a record-42 actual vertex contains a dotted pair")
        options = _placement_options(support)
        passing: set[LabelState] = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in options:
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(_restriction_finite(row, labels) for row in restrictions):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            LabelStage(
                support_index=support_index,
                support=support,
                cluster_union=cluster,
                restriction_size_histogram=tuple(
                    sorted(Counter(map(len, restrictions)).items())
                ),
                raw_extension_count=(
                    (stages[-1].passing_state_count if stages else 1) * len(options)
                ),
                compatible_extension_count=compatible,
                passing_state_count=len(states),
                passing_state_sha256=_label_state_digest(states),
            )
        )
        remaining.remove(support_index)
    audit = FamilyLabelAudit(
        family_index=family_index,
        missing_fives=missing_fives,
        rank_nine_vertex_count=len(rank_nine),
        stages=tuple(stages),
    )
    if not audit.compact_no_go:
        raise AssertionError(f"record-42 a=0 family {family_index} survived labels")
    return audit


@dataclass(frozen=True, slots=True)
class A0FamilyAudit:
    family_index: int
    missing_fives: tuple[int, ...]
    minimal_transversal_histogram: tuple[tuple[int, int], ...]
    facet_complement_sha256: str
    packed_global_audit: tuple[tuple[str, object], ...]

    @property
    def global_audit(self) -> dict[str, object]:
        return dict(self.packed_global_audit)

    @property
    def topology_valid(self) -> bool:
        return bool(self.global_audit["passes_connected_homology_audit"])

    def to_record(self) -> dict[str, object]:
        return {
            "family_index": self.family_index,
            "missing_fives": [list(vertices(value)) for value in self.missing_fives],
            "minimal_transversal_size_histogram": dict(
                self.minimal_transversal_histogram
            ),
            "facet_complement_sha256": self.facet_complement_sha256,
            "topology_valid": self.topology_valid,
            "global_selection_audit": self.global_audit,
        }


@dataclass(frozen=True, slots=True)
class A0Audit:
    active_row_configuration_count: int
    graph_cover_prefilter_path_count: int
    prefilter_repeated_row_path_count: int
    small_transversal_failure_count: int
    blocker_valid_path_count: int
    facet_count_histogram: tuple[tuple[int, int], ...]
    target_path_count: int
    exact_target_family_count: int
    target_family_sha256: str
    family_audits: tuple[A0FamilyAudit, ...]
    topology_valid_family_indices: tuple[int, ...]
    topology_orbit_count: int
    topology_orbit_key_sha256: str
    label_audits: tuple[FamilyLabelAudit, ...]
    label_stage_transcript_sha256: str

    @property
    def verified(self) -> bool:
        return (
            self.active_row_configuration_count == 82
            and self.graph_cover_prefilter_path_count == 1_098_455
            and self.prefilter_repeated_row_path_count == 28_805
            and self.small_transversal_failure_count == 0
            and self.blocker_valid_path_count == 1_069_650
            and self.facet_count_histogram == EXPECTED_A0_FACET_HISTOGRAM
            and self.target_path_count == 35
            and self.exact_target_family_count == 15
            and self.target_family_sha256 == EXPECTED_TARGET_FAMILY_SHA256
            and len(self.family_audits) == 15
            and self.topology_valid_family_indices
            == EXPECTED_TOPOLOGY_VALID_FAMILY_INDICES
            and self.topology_orbit_count == 1
            and self.topology_orbit_key_sha256
            == EXPECTED_TOPOLOGY_ORBIT_KEY_SHA256
            and len(self.label_audits) == 15
            and all(audit.compact_no_go for audit in self.label_audits)
            and self.label_stage_transcript_sha256
            == EXPECTED_LABEL_STAGE_TRANSCRIPT_SHA256
        )

    def to_record(self) -> dict[str, object]:
        return {
            "enumeration_scope": {
                "active_row_configuration_count": (
                    self.active_row_configuration_count
                ),
                "active_row_configuration_definition": (
                    "all unordered multisets of five pairwise-disjoint "
                    "independent subsets of the active missing-edge C5; "
                    "there are 11 possible row masks and 82 configurations"
                ),
                "exact_graph_cover_prefilter_paths": (
                    self.graph_cover_prefilter_path_count
                ),
                "prefilter_paths_with_repeated_rows": (
                    self.prefilter_repeated_row_path_count
                ),
                "small_transversal_failure_count": (
                    self.small_transversal_failure_count
                ),
                "blocker_valid_paths": self.blocker_valid_path_count,
                "losslessness": (
                    "strict singleton link intervals put every active vertex "
                    "in at most one A5; graph independence and alpha(C5)=2 "
                    "give the 82 active configurations; isolate columns are "
                    "then enumerated by their complete S9 multiset"
                ),
            },
            "facet_count_histogram_on_blocker_valid_paths": dict(
                self.facet_count_histogram
            ),
            "scheduled_facet_count": 102,
            "scheduled_target_path_count": self.target_path_count,
            "exact_unordered_target_family_count": self.exact_target_family_count,
            "target_family_sha256": self.target_family_sha256,
            "target_family_audits": [audit.to_record() for audit in self.family_audits],
            "topology_valid_family_indices": list(
                self.topology_valid_family_indices
            ),
            "topology_valid_family_count": len(self.topology_valid_family_indices),
            "exact_AutC5_x_S9_x_S5_orbit_count": self.topology_orbit_count,
            "topology_orbit_key_sha256": self.topology_orbit_key_sha256,
            "orbit_quotient_proof": (
                "all 10 active-cycle automorphisms and all 120 row "
                "permutations are explicit; for each choice, equality modulo "
                "S9 is equivalent to equality of the nine isolate-column "
                "multisets"
            ),
            "representative_missing_fives": [
                list(face) for face in A0_REPRESENTATIVE_FIVES
            ],
            "unrestricted_label_audits_for_all_15_targets": [
                audit.to_record() for audit in self.label_audits
            ],
            "label_stage_transcript_sha256": (
                self.label_stage_transcript_sha256
            ),
            "verified": self.verified,
        }


def _pack_global_audit(record: Mapping[str, object]) -> tuple[tuple[str, object], ...]:
    keys = (
        "facet_count",
        "face_vector",
        "ridge_degree_histogram",
        "all_ridge_degrees_two",
        "dual_component_count",
        "dual_graph_connected",
        "covers_all_independent_triples",
        "missing_four_count",
        "largest_minimal_nonface_order",
        "matches_dehn_sommerville_schedule",
        "strict_local_g_constraints_hold",
        "local_macaulay_constraints_hold",
        "all_links_have_spherical_euler_characteristic",
        "all_links_are_mod2_homology_spheres",
        "passes_selection_model",
        "passes_connected_homology_audit",
    )
    return tuple((key, record[key]) for key in keys)


def _audit_a0(
    support_core_masks: Sequence[int], schedules: Sequence[Mapping[str, object]]
) -> A0Audit:
    graph_covers = _minimal_graph_covers()
    if graph_covers != EXPECTED_MINIMAL_GRAPH_COVERS:
        raise AssertionError("record-42 minimal graph covers changed")
    five_covers = graph_vertex_cover_masks(TARGET_GRAPH6, 5)
    all_five = (1 << len(five_covers)) - 1
    five_hit_words: dict[int, int] = {}
    prefilter_paths = 0
    repeated = 0
    blocker_paths = 0
    small_transversal_failures = 0
    histogram: Counter[int] = Counter()
    target_paths: list[tuple[int, ...]] = []

    active_configurations = _active_row_configurations()
    if tuple(_local_a5_interval((vertex,)) for vertex in ACTIVE_VERTICES) != (
        (0, 1),
    ) * 5:
        raise AssertionError("record-42 active singleton A5 intervals changed")
    if tuple(_local_a5_interval((vertex,)) for vertex in ISOLATED_VERTICES) != (
        (0, 4),
    ) * 9:
        raise AssertionError("record-42 isolate singleton A5 intervals changed")

    for active_rows in active_configurations:
        missed_rows = tuple(
            sum(
                1 << row
                for row, active in enumerate(active_rows)
                if not active & cover
            )
            for cover in graph_covers
        )
        allowed = tuple(
            column
            for column in range(1 << ROW_COUNT)
            if column.bit_count() <= 4
            and all(column & missed != missed for missed in missed_rows)
        )
        row_sums = tuple(5 - row.bit_count() for row in active_rows)
        for distribution in _column_distributions(allowed, row_sums):
            prefilter_paths += 1
            fives = _faces_from_distribution(active_rows, allowed, distribution)
            if len(set(fives)) != ROW_COUNT:
                repeated += 1
                continue
            if _has_small_transversal(fives, graph_covers):
                small_transversal_failures += 1
                continue
            blocker_paths += 1
            word = all_five
            for five in fives:
                if five not in five_hit_words:
                    five_hit_words[five] = _hit_word(five, five_covers)
                word &= five_hit_words[five]
            facet_count = word.bit_count()
            histogram[facet_count] += 1
            if facet_count == 102:
                target_paths.append(fives)

    if small_transversal_failures:
        raise AssertionError("record-42 graph-cover prefilter ceased to be exact")
    target_families = tuple(sorted(set(target_paths)))
    family_audits = []
    label_audits = []
    topology_keys = []
    topology_indices = []
    for family_index, family in enumerate(target_families):
        transversals = _minimal_transversals(EDGE_MASKS + family)
        transversal_histogram = tuple(
            sorted(Counter(value.bit_count() for value in transversals).items())
        )
        five_transversals = tuple(
            value for value in transversals if value.bit_count() == 5
        )
        global_audit = audit_global_facet_selection(
            TARGET_GRAPH6, support_core_masks, schedules, five_transversals
        )
        packed = _pack_global_audit(global_audit.to_record())
        family_audit = A0FamilyAudit(
            family_index=family_index,
            missing_fives=family,
            minimal_transversal_histogram=transversal_histogram,
            facet_complement_sha256=sha256(
                b"".join(value.to_bytes(2, "little") for value in five_transversals)
            ).hexdigest(),
            packed_global_audit=packed,
        )
        family_audits.append(family_audit)
        if family_audit.topology_valid:
            topology_indices.append(family_index)
            topology_keys.append(_orbit_key(family))
        label_audits.append(_audit_family_labels(family_index, family))

    representative_masks = tuple(mask(face) for face in A0_REPRESENTATIVE_FIVES)
    representative_index = target_families.index(tuple(sorted(representative_masks)))
    if representative_index not in topology_indices:
        raise AssertionError("the specified record-42 a=0 representative is not valid")
    orbit_keys = tuple(sorted(set(topology_keys)))
    audit = A0Audit(
        active_row_configuration_count=len(active_configurations),
        graph_cover_prefilter_path_count=prefilter_paths,
        prefilter_repeated_row_path_count=repeated,
        small_transversal_failure_count=small_transversal_failures,
        blocker_valid_path_count=blocker_paths,
        facet_count_histogram=tuple(sorted(histogram.items())),
        target_path_count=len(target_paths),
        exact_target_family_count=len(target_families),
        target_family_sha256=_canonical_digest(
            [[list(vertices(face)) for face in family] for family in target_families]
        ),
        family_audits=tuple(family_audits),
        topology_valid_family_indices=tuple(topology_indices),
        topology_orbit_count=len(orbit_keys),
        topology_orbit_key_sha256=_canonical_digest(orbit_keys),
        label_audits=tuple(label_audits),
        label_stage_transcript_sha256=_canonical_digest(
            [audit.to_record() for audit in label_audits]
        ),
    )
    if not audit.verified:
        raise AssertionError("record-42 a=0 exhaustive census changed")
    return audit


@dataclass(frozen=True, slots=True)
class Record42AllCompletionsAudit:
    decoded_edges: tuple[tuple[int, int], ...]
    schedule: ScheduleAudit
    positive: PositiveBranchAudit
    a0: A0Audit

    @property
    def verified(self) -> bool:
        return (
            self.decoded_edges == TARGET_MISSING_EDGES
            and self.schedule.verified
            and self.positive.verified
            and self.a0.verified
        )

    def to_record(self) -> dict[str, object]:
        return {
            "record_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
            "decoded_missing_edges": [list(edge) for edge in self.decoded_edges],
            "schedule_audit": self.schedule.to_record(),
            "positive_branch_audit": self.positive.to_record(),
            "a0_audit": self.a0.to_record(),
            "all_incidence_branches_exhausted": self.verified,
            "verified": self.verified,
        }


@lru_cache(maxsize=1)
def audit_record42_all_completions(
    support_record_json: str,
) -> Record42AllCompletionsAudit:
    record = json.loads(support_record_json)
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("record-42 support record changed")
    schedules = record.get("dehn_sommerville_schedule")
    support_core = record.get("ridge_support_core_masks")
    if not isinstance(schedules, list) or not isinstance(support_core, list):
        raise ValueError("record-42 support record is incomplete")
    rows = tuple(
        (
            int(schedule["missing_four_count"]),
            int(schedule["missing_five_count_offset"]),
            int(schedule["f_vector"][-1]),
            summed_vertex_link_g4(schedule["f_vector"]),
        )
        for schedule in schedules
    )
    schedule_audit = ScheduleAudit(
        rows=rows,
        surviving_a=tuple(a for a, _b, _f, total in rows if total >= 14),
        excluded_a=tuple(a for a, _b, _f, total in rows if total < 14),
    )
    audit = Record42AllCompletionsAudit(
        decoded_edges=decode_graph6(TARGET_GRAPH6)[1],
        schedule=schedule_audit,
        positive=_audit_positive_branches(),
        a0=_audit_a0(support_core, schedules),
    )
    if not audit.verified:
        raise AssertionError("record-42 whole-shadow audit changed")
    return audit


__all__ = [
    "A0_REPRESENTATIVE_FIVES",
    "ACTIVE_VERTICES",
    "EXPECTED_A0_FACET_HISTOGRAM",
    "EXPECTED_A2_SHADOW_HISTOGRAM",
    "EXPECTED_SCHEDULES",
    "PRODUCT_MISSING_FIVE",
    "PRODUCT_MISSING_FOUR",
    "Record42AllCompletionsAudit",
    "TARGET_GRAPH6",
    "TARGET_MISSING_EDGES",
    "TARGET_RECORD_INDEX",
    "audit_record42_all_completions",
]
