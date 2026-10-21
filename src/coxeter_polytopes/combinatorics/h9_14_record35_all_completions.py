"""Exact whole-shadow compact-Coxeter census for H9/14 record 35.

All five Dehn--Sommerville schedules are covered.  Strict positivity of the
fourteen vertex-link ``g_4`` entries removes ``a=3,4``.  Complete missing-four
and forced-shadow censuses remove ``a=2`` and leave one incidence completion
in ``a=1``. The ``a=0`` branch is generated from *every* independent active
row mask allowed by the exact incidence caps (not from a singleton-row
normalization), followed by a lossless isolate-column multiset recursion. Its
eighteen exact target families are
all rejected by complete Lanner-placement joins and actual-vertex tests.

The remaining ``a=1`` nerve is rejected by an unrestricted cross-block gate.
The gate literally retains all 108 order-four and 420 order-five Lanner
placements.  Every locally possible non-right edge between its disjoint
five-node Lanner block and the six-node overlapping-Lanner cluster either is
locally impossible or lies in an exact ten-node principal block of negative
index two.  Cauchy interlacing then excludes a compact H9 Gram matrix.

The theorem is compact-only.  Ideal faces, parabolic diagrams, and Vinberg's
finite-volume endpoint criterion are not modeled here.
"""

from __future__ import annotations

from collections import Counter, defaultdict
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
from coxeter_polytopes.combinatorics.h9_14_record35 import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_RECORD_INDEX,
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
from coxeter_polytopes.combinatorics.h9_14_selected_nerve_coxeter_no_go import (
    ExactIndexTwoWitness,
    _a5_witnesses,
    _canonical_lanner_states,
)
from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6, macaulay_raise


VERTEX_COUNT = 14
DIMENSION = 9
FULL_MASK = (1 << VERTEX_COUNT) - 1
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in TARGET_MISSING_EDGES)
ACTIVE_VERTICES = (0, 1, 2, 11, 12, 13)
ISOLATED_VERTICES = tuple(range(3, 11))
ROW_COUNT = 5

EXPECTED_SCHEDULES = (
    (0, 5, 102, 30),
    (1, 11, 100, 24),
    (2, 17, 98, 18),
    (3, 23, 96, 12),
    (4, 29, 94, 6),
)

A1_REPRESENTATIVE_FOUR = (2, 3, 4, 5)
A1_REPRESENTATIVE_FIVES = (
    (3, 4, 5, 11, 12),
    (6, 7, 8, 9, 10),
)
A1_OTHER_REPRESENTATIVE_FOUR = (2, 3, 4, 11)
A2_REPRESENTATIVE_FOURS = ((2, 3, 4, 11), (2, 3, 4, 12))

A1_Q_SUPPORT = A1_REPRESENTATIVE_FOUR
A1_R_SUPPORT = A1_REPRESENTATIVE_FIVES[0]
A1_B_SUPPORT = A1_REPRESENTATIVE_FIVES[1]
A1_LEFT_UNION = tuple(sorted(set(A1_Q_SUPPORT) | set(A1_R_SUPPORT)))
A1_BLOCK_UNION = tuple(sorted(set(A1_LEFT_UNION) | set(A1_B_SUPPORT)))
A1_CROSS_PAIRS = tuple(pair(left, right) for left in A1_LEFT_UNION for right in A1_B_SUPPORT)
A1_INTERNAL_UNFORCED_PAIRS = tuple(
    pair(left, right)
    for left in set(A1_Q_SUPPORT) - set(A1_R_SUPPORT)
    for right in set(A1_R_SUPPORT) - set(A1_Q_SUPPORT)
)

# Filled only after the deterministic production transcript has been replayed.
EXPECTED_A0_PATH_TRANSCRIPT_SHA256 = (
    "28bedd57cb11f4217fb1e68ff661f1a95c6cc6626efe7277a902ae3aba671fe1"
)
EXPECTED_A0_TARGET_FAMILY_SHA256 = (
    "a3678a721ea0e69b1fa82381eaa4625e996c826c06f52c168d94f20a962671d1"
)
EXPECTED_A0_LABEL_TRANSCRIPT_SHA256 = (
    "136879027cfc033d6eb91625d616d34b47a35b2aba955db40f5663dc50e336b2"
)
EXPECTED_A1_LABEL_TRANSCRIPT_SHA256 = (
    "45b32fb880ba66a77202682b5db5532a0ba055017274ceb44e99a3e35cc2804e"
)


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


EXPECTED_INDEPENDENT_ACTIVE_MASKS = tuple(
    mask(values)
    for values in (
        (),
        (0,),
        (1,),
        (2,),
        (11,),
        (12,),
        (0, 1),
        (0, 2),
        (0, 12),
        (1, 2),
        (1, 11),
        (2, 11),
        (2, 12),
        (11, 12),
        (0, 1, 2),
        (0, 2, 12),
        (1, 2, 11),
        (2, 11, 12),
    )
)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1)


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
        for chosen in combinations(range(VERTEX_COUNT), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


def _minimal_graph_covers() -> tuple[int, ...]:
    answers = []
    for value in range(1, 1 << VERTEX_COUNT):
        if not all(value & edge for edge in EDGE_MASKS):
            continue
        if all(
            not all((value ^ (1 << vertex)) & edge for edge in EDGE_MASKS)
            for vertex in vertices(value)
        ):
            answers.append(value)
    return tuple(answers)


EXPECTED_MINIMAL_GRAPH_COVERS = (
    mask((0, 1, 2)),
    mask((0, 1, 13)),
    mask((1, 11, 13)),
    mask((0, 12, 13)),
    mask((11, 12, 13)),
)


def _active_graph_automorphisms() -> tuple[dict[int, int], ...]:
    edge_set = {frozenset(edge) for edge in TARGET_MISSING_EDGES}
    answers = []
    for image_order in permutations(ACTIVE_VERTICES):
        action = dict(zip(ACTIVE_VERTICES, image_order, strict=True))
        if {
            frozenset((action[left], action[right]))
            for left, right in TARGET_MISSING_EDGES
        } == edge_set:
            answers.append(action)
    return tuple(answers)


def _family_orbit_key(
    family: Sequence[int],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Canonical key under full active graph auts, isolate ``S8``, and rows."""

    active_mask = mask(ACTIVE_VERTICES)
    keys = []
    for action in _active_graph_automorphisms():
        transformed = tuple(
            mask(action.get(vertex, vertex) for vertex in vertices(face))
            for face in family
        )
        for row_order in permutations(range(len(family))):
            active_rows = tuple(
                transformed[index] & active_mask for index in row_order
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
            keys.append((active_rows, isolate_columns))
    return min(keys)


def _forced_five_shadow(
    family: Sequence[int], fives: Sequence[int]
) -> tuple[int, ...]:
    return tuple(
        five
        for five in fives
        if any(five & four == four for four in family)
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


def _shadow_gate(
    family: Sequence[int],
    shadow: Sequence[int],
    base_g3: Sequence[int],
    base_g4: Sequence[int],
    total_fives: int,
) -> tuple[bool, bool, int, tuple[int, ...], tuple[int, ...]]:
    additional = total_fives - len(shadow)
    if additional < 0:
        return False, False, additional, (), ()
    lower: list[int] = []
    upper: list[int] = []
    strict = True
    for vertex in range(VERTEX_COUNT):
        a4_incidence = sum(bool(face & (1 << vertex)) for face in family)
        previous = base_g3[vertex] - a4_incidence
        forced_incidence = sum(bool(face & (1 << vertex)) for face in shadow)
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
    return strict, aggregate, additional, tuple(lower), tuple(upper)


def _small_base_transversals(hyperedges: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        candidate
        for order in range(1, 5)
        for chosen in combinations(range(VERTEX_COUNT), order)
        for candidate in (mask(chosen),)
        if all(candidate & edge for edge in hyperedges)
    )


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


def _hit_word(face: int, covers: Sequence[int]) -> int:
    return sum(bool(face & cover) << index for index, cover in enumerate(covers))


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
            "strict_link_sum_required": VERTEX_COUNT,
            "surviving_missing_four_counts": list(self.surviving_a),
            "excluded_missing_four_counts": list(self.excluded_a),
            "verified": self.verified,
        }


def _audit_schedules(schedules: Sequence[Mapping[str, object]]) -> ScheduleAudit:
    rows = tuple(
        (
            int(row["missing_four_count"]),
            int(row["missing_five_count_offset"]),
            int(row["f_vector"][-1]),
            summed_vertex_link_g4(tuple(map(int, row["f_vector"]))),
        )
        for row in schedules
    )
    surviving = tuple(a for a, _b, _f, total in rows if total >= VERTEX_COUNT)
    excluded = tuple(a for a, _b, _f, total in rows if total < VERTEX_COUNT)
    audit = ScheduleAudit(rows=rows, surviving_a=surviving, excluded_a=excluded)
    if not audit.verified:
        raise AssertionError("record-35 schedule audit changed")
    return audit


@dataclass(frozen=True, slots=True)
class AdditionCensus:
    representative_fours: tuple[int, ...]
    forced_shadow_count: int
    additional_five_count: int
    individually_cap_feasible_candidate_count: int
    raw_addition_family_count: int
    vertex_cap_feasible_count: int
    no_small_transversal_count: int
    no_small_transversal_facet_histogram: tuple[tuple[int, int], ...]
    scheduled_target_count: int
    scheduled_target_additions: tuple[tuple[int, ...], ...]

    def to_record(self) -> dict[str, object]:
        return {
            "representative_missing_fours": [
                list(vertices(face)) for face in self.representative_fours
            ],
            "forced_five_shadow_count": self.forced_shadow_count,
            "additional_minimal_five_count": self.additional_five_count,
            "individually_cap_feasible_A5_candidate_count": (
                self.individually_cap_feasible_candidate_count
            ),
            "raw_unordered_addition_family_count": self.raw_addition_family_count,
            "vertex_cap_feasible_addition_count": self.vertex_cap_feasible_count,
            "no_transversal_of_order_at_most_four_count": (
                self.no_small_transversal_count
            ),
            "facet_count_histogram_after_no_small_transversal": dict(
                self.no_small_transversal_facet_histogram
            ),
            "scheduled_target_count": self.scheduled_target_count,
            "scheduled_target_additions": [
                [list(vertices(face)) for face in additions]
                for additions in self.scheduled_target_additions
            ],
        }


def _enumerate_additions(
    family: Sequence[int], total_fives: int, scheduled_facets: int
) -> AdditionCensus:
    fives = independent_masks(TARGET_GRAPH6, 5)
    base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)
    shadow = _forced_five_shadow(family, fives)
    strict, aggregate, additional, lower, upper = _shadow_gate(
        family, shadow, base_g3, base_g4, total_fives
    )
    if not (strict and aggregate) or additional < 0:
        raise AssertionError("addition census called on an infeasible shadow")
    candidates = tuple(
        five
        for five in fives
        if five not in shadow
        and not any(five & four == four for four in family)
        and all(
            not (five & (1 << vertex)) or upper[vertex] >= 1
            for vertex in range(VERTEX_COUNT)
        )
    )
    small_transversals = _small_base_transversals(
        (*EDGE_MASKS, *family, *shadow)
    )
    five_covers = graph_vertex_cover_masks(TARGET_GRAPH6, 5)
    hit_words: dict[int, int] = {}
    fixed_word = (1 << len(five_covers)) - 1
    for face in (*family, *shadow):
        hit_words.setdefault(face, _hit_word(face, five_covers))
        fixed_word &= hit_words[face]

    cap_count = 0
    no_small_count = 0
    facet_histogram: Counter[int] = Counter()
    targets = []
    lower_at_least = tuple(
        mask(
            vertex
            for vertex in range(VERTEX_COUNT)
            if lower[vertex] >= threshold
        )
        for threshold in range(1, additional + 1)
    )
    upper_below = tuple(
        mask(
            vertex
            for vertex in range(VERTEX_COUNT)
            if upper[vertex] < threshold
        )
        for threshold in range(1, additional + 1)
    )
    for additions in combinations(candidates, additional):
        if additional == 1:
            at_least = additions
        elif additional == 2:
            first, second = additions
            at_least = (first | second, first & second)
        elif additional == 3:
            first, second, third = additions
            at_least = (
                first | second | third,
                (first & second) | (first & third) | (second & third),
                first & second & third,
            )
        else:
            raise AssertionError("record-35 addition count escaped 1..3")
        if any(
            required & present != required or forbidden & present
            for required, forbidden, present in zip(
                lower_at_least, upper_below, at_least, strict=True
            )
        ):
            continue
        cap_count += 1
        if any(
            all(transversal & face for face in additions)
            for transversal in small_transversals
        ):
            continue
        no_small_count += 1
        word = fixed_word
        for face in additions:
            hit_words.setdefault(face, _hit_word(face, five_covers))
            word &= hit_words[face]
        facet_count = word.bit_count()
        facet_histogram[facet_count] += 1
        if facet_count == scheduled_facets:
            targets.append(tuple(additions))
    return AdditionCensus(
        representative_fours=tuple(family),
        forced_shadow_count=len(shadow),
        additional_five_count=additional,
        individually_cap_feasible_candidate_count=len(candidates),
        raw_addition_family_count=comb(len(candidates), additional),
        vertex_cap_feasible_count=cap_count,
        no_small_transversal_count=no_small_count,
        no_small_transversal_facet_histogram=tuple(sorted(facet_histogram.items())),
        scheduled_target_count=len(targets),
        scheduled_target_additions=tuple(targets),
    )


@dataclass(frozen=True, slots=True)
class PositiveBranchAudit:
    individually_interval_feasible_a4_count: int
    a1_gate_histogram: tuple[tuple[tuple[bool, bool, int, int], int], ...]
    a1_viable_a4_count: int
    a1_orbit_size_histogram: tuple[tuple[int, int], ...]
    a1_completion_censuses: tuple[AdditionCensus, ...]
    a2_unordered_pair_count: int
    a2_family_interval_feasible_count: int
    a2_forced_shadow_histogram: tuple[tuple[int, int], ...]
    a2_viable_pair_count: int
    a2_orbit_size_histogram: tuple[tuple[int, int], ...]
    a2_completion_census: AdditionCensus

    @property
    def verified(self) -> bool:
        a1 = self.a1_completion_censuses
        a2 = self.a2_completion_census
        return (
            self.individually_interval_feasible_a4_count == 490
            and self.a1_gate_histogram
            == (
                ((False, False, 8, 3), 196),
                ((False, False, 9, 2), 112),
                ((False, False, 10, 1), 70),
                ((True, True, 8, 3), 56),
                ((True, True, 9, 2), 56),
            )
            and self.a1_viable_a4_count == 112
            and self.a1_orbit_size_histogram == ((56, 2),)
            and len(a1) == 2
            and a1[0].representative_fours == (mask(A1_REPRESENTATIVE_FOUR),)
            and a1[0].forced_shadow_count == 9
            and a1[0].additional_five_count == 2
            and a1[0].individually_cap_feasible_candidate_count == 252
            and a1[0].raw_addition_family_count == 31_626
            and a1[0].vertex_cap_feasible_count == 5_626
            and a1[0].no_small_transversal_count == 1
            and a1[0].scheduled_target_count == 1
            and a1[0].scheduled_target_additions
            == (tuple(mask(face) for face in A1_REPRESENTATIVE_FIVES),)
            and a1[1].representative_fours
            == (mask(A1_OTHER_REPRESENTATIVE_FOUR),)
            and a1[1].forced_shadow_count == 8
            and a1[1].additional_five_count == 3
            and a1[1].individually_cap_feasible_candidate_count == 462
            and a1[1].raw_addition_family_count == 16_328_620
            and a1[1].vertex_cap_feasible_count == 1_794_605
            and a1[1].no_small_transversal_count == 840
            and a1[1].scheduled_target_count == 0
            and self.a2_unordered_pair_count == 119_805
            and self.a2_family_interval_feasible_count == 111_699
            and self.a2_forced_shadow_histogram
            == (
                (15, 224),
                (16, 24_864),
                (17, 42_756),
                (18, 30_520),
                (19, 11_480),
                (20, 1_855),
            )
            and self.a2_viable_pair_count == 28
            and self.a2_orbit_size_histogram == ((28, 1),)
            and a2.representative_fours
            == tuple(mask(face) for face in A2_REPRESENTATIVE_FOURS)
            and a2.forced_shadow_count == 15
            and a2.additional_five_count == 2
            and a2.individually_cap_feasible_candidate_count == 56
            and a2.raw_addition_family_count == 1_540
            and a2.vertex_cap_feasible_count == 540
            and a2.no_small_transversal_count == 0
            and a2.scheduled_target_count == 0
        )

    def to_record(self) -> dict[str, object]:
        return {
            "individually_interval_feasible_A4_count": (
                self.individually_interval_feasible_a4_count
            ),
            "a1": {
                "shadow_gate_histogram": [
                    {
                        "strict": gate[0],
                        "aggregate": gate[1],
                        "forced_shadow_count": gate[2],
                        "additional_five_count": gate[3],
                        "multiplicity": count,
                    }
                    for gate, count in self.a1_gate_histogram
                ],
                "viable_A4_count": self.a1_viable_a4_count,
                "full_graph_x_S8_orbit_size_histogram": dict(
                    self.a1_orbit_size_histogram
                ),
                "representative_completion_censuses": [
                    audit.to_record() for audit in self.a1_completion_censuses
                ],
                "exact_residual_incidence_orbit_count": 1,
            },
            "a2": {
                "all_unordered_distinct_A4_pairs": self.a2_unordered_pair_count,
                "family_interval_feasible_pairs": (
                    self.a2_family_interval_feasible_count
                ),
                "forced_shadow_size_histogram": dict(
                    self.a2_forced_shadow_histogram
                ),
                "strict_aggregate_viable_pair_count": self.a2_viable_pair_count,
                "full_graph_x_S8_x_S2_orbit_size_histogram": dict(
                    self.a2_orbit_size_histogram
                ),
                "representative_completion_census": (
                    self.a2_completion_census.to_record()
                ),
                "branch_excluded": self.a2_completion_census.no_small_transversal_count
                == 0,
            },
            "verified": self.verified,
        }


def _audit_positive_branches() -> PositiveBranchAudit:
    local = analyze_local_link_graph(TARGET_GRAPH6)
    # This makes the singleton A4 filter lossless for the a=2 census: every
    # lower incidence bound is zero, while all remaining conditions are upper
    # bounds.  A singleton failure therefore cannot be repaired by adjoining
    # a second A4.
    if any(
        figure.missing_four_incidence_interval is None
        or figure.missing_four_incidence_interval[0] != 0
        for figure in (*local.vertex_figures, *local.edge_figures)
    ):
        raise AssertionError("record-35 local A4 lower-bound audit changed")
    fours = tuple(
        face
        for face in independent_masks(TARGET_GRAPH6, 4)
        if _a4_interval_feasible(face, local)
    )
    fives = independent_masks(TARGET_GRAPH6, 5)
    base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)

    a1_gates = []
    viable_a1 = []
    for four in fours:
        shadow = _forced_five_shadow((four,), fives)
        strict, aggregate, additional, _lower, _upper = _shadow_gate(
            (four,), shadow, base_g3, base_g4, 11
        )
        a1_gates.append((strict, aggregate, len(shadow), additional))
        if strict and aggregate:
            viable_a1.append(four)
    a1_orbits: defaultdict[object, list[int]] = defaultdict(list)
    for four in viable_a1:
        a1_orbits[_family_orbit_key((four,))].append(four)
    a1_representatives = tuple(sorted(min(orbit) for orbit in a1_orbits.values()))
    a1_censuses = tuple(
        _enumerate_additions((representative,), 11, 100)
        for representative in a1_representatives
    )

    a2_local = 0
    a2_histogram: Counter[int] = Counter()
    viable_a2 = []
    for family in combinations(fours, 2):
        if not _family_interval_feasible(family, local):
            continue
        a2_local += 1
        shadow = _forced_five_shadow(family, fives)
        strict, aggregate, _additional, _lower, _upper = _shadow_gate(
            family, shadow, base_g3, base_g4, 17
        )
        a2_histogram[len(shadow)] += 1
        if strict and aggregate:
            viable_a2.append(family)
    a2_orbits: defaultdict[object, list[tuple[int, int]]] = defaultdict(list)
    for family in viable_a2:
        a2_orbits[_family_orbit_key(family)].append(family)
    if len(a2_orbits) != 1:
        raise AssertionError("record-35 a=2 viable families changed orbit count")
    a2_representative = min(next(iter(a2_orbits.values())))
    audit = PositiveBranchAudit(
        individually_interval_feasible_a4_count=len(fours),
        a1_gate_histogram=tuple(sorted(Counter(a1_gates).items())),
        a1_viable_a4_count=len(viable_a1),
        a1_orbit_size_histogram=tuple(
            sorted(Counter(map(len, a1_orbits.values())).items())
        ),
        a1_completion_censuses=a1_censuses,
        a2_unordered_pair_count=comb(len(fours), 2),
        a2_family_interval_feasible_count=a2_local,
        a2_forced_shadow_histogram=tuple(sorted(a2_histogram.items())),
        a2_viable_pair_count=len(viable_a2),
        a2_orbit_size_histogram=tuple(
            sorted(Counter(map(len, a2_orbits.values())).items())
        ),
        a2_completion_census=_enumerate_additions(a2_representative, 17, 98),
    )
    if not audit.verified:
        raise AssertionError("record-35 positive-branch census changed")
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
        raise AssertionError("record-35 isolate-column count changed")
    faces = list(active_rows)
    for isolate, column in zip(ISOLATED_VERTICES, columns, strict=True):
        for row in range(ROW_COUNT):
            if column & (1 << row):
                faces[row] |= 1 << isolate
    return tuple(sorted(faces))


def _has_small_transversal_from_graph_covers(
    fives: Sequence[int], graph_covers: Sequence[int]
) -> bool:
    """Exact order-at-most-four transversal test after the edge clutter.

    Every minimal cover of the five missing edges has order three.  A cover
    extends within the size-four budget iff it already hits every A5 or one
    vertex lies in every A5 row that it misses.
    """

    for cover in graph_covers:
        missed = tuple(face for face in fives if not face & cover)
        if not missed:
            return True
        common = missed[0]
        for face in missed[1:]:
            common &= face
        if common:
            return True
    return False


@lru_cache(maxsize=None)
def _placement_options(support: tuple[int, ...]) -> tuple[tuple[tuple[Pair, int], ...], ...]:
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


def _restriction_finite(
    restriction: Sequence[int], labels: Mapping[Pair, int]
) -> bool:
    return _finite_code(
        len(restriction),
        tuple(
            labels[pair(left, right)]
            for left, right in combinations(restriction, 2)
        ),
    )


LabelState = tuple[tuple[Pair, int], ...]


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
    restriction_size_histogram: tuple[tuple[int, int], ...]
    raw_extension_count: int
    compatible_extension_count: int
    passing_state_count: int
    passing_state_sha256: str

    def to_record(self) -> dict[str, object]:
        return {
            "support_index": self.support_index,
            "support": list(self.support),
            "complete_lanner_placement_count": len(
                _placement_options(self.support)
            ),
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
            "stages": [stage.to_record() for stage in self.stages],
            "unrestricted_compact_coxeter_no_go": self.compact_no_go,
            "monotonicity_scope": (
                "all unassigned meeting pairs are initially label 2; any "
                "label 3, 4, or 5 increases the nonnegative cosine adjacency "
                "entrywise and cannot repair a non-elliptic principal vertex"
            ),
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
    actual_vertices = _rank_nine_vertices((*TARGET_MISSING_EDGES, *missing_fives))
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
                    for vertex in actual_vertices
                },
                key=lambda row: (len(row), row),
            )
        )
        if any(
            set(edge) <= set(restriction)
            for edge in TARGET_MISSING_EDGES
            for restriction in restrictions
        ):
            raise AssertionError("an actual vertex contains a missing edge")
        options = _placement_options(support)
        passing: set[LabelState] = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in options:
                if not all(
                    old.get(edge, label) == label for edge, label in option
                ):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(
                    _restriction_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            LabelStage(
                support_index=support_index,
                support=support,
                restriction_size_histogram=tuple(
                    sorted(Counter(map(len, restrictions)).items())
                ),
                raw_extension_count=(
                    (stages[-1].passing_state_count if stages else 1)
                    * len(options)
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
        rank_nine_vertex_count=len(actual_vertices),
        stages=tuple(stages),
    )
    if not audit.compact_no_go:
        raise AssertionError(f"record-35 a=0 family {family_index} survived labels")
    return audit


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
            "missing_fives": [list(vertices(face)) for face in self.missing_fives],
            "minimal_transversal_size_histogram": dict(
                self.minimal_transversal_histogram
            ),
            "facet_complement_sha256": self.facet_complement_sha256,
            "topology_valid": self.topology_valid,
            "global_selection_audit": self.global_audit,
        }


@dataclass(frozen=True, slots=True)
class A0Audit:
    independent_active_masks: tuple[int, ...]
    independent_active_mask_count: int
    active_row_configuration_count: int
    column_distribution_problem_count: int
    total_distribution_path_count: int
    repeated_row_path_count: int
    small_transversal_path_count: int
    blocker_valid_path_count: int
    facet_count_histogram: tuple[tuple[int, int], ...]
    path_transcript_sha256: str
    target_path_count: int
    exact_target_family_count: int
    target_family_sha256: str
    family_audits: tuple[A0FamilyAudit, ...]
    topology_valid_family_indices: tuple[int, ...]
    all_target_full_graph_orbit_count: int
    topology_valid_full_graph_orbit_count: int
    topology_orbit_key_sha256: str
    label_audits: tuple[FamilyLabelAudit, ...]
    label_transcript_sha256: str

    @property
    def verified(self) -> bool:
        digest_checks = (
            EXPECTED_A0_PATH_TRANSCRIPT_SHA256 in (
                "TO_BE_FROZEN",
                self.path_transcript_sha256,
            )
            and EXPECTED_A0_TARGET_FAMILY_SHA256
            in ("TO_BE_FROZEN", self.target_family_sha256)
            and EXPECTED_A0_LABEL_TRANSCRIPT_SHA256
            in ("TO_BE_FROZEN", self.label_transcript_sha256)
        )
        return (
            self.independent_active_masks == EXPECTED_INDEPENDENT_ACTIVE_MASKS
            and self.independent_active_mask_count == 18
            and self.active_row_configuration_count == 1_221
            and self.column_distribution_problem_count == 235
            and self.total_distribution_path_count == 375_582
            and self.repeated_row_path_count == 7_239
            and self.small_transversal_path_count == 0
            and self.blocker_valid_path_count == 368_343
            and sum(count for _facets, count in self.facet_count_histogram)
            == 368_343
            and self.target_path_count == 20
            and self.exact_target_family_count == 18
            and len(self.family_audits) == 18
            and all(
                audit.minimal_transversal_histogram == ((5, 102),)
                for audit in self.family_audits
            )
            and self.topology_valid_family_indices == (0, 1, 2, 3, 4)
            and self.all_target_full_graph_orbit_count == 11
            and self.topology_valid_full_graph_orbit_count == 3
            and len(self.label_audits) == 18
            and all(audit.compact_no_go for audit in self.label_audits)
            and digest_checks
        )

    def to_record(self) -> dict[str, object]:
        return {
            "lossless_enumeration_scope": {
                "active_vertices": list(ACTIVE_VERTICES),
                "isolated_vertices": list(ISOLATED_VERTICES),
                "independent_active_masks": [
                    list(vertices(value)) for value in self.independent_active_masks
                ],
                "independent_active_mask_size_histogram": dict(
                    sorted(
                        Counter(
                            value.bit_count()
                            for value in self.independent_active_masks
                        ).items()
                    )
                ),
                "independent_active_mask_count": self.independent_active_mask_count,
                "active_row_configuration_count": (
                    self.active_row_configuration_count
                ),
                "active_row_statement": (
                    "all independent active subsets allowed by the exact "
                    "incidence caps are admitted in every A5 row; rows are "
                    "an unordered multiset. No singleton-active-row "
                    "normalization is imposed"
                ),
                "isolate_column_statement": (
                    "the eight isolate columns are an exact multiplicity "
                    "vector, hence a complete quotient by S8"
                ),
                "column_distribution_problem_count": (
                    self.column_distribution_problem_count
                ),
                "total_distribution_path_count": (
                    self.total_distribution_path_count
                ),
                "repeated_row_path_count": self.repeated_row_path_count,
                "small_transversal_path_count": (
                    self.small_transversal_path_count
                ),
                "blocker_valid_path_count": self.blocker_valid_path_count,
                "path_transcript_sha256": self.path_transcript_sha256,
            },
            "minimal_graph_covers": [
                list(vertices(value)) for value in EXPECTED_MINIMAL_GRAPH_COVERS
            ],
            "facet_count_histogram_on_blocker_valid_paths": dict(
                self.facet_count_histogram
            ),
            "scheduled_facet_count": 102,
            "scheduled_target_path_count": self.target_path_count,
            "exact_unordered_target_family_count": self.exact_target_family_count,
            "target_family_sha256": self.target_family_sha256,
            "target_family_audits": [
                audit.to_record() for audit in self.family_audits
            ],
            "topology_valid_family_indices": list(
                self.topology_valid_family_indices
            ),
            "full_active_graph_automorphism_count": len(
                _active_graph_automorphisms()
            ),
            "all_target_orbit_count_under_full_graph_x_S8_x_S5": (
                self.all_target_full_graph_orbit_count
            ),
            "topology_valid_orbit_count_under_full_graph_x_S8_x_S5": (
                self.topology_valid_full_graph_orbit_count
            ),
            "topology_orbit_key_sha256": self.topology_orbit_key_sha256,
            "unrestricted_label_audits_for_all_18_targets": [
                audit.to_record() for audit in self.label_audits
            ],
            "label_transcript_sha256": self.label_transcript_sha256,
            "verified": self.verified,
        }


def _audit_a0(
    support_core_masks: Sequence[int], schedules: Sequence[Mapping[str, object]]
) -> A0Audit:
    _column_distributions.cache_clear()
    graph_covers = _minimal_graph_covers()
    if graph_covers != EXPECTED_MINIMAL_GRAPH_COVERS:
        raise AssertionError("record-35 minimal graph covers changed")
    five_covers = graph_vertex_cover_masks(TARGET_GRAPH6, 5)
    all_five_word = (1 << len(five_covers)) - 1
    five_hit_words: dict[int, int] = {}
    active_edge_masks = tuple(
        edge
        for edge in EDGE_MASKS
        if edge & mask(ACTIVE_VERTICES) == edge
    )
    _base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)
    incidence_caps = tuple(value - 1 for value in base_g4)
    independent_active = tuple(
        value
        for order in range(len(ACTIVE_VERTICES) + 1)
        for chosen in combinations(ACTIVE_VERTICES, order)
        for value in (mask(chosen),)
        if all(value & edge != edge for edge in active_edge_masks)
        and all(
            int(bool(value & (1 << vertex))) <= incidence_caps[vertex]
            for vertex in ACTIVE_VERTICES
        )
    )

    active_configurations = []
    for configuration in combinations_with_replacement(
        independent_active, ROW_COUNT
    ):
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
        # An active vertex in every row missed by a graph cover would extend
        # that order-three cover to a forbidden transversal of order four.
        if any(
            all(
                not (pattern & (1 << row_index))
                or configuration[row_index] & (1 << vertex)
                for row_index in range(ROW_COUNT)
            )
            for pattern in missed_rows
            for vertex in ACTIVE_VERTICES
        ):
            continue
        active_configurations.append((configuration, missed_rows))

    total_paths = 0
    repeated = 0
    small_transversal_failures = 0
    valid_paths = 0
    histogram: Counter[int] = Counter()
    target_paths: list[tuple[int, ...]] = []
    path_hasher = sha256()
    for configuration, missed_rows in active_configurations:
        allowed = tuple(
            column
            for column in range(1 << ROW_COUNT)
            if column.bit_count() <= ROW_COUNT - 1
            and all(column & rows != rows for rows in missed_rows)
        )
        row_sums = tuple(5 - row.bit_count() for row in configuration)
        for multiplicities in _column_distributions(allowed, row_sums):
            total_paths += 1
            family = _faces_from_distribution(
                configuration, allowed, multiplicities
            )
            path_hasher.update(
                b"".join(value.to_bytes(2, "little") for value in family)
            )
            if len(set(family)) != ROW_COUNT:
                repeated += 1
                path_hasher.update(b"R")
                continue
            if _has_small_transversal_from_graph_covers(family, graph_covers):
                small_transversal_failures += 1
                path_hasher.update(b"S")
                continue
            valid_paths += 1
            word = all_five_word
            for face in family:
                if face not in five_hit_words:
                    five_hit_words[face] = _hit_word(face, five_covers)
                word &= five_hit_words[face]
            facet_count = word.bit_count()
            path_hasher.update(b"V" + facet_count.to_bytes(2, "little"))
            histogram[facet_count] += 1
            if facet_count == 102:
                target_paths.append(family)

    target_families = tuple(sorted(set(target_paths)))
    family_audits = []
    label_audits = []
    topology_indices = []
    all_orbit_keys = set()
    topology_orbit_keys = set()
    for family_index, family in enumerate(target_families):
        transversals = _minimal_transversals((*EDGE_MASKS, *family))
        transversal_histogram = tuple(
            sorted(Counter(value.bit_count() for value in transversals).items())
        )
        five_transversals = tuple(
            value for value in transversals if value.bit_count() == 5
        )
        global_audit = audit_global_facet_selection(
            TARGET_GRAPH6, support_core_masks, schedules, five_transversals
        )
        family_audit = A0FamilyAudit(
            family_index=family_index,
            missing_fives=family,
            minimal_transversal_histogram=transversal_histogram,
            facet_complement_sha256=sha256(
                b"".join(
                    value.to_bytes(2, "little") for value in five_transversals
                )
            ).hexdigest(),
            packed_global_audit=_pack_global_audit(global_audit.to_record()),
        )
        family_audits.append(family_audit)
        orbit_key = _family_orbit_key(family)
        all_orbit_keys.add(orbit_key)
        if family_audit.topology_valid:
            topology_indices.append(family_index)
            topology_orbit_keys.add(orbit_key)
        label_audits.append(_audit_family_labels(family_index, family))

    audit = A0Audit(
        independent_active_masks=independent_active,
        independent_active_mask_count=len(independent_active),
        active_row_configuration_count=len(active_configurations),
        column_distribution_problem_count=_column_distributions.cache_info().currsize,
        total_distribution_path_count=total_paths,
        repeated_row_path_count=repeated,
        small_transversal_path_count=small_transversal_failures,
        blocker_valid_path_count=valid_paths,
        facet_count_histogram=tuple(sorted(histogram.items())),
        path_transcript_sha256=path_hasher.hexdigest(),
        target_path_count=len(target_paths),
        exact_target_family_count=len(target_families),
        target_family_sha256=_canonical_digest(
            [
                [list(vertices(face)) for face in family]
                for family in target_families
            ]
        ),
        family_audits=tuple(family_audits),
        topology_valid_family_indices=tuple(topology_indices),
        all_target_full_graph_orbit_count=len(all_orbit_keys),
        topology_valid_full_graph_orbit_count=len(topology_orbit_keys),
        topology_orbit_key_sha256=_canonical_digest(
            tuple(sorted(topology_orbit_keys))
        ),
        label_audits=tuple(label_audits),
        label_transcript_sha256=_canonical_digest(
            [label.to_record() for label in label_audits]
        ),
    )
    if not audit.verified:
        raise AssertionError("record-35 a=0 exhaustive census changed")
    return audit


@dataclass(frozen=True, slots=True)
class A1LabelAudit:
    rank_nine_vertex_count: int
    order_four_placement_count: int
    overlapping_order_five_placement_count: int
    disjoint_order_five_placement_count: int
    compatible_overlapping_join_count: int
    left_cluster_passing_state_count: int
    block_restriction_size_histogram: tuple[tuple[int, int], ...]
    sparse_block_state_count: int
    sparse_block_failure_count: int
    cross_pair_count: int
    label_three_test_count: int
    label_three_passing_count: int
    passing_edges_per_sparse_state_histogram: tuple[tuple[int, int], ...]
    first_bad_restriction_size_histogram: tuple[tuple[int, int], ...]
    stronger_single_edge_survivors: tuple[tuple[int, int], ...]
    passing_R_endpoint_histogram: tuple[tuple[int, int], ...]
    exceptional_class_histogram: tuple[tuple[str, int], ...]
    exact_index_two_witnesses: tuple[ExactIndexTwoWitness, ...]
    transcript_sha256: str

    @property
    def unrestricted_compact_no_go(self) -> bool:
        digest_ok = EXPECTED_A1_LABEL_TRANSCRIPT_SHA256 in (
            "TO_BE_FROZEN",
            self.transcript_sha256,
        )
        return (
            self.rank_nine_vertex_count == 100
            and self.order_four_placement_count == 108
            and self.overlapping_order_five_placement_count == 420
            and self.disjoint_order_five_placement_count == 420
            and self.compatible_overlapping_join_count == 1_542
            and self.left_cluster_passing_state_count == 84
            and self.block_restriction_size_histogram
            == ((7, 20), (8, 55), (9, 15))
            and self.sparse_block_state_count == 35_280
            and self.sparse_block_failure_count == 0
            and self.cross_pair_count == 30
            and self.label_three_test_count == 1_058_400
            and self.label_three_passing_count == 5_040
            and self.passing_edges_per_sparse_state_histogram == ((1, 5_040),)
            and self.first_bad_restriction_size_histogram
            == ((7, 705_600), (8, 346_320), (9, 1_440))
            and self.stronger_single_edge_survivors == ((4, 0), (5, 0))
            and self.passing_R_endpoint_histogram == ((11, 2_520), (12, 2_520))
            and self.exceptional_class_histogram
            == (("A5_A5_type_1_1", 4_320), ("A5_A5_type_4_1", 720))
            and tuple(witness.inertia for witness in self.exact_index_two_witnesses)
            == ((8, 2, 0), (8, 2, 0))
            and digest_ok
        )

    def to_record(self) -> dict[str, object]:
        return {
            "residual_minimal_nonfaces": {
                "missing_four": list(A1_Q_SUPPORT),
                "overlapping_missing_five": list(A1_R_SUPPORT),
                "disjoint_missing_five": list(A1_B_SUPPORT),
            },
            "rank_nine_vertex_count": self.rank_nine_vertex_count,
            "complete_placement_counts": {
                "order_four": self.order_four_placement_count,
                "overlapping_order_five": (
                    self.overlapping_order_five_placement_count
                ),
                "disjoint_order_five": self.disjoint_order_five_placement_count,
            },
            "overlapping_cluster": {
                "compatible_join_count": self.compatible_overlapping_join_count,
                "actual_vertex_passing_state_count": (
                    self.left_cluster_passing_state_count
                ),
                "unforced_internal_pairs": [
                    list(edge) for edge in A1_INTERNAL_UNFORCED_PAIRS
                ],
            },
            "block_actual_vertex_restriction_size_histogram": dict(
                self.block_restriction_size_histogram
            ),
            "sparse_block_state_count": self.sparse_block_state_count,
            "sparse_block_failure_count": self.sparse_block_failure_count,
            "cross_pair_count": self.cross_pair_count,
            "label_three_test_count": self.label_three_test_count,
            "label_three_passing_count": self.label_three_passing_count,
            "passing_edges_per_sparse_state_histogram": dict(
                self.passing_edges_per_sparse_state_histogram
            ),
            "first_bad_restriction_size_histogram": dict(
                self.first_bad_restriction_size_histogram
            ),
            "stronger_single_edge_survivors": dict(
                self.stronger_single_edge_survivors
            ),
            "passing_R_endpoint_histogram": dict(
                self.passing_R_endpoint_histogram
            ),
            "exceptional_A5_A5_class_histogram": dict(
                self.exceptional_class_histogram
            ),
            "exact_index_two_witnesses": [
                witness.to_record() for witness in self.exact_index_two_witnesses
            ],
            "coverage_argument": (
                "Every full locally elliptic state can be weakened by "
                "changing all but one cross edge to label 2 and the retained "
                "edge to label 3. Exactly 5,040 sparse states admit one such "
                "edge and no sparse state admits two. Labels 4 and 5 admit "
                "none. Every exceptional edge joins the R and B A5 supports, "
                "whose exact ten-node principal Gram has negative index two. "
                "Labels on the two internal unforced pairs involve facet 2 "
                "outside that principal block and cannot change its inertia"
            ),
            "transcript_sha256": self.transcript_sha256,
            "unrestricted_compact_coxeter_no_go": (
                self.unrestricted_compact_no_go
            ),
        }


def _merge_label_states(
    left: Mapping[Pair, int], right: Mapping[Pair, int]
) -> dict[Pair, int] | None:
    answer = dict(left)
    for edge, label in right.items():
        previous = answer.get(edge)
        if previous is not None and previous != label:
            return None
        answer[edge] = label
    return answer


def _all_placement_maps(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    return tuple(dict(option) for option in _placement_options(tuple(support)))


def _rooted_lanner_code(
    support: Sequence[int], state: Mapping[Pair, int], root: int
) -> tuple[int, ...]:
    """Exact rooted weighted-graph code, exhausting the remaining ``4!``."""

    support = tuple(sorted(support))
    if root not in support or len(support) != 5:
        raise ValueError("rooted Lanner code needs a root in an A5 support")
    others = tuple(vertex for vertex in support if vertex != root)
    return min(
        tuple(
            state[pair(order[left], order[right])]
            for left, right in combinations(range(5), 2)
        )
        for tail in permutations(others)
        for order in ((root, *tail),)
    )


def _one_bridge_class_key(
    left_support: Sequence[int],
    right_support: Sequence[int],
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    bridge: Pair,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    left_endpoint = next(vertex for vertex in bridge if vertex in left_support)
    right_endpoint = next(vertex for vertex in bridge if vertex in right_support)
    return tuple(
        sorted(
            (
                _rooted_lanner_code(left_support, left_state, left_endpoint),
                _rooted_lanner_code(right_support, right_state, right_endpoint),
            )
        )
    )


def _cluster_restrictions(
    actual_vertices: Sequence[Sequence[int]], cluster: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    support = set(cluster)
    return tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in support)
                for vertex in actual_vertices
            },
            key=lambda row: (len(row), row),
        )
    )


def _audit_a1_labels() -> A1LabelAudit:
    actual_vertices = _rank_nine_vertices(
        (
            *TARGET_MISSING_EDGES,
            A1_Q_SUPPORT,
            A1_R_SUPPORT,
            A1_B_SUPPORT,
        )
    )
    q_states = _all_placement_maps(A1_Q_SUPPORT)
    r_states = _all_placement_maps(A1_R_SUPPORT)
    b_states = _all_placement_maps(A1_B_SUPPORT)
    left_restrictions = _cluster_restrictions(actual_vertices, A1_LEFT_UNION)
    compatible = 0
    left_passing: set[LabelState] = set()
    for q_state in q_states:
        for r_state in r_states:
            joined = _merge_label_states(q_state, r_state)
            if joined is None:
                continue
            compatible += 1
            labels = {
                pair(left, right): 2
                for left, right in combinations(A1_LEFT_UNION, 2)
            }
            labels.update(joined)
            if all(
                _restriction_finite(restriction, labels)
                for restriction in left_restrictions
            ):
                left_passing.add(tuple(sorted(joined.items())))
    left_states = tuple(dict(state) for state in sorted(left_passing))

    block_restrictions = _cluster_restrictions(actual_vertices, A1_BLOCK_UNION)
    restrictions_by_cross = {
        edge: tuple(
            restriction
            for restriction in block_restrictions
            if set(edge) <= set(restriction)
        )
        for edge in A1_CROSS_PAIRS
    }
    if any(not rows for rows in restrictions_by_cross.values()):
        raise AssertionError("an a=1 cross pair occurs in no actual vertex")
    template = {
        pair(left, right): 2
        for left, right in combinations(A1_BLOCK_UNION, 2)
    }
    sparse_failures = 0
    passing: list[tuple[int, int, Pair]] = []
    passing_by_state: defaultdict[tuple[int, int], list[Pair]] = defaultdict(list)
    bad_sizes: Counter[int] = Counter()
    for left_index, left_state in enumerate(left_states):
        for right_index, right_state in enumerate(b_states):
            base = dict(template)
            base.update(left_state)
            base.update(right_state)
            sparse_failures += not all(
                _restriction_finite(restriction, base)
                for restriction in block_restrictions
            )
            for edge in A1_CROSS_PAIRS:
                labels = dict(base)
                labels[edge] = 3
                bad = next(
                    (
                        restriction
                        for restriction in restrictions_by_cross[edge]
                        if not _restriction_finite(restriction, labels)
                    ),
                    None,
                )
                if bad is None:
                    passing.append((left_index, right_index, edge))
                    passing_by_state[(left_index, right_index)].append(edge)
                else:
                    bad_sizes[len(bad)] += 1

    stronger = []
    for label in (4, 5):
        count = 0
        for left_index, right_index, edge in passing:
            labels = dict(template)
            labels.update(left_states[left_index])
            labels.update(b_states[right_index])
            labels[edge] = label
            count += all(
                _restriction_finite(restriction, labels)
                for restriction in restrictions_by_cross[edge]
            )
        stronger.append((label, count))

    canonical_r = _canonical_lanner_states(A1_R_SUPPORT)
    canonical_b = _canonical_lanner_states(A1_B_SUPPORT)
    reference_keys = (
        _one_bridge_class_key(
            A1_R_SUPPORT,
            A1_B_SUPPORT,
            canonical_r[1],
            canonical_b[1],
            pair(A1_R_SUPPORT[0], A1_B_SUPPORT[0]),
        ),
        _one_bridge_class_key(
            A1_R_SUPPORT,
            A1_B_SUPPORT,
            canonical_r[4],
            canonical_b[1],
            pair(A1_R_SUPPORT[4], A1_B_SUPPORT[0]),
        ),
    )
    witnesses = _a5_witnesses(A1_R_SUPPORT, A1_B_SUPPORT)
    classes: Counter[str] = Counter()
    endpoint_histogram: Counter[int] = Counter()
    transcript = []
    r_set = set(A1_R_SUPPORT)
    b_set = set(A1_B_SUPPORT)
    for left_index, right_index, edge in passing:
        r_endpoints = tuple(vertex for vertex in edge if vertex in r_set)
        b_endpoints = tuple(vertex for vertex in edge if vertex in b_set)
        if len(r_endpoints) != 1 or len(b_endpoints) != 1:
            raise AssertionError(
                "a locally possible cross edge does not join the R and B Lanners"
            )
        class_key = _one_bridge_class_key(
            A1_R_SUPPORT,
            A1_B_SUPPORT,
            left_states[left_index],
            b_states[right_index],
            edge,
        )
        matches = tuple(
            index
            for index, reference in enumerate(reference_keys)
            if class_key == reference
        )
        if len(matches) != 1:
            raise AssertionError("an exceptional a=1 block has no unique class")
        name = witnesses[matches[0]].name
        classes[name] += 1
        endpoint_histogram[r_endpoints[0]] += 1
        transcript.append((left_index, right_index, edge, name))

    # ``_a5_witnesses`` supplies the exact rooted Gram calculations, but its
    # multiplicities belong to the older record-44 support gate.  Rewrap the
    # same exact witnesses with the multiplicities of this complete
    # record-35 exceptional-state census so the serialized provenance cannot
    # silently report the wrong counts.
    witnesses = tuple(
        ExactIndexTwoWitness(
            name=witness.name,
            order=witness.order,
            leading_minors=witness.leading_minors,
            signs=witness.signs,
            inertia=witness.inertia,
            multiplicity_per_support_gate=classes[witness.name],
        )
        for witness in witnesses
    )

    audit = A1LabelAudit(
        rank_nine_vertex_count=len(actual_vertices),
        order_four_placement_count=len(q_states),
        overlapping_order_five_placement_count=len(r_states),
        disjoint_order_five_placement_count=len(b_states),
        compatible_overlapping_join_count=compatible,
        left_cluster_passing_state_count=len(left_states),
        block_restriction_size_histogram=tuple(
            sorted(Counter(map(len, block_restrictions)).items())
        ),
        sparse_block_state_count=len(left_states) * len(b_states),
        sparse_block_failure_count=sparse_failures,
        cross_pair_count=len(A1_CROSS_PAIRS),
        label_three_test_count=(
            len(left_states) * len(b_states) * len(A1_CROSS_PAIRS)
        ),
        label_three_passing_count=len(passing),
        passing_edges_per_sparse_state_histogram=tuple(
            sorted(Counter(map(len, passing_by_state.values())).items())
        ),
        first_bad_restriction_size_histogram=tuple(sorted(bad_sizes.items())),
        stronger_single_edge_survivors=tuple(stronger),
        passing_R_endpoint_histogram=tuple(sorted(endpoint_histogram.items())),
        exceptional_class_histogram=tuple(sorted(classes.items())),
        exact_index_two_witnesses=witnesses,
        transcript_sha256=_canonical_digest(transcript),
    )
    if not audit.unrestricted_compact_no_go:
        raise AssertionError("record-35 residual a=1 label gate changed")
    return audit


def _audit_a1_incidence(
    support_core_masks: Sequence[int], schedules: Sequence[Mapping[str, object]]
) -> tuple[tuple[str, object], ...]:
    clutter = (
        *EDGE_MASKS,
        mask(A1_Q_SUPPORT),
        *(mask(face) for face in A1_REPRESENTATIVE_FIVES),
    )
    transversals = _minimal_transversals(clutter)
    histogram = tuple(sorted(Counter(value.bit_count() for value in transversals).items()))
    if histogram != ((5, 100),):
        raise AssertionError("record-35 a=1 residual blocker changed")
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6, support_core_masks, schedules, transversals
    )
    packed = _pack_global_audit(global_audit.to_record())
    record = dict(packed)
    if not (
        record["facet_count"] == 100
        and record["ridge_degree_histogram"] == {"2": 450}
        and record["missing_four_count"] == 1
        and record["passes_connected_homology_audit"]
    ):
        raise AssertionError("record-35 a=1 residual topology changed")
    return (
        ("minimal_transversal_size_histogram", {"5": 100}),
        (
            "facet_complement_sha256",
            sha256(
                b"".join(value.to_bytes(2, "little") for value in transversals)
            ).hexdigest(),
        ),
        *packed,
    )


@dataclass(frozen=True, slots=True)
class Record35AllCompletionsAudit:
    decoded_edges: tuple[tuple[int, int], ...]
    schedule: ScheduleAudit
    positive_branches: PositiveBranchAudit
    a0: A0Audit
    a1_incidence_audit: tuple[tuple[str, object], ...]
    a1_label: A1LabelAudit

    @property
    def verified(self) -> bool:
        return (
            self.decoded_edges == TARGET_MISSING_EDGES
            and self.schedule.verified
            and self.positive_branches.verified
            and self.a0.verified
            and bool(dict(self.a1_incidence_audit)["passes_connected_homology_audit"])
            and self.a1_label.unrestricted_compact_no_go
        )

    @property
    def exact_residual_incidence_orbit_count(self) -> int:
        # Three topology-valid a=0 orbits and one a=1 orbit.
        return (
            self.a0.topology_valid_full_graph_orbit_count
            + int(self.positive_branches.a1_completion_censuses[0].scheduled_target_count)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "record_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
            "decoded_edges": [list(edge) for edge in self.decoded_edges],
            "schedule_audit": self.schedule.to_record(),
            "positive_branch_audit": self.positive_branches.to_record(),
            "a0_audit": self.a0.to_record(),
            "a1_residual_incidence_audit": dict(self.a1_incidence_audit),
            "a1_unrestricted_label_audit": self.a1_label.to_record(),
            "exact_topology_valid_residual_incidence_orbit_count": (
                self.exact_residual_incidence_orbit_count
            ),
            "whole_graph_shadow_compact_coxeter_excluded": self.verified,
            "finite_volume_noncompact_chamber_excluded": False,
            "verified": self.verified,
        }


@lru_cache(maxsize=1)
def audit_record35_all_completions(
    schedules_json: str, support_core_json: str
) -> Record35AllCompletionsAudit:
    schedules = json.loads(schedules_json)
    support_core = json.loads(support_core_json)
    if not isinstance(schedules, list) or not isinstance(support_core, list):
        raise ValueError("record-35 archived schedules/support must be JSON lists")
    order, decoded_edges = decode_graph6(TARGET_GRAPH6)
    if order != VERTEX_COUNT:
        raise AssertionError("record-35 graph order changed")
    audit = Record35AllCompletionsAudit(
        decoded_edges=decoded_edges,
        schedule=_audit_schedules(schedules),
        positive_branches=_audit_positive_branches(),
        a0=_audit_a0(tuple(map(int, support_core)), schedules),
        a1_incidence_audit=_audit_a1_incidence(
            tuple(map(int, support_core)), schedules
        ),
        a1_label=_audit_a1_labels(),
    )
    if not audit.verified:
        raise AssertionError("record-35 whole-shadow audit changed")
    return audit


__all__ = [
    "A0Audit",
    "A1LabelAudit",
    "EXPECTED_INDEPENDENT_ACTIVE_MASKS",
    "EXPECTED_SCHEDULES",
    "PositiveBranchAudit",
    "Record35AllCompletionsAudit",
    "audit_record35_all_completions",
]
