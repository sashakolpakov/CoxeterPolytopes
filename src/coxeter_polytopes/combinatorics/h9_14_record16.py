"""Analytic ridge obstruction for the ``a=0`` branch of H9/14 record 16.

The missing-edge graph is ``K_4-e`` on ``{0,1,12,13}`` together with ten
isolates.  In the ``a=0`` Dehn--Sommerville schedule the remaining minimal
nonfaces are eight five-sets on ``W={0,...,11}``, and there must be 110
facets.  Blocker duality splits their five-element complements into two
families.  The ridge equations then reduce those two families to a complete
bipartite graph and a small two-graph.  A Cauchy--Schwarz bound gives the
contradiction without SAT, SMT, or a missing-five census.

This module checks the finite incidence reduction and all numerical cases.
The accompanying report supplies the short mathematical proof of the two
structural lemmas used by the arithmetic replay.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 16
TARGET_GRAPH6 = "M???????????o?o?_"
TARGET_VERTEX_COUNT = 14
TARGET_A0_FACET_COUNT = 110
TARGET_A0_MISSING_FIVE_COUNT = 8

ACTIVE_VERTICES = (0, 1)
ISOLATED_VERTICES = tuple(range(2, 12))
WORK_VERTICES = ACTIVE_VERTICES + ISOLATED_VERTICES
SPECIAL_VERTICES = (12, 13)
TARGET_MISSING_EDGES = (
    (0, 12),
    (1, 12),
    (0, 13),
    (1, 13),
    (12, 13),
)
EXPECTED_MINIMAL_GRAPH_COVERS = (
    (0, 1, 12),
    (0, 1, 13),
    (12, 13),
)


def mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(
        vertex
        for vertex in range(TARGET_VERTEX_COUNT)
        if value & (1 << vertex)
    )


EDGE_MASKS = tuple(mask(edge) for edge in TARGET_MISSING_EDGES)


def _is_graph_cover(candidate: int) -> bool:
    return all(candidate & edge for edge in EDGE_MASKS)


def minimal_graph_covers() -> tuple[int, ...]:
    """Enumerate the inclusion-minimal covers of the five-edge shadow."""

    result = []
    for candidate in range(1 << TARGET_VERTEX_COUNT):
        if not _is_graph_cover(candidate):
            continue
        if all(
            not _is_graph_cover(candidate ^ (1 << vertex))
            for vertex in vertices(candidate)
        ):
            result.append(candidate)
    return tuple(result)


def _base_vertex_link_g3_g4() -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Compute link ``g3,g4`` before adding four- or five-nonfaces."""

    independent: dict[int, tuple[int, ...]] = {}
    for size in range(2, 6):
        values = []
        for chosen in combinations(range(TARGET_VERTEX_COUNT), size):
            candidate = mask(chosen)
            if not any(candidate & edge == edge for edge in EDGE_MASKS):
                values.append(candidate)
        independent[size] = tuple(values)
    g3_values = []
    g4_values = []
    for vertex in range(TARGET_VERTEX_COUNT):
        counts = {
            size: sum(
                bool(face & (1 << vertex)) for face in independent[size]
            )
            for size in independent
        }
        f0, f1, f2, f3 = (counts[size] for size in range(2, 6))
        g3_values.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4_values.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    return tuple(g3_values), tuple(g4_values)


Y_KEYS = tuple(mask(triple) for triple in combinations(WORK_VERTICES, 3))
Z_KEYS = tuple(mask(pair) for pair in combinations(ISOLATED_VERTICES, 2))


@dataclass(frozen=True)
class BlockerTypeAudit:
    graph_cover_five_sets: int
    y_type_five_sets: int
    forbidden_y_01i_five_sets: int
    z_type_five_sets: int
    paired_z_variables: int
    every_graph_cover_classified_once: bool

    @property
    def verified(self) -> bool:
        return (
            self.graph_cover_five_sets == 310
            and self.y_type_five_sets == 220
            and self.forbidden_y_01i_five_sets == 10
            and self.z_type_five_sets == 90
            and self.paired_z_variables == 45
            and self.every_graph_cover_classified_once
        )

    def to_record(self) -> dict[str, int | bool]:
        return {
            "graph_cover_five_sets": self.graph_cover_five_sets,
            "y_type_five_sets": self.y_type_five_sets,
            "forbidden_y_01i_five_sets": self.forbidden_y_01i_five_sets,
            "z_type_five_sets": self.z_type_five_sets,
            "paired_z_variables": self.paired_z_variables,
            "every_graph_cover_classified_once": (
                self.every_graph_cover_classified_once
            ),
            "verified": self.verified,
        }


def audit_blocker_types() -> BlockerTypeAudit:
    """Classify all five-element covers into the symbolic ``y,z`` types."""

    special_mask = mask(SPECIAL_VERTICES)
    active_mask = mask(ACTIVE_VERTICES)
    isolate_mask = mask(ISOLATED_VERTICES)
    y_count = 0
    forbidden_y = 0
    z_count = 0
    cover_count = 0
    classified = True
    z_pairs: Counter[int] = Counter()
    for chosen in combinations(range(TARGET_VERTEX_COUNT), 5):
        candidate = mask(chosen)
        if not _is_graph_cover(candidate):
            continue
        cover_count += 1
        is_y = candidate & special_mask == special_mask
        one_special = (candidate & special_mask).bit_count() == 1
        is_z = (
            one_special
            and candidate & active_mask == active_mask
            and (candidate & isolate_mask).bit_count() == 2
        )
        if is_y == is_z:
            classified = False
            continue
        if is_y:
            triple = candidate ^ special_mask
            if triple.bit_count() != 3 or triple & ~mask(WORK_VERTICES):
                classified = False
                continue
            y_count += 1
            forbidden_y += triple & active_mask == active_mask
        else:
            pair = candidate & isolate_mask
            z_pairs[pair] += 1
            z_count += 1
    paired = all(multiplicity == 2 for multiplicity in z_pairs.values())
    audit = BlockerTypeAudit(
        graph_cover_five_sets=cover_count,
        y_type_five_sets=y_count,
        forbidden_y_01i_five_sets=forbidden_y,
        z_type_five_sets=z_count,
        paired_z_variables=len(z_pairs),
        every_graph_cover_classified_once=classified and paired,
    )
    if not audit.verified:
        raise AssertionError("record-16 blocker type reduction changed")
    return audit


def _candidate_facet_terms() -> tuple[tuple[int, tuple[str, int]], ...]:
    """Return every symbolic a=0 facet and its ``y`` or ``z`` variable.

    A ``y_T`` facet is ``{12,13} union T``.  Each ``z_S`` variable supplies
    the two facets ``{0,1,12} union S`` and ``{0,1,13} union S``.
    """

    special = mask(SPECIAL_VERTICES)
    active = mask(ACTIVE_VERTICES)
    terms: list[tuple[int, tuple[str, int]]] = []
    for triple in Y_KEYS:
        terms.append((special | triple, ("y", triple)))
    for pair in Z_KEYS:
        for special_vertex in SPECIAL_VERTICES:
            terms.append(
                (
                    active | (1 << special_vertex) | pair,
                    ("z", pair),
                )
            )
    return tuple(terms)


@dataclass(frozen=True)
class RidgeSystemAudit:
    total_six_sets: int
    nonempty_ridge_groups: int
    two_special_groups: int
    two_special_groups_with_doubled_z: int
    one_special_three_isolate_groups: int
    y_facet_six_set_incidences: int
    z_facet_six_set_incidences: int

    @property
    def verified(self) -> bool:
        return (
            self.total_six_sets == 3_003
            and self.nonempty_ridge_groups == 735
            and self.two_special_groups == 495
            and self.two_special_groups_with_doubled_z == 45
            and self.one_special_three_isolate_groups == 240
            and self.y_facet_six_set_incidences == 1_980
            and self.z_facet_six_set_incidences == 810
        )

    def to_record(self) -> dict[str, int | bool]:
        return {
            "total_six_sets": self.total_six_sets,
            "nonempty_ridge_groups": self.nonempty_ridge_groups,
            "two_special_groups": self.two_special_groups,
            "two_special_groups_with_doubled_z": (
                self.two_special_groups_with_doubled_z
            ),
            "one_special_three_isolate_groups": (
                self.one_special_three_isolate_groups
            ),
            "y_facet_six_set_incidences": self.y_facet_six_set_incidences,
            "z_facet_six_set_incidences": self.z_facet_six_set_incidences,
            "verified": self.verified,
        }


def audit_symbolic_ridge_system() -> RidgeSystemAudit:
    """Mechanically reconstruct the complete symbolic ridge system.

    The output verifies that no six-set equation is omitted from the proof:

    * ``{12,13} union Q`` has the four ``y`` variables indexed by the
      triples of ``Q`` and, when ``Q={0,1,i,j}``, two copies of ``z_ij``;
    * ``{p,0,1,i,j,k}`` has ``z_ij+z_ik+z_jk``;
    * every other six-set contains no candidate facet.
    """

    terms = _candidate_facet_terms()
    special_mask = mask(SPECIAL_VERTICES)
    active_mask = mask(ACTIVE_VERTICES)
    isolate_mask = mask(ISOLATED_VERTICES)
    nonempty = 0
    two_special = 0
    doubled_z = 0
    one_special = 0
    y_incidences = 0
    z_incidences = 0

    for chosen in combinations(range(TARGET_VERTEX_COUNT), 6):
        six_set = mask(chosen)
        observed = Counter(
            key for facet, key in terms if facet & six_set == facet
        )
        if observed:
            nonempty += 1
        y_incidences += sum(
            coefficient
            for (kind, _), coefficient in observed.items()
            if kind == "y"
        )
        z_incidences += sum(
            coefficient
            for (kind, _), coefficient in observed.items()
            if kind == "z"
        )

        special_count = (six_set & special_mask).bit_count()
        if special_count == 2:
            four_set = six_set ^ special_mask
            expected = Counter(
                ("y", mask(triple))
                for triple in combinations(vertices(four_set), 3)
            )
            if four_set & active_mask == active_mask:
                isolate_pair = four_set & isolate_mask
                if isolate_pair.bit_count() != 2:
                    raise AssertionError("malformed active four-set")
                expected[("z", isolate_pair)] += 2
                doubled_z += 1
            if observed != expected:
                raise AssertionError("two-special ridge reduction changed")
            two_special += 1
        elif special_count == 1 and six_set & active_mask == active_mask:
            isolate_triple = six_set & isolate_mask
            if isolate_triple.bit_count() == 3:
                expected = Counter(
                    ("z", mask(pair))
                    for pair in combinations(vertices(isolate_triple), 2)
                )
                if observed != expected:
                    raise AssertionError(
                        "one-special ridge reduction changed"
                    )
                one_special += 1
            elif observed:
                raise AssertionError("unexpected one-special ridge group")
        elif observed:
            raise AssertionError("unexpected symbolic ridge group")

    audit = RidgeSystemAudit(
        total_six_sets=3_003,
        nonempty_ridge_groups=nonempty,
        two_special_groups=two_special,
        two_special_groups_with_doubled_z=doubled_z,
        one_special_three_isolate_groups=one_special,
        y_facet_six_set_incidences=y_incidences,
        z_facet_six_set_incidences=z_incidences,
    )
    if not audit.verified:
        raise AssertionError("record-16 symbolic ridge audit changed")
    return audit


@dataclass(frozen=True)
class CutCase:
    smaller_cut_side: int
    larger_cut_side: int
    z_edge_count: int
    required_y_count: int
    cauchy_y_upper_bound: int
    contradiction_gap: int

    @property
    def excluded(self) -> bool:
        return self.contradiction_gap > 0

    def to_record(self) -> dict[str, int | bool]:
        return {
            "smaller_cut_side": self.smaller_cut_side,
            "larger_cut_side": self.larger_cut_side,
            "z_edge_count": self.z_edge_count,
            "required_y_count": self.required_y_count,
            "cauchy_y_upper_bound": self.cauchy_y_upper_bound,
            "contradiction_gap": self.contradiction_gap,
            "excluded": self.excluded,
        }


def cut_cases() -> tuple[CutCase, ...]:
    """Return all complete-bipartite ``z`` cases up to swapping sides."""

    result = []
    for smaller in range(6):
        larger = len(ISOLATED_VERTICES) - smaller
        z_edges = smaller * larger
        required_y = TARGET_A0_FACET_COUNT - 2 * z_edges
        upper = 9 * larger
        result.append(
            CutCase(
                smaller_cut_side=smaller,
                larger_cut_side=larger,
                z_edge_count=z_edges,
                required_y_count=required_y,
                cauchy_y_upper_bound=upper,
                contradiction_gap=required_y - upper,
            )
        )
    return tuple(result)


def _truth_tables() -> dict[str, object]:
    """Replay the local Boolean implications used in the proof."""

    z_triangle = tuple(
        {
            "edge_mask": value,
            "edge_count": value.bit_count(),
            "ridge_allowed": value.bit_count() in (0, 2),
        }
        for value in range(8)
    )
    x_triangle = tuple(
        {
            "x_edge_count": edge_count,
            "y_on_isolate_triple": edge_count % 2,
            "ridge_degree": edge_count + edge_count % 2,
            "ridge_allowed": edge_count + edge_count % 2 in (0, 2),
        }
        for edge_count in range(4)
    )
    q_01ij = tuple(
        {
            "z": z_value,
            "y_0ij": y0,
            "y_1ij": y1,
            "ridge_degree": 2 * z_value + y0 + y1,
            "ridge_allowed": 2 * z_value + y0 + y1 in (0, 2),
        }
        for z_value in range(2)
        for y0 in range(2)
        for y1 in range(2)
    )
    return {
        "z_triangle_edge_counts_allowed": [0, 2],
        "z_triangle_truth_table": list(z_triangle),
        "x_triangle_truth_table": list(x_triangle),
        "q_01ij_truth_table": list(q_01ij),
        "induced_2k2_y_count_on_four_set": 4,
        "induced_2k2_ridge_allowed": False,
        "two_cut_side_edges_induce_2k2": True,
    }


def _z_cut_replay() -> dict[str, object]:
    """Exhaust all labeled solutions of the ``z`` triple equations.

    Choose the nine edges from isolate 0 to the other isolates.  For every
    pair ``i,j>0``, the triple equation uniquely forces
    ``z_ij = z_0i xor z_0j``.  Thus the 512 seeds exhaust, rather than sample,
    the solution set.
    """

    local_vertices = tuple(range(len(ISOLATED_VERTICES)))
    local_pairs = tuple(combinations(local_vertices, 2))
    pair_index = {pair: index for index, pair in enumerate(local_pairs)}
    graphs: set[int] = set()
    edge_histogram: Counter[int] = Counter()
    for seed in range(1 << (len(local_vertices) - 1)):
        adjacency_to_zero = {
            vertex: bool(seed & (1 << (vertex - 1)))
            for vertex in local_vertices[1:]
        }
        graph = 0
        for first, second in local_pairs:
            if first == 0:
                edge = adjacency_to_zero[second]
            else:
                edge = (
                    adjacency_to_zero[first]
                    ^ adjacency_to_zero[second]
                )
            if edge:
                graph |= 1 << pair_index[(first, second)]
        for triple in combinations(local_vertices, 3):
            edge_count = sum(
                bool(graph & (1 << pair_index[pair]))
                for pair in combinations(triple, 2)
            )
            if edge_count not in (0, 2):
                raise AssertionError("forced z graph violates a triple")
        opposite_zero = {
            vertex
            for vertex, adjacent in adjacency_to_zero.items()
            if adjacent
        }
        zero_side = set(local_vertices) - opposite_zero
        expected_edges = {
            tuple(sorted((first, second)))
            for first in zero_side
            for second in opposite_zero
        }
        observed_edges = {
            pair
            for pair, index in pair_index.items()
            if graph & (1 << index)
        }
        if observed_edges != expected_edges:
            raise AssertionError("forced z graph is not the expected cut")
        graphs.add(graph)
        edge_histogram[len(observed_edges)] += 1
    if len(graphs) != 512 or sum(edge_histogram.values()) != 512:
        raise AssertionError("z cut replay did not exhaust 512 seeds")
    return {
        "seed_edge_count": 9,
        "seed_count": 512,
        "distinct_solution_count": len(graphs),
        "all_solutions_complete_bipartite": True,
        "edge_count_histogram": {
            str(edge_count): edge_histogram[edge_count]
            for edge_count in sorted(edge_histogram)
        },
    }


@dataclass(frozen=True)
class Record16A0Audit:
    graph6: str
    missing_edges: tuple[tuple[int, int], ...]
    minimal_graph_cover_masks: tuple[int, ...]
    base_vertex_link_g3: tuple[int, ...]
    base_vertex_link_g4: tuple[int, ...]
    missing_five_count: int
    target_facet_count: int
    blocker_types: BlockerTypeAudit
    ridge_system: RidgeSystemAudit
    truth_tables: Mapping[str, object]
    z_cut_replay: Mapping[str, object]
    cases: tuple[CutCase, ...]

    @property
    def verified(self) -> bool:
        return (
            self.graph6 == TARGET_GRAPH6
            and self.missing_edges == TARGET_MISSING_EDGES
            and self.minimal_graph_cover_masks
            == tuple(mask(cover) for cover in EXPECTED_MINIMAL_GRAPH_COVERS)
            and self.base_vertex_link_g3
            == (4, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1)
            and self.base_vertex_link_g4
            == (5, 5, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 1, 1)
            and self.missing_five_count == TARGET_A0_MISSING_FIVE_COUNT
            and self.target_facet_count == TARGET_A0_FACET_COUNT
            and self.blocker_types.verified
            and self.ridge_system.verified
            and self.z_cut_replay.get("distinct_solution_count") == 512
            and self.z_cut_replay.get("all_solutions_complete_bipartite")
            is True
            and [case.required_y_count for case in self.cases]
            == [110, 92, 78, 68, 62, 60]
            and [case.cauchy_y_upper_bound for case in self.cases]
            == [90, 81, 72, 63, 54, 45]
            and all(case.excluded for case in self.cases)
        )

    def to_record(self) -> dict[str, object]:
        return {
            "graph6": self.graph6,
            "missing_edges": [list(edge) for edge in self.missing_edges],
            "minimal_graph_cover_masks": list(self.minimal_graph_cover_masks),
            "minimal_graph_cover_vertices": [
                list(vertices(value)) for value in self.minimal_graph_cover_masks
            ],
            "base_vertex_link_g3": list(self.base_vertex_link_g3),
            "base_vertex_link_g4": list(self.base_vertex_link_g4),
            "a0_support_consequence": (
                "strict vertex-link g4 at vertices 12 and 13 forces every "
                "additional minimal five-face onto W={0,...,11}"
            ),
            "missing_four_count": 0,
            "missing_five_count": self.missing_five_count,
            "target_facet_count": self.target_facet_count,
            "symbolic_facet_variables": {
                "y_triples": len(Y_KEYS),
                "z_isolate_pairs": len(Z_KEYS),
                "facets_per_y_variable": 1,
                "facets_per_z_variable": 2,
                "facet_count_equation": "|y| + 2|z| = 110",
            },
            "blocker_type_audit": self.blocker_types.to_record(),
            "ridge_system": self.ridge_system.to_record(),
            "truth_tables": dict(self.truth_tables),
            "z_cut_exhaustive_replay": dict(self.z_cut_replay),
            "two_graph_count_formula": "|y| = 12e - sum(d_i^2)",
            "cauchy_bound": (
                "s|y| <= 12es - 4e^2 = "
                "9s^2 - (2e-3s)^2, hence |y| <= 9s"
            ),
            "cut_cases": [case.to_record() for case in self.cases],
            "conclusion": (
                "every possible complete-bipartite z cut contradicts the "
                "required facet count; the a=0 schedule is impossible"
            ),
            "verified": self.verified,
        }


def audit_record16_a0(
    graph6: str,
    schedules: Sequence[Mapping[str, object]],
) -> Record16A0Audit:
    """Build the exact, solver-independent audit of the ``a=0`` branch."""

    order, edges = decode_graph6(graph6)
    if order != TARGET_VERTEX_COUNT or tuple(edges) != TARGET_MISSING_EDGES:
        raise ValueError("record-16 missing-edge graph changed")
    a0 = [
        schedule
        for schedule in schedules
        if schedule.get("missing_four_count") == 0
    ]
    if len(a0) != 1:
        raise ValueError("record-16 source does not have one a=0 schedule")
    vector = a0[0].get("f_vector")
    if not isinstance(vector, (list, tuple)) or len(vector) != 9:
        raise ValueError("record-16 a=0 schedule has a malformed f-vector")
    if (
        a0[0].get("missing_five_count_offset")
        != TARGET_A0_MISSING_FIVE_COUNT
        or int(vector[-1]) != TARGET_A0_FACET_COUNT
    ):
        raise ValueError("record-16 a=0 schedule changed")

    base_g3, base_g4 = _base_vertex_link_g3_g4()
    audit = Record16A0Audit(
        graph6=graph6,
        missing_edges=tuple(edges),
        minimal_graph_cover_masks=minimal_graph_covers(),
        base_vertex_link_g3=base_g3,
        base_vertex_link_g4=base_g4,
        missing_five_count=TARGET_A0_MISSING_FIVE_COUNT,
        target_facet_count=TARGET_A0_FACET_COUNT,
        blocker_types=audit_blocker_types(),
        ridge_system=audit_symbolic_ridge_system(),
        truth_tables=_truth_tables(),
        z_cut_replay=_z_cut_replay(),
        cases=cut_cases(),
    )
    if not audit.verified:
        raise AssertionError("record-16 a=0 analytic certificate changed")
    return audit
