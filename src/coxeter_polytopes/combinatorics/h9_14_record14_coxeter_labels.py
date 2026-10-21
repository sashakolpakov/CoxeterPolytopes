"""Unrestricted compact-Coxeter no-go for one record-14 completion.

The selected ordinary-polytopal record-14 nerve contains two disjoint
five-node minimal nonfaces.  In a compact Coxeter realization these support
Lanner diagrams, while every actual rank-nine vertex is elliptic.  This
module quotients all 420 labelled order-five Lanner placements on each
support by a verified ``S3 x S2 x S3`` nerve symmetry and applies the exact
local finite-diagram gate.  The residual combined Gram blocks all have
negative index two.

The conclusion concerns compact realizations of this selected completion
only.  It neither excludes other record-14 completions nor tests Vinberg's
finite-volume criterion for noncompact chambers.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations, product
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_disjoint_lanner_gate import (
    DisjointLannerCrossGateAudit,
    audit_disjoint_lanner_cross_gate,
)
from coxeter_polytopes.combinatorics.h9_14_record14 import (
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    record14_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_selected_nerve_coxeter_no_go import (
    ExactIndexTwoWitness,
    _a5_witnesses,
    _isomorphic_blocks,
    _labelled_block_graph,
)


FACET_COUNT = 14
FULL_MASK = (1 << FACET_COUNT) - 1
LEFT_SUPPORT = TARGET_MISSING_FIVES[2]
RIGHT_SUPPORT = TARGET_MISSING_FIVES[1]
SYMMETRY_BLOCKS = ((3, 4, 5), (6, 7), (9, 10, 11))
EXPECTED_QUOTIENT_PASS_HISTOGRAM = (
    (2, 1),
    (5, 1),
    (6, 2),
    (8, 2),
    (9, 1),
    (10, 1),
    (11, 1),
    (12, 2),
    (14, 3),
    (15, 1),
    (16, 1),
    (18, 2),
    (20, 1),
    (26, 2),
    (30, 1),
    (32, 2),
    (56, 1),
)


def _actual_vertices(complements: Iterable[int]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        sorted(
            tuple(
                facet
                for facet in range(FACET_COUNT)
                if (FULL_MASK ^ complement) >> facet & 1
            )
            for complement in complements
        )
    )


def _all_lanner_states(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    edges = tuple(combinations(sorted(support), 2))
    return tuple(
        dict(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(tuple(support)))
    )


def _symmetry_group() -> tuple[dict[int, int], ...]:
    result = []
    for images in product(
        *(tuple(permutations(block)) for block in SYMMETRY_BLOCKS)
    ):
        permutation = {vertex: vertex for vertex in range(FACET_COUNT)}
        for block, image in zip(SYMMETRY_BLOCKS, images, strict=True):
            permutation.update(zip(block, image, strict=True))
        result.append(permutation)
    return tuple(result)


def _state_key(
    state: Mapping[Pair, int], support: Sequence[int]
) -> tuple[int, ...]:
    return tuple(
        state[pair(left, right)]
        for left, right in combinations(sorted(support), 2)
    )


def _move_state(
    state: Mapping[Pair, int], permutation: Mapping[int, int]
) -> dict[Pair, int]:
    return {
        pair(permutation[left], permutation[right]): label
        for (left, right), label in state.items()
    }


def _state_orbit_representatives(
    support: Sequence[int], group: Sequence[Mapping[int, int]]
) -> tuple[tuple[dict[Pair, int], ...], tuple[tuple[int, int], ...]]:
    states = _all_lanner_states(support)
    seen: set[tuple[int, ...]] = set()
    representatives = []
    orbit_sizes = []
    for state in states:
        key = _state_key(state, support)
        if key in seen:
            continue
        orbit = {
            _state_key(_move_state(state, permutation), support)
            for permutation in group
        }
        seen.update(orbit)
        representatives.append(state)
        orbit_sizes.append(len(orbit))
    if len(seen) != len(states):
        raise AssertionError("record-14 Lanner symmetry quotient is incomplete")
    return (
        tuple(representatives),
        tuple(sorted(Counter(orbit_sizes).items())),
    )


def _moved_nonfaces(
    permutation: Mapping[int, int],
) -> set[frozenset[int]]:
    rows = (*TARGET_MISSING_EDGES, *TARGET_MISSING_FOURS, *TARGET_MISSING_FIVES)
    return {
        frozenset(permutation[vertex] for vertex in row) for row in rows
    }


def _symmetry_verified(
    group: Sequence[Mapping[int, int]], vertices: Sequence[Sequence[int]]
) -> bool:
    rows = (*TARGET_MISSING_EDGES, *TARGET_MISSING_FOURS, *TARGET_MISSING_FIVES)
    expected_rows = {frozenset(row) for row in rows}
    expected_vertices = {frozenset(vertex) for vertex in vertices}
    return all(
        _moved_nonfaces(permutation) == expected_rows
        and {
            frozenset(permutation[vertex] for vertex in actual_vertex)
            for actual_vertex in vertices
        }
        == expected_vertices
        for permutation in group
    )


def _restriction_family(
    vertices: Sequence[Sequence[int]], union: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    union_set = set(union)
    return tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in union_set)
                for vertex in vertices
            },
            key=lambda item: (len(item), item),
        )
    )


def _labels(
    union: Sequence[int],
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    cross: Mapping[Pair, int] | None = None,
) -> dict[Pair, int]:
    result = {pair(a, b): 2 for a, b in combinations(union, 2)}
    result.update(left_state)
    result.update(right_state)
    if cross:
        result.update(cross)
    return result


def _affected_restrictions_are_finite(
    restrictions: Sequence[Sequence[int]],
    labels: Mapping[Pair, int],
    edges: Sequence[Pair],
) -> bool:
    return all(
        is_finite_coxeter_labeling(restriction, labels)
        for restriction in restrictions
        if any(set(edge) <= set(restriction) for edge in edges)
    )


def _triple_key(
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    edge: Pair,
    permutation: Mapping[int, int],
) -> tuple[tuple[int, ...], tuple[int, ...], Pair]:
    return (
        _state_key(_move_state(left_state, permutation), LEFT_SUPPORT),
        _state_key(_move_state(right_state, permutation), RIGHT_SUPPORT),
        pair(permutation[edge[0]], permutation[edge[1]]),
    )


@dataclass(frozen=True, slots=True)
class Record14CompactCoxeterAudit:
    actual_vertex_count: int
    symmetry_group_order: int
    symmetry_preserves_nerve: bool
    left_full_state_count: int
    right_full_state_count: int
    left_orbit_representative_count: int
    right_orbit_representative_count: int
    left_state_orbit_size_histogram: tuple[tuple[int, int], ...]
    right_state_orbit_size_histogram: tuple[tuple[int, int], ...]
    quotient_gate: DisjointLannerCrossGateAudit
    baseline_failure_count: int
    full_single_edge_trial_count: int
    quotient_passing_trial_count: int
    passing_triple_orbit_count: int
    passing_triple_orbit_size_histogram: tuple[tuple[int, int], ...]
    full_single_edge_passing_count: int
    stronger_label_survivors: tuple[tuple[int, int], ...]
    full_raw_double_edge_trial_count: int
    quotient_residual_double_edge_trial_count: int
    exceptional_orbit_class_histogram: tuple[tuple[str, int], ...]
    exceptional_full_class_histogram: tuple[tuple[str, int], ...]
    exact_index_two_witnesses: tuple[ExactIndexTwoWitness, ...]
    unrestricted_compact_coxeter_no_go: bool

    def to_record(self) -> dict[str, object]:
        return {
            "rank_nine_actual_vertex_count": self.actual_vertex_count,
            "disjoint_lanner_supports": [list(LEFT_SUPPORT), list(RIGHT_SUPPORT)],
            "nerve_symmetry_blocks": [list(block) for block in SYMMETRY_BLOCKS],
            "symmetry_group_order": self.symmetry_group_order,
            "symmetry_preserves_nonfaces_and_vertices": self.symmetry_preserves_nerve,
            "left_full_lanner_state_count": self.left_full_state_count,
            "right_full_lanner_state_count": self.right_full_state_count,
            "left_orbit_representative_count": self.left_orbit_representative_count,
            "right_orbit_representative_count": self.right_orbit_representative_count,
            "left_state_orbit_size_histogram": dict(
                self.left_state_orbit_size_histogram
            ),
            "right_state_orbit_size_histogram": dict(
                self.right_state_orbit_size_histogram
            ),
            "quotient_gate": self.quotient_gate.to_record(),
            "baseline_failure_count": self.baseline_failure_count,
            "full_single_edge_trial_count": self.full_single_edge_trial_count,
            "quotient_passing_trial_count": self.quotient_passing_trial_count,
            "passing_triple_orbit_count": self.passing_triple_orbit_count,
            "passing_triple_orbit_size_histogram": dict(
                self.passing_triple_orbit_size_histogram
            ),
            "full_single_edge_passing_count": self.full_single_edge_passing_count,
            "stronger_label_survivors": dict(self.stronger_label_survivors),
            "full_raw_double_edge_trial_count": self.full_raw_double_edge_trial_count,
            "quotient_residual_double_edge_trial_count": (
                self.quotient_residual_double_edge_trial_count
            ),
            "exceptional_orbit_class_histogram": dict(
                self.exceptional_orbit_class_histogram
            ),
            "exceptional_full_class_histogram": dict(
                self.exceptional_full_class_histogram
            ),
            "exact_index_two_witnesses": [
                witness.to_record() for witness in self.exact_index_two_witnesses
            ],
            "ambient_compact_gram_inertia": [9, 1, 4],
            "unrestricted_compact_coxeter_no_go": (
                self.unrestricted_compact_coxeter_no_go
            ),
            "selected_nerve_only": True,
            "finite_volume_noncompact_checked": False,
        }


def audit_record14_compact_coxeter_labels() -> Record14CompactCoxeterAudit:
    vertices = _actual_vertices(record14_facet_complements())
    group = _symmetry_group()
    if len(group) != 72 or not _symmetry_verified(group, vertices):
        raise AssertionError("record-14 nerve symmetry changed")

    left_full = _all_lanner_states(LEFT_SUPPORT)
    right_full = _all_lanner_states(RIGHT_SUPPORT)
    left_states, left_orbits = _state_orbit_representatives(LEFT_SUPPORT, group)
    right_states, right_orbits = _state_orbit_representatives(RIGHT_SUPPORT, group)
    if (
        len(left_full),
        len(right_full),
        len(left_states),
        len(right_states),
        left_orbits,
        right_orbits,
    ) != (420, 420, 39, 73, ((6, 8), (12, 31)), ((3, 6), (6, 67))):
        raise AssertionError("record-14 Lanner orbit census changed")

    gate = audit_disjoint_lanner_cross_gate(
        vertices,
        LEFT_SUPPORT,
        RIGHT_SUPPORT,
        left_states,
        right_states,
    )
    pass_histogram = tuple(
        sorted(Counter(count for _edge, count in gate.passing_count_by_cross_pair).items())
    )
    if (
        gate.unique_actual_restriction_size_histogram != ((7, 47), (8, 25))
        or gate.label_three_test_count != 71_175
        or len(gate.passing_trials) != 420
        or pass_histogram != EXPECTED_QUOTIENT_PASS_HISTOGRAM
    ):
        raise AssertionError("record-14 quotient cross-edge census changed")

    union = tuple(sorted((*LEFT_SUPPORT, *RIGHT_SUPPORT)))
    restrictions = _restriction_family(vertices, union)
    baseline_failures = sum(
        not all(
            is_finite_coxeter_labeling(restriction, labels)
            for restriction in restrictions
        )
        for left_state in left_states
        for right_state in right_states
        for labels in (_labels(union, left_state, right_state),)
    )
    stronger = []
    for label in (4, 5):
        passing = sum(
            _affected_restrictions_are_finite(
                restrictions,
                _labels(
                    union,
                    left_states[left_index],
                    right_states[right_index],
                    {edge: label},
                ),
                (edge,),
            )
            for left_index, right_index, edge in gate.passing_trials
        )
        stronger.append((label, passing))

    by_state: defaultdict[tuple[int, int], list[Pair]] = defaultdict(list)
    for left_index, right_index, edge in gate.passing_trials:
        by_state[(left_index, right_index)].append(edge)
    residual_double = sum(
        len(tuple(combinations(set(edges), 2))) for edges in by_state.values()
    )

    orbit_representatives: dict[
        tuple[tuple[int, ...], tuple[int, ...], Pair],
        tuple[int, int, Pair, int],
    ] = {}
    for left_index, right_index, edge in gate.passing_trials:
        orbit = {
            _triple_key(
                left_states[left_index],
                right_states[right_index],
                edge,
                permutation,
            )
            for permutation in group
        }
        orbit_representatives.setdefault(
            min(orbit), (left_index, right_index, edge, len(orbit))
        )

    witnesses = _a5_witnesses(LEFT_SUPPORT, RIGHT_SUPPORT)
    # The witness helper performs the independent exact-minor assertion.
    # Build references from canonical states to classify every residual
    # weighted block up to isomorphism.
    from coxeter_polytopes.combinatorics.h9_14_selected_nerve_coxeter_no_go import (
        _canonical_lanner_states,
    )

    canonical_left = _canonical_lanner_states(LEFT_SUPPORT)
    canonical_right = _canonical_lanner_states(RIGHT_SUPPORT)
    references = (
        _labelled_block_graph(
            LEFT_SUPPORT,
            RIGHT_SUPPORT,
            canonical_left[1],
            canonical_right[1],
            pair(LEFT_SUPPORT[0], RIGHT_SUPPORT[0]),
        ),
        _labelled_block_graph(
            LEFT_SUPPORT,
            RIGHT_SUPPORT,
            canonical_left[4],
            canonical_right[1],
            pair(LEFT_SUPPORT[4], RIGHT_SUPPORT[0]),
        ),
    )
    orbit_classes: Counter[str] = Counter()
    full_classes: Counter[str] = Counter()
    names = (witnesses[0].name, witnesses[1].name)
    for left_index, right_index, edge, orbit_size in orbit_representatives.values():
        graph = _labelled_block_graph(
            LEFT_SUPPORT,
            RIGHT_SUPPORT,
            left_states[left_index],
            right_states[right_index],
            edge,
        )
        matches = tuple(
            index
            for index, reference in enumerate(references)
            if _isomorphic_blocks(graph, reference)
        )
        if len(matches) != 1:
            raise AssertionError("record-14 exceptional block classification changed")
        name = names[matches[0]]
        orbit_classes[name] += 1
        full_classes[name] += orbit_size

    full_single_trials = len(left_full) * len(right_full) * 25
    full_double_trials = len(left_full) * len(right_full) * 300
    orbit_size_histogram = tuple(
        sorted(Counter(value[3] for value in orbit_representatives.values()).items())
    )
    full_passing = sum(value[3] for value in orbit_representatives.values())
    verified = (
        len(vertices) == 106
        and baseline_failures == 0
        and tuple(stronger) == ((4, 0), (5, 0))
        and residual_double == 0
        and len(orbit_representatives) == 420
        and orbit_size_histogram == ((36, 40), (72, 380))
        and full_single_trials == 4_410_000
        and full_passing == 28_800
        and full_double_trials == 52_920_000
        and orbit_classes
        == {witnesses[0].name: 200, witnesses[1].name: 220}
        and full_classes
        == {witnesses[0].name: 14_400, witnesses[1].name: 14_400}
        and all(witness.inertia == (8, 2, 0) for witness in witnesses)
    )
    if not verified:
        raise AssertionError("record-14 compact-Coxeter no-go audit failed")

    return Record14CompactCoxeterAudit(
        actual_vertex_count=len(vertices),
        symmetry_group_order=len(group),
        symmetry_preserves_nerve=True,
        left_full_state_count=len(left_full),
        right_full_state_count=len(right_full),
        left_orbit_representative_count=len(left_states),
        right_orbit_representative_count=len(right_states),
        left_state_orbit_size_histogram=left_orbits,
        right_state_orbit_size_histogram=right_orbits,
        quotient_gate=gate,
        baseline_failure_count=baseline_failures,
        full_single_edge_trial_count=full_single_trials,
        quotient_passing_trial_count=len(gate.passing_trials),
        passing_triple_orbit_count=len(orbit_representatives),
        passing_triple_orbit_size_histogram=orbit_size_histogram,
        full_single_edge_passing_count=full_passing,
        stronger_label_survivors=tuple(stronger),
        full_raw_double_edge_trial_count=full_double_trials,
        quotient_residual_double_edge_trial_count=residual_double,
        exceptional_orbit_class_histogram=tuple(sorted(orbit_classes.items())),
        exceptional_full_class_histogram=tuple(sorted(full_classes.items())),
        exact_index_two_witnesses=witnesses,
        unrestricted_compact_coxeter_no_go=True,
    )


__all__ = [
    "LEFT_SUPPORT",
    "RIGHT_SUPPORT",
    "Record14CompactCoxeterAudit",
    "audit_record14_compact_coxeter_labels",
]
