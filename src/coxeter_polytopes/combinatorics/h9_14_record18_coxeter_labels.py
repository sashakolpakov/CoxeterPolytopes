"""Unrestricted compact-Coxeter no-go for the selected record-18 nerve.

Two disjoint order-five minimal nonfaces must support Lanner diagrams in a
compact Coxeter realization.  This module exhausts their labelled placements
modulo a verified eight-element subgroup of the selected nerve symmetry and
uses exact actual-vertex ellipticity plus principal-block inertia.

Only this selected abstract nerve candidate is excluded, and only in the
compact case.  No finite-volume noncompact conclusion is made.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations
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
from coxeter_polytopes.combinatorics.h9_14_record18 import (
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    record18_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_selected_nerve_coxeter_no_go import (
    ExactIndexTwoWitness,
    _a5_witnesses,
    _canonical_lanner_states,
)


FACET_COUNT = 14
FULL_MASK = (1 << FACET_COUNT) - 1
LEFT_SUPPORT = TARGET_MISSING_FIVES[1]
RIGHT_SUPPORT = TARGET_MISSING_FIVES[7]


def _transposition(first: int, second: int) -> tuple[int, ...]:
    result = list(range(FACET_COUNT))
    result[first], result[second] = result[second], result[first]
    return tuple(result)


SYMMETRY_GENERATORS = (
    _transposition(5, 6),
    _transposition(7, 8),
    (0, 1, 3, 2, 9, 7, 8, 5, 6, 4, 11, 10, 12, 13),
)


def _compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(FACET_COUNT))


def _symmetry_group() -> tuple[dict[int, int], ...]:
    identity = tuple(range(FACET_COUNT))
    generated = {identity}
    frontier = [identity]
    while frontier:
        element = frontier.pop()
        for generator in SYMMETRY_GENERATORS:
            for composite in (
                _compose(element, generator),
                _compose(generator, element),
            ):
                if composite not in generated:
                    generated.add(composite)
                    frontier.append(composite)
    return tuple(
        dict(enumerate(permutation)) for permutation in sorted(generated)
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


def _all_states(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    edges = tuple(combinations(sorted(support), 2))
    return tuple(
        dict(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(tuple(support)))
    )


def _state_key(
    state: Mapping[Pair, int], support: Sequence[int]
) -> tuple[int, ...]:
    return tuple(
        state[pair(first, second)]
        for first, second in combinations(sorted(support), 2)
    )


def _move_state(
    state: Mapping[Pair, int], permutation: Mapping[int, int]
) -> dict[Pair, int]:
    return {
        pair(permutation[first], permutation[second]): label
        for (first, second), label in state.items()
    }


def _left_orbit_representatives(
    group: Sequence[Mapping[int, int]],
) -> tuple[tuple[dict[Pair, int], ...], tuple[tuple[int, int], ...]]:
    states = _all_states(LEFT_SUPPORT)
    seen = set()
    representatives = []
    sizes = []
    for state in states:
        if _state_key(state, LEFT_SUPPORT) in seen:
            continue
        orbit = {
            _state_key(_move_state(state, permutation), LEFT_SUPPORT)
            for permutation in group
        }
        seen.update(orbit)
        representatives.append(state)
        sizes.append(len(orbit))
    if len(seen) != len(states):
        raise AssertionError("record-18 left Lanner quotient is incomplete")
    return tuple(representatives), tuple(sorted(Counter(sizes).items()))


def _symmetry_verified(
    group: Sequence[Mapping[int, int]], vertices: Sequence[Sequence[int]]
) -> bool:
    nonfaces = {
        frozenset(face) for face in (*TARGET_MISSING_EDGES, *TARGET_MISSING_FIVES)
    }
    actual_vertices = {frozenset(vertex) for vertex in vertices}
    return all(
        {
            frozenset(permutation[vertex] for vertex in face)
            for face in nonfaces
        }
        == nonfaces
        and {
            frozenset(permutation[vertex] for vertex in actual)
            for actual in vertices
        }
        == actual_vertices
        and {permutation[vertex] for vertex in LEFT_SUPPORT} == set(LEFT_SUPPORT)
        and {permutation[vertex] for vertex in RIGHT_SUPPORT} == set(RIGHT_SUPPORT)
        for permutation in group
    )


def _restrictions(
    vertices: Sequence[Sequence[int]], union: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    selected = set(union)
    return tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in selected)
                for vertex in vertices
            },
            key=lambda item: (len(item), item),
        )
    )


def _labels(
    union: Sequence[int],
    left: Mapping[Pair, int],
    right: Mapping[Pair, int],
    cross: Mapping[Pair, int] | None = None,
) -> dict[Pair, int]:
    result = {pair(a, b): 2 for a, b in combinations(union, 2)}
    result.update(left)
    result.update(right)
    if cross:
        result.update(cross)
    return result


def _affected_finite(
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
    left: Mapping[Pair, int],
    right: Mapping[Pair, int],
    edge: Pair,
    permutation: Mapping[int, int],
) -> tuple[tuple[int, ...], tuple[int, ...], Pair]:
    return (
        _state_key(_move_state(left, permutation), LEFT_SUPPORT),
        _state_key(_move_state(right, permutation), RIGHT_SUPPORT),
        pair(permutation[edge[0]], permutation[edge[1]]),
    )


def _rooted_state_code(
    state: Mapping[Pair, int], support: Sequence[int], root: int
) -> tuple[int, ...]:
    """Canonical exact code of a labelled Lanner diagram with one root."""

    remainder = tuple(vertex for vertex in support if vertex != root)
    return min(
        tuple(
            state[pair(order[left], order[right])]
            for left, right in combinations(range(5), 2)
        )
        for image in permutations(remainder)
        for order in ((root, *image),)
    )


def _single_cross_block_code(
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    edge: Pair,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    left_root = edge[0] if edge[0] in LEFT_SUPPORT else edge[1]
    right_root = edge[1] if edge[1] in RIGHT_SUPPORT else edge[0]
    codes = (
        _rooted_state_code(left_state, LEFT_SUPPORT, left_root),
        _rooted_state_code(right_state, RIGHT_SUPPORT, right_root),
    )
    return tuple(sorted(codes))  # type: ignore[return-value]


@dataclass(frozen=True, slots=True)
class Record18CompactCoxeterAudit:
    actual_vertex_count: int
    symmetry_group_order: int
    symmetry_preserves_selected_nerve: bool
    full_left_state_count: int
    full_right_state_count: int
    left_orbit_representative_count: int
    left_state_orbit_size_histogram: tuple[tuple[int, int], ...]
    quotient_gate: DisjointLannerCrossGateAudit
    baseline_failure_count: int
    quotient_single_edge_trial_count: int
    full_single_edge_trial_count: int
    quotient_single_edge_passing_count: int
    passing_triple_orbit_count: int
    passing_triple_orbit_size_histogram: tuple[tuple[int, int], ...]
    full_single_edge_passing_count: int
    stronger_label_survivors: tuple[tuple[int, int], ...]
    quotient_raw_double_edge_trial_count: int
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
            "symmetry_generators": [list(value) for value in SYMMETRY_GENERATORS],
            "symmetry_group_order": self.symmetry_group_order,
            "symmetry_preserves_selected_nerve": self.symmetry_preserves_selected_nerve,
            "full_lanner_state_counts": [
                self.full_left_state_count,
                self.full_right_state_count,
            ],
            "left_orbit_representative_count": self.left_orbit_representative_count,
            "left_state_orbit_size_histogram": dict(
                self.left_state_orbit_size_histogram
            ),
            "right_states_retained_literally": self.full_right_state_count,
            "quotient_gate": self.quotient_gate.to_record(),
            "baseline_failure_count": self.baseline_failure_count,
            "quotient_single_edge_trial_count": self.quotient_single_edge_trial_count,
            "full_single_edge_trial_count": self.full_single_edge_trial_count,
            "quotient_single_edge_passing_count": self.quotient_single_edge_passing_count,
            "passing_triple_orbit_count": self.passing_triple_orbit_count,
            "passing_triple_orbit_size_histogram": dict(
                self.passing_triple_orbit_size_histogram
            ),
            "full_single_edge_passing_count": self.full_single_edge_passing_count,
            "stronger_label_survivors": dict(self.stronger_label_survivors),
            "quotient_raw_double_edge_trial_count": (
                self.quotient_raw_double_edge_trial_count
            ),
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


@lru_cache(maxsize=1)
def audit_record18_compact_coxeter_labels() -> Record18CompactCoxeterAudit:
    vertices = _actual_vertices(record18_facet_complements())
    group = _symmetry_group()
    if len(group) != 8 or not _symmetry_verified(group, vertices):
        raise AssertionError("record-18 selected-nerve symmetry changed")
    left_full = _all_states(LEFT_SUPPORT)
    right_states = _all_states(RIGHT_SUPPORT)
    left_states, left_orbits = _left_orbit_representatives(group)
    if (
        len(left_full),
        len(right_states),
        len(left_states),
        left_orbits,
    ) != (420, 420, 57, ((4, 9), (8, 48))):
        raise AssertionError("record-18 Lanner-state quotient changed")

    gate = audit_disjoint_lanner_cross_gate(
        vertices, LEFT_SUPPORT, RIGHT_SUPPORT, left_states, right_states
    )
    pass_histogram = tuple(
        sorted(Counter(value for _edge, value in gate.passing_count_by_cross_pair).items())
    )
    if (
        gate.unique_actual_restriction_size_histogram
        != ((6, 16), (7, 44), (8, 25))
        or gate.label_three_test_count != 598_500
        or len(gate.passing_trials) != 3_780
        or pass_histogram != ((84, 10), (132, 5), (180, 5), (276, 5))
        or gate.first_bad_restriction_size_histogram
        != ((6, 310_276), (7, 281_322), (8, 3_122))
    ):
        raise AssertionError("record-18 quotient cross-edge census changed")

    union = tuple(sorted((*LEFT_SUPPORT, *RIGHT_SUPPORT)))
    restrictions = _restrictions(vertices, union)
    # No actual vertex contains either minimal nonface.  With all cross pairs
    # right-angled, every restriction is therefore an orthogonal sum of
    # proper (hence elliptic) subdiagrams of the two Lanner diagrams.
    if any(
        set(LEFT_SUPPORT) <= set(restriction)
        or set(RIGHT_SUPPORT) <= set(restriction)
        for restriction in restrictions
    ):
        raise AssertionError("an actual vertex contains a Lanner support")
    baseline_failures = 0
    stronger = tuple(
        (
            label,
            sum(
                _affected_finite(
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
            ),
        )
        for label in (4, 5)
    )
    by_state: defaultdict[tuple[int, int], list[Pair]] = defaultdict(list)
    for left_index, right_index, edge in gate.passing_trials:
        by_state[(left_index, right_index)].append(edge)
    residual_double = sum(
        len(tuple(combinations(set(edges), 2))) for edges in by_state.values()
    )

    triple_orbits: dict[
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
        triple_orbits.setdefault(
            min(orbit), (left_index, right_index, edge, len(orbit))
        )

    witnesses = _a5_witnesses(LEFT_SUPPORT, RIGHT_SUPPORT)
    canonical_left = _canonical_lanner_states(LEFT_SUPPORT)
    canonical_right = _canonical_lanner_states(RIGHT_SUPPORT)
    references = (
        _single_cross_block_code(
            canonical_left[1],
            canonical_right[1],
            pair(LEFT_SUPPORT[0], RIGHT_SUPPORT[0]),
        ),
        _single_cross_block_code(
            canonical_left[4],
            canonical_right[1],
            pair(LEFT_SUPPORT[4], RIGHT_SUPPORT[0]),
        ),
    )
    names = (witnesses[0].name, witnesses[1].name)
    orbit_classes: Counter[str] = Counter()
    full_classes: Counter[str] = Counter()
    for left_index, right_index, edge, orbit_size in triple_orbits.values():
        block_code = _single_cross_block_code(
            left_states[left_index],
            right_states[right_index],
            edge,
        )
        matches = tuple(
            index
            for index, reference in enumerate(references)
            if block_code == reference
        )
        if len(matches) != 1:
            raise AssertionError("record-18 exceptional block classification changed")
        name = names[matches[0]]
        orbit_classes[name] += 1
        full_classes[name] += orbit_size

    triple_orbit_histogram = tuple(
        sorted(Counter(value[3] for value in triple_orbits.values()).items())
    )
    full_passing = sum(value[3] for value in triple_orbits.values())
    quotient_raw_double = len(left_states) * len(right_states) * 300
    full_raw_double = len(left_full) * len(right_states) * 300
    verified = (
        len(vertices) == 110
        and baseline_failures == 0
        and stronger == ((4, 0), (5, 0))
        and residual_double == 0
        and len(triple_orbits) == 3_660
        and triple_orbit_histogram == ((4, 120), (8, 3_540))
        and full_passing == 28_800
        and orbit_classes == {names[0]: 1_800, names[1]: 1_860}
        and full_classes == {names[0]: 14_400, names[1]: 14_400}
        and all(witness.inertia == (8, 2, 0) for witness in witnesses)
    )
    if not verified:
        raise AssertionError("record-18 compact-Coxeter audit failed")
    return Record18CompactCoxeterAudit(
        actual_vertex_count=len(vertices),
        symmetry_group_order=len(group),
        symmetry_preserves_selected_nerve=True,
        full_left_state_count=len(left_full),
        full_right_state_count=len(right_states),
        left_orbit_representative_count=len(left_states),
        left_state_orbit_size_histogram=left_orbits,
        quotient_gate=gate,
        baseline_failure_count=baseline_failures,
        quotient_single_edge_trial_count=gate.label_three_test_count,
        full_single_edge_trial_count=len(left_full) * len(right_states) * 25,
        quotient_single_edge_passing_count=len(gate.passing_trials),
        passing_triple_orbit_count=len(triple_orbits),
        passing_triple_orbit_size_histogram=triple_orbit_histogram,
        full_single_edge_passing_count=full_passing,
        stronger_label_survivors=stronger,
        quotient_raw_double_edge_trial_count=quotient_raw_double,
        full_raw_double_edge_trial_count=full_raw_double,
        quotient_residual_double_edge_trial_count=residual_double,
        exceptional_orbit_class_histogram=tuple(sorted(orbit_classes.items())),
        exceptional_full_class_histogram=tuple(sorted(full_classes.items())),
        exact_index_two_witnesses=witnesses,
        unrestricted_compact_coxeter_no_go=True,
    )


__all__ = [
    "LEFT_SUPPORT",
    "RIGHT_SUPPORT",
    "Record18CompactCoxeterAudit",
    "audit_record18_compact_coxeter_labels",
]
