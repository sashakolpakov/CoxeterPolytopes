"""Unrestricted compact-Coxeter no-go for the selected record-31 nerve.

Two disjoint order-five minimal nonfaces are used as Lanner supports.  Four
verified clone transpositions reduce their complete labelled-placement
census to 110 states on each support.  Exact actual-vertex ellipticity then
leaves four one-cross-edge block classes.  The two large classes have a
principal Gram block of negative index two; the two small classes cannot be
extended consistently to Lanner diagrams on the other seven missing
five-faces.  Multiple cross edges fail the actual-vertex test.

The theorem concerns only this selected nerve and only compact Coxeter
realizations.  It does not exclude other completions of graph-shadow 31 or a
finite-volume noncompact realization with ideal faces.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations
from typing import Mapping, Sequence

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
from coxeter_polytopes.combinatorics.h9_14_record31 import (
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    target_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_selected_nerve_coxeter_no_go import (
    ExactIndexTwoWitness,
    _a5_witnesses,
    _canonical_lanner_states,
)


FACET_COUNT = 14
FULL_MASK = (1 << FACET_COUNT) - 1
LEFT_SUPPORT = TARGET_MISSING_FIVES[0]
RIGHT_SUPPORT = TARGET_MISSING_FIVES[1]
RESIDUAL_SUPPORTS = TARGET_MISSING_FIVES[2:]


def _transposition(first: int, second: int) -> tuple[int, ...]:
    result = list(range(FACET_COUNT))
    result[first], result[second] = result[second], result[first]
    return tuple(result)


SYMMETRY_GENERATORS = (
    _transposition(3, 4),
    _transposition(5, 6),
    _transposition(7, 8),
    _transposition(9, 10),
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
            composite = _compose(element, generator)
            if composite not in generated:
                generated.add(composite)
                frontier.append(composite)
    return tuple(dict(enumerate(value)) for value in sorted(generated))


def _actual_vertices() -> tuple[tuple[int, ...], ...]:
    return tuple(
        sorted(
            tuple(
                facet
                for facet in range(FACET_COUNT)
                if (FULL_MASK ^ complement) >> facet & 1
            )
            for complement in target_facet_complements()
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


def _orbit_representatives(
    support: Sequence[int], group: Sequence[Mapping[int, int]]
) -> tuple[tuple[dict[Pair, int], ...], tuple[tuple[int, int], ...]]:
    states = _all_states(support)
    seen: set[tuple[int, ...]] = set()
    representatives = []
    sizes = []
    for state in states:
        if _state_key(state, support) in seen:
            continue
        orbit = {
            _state_key(_move_state(state, permutation), support)
            for permutation in group
        }
        seen.update(orbit)
        representatives.append(state)
        sizes.append(len(orbit))
    if len(seen) != len(states):
        raise AssertionError("record-31 Lanner quotient is incomplete")
    return tuple(representatives), tuple(sorted(Counter(sizes).items()))


def _symmetry_verified(
    group: Sequence[Mapping[int, int]], vertices: Sequence[Sequence[int]]
) -> bool:
    nonfaces = {
        frozenset(face) for face in (*TARGET_MISSING_EDGES, *TARGET_MISSING_FIVES)
    }
    actual = {frozenset(vertex) for vertex in vertices}
    return all(
        {
            frozenset(permutation[vertex] for vertex in face)
            for face in nonfaces
        }
        == nonfaces
        and {
            frozenset(permutation[vertex] for vertex in vertex_set)
            for vertex_set in actual
        }
        == actual
        and {permutation[vertex] for vertex in LEFT_SUPPORT}
        == set(LEFT_SUPPORT)
        and {permutation[vertex] for vertex in RIGHT_SUPPORT}
        == set(RIGHT_SUPPORT)
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


def _rooted_state_code(
    state: Mapping[Pair, int], support: Sequence[int], root: int
) -> tuple[int, ...]:
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
    return tuple(
        sorted(
            (
                _rooted_state_code(left_state, LEFT_SUPPORT, left_root),
                _rooted_state_code(right_state, RIGHT_SUPPORT, right_root),
            )
        )
    )  # type: ignore[return-value]


def _exceptional_references() -> tuple[
    tuple[str, tuple[tuple[int, ...], tuple[int, ...]]], ...
]:
    left = _canonical_lanner_states(LEFT_SUPPORT)
    right = _canonical_lanner_states(RIGHT_SUPPORT)
    specifications = (
        ("A5_A5_type_0_1", 0, 0, 1, 0),
        ("A5_A5_type_1_1", 1, 0, 1, 0),
        ("A5_A5_type_1_1_opposite_roots", 1, 4, 1, 0),
        ("A5_A5_type_4_1", 4, 4, 1, 0),
    )
    result = []
    for name, left_type, left_root, right_type, right_root in specifications:
        result.append(
            (
                name,
                _single_cross_block_code(
                    left[left_type],
                    right[right_type],
                    pair(LEFT_SUPPORT[left_root], RIGHT_SUPPORT[right_root]),
                ),
            )
        )
    if len({code for _name, code in result}) != len(result):
        raise AssertionError("record-31 rooted block references collide")
    return tuple(result)


def _compatible_state(
    state: Mapping[Pair, int], fixed: Mapping[Pair, int]
) -> bool:
    return all(edge not in fixed or fixed[edge] == label for edge, label in state.items())


def _residual_lanner_completion_search(
    fixed: dict[Pair, int],
    state_domains: Sequence[Sequence[Mapping[Pair, int]]],
) -> tuple[bool, int]:
    """Exact MRV backtracking over the other seven Lanner supports."""

    node_count = 0

    def visit(remaining: tuple[int, ...]) -> bool:
        nonlocal node_count
        node_count += 1
        if not remaining:
            return True
        domains = []
        for index in remaining:
            domain = tuple(
                state
                for state in state_domains[index]
                if _compatible_state(state, fixed)
            )
            if not domain:
                return False
            domains.append((len(domain), index, domain))
        _size, chosen, domain = min(domains, key=lambda item: (item[0], item[1]))
        rest = tuple(index for index in remaining if index != chosen)
        for state in domain:
            added = {edge: label for edge, label in state.items() if edge not in fixed}
            fixed.update(added)
            if visit(rest):
                return True
            for edge in added:
                fixed.pop(edge)
        return False

    return visit(tuple(range(len(state_domains)))), node_count


@dataclass(frozen=True, slots=True)
class Record31CompactCoxeterAudit:
    actual_vertex_count: int
    symmetry_group_order: int
    symmetry_preserves_selected_nerve: bool
    full_left_state_count: int
    full_right_state_count: int
    left_orbit_representative_count: int
    right_orbit_representative_count: int
    left_state_orbit_size_histogram: tuple[tuple[int, int], ...]
    right_state_orbit_size_histogram: tuple[tuple[int, int], ...]
    quotient_gate: DisjointLannerCrossGateAudit
    stronger_label_survivors: tuple[tuple[int, int], ...]
    quotient_raw_double_edge_trial_count: int
    full_raw_double_edge_trial_count: int
    quotient_residual_double_edge_trial_count: int
    quotient_residual_double_edge_passing_count: int
    exceptional_class_histogram: tuple[tuple[str, int], ...]
    index_two_single_edge_trial_count: int
    residual_lanner_csp_trial_count: int
    residual_lanner_csp_completing_count: int
    residual_lanner_csp_node_histogram: tuple[tuple[int, int], ...]
    exact_index_two_witnesses: tuple[ExactIndexTwoWitness, ...]
    unrestricted_compact_coxeter_no_go: bool

    def to_record(self) -> dict[str, object]:
        gate = self.quotient_gate
        gate_summary = {
            "left_support": list(gate.left_support),
            "right_support": list(gate.right_support),
            "left_state_count": gate.left_state_count,
            "right_state_count": gate.right_state_count,
            "sparse_state_count": gate.sparse_state_count,
            "cross_pairs": [list(edge) for edge in gate.cross_pairs],
            "unique_actual_restriction_size_histogram": dict(
                gate.unique_actual_restriction_size_histogram
            ),
            "tested_label": 3,
            "label_three_test_count": gate.label_three_test_count,
            "passing_count": len(gate.passing_trials),
            "passing_count_by_cross_pair": [
                {"pair": list(edge), "passing_count": count}
                for edge, count in gate.passing_count_by_cross_pair
            ],
            "first_bad_restriction_size_histogram": dict(
                gate.first_bad_restriction_size_histogram
            ),
        }
        return {
            "rank_nine_actual_vertex_count": self.actual_vertex_count,
            "disjoint_lanner_supports": [
                list(LEFT_SUPPORT),
                list(RIGHT_SUPPORT),
            ],
            "residual_lanner_supports": [list(face) for face in RESIDUAL_SUPPORTS],
            "symmetry_generators": [list(value) for value in SYMMETRY_GENERATORS],
            "symmetry_group_order": self.symmetry_group_order,
            "symmetry_preserves_selected_nerve": self.symmetry_preserves_selected_nerve,
            "full_lanner_state_counts": [
                self.full_left_state_count,
                self.full_right_state_count,
            ],
            "orbit_representative_counts": [
                self.left_orbit_representative_count,
                self.right_orbit_representative_count,
            ],
            "left_state_orbit_size_histogram": dict(
                self.left_state_orbit_size_histogram
            ),
            "right_state_orbit_size_histogram": dict(
                self.right_state_orbit_size_histogram
            ),
            "quotient_gate": gate_summary,
            "stronger_label_survivors": dict(self.stronger_label_survivors),
            "quotient_raw_double_edge_trial_count": (
                self.quotient_raw_double_edge_trial_count
            ),
            "full_raw_double_edge_trial_count": self.full_raw_double_edge_trial_count,
            "quotient_residual_double_edge_trial_count": (
                self.quotient_residual_double_edge_trial_count
            ),
            "quotient_residual_double_edge_passing_count": (
                self.quotient_residual_double_edge_passing_count
            ),
            "exceptional_class_histogram": dict(self.exceptional_class_histogram),
            "index_two_single_edge_trial_count": self.index_two_single_edge_trial_count,
            "residual_lanner_csp_trial_count": self.residual_lanner_csp_trial_count,
            "residual_lanner_csp_completing_count": (
                self.residual_lanner_csp_completing_count
            ),
            "residual_lanner_csp_node_histogram": dict(
                self.residual_lanner_csp_node_histogram
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
def audit_record31_compact_coxeter_labels() -> Record31CompactCoxeterAudit:
    vertices = _actual_vertices()
    group = _symmetry_group()
    if len(group) != 16 or not _symmetry_verified(group, vertices):
        raise AssertionError("record-31 selected-nerve symmetry changed")

    left_full = _all_states(LEFT_SUPPORT)
    right_full = _all_states(RIGHT_SUPPORT)
    left_states, left_orbits = _orbit_representatives(LEFT_SUPPORT, group)
    right_states, right_orbits = _orbit_representatives(RIGHT_SUPPORT, group)
    expected_orbits = ((2, 10), (4, 100))
    if (
        len(left_full),
        len(right_full),
        len(left_states),
        len(right_states),
        left_orbits,
        right_orbits,
    ) != (420, 420, 110, 110, expected_orbits, expected_orbits):
        raise AssertionError("record-31 Lanner-state quotient changed")

    gate = audit_disjoint_lanner_cross_gate(
        vertices,
        LEFT_SUPPORT,
        RIGHT_SUPPORT,
        left_states,
        right_states,
    )
    pass_histogram = tuple(
        sorted(
            Counter(
                count for _edge, count in gate.passing_count_by_cross_pair
            ).items()
        )
    )
    if (
        gate.unique_actual_restriction_size_histogram
        != ((6, 28), (7, 20), (8, 21))
        or gate.label_three_test_count != 302_500
        or len(gate.passing_trials) != 1_932
        or pass_histogram
        != (
            (16, 1),
            (23, 2),
            (33, 1),
            (38, 1),
            (49, 2),
            (54, 1),
            (56, 3),
            (69, 2),
            (72, 1),
            (79, 2),
            (84, 1),
            (102, 1),
            (117, 1),
            (118, 1),
            (120, 1),
            (136, 1),
            (137, 2),
            (160, 1),
        )
        or gate.first_bad_restriction_size_histogram
        != ((6, 214_866), (7, 17_754), (8, 67_948))
    ):
        raise AssertionError("record-31 quotient cross-edge census changed")

    union = tuple(sorted((*LEFT_SUPPORT, *RIGHT_SUPPORT)))
    restrictions = _restrictions(vertices, union)
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

    by_state: defaultdict[tuple[int, int], set[Pair]] = defaultdict(set)
    for left_index, right_index, edge in gate.passing_trials:
        by_state[(left_index, right_index)].add(edge)
    residual_pairs = tuple(
        (left_index, right_index, edges)
        for (left_index, right_index), passing_edges in sorted(by_state.items())
        for edges in combinations(sorted(passing_edges), 2)
    )
    double_passing = sum(
        _affected_finite(
            restrictions,
            _labels(
                union,
                left_states[left_index],
                right_states[right_index],
                {edge: 3 for edge in edges},
            ),
            edges,
        )
        for left_index, right_index, edges in residual_pairs
    )

    references = _exceptional_references()
    classes: Counter[str] = Counter()
    classified_trials = []
    for left_index, right_index, edge in gate.passing_trials:
        code = _single_cross_block_code(
            left_states[left_index], right_states[right_index], edge
        )
        matches = tuple(name for name, reference in references if code == reference)
        if len(matches) != 1:
            raise AssertionError("record-31 exceptional block classification changed")
        name = matches[0]
        classes[name] += 1
        classified_trials.append((left_index, right_index, edge, name))

    residual_names = {
        "A5_A5_type_0_1",
        "A5_A5_type_1_1_opposite_roots",
    }
    state_domains = tuple(_all_states(support) for support in RESIDUAL_SUPPORTS)
    csp_nodes = Counter()
    csp_completing = 0
    csp_trials = 0
    for left_index, right_index, edge, name in classified_trials:
        if name not in residual_names:
            continue
        csp_trials += 1
        fixed = _labels(
            union,
            left_states[left_index],
            right_states[right_index],
            {edge: 3},
        )
        completing, nodes = _residual_lanner_completion_search(fixed, state_domains)
        csp_nodes[nodes] += 1
        csp_completing += completing

    witnesses = _a5_witnesses(LEFT_SUPPORT, RIGHT_SUPPORT)
    index_two_count = classes["A5_A5_type_1_1"] + classes["A5_A5_type_4_1"]
    expected_classes = {
        "A5_A5_type_0_1": 36,
        "A5_A5_type_1_1": 900,
        "A5_A5_type_1_1_opposite_roots": 36,
        "A5_A5_type_4_1": 960,
    }
    quotient_raw_double = len(left_states) * len(right_states) * 300
    full_raw_double = len(left_full) * len(right_full) * 300
    verified = (
        len(vertices) == 122
        and stronger == ((4, 0), (5, 0))
        and len(residual_pairs) == 36
        and double_passing == 0
        and classes == expected_classes
        and index_two_count == 1_860
        and csp_trials == 72
        and csp_completing == 0
        and all(witness.inertia == (8, 2, 0) for witness in witnesses)
    )
    if not verified:
        raise AssertionError("record-31 compact-Coxeter audit failed")
    return Record31CompactCoxeterAudit(
        actual_vertex_count=len(vertices),
        symmetry_group_order=len(group),
        symmetry_preserves_selected_nerve=True,
        full_left_state_count=len(left_full),
        full_right_state_count=len(right_full),
        left_orbit_representative_count=len(left_states),
        right_orbit_representative_count=len(right_states),
        left_state_orbit_size_histogram=left_orbits,
        right_state_orbit_size_histogram=right_orbits,
        quotient_gate=gate,
        stronger_label_survivors=stronger,
        quotient_raw_double_edge_trial_count=quotient_raw_double,
        full_raw_double_edge_trial_count=full_raw_double,
        quotient_residual_double_edge_trial_count=len(residual_pairs),
        quotient_residual_double_edge_passing_count=double_passing,
        exceptional_class_histogram=tuple(sorted(classes.items())),
        index_two_single_edge_trial_count=index_two_count,
        residual_lanner_csp_trial_count=csp_trials,
        residual_lanner_csp_completing_count=csp_completing,
        residual_lanner_csp_node_histogram=tuple(sorted(csp_nodes.items())),
        exact_index_two_witnesses=witnesses,
        unrestricted_compact_coxeter_no_go=True,
    )


__all__ = [
    "LEFT_SUPPORT",
    "RIGHT_SUPPORT",
    "Record31CompactCoxeterAudit",
    "audit_record31_compact_coxeter_labels",
]
