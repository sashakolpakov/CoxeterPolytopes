"""Unrestricted compact-Coxeter audits for six displayed H9/14 nerves.

The selected ordinary-polytopal nerves 30, 32, 42, 43, 44, and 47 each
contain two disjoint Lanner supports.  Exact actual-vertex enumeration limits
the possible finite non-right edges between the supports.  If there is no
such edge, their orthogonal principal Gram has negative index two.  Every
locally possible exceptional single-edge block also has exact negative index
two, while two exceptional edges are locally impossible.  Cauchy interlacing
therefore excludes ambient Gram inertia ``(9,1,4)``.

The decisive objects are principal Gram blocks.  Labels involving nodes
outside a chosen block cannot change its inertia.  Thus it is enough to cover
every entry between the two Lanner supports inside that block: a forbidden
entry gives a monotonicity-stable bad actual-vertex subdiagram, no entries give
an orthogonal Lanner direct sum of negative index two, and the explicitly
surviving one-entry blocks also have exact negative index two.  Full support
symmetry orders are verified on the complete actual-vertex family before any
canonical Lanner placement list is used.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations
from typing import Iterable, Mapping, Sequence

import networkx as nx
import sympy as sp

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    lanner_isomorphism_codes,
    lanner_placement_codes,
    pair,
    representative_lanner_labels,
)
from coxeter_polytopes.combinatorics.h9_14_disjoint_lanner_gate import (
    DisjointLannerCrossGateAudit,
    audit_disjoint_lanner_cross_gate,
)
from coxeter_polytopes.combinatorics.h9_14_record32_coxeter_labels import (
    audit_record32_coxeter_labels,
)
from coxeter_polytopes.combinatorics.h9_14_record43 import (
    TARGET_MISSING_FIVES as RECORD43_FIVES,
    record43_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record44 import (
    TARGET_MISSING_FIVES as RECORD44_FIVES,
    record44_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record47 import (
    TARGET_MISSING_FIVES as RECORD47_FIVES,
    record47_facet_complements,
)


FACET_COUNT = 14
FULL_MASK = (1 << FACET_COUNT) - 1
RECORD30_NONFACES = (
    (0, 11),
    (1, 12),
    (2, 13),
    (3, 4, 5, 6),
    (7, 8, 9, 10),
)
RECORD42_NONFACES = (
    (0, 11),
    (0, 12),
    (1, 12),
    (1, 13),
    (11, 13),
    (2, 3, 4, 5),
    (6, 7, 8, 9, 10),
)


def _vertices_from_nonfaces(
    nonfaces: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    supports = tuple(set(face) for face in nonfaces)
    return tuple(
        subset
        for subset in combinations(range(FACET_COUNT), 9)
        if not any(support <= set(subset) for support in supports)
    )


def _vertices_from_complements(complements: Iterable[int]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        sorted(
            tuple(vertex for vertex in range(FACET_COUNT) if (FULL_MASK ^ mask) >> vertex & 1)
            for mask in complements
        )
    )


def _canonical_lanner_states(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    return tuple(
        representative_lanner_labels(support, type_index)
        for type_index in range(len(lanner_isomorphism_codes(len(tuple(support)))))
    )


def _all_lanner_states(support: Sequence[int]) -> tuple[dict[Pair, int], ...]:
    ordered = tuple(sorted(support))
    edges = tuple(combinations(ordered, 2))
    return tuple(
        dict(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(ordered))
    )


def _restricted_finite(
    vertices: Sequence[Sequence[int]],
    left: Sequence[int],
    right: Sequence[int],
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    cross_labels: Mapping[Pair, int],
) -> bool:
    # Local import keeps the reusable generic gate the only shared enumerator.
    from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
        is_finite_coxeter_labeling,
    )

    union = tuple(sorted((*left, *right)))
    labels = {pair(a, b): 2 for a, b in combinations(union, 2)}
    labels.update(left_state)
    labels.update(right_state)
    labels.update(cross_labels)
    affected_edges = tuple(cross_labels)
    restrictions = {
        tuple(facet for facet in vertex if facet in union)
        for vertex in vertices
        if any(set(edge) <= set(vertex) for edge in affected_edges)
    }
    return all(is_finite_coxeter_labeling(restriction, labels) for restriction in restrictions)


def _stronger_single_edge_survivors(
    gate: DisjointLannerCrossGateAudit,
    vertices: Sequence[Sequence[int]],
    left_states: Sequence[Mapping[Pair, int]],
    right_states: Sequence[Mapping[Pair, int]],
) -> tuple[tuple[int, int], ...]:
    counts = []
    for label in (4, 5):
        passing = sum(
            _restricted_finite(
                vertices,
                gate.left_support,
                gate.right_support,
                left_states[left_index],
                right_states[right_index],
                {edge: label},
            )
            for left_index, right_index, edge in gate.passing_trials
        )
        counts.append((label, passing))
    return tuple(counts)


def _double_edge_gate(
    gate: DisjointLannerCrossGateAudit,
    vertices: Sequence[Sequence[int]],
    left_states: Sequence[Mapping[Pair, int]],
    right_states: Sequence[Mapping[Pair, int]],
) -> tuple[int, int, int]:
    """Return raw trials, monotonicity-residual trials, residual survivors."""

    by_state: defaultdict[tuple[int, int], list[Pair]] = defaultdict(list)
    for left_index, right_index, edge in gate.passing_trials:
        by_state[(left_index, right_index)].append(edge)
    residual = 0
    passing = 0
    for (left_index, right_index), edges in by_state.items():
        for first, second in combinations(sorted(set(edges)), 2):
            residual += 1
            passing += _restricted_finite(
                vertices,
                gate.left_support,
                gate.right_support,
                left_states[left_index],
                right_states[right_index],
                {first: 3, second: 3},
            )
    raw = (
        len(left_states)
        * len(right_states)
        * len(tuple(combinations(gate.cross_pairs, 2)))
    )
    return raw, residual, passing


def _labelled_block_graph(
    left: Sequence[int],
    right: Sequence[int],
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    edge: Pair,
) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from((vertex, {"side": "L"}) for vertex in left)
    graph.add_nodes_from((vertex, {"side": "R"}) for vertex in right)
    labels = {pair(a, b): 2 for a, b in combinations(tuple((*left, *right)), 2)}
    labels.update(left_state)
    labels.update(right_state)
    labels[edge] = 3
    for a, b in combinations(tuple((*left, *right)), 2):
        graph.add_edge(a, b, label=labels[pair(a, b)])
    return graph


def _isomorphic_blocks(first: nx.Graph, second: nx.Graph) -> bool:
    return nx.is_isomorphic(
        first,
        second,
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("label", None),
    )


def _twice_gram(
    order: Sequence[int], labels: Mapping[Pair, int]
) -> sp.Matrix:
    sqrt2 = sp.sqrt(2)
    sqrt5 = sp.sqrt(5)
    value = {2: sp.Integer(0), 3: -sp.Integer(1), 4: -sqrt2, 5: -(1 + sqrt5) / 2}
    return sp.Matrix(
        [
            [
                sp.Integer(2)
                if left == right
                else value[labels[pair(left, right)]]
                for right in order
            ]
            for left in order
        ]
    )


def _leading_minors(order: Sequence[int], labels: Mapping[Pair, int]) -> tuple[sp.Expr, ...]:
    matrix = _twice_gram(order, labels)
    return tuple(sp.factor(matrix[:size, :size].det()) for size in range(1, len(order) + 1))


def _exact_linear_radical_sign(value: sp.Expr) -> str:
    """Prove a sign by rational enclosures, with no floating-point fallback."""

    bounds = {
        sp.sqrt(2): (sp.Rational(7, 5), sp.Rational(3, 2)),
        sp.sqrt(5): (sp.Rational(11, 5), sp.Rational(9, 4)),
        sp.sqrt(10): (sp.Rational(79, 25), sp.Rational(19, 6)),
    }
    for radicand, (lower, upper) in zip((2, 5, 10), bounds.values(), strict=True):
        if not lower**2 < radicand < upper**2:
            raise AssertionError("a rational radical enclosure changed")

    remainder = sp.expand(value)
    lower_bound = sp.Integer(0)
    upper_bound = sp.Integer(0)
    for radical, (lower, upper) in bounds.items():
        coefficient = sp.expand(remainder).coeff(radical)
        if coefficient.is_Rational is not True:
            raise AssertionError("a leading minor is not linear in the radicals")
        remainder = sp.expand(remainder - coefficient * radical)
        if coefficient >= 0:
            lower_bound += coefficient * lower
            upper_bound += coefficient * upper
        else:
            lower_bound += coefficient * upper
            upper_bound += coefficient * lower
    if remainder.is_Rational is not True:
        raise AssertionError("a leading minor has an unsupported radical term")
    lower_bound += remainder
    upper_bound += remainder
    if lower_bound > 0:
        return "+"
    if upper_bound < 0:
        return "-"
    raise AssertionError(
        f"rational radical bounds do not decide sign: {value}, "
        f"interval=({lower_bound}, {upper_bound})"
    )


def _exact_index_two_witness(
    name: str,
    order: Sequence[int],
    labels: Mapping[Pair, int],
    expected_minors: Sequence[sp.Expr],
    multiplicity: int,
) -> "ExactIndexTwoWitness":
    observed = _leading_minors(order, labels)
    expected = tuple(sp.sympify(value) for value in expected_minors)
    if any(
        sp.simplify(actual - wanted) != 0
        for actual, wanted in zip(observed, expected, strict=True)
    ):
        raise AssertionError(f"{name} exact leading minors changed")
    signs = "".join(_exact_linear_radical_sign(value) for value in expected)
    augmented_signs = "+" + signs  # Delta_0 = 1.
    negative_index = sum(
        left != right
        for left, right in zip(
            augmented_signs[:-1], augmented_signs[1:], strict=True
        )
    )
    inertia = (len(order) - negative_index, negative_index, 0)
    if negative_index != 2:
        raise AssertionError(f"{name} no longer has exact negative index two")
    return ExactIndexTwoWitness(
        name=name,
        order=tuple(order),
        leading_minors=tuple(sp.sstr(value) for value in expected),
        signs=signs,
        inertia=inertia,
        multiplicity_per_support_gate=multiplicity,
    )


@dataclass(frozen=True, slots=True)
class ExactIndexTwoWitness:
    name: str
    order: tuple[int, ...]
    leading_minors: tuple[str, ...]
    signs: str
    inertia: tuple[int, int, int]
    multiplicity_per_support_gate: int

    def to_record(self) -> dict[str, object]:
        return {
            "name": self.name,
            "twice_gram_facet_order": list(self.order),
            "leading_principal_minors": list(self.leading_minors),
            "signs": self.signs,
            "inertia": list(self.inertia),
            "multiplicity_per_support_gate": self.multiplicity_per_support_gate,
        }


@dataclass(frozen=True, slots=True)
class SelectedNerveAudit:
    record_index: int
    vertex_count: int
    support_gates: tuple[DisjointLannerCrossGateAudit, ...]
    stronger_single_edge_survivors: tuple[tuple[int, int], ...]
    raw_double_edge_trial_count: int
    monotonicity_residual_double_trial_count: int
    double_edge_passing_count: int
    exceptional_class_histogram: tuple[tuple[str, int], ...]
    exact_index_two_witnesses: tuple[ExactIndexTwoWitness, ...]
    normalizing_support_symmetry_orders: tuple[int, ...]
    unrestricted_compact_coxeter_no_go: bool

    def to_record(self) -> dict[str, object]:
        return {
            "record_index": self.record_index,
            "rank_nine_vertex_count": self.vertex_count,
            "support_gates": [gate.to_record() for gate in self.support_gates],
            "stronger_single_edge_survivors": {
                str(label): count for label, count in self.stronger_single_edge_survivors
            },
            "raw_double_edge_trial_count": self.raw_double_edge_trial_count,
            "monotonicity_residual_double_trial_count": (
                self.monotonicity_residual_double_trial_count
            ),
            "double_edge_passing_count": self.double_edge_passing_count,
            "exceptional_class_histogram": {
                name: count for name, count in self.exceptional_class_histogram
            },
            "exact_index_two_witnesses": [
                witness.to_record() for witness in self.exact_index_two_witnesses
            ],
            "normalizing_support_symmetry_orders": list(
                self.normalizing_support_symmetry_orders
            ),
            "unrestricted_compact_coxeter_no_go": (
                self.unrestricted_compact_coxeter_no_go
            ),
        }


def _pointwise_outside_support_symmetry_order(
    vertices: Sequence[Sequence[int]], support: Sequence[int]
) -> int:
    """Count support permutations preserving the full maximal-face family."""

    block = tuple(support)
    target = {frozenset(vertex) for vertex in vertices}
    order = 0
    for image_block in permutations(block):
        image = dict(zip(block, image_block, strict=True))
        transported = {
            frozenset(image.get(vertex, vertex) for vertex in face)
            for face in target
        }
        order += transported == target
    return order


def _full_block_labels(
    left: Sequence[int],
    right: Sequence[int],
    left_state: Mapping[Pair, int],
    right_state: Mapping[Pair, int],
    edge: Pair,
) -> dict[Pair, int]:
    union = tuple((*left, *right))
    labels = {pair(a, b): 2 for a, b in combinations(union, 2)}
    labels.update(left_state)
    labels.update(right_state)
    labels[edge] = 3
    return labels


def _record30_audit() -> SelectedNerveAudit:
    vertices = _vertices_from_nonfaces(RECORD30_NONFACES)
    left, right = RECORD30_NONFACES[-2:]
    left_states = _canonical_lanner_states(left)
    right_states = _canonical_lanner_states(right)
    symmetry_orders = (
        _pointwise_outside_support_symmetry_order(vertices, left),
        _pointwise_outside_support_symmetry_order(vertices, right),
    )
    gate = audit_disjoint_lanner_cross_gate(vertices, left, right, left_states, right_states)
    verified = (
        len(vertices) == 128
        and symmetry_orders == (24, 24)
        and gate.unique_actual_restriction_size_histogram == ((6, 16),)
        and gate.label_three_test_count == 1_296
        and gate.all_cross_edges_forbidden
    )
    if not verified:
        raise AssertionError("record-30 unrestricted support audit changed")
    return SelectedNerveAudit(
        record_index=30,
        vertex_count=len(vertices),
        support_gates=(gate,),
        stronger_single_edge_survivors=((4, 0), (5, 0)),
        raw_double_edge_trial_count=(
            len(left_states)
            * len(right_states)
            * len(tuple(combinations(gate.cross_pairs, 2)))
        ),
        monotonicity_residual_double_trial_count=0,
        double_edge_passing_count=0,
        exceptional_class_histogram=(),
        exact_index_two_witnesses=(),
        normalizing_support_symmetry_orders=symmetry_orders,
        unrestricted_compact_coxeter_no_go=True,
    )


def _record42_witness(
    left: Sequence[int], right: Sequence[int]
) -> ExactIndexTwoWitness:
    left_state = _canonical_lanner_states(left)[4]
    right_state = _canonical_lanner_states(right)[1]
    labels = _full_block_labels(left, right, left_state, right_state, pair(2, 6))
    order = (9, 3, 7, 4, 10, 5, 6, 8, 2)
    expected = (
        sp.Integer(2), sp.Integer(4), sp.Integer(6), sp.Integer(9),
        sp.Integer(12), sp.Integer(8), sp.Integer(10),
        8 - 4 * sp.sqrt(5),
        -8 * sp.sqrt(2) - 11 + 5 * sp.sqrt(5) + 4 * sp.sqrt(10),
    )
    return _exact_index_two_witness(
        name="record42_unique_single_edge_class",
        order=order,
        labels=labels,
        expected_minors=expected,
        multiplicity=2,
    )


def _record42_audit() -> SelectedNerveAudit:
    vertices = _vertices_from_nonfaces(RECORD42_NONFACES)
    left, right = RECORD42_NONFACES[-2:]
    left_states = _canonical_lanner_states(left)
    right_states = _canonical_lanner_states(right)
    symmetry_orders = (
        _pointwise_outside_support_symmetry_order(vertices, left),
        _pointwise_outside_support_symmetry_order(vertices, right),
    )
    gate = audit_disjoint_lanner_cross_gate(vertices, left, right, left_states, right_states)
    expected_trials = ((4, 1, pair(2, 6)), (4, 1, pair(4, 6)))
    stronger = _stronger_single_edge_survivors(gate, vertices, left_states, right_states)
    raw_double, residual_double, double_passing = _double_edge_gate(
        gate, vertices, left_states, right_states
    )
    # The two exceptional blocks are isomorphic by the displayed automorphism
    # of the canonical order-four Lanner diagram.
    q_auto = {2: 4, 3: 5, 4: 2, 5: 3}
    q_labels = left_states[4]
    auto_ok = all(
        q_labels[pair(a, b)] == q_labels[pair(q_auto[a], q_auto[b])]
        for a, b in combinations(left, 2)
    )
    witness = _record42_witness(left, right)
    verified = (
        len(vertices) == 100
        and symmetry_orders == (24, 120)
        and gate.label_three_test_count == 900
        and gate.passing_trials == expected_trials
        and stronger == ((4, 0), (5, 0))
        and raw_double == 8_550
        and residual_double == 1
        and double_passing == 0
        and auto_ok
        and witness.inertia == (7, 2, 0)
    )
    if not verified:
        raise AssertionError("record-42 unrestricted support audit changed")
    return SelectedNerveAudit(
        record_index=42,
        vertex_count=len(vertices),
        support_gates=(gate,),
        stronger_single_edge_survivors=stronger,
        raw_double_edge_trial_count=raw_double,
        monotonicity_residual_double_trial_count=residual_double,
        double_edge_passing_count=double_passing,
        exceptional_class_histogram=((witness.name, 2),),
        exact_index_two_witnesses=(witness,),
        normalizing_support_symmetry_orders=symmetry_orders,
        unrestricted_compact_coxeter_no_go=True,
    )


def _a5_witnesses(left: Sequence[int], right: Sequence[int]) -> tuple[ExactIndexTwoWitness, ...]:
    left_states = _canonical_lanner_states(left)
    right_states = _canonical_lanner_states(right)
    # The exact formulas are verified on record 44; every other support gate
    # is a relabeling of one of these two complete weighted blocks.
    labels_a = _full_block_labels(left, right, left_states[1], right_states[1], pair(left[0], right[0]))
    order_a = tuple((*left, *right))
    expected_a = (
        2, 4, 8, 6 - 2 * sp.sqrt(5), 4 - 2 * sp.sqrt(5),
        (9 - 5 * sp.sqrt(5)) / 2, 9 - 5 * sp.sqrt(5),
        18 - 10 * sp.sqrt(5), 26 - 12 * sp.sqrt(5),
        (25 - 11 * sp.sqrt(5)) / 2,
    )
    labels_b = _full_block_labels(left, right, left_states[4], right_states[1], pair(left[4], right[0]))
    order_b = tuple(reversed(tuple((*left, *right))))
    expected_b = (
        2, 4, 5 - sp.sqrt(5), (7 - 3 * sp.sqrt(5)) / 2,
        4 - 2 * sp.sqrt(5), (9 - 5 * sp.sqrt(5)) / 2,
        9 - 5 * sp.sqrt(5), 9 - 5 * sp.sqrt(5), 1 - sp.sqrt(5),
        (-23 - 16 * sp.sqrt(2) + 11 * sp.sqrt(5) + 8 * sp.sqrt(10)) / 2,
    )
    return (
        _exact_index_two_witness(
            name="A5_A5_type_1_1",
            order=order_a,
            labels=labels_a,
            expected_minors=expected_a,
            multiplicity=120,
        ),
        _exact_index_two_witness(
            name="A5_A5_type_4_1",
            order=order_b,
            labels=labels_b,
            expected_minors=expected_b,
            multiplicity=180,
        ),
    )


def _multiwedge_record_audit(
    record_index: int,
    missing_fives: Sequence[Sequence[int]],
    complements: Iterable[int],
) -> SelectedNerveAudit:
    vertices = _vertices_from_complements(complements)
    left = tuple(missing_fives[0])
    left_states = _canonical_lanner_states(left)
    symmetry_order = _pointwise_outside_support_symmetry_order(vertices, left)
    gates = []
    stronger_total = Counter({4: 0, 5: 0})
    raw_double = residual_double = double_passing = 0
    class_total = Counter()
    witnesses: tuple[ExactIndexTwoWitness, ...] | None = None
    for right_raw in missing_fives[1:]:
        right = tuple(right_raw)
        right_states = _all_lanner_states(right)
        gate = audit_disjoint_lanner_cross_gate(
            vertices, left, right, left_states, right_states
        )
        stronger = _stronger_single_edge_survivors(
            gate, vertices, left_states, right_states
        )
        raw, residual, passing = _double_edge_gate(
            gate, vertices, left_states, right_states
        )
        if witnesses is None:
            witnesses = _a5_witnesses(left, right)
        reference_graphs = []
        for witness_index, (left_type, right_type, left_endpoint) in enumerate(
            ((1, 1, left[0]), (4, 1, left[4]))
        ):
            reference_graphs.append(
                _labelled_block_graph(
                    left,
                    right,
                    left_states[left_type],
                    _canonical_lanner_states(right)[right_type],
                    pair(left_endpoint, right[0]),
                )
            )
        class_hist = Counter()
        for left_index, right_index, edge in gate.passing_trials:
            graph = _labelled_block_graph(
                left, right, left_states[left_index], right_states[right_index], edge
            )
            matches = [
                index for index, reference in enumerate(reference_graphs)
                if _isomorphic_blocks(graph, reference)
            ]
            if len(matches) != 1:
                raise AssertionError("an exceptional A5 block has no unique class")
            class_hist[witnesses[matches[0]].name] += 1
        gates.append(gate)
        stronger_total.update(dict(stronger))
        raw_double += raw
        residual_double += residual
        double_passing += passing
        class_total.update(class_hist)
        expected_restrictions = {
            43: (
                ((6, 15), (7, 20), (8, 25)),
                ((7, 35), (8, 25)),
                ((6, 15), (7, 20), (8, 25)),
            ),
            44: (((7, 20), (8, 25)),) * 3,
            47: (((7, 35), (8, 25)),) * 3,
        }[record_index][len(gates) - 1]
        if (
            gate.unique_actual_restriction_size_histogram != expected_restrictions
            or gate.label_three_test_count != 52_500
            or len(gate.passing_trials) != 300
            or tuple(sorted(Counter(count for _edge, count in gate.passing_count_by_cross_pair).items()))
            != ((0, 15), (24, 5), (36, 5))
            or stronger != ((4, 0), (5, 0))
            or raw != 630_000
            or passing != 0
            or class_hist != {
                "A5_A5_type_1_1": 120,
                "A5_A5_type_4_1": 180,
            }
        ):
            raise AssertionError(f"record-{record_index} support gate changed")
    if witnesses is None:
        raise AssertionError("multiwedge record has no residual A5 support")
    verified = (
        len(vertices) == 90
        and symmetry_order == 120
        and len(gates) == 3
        and stronger_total == {4: 0, 5: 0}
        and raw_double == 1_890_000
        and double_passing == 0
        and class_total == {
            "A5_A5_type_1_1": 360,
            "A5_A5_type_4_1": 540,
        }
        and all(witness.inertia == (8, 2, 0) for witness in witnesses)
    )
    if not verified:
        raise AssertionError(f"record-{record_index} unrestricted support audit changed")
    return SelectedNerveAudit(
        record_index=record_index,
        vertex_count=len(vertices),
        support_gates=tuple(gates),
        stronger_single_edge_survivors=tuple(sorted(stronger_total.items())),
        raw_double_edge_trial_count=raw_double,
        monotonicity_residual_double_trial_count=residual_double,
        double_edge_passing_count=double_passing,
        exceptional_class_histogram=tuple(sorted(class_total.items())),
        exact_index_two_witnesses=witnesses,
        normalizing_support_symmetry_orders=(symmetry_order,),
        unrestricted_compact_coxeter_no_go=True,
    )


@lru_cache(maxsize=1)
def audit_selected_nerve_coxeter_no_gos() -> tuple[SelectedNerveAudit, ...]:
    records = (
        _record30_audit(),
        # Record 32 uses all 546 overlapping-A4 joins before its cross gate.
        SelectedNerveAudit(
            record_index=32,
            vertex_count=110,
            support_gates=(),
            stronger_single_edge_survivors=((4, 0), (5, 0)),
            raw_double_edge_trial_count=0,
            monotonicity_residual_double_trial_count=0,
            double_edge_passing_count=0,
            exceptional_class_histogram=(),
            exact_index_two_witnesses=(),
            normalizing_support_symmetry_orders=(),
            unrestricted_compact_coxeter_no_go=(
                audit_record32_coxeter_labels().unrestricted_compact_coxeter_no_go
            ),
        ),
        _record42_audit(),
        _multiwedge_record_audit(43, RECORD43_FIVES, record43_facet_complements()),
        _multiwedge_record_audit(44, RECORD44_FIVES, record44_facet_complements()),
        _multiwedge_record_audit(47, RECORD47_FIVES, record47_facet_complements()),
    )
    if tuple(record.record_index for record in records) != (30, 32, 42, 43, 44, 47):
        raise AssertionError("selected-nerve record order changed")
    flags = {
        record.record_index: record.unrestricted_compact_coxeter_no_go
        for record in records
    }
    if flags != {30: True, 32: True, 42: True, 43: True, 44: True, 47: True}:
        raise AssertionError("the unrestricted selected-nerve scope changed")
    return records


__all__ = [
    "ExactIndexTwoWitness",
    "SelectedNerveAudit",
    "audit_selected_nerve_coxeter_no_gos",
]
