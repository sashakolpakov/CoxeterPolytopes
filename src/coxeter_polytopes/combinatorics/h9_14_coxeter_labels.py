"""Exact Coxeter-label gates for compact ``H^9`` fourteen-facet nerves.

Burcroff's no-multi-multiple-edge lemma reduces every finite dihedral label
in the target problem to ``2, 3, 4, 5``.  This module treats that full domain.
It is combinatorial: dotted pairs are recorded as missing edges, but their
continuous Gram entries and the global Lorentzian Gram equations are not
solved here.

The two local inputs are classical classifications.

* In the Burcroff-reduced alphabet, a vertex link is elliptic precisely when
  every connected component of its non-right labelled graph is one of
  ``A, B, D, E, F, H`` or ``I_2(m)`` for ``m=3,4,5``.
* Lanner diagrams of orders four and five have respectively nine and five
  isomorphism types.  Their complete labelled placement orbits contain 108
  and 420 elements.

Thus an accepted labelling proves every rank-nine vertex diagram finite and
every non-edge minimal nonface Lanner.  It does *not* prove that choices for
the dotted entries extend to one Gram matrix of inertia ``(9, 1, 4)``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations, product
from typing import Iterable, Iterator, Mapping, Sequence


Pair = tuple[int, int]
LabelMap = Mapping[Pair, int]

FINITE_LABELS = (2, 3, 4, 5)
NONRIGHT_LABELS = (3, 4, 5)


def pair(left: int, right: int) -> Pair:
    """Return a canonically ordered pair of distinct facet labels."""

    if left == right:
        raise ValueError("a Coxeter pair needs two distinct facets")
    return (left, right) if left < right else (right, left)


def _path_labels(
    component: frozenset[int], adjacency: Mapping[int, tuple[tuple[int, int], ...]]
) -> tuple[int, ...]:
    endpoint = min(
        vertex
        for vertex in component
        if len(adjacency[vertex]) == 1
    )
    labels: list[int] = []
    previous: int | None = None
    current = endpoint
    while True:
        following = [
            (vertex, label)
            for vertex, label in adjacency[current]
            if vertex != previous
        ]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("a path traversal branched")
        next_vertex, label = following[0]
        labels.append(label)
        previous, current = current, next_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def _branch_arm_lengths(
    component: frozenset[int], adjacency: Mapping[int, tuple[tuple[int, int], ...]]
) -> tuple[int, ...]:
    branch = next(
        vertex for vertex in component if len(adjacency[vertex]) == 3
    )
    lengths: list[int] = []
    for neighbor, _label in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(
                vertex
                for vertex, _ in adjacency[current]
                if vertex != previous
            )
            previous, current = current, following
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def _finite_connected_component(
    component: frozenset[int], labels: LabelMap
) -> bool:
    """Recognize one irreducible finite Coxeter diagram exactly."""

    order = len(component)
    if order == 1:
        return True

    adjacency: dict[int, tuple[tuple[int, int], ...]] = {}
    for vertex in component:
        adjacency[vertex] = tuple(
            (other, labels[pair(vertex, other)])
            for other in component
            if other != vertex and labels[pair(vertex, other)] > 2
        )
    edge_count = sum(map(len, adjacency.values())) // 2
    if edge_count != order - 1:
        return False
    degrees = Counter(len(adjacency[vertex]) for vertex in component)

    if max(degrees) <= 2:
        sequence = _path_labels(component, adjacency)
        if order == 2:
            return sequence[0] in NONRIGHT_LABELS  # A2, B2, or I2(5).
        if all(label == 3 for label in sequence):
            return True  # A_r.
        if sequence == (3,) * (order - 2) + (4,):
            return True  # B_r, with the multiple edge at an endpoint.
        if order == 3 and sequence == (3, 5):
            return True  # H3.
        if order == 4 and sequence == (3, 3, 5):
            return True  # H4.
        if order == 4 and sequence == (3, 4, 3):
            return True  # F4.
        return False

    if degrees[3] != 1 or any(degree > 3 for degree in degrees):
        return False
    if any(
        labels[pair(left, right)] != 3
        for left in component
        for right in component
        if left < right and labels[pair(left, right)] > 2
    ):
        return False
    arms = _branch_arm_lengths(component, adjacency)
    return arms == (1, 1, order - 3) or arms in {
        (1, 2, 2),  # E6.
        (1, 2, 3),  # E7.
        (1, 2, 4),  # E8.
    }


def is_finite_coxeter_labeling(
    facets: Iterable[int], labels: LabelMap
) -> bool:
    """Whether the induced H9/14-reduced labelled diagram is finite/elliptic.

    The function is complete for Burcroff's target-specific alphabet
    ``{2,3,4,5}``; it intentionally rejects labels outside that alphabet even
    though ``I_2(m)`` is finite for arbitrary ``m``.  No floating-point Gram
    test is used.
    """

    vertices = tuple(sorted(set(facets)))
    if any(
        pair(left, right) not in labels
        or labels[pair(left, right)] not in FINITE_LABELS
        for left, right in combinations(vertices, 2)
    ):
        return False

    unseen = set(vertices)
    components: list[frozenset[int]] = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        frontier = [start]
        while frontier:
            vertex = frontier.pop()
            neighbors = {
                other
                for other in unseen
                if labels[pair(vertex, other)] > 2
            }
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        components.append(frozenset(component))
    return all(_finite_connected_component(component, labels) for component in components)


# Representatives for the classical nine order-four and five order-five
# Lanner diagrams.  Unlisted pairs have label two.  The numbering is local to
# this certificate and makes no assertion about numbering in another table.
_LANNER_REPRESENTATIVES: dict[int, tuple[tuple[tuple[Pair, int], ...], ...]] = {
    4: (
        (((0, 1), 3), ((0, 2), 3), ((0, 3), 5)),
        (((0, 1), 5), ((0, 2), 3), ((1, 3), 3)),
        (((0, 1), 3), ((0, 2), 4), ((1, 3), 5)),
        (((0, 1), 3), ((0, 2), 5), ((1, 3), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 3), 4)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 3), 5)),
        (((0, 1), 3), ((0, 2), 4), ((1, 3), 4), ((2, 3), 3)),
        (((0, 1), 3), ((0, 2), 4), ((1, 3), 5), ((2, 3), 3)),
        (((0, 1), 3), ((0, 2), 5), ((1, 3), 5), ((2, 3), 3)),
    ),
    5: (
        (((0, 1), 3), ((0, 2), 3), ((0, 3), 3), ((1, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 4), ((2, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 5), ((2, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 3), ((3, 4), 4)),
    ),
}


def _representative_code(order: int, edges: Sequence[tuple[Pair, int]]) -> tuple[int, ...]:
    labels = {pair(left, right): 2 for left, right in combinations(range(order), 2)}
    labels.update(edges)
    return tuple(labels[pair(left, right)] for left, right in combinations(range(order), 2))


def _permuted_code(
    order: int, code: tuple[int, ...], permutation: tuple[int, ...]
) -> tuple[int, ...]:
    source_pairs = tuple(combinations(range(order), 2))
    source = dict(zip(source_pairs, code, strict=True))
    target = {pair(permutation[left], permutation[right]): label for (left, right), label in source.items()}
    return tuple(target[pair(left, right)] for left, right in source_pairs)


@lru_cache(maxsize=None)
def lanner_placement_codes(order: int) -> tuple[tuple[int, ...], ...]:
    """Return every labelled placement of every Lanner type of this order."""

    if order not in _LANNER_REPRESENTATIVES:
        raise ValueError("the H9/14 gate uses Lanner diagrams of order four or five")
    placements = {
        _permuted_code(order, _representative_code(order, representative), permutation)
        for representative in _LANNER_REPRESENTATIVES[order]
        for permutation in permutations(range(order))
    }
    expected = {4: 108, 5: 420}[order]
    if len(placements) != expected:
        raise AssertionError(f"the order-{order} Lanner placement census changed")
    return tuple(sorted(placements))


@lru_cache(maxsize=None)
def lanner_isomorphism_codes(order: int) -> tuple[tuple[int, ...], ...]:
    """Return one canonical code for each Lanner isomorphism type."""

    representatives = {
        min(
            _permuted_code(order, code, permutation)
            for permutation in permutations(range(order))
        )
        for code in lanner_placement_codes(order)
    }
    expected = {4: 9, 5: 5}[order]
    if len(representatives) != expected:
        raise AssertionError(f"the order-{order} Lanner isomorphism census changed")
    return tuple(sorted(representatives))


def induced_label_code(facets: Iterable[int], labels: LabelMap) -> tuple[int, ...]:
    vertices = tuple(sorted(facets))
    return tuple(labels[pair(left, right)] for left, right in combinations(vertices, 2))


def is_lanner_coxeter_labeling(facets: Iterable[int], labels: LabelMap) -> bool:
    vertices = tuple(sorted(set(facets)))
    if len(vertices) not in (4, 5):
        return False
    if any(pair(left, right) not in labels for left, right in combinations(vertices, 2)):
        return False
    return induced_label_code(vertices, labels) in lanner_placement_codes(len(vertices))


def representative_lanner_labels(
    facets: Iterable[int], type_index: int = 0
) -> dict[Pair, int]:
    """Place one canonical Lanner type on the specified four/five facets."""

    vertices = tuple(sorted(set(facets)))
    types = lanner_isomorphism_codes(len(vertices))
    if not 0 <= type_index < len(types):
        raise IndexError("Lanner type index is out of range")
    return {
        pair(vertices[left], vertices[right]): label
        for (left, right), label in zip(
            combinations(range(len(vertices)), 2), types[type_index], strict=True
        )
    }


@dataclass(frozen=True, slots=True)
class H914CoxeterLabelAudit:
    facet_count: int
    vertex_count: int
    meeting_pair_count: int
    dotted_pair_count: int
    finite_label_histogram: tuple[tuple[int, int], ...]
    minimal_nonface_order_histogram: tuple[tuple[int, int], ...]
    bad_vertex_indices: tuple[int, ...]
    bad_lanner_nonfaces: tuple[tuple[int, ...], ...]
    missing_or_extra_label_pairs: tuple[Pair, ...]

    @property
    def passes_discrete_gate(self) -> bool:
        return (
            not self.bad_vertex_indices
            and not self.bad_lanner_nonfaces
            and not self.missing_or_extra_label_pairs
        )


@dataclass(frozen=True, slots=True)
class H914ProductLabelRepresentative:
    """One symmetry representative in a product nerve's orthogonal family.

    ``labels`` lists all meeting pairs, including those labelled two, so the
    object can be passed directly to the local audit or a later Gram solver.
    The symmetry quotient is by the full minimal-nonface-hypergraph
    automorphism group: internal permutations of simplex blocks, interchange
    of the two equal four-blocks in record 30, and the dihedral action on the
    pentagon block in record 42.
    """

    family: str
    lanner_type_indices: tuple[int, ...]
    cycle_color_code: tuple[int, ...]
    labels: tuple[tuple[int, int, int], ...]

    def label_map(self) -> dict[Pair, int]:
        return {(left, right): label for left, right, label in self.labels}


@lru_cache(maxsize=None)
def five_cycle_color_orbit_codes() -> tuple[tuple[int, ...], ...]:
    """One code for each four-coloring of a 5-cycle modulo ``D_5``."""

    representatives = set()
    for raw in product(FINITE_LABELS, repeat=5):
        rotations = [raw[index:] + raw[:index] for index in range(5)]
        reverse = tuple(reversed(raw))
        rotations.extend(reverse[index:] + reverse[:index] for index in range(5))
        representatives.add(min(rotations))
    if len(representatives) != 136:
        raise AssertionError("the four-color length-five bracelet census changed")
    return tuple(sorted(representatives))


def minimal_nonface_component_partition(
    minimal_nonfaces: Sequence[Iterable[int]],
) -> tuple[tuple[int, ...], ...]:
    """Return components of the minimal-nonface intersection hypergraph.

    Every order-two minimal nonface is a connected dotted Coxeter diagram and
    every higher minimal nonface is a connected Lanner diagram.  Overlapping
    supports therefore lie in one connected component of the nonzero Gram
    graph before any additional non-right meeting edge is added.
    """

    supports = [set(face) for face in minimal_nonfaces]
    if not supports or any(len(face) < 2 for face in supports):
        raise ValueError("minimal nonfaces must have order at least two")
    components: list[set[int]] = []
    for support in supports:
        touching = [index for index, component in enumerate(components) if component & support]
        merged = set(support)
        for index in reversed(touching):
            merged.update(components.pop(index))
        components.append(merged)
    return tuple(sorted(tuple(sorted(component)) for component in components))


def minimum_cross_component_nonright_edges(
    minimal_nonfaces: Sequence[Iterable[int]],
) -> int:
    """Necessary cross-component edge count for Lorentzian negative index one.

    Each minimal-nonface component contains an indefinite principal block.
    If the nonzero Gram graph does not connect all such components, its block
    decomposition has at least two negative directions.  A graph on ``c``
    initial components needs at least ``c-1`` added edges to become connected.
    """

    return len(minimal_nonface_component_partition(minimal_nonfaces)) - 1


def _cycle_edge_order(vertices: frozenset[int], edges: set[Pair]) -> tuple[Pair, ...]:
    adjacency = {
        vertex: tuple(
            sorted(other for other in vertices if other != vertex and pair(vertex, other) in edges)
        )
        for vertex in vertices
    }
    if len(vertices) != 5 or any(len(neighbors) != 2 for neighbors in adjacency.values()):
        raise ValueError("the polygon meeting graph must be a five-cycle")
    start = min(vertices)
    previous, current = start, min(adjacency[start])
    cycle_vertices = [start, current]
    while len(cycle_vertices) < 5:
        following = next(vertex for vertex in adjacency[current] if vertex != previous)
        if following == start:
            raise ValueError("the polygon meeting graph closed too early")
        cycle_vertices.append(following)
        previous, current = current, following
    if start not in adjacency[current]:
        raise ValueError("the polygon meeting graph is disconnected")
    return tuple(
        pair(cycle_vertices[index], cycle_vertices[(index + 1) % 5])
        for index in range(5)
    )


def iter_product_orthogonal_label_representatives(
    minimal_nonfaces: Sequence[Iterable[int]],
) -> Iterator[H914ProductLabelRepresentative]:
    """Enumerate the exact product-family symmetry representatives.

    This recognizes the minimal-nonface clutters of records 30 and 42.  All
    meeting pairs outside the named factor blocks receive label two.  The
    resulting iterator has length 45 or 6120 respectively.  It is a complete
    orbit enumeration *within this orthogonal-join ansatz*, not among all
    Coxeter labelings of the nerve.
    """

    nonfaces = tuple(tuple(sorted(set(face))) for face in minimal_nonfaces)
    by_order: dict[int, list[tuple[int, ...]]] = {}
    for face in nonfaces:
        by_order.setdefault(len(face), []).append(face)
    dotted = {pair(*face) for face in by_order.get(2, [])}
    meeting = set(combinations(range(14), 2)) - dotted

    if Counter(map(len, nonfaces)) == {2: 3, 4: 2}:
        four_blocks = tuple(sorted(by_order[4]))
        if len(set().union(*map(set, nonfaces))) != 14 or any(
            set(left) & set(right)
            for left, right in combinations(nonfaces, 2)
        ):
            raise ValueError("the record-30 product blocks must partition 14 facets")
        for first_type, second_type in combinations_with_replacement(range(9), 2):
            labels = {edge: 2 for edge in meeting}
            labels.update(representative_lanner_labels(four_blocks[0], first_type))
            labels.update(representative_lanner_labels(four_blocks[1], second_type))
            yield H914ProductLabelRepresentative(
                family="record30_Delta1^3_x_Delta3^2",
                lanner_type_indices=(first_type, second_type),
                cycle_color_code=(),
                labels=tuple(
                    (left, right, label)
                    for (left, right), label in sorted(labels.items())
                ),
            )
        return

    if Counter(map(len, nonfaces)) == {2: 5, 4: 1, 5: 1}:
        polygon_vertices = frozenset(vertex for edge in by_order[2] for vertex in edge)
        four_block = by_order[4][0]
        five_block = by_order[5][0]
        blocks = [polygon_vertices, set(four_block), set(five_block)]
        if set().union(*blocks) != set(range(14)) or any(
            set(left) & set(right) for left, right in combinations(blocks, 2)
        ):
            raise ValueError("the record-42 product blocks must partition 14 facets")
        polygon_meeting = {
            edge
            for edge in combinations(sorted(polygon_vertices), 2)
            if edge in meeting
        }
        cycle_edges = _cycle_edge_order(polygon_vertices, polygon_meeting)
        for four_type in range(9):
            for five_type in range(5):
                for colors in five_cycle_color_orbit_codes():
                    labels = {edge: 2 for edge in meeting}
                    labels.update(representative_lanner_labels(four_block, four_type))
                    labels.update(representative_lanner_labels(five_block, five_type))
                    labels.update(zip(cycle_edges, colors, strict=True))
                    yield H914ProductLabelRepresentative(
                        family="record42_C5_x_Delta3_x_Delta4",
                        lanner_type_indices=(four_type, five_type),
                        cycle_color_code=colors,
                        labels=tuple(
                            (left, right, label)
                            for (left, right), label in sorted(labels.items())
                        ),
                    )
        return

    raise ValueError("minimal nonfaces are neither product record 30 nor record 42")


def audit_h9_14_coxeter_labels(
    vertex_facets: Sequence[Iterable[int]],
    minimal_nonfaces: Sequence[Iterable[int]],
    labels: LabelMap,
) -> H914CoxeterLabelAudit:
    """Replay the complete local finite/Lanner label gate.

    ``vertex_facets`` must contain every rank-nine vertex set.  Minimal
    nonfaces of order two are treated as dotted pairs; those of orders four
    and five must carry one of the complete Lanner placements.
    """

    vertices = tuple(tuple(sorted(set(vertex))) for vertex in vertex_facets)
    if not vertices or any(len(vertex) != 9 for vertex in vertices):
        raise ValueError("the H9/14 gate requires the complete rank-nine vertex list")
    ground = tuple(sorted(set().union(*(set(vertex) for vertex in vertices))))
    if len(ground) != 14:
        raise ValueError("the H9/14 gate requires fourteen facets")

    nonfaces = tuple(tuple(sorted(set(face))) for face in minimal_nonfaces)
    if any(len(face) not in (2, 4, 5) for face in nonfaces):
        raise ValueError("Burcroff's H9/14 reduction permits only orders 2, 4, and 5")
    meeting_pairs = {
        pair(left, right)
        for vertex in vertices
        for left, right in combinations(vertex, 2)
    }
    dotted_pairs = {pair(*face) for face in nonfaces if len(face) == 2}
    all_pairs = {pair(left, right) for left, right in combinations(ground, 2)}
    if meeting_pairs | dotted_pairs != all_pairs or meeting_pairs & dotted_pairs:
        raise ValueError("the vertex list and dotted minimal nonfaces disagree on pairs")

    supplied_pairs = set(labels)
    missing_or_extra = tuple(sorted(meeting_pairs ^ supplied_pairs))
    bad_domain_pairs = {
        edge for edge in supplied_pairs & meeting_pairs if labels[edge] not in FINITE_LABELS
    }
    missing_or_extra = tuple(sorted(set(missing_or_extra) | bad_domain_pairs))
    bad_vertices = tuple(
        index
        for index, vertex in enumerate(vertices)
        if not is_finite_coxeter_labeling(vertex, labels)
    )
    bad_lanners = tuple(
        face
        for face in nonfaces
        if len(face) > 2 and not is_lanner_coxeter_labeling(face, labels)
    )
    return H914CoxeterLabelAudit(
        facet_count=14,
        vertex_count=len(vertices),
        meeting_pair_count=len(meeting_pairs),
        dotted_pair_count=len(dotted_pairs),
        finite_label_histogram=tuple(sorted(Counter(labels.values()).items())),
        minimal_nonface_order_histogram=tuple(
            sorted(Counter(map(len, nonfaces)).items())
        ),
        bad_vertex_indices=bad_vertices,
        bad_lanner_nonfaces=bad_lanners,
        missing_or_extra_label_pairs=missing_or_extra,
    )


__all__ = [
    "FINITE_LABELS",
    "H914CoxeterLabelAudit",
    "H914ProductLabelRepresentative",
    "NONRIGHT_LABELS",
    "Pair",
    "audit_h9_14_coxeter_labels",
    "five_cycle_color_orbit_codes",
    "induced_label_code",
    "iter_product_orthogonal_label_representatives",
    "is_finite_coxeter_labeling",
    "is_lanner_coxeter_labeling",
    "lanner_isomorphism_codes",
    "lanner_placement_codes",
    "minimal_nonface_component_partition",
    "minimum_cross_component_nonright_edges",
    "pair",
    "representative_lanner_labels",
]
