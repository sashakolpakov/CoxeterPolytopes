"""Exact sparse Coxeter-label census for H9/14 records 43, 44, and 47.

The three ordinary-polytopal multiwedge nerves each have four minimal
nonfaces of order five.  This module exhausts the following deliberately
restricted, but useful, discrete ansatz:

* every order-five minimal nonface carries one of the 420 labelled placements
  of the five classical order-five Lanner diagrams with labels in
  ``{2, 3, 4, 5}``;
* placements agree on every pair belonging to two missing five-faces; and
* every meeting pair outside those four supports has label two.

The first missing five-face is the ``boundary(Delta_4)`` join factor in each
multiwedge construction.  It is disjoint from the other three missing
five-faces.  The 90 vertices consequently factor as its five four-node
deletions times 18 residual vertices.  This makes the exhaustive calculation
small: join only the three overlapping placement tables, test all 18 residual
vertices, and multiply by the 420 independent placements on the join factor.
The Cartesian-product identity is checked against all 90 archived vertices,
and one final representative is replayed through the general 90-vertex
label verifier.

This factorization also gives a global no-go for the sparse ansatz.  All
cross-factor pairs have label two, so the normalized Gram matrix splits off
the isolated Lanner block, which has one negative direction.  The residual
block contains another Lanner principal submatrix and hence has at least one
negative direction by interlacing.  Thus the full negative index is at least
two, not the Lorentzian value one.  The census supplies exact local label
witnesses, not hyperbolic Coxeter examples.
"""

from __future__ import annotations

from collections import Counter, deque
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    audit_h9_14_coxeter_labels,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    minimal_nonface_component_partition,
    minimum_cross_component_nonright_edges,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_multiwedge import (
    RECORD43,
    RECORD44,
    RECORD47,
)
from coxeter_polytopes.combinatorics.h9_14_record43 import (
    record43_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record44 import (
    record44_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record47 import (
    record47_facet_complements,
)


FACET_COUNT = 14
FULL_MASK = (1 << FACET_COUNT) - 1
ALL_PAIRS = tuple(combinations(range(FACET_COUNT), 2))
LOCAL_FIVE_PAIRS = tuple(combinations(range(5), 2))

Permutation = tuple[int, ...]
LabelState = tuple[tuple[Pair, int], ...]


def _permutation(*cycles: tuple[int, ...]) -> Permutation:
    image = list(range(FACET_COUNT))
    for cycle in cycles:
        if len(cycle) < 2:
            raise ValueError("a permutation cycle must have at least two entries")
        for source, target in zip(cycle, cycle[1:] + cycle[:1], strict=True):
            image[source] = target
    if len(set(image)) != FACET_COUNT:
        raise ValueError("cycles do not define a permutation")
    return tuple(image)


def _adjacent_transpositions(vertices: tuple[int, ...]) -> tuple[Permutation, ...]:
    return tuple(
        _permutation((vertices[index], vertices[index + 1]))
        for index in range(len(vertices) - 1)
    )


# These generators act on the full minimal-nonface hypergraph, not merely on
# its missing-edge graph.  The latter has spurious symmetric isolated nodes.
_AUTOMORPHISM_GENERATORS: dict[int, tuple[Permutation, ...]] = {
    43: (
        *_adjacent_transpositions((3, 4, 5, 6, 7)),
        *_adjacent_transpositions((8, 9, 10)),
        _permutation((0, 13), (1, 12), (2, 11)),
    ),
    44: (
        *_adjacent_transpositions((2, 3, 4, 5, 6)),
        *_adjacent_transpositions((7, 8, 9, 10)),
        _permutation((0, 13), (1, 12)),
    ),
    47: (
        *_adjacent_transpositions((3, 4, 5, 6, 7)),
        *_adjacent_transpositions((8, 9, 10)),
        _permutation((0, 12), (11, 1)),
        _permutation((12, 13), (1, 2)),
    ),
}


@dataclass(frozen=True, slots=True)
class _SparseSpecification:
    record_index: int
    graph6: str
    missing_edges: tuple[tuple[int, int], ...]
    missing_fives: tuple[tuple[int, ...], ...]
    simplex_join_vertices: tuple[int, ...]
    facet_complements: Callable[[], tuple[int, ...]]


_SPECIFICATIONS = {
    43: _SparseSpecification(
        record_index=43,
        graph6=RECORD43.graph6,
        missing_edges=RECORD43.missing_edges,
        missing_fives=RECORD43.missing_fives,
        simplex_join_vertices=RECORD43.simplex_join_vertices,
        facet_complements=record43_facet_complements,
    ),
    44: _SparseSpecification(
        record_index=44,
        graph6=RECORD44.graph6,
        missing_edges=RECORD44.missing_edges,
        missing_fives=RECORD44.missing_fives,
        simplex_join_vertices=RECORD44.simplex_join_vertices,
        facet_complements=record44_facet_complements,
    ),
    47: _SparseSpecification(
        record_index=47,
        graph6=RECORD47.graph6,
        missing_edges=RECORD47.missing_edges,
        missing_fives=RECORD47.missing_fives,
        simplex_join_vertices=RECORD47.simplex_join_vertices,
        facet_complements=record47_facet_complements,
    ),
}


_EXPECTED = {
    43: {
        "cluster_join_counts": (420, 1_332, 4_212),
        "passing_cluster_count": 30,
        "full_labelled_count": 12_600,
        "automorphism_group_order": 1_440,
        "orbit_count": 15,
        "orbit_size_histogram": ((360, 3), (720, 8), (1_440, 4)),
    },
    44: {
        "cluster_join_counts": (420, 1_332, 4_620),
        "passing_cluster_count": 144,
        "full_labelled_count": 60_480,
        "automorphism_group_order": 5_760,
        "orbit_count": 25,
        "orbit_size_histogram": ((1_440, 12), (2_880, 11), (5_760, 2)),
    },
    47: {
        "cluster_join_counts": (420, 1_332, 222),
        "passing_cluster_count": 18,
        "full_labelled_count": 7_560,
        "automorphism_group_order": 4_320,
        "orbit_count": 5,
        "orbit_size_histogram": ((1_080, 3), (2_160, 2)),
    },
}


def _vertices(mask: int) -> tuple[int, ...]:
    if mask < 0 or mask & ~FULL_MASK:
        raise ValueError("facet mask lies outside the H9/14 ground set")
    return tuple(vertex for vertex in range(FACET_COUNT) if mask & (1 << vertex))


@lru_cache(maxsize=None)
def _placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    if len(support) != 5 or len(set(support)) != 5:
        raise ValueError("an order-five Lanner support must contain five facets")
    support_pairs = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(support_pairs, code, strict=True))
        for code in lanner_placement_codes(5)
    )


def _merge_states(left: LabelState, right: LabelState) -> LabelState | None:
    merged = dict(left)
    for edge, label in right:
        old = merged.get(edge)
        if old is not None and old != label:
            return None
        merged[edge] = label
    return tuple(sorted(merged.items()))


def _join_placement_tables(
    supports: Sequence[tuple[int, ...]],
) -> tuple[tuple[LabelState, ...], tuple[int, ...]]:
    """Join placement tables exactly, retaining one state per union labelling.

    Restriction of a union labelling to any support recovers its unique table
    row.  Deduplicating union states therefore loses no placement tuple.
    """

    states: tuple[LabelState, ...] = ((),)
    layer_counts: list[int] = []
    for support in supports:
        next_states = {
            merged
            for state in states
            for option in _placement_options(tuple(support))
            if (merged := _merge_states(state, option)) is not None
        }
        states = tuple(sorted(next_states))
        layer_counts.append(len(states))
    return states, tuple(layer_counts)


@lru_cache(maxsize=None)
def _finite_code(code: tuple[int, ...], order: int) -> bool:
    local_labels = {
        edge: label
        for edge, label in zip(combinations(range(order), 2), code, strict=True)
    }
    return is_finite_coxeter_labeling(range(order), local_labels)


def _induced_code(facets: Sequence[int], labels: Mapping[Pair, int]) -> tuple[int, ...]:
    vertices = tuple(sorted(facets))
    return tuple(labels[pair(left, right)] for left, right in combinations(vertices, 2))


def _all_deletions_are_finite(
    support: tuple[int, ...], option: LabelState
) -> bool:
    labels = dict(option)
    return all(
        _finite_code(
            _induced_code(
                tuple(vertex for vertex in support if vertex != deleted), labels
            ),
            len(support) - 1,
        )
        for deleted in support
    )


def _compose(left: Permutation, right: Permutation) -> Permutation:
    return tuple(left[right[index]] for index in range(FACET_COUNT))


def _generated_group_order(generators: Sequence[Permutation]) -> int:
    identity = tuple(range(FACET_COUNT))
    group = {identity}
    frontier = deque((identity,))
    while frontier:
        current = frontier.popleft()
        for generator in generators:
            image = _compose(generator, current)
            if image not in group:
                group.add(image)
                frontier.append(image)
    return len(group)


def _validate_automorphisms(
    generators: Sequence[Permutation],
    vertices: Sequence[tuple[int, ...]],
    nonfaces: Sequence[tuple[int, ...]],
) -> None:
    vertex_set = {frozenset(vertex) for vertex in vertices}
    nonface_set = {frozenset(face) for face in nonfaces}
    for generator in generators:
        if {
            frozenset(generator[vertex] for vertex in face) for face in nonfaces
        } != nonface_set:
            raise AssertionError("a declared symmetry does not preserve minimal nonfaces")
        if {
            frozenset(generator[vertex] for vertex in face) for face in vertices
        } != vertex_set:
            raise AssertionError("a declared symmetry does not preserve all vertices")


def _code_action_maps(
    meeting_pairs: tuple[Pair, ...], generators: Sequence[Permutation]
) -> tuple[tuple[int, ...], ...]:
    position = {edge: index for index, edge in enumerate(meeting_pairs)}
    result = []
    for generator in generators:
        result.append(
            tuple(
                position[pair(generator[left], generator[right])]
                for left, right in meeting_pairs
            )
        )
    return tuple(result)


def _act_on_code(code: tuple[int, ...], action: tuple[int, ...]) -> tuple[int, ...]:
    image = [0] * len(code)
    for source, target in enumerate(action):
        image[target] = code[source]
    return tuple(image)


def _orbit_census(
    codes: set[tuple[int, ...]], actions: Sequence[tuple[int, ...]]
) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, int], ...]]:
    remaining = set(codes)
    sizes: Counter[int] = Counter()
    representatives: list[tuple[int, ...]] = []
    while remaining:
        seed = min(remaining)
        orbit = {seed}
        frontier = deque((seed,))
        while frontier:
            current = frontier.popleft()
            for action in actions:
                image = _act_on_code(current, action)
                if image not in codes:
                    raise AssertionError("the sparse family is not symmetry invariant")
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        sizes[len(orbit)] += 1
        representatives.append(seed)
        remaining.difference_update(orbit)
    return tuple(representatives), tuple(sorted(sizes.items()))


@dataclass(frozen=True, slots=True)
class H914SparseLabelAudit:
    """Complete exact census for one sparse multiwedge label ansatz."""

    record_index: int
    graph6: str
    isolated_simplex_lanner_face: tuple[int, ...]
    overlapping_lanner_faces: tuple[tuple[int, ...], ...]
    vertex_count: int
    residual_vertex_count: int
    meeting_pair_count: int
    dotted_pair_count: int
    lanner_placements_per_five_face: int
    cluster_join_counts: tuple[int, ...]
    compatible_four_lanner_assignment_count: int
    passing_cluster_assignment_count: int
    full_labelled_assignment_count: int
    nerve_automorphism_group_order: int
    symmetry_orbit_count: int
    symmetry_orbit_size_histogram: tuple[tuple[int, int], ...]
    symmetry_orbit_representative_codes: tuple[tuple[int, ...], ...]
    representative_labels: tuple[tuple[Pair, int], ...]
    representative_nonright_labels: tuple[tuple[Pair, int], ...]
    representative_finite_label_histogram: tuple[tuple[int, int], ...]
    representative_all_ninety_vertices_finite: bool
    representative_all_four_missing_fives_lanner: bool
    orthogonal_negative_index_lower_bound: int
    minimal_nonface_components: tuple[tuple[int, ...], ...]
    unrestricted_minimum_cross_component_nonright_edges: int

    @property
    def representative_label_map(self) -> dict[Pair, int]:
        return dict(self.representative_labels)

    def orbit_representative_label_maps(self) -> tuple[dict[Pair, int], ...]:
        """Return every sparse symmetry representative for later Gram gates."""

        meeting_pairs = tuple(edge for edge, _label in self.representative_labels)
        return tuple(
            dict(zip(meeting_pairs, code, strict=True))
            for code in self.symmetry_orbit_representative_codes
        )

    @property
    def verified(self) -> bool:
        expected = _EXPECTED[self.record_index]
        return (
            self.vertex_count == 90
            and self.residual_vertex_count == 18
            and self.meeting_pair_count == 85
            and self.dotted_pair_count == 6
            and self.lanner_placements_per_five_face == 420
            and self.cluster_join_counts == expected["cluster_join_counts"]
            and self.compatible_four_lanner_assignment_count
            == 420 * self.cluster_join_counts[-1]
            and self.passing_cluster_assignment_count
            == expected["passing_cluster_count"]
            and self.full_labelled_assignment_count
            == expected["full_labelled_count"]
            and self.full_labelled_assignment_count
            == 420 * self.passing_cluster_assignment_count
            and self.nerve_automorphism_group_order
            == expected["automorphism_group_order"]
            and self.symmetry_orbit_count == expected["orbit_count"]
            and len(self.symmetry_orbit_representative_codes)
            == self.symmetry_orbit_count
            and self.symmetry_orbit_size_histogram
            == expected["orbit_size_histogram"]
            and self.representative_all_ninety_vertices_finite
            and self.representative_all_four_missing_fives_lanner
            and self.orthogonal_negative_index_lower_bound == 2
            and len(self.minimal_nonface_components) == 2
            and self.unrestricted_minimum_cross_component_nonright_edges == 1
        )

    def to_record(self) -> dict[str, object]:
        return {
            "record_index": self.record_index,
            "graph6": self.graph6,
            "sparse_ansatz": {
                "lanner_supports": [
                    list(self.isolated_simplex_lanner_face),
                    *[list(face) for face in self.overlapping_lanner_faces],
                ],
                "default_label_on_every_other_meeting_pair": 2,
                "isolated_simplex_lanner_face": list(
                    self.isolated_simplex_lanner_face
                ),
                "isolated_factor": "boundary of Delta_4",
            },
            "vertex_count": self.vertex_count,
            "residual_vertex_count": self.residual_vertex_count,
            "meeting_pair_count": self.meeting_pair_count,
            "dotted_pair_count": self.dotted_pair_count,
            "lanner_placements_per_five_face": (
                self.lanner_placements_per_five_face
            ),
            "cluster_join_counts": list(self.cluster_join_counts),
            "compatible_four_lanner_assignment_count": (
                self.compatible_four_lanner_assignment_count
            ),
            "passing_cluster_assignment_count": (
                self.passing_cluster_assignment_count
            ),
            "full_labelled_assignment_count": self.full_labelled_assignment_count,
            "nerve_automorphism_group_order": (
                self.nerve_automorphism_group_order
            ),
            "symmetry_orbit_count": self.symmetry_orbit_count,
            "symmetry_orbit_size_histogram": {
                str(size): multiplicity
                for size, multiplicity in self.symmetry_orbit_size_histogram
            },
            "symmetry_orbit_representative_nonright_labels": [
                [
                    [left, right, label]
                    for (left, right), label in zip(
                        (edge for edge, _ in self.representative_labels),
                        code,
                        strict=True,
                    )
                    if label > 2
                ]
                for code in self.symmetry_orbit_representative_codes
            ],
            "representative": {
                "labels": [
                    [left, right, label]
                    for (left, right), label in self.representative_labels
                ],
                "nonright_labels": [
                    [left, right, label]
                    for (left, right), label in self.representative_nonright_labels
                ],
                "finite_label_histogram": {
                    str(label): count
                    for label, count in self.representative_finite_label_histogram
                },
                "all_ninety_vertices_finite": (
                    self.representative_all_ninety_vertices_finite
                ),
                "all_four_missing_fives_lanner": (
                    self.representative_all_four_missing_fives_lanner
                ),
            },
            "global_sparse_no_go": {
                "orthogonal_negative_index_lower_bound": (
                    self.orthogonal_negative_index_lower_bound
                ),
                "reason": (
                    "the isolated simplex-factor Lanner block contributes one "
                    "negative direction and the residual block contains another "
                    "Lanner principal submatrix"
                ),
                "rules_out_lorentzian_negative_index_one": True,
            },
            "unrestricted_forced_cross_edge_gate": {
                "minimal_nonface_components": [
                    list(component) for component in self.minimal_nonface_components
                ],
                "component_count": len(self.minimal_nonface_components),
                "minimum_nonright_meeting_edges_between_components": (
                    self.unrestricted_minimum_cross_component_nonright_edges
                ),
                "scope": (
                    "necessary for every global labeling, not only the sparse ansatz"
                ),
            },
            "verified": self.verified,
        }


@lru_cache(maxsize=3)
def audit_sparse_h9_14_labels(record_index: int) -> H914SparseLabelAudit:
    """Exhaust and symmetry-quotient the sparse ansatz for one record."""

    try:
        specification = _SPECIFICATIONS[record_index]
    except KeyError as error:
        raise ValueError("the sparse census covers records 43, 44, and 47") from error

    isolated = tuple(specification.simplex_join_vertices)
    if isolated != specification.missing_fives[0]:
        raise AssertionError("the first A5 is no longer the simplex join factor")
    cluster = tuple(specification.missing_fives[1:])
    if any(set(isolated) & set(face) for face in cluster):
        raise AssertionError("the simplex-factor A5 is no longer isolated")

    nonfaces = tuple(specification.missing_edges) + tuple(
        specification.missing_fives
    )
    nonface_components = minimal_nonface_component_partition(nonfaces)
    minimum_cross_edges = minimum_cross_component_nonright_edges(nonfaces)
    if len(nonface_components) != 2 or minimum_cross_edges != 1:
        raise AssertionError("the multiwedge minimal-nonface components changed")
    dotted = {pair(*edge) for edge in specification.missing_edges}
    meeting_pairs = tuple(edge for edge in ALL_PAIRS if edge not in dotted)
    facet_complements = tuple(specification.facet_complements())
    vertices = tuple(
        sorted(_vertices(FULL_MASK ^ complement) for complement in facet_complements)
    )
    if len(vertices) != 90:
        raise AssertionError("the archived multiwedge no longer has 90 vertices")

    isolated_deletions = tuple(
        tuple(vertex for vertex in isolated if vertex != deleted)
        for deleted in isolated
    )
    residual_vertices = tuple(
        sorted({
            tuple(vertex for vertex in maximal if vertex not in isolated)
            for maximal in vertices
        })
    )
    factored_vertices = {
        tuple(sorted((*left, *right)))
        for left in isolated_deletions
        for right in residual_vertices
    }
    if len(residual_vertices) != 18 or factored_vertices != set(vertices):
        raise AssertionError("the 5 by 18 vertex factorization did not reproduce")

    isolated_options = _placement_options(isolated)
    if len(isolated_options) != 420 or not all(
        _all_deletions_are_finite(isolated, option)
        for option in isolated_options
    ):
        raise AssertionError("an isolated Lanner placement has a bad deletion")

    cluster_states, join_counts = _join_placement_tables(cluster)
    passing_cluster_states: list[LabelState] = []
    for state in cluster_states:
        state_labels = {edge: 2 for edge in meeting_pairs}
        state_labels.update(state)
        if all(
            _finite_code(_induced_code(vertex, state_labels), len(vertex))
            for vertex in residual_vertices
        ):
            passing_cluster_states.append(state)
    passing_cluster_states = sorted(passing_cluster_states)

    full_codes: set[tuple[int, ...]] = set()
    for isolated_state in isolated_options:
        for cluster_state in passing_cluster_states:
            labels = {edge: 2 for edge in meeting_pairs}
            labels.update(isolated_state)
            labels.update(cluster_state)
            full_codes.add(tuple(labels[edge] for edge in meeting_pairs))
    expected_full_count = len(isolated_options) * len(passing_cluster_states)
    if len(full_codes) != expected_full_count:
        raise AssertionError("two sparse placement tuples produced one labelling")
    representative_code = min(full_codes)
    representative_labels = tuple(zip(meeting_pairs, representative_code, strict=True))
    representative_map = dict(representative_labels)
    representative_audit = audit_h9_14_coxeter_labels(
        vertices, nonfaces, representative_map
    )
    if not representative_audit.passes_discrete_gate:
        raise AssertionError("the deterministic sparse representative failed replay")

    generators = _AUTOMORPHISM_GENERATORS[record_index]
    _validate_automorphisms(generators, vertices, nonfaces)
    group_order = _generated_group_order(generators)
    actions = _code_action_maps(meeting_pairs, generators)
    orbit_representatives, orbit_histogram = _orbit_census(full_codes, actions)
    orbit_count = sum(multiplicity for _, multiplicity in orbit_histogram)

    result = H914SparseLabelAudit(
        record_index=record_index,
        graph6=specification.graph6,
        isolated_simplex_lanner_face=isolated,
        overlapping_lanner_faces=cluster,
        vertex_count=len(vertices),
        residual_vertex_count=len(residual_vertices),
        meeting_pair_count=len(meeting_pairs),
        dotted_pair_count=len(dotted),
        lanner_placements_per_five_face=len(isolated_options),
        cluster_join_counts=join_counts,
        compatible_four_lanner_assignment_count=(
            len(isolated_options) * len(cluster_states)
        ),
        passing_cluster_assignment_count=len(passing_cluster_states),
        full_labelled_assignment_count=len(full_codes),
        nerve_automorphism_group_order=group_order,
        symmetry_orbit_count=orbit_count,
        symmetry_orbit_size_histogram=orbit_histogram,
        symmetry_orbit_representative_codes=orbit_representatives,
        representative_labels=representative_labels,
        representative_nonright_labels=tuple(
            (edge, label) for edge, label in representative_labels if label > 2
        ),
        representative_finite_label_histogram=tuple(
            sorted(Counter(representative_code).items())
        ),
        representative_all_ninety_vertices_finite=(
            not representative_audit.bad_vertex_indices
        ),
        representative_all_four_missing_fives_lanner=(
            not representative_audit.bad_lanner_nonfaces
        ),
        orthogonal_negative_index_lower_bound=2,
        minimal_nonface_components=nonface_components,
        unrestricted_minimum_cross_component_nonright_edges=minimum_cross_edges,
    )
    if not result.verified:
        raise AssertionError(f"record-{record_index} sparse census changed")
    return result


def audit_all_sparse_h9_14_labels() -> tuple[H914SparseLabelAudit, ...]:
    """Return the deterministic sparse censuses for records 43, 44, and 47."""

    return tuple(audit_sparse_h9_14_labels(index) for index in (43, 44, 47))


__all__ = [
    "H914SparseLabelAudit",
    "audit_all_sparse_h9_14_labels",
    "audit_sparse_h9_14_labels",
]
