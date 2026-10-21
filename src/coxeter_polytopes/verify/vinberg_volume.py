"""Vinberg's exact diagram test for compactness and finite covolume.

This module implements the edge-endpoint criterion in Vinberg's Proposition
4.2 (in the diagram form recorded as Propositions 1.9--1.10 of CoxIter).
It is deliberately gated on a *complete* finite list of simple roots because
its return value describes the chamber represented by the supplied diagram,
including its final vertices and cusps.  For an incomplete but genuinely
oriented wall prefix, failure of the criterion is inconclusive: later walls
may cut off an infinite-volume escape.  Success does give a one-sided finite-
volume outer bound by monotonicity, and a method-valid Vinberg prefix may then
be promoted to a completeness certificate by the separate stopping theorem.

For a Coxeter polyhedron in hyperbolic dimension ``n``:

* elliptic subdiagrams on ``n - 1`` nodes are its edges;
* elliptic subdiagrams on ``n`` nodes are its finite vertices; and
* componentwise affine parabolic subdiagrams of Coxeter rank ``n - 1``
  are its ideal-vertex neighbourhoods.

Every edge must have exactly two endpoint extensions.  Compactness permits
only elliptic endpoint extensions; finite covolume also permits parabolic
ones.  A reducible parabolic neighbourhood is counted once, not once per
added node: if it has ``c`` affine components, it contains ``n - 1 + c``
nodes and an incident elliptic edge is obtained by deleting one node from
each component.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from typing import Literal

from coxeter_polytopes.exact.diagram import (
    connected_components,
    induced_subgram,
    is_affine_parabolic_component,
    is_parabolic,
    parabolic_rank,
)
from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import Inertia, Matrix, inertia


Decision = Literal[
    "compact",
    "finite_volume_noncompact",
    "infinite_volume",
    "undecidable_incomplete_simple_roots",
]
EndpointKind = Literal["finite", "ideal"]


@dataclass(frozen=True, slots=True)
class VinbergEndpoint:
    """One finite vertex or full ideal-vertex parabolic neighbourhood."""

    kind: EndpointKind
    facets: tuple[int, ...]
    affine_components: tuple[tuple[int, ...], ...] = ()


@dataclass(frozen=True, slots=True)
class VinbergEdgeAudit:
    """The endpoint extensions of one elliptic rank-``n-1`` edge diagram."""

    facets: tuple[int, ...]
    finite_endpoints: tuple[tuple[int, ...], ...]
    ideal_endpoints: tuple[tuple[int, ...], ...]

    @property
    def endpoint_count(self) -> int:
        return len(self.finite_endpoints) + len(self.ideal_endpoints)


@dataclass(frozen=True, slots=True)
class VinbergParabolicNeighborhoodAudit:
    """Vinberg's full zero-neighbourhood test for one connected affine diagram."""

    affine_component: tuple[int, ...]
    zero_neighborhood: tuple[int, ...]
    neighborhood_components: tuple[tuple[int, ...], ...]
    is_rank_dimension_minus_one_parabolic: bool


@dataclass(frozen=True, slots=True)
class VinbergVolumeVerification:
    """Exact, reproducible output of Vinberg's volume criterion."""

    dimension: int
    facet_count: int
    simple_roots_complete: bool
    decision: Decision
    finite_volume: bool | None
    compact: bool | None
    finite_vertices: tuple[tuple[int, ...], ...]
    ideal_vertices: tuple[VinbergEndpoint, ...]
    parabolic_neighborhood_audits: tuple[VinbergParabolicNeighborhoodAudit, ...]
    edge_audits: tuple[VinbergEdgeAudit, ...]
    compact_bad_edges: tuple[tuple[int, ...], ...]
    finite_volume_bad_edges: tuple[tuple[int, ...], ...]

    @property
    def elliptic_edge_count(self) -> int:
        return len(self.edge_audits)

    @property
    def finite_vertex_count(self) -> int:
        return len(self.finite_vertices)

    @property
    def ideal_vertex_count(self) -> int:
        return len(self.ideal_vertices)


def _coerce_wall_gram(wall_gram: Matrix) -> list[list[Golden]]:
    """Coerce and validate the global Coxeter sign convention."""

    size = len(wall_gram)
    matrix = induced_subgram(wall_gram, range(size))
    if any(matrix[index][index].sign() <= 0 for index in range(size)):
        raise ValueError("wall normals must have positive square")
    if any(
        matrix[left][right].sign() > 0
        for left in range(size)
        for right in range(left)
    ):
        raise ValueError("an acute Coxeter wall Gram matrix has nonpositive cross terms")
    return matrix


def _elliptic_subdiagrams(
    matrix: list[list[Golden]], order: int
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        vertices
        for vertices in combinations(range(len(matrix)), order)
        if _is_positive_definite(_subgram(matrix, vertices))
    )


def _subgram(
    matrix: list[list[Golden]], vertices: tuple[int, ...]
) -> list[list[Golden]]:
    """Fast principal submatrix after the global Gram has been validated."""

    return [[matrix[left][right] for right in vertices] for left in vertices]


def _is_positive_definite(matrix: list[list[Golden]]) -> bool:
    """Exact Sylvester test, with a fraction-free rational fast path."""

    simply_laced = _finite_simply_laced_cartan(matrix)
    if simply_laced is not None:
        return simply_laced
    if any(
        (
            matrix[left][left] * matrix[right][right]
            - matrix[left][right] * matrix[left][right]
        ).sign()
        <= 0
        for left in range(len(matrix))
        for right in range(left)
    ):
        return False
    if all(entry.is_rational for row in matrix for entry in row):
        work = [[entry.a for entry in row] for row in matrix]
        previous_pivot = Fraction(1)
        for pivot_index in range(len(work)):
            pivot = work[pivot_index][pivot_index]
            if pivot <= 0:
                return False
            if pivot_index == len(work) - 1:
                return True
            for row in range(pivot_index + 1, len(work)):
                for column in range(row, len(work)):
                    value = (
                        work[row][column] * pivot
                        - work[row][pivot_index] * work[pivot_index][column]
                    ) / previous_pivot
                    work[row][column] = value
                    work[column][row] = value
            previous_pivot = pivot
        return True
    return inertia(matrix) == Inertia(len(matrix), 0, 0)


def _finite_simply_laced_cartan(matrix: list[list[Golden]]) -> bool | None:
    """Recognize finite ADE diagrams, or decline for non-Cartan input.

    Complete K3 curve diagrams are large but have diagonal two and integral
    cross terms.  Any cross term at most ``-2`` already makes its two-node
    principal minor nonpositive.  Otherwise the matrix is ``2I-A`` for a
    simple graph, and positivity is exactly the finite ADE tree condition.
    Returning ``None`` selects the general exact-arithmetic fallback.
    """

    size = len(matrix)
    if any(matrix[index][index] != 2 for index in range(size)):
        return None
    if any(
        not matrix[left][right].is_rational
        or matrix[left][right].a.denominator != 1
        for left in range(size)
        for right in range(left)
    ):
        return None
    if any(
        matrix[left][right].a <= -2
        for left in range(size)
        for right in range(left)
    ):
        return False
    if any(
        matrix[left][right].a not in (0, -1)
        for left in range(size)
        for right in range(left)
    ):
        return None

    adjacency = [set() for _ in range(size)]
    for left in range(size):
        for right in range(left):
            if matrix[left][right] == -1:
                adjacency[left].add(right)
                adjacency[right].add(left)

    unseen = set(range(size))
    while unseen:
        start = min(unseen)
        component = {start}
        frontier = [start]
        unseen.remove(start)
        while frontier:
            vertex = frontier.pop()
            new = adjacency[vertex] & unseen
            unseen.difference_update(new)
            component.update(new)
            frontier.extend(new)

        degree_sum = sum(len(adjacency[vertex] & component) for vertex in component)
        if degree_sum != 2 * (len(component) - 1):
            return False
        degrees = {
            vertex: len(adjacency[vertex] & component) for vertex in component
        }
        if max(degrees.values(), default=0) <= 2:
            continue  # A_r (including A_1).
        branch_vertices = [vertex for vertex, degree in degrees.items() if degree == 3]
        if len(branch_vertices) != 1 or any(degree > 3 for degree in degrees.values()):
            return False
        branch = branch_vertices[0]
        arm_lengths: list[int] = []
        for neighbor in adjacency[branch] & component:
            length = 1
            previous, current = branch, neighbor
            while degrees[current] == 2:
                following = next(
                    vertex
                    for vertex in adjacency[current] & component
                    if vertex != previous
                )
                previous, current = current, following
                length += 1
            arm_lengths.append(length)
        arms = tuple(sorted(arm_lengths))
        if not (
            (arms[0:2] == (1, 1))
            or arms in {(1, 2, 2), (1, 2, 3), (1, 2, 4)}
        ):
            return False
    return True


def _connected_affine_subdiagrams(
    matrix: list[list[Golden]], dimension: int
) -> tuple[tuple[int, ...], ...]:
    """Enumerate connected affine diagrams, pruning non-elliptic branches.

    Every proper subdiagram of an irreducible affine Coxeter diagram is
    elliptic.  We therefore grow connected elliptic subsets along diagram
    edges and stop a branch as soon as it becomes affine or indefinite.  This
    is substantially smaller than scanning the full power set on sparse
    high-dimensional diagrams.
    """

    affine: list[tuple[int, ...]] = []
    candidates = {(vertex,) for vertex in range(len(matrix))}
    # Invalid affine diagrams can themselves exceed the cusp rank.  They must
    # be found and rejected rather than excluded by the desired conclusion.
    # The ambient Lorentzian Gram rank bounds a connected affine diagram to
    # order at most n+2; the extra two layers are cheap after elliptic pruning.
    maximum_order = min(len(matrix), dimension + 2)
    for _order in range(1, maximum_order + 1):
        next_candidates: set[tuple[int, ...]] = set()
        for facets in sorted(candidates):
            subgram = _subgram(matrix, facets)
            if len(facets) >= 2 and is_affine_parabolic_component(subgram):
                affine.append(facets)
                continue
            if not _is_positive_definite(subgram):
                continue
            facet_set = set(facets)
            adjacent = {
                vertex
                for node in facets
                for vertex in range(len(matrix))
                if vertex not in facet_set and matrix[node][vertex] != 0
            }
            next_candidates.update(
                tuple(sorted((*facets, vertex))) for vertex in adjacent
            )
        candidates = next_candidates
        if not candidates:
            break
    return tuple(affine)


def _parabolic_vertex_neighbourhoods(
    matrix: list[list[Golden]], dimension: int
) -> tuple[
    tuple[VinbergEndpoint, ...],
    tuple[VinbergParabolicNeighborhoodAudit, ...],
]:
    """Apply Vinberg's full-neighbourhood test at every affine component.

    If ``J`` is connected affine, let ``Z(J)`` be all diagram nodes outside
    ``J`` orthogonal to every node of ``J``.  The facets through the putative
    cusp are exactly ``N(J) = J union Z(J)``.  Finite covolume requires this
    *full* neighbourhood, not an arbitrary parabolic subset, to be a
    componentwise affine diagram of rank ``n - 1``.  If it has ``c``
    components this is equivalently ``|N(J)| = n - 1 + c``.
    """

    target_rank = dimension - 1
    audits: list[VinbergParabolicNeighborhoodAudit] = []
    endpoints_by_facets: dict[tuple[int, ...], VinbergEndpoint] = {}
    for affine_component in _connected_affine_subdiagrams(matrix, dimension):
        component_set = set(affine_component)
        zero_nodes = tuple(
            vertex
            for vertex in range(len(matrix))
            if vertex not in component_set
            and all(matrix[vertex][node] == 0 for node in affine_component)
        )
        neighborhood = tuple(sorted((*affine_component, *zero_nodes)))
        neighborhood_gram = _subgram(matrix, neighborhood)
        local_components = connected_components(neighborhood_gram)
        components = tuple(
            tuple(neighborhood[local] for local in component)
            for component in local_components
        )
        is_full_parabolic = (
            len(neighborhood) == target_rank + len(components)
            and is_parabolic(neighborhood_gram)
            and parabolic_rank(neighborhood_gram) == target_rank
        )
        audits.append(
            VinbergParabolicNeighborhoodAudit(
                affine_component,
                neighborhood,
                components,
                is_full_parabolic,
            )
        )
        if is_full_parabolic:
            endpoints_by_facets.setdefault(
                neighborhood,
                VinbergEndpoint("ideal", neighborhood, components),
            )
    return tuple(endpoints_by_facets.values()), tuple(audits)


def _is_incident_ideal_extension(
    edge: tuple[int, ...], endpoint: VinbergEndpoint
) -> bool:
    """Whether ``edge`` reaches the full parabolic endpoint neighbourhood.

    At a reducible cusp, containment alone is supplemented by the exact local
    incidence rule: one node must be omitted from every affine component.
    Thus the entire parabolic neighbourhood is one extension of the edge.
    """

    edge_set = set(edge)
    if not edge_set.issubset(endpoint.facets):
        return False
    return all(len(set(component) - edge_set) == 1 for component in endpoint.affine_components)


def verify_vinberg_volume_criterion(
    wall_gram: Matrix,
    dimension: int,
    *,
    simple_roots_complete: bool,
) -> VinbergVolumeVerification:
    """Apply Vinberg's exact compactness/finite-covolume diagram criterion.

    ``simple_roots_complete`` is a mandatory hard gate.  When it is false the
    function returns an undecidable result without enumerating the prefix's
    apparent vertices or edges.  When it is true, ``wall_gram`` must be the
    Gram matrix of all irredundant simple walls of an acute Coxeter chamber in
    hyperbolic dimension ``dimension``.
    """

    if not isinstance(dimension, int) or dimension < 2:
        raise ValueError("hyperbolic dimension must be an integer at least two")
    if not isinstance(simple_roots_complete, bool):
        raise TypeError("simple_roots_complete must be a bool")

    matrix = _coerce_wall_gram(wall_gram)
    facet_count = len(matrix)
    if not simple_roots_complete:
        return VinbergVolumeVerification(
            dimension=dimension,
            facet_count=facet_count,
            simple_roots_complete=False,
            decision="undecidable_incomplete_simple_roots",
            finite_volume=None,
            compact=None,
            finite_vertices=(),
            ideal_vertices=(),
            parabolic_neighborhood_audits=(),
            edge_audits=(),
            compact_bad_edges=(),
            finite_volume_bad_edges=(),
        )

    if facet_count and len(connected_components(matrix)) != 1:
        raise ValueError("a hyperbolic Coxeter wall Gram matrix must be indecomposable")

    expected_inertia = Inertia(dimension, 1, facet_count - dimension - 1)
    if expected_inertia.zero < 0 or inertia(matrix) != expected_inertia:
        raise ValueError(
            "a complete hyperbolic wall Gram matrix must have inertia "
            f"({dimension}, 1, {max(expected_inertia.zero, 0)})"
        )

    edges = _elliptic_subdiagrams(matrix, dimension - 1)
    finite_vertices = _elliptic_subdiagrams(matrix, dimension)
    ideal_vertices, parabolic_neighborhood_audits = _parabolic_vertex_neighbourhoods(
        matrix, dimension
    )

    edge_audits: list[VinbergEdgeAudit] = []
    for edge in edges:
        edge_set = set(edge)
        finite_endpoints = tuple(
            vertex for vertex in finite_vertices if edge_set.issubset(vertex)
        )
        ideal_endpoints = tuple(
            endpoint.facets
            for endpoint in ideal_vertices
            if _is_incident_ideal_extension(edge, endpoint)
        )
        edge_audits.append(
            VinbergEdgeAudit(edge, finite_endpoints, ideal_endpoints)
        )

    compact_bad_edges = tuple(
        audit.facets for audit in edge_audits if len(audit.finite_endpoints) != 2
    )
    finite_volume_bad_edges = tuple(
        audit.facets for audit in edge_audits if audit.endpoint_count != 2
    )
    has_any_vertex = bool(finite_vertices or ideal_vertices)
    bad_parabolic_neighborhood = any(
        not audit.is_rank_dimension_minus_one_parabolic
        for audit in parabolic_neighborhood_audits
    )
    finite_volume = (
        bool(edges)
        and has_any_vertex
        and not finite_volume_bad_edges
        and not bad_parabolic_neighborhood
    )
    compact = (
        bool(edges)
        and bool(finite_vertices)
        and not parabolic_neighborhood_audits
        and not compact_bad_edges
    )
    if compact:
        decision: Decision = "compact"
    elif finite_volume:
        decision = "finite_volume_noncompact"
    else:
        decision = "infinite_volume"

    return VinbergVolumeVerification(
        dimension=dimension,
        facet_count=facet_count,
        simple_roots_complete=True,
        decision=decision,
        finite_volume=finite_volume,
        compact=compact,
        finite_vertices=finite_vertices,
        ideal_vertices=ideal_vertices,
        parabolic_neighborhood_audits=parabolic_neighborhood_audits,
        edge_audits=tuple(edge_audits),
        compact_bad_edges=compact_bad_edges,
        finite_volume_bad_edges=finite_volume_bad_edges,
    )
