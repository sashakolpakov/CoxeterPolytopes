r"""Sparse exact Vinberg diagnostics over :math:`\mathbf Q(\sqrt2)`.

The generic volume verifier deliberately uses a simple combinations loop.
That is ideal for the small archived chambers, but it is not viable for the
hundreds of walls occurring in the intrinsic silver K3 search.  This module
enumerates the same elliptic edge diagrams by hereditary positive-definite
backtracking.  Pairwise-positive two-by-two minors are represented by bit
sets, and an exact silver inertia test is still applied to every surviving
extension.

The result is intentionally a *prefix* diagnostic.  Compact stopping requires
two finite endpoints on every edge and no affine component.  Finite-volume
stopping also admits ideal endpoints, but only after Vinberg's full
zero-neighbourhood test, including reducible affine cusp links.  Either
success is one-sided and must be joined to an exactly replayed Vinberg prefix;
failure only identifies currently open edges or contaminated neighbourhoods,
which later walls may cut off.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Callable, Sequence, TypeAlias

from coxeter_polytopes.exact.matrix import Inertia
from coxeter_polytopes.exact.silver import SQRT2, Silver, silver_inertia
from coxeter_polytopes.verify.vinberg_volume import (
    VinbergEdgeAudit,
    VinbergEndpoint,
    VinbergParabolicNeighborhoodAudit,
)


SilverMatrix: TypeAlias = Sequence[Sequence[Silver]]
Progress: TypeAlias = Callable[[str], None]


def _subgram(
    matrix: SilverMatrix, vertices: Sequence[int]
) -> tuple[tuple[Silver, ...], ...]:
    return tuple(
        tuple(matrix[left][right] for right in vertices) for left in vertices
    )


def _is_positive_definite(
    matrix: SilverMatrix, vertices: tuple[int, ...]
) -> bool:
    classified = _finite_34_cartan(matrix, vertices)
    if classified is not None:
        return classified
    return silver_inertia(_subgram(matrix, vertices), embedding=0) == Inertia(
        len(vertices), 0, 0
    )


def _finite_34_cartan(
    matrix: SilverMatrix, vertices: tuple[int, ...]
) -> bool | None:
    """Recognize finite Coxeter diagrams using only labels three and four.

    Norm-two silver roots have diagonal two.  On every pair occurring in the
    present Vinberg prefixes, a positive two-node subdiagram has pairing
    ``0``, ``-1``, or ``-sqrt(2)``.  In that case finite type is a purely
    graph-theoretic ADE/B/F classification.  Any other value selects the
    general exact-inertia fallback rather than being guessed.
    """

    if any(matrix[node][node] != 2 for node in vertices):
        return None
    adjacency = {node: set() for node in vertices}
    four_edges: set[frozenset[int]] = set()
    for offset, left in enumerate(vertices):
        for right in vertices[:offset]:
            value = matrix[left][right]
            if value == 0:
                continue
            if value == -1:
                adjacency[left].add(right)
                adjacency[right].add(left)
                continue
            if value == -SQRT2:
                adjacency[left].add(right)
                adjacency[right].add(left)
                four_edges.add(frozenset((left, right)))
                continue
            return None

    unseen = set(vertices)
    while unseen:
        start = min(unseen)
        component = {start}
        frontier = [start]
        unseen.remove(start)
        while frontier:
            node = frontier.pop()
            new = adjacency[node] & unseen
            unseen.difference_update(new)
            component.update(new)
            frontier.extend(new)
        edge_count = sum(len(adjacency[node] & component) for node in component) // 2
        if edge_count != len(component) - 1:
            return False
        degrees = {node: len(adjacency[node] & component) for node in component}
        local_fours = {
            edge for edge in four_edges if edge.issubset(component)
        }
        if not local_fours:
            if max(degrees.values(), default=0) <= 2:
                continue  # A_r, including A_1.
            branches = [node for node, degree in degrees.items() if degree == 3]
            if len(branches) != 1 or any(degree > 3 for degree in degrees.values()):
                return False
            branch = branches[0]
            arms: list[int] = []
            for neighbor in adjacency[branch] & component:
                length = 1
                previous, current = branch, neighbor
                while degrees[current] == 2:
                    following = next(
                        node
                        for node in adjacency[current] & component
                        if node != previous
                    )
                    previous, current = current, following
                    length += 1
                arms.append(length)
            arm_tuple = tuple(sorted(arms))
            if not (
                arm_tuple[:2] == (1, 1)
                or arm_tuple in {(1, 2, 2), (1, 2, 3), (1, 2, 4)}
            ):
                return False
            continue

        if len(local_fours) != 1 or max(degrees.values(), default=0) > 2:
            return False
        four = next(iter(local_fours))
        endpoints = {node for node, degree in degrees.items() if degree <= 1}
        if four & endpoints:
            continue  # B_r/C_r, including I_2(4).
        if len(component) == 4:
            continue  # The central edge of F_4.
        return False
    return True


def _least_bit_index(bits: int) -> int:
    return (bits & -bits).bit_length() - 1


@dataclass(frozen=True, slots=True)
class SilverSparseVinbergAudit:
    """Exact sparse Vinberg edge, finite-vertex, and ideal-cusp census."""

    dimension: int
    facet_count: int
    backtracking_positive_prefix_visits: tuple[int, ...]
    elliptic_edges: tuple[tuple[int, ...], ...]
    finite_vertices: tuple[tuple[int, ...], ...]
    finite_endpoint_indices_by_edge: tuple[tuple[int, ...], ...]
    connected_affine_subdiagrams: tuple[tuple[int, ...], ...]
    ideal_vertices: tuple[VinbergEndpoint, ...]
    parabolic_neighborhood_audits: tuple[
        VinbergParabolicNeighborhoodAudit, ...
    ]
    ideal_endpoint_indices_by_edge: tuple[tuple[int, ...], ...]
    edge_audits: tuple[VinbergEdgeAudit, ...]
    invalid_finite_angle_pairs: tuple[tuple[int, int, Silver], ...]

    @property
    def finite_endpoint_histogram(self) -> tuple[tuple[int, int], ...]:
        return tuple(
            sorted(Counter(map(len, self.finite_endpoint_indices_by_edge)).items())
        )

    @property
    def ideal_endpoint_histogram(self) -> tuple[tuple[int, int], ...]:
        return tuple(
            sorted(Counter(map(len, self.ideal_endpoint_indices_by_edge)).items())
        )

    @property
    def finite_plus_ideal_endpoint_histogram(self) -> tuple[tuple[int, int], ...]:
        return tuple(
            sorted(Counter(audit.endpoint_count for audit in self.edge_audits).items())
        )

    @property
    def open_edges(self) -> tuple[tuple[int, ...], ...]:
        """Legacy compact-open edges, counting finite endpoints only."""

        return tuple(
            edge
            for edge, endpoints in zip(
                self.elliptic_edges,
                self.finite_endpoint_indices_by_edge,
                strict=True,
            )
            if len(endpoints) != 2
        )

    @property
    def finite_volume_open_edges(self) -> tuple[tuple[int, ...], ...]:
        return tuple(
            audit.facets for audit in self.edge_audits if audit.endpoint_count != 2
        )

    @property
    def bad_parabolic_neighborhoods(
        self,
    ) -> tuple[VinbergParabolicNeighborhoodAudit, ...]:
        return tuple(
            audit
            for audit in self.parabolic_neighborhood_audits
            if not audit.is_rank_dimension_minus_one_parabolic
        )

    @property
    def compact_stopping_available(self) -> bool:
        return (
            bool(self.elliptic_edges)
            and bool(self.finite_vertices)
            and not self.open_edges
            and not self.connected_affine_subdiagrams
            and not self.invalid_finite_angle_pairs
        )

    @property
    def finite_volume_stopping_available(self) -> bool:
        return (
            bool(self.elliptic_edges)
            and bool(self.finite_vertices or self.ideal_vertices)
            and not self.finite_volume_open_edges
            and not self.bad_parabolic_neighborhoods
            and not self.invalid_finite_angle_pairs
        )


class _SilverDiagramSearch:
    def __init__(
        self,
        matrix: SilverMatrix,
        dimension: int,
        progress: Progress | None = None,
    ) -> None:
        self.matrix = tuple(tuple(Silver.coerce(value) for value in row) for row in matrix)
        self.dimension = dimension
        self.size = len(self.matrix)
        self.progress = progress
        if dimension < 2:
            raise ValueError("hyperbolic dimension must be at least two")
        if any(len(row) != self.size for row in self.matrix):
            raise ValueError("wall Gram matrix must be square")
        if any(
            self.matrix[left][right] != self.matrix[right][left]
            for left in range(self.size)
            for right in range(left)
        ):
            raise ValueError("wall Gram matrix must be symmetric")
        if any(self.matrix[index][index].sign(0) <= 0 for index in range(self.size)):
            raise ValueError("wall normals must have positive square")
        if any(self.matrix[index][index] != 2 for index in range(self.size)):
            raise ValueError("silver prefix diagnostics require norm-two walls")
        if any(
            self.matrix[left][right].sign(0) > 0
            for left in range(self.size)
            for right in range(left)
        ):
            raise ValueError("wall Gram matrix is not acute")

        self.compatible = [0] * self.size
        self.nonorthogonal = [0] * self.size
        self.invalid_finite_angle_pairs: list[tuple[int, int, Silver]] = []
        for left in range(self.size):
            diagonal_left = self.matrix[left][left]
            for right in range(left):
                pairing = self.matrix[left][right]
                if pairing:
                    self.nonorthogonal[left] |= 1 << right
                    self.nonorthogonal[right] |= 1 << left
                minor = (
                    diagonal_left * self.matrix[right][right] - pairing * pairing
                )
                if minor.sign(0) > 0:
                    self.compatible[left] |= 1 << right
                    self.compatible[right] |= 1 << left
                    if pairing not in (0, -1, -SQRT2):
                        self.invalid_finite_angle_pairs.append(
                            (right, left, pairing)
                        )
            if self.progress is not None and (left + 1) % 50 == 0:
                self.progress(f"pair graph rows {left + 1}/{self.size}")
        self.all_bits = (1 << self.size) - 1
        self.orthogonal = [
            self.all_bits & ~(1 << vertex) & ~self.nonorthogonal[vertex]
            for vertex in range(self.size)
        ]
        self._positive_cache: dict[tuple[int, ...], bool] = {}
        self._affine_component_cache: dict[tuple[int, ...], bool] = {}

    def positive(self, vertices: tuple[int, ...]) -> bool:
        cached = self._positive_cache.get(vertices)
        if cached is None:
            cached = _is_positive_definite(self.matrix, vertices)
            self._positive_cache[vertices] = cached
        return cached

    def elliptic(self, order: int) -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...]]:
        """Enumerate every positive-definite principal subset of fixed order."""

        if order < 1 or order > self.size:
            return (), ()
        all_bits = (1 << self.size) - 1
        output: list[tuple[int, ...]] = []
        counts = [0] * (order + 1)

        def visit(prefix: tuple[int, ...], candidates: int) -> None:
            needed = order - len(prefix)
            if candidates.bit_count() < needed:
                return
            remaining = candidates
            while remaining:
                vertex = _least_bit_index(remaining)
                remaining &= ~(1 << vertex)
                if remaining.bit_count() + 1 < needed:
                    break
                extended = (*prefix, vertex)
                if len(extended) >= 2 and not self.positive(extended):
                    continue
                counts[len(extended)] += 1
                if len(extended) == order:
                    output.append(extended)
                    continue
                # ``remaining`` already contains only vertices larger than
                # ``vertex`` and compatible with every earlier prefix node.
                visit(extended, remaining & self.compatible[vertex])
                if self.progress is not None and not prefix:
                    self.progress(
                        f"elliptic roots through {vertex + 1}/{self.size}; "
                        f"rank-{order} total={len(output)}"
                    )

        visit((), all_bits)
        return tuple(output), tuple(counts[1:])

    def finite_endpoints(
        self, edges: Sequence[tuple[int, ...]]
    ) -> tuple[tuple[tuple[int, ...], ...], tuple[tuple[int, ...], ...]]:
        """Extend every elliptic rank ``n-1`` edge by one finite node."""

        all_bits = (1 << self.size) - 1
        vertices: set[tuple[int, ...]] = set()
        endpoint_rows: list[tuple[int, ...]] = []
        for edge_index, edge in enumerate(edges):
            candidates = all_bits
            edge_bits = 0
            for node in edge:
                candidates &= self.compatible[node]
                edge_bits |= 1 << node
            candidates &= ~edge_bits
            endpoints: list[int] = []
            while candidates:
                node = _least_bit_index(candidates)
                candidates &= ~(1 << node)
                vertex = tuple(sorted((*edge, node)))
                if self.positive(vertex):
                    endpoints.append(node)
                    vertices.add(vertex)
            endpoint_rows.append(tuple(endpoints))
            if self.progress is not None and (edge_index + 1) % 1000 == 0:
                self.progress(
                    f"endpoint rows {edge_index + 1}/{len(edges)}; "
                    f"finite vertices={len(vertices)}"
                )
        return tuple(sorted(vertices)), tuple(endpoint_rows)

    def connected_affine(self) -> tuple[tuple[int, ...], ...]:
        """Enumerate irreducible affine diagrams by elliptic growth."""

        current = {(vertex,) for vertex in range(self.size)}
        affine: set[tuple[int, ...]] = set()
        maximum_order = min(self.size, self.dimension + 2)
        for order in range(1, maximum_order + 1):
            following: set[tuple[int, ...]] = set()
            for vertices in sorted(current):
                if order >= 2:
                    # Almost every tuple in the large silver prefixes is a
                    # finite ADE/B/F diagram.  Classify those from their
                    # labelled graph and reserve exact field elimination for
                    # the much smaller non-finite residue.  This preserves
                    # the exact test while avoiding millions of redundant
                    # inertia computations during connected growth.
                    finite_34 = _finite_34_cartan(self.matrix, vertices)
                    if finite_34 is not True:
                        gram = _subgram(self.matrix, vertices)
                        observed = silver_inertia(gram, embedding=0)
                        if observed == Inertia(order - 1, 0, 1) and all(
                            self.positive(
                                vertices[:index] + vertices[index + 1 :]
                            )
                            for index in range(order)
                        ):
                            affine.add(vertices)
                            continue
                        if observed != Inertia(order, 0, 0):
                            continue
                if order == maximum_order:
                    continue
                vertices_bits = sum(1 << vertex for vertex in vertices)
                adjacent = 0
                for vertex in vertices:
                    adjacent |= self.nonorthogonal[vertex]
                adjacent &= ~vertices_bits
                if order >= 2:
                    # Every proper subdiagram of an irreducible affine
                    # diagram is elliptic.  In particular a new node must
                    # form a positive two-node diagram with every old node.
                    # Keep the unrestricted order-one growth so rank-two
                    # affine diagrams are still detected, then impose this
                    # necessary condition on all later extensions.
                    for vertex in vertices:
                        adjacent &= self.compatible[vertex]
                while adjacent:
                    vertex = _least_bit_index(adjacent)
                    adjacent &= ~(1 << vertex)
                    following.add(tuple(sorted((*vertices, vertex))))
            current = following
            if self.progress is not None:
                self.progress(
                    f"affine growth order {order}; next={len(current)}, "
                    f"affine={len(affine)}"
                )
            if not current:
                break
        return tuple(sorted(affine))

    def components_from_bits(self, vertices_bits: int) -> tuple[tuple[int, ...], ...]:
        """Return exact diagram components using bit-set adjacency."""

        remaining = vertices_bits
        output: list[tuple[int, ...]] = []
        while remaining:
            first_bit = remaining & -remaining
            remaining &= ~first_bit
            component_bits = 0
            frontier = first_bit
            while frontier:
                node_bit = frontier & -frontier
                frontier &= ~node_bit
                node = node_bit.bit_length() - 1
                component_bits |= node_bit
                neighbours = self.nonorthogonal[node] & remaining
                remaining &= ~neighbours
                frontier |= neighbours
            output.append(
                tuple(
                    vertex
                    for vertex in range(self.size)
                    if component_bits & (1 << vertex)
                )
            )
        return tuple(output)

    def is_affine_component(self, vertices: tuple[int, ...]) -> bool:
        """Silver-field analogue of ``is_affine_parabolic_component``."""

        cached = self._affine_component_cache.get(vertices)
        if cached is not None:
            return cached
        vertices_bits = sum(1 << vertex for vertex in vertices)
        cached = (
            len(vertices) >= 2
            and len(self.components_from_bits(vertices_bits)) == 1
            and silver_inertia(_subgram(self.matrix, vertices), embedding=0)
            == Inertia(len(vertices) - 1, 0, 1)
            and all(
                self.positive(vertices[:index] + vertices[index + 1 :])
                for index in range(len(vertices))
            )
        )
        self._affine_component_cache[vertices] = cached
        return cached

    def parabolic_vertex_neighborhoods(
        self,
        affine_subdiagrams: Sequence[tuple[int, ...]],
    ) -> tuple[
        tuple[VinbergEndpoint, ...],
        tuple[VinbergParabolicNeighborhoodAudit, ...],
    ]:
        """Apply the full zero-neighbourhood cusp test with bit intersections.

        This is exactly the generic ``vinberg_volume`` semantics.  For each
        connected affine ``J``, all walls orthogonal to every node of ``J``
        are included; no favourable parabolic subset may be selected from
        that full neighbourhood.
        """

        target_rank = self.dimension - 1
        audits: list[VinbergParabolicNeighborhoodAudit] = []
        endpoints_by_facets: dict[tuple[int, ...], VinbergEndpoint] = {}
        for index, affine_component in enumerate(affine_subdiagrams):
            component_bits = sum(1 << vertex for vertex in affine_component)
            zero_bits = self.all_bits & ~component_bits
            for vertex in affine_component:
                zero_bits &= self.orthogonal[vertex]
            neighborhood_bits = component_bits | zero_bits
            neighborhood = tuple(
                vertex
                for vertex in range(self.size)
                if neighborhood_bits & (1 << vertex)
            )
            components = self.components_from_bits(neighborhood_bits)
            correct_size = len(neighborhood) == target_rank + len(components)
            is_full_parabolic = correct_size and all(
                self.is_affine_component(component) for component in components
            )
            audit = VinbergParabolicNeighborhoodAudit(
                affine_component=affine_component,
                zero_neighborhood=neighborhood,
                neighborhood_components=components,
                is_rank_dimension_minus_one_parabolic=is_full_parabolic,
            )
            audits.append(audit)
            if is_full_parabolic:
                endpoints_by_facets.setdefault(
                    neighborhood,
                    VinbergEndpoint("ideal", neighborhood, components),
                )
            if self.progress is not None and (index + 1) % 1000 == 0:
                self.progress(
                    f"full affine neighbourhoods {index + 1}/"
                    f"{len(affine_subdiagrams)}; ideal="
                    f"{len(endpoints_by_facets)}"
                )
        return (
            tuple(sorted(endpoints_by_facets.values(), key=lambda row: row.facets)),
            tuple(audits),
        )

    def ideal_endpoint_rows(
        self,
        edges: Sequence[tuple[int, ...]],
        ideal_vertices: Sequence[VinbergEndpoint],
    ) -> tuple[tuple[int, ...], ...]:
        """Index ideal endpoints by incident edge without scanning all pairs."""

        edge_index = {edge: index for index, edge in enumerate(edges)}
        rows: list[set[int]] = [set() for _ in edges]
        for endpoint_index, endpoint in enumerate(ideal_vertices):
            facet_bits = sum(1 << vertex for vertex in endpoint.facets)

            def omit_one_from_each(component_index: int, omitted_bits: int) -> None:
                if component_index == len(endpoint.affine_components):
                    edge_bits = facet_bits & ~omitted_bits
                    edge = tuple(
                        vertex
                        for vertex in range(self.size)
                        if edge_bits & (1 << vertex)
                    )
                    if len(edge) != self.dimension - 1:
                        raise AssertionError(
                            "rank-n-1 ideal edge construction changed"
                        )
                    found = edge_index.get(edge)
                    if found is None:
                        raise AssertionError(
                            "a full parabolic endpoint generated a non-elliptic edge"
                        )
                    rows[found].add(endpoint_index)
                    return
                for vertex in endpoint.affine_components[component_index]:
                    omit_one_from_each(
                        component_index + 1, omitted_bits | (1 << vertex)
                    )

            omit_one_from_each(0, 0)
        return tuple(tuple(sorted(row)) for row in rows)


def audit_silver_vinberg_wall_gram(
    wall_gram: SilverMatrix,
    dimension: int,
    *,
    progress: Progress | None = None,
) -> SilverSparseVinbergAudit:
    """Run the scalable exact compact/finite-volume stopping diagnostic.

    No completeness flag is accepted here.  A caller must separately prove
    that the walls are an exact, correctly ordered Vinberg prefix before
    using either successful stopping flag as a completeness theorem.  A
    failed prefix audit remains inconclusive because later walls can close
    apparent open ends.
    """

    search = _SilverDiagramSearch(wall_gram, dimension, progress)
    edges, prefix_counts = search.elliptic(dimension - 1)
    finite_vertices, finite_endpoint_rows = search.finite_endpoints(edges)
    affine = search.connected_affine()
    ideal_vertices, neighborhood_audits = search.parabolic_vertex_neighborhoods(
        affine
    )
    ideal_endpoint_rows = search.ideal_endpoint_rows(edges, ideal_vertices)
    edge_audits = tuple(
        VinbergEdgeAudit(
            facets=edge,
            finite_endpoints=tuple(
                tuple(sorted((*edge, node)))
                for node in finite_endpoint_rows[index]
            ),
            ideal_endpoints=tuple(
                ideal_vertices[endpoint_index].facets
                for endpoint_index in ideal_endpoint_rows[index]
            ),
        )
        for index, edge in enumerate(edges)
    )
    return SilverSparseVinbergAudit(
        dimension=dimension,
        facet_count=len(wall_gram),
        backtracking_positive_prefix_visits=prefix_counts,
        elliptic_edges=edges,
        finite_vertices=finite_vertices,
        finite_endpoint_indices_by_edge=finite_endpoint_rows,
        connected_affine_subdiagrams=affine,
        ideal_vertices=ideal_vertices,
        parabolic_neighborhood_audits=neighborhood_audits,
        ideal_endpoint_indices_by_edge=ideal_endpoint_rows,
        edge_audits=edge_audits,
        invalid_finite_angle_pairs=tuple(search.invalid_finite_angle_pairs),
    )


__all__ = ["SilverSparseVinbergAudit", "audit_silver_vinberg_wall_gram"]
