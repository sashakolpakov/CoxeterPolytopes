"""Exact global facet-selection model for the H9/14 nerve problem.

The second graph-shadow layer retains five-set complements of possible
nine-vertex facets.  This module turns that support into a finite Boolean
decision problem.  Its constraints are all necessary for the boundary of a
simplicial 9-polytope whose missing faces support Lanner diagrams:

* every selected complement is in the certified ridge-support core;
* every eight-face is contained in zero or two selected facets;
* every independent triple is a face;
* all face indicators are *defined*, in both directions, by the selected
  facets and obey one complete Dehn--Sommerville schedule; and
* no inclusion-minimal nonface has more than five vertices.

The last condition is important.  Merely finding a closed pseudomanifold
with the right number of facets is not enough: its large minimal nonfaces
cannot be Lanner diagrams.

``solve_facet_selection`` uses Z3 only as a search backend.  Importing this
module and independently auditing a SAT witness require no optional package.
An UNSAT return is an exact backend result, but it is not by itself the
independent proof certificate required for the final frontier theorem.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from itertools import combinations, product
from math import comb
from time import monotonic
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.missing_edges import (
    TARGET_POLYTOPE_DIMENSION,
    TARGET_VERTEX_COUNT,
    decode_graph6,
    macaulay_raise,
)
from coxeter_polytopes.combinatorics.nerve_support import (
    FACET_COMPLEMENT_SIZE,
)
from coxeter_polytopes.combinatorics.local_link import (
    analyze_local_link_graph,
)
from coxeter_polytopes.combinatorics.simplicial_homology import (
    first_mod2_homology_sphere_link_failure,
)


FULL_VERTEX_MASK = (1 << TARGET_VERTEX_COUNT) - 1
MAX_LANNER_ORDER = 5


def _mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        result |= 1 << vertex
    return result


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(
        vertex
        for vertex in range(TARGET_VERTEX_COUNT)
        if mask & (1 << vertex)
    )


def _edge_masks(graph6: str) -> tuple[int, ...]:
    vertex_count, edges = decode_graph6(graph6)
    if vertex_count != TARGET_VERTEX_COUNT:
        raise ValueError("facet completion is specific to 14 vertices")
    return tuple((1 << first) | (1 << second) for first, second in edges)


def independent_masks_by_size(graph6: str) -> tuple[tuple[int, ...], ...]:
    """Enumerate all subsets avoiding the missing edges, grouped by size."""

    edge_masks = _edge_masks(graph6)
    result: list[tuple[int, ...]] = []
    for size in range(TARGET_VERTEX_COUNT + 1):
        masks = []
        for vertices in combinations(range(TARGET_VERTEX_COUNT), size):
            mask = _mask(vertices)
            if not any(mask & edge == edge for edge in edge_masks):
                masks.append(mask)
        result.append(tuple(masks))
    return tuple(result)


def generated_face_table(
    facet_complement_masks: Iterable[int],
) -> tuple[bool, ...]:
    """Return the complete face table of the generated pure complex."""

    selected = tuple(sorted(set(facet_complement_masks)))
    table = [False] * (1 << TARGET_VERTEX_COUNT)
    table[0] = True
    for complement in selected:
        if complement < 0 or complement & ~FULL_VERTEX_MASK:
            raise ValueError("facet complement lies outside the ground set")
        if complement.bit_count() != FACET_COMPLEMENT_SIZE:
            raise ValueError("facet complements must be five-sets")
        facet = FULL_VERTEX_MASK ^ complement
        subset = facet
        while True:
            table[subset] = True
            if subset == 0:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def explicit_product_facet_selection(graph6: str) -> tuple[int, ...] | None:
    """Return one of the two exact product-fixture completions when matched.

    These witnesses calibrate the global solver and are useful even when a
    bounded generic backend times out.  They are ordinary polytopal spheres,
    not hyperbolic Coxeter examples.
    """

    vertex_count, raw_edges = decode_graph6(graph6)
    if vertex_count != TARGET_VERTEX_COUNT:
        return None
    edges = {tuple(sorted(edge)) for edge in raw_edges}
    degrees = Counter(vertex for edge in edges for vertex in edge)
    active = tuple(sorted(vertex for vertex in range(vertex_count) if degrees[vertex]))
    isolated = tuple(vertex for vertex in range(vertex_count) if not degrees[vertex])

    facets: list[int] = []
    # Delta_1^3 x Delta_3^2: H is three pairwise disjoint missing edges.
    if (
        len(edges) == 3
        and len(active) == 6
        and all(degrees[vertex] == 1 for vertex in active)
        and len(isolated) == 8
    ):
        pairs = tuple(sorted(edges))
        first_block = isolated[:4]
        second_block = isolated[4:]
        for choices in product(*pairs):
            for omitted_first in first_block:
                for omitted_second in second_block:
                    facet_vertices = set(choices)
                    facet_vertices.update(set(first_block) - {omitted_first})
                    facet_vertices.update(set(second_block) - {omitted_second})
                    facets.append(_mask(facet_vertices))

    # C_5 x Delta_3 x Delta_4: H is a five-cycle and nine isolates.
    elif (
        len(edges) == 5
        and len(active) == 5
        and all(degrees[vertex] == 2 for vertex in active)
        and len(isolated) == 9
    ):
        # Degree two alone could describe a disconnected cycle union, but on
        # five active vertices it is necessarily one C5 (no simple 2-cycle).
        allowed_pairs = tuple(
            pair for pair in combinations(active, 2) if pair not in edges
        )
        if len(allowed_pairs) != 5:
            return None
        first_block = isolated[:4]
        second_block = isolated[4:]
        for allowed_pair in allowed_pairs:
            for omitted_first in first_block:
                for omitted_second in second_block:
                    facet_vertices = set(allowed_pair)
                    facet_vertices.update(set(first_block) - {omitted_first})
                    facet_vertices.update(set(second_block) - {omitted_second})
                    facets.append(_mask(facet_vertices))
    else:
        return None

    return tuple(sorted(FULL_VERTEX_MASK ^ facet for facet in facets))


def forced_small_transversal_obstruction(
    graph6: str, schedules: Sequence[Mapping[str, object]]
) -> dict[str, int] | None:
    """Certify the special ``A4=A5=empty`` transversal contradiction.

    If every admissible schedule has no missing four- or five-faces, the
    missing-edge graph is the complete minimal-nonface clutter.  Facet
    complements are then exactly its inclusion-minimal vertex covers.  A
    graph with at most four edges has such a cover of size at most four,
    contradicting the required complement size five.
    """

    if not schedules:
        return None
    if any(
        schedule.get("missing_four_count") != 0
        or schedule.get("missing_five_count_offset") != 0
        for schedule in schedules
    ):
        return None
    _, edges = decode_graph6(graph6)
    edge_count = len(edges)
    if edge_count >= FACET_COMPLEMENT_SIZE:
        return None
    return {
        "missing_four_count": 0,
        "missing_five_count": 0,
        "missing_edge_count": edge_count,
        "minimal_vertex_cover_upper_bound": edge_count,
        "required_facet_complement_size": FACET_COMPLEMENT_SIZE,
    }


def minimal_nonfaces_from_table(
    face_table: Sequence[bool],
) -> dict[int, tuple[int, ...]]:
    """Enumerate all inclusion-minimal nonfaces from a full face table."""

    if len(face_table) != 1 << TARGET_VERTEX_COUNT:
        raise ValueError("face table has the wrong ground-set size")
    result: dict[int, list[int]] = {}
    for subset in range(1, 1 << TARGET_VERTEX_COUNT):
        if face_table[subset]:
            continue
        if all(
            face_table[subset ^ (1 << vertex)]
            for vertex in _vertices(subset)
        ):
            result.setdefault(subset.bit_count(), []).append(subset)
    return {size: tuple(values) for size, values in sorted(result.items())}


def _ridge_groups(
    facet_complement_masks: Sequence[int],
) -> dict[int, tuple[int, ...]]:
    groups: dict[int, list[int]] = {}
    for complement in facet_complement_masks:
        outside = FULL_VERTEX_MASK ^ complement
        while outside:
            bit = outside & -outside
            groups.setdefault(complement | bit, []).append(complement)
            outside ^= bit
    return {
        six_set: tuple(sorted(group))
        for six_set, group in groups.items()
    }


def _dual_components(
    facet_complement_masks: Sequence[int],
) -> tuple[tuple[int, ...], ...]:
    selected = tuple(sorted(facet_complement_masks))
    if not selected:
        return ()
    adjacency = {complement: set() for complement in selected}
    for group in _ridge_groups(selected).values():
        if len(group) == 2:
            first, second = group
            adjacency[first].add(second)
            adjacency[second].add(first)
    unseen = set(selected)
    components: list[tuple[int, ...]] = []
    while unseen:
        root = min(unseen)
        seen = {root}
        queue = deque([root])
        while queue:
            current = queue.popleft()
            for neighbor in adjacency[current] - seen:
                seen.add(neighbor)
                queue.append(neighbor)
        unseen.difference_update(seen)
        components.append(tuple(sorted(seen)))
    return tuple(components)


def _first_local_euler_failure(
    face_table: Sequence[bool],
) -> tuple[int, int, int] | None:
    """Return ``(face, observed, expected)`` for the first bad link.

    The link of a k-vertex face in an 8-dimensional sphere is a sphere of
    dimension ``8-k``.  We check its Euler characteristic exactly for every
    face of at most eight vertices.  This is necessary, but not sufficient,
    for the links to be homology spheres.
    """

    for face in range(1 << TARGET_VERTEX_COUNT):
        if not face_table[face] or face.bit_count() > 8:
            continue
        remaining = FULL_VERTEX_MASK ^ face
        link_euler = 0
        subset = remaining
        while subset:
            if face_table[face | subset]:
                link_euler += 1 if subset.bit_count() % 2 else -1
            subset = (subset - 1) & remaining
        link_dimension = 8 - face.bit_count()
        expected = 1 + (-1) ** link_dimension
        if link_euler != expected:
            return face, link_euler, expected
    return None


def _link_g_entries_from_table(
    face_table: Sequence[bool], face_mask: int
) -> tuple[int, tuple[int, ...]]:
    """Return ``(link_vertex_count, g_1,...,g_floor(D/2))``."""

    face_size = face_mask.bit_count()
    polytope_dimension = TARGET_POLYTOPE_DIMENSION - face_size
    remaining_vertices = _vertices(FULL_VERTEX_MASK ^ face_mask)
    link_f = {-1: 1}
    for dimension in range(polytope_dimension):
        order = dimension + 1
        link_f[dimension] = sum(
            face_table[face_mask | _mask(extension)]
            for extension in combinations(remaining_vertices, order)
        )
    h = []
    for index in range(polytope_dimension // 2 + 1):
        h.append(
            sum(
                (-1) ** (index - face_order)
                * comb(
                    polytope_dimension - face_order,
                    polytope_dimension - index,
                )
                * link_f[face_order - 1]
                for face_order in range(index + 1)
            )
        )
    return link_f[0], tuple(h[index] - h[index - 1] for index in range(1, len(h)))


def _first_strict_local_g_failure(
    graph6: str, face_table: Sequence[bool]
) -> tuple[int, int, int, int] | None:
    """Return ``(face, order, observed_g, minimum_g)`` on failure."""

    independent = independent_masks_by_size(graph6)
    strict_link_orders = {1: (2, 3, 4), 2: (2, 3), 3: (2,)}
    for face_size, strict_orders in strict_link_orders.items():
        polytope_dimension = TARGET_POLYTOPE_DIMENSION - face_size
        for face_mask in independent[face_size]:
            if not face_table[face_mask]:
                return face_mask, 0, -1, 0
            link_vertices, g_entries = _link_g_entries_from_table(
                face_table, face_mask
            )
            # A simplex face figure is impossible here: its boundary has one
            # minimal nonface of order D+1>5, contradicting link inheritance
            # of the Lanner bound.
            if link_vertices < polytope_dimension + 2:
                return face_mask, 0, link_vertices, polytope_dimension + 2
            for order, observed in enumerate(g_entries, start=1):
                minimum = 1 if order in strict_orders else 0
                if observed < minimum:
                    return face_mask, order, observed, minimum
    return None


def _first_local_macaulay_failure(
    graph6: str, face_table: Sequence[bool]
) -> tuple[int, int, int, int] | None:
    """Return ``(face, next_order, observed, upper_bound)`` on failure."""

    independent = independent_masks_by_size(graph6)
    for face_size in (1, 2, 3):
        for face_mask in independent[face_size]:
            _, g_entries = _link_g_entries_from_table(face_table, face_mask)
            for previous_order in range(1, len(g_entries)):
                previous = g_entries[previous_order - 1]
                observed = g_entries[previous_order]
                if previous < 0:
                    return face_mask, previous_order + 1, observed, -1
                upper = macaulay_raise(previous, previous_order)
                if observed > upper:
                    return face_mask, previous_order + 1, observed, upper
    return None


def _missing_four_link_intervals_hold(
    graph6: str,
    missing_four_masks: Sequence[int],
) -> bool:
    local = analyze_local_link_graph(graph6)
    if not local.passes_edge_and_triple_feasibility:
        return False
    for figure in local.vertex_figures:
        interval = figure.missing_four_incidence_interval
        if interval is None:
            return False
        incidence = sum(
            bool(subset & (1 << figure.vertex)) for subset in missing_four_masks
        )
        if not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.edge_figures:
        interval = figure.missing_four_incidence_interval
        if interval is None:
            return False
        edge_mask = _mask(figure.edge)
        incidence = sum(
            subset & edge_mask == edge_mask for subset in missing_four_masks
        )
        if not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.triple_figures:
        triple_mask = _mask(figure.triple)
        incidence = sum(
            subset & triple_mask == triple_mask for subset in missing_four_masks
        )
        if incidence > figure.missing_four_incidence_max:
            return False
    return True


@dataclass(frozen=True)
class GlobalFacetSelectionAudit:
    """Independent exact audit of one selected facet collection."""

    graph6: str
    selected_facet_complement_masks: tuple[int, ...]
    selected_inside_support_core: bool
    facet_count: int
    face_vector: tuple[int, ...]
    ridge_degree_histogram: dict[int, int]
    all_ridge_degrees_two: bool
    dual_component_count: int
    dual_graph_connected: bool
    covers_all_independent_triples: bool
    missing_four_count: int
    minimal_nonfaces_by_size: dict[int, tuple[int, ...]]
    largest_minimal_nonface_order: int
    lanner_order_bound_holds: bool
    matching_schedule_index: int | None
    matches_dehn_sommerville_schedule: bool
    missing_four_link_intervals_hold: bool
    first_strict_local_g_failure: tuple[int, int, int, int] | None
    strict_local_g_constraints_hold: bool
    first_local_macaulay_failure: tuple[int, int, int, int] | None
    local_macaulay_constraints_hold: bool
    first_local_euler_failure: tuple[int, int, int] | None
    all_links_have_spherical_euler_characteristic: bool
    first_mod2_homology_failure: tuple[int, int, tuple[int, ...]] | None
    all_links_are_mod2_homology_spheres: bool

    @property
    def passes_selection_model(self) -> bool:
        """Whether the witness satisfies every constraint sent to Z3."""

        return (
            self.selected_inside_support_core
            and self.all_ridge_degrees_two
            and self.covers_all_independent_triples
            and self.lanner_order_bound_holds
            and self.matches_dehn_sommerville_schedule
            and self.missing_four_link_intervals_hold
            and self.strict_local_g_constraints_hold
            and self.local_macaulay_constraints_hold
        )

    @property
    def passes_connected_euler_audit(self) -> bool:
        """Stronger necessary manifold audit, not part of base SAT model."""

        return (
            self.passes_selection_model
            and self.dual_graph_connected
            and self.all_links_have_spherical_euler_characteristic
        )

    @property
    def passes_connected_homology_audit(self) -> bool:
        """Whether every face link passes exact mod-two homology."""

        return (
            self.passes_connected_euler_audit
            and self.all_links_are_mod2_homology_spheres
        )

    def to_record(self) -> dict[str, object]:
        return {
            "graph6": self.graph6,
            "selected_facet_complement_masks": list(
                self.selected_facet_complement_masks
            ),
            "selected_inside_support_core": self.selected_inside_support_core,
            "facet_count": self.facet_count,
            "face_vector": list(self.face_vector),
            "ridge_degree_histogram": {
                str(degree): count
                for degree, count in sorted(self.ridge_degree_histogram.items())
            },
            "all_ridge_degrees_two": self.all_ridge_degrees_two,
            "dual_component_count": self.dual_component_count,
            "dual_graph_connected": self.dual_graph_connected,
            "covers_all_independent_triples": self.covers_all_independent_triples,
            "missing_four_count": self.missing_four_count,
            "minimal_nonfaces_by_size": {
                str(size): list(values)
                for size, values in self.minimal_nonfaces_by_size.items()
            },
            "largest_minimal_nonface_order": self.largest_minimal_nonface_order,
            "lanner_order_bound_holds": self.lanner_order_bound_holds,
            "matching_schedule_index": self.matching_schedule_index,
            "matches_dehn_sommerville_schedule": self.matches_dehn_sommerville_schedule,
            "missing_four_link_intervals_hold": (
                self.missing_four_link_intervals_hold
            ),
            "first_strict_local_g_failure": (
                list(self.first_strict_local_g_failure)
                if self.first_strict_local_g_failure is not None
                else None
            ),
            "strict_local_g_constraints_hold": self.strict_local_g_constraints_hold,
            "first_local_macaulay_failure": (
                list(self.first_local_macaulay_failure)
                if self.first_local_macaulay_failure is not None
                else None
            ),
            "local_macaulay_constraints_hold": self.local_macaulay_constraints_hold,
            "first_local_euler_failure": (
                list(self.first_local_euler_failure)
                if self.first_local_euler_failure is not None
                else None
            ),
            "all_links_have_spherical_euler_characteristic": (
                self.all_links_have_spherical_euler_characteristic
            ),
            "first_mod2_homology_failure": (
                [
                    self.first_mod2_homology_failure[0],
                    self.first_mod2_homology_failure[1],
                    list(self.first_mod2_homology_failure[2]),
                ]
                if self.first_mod2_homology_failure is not None
                else None
            ),
            "all_links_are_mod2_homology_spheres": (
                self.all_links_are_mod2_homology_spheres
            ),
            "passes_selection_model": self.passes_selection_model,
            "passes_connected_euler_audit": self.passes_connected_euler_audit,
            "passes_connected_homology_audit": (
                self.passes_connected_homology_audit
            ),
        }


def audit_global_facet_selection(
    graph6: str,
    support_core_masks: Sequence[int],
    schedules: Sequence[Mapping[str, object]],
    selected_facet_complement_masks: Iterable[int],
) -> GlobalFacetSelectionAudit:
    """Verify a SAT witness without using the search backend."""

    selected = tuple(sorted(set(selected_facet_complement_masks)))
    if not selected:
        raise ValueError("a facet selection witness cannot be empty")
    support = set(support_core_masks)
    inside_support = all(complement in support for complement in selected)
    face_table = generated_face_table(selected)
    face_vector = tuple(
        sum(
            face_table[_mask(vertices)]
            for vertices in combinations(range(TARGET_VERTEX_COUNT), size)
        )
        for size in range(1, TARGET_POLYTOPE_DIMENSION + 1)
    )

    ridge_histogram = Counter(
        len(group) for group in _ridge_groups(selected).values()
    )
    all_ridges_two = bool(ridge_histogram) and set(ridge_histogram) == {2}
    components = _dual_components(selected)

    edge_masks = _edge_masks(graph6)
    covers_triples = True
    for vertices in combinations(range(TARGET_VERTEX_COUNT), 3):
        triple = _mask(vertices)
        if any(triple & edge == edge for edge in edge_masks):
            continue
        if not face_table[triple]:
            covers_triples = False
            break

    minimal_nonfaces = minimal_nonfaces_from_table(face_table)
    largest_order = max(minimal_nonfaces, default=0)
    missing_four_count = len(minimal_nonfaces.get(4, ()))
    matching_schedule_index: int | None = None
    for index, schedule in enumerate(schedules):
        raw_vector = schedule.get("f_vector")
        if (
            schedule.get("missing_four_count") == missing_four_count
            and isinstance(raw_vector, (list, tuple))
            and tuple(raw_vector) == face_vector
        ):
            matching_schedule_index = index
            break

    local_euler_failure = _first_local_euler_failure(face_table)
    missing_four_masks = minimal_nonfaces.get(4, ())
    link_intervals_hold = _missing_four_link_intervals_hold(
        graph6, missing_four_masks
    )
    strict_local_g_failure = _first_strict_local_g_failure(graph6, face_table)
    local_macaulay_failure = _first_local_macaulay_failure(graph6, face_table)
    homology_failure = first_mod2_homology_sphere_link_failure(
        face_table, TARGET_POLYTOPE_DIMENSION - 1
    )
    serialized_homology_failure = (
        (
            homology_failure.face_mask,
            homology_failure.expected_sphere_dimension,
            homology_failure.observed_reduced_betti,
        )
        if homology_failure is not None
        else None
    )
    return GlobalFacetSelectionAudit(
        graph6=graph6,
        selected_facet_complement_masks=selected,
        selected_inside_support_core=inside_support,
        facet_count=len(selected),
        face_vector=face_vector,
        ridge_degree_histogram=dict(sorted(ridge_histogram.items())),
        all_ridge_degrees_two=all_ridges_two,
        dual_component_count=len(components),
        dual_graph_connected=len(components) == 1,
        covers_all_independent_triples=covers_triples,
        missing_four_count=missing_four_count,
        minimal_nonfaces_by_size=minimal_nonfaces,
        largest_minimal_nonface_order=largest_order,
        lanner_order_bound_holds=largest_order <= MAX_LANNER_ORDER,
        matching_schedule_index=matching_schedule_index,
        matches_dehn_sommerville_schedule=matching_schedule_index is not None,
        missing_four_link_intervals_hold=link_intervals_hold,
        first_strict_local_g_failure=strict_local_g_failure,
        strict_local_g_constraints_hold=strict_local_g_failure is None,
        first_local_macaulay_failure=local_macaulay_failure,
        local_macaulay_constraints_hold=local_macaulay_failure is None,
        first_local_euler_failure=local_euler_failure,
        all_links_have_spherical_euler_characteristic=(
            local_euler_failure is None
        ),
        first_mod2_homology_failure=serialized_homology_failure,
        all_links_are_mod2_homology_spheres=(homology_failure is None),
    )


@dataclass(frozen=True)
class FacetSelectionSolveResult:
    """One finite Boolean-solver disposition."""

    status: str
    selected_facet_complement_masks: tuple[int, ...]
    backend: str
    backend_version: str
    elapsed_seconds: float
    assertion_count: int
    reason_unknown: str | None

    def to_record(self) -> dict[str, object]:
        return {
            "status": self.status,
            "selected_facet_complement_masks": list(
                self.selected_facet_complement_masks
            ),
            "backend": self.backend,
            "backend_version": self.backend_version,
            "elapsed_seconds": self.elapsed_seconds,
            "assertion_count": self.assertion_count,
            "reason_unknown": self.reason_unknown,
        }


@dataclass(frozen=True)
class FacetSelectionConstraintSystem:
    """A composable Z3 instance of the exact facet-selection model.

    The fields deliberately expose only the primary facet literals in
    addition to the solver.  Auxiliary face, ridge, and schedule variables
    remain implementation details, while callers such as the joint Gale gate
    can couple further *necessary* constraints to the same primary literals.
    """

    solver: object
    core: tuple[int, ...]
    variables: tuple[object, ...]
    facet_literals: Mapping[int, object]


def _passes_nonsimplex_local_precheck(graph6: str) -> bool:
    """Return the exact local-link guard shared by build and solve paths."""

    local_link_audit = analyze_local_link_graph(graph6)
    return (
        local_link_audit.passes_edge_and_triple_feasibility
        and all(figure.g1 > 0 for figure in local_link_audit.vertex_figures)
        and all(figure.g1 > 0 for figure in local_link_audit.edge_figures)
        and all(
            figure.candidate_link_vertex_count >= 8
            for figure in local_link_audit.triple_figures
        )
    )


def build_facet_selection_constraint_system(
    graph6: str,
    support_core_masks: Sequence[int],
    schedules: Sequence[Mapping[str, object]],
    *,
    timeout_ms: int | None = None,
) -> FacetSelectionConstraintSystem:
    """Build, but do not solve, the exact global facet-selection model.

    This is the single production constraint builder used by
    :func:`solve_facet_selection` and by compositional necessary-condition
    experiments.  The input must pass the independently exact local-link
    precheck; standalone solving turns a failed precheck directly into an
    exact UNSAT disposition before calling this function.
    """

    if not _passes_nonsimplex_local_precheck(graph6):
        raise ValueError(
            "facet constraint system fails the exact local-link precheck"
        )
    local_link_audit = analyze_local_link_graph(graph6)

    try:
        import z3  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - depends on local extras
        raise RuntimeError(
            "solve_facet_selection requires the optional z3-solver package"
        ) from exc

    if not support_core_masks:
        raise ValueError("support core cannot be empty")
    core = tuple(sorted(set(support_core_masks)))
    if tuple(support_core_masks) != core:
        raise ValueError("support core must be sorted and duplicate-free")
    if not schedules:
        raise ValueError("at least one Dehn--Sommerville schedule is required")

    solver = z3.Solver()
    if timeout_ms is not None:
        if timeout_ms <= 0:
            raise ValueError("timeout must be positive")
        solver.set(timeout=timeout_ms)

    variables = tuple(z3.Bool(f"x_{index}") for index in range(len(core)))
    core_index = {complement: index for index, complement in enumerate(core)}

    # A six-set contains either zero or exactly two selected complements.
    ridge_groups: dict[int, list[int]] = {}
    for index, complement in enumerate(core):
        outside = FULL_VERTEX_MASK ^ complement
        while outside:
            bit = outside & -outside
            ridge_groups.setdefault(complement | bit, []).append(index)
            outside ^= bit
    for ridge_index, indices in enumerate(ridge_groups.values()):
        active = z3.Bool(f"ridge_active_{ridge_index}")
        solver.add(
            z3.PbEq(
                [(variables[index], 1) for index in indices]
                + [(active, -2)],
                0,
            )
        )

    independent = independent_masks_by_size(graph6)

    # Build containment lists by enumerating subsets of each possible facet.
    cover_indices: dict[int, list[int]] = {}
    for index, complement in enumerate(core):
        facet_vertices = _vertices(FULL_VERTEX_MASK ^ complement)
        for size in range(1, TARGET_POLYTOPE_DIMENSION):
            for vertices in combinations(facet_vertices, size):
                cover_indices.setdefault(_mask(vertices), []).append(index)

    face: dict[int, object] = {0: z3.BoolVal(True)}
    face_counts: dict[int, object] = {}
    for size in range(1, TARGET_POLYTOPE_DIMENSION + 1):
        size_faces = []
        for subset in independent[size]:
            covering = [variables[index] for index in cover_indices.get(subset, ())]
            if size <= 3:
                value = z3.BoolVal(True)
                # The no-missing-triples theorem says every such independent
                # set is a face, not merely that it is combinatorially allowed.
                solver.add(z3.Or(covering))
            elif size == TARGET_POLYTOPE_DIMENSION:
                complement = FULL_VERTEX_MASK ^ subset
                index = core_index.get(complement)
                value = variables[index] if index is not None else z3.BoolVal(False)
            else:
                value = z3.Bool(f"face_{subset}")
                solver.add(value == z3.Or(covering))
            face[subset] = value
            size_faces.append(value)
        face_counts[size] = z3.Sum(
            [z3.If(value, 1, 0) for value in size_faces]
        )

    # Select one whole DS schedule.  The archived schedules form one
    # consecutive affine family in the missing-four count.  We verify that
    # fact before using its much smaller linear encoding.
    parsed_schedules: list[tuple[int, tuple[int, ...]]] = []
    for schedule in schedules:
        missing_four_count = schedule.get("missing_four_count")
        raw_vector = schedule.get("f_vector")
        if not isinstance(missing_four_count, int):
            raise ValueError("DS schedule lacks an integer missing-four count")
        if not isinstance(raw_vector, (list, tuple)) or len(raw_vector) != 9:
            raise ValueError("DS schedule has an invalid f-vector")
        parsed_schedules.append(
            (missing_four_count, tuple(int(entry) for entry in raw_vector))
        )
    parsed_schedules.sort()
    expected_counts = list(
        range(parsed_schedules[0][0], parsed_schedules[-1][0] + 1)
    )
    if [count for count, _ in parsed_schedules] != expected_counts:
        raise ValueError("DS schedules must cover one consecutive interval")
    base_count, base_vector = parsed_schedules[0]
    if len(parsed_schedules) > 1:
        slopes = tuple(
            parsed_schedules[1][1][index] - base_vector[index]
            for index in range(9)
        )
        for offset, (_, vector) in enumerate(parsed_schedules):
            if vector != tuple(
                base_vector[index] + offset * slopes[index]
                for index in range(9)
            ):
                raise ValueError("DS schedule family is not affine")
    else:
        slopes = (0,) * 9
    missing_four = z3.Int("missing_four_count")
    solver.add(
        missing_four >= parsed_schedules[0][0],
        missing_four <= parsed_schedules[-1][0],
        missing_four == len(independent[4]) - face_counts[4],
    )
    for size in range(4, TARGET_POLYTOPE_DIMENSION + 1):
        solver.add(
            face_counts[size]
            == base_vector[size - 1]
            + (missing_four - base_count) * slopes[size - 1]
        )

    # Apply the graph-derived Macaulay intervals directly to the missing
    # four-face indicators.  Since every independent triple is a face, an
    # independent four-set A is a missing four-face exactly when face_A is
    # false.
    missing_four_indicators = {
        subset: z3.Not(face[subset]) for subset in independent[4]
    }
    for figure in local_link_audit.vertex_figures:
        interval = figure.missing_four_incidence_interval
        if interval is None:  # defensive; the graph precheck handled this
            raise AssertionError("surviving graph has an empty vertex interval")
        incidence = z3.Sum(
            [
                z3.If(indicator, 1, 0)
                for subset, indicator in missing_four_indicators.items()
                if subset & (1 << figure.vertex)
            ]
        )
        solver.add(incidence >= interval[0], incidence <= interval[1])
    for figure in local_link_audit.edge_figures:
        interval = figure.missing_four_incidence_interval
        if interval is None:
            raise AssertionError("surviving graph has an empty edge interval")
        edge_mask = _mask(figure.edge)
        incidence = z3.Sum(
            [
                z3.If(indicator, 1, 0)
                for subset, indicator in missing_four_indicators.items()
                if subset & edge_mask == edge_mask
            ]
        )
        solver.add(incidence >= interval[0], incidence <= interval[1])
    for figure in local_link_audit.triple_figures:
        triple_mask = _mask(figure.triple)
        incidence = z3.Sum(
            [
                z3.If(indicator, 1, 0)
                for subset, indicator in missing_four_indicators.items()
                if subset & triple_mask == triple_mask
            ]
        )
        solver.add(incidence <= figure.missing_four_incidence_max)

    def link_h_entry(
        face_mask: int,
        polytope_dimension: int,
        index: int,
    ):
        """Return one exact h-entry of the link of ``face_mask``."""

        face_size = face_mask.bit_count()
        terms = []
        for face_order in range(index + 1):
            coefficient = (
                (-1) ** (index - face_order)
                * comb(
                    polytope_dimension - face_order,
                    polytope_dimension - index,
                )
            )
            if face_order == 0:
                terms.append(coefficient)
                continue
            containing = [
                value
                for superset, value in face.items()
                if superset.bit_count() == face_size + face_order
                and superset & face_mask == face_mask
            ]
            terms.append(
                coefficient
                * z3.Sum([z3.If(value, 1, 0) for value in containing])
            )
        return z3.Sum(terms)

    # Low-dimensional links are themselves boundaries of simplicial convex
    # polytopes.  Their g-vectors are nonnegative.  Murai--Nevo makes the
    # indicated inequalities strict for every nonsimplex link: equality
    # would force a missing link face of order >5, impossible because every
    # link missing face contains a global minimal nonface of no larger order.
    strict_link_orders = {1: (2, 3, 4), 2: (2, 3), 3: (2,)}
    for face_size, strict_orders in strict_link_orders.items():
        link_dimension = TARGET_POLYTOPE_DIMENSION - face_size
        for face_mask in independent[face_size]:
            h_entries = [
                link_h_entry(face_mask, link_dimension, index)
                for index in range(link_dimension // 2 + 1)
            ]
            g_entries = [
                h_entries[index] - h_entries[index - 1]
                for index in range(1, len(h_entries))
            ]
            link_vertex_count = h_entries[1] + link_dimension
            # The D-simplex branch is itself forbidden: its boundary has a
            # missing face of order D+1>5.  Every figure here is nonsimplex.
            solver.add(link_vertex_count >= link_dimension + 2)
            for order, g_entry in enumerate(g_entries, start=1):
                if order in strict_orders:
                    solver.add(g_entry >= 1)
                else:
                    solver.add(g_entry >= 0)
            maximum_previous = 4
            for previous_order in range(1, len(g_entries)):
                previous = g_entries[previous_order - 1]
                next_entry = g_entries[previous_order]
                solver.add(
                    z3.Or(
                        *[
                            z3.And(
                                previous == value,
                                next_entry
                                <= macaulay_raise(value, previous_order),
                            )
                            for value in range(maximum_previous + 1)
                        ]
                    )
                )
                maximum_previous = macaulay_raise(
                    maximum_previous, previous_order
                )

    # Every minimal nonface of a compact Coxeter polytope is Lanner and has
    # order at most five.  For each larger independent S, either S is a face
    # or at least one immediate subset is already a nonface.  Nonindependent
    # sets contain a missing edge and need no additional clause.
    for size in range(MAX_LANNER_ORDER + 1, TARGET_VERTEX_COUNT + 1):
        for subset in independent[size]:
            subset_face = face.get(subset, z3.BoolVal(False))
            predecessor_nonfaces = [
                z3.Not(face.get(subset ^ (1 << vertex), z3.BoolVal(False)))
                for vertex in _vertices(subset)
            ]
            solver.add(z3.Or(subset_face, *predecessor_nonfaces))

    return FacetSelectionConstraintSystem(
        solver=solver,
        core=core,
        variables=variables,
        facet_literals={
            complement: variables[index]
            for index, complement in enumerate(core)
        },
    )


def solve_facet_selection(
    graph6: str,
    support_core_masks: Sequence[int],
    schedules: Sequence[Mapping[str, object]],
    *,
    timeout_ms: int | None = None,
) -> FacetSelectionSolveResult:
    """Decide the exact global facet-selection model with Z3.

    ``timeout_ms=None`` requests an unbounded exact decision.  ``unknown`` is
    always preserved as unknown and is never interpreted as nonexistence.
    """

    start = monotonic()
    # Keep the nonsimplex checks explicit here as a defense-in-depth guard:
    # the boundary of a D-simplex has a missing face of order D+1>5.
    if not _passes_nonsimplex_local_precheck(graph6):
        return FacetSelectionSolveResult(
            status="unsat",
            selected_facet_complement_masks=(),
            backend="exact_local_link_precheck",
            backend_version="1.0",
            elapsed_seconds=monotonic() - start,
            assertion_count=0,
            reason_unknown=None,
        )

    system = build_facet_selection_constraint_system(
        graph6,
        support_core_masks,
        schedules,
        timeout_ms=timeout_ms,
    )
    try:
        import z3  # type: ignore[import-not-found]
    except ImportError as exc:  # pragma: no cover - depends on local extras
        raise RuntimeError(
            "solve_facet_selection requires the optional z3-solver package"
        ) from exc
    solver = system.solver
    assertion_count = len(solver.assertions())
    disposition = solver.check()
    elapsed = monotonic() - start
    if disposition == z3.sat:
        model = solver.model()
        selected = tuple(
            complement
            for index, complement in enumerate(system.core)
            if z3.is_true(
                model.eval(system.variables[index], model_completion=True)
            )
        )
        status = "sat"
        reason_unknown = None
    elif disposition == z3.unsat:
        selected = ()
        status = "unsat"
        reason_unknown = None
    else:
        selected = ()
        status = "unknown"
        reason_unknown = solver.reason_unknown()

    return FacetSelectionSolveResult(
        status=status,
        selected_facet_complement_masks=selected,
        backend="z3",
        backend_version=z3.get_version_string(),
        elapsed_seconds=elapsed,
        assertion_count=assertion_count,
        reason_unknown=reason_unknown,
    )
