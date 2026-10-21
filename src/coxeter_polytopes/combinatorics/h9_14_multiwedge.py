"""Exact multiwedge factorizations for H9/14 records 43, 44, and 47.

Each fourteen-vertex nerve is reconstructed from a small polytopal sphere by
one simplicial multiwedge and a join with the boundary of a four-simplex.
The finite audit works entirely at the level of minimal nonfaces.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable

from coxeter_polytopes.combinatorics.h9_14_record43 import (
    TARGET_GRAPH6 as RECORD43_GRAPH6,
    TARGET_MISSING_EDGES as RECORD43_EDGES,
    TARGET_MISSING_FIVES as RECORD43_FIVES,
)
from coxeter_polytopes.combinatorics.h9_14_record44 import (
    TARGET_GRAPH6 as RECORD44_GRAPH6,
    TARGET_MISSING_EDGES as RECORD44_EDGES,
    TARGET_MISSING_FIVES as RECORD44_FIVES,
)
from coxeter_polytopes.combinatorics.h9_14_record47 import (
    TARGET_GRAPH6 as RECORD47_GRAPH6,
    TARGET_MISSING_EDGES as RECORD47_EDGES,
    TARGET_MISSING_FIVES as RECORD47_FIVES,
)


CLONE_VERTEX = 14
SIMPLEX_JOIN_DIMENSION = 4
TARGET_POLYTOPE_DIMENSION = 9
TARGET_FACET_COUNT = 14


@dataclass(frozen=True)
class StackingStep:
    new_vertex: int
    replaced_facet: tuple[int, ...]

    def to_record(self) -> dict[str, object]:
        return {
            "new_vertex": self.new_vertex,
            "replaced_facet": list(self.replaced_facet),
        }


@dataclass(frozen=True)
class MultiwedgeSpecification:
    record_index: int
    graph6: str
    missing_edges: tuple[tuple[int, int], ...]
    missing_fives: tuple[tuple[int, ...], ...]
    simplex_join_vertices: tuple[int, ...]
    clone_vertices: tuple[int, ...]
    base_vertices: tuple[int, ...]
    initial_base_facets: tuple[tuple[int, ...], ...]
    stacking_steps: tuple[StackingStep, ...]
    base_description: str


def _tetrahedron_boundary(vertices: tuple[int, int, int, int]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(face) for face in combinations(vertices, 3))


RECORD43 = MultiwedgeSpecification(
    record_index=43,
    graph6=RECORD43_GRAPH6,
    missing_edges=RECORD43_EDGES,
    missing_fives=RECORD43_FIVES,
    simplex_join_vertices=(3, 4, 5, 6, 7),
    clone_vertices=(8, 9, 10),
    base_vertices=(0, 1, 2, 11, 12, 13, CLONE_VERTEX),
    initial_base_facets=_tetrahedron_boundary((1, 2, 11, CLONE_VERTEX)),
    stacking_steps=(
        StackingStep(12, (2, 11, CLONE_VERTEX)),
        StackingStep(0, (1, 2, CLONE_VERTEX)),
        StackingStep(13, (11, 12, CLONE_VERTEX)),
    ),
    base_description="boundary of an explicit three-times stacked tetrahedron",
)

RECORD44 = MultiwedgeSpecification(
    record_index=44,
    graph6=RECORD44_GRAPH6,
    missing_edges=RECORD44_EDGES,
    missing_fives=RECORD44_FIVES,
    simplex_join_vertices=(2, 3, 4, 5, 6),
    clone_vertices=(7, 8, 9, 10),
    base_vertices=(0, 1, 11, 12, 13, CLONE_VERTEX),
    initial_base_facets=(
        (CLONE_VERTEX, 0),
        (0, 1),
        (1, 11),
        (11, 12),
        (12, 13),
        (13, CLONE_VERTEX),
    ),
    stacking_steps=(),
    base_description="boundary cycle of a hexagon",
)

RECORD47 = MultiwedgeSpecification(
    record_index=47,
    graph6=RECORD47_GRAPH6,
    missing_edges=RECORD47_EDGES,
    missing_fives=RECORD47_FIVES,
    simplex_join_vertices=(3, 4, 5, 6, 7),
    clone_vertices=(8, 9, 10),
    base_vertices=(0, 1, 2, 11, 12, 13, CLONE_VERTEX),
    initial_base_facets=_tetrahedron_boundary((1, 2, 11, CLONE_VERTEX)),
    stacking_steps=(
        StackingStep(0, (1, 2, CLONE_VERTEX)),
        StackingStep(12, (2, 11, CLONE_VERTEX)),
        StackingStep(13, (1, 11, CLONE_VERTEX)),
    ),
    base_description="boundary of an explicit three-times stacked tetrahedron",
)

SPECIFICATIONS = (RECORD43, RECORD44, RECORD47)


def _stacked_facets(specification: MultiwedgeSpecification) -> tuple[frozenset[int], ...]:
    facets = {frozenset(facet) for facet in specification.initial_base_facets}
    for step in specification.stacking_steps:
        replaced = frozenset(step.replaced_facet)
        if replaced not in facets:
            raise AssertionError(
                f"record {specification.record_index} stacking facet disappeared"
            )
        facets.remove(replaced)
        for ridge in combinations(sorted(replaced), len(replaced) - 1):
            facets.add(frozenset((*ridge, step.new_vertex)))
    return tuple(sorted(facets, key=lambda face: (len(face), tuple(sorted(face)))))


def _all_faces(facets: Iterable[frozenset[int]]) -> set[frozenset[int]]:
    faces = {frozenset()}
    for facet in facets:
        ordered = tuple(sorted(facet))
        for size in range(1, len(ordered) + 1):
            faces.update(frozenset(face) for face in combinations(ordered, size))
    return faces


def _minimal_nonfaces(
    vertices: tuple[int, ...], facets: tuple[frozenset[int], ...]
) -> tuple[frozenset[int], ...]:
    faces = _all_faces(facets)
    result = []
    for size in range(1, len(vertices) + 1):
        for subset in combinations(vertices, size):
            candidate = frozenset(subset)
            if candidate in faces:
                continue
            if all(candidate - {vertex} in faces for vertex in candidate):
                result.append(candidate)
    return tuple(sorted(result, key=lambda face: (len(face), tuple(sorted(face)))))


def _face_vector(facets: tuple[frozenset[int], ...]) -> tuple[int, ...]:
    faces = _all_faces(facets)
    dimension = max(map(len, facets))
    return tuple(sum(len(face) == size for face in faces) for size in range(1, dimension + 1))


def _expand_nonface(
    nonface: frozenset[int], clone_vertices: tuple[int, ...]
) -> frozenset[int]:
    expanded = set()
    for vertex in nonface:
        if vertex == CLONE_VERTEX:
            expanded.update(clone_vertices)
        else:
            expanded.add(vertex)
    return frozenset(expanded)


def _group_by_size(
    faces: Iterable[frozenset[int]],
) -> tuple[tuple[int, tuple[tuple[int, ...], ...]], ...]:
    groups: dict[int, list[tuple[int, ...]]] = {}
    for face in faces:
        groups.setdefault(len(face), []).append(tuple(sorted(face)))
    return tuple(
        (size, tuple(sorted(entries)))
        for size, entries in sorted(groups.items())
    )


@dataclass(frozen=True)
class MultiwedgePolytopalityAudit:
    record_index: int
    graph6: str
    base_description: str
    base_vertices: tuple[int, ...]
    base_facets: tuple[tuple[int, ...], ...]
    base_face_vector: tuple[int, ...]
    base_minimal_nonfaces_by_size: tuple[
        tuple[int, tuple[tuple[int, ...], ...]], ...
    ]
    stacking_steps: tuple[StackingStep, ...]
    clone_vertex: int
    clone_vertices: tuple[int, ...]
    simplex_join_vertices: tuple[int, ...]
    expanded_minimal_nonfaces_by_size: tuple[
        tuple[int, tuple[tuple[int, ...], ...]], ...
    ]
    target_minimal_nonfaces_by_size: tuple[
        tuple[int, tuple[tuple[int, ...], ...]], ...
    ]
    base_polytope_dimension: int
    wedge_count: int
    multiwedge_factor_dimension: int
    product_dimension: int
    product_facet_count: int
    exact_nonface_match: bool

    @property
    def verified(self) -> bool:
        expected_base_vector = (
            (6, 6) if self.record_index == 44 else (7, 15, 10)
        )
        return (
            self.record_index in (43, 44, 47)
            and self.base_face_vector == expected_base_vector
            and self.exact_nonface_match
            and self.multiwedge_factor_dimension == 5
            and self.product_dimension == TARGET_POLYTOPE_DIMENSION
            and self.product_facet_count == TARGET_FACET_COUNT
        )

    def to_record(self) -> dict[str, object]:
        return {
            "record_index": self.record_index,
            "graph6": self.graph6,
            "base_description": self.base_description,
            "base_vertices": list(self.base_vertices),
            "base_facets": [list(face) for face in self.base_facets],
            "base_face_vector": list(self.base_face_vector),
            "base_minimal_nonfaces_by_size": {
                str(size): [list(face) for face in faces]
                for size, faces in self.base_minimal_nonfaces_by_size
            },
            "stacking_steps": [step.to_record() for step in self.stacking_steps],
            "clone_vertex": self.clone_vertex,
            "clone_vertices": list(self.clone_vertices),
            "clone_multiplicity": len(self.clone_vertices),
            "simplicial_wedge_count": self.wedge_count,
            "simplex_join_vertices": list(self.simplex_join_vertices),
            "expanded_minimal_nonfaces_by_size": {
                str(size): [list(face) for face in faces]
                for size, faces in self.expanded_minimal_nonfaces_by_size
            },
            "target_minimal_nonfaces_by_size": {
                str(size): [list(face) for face in faces]
                for size, faces in self.target_minimal_nonfaces_by_size
            },
            "exact_nonface_match": self.exact_nonface_match,
            "base_polytope_dimension": self.base_polytope_dimension,
            "multiwedge_factor_dimension": self.multiwedge_factor_dimension,
            "product_with_simplex_dimension": self.product_dimension,
            "product_facet_count": self.product_facet_count,
            "ordinary_polytopal": self.verified,
        }


def audit_specification(
    specification: MultiwedgeSpecification,
) -> MultiwedgePolytopalityAudit:
    facets = _stacked_facets(specification)
    base_nonfaces = _minimal_nonfaces(specification.base_vertices, facets)
    expanded = tuple(
        _expand_nonface(nonface, specification.clone_vertices)
        for nonface in base_nonfaces
    ) + (frozenset(specification.simplex_join_vertices),)
    target = tuple(frozenset(edge) for edge in specification.missing_edges) + tuple(
        frozenset(face) for face in specification.missing_fives
    )
    base_dimension = max(map(len, facets))
    wedge_count = len(specification.clone_vertices) - 1
    audit = MultiwedgePolytopalityAudit(
        record_index=specification.record_index,
        graph6=specification.graph6,
        base_description=specification.base_description,
        base_vertices=specification.base_vertices,
        base_facets=tuple(tuple(sorted(face)) for face in facets),
        base_face_vector=_face_vector(facets),
        base_minimal_nonfaces_by_size=_group_by_size(base_nonfaces),
        stacking_steps=specification.stacking_steps,
        clone_vertex=CLONE_VERTEX,
        clone_vertices=specification.clone_vertices,
        simplex_join_vertices=specification.simplex_join_vertices,
        expanded_minimal_nonfaces_by_size=_group_by_size(expanded),
        target_minimal_nonfaces_by_size=_group_by_size(target),
        base_polytope_dimension=base_dimension,
        wedge_count=wedge_count,
        multiwedge_factor_dimension=base_dimension + wedge_count,
        product_dimension=(
            base_dimension + wedge_count + SIMPLEX_JOIN_DIMENSION
        ),
        product_facet_count=(
            len(specification.base_vertices)
            + wedge_count
            + SIMPLEX_JOIN_DIMENSION
            + 1
        ),
        exact_nonface_match=set(expanded) == set(target),
    )
    if not audit.verified:
        raise AssertionError(
            f"record {specification.record_index} multiwedge audit failed"
        )
    return audit


def audit_all() -> tuple[MultiwedgePolytopalityAudit, ...]:
    return tuple(audit_specification(specification) for specification in SPECIFICATIONS)


__all__ = [
    "MultiwedgePolytopalityAudit",
    "RECORD43",
    "RECORD44",
    "RECORD47",
    "SPECIFICATIONS",
    "audit_all",
    "audit_specification",
]
