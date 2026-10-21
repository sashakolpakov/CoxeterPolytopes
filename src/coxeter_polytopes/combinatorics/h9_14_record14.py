"""Exact ordinary-polytopal residual for H9/14 graph-shadow record 14.

The module freezes one incidence completion found by the whole-shadow
``a=0,...,4`` reduction.  It checks the complete blocker, the exact facet
family, and an ordinary-polytopality construction: the completion is a
five-fold simplicial multiwedge of an explicitly realizable simplicial
four-polytope on nine vertices.  The base realization is certified by an
integral rank-four Gale matrix.

Nothing here supplies Coxeter labels, a Lorentzian Gram matrix, or a finite-
volume hyperbolic chamber.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable

from coxeter_polytopes.combinatorics.facet_completion import generated_face_table
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


TARGET_RECORD_INDEX = 14
TARGET_GRAPH6 = "M???????????o?w??"
TARGET_VERTEX_COUNT = 14
TARGET_DIMENSION = 9
FULL_MASK = (1 << TARGET_VERTEX_COUNT) - 1
TARGET_F_VECTOR = (14, 86, 310, 733, 1192, 1344, 1020, 477, 106)

TARGET_MISSING_EDGES = (
    (0, 12),
    (1, 12),
    (0, 13),
    (1, 13),
    (2, 13),
)
TARGET_MISSING_FOURS = (
    (1, 3, 4, 5),
    (1, 9, 10, 11),
)
TARGET_MISSING_FIVES = (
    (0, 2, 6, 7, 8),
    (2, 8, 9, 10, 11),
    (3, 4, 5, 6, 7),
    (8, 9, 10, 11, 12),
)
EXPECTED_MINIMAL_GRAPH_COVERS = (
    (0, 1, 2),
    (0, 1, 13),
    (12, 13),
)

CLONE_BLOCKS = ((3, 4, 5), (6, 7), (9, 10, 11))
BASE_VERTICES = (0, 1, 2, 3, 6, 8, 9, 12, 13)
BASE_DIMENSION = 4
BASE_GALE_MATRIX = (
    (1, 0, 0, 0, 31, -53, -52, -31, 47),
    (0, 1, 0, 0, 21, -39, -32, -64, 60),
    (0, 0, 1, 0, 23, -21, -22, -50, 45),
    (0, 0, 0, 1, 50, -20, -44, -10, 15),
)
BASE_POSITIVE_DEPENDENCE = (
    107_156,
    652_142,
    43_267,
    470_205,
    32_655,
    42_227,
    58_450,
    128_091,
    172_953,
)


def mask(vertices: Iterable[int]) -> int:
    result = 0
    for vertex in vertices:
        if not 0 <= vertex < TARGET_VERTEX_COUNT:
            raise ValueError("record-14 vertex lies outside [14]")
        result |= 1 << vertex
    return result


def vertices(value: int) -> tuple[int, ...]:
    if value < 0 or value & ~FULL_MASK:
        raise ValueError("record-14 mask lies outside [14]")
    return tuple(vertex for vertex in range(TARGET_VERTEX_COUNT) if value & (1 << vertex))


EDGE_MASKS = tuple(mask(edge) for edge in TARGET_MISSING_EDGES)


def record14_clutter() -> tuple[int, ...]:
    return (
        EDGE_MASKS
        + tuple(mask(face) for face in TARGET_MISSING_FOURS)
        + tuple(mask(face) for face in TARGET_MISSING_FIVES)
    )


def is_transversal(candidate: int, clutter: tuple[int, ...]) -> bool:
    return all(candidate & nonface for nonface in clutter)


def is_minimal_transversal(candidate: int, clutter: tuple[int, ...]) -> bool:
    if not is_transversal(candidate, clutter):
        return False
    return all(
        not is_transversal(candidate ^ (1 << vertex), clutter)
        for vertex in vertices(candidate)
    )


def record14_facet_complements() -> tuple[int, ...]:
    clutter = record14_clutter()
    result = tuple(
        candidate
        for candidate in range(1, FULL_MASK + 1)
        if is_minimal_transversal(candidate, clutter)
    )
    if Counter(value.bit_count() for value in result) != {5: 106}:
        raise AssertionError("record-14 blocker is not pure with 106 five-sets")
    return result


def minimal_graph_covers() -> tuple[int, ...]:
    return tuple(
        candidate
        for candidate in range(1, FULL_MASK + 1)
        if is_minimal_transversal(candidate, EDGE_MASKS)
    )


def _all_faces(facets: Iterable[frozenset[int]]) -> set[frozenset[int]]:
    result = {frozenset()}
    for facet in facets:
        ordered = tuple(sorted(facet))
        for size in range(1, len(ordered) + 1):
            result.update(frozenset(face) for face in combinations(ordered, size))
    return result


def _minimal_nonfaces(
    ground: tuple[int, ...], facets: tuple[frozenset[int], ...]
) -> tuple[frozenset[int], ...]:
    faces = _all_faces(facets)
    result = []
    for size in range(1, len(ground) + 1):
        for chosen in combinations(ground, size):
            candidate = frozenset(chosen)
            if candidate in faces:
                continue
            if all(candidate - {vertex} in faces for vertex in candidate):
                result.append(candidate)
    return tuple(result)


def _determinant_4(columns: tuple[int, int, int, int]) -> int:
    matrix = [
        [BASE_GALE_MATRIX[row][column] for column in columns]
        for row in range(4)
    ]
    sign = 1
    previous = 1
    for pivot in range(3):
        if matrix[pivot][pivot] == 0:
            replacement = next(
                (
                    row
                    for row in range(pivot + 1, 4)
                    if matrix[row][pivot] != 0
                ),
                None,
            )
            if replacement is None:
                return 0
            matrix[pivot], matrix[replacement] = (
                matrix[replacement],
                matrix[pivot],
            )
            sign = -sign
        for row in range(pivot + 1, 4):
            for column in range(pivot + 1, 4):
                matrix[row][column] = (
                    matrix[row][column] * matrix[pivot][pivot]
                    - matrix[row][pivot] * matrix[pivot][column]
                ) // previous
        previous = matrix[pivot][pivot]
    return sign * matrix[3][3]


def _positive_circuit_supports() -> tuple[frozenset[int], ...]:
    result = []
    for support in combinations(range(len(BASE_VERTICES)), 5):
        coefficients = tuple(
            (-1 if index % 2 else 1)
            * _determinant_4(support[:index] + support[index + 1 :])
            for index in range(5)
        )
        if all(value > 0 for value in coefficients) or all(
            value < 0 for value in coefficients
        ):
            result.append(frozenset(support))
    return tuple(result)


def _collapsed_nonfaces() -> tuple[frozenset[int], ...]:
    representative = {
        vertex: block[0] for block in CLONE_BLOCKS for vertex in block
    }
    result = []
    for nonface in record14_clutter():
        original = set(vertices(nonface))
        for block in CLONE_BLOCKS:
            intersection = original.intersection(block)
            if intersection and intersection != set(block):
                raise AssertionError("record-14 clone block is only partly present")
        result.append(
            frozenset(representative.get(vertex, vertex) for vertex in original)
        )
    return tuple(result)


def _expand_base_nonface(nonface: frozenset[int]) -> frozenset[int]:
    blocks = {block[0]: block for block in CLONE_BLOCKS}
    result: set[int] = set()
    for vertex in nonface:
        result.update(blocks.get(vertex, (vertex,)))
    return frozenset(result)


@dataclass(frozen=True)
class Record14PolytopalityAudit:
    decoded_edges: tuple[tuple[int, int], ...]
    minimal_graph_covers: tuple[tuple[int, ...], ...]
    blocker_histogram: tuple[tuple[int, int], ...]
    facet_count: int
    face_vector: tuple[int, ...]
    ridge_degree_histogram: tuple[tuple[int, int], ...]
    base_vertex_count: int
    base_f_vector: tuple[int, ...]
    base_gale_rank: int
    base_nonzero_minor_count: int
    base_positive_circuit_count: int
    base_positive_circuits_match_facets: bool
    base_positive_dependence_verified: bool
    every_base_element_in_and_out_of_positive_circuit: bool
    clone_block_sizes: tuple[int, ...]
    exact_multiwedge_nonface_match: bool
    target_dimension: int

    @property
    def verified(self) -> bool:
        return (
            self.decoded_edges == TARGET_MISSING_EDGES
            and self.minimal_graph_covers == EXPECTED_MINIMAL_GRAPH_COVERS
            and self.blocker_histogram == ((5, 106),)
            and self.facet_count == 106
            and self.face_vector == TARGET_F_VECTOR
            and self.ridge_degree_histogram == ((2, 477),)
            and self.base_vertex_count == 9
            and self.base_f_vector == (9, 28, 38, 19)
            and self.base_gale_rank == 4
            and self.base_nonzero_minor_count == 126
            and self.base_positive_circuit_count == 19
            and self.base_positive_circuits_match_facets
            and self.base_positive_dependence_verified
            and self.every_base_element_in_and_out_of_positive_circuit
            and self.clone_block_sizes == (3, 2, 3)
            and self.exact_multiwedge_nonface_match
            and self.target_dimension == TARGET_DIMENSION
        )

    def to_record(self) -> dict[str, object]:
        return {
            "decoded_edges": [list(edge) for edge in self.decoded_edges],
            "minimal_graph_covers": [list(cover) for cover in self.minimal_graph_covers],
            "minimal_transversal_histogram": dict(self.blocker_histogram),
            "facet_count": self.facet_count,
            "face_vector": list(self.face_vector),
            "ridge_degree_histogram": dict(self.ridge_degree_histogram),
            "base_vertex_count": self.base_vertex_count,
            "base_vertices": list(BASE_VERTICES),
            "base_f_vector": list(self.base_f_vector),
            "base_gale_rank": self.base_gale_rank,
            "base_gale_matrix": [list(row) for row in BASE_GALE_MATRIX],
            "base_nonzero_minor_count": self.base_nonzero_minor_count,
            "base_positive_circuit_count": self.base_positive_circuit_count,
            "base_positive_circuits_match_facets": (
                self.base_positive_circuits_match_facets
            ),
            "base_positive_dependence": list(BASE_POSITIVE_DEPENDENCE),
            "base_positive_dependence_verified": (
                self.base_positive_dependence_verified
            ),
            "every_base_element_in_and_out_of_positive_circuit": (
                self.every_base_element_in_and_out_of_positive_circuit
            ),
            "clone_blocks": [list(block) for block in CLONE_BLOCKS],
            "clone_block_sizes": list(self.clone_block_sizes),
            "simplicial_wedge_count": sum(size - 1 for size in self.clone_block_sizes),
            "exact_multiwedge_nonface_match": self.exact_multiwedge_nonface_match,
            "target_dimension": self.target_dimension,
            "ordinary_polytopal": self.verified,
            "verified": self.verified,
        }


def audit_record14_polytopality() -> Record14PolytopalityAudit:
    order, decoded_edges = decode_graph6(TARGET_GRAPH6)
    if order != TARGET_VERTEX_COUNT:
        raise AssertionError("record-14 graph order changed")
    graph_covers = tuple(vertices(value) for value in minimal_graph_covers())
    facet_complements = record14_facet_complements()
    face_table = generated_face_table(facet_complements)
    face_vector = tuple(
        sum(
            face_table[mask(chosen)]
            for chosen in combinations(range(TARGET_VERTEX_COUNT), size)
        )
        for size in range(1, TARGET_DIMENSION + 1)
    )
    ridge_histogram = Counter(
        sum(complement & mask(chosen) == complement for complement in facet_complements)
        for chosen in combinations(range(TARGET_VERTEX_COUNT), 6)
    )
    ridge_histogram.pop(0, None)

    base_nonfaces = _collapsed_nonfaces()
    base_faces = tuple(
        frozenset(chosen)
        for size in range(1, len(BASE_VERTICES) + 1)
        for chosen in combinations(BASE_VERTICES, size)
        if not any(nonface <= frozenset(chosen) for nonface in base_nonfaces)
    )
    base_facets = tuple(
        face
        for face in base_faces
        if not any(face < other for other in base_faces)
    )
    if set(_minimal_nonfaces(BASE_VERTICES, base_facets)) != set(base_nonfaces):
        raise AssertionError("record-14 base nonface reconstruction changed")
    base_f_vector = tuple(
        sum(len(face) == size for face in _all_faces(base_facets))
        for size in range(1, BASE_DIMENSION + 1)
    )

    minors = tuple(_determinant_4(basis) for basis in combinations(range(9), 4))
    positive_circuits = _positive_circuit_supports()
    base_index = {vertex: index for index, vertex in enumerate(BASE_VERTICES)}
    expected_circuits = {
        frozenset(
            index
            for vertex, index in base_index.items()
            if vertex not in facet
        )
        for facet in base_facets
    }
    covered = set().union(*positive_circuits)
    avoided = set().union(
        *(set(range(9)) - set(circuit) for circuit in positive_circuits)
    )
    dependence_verified = all(value > 0 for value in BASE_POSITIVE_DEPENDENCE) and all(
        sum(
            BASE_GALE_MATRIX[row][column]
            * BASE_POSITIVE_DEPENDENCE[column]
            for column in range(9)
        )
        == 0
        for row in range(4)
    )
    expanded = {_expand_base_nonface(nonface) for nonface in base_nonfaces}
    target_nonfaces = {frozenset(vertices(value)) for value in record14_clutter()}

    audit = Record14PolytopalityAudit(
        decoded_edges=decoded_edges,
        minimal_graph_covers=graph_covers,
        blocker_histogram=tuple(
            sorted(Counter(value.bit_count() for value in facet_complements).items())
        ),
        facet_count=len(facet_complements),
        face_vector=face_vector,
        ridge_degree_histogram=tuple(sorted(ridge_histogram.items())),
        base_vertex_count=len(BASE_VERTICES),
        base_f_vector=base_f_vector,
        base_gale_rank=4,
        base_nonzero_minor_count=sum(value != 0 for value in minors),
        base_positive_circuit_count=len(positive_circuits),
        base_positive_circuits_match_facets=set(positive_circuits)
        == expected_circuits,
        base_positive_dependence_verified=dependence_verified,
        every_base_element_in_and_out_of_positive_circuit=(
            covered == set(range(9)) and avoided == set(range(9))
        ),
        clone_block_sizes=tuple(len(block) for block in CLONE_BLOCKS),
        exact_multiwedge_nonface_match=expanded == target_nonfaces,
        target_dimension=BASE_DIMENSION
        + sum(len(block) - 1 for block in CLONE_BLOCKS),
    )
    if not audit.verified:
        raise AssertionError(
            f"record-14 ordinary-polytopality audit failed: {audit!r}"
        )
    return audit


__all__ = [
    "BASE_GALE_MATRIX",
    "BASE_POSITIVE_DEPENDENCE",
    "CLONE_BLOCKS",
    "TARGET_F_VECTOR",
    "TARGET_GRAPH6",
    "TARGET_MISSING_EDGES",
    "TARGET_MISSING_FOURS",
    "TARGET_MISSING_FIVES",
    "TARGET_RECORD_INDEX",
    "audit_record14_polytopality",
    "record14_clutter",
    "record14_facet_complements",
]
