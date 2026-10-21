#!/usr/bin/env python3
"""Independent coordinate and nonface replay of the three multiwedge upgrades.

The checker imports none of the record-specific or multiwedge construction
modules.  Exact integer coordinates certify the small convex bases; the rest
is reconstructed from hard-coded clone and join blocks.
"""

from __future__ import annotations

from itertools import combinations
import json


U = 14
TARGET_FACE_POLYNOMIAL = (1, 14, 85, 300, 690, 1088, 1190, 880, 405, 90)
RESIDUAL_FACE_POLYNOMIAL = (1, 9, 30, 50, 45, 18)


DATA = {
    43: {
        "coordinates": {
            U: (0, 0, 0),
            1: (30, 0, 0),
            2: (0, 30, 0),
            11: (0, 0, 30),
            12: (-2, 10, 10),
            0: (10, 10, -2),
            13: (-1, 3, 13),
        },
        "join": (3, 4, 5, 6, 7),
        "clones": (8, 9, 10),
        "edges": ((0, 11), (0, 12), (1, 12), (0, 13), (1, 13), (2, 13)),
        "fives": (
            (3, 4, 5, 6, 7),
            (1, 2, 8, 9, 10),
            (2, 8, 9, 10, 11),
            (8, 9, 10, 11, 12),
        ),
    },
    44: {
        "coordinates": {
            U: (2, 0),
            0: (1, 2),
            1: (-1, 2),
            11: (-2, 0),
            12: (-1, -2),
            13: (1, -2),
        },
        "join": (2, 3, 4, 5, 6),
        "clones": (7, 8, 9, 10),
        "edges": ((0, 11), (0, 12), (1, 12), (0, 13), (1, 13), (11, 13)),
        "fives": (
            (2, 3, 4, 5, 6),
            (1, 7, 8, 9, 10),
            (7, 8, 9, 10, 11),
            (7, 8, 9, 10, 12),
        ),
    },
    47: {
        "coordinates": {
            U: (0, 0, 0),
            1: (30, 0, 0),
            2: (0, 30, 0),
            11: (0, 0, 30),
            0: (10, 10, -2),
            12: (-2, 10, 10),
            13: (10, -2, 10),
        },
        "join": (3, 4, 5, 6, 7),
        "clones": (8, 9, 10),
        "edges": ((0, 11), (0, 12), (1, 12), (0, 13), (2, 13), (12, 13)),
        "fives": (
            (3, 4, 5, 6, 7),
            (1, 2, 8, 9, 10),
            (1, 8, 9, 10, 11),
            (2, 8, 9, 10, 11),
        ),
    },
}


def _cross_2(a, b, c) -> int:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (
        c[0] - a[0]
    )


def _cross_3(a, b, c) -> tuple[int, int, int]:
    u = tuple(b[index] - a[index] for index in range(3))
    v = tuple(c[index] - a[index] for index in range(3))
    return (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )


def _supporting_facets(coordinates) -> tuple[frozenset[int], ...]:
    labels = tuple(sorted(coordinates))
    dimension = len(next(iter(coordinates.values())))
    facets = []
    for labels_on_plane in combinations(labels, dimension):
        points = tuple(coordinates[label] for label in labels_on_plane)
        if dimension == 2:
            signs = [
                _cross_2(points[0], points[1], coordinates[label])
                for label in labels
                if label not in labels_on_plane
            ]
            nondegenerate = points[0] != points[1]
        elif dimension == 3:
            normal = _cross_3(*points)
            signs = [
                sum(
                    normal[index]
                    * (coordinates[label][index] - points[0][index])
                    for index in range(3)
                )
                for label in labels
                if label not in labels_on_plane
            ]
            nondegenerate = normal != (0, 0, 0)
        else:  # pragma: no cover - fixed exact data
            raise ValueError("independent base audit expects dimension two or three")
        if nondegenerate and all(sign != 0 for sign in signs) and (
            all(sign > 0 for sign in signs) or all(sign < 0 for sign in signs)
        ):
            facets.append(frozenset(labels_on_plane))
    return tuple(sorted(facets, key=lambda face: tuple(sorted(face))))


def _all_faces(facets) -> set[frozenset[int]]:
    faces = {frozenset()}
    for facet in facets:
        ordered = tuple(sorted(facet))
        for size in range(1, len(ordered) + 1):
            faces.update(frozenset(face) for face in combinations(ordered, size))
    return faces


def _minimal_nonfaces(vertices, facets) -> tuple[frozenset[int], ...]:
    faces = _all_faces(facets)
    result = []
    for size in range(1, len(vertices) + 1):
        for subset in combinations(vertices, size):
            candidate = frozenset(subset)
            if candidate not in faces and all(
                candidate - {vertex} in faces for vertex in candidate
            ):
                result.append(candidate)
    return tuple(result)


def _face_polynomial(vertices, minimal_nonfaces) -> tuple[int, ...]:
    counts = [0] * (len(vertices) + 1)
    for size in range(len(vertices) + 1):
        for subset in combinations(vertices, size):
            candidate = frozenset(subset)
            if not any(nonface <= candidate for nonface in minimal_nonfaces):
                counts[size] += 1
    while counts and counts[-1] == 0:
        counts.pop()
    return tuple(counts)


def replay() -> dict[str, object]:
    records = []
    for index, data in sorted(DATA.items()):
        coordinates = data["coordinates"]
        base_facets = _supporting_facets(coordinates)
        base_nonfaces = _minimal_nonfaces(tuple(sorted(coordinates)), base_facets)
        expanded = []
        for nonface in base_nonfaces:
            face = set()
            for vertex in nonface:
                face.update(data["clones"] if vertex == U else (vertex,))
            expanded.append(frozenset(face))
        expanded.append(frozenset(data["join"]))
        target = tuple(frozenset(edge) for edge in data["edges"]) + tuple(
            frozenset(face) for face in data["fives"]
        )
        if set(expanded) != set(target):
            raise AssertionError(f"record {index} independent nonface expansion failed")

        residual_vertices = tuple(
            vertex for vertex in range(14) if vertex not in data["join"]
        )
        residual_nonfaces = tuple(
            nonface for nonface in target if nonface != frozenset(data["join"])
        )
        residual = _face_polynomial(residual_vertices, residual_nonfaces)
        full = _face_polynomial(tuple(range(14)), target)
        if residual != RESIDUAL_FACE_POLYNOMIAL or full != TARGET_FACE_POLYNOMIAL:
            raise AssertionError(f"record {index} independent face polynomial failed")
        records.append(
            {
                "record_index": index,
                "base_dimension": len(next(iter(coordinates.values()))),
                "base_vertex_count": len(coordinates),
                "base_facet_count": len(base_facets),
                "base_facets": [
                    list(sorted(facet)) for facet in base_facets
                ],
                "base_minimal_nonface_count": len(base_nonfaces),
                "exact_nonface_expansion": True,
                "residual_face_polynomial": list(residual),
                "full_face_polynomial": list(full),
            }
        )
    return {
        "records": records,
        "ordinary_polytopality_follows_from_facet_wedge_lemma": True,
    }


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2, sort_keys=True))
