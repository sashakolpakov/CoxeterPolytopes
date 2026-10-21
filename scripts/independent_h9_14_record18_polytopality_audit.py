#!/usr/bin/env python3
"""Hostile stdlib replay of the record-18 rational-polytopality archive.

This script deliberately does not import the proof implementation.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record18_polytopality.json"
EXPECTED_DIGEST = (
    "6dae56855a8380ee8cffd8b3e2572df90164f98e483e717f9323de8289d4cb90"
)
EXPECTED_SOURCE_HASH = (
    "9de308c9a4dc54b84daaa2a8e4c38ca062da4dee28f352a41a2beceabaff165e"
)
EXPECTED_EDGES = ((0, 12), (1, 12), (2, 13), (3, 13), (12, 13))
EXPECTED_FIVES = (
    (0, 2, 4, 5, 6),
    (0, 5, 6, 7, 8),
    (0, 3, 7, 8, 9),
    (3, 7, 8, 9, 10),
    (2, 4, 5, 6, 11),
    (1, 2, 4, 10, 11),
    (1, 3, 9, 10, 11),
    (1, 4, 9, 10, 11),
)
EXPECTED_BASE_VECTOR = (12, 61, 175, 310, 341, 217, 62)


def det(matrix: list[list[Fraction]]) -> Fraction:
    rows = [row[:] for row in matrix]
    result = Fraction(1)
    for column in range(len(rows)):
        pivot = next(
            (row for row in range(column, len(rows)) if rows[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            result = -result
        value = rows[column][column]
        result *= value
        for row in range(column + 1, len(rows)):
            scale = rows[row][column] / value
            for later in range(column + 1, len(rows)):
                rows[row][later] -= scale * rows[column][later]
    return result


def support_facets(
    labels: tuple[int, ...], coordinates: dict[int, tuple[Fraction, ...]]
) -> set[frozenset[int]]:
    supports: set[frozenset[int]] = set()
    for chosen in combinations(labels, 7):
        outside = tuple(vertex for vertex in labels if vertex not in chosen)
        values = []
        for vertex in outside:
            rows = [
                [Fraction(1), *coordinates[item]]
                for item in (*chosen, vertex)
            ]
            values.append(det(rows))
        nonzero = [value for value in values if value]
        if nonzero and (
            all(value > 0 for value in nonzero)
            or all(value < 0 for value in nonzero)
        ):
            supports.add(
                frozenset(chosen)
                | frozenset(
                    vertex
                    for vertex, value in zip(outside, values)
                    if not value
                )
            )
    maximal = {
        face for face in supports if not any(face < other for other in supports)
    }
    if any(len(face) != 7 for face in maximal):
        raise AssertionError(
            "hostile record-18 replay found a nonsimplicial supporting facet"
        )
    return maximal


def cyclic_facets(order: tuple[int, ...]) -> set[frozenset[int]]:
    result = set()
    for positions in combinations(range(10), 7):
        omitted = [position for position in range(10) if position not in positions]
        if all(
            sum(first < position < second for position in positions) % 2 == 0
            for first, second in combinations(omitted, 2)
        ):
            result.add(frozenset(order[position] for position in positions))
    return result


def stellar(
    facets: set[frozenset[int]], face: tuple[int, ...], vertex: int
) -> set[frozenset[int]]:
    carrier = frozenset(face)
    star = {facet for facet in facets if carrier <= facet}
    if not star or any(vertex in facet for facet in facets):
        raise AssertionError("hostile record-18 replay lost a stellar carrier")
    return (facets - star) | {
        (facet - {removed}) | {vertex}
        for facet in star
        for removed in carrier
    }


def radial_beyond(
    coordinates: dict[int, tuple[Fraction, ...]],
    labels: tuple[int, ...],
    face: tuple[int, ...],
    epsilon: Fraction,
) -> tuple[Fraction, ...]:
    barycentre = tuple(
        sum((coordinates[vertex][axis] for vertex in face), Fraction(0))
        / len(face)
        for axis in range(7)
    )
    interior = tuple(
        sum((coordinates[vertex][axis] for vertex in labels), Fraction(0))
        / len(labels)
        for axis in range(7)
    )
    return tuple(
        barycentre[axis] + epsilon * (barycentre[axis] - interior[axis])
        for axis in range(7)
    )


def all_faces(facets: set[frozenset[int]]) -> set[frozenset[int]]:
    result = {frozenset()}
    for facet in facets:
        for size in range(1, 8):
            result.update(frozenset(face) for face in combinations(facet, size))
    return result


def minimal_nonfaces(
    labels: tuple[int, ...], facets: set[frozenset[int]]
) -> tuple[frozenset[int], ...]:
    faces = all_faces(facets)
    result = []
    for size in range(1, len(labels) + 1):
        for chosen in combinations(labels, size):
            candidate = frozenset(chosen)
            if candidate not in faces and all(
                candidate - {vertex} in faces for vertex in candidate
            ):
                result.append(candidate)
    return tuple(sorted(result, key=lambda face: (len(face), tuple(sorted(face)))))


def main() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = dict(document)
    claimed = payload.pop("certificate_sha256")
    digest = sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if claimed != digest or digest != EXPECTED_DIGEST:
        raise AssertionError("hostile record-18 replay rejected the archive digest")
    proof = document["proof_source"]
    source_hash = sha256((ROOT / proof["path"]).read_bytes()).hexdigest()
    if (
        source_hash != proof["sha256"]
        or source_hash != EXPECTED_SOURCE_HASH
    ):
        raise AssertionError("hostile record-18 replay rejected the source hash")
    target = document["target"]
    if (
        tuple(map(tuple, target["missing_edges"])) != EXPECTED_EDGES
        or tuple(map(tuple, target["missing_fives"])) != EXPECTED_FIVES
    ):
        raise AssertionError("hostile record-18 replay target clutter changed")

    audit = document["audit"]
    order = tuple(audit["cyclic_order"])
    if order != (1, 0, 10, 5, 9, 2, 3, 4, 7, 11):
        raise AssertionError("hostile record-18 replay cyclic order changed")
    labels = tuple(sorted(int(label) for label in audit["base_coordinates"]))
    coordinates = {
        int(label): tuple(Fraction(value) for value in row)
        for label, row in audit["base_coordinates"].items()
    }
    for parameter, label in enumerate(order):
        expected = tuple(Fraction(parameter**power) for power in range(1, 8))
        if coordinates[label] != expected:
            raise AssertionError("hostile record-18 replay moment curve changed")

    subdivisions = audit["stellar_subdivisions"]
    facets = cyclic_facets(order)
    if len(facets) != 40:
        raise AssertionError("hostile record-18 replay cyclic facets changed")
    current_labels = order
    for step in subdivisions:
        face = tuple(step["face"])
        vertex = int(step["new_vertex"])
        epsilon = Fraction(step["radial_beyond_epsilon"])
        if radial_beyond(coordinates, current_labels, face, epsilon) != coordinates[vertex]:
            raise AssertionError("hostile record-18 replay radial placement changed")
        facets = stellar(facets, face, vertex)
        if len(facets) != int(step["facet_count"]):
            raise AssertionError("hostile record-18 replay stellar facets changed")
        current_labels = (*current_labels, vertex)

    geometric = support_facets(labels, coordinates)
    if geometric != facets or len(geometric) != 62:
        raise AssertionError("hostile record-18 coordinate supports differ")
    faces = all_faces(facets)
    vector = tuple(
        sum(len(face) == size for face in faces) for size in range(1, 8)
    )
    if vector != EXPECTED_BASE_VECTOR or list(vector) != audit["base_face_vector"]:
        raise AssertionError("hostile record-18 base face vector changed")
    nonfaces = minimal_nonfaces(labels, facets)
    archived_nonfaces = tuple(
        frozenset(face) for face in audit["base_minimal_nonfaces"]
    )
    if nonfaces != archived_nonfaces:
        raise AssertionError("hostile record-18 base minimal nonfaces changed")
    blocks = {int(key): tuple(value) for key, value in audit["clone_blocks"].items()}
    expanded = []
    for nonface in nonfaces:
        value: set[int] = set()
        for vertex in nonface:
            value.update(blocks.get(vertex, (vertex,)))
        expanded.append(frozenset(value))
    expanded_tuple = tuple(
        sorted(expanded, key=lambda face: (len(face), tuple(sorted(face))))
    )
    expected = tuple(
        sorted(
            (frozenset(face) for face in EXPECTED_EDGES + EXPECTED_FIVES),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    if (
        expanded_tuple != expected
        or [list(sorted(face)) for face in expanded_tuple]
        != audit["expanded_minimal_nonfaces"]
    ):
        raise AssertionError("hostile record-18 inverse multiwedge mismatch")
    print(
        json.dumps(
            {
                "verified": True,
                "certificate_sha256": digest,
                "base_facets": len(facets),
                "base_face_vector": vector,
            }
        )
    )


if __name__ == "__main__":
    main()
