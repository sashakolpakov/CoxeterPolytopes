#!/usr/bin/env python3
"""Hostile stdlib replay of the record-17 rational-polytopality archive.

This script deliberately does not import the proof implementation.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record17_polytopality.json"
EXPECTED_DIGEST = "d17332be82e6bddb0ebfe320d678238cc0fa182f2a1b03d81931217826c72e31"
EXPECTED_SOURCE_HASH = "00722328898ff4b0fe7cdf1efade2264300680cf8abccdbbc55388215efd9aae"
EXPECTED_EDGES = ((0, 12), (1, 12), (0, 13), (2, 13), (12, 13))
EXPECTED_FIVES = (
    (0, 3, 4, 5, 6), (1, 3, 4, 5, 6), (2, 3, 4, 5, 6),
    (0, 5, 6, 7, 8), (0, 7, 8, 9, 10), (1, 2, 3, 4, 11),
    (1, 2, 9, 10, 11), (7, 8, 9, 10, 11),
)


def det(matrix: list[list[Fraction]]) -> Fraction:
    rows = [row[:] for row in matrix]
    result = Fraction(1)
    for column in range(len(rows)):
        pivot = next((r for r in range(column, len(rows)) if rows[r][column]), None)
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


def support_facets(labels: tuple[int, ...], coordinates: dict[int, tuple[Fraction, ...]]) -> set[frozenset[int]]:
    supports: set[frozenset[int]] = set()
    for chosen in combinations(labels, 5):
        outside = tuple(v for v in labels if v not in chosen)
        values = []
        for vertex in outside:
            rows = [[Fraction(1), *coordinates[v]] for v in (*chosen, vertex)]
            values.append(det(rows))
        nonzero = [value for value in values if value]
        if nonzero and (all(value > 0 for value in nonzero) or all(value < 0 for value in nonzero)):
            supports.add(frozenset(chosen) | frozenset(v for v, value in zip(outside, values) if not value))
    maximal = {face for face in supports if not any(face < other for other in supports)}
    if any(len(face) != 5 for face in maximal):
        raise AssertionError("hostile replay found a nonsimplicial supporting facet")
    return maximal


def cyclic_facets(order: tuple[int, ...]) -> set[frozenset[int]]:
    result = set()
    for positions in combinations(range(8), 5):
        omitted = [p for p in range(8) if p not in positions]
        if all(sum(a < p < b for p in positions) % 2 == 0 for a, b in combinations(omitted, 2)):
            result.add(frozenset(order[p] for p in positions))
    return result


def stellar(facets: set[frozenset[int]], face: tuple[int, ...], vertex: int) -> set[frozenset[int]]:
    carrier = frozenset(face)
    star = {facet for facet in facets if carrier <= facet}
    if not star:
        raise AssertionError("hostile replay lost a stellar carrier")
    return (facets - star) | {(facet - {removed}) | {vertex} for facet in star for removed in carrier}


def all_faces(facets: set[frozenset[int]]) -> set[frozenset[int]]:
    result = {frozenset()}
    for facet in facets:
        for size in range(1, 6):
            result.update(frozenset(face) for face in combinations(facet, size))
    return result


def minimal_nonfaces(labels: tuple[int, ...], facets: set[frozenset[int]]) -> tuple[frozenset[int], ...]:
    faces = all_faces(facets)
    result = []
    for size in range(1, len(labels) + 1):
        for chosen in combinations(labels, size):
            candidate = frozenset(chosen)
            if candidate not in faces and all(candidate - {v} in faces for v in candidate):
                result.append(candidate)
    return tuple(sorted(result, key=lambda face: (len(face), tuple(sorted(face)))))


def main() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = dict(document)
    claimed = payload.pop("certificate_sha256")
    digest = sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if claimed != digest or digest != EXPECTED_DIGEST:
        raise AssertionError("hostile replay rejected the archive digest")
    proof = document["proof_source"]
    source_hash = sha256((ROOT / proof["path"]).read_bytes()).hexdigest()
    if source_hash != proof["sha256"] or source_hash != EXPECTED_SOURCE_HASH:
        raise AssertionError("hostile replay rejected the proof-source hash")
    target = document["target"]
    if tuple(map(tuple, target["missing_edges"])) != EXPECTED_EDGES or tuple(map(tuple, target["missing_fives"])) != EXPECTED_FIVES:
        raise AssertionError("hostile replay target clutter changed")

    audit = document["audit"]
    order = tuple(audit["cyclic_order"])
    subdivisions = audit["stellar_subdivisions"]
    facets = cyclic_facets(order)
    if len(facets) != 20:
        raise AssertionError("hostile replay cyclic facet count changed")
    for step in subdivisions:
        facets = stellar(facets, tuple(step["face"]), int(step["new_vertex"]))
        if len(facets) != int(step["facet_count"]):
            raise AssertionError("hostile replay stellar facet count changed")
    labels = tuple(sorted(int(label) for label in audit["base_coordinates"]))
    coordinates = {
        int(label): tuple(Fraction(value) for value in row)
        for label, row in audit["base_coordinates"].items()
    }
    geometric = support_facets(labels, coordinates)
    if geometric != facets or len(geometric) != 32:
        raise AssertionError("hostile replay coordinate supports differ")
    faces = all_faces(facets)
    vector = tuple(sum(len(face) == size for face in faces) for size in range(1, 6))
    if vector != (10, 40, 80, 80, 32) or list(vector) != audit["base_face_vector"]:
        raise AssertionError("hostile replay base face vector changed")
    nonfaces = minimal_nonfaces(labels, facets)
    archived_nonfaces = tuple(frozenset(face) for face in audit["base_minimal_nonfaces"])
    if nonfaces != archived_nonfaces:
        raise AssertionError("hostile replay base minimal nonfaces changed")
    blocks = {int(k): tuple(v) for k, v in audit["clone_blocks"].items()}
    expanded = []
    for nonface in nonfaces:
        value: set[int] = set()
        for vertex in nonface:
            value.update(blocks.get(vertex, (vertex,)))
        expanded.append(frozenset(value))
    expanded = tuple(sorted(expanded, key=lambda face: (len(face), tuple(sorted(face)))))
    expected = tuple(sorted((frozenset(face) for face in EXPECTED_EDGES + EXPECTED_FIVES), key=lambda face: (len(face), tuple(sorted(face)))))
    if expanded != expected or [list(sorted(face)) for face in expanded] != audit["expanded_minimal_nonfaces"]:
        raise AssertionError("hostile replay inverse multiwedge mismatch")
    print(json.dumps({"verified": True, "certificate_sha256": digest, "base_facets": len(facets), "base_face_vector": vector}))


if __name__ == "__main__":
    main()

