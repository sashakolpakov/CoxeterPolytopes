#!/usr/bin/env python3
"""Hostile standard-library replay of record-3 ordinary polytopality."""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record3_polytopality.json"
EXPECTED_CERTIFICATE_SHA256 = "500c0874643b21cc5f6c0eadaf278695d9516db1832a895afb84042fa1da4685"
EXPECTED_PROOF_SOURCE_SHA256 = "c56890b3049f5aa379e6a9caf2033a9558c09699647bdfb741916657882820a7"
EXPECTED_SUPPORT_SOURCE_SHA256 = "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
EXPECTED_GRAPH6 = "M???????????_?O??"
EXPECTED_EDGES = ((0, 12), (1, 13))
EXPECTED_FOURS = (
    (1, 5, 10, 12), (3, 4, 9, 12), (3, 4, 9, 13), (5, 7, 10, 11),
)
EXPECTED_FIVES = (
    (0, 1, 2, 6, 8),
    (0, 2, 4, 6, 8),
    (0, 3, 4, 6, 9),
    (1, 2, 5, 8, 10),
    (2, 4, 6, 8, 13),
    (2, 6, 7, 8, 11),
    (3, 4, 7, 9, 11),
    (3, 5, 9, 10, 12),
    (3, 5, 9, 10, 13),
)
EXPECTED_BASE_VECTOR = (10, 42, 88, 90, 36)
EXPECTED_TARGET_VECTOR = (14, 89, 340, 866, 1528, 1862, 1504, 729, 162)
FULL = (1 << 14) - 1


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decode_graph6(value: str) -> tuple[int, tuple[tuple[int, int], ...]]:
    codes = [ord(character) - 63 for character in value]
    order = codes[0]
    count = comb(order, 2)
    bits = [(code >> shift) & 1 for code in codes[1:] for shift in range(5, -1, -1)]
    if len(codes) != 1 + (count + 5) // 6 or any(bits[count:]):
        raise AssertionError("hostile record-3 graph6 syntax failed")
    edges = []
    index = 0
    for upper in range(1, order):
        for lower in range(upper):
            if bits[index]:
                edges.append((lower, upper))
            index += 1
    return order, tuple(edges)


def determinant(matrix: list[list[Fraction]]) -> Fraction:
    rows = [row[:] for row in matrix]
    result = Fraction(1)
    for column in range(len(rows)):
        pivot = next((row for row in range(column, len(rows)) if rows[row][column]), None)
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


def supporting_facets(labels: tuple[int, ...], coordinates: dict[int, tuple[Fraction, ...]]) -> set[frozenset[int]]:
    supports: set[frozenset[int]] = set()
    for chosen in combinations(labels, 5):
        outside = tuple(label for label in labels if label not in chosen)
        values = [
            determinant([[Fraction(1), *coordinates[item]] for item in (*chosen, label)])
            for label in outside
        ]
        nonzero = [value for value in values if value]
        if nonzero and (all(value > 0 for value in nonzero) or all(value < 0 for value in nonzero)):
            supports.add(
                frozenset(chosen)
                | frozenset(label for label, value in zip(outside, values) if not value)
            )
    maximal = {face for face in supports if not any(face < other for other in supports)}
    if any(len(face) != 5 for face in maximal):
        raise AssertionError("hostile record-3 base is nonsimplicial")
    return maximal


def cyclic_facets() -> set[frozenset[int]]:
    result = set()
    for chosen in combinations(range(8), 5):
        omitted = tuple(position for position in range(8) if position not in chosen)
        if all(
            sum(first < selected < second for selected in chosen) % 2 == 0
            for first, second in combinations(omitted, 2)
        ):
            result.add(frozenset(chosen))
    return result


def stellar(facets: set[frozenset[int]], face: tuple[int, ...], vertex: int) -> set[frozenset[int]]:
    carrier = frozenset(face)
    star = {facet for facet in facets if carrier <= facet}
    if not star:
        raise AssertionError("hostile record-3 stellar carrier missing")
    return (facets - star) | {
        (facet - {removed}) | {vertex} for facet in star for removed in carrier
    }


def radial(coordinates: dict[int, tuple[Fraction, ...]], labels: tuple[int, ...], face: tuple[int, ...], epsilon: Fraction) -> tuple[Fraction, ...]:
    barycentre = tuple(
        sum((coordinates[vertex][axis] for vertex in face), Fraction(0)) / len(face)
        for axis in range(5)
    )
    interior = tuple(
        sum((coordinates[vertex][axis] for vertex in labels), Fraction(0)) / len(labels)
        for axis in range(5)
    )
    return tuple(
        barycentre[axis] + epsilon * (barycentre[axis] - interior[axis])
        for axis in range(5)
    )


def all_faces(facets: set[frozenset[int]]) -> set[frozenset[int]]:
    result = {frozenset()}
    for facet in facets:
        for size in range(1, 6):
            result.update(frozenset(face) for face in combinations(facet, size))
    return result


def minimal_nonfaces(labels: tuple[int, ...], faces: set[frozenset[int]]) -> tuple[frozenset[int], ...]:
    result = []
    for size in range(1, len(labels) + 1):
        for chosen in combinations(labels, size):
            candidate = frozenset(chosen)
            if candidate not in faces and all(candidate - {vertex} in faces for vertex in candidate):
                result.append(candidate)
    return tuple(sorted(result, key=lambda face: (len(face), tuple(sorted(face)))))


def mask(items: frozenset[int] | tuple[int, ...]) -> int:
    return sum(1 << vertex for vertex in items)


def blocker(clutter: tuple[int, ...]) -> tuple[int, ...]:
    def hits(candidate: int) -> bool:
        return all(candidate & row for row in clutter)
    return tuple(
        candidate for candidate in range(1, 1 << 14)
        if hits(candidate)
        and all(
            not hits(candidate ^ (1 << vertex))
            for vertex in range(14) if candidate & (1 << vertex)
        )
    )


def face_table(facets: tuple[int, ...]) -> tuple[bool, ...]:
    table = [False] * (1 << 14)
    table[0] = True
    for complement in facets:
        facet = FULL ^ complement
        subset = facet
        while True:
            table[subset] = True
            if not subset:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def table_minimal_nonfaces(table: tuple[bool, ...]) -> tuple[int, ...]:
    return tuple(
        candidate for candidate in range(1, 1 << 14)
        if not table[candidate]
        and all(
            table[candidate ^ (1 << vertex)]
            for vertex in range(14) if candidate & (1 << vertex)
        )
    )


def ridge_audit(facets: tuple[int, ...]) -> tuple[Counter[int], int]:
    groups: dict[int, list[int]] = {}
    for complement in facets:
        outside = FULL ^ complement
        while outside:
            bit = outside & -outside
            groups.setdefault(complement | bit, []).append(complement)
            outside ^= bit
    histogram = Counter(len(group) for group in groups.values())
    adjacency = {facet: set() for facet in facets}
    for group in groups.values():
        if len(group) == 2:
            first, second = group
            adjacency[first].add(second)
            adjacency[second].add(first)
    seen = {facets[0]}
    queue = deque([facets[0]])
    while queue:
        current = queue.popleft()
        for neighbor in adjacency[current] - seen:
            seen.add(neighbor)
            queue.append(neighbor)
    return histogram, len(seen)


def main() -> None:
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = dict(document)
    claimed = payload.pop("certificate_sha256")
    digest = sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if claimed != digest or digest != EXPECTED_CERTIFICATE_SHA256:
        raise AssertionError("hostile record-3 archive digest failed")
    sources = document["sources"]
    proof_path = ROOT / sources["proof"]["path"]
    support_path = ROOT / sources["ridge_support"]["path"]
    if (
        file_hash(proof_path) != EXPECTED_PROOF_SOURCE_SHA256
        or file_hash(proof_path) != sources["proof"]["sha256"]
        or file_hash(support_path) != EXPECTED_SUPPORT_SOURCE_SHA256
        or file_hash(support_path) != sources["ridge_support"]["sha256"]
    ):
        raise AssertionError("hostile record-3 source closure failed")

    target = document["target"]
    if (
        target["second_layer_survivor_index"] != 3
        or target["graph6"] != EXPECTED_GRAPH6
        or tuple(map(tuple, target["missing_edges"])) != EXPECTED_EDGES
        or tuple(map(tuple, target["missing_fours"])) != EXPECTED_FOURS
        or tuple(map(tuple, target["missing_fives"])) != EXPECTED_FIVES
        or decode_graph6(EXPECTED_GRAPH6) != (14, EXPECTED_EDGES)
    ):
        raise AssertionError("hostile record-3 target failed")

    audit = document["polytopality_audit"]
    coordinates = {
        int(label): tuple(Fraction(value) for value in row)
        for label, row in audit["base_coordinates"].items()
    }
    for parameter in range(8):
        expected = tuple(Fraction(parameter**power) for power in range(1, 6))
        if coordinates[parameter] != expected:
            raise AssertionError("hostile record-3 moment curve failed")
    facets = cyclic_facets()
    current = tuple(range(8))
    for step, expected_count in zip(audit["stellar_subdivisions"], (28, 36)):
        carrier = tuple(step["face"])
        vertex = int(step["new_vertex"])
        epsilon = Fraction(step["radial_beyond_epsilon"])
        if radial(coordinates, current, carrier, epsilon) != coordinates[vertex]:
            raise AssertionError("hostile record-3 radial formula failed")
        facets = stellar(facets, carrier, vertex)
        if len(facets) != expected_count or len(facets) != step["facet_count"]:
            raise AssertionError("hostile record-3 stellar count failed")
        current = (*current, vertex)
    if supporting_facets(tuple(range(10)), coordinates) != facets:
        raise AssertionError("hostile record-3 support determinants failed")
    base_faces = all_faces(facets)
    base_vector = tuple(sum(len(face) == size for face in base_faces) for size in range(1, 6))
    if base_vector != EXPECTED_BASE_VECTOR or list(base_vector) != audit["base_face_vector"]:
        raise AssertionError("hostile record-3 base vector failed")
    base_nonfaces = minimal_nonfaces(tuple(range(10)), base_faces)
    if [list(sorted(face)) for face in base_nonfaces] != audit["base_minimal_nonfaces"]:
        raise AssertionError("hostile record-3 base clutter failed")

    blocks = {int(key): tuple(value) for key, value in audit["base_blocks"].items()}
    if sorted(sum((list(block) for block in blocks.values()), [])) != list(range(14)):
        raise AssertionError("hostile record-3 block partition failed")
    expanded = sorted(
        (
            frozenset(vertex for base in nonface for vertex in blocks[base])
            for nonface in base_nonfaces
        ),
        key=lambda face: (len(face), tuple(sorted(face))),
    )
    expected_nonfaces = sorted(
        (frozenset(face) for face in EXPECTED_EDGES + EXPECTED_FOURS + EXPECTED_FIVES),
        key=lambda face: (len(face), tuple(sorted(face))),
    )
    if expanded != expected_nonfaces:
        raise AssertionError("hostile record-3 expanded clutter failed")
    clutter = tuple(mask(face) for face in expected_nonfaces)
    facet_complements = blocker(clutter)
    if Counter(value.bit_count() for value in facet_complements) != {5: 162}:
        raise AssertionError("hostile record-3 blocker purity failed")
    if list(facet_complements) != target["facet_complement_masks"]:
        raise AssertionError("hostile record-3 blocker list failed")
    table = face_table(facet_complements)
    vector = tuple(
        sum(table[mask(frozenset(face))] for face in combinations(range(14), size))
        for size in range(1, 10)
    )
    if vector != EXPECTED_TARGET_VECTOR or list(vector) != target["face_vector"]:
        raise AssertionError("hostile record-3 target vector failed")
    if table_minimal_nonfaces(table) != tuple(sorted(clutter)):
        raise AssertionError("hostile record-3 target clutter failed")
    ridge_histogram, connected = ridge_audit(facet_complements)
    if ridge_histogram != {2: 729} or connected != 162:
        raise AssertionError("hostile record-3 ridge audit failed")

    support = json.loads(support_path.read_text(encoding="utf-8"))
    record = [item for item in support["records"] if item["survives_second_layer"]][3]
    if (
        record["graph6"] != EXPECTED_GRAPH6
        or not set(facet_complements) <= set(record["ridge_support_core_masks"])
        or record["dehn_sommerville_schedule"][4]["f_vector"] != list(EXPECTED_TARGET_VECTOR)
    ):
        raise AssertionError("hostile record-3 support archive failed")
    print(json.dumps({
        "verified": True,
        "base_facets": 36,
        "target_facets": 162,
        "target_face_vector": list(vector),
        "certificate_sha256": digest,
        "ordinary_polytopal_only": True,
    }, sort_keys=True))


if __name__ == "__main__":
    main()
