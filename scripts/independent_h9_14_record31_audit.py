#!/usr/bin/env python3
"""Hostile stdlib replay of the record-31 ordinary-polytopality proof."""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record31_polytopality.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "12b01aec7f3de0bb4e71040291753e8362787b9c7063704877bd24693884a799"
)
EXPECTED_PROOF_SOURCE_SHA256 = (
    "ad30ec2b567c865165d805f1f2375211cf03aec766b574e8fe2fcf6cb987f890"
)
EXPECTED_SUPPORT_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_GRAPH6 = "M?????????O?O?o??"
EXPECTED_EDGES = ((0, 11), (1, 12), (0, 13), (1, 13))
EXPECTED_FIVES = (
    (0, 5, 6, 9, 10),
    (2, 3, 4, 7, 8),
    (2, 5, 6, 9, 10),
    (2, 5, 6, 11, 12),
    (2, 7, 8, 11, 12),
    (3, 4, 7, 8, 11),
    (3, 4, 7, 8, 13),
    (3, 4, 9, 10, 13),
    (5, 6, 9, 10, 13),
)
EXPECTED_BASE_VECTOR = (10, 41, 84, 85, 34)
EXPECTED_TARGET_VECTOR = (14, 87, 319, 771, 1287, 1491, 1158, 549, 122)
FULL_MASK = (1 << 14) - 1


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def det(matrix: list[list[Fraction]]) -> Fraction:
    rows = [row[:] for row in matrix]
    result = Fraction(1)
    for column in range(len(rows)):
        pivot = next(
            (row for row in range(column, len(rows)) if rows[row][column]), None
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
    for chosen in combinations(labels, 5):
        outside = tuple(vertex for vertex in labels if vertex not in chosen)
        values = []
        for vertex in outside:
            rows = [
                [Fraction(1), *coordinates[item]] for item in (*chosen, vertex)
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
    if any(len(face) != 5 for face in maximal):
        raise AssertionError("hostile record-31 replay found a nonsimplicial facet")
    return maximal


def cyclic_facets(order: tuple[int, ...]) -> set[frozenset[int]]:
    result = set()
    for positions in combinations(range(8), 5):
        omitted = tuple(position for position in range(8) if position not in positions)
        if all(
            sum(first < selected < second for selected in positions) % 2 == 0
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
        raise AssertionError("hostile record-31 replay lost a stellar carrier")
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
        for axis in range(5)
    )
    interior = tuple(
        sum((coordinates[vertex][axis] for vertex in labels), Fraction(0))
        / len(labels)
        for axis in range(5)
    )
    return tuple(
        barycentre[axis] + epsilon * (barycentre[axis] - interior[axis])
        for axis in range(5)
    )


def all_abstract_faces(
    facets: set[frozenset[int]], dimension: int
) -> set[frozenset[int]]:
    result = {frozenset()}
    for facet in facets:
        for size in range(1, dimension + 1):
            result.update(frozenset(face) for face in combinations(facet, size))
    return result


def minimal_nonfaces(
    labels: tuple[int, ...], faces: set[frozenset[int]]
) -> tuple[frozenset[int], ...]:
    result = []
    for size in range(1, len(labels) + 1):
        for chosen in combinations(labels, size):
            candidate = frozenset(chosen)
            if candidate not in faces and all(
                candidate - {vertex} in faces for vertex in candidate
            ):
                result.append(candidate)
    return tuple(sorted(result, key=lambda face: (len(face), tuple(sorted(face)))))


def mask(vertices: tuple[int, ...] | frozenset[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def blocker(clutter: tuple[int, ...]) -> tuple[int, ...]:
    def hits(candidate: int) -> bool:
        return all(candidate & nonface for nonface in clutter)

    result = []
    for candidate in range(1, 1 << 14):
        if not hits(candidate):
            continue
        if all(
            not hits(candidate ^ (1 << vertex))
            for vertex in range(14)
            if candidate & (1 << vertex)
        ):
            result.append(candidate)
    return tuple(result)


def target_face_table(facet_complements: tuple[int, ...]) -> tuple[bool, ...]:
    table = [False] * (1 << 14)
    table[0] = True
    for complement in facet_complements:
        facet = FULL_MASK ^ complement
        subset = facet
        while True:
            table[subset] = True
            if not subset:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def target_minimal_nonfaces(table: tuple[bool, ...]) -> tuple[int, ...]:
    result = []
    for candidate in range(1, 1 << 14):
        if table[candidate]:
            continue
        if all(
            table[candidate ^ (1 << vertex)]
            for vertex in range(14)
            if candidate & (1 << vertex)
        ):
            result.append(candidate)
    return tuple(result)


def ridge_and_connectivity(facets: tuple[int, ...]) -> tuple[Counter[int], int]:
    groups: dict[int, list[int]] = {}
    for complement in facets:
        outside = FULL_MASK ^ complement
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
    digest = sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if claimed != digest or digest != EXPECTED_CERTIFICATE_SHA256:
        raise AssertionError("hostile record-31 replay rejected archive digest")

    sources = document["sources"]
    proof_path = ROOT / sources["proof"]["path"]
    support_path = ROOT / sources["ridge_support"]["path"]
    if (
        file_hash(proof_path) != sources["proof"]["sha256"]
        or file_hash(proof_path) != EXPECTED_PROOF_SOURCE_SHA256
        or file_hash(support_path) != sources["ridge_support"]["sha256"]
        or file_hash(support_path) != EXPECTED_SUPPORT_SOURCE_SHA256
    ):
        raise AssertionError("hostile record-31 replay rejected source closure")

    target = document["target"]
    if (
        target["second_layer_survivor_index"] != 31
        or target["graph6"] != EXPECTED_GRAPH6
        or tuple(map(tuple, target["missing_edges"])) != EXPECTED_EDGES
        or tuple(map(tuple, target["missing_fives"])) != EXPECTED_FIVES
    ):
        raise AssertionError("hostile record-31 target changed")

    audit = document["polytopality_audit"]
    order = tuple(audit["cyclic_order"])
    if order != (2, 14, 15, 16, 17, 11, 13, 12):
        raise AssertionError("hostile record-31 cyclic order changed")
    labels = tuple(sorted(int(label) for label in audit["base_coordinates"]))
    coordinates = {
        int(label): tuple(Fraction(value) for value in row)
        for label, row in audit["base_coordinates"].items()
    }
    for parameter, label in enumerate(order):
        expected = tuple(Fraction(parameter**power) for power in range(1, 6))
        if coordinates[label] != expected:
            raise AssertionError("hostile record-31 moment curve changed")

    facets = cyclic_facets(order)
    if len(facets) != 20:
        raise AssertionError("hostile record-31 cyclic facets changed")
    current = order
    for step in audit["stellar_subdivisions"]:
        face = tuple(step["face"])
        vertex = int(step["new_vertex"])
        epsilon = Fraction(step["radial_beyond_epsilon"])
        if radial_beyond(coordinates, current, face, epsilon) != coordinates[vertex]:
            raise AssertionError("hostile record-31 radial placement changed")
        facets = stellar(facets, face, vertex)
        if len(facets) != int(step["facet_count"]):
            raise AssertionError("hostile record-31 stellar facets changed")
        current = (*current, vertex)

    geometric = support_facets(labels, coordinates)
    if geometric != facets or len(geometric) != 34:
        raise AssertionError("hostile record-31 coordinate supports differ")
    base_faces = all_abstract_faces(facets, 5)
    base_vector = tuple(
        sum(len(face) == size for face in base_faces) for size in range(1, 6)
    )
    if base_vector != EXPECTED_BASE_VECTOR or list(base_vector) != audit["base_face_vector"]:
        raise AssertionError("hostile record-31 base face vector changed")
    base_nonfaces = minimal_nonfaces(labels, base_faces)
    if [list(sorted(face)) for face in base_nonfaces] != audit["base_minimal_nonfaces"]:
        raise AssertionError("hostile record-31 base minimal nonfaces changed")

    blocks = {int(key): tuple(value) for key, value in audit["clone_blocks"].items()}
    expanded = []
    for nonface in base_nonfaces:
        value: set[int] = set()
        for vertex in nonface:
            value.update(blocks.get(vertex, (vertex,)))
        expanded.append(frozenset(value))
    expanded = sorted(expanded, key=lambda face: (len(face), tuple(sorted(face))))
    expected_nonfaces = sorted(
        (frozenset(face) for face in EXPECTED_EDGES + EXPECTED_FIVES),
        key=lambda face: (len(face), tuple(sorted(face))),
    )
    if expanded != expected_nonfaces:
        raise AssertionError("hostile record-31 multiwedge clutter differs")

    clutter = tuple(mask(face) for face in expected_nonfaces)
    facet_complements = blocker(clutter)
    if (
        Counter(value.bit_count() for value in facet_complements) != {5: 122}
        or list(facet_complements) != target["facet_complement_masks"]
    ):
        raise AssertionError("hostile record-31 blocker changed")
    table = target_face_table(facet_complements)
    vector = tuple(
        sum(table[mask(frozenset(face))] for face in combinations(range(14), size))
        for size in range(1, 10)
    )
    if vector != EXPECTED_TARGET_VECTOR or list(vector) != target["face_vector"]:
        raise AssertionError("hostile record-31 target face vector changed")
    if target_minimal_nonfaces(table) != tuple(sorted(clutter)):
        raise AssertionError("hostile record-31 target minimal nonfaces changed")
    ridge_histogram, connected = ridge_and_connectivity(facet_complements)
    if ridge_histogram != {2: 549} or connected != 122:
        raise AssertionError("hostile record-31 ridge/connectedness audit failed")

    support = json.loads(support_path.read_text(encoding="utf-8"))
    survivors = [record for record in support["records"] if record["survives_second_layer"]]
    source_record = survivors[31]
    if (
        source_record["graph6"] != EXPECTED_GRAPH6
        or not set(facet_complements) <= set(source_record["ridge_support_core_masks"])
        or source_record["dehn_sommerville_schedule"][0]["f_vector"]
        != list(EXPECTED_TARGET_VECTOR)
    ):
        raise AssertionError("hostile record-31 frozen-support audit failed")

    print(
        json.dumps(
            {
                "verified": True,
                "base_facets": 34,
                "target_facets": 122,
                "target_face_vector": list(vector),
                "certificate_sha256": digest,
                "ordinary_polytopal_only": True,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
