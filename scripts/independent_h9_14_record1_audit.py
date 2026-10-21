#!/usr/bin/env python3
"""Hostile stdlib replay of the record-1 ordinary-polytopality proof."""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record1_polytopality.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "5684d466a7bca45828da56aae3bd1ea91c5390e08f966c98bccf71b427f2b976"
)
EXPECTED_PROOF_SOURCE_SHA256 = (
    "d11365f9794a1b2fd450e08f60a5f334b4203b2e48a69ed44578555e3de9a6a3"
)
EXPECTED_SUPPORT_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_GRAPH6 = "M?????????????w??"
EXPECTED_EDGES = ((0, 13), (1, 13), (2, 13))
EXPECTED_FOURS = (
    (0, 1, 3, 7),
    (0, 2, 6, 11),
    (0, 3, 4, 9),
    (1, 7, 8, 12),
)
EXPECTED_FIVES = (
    (0, 3, 6, 7, 11),
    (1, 3, 5, 7, 10),
    (1, 4, 5, 9, 10),
    (1, 5, 8, 10, 12),
    (2, 4, 5, 9, 10),
    (2, 4, 6, 9, 11),
    (2, 5, 8, 10, 12),
    (2, 6, 8, 11, 12),
    (3, 4, 5, 9, 10),
    (6, 7, 8, 11, 12),
)
EXPECTED_BASE_VECTOR = (10, 42, 88, 90, 36)
EXPECTED_TARGET_VECTOR = (14, 88, 331, 831, 1451, 1757, 1414, 684, 152)
EXPECTED_BLOCKS = {1: (1, 10), 3: (3, 11), 4: (4, 12), 8: (8, 13)}
EXPECTED_RELABEL = {
    0: 3,
    1: 4,
    2: 0,
    3: 5,
    4: 6,
    5: 1,
    6: 2,
    7: 7,
    8: 8,
    9: 13,
    10: 9,
    11: 10,
    12: 11,
    13: 12,
}
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
        raise AssertionError("hostile record-1 replay found a nonsimplicial facet")
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
        raise AssertionError("hostile record-1 replay lost a stellar carrier")
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
        raise AssertionError("hostile record-1 replay rejected archive digest")

    sources = document["sources"]
    proof_path = ROOT / sources["proof"]["path"]
    support_path = ROOT / sources["ridge_support"]["path"]
    if (
        file_hash(proof_path) != sources["proof"]["sha256"]
        or file_hash(proof_path) != EXPECTED_PROOF_SOURCE_SHA256
        or file_hash(support_path) != sources["ridge_support"]["sha256"]
        or file_hash(support_path) != EXPECTED_SUPPORT_SOURCE_SHA256
    ):
        raise AssertionError("hostile record-1 replay rejected source closure")

    target = document["target"]
    if (
        target["second_layer_survivor_index"] != 1
        or target["graph6"] != EXPECTED_GRAPH6
        or tuple(map(tuple, target["missing_edges"])) != EXPECTED_EDGES
        or tuple(map(tuple, target["missing_fours"])) != EXPECTED_FOURS
        or tuple(map(tuple, target["missing_fives"])) != EXPECTED_FIVES
        or target["missing_four_schedule"] != 4
    ):
        raise AssertionError("hostile record-1 target changed")

    audit = document["polytopality_audit"]
    order = tuple(audit["cyclic_order"])
    if order != tuple(range(8)):
        raise AssertionError("hostile record-1 cyclic order changed")
    labels = tuple(sorted(int(label) for label in audit["base_coordinates"]))
    if labels != tuple(range(10)):
        raise AssertionError("hostile record-1 base labels changed")
    coordinates = {
        int(label): tuple(Fraction(value) for value in row)
        for label, row in audit["base_coordinates"].items()
    }
    for parameter, label in enumerate(order):
        expected = tuple(Fraction(parameter**power) for power in range(1, 6))
        if coordinates[label] != expected:
            raise AssertionError("hostile record-1 moment curve changed")

    facets = cyclic_facets(order)
    if len(facets) != 20:
        raise AssertionError("hostile record-1 cyclic facets changed")
    current = order
    expected_steps = (
        ((0, 1, 2), 8, Fraction(1, 16), 30),
        ((0, 1, 3), 9, Fraction(1, 16), 36),
    )
    if len(audit["stellar_subdivisions"]) != len(expected_steps):
        raise AssertionError("hostile record-1 stellar schedule length changed")
    for step, expected_step in zip(audit["stellar_subdivisions"], expected_steps):
        face = tuple(step["face"])
        vertex = int(step["new_vertex"])
        epsilon = Fraction(step["radial_beyond_epsilon"])
        count = int(step["facet_count"])
        if (face, vertex, epsilon, count) != expected_step:
            raise AssertionError("hostile record-1 stellar schedule changed")
        if radial_beyond(coordinates, current, face, epsilon) != coordinates[vertex]:
            raise AssertionError("hostile record-1 radial placement changed")
        facets = stellar(facets, face, vertex)
        if len(facets) != count:
            raise AssertionError("hostile record-1 stellar facets changed")
        current = (*current, vertex)

    geometric = support_facets(labels, coordinates)
    if geometric != facets or len(geometric) != 36:
        raise AssertionError("hostile record-1 coordinate supports differ")
    base_faces = all_abstract_faces(facets, 5)
    base_vector = tuple(
        sum(len(face) == size for face in base_faces) for size in range(1, 6)
    )
    if base_vector != EXPECTED_BASE_VECTOR or list(base_vector) != audit["base_face_vector"]:
        raise AssertionError("hostile record-1 base face vector changed")
    base_nonfaces = minimal_nonfaces(labels, base_faces)
    if [list(sorted(face)) for face in base_nonfaces] != audit["base_minimal_nonfaces"]:
        raise AssertionError("hostile record-1 base minimal nonfaces changed")

    blocks = {
        int(key): tuple(value)
        for key, value in audit["construction_clone_blocks"].items()
    }
    relabel = {
        int(key): int(value)
        for key, value in audit["construction_to_target"].items()
    }
    if blocks != EXPECTED_BLOCKS or relabel != EXPECTED_RELABEL:
        raise AssertionError("hostile record-1 wedge data changed")
    expanded = []
    for nonface in base_nonfaces:
        construction: set[int] = set()
        for vertex in nonface:
            construction.update(blocks.get(vertex, (vertex,)))
        expanded.append(frozenset(relabel[vertex] for vertex in construction))
    expanded = sorted(expanded, key=lambda face: (len(face), tuple(sorted(face))))
    expected_nonfaces = sorted(
        (
            frozenset(face)
            for face in EXPECTED_EDGES + EXPECTED_FOURS + EXPECTED_FIVES
        ),
        key=lambda face: (len(face), tuple(sorted(face))),
    )
    if expanded != expected_nonfaces:
        raise AssertionError("hostile record-1 multiwedge clutter differs")

    clutter = tuple(mask(face) for face in expected_nonfaces)
    facet_complements = blocker(clutter)
    if (
        Counter(value.bit_count() for value in facet_complements) != {5: 152}
        or list(facet_complements) != target["facet_complement_masks"]
    ):
        raise AssertionError("hostile record-1 blocker changed")
    table = target_face_table(facet_complements)
    vector = tuple(
        sum(table[mask(frozenset(face))] for face in combinations(range(14), size))
        for size in range(1, 10)
    )
    if vector != EXPECTED_TARGET_VECTOR or list(vector) != target["face_vector"]:
        raise AssertionError("hostile record-1 target face vector changed")
    if target_minimal_nonfaces(table) != tuple(sorted(clutter)):
        raise AssertionError("hostile record-1 target minimal nonfaces changed")
    ridge_histogram, connected = ridge_and_connectivity(facet_complements)
    if ridge_histogram != {2: 684} or connected != 152:
        raise AssertionError("hostile record-1 ridge/connectedness audit failed")

    support = json.loads(support_path.read_text(encoding="utf-8"))
    survivors = [record for record in support["records"] if record["survives_second_layer"]]
    source_record = survivors[1]
    if (
        source_record["graph6"] != EXPECTED_GRAPH6
        or not set(facet_complements) <= set(source_record["ridge_support_core_masks"])
        or source_record["dehn_sommerville_schedule"][4]["f_vector"]
        != list(EXPECTED_TARGET_VECTOR)
        or source_record["dehn_sommerville_schedule"][4]["missing_four_count"] != 4
    ):
        raise AssertionError("hostile record-1 frozen-support audit failed")

    print(
        json.dumps(
            {
                "verified": True,
                "base_facets": 36,
                "target_facets": 152,
                "target_face_vector": list(vector),
                "certificate_sha256": digest,
                "ordinary_polytopal_only": True,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
