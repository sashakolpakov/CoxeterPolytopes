#!/usr/bin/env python3
"""Hostile replay of the record-4 rational ordinary-polytopality proof.

No record-4 construction or certificate code is imported.  This script
independently reconstructs both rational stellar subdivisions, the four
wedges, the complete blocker, and every global face/link audit.
"""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
    generated_face_table,
)


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record4_polytopality.json"
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
GRAPH6 = "M???????????_?o??"
EDGES = ((0, 12), (0, 13), (1, 13))
FOURS = ((1, 2, 4, 9), (1, 2, 5, 10))
FIVES = (
    (1, 3, 5, 8, 10),
    (2, 4, 6, 9, 11),
    (2, 4, 7, 9, 12),
    (2, 5, 7, 10, 12),
    (3, 5, 7, 8, 10),
    (3, 5, 8, 10, 13),
    (3, 6, 7, 8, 11),
    (3, 6, 8, 11, 13),
    (4, 6, 7, 9, 11),
    (4, 6, 9, 11, 13),
)
EXPECTED_BASE_VECTOR = (10, 42, 88, 90, 36)
EXPECTED_TARGET_VECTOR = (14, 88, 330, 823, 1424, 1708, 1364, 657, 146)
EXPECTED_FACET_SHA256 = "1bbd853f610d450b5e379ed9de03d97ee53ac3718da94fb585aca7b6192115a9"
EXPECTED_FACE_SHA256 = "ce13e1bcb9207ce4e18e50dc943148d2767d9e8f9e9227bd3cdd36b48fa518d0"
EXPECTED_SUMMARY_SHA256 = "1d34dcf49e43151af6ba22dde78593468b76ce459cdcaca3dd6f46e51a7ba687"
FULL = (1 << 14) - 1

COORDINATES = {
    **{
        label: tuple(Fraction(label**power) for power in range(1, 6))
        for label in range(8)
    },
    8: (
        Fraction(5, 16), Fraction(-9, 16), Fraction(-179, 32),
        Fraction(-36), Fraction(-7235, 32),
    ),
    9: (
        Fraction(2705, 768), Fraction(57745, 2304), Fraction(271585, 1536),
        Fraction(358073, 288), Fraction(13407505, 1536),
    ),
}
CLONES = {3: (3, 10), 4: (4, 11), 5: (5, 12), 6: (6, 13)}
RELABEL = {
    0: 0, 1: 13, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 12,
    8: 7, 9: 1, 10: 8, 11: 9, 12: 10, 13: 11,
}


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def msk(values) -> int:
    return sum(1 << value for value in values)


def verts(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if value & (1 << vertex))


def decode(word: str) -> tuple[int, tuple[tuple[int, int], ...]]:
    codes = [ord(character) - 63 for character in word]
    order = codes[0]
    count = order * (order - 1) // 2
    bits = [
        (code >> shift) & 1
        for code in codes[1:]
        for shift in range(5, -1, -1)
    ]
    if len(codes) != 1 + (count + 5) // 6 or any(bits[count:]):
        raise AssertionError("hostile record-4 graph6 framing failed")
    edges = []
    position = 0
    for upper in range(1, order):
        for lower in range(upper):
            if bits[position]:
                edges.append((lower, upper))
            position += 1
    return order, tuple(edges)


def determinant(matrix: list[list[Fraction]]) -> Fraction:
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


def supporting_facets(labels: tuple[int, ...]) -> set[frozenset[int]]:
    supports: set[frozenset[int]] = set()
    for chosen in combinations(labels, 5):
        outside = tuple(vertex for vertex in labels if vertex not in chosen)
        values = [
            determinant(
                [[Fraction(1), *COORDINATES[item]] for item in (*chosen, vertex)]
            )
            for vertex in outside
        ]
        nonzero = [value for value in values if value]
        if nonzero and (
            all(value > 0 for value in nonzero)
            or all(value < 0 for value in nonzero)
        ):
            supports.add(
                frozenset(chosen)
                | frozenset(
                    vertex for vertex, value in zip(outside, values) if not value
                )
            )
    maximal = {face for face in supports if not any(face < other for other in supports)}
    if any(len(face) != 5 for face in maximal):
        raise AssertionError("hostile record-4 rational hull is nonsimplicial")
    return maximal


def cyclic_facets() -> set[frozenset[int]]:
    result = set()
    for positions in combinations(range(8), 5):
        omitted = tuple(position for position in range(8) if position not in positions)
        if all(
            sum(first < selected < second for selected in positions) % 2 == 0
            for first, second in combinations(omitted, 2)
        ):
            result.add(frozenset(positions))
    return result


def stellar(
    facets: set[frozenset[int]], carrier: tuple[int, ...], vertex: int
) -> set[frozenset[int]]:
    face = frozenset(carrier)
    star = {facet for facet in facets if face <= facet}
    if not star or any(vertex in facet for facet in facets):
        raise AssertionError("hostile record-4 stellar carrier failed")
    return (facets - star) | {
        (facet - {removed}) | {vertex}
        for facet in star
        for removed in face
    }


def radial(
    labels: tuple[int, ...], carrier: tuple[int, ...], epsilon: Fraction
) -> tuple[Fraction, ...]:
    barycentre = tuple(
        sum((COORDINATES[vertex][axis] for vertex in carrier), Fraction())
        / len(carrier)
        for axis in range(5)
    )
    interior = tuple(
        sum((COORDINATES[vertex][axis] for vertex in labels), Fraction())
        / len(labels)
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
            result.update(frozenset(value) for value in combinations(facet, size))
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


def blocker(clutter: tuple[int, ...]) -> tuple[int, ...]:
    def hits(candidate: int) -> bool:
        return all(candidate & edge for edge in clutter)

    return tuple(
        candidate
        for candidate in range(1, 1 << 14)
        if hits(candidate)
        and all(not hits(candidate ^ (1 << vertex)) for vertex in verts(candidate))
    )


def ridge_replay(facets: tuple[int, ...]) -> tuple[Counter[int], int]:
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
    queue = deque(seen)
    while queue:
        current = queue.popleft()
        for neighbor in adjacency[current] - seen:
            seen.add(neighbor)
            queue.append(neighbor)
    return histogram, len(seen)


def packed_face_digest(facets: tuple[int, ...]) -> str:
    table = generated_face_table(facets)
    packed = bytes(
        sum(bool(table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(table), 8)
    )
    return sha256(packed).hexdigest()


def replay() -> dict[str, object]:
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = dict(archive)
    claimed = payload.pop("certificate_sha256")
    observed = sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if claimed != observed:
        raise AssertionError("hostile record-4 archive digest failed")
    if archive.get("theorem_status") != "certified_explicit_ordinary_polytopal":
        raise AssertionError("hostile record-4 theorem status changed")
    if archive["sources"]["hostile_replay"]["sha256"] != file_hash(Path(__file__)):
        raise AssertionError("hostile record-4 checker is not source-pinned")
    if decode(GRAPH6) != (14, EDGES):
        raise AssertionError("hostile record-4 graph6 decode changed")

    cyclic = cyclic_facets()
    first = stellar(cyclic, (0, 1), 8)
    second = stellar(first, (0, 7), 9)
    if (len(cyclic), len(first), len(second)) != (20, 28, 36):
        raise AssertionError("hostile record-4 stellar counts changed")
    if radial(tuple(range(8)), (0, 1), Fraction(1, 16)) != COORDINATES[8]:
        raise AssertionError("hostile record-4 first radial identity failed")
    if radial(tuple(range(9)), (0, 7), Fraction(1, 16)) != COORDINATES[9]:
        raise AssertionError("hostile record-4 second radial identity failed")
    if supporting_facets(tuple(range(9))) != first:
        raise AssertionError("hostile record-4 first rational hull mismatch")
    if supporting_facets(tuple(range(10))) != second:
        raise AssertionError("hostile record-4 second rational hull mismatch")

    faces = all_faces(second)
    base_vector = tuple(sum(len(face) == size for face in faces) for size in range(1, 6))
    if base_vector != EXPECTED_BASE_VECTOR:
        raise AssertionError("hostile record-4 base vector changed")
    base_nonfaces = minimal_nonfaces(tuple(range(10)), faces)
    expanded = []
    for nonface in base_nonfaces:
        value = set()
        for vertex in nonface:
            value.update(CLONES.get(vertex, (vertex,)))
        expanded.append(frozenset(RELABEL[vertex] for vertex in value))
    target = tuple(
        sorted(
            (frozenset(value) for value in EDGES + FOURS + FIVES),
            key=lambda face: (len(face), tuple(sorted(face))),
        )
    )
    if tuple(sorted(expanded, key=lambda face: (len(face), tuple(sorted(face))))) != target:
        raise AssertionError("hostile record-4 wedge expansion mismatch")

    facets = blocker(tuple(msk(value) for value in EDGES + FOURS + FIVES))
    if Counter(value.bit_count() for value in facets) != {5: 146}:
        raise AssertionError("hostile record-4 blocker is not pure rank five")
    if tuple(archive["target"]["facet_complement_masks"]) != facets:
        raise AssertionError("hostile record-4 blocker differs from archive")
    facet_digest = sha256(
        b"".join(value.to_bytes(2, "little") for value in facets)
    ).hexdigest()
    if EXPECTED_FACET_SHA256 != "TO_BE_FROZEN" and facet_digest != EXPECTED_FACET_SHA256:
        raise AssertionError("hostile record-4 facet digest changed")
    ridge_histogram, connected = ridge_replay(facets)
    if ridge_histogram != {2: 657} or connected != 146:
        raise AssertionError("hostile record-4 ridge replay failed")

    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    source = [record for record in support["records"] if record["survives_second_layer"]][4]
    if source["graph6"] != GRAPH6:
        raise AssertionError("hostile record-4 support index changed")
    audit = audit_global_facet_selection(
        GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        facets,
    )
    if not audit.passes_connected_homology_audit or audit.face_vector != EXPECTED_TARGET_VECTOR:
        raise AssertionError("hostile record-4 full face/link audit failed")
    face_digest = packed_face_digest(facets)
    if EXPECTED_FACE_SHA256 != "TO_BE_FROZEN" and face_digest != EXPECTED_FACE_SHA256:
        raise AssertionError("hostile record-4 face digest changed")

    summary = {
        "certificate_sha256": claimed,
        "facet_complement_sha256": facet_digest,
        "packed_face_table_sha256": face_digest,
        "stellar_facet_counts": [20, 28, 36],
        "base_face_vector": list(base_vector),
        "minimal_nonface_size_histogram": {
            str(size): count for size, count in Counter(map(len, target)).items()
        },
        "blocker_order_histogram": {"5": 146},
        "ridge_degree_histogram": {"2": 657},
        "dual_graph_vertex_count": connected,
        "target_face_vector": list(audit.face_vector),
        "all_links_mod2_homology_spheres": audit.all_links_are_mod2_homology_spheres,
        "ordinary_rational_polytopal": True,
    }
    digest_payload = dict(summary)
    digest_payload.pop("certificate_sha256")
    summary_digest = sha256(
        json.dumps(digest_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    if EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN" and summary_digest != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("hostile record-4 summary digest changed")
    summary["summary_sha256"] = summary_digest
    return summary


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2, sort_keys=True))
