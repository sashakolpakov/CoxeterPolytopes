"""Independent replay helpers for the record-34 ``a=0`` census.

The exhaustive generator is ``scripts/certify_h9_14_record34_a0.cpp``.
This module deliberately does not share its implementation: it reconstructs
the blocker, face vector, minimal transversals, ridge degrees, and symmetry
orbits directly from the archived six missing five-faces.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


VERTEX_COUNT = 14
A_VERTICES = (0, 11, 13)
B_VERTICES = (1, 12)
ACTIVE_VERTICES = frozenset(A_VERTICES + B_VERTICES)
ISOLATED_VERTICES = tuple(range(2, 11))
MISSING_EDGES = ((0, 11), (0, 13), (11, 13), (1, 12))
TARGET_F_VECTOR = (14, 87, 318, 762, 1254, 1428, 1092, 513, 114)


def mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


EDGE_MASKS = tuple(mask(edge) for edge in MISSING_EDGES)
ALL_VERTICES_MASK = (1 << VERTEX_COUNT) - 1


def vertices(value: int) -> tuple[int, ...]:
    return tuple(
        vertex for vertex in range(VERTEX_COUNT) if value & (1 << vertex)
    )


def is_transversal(value: int, missing_fives: Sequence[int]) -> bool:
    return all(value & edge for edge in EDGE_MASKS) and all(
        value & face for face in missing_fives
    )


def minimal_transversals(missing_fives: Sequence[int]) -> tuple[int, ...]:
    values = []
    for value in range(1 << VERTEX_COUNT):
        if not is_transversal(value, missing_fives):
            continue
        if all(
            not is_transversal(value ^ (1 << vertex), missing_fives)
            for vertex in vertices(value)
        ):
            values.append(value)
    return tuple(values)


def face_vector(missing_fives: Sequence[int]) -> tuple[int, ...]:
    counts = Counter()
    for value in range(1 << VERTEX_COUNT):
        if any(value & edge == edge for edge in EDGE_MASKS):
            continue
        if any(value & face == face for face in missing_fives):
            continue
        size = value.bit_count()
        if size:
            counts[size] += 1
    return tuple(counts[size] for size in range(1, 10))


def ridge_audit(
    missing_fives: Sequence[int],
) -> tuple[dict[int, int], int]:
    """Return positive ridge-degree histogram and degree-two dual components."""

    transversals = tuple(
        value
        for value in range(1 << VERTEX_COUNT)
        if value.bit_count() == 5 and is_transversal(value, missing_fives)
    )
    indices = {value: index for index, value in enumerate(transversals)}
    parent = list(range(len(transversals)))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(first: int, second: int) -> None:
        first = find(first)
        second = find(second)
        if first != second:
            parent[second] = first

    histogram: Counter[int] = Counter()
    for chosen in combinations(range(VERTEX_COUNT), 6):
        six_set = mask(chosen)
        incident = []
        for vertex in chosen:
            candidate = six_set ^ (1 << vertex)
            if candidate in indices:
                incident.append(indices[candidate])
        if incident:
            histogram[len(incident)] += 1
        if len(incident) == 2:
            union(*incident)

    components = len({find(index) for index in range(len(transversals))})
    return dict(sorted(histogram.items())), components


def canonical_key(missing_fives: Sequence[int]) -> tuple[int, ...]:
    """Canonicalize under ``S3(A) x S2(B) x S9(I) x S6(rows)``."""

    rows = tuple(frozenset(vertices(face)) for face in missing_fives)
    best: tuple[int, ...] | None = None
    for a_image in permutations(A_VERTICES):
        for b_image in permutations(B_VERTICES):
            active_map = dict(zip(A_VERTICES, a_image, strict=True))
            active_map.update(zip(B_VERTICES, b_image, strict=True))
            labels = []
            isolate_parts = []
            for row in rows:
                active = tuple(row & ACTIVE_VERTICES)
                if len(active) > 1:
                    raise ValueError("an A5 row has more than one active vertex")
                labels.append(-1 if not active else active_map[active[0]])
                isolate_parts.append(row - ACTIVE_VERTICES)
            for row_order in permutations(range(6)):
                ordered_labels = tuple(labels[row] for row in row_order)
                patterns = tuple(
                    sorted(
                        sum(
                            int(vertex in isolate_parts[row]) << position
                            for position, row in enumerate(row_order)
                        )
                        for vertex in ISOLATED_VERTICES
                    )
                )
                candidate = ordered_labels + patterns
                if best is None or candidate < best:
                    best = candidate
    if best is None:
        raise AssertionError("canonicalization unexpectedly found no image")
    return best


def certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return sha256(encoded).hexdigest()


def verify_manifest(path: str | Path) -> dict[str, object]:
    path = Path(path)
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema_version") != 1:
        raise ValueError("unexpected record-34 a=0 schema")
    if tuple(map(tuple, document["missing_edges"])) != MISSING_EDGES:
        raise ValueError("record-34 missing graph changed")
    if tuple(document["target_f_vector"]) != TARGET_F_VECTOR:
        raise ValueError("record-34 target f-vector changed")

    strata = document["strata"]
    aggregate_keys = (
        "genbg_representatives",
        "active_colourings",
        "distinct_row_states",
        "no_small_transversal_states",
        "target_114_states",
        "pure_target_states",
        "ridge_valid_states",
    )
    observed_totals = {
        key: sum(int(record[key]) for record in strata) for key in aggregate_keys
    }
    if observed_totals != document["totals"]:
        raise ValueError("stratum totals do not add to the archived aggregate")

    records = document["target_states"]
    if len(records) != 18:
        raise ValueError("the target-state list is not the complete 18-state list")
    keys = []
    for record in records:
        faces = tuple(mask(face) for face in record["missing_fives"])
        if len(faces) != 6 or len(set(faces)) != 6:
            raise ValueError("target state does not have six distinct A5 rows")
        if any(face.bit_count() != 5 for face in faces):
            raise ValueError("target state contains a wrongly sized A5")
        if any(any(face & edge == edge for edge in EDGE_MASKS) for face in faces):
            raise ValueError("target A5 contains a missing graph edge")
        if any((face & mask(ACTIVE_VERTICES)).bit_count() > 1 for face in faces):
            raise ValueError("target A5 violates the one-active-vertex gate")
        incidences = Counter(
            vertex for face in faces for vertex in vertices(face)
        )
        if any(incidences[vertex] > 1 for vertex in A_VERTICES):
            raise ValueError("an A vertex occurs in more than one A5")
        if any(incidences[vertex] > 2 for vertex in B_VERTICES):
            raise ValueError("a B vertex occurs in more than two A5s")
        if any(incidences[vertex] > 5 for vertex in ISOLATED_VERTICES):
            raise ValueError("an isolate occurs in all six A5s")

        blocker = minimal_transversals(faces)
        if len(blocker) != 114 or {value.bit_count() for value in blocker} != {5}:
            raise ValueError("target blocker is not pure with 114 order-five members")
        if face_vector(faces) != TARGET_F_VECTOR:
            raise ValueError("target face vector does not replay")
        histogram, components = ridge_audit(faces)
        expected_histogram = {
            int(degree): int(count)
            for degree, count in record["ridge_degree_histogram"].items()
        }
        if histogram != expected_histogram:
            raise ValueError("target ridge histogram does not replay")
        if components != record["degree_two_dual_components"]:
            raise ValueError("target degree-two dual components do not replay")
        if set(histogram) <= {2}:
            raise ValueError("a purported target state unexpectedly passes ridges")
        keys.append(canonical_key(faces))

    key_counts = Counter(keys)
    expected_multiplicities = sorted(document["canonical_orbit_multiplicities"])
    if sorted(key_counts.values()) != expected_multiplicities:
        raise ValueError("the target-state symmetry quotient changed")
    if len(key_counts) != document["canonical_orbit_count"]:
        raise ValueError("the target-state orbit count changed")

    source = path.parents[2] / document["enumerator_source"]
    source_digest = sha256(source.read_bytes()).hexdigest()
    if source_digest != document["enumerator_source_sha256"]:
        raise ValueError("the C++ enumerator source digest changed")
    if certificate_digest(document) != document["certificate_sha256"]:
        raise ValueError("the record-34 a=0 certificate digest changed")
    return document

