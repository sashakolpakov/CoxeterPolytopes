#!/usr/bin/env python3
"""Independent exact replay of the record-14 polytopal completion."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json


GRAPH6 = "M???????????o?w??"
EDGES = ((0, 12), (1, 12), (0, 13), (1, 13), (2, 13))
A4 = ((1, 3, 4, 5), (1, 9, 10, 11))
A5 = (
    (0, 2, 6, 7, 8),
    (2, 8, 9, 10, 11),
    (3, 4, 5, 6, 7),
    (8, 9, 10, 11, 12),
)
CLONES = ((3, 4, 5), (6, 7), (9, 10, 11))
BASE = (0, 1, 2, 3, 6, 8, 9, 12, 13)
GALE = (
    (1, 0, 0, 0, 31, -53, -52, -31, 47),
    (0, 1, 0, 0, 21, -39, -32, -64, 60),
    (0, 0, 1, 0, 23, -21, -22, -50, 45),
    (0, 0, 0, 1, 50, -20, -44, -10, 15),
)
POSITIVE_KERNEL = (
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
EXPECTED_SUMMARY_SHA256 = (
    "149d05f32c19e056e33d750c228784c28d2d61d06952c3a61018b509195e2ff4"
)


def bitset(values) -> int:
    return sum(1 << value for value in values)


def elements(value: int, order: int = 14) -> tuple[int, ...]:
    return tuple(index for index in range(order) if value & (1 << index))


def minimal_transversals(rows: tuple[int, ...], order: int) -> tuple[int, ...]:
    result = []
    for candidate in range(1, 1 << order):
        if not all(candidate & row for row in rows):
            continue
        if all(
            not all((candidate ^ (1 << vertex)) & row for row in rows)
            for vertex in elements(candidate, order)
        ):
            result.append(candidate)
    return tuple(result)


def determinant(columns: tuple[int, int, int, int]) -> int:
    total = 0
    for ordering in permutations(range(4)):
        inversions = sum(
            ordering[left] > ordering[right]
            for left in range(4)
            for right in range(left + 1, 4)
        )
        product = 1
        for row, selected_row in enumerate(ordering):
            product *= GALE[selected_row][columns[row]]
        total += (-1 if inversions % 2 else 1) * product
    return total


def audit() -> dict[str, object]:
    clutter = tuple(bitset(row) for row in (*EDGES, *A4, *A5))
    facets = minimal_transversals(clutter, 14)
    if Counter(value.bit_count() for value in facets) != {5: 106}:
        raise AssertionError("independent record-14 blocker changed")

    faces = tuple(
        candidate
        for candidate in range(1 << 14)
        if not any(candidate & row == row for row in clutter)
    )
    f_vector = tuple(
        sum(candidate.bit_count() == size for candidate in faces)
        for size in range(1, 10)
    )
    if f_vector != (14, 86, 310, 733, 1192, 1344, 1020, 477, 106):
        raise AssertionError("independent record-14 face vector changed")
    ridge_histogram = Counter(
        sum(facet & bitset(six) == facet for facet in facets)
        for six in combinations(range(14), 6)
    )
    ridge_histogram.pop(0, None)
    if ridge_histogram != {2: 477}:
        raise AssertionError("independent record-14 ridge degrees changed")

    representative = {
        vertex: block[0] for block in CLONES for vertex in block
    }
    base_rows = {
        frozenset(representative.get(vertex, vertex) for vertex in row)
        for row in (*EDGES, *A4, *A5)
    }
    base_faces = tuple(
        frozenset(chosen)
        for size in range(1, 10)
        for chosen in combinations(BASE, size)
        if not any(row <= frozenset(chosen) for row in base_rows)
    )
    base_facets = tuple(
        face for face in base_faces if not any(face < other for other in base_faces)
    )
    base_f_vector = tuple(
        sum(len(face) == size for face in base_faces) for size in range(1, 5)
    )
    if len(base_facets) != 19 or base_f_vector != (9, 28, 38, 19):
        raise AssertionError("independent record-14 base complex changed")

    minors = {
        basis: determinant(basis) for basis in combinations(range(9), 4)
    }
    if any(value == 0 for value in minors.values()):
        raise AssertionError("record-14 Gale matrix ceased to be uniform")
    positive_circuits = set()
    for support in combinations(range(9), 5):
        coefficients = tuple(
            (-1 if index % 2 else 1)
            * minors[support[:index] + support[index + 1 :]]
            for index in range(5)
        )
        if all(value > 0 for value in coefficients) or all(
            value < 0 for value in coefficients
        ):
            positive_circuits.add(frozenset(support))
    base_index = {vertex: index for index, vertex in enumerate(BASE)}
    expected_circuits = {
        frozenset(
            index
            for vertex, index in base_index.items()
            if vertex not in facet
        )
        for facet in base_facets
    }
    if positive_circuits != expected_circuits:
        raise AssertionError("record-14 Gale positive circuits changed")
    if not all(value > 0 for value in POSITIVE_KERNEL):
        raise AssertionError("record-14 Gale dependence is not positive")
    kernel_product = tuple(
        sum(GALE[row][column] * POSITIVE_KERNEL[column] for column in range(9))
        for row in range(4)
    )
    if kernel_product != (0, 0, 0, 0):
        raise AssertionError("record-14 positive Gale dependence changed")

    blocks = {block[0]: block for block in CLONES}
    expanded = {
        frozenset(
            expanded_vertex
            for vertex in row
            for expanded_vertex in blocks.get(vertex, (vertex,))
        )
        for row in base_rows
    }
    if expanded != {frozenset(row) for row in (*EDGES, *A4, *A5)}:
        raise AssertionError("independent record-14 multiwedge expansion changed")

    result: dict[str, object] = {
        "graph6": GRAPH6,
        "minimal_nonface_histogram": {"2": 5, "4": 2, "5": 4},
        "facet_count": len(facets),
        "facet_complement_sha256": sha256(
            ("\n".join(map(str, facets)) + "\n").encode()
        ).hexdigest(),
        "f_vector": list(f_vector),
        "ridge_degree_histogram": {str(key): value for key, value in ridge_histogram.items()},
        "base_f_vector": list(base_f_vector),
        "base_gale_nonzero_minors": len(minors),
        "base_positive_circuit_count": len(positive_circuits),
        "base_positive_circuits_match_facets": True,
        "positive_kernel_product": list(kernel_product),
        "clone_block_sizes": [len(block) for block in CLONES],
        "multiwedge_dimension": 4 + sum(len(block) - 1 for block in CLONES),
        "ordinary_polytopal": True,
    }
    digest = sha256(
        json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    result["summary_sha256"] = digest
    if EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("record-14 independent summary digest changed")
    return result


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
