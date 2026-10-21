#!/usr/bin/env python3
"""Exact replay of the I_(1,9) Hamming P-first effectivity obstruction."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import subprocess
from collections import Counter
from fractions import Fraction
from pathlib import Path
from typing import Iterable


Vector = tuple[int, ...]
Matrix = list[list[int]]
DIMENSION = 10
METRIC = (1,) + (-1,) * 9

HAMMING_GENERATORS: tuple[tuple[int, ...], ...] = (
    (0, 0, 0, 1, 0, 1, 1, 1),
    (0, 0, 1, 0, 1, 0, 1, 1),
    (0, 1, 0, 0, 1, 1, 0, 1),
    (1, 0, 0, 0, 1, 1, 1, 0),
)


def dot(x: Vector, y: Vector) -> int:
    return sum(sign * a * b for sign, a, b in zip(METRIC, x, y))


def add(x: Vector, y: Vector) -> Vector:
    return tuple(a + b for a, b in zip(x, y))


def subtract(x: Vector, y: Vector) -> Vector:
    return tuple(a - b for a, b in zip(x, y))


def scale(coefficient: int, x: Vector) -> Vector:
    return tuple(coefficient * entry for entry in x)


def bareiss_determinant(matrix: list[list[int]]) -> int:
    work = [row[:] for row in matrix]
    n = len(work)
    assert n and all(len(row) == n for row in work)
    sign = 1
    previous = 1
    for pivot_index in range(n - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (
                    row
                    for row in range(pivot_index + 1, n)
                    if work[row][pivot_index]
                ),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, n):
            for column in range(pivot_index + 1, n):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                assert numerator % previous == 0
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def rational_rank(rows: Iterable[Vector]) -> int:
    matrix = [list(map(Fraction, row)) for row in rows]
    if not matrix:
        return 0
    rank = 0
    columns = len(matrix[0])
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        value = matrix[rank][column]
        matrix[rank] = [entry / value for entry in matrix[rank]]
        for row in range(len(matrix)):
            if row != rank and matrix[row][column]:
                coefficient = matrix[row][column]
                matrix[row] = [
                    entry - coefficient * pivot_entry
                    for entry, pivot_entry in zip(matrix[row], matrix[rank])
                ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return [
        [
            sum(left[row][middle] * right[middle][column] for middle in range(DIMENSION))
            for column in range(DIMENSION)
        ]
        for row in range(DIMENSION)
    ]


def transpose(matrix: Matrix) -> Matrix:
    return [list(row) for row in zip(*matrix)]


def reflection_matrix(root: Vector) -> Matrix:
    square = dot(root, root)
    assert square < 0
    columns: list[Vector] = []
    for column in range(DIMENSION):
        basis = tuple(1 if row == column else 0 for row in range(DIMENSION))
        coefficient = Fraction(-2 * dot(basis, root), square)
        image = tuple(Fraction(basis[row]) + coefficient * root[row] for row in range(DIMENSION))
        assert all(entry.denominator == 1 for entry in image)
        columns.append(tuple(int(entry) for entry in image))
    return [list(row) for row in zip(*columns)]


def verify_reflection(root: Vector) -> Matrix:
    matrix = reflection_matrix(root)
    identity = [
        [1 if row == column else 0 for column in range(DIMENSION)]
        for row in range(DIMENSION)
    ]
    metric_matrix = [
        [METRIC[row] if row == column else 0 for column in range(DIMENSION)]
        for row in range(DIMENSION)
    ]
    assert matmul(matrix, matrix) == identity
    assert matmul(matmul(transpose(matrix), metric_matrix), matrix) == metric_matrix
    return matrix


def hamming_code() -> tuple[tuple[int, ...], ...]:
    words: set[tuple[int, ...]] = set()
    for mask in range(16):
        word = [0] * 8
        for generator_index, generator in enumerate(HAMMING_GENERATORS):
            if (mask >> generator_index) & 1:
                word = [left ^ right for left, right in zip(word, generator)]
        words.add(tuple(word))
    result = tuple(sorted(words))
    assert len(result) == 16
    assert Counter(sum(word) for word in result) == {0: 1, 4: 14, 8: 1}
    assert min(
        sum(left ^ right for left, right in zip(first, second))
        for index, first in enumerate(result)
        for second in result[index + 1 :]
    ) == 4
    return result


def packet() -> dict[str, object]:
    canonical: Vector = (-3,) + (1,) * 9
    anticanonical: Vector = (3,) + (-1,) * 9
    exceptional_9: Vector = (0,) * 9 + (1,)

    roots: list[Vector] = []
    for pair in range(4):
        root = [0] * DIMENSION
        root[1 + 2 * pair] = 1
        root[2 + 2 * pair] = -1
        roots.append(tuple(root))
    for omitted_pair in (3, 2, 1, 0):
        root = [0] * DIMENSION
        root[0] = 2
        for pair in range(4):
            if pair != omitted_pair:
                root[1 + 2 * pair] = -1
                root[2 + 2 * pair] = -1
        roots.append(tuple(root))
    assert all(dot(root, root) == -2 for root in roots)
    assert all(dot(anticanonical, root) == 0 for root in roots)
    assert [
        [dot(first, second) for second in roots] for first in roots
    ] == [
        [-2 if row == column else 0 for column in range(8)] for row in range(8)
    ]

    opposite_roots = [subtract(anticanonical, root) for root in roots]
    code = hamming_code()
    exceptional: dict[str, Vector] = {}
    for word in code:
        weight = sum(word)
        entries: list[int] = []
        for coordinate in range(DIMENSION):
            value = (
                Fraction(exceptional_9[coordinate])
                - Fraction(
                    sum(word[index] * roots[index][coordinate] for index in range(8)),
                    2,
                )
                + Fraction(weight * anticanonical[coordinate], 4)
            )
            assert value.denominator == 1
            entries.append(int(value))
        exceptional["".join(map(str, word))] = tuple(entries)

    wall_records: list[dict[str, object]] = []
    for index, root in enumerate(roots, 1):
        wall_records.append(
            {"name": f"A{index}", "kind": "minus_two", "vector": root}
        )
    for index, root in enumerate(opposite_roots, 1):
        wall_records.append(
            {"name": f"Beta{index}", "kind": "minus_two", "vector": root}
        )
    for label, root in exceptional.items():
        wall_records.append(
            {"name": f"E_{label}", "kind": "minus_one", "vector": root}
        )
    walls = [record["vector"] for record in wall_records]
    assert all(isinstance(wall, tuple) for wall in walls)

    strict_interior: Vector = (5, -2, -1, -2, -1, -2, -1, -2, -1, -1)
    assert dot(strict_interior, strict_interior) == 4
    assert [dot(strict_interior, wall) for wall in walls] == [1] * 32

    assert all(entry % 2 == 1 for entry in canonical)
    assert dot(canonical, canonical) == 0 == 10 - DIMENSION
    assert anticanonical == scale(-1, canonical)
    assert math.gcd(*map(abs, anticanonical)) == 1
    assert dot(anticanonical, anticanonical) == 0
    assert [dot(anticanonical, wall) for wall in walls] == [0] * 16 + [1] * 16
    assert [dot(canonical, wall) for wall in walls[16:]] == [-1] * 16

    for index in range(8):
        assert add(roots[index], opposite_roots[index]) == anticanonical
    assert all(
        dot(roots[first], roots[second]) == 0
        and dot(roots[first], opposite_roots[second]) == 0
        and dot(opposite_roots[first], opposite_roots[second]) == 0
        for first in range(8)
        for second in range(8)
        if first != second
    )

    for label, curve in exceptional.items():
        word = tuple(map(int, label))
        assert dot(curve, curve) == -1
        assert dot(anticanonical, curve) == 1
        assert [dot(curve, root) for root in roots] == list(word)
        assert [dot(curve, root) for root in opposite_roots] == [
            1 - entry for entry in word
        ]

    coxeter_histogram: Counter[str] = Counter()
    for first_index, first in enumerate(walls):
        first_norm = -dot(first, first)
        for second in walls[first_index + 1 :]:
            second_norm = -dot(second, second)
            pairing = dot(first, second)
            assert pairing >= 0
            if pairing == 0:
                order = "2"
            elif 2 * pairing * pairing == first_norm * second_norm:
                order = "4"
            elif pairing * pairing == first_norm * second_norm:
                order = "infinity"
            else:
                raise AssertionError(
                    f"non-Coxeter pair ({first_norm}, {second_norm}, {pairing})"
                )
            coxeter_histogram[order] += 1
    assert coxeter_histogram == {"2": 352, "4": 128, "infinity": 16}

    reflection_matrices = {
        record["name"]: verify_reflection(record["vector"]) for record in wall_records
    }
    reflection_bytes = json.dumps(
        reflection_matrices, sort_keys=True, separators=(",", ":")
    ).encode()

    spanning_names = (
        "A1",
        "A2",
        "A3",
        "A4",
        "E_00000000",
        "E_00010111",
        "E_00101011",
        "E_01001101",
        "E_10001110",
        "E_00111100",
    )
    by_name = {record["name"]: record["vector"] for record in wall_records}
    spanning_determinant = bareiss_determinant(
        [list(by_name[name]) for name in spanning_names]
    )
    assert spanning_determinant == -1

    distinguished_word = "00111100"
    missing_root = subtract(exceptional[distinguished_word], exceptional["00000000"])
    assert missing_root == (1, 0, 0, 0, 0, -1, 0, -1, 0, -1)
    assert dot(missing_root, missing_root) == -2
    assert dot(anticanonical, missing_root) == 0
    assert dot(strict_interior, missing_root) == 0
    missing_reflection = verify_reflection(missing_root)

    return {
        "canonical": canonical,
        "anticanonical": anticanonical,
        "orthogonal_roots": roots,
        "opposite_roots": opposite_roots,
        "hamming_code": code,
        "exceptional_walls": exceptional,
        "walls": wall_records,
        "strict_interior": strict_interior,
        "coxeter_order_histogram": dict(sorted(coxeter_histogram.items())),
        "wall_reflection_matrices_sha256": hashlib.sha256(reflection_bytes).hexdigest(),
        "spanning_names": spanning_names,
        "spanning_determinant": spanning_determinant,
        "missing_root": missing_root,
        "missing_root_reflection_matrix": missing_reflection,
    }


def normalize_ray(entries: list[Fraction]) -> Vector:
    denominator_lcm = 1
    for entry in entries:
        denominator_lcm = math.lcm(denominator_lcm, entry.denominator)
    integral = [int(entry * denominator_lcm) for entry in entries]
    divisor = math.gcd(*map(abs, integral))
    assert divisor
    primitive = tuple(entry // divisor for entry in integral)
    if primitive[0] < 0:
        primitive = scale(-1, primitive)
    assert primitive[0] > 0
    return primitive


def exact_polyhedron_statistics(walls: list[Vector]) -> dict[str, object]:
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        fallback = Path("/opt/homebrew/bin/cddexec_gmp")
        if fallback.exists():
            executable = str(fallback)
    if executable is None:
        raise RuntimeError("exact finite-volume replay requires cddexec_gmp")

    lines = ["H-representation", "begin", f"{len(walls)} {DIMENSION + 1} rational"]
    for wall in walls:
        coefficients = [wall[0], *(-entry for entry in wall[1:])]
        lines.append("0 " + " ".join(map(str, coefficients)))
    lines.append("end")
    completed = subprocess.run(
        [executable, "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    output = completed.stdout.splitlines()
    representation = output.index("V-representation")
    begin = output.index("begin", representation)
    count = int(output[begin + 1].split()[0])
    rays: list[Vector] = []
    for line in output[begin + 2 : begin + 2 + count]:
        fields = line.split()
        assert fields[0] == "0"
        rays.append(normalize_ray([Fraction(field) for field in fields[1:]]))
    assert len(rays) == count
    assert len(set(rays)) == count
    rays = sorted(rays)

    square_histogram = Counter(dot(ray, ray) for ray in rays)
    assert square_histogram == {0: 226, 1: 1040, 2: 1040}
    assert all(ray[0] > 0 and dot(ray, ray) >= 0 for ray in rays)
    active_histogram: Counter[str] = Counter()
    active_counts = [0] * len(walls)
    active_ranks: list[int] = []
    active_ray_lists: list[list[Vector]] = [[] for _ in walls]
    for ray in rays:
        active = [index for index, wall in enumerate(walls) if dot(ray, wall) == 0]
        kind = "ideal" if dot(ray, ray) == 0 else "finite"
        active_histogram[f"{kind}_{len(active)}"] += 1
        for index in active:
            active_counts[index] += 1
            active_ray_lists[index].append(ray)
    for active_rays in active_ray_lists:
        active_ranks.append(rational_rank(active_rays))
    assert active_histogram == {"finite_9": 2080, "ideal_12": 224, "ideal_16": 2}
    assert active_counts == [670] * 32
    assert active_ranks == [9] * 32

    canonical_rays = "".join(
        " ".join(map(str, ray)) + "\n" for ray in rays
    ).encode()
    ray_hash = hashlib.sha256(canonical_rays).hexdigest()
    assert ray_hash == "33ad7a1bb8df2dda4c4fd82b7cea62577f9f64d1f9972222707923104e248b09"
    assert max(abs(entry) for ray in rays for entry in ray) == 15
    return {
        "extreme_ray_count": count,
        "primitive_ray_square_histogram": {
            str(square): multiplicity
            for square, multiplicity in sorted(square_histogram.items())
        },
        "active_histogram": dict(sorted(active_histogram.items())),
        "facet_active_ray_counts": active_counts,
        "facet_active_span_ranks": active_ranks,
        "maximum_absolute_primitive_coordinate": 15,
        "sorted_primitive_ray_sha256": ray_hash,
        "finite_volume_reason": (
            "all exact extreme rays are future timelike or future isotropic; "
            "there are no hyperideal directions"
        ),
    }


def payload() -> dict[str, object]:
    built = packet()
    walls = [record["vector"] for record in built["walls"]]
    polyhedron = exact_polyhedron_statistics(walls)
    raw: dict[str, object] = {
        "schema": "rational_pfirst_k_escape_hamming_effectivity_no_go/v1",
        "lattice": "I_(1,9)",
        "metric_diagonal": list(METRIC),
        "rank": DIMENSION,
        "canonical": list(built["canonical"]),
        "canonical_square": dot(built["canonical"], built["canonical"]),
        "anticanonical_cusp": list(built["anticanonical"]),
        "hamming_generators": [list(row) for row in HAMMING_GENERATORS],
        "hamming_codewords": ["".join(map(str, word)) for word in built["hamming_code"]],
        "orthogonal_minus_two_roots": [
            list(root) for root in built["orthogonal_roots"]
        ],
        "opposite_minus_two_roots": [
            list(root) for root in built["opposite_roots"]
        ],
        "minus_one_walls": {
            label: list(root) for label, root in built["exceptional_walls"].items()
        },
        "strict_interior": list(built["strict_interior"]),
        "strict_interior_square": dot(built["strict_interior"], built["strict_interior"]),
        "strict_wall_slacks": [1] * 32,
        "wall_count": 32,
        "wall_kind_histogram": {"minus_one": 16, "minus_two": 16},
        "coxeter_order_histogram": built["coxeter_order_histogram"],
        "wall_reflection_matrices_sha256": built[
            "wall_reflection_matrices_sha256"
        ],
        "unimodular_wall_subset": list(built["spanning_names"]),
        "unimodular_wall_subset_determinant": built["spanning_determinant"],
        "polyhedron": polyhedron,
        "canonical_escape": {
            "minus_K_is_ideal_cusp": True,
            "minus_K_pairings_with_minus_two_walls": [0] * 16,
            "minus_K_pairings_with_minus_one_walls": [1] * 16,
            "K_pairing_with_every_minus_one_wall": -1,
            "K_is_not_in_P": True,
        },
        "simultaneous_effectivity_obstruction": {
            "anticanonical_decompositions": 8,
            "decomposition_formula": "A_i + Beta_i = -K",
            "cross_pairings_between_distinct_decompositions": 0,
            "two_distinct_decompositions_give_basepoint_free_pencil": True,
            "minimum_euler_contribution_per_reducible_fiber": 2,
            "minimum_euler_contribution_of_eight_fibers": 16,
            "rational_surface_c2": 12,
            "conclusion_over_C": (
                "the 16 minus-two wall classes cannot all be represented "
                "simultaneously by distinct irreducible curves"
            ),
        },
        "explicit_missing_reflection": {
            "distinguished_codeword": "00111100",
            "formula": "r = E_00111100 - E_00000000",
            "vector": list(built["missing_root"]),
            "square": dot(built["missing_root"], built["missing_root"]),
            "anticanonical_pairing": dot(
                built["anticanonical"], built["missing_root"]
            ),
            "strict_interior_pairing": dot(
                built["strict_interior"], built["missing_root"]
            ),
            "reflection_matrix": built["missing_root_reflection_matrix"],
            "proper_overgroup": True,
            "reason": (
                "the integral missing mirror contains a strict interior point "
                "of the fundamental chamber"
            ),
        },
        "scope": {
            "K_escape_plus_simultaneous_effectivity_counterexample": "PROVED",
            "relative_reflection_maximality_of_this_packet": "DISPROVED",
            "maximal_plus_K_escape_nonrealizability_counterexample": "OPEN",
            "universal_polytope_to_surface_existence": "OPEN",
            "necessary_maximality_consequence": (
                "every host-integral reflecting mirror is disjoint from the "
                "interior of a reflection-maximal chamber; equivalently each "
                "cusp link is host-reflection-saturated"
            ),
        },
    }
    canonical = json.dumps(raw, sort_keys=True, separators=(",", ":")).encode()
    raw["payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    return raw


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    result = payload()
    default = Path(
        "data/enumerations/rational_pfirst_k_escape_hamming_effectivity_no_go.json"
    )
    target = args.output or default
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        assert target.read_text() == rendered
    elif args.output:
        target.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
