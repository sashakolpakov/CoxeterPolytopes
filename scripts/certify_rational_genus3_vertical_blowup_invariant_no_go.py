#!/usr/bin/env python3
"""Exact certificate for the genus-three vertical-blowup invariant no-go."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from typing import Sequence


Vector = tuple[int, ...]


def pair(left: Sequence[int], right: Sequence[int]) -> int:
    """Negative exceptional form in the e_i coordinates."""

    return -sum(a * b for a, b in zip(left, right))


def fibration_row(core_rho: int, blowups: int) -> dict[str, int]:
    if core_rho < 1:
        raise ValueError("a rational surface has positive Picard rank")
    if blowups < 0:
        raise ValueError("the blowup count must be nonnegative")
    rho = core_rho + blowups
    k_f_square = 26 - rho
    e_f = rho + 10
    return {
        "core_rho": core_rho,
        "vertical_blowups": blowups,
        "rho": rho,
        "surface_K_square": 10 - rho,
        "surface_euler": rho + 2,
        "relative_K_square": k_f_square,
        "relative_euler": e_f,
        "relative_noether_sum": k_f_square + e_f,
        "hodge_degree": 3,
    }


def tail_vectors(length: int) -> tuple[Vector, ...]:
    if length < 1:
        raise ValueError("the tail length must be positive")
    roots: list[Vector] = []
    for index in range(length - 1):
        root = [0] * length
        root[index] = 1
        root[index + 1] = -1
        roots.append(tuple(root))
    endpoint = [0] * length
    endpoint[-1] = 1
    roots.append(tuple(endpoint))
    return tuple(roots)


def reflection(root: Vector, vector: Vector) -> Vector:
    square = pair(root, root)
    coefficient = Fraction(2 * pair(vector, root), square)
    result = tuple(
        Fraction(entry) - coefficient * root_entry
        for entry, root_entry in zip(vector, root)
    )
    if any(entry.denominator != 1 for entry in result):
        raise AssertionError("tail reflection ceased to be integral")
    return tuple(int(entry) for entry in result)


def orbit(seed: Vector, roots: Sequence[Vector]) -> set[Vector]:
    found = {seed}
    queue = [seed]
    while queue:
        vector = queue.pop()
        for root in roots:
            image = reflection(root, vector)
            if image not in found:
                found.add(image)
                queue.append(image)
    return found


def rank(vectors: Sequence[Vector]) -> int:
    if not vectors:
        return 0
    matrix = [[Fraction(entry) for entry in vector] for vector in vectors]
    row = 0
    columns = len(matrix[0])
    for column in range(columns):
        pivot = next(
            (index for index in range(row, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        value = matrix[row][column]
        matrix[row] = [entry / value for entry in matrix[row]]
        for index in range(len(matrix)):
            if index == row or not matrix[index][column]:
                continue
            value = matrix[index][column]
            matrix[index] = [
                entry - value * pivot_entry
                for entry, pivot_entry in zip(matrix[index], matrix[row])
            ]
        row += 1
        if row == len(matrix):
            break
    return row


def tail_row(length: int) -> dict[str, object]:
    roots = tail_vectors(length)
    gram = [[pair(left, right) for right in roots] for left in roots]
    seed = (-1, *((0,) * (length - 1)))
    seed_orbit = orbit(seed, roots)
    expected_orbit = {
        tuple(sign if coordinate == index else 0 for coordinate in range(length))
        for index in range(length)
        for sign in (-1, 1)
    }
    if seed_orbit != expected_orbit:
        raise AssertionError((length, seed_orbit, expected_orbit))
    return {
        "length": length,
        "surface_root_squares": [pair(root, root) for root in roots],
        "surface_gram": gram,
        "host_root_squares": [2 * pair(root, root) for root in roots],
        "canonical_host_colors": [0] * (length - 1) + [-2],
        "mixed_projection": list(seed),
        "mixed_orbit_size": len(seed_orbit),
        "mixed_orbit_rank": rank(sorted(seed_orbit)),
        "signed_coordinate_orbit": [list(vector) for vector in sorted(seed_orbit)],
    }


def build_certificate() -> dict[str, object]:
    core_ranks = (6, 11, 18, 26)
    blowup_counts = (0, 1, 2, 5, 17, 64)
    fibration = [
        fibration_row(core_rho, blowups)
        for core_rho in core_ranks
        for blowups in blowup_counts
    ]
    if any(row["relative_noether_sum"] != 36 for row in fibration):
        raise AssertionError("relative Noether identity changed")
    tails = [tail_row(length) for length in (1, 2, 3, 5, 11, 37)]
    return {
        "theorem": "genus-three vertical blowup invariant no-go",
        "genus": 3,
        "base_genus": 0,
        "hodge_degree": 3,
        "relative_noether_constant": 36,
        "fibration_samples": fibration,
        "tail_samples": tails,
        "stable_data_bounds_tail": False,
        "relative_minimal_data_bounds_tail": False,
        "raw_relative_K_nef_after_vertical_blowup": False,
        "finite_volume_orbit_span_bounds_tail": False,
        "complete_nef_counterfamily_claimed": False,
        "reflection_maximal_counterfamily_claimed": False,
        "remaining_inputs": [
            "absolute total marked-facet degree",
            "finite equivariant contraction orbit or descent",
        ],
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pretty", action="store_true")
    arguments = parser.parse_args()
    payload = build_certificate()
    envelope = {"payload": payload, "payload_sha256": digest(payload)}
    print(json.dumps(envelope, indent=2 if arguments.pretty else None, sort_keys=True))


if __name__ == "__main__":
    main()
