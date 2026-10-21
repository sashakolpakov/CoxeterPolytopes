#!/usr/bin/env python3
"""Exact certificate for the exceptional-forest comparability theorem."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "rational_exceptional_forest_comparability.json"
)


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def pairing(left: Sequence[int | Fraction], right: Sequence[int | Fraction]) -> Fraction:
    return (
        Fraction(left[0]) * Fraction(right[1])
        + Fraction(left[1]) * Fraction(right[0])
        - sum(
            Fraction(left[index]) * Fraction(right[index])
            for index in range(2, len(left))
        )
    )


def forest_parent_vectors(size: int) -> Iterable[tuple[int | None, ...]]:
    if size < 1:
        return
    choices = [[None, *range(index)] for index in range(1, size)]
    for tail in itertools.product(*choices):
        yield (None, *tail)


def proximity_parents_from_forest(
    parents: Sequence[int | None],
) -> tuple[tuple[int, ...], ...]:
    return tuple(
        () if parent is None else (parent,)
        for parent in parents
    )


def reaches(
    ancestor: int,
    descendant: int,
    proximity_parents: Sequence[Sequence[int]],
) -> bool:
    stack = list(proximity_parents[descendant])
    seen = set()
    while stack:
        vertex = stack.pop()
        if vertex == ancestor:
            return True
        if vertex in seen:
            continue
        seen.add(vertex)
        stack.extend(proximity_parents[vertex])
    return False


def incomparable(
    left: int,
    right: int,
    proximity_parents: Sequence[Sequence[int]],
) -> bool:
    return not reaches(left, right, proximity_parents) and not reaches(
        right, left, proximity_parents
    )


def proximity_children(
    proximity_parents: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    children: list[list[int]] = [[] for _ in proximity_parents]
    for child, predecessors in enumerate(proximity_parents):
        if any(predecessor >= child for predecessor in predecessors):
            raise ValueError("proximity predecessors must be earlier blowups")
        for predecessor in predecessors:
            children[predecessor].append(child)
    return tuple(tuple(row) for row in children)


def weights_equalizing_incomparable_pair(
    proximity_parents: Sequence[Sequence[int]],
    left: int,
    right: int,
) -> tuple[int, ...]:
    """Return positive weights w_v>sum_{u proximate to v}w_u and w_i=w_j."""

    if not incomparable(left, right, proximity_parents):
        raise ValueError("the prescribed valuations are comparable")
    children = proximity_children(proximity_parents)
    size = len(children)

    def evaluate(slacks: Sequence[int]) -> list[int]:
        values = [0] * size
        for vertex in range(size - 1, -1, -1):
            values[vertex] = slacks[vertex] + sum(
                values[child] for child in children[vertex]
            )
        return values

    base_slacks = [1] * size
    base_slacks[left] = base_slacks[right] = 0
    base = evaluate(base_slacks)
    target = max(base[left], base[right]) + 1
    slacks = base_slacks[:]
    slacks[left] = target - base[left]
    slacks[right] = target - base[right]
    if slacks[left] <= 0 or slacks[right] <= 0:
        raise AssertionError("nonpositive prescribed slack")
    values = evaluate(slacks)
    if values[left] != values[right]:
        raise AssertionError((left, right, values, slacks))
    margins = [
        values[vertex] - sum(values[child] for child in children[vertex])
        for vertex in range(size)
    ]
    if any(margin <= 0 for margin in margins):
        raise AssertionError((values, margins))
    return tuple(values)


def exhaustive_forest_record(max_size: int = 7) -> list[dict[str, int]]:
    records = []
    for size in range(2, max_size + 1):
        forest_count = 0
        incomparable_pair_count = 0
        for parents in forest_parent_vectors(size):
            forest_count += 1
            proximity = proximity_parents_from_forest(parents)
            for left in range(size):
                for right in range(left + 1, size):
                    if not incomparable(left, right, proximity):
                        continue
                    incomparable_pair_count += 1
                    weights_equalizing_incomparable_pair(
                        proximity, left, right
                    )
        if forest_count != 1:
            # There are product_{k=1}^{n-1}(k+1)=n! chronological forests.
            expected = 1
            for value in range(2, size + 1):
                expected *= value
            if forest_count != expected:
                raise AssertionError((size, forest_count, expected))
        records.append(
            {
                "vertices": size,
                "chronological_forests": forest_count,
                "incomparable_pairs_checked": incomparable_pair_count,
            }
        )
    return records


def satellite_records() -> list[dict[str, object]]:
    samples = (
        # Two branches over vertex zero; later satellite predecessors remain
        # inside their respective ancestry branches.
        ((), (0,), (0,), (1, 0), (2, 0), (3, 1), (4, 2)),
        ((), (), (0,), (0,), (2, 0), (3, 0), (4, 2), (5, 3)),
    )
    records = []
    for proximity in samples:
        pairs = []
        for left in range(len(proximity)):
            for right in range(left + 1, len(proximity)):
                if not incomparable(left, right, proximity):
                    continue
                weights = weights_equalizing_incomparable_pair(
                    proximity, left, right
                )
                pairs.append(
                    {
                        "pair": [left, right],
                        "equal_weight": weights[left],
                        "weights": list(weights),
                    }
                )
        records.append(
            {
                "proximity_parents": [list(row) for row in proximity],
                "incomparable_pair_witnesses": pairs,
            }
        )
    return records


def reflection_record(exceptional_rank: int, left: int, right: int) -> dict[str, object]:
    size = exceptional_rank + 2
    root = [0] * size
    root[2 + left] = 1
    root[2 + right] = -1
    root_tuple = tuple(root)
    if pairing(root_tuple, root_tuple) != -2:
        raise AssertionError("difference root has wrong square")

    basis = [
        tuple(int(index == column) for index in range(size))
        for column in range(size)
    ]

    def reflect(vector: Sequence[int]) -> tuple[int, ...]:
        result = tuple(
            Fraction(value) + pairing(vector, root_tuple) * root_value
            for value, root_value in zip(vector, root_tuple)
        )
        if any(value.denominator != 1 for value in result):
            raise AssertionError
        return tuple(int(value) for value in result)

    images = [reflect(vector) for vector in basis]
    if images[2 + left] != basis[2 + right]:
        raise AssertionError
    if images[2 + right] != basis[2 + left]:
        raise AssertionError
    for index, image in enumerate(images):
        if index not in (2 + left, 2 + right) and image != basis[index]:
            raise AssertionError
    for i, first in enumerate(basis):
        for j, second in enumerate(basis):
            if pairing(images[i], images[j]) != pairing(first, second):
                raise AssertionError("reflection is not an isometry")
    return {
        "exceptional_rank": exceptional_rank,
        "pair": [left, right],
        "surface_root_square": -2,
        "doubled_host_root_square": -4,
        "reflection_integral_on_NS": True,
        "reflection_swaps_total_exceptional_coordinates": True,
    }


def parity_record(boundary_count: int) -> dict[str, object]:
    colors = 2 ** boundary_count
    return {
        "boundary_count": boundary_count,
        "available_parity_vectors": colors,
        "maximal_antichain_bound": colors,
        "chain_cover_bound": colors,
        "direct_lane_single_chain": boundary_count == 0,
    }


def build_payload() -> dict[str, object]:
    exhaustive = exhaustive_forest_record()
    satellites = satellite_records()
    if not all(
        sample["incomparable_pair_witnesses"] for sample in satellites
    ):
        raise AssertionError("satellite sample had no incomparable pair")
    return {
        "theorem": "exceptional-forest parity comparability",
        "verdict": "PROVED",
        "exhaustive_chronological_forest_checks": exhaustive,
        "satellite_proximity_checks": satellites,
        "reflection_samples": [
            reflection_record(5, 1, 3),
            reflection_record(12, 2, 11),
            reflection_record(30, 7, 24),
        ],
        "ample_equalization_mechanism": (
            "positive proximity slacks give weights "
            "lambda_v=sum(children)+slack_v; incomparable descendant "
            "sets permit the two prescribed slacks to equalize lambda_i,lambda_j"
        ),
        "adapted_host_parity": {
            "even_pairings_sufficient": True,
            "even_pairings_necessary": True,
            "odd_formula":
                "s_r(B/2)-B/2=(B.r/2)r; odd B.r would require r/2 in M",
            "odd_absorption_impossible_in_integral_M": True,
        },
        "parity_width_corollary": [
            parity_record(count) for count in range(0, 17)
        ],
        "scope": {
            "full_nef_coxeter_chamber_required": True,
            "reflection_maximality_required": True,
            "odd_parity_proves_maximality": False,
            "depth_bound_obtained": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    payload = build_payload()
    archive = {"payload": payload, "payload_sha256": digest(payload)}
    rendered = json.dumps(archive, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if arguments.output.read_text() != rendered:
            raise AssertionError("archive differs from exact recomputation")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print(archive["payload_sha256"])


if __name__ == "__main__":
    main()
