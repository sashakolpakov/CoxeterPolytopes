#!/usr/bin/env python3
"""Extend the frozen exact silver prime-seven D24+ census to H=500.

The dense equal-height packets are checked with integer coefficient matrices.
If every coordinate is written over a common positive denominator, a
Q(sqrt(2)) pairing is represented by two integer matrices.  This is exact;
the explicit overflow bounds below are checked before every NumPy product.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from hashlib import sha256
from itertools import groupby
from math import lcm
from pathlib import Path
from typing import Sequence

import numpy as np

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_p7_d24_plus import (
    H300_REPLACEMENT_WITNESS_VECTOR,
    H500_REPLACEMENT_WITNESS_VECTOR,
    P7_TIMELIKE,
    SILVER_P7_D24_SEED,
    SilverP7D24Candidate,
    certify_h500_replacement_hyperideal_ray,
    enumerate_p7_d24_control_levels,
    enumerate_p7_d24_level,
    silver_d24_height_zero_simple_roots,
    silver_p7_d24_pair,
)


ROOT = Path(__file__).resolve().parents[1]
H300_ARTIFACT = ROOT / "data/enumerations/silver_p7_d24_plus_h300.json"
H300_ARTIFACT_SHA256 = "4ea87341148489598dc85b3cc8271965a036e5ffba2f24052784705e4af614cf"
INT64_MAX = np.iinfo(np.int64).max
MAX_EXACT_SIGN_COEFFICIENT = 2_000_000_000


@dataclass(slots=True)
class ExactMatrixAudit:
    pairing_count: int = 0
    same_height_unordered_pair_count: int = 0
    maximum_proved_intermediate_bound: int = 0
    maximum_abs_pair_coefficient: int = 0
    zero_boundary_scalar_replay_count: int = 0
    strict_scalar_replay_count: int = 0


MATRIX_AUDIT = ExactMatrixAudit()


def _silver(value: dict[str, str]) -> Silver:
    return Silver.from_json(value)


def _vector(values: Sequence[dict[str, str]]) -> tuple[Silver, ...]:
    return tuple(_silver(value) for value in values)


def _value(value: Silver) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in (0, 1)],
    }


def _serialize_vector(vector: Sequence[Silver]) -> list[dict[str, str]]:
    return [value.to_json() for value in vector]


def _root(root: SilverP7D24Candidate) -> dict[str, object]:
    return {
        "norm": _value(root.norm),
        "height": _value(root.height),
        "control": _value(root.control),
        "target": _value(root.target),
        "coset": root.coset,
        "vector_orthonormal": _serialize_vector(root.vector),
    }


def _root_stream_sha256(roots: Sequence[SilverP7D24Candidate]) -> str:
    digest = sha256()
    for root in roots:
        digest.update(
            json.dumps(_root(root), sort_keys=True, separators=(",", ":")).encode()
        )
        digest.update(b"\n")
    return digest.hexdigest()


def _frozen_h300() -> tuple[dict[str, object], tuple[tuple[Silver, ...], ...]]:
    if sha256(H300_ARTIFACT.read_bytes()).hexdigest() != H300_ARTIFACT_SHA256:
        raise AssertionError("the pinned H=300 input artifact changed")
    payload = json.loads(H300_ARTIFACT.read_text())
    roots = tuple(
        _vector(root["vector_orthonormal"])
        for packet in payload["enumeration"]["packets"]
        for root in packet["accepted_roots"]
    )
    if len(roots) != 324:
        raise AssertionError("the frozen H=300 accepted-root census changed")
    return payload, (*silver_d24_height_zero_simple_roots(), *roots)


def _coefficient_matrices(
    vectors: Sequence[Sequence[Silver]],
) -> tuple[np.ndarray, np.ndarray, int]:
    denominator = 1
    for vector in vectors:
        for value in vector:
            denominator = lcm(
                denominator, value.a.denominator, value.b.denominator
            )
    rational = np.asarray(
        [
            [
                int(value.a.numerator * (denominator // value.a.denominator))
                for value in vector
            ]
            for vector in vectors
        ],
        dtype=np.int64,
    )
    radical = np.asarray(
        [
            [
                int(value.b.numerator * (denominator // value.b.denominator))
                for value in vector
            ]
            for vector in vectors
        ],
        dtype=np.int64,
    )
    return rational, radical, denominator


def _maximum_abs(matrix: np.ndarray) -> int:
    return int(np.max(np.abs(matrix), initial=0))


def _pair_coefficients(
    left_a: np.ndarray,
    left_b: np.ndarray,
    right_a: np.ndarray,
    right_b: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return exact numerator coefficients for all Lorentzian pairings."""

    la_spatial = _maximum_abs(left_a[:, 1:])
    lb_spatial = _maximum_abs(left_b[:, 1:])
    ra_spatial = _maximum_abs(right_a[:, 1:])
    rb_spatial = _maximum_abs(right_b[:, 1:])
    la_control = _maximum_abs(left_a[:, :1])
    lb_control = _maximum_abs(left_b[:, :1])
    ra_control = _maximum_abs(right_a[:, :1])
    rb_control = _maximum_abs(right_b[:, :1])
    spatial_a_bound = 24 * (
        la_spatial * ra_spatial + 2 * lb_spatial * rb_spatial
    )
    spatial_b_bound = 24 * (
        la_spatial * rb_spatial + lb_spatial * ra_spatial
    )
    control_a_bound = la_control * ra_control + 2 * lb_control * rb_control
    control_b_bound = la_control * rb_control + lb_control * ra_control
    output_bound = max(
        spatial_a_bound + control_a_bound + 4 * control_b_bound,
        spatial_b_bound + control_b_bound + 2 * control_a_bound,
    )
    MATRIX_AUDIT.maximum_proved_intermediate_bound = max(
        MATRIX_AUDIT.maximum_proved_intermediate_bound, output_bound
    )
    if output_bound > INT64_MAX:
        raise OverflowError("integer pairing formula can overflow int64")

    spatial_a = left_a[:, 1:] @ right_a[:, 1:].T
    spatial_a += 2 * (left_b[:, 1:] @ right_b[:, 1:].T)
    spatial_b = left_a[:, 1:] @ right_b[:, 1:].T
    spatial_b += left_b[:, 1:] @ right_a[:, 1:].T

    la = left_a[:, :1]
    lb = left_b[:, :1]
    ra = right_a[:, :1].T
    rb = right_b[:, :1].T
    control_a = la * ra + 2 * lb * rb
    control_b = la * rb + lb * ra
    # -(-1+2sqrt(2))*(control_a+control_b*sqrt(2)).
    output_a = spatial_a + control_a - 4 * control_b
    output_b = spatial_b + control_b - 2 * control_a
    maximum = max(
        int(np.max(np.abs(output_a), initial=0)),
        int(np.max(np.abs(output_b), initial=0)),
    )
    if maximum > MAX_EXACT_SIGN_COEFFICIENT:
        raise OverflowError("pair coefficient too large for exact int64 sign test")
    MATRIX_AUDIT.pairing_count += output_a.size
    MATRIX_AUDIT.maximum_abs_pair_coefficient = max(
        MATRIX_AUDIT.maximum_abs_pair_coefficient, maximum
    )
    return output_a, output_b


def _positive_mask(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Exact identity-embedding sign test for integer a+b*sqrt(2)."""

    field_norm = a * a - 2 * b * b
    return (
        ((a > 0) & (b >= 0))
        | ((a > 0) & (b < 0) & (field_norm > 0))
        | ((a < 0) & (b > 0) & (field_norm < 0))
        | ((a == 0) & (b > 0))
    )


def _survivor_indices(
    candidates: Sequence[SilverP7D24Candidate],
    walls: Sequence[Sequence[Silver]],
    *,
    block_size: int = 256,
) -> tuple[int, ...]:
    candidate_a, candidate_b, candidate_denominator = _coefficient_matrices(
        tuple(root.vector for root in candidates)
    )
    wall_a, wall_b, wall_denominator = _coefficient_matrices(walls)
    if candidate_denominator <= 0 or wall_denominator <= 0:
        raise AssertionError("coordinate denominators must be positive")
    survivors = []
    boundary_pairs: list[tuple[int, int]] = []
    for start in range(0, len(candidates), block_size):
        stop = min(start + block_size, len(candidates))
        a, b = _pair_coefficients(
            candidate_a[start:stop], candidate_b[start:stop], wall_a, wall_b
        )
        mask = ~np.any(_positive_mask(a, b), axis=1)
        local_survivors = np.flatnonzero(mask)
        survivors.extend(start + index for index in local_survivors)
        zero = (a[local_survivors] == 0) & (b[local_survivors] == 0)
        boundary_pairs.extend(
            (start + int(local_index), int(wall_index))
            for survivor_row, local_index in enumerate(local_survivors)
            for wall_index in np.flatnonzero(zero[survivor_row])
        )
    # Boundary replay uses the independent scalar field implementation.  The
    # strict signs themselves already come from exact integer arithmetic, not
    # floating-point screening.
    if any(
        silver_p7_d24_pair(candidates[index].vector, walls[wall_index])
        for index, wall_index in boundary_pairs
    ):
        raise AssertionError("integer-matrix boundary disagrees with scalar replay")
    MATRIX_AUDIT.zero_boundary_scalar_replay_count += len(boundary_pairs)
    # Also replay two strict inequalities for every retained root so a bad
    # coordinate conversion cannot hide behind an empty boundary set.
    replay_walls = tuple(dict.fromkeys((0, len(walls) // 2, len(walls) - 1)))
    if any(
        silver_p7_d24_pair(candidates[index].vector, walls[wall_index]).sign(0) > 0
        for index in survivors
        for wall_index in replay_walls
    ):
        raise AssertionError("integer-matrix filtering disagrees with scalar replay")
    MATRIX_AUDIT.strict_scalar_replay_count += len(survivors) * len(replay_walls)
    return tuple(survivors)


def _assert_same_height_acute(
    roots: Sequence[SilverP7D24Candidate], *, block_size: int = 256
) -> None:
    rational, radical, denominator = _coefficient_matrices(
        tuple(root.vector for root in roots)
    )
    if denominator <= 0:
        raise AssertionError("coordinate denominator must be positive")
    MATRIX_AUDIT.same_height_unordered_pair_count += len(roots) * (len(roots) - 1) // 2
    for start in range(0, len(roots), block_size):
        stop = min(start + block_size, len(roots))
        a, b = _pair_coefficients(
            rational[start:stop], radical[start:stop], rational, radical
        )
        positive = _positive_mask(a, b)
        for local_index in range(stop - start):
            positive[local_index, : start + local_index + 1] = False
        locations = np.argwhere(positive)
        if locations.size:
            left = start + int(locations[0, 0])
            right = int(locations[0, 1])
            exact = silver_p7_d24_pair(roots[left].vector, roots[right].vector)
            raise AssertionError(
                f"same-height roots {left},{right} have positive pairing {exact}"
            )


def build_audit() -> dict[str, object]:
    global MATRIX_AUDIT
    MATRIX_AUDIT = ExactMatrixAudit()
    h300, frozen_walls = _frozen_h300()
    walls = list(frozen_walls)
    levels = tuple(
        level
        for level in enumerate_p7_d24_control_levels(500)
        if (level.height - Silver(300)).sign(0) > 0
    )
    print(f"new_control_levels={len(levels)} initial_wall_count={len(walls)}", flush=True)
    if len(levels) != 16:
        raise AssertionError("the H=(300,500] control census changed")
    accepted: list[SilverP7D24Candidate] = []
    packet_data = []
    all_enumerations = []
    for height, iterator in groupby(levels, key=lambda level: level.height):
        packet_levels = tuple(iterator)
        enumerations = tuple(enumerate_p7_d24_level(level) for level in packet_levels)
        all_enumerations.extend(enumerations)
        candidates = tuple(
            root for enumeration in enumerations for root in enumeration.candidates
        )
        print(
            f"height={height} levels={len(packet_levels)} candidates={len(candidates)}",
            flush=True,
        )
        indices = _survivor_indices(candidates, walls)
        survivors = tuple(candidates[index] for index in indices)
        print(f"  survivors={len(survivors)} checking_equal_height", flush=True)
        _assert_same_height_acute(survivors)
        walls.extend(root.vector for root in survivors)
        accepted.extend(survivors)
        packet_data.append((height, enumerations, candidates, survivors))
        print(f"  wall_count={len(walls)}", flush=True)
    print(
        f"new_accepted={len(accepted)} final_wall_count={len(walls)}",
        flush=True,
    )
    if len(accepted) != 4026 or len(walls) != 4374:
        raise AssertionError("the exact H=500 wall census changed")

    h300_witness_new_pairings = tuple(
        silver_p7_d24_pair(H300_REPLACEMENT_WITNESS_VECTOR, root.vector)
        for root in accepted
    )
    if sum(value.sign(0) > 0 for value in h300_witness_new_pairings) != 74:
        raise AssertionError("the H=300 replacement-ray cut census changed")

    replacement = certify_h500_replacement_hyperideal_ray(tuple(walls))
    if replacement.vector != H500_REPLACEMENT_WITNESS_VECTOR:
        raise AssertionError("the H=500 replacement vector changed")
    coordinate = replacement.vector[1]
    bounds = []
    for wall_index, wall in enumerate(walls[24:], start=24):
        partial_sum = sum(wall[1:6], start=Silver(0))
        if partial_sum.sign(0) > 0:
            bounds.append((P7_TIMELIKE * wall[0] / partial_sum, wall_index))
    minimum_bound, limiting_wall_index = min(bounds)
    if minimum_bound != coordinate or limiting_wall_index != 4361:
        raise AssertionError("the exact five-coordinate limiter changed")
    new_root_by_wall = {
        348 + index: root for index, root in enumerate(accepted)
    }
    limiter = new_root_by_wall[limiting_wall_index]

    extension_packets = []
    for height, enumerations, candidates, survivors in packet_data:
        extension_packets.append(
            {
                "height": _value(height),
                "control_level_count": len(enumerations),
                "raw_reflective_count": sum(
                    item.reflection_integral_count for item in enumerations
                ),
                "imprimitive_sqrt2_count": sum(
                    item.imprimitive_sqrt2_count for item in enumerations
                ),
                "primitive_candidate_count": len(candidates),
                "accepted_count": len(survivors),
                "accepted_spinor_count": sum(
                    root.coset != "D24" for root in survivors
                ),
                "accepted_root_stream_sha256": _root_stream_sha256(survivors),
                "accepted_roots": [_root(root) for root in survivors],
            }
        )

    old_enumeration = h300["enumeration"]
    return {
        "schema_version": "0.1.0",
        "field": "Q(sqrt(2)), ring of integers Z[sqrt(2)]",
        "lattice": "<-(2*sqrt(2)-1)> + D24+",
        "hyperbolic_dimension": 24,
        "height_bound": 500,
        "seed_id": SILVER_P7_D24_SEED.seed_id,
        "frozen_h300_lineage": {
            "artifact": str(H300_ARTIFACT.relative_to(ROOT)),
            "sha256": H300_ARTIFACT_SHA256,
            "accepted_root_count": old_enumeration["accepted_root_count"],
            "final_wall_count": old_enumeration["final_wall_count"],
        },
        "canonical_mirror_norms": [
            _value(norm) for norm in SILVER_P7_D24_SEED.canonical_mirror_norms
        ],
        "complete_mirror_norm_proof": SILVER_P7_D24_SEED.complete_norm_class_proof,
        "enumeration": {
            "control_level_count": old_enumeration["control_level_count"] + len(levels),
            "height_packet_count": old_enumeration["height_packet_count"] + len(packet_data),
            "accepted_root_count": old_enumeration["accepted_root_count"] + len(accepted),
            "final_wall_count": len(walls),
            "extension_control_level_count": len(levels),
            "extension_accepted_root_count": len(accepted),
            "extension_accepted_root_stream_sha256": _root_stream_sha256(accepted),
            "extension_control_levels": [
                {
                    "norm": _value(item.level.norm),
                    "height": _value(item.level.height),
                    "control": _value(item.level.control),
                    "spatial_square_target": _value(item.level.target),
                    "reduced_target_after_p7_divisibility": (
                        None if item.reduced_target is None else _value(item.reduced_target)
                    ),
                    "glue_pattern_count": item.glue_pattern_count,
                    "reflection_integral_count_before_primitivity": item.reflection_integral_count,
                    "imprimitive_sqrt2_count": item.imprimitive_sqrt2_count,
                    "primitive_candidate_count": len(item.candidates),
                }
                for item in all_enumerations
            ],
            "extension_packets": extension_packets,
            "new_nonempty_packet_summary": [
                {
                    "height": packet["height"],
                    "raw_reflective_count": packet["raw_reflective_count"],
                    "imprimitive_sqrt2_count": packet["imprimitive_sqrt2_count"],
                    "primitive_candidate_count": packet["primitive_candidate_count"],
                    "accepted_count": packet["accepted_count"],
                    "accepted_spinor_count": packet["accepted_spinor_count"],
                }
                for packet in extension_packets
                if packet["accepted_count"]
            ],
        },
        "exact_dense_packet_accelerator": {
            "method": (
                "Coordinates are lifted over positive common denominators to "
                "integer coefficient matrices A+B*sqrt(2). Lorentz pairings "
                "are computed by exact int64 matrix products. A symbolic "
                "absolute-value bound is checked before every product and "
                "combination. Signs use the exact comparison A^2 versus "
                "2*B^2; no floating-point value enters selection. All zero "
                "boundaries and three strict inequalities per survivor are "
                "replayed through the independent scalar Silver arithmetic."
            ),
            "pairing_matrix_entry_count": MATRIX_AUDIT.pairing_count,
            "same_height_unordered_pair_count": MATRIX_AUDIT.same_height_unordered_pair_count,
            "maximum_proved_intermediate_absolute_bound": MATRIX_AUDIT.maximum_proved_intermediate_bound,
            "int64_max": int(INT64_MAX),
            "maximum_abs_pair_coefficient": MATRIX_AUDIT.maximum_abs_pair_coefficient,
            "exact_sign_square_safety_limit": MAX_EXACT_SIGN_COEFFICIENT,
            "zero_boundary_scalar_replay_count": MATRIX_AUDIT.zero_boundary_scalar_replay_count,
            "strict_scalar_replay_count": MATRIX_AUDIT.strict_scalar_replay_count,
            "simultaneous_height_semantics": True,
        },
        "h300_replacement_ray_cut_audit": {
            "vector": _serialize_vector(H300_REPLACEMENT_WITNESS_VECTOR),
            "new_wall_count": len(h300_witness_new_pairings),
            "negative_pairing_count": sum(
                value.sign(0) < 0 for value in h300_witness_new_pairings
            ),
            "zero_pairing_count": sum(not value for value in h300_witness_new_pairings),
            "positive_pairing_count": sum(
                value.sign(0) > 0 for value in h300_witness_new_pairings
            ),
            "minimum_pairing": _value(min(h300_witness_new_pairings)),
            "maximum_pairing": _value(max(h300_witness_new_pairings)),
        },
        "h500_replacement_hyperideal_ray": {
            "construction": (
                "v=(1;c,c,c,c,c,0^19), where c is the minimum of "
                "(2*sqrt(2)-1)*x/(y1+...+y5) over positive-height walls "
                "with positive denominator"
            ),
            "coordinate": _value(coordinate),
            "vector": _serialize_vector(replacement.vector),
            "norm": _value(replacement.norm),
            "wall_count": replacement.wall_count,
            "negative_pairing_count": sum(
                value.sign(0) < 0 for value in replacement.pairings
            ),
            "zero_pairing_count": sum(not value for value in replacement.pairings),
            "positive_pairing_count": sum(
                value.sign(0) > 0 for value in replacement.pairings
            ),
            "minimum_pairing": _value(min(replacement.pairings)),
            "maximum_pairing": _value(max(replacement.pairings)),
            "active_facets": list(replacement.active_facets),
            "active_rank": replacement.active_rank,
            "limiting_wall_index": limiting_wall_index,
            "limiting_wall": _root(limiter),
            "all_wall_pairings_nonpositive": True,
            "finite_volume": False,
            "compact": False,
        },
        "geometry_status": {
            "h300_replacement_ray_survives": False,
            "h500_finite_volume": False,
            "h500_compact": False,
            "eventual_chamber_status": "unresolved",
        },
        "discovery": None,
        "scope": (
            "The independently frozen exact H=300 prefix is extended through "
            "algebraic height 500 using all 16 newly admissible controls, all "
            "D24 root/spinor glue classes, the complete four canonical mirror "
            "norms, and simultaneous-height Vinberg selection. The resulting "
            "4,374-wall bounded prefix has an explicit exact rank-24 "
            "hyperideal ray, hence this prefix has infinite volume. Walls "
            "above height 500 may cut it, so the eventual chamber remains "
            "unresolved and no global nonreflectivity claim is made."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/silver_p7_d24_plus_h500.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "control_level_count": payload["enumeration"]["control_level_count"],
                "accepted_root_count": payload["enumeration"]["accepted_root_count"],
                "final_wall_count": payload["enumeration"]["final_wall_count"],
                "new_nonempty_packet_counts": [
                    packet["accepted_count"]
                    for packet in payload["enumeration"]["new_nonempty_packet_summary"]
                ],
                "replacement_norm": payload["h500_replacement_hyperideal_ray"]["norm"],
                "eventual_chamber_status": payload["geometry_status"]["eventual_chamber_status"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
