#!/usr/bin/env python3
"""Freeze the three exact split-prime-13 ``D24+`` prefixes through H=300."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_bna1 import _functional_value, _row_rank
from coxeter_polytopes.search.cubic7_d24_plus_coefficients import (
    analyze_coefficient_prefix,
    d24_coefficient_pair,
)
from coxeter_polytopes.search.cubic7_d24_plus_split13 import (
    P13_D24_SEEDS,
    enumerate_split13_control_levels,
    run_split13_prefix,
    split13_height_zero_patterns,
    split13_prime,
)
from coxeter_polytopes.search.cubic7_d24_plus_unit import (
    d24_height_zero_simple_roots_orthonormal,
)
from coxeter_polytopes.search.cubic7_dn import (
    d24_plus_basis_audit,
    d24_plus_membership,
)


HEIGHT_BOUND = 300
EXPECTED = {
    "cubic7-H24-D24+-a-theta-p13-r7": {
        "controls": 46,
        "accepted": 5,
        "packets": (
            ("19 + 54*theta + 24*theta^2", 1, 1),
            ("38 + 108*theta + 48*theta^2", 5, 4),
        ),
        "snapshots": ((25, 25, 1, 0, 24), (29, 895, 3, 0, 892)),
    },
    "cubic7-H24-D24+-a-theta-p13-r10": {
        "controls": 46,
        "accepted": 1,
        "packets": (("33 + 94*theta + 42*theta^2", 1, 1),),
        "snapshots": ((25, 25, 2, 0, 23),),
    },
    "cubic7-H24-D24+-a-theta-p13-r8": {
        "controls": 47,
        "accepted": 1,
        "packets": (("6 + 18*theta + 8*theta^2", 1, 1),),
        "snapshots": ((25, 25, 2, 0, 23),),
    },
}


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def _vector(vector: tuple[Cubic7, ...]) -> list[dict[str, str]]:
    return [value.to_json() for value in vector]


def _functional(seed, wall):
    return (-seed.timelike_coefficient * wall[0], *wall[1:])


def _residue_root(prime: Cubic7) -> int:
    roots = tuple(
        root
        for root in (7, 8, 10)
        if (
            int(prime.a)
            + int(prime.b) * root
            + int(prime.c) * root * root
        )
        % 13
        == 0
    )
    if len(roots) != 1:
        raise AssertionError("a split-prime generator lost its unique residue root")
    return roots[0]


def _seed_audit(seed) -> dict[str, object]:
    prime = split13_prime(seed)
    expected = EXPECTED[seed.seed_id]
    levels = enumerate_split13_control_levels(seed, HEIGHT_BOUND)
    search = run_split13_prefix(seed, HEIGHT_BOUND)
    geometry = analyze_coefficient_prefix(search)
    nonempty = tuple(packet for packet in search.packets if packet.accepted_roots)
    packet_census = tuple(
        (str(packet.height), packet.candidate_count, len(packet.accepted_roots))
        for packet in nonempty
    )
    snapshots = tuple(
        (
            snapshot.wall_count,
            snapshot.extreme_ray_count,
            snapshot.finite_vertex_count,
            snapshot.cusp_count,
            snapshot.hyperideal_ray_count,
        )
        for snapshot in geometry.snapshots
    )
    if len(levels) != expected["controls"] or len(search.packets) != len(levels):
        raise AssertionError("the split-prime H=300 control census changed")
    if len(search.accepted_roots) != expected["accepted"]:
        raise AssertionError("the split-prime H=300 accepted-root census changed")
    if packet_census != expected["packets"] or snapshots != expected["snapshots"]:
        raise AssertionError("the split-prime H=300 packet geometry changed")

    height_zero = split13_height_zero_patterns(seed)
    if len(height_zero) != 1:
        raise AssertionError("a split-prime seed acquired a new height-zero orbit")
    hz_norm, hz_coset, hz_spatial = height_zero[0]
    if (
        hz_norm != Cubic7(2)
        or hz_coset != "D"
        or hz_spatial[:2] != (Cubic7(1), Cubic7(1))
        or any(hz_spatial[2:])
    ):
        raise AssertionError("the height-zero orbit is not the classical D24 root")

    walls = (
        *d24_height_zero_simple_roots_orthonormal(),
        *(root.vector for root in search.accepted_roots),
    )
    functionals = tuple(_functional(seed, wall) for wall in walls)
    hyperideal = next(ray for ray in geometry.rays if ray.kind == "hyperideal")
    witness_pairings = tuple(
        _functional_value(functional, hyperideal.vector) for functional in functionals
    )
    active = tuple(index for index, value in enumerate(witness_pairings) if not value)
    active_rank = _row_rank(tuple(functionals[index] for index in active), 25)
    if any(value.sign(0) > 0 for value in witness_pairings):
        raise AssertionError("the archived hyperideal witness violates a wall")
    if active_rank != 24 or hyperideal.norm.sign(0) <= 0:
        raise AssertionError("the archived witness is not an extreme hyperideal ray")

    same_height = []
    earlier_walls = list(d24_height_zero_simple_roots_orthonormal())
    for packet in nonempty:
        roots = packet.accepted_roots
        within = tuple(
            d24_coefficient_pair(seed, left.vector, right.vector)
            for index, left in enumerate(roots)
            for right in roots[index + 1 :]
        )
        earlier = tuple(
            d24_coefficient_pair(seed, root.vector, wall)
            for root in roots
            for wall in earlier_walls
        )
        if any(value.sign(0) > 0 for value in (*within, *earlier)):
            raise AssertionError("the archived simultaneous-height selection is not acute")
        same_height.append(
            {
                "height": _value(packet.height),
                "within_packet": {
                    "negative": sum(value.sign(0) < 0 for value in within),
                    "zero": sum(not value for value in within),
                    "positive": 0,
                },
                "against_earlier_walls": {
                    "negative": sum(value.sign(0) < 0 for value in earlier),
                    "zero": sum(not value for value in earlier),
                    "positive": 0,
                },
                "passed": True,
            }
        )
        earlier_walls.extend(root.vector for root in roots)

    p_divisible_roots = tuple(
        root
        for root in search.accepted_roots
        if root.norm in (prime, 2 * prime)
    )
    for root in p_divisible_roots:
        reduced = tuple(value / prime for value in root.vector[1:])
        if not all(value.is_integral for value in reduced):
            raise AssertionError("a p-divisible root escaped the spatial divisibility gate")
        if d24_plus_membership(reduced) not in ("D", "spinor"):
            raise AssertionError("a p-divisible root reduced outside D24+")

    norm_counts = Counter(level.norm for level in levels)
    accepted_norm_counts = Counter(root.norm for root in search.accepted_roots)
    return {
        "seed_id": seed.seed_id,
        "support": seed.support,
        "prime_generator": _value(prime),
        "prime_norm": int(prime.norm()),
        "residue_root_mod_13": _residue_root(prime),
        "timelike_coefficient": _value(seed.timelike_coefficient),
        "canonical_mirror_norms": [_value(norm) for norm in seed.canonical_mirror_norms],
        "complete_local_norm_proof": seed.complete_norm_class_proof,
        "p_divisible_spatial_gate": {
            "statement": (
                "For q=p or 2p, reflection integrality is equivalent to y in pD24+; "
                "write y=pz and enumerate z in D24+ at target/p^2."
            ),
            "p_divisible_control_count": sum(
                level.norm in (prime, 2 * prime) for level in levels
            ),
            "p_divisible_accepted_root_count": len(p_divisible_roots),
            "every_accepted_p_divisible_spatial_vector_reduces_to_D24_plus": True,
        },
        "height_zero_audit": {
            "dominant_reflective_orbit_count": 1,
            "norm": _value(hz_norm),
            "coset": hz_coset,
            "spatial_vector": _vector(hz_spatial),
            "complete_chamber": "D24 with 24 classical simple roots",
        },
        "control_level_count": len(levels),
        "control_counts_by_norm": [
            {"norm": _value(norm), "count": norm_counts[norm]}
            for norm in seed.canonical_mirror_norms
        ],
        "control_levels": [
            {
                "height": _value(level.height),
                "norm": _value(level.norm),
                "control": _value(level.control),
                "spatial_square_target": _value(level.target),
                "p_reduced_target": (
                    _value(level.target / (prime * prime))
                    if level.norm in (prime, 2 * prime)
                    and (level.target / (prime * prime)).is_integral
                    else None
                ),
            }
            for level in levels
        ],
        "height_packet_count": len(search.packets),
        "height_packets": [
            {
                "height": _value(packet.height),
                "candidate_count": packet.candidate_count,
                "accepted_count": len(packet.accepted_roots),
                "accepted_cosets": {
                    coset: sum(root.coset == coset for root in packet.accepted_roots)
                    for coset in sorted({root.coset for root in packet.accepted_roots})
                },
            }
            for packet in search.packets
        ],
        "nonempty_packet_count": len(nonempty),
        "candidate_root_count": sum(packet.candidate_count for packet in search.packets),
        "accepted_root_count": len(search.accepted_roots),
        "accepted_norm_counts": [
            {"norm": _value(norm), "count": accepted_norm_counts[norm]}
            for norm in seed.canonical_mirror_norms
            if accepted_norm_counts[norm]
        ],
        "spinor_accepted_count": sum(root.coset != "D" for root in search.accepted_roots),
        "accepted_roots": [
            {
                "height": _value(root.height),
                "norm": _value(root.norm),
                "control": _value(root.control),
                "spatial_square_target": _value(root.target),
                "coset": root.coset,
                "vector_orthonormal": _vector(root.vector),
            }
            for root in search.accepted_roots
        ],
        "same_height_audit": same_height,
        "geometry": {
            "engine": "exact double description; adjacency by common-active rank 23",
            "snapshots": [
                {
                    "height": _value(snapshot.height),
                    "wall_count": snapshot.wall_count,
                    "extreme_ray_count": snapshot.extreme_ray_count,
                    "finite_vertex_count": snapshot.finite_vertex_count,
                    "cusp_count": snapshot.cusp_count,
                    "hyperideal_ray_count": snapshot.hyperideal_ray_count,
                }
                for snapshot in geometry.snapshots
            ],
            "final": {
                "wall_count": len(walls),
                "extreme_ray_count": len(geometry.rays),
                "finite_vertex_count": geometry.finite_vertex_count,
                "cusp_count": geometry.cusp_count,
                "hyperideal_ray_count": geometry.hyperideal_ray_count,
                "finite_volume": geometry.finite_volume,
                "compact": geometry.compact,
                "ray_stream_sha256": geometry.ray_stream_sha256,
            },
            "hyperideal_witness": {
                "vector": _vector(hyperideal.vector),
                "norm": _value(hyperideal.norm),
                "active_facets": list(active),
                "active_rank": active_rank,
                "wall_pairing_census": {
                    "negative": sum(value.sign(0) < 0 for value in witness_pairings),
                    "zero": sum(not value for value in witness_pairings),
                    "positive": 0,
                },
            },
        },
        "frontier": {
            "last_enumerated_height": _value(levels[-1].height),
            "last_enumerated_norm": _value(levels[-1].norm),
            "last_enumerated_control": _value(levels[-1].control),
            "last_enumerated_target": _value(levels[-1].target),
            "all_controls_at_height_at_most_300_exhausted": True,
        },
    }


def build_audit() -> dict[str, object]:
    if not d24_plus_basis_audit():
        raise AssertionError("the D24+ basis lost even unimodularity")
    seeds = [_seed_audit(seed) for seed in P13_D24_SEEDS]
    return {
        "schema_version": "0.1.0",
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "hyperbolic_dimension": 24,
        "spatial_lattice": "D24+ (even unimodular Niemeier lattice)",
        "height_bound": HEIGHT_BOUND,
        "seed_order": [seed.seed_id for seed in P13_D24_SEEDS],
        "galois_scope": (
            "The three prime ideals are Galois-conjugate, but their totally-positive "
            "generators and theta*p timelike coefficients are retained separately "
            "because embedding zero is the distinguished Lorentzian place."
        ),
        "seeds": seeds,
        "comparison": {
            "finite_volume_survivors": [],
            "compact_survivors": [],
            "fewest_final_hyperideal_rays": min(
                seed["geometry"]["final"]["hyperideal_ray_count"] for seed in seeds
            ),
            "best_bounded_prefixes": [
                seed["seed_id"]
                for seed in seeds
                if seed["geometry"]["final"]["hyperideal_ray_count"]
                == min(
                    other["geometry"]["final"]["hyperideal_ray_count"]
                    for other in seeds
                )
            ],
        },
        "discovery": None,
        "scope": (
            "For each of the three distinguished split-prime-13 directions, every "
            "canonical mirror norm, compatible control, D24 root coset, and all "
            "seven nonzero D24+ spinor residues are exhausted through exact Vinberg "
            "height H=300. No global nonreflectivity conclusion is claimed."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/cubic7_split13_d24_plus_h300.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "height_bound": payload["height_bound"],
                "comparison": payload["comparison"],
                "seeds": [
                    {
                        "seed_id": seed["seed_id"],
                        "controls": seed["control_level_count"],
                        "accepted_roots": seed["accepted_root_count"],
                        "geometry": seed["geometry"]["final"],
                    }
                    for seed in payload["seeds"]
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
