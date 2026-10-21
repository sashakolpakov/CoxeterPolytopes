#!/usr/bin/env python3
"""Freeze the exact H=1000 prefix for ``<-theta(2-theta)> + D24+``.

Floating point linear programming is used only to propose an active set for
one ray.  The archived certificate is then rebuilt and checked entirely in
``Q(theta)``: all wall pairings, the active rank, and the Lorentz norm are
exact.  No numerical conclusion is trusted by the audit.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.cubic7 import CUBIC7_ONE, CUBIC7_ZERO, Cubic7
from coxeter_polytopes.search.cubic7_bna1 import (
    _functional_value,
    _row_rank,
)
from coxeter_polytopes.search.cubic7_d24_plus_coefficients import (
    DELTA7,
    P7_D24_SEED,
    analyze_coefficient_prefix,
    d24_coefficient_pair,
    d24_height_zero_simple_roots_orthonormal,
    enumerate_coefficient_control_levels,
    run_coefficient_prefix,
)


HEIGHT_BOUND = 1000
CERTIFIED_WITNESS = (
    Cubic7(Fraction(23, 14), Fraction(5, 14), Fraction(-3, 14)),
    CUBIC7_ONE,
    CUBIC7_ONE,
    CUBIC7_ONE,
    *(CUBIC7_ZERO for _ in range(21)),
)
CERTIFIED_ACTIVE_FACETS = (
    0,
    1,
    *range(3, 24),
    62,
    142,
)


def _value(value: Cubic7) -> dict[str, object]:
    return {
        "exact": value.to_json(),
        "display": str(value),
        "signs": [value.sign(place) for place in range(3)],
    }


def _vector(vector: tuple[Cubic7, ...]) -> list[dict[str, str]]:
    return [value.to_json() for value in vector]


def _functional(wall: tuple[Cubic7, ...]) -> tuple[Cubic7, ...]:
    return (-P7_D24_SEED.timelike_coefficient * wall[0], *wall[1:])


def build_audit() -> dict[str, object]:
    seed = P7_D24_SEED
    search = run_coefficient_prefix(seed, HEIGHT_BOUND)
    levels = enumerate_coefficient_control_levels(seed, HEIGHT_BOUND)
    nonempty = tuple(packet for packet in search.packets if packet.accepted_roots)
    if len(levels) != 94 or len(search.packets) != 94:
        raise AssertionError("the frozen H=1000 control census changed")
    if tuple(len(packet.accepted_roots) for packet in nonempty) != (1, 77, 77):
        raise AssertionError("the frozen H=1000 packet census changed")
    if len(search.accepted_roots) != 155:
        raise AssertionError("the frozen H=1000 wall census changed")

    height_zero = d24_height_zero_simple_roots_orthonormal()
    walls = (*height_zero, *(root.vector for root in search.accepted_roots))
    functionals = tuple(_functional(wall) for wall in walls)

    old_search = run_coefficient_prefix(seed, 300)
    old_geometry = analyze_coefficient_prefix(old_search)
    old_hyperideal = tuple(
        ray.vector for ray in old_geometry.rays if ray.kind == "hyperideal"
    )
    if len(old_hyperideal) != 23:
        raise AssertionError("the H=300 hyperideal-ray baseline changed")
    packet2 = tuple(root.vector for root in nonempty[1].accepted_roots)
    packet3 = tuple(root.vector for root in nonempty[2].accepted_roots)
    old_ray_audit = []
    for index, ray in enumerate(old_hyperideal):
        values2 = tuple(d24_coefficient_pair(seed, ray, wall) for wall in packet2)
        values3 = tuple(d24_coefficient_pair(seed, ray, wall) for wall in packet3)
        old_ray_audit.append(
            {
                "old_hyperideal_ray_index": index,
                "vector": _vector(ray),
                "norm": _value(d24_coefficient_pair(seed, ray, ray)),
                "packet_2_positive_pairing_count": sum(
                    value.sign(0) > 0 for value in values2
                ),
                "packet_2_zero_pairing_count": sum(not value for value in values2),
                "packet_2_max_pairing": _value(max(values2)),
                "survives_packet_2": all(value.sign(0) <= 0 for value in values2),
                "packet_3_positive_pairing_count": sum(
                    value.sign(0) > 0 for value in values3
                ),
                "packet_3_zero_pairing_count": sum(not value for value in values3),
                "packet_3_max_pairing": _value(max(values3)),
                "survives_both_new_packets": all(
                    value.sign(0) <= 0 for value in (*values2, *values3)
                ),
            }
        )

    witness = CERTIFIED_WITNESS
    pairings = tuple(_functional_value(functional, witness) for functional in functionals)
    active = tuple(index for index, value in enumerate(pairings) if not value)
    norm = d24_coefficient_pair(seed, witness, witness)
    active_rank = _row_rank(tuple(functionals[index] for index in active), 25)
    if any(value.sign(0) > 0 for value in pairings):
        raise AssertionError("the replacement witness violates a frozen wall")
    if active_rank != 24 or norm.sign(0) <= 0:
        raise AssertionError("the replacement witness is not an extreme hyperideal ray")
    if active != CERTIFIED_ACTIVE_FACETS:
        raise AssertionError("the replacement witness active set changed")

    same_height = []
    earlier_walls = list(height_zero)
    for packet in nonempty:
        roots = packet.accepted_roots
        within_values = tuple(
            d24_coefficient_pair(seed, left.vector, right.vector)
            for index, left in enumerate(roots)
            for right in roots[index + 1 :]
        )
        earlier_values = tuple(
            d24_coefficient_pair(seed, root.vector, wall)
            for root in roots
            for wall in earlier_walls
        )
        same_height.append(
            {
                "height": _value(packet.height),
                "accepted_count": len(roots),
                "pair_count": len(roots) * (len(roots) - 1) // 2,
                "negative_pairing_count": sum(value.sign(0) < 0 for value in within_values),
                "zero_pairing_count": sum(not value for value in within_values),
                "positive_pairing_count": sum(value.sign(0) > 0 for value in within_values),
                "maximum_pairing": _value(max(within_values)) if within_values else None,
                "earlier_wall_pairing_census": {
                    "negative": sum(value.sign(0) < 0 for value in earlier_values),
                    "zero": sum(not value for value in earlier_values),
                    "positive": sum(value.sign(0) > 0 for value in earlier_values),
                },
                "passed": all(value.sign(0) <= 0 for value in (*within_values, *earlier_values)),
            }
        )
        earlier_walls.extend(root.vector for root in roots)

    return {
        "schema_version": "0.1.0",
        "seed_id": seed.seed_id,
        "field": "Q(theta), theta^3+theta^2-2theta-1=0",
        "lattice": "<-theta(2-theta)> + D24+",
        "hyperbolic_dimension": 24,
        "spatial_lattice": "D24+ (even unimodular Niemeier lattice)",
        "height_bound": HEIGHT_BOUND,
        "timelike_coefficient": _value(seed.timelike_coefficient),
        "ideal_support": seed.support,
        "canonical_mirror_norms": [_value(norm) for norm in seed.canonical_mirror_norms],
        "complete_mirror_norm_proof": seed.complete_norm_class_proof,
        "delta_divisible_level_reduction": {
            "delta": _value(DELTA7),
            "statement": (
                "For q=delta or 2*delta, self-duality forces y in delta*D24+. "
                "The exact enumerator writes y=delta*z and searches D24+ at "
                "target/delta^2, including the root and all seven spinor residues."
            ),
        },
        "control_level_count": len(levels),
        "height_packet_count": len(search.packets),
        "control_levels": [
            {
                "height": _value(level.height),
                "norm": _value(level.norm),
                "control": _value(level.control),
                "spatial_square_target": _value(level.target),
            }
            for level in levels
        ],
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
        "accepted_root_count": len(search.accepted_roots),
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
        "old_h300_hyperideal_ray_audit": {
            "old_ray_count": len(old_hyperideal),
            "survives_packet_2_count": sum(row["survives_packet_2"] for row in old_ray_audit),
            "survives_both_new_packets_count": sum(
                row["survives_both_new_packets"] for row in old_ray_audit
            ),
            "rays": old_ray_audit,
        },
        "certified_hyperideal_extreme_ray": {
            "vector": _vector(witness),
            "norm": _value(norm),
            "active_facets": list(active),
            "active_rank": active_rank,
            "wall_pairing_census": {
                "negative": sum(value.sign(0) < 0 for value in pairings),
                "zero": sum(not value for value in pairings),
                "positive": sum(value.sign(0) > 0 for value in pairings),
            },
            "wall_count": len(walls),
            "proposal_method": (
                "deterministic floating-point LP proposed an active set only; "
                "the displayed ray was then reconstructed from an exact rank-24 "
                "kernel, and this replay rechecks its exact active rank, norm, "
                "and every wall pairing"
            ),
            "proposal_label": "positive third spatial coordinate",
            "certificate_is_entirely_exact": True,
        },
        "geometry_conclusion": {
            "wall_count": len(walls),
            "finite_volume": False,
            "compact": False,
            "reason": "an exact rank-24 hyperideal extreme ray satisfies every H<=1000 wall",
            "full_double_description_claimed": False,
        },
        "frontier": {
            "last_enumerated_height": _value(levels[-1].height),
            "last_enumerated_norm": _value(levels[-1].norm),
            "last_enumerated_control": _value(levels[-1].control),
            "last_enumerated_target": _value(levels[-1].target),
            "all_controls_at_height_at_most_1000_exhausted": True,
        },
        "discovery": None,
        "scope": (
            "Every canonical mirror-norm class, compatible control, D24 root coset, "
            "and all seven nonzero D24+ spinor glue residues are exhausted through "
            "exact Vinberg height H=1000. The certified open prefix does not imply "
            "global nonreflectivity."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/enumerations/cubic7_p7_d24_plus_h1000.json"),
    )
    args = parser.parse_args()
    payload = build_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "controls": payload["control_level_count"],
                "accepted_roots": payload["accepted_root_count"],
                "old_ray_survival": {
                    key: payload["old_h300_hyperideal_ray_audit"][key]
                    for key in (
                        "old_ray_count",
                        "survives_packet_2_count",
                        "survives_both_new_packets_count",
                    )
                },
                "certificate": payload["certified_hyperideal_extreme_ray"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
