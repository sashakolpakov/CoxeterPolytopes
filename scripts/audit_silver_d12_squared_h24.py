#!/usr/bin/env python3
"""Freeze three exact silver N(D12^2) Vinberg prefixes."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_d12_squared import (
    SQRT2_D12_SQUARED_SEED,
    TWO_SQRT2_D12_SQUARED_SEED,
    UNIT_D12_SQUARED_SEED,
    d12_squared_basis_audit,
    d12_squared_seed_pair,
    first_packet_geometry,
    height_zero_seed_dominant_roots,
    run_d12_squared_seed_prefix,
    strict_hyperideal_witness,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/silver_d12_squared_h24.json"

FROZEN_SEARCHES = (
    (UNIT_D12_SQUARED_SEED, 30, (2, 8, 89), (0, 4, 33)),
    (SQRT2_D12_SQUARED_SEED, 40, (2, 17), (0, 5)),
    (TWO_SQRT2_D12_SQUARED_SEED, 27, (4, 8, 15), (0, 2, 5)),
)


def silver(value) -> dict[str, str]:
    return Silver.coerce(value).to_json()


def vector(values) -> list[dict[str, str]]:
    return [silver(value) for value in values]


def weight_payload(weight) -> dict:
    return {
        "square": silver(weight.square),
        "weights": vector(weight.weights),
        "coordinates": vector(weight.coordinates),
        "residue": [[*bits] for bits in weight.residue],
    }


def root_payload(root) -> dict:
    return {
        "norm": silver(root.norm),
        "height": silver(root.height),
        "control": silver(root.control),
        "target": silver(root.target),
        "left": weight_payload(root.left),
        "right": weight_payload(root.right),
        "glue_nontrivial": root.glue_nontrivial,
    }


def seed_payload(seed, height_bound, expected_nonempty, expected_glue) -> dict:
    prefix = run_d12_squared_seed_prefix(seed, height_bound)
    nonempty = tuple(packet for packet in prefix.packets if packet.roots)
    counts = tuple(len(packet.roots) for packet in nonempty)
    glue_counts = tuple(
        sum(root.glue_nontrivial for root in packet.roots) for packet in nonempty
    )
    if counts != expected_nonempty or glue_counts != expected_glue:
        raise ArithmeticError(
            f"the frozen {seed.seed_id} packets changed: {counts}, {glue_counts}"
        )
    if any(
        d12_squared_seed_pair(seed, left, right).sign(0) > 0
        for index, left in enumerate(prefix.accepted_roots)
        for right in prefix.accepted_roots[index + 1 :]
    ):
        raise ArithmeticError(f"the accepted {seed.seed_id} walls are not acute")
    if any(
        d12_squared_seed_pair(seed, root, root) != root.norm
        for root in prefix.accepted_roots
    ):
        raise ArithmeticError(f"an accepted {seed.seed_id} wall has wrong norm")
    witness = strict_hyperideal_witness(seed)
    if witness.height_bound != height_bound:
        raise ArithmeticError("the witness and frozen prefix bounds disagree")
    if not witness.hyperideal or not witness.all_positive_walls_strict:
        raise ArithmeticError(f"the {seed.seed_id} strict witness failed")
    roots = [root_payload(root) for root in prefix.accepted_roots]
    root_digest = sha256(
        json.dumps(roots, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "seed_id": seed.seed_id,
        "timelike_coefficient": silver(seed.timelike_coefficient),
        "height_bound": height_bound,
        "canonical_mirror_norms": [silver(norm) for norm in seed.canonical_mirror_norms],
        "complete_norm_class_proof": seed.complete_norm_class_proof,
        "affine_weight_cap": silver(seed.affine_weight_cap),
        "height_zero_dominant_orbits": len(height_zero_seed_dominant_roots(seed)),
        "controls": [
            {
                "norm": silver(level.norm),
                "height": silver(level.height),
                "control": silver(level.control),
                "target": silver(level.target),
            }
            for level in prefix.controls
        ],
        "packets": [
            {
                "height": silver(packet.height),
                "control_count": len(packet.controls),
                "candidate_count": packet.candidate_count,
                "rejected_count": packet.rejected_count,
                "accepted_count": len(packet.roots),
                "root_lattice_count": sum(not root.glue_nontrivial for root in packet.roots),
                "cross_glue_count": sum(root.glue_nontrivial for root in packet.roots),
            }
            for packet in prefix.packets
        ],
        "accepted_roots": roots,
        "accepted_roots_sha256": root_digest,
        "all_accepted_roots_have_declared_norm": True,
        "all_accepted_pairings_nonpositive": True,
        "strict_hyperideal_witness": {
            "block": witness.block,
            "dynkin_index_zero_based": witness.dynkin_index,
            "weight": str(witness.weight),
            "norm": silver(witness.norm),
            "height_zero_nonpositive": witness.height_zero_nonpositive,
            "all_positive_walls_strict": witness.all_positive_walls_strict,
            "hyperideal": witness.hyperideal,
            "pairing_histogram": dict(
                sorted(Counter(str(value) for value in witness.positive_wall_pairings).items())
            ),
        },
    }


def main() -> None:
    if not d12_squared_basis_audit():
        raise ArithmeticError("the D12^2 Niemeier basis audit failed")
    searches = [seed_payload(*row) for row in FROZEN_SEARCHES]
    unit_geometry = first_packet_geometry()
    payload = {
        "scope": {
            "field": "Q(sqrt(2))",
            "hyperbolic_dimension": 24,
            "spatial_lattice": (
                "Niemeier N(D12^2), exact vector-spinor graph glue "
                "{(0,0),(v,s),(s,v),(c,c)}"
            ),
            "claim": (
                "Three complete primitive chronological prefixes for the unit, "
                "sqrt(2), and 2sqrt(2) coefficients. Each has an explicit strict "
                "hyperideal ray. These are bounded screens, not global "
                "nonreflectivity theorems."
            ),
        },
        "glue": {
            "root_lattice_discriminant": 16,
            "glue_index": 4,
            "code": ["(0,0)", "(v,s)", "(s,v)", "(c,c)"],
            "diagonal_v_v_excluded_because_it_would_add_roots": True,
            "explicit_basis_even_unimodular": True,
            "height_zero_root_system": "D12^2",
            "height_zero_root_count": 528,
        },
        "unit_first_affine_product_geometry": {
            "extreme_ray_count": unit_geometry.extreme_ray_count,
            "finite_count": unit_geometry.finite_count,
            "cusp_count": unit_geometry.cusp_count,
            "hyperideal_count": unit_geometry.hyperideal_count,
        },
        "searches": searches,
        "outcome": (
            "No compact or finite-volume example in any frozen prefix. The "
            "2sqrt(2) seed is the most wall-rich early branch (27 accepted walls), "
            "but its sparse strict hyperideal witness survives through H=27; the "
            "sqrt(2) witness survives three further empty shells through H=40."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "sha256": sha256(OUTPUT.read_bytes()).hexdigest(),
                "accepted_counts": {
                    search["seed_id"]: len(search["accepted_roots"])
                    for search in searches
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
