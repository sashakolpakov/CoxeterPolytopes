#!/usr/bin/env python3
"""Print the finite rank-five/rank-six AlVin RC7 seed manifest."""

from __future__ import annotations

import argparse
import json

from coxeter_polytopes.search.cubic7_seed_manifest import (
    initial_cubic7_diagonal_seeds,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--indent", type=int, default=2)
    parser.add_argument("--binary", default="./alvin")
    parser.add_argument("--without-control", action="store_true")
    arguments = parser.parse_args()

    seeds = initial_cubic7_diagonal_seeds(
        include_control=not arguments.without_control
    )
    record = {
        "schema": "coxeter.cubic7.diagonal_seed_manifest.v1",
        "field": "Q(2*cos(2*pi/7))",
        "power_basis": ["1", "theta", "theta^2"],
        "cutoff": {
            "kind": "squarefree ideal-radical norm",
            "maximum": 128,
            "atom_palette": "(2), ramified p7, and the three primes above 13",
            "included_support": "products of distinct atoms within the cutoff",
            "deferred_support": "inert primes >=3 and split primes >=29",
        },
        "equivalence_reduction": (
            "permutations of equal unit coordinates, global similarity, "
            "and multiplication of t by unit squares"
        ),
        "ranks": [5, 6],
        "seeds": [
            {
                "id": seed.seed_id,
                "coefficient_power_basis": seed.coefficient.to_json(),
                "coefficient_rc7_basis": list(seed.rc7_coefficient),
                "field_norm": seed.field_norm,
                "radical_norm": seed.radical_norm,
                "support": seed.support,
                "priority": seed.priority,
                "estimated_cost": seed.estimated_cost,
                "rationale": seed.rationale,
                "calibration": seed.calibration,
                "jobs": [
                    {
                        "rank": rank,
                        "qf": seed.alvin_qf(rank),
                        "command": seed.alvin_command(rank, arguments.binary),
                    }
                    for rank in (5, 6)
                ],
            }
            for seed in seeds
        ],
    }
    print(json.dumps(record, indent=arguments.indent, sort_keys=True))


if __name__ == "__main__":
    main()
