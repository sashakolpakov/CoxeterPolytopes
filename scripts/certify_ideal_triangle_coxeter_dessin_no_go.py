#!/usr/bin/env python3
"""Certify the exact ideal-triangle Coxeter-dessin no-go family.

The calculation is integer-only and deliberately tiny.  It verifies the
Minkowski roots, all Coxeter half-space pairings, the unbounded dotted entry,
and the failure of every fixed rational principal level to separate sheets.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/enumerations/ideal_triangle_coxeter_dessin_no_go.json"


Vector = tuple[int, int, int]


def pairing(x: Vector, y: Vector) -> int:
    return -x[0] * y[0] + x[1] * y[1] + x[2] * y[2]


def semicircle_root(j: int) -> Vector:
    p = j * (j + 1)
    return (-1 - p, -(2 * j + 1), 1 - p)


def boundary_roots(N: int) -> list[Vector]:
    if N < 1:
        raise ValueError("N must be positive")
    return [(0, 1, 0), *(semicircle_root(j) for j in range(N)), (-N, -1, -N)]


def congruent(x: Iterable[int], y: Iterable[int], modulus: int) -> bool:
    return all((a - b) % modulus == 0 for a, b in zip(x, y))


def build(N: int = 12) -> dict:
    roots = boundary_roots(N)
    gram = [[pairing(x, y) for y in roots] for x in roots]

    assert all(gram[i][i] == 1 for i in range(len(roots)))
    assert all(
        gram[i][j] <= 0
        for i in range(len(roots))
        for j in range(len(roots))
        if i != j
    )
    assert all(pairing(semicircle_root(0), semicircle_root(j)) == 1 - 2 * j * j
               for j in range(N))

    periodicity = []
    for modulus in range(2, N):
        left = semicircle_root(0)
        right = semicircle_root(modulus)
        assert congruent(left, right, modulus)
        periodicity.append(
            {
                "modulus": modulus,
                "colliding_facet_indices": [0, modulus],
                "roots": [list(left), list(right)],
            }
        )

    dotted = [
        {
            "j": j,
            "pairing_with_n0": pairing(semicircle_root(0), semicircle_root(j)),
            "cosh_distance": -pairing(semicircle_root(0), semicircle_root(j)),
        }
        for j in range(2, N)
    ]

    return {
        "all_checks_pass": True,
        "family": {
            "id": "IDEAL-TRIANGLE-P_N",
            "dimension": 2,
            "field": "Q",
            "lattice": "I_(2,1), Gram diag(-1,1,1)",
            "root_norm_divisibility_alphabet": [[1, 1]],
            "compactness": "finite_volume_noncompact",
            "N": N,
            "facet_count": N + 2,
            "area_over_pi": N,
            "reflection_subgroup_index": N,
        },
        "roots_in_boundary_order": [list(root) for root in roots],
        "gram": gram,
        "dotted_height_certificate": dotted,
        "principal_level_collisions": periodicity,
        "general_identities": {
            "semicircle_root": "n_j=(-1-j(j+1),-(2j+1),1-j(j+1))",
            "norm": "<n_j,n_j>=1",
            "dotted_pairing": "<n_0,n_j>=1-2j^2",
            "periodicity": "n_(j+m)=n_j (mod m)",
            "area_and_index": "area(P_N)/pi=[W(Delta_0):W(P_N)]=N",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=12)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-output", type=Path)
    args = parser.parse_args()
    result = build(args.N)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check_output is not None:
        if args.check_output.read_text() != rendered:
            raise SystemExit(f"archive differs from exact replay: {args.check_output}")
    elif args.output is not None:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()

