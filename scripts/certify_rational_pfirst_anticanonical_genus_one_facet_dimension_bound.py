#!/usr/bin/env python3
"""Certify the sharp anticanonical genus-one facet rank theorem.

The theorem itself is an elementary unimodularity/adjunction argument in the
companion report.  This script pins and rechecks the two exact endpoint
carriers X_10 and X_11 from the independent rational reservoir archive.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/rational_pfirst_anticanonical_genus_one_facet_dimension_bound.json"
RESERVOIR = ROOT / "data/enumerations/other_algebraic_reservoirs_certificate.json"
RESERVOIR_SHA256 = "bb02c9fdf14d3370e97ffd2bb89a3a7801f8bd578fe73972810b2ddc4e58882e"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dot(x: Sequence[int], y: Sequence[int]) -> int:
    return x[0] * y[0] - sum(a * b for a, b in zip(x[1:], y[1:]))


def reflection_integral(root: Sequence[int]) -> bool:
    square = dot(root, root)
    for j in range(len(root)):
        basis = [0] * len(root)
        basis[j] = 1
        coefficient = Fraction(-2 * dot(basis, root), square)
        if any((Fraction(basis[i]) + coefficient * root[i]).denominator != 1 for i in range(len(root))):
            return False
    return True


def endpoint(source: dict[str, object], n: int) -> dict[str, object]:
    case = source[f"p{n}_rational_surface"]
    labels = case["simple_curve_labels"]
    vectors = [list(v) for v in zip(*case["simple_curve_columns"])]
    f = vectors[-1]
    expected_f = [3] + [-1] * n
    k = [-z for z in expected_f]
    if f != expected_f or labels[-1] != f"3H-E1-...-E{n}=-K_X":
        raise ArithmeticError(f"P{n} anticanonical facet changed")
    square = dot(f, f)
    canonical_pairing = dot(k, f)
    arithmetic_genus = 1 + (square + canonical_pairing) // 2
    if square != 9 - n or canonical_pairing != n - 9 or arithmetic_genus != 1:
        raise ArithmeticError(f"P{n} adjunction data changed")
    if not reflection_integral(f):
        raise ArithmeticError(f"P{n} anticanonical reflection ceased to be integral")
    cone = case["exact_cone_extreme_rays"]
    if cone["hyperideal"] != 0 or cone["ideal"] != 2:
        raise ArithmeticError(f"P{n} volume certificate changed")
    if not case["all_checks_pass"] or not all(case["checks"].values()):
        raise ArithmeticError(f"P{n} carrier certificate failed")
    return {
        "n": n,
        "rho": n + 1,
        "lattice": case["lattice"],
        "surface": case["surface"],
        "known_chamber": case["known_chamber"],
        "facet_count": len(labels),
        "anticanonical_facet": f,
        "F_square": square,
        "K_dot_F": canonical_pairing,
        "arithmetic_genus": arithmetic_genus,
        "smooth_genus_one_in_source_construction": True,
        "reflection_integral": True,
        "complete_nef_cone_certified": True,
        "finite_vertices": case["known_exact_face_data"]["finite_vertices"],
        "ideal_vertices": case["known_exact_face_data"]["cusps"],
        "hyperideal_rays": cone["hyperideal"],
    }


def build() -> dict[str, object]:
    if sha256(RESERVOIR) != RESERVOIR_SHA256:
        raise ArithmeticError("frozen rational reservoir archive changed")
    source = json.loads(RESERVOIR.read_text())
    endpoints = [endpoint(source, 10), endpoint(source, 11)]
    if [row["F_square"] for row in endpoints] != [-1, -2]:
        raise ArithmeticError("the two allowed norms are not both attained")
    return {
        "schema": "rational-pfirst-anticanonical-genus-one-facet-dimension-bound-v1",
        "status": "PASS_SHARP_ANTICANONICAL_FACET_BOUND_N_LE_11",
        "all_checks_pass": True,
        "theorem": {
            "hypotheses": [
                "S smooth projective rational with rho(S)>=3",
                "F=-K_S is a reduced irreducible negative curve defining a genuine full-nef facet",
                "reflection in F-perp preserves NS(S)",
            ],
            "deduction": [
                "NS(S) is unimodular and K_S is primitive",
                "F^2 is a negative divisor of 2",
                "F^2=K_S^2 is -1 or -2",
                "rho(S)=11 or 12",
                "hyperbolic nef dimension n=rho-1 is 10 or 11",
            ],
            "bound": "n<=11",
            "sharp": True,
            "classification_inputs": [],
        },
        "exact_returns": endpoints,
        "source_archive": str(RESERVOIR.relative_to(ROOT)),
        "source_archive_sha256": RESERVOIR_SHA256,
        "scope": {
            "facet_is_literal_curve_class_not_only_same_hyperplane": True,
            "genus_means_arithmetic_genus_one; endpoint cubics are smooth": True,
            "not_claimed": [
                "every rational full-nef chamber has an anticanonical facet",
                "a bound for facets proportional to -K with nonprimitive curve class",
                "a bound for adapted (-4) boundary facets B=2b",
                "classification of all chambers in dimensions 10 and 11",
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    payload = build()
    if args.check_archive:
        if payload != json.loads(ARCHIVE.read_text()):
            raise SystemExit("live anticanonical-facet certificate differs from archive")
        print("archived anticanonical-facet certificate matches")
    elif args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(payload["status"])
        for endpoint_row in payload["exact_returns"]:
            print(endpoint_row["known_chamber"], endpoint_row["F_square"])


if __name__ == "__main__":
    main()
