#!/usr/bin/env python3
"""Exact arithmetic replay for the wrong-log orbifold rank gate."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/rational_wrong_log_orbifold_rank_gate.json"

PRIMARY_SOURCE_PINS = {
    "arxiv_id": "2607.18699v1",
    "theorem": "Iwai-Jinnouchi-Zhang Theorem 1.3",
    "pdf_sha256": "2739f4489b0a4cf0722141ea1954f871fae65269c8aec63c02f66bfad398c501",
    "source_tar_sha256": "a18775cdb3f99507da702375c03e983760c894ab9ddd36eaebc9ee96910e8442",
    "tex_sha256": "6cb03d24e199e4524ed0e2cec4e7bce6371c2aee09f910bde9deb4a297877517",
}


def canonical_degree(self_square: int) -> int:
    """K.C for a smooth rational curve."""
    return -2 - self_square


def natural_wrong_log_pairing(
    self_square: int,
    boundary_degree: int,
    other_boundary_degree: int = 0,
) -> Fraction:
    """D.C for D=-(K+B/2).

    For a boundary component, boundary_degree is its own -4 contribution and
    other_boundary_degree is the sum of intersections with the other B_j.
    """
    if self_square == -4:
        if boundary_degree != -4:
            raise ValueError("a boundary component must include its own -4")
        return -Fraction(
            canonical_degree(self_square)
            + Fraction(boundary_degree + other_boundary_degree, 2)
        )
    return -Fraction(
        canonical_degree(self_square) + Fraction(boundary_degree, 2)
    )


def d_square(picard_rank: int, boundary_count: int) -> int:
    return 10 - picard_rank + boundary_count


def orbifold_c2_times_four(picard_rank: int, boundary_count: int) -> int:
    """Four times c2_orb after contracting b disjoint -4 curves."""
    return 4 * (picard_rank + 2) - 7 * boundary_count


def exact_picard_ceiling() -> int:
    """Eliminate b symbolically, without an enumeration horizon.

    D^2>=0 gives b>=rho-10, while c2_orb>=0 gives
    7b<=4rho+8.  Hence 7(rho-10)<=4rho+8.
    """
    return (7 * 10 + 8) // (7 - 4)


def admissible_rank_pairs(limit: int = 200) -> list[tuple[int, int]]:
    """Bounded calibration of integer pairs allowed by both inequalities."""
    return [
        (rho, b)
        for rho in range(3, limit + 1)
        for b in range(1, limit + 1)
        if d_square(rho, b) >= 0
        and orbifold_c2_times_four(rho, b) >= 0
    ]


def _bilinear(
    left: tuple[Fraction, ...],
    gram: tuple[tuple[int, ...], ...],
    right: tuple[Fraction, ...],
) -> Fraction:
    return sum(
        left[i] * gram[i][j] * right[j]
        for i in range(len(left))
        for j in range(len(right))
    )


def f4_control() -> dict[str, Any]:
    gram = (
        (-4, 1, 0),
        (1, 0, 0),
        (0, 0, -1),
    )
    B = (Fraction(1), Fraction(0), Fraction(0))
    E = (Fraction(0), Fraction(0), Fraction(1))
    G = (Fraction(0), Fraction(1), Fraction(-1))
    K = (Fraction(-2), Fraction(-6), Fraction(1))
    log_class = tuple(K[i] + B[i] / 2 for i in range(3))
    pairings = {
        "B": _bilinear(log_class, gram, B),
        "E": _bilinear(log_class, gram, E),
        "G": _bilinear(log_class, gram, G),
    }
    return {
        "rho": 3,
        "boundary_count": 1,
        "K_plus_half_B_pairings": {
            key: str(value) for key, value in pairings.items()
        },
        "D_square": str(_bilinear(log_class, gram, log_class)),
        "orbifold_c2": str(
            Fraction(orbifold_c2_times_four(3, 1), 4)
        ),
        "passes": pairings
        == {"B": Fraction(0), "E": Fraction(-1), "G": Fraction(-1, 2)}
        and _bilinear(log_class, gram, log_class) == 8,
    }


def ruled_curve(
    bidegree: tuple[int, int],
    private_blowups: int,
) -> dict[str, int | list[int]]:
    a, b = bidegree
    return {
        "bidegree": [a, b],
        "private_blowups": private_blowups,
        "self_square": 2 * a * b - private_blowups,
        "canonical_degree": -2 * (a + b) + private_blowups,
    }


def base_intersection(
    left: tuple[int, int], right: tuple[int, int]
) -> int:
    return left[0] * right[1] + left[1] * right[0]


def hostile_geometric_controls() -> list[dict[str, Any]]:
    rows = [
        {
            "id": "H1_boundary_boundary",
            "left": ruled_curve((1, 1), 6),
            "right": ruled_curve((1, 1), 6),
            "intersection": base_intersection((1, 1), (1, 1)),
            "expected": [-4, -4, 2, 2, 2],
        },
        {
            "id": "H2_boundary_root",
            "left": ruled_curve((1, 1), 6),
            "right": ruled_curve((1, 1), 4),
            "intersection": base_intersection((1, 1), (1, 1)),
            "expected": [-4, -2, 2, 2, 0],
        },
        {
            "id": "H3_boundary_exceptional",
            "left": ruled_curve((1, 1), 6),
            "right": ruled_curve((1, 2), 5),
            "intersection": base_intersection((1, 1), (1, 2)),
            "expected": [-4, -1, 3, 2, -1],
        },
    ]
    for row in rows:
        actual = [
            row["left"]["self_square"],
            row["right"]["self_square"],
            row["intersection"],
            row["left"]["canonical_degree"],
            row["right"]["canonical_degree"],
        ]
        row["actual"] = actual
        row["passes"] = actual == row["expected"]
    boundary_fiber = ruled_curve((1, 0), 4)
    exceptional_fiber = ruled_curve((0, 1), 1)
    component_degrees = [base_intersection((1, 0), (0, 1))] * 3
    total_degree = sum(component_degrees)
    square_sum = sum(value * value for value in component_degrees)
    coefficient_lower_bound = Fraction(
        Fraction(total_degree, 2) - 1,
        1 - Fraction(square_sum, 4),
    )
    rows.append(
        {
            "id": "H3_distributed_boundary_exceptional",
            "boundary_components": [boundary_fiber.copy() for _ in range(3)],
            "exceptional": exceptional_fiber,
            "boundary_pairwise_intersection": base_intersection(
                (1, 0), (1, 0)
            ),
            "component_boundary_degrees": component_degrees,
            "aggregate_boundary_degree": total_degree,
            "sum_squared_degrees": square_sum,
            "forced_exceptional_coefficient_lower_bound": str(
                coefficient_lower_bound
            ),
            "passes": boundary_fiber["self_square"] == -4
            and boundary_fiber["canonical_degree"] == 2
            and exceptional_fiber["self_square"] == -1
            and exceptional_fiber["canonical_degree"] == -1
            and total_degree == 3
            and square_sum == 3
            and coefficient_lower_bound == 2,
        }
    )
    return rows


def build_certificate() -> dict[str, Any]:
    maximum_rank = exact_picard_ceiling()
    pairs = admissible_rank_pairs(limit=maximum_rank)
    extremal_pairs = [[rho, b] for rho, b in pairs if rho == maximum_rank]
    pairing_table = {
        "boundary_disjoint": str(natural_wrong_log_pairing(-4, -4, 0)),
        "boundary_meets_boundary_twice": str(
            natural_wrong_log_pairing(-4, -4, 2)
        ),
        "root_boundary_degree_0": str(natural_wrong_log_pairing(-2, 0)),
        "root_boundary_degree_2": str(natural_wrong_log_pairing(-2, 2)),
        "exceptional_boundary_degree_0": str(
            natural_wrong_log_pairing(-1, 0)
        ),
        "exceptional_boundary_degree_1": str(
            natural_wrong_log_pairing(-1, 1)
        ),
        "exceptional_boundary_degree_2": str(
            natural_wrong_log_pairing(-1, 2)
        ),
        "exceptional_boundary_degree_3": str(
            natural_wrong_log_pairing(-1, 3)
        ),
    }
    controls = hostile_geometric_controls()
    payload = {
        "schema_version": 2,
        "primary_source_pins": PRIMARY_SOURCE_PINS,
        "elementary_pairing_table_D_equals_minus_K_plus_half_B": pairing_table,
        "rank_sieve": {
            "conditions": [
                "10-rho+b >= 0",
                "4*(rho+2)-7*b >= 0",
            ],
            "maximum_picard_rank": maximum_rank,
            "maximum_hyperbolic_dimension": maximum_rank - 1,
            "extremal_integer_pairs": extremal_pairs,
            "symbolic_elimination": (
                "7*(rho-10)<=4*rho+8, hence 3*rho<=78"
            ),
        },
        "supported_boundary_obstructions": {
            "H1_at_intersection_2": "a_i+a_j>=2 contradicts a_i,a_j<1",
            "H2_at_intersection_2": "a_B>=1 contradicts a_B<1",
            "H3_aggregate_degree_at_least_3": (
                "(1-Q/4)*a_E>=T/2-1; T>=3 forces a_E>=2 or impossibility"
            ),
        },
        "F4_control": f4_control(),
        "geometric_hostile_controls": controls,
    }
    payload["all_checks_pass"] = (
        pairing_table
        == {
            "boundary_disjoint": "0",
            "boundary_meets_boundary_twice": "-1",
            "root_boundary_degree_0": "0",
            "root_boundary_degree_2": "-1",
            "exceptional_boundary_degree_0": "1",
            "exceptional_boundary_degree_1": "1/2",
            "exceptional_boundary_degree_2": "0",
            "exceptional_boundary_degree_3": "-1/2",
        }
        and maximum_rank == 26
        and extremal_pairs == [[26, 16]]
        and payload["F4_control"]["passes"]
        and all(row["passes"] for row in controls)
    )
    return payload


def certify(payload: dict[str, Any]) -> None:
    if not payload["all_checks_pass"]:
        raise AssertionError("wrong-log orbifold certificate failed")
    if payload["rank_sieve"]["maximum_picard_rank"] != 26:
        raise AssertionError("unexpected Picard ceiling")
    if payload["rank_sieve"]["maximum_hyperbolic_dimension"] != 25:
        raise AssertionError("unexpected dimension ceiling")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()

    payload = build_certificate()
    certify(payload)
    if args.check_archive:
        archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
        if archived != payload:
            raise AssertionError(f"archive mismatch: {ARCHIVE}")
    if args.print_json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print("wrong-log orbifold rank gate: PASS (rho<=26, n<=25)")


if __name__ == "__main__":
    main()
