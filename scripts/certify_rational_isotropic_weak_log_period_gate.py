#!/usr/bin/env python3
"""Exact replay for the rational isotropic weak-log period gate.

The producer uses only integer arithmetic.  It checks the marked
``I_(1,9)`` blowup lattice, the Riemann--Roch polynomial for ``rF``, and the
cohomology recursion dictated by a normal character of a specified torsion
order.  The geometric facts about elliptic degree-zero line bundles and the
choice of nine points on a smooth cubic are proved in the accompanying
report; a finite computation is not presented as a substitute for them.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import gcd
from pathlib import Path
from typing import Iterable, Sequence


Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]

BASIS = ("H",) + tuple(f"E_{index}" for index in range(1, 10))
GRAM: Matrix = tuple(
    tuple(
        (1 if row == column == 0 else -1 if row == column else 0)
        for column in range(10)
    )
    for row in range(10)
)
CANONICAL: Vector = (-3,) + (1,) * 9
FIBRE: Vector = (3,) + (-1,) * 9


def matrix_vector(matrix: Sequence[Sequence[int]], vector: Sequence[int]) -> Vector:
    if any(len(row) != len(vector) for row in matrix):
        raise ValueError("matrix and vector dimensions differ")
    return tuple(sum(entry * coordinate for entry, coordinate in zip(row, vector)) for row in matrix)


def dot(left: Sequence[int], right: Sequence[int], gram: Matrix = GRAM) -> int:
    if len(left) != len(right) or len(gram) != len(left):
        raise ValueError("intersection dimensions differ")
    return sum(a * b for a, b in zip(left, matrix_vector(gram, right)))


def determinant_diagonal(matrix: Matrix) -> int:
    if any(
        matrix[row][column] != 0
        for row in range(len(matrix))
        for column in range(len(matrix))
        if row != column
    ):
        raise ValueError("this exact helper accepts only diagonal matrices")
    answer = 1
    for index in range(len(matrix)):
        answer *= matrix[index][index]
    return answer


def exceptional(index: int) -> Vector:
    if not 1 <= index <= 9:
        raise ValueError("exceptional index must lie in 1,...,9")
    return tuple(1 if coordinate == index else 0 for coordinate in range(10))


def riemann_roch_euler(r: int) -> int:
    """Return chi(O_S(rF)) from the exact surface Riemann--Roch formula."""

    if r < 0:
        raise ValueError("r must be nonnegative")
    r_fibre = tuple(r * coordinate for coordinate in FIBRE)
    r_fibre_minus_k = tuple(
        coordinate - canonical
        for coordinate, canonical in zip(r_fibre, CANONICAL)
    )
    numerator = dot(r_fibre, r_fibre_minus_k)
    if numerator % 2:
        raise AssertionError("Riemann--Roch numerator is not even")
    return 1 + numerator // 2


def torsion_recurrence(order: int, horizon: int) -> list[dict[str, object]]:
    """Replay the exact h^0/h^1 recurrence for a character of exact order m."""

    if order <= 0:
        raise ValueError("torsion order must be positive")
    if horizon < 0:
        raise ValueError("horizon must be nonnegative")
    rows: list[dict[str, object]] = []
    h0 = 1
    for r in range(horizon + 1):
        if r:
            normal_power_trivial = r % order == 0
            h0 += int(normal_power_trivial)
        else:
            normal_power_trivial = True
        expected_h0 = r // order + 1
        if h0 != expected_h0 or riemann_roch_euler(r) != 1:
            raise AssertionError("torsion recurrence failed")
        h2 = 0
        h1 = h0 - 1
        rows.append(
            {
                "r": r,
                "normal_power_trivial": normal_power_trivial,
                "elliptic_restriction_h0": (
                    None if r == 0 else int(normal_power_trivial)
                ),
                "restriction_map_rank": (
                    None if r == 0 else int(normal_power_trivial)
                ),
                "h0": h0,
                "h1": h1,
                "h2": h2,
                "chi": h0 - h1 + h2,
                "closed_formula_h0": expected_h0,
            }
        )
    return rows


def nontorsion_recurrence(horizon: int) -> list[dict[str, object]]:
    """Replay the recurrence when every positive normal power is nontrivial."""

    if horizon < 0:
        raise ValueError("horizon must be nonnegative")
    rows: list[dict[str, object]] = []
    for r in range(horizon + 1):
        chi = riemann_roch_euler(r)
        row = {
            "r": r,
            "normal_power_trivial": r == 0,
            "elliptic_restriction_h0": None if r == 0 else 0,
            "restriction_map_rank": None if r == 0 else 0,
            "h0": 1,
            "h1": 0,
            "h2": 0,
            "chi": chi,
        }
        if chi != 1:
            raise AssertionError("non-torsion recurrence failed")
        rows.append(row)
    return rows


def lattice_record() -> dict[str, object]:
    exceptional_pairings = {
        f"F.E_{index}": dot(FIBRE, exceptional(index))
        for index in range(1, 10)
    }
    record: dict[str, object] = {
        "basis": list(BASIS),
        "gram": [list(row) for row in GRAM],
        "signature": [1, 9],
        "determinant": determinant_diagonal(GRAM),
        "canonical": list(CANONICAL),
        "F": list(FIBRE),
        "F_primitive": gcd(*map(abs, FIBRE)) == 1,
        "F_equals_minus_K": all(
            fibre == -canonical
            for fibre, canonical in zip(FIBRE, CANONICAL)
        ),
        "K_square": dot(CANONICAL, CANONICAL),
        "F_square": dot(FIBRE, FIBRE),
        "K_dot_F": dot(CANONICAL, FIBRE),
        "adjunction_two_g_minus_two": dot(
            FIBRE,
            tuple(fibre + canonical for fibre, canonical in zip(FIBRE, CANONICAL)),
        ),
        "normal_degree": dot(FIBRE, FIBRE),
        "exceptional_pairings": exceptional_pairings,
    }
    if record["determinant"] != -1:
        raise AssertionError("I_(1,9) determinant failed")
    if not record["F_equals_minus_K"]:
        raise AssertionError("F=-K failed")
    if not record["F_primitive"]:
        raise AssertionError("F is not primitive")
    for key in ("K_square", "F_square", "K_dot_F", "adjunction_two_g_minus_two", "normal_degree"):
        if record[key] != 0:
            raise AssertionError(f"{key} failed")
    if set(exceptional_pairings.values()) != {1}:
        raise AssertionError("F.E_i=1 failed")
    return record


def normalize_orders(orders: Iterable[int]) -> tuple[int, ...]:
    normalized = tuple(dict.fromkeys(int(order) for order in orders))
    if not normalized or any(order <= 0 for order in normalized):
        raise ValueError("at least one positive torsion order is required")
    return normalized


def build_certificate(
    horizon: int = 64,
    orders: Iterable[int] = (1, 2, 3, 4, 7, 12),
) -> dict[str, object]:
    if horizon < 1:
        raise ValueError("horizon must be at least one")
    normalized_orders = normalize_orders(orders)
    return {
        "schema": "coxeter.rational_isotropic_weak_log_period_gate.v1",
        "theorem": {
            "ground_field": "C",
            "surface": "smooth projective rational surface",
            "divisor": (
                "smooth connected effective genus-one Cartier divisor F, "
                "primitive in Pic(S)=NS(S)"
            ),
            "numerical_hypotheses": ["F^2=0", "K.F=0", "F nef"],
            "normal_character": "N=O_F(F) in Pic^0(F)",
            "equivalence": [
                "F semiample",
                "some |rF| moves",
                "N torsion",
            ],
            "exact_torsion_formula": "h0(rF)=floor(r/m)+1; h1(rF)=floor(r/m)",
            "non_torsion_formula": "h0(rF)=1; h1(rF)=0 for every r>=0",
        },
        "blowup_lattice": lattice_record(),
        "riemann_roch": {
            "formula": "chi(O(rF))=1+rF.(rF-K)/2",
            "values_through_horizon": [
                riemann_roch_euler(r) for r in range(horizon + 1)
            ],
            "H2_vanishing_argument": (
                "effectivity of K-rF together with rF would make K effective, "
                "contradicting p_g=0"
            ),
        },
        "torsion_profiles": {
            str(order): torsion_recurrence(order, horizon)
            for order in normalized_orders
        },
        "nontorsion_profile": nontorsion_recurrence(horizon),
        "same_marked_lattice_examples": {
            "common": {
                "surface": "Blow up nine distinct points of a smooth plane cubic C",
                "NS": "I_(1,9)",
                "K": "-3H+sum(E_i)",
                "F": "3H-sum(E_i)=-K",
                "normal_character": "O_C(3)(-p_1-...-p_9)",
                "F_nef": True,
                "nef_proof": "Bezout plus F.E_i=1 and F^2=0",
            },
            "trivial": {
                "points": "transverse base locus of a cubic pencil containing C",
                "normal_order": 1,
                "first_moving_multiple": 1,
                "F_semiample": True,
            },
            "nontorsion": {
                "points": "choose p_9 outside the countable torsion preimage",
                "normal_order": "infinite",
                "first_moving_multiple": None,
                "F_semiample": False,
            },
        },
        "scope": {
            "smooth_genus_one_representative_assumed": True,
            "every_isotropic_weak_log_class_covered": False,
            "same_marked_NS_K_F_and_nefness_proved": True,
            "same_complete_nef_cone_claimed": False,
            "prescribed_Coxeter_polytope_realized": False,
            "complete_nef_chamber_returned": False,
            "all_physical_facets_returned": False,
            "finite_hyperbolic_volume_claimed": False,
            "reflection_saturation_claimed": False,
            "reflection_maximality_claimed": False,
            "raw_polytope_exhaustivity_claimed": False,
            "computer_replay_replaces_geometric_proof": False,
        },
        "replay_horizon": horizon,
        "replayed_torsion_orders": list(normalized_orders),
    }


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def envelope(
    horizon: int = 64,
    orders: Iterable[int] = (1, 2, 3, 4, 7, 12),
) -> dict[str, object]:
    payload = build_certificate(horizon=horizon, orders=orders)
    return {"payload": payload, "payload_sha256": digest(payload)}


def encoded_certificate(
    horizon: int = 64,
    orders: Iterable[int] = (1, 2, 3, 4, 7, 12),
) -> str:
    return json.dumps(
        envelope(horizon=horizon, orders=orders), indent=2, sort_keys=True
    ) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--horizon", type=int, default=64)
    parser.add_argument("--orders", type=int, nargs="*", default=[1, 2, 3, 4, 7, 12])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    arguments = parser.parse_args()
    encoded = encoded_certificate(
        horizon=arguments.horizon,
        orders=arguments.orders,
    )
    if arguments.check is not None:
        if arguments.check.read_text(encoding="utf-8") != encoded:
            raise SystemExit(f"archive mismatch: {arguments.check}")
        print(f"PASS archive replay: {arguments.check}")
        return
    if arguments.output is not None:
        arguments.output.write_text(encoded, encoding="utf-8")
        print(arguments.output)
        return
    print(encoded, end="")


if __name__ == "__main__":
    main()
