#!/usr/bin/env python3
"""Exact coordinate and weighted-group scout for Roulleau's 59-wall rank-12 stock."""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path

import sympy as sp

from scripts.search_k3_high_rank_involution_folds import (
    _group_checks,
    _matrix_payload,
    _weighted_graph_automorphisms,
)
from scripts.search_k3_nef_subgroup_folds import _order
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_r12_ud4x2a1x2_59_scout.json"
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
SOURCE_TEX_SHA256 = "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
RANK = 12
WALL_COUNT = 59


def _vector(**coefficients: int) -> sp.Matrix:
    row = sp.zeros(1, RANK)
    for name, coefficient in coefficients.items():
        if not name.startswith("A"):
            raise ValueError(name)
        row[0, int(name[1:]) - 1] = coefficient
    return row


def _basis_gram() -> sp.Matrix:
    gram = -2 * sp.eye(RANK)
    ordinary_edges = (
        (12, 1),
        (12, 2),
        (12, 3),
        (12, 4),
        (12, 5),
        (12, 9),
        (8, 6),
        (8, 7),
        (9, 8),
        (8, 10),
        (10, 11),
    )
    for left, right in ordinary_edges:
        gram[left - 1, right - 1] = gram[right - 1, left - 1] = 1
    return gram


def coordinate_model() -> tuple[sp.Matrix, list[dict[str, object]], sp.Matrix]:
    basis_gram = _basis_gram()
    basis = [sp.eye(RANK).row(index) for index in range(RANK)]
    named: dict[str, sp.Matrix] = {f"A{index + 1}": row for index, row in enumerate(basis)}
    named["A13"] = sp.Matrix([[-1, -1, -1, -1, -1, 3, 3, 6, 2, 4, 2, -2]])
    fiber = named["A6"] + named["A7"] + 2 * named["A8"] + named["A9"] + named["A10"]
    for index in range(1, 6):
        named[f"A{13 + index}"] = fiber - named[f"A{index}"]
    named["A19"] = fiber - named["A13"]

    walls: list[dict[str, object]] = [
        {"name": f"A{index}", "kind": "A", "coordinates": named[f"A{index}"]}
        for index in range(1, 20)
    ]
    pairs = tuple(combinations(range(1, 6), 2))
    common = 6 * named["A8"] + 3 * named["A9"] + 4 * named["A10"] + 2 * named["A11"]
    for kind in ("B1", "B2", "C1", "C2"):
        for left, right in pairs:
            complement = sorted(set(range(1, 6)) - {left, right})
            if kind == "B1":
                row = (
                    -sum((named[f"A{i}"] for i in complement), sp.zeros(1, RANK))
                    + 3 * named["A6"] + 2 * named["A7"] + common
                )
            elif kind == "B2":
                row = (
                    -sum((named[f"A{i}"] for i in complement), sp.zeros(1, RANK))
                    + 2 * named["A6"] + 3 * named["A7"] + common
                )
            elif kind == "C1":
                row = (
                    -named[f"A{left}"] - named[f"A{right}"] - named["A13"]
                    + 2 * named["A6"] + 3 * named["A7"] + common
                )
            else:
                row = (
                    -named[f"A{left}"] - named[f"A{right}"] - named["A13"]
                    + 3 * named["A6"] + 2 * named["A7"] + common
                )
            walls.append(
                {
                    "name": f"{kind}_{left}{right}",
                    "kind": kind,
                    "pair": [left, right],
                    "coordinates": row,
                }
            )
    wall_rows = sp.Matrix.vstack(*(record["coordinates"] for record in walls))
    wall_gram = wall_rows * basis_gram * wall_rows.T

    d2 = sp.Matrix([[0, 0, 0, 0, 0, 1, 1, 3, 2, 2, 1, 1]])
    formula_checks = {
        "A13_printed_coordinate": named["A13"]
        == sp.Matrix([[-1, -1, -1, -1, -1, 3, 3, 6, 2, 4, 2, -2]]),
        "fiber_pair_relations": all(
            named[f"A{index}"] + named[f"A{13 + index}"] == fiber
            for index in range(1, 6)
        ) and named["A13"] + named["A19"] == fiber,
        "B_plus_C_is_2D2": all(
            walls[19 + offset]["coordinates"] + walls[39 + offset]["coordinates"] == 2 * d2
            and walls[29 + offset]["coordinates"] + walls[49 + offset]["coordinates"] == 2 * d2
            for offset in range(10)
        ),
    }
    if not all(formula_checks.values()):
        raise ArithmeticError(f"printed formula reconstruction failed: {formula_checks}")
    return basis_gram, walls, wall_gram


def build_scout() -> dict[str, object]:
    import hashlib

    if hashlib.sha256(SOURCE_TEX.read_bytes()).hexdigest() != SOURCE_TEX_SHA256:
        raise ArithmeticError("pinned Roulleau source hash changed")
    basis_gram, walls, wall_gram = coordinate_model()
    coordinate_checks = {
        "wall_count": len(walls) == WALL_COUNT,
        "basis_determinant": int(basis_gram.det()) == -64,
        "basis_inertia": exact_inertia(basis_gram) == (1, 11, 0),
        "wall_rank": wall_gram.rank() == RANK,
        "wall_inertia": exact_inertia(wall_gram) == (1, 11, 47),
        "all_roots_minus_two": all(wall_gram[index, index] == -2 for index in range(WALL_COUNT)),
        "all_coordinates_integral": all(
            value.is_Integer for record in walls for value in record["coordinates"]
        ),
    }
    if not all(coordinate_checks.values()):
        raise ArithmeticError(f"coordinate checks failed: {coordinate_checks}")
    group, search_nodes = _weighted_graph_automorphisms(wall_gram)
    group_checks = _group_checks(group)
    order_histogram = {
        str(order): count
        for order, count in sorted(Counter(_order(element) for element in group).items())
    }
    expected_histogram = {
        "1": 1,
        "2": 303,
        "3": 80,
        "4": 720,
        "5": 144,
        "6": 1200,
        "10": 432,
    }
    if not (
        all(group_checks.values())
        and len(group) == 2880
        and order_histogram == expected_histogram
    ):
        raise ArithmeticError("weighted automorphism denominator changed")
    return {
        "schema": "coxeter.roulleau_r12_ud4x2a1x2_59_scout.v1",
        "source_tex_sha256": SOURCE_TEX_SHA256,
        "source_lines": [8804, 9033],
        "lattice": "U+D4(-1)^2+A1(-1)^2",
        "basis_gram": _matrix_payload(basis_gram),
        "walls": [
            {
                **{key: value for key, value in record.items() if key != "coordinates"},
                "coordinates": _matrix_payload(record["coordinates"])[0],
            }
            for record in walls
        ],
        "wall_gram": _matrix_payload(wall_gram),
        "coordinate_checks": coordinate_checks,
        "weighted_group": {
            "order": len(group),
            "element_order_histogram": order_histogram,
            "histogram_recognition": "matches S6 x C2 x C2; GAP structure proof deferred to full certifier",
            "search_nodes": search_nodes,
            "checks": group_checks,
            "elements_zero_based": [list(element) for element in group],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_scout()
    if args.check_archive:
        if payload != json.loads(args.output.read_text()):
            raise SystemExit("recomputed scout differs from archive")
    else:
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
