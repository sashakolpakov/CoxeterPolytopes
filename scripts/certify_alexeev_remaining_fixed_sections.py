#!/usr/bin/env python3
"""Certify all high fixed sections in three remaining Alexeev stocks.

The input consists of the complete primitive-root Gram matrices for the
reflective 2-elementary lattices ``(18,4,0)``, ``(15,11,1)``, and
``(15,7,1)`` in Alexeev's arXiv ancillary file.  For each stock this script

* reconstructs the full automorphism group of the *weighted* wall Gram;
* asks GAP for every subgroup conjugacy class and freezes representatives;
* proves that the list of proper fixed spaces above H9 is complete; and
* constructs every resulting strict, integral Coxeter section exactly.

This is the construction certificate.  Finite volume is deliberately left
to the companion literal Vinberg edge/cusp certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
import scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries as inventory_tools
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    SOURCE,
    SOURCE_SHA256,
    _ambient_reconstruction,
)
from scripts.screen_alexeev_20_2_1_15_9_1_fixed_sections import _screen_case


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_remaining_high_fixed_section_construction.json"
)
SCHEMA = "coxeter.alexeev_remaining_high_fixed_sections.construction.v1"


CONFIGS: dict[tuple[int, int, int], dict[str, Any]] = {
    (18, 4, 0): {
        "wall_count": 48,
        "rank": 18,
        "determinant": -16,
        "diagram_group_name": "S2 wreath (S4 x S4)",
        "diagram_group_order": 1152,
        "geometric_context": (
            "Alexeev's reflective (18,4,0) K3 Picard lattice; the full "
            "weighted diagram group is used as the search denominator"
        ),
    },
    (15, 11, 1): {
        "wall_count": 66,
        "rank": 15,
        "determinant": 2048,
        "diagram_group_name": "S2 wreath (S3 x S3)",
        "diagram_group_order": 72,
        "geometric_context": (
            "Alexeev's reflective (15,11,1) lattice; the full weighted "
            "diagram group is used as the search denominator"
        ),
    },
    (15, 7, 1): {
        "wall_count": 60,
        "rank": 15,
        "determinant": 128,
        "diagram_group_name": "S5",
        "diagram_group_order": 120,
        "geometric_context": (
            "Alexeev's reflective (15,7,1) K3 Picard lattice; the full "
            "weighted diagram group is used as the search denominator"
        ),
    },
}


EXPECTED_HIGH_CLASSES = {
    (18, 4, 0): [
        # class, order, conjugacy-class size, fixed rank, dimension, cyclic
        (3, 2, 36, 12, 11, True),
        (4, 2, 24, 11, 10, True),
        (6, 2, 12, 14, 13, True),
        (18, 4, 36, 11, 10, False),
    ],
    (15, 11, 1): [(2, 2, 6, 11, 10, True)],
    (15, 7, 1): [(1, 2, 10, 11, 10, True)],
}


EXPECTED_SECTIONS = {
    "A18_4_0_C03": {
        "dimension": 11,
        "orbits": 30,
        "facets": 28,
        "redundant": 2,
        "determinant": -256,
        "smith": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4],
        "pair_counts": {"2": 264, "3": 16, "4": 12, "infinity": 46, "ultraparallel": 40},
        "facet_gram_sha256": "1324606e59c1f994ab871fc6d354d22b405d96a8929e2d4d8e71b5ac4df469c1",
    },
    "A18_4_0_C04": {
        "dimension": 10,
        "orbits": 31,
        "facets": 28,
        "redundant": 3,
        "determinant": 1024,
        "smith": [1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 4],
        "pair_counts": {"2": 251, "3": 24, "4": 20, "infinity": 61, "ultraparallel": 22},
        "facet_gram_sha256": "5fd6cf9a5380ca9c0c34495f4d1009b9b713006ae10879daad6f3e1b1e0e1bd5",
    },
    "A18_4_0_C06": {
        "dimension": 13,
        "orbits": 31,
        "facets": 31,
        "redundant": 0,
        "determinant": -64,
        "smith": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4, 4],
        "pair_counts": {"2": 327, "3": 44, "4": 28, "infinity": 30, "ultraparallel": 36},
        "facet_gram_sha256": "d5d459ccaa4cebdbd9dfa42ca16c66446c39537eda213055bdc4fb91603dc379",
    },
    "A18_4_0_C18": {
        "dimension": 10,
        "orbits": 22,
        "facets": 21,
        "redundant": 1,
        "determinant": 128,
        "smith": [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 8],
        "pair_counts": {"2": 155, "3": 12, "4": 6, "infinity": 15, "ultraparallel": 22},
        "facet_gram_sha256": "ca4aef07f8579548edd87326d8488993c889cbf3c7cdb32909f4fff542ed511e",
    },
    "A15_11_1_C02": {
        "dimension": 10,
        "orbits": 40,
        "facets": 21,
        "redundant": 19,
        "determinant": 8192,
        "smith": [2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4],
        "pair_counts": {"2": 161, "3": 12, "4": 21, "infinity": 13, "ultraparallel": 3},
        "facet_gram_sha256": "467a01f91899c1ebf5bab5d6aaa4d5a240bf51323a23413f32d1034f1f187103",
    },
    "A15_7_1_C01": {
        "dimension": 10,
        "orbits": 40,
        "facets": 27,
        "redundant": 13,
        "determinant": 512,
        "smith": [1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4],
        "pair_counts": {"2": 235, "3": 15, "4": 14, "infinity": 48, "ultraparallel": 39},
        "facet_gram_sha256": "2b925864afa0e0a586dfc0a2d819ab491b8dc198c3b5b961f44435a6f23c5973",
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _high_signature(example: dict[str, Any]) -> list[tuple[Any, ...]]:
    return [
        (
            int(row["class_index"]),
            int(row["order"]),
            int(row["conjugacy_class_size"]),
            int(row["fixed_rank"]),
            int(row["hyperbolic_dimension"]),
            bool(row["cyclic"]),
        )
        for row in example["subgroup_classes"]
        if row["proper_fixed_section_above_h9"]
    ]


def _all_wall_sum_interior(wall_gram: sp.Matrix, rank: int) -> sp.Matrix:
    ambient_gram, wall_rows = _ambient_reconstruction(wall_gram, rank)
    interior = sum(
        (wall_rows.row(index) for index in range(wall_rows.rows)),
        sp.zeros(1, rank),
    )
    square = (interior * ambient_gram * interior.T)[0]
    minimum = min(interior * ambient_gram * wall_rows.T)
    if square != 7360 or minimum != 16:
        raise ArithmeticError("the exact (15,7,1) all-wall interior changed")
    return interior


def _audit_section(case: dict[str, Any]) -> dict[str, bool]:
    expected = EXPECTED_SECTIONS[case["id"]]
    lattice = case["saturated_fixed_lattice"]
    reduction = case["exact_reduction"]
    checks = {
        "dimension": case["case"]["hyperbolic_dimension"] == expected["dimension"],
        "wall_orbits": case["case"]["wall_orbit_count"] == expected["orbits"],
        "strict_facets": reduction["active_hyperbolic_facet_count"] == expected["facets"],
        "exact_redundancies": reduction["conic_redundant_count"] == expected["redundant"],
        "no_unresolved_facets": reduction["unresolved_orbit_indices"] == [],
        "fixed_lattice": (
            lattice["determinant"] == expected["determinant"]
            and lattice["smith_diagonal_absolute"] == expected["smith"]
            and lattice["inertia"] == [1, expected["dimension"], 0]
        ),
        "coxeter_pairs": (
            case["coxeter_gate"]["all_pairs_coxeter"] is True
            and case["coxeter_gate"]["coxeter_edge_counts"]
            == expected["pair_counts"]
        ),
        "integral_reflections": case["restricted_mirrors"]
        ["all_active_reflections_integral_on_saturated_fixed_lattice"]
        is True,
        "facet_gram_hash": _payload_sha256(
            case["restricted_mirrors"]["positive_wall_gram"]
        )
        == expected["facet_gram_sha256"],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"section replay changed for {case['id']}: {checks}")
    return checks


def build(source: Path = SOURCE) -> dict[str, Any]:
    if _sha256(source) != SOURCE_SHA256:
        raise ArithmeticError("Alexeev ancillary source hash changed")
    source_text = source.read_text(encoding="utf-8")

    # _inventory_one is the already tested generic exact group/subgroup
    # engine.  Swap its small configuration table only during this replay.
    original_configs = inventory_tools.CONFIGS
    inventory_tools.CONFIGS = CONFIGS
    try:
        examples = [
            inventory_tools._inventory_one(source_text, key) for key in CONFIGS
        ]
    finally:
        inventory_tools.CONFIGS = original_configs

    for example in examples:
        key = tuple(map(int, example["key"]))
        if _high_signature(example) != EXPECTED_HIGH_CLASSES[key]:
            raise ArithmeticError(f"high subgroup queue changed for {key}")

    cases = []
    for example in examples:
        key = tuple(map(int, example["key"]))
        wall_gram = _extract_matrix(source_text, key)
        source_interior = (
            _all_wall_sum_interior(wall_gram, int(example["ambient"]["rank"]))
            if key == (15, 7, 1)
            else None
        )
        for class_record in example["subgroup_classes"]:
            if not class_record["proper_fixed_section_above_h9"]:
                continue
            case = _screen_case(
                key=key,
                example=example,
                class_record=class_record,
                wall_gram=wall_gram,
                source_interior=source_interior,
            )
            case["construction_checks"] = _audit_section(case)
            cases.append(case)
            print(
                json.dumps(
                    {
                        "id": case["id"],
                        "dimension": case["case"]["hyperbolic_dimension"],
                        "facets": case["exact_reduction"][
                            "active_hyperbolic_facet_count"
                        ],
                        "decision": case["decision"]["classification"],
                    }
                ),
                flush=True,
            )

    if {case["id"] for case in cases} != set(EXPECTED_SECTIONS):
        raise ArithmeticError("the complete high section identifier set changed")
    return {
        "schema": SCHEMA,
        "source": {
            "path": str(source),
            "sha256": SOURCE_SHA256,
            "arxiv": "2209.09110v4",
            "attachment": "all_cases_in_2el_paper.txt",
        },
        "scope": {
            "ambient_keys": [list(key) for key in CONFIGS],
            "minimum_fixed_rank": 11,
            "full_weighted_diagram_groups": True,
            "all_subgroup_conjugacy_classes": True,
            "geometric_automorphism_membership_claimed": False,
            "finite_volume_claimed": False,
        },
        "inventory": {
            "examples": examples,
            "weighted_group_orders": {
                "18,4,0": 1152,
                "15,11,1": 72,
                "15,7,1": 120,
            },
            "subgroup_class_counts": {
                "18,4,0": 221,
                "15,11,1": 26,
                "15,7,1": 19,
            },
            "proper_high_class_count": 6,
        },
        "cases": cases,
        "summary": {
            "case_count": len(cases),
            "dimension_histogram": {
                str(key): value
                for key, value in sorted(
                    Counter(
                        case["case"]["hyperbolic_dimension"] for case in cases
                    ).items()
                )
            },
            "facet_count_histogram": {
                str(key): value
                for key, value in sorted(
                    Counter(
                        case["exact_reduction"]["active_hyperbolic_facet_count"]
                        for case in cases
                    ).items()
                )
            },
            "all_six_strict_integral_coxeter_candidates": all(
                all(case["construction_checks"].values()) for case in cases
            ),
        },
        "all_checks_pass": all(
            all(case["construction_checks"].values()) for case in cases
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build(args.source)
    if args.check_archive:
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt remaining-stock construction changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
