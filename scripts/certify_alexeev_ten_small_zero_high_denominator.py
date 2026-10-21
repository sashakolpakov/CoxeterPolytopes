#!/usr/bin/env python3
"""Close the full subgroup denominator for ten small Alexeev chambers.

These are the ten explicit ancillary matrices not covered by the previously
frozen five medium chambers, the dual 118-wall pair, or the exceptional
210-wall Borcherds chamber.  The exact complete weighted wall-Gram group and
every subgroup conjugacy class are reconstructed for each stock.  No proper
subgroup fixes rank 11 or more.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries as tools
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    SOURCE,
    SOURCE_SHA256,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_ten_small_zero_high_denominator.json"
)
SCHEMA = "coxeter.alexeev_ten_small_zero_high_denominator.v1"


CONFIGS: dict[tuple[int, int, int], dict[str, Any]] = {
    (11, 11, 1): {
        "wall_count": 12,
        "rank": 11,
        "determinant": 2048,
        "diagram_group_name": "trivial",
        "diagram_group_order": 1,
        "geometric_context": "Alexeev complete reflective 2-elementary chamber",
    },
    (12, 10, 1): {
        "wall_count": 14,
        "rank": 12,
        "determinant": -1024,
        "diagram_group_name": "trivial",
        "diagram_group_order": 1,
        "geometric_context": "Alexeev complete reflective 2-elementary chamber",
    },
    (13, 9, 1): {
        "wall_count": 19,
        "rank": 13,
        "determinant": 512,
        "diagram_group_name": "S3",
        "diagram_group_order": 6,
        "geometric_context": "Alexeev complete reflective 2-elementary chamber",
    },
    (14, 8, 0): {
        "wall_count": 17,
        "rank": 14,
        "determinant": -256,
        "diagram_group_name": "S3",
        "diagram_group_order": 6,
        "geometric_context": "Alexeev complete reflective 2-elementary chamber",
    },
    (14, 8, 1): {
        "wall_count": 24,
        "rank": 14,
        "determinant": -256,
        "diagram_group_name": "S4",
        "diagram_group_order": 24,
        "geometric_context": "Alexeev complete reflective 2-elementary chamber",
    },
    (13, 11, 1): {
        "wall_count": 16,
        "rank": 13,
        "determinant": 2048,
        "diagram_group_name": "C2",
        "diagram_group_order": 2,
        "geometric_context": "Alexeev complete reflective 2-elementary chamber",
    },
    (14, 10, 1): {
        "wall_count": 21,
        "rank": 14,
        "determinant": -1024,
        "diagram_group_name": "D6 (order 12)",
        "diagram_group_order": 12,
        "geometric_context": "Alexeev complete reflective 2-elementary chamber",
    },
    (14, 12, 1): {
        "wall_count": 17,
        "rank": 14,
        "determinant": -4096,
        "diagram_group_name": "C2",
        "diagram_group_order": 2,
        "geometric_context": "Alexeev complete reflective 2-elementary chamber",
    },
    (15, 13, 1): {
        "wall_count": 27,
        "rank": 15,
        "determinant": 8192,
        "diagram_group_name": "S3",
        "diagram_group_order": 6,
        "geometric_context": "Alexeev complete reflective 2-elementary chamber",
    },
    (16, 14, 1): {
        "wall_count": 29,
        "rank": 16,
        "determinant": -16384,
        "diagram_group_name": "D6 (order 12)",
        "diagram_group_order": 12,
        "geometric_context": "Alexeev complete reflective 2-elementary chamber",
    },
}


EXPECTED_CLASS_SIGNATURES = {
    (11, 11, 1): [(0, 1, 1, 11, True)],
    (12, 10, 1): [(0, 1, 1, 12, True)],
    (13, 9, 1): [
        (0, 1, 1, 13, True), (1, 2, 3, 9, True),
        (2, 3, 1, 5, True), (3, 6, 1, 5, False),
    ],
    (14, 8, 0): [
        (0, 1, 1, 14, True), (1, 2, 3, 10, True),
        (2, 3, 1, 6, True), (3, 6, 1, 6, False),
    ],
    (14, 8, 1): [
        (0, 1, 1, 14, True), (1, 2, 3, 10, True),
        (2, 2, 6, 10, True), (3, 3, 4, 6, True),
        (4, 4, 1, 8, False), (5, 4, 3, 8, False),
        (6, 4, 3, 6, True), (7, 6, 4, 6, False),
        (8, 8, 3, 6, False), (9, 12, 1, 4, False),
        (10, 24, 1, 4, False),
    ],
    (13, 11, 1): [(0, 1, 1, 13, True), (1, 2, 1, 7, True)],
    (14, 10, 1): [
        (0, 1, 1, 14, True), (1, 2, 1, 10, True),
        (2, 2, 3, 10, True), (3, 2, 3, 8, True),
        (4, 3, 1, 6, True), (5, 4, 3, 7, False),
        (6, 6, 1, 4, True), (7, 6, 1, 6, False),
        (8, 6, 1, 4, False), (9, 12, 1, 4, False),
    ],
    (14, 12, 1): [(0, 1, 1, 14, True), (1, 2, 1, 8, True)],
    (15, 13, 1): [
        (0, 1, 1, 15, True), (1, 2, 3, 9, True),
        (2, 3, 1, 5, True), (3, 6, 1, 4, False),
    ],
    (16, 14, 1): [
        (0, 1, 1, 16, True), (1, 2, 3, 10, True),
        (2, 2, 1, 9, True), (3, 2, 3, 9, True),
        (4, 3, 1, 6, True), (5, 4, 3, 6, False),
        (6, 6, 1, 5, False), (7, 6, 1, 3, True),
        (8, 6, 1, 4, False), (9, 12, 1, 3, False),
    ],
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _signature(example: dict[str, Any]) -> list[tuple[int, int, int, int, bool]]:
    return [
        (
            int(row["class_index"]),
            int(row["order"]),
            int(row["conjugacy_class_size"]),
            int(row["fixed_rank"]),
            bool(row["cyclic"]),
        )
        for row in example["subgroup_classes"]
    ]


def build(source: Path = SOURCE) -> dict[str, Any]:
    if _sha256(source) != SOURCE_SHA256:
        raise ArithmeticError("Alexeev ancillary source hash changed")
    source_text = source.read_text(encoding="utf-8")
    original_configs = tools.CONFIGS
    tools.CONFIGS = CONFIGS
    try:
        examples = [tools._inventory_one(source_text, key) for key in CONFIGS]
    finally:
        tools.CONFIGS = original_configs

    checks = {}
    proper_fixed_rank_histogram: Counter[int] = Counter()
    for example in examples:
        key = tuple(map(int, example["key"]))
        signature_matches = _signature(example) == EXPECTED_CLASS_SIGNATURES[key]
        no_proper_high = not any(
            row["order"] > 1 and row["fixed_rank"] >= 11
            for row in example["subgroup_classes"]
        )
        every_class_stored = (
            len(example["subgroup_classes"])
            == example["summary"]["subgroup_conjugacy_class_count"]
        )
        group_complete = (
            len(example["weighted_diagram_group"]["elements"])
            == CONFIGS[key]["diagram_group_order"]
            and example["weighted_diagram_group"]
            ["all_actions_integral_unimodular_Lorentz"]
            is True
        )
        checks[",".join(map(str, key))] = {
            "exact_class_signature": signature_matches,
            "every_subgroup_class_stored": every_class_stored,
            "full_weighted_group_reconstructed": group_complete,
            "zero_proper_fixed_rank_at_least_11": no_proper_high,
        }
        if not all(checks[",".join(map(str, key))].values()):
            raise ArithmeticError(f"small-stock denominator changed for {key}")
        proper_fixed_rank_histogram.update(
            int(row["fixed_rank"])
            for row in example["subgroup_classes"]
            if row["order"] > 1
        )

    return {
        "schema": SCHEMA,
        "source": {
            "path": str(source),
            "sha256": SOURCE_SHA256,
            "arxiv": "2209.09110v4",
            "attachment": "all_cases_in_2el_paper.txt",
        },
        "scope": {
            "matrix_keys": [list(key) for key in CONFIGS],
            "full_weighted_wall_gram_automorphism_group": True,
            "all_subgroups_modulo_conjugacy": True,
            "proper_fixed_rank_threshold": 11,
            "excluded_exceptional_matrix": [22, 2, 0],
            "exception_reason": (
                "its order-241920 group has a separately frozen all-cyclic "
                "and targeted high-V4 certificate, not an all-subgroup census"
            ),
        },
        "examples": examples,
        "checks_by_key": checks,
        "summary": {
            "stock_count": len(examples),
            "total_subgroup_conjugacy_class_count": sum(
                len(example["subgroup_classes"]) for example in examples
            ),
            "proper_subgroup_fixed_rank_histogram": {
                str(rank): count
                for rank, count in sorted(proper_fixed_rank_histogram.items())
            },
            "maximum_proper_fixed_rank": max(proper_fixed_rank_histogram),
            "proper_fixed_rank_at_least_11_count": 0,
            "all_ten_stocks_closed_negatively_above_H9": True,
        },
        "all_checks_pass": all(
            all(row.values()) for row in checks.values()
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
            raise ArithmeticError("rebuilt ten-small-stock denominator changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
