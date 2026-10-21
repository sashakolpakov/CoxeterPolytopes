#!/usr/bin/env python3
"""Compose the exact fishing denominator for all 18 Alexeev matrices.

Seventeen matrices have a complete weighted-diagram group and every subgroup
conjugacy class enumerated, split across four explicit source certificates.
The exceptional 210-wall ``(22,2,0)`` stock is kept logically separate: all
28 element classes (hence all cyclic fixed sections) are closed, and a
Reynolds/closure argument proves that the unique possible high-rank
noncyclic class is the already certified ``V4``.  We do *not* claim an
all-subgroup conjugacy census for that order-241920 group.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Any

from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    SOURCE,
    SOURCE_SHA256,
)


ROOT = Path(__file__).resolve().parents[1]
TEN_SMALL = ROOT / "data/enumerations/alexeev_ten_small_zero_high_denominator.json"
TWO_MEDIUM = (
    ROOT / "data/enumerations/alexeev_20_2_1_15_9_1_wall_symmetry_inventory.json"
)
THREE_REMAINING = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json"
)
TWO_DUAL_118 = (
    ROOT / "data/enumerations/alexeev_118_dual_h11_29_construction.json"
)
EXCEPTIONAL_CYCLIC = (
    ROOT / "data/enumerations/borcherds_char2_all_cyclic_element_completion.json"
)
EXCEPTIONAL_V4 = (
    ROOT / "data/enumerations/borcherds_char2_v4_h10_34_certificate.json"
)
OUTPUT = ROOT / "data/enumerations/alexeev_18_stock_denominator.json"
SCHEMA = "coxeter.alexeev_18_stock_denominator.v1"


EXPECTED_ALL_KEYS = {
    (11, 11, 1), (12, 10, 1), (13, 9, 1), (14, 8, 0), (14, 8, 1),
    (15, 7, 1), (16, 6, 1), (18, 4, 0), (20, 2, 1), (13, 11, 1),
    (14, 10, 1), (15, 9, 1), (22, 2, 0), (14, 12, 1), (15, 11, 1),
    (15, 13, 1), (16, 12, 1), (16, 14, 1),
}
EXCEPTIONAL_KEY = (22, 2, 0)
EXPECTED_STANDARD_ROWS = {
    (11, 11, 1): (12, 1, 1, 0),
    (12, 10, 1): (14, 1, 1, 0),
    (13, 9, 1): (19, 6, 4, 0),
    (13, 11, 1): (16, 2, 2, 0),
    (14, 8, 0): (17, 6, 4, 0),
    (14, 8, 1): (24, 24, 11, 0),
    (14, 10, 1): (21, 12, 10, 0),
    (14, 12, 1): (17, 2, 2, 0),
    (15, 7, 1): (60, 120, 19, 1),
    (15, 9, 1): (98, 336, 23, 1),
    (15, 11, 1): (66, 72, 26, 1),
    (15, 13, 1): (27, 6, 4, 0),
    (16, 6, 1): (118, 720, 56, 1),
    (16, 12, 1): (118, 720, 56, 1),
    (16, 14, 1): (29, 12, 10, 0),
    (18, 4, 0): (48, 1152, 221, 4),
    (20, 2, 1): (50, 120, 19, 2),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _source_keys(text: str) -> set[tuple[int, int, int]]:
    return {
        tuple(map(int, match.groups()))
        for match in re.finditer(
            r"coxeter_matrix\[\((\d+),\s*(\d+),\s*(\d+)\)\]", text
        )
    }


def _standard_row(
    *,
    key: tuple[int, int, int],
    walls: int,
    group_order: int,
    subgroup_classes: int,
    high_classes: int,
    source_name: str,
) -> dict[str, Any]:
    return {
        "key": list(key),
        "wall_count": int(walls),
        "weighted_diagram_group_order": int(group_order),
        "subgroup_conjugacy_class_count": int(subgroup_classes),
        "proper_fixed_rank_at_least_11_class_count": int(high_classes),
        "all_subgroup_conjugacy_classes_enumerated": True,
        "source_certificate": source_name,
    }


def build() -> dict[str, Any]:
    source_text = SOURCE.read_text(encoding="utf-8")
    if _sha256(SOURCE) != SOURCE_SHA256:
        raise ArithmeticError("Alexeev ancillary source hash changed")
    explicit_keys = _source_keys(source_text)
    if explicit_keys != EXPECTED_ALL_KEYS:
        raise ArithmeticError("the set of 18 explicit Alexeev matrices changed")

    paths = [
        TEN_SMALL,
        TWO_MEDIUM,
        THREE_REMAINING,
        TWO_DUAL_118,
        EXCEPTIONAL_CYCLIC,
        EXCEPTIONAL_V4,
    ]
    ten, two, three, dual, cyclic, v4 = map(_load, paths)
    if not (
        ten["all_checks_pass"] is True
        and ten["summary"]["all_ten_stocks_closed_negatively_above_H9"] is True
        and three["all_checks_pass"] is True
        and dual["all_checks_pass"] is True
        and cyclic["checks"]["all_element_classes_exhausted"] is True
        and v4["all_checks_pass"] is True
        and v4["checks"]["complete_high_subgroup_filter"] is True
    ):
        raise ArithmeticError("an upstream denominator certificate is incomplete")

    rows = []
    for example in ten["examples"]:
        rows.append(
            _standard_row(
                key=tuple(example["key"]),
                walls=example["walls"]["count"],
                group_order=example["weighted_diagram_group"]["order"],
                subgroup_classes=example["summary"]
                ["subgroup_conjugacy_class_count"],
                high_classes=example["summary"]
                ["proper_fixed_rank_at_least_11_count"],
                source_name=TEN_SMALL.name,
            )
        )
    for example in two["examples"]:
        rows.append(
            _standard_row(
                key=tuple(example["key"]),
                walls=example["walls"]["count"],
                group_order=example["weighted_diagram_group"]["order"],
                subgroup_classes=example["summary"]
                ["subgroup_conjugacy_class_count"],
                high_classes=example["summary"]
                ["proper_fixed_rank_at_least_11_count"],
                source_name=TWO_MEDIUM.name,
            )
        )
    for example in three["inventory"]["examples"]:
        rows.append(
            _standard_row(
                key=tuple(example["key"]),
                walls=example["walls"]["count"],
                group_order=example["weighted_diagram_group"]["order"],
                subgroup_classes=example["summary"]
                ["subgroup_conjugacy_class_count"],
                high_classes=example["summary"]
                ["proper_fixed_rank_at_least_11_count"],
                source_name=THREE_REMAINING.name,
            )
        )
    for stock in dual["stocks"]:
        rows.append(
            _standard_row(
                key=tuple(stock["key"]),
                walls=stock["source_wall_count"],
                group_order=stock["weighted_wall_group"]["order"],
                subgroup_classes=stock["subgroup_inventory"]
                ["gap_subgroup_conjugacy_class_count"],
                high_classes=stock["subgroup_inventory"]
                ["proper_fixed_rank_at_least_11_class_count"],
                source_name=TWO_DUAL_118.name,
            )
        )

    rows.sort(key=lambda row: tuple(row["key"]))
    full_keys = {tuple(row["key"]) for row in rows}
    expected_full = EXPECTED_ALL_KEYS - {EXCEPTIONAL_KEY}
    actual_standard_rows = {
        tuple(row["key"]): (
            row["wall_count"],
            row["weighted_diagram_group_order"],
            row["subgroup_conjugacy_class_count"],
            row["proper_fixed_rank_at_least_11_class_count"],
        )
        for row in rows
    }
    if not (
        len(rows) == 17
        and full_keys == expected_full
        and len(full_keys) == len(rows)
        and actual_standard_rows == EXPECTED_STANDARD_ROWS
        and sum(row["subgroup_conjugacy_class_count"] for row in rows) == 469
        and sum(row["proper_fixed_rank_at_least_11_class_count"] for row in rows)
        == 11
    ):
        raise ArithmeticError("the 17-stock all-subgroup composition changed")

    exceptional = {
        "key": list(EXCEPTIONAL_KEY),
        "wall_count": cyclic["ambient"]["wall_count"],
        "weighted_diagram_group_order": v4["subgroup_class_proof"]
        ["full_group_order"],
        "element_conjugacy_class_count": len(cyclic["classes"]),
        "all_element_classes_enumerated": cyclic["checks"]
        ["all_element_classes_exhausted"],
        "nonidentity_cyclic_high_class_names": cyclic["rank_gate"]
        ["nonidentity_high_class_names"],
        "nonidentity_cyclic_high_class_count": cyclic["rank_gate"]
        ["nonidentity_high_class_count"],
        "complete_high_rank_noncyclic_filter": v4["checks"]
        ["complete_high_subgroup_filter"],
        "unique_high_noncyclic_class": "V4 with fusion 2a+2c+2c",
        "all_subgroup_conjugacy_classes_enumerated": False,
        "scope_warning": (
            "all 28 element classes and the complete high-rank noncyclic "
            "filter are proved; no all-subgroup conjugacy census is claimed"
        ),
        "source_certificates": [EXCEPTIONAL_CYCLIC.name, EXCEPTIONAL_V4.name],
    }
    if not (
        exceptional["wall_count"] == 210
        and exceptional["weighted_diagram_group_order"] == 241920
        and exceptional["element_conjugacy_class_count"] == 28
        and exceptional["nonidentity_cyclic_high_class_count"] == 4
        and exceptional["complete_high_rank_noncyclic_filter"] is True
    ):
        raise ArithmeticError("the exceptional 210-wall scope changed")

    return {
        "schema": SCHEMA,
        "primary_source": {
            "path": str(SOURCE),
            "sha256": SOURCE_SHA256,
            "explicit_matrix_count": len(explicit_keys),
            "explicit_matrix_keys": [list(key) for key in sorted(explicit_keys)],
        },
        "source_certificates": {
            path.name: {
                "path": str(path.relative_to(ROOT)),
                "sha256": _sha256(path),
            }
            for path in paths
        },
        "seventeen_full_all_subgroup_stocks": rows,
        "exceptional_210_wall_stock": exceptional,
        "summary": {
            "explicit_Alexeev_matrix_count": 18,
            "all_subgroup_closed_stock_count": 17,
            "all_subgroup_closed_total_conjugacy_classes": 469,
            "all_subgroup_closed_proper_high_class_count": 11,
            "exceptional_stock_count": 1,
            "exceptional_element_class_count": 28,
            "exceptional_cyclic_high_class_count": 4,
            "exceptional_high_noncyclic_class_count": 1,
            "all_18_high_rank_fishing_denominators_closed": True,
            "all_18_full_subgroup_conjugacy_denominators_closed": False,
        },
        "claim": (
            "All subgroup conjugacy classes are exhausted for 17 of 18 "
            "explicit Alexeev matrices. For (22,2,0), all cyclic classes and "
            "a complete theorem-level high-noncyclic filter close the H10+ "
            "fishing denominator, without claiming all subgroup classes."
        ),
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build()
    if args.check_archive:
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt 18-stock denominator changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
