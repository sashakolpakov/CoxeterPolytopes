#!/usr/bin/env python3
"""Verify the production record-18 a=3 addition scout archive."""

from __future__ import annotations

from collections import Counter
import json
from math import comb
from pathlib import Path
from typing import Mapping


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record18_a3_addition_scout.json"
POSITIVE = ROOT / "data/enumerations/h9_14_record18_positive_A4_scout.json"


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record18_a3_addition_scout":
        raise ValueError("wrong record-18 a3 scout kind")
    if (
        document.get("theorem_certificate") is not False
        or document.get("independent_replay_complete") is not False
        or document.get("compact_only") is not True
        or document.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise ValueError("record-18 a3 scout scope was over-promoted")
    rows = document.get("orbit_audits")
    if not isinstance(rows, list) or len(rows) != 40:
        raise ValueError("record-18 a3 orbit audit list changed")
    positive = json.loads(POSITIVE.read_text(encoding="utf-8"))["branches"][2][
        "orbits"
    ]
    if [row["orbit_index"] for row in rows] != list(range(40)):
        raise ValueError("record-18 a3 orbit indices changed")
    if [row["positive_orbit_key"] for row in rows] != [
        row["key"] for row in positive
    ]:
        raise ValueError("record-18 a3 orbit keys are not linked to the positive census")
    if [row["representative_A4s"] for row in rows] != [
        row["representative"] for row in positive
    ]:
        raise ValueError("record-18 a3 representatives changed")
    if Counter(
        (row["forced_A5_shadow_count"], row["additional_A5_count"])
        for row in rows
    ) != {(24, 2): 29, (23, 3): 10, (22, 4): 1}:
        raise ValueError("record-18 a3 shadow/addition histogram changed")
    for row in rows:
        if row["forced_A5_shadow_count"] + row["additional_A5_count"] != 26:
            raise ValueError("record-18 a3 scheduled A5 count changed")
        if row["raw_unordered_addition_count"] != comb(
            row["candidate_count"], row["additional_A5_count"]
        ):
            raise ValueError("record-18 a3 raw addition count changed")
        if sum(row["facet_count_histogram"].values()) != row[
            "no_small_transversal_addition_count"
        ]:
            raise ValueError("record-18 a3 facet histogram does not sum")
        if sum(row["blocker_size_histogram"].values()) != row[
            "no_small_transversal_addition_count"
        ]:
            raise ValueError("record-18 a3 blocker histogram does not sum")
        if row["scheduled_pure_target_count"] != 0:
            raise ValueError("record-18 a3 archive contains an undeclared target")
    expected = {
        "aggregate_raw_addition_count": 14_933_639_143,
        "aggregate_interval_feasible_addition_count": 2_392_528_542,
        "aggregate_no_small_transversal_addition_count": 6_011,
        "scheduled_pure_104_facet_target_count": 0,
        "a3_branch_empty_before_topology_and_labels": True,
    }
    if any(document.get(key) != value for key, value in expected.items()):
        raise ValueError("record-18 a3 aggregate transcript changed")
    if sum(row["raw_unordered_addition_count"] for row in rows) != expected[
        "aggregate_raw_addition_count"
    ] or sum(row["vertex_interval_feasible_addition_count"] for row in rows) != expected[
        "aggregate_interval_feasible_addition_count"
    ] or sum(row["no_small_transversal_addition_count"] for row in rows) != expected[
        "aggregate_no_small_transversal_addition_count"
    ]:
        raise ValueError("record-18 a3 aggregate sums are inconsistent")
    return {
        "record_index": 18,
        "a3_positive_orbit_count": 40,
        "production_a3_branch_empty": True,
        "independent_replay_complete": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def main() -> None:
    print(json.dumps(verify(json.loads(ARCHIVE.read_text(encoding="utf-8"))), indent=2))


if __name__ == "__main__":
    main()
