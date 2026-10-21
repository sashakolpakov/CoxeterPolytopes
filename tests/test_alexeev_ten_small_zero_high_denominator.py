from __future__ import annotations

import json

from scripts.certify_alexeev_ten_small_zero_high_denominator import (
    CONFIGS,
    EXPECTED_CLASS_SIGNATURES,
    OUTPUT,
)


def test_ten_complete_small_stocks_have_no_proper_high_fixed_section() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert payload["source"]["sha256"] == (
        "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
    )
    assert payload["summary"] == {
        "all_ten_stocks_closed_negatively_above_H9": True,
        "maximum_proper_fixed_rank": 10,
        "proper_fixed_rank_at_least_11_count": 0,
        "proper_subgroup_fixed_rank_histogram": {
            "3": 2,
            "4": 7,
            "5": 4,
            "6": 10,
            "7": 2,
            "8": 4,
            "9": 4,
            "10": 6,
        },
        "stock_count": 10,
        "total_subgroup_conjugacy_class_count": 49,
    }
    by_key = {tuple(row["key"]): row for row in payload["examples"]}
    assert set(by_key) == set(CONFIGS) == set(EXPECTED_CLASS_SIGNATURES)
    for key, row in by_key.items():
        group = row["weighted_diagram_group"]
        assert group["order"] == CONFIGS[key]["diagram_group_order"]
        assert len(group["elements"]) == group["order"]
        assert group["all_actions_integral_unimodular_Lorentz"] is True
        signature = [
            (
                subgroup["class_index"],
                subgroup["order"],
                subgroup["conjugacy_class_size"],
                subgroup["fixed_rank"],
                subgroup["cyclic"],
            )
            for subgroup in row["subgroup_classes"]
        ]
        assert signature == EXPECTED_CLASS_SIGNATURES[key]
        assert all(
            subgroup["order"] == 1 or subgroup["fixed_rank"] <= 10
            for subgroup in row["subgroup_classes"]
        )
        assert row["summary"]["proper_fixed_rank_at_least_11_count"] == 0


def test_trivial_weighted_groups_are_frozen_without_gap_degree_ambiguity() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    by_key = {tuple(row["key"]): row for row in payload["examples"]}
    for key in ((11, 11, 1), (12, 10, 1)):
        row = by_key[key]
        assert row["weighted_diagram_group"]["order"] == 1
        assert len(row["subgroup_classes"]) == 1
        identity = row["subgroup_classes"][0]
        assert identity["class_index"] == 0
        assert identity["order"] == 1
        assert identity["conjugacy_class_size"] == 1
        assert identity["representative_generators"] == []
