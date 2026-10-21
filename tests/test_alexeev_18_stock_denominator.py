from __future__ import annotations

import json

from scripts.certify_alexeev_18_stock_denominator import (
    EXCEPTIONAL_KEY,
    EXPECTED_ALL_KEYS,
    EXPECTED_STANDARD_ROWS,
    OUTPUT,
)


def test_seventeen_all_subgroup_denominators_and_exceptional_scope() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert payload["primary_source"]["sha256"] == (
        "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
    )
    assert payload["summary"] == {
        "all_18_full_subgroup_conjugacy_denominators_closed": False,
        "all_18_high_rank_fishing_denominators_closed": True,
        "all_subgroup_closed_proper_high_class_count": 11,
        "all_subgroup_closed_stock_count": 17,
        "all_subgroup_closed_total_conjugacy_classes": 469,
        "exceptional_cyclic_high_class_count": 4,
        "exceptional_element_class_count": 28,
        "exceptional_high_noncyclic_class_count": 1,
        "exceptional_stock_count": 1,
        "explicit_Alexeev_matrix_count": 18,
    }

    rows = payload["seventeen_full_all_subgroup_stocks"]
    assert len(rows) == 17
    assert {tuple(row["key"]) for row in rows} == (
        EXPECTED_ALL_KEYS - {EXCEPTIONAL_KEY}
    )
    assert all(row["all_subgroup_conjugacy_classes_enumerated"] for row in rows)
    assert sum(row["subgroup_conjugacy_class_count"] for row in rows) == 469
    assert sum(
        row["proper_fixed_rank_at_least_11_class_count"] for row in rows
    ) == 11
    assert {
        tuple(row["key"]): (
            row["wall_count"],
            row["weighted_diagram_group_order"],
            row["subgroup_conjugacy_class_count"],
            row["proper_fixed_rank_at_least_11_class_count"],
        )
        for row in rows
    } == EXPECTED_STANDARD_ROWS

    exceptional = payload["exceptional_210_wall_stock"]
    assert tuple(exceptional["key"]) == EXCEPTIONAL_KEY
    assert exceptional["wall_count"] == 210
    assert exceptional["weighted_diagram_group_order"] == 241_920
    assert exceptional["element_conjugacy_class_count"] == 28
    assert exceptional["all_element_classes_enumerated"] is True
    assert exceptional["nonidentity_cyclic_high_class_names"] == [
        "2a",
        "2c",
        "2d",
        "3b",
    ]
    assert exceptional["complete_high_rank_noncyclic_filter"] is True
    assert exceptional["unique_high_noncyclic_class"] == (
        "V4 with fusion 2a+2c+2c"
    )
    assert exceptional["all_subgroup_conjugacy_classes_enumerated"] is False


def test_ten_small_stocks_supply_the_zero_high_rank_block() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rows = payload["seventeen_full_all_subgroup_stocks"]
    zero_rows = [
        row
        for row in rows
        if row["source_certificate"]
        == "alexeev_ten_small_zero_high_denominator.json"
    ]
    assert len(zero_rows) == 10
    assert sum(row["subgroup_conjugacy_class_count"] for row in zero_rows) == 49
    assert max(row["proper_fixed_rank_at_least_11_class_count"] for row in zero_rows) == 0
