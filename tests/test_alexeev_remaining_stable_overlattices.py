from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_alexeev_remaining_stable_overlattices import (
    OUTPUT,
    SCHEMA,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def _archive() -> dict:
    return json.loads(OUTPUT.read_text())


def test_stable_overlattice_archive_rebuilds_exactly() -> None:
    assert build_certificate() == _archive()


def test_all_six_discriminant_module_denominators_are_closed() -> None:
    payload = _archive()
    assert payload["schema"] == SCHEMA
    expected = {
        "A18_4_0_C03": ([1, 2, 2, 2, 4], [1, 2], 4, -16),
        "A18_4_0_C04": ([1, 2], [1], 2, 256),
        "A18_4_0_C06": ([1, 4], [1, 4], 4, -4),
        "A18_4_0_C18": ([1, 2], [1, 2], 2, 32),
        "A15_11_1_C02": ([1, 4], [1, 4], 4, 512),
        "A15_7_1_C01": ([1, 4], [1, 4], 4, 32),
    }
    for case in payload["cases"]:
        integral, even, index, determinant = expected[case["id"]]
        audit = case["exhaustive_discriminant_module_audit"]
        maximal = case["unique_maximal_reflection_stable_overlattice"]
        assert audit["stable_integral_subgroup_orders"] == integral
        assert audit["stable_even_subgroup_orders"] == even
        assert audit["maximal_stable_subgroup_positions_zero_based"]
        assert maximal["unique"]
        assert maximal["index_over_original_lattice"] == index
        assert maximal["determinant"] == determinant
        assert all(maximal["checks"].values())


def test_literal_normalizations_and_maximal_group_indices() -> None:
    payload = _archive()
    cases = {case["id"]: case for case in payload["cases"]}
    assert all(
        cases["A18_4_0_C04"]["normalization_to_strongly_squarefree_class"][
            "checks"
        ].values()
    )
    assert all(
        cases["A15_11_1_C02"]["normalization_to_strongly_squarefree_class"][
            "checks"
        ].values()
    )
    assert payload["summary"]["exact_numerical_indices"] == {
        "A15_11_1_C02": 6720,
        "A15_7_1_C01": 2016,
        "A18_4_0_C04": 130560,
        "A18_4_0_C18": 144,
    }
    assert payload["summary"]["zero_characteristic_index_values_open"] == [
        "A18_4_0_C03",
        "A18_4_0_C06",
    ]
    assert payload["summary"]["new_arithmetic_commensurability_class_count"] == 0
