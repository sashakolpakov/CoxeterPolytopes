from __future__ import annotations

import json

from scripts.certify_alexeev_118_dual_h11_29_stable_overlattices import (
    OUTPUT,
    SCHEMA,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text())


def test_both_actual_lattices_have_closed_stable_overlattice_denominators() -> None:
    payload = _archive()
    assert payload["schema"] == SCHEMA
    cases = {case["id"]: case for case in payload["cases"]}
    assert set(cases) == {"A16_6_1_C01", "A16_12_1_C01"}
    for case in cases.values():
        audit = case["exhaustive_discriminant_module_audit"]
        maximal = case["unique_maximal_reflection_stable_overlattice"]
        assert audit["stable_integral_subgroup_orders"] == [1, 2, 4, 8]
        assert audit["stable_even_subgroup_orders"] == [1, 4]
        assert maximal["index_over_original_lattice"] == 8
        assert all(maximal["checks"].values())
    assert cases["A16_6_1_C01"]["unique_maximal_reflection_stable_overlattice"][
        "determinant"
    ] == -4
    assert cases["A16_12_1_C01"]["unique_maximal_reflection_stable_overlattice"][
        "determinant"
    ] == -256


def test_two_stocks_normalize_to_one_literal_integral_reflection_system() -> None:
    payload = _archive()
    common = payload["common_normalized_reflection_system"]
    assert all(common["checks"].values())
    assert sorted(common["source_to_dual_facet_permutation_one_based"]) == list(
        range(1, 30)
    )
    assert len(common["source_root_factors_before_A"]) == 29
    assert payload["summary"]["one_common_normalized_integral_reflection_system"]
    lineage = payload["maximal_reflection_lineage"]
    assert lineage["entry_one_based"] == 5
    assert lineage["full_chamber_facets"] == 14
    assert lineage["full_chamber_finite_vertices"] == 56
    assert lineage["full_chamber_cusps"] == 5
    assert all(lineage["checks"].values())


def test_same_header_collision_is_resolved_by_discriminant_form() -> None:
    collision = _archive()["same_header_collision_with_A18_4_0_C03"]
    assert all(collision["checks"].values())
    assert "3/4" in collision["A18_C03_discriminant_q_mod_2_histogram"]
    assert "3/4" not in collision["A16_6_discriminant_q_mod_2_histogram"]
    assert _archive()["summary"]["new_arithmetic_commensurability_class_count"] == 0
