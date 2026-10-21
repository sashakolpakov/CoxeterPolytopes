from __future__ import annotations

import json

from scripts.certify_alexeev_remaining_fixed_sections_lineage import (
    EXPECTED_DUPLICATE_PERMUTATION,
    EXPECTED_EULER,
    OUTPUT,
)


def test_six_case_bounded_lineage_certificate() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert {
        row["id"]: row["solomon"]["rational_euler_characteristic"]
        for row in payload["cases"]
    } == EXPECTED_EULER
    assert payload["prior_ledger_exact_duplicate"][
        "source_to_target_facet_permutation_one_based"
    ] == EXPECTED_DUPLICATE_PERMUTATION
    assert payload["within_batch_same_header_comparisons"] == [
        {
            "same_signed_normalized_gram": False,
            "source": "A18_4_0_C18",
            "source_to_target_facet_permutation_zero_based": None,
            "target": "A15_11_1_C02",
        }
    ]


def test_six_case_coxiter_and_priority_scope() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    coxiter = payload["coxiter_bounded_comparison"]
    assert coxiter["targets"]["H10_facets_28"][
        "doubling_header_match_count"
    ] == 2
    assert coxiter["no_exact_relation_graph_match"] is True
    assert payload["scope"]["historical_ceiling"]["doi"] == (
        "https://doi.org/10.2140/gt.2006.10.737"
    )
    stable = payload["reflection_stable_overlattice_status"]
    assert stable["permanent_exhaustive_certificate"] == (
        "data/enumerations/"
        "alexeev_remaining_high_fixed_section_stable_overlattices.json"
    )
    assert stable["schema"] == (
        "coxeter.alexeev_remaining_high_fixed_sections."
        "stable_overlattices.v1"
    )
    assert stable["status"] == "closed for all six displayed reflection groups"
    assert stable["summary"]["exact_numerical_indices"] == {
        "A15_11_1_C02": 6720,
        "A15_7_1_C01": 2016,
        "A18_4_0_C04": 130560,
        "A18_4_0_C18": 144,
    }
    assert stable["summary"]["zero_characteristic_index_values_open"] == [
        "A18_4_0_C03",
        "A18_4_0_C06",
    ]


def test_six_unique_stable_overlattices_and_full_chambers() -> None:
    payload = json.loads(OUTPUT.read_text(encoding="utf-8"))
    expected = {
        "A18_4_0_C03": (4, -16, "List11", 15, None),
        "A18_4_0_C04": (2, 256, "List10", 11, 130560),
        "A18_4_0_C06": (4, -4, "List13", 4, None),
        "A18_4_0_C18": (2, 32, "List10", 20, 144),
        "A15_11_1_C02": (4, 512, "List10", 6, 6720),
        "A15_7_1_C01": (4, 32, "List10", 20, 2016),
    }
    for case in payload["cases"]:
        index, determinant, census, entry, group_index = expected[case["id"]]
        stable = case["reflection_stable_overlattice"]
        chamber = case["containing_full_reflection_chamber"]
        assert stable["unique_maximal"] is True
        assert stable["index_over_fixed_lattice"] == index
        assert stable["determinant"] == determinant
        assert chamber["scharlau_kirschmer_list"] == census
        assert chamber["entry_one_based"] == entry
        assert chamber["section_to_full_index"] == group_index
