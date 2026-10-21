from __future__ import annotations

import json

from scripts.certify_triangular_looijenga_h4_worldwide_lineage import OUTPUT


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_six_and_seven_facet_types_are_not_new_abstract_polytopes() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    conclusion = payload["main_conclusion"]
    assert conclusion["T6a_T6b_T7_new_abstract_polytopes"] is False
    assert conclusion["T6a_T6b_T7_new_AG_carriers_in_this_project"] is True

    lineage = payload["lineage_decisions"]
    assert lineage["T6a"]["exact_source_family"] == (
        "H4 pyramid family 07_02 at (k,l)=(2,4)"
    )
    assert lineage["T6b"]["exact_source_family"] == (
        "H4 pyramid family 07_01 at (k,l)=(2,4)"
    )
    assert lineage["T7"]["independent_complete_census_match"] == (
        "Ma--Zheng P_{2,1}"
    )
    assert all(
        lineage[type_id]["abstract_polytope_status"] == "previously published"
        for type_id in ("T6a", "T6b", "T7")
    )


def test_ma_zheng_complete_census_has_one_exact_t7_match() -> None:
    audit = _payload()["sources"]["ma_zheng"]
    assert audit["classification_complete_for_finite_volume_H4_seven_facets"] is True
    assert audit["result_row_count"] == 331
    assert audit["distinct_result_id_count"] == 331
    assert audit["exact_unlabelled_match_count"] == 1
    assert [match["id"] for match in audit["matches"]] == ["P_{2,1}"]
    assert audit["matches"][0][
        "lexicographically_first_candidate_to_table_permutation"
    ] == [1, 4, 5, 0, 3, 2, 6]
    assert audit["published_invariants"] == {
        "arithmetic": True,
        "cusp_count": 3,
        "volume": "pi^2/48",
    }


def test_exact_volumes_indices_and_mod_two_preimages() -> None:
    payload = _payload()
    p4 = payload["vinberg_P4"]
    assert p4["volume"]["volume_over_pi_squared"] == "1/1440"
    assert p4["mod_two_image"] == {"coxeter_type": "W(A4)=S5", "order": 120}

    expected = {
        "T6a": ("1/144", 10, "W(A2+A1)", 12),
        "T6b": ("1/72", 20, "W(A2)", 6),
        "T7": ("1/48", 30, "W(A1+A1)", 4),
        "T8": ("1/24", 60, "W(A1)", 2),
    }
    for type_id, (volume, index, image_type, image_order) in expected.items():
        record = payload["types"][type_id]
        assert record["volume"]["volume_over_pi_squared"] == volume
        assert record["index_in_vinberg_P4"] == index
        assert record["mod_two_image"] == {
            "coxeter_type": image_type,
            "index_in_W_A4_equals_S5": index,
            "order": image_order,
        }
        assert record["equals_full_preimage_of_mod_two_image"] is True
        assert record["arithmetic"] is True


def test_pinned_coxiter_matches_and_bounded_t8_exclusion() -> None:
    audit = _payload()["coxiter_bounded_audit"]
    assert audit["ordinary_file_count"] == 1_214
    assert audit["ordinary_H4_facet_histogram"] == {
        "5": 64,
        "6": 75,
        "7": 42,
        "8": 1,
    }
    ordinary = audit["ordinary_exact_matches"]
    assert [record["path"] for record in ordinary["T6a"]] == [
        "graphs/Tumarkin pyramids n+2/4-tum04_07_02-k=2-l=4.coxiter"
    ]
    assert [record["path"] for record in ordinary["T6b"]] == [
        "graphs/Tumarkin pyramids n+2/4-tum04_07_01-k=2-l=4.coxiter"
    ]
    assert [record["path"] for record in ordinary["T7"]] == [
        "graphs/doubly truncated orthoscheme/4-imohf85_13.coxiter"
    ]
    assert ordinary["T8"] == []
    assert audit["doubling_file_count"] == 2_959
    assert all(matches == [] for matches in audit["doubling_exact_matches"].values())
    assert audit["scope"] == "bounded pinned corpus, not a worldwide priority census"


def test_t8_priority_is_left_open_but_group_lineage_is_exact() -> None:
    payload = _payload()
    decision = payload["lineage_decisions"]["T8"]
    assert decision["worldwide_novelty_claimed"] is False
    assert decision[
        "complete_noncompact_H4_eight_facet_classification_available"
    ] is False
    assert decision["exact_group_theoretic_lineage"] == (
        "index-60 full preimage of W(A1) under W(P4)->W(A4)=S5 modulo two"
    )
    assert payload["main_conclusion"]["T8_worldwide_abstract_priority_resolved"] is False
