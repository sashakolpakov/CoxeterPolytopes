from __future__ import annotations

import json

from scripts.certify_alexeev_parked_rank10_h9_lineage_overclosures import (
    IMHOF_MATCH,
    OUTPUT,
    SCHEMA,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == _archive()


def test_five_unique_maximal_stable_overclosures_and_named_classes() -> None:
    payload = _archive()
    assert payload["schema"] == SCHEMA
    assert payload["all_checks_pass"]
    assert all(payload["checks"].values())
    expected = {
        "A14_8_0_C01": (4, -64, [10, 6, 0], "U(2) + D4(-1)^2", 23),
        "A14_8_1_C02": (8, -16, [10, 4, 0], "U + D4(-1)^2", 23),
        "A14_10_1_C01": (4, -256, [10, 8, 0], "U + E8(-2)", 7),
        "A14_10_1_C02": (8, -64, [10, 6, 0], "U(2) + D4(-1)^2", 23),
        "A16_14_1_C01": (32, -16, [10, 4, 0], "U + D4(-1)^2", 23),
    }
    for case in payload["cases"]:
        maximal = case["unique_maximal_reflection_stable_overclosure"]
        normalized = case["normalization_to_strongly_squarefree_class"]
        index, determinant, invariants, isometry_type, sk_entry = expected[case["id"]]
        assert maximal["unique"]
        assert maximal["index_over_original_lattice"] == index
        assert maximal["determinant"] == determinant
        assert maximal["two_elementary_invariants_r_a_delta"] == invariants
        assert maximal["isometry_type"] == isometry_type
        assert all(maximal["checks"].values())
        assert normalized["scharlau_kirschmer_entry_one_based"] == sk_entry
        assert normalized["section_reflection_group_has_finite_index_in_full_group"]
        assert normalized["numerical_index"] is None
        assert all(normalized["checks"].values())
        assert case["rational_similarity_to_odd_unimodular_I_9_1"][
            "rationally_similar"
        ]


def test_bounded_coxiter_collision_is_exactly_im_hof() -> None:
    payload = _archive()
    comparison = payload["coxiter_bounded_comparison"]
    assert comparison["ordinary_coxiter_file_count"] == 1214
    assert comparison["doubling_coxiter_file_count"] == 2959
    rows = comparison["comparisons"]
    assert rows["A14_8_0_C01"]["exact_labelled_graph_match_paths"] == [
        IMHOF_MATCH
    ]
    assert rows["A14_8_0_C01"]["exact_match"][
        "candidate_to_reference_wall_permutation_one_based"
    ] == [6, 11, 10, 9, 8, 1, 2, 3, 4, 5, 7, 12]
    assert all(
        not rows[identifier]["exact_labelled_graph_match_paths"]
        for identifier in (
            "A14_8_1_C02",
            "A14_10_1_C01",
            "A14_10_1_C02",
            "A16_14_1_C01",
        )
    )
    assert rows["A14_8_1_C02"]["same_header_record_count"] == 0
    assert rows["A14_10_1_C01"]["same_header_record_count"] == 8
    assert rows["A14_10_1_C02"]["same_header_record_count"] == 8
    assert rows["A16_14_1_C01"]["same_header_record_count"] == 8


def test_parallel_redoubling_scope_and_named_stock_exclusions() -> None:
    payload = _archive()
    cases = {case["id"]: case for case in payload["cases"]}
    expected_pairs = {
        "A14_8_0_C01": [],
        "A14_8_1_C02": [[15, 16]],
        "A14_10_1_C01": [],
        "A14_10_1_C02": [[6, 15]],
        "A16_14_1_C01": [[1, 14], [8, 15]],
    }
    for identifier, pairs in expected_pairs.items():
        redoubling = cases[identifier]["redoubling"]
        assert redoubling["disjoint_doubling_wall_pairs_one_based"] == pairs
        assert redoubling["redoublable_by_Allcock_definition"] is bool(pairs)
        assert redoubling[
            "Allcock_Corollary_4_2_reflection_subgroups_of_every_positive_index"
        ] is bool(pairs)
        assert not redoubling[
            "Allcock_Theorem_4_3_exponential_count_applies_from_certified_walls"
        ]

    lineage = payload["primary_source_lineage"]
    assert lineage["Allcock_2006"]["Figure_4_1_H9_Q"][
        "facet_count_from_vector_figure"
    ] == 14
    assert lineage["Bredon_ADEG"][
        "all_five_outside_ADEG_because_each_uses_label_4"
    ]
    assert lineage["Tumarkin"]["none_of_five_literal_n_plus_2_by_facet_count"]
    assert lineage["Vinberg_Kaplinskaja_P9"][
        "A14_10_1_C01_scaled_dual_class_is_P9_entry7"
    ]
    assert payload["summary"]["new_arithmetic_commensurability_class_count"] == 0
    assert payload["summary"]["new_maximal_reflection_group_count"] == 0


def test_conditional_k3_interpretation_is_classical_and_scoped() -> None:
    interpretation = _archive()["conditional_algebro_geometric_interpretation"]
    assert "do not construct" in interpretation["scope_caveat"]
    assert "no explicit polynomial normal form" in interpretation[
        "equation_family_status"
    ]
    assert all(interpretation["source_marker_checks"].values())
    expected = {
        "(10, 4, 0)": ("U + D4(-1)^2", "genus-4", 3),
        "(10, 6, 0)": ("U(2) + D4(-1)^2", "genus-3", 2),
        "(10, 8, 0)": ("U + E8(-2)", "two disjoint smooth elliptic curves", None),
    }
    for invariants, (lattice, locus, rational_count) in expected.items():
        row = interpretation["types"][invariants]
        assert row["lattice"] == lattice
        assert row["nonsymplectic_K3_moduli_dimension"] == 10
        assert locus in row["fixed_locus"]
        if rational_count is None:
            assert row["fixed_locus_exceptional_case"]
            assert "two smooth fibers" in row["quotient_model"]
        else:
            assert not row["fixed_locus_exceptional_case"]
            assert row["additional_rational_curve_count"] == rational_count
