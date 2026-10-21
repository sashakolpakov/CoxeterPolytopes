from __future__ import annotations

import json

from scripts.certify_shimada_disc16_d5_h10_lineage import OUTPUT


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_saturated_fixed_lattice_is_twice_odd_unimodular() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    lattice = payload["saturated_fixed_lattice"]
    assert lattice["rank"] == 11
    assert lattice["inertia"] == [1, 10, 0]
    assert lattice["determinant"] == 2**11
    assert lattice["smith_diagonal"] == [2] * 11
    assert lattice["target_lattice"] == "2*I_{1,10}"
    assert lattice["half_scaled_determinant"] == 1
    assert lattice["half_scaled_odd_unimodular"] is True
    assert lattice["explicit_isometry_determinant"] == -1
    assert lattice["integrally_isometric"] is True


def test_all_thirty_mirrors_are_integral_and_span_the_lattice() -> None:
    mirrors = _payload()["primitive_mirrors"]
    assert mirrors["facet_count"] == 30
    assert mirrors["fixed_lattice_square_divisibility_histogram"] == {
        "-4,2": 12,
        "-2,2": 18,
    }
    assert mirrors["half_scaled_square_divisibility_histogram"] == {
        "-2,1": 12,
        "-1,1": 18,
    }
    assert mirrors["all_reflections_integral_involutive_isometries"] is True
    assert mirrors["target_root_matrix_rank"] == 11
    assert mirrors["target_root_span_smith_diagonal"] == [1] * 11
    assert mirrors["roots_generate_entire_fixed_lattice"] is True
    assert len(mirrors["records"]) == 30
    assert all(
        record["reflection_integral_involutive_isometry"]
        for record in mirrors["records"]
    )


def test_maximal_chamber_matches_vinberg_alexeev_sk_and_coxiter() -> None:
    maximal = _payload()["known_maximal_reflection_chamber"]
    assert maximal["name"] == "Vinberg P10"
    assert maximal["facet_count"] == 12
    assert maximal["finite_vertex_count"] == 17
    assert maximal["ideal_vertex_count"] == 2
    assert maximal["primitive_square_divisibility_histogram"] == {
        "-4,2": 10,
        "-2,2": 2,
    }
    assert maximal["all_reflections_integral"] is True
    assert maximal["coxiter"]["p10_root_to_coxiter_zero_based_permutation"] == list(
        range(12)
    )
    assert maximal["coxiter"]["entrywise_relation_matrix_match"] is True
    assert maximal["alexeev_11_11_1"][
        "p10_root_to_alexeev_zero_based_permutation"
    ] == [2, 3, 4, 5, 6, 7, 8, 10, 0, 11, 1, 9]
    assert maximal["alexeev_11_11_1"]["entrywise_weighted_gram_match"] is True
    assert maximal["scharlau_kirschmer_list10_entry1"]["genus_symbol"] == (
        "I_{10,1}()"
    )
    assert maximal["scharlau_kirschmer_list10_entry1"]["cusp_count"] == 2


def test_solomon_index_volume_and_face_vector() -> None:
    solomon = _payload()["solomon_index"]
    assert solomon["candidate_30_wall_group"][
        "rational_euler_characteristic"
    ] == "-527/483840"
    assert solomon["vinberg_p10_group"]["rational_euler_characteristic"] == (
        "-1/11147673600"
    )
    assert solomon["index"] == 12_142_080
    assert solomon["index_prime_factorization"] == {
        "2": 9,
        "3": 2,
        "5": 1,
        "17": 1,
        "31": 1,
    }
    assert solomon["hyperbolic_covolume"] == "pi^5 * 527/14288400"
    assert solomon["full_face_vector"] == [
        1_540,
        10_030,
        25_206,
        34_857,
        30_452,
        17_784,
        7_072,
        1_890,
        321,
        30,
        1,
    ]


def test_bounded_corpus_exclusion_is_not_a_priority_claim() -> None:
    payload = _payload()
    bounded = payload["bounded_weighted_diagram_audit"]
    snapshot = bounded["coxiter_snapshot"]
    assert snapshot["ordinary_file_count"] == 1_214
    assert snapshot["ordinary_H10_facet_histogram"] == {
        "11": 14,
        "12": 8,
        "13": 4,
        "14": 1,
    }
    assert snapshot["ordinary_H10_30_candidates"] == []
    assert snapshot["ordinary_exact_weighted_diagram_matches"] == []
    assert snapshot["doubling_file_count"] == 2_959
    assert snapshot["doubling_H10_facet_histogram"] == {
        "12": 3,
        "13": 4,
        "14": 3,
        "15": 1,
        "16": 3,
        "19": 1,
        "20": 3,
        "28": 2,
        "44": 2,
    }
    assert snapshot["doubling_H10_30_candidates"] == []
    assert snapshot["doubling_exact_weighted_diagram_matches"] == []
    assert bounded["does_not_prove_worldwide_absence"] is True
    assert payload["lineage_decision"]["worldwide_priority_claimed"] is False

