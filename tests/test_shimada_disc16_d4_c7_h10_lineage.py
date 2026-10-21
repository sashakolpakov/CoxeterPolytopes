from __future__ import annotations

import json

from scripts.certify_shimada_disc16_d4_c7_h10_lineage import OUTPUT


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_fixed_lattice_and_root_span_are_not_conflated_with_overlattice() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    assert payload["candidate"][
        "source_chamber_archived_group_equals_full_weighted_wall_group"
    ] is True
    fixed = payload["saturated_divisor_fixed_lattice"]
    assert fixed["rank"] == 11
    assert fixed["inertia"] == [1, 10, 0]
    assert fixed["determinant"] == 512
    assert fixed["smith_diagonal"] == [1, 1, 1, 1, 2, 2, 2, 2, 2, 4, 4]
    assert fixed["discriminant_group"] == "(Z/2)^5 + (Z/4)^2"
    mirrors = payload["primitive_candidate_mirrors"]
    assert mirrors["facet_count"] == 27
    assert mirrors["fixed_lattice_square_divisibility_histogram"] == {
        "-8,4": 3,
        "-4,2": 8,
        "-2,1": 6,
        "-2,2": 10,
    }
    assert mirrors["fixed_root_matrix_rank"] == 11
    assert mirrors["fixed_root_span_smith_diagonal"] == [1] * 11
    assert mirrors["roots_generate_entire_fixed_lattice"] is True


def test_discriminant_exhaustion_gives_unique_index_four_stable_overlattice() -> None:
    stable = _payload()["unique_reflection_stable_integral_overlattice"]
    audit = stable["discriminant_exhaustion"]
    assert audit["discriminant_group_order"] == 512
    assert audit["self_integral_coset_count"] == 256
    assert audit["invariant_integral_subgroup_count"] == 2
    assert audit["invariant_integral_subgroup_orders"] == [1, 4]
    assert audit["only_trivial_and_unique_order_four_subgroups_survive"] is True
    assert audit["unique_nontrivial_subgroup_even"] is True
    assert stable["index_over_L"] == 4
    assert stable["determinant"] == 32
    assert stable["smith_diagonal"] == [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2]
    assert stable["two_elementary_invariants"] == {"r": 11, "a": 5, "delta": 1}
    assert stable["delta_one_witness_square"] == "-1/2"
    assert stable["standard_lattice"] == "U + D6(-1) + A1(-1)^3"
    assert stable["explicit_isometry_determinant"] == 1
    assert stable["all_27_reflections_integral"] is True
    assert stable["candidate_root_content_histogram"] == {"1": 24, "2": 3}


def test_full_reflection_chamber_matches_scharlau_kirschmer_entry20() -> None:
    maximal = _payload()["known_full_reflection_chamber"]
    assert maximal["scharlau_kirschmer_list10_entry"] == 20
    assert maximal["scharlau_kirschmer_record"] == {
        "determinant": 32,
        "genus_symbol": "II_{10,1}([2^{5}]_{1})",
        "finite_face_counts_facets_to_vertices": [
            14,
            89,
            340,
            868,
            1_554,
            1_986,
            1_800,
            1_114,
            424,
            70,
        ],
        "cusp_count": 5,
    }
    assert maximal["height_zero_type"] == "B6 + B3 + A1"
    assert maximal["facet_count"] == 14
    assert maximal["finite_vertex_count"] == 70
    assert maximal["ideal_vertex_count"] == 5
    assert maximal["hyperideal_ray_count"] == 0
    assert maximal["primitive_square_divisibility_histogram"] == {
        "-4,2": 4,
        "-2,1": 8,
        "-2,2": 2,
    }
    assert maximal["relation_histogram"] == {"0": 2, "2": 76, "3": 9, "4": 4}
    assert maximal["full_integral_root_reflection_group"] is True


def test_solomon_index_volume_and_both_face_vectors() -> None:
    solomon = _payload()["solomon_index"]
    assert solomon["candidate_27_wall_group"][
        "rational_euler_characteristic"
    ] == "-527/1658880"
    assert solomon["full_14_wall_group"][
        "rational_euler_characteristic"
    ] == "-527/3344302080"
    assert solomon["index"] == 2_016
    assert solomon["index_prime_factorization"] == {"2": 5, "3": 2, "7": 1}
    assert solomon["candidate_hyperbolic_covolume"] == "pi^5 * 527/48988800"
    assert solomon["candidate_full_face_vector"] == [
        962,
        6_221,
        15_784,
        22_234,
        19_934,
        12_061,
        5_034,
        1_435,
        264,
        27,
        1,
    ]
    assert solomon["maximal_full_face_vector"] == [
        75,
        424,
        1_114,
        1_800,
        1_986,
        1_554,
        868,
        340,
        89,
        14,
        1,
    ]


def test_bounded_corpus_exclusion_does_not_become_worldwide_priority() -> None:
    payload = _payload()
    snapshot = payload["bounded_weighted_diagram_audit"]["coxiter_snapshot"]
    assert snapshot["ordinary_file_count"] == 1_214
    assert snapshot["ordinary_H10_27_candidates"] == []
    assert snapshot["doubling_file_count"] == 2_959
    assert snapshot["doubling_H10_27_candidates"] == []
    assert payload["bounded_weighted_diagram_audit"][
        "does_not_prove_worldwide_absence"
    ] is True
    decision = payload["lineage_decision"]
    assert decision["new_positive_algebraic_extraction"] is True
    assert decision["new_reflection_stable_lattice"] is False
    assert decision["new_maximal_reflection_group"] is False
    assert decision["new_arithmetic_commensurability_class"] is False
    assert decision["worldwide_priority_claimed"] is False
