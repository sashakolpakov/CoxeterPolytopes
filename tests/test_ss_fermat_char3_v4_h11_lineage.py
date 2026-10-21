from __future__ import annotations

import json

from scripts.certify_ss_fermat_char3_v4_h11_lineage import OUTPUT


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_v4_h11_stable_lattice_is_sk_entry20() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    stable = payload["minimal_reflection_stable_lattice"]
    assert stable["index_over_divisor_fixed_lattice"] == 4
    assert stable["determinant_absolute"] == 48
    assert stable["smith_diagonal_absolute"] == [1] * 8 + [2, 2, 2, 6]

    maximal = payload["known_maximal_integral_reflection_chamber"]
    assert maximal["source"]["entry_one_based"] == 20
    assert maximal["source"]["genus_symbol"] == "II_{11,1}(2^{4} 3^{-1})"
    assert maximal["explicit_isometry_determinant"] == -1
    assert maximal["source_gram_sha256"] == (
        "c854fc82d4602c0e6cb3752d76e0a195e51194532754feaab85d2b92b99e455f"
    )
    assert maximal["facet_count"] == 22
    assert maximal["finite_vertex_count"] == 543
    assert maximal["cusp_count"] == 22
    assert maximal["full_face_vector_f0_through_polytope"] == [
        565,
        3564,
        9973,
        16851,
        19295,
        15660,
        9126,
        3789,
        1089,
        204,
        22,
        1,
    ]


def test_all_328_mirrors_lie_in_the_known_full_reflection_group() -> None:
    payload = _payload()
    subgroup = payload["induced_328_wall_reflection_subgroup"]
    assert subgroup["facet_count"] == 328
    assert subgroup["primitive_square_divisibility_histogram"] == {
        "-6,3": 276,
        "-4,2": 24,
        "-2,1": 28,
    }
    assert subgroup["all_reflections_integral_involutive_isometries"] is True
    mirrors = subgroup["mirrors_in_stable_and_sk_bases"]
    assert len(mirrors) == 328
    assert len({tuple(row["stable_lattice_coordinates"]) for row in mirrors}) == 328
    assert len({tuple(row["sk_entry20_coordinates"]) for row in mirrors}) == 328
    assert all(
        row["integral_involutive_isometry_in_both_bases"] for row in mirrors
    )
    assert "0/0" in subgroup["numeric_index_status"]


def test_near_miss_and_rational_lineage_are_separated() -> None:
    payload = _payload()
    near = payload["near_miss_decomposition"]
    assert near["candidate"] == "U + D8(-1) + A2(-2)"
    assert near["same_rank_determinant_and_smith_diagonal"] is True
    assert near["stable_discriminant_q_mod_2Z_histogram"] == {
        "0": 10,
        "1": 6,
        "1/3": 12,
        "4/3": 20,
    }
    assert near["candidate_discriminant_q_mod_2Z_histogram"] == {
        "0": 6,
        "1": 10,
        "2/3": 12,
        "5/3": 20,
    }
    assert near["integrally_isometric"] is False

    comparisons = payload["rational_and_census_comparisons"]
    assert comparisons["mcleod_f3_H11"]["exact_rational_comparison"][
        "rationally_similar"
    ] is True
    assert comparisons["alexeev_2_elementary_census"]["applicable"] is False


def test_bounded_corpus_and_lineage_decision() -> None:
    payload = _payload()
    bounded = payload["bounded_weighted_diagram_audit"]
    assert bounded["literal_328_facet_record_present"] is False
    snapshot = bounded["coxiter_snapshot"]
    assert snapshot["ordinary_H11_facet_histogram"] == {
        "12": 14,
        "13": 6,
        "15": 2,
    }
    assert snapshot["doubling_H11_facet_histogram"] == {
        "15": 2,
        "16": 2,
        "18": 2,
        "22": 2,
        "30": 2,
        "46": 2,
    }
    assert payload["lineage_decision"] == {
        "literal_328_wall_subgroup_chamber_priority": (
            "open after bounded corpus exclusion"
        ),
        "new_arithmetic_commensurability_class": False,
        "new_maximal_reflection_group": False,
        "new_stable_integral_reflection_lattice": False,
        "positive_algebraic_construction": True,
    }
