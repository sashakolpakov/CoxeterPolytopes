import json
from pathlib import Path

from scripts.certify_k3_high_rank_full_root_classification import (
    EXPECTED_CASES,
    EXPECTED_RANK_COUNTS,
    OUTPUT_PATH,
    build_artifact,
)


ROOT = Path(__file__).resolve().parents[1]


def _artifact() -> dict:
    return json.loads(OUTPUT_PATH.read_text())


def test_frozen_synthesis_rebuilds_exactly() -> None:
    frozen = _artifact()
    assert build_artifact(frozen["source_check"]) == frozen


def test_ambient_class_does_not_assume_volume_and_fv_locus_has_nineteen_types() -> None:
    frozen = _artifact()
    classification = frozen["classification"]
    assert frozen["schema"].endswith(".v3")
    assert frozen["surface_free_class"]["name"] == "K3Root[11,20]"
    assert frozen["surface_free_class"]["surface_language_used_in_definition"] is False
    assert frozen["surface_free_class"]["finite_volume_or_finite_index_used_in_definition"] is False
    assert all(
        "finite index" not in condition and "finite volume" not in condition
        for condition in frozen["surface_free_class"]["definition"]
    )
    assert frozen["finite_volume_sublocus"]["name"] == "K3RootFV[11,20]"
    assert "not an ambient hypothesis" in frozen["finite_volume_sublocus"]["finite_index_equivalence"]
    assert classification["classified_locus"] == "finite_volume_sublocus_only"
    assert classification["ambient_class_exhausted"] is False
    assert classification["finite_volume_rank_counts"] == {
        str(rank): count for rank, count in EXPECTED_RANK_COUNTS.items()
    }
    assert classification["rank_11_through_19_finite_volume_type_count"] == 19
    assert classification["rank_20_finite_volume_type_count"] == 0
    assert classification["finite_volume_hyperbolic_dimension_range_realized"] == [10, 18]
    assert classification["maximum_finite_volume_hyperbolic_dimension"] == 18
    assert len(frozen["cases"]) == len(EXPECTED_CASES) == 19


def test_cases_are_pairwise_distinct_and_have_unique_k3_embeddings() -> None:
    frozen = _artifact()
    distinction_keys = []
    for case in frozen["cases"]:
        distinction_keys.append(tuple(case["pairwise_distinction_key"]))
        inequality = case["nikulin_stable_unique_embedding_inequality"]
        assert inequality["passes"] is True
        assert inequality["left"] <= inequality["right"]
        assert case["exact_picard_moduli_dimension"] == 20 - case["picard_rank"]
        assert case["contains_U"] is True
        assert case["complete_root_chamber"] is True
    assert len(set(distinction_keys)) == 19
    assert sum(case["facet_count"] for case in frozen["cases"]) == 477


def test_u_r_c_v_and_volume_claim_boundaries_are_explicit() -> None:
    frozen = _artifact()
    theorem = frozen["transfer_theorem"]
    assert theorem["U"]["status"] == "closed"
    assert theorem["R"]["status"].startswith("closed")
    assert theorem["C"]["status"] == "closed"
    assert theorem["V"]["status"] == "closed_literally_and_structurally_for_all_19"
    assert theorem["V"]["literal_vinberg_replay_count"] == 19
    assert theorem["V"]["structural_only_at_synthesis_level_count"] == 0
    assert theorem["V"]["noncompact_all_19"] is True
    assert "not assumed in the ambient class" in theorem["V"]["proof"]
    assert theorem["V"]["literal_audit_totals"] == {
        "bad_full_parabolic_neighborhood_count": 0,
        "connected_affine_audit_count": 4255,
        "elliptic_edge_count": 117845,
        "facet_count": 477,
        "finite_vertex_count": 14098,
        "finite_volume_bad_edge_count": 0,
        "ideal_vertex_count": 1156,
    }
    assert theorem["strict_m3c"]["status"] == "open"
    assert len(theorem["strict_m3c"]["reasons"]) == 3
    assert any("IC0" in reason for reason in theorem["strict_m3c"]["reasons"])

    literal = [
        case for case in frozen["cases"]
        if case["literal_vinberg_status"]["available"]
    ]
    structural_only = [
        case for case in frozen["cases"]
        if not case["literal_vinberg_status"]["available"]
    ]
    assert [case["picard_rank"] for case in literal] == [
        11, 11, 11, 11,
        12, 12, 12, 12,
        13, 13, 13,
        14, 14, 14,
        15, 16, 17, 18, 19,
    ]
    assert structural_only == []
    assert all(
        case["literal_vinberg_status"]["decision"]
        == "finite_volume_noncompact"
        for case in literal
    )
    assert all(
        case["exact_gram_replay_in_this_synthesis"]["available"]
        for case in frozen["cases"]
    )


def test_oriented_arithmetic_period_layer_is_exact_without_group_conflation() -> None:
    frozen = _artifact()
    oriented = frozen["oriented_arithmetic_period_classification"]
    assert oriented["status"] == "closed_for_all_19"
    assert oriented["primitive_embedding_orbits_under_O_plus_per_case"] == 1
    assert (
        oriented[
            "complete_chamber_orbits_in_a_fixed_component_under_O_plus_per_case"
        ]
        == 1
    )
    assert oriented["pointwise_N_polarized_T_period_image"] == "tilde O^+(T)"
    assert oriented["unmarked_ambient_chamber_normalizer"] == "Gamma^+_{j,C}"
    assert oriented["unmarked_chamber_normalizer_T_period_image"] == "O^+(T)"
    assert "not an identification" in oriented["normalizer_warning"]
    assert "not independently certified" in oriented["geometric_monodromy_status"]
    assert oriented["finite_discriminant_quotients"]["maximum_order"] == 40320
    for case in frozen["cases"]:
        arithmetic = case["oriented_arithmetic_period_data"]
        assert arithmetic["primitive_embedding_orbits_under_O_plus"] == 1
        assert arithmetic[
            "complete_chamber_orbits_in_one_component_under_O_plus"
        ] == 1
        assert arithmetic["unmarked_ambient_chamber_normalizer"] == "Gamma^+_{j,C}"
        assert arithmetic["unmarked_chamber_normalizer_T_period_image"] == "O^+(T)"


def test_rank11_22_case_uses_the_independent_standard_lattice_reconstruction() -> None:
    frozen = _artifact()
    case = next(case for case in frozen["cases"] if case["key"] == "R11_UD4D4A1")
    gram = case["exact_gram_replay_in_this_synthesis"]
    literal = case["literal_vinberg_status"]
    assert gram["matrix_provenance"]["kind"] == (
        "independent_standard_lattice_norm_two_vinberg_reconstruction"
    )
    assert gram["expected_signature"] == [1, 10, 11]
    assert literal["elliptic_edge_count"] == 1343
    assert literal["finite_vertex_count"] == 185
    assert literal["ideal_vertex_count"] == 25
    assert literal["connected_affine_audit_count"] == 72
    assert literal["finite_volume_bad_edge_count"] == 0
    assert literal["bad_full_parabolic_neighborhood_count"] == 0


def test_new_frozen_package_hashes_are_pinned_in_the_synthesis() -> None:
    sources = _artifact()["source_artifacts"]
    assert sources["rank11_14_literal_vinberg"]["sha256"] == (
        "f97b713d76fbbd7e18d93c018b7055164da7768f9f21355f337e8aff085a81cb"
    )
    assert sources["oriented_arithmetic_period_groups"]["sha256"] == (
        "c9d886a568abd636c921f4f8c478a35f4d3f5a94b1adbaa3f44493cabb54b0a3"
    )


def test_moduli_dimensions_and_novelty_are_not_overclaimed() -> None:
    frozen = _artifact()
    assert frozen["moduli"]["irreducible_family_count"] == 19
    assert "full coarse N-polarized" in frozen["moduli"]["irreducibility_scope"]
    assert "exact-Picard locus is asserted only nonempty" in frozen["moduli"]["irreducibility_scope"]
    assert frozen["moduli"]["dimensions_by_rank"] == {
        str(rank): 20 - rank for rank in range(11, 20)
    }
    assert "arithmetic period/change-of-marking" in frozen["moduli"][
        "orientation_and_monodromy_boundary"
    ]
    assert frozen["moduli"][
        "primitive_embedding_orbits_under_oriented_k3_lattice_isometry"
    ] == 19
    assert frozen["classification"]["worldwide_novelty"] is False
    assert all(
        case["novelty"] == "published_known_full_nef_Coxeter_polytope"
        for case in frozen["cases"]
    )
