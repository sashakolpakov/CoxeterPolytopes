import json

from scripts.certify_k3_high_rank_oriented_monodromy import (
    EXPECTED_DISCRIMINANT_ORTHOGONAL_ORDERS,
    OUTPUT_PATH,
    build_artifact,
)


def _artifact() -> dict:
    return json.loads(OUTPUT_PATH.read_text())


def test_frozen_oriented_period_group_manifest_rebuilds_exactly() -> None:
    frozen = _artifact()
    source_check = frozen["source_artifacts"]["source_check"]
    assert source_check["all_checks_pass"] is True
    assert build_artifact(source_check) == frozen


def test_all_nineteen_full_O_embedding_orbits_remain_single_O_plus_orbits() -> None:
    frozen = _artifact()
    assert len(frozen["rows"]) == 19
    reflection = frozen["u_positive_reflection"]
    assert reflection["h_square"] == 2
    assert reflection["reflection_matrix"] == [[0, -1], [-1, 0]]
    assert reflection["reflection_determinant"] == -1
    for row in frozen["rows"]:
        assert row["N_splits_U"] is True
        assert row["T_splits_U"] is True
        assert row["primitive_embedding_orbits"]["under_full_O_Lambda_K3"] == 1
        assert row["primitive_embedding_orbits"]["under_O_plus_Lambda_K3"] == 1


def test_positive_cone_and_complete_chamber_orbits_keep_marking_conventions() -> None:
    frozen = _artifact()
    for row in frozen["rows"]:
        cone = row["positive_cone_components"]
        assert cone == {
            "number": 2,
            "orbits_under_O_N": 1,
            "orbits_under_O_plus_N": 2,
            "stabilizer_of_either_component": "O^+(N)",
        }
        chamber = row["complete_root_chambers_in_one_positive_component"]
        assert chamber["raw_cardinality"] == "countably infinite W_2(N)-torsor"
        assert chamber["W_2_action"] == "simply_transitive"
        assert chamber["orbits_under_W_2_N"] == 1
        assert chamber["orbits_under_O_plus_N"] == 1
        assert chamber["W_2_discriminant_action"] == "trivial"


def test_pointwise_and_unmarked_period_groups_have_exact_distinct_images() -> None:
    frozen = _artifact()
    assert frozen["schema"].endswith("oriented_period_groups.v2")
    actual_orders = {}
    for row in frozen["rows"]:
        key = row["key"]
        arithmetic = row["oriented_period_change_of_marking_data"]
        assert arithmetic["pointwise_N_polarized_T_period_image"].startswith(
            "tilde O^+(T)"
        )
        assert arithmetic["pointwise_discriminant_image"] == "identity"
        assert arithmetic["unmarked_ambient_chamber_normalizer"] == "Gamma^+_{j,C}"
        assert arithmetic["unmarked_chamber_normalizer_T_period_image"] == "O^+(T)"
        assert arithmetic["unmarked_T_discriminant_image"] == "all O(q_T)"
        assert arithmetic["geometric_monodromy_identification"].startswith(
            "conditional_only"
        )
        order = row["discriminant_form_orthogonal_group_order"]
        assert arithmetic["arithmetic_group_index"] == order
        assert arithmetic["coarse_period_quotient_degree"].startswith("not asserted")
        assert row["discriminant_representations"]["stable_index_T"] == order
        actual_orders[key] = order
    assert actual_orders == EXPECTED_DISCRIMINANT_ORTHOGONAL_ORDERS
    assert frozen["finite_discriminant_quotients"]["maximum_order"] == 40320


def test_no_group_conflation_or_unsupported_finite_presentations() -> None:
    frozen = _artifact()
    boundary = frozen["presentation_boundary"]
    assert boundary["explicit_finite_presentations_supplied"] is False
    assert "do not derive" in boundary["reason"]
    assert len(boundary["not_conflated"]) == 6
    for row in frozen["rows"]:
        assert "comparison-only" in row["separation_warning"]
        assert row["generic_surface_automorphism_group_source_label"]


def test_q_form_provenance_literal_V_and_strict_residual_are_honest() -> None:
    frozen = _artifact()
    sources = frozen["source_artifacts"]
    assert sources["dependency_graph"].startswith("acyclic upstream-only")
    assert "classification_manifest" not in sources
    assert sources["roulleau_19_stock_inventory"].endswith(
        "roulleau_rank11_14_source_stock_inventory.json"
    )
    assert sources["rank15_19_chamber_manifest"].endswith(
        "k3_high_rank_nef_chambers.json"
    )
    assert sources["rank15_19_literal_vinberg_manifest"].endswith(
        "k3_high_rank_nef_vinberg_volume.json"
    )

    provenance = frozen["finite_quadratic_form_provenance"]
    assert provenance["ADE_q_forms"].startswith("table-derived")
    assert "independently recomputed" in provenance["smith_invariants"]
    assert "not a second inverse-Gram derivation" in provenance["orthogonal_group_orders"]
    assert all("not independently derived" in row["q_form_provenance"] for row in frozen["rows"])

    nikulin = sources["nikulin_surjectivity_input"]
    assert nikulin["provenance"] == "bibliographic_primary_web_input_not_locally_hash_pinned"
    assert nikulin["primary_web"].startswith("https://www.mathnet.ru/")

    literal = frozen["literal_vinberg_status"]
    assert literal["all_19_closed"] is True
    assert literal["remaining_literal_V_gap_count"] == 0
    assert literal["ranks_11_14_case_count"] + literal["ranks_15_19_case_count"] == 19

    residual = frozen["strict_M3c_residual"]
    assert residual["status"] == "open"
    assert "closed for all 19" in residual["literal_V_status"]
    assert len(residual["reasons"]) == 3
