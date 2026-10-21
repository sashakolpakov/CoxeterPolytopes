from __future__ import annotations

import json

from scripts.certify_kondo_rho16_h14_bh_order4_lift import (
    FOLD,
    OUTPUT,
    SEARCH_CENSUS,
    SOURCE,
    _archive_pair,
    build_certificate,
)


CERTIFICATE = json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_order4_lift_certificate_replays_offline() -> None:
    assert build_certificate() == CERTIFICATE
    assert CERTIFICATE["all_checks_pass"] is True
    assert CERTIFICATE["summary"] == {
        "automorphism_order": 4,
        "compact": 0,
        "effective_database_rows": 2,
        "finite_volume_noncompact": 1,
        "full_group_enrichments": 1,
        "genuine_K3_automorphism_realizations": 2,
        "hyperbolic_dimension": 14,
        "new_abstract_diagram_claims": 0,
        "picard_action_order": 2,
    }


def test_order_two_no_go_is_retained_and_higher_order_escape_closes() -> None:
    boundary = CERTIFICATE["order_two_boundary_reconciliation"]
    assert "not induced by an order-two" in boundary["prior_theorem_retained"]
    assert "sigma has order four" in boundary["new_mechanism"]
    assert boundary["contradiction"] is False
    assert boundary["previous_special_higher_order_period_question_closed"]


def test_full_group_enrichment_and_pure_lineage_are_exact() -> None:
    bh = CERTIFICATE["brandhorst_hofmann"]
    pure = bh["pure_cyclic_row"]
    full = bh["full_group_enrichment"]
    assert (pure["row_id"], pure["group"], pure["global_automorphism_order"]) == (
        "0.4.2.2",
        "C4",
        4,
    )
    assert (
        full["row_id"],
        full["group"],
        full["selected_element"],
        full["selected_element_order"],
        full["picard_image"],
        full["kernel_of_group_to_O_NS"],
    ) == ("1.4.2.1", "C4 x C2", "s*g", 4, "C2 x C2", "<g^2> = C2")
    assert full["selected_cyclic_subgroup_is_unique_pure_row"] == "0.4.2.2"
    assert full["common_fixed_lattice_rank_signature_determinant"] == [
        8,
        [1, 7, 0],
        -256,
    ]


def test_archived_picard_pair_and_complete_graph_class_are_pinned() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    fold = json.loads(FOLD.read_text(encoding="utf-8"))
    pair = _archive_pair(source, fold)
    assert pair["fixed_gram"].rows == 15
    assert pair["fixed_gram"].det() == 4
    assert pair["coinvariant_gram"].tolist() == [[-4]]
    assert pair["coinvariant_action"].tolist() == [[-1]]
    assert pair["graph_fixed_rank_classes"] == {
        "rho16_involution_1": 15,
        "rho16_involution_2": 8,
        "rho16_involution_3": 9,
    }
    exact = CERTIFICATE["exact_picard_action_conjugacy"]
    assert exact["fixed_lattices"]["nikulin_lifting"].startswith("rank 15")
    assert exact["primitive_glue"]["index"] == 2
    assert exact["primitive_glue"]["identity_maps_carry_both_pairs"]


def test_vinberg_gate_is_finite_volume_not_signature_only() -> None:
    polytope = CERTIFICATE["hyperbolic_polytope"]
    assert (
        polytope["id"],
        polytope["dimension"],
        polytope["facets"],
        polytope["classification"],
    ) == (
        "R16-UE8D6-C2-H14-17",
        14,
        17,
        "finite_volume_noncompact",
    )
    assert (
        polytope["finite_vertices"],
        polytope["cusps"],
        polytope["hyperideal_vertices"],
        polytope["elliptic_edges"],
        polytope["bad_edges"],
        polytope["parabolic_neighborhood_audits"],
    ) == (89, 5, 0, 704, 0, 8)
    assert polytope["complete_simple_roots"]


def test_finite_search_census_and_novelty_scope_are_precise() -> None:
    assert CERTIFICATE["finite_database_search_census"] == SEARCH_CENSUS
    assert SEARCH_CENSUS["pure_exact_rho16_target_hits"] == ["0.4.2.2"]
    assert SEARCH_CENSUS["full_non_order2_exact_rho16_target_hits"] == [
        "1.4.2.1"
    ]
    assert SEARCH_CENSUS["full_order2_exact_rho16_target_hits"] == []
    scope = CERTIFICATE["novelty_scope"]
    assert scope["new_project_level_genuine_K3_complete_nef_realization"]
    assert scope["new_project_level_full_group_enrichment"]
    assert scope["known_polytope_new_realization"]
    assert not scope["new_abstract_coxeter_polytope_claimed"]
    assert not scope["new_dimension_claimed"]
    assert not scope["worldwide_priority_claimed"]
