import json

from scripts.certify_k3_rho16_h14_fold_identification import (
    CERTIFICATE_PATH,
    build_certificate,
    verify_certificate,
)


def test_rho16_h14_fold_identification_replays() -> None:
    archived = json.loads(CERTIFICATE_PATH.read_text())
    assert archived == build_certificate()
    assert all(verify_certificate(archived).values())
    assert archived["all_checks_pass"] is True


def test_proposed_action_is_kondos_rejected_graph_reflection() -> None:
    payload = build_certificate()
    action = payload["proposed_action"]
    assert action["weighted_graph_automorphism_group"] == "C2 x C2"
    assert action["weighted_graph_automorphism_count"] == 4
    assert action["fixed_wall_count"] == 15
    assert action["unique_nonidentity_action_fixing_15_walls"] is True
    lattice = payload["lattice_action"]
    assert lattice["discriminant_group"] == "(Z/2)^2"
    assert lattice["discriminant_action"] == [0, 2, 1, 3]
    assert lattice["discriminant_action_is_nontrivial"] is True
    assert payload["algebraic_realization"]["realized_K3_involution"] is False


def test_fold_is_exactly_published_vinberg_h14_chamber() -> None:
    payload = build_certificate()
    fold = payload["folded_diagram"]
    assert fold["hyperbolic_dimension"] == 14
    assert fold["facet_count"] == 17
    assert fold["orbit_gram_inertia"] == [1, 14, 2]
    assert fold["edge_label_histogram"] == {"3": 14, "4": 2, "infinity": 1}
    published = payload["published_identification"]
    assert published["name"] == "14-vinb85.coxiter"
    assert published["exact_edge_labeled_graph_match"] is True
    assert published["published_total_vertex_count"] == 94
    assert published["published_finite_vertex_count"] == 89
    assert published["published_ideal_vertex_count"] == 5
    assert published["published_edge_count"] == 704
    assert payload["conclusion"]["new_abstract_hyperbolic_coxeter_polytope"] is False
    assert payload["conclusion"]["new_algebro_geometric_quotient_example"] is False
