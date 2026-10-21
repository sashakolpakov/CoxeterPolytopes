import json

from scripts.certify_k3_high_rank_all_involution_folds import (
    CERTIFICATE_PATH,
    build_certificate,
    verify_certificate,
)


def _archived() -> dict:
    return json.loads(CERTIFICATE_PATH.read_text())


def test_exhaustive_k3_involution_fold_certificate_replays() -> None:
    archived = _archived()
    assert archived == build_certificate()
    assert all(verify_certificate(archived).values())
    assert archived["all_checks_pass"] is True
    assert archived["totals"] == {
        "weighted_graph_automorphisms": 8,
        "nonidentity_involutions": 4,
        "direct_finite_volume_fold_certificates": 4,
        "dimension_gt_9_target_actions": [
            "rho19_involution_1",
            "rho19_involution_2",
            "rho19_involution_3",
        ],
        "dimension_gt_9_geometric_conjugacy_classes": 1,
        "certified_new_folds": 0,
    }


def test_weighted_graph_involution_census_is_complete() -> None:
    cases = {case["ambient_rho"]: case for case in _archived()["cases"]}
    assert cases[18]["weighted_graph_automorphism_group"] == "C2"
    assert cases[18]["automorphism_count"] == 2
    assert cases[18]["element_order_histogram"] == {"1": 1, "2": 1}
    assert cases[18]["nonidentity_involution_count"] == 1
    assert cases[18]["nonidentity_involution_conjugacy_class_count"] == 1

    assert cases[19]["weighted_graph_automorphism_group"] == "S3"
    assert cases[19]["automorphism_count"] == 6
    assert cases[19]["element_order_histogram"] == {"1": 1, "2": 3, "3": 2}
    assert cases[19]["nonidentity_involution_count"] == 3
    assert cases[19]["nonidentity_involution_conjugacy_class_count"] == 1
    assert all(case["all_checks_pass"] for case in cases.values())


def test_every_fold_passes_exact_coxeter_reflection_and_volume_gates() -> None:
    folds = [fold for case in _archived()["cases"] for fold in case["folds"]]
    assert len(folds) == 4
    assert all(fold["finite_volume"] is True for fold in folds)
    assert all(fold["compact"] is False for fold in folds)
    assert all(
        fold["gate_checks"]["all_facet_pairs_are_exact_coxeter_relations"]
        and fold["gate_checks"]["all_folded_facet_reflections_integral"]
        and fold["gate_checks"]["all_retained_walls_irredundant"]
        and fold["gate_checks"]["vinberg_coxiter_finite_volume"]
        and fold["gate_checks"]["direct_full_audit_passes"]
        for fold in folds
    )
    assert all(
        row["role"] in {
            "disjoint_reflective_facet",
            "redundant_timelike_inequality",
        }
        for fold in folds
        for row in fold["two_curve_orbit_audit"]
    )

    rank18 = [fold for fold in folds if fold["ambient_rho"] == 18]
    rank19 = [fold for fold in folds if fold["ambient_rho"] == 19]
    assert len(rank18) == 1 and rank18[0]["hyperbolic_dimension"] == 9
    assert rank18[0]["passes_dimension_gt_9_target"] is False
    assert rank18[0]["finite_vertex_count"] == 8
    assert rank18[0]["ideal_vertex_count"] == 2
    assert rank18[0]["literal_coxiter_prop_1_10"]["extension_histogram"] == {
        "2": 45
    }

    assert len(rank19) == 3
    assert all(fold["hyperbolic_dimension"] == 10 for fold in rank19)
    assert all(fold["passes_dimension_gt_9_target"] is True for fold in rank19)
    assert all(fold["finite_vertex_count"] == 34 for fold in rank19)
    assert all(fold["ideal_vertex_count"] == 4 for fold in rank19)
    assert all(
        fold["literal_coxiter_prop_1_10"]["extension_histogram"] == {"2": 208}
        for fold in rank19
    )
    assert {
        tuple(fold["redundant_timelike_orbits"][0]) for fold in rank19
    } == {(22, 23), (22, 24), (23, 24)}


def test_realization_and_novelty_are_not_inferred_from_graph_symmetry() -> None:
    archived = _archived()
    folds = [fold for case in archived["cases"] for fold in case["folds"]]
    assert all(
        fold["algebraic_realization"]["status"]
        == "published_surface_automorphism"
        for fold in folds
    )
    assert sum(fold["equals_archived_published_action"] for fold in folds) == 2
    assert all(fold["gate_checks"]["exactly_conjugate_to_published_action"] for fold in folds)
    assert all(fold["gate_checks"]["orbit_gram_transport_matches"] for fold in folds)
    assert all(fold["certified_new_fold"] is False for fold in folds)
    assert archived["new_hyperbolic_example"] is False
