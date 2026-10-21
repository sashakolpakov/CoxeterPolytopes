import json
from pathlib import Path

from scripts.scout_roulleau_rank12_graph_folds import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json"
CERTIFICATE = json.loads(ARTIFACT.read_text())


def test_two_complete_weighted_subgroup_denominators() -> None:
    assert CERTIFICATE["all_checks_pass"] is True
    assert CERTIFICATE["summary"] == {
        "bounded_literal_diagram_priority_candidate_count": 2,
        "finite_volume_noncompact_type_count": 3,
        "general_surface_geometric_action_count": 0,
        "known_abstract_full_chamber_count": 1,
        "new_arithmetic_commensurability_class_count": 0,
        "proper_high_coxeter_survivor_class_count": 3,
        "stock_count": 2,
        "subgroup_conjugacy_class_count": 12,
    }
    first, second = CERTIFICATE["stocks"]
    assert (first["weighted_diagram_group"]["structure"], first["subgroup_class_count"]) == (
        "C2",
        2,
    )
    assert (
        second["weighted_diagram_group"]["structure"],
        second["subgroup_class_count"],
    ) == ("D12 = S3 x C2", 10)
    assert first["proper_high_coxeter_survivor_class_indices"] == [2]
    assert second["proper_high_coxeter_survivor_class_indices"] == [2, 3]


def test_three_exact_finite_volume_noncompact_catches() -> None:
    cases = [case for stock in CERTIFICATE["stocks"] for case in stock["survivors"]]
    assert [
        (
            case["geometry"]["hyperbolic_dimension"],
            case["geometry"]["facet_count"],
            case["geometry"]["finite_vertex_count"],
            case["geometry"]["ideal_vertex_count"],
            case["geometry"]["vinberg_volume_criterion"]["elliptic_edge_count"],
        )
        for case in cases
    ] == [
        (10, 12, 17, 2, 99),
        (10, 15, 42, 10, 294),
        (10, 16, 70, 6, 426),
    ]
    for case in cases:
        vinberg = case["geometry"]["vinberg_volume_criterion"]
        assert vinberg["decision"] == "finite_volume_noncompact"
        assert vinberg["finite_volume_bad_edges"] == []
        assert case["geometry"]["all_checks_pass"] is True
        assert case["fixed_lattice"]["all_reflections_integral"] is True


def test_full_reflection_lineage_and_priority_are_separated() -> None:
    cases = [case for stock in CERTIFICATE["stocks"] for case in stock["survivors"]]
    assert [case["maximal_reflection_lineage"]["reflection_subgroup_index"] for case in cases] == [
        1,
        1054,
        992,
    ]
    assert [case["fixed_lattice"]["determinant"] for case in cases] == [4, 16, 16]
    assert [len(case["fixed_lattice"]["maximal_stable_overlattices"]) for case in cases] == [
        1,
        3,
        1,
    ]
    assert cases[0]["priority_status"] == "known_abstract_full_chamber"
    assert all(
        case["priority_status"]
        == "literal_diagram_unmatched_in_pinned_corpora_worldwide_priority_open"
        and not case["coxiter_lineage"]["ordinary_exact_matches"]
        and not case["coxiter_lineage"]["doubling_exact_matches"]
        for case in cases[1:]
    )
    assert cases[1]["prior_project_ledger_comparison"] == {
        "candidate_artifact": "data/enumerations/roulleau_r14_ue8a1x4_all_subgroup_folds.json",
        "candidate_relation_matrix_sha256": "caddc8bb3c5592b38c6a015bd6eeba06b84dd38799ac097d0c38c6bf3dcc1a0d",
        "candidate_id": "R14-E8A1x4-V4-H10-15",
        "exact_relation_graph_isomorphic": False,
        "status": "distinct_literal_H10_15_diagram",
    }
    assert cases[2]["prior_project_ledger_comparison"] == {
        "same_dimension_and_facet_prior_ledger_candidate_count": 0
    }
    assert all(
        case["neron_severi_action"]["surface_realization_status"].startswith(
            "not_realized_on_the_general_K3"
        )
        for case in cases
    )


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == CERTIFICATE
