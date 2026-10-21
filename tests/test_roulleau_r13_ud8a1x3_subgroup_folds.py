import json
from pathlib import Path

from scripts.certify_roulleau_r13_ud8a1x3_subgroup_folds import (
    _coordinate_model,
    _first_nineteen_gram,
    build_certificate,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/roulleau_r13_ud8a1x3_all_subgroup_folds.json"
CERTIFICATE = json.loads(ARTIFACT.read_text())


def test_coordinate_formula_source_reconstruction() -> None:
    basis_gram, walls, wall_gram, labels = _coordinate_model()
    source = CERTIFICATE["source_reconstruction"]
    assert basis_gram.det() == 32
    assert exact_inertia(basis_gram) == (1, 12, 0)
    assert walls.shape == (13, 39)
    assert all(value.is_Integer for value in walls)
    assert wall_gram.rank() == 13
    assert exact_inertia(wall_gram) == (1, 12, 26)
    assert wall_gram[:19, :19] == _first_nineteen_gram()
    assert all(wall_gram[index, index] == -2 for index in range(39))
    assert labels == [
        (1, 2), (1, 3), (1, 4), (1, 5), (2, 3),
        (2, 4), (2, 5), (3, 4), (3, 5), (4, 5),
    ]
    assert "exactly consistent" in source["printed_graph_status"]
    assert "AJ/BJ formulas" in source["printed_graph_status"]
    assert all(source["checks"].values())


def test_complete_s5_times_c2_subgroup_denominator() -> None:
    assert CERTIFICATE["schema"] == "coxeter.roulleau_r13_ud8a1x3.all_subgroup_folds.v1"
    assert CERTIFICATE["all_checks_pass"] is True
    group = CERTIFICATE["weighted_wall_group"]
    census = CERTIFICATE["subgroup_census"]
    assert group["order"] == 240
    assert group["recognition"] == "S5 x C2"
    assert group["id_group"] == [240, 189]
    assert group["central_C2_fixed_rank"] == 5
    assert census["conjugacy_class_count"] == 57
    assert census["proper_fixed_rank_at_least_11_count"] == 5
    high = [
        (row["class_index_one_based"], row["subgroup_order"], row["fixed_rank"], row["status"])
        for row in census["classes"]
        if row["fixed_rank"] >= 11 and row["subgroup_order"] > 1
    ]
    assert high == [
        (4, 2, 12, "exact_coxeter_candidate"),
        (6, 2, 11, "isotropic_orbit_normal_reject"),
        (7, 3, 11, "exact_coxeter_candidate"),
        (12, 4, 11, "isotropic_orbit_normal_reject"),
        (22, 6, 11, "exact_coxeter_candidate"),
    ]


def test_two_literal_vinberg_positive_types() -> None:
    expected = [
        (11, 31, 3948, 552, 41, 121),
        (10, 23, 1086, 170, 17, 42),
    ]
    observed = []
    for case in CERTIFICATE["positive_types"]:
        vinberg = case["geometry"]["vinberg"]
        observed.append(
            (
                case["hyperbolic_dimension"],
                case["facet_count"],
                vinberg["elliptic_edge_count"],
                vinberg["finite_vertex_count"],
                vinberg["ideal_vertex_count"],
                vinberg["connected_affine_component_count"],
            )
        )
        assert vinberg["finite_volume"] is True
        assert vinberg["compact"] is False
        assert vinberg["bad_edge_count"] == 0
        assert vinberg["bad_parabolic_neighborhood_count"] == 0
        assert all(vinberg["checks"].values())
        assert case["folded_walls"]["strict_facet_witnesses"]["all_orbit_inequalities_are_genuine_facets"]
    assert observed == expected


def test_c3_and_s3_are_the_same_integral_h10_chamber() -> None:
    identification = CERTIFICATE["C3_S3_integral_chamber_identification"]
    assert identification["source_class"] == 7
    assert identification["target_class"] == 22
    assert identification["determinant"] == 1
    assert identification["facet_permutation_zero_based"] == [
        1, 0, 2, 3, 4, 5, 6, 7, 8, 9, 11, 10, 12, 13, 14,
        17, 18, 21, 22, 19, 20, 15, 16,
    ]
    assert all(identification["checks"].values())


def test_stable_overclosures_and_bounded_lineage() -> None:
    h11, h10 = CERTIFICATE["positive_types"]
    assert h11["fixed_lattice"]["determinant"] == -32
    assert h10["fixed_lattice"]["determinant"] == 24
    h11_maxima = h11["stable_overlattices"]["maximal_reflection_stable_overclosures"]
    h10_maxima = h10["stable_overlattices"]["maximal_reflection_stable_overclosures"]
    assert [(row["index"], row["determinant"], row["even"]) for row in h11_maxima] == [
        (4, -2, False), (4, -2, False), (4, -2, False),
    ]
    assert [(row["index"], row["determinant"], row["even"]) for row in h10_maxima] == [
        (2, 6, False), (2, 6, True), (2, 6, True),
    ]
    assert [row["Scharlau_Kirschmer_target"]["entry_one_based"] for row in h11_maxima] == [2, 2, 2]
    assert [row["Scharlau_Kirschmer_target"]["entry_one_based"] for row in h10_maxima] == [8, 9, 9]
    assert [row["section_to_full_reflection_group_index"] for row in h10_maxima] == [None, 2046, 2046]
    assert all(
        row["exact_full_chamber_bridge"]["source_type_id"] == "R13_UE8A1x3_C3_S3"
        and all(row["exact_full_chamber_bridge"]["checks"].values())
        for row in h10_maxima
        if row["even"]
    )
    assert all(
        not case["lineage"]["CoxIter"]["ordinary_exact_matches"]
        and not case["lineage"]["CoxIter"]["doubling_exact_matches"]
        and case["lineage"]["new_arithmetic_commensurability_class"] is False
        for case in CERTIFICATE["positive_types"]
    )


def test_no_survivor_is_promoted_to_a_general_surface_action() -> None:
    assert all(
        not case["surface_geometry"]["selected_nonidentity_subgroup_realized_on_general_surface"]
        and case["surface_geometry"]["surface_automorphism_not_needed_for_polytope_proof"]
        for case in CERTIFICATE["positive_types"]
    )


def test_full_certificate_replays() -> None:
    assert build_certificate() == CERTIFICATE
