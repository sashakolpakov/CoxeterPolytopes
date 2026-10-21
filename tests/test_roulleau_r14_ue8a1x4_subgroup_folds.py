import json
from pathlib import Path

from scripts.certify_roulleau_r14_ue8a1x4_subgroup_folds import (
    _coordinate_model,
    _printed_tikz_gram,
    build_certificate,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/roulleau_r14_ue8a1x4_all_subgroup_folds.json"
CERTIFICATE = json.loads(ARTIFACT.read_text())


def test_coordinate_lane_and_printed_cross_graph_inconsistency_are_frozen() -> None:
    basis_gram, walls, coordinate_gram = _coordinate_model()
    literal = _printed_tikz_gram(thick_precedence=False)
    thick_first = _printed_tikz_gram(thick_precedence=True)
    source = CERTIFICATE["source_reconstruction"]

    assert basis_gram.det() == -16
    assert walls.shape == (14, 27)
    assert coordinate_gram.rank() == 14
    assert exact_inertia(coordinate_gram) == (1, 13, 13)
    assert all(coordinate_gram[index, index] == -2 for index in range(27))

    # Never silently replace the coordinate/equivalence reconstruction with
    # the impossible final TikZ cross-incidence labels in the paper.
    assert literal.rank() == 24
    assert exact_inertia(literal) == (6, 18, 3)
    assert thick_first.rank() == 22
    assert exact_inertia(thick_first) == (5, 17, 5)
    assert source["source_inconsistency"]["coordinate_vs_last_draw_lower_triangle_difference_count"] == 25
    assert "coordinate/equivalence lane is authoritative" in source["source_inconsistency"]["conclusion"]
    assert all(source["checks"].values())


def test_complete_s4_times_c2_subgroup_denominator_and_four_diagrams() -> None:
    assert CERTIFICATE["schema"] == "coxeter.roulleau_r14_ue8a1x4.all_subgroup_folds.v1"
    assert CERTIFICATE["all_checks_pass"] is True
    group = CERTIFICATE["weighted_wall_group"]
    census = CERTIFICATE["subgroup_census"]
    assert group["order"] == 48
    assert group["recognition"] == "S4 x C2"
    assert census["subgroup_count"] == 98
    assert census["conjugacy_class_count"] == 33
    assert census["high_class_to_distinct_diagram_representative"] == {
        "2": 2,
        "3": 3,
        "7": 7,
        "8": 3,
        "11": 11,
        "12": 11,
        "17": 7,
        "20": 11,
        "28": 11,
        "30": 11,
    }
    assert [(case["geometry"]["hyperbolic_dimension"], case["geometry"]["facet_count"]) for case in CERTIFICATE["cases"]] == [
        (12, 23),
        (11, 19),
        (11, 19),
        (10, 15),
    ]


def test_all_four_literal_vinberg_volume_certificates() -> None:
    expected = [
        (2308, 320, 17),
        (795, 114, 13),
        (833, 128, 7),
        (230, 37, 5),
    ]
    observed = []
    for case in CERTIFICATE["cases"]:
        vinberg = case["geometry"]["vinberg"]
        observed.append(
            (
                vinberg["elliptic_edge_count"],
                vinberg["finite_vertex_count"],
                vinberg["ideal_vertex_count"],
            )
        )
        assert vinberg["finite_volume"] is True
        assert vinberg["compact"] is False
        assert vinberg["bad_edge_count"] == 0
        assert vinberg["bad_parabolic_neighborhood_count"] == 0
        assert all(vinberg["checks"].values())
        assert case["folded_walls"]["strict_facet_witnesses"]["all_orbit_inequalities_are_genuine_facets"]
    assert observed == expected


def test_stable_overclosures_and_bounded_lineage() -> None:
    assert [case["fixed_lattice"]["determinant"] for case in CERTIFICATE["cases"]] == [16, -16, -12, 8]
    assert [case["stable_overlattices"]["unique_maximal_reflection_stable_overclosure"]["index"] for case in CERTIFICATE["cases"]] == [4, 4, 2, 2]
    assert [case["lineage"]["section_to_full_reflection_group_index"] for case in CERTIFICATE["cases"]] == [8320, None, None, 1056]
    assert all(
        not case["lineage"]["CoxIter"]["ordinary_exact_matches"]
        and not case["lineage"]["CoxIter"]["doubling_exact_matches"]
        and case["lineage"]["new_arithmetic_commensurability_class"] is False
        for case in CERTIFICATE["cases"]
    )
    collision = CERTIFICATE["cases"][2]["lineage"]["original_fixed_lattice_SK_collision"]
    assert collision["entry_one_based"] == 12
    assert collision["finite_special_counts_and_cusps_match_exactly"] is True


def test_algebraic_carrier_does_not_assert_a_general_surface_action() -> None:
    assert all(
        case["surface_geometry"]["ambient_object"].startswith("complete nef cone of a K3 surface")
        and case["surface_geometry"]["surface_automorphism_not_needed_for_polytope_proof"]
        and not case["surface_geometry"]["selected_nonidentity_graph_subgroup_realized_on_general_surface"]
        for case in CERTIFICATE["cases"]
    )


def test_full_certificate_replays() -> None:
    assert build_certificate() == CERTIFICATE
