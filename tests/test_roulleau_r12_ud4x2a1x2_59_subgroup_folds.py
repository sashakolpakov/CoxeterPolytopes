import json

from scripts.certify_roulleau_r12_ud4x2a1x2_59_subgroup_folds import (
    OUTPUT,
    build_certificate,
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text())


def test_complete_subgroup_denominator_and_two_high_classes() -> None:
    payload = _archive()
    census = payload["subgroup_census"]
    assert payload["all_checks_pass"] is True
    assert payload["weighted_wall_group"]["order"] == 2880
    assert payload["weighted_wall_group"]["structure"] == "C2 x C2 x S6"
    assert census["class_count"] == 914
    assert census["proper_high_class_indices"] == [4, 12]
    assert census["exploratory_all_element_GAP_aliases"] == {"4": 3, "12": 9}
    high = [
        (row["class_index"], row["order"], row["conjugacy_class_size"], row["fixed_rank"])
        for row in census["classes"]
        if row["proper_fixed_section_above_H9"]
    ]
    assert high == [(4, 2, 1, 11), (12, 2, 15, 11)]


def test_two_literal_vinberg_positive_examples() -> None:
    expected = [
        ("R12-D4x2A1x2-p67-H10-38", 38, 4518, 526, 141, 307),
        ("R12-D4x2A1x2-S6transposition-H10-45", 45, 7332, 1048, 93, 302),
    ]
    observed = []
    for case in _archive()["cases"]:
        vinberg = case["vinberg_finite_volume_audit"]
        observed.append(
            (
                case["id"],
                len(case["normalized_relation_matrix"]),
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
        assert case["strict_facet_witnesses"][
            "all_orbit_inequalities_are_genuine_facets"
        ]
        assert all(case["exact_coxeter_gate"].values())
    assert observed == expected


def test_stable_p10_lineage_and_exact_indices() -> None:
    expected = [(15, 2276640, "-527/2580480"), (3, 2023680, "-527/2903040")]
    observed = []
    for case in _archive()["cases"]:
        fixed = case["fixed_lattice"]
        lineage = case["maximal_reflection_lineage"]
        maxima = fixed["maximal_stable_overlattices"]
        assert fixed["determinant"] == 64
        assert fixed["smith_invariants"] == [1] * 6 + [2] * 4 + [4]
        assert all(
            row["index_over_fixed_lattice"] == 8
            and row["determinant"] == 1
            and row["parity"] == "odd"
            for row in maxima
        )
        observed.append(
            (
                len(maxima),
                lineage["index_in_each_P10_reflection_group"],
                lineage["section_euler_characteristic"],
            )
        )
        assert not case["bounded_lineage"]["coxiter_scan"]["ordinary_exact_matches"]
        assert not case["bounded_lineage"]["coxiter_scan"]["doubling_exact_matches"]
    assert observed == expected
    p10 = _archive()["P10_denominator"]
    assert p10["solomon"]["rational_euler_characteristic"] == "-1/11147673600"
    assert p10["checks"]["componentwise_helper_equals_direct"] is True


def test_neither_action_is_promoted_to_the_general_surface() -> None:
    payload = _archive()
    assert payload["independent_surface_action_certificate"]["all_checks_pass"]
    assert payload["weighted_wall_group"]["first_D2_deck_fixed_rank"] == 4
    assert payload["weighted_wall_group"]["second_D8_deck_NS_action"] == "identity"
    p67, transposition = payload["cases"]
    assert p67["geometric_realization"]["p67_exact_permutation_match"] is True
    assert transposition["geometric_realization"]["p67_exact_permutation_match"] is False
    assert all(
        not case["geometric_realization"]["is_general_surface_automorphism_action"]
        and "special_period_open" in case["geometric_realization"]["status"]
        and case["geometric_realization"][
            "selected_action_excluded_from_general_surface_image"
        ]
        for case in payload["cases"]
    )


def test_full_certificate_replays() -> None:
    assert build_certificate() == _archive()
