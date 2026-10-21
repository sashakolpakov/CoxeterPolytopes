from __future__ import annotations

from fractions import Fraction
import json

from scripts.certify_alexeev_remaining_fixed_sections import (
    EXPECTED_HIGH_CLASSES,
    EXPECTED_SECTIONS,
    OUTPUT as CONSTRUCTION,
    _payload_sha256,
)
from scripts.certify_alexeev_remaining_fixed_sections_vinberg import (
    EXPECTED as EXPECTED_VOLUME,
    OUTPUT as VINBERG,
)
from scripts.certify_alexeev_remaining_fixed_sections_lineage import (
    EXPECTED_DUPLICATE_PERMUTATION,
    EXPECTED_EULER,
    OUTPUT as LINEAGE,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _euler_from_connected_components,
    _relation_matrix,
)
from scripts.certify_shimada_disc16_d4_duplicate_survivors import _match
from scripts.certify_shimada_disc16_d5_h10_lineage import (
    _connected_elliptic_masks,
)


D4_C007 = (
    CONSTRUCTION.parents[0]
    / "shimada_disc16_atlas_archived_high_sections"
    / "D0004_C007.json"
)


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_all_three_weighted_denominators_and_high_classes_are_complete() -> None:
    payload = _load(CONSTRUCTION)
    assert payload["all_checks_pass"] is True
    assert payload["source"]["sha256"] == (
        "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
    )
    assert payload["inventory"]["weighted_group_orders"] == {
        "18,4,0": 1152,
        "15,11,1": 72,
        "15,7,1": 120,
    }
    assert payload["inventory"]["subgroup_class_counts"] == {
        "18,4,0": 221,
        "15,11,1": 26,
        "15,7,1": 19,
    }
    by_key = {
        tuple(example["key"]): example
        for example in payload["inventory"]["examples"]
    }
    assert set(by_key) == set(EXPECTED_HIGH_CLASSES)
    for key, expected in EXPECTED_HIGH_CLASSES.items():
        high = [
            (
                row["class_index"],
                row["order"],
                row["conjugacy_class_size"],
                row["fixed_rank"],
                row["hyperbolic_dimension"],
                row["cyclic"],
            )
            for row in by_key[key]["subgroup_classes"]
            if row["proper_fixed_section_above_h9"]
        ]
        assert high == expected
        group = by_key[key]["weighted_diagram_group"]
        assert len(group["elements"]) == group["order"]
        assert group["all_actions_integral_unimodular_Lorentz"] is True


def test_all_six_sections_are_strict_integral_coxeter_chambers() -> None:
    payload = _load(CONSTRUCTION)
    by_id = {case["id"]: case for case in payload["cases"]}
    assert set(by_id) == set(EXPECTED_SECTIONS)
    for identifier, expected in EXPECTED_SECTIONS.items():
        case = by_id[identifier]
        lattice = case["saturated_fixed_lattice"]
        reduction = case["exact_reduction"]
        assert case["case"]["hyperbolic_dimension"] == expected["dimension"]
        assert case["case"]["wall_orbit_count"] == expected["orbits"]
        assert reduction["active_hyperbolic_facet_count"] == expected["facets"]
        assert reduction["conic_redundant_count"] == expected["redundant"]
        assert reduction["unresolved_orbit_indices"] == []
        assert len(reduction["strict_facet_witnesses"]) == expected["facets"]
        assert lattice["determinant"] == expected["determinant"]
        assert lattice["smith_diagonal_absolute"] == expected["smith"]
        assert lattice["inertia"] == [1, expected["dimension"], 0]
        assert case["coxeter_gate"]["coxeter_edge_counts"] == expected[
            "pair_counts"
        ]
        assert case["coxeter_gate"]["all_pairs_coxeter"] is True
        assert case["restricted_mirrors"][
            "all_active_reflections_integral_on_saturated_fixed_lattice"
        ] is True
        assert _payload_sha256(
            case["restricted_mirrors"]["positive_wall_gram"]
        ) == expected["facet_gram_sha256"]
    assert by_id["A15_7_1_C01"]["strict_fixed_interior"] == {
        **by_id["A15_7_1_C01"]["strict_fixed_interior"],
        "source_method": "caller-supplied exact strict timelike interior",
        "source_square": "7360",
    }


def test_literal_vinberg_gate_proves_six_finite_volume_noncompact_cases() -> None:
    construction = {case["id"]: case for case in _load(CONSTRUCTION)["cases"]}
    payload = _load(VINBERG)
    by_id = {case["id"]: case for case in payload["cases"]}
    assert payload["summary"] == {
        "bad_edge_count": 0,
        "bad_parabolic_neighborhood_count": 0,
        "case_count": 6,
        "finite_volume_noncompact_count": 6,
    }
    assert set(by_id) == set(EXPECTED_VOLUME)
    for identifier, expected in EXPECTED_VOLUME.items():
        row = by_id[identifier]
        source = construction[identifier]
        assert row["dimension"] == expected["dimension"]
        assert row["facet_count"] == expected["facets"]
        assert row["finite_special_subset_counts_by_rank"] == expected[
            "finite_special_counts"
        ]
        assert row["elliptic_edge_count"] == expected["finite_special_counts"][-2]
        assert row["finite_vertex_count"] == expected["finite_special_counts"][-1]
        assert row["connected_affine_component_count"] == expected[
            "affine_components"
        ]
        assert row["ideal_vertex_count"] == expected["cusps"]
        assert row["edge_endpoint_type_histogram"] == expected[
            "endpoint_histogram"
        ]
        assert row["bad_edge_count"] == 0
        assert row["bad_parabolic_neighborhood_count"] == 0
        assert row["finite_volume"] is True
        assert row["compact"] is False
        assert source["coxeter_gate"]["all_pairs_coxeter"] is True


def test_h10_solomon_characteristics_and_the_one_exact_old_duplicate() -> None:
    cases = {case["id"]: case for case in _load(CONSTRUCTION)["cases"]}
    expected_chi = {
        "A18_4_0_C04": Fraction(-527, 1_451_520),
        "A18_4_0_C18": Fraction(-527, 23_224_320),
        "A15_11_1_C02": Fraction(-527, 4_976_640),
        "A15_7_1_C01": Fraction(-527, 1_658_880),
    }
    for identifier, expected in expected_chi.items():
        gram = cases[identifier]["restricted_mirrors"]["positive_wall_gram"]
        relations = _relation_matrix(gram)
        connected, _ = _connected_elliptic_masks(gram)
        euler = _euler_from_connected_components(relations, connected)
        assert Fraction(euler["rational_euler_characteristic"]) == expected

    old = _load(D4_C007)["coxeter_gate"]["positive_wall_gram"]
    new = cases["A15_7_1_C01"]["restricted_mirrors"]["positive_wall_gram"]
    comparison = _match(new, old)
    assert comparison["same_hyperbolic_chamber_up_to_isometry"] is True
    assert comparison["source_to_target_facet_permutation_one_based"] == [
        3, 1, 9, 2, 4, 7, 12, 17, 5, 8, 13, 18, 11, 16, 10, 21, 14, 6,
        15, 23, 19, 20, 22, 24, 25, 26, 27,
    ]


def test_bounded_lineage_scope_and_exact_comparisons_are_frozen() -> None:
    payload = _load(LINEAGE)
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert {
        row["id"]: row["solomon"]["rational_euler_characteristic"]
        for row in payload["cases"]
    } == EXPECTED_EULER
    assert payload["prior_ledger_exact_duplicate"][
        "source_to_target_facet_permutation_one_based"
    ] == EXPECTED_DUPLICATE_PERMUTATION
    assert payload["within_batch_same_header_comparisons"] == [
        {
            "same_signed_normalized_gram": False,
            "source": "A18_4_0_C18",
            "source_to_target_facet_permutation_zero_based": None,
            "target": "A15_11_1_C02",
        }
    ]
    coxiter = payload["coxiter_bounded_comparison"]
    assert coxiter["targets"]["H10_facets_28"][
        "doubling_header_match_count"
    ] == 2
    assert coxiter["no_exact_relation_graph_match"] is True
    stable = payload["reflection_stable_overlattice_status"]
    assert stable["status"] == "closed for all six displayed reflection groups"
    assert stable["summary"]["unique_maximal_stable_overlattice_count"] == 6
    assert stable["summary"]["new_arithmetic_commensurability_class_count"] == 0
    assert stable["summary"]["exact_numerical_indices"] == {
        "A15_11_1_C02": 6720,
        "A15_7_1_C01": 2016,
        "A18_4_0_C04": 130560,
        "A18_4_0_C18": 144,
    }
