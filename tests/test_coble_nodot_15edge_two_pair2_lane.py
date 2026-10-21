from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json

from scripts.certify_coble_nodot_15edge_two_pair2_lane import OUTPUT, build


ARCHIVE_SHA256 = "85a2a347b84954b20e3f796958a244755fd994a76508fbfa3bad84f077c15167"


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_complete_exact_unordered_edge_pair_denominator() -> None:
    assert sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert payload["scope"] == {
        "all_Lorentzian_corank_two_determinants_retained": True,
        "all_even_unimodular_E10_extensions_retained": True,
        "all_other_nonzero_pairings_equal_one": True,
        "all_unordered_marked_edge_pairs_of_every_unlabelled_support_tested": True,
        "connected_ordinary_support": True,
        "exact_rank_on_every_marked_trial_without_float_prefilter": True,
        "not_the_full_D_1_split_nodot_denominator": True,
        "ordinary_support_edge_count": 15,
        "ordinary_support_maximum_degree": 3,
        "pairing_alphabet": [0, 1, 2],
        "pairing_two_edge_count": 2,
        "vertices": 12,
        "weighted_trials_are_edge_pair_markings_not_vertex_labellings": True,
    }
    denominator = payload["denominator"]
    assert denominator["underlying_graph_count"] == 5_741
    assert denominator["weighted_trial_count"] == 602_805
    assert denominator["exact_rank_histogram"] == {
        "10": 83,
        "11": 7_362,
        "12": 595_360,
    }
    assert denominator["rank_ten_exact_inertia_histogram"] == {"9,1,2": 83}
    assert denominator["exact_Lorentzian_corank_two_marked_trial_count"] == 83
    assert denominator["weighted_diagram_isometry_class_count"] == 70
    assert Counter(
        record["marked_trial_member_count"] for record in payload["class_records"]
    ) == {1: 60, 2: 8, 3: 1, 4: 1}


def test_determinant_and_all_e10_saturations() -> None:
    payload = _archive()
    denominator = payload["denominator"]
    assert denominator["radical_quotient_determinant_histogram"] == {
        "-1": 51,
        "-12": 1,
        "-13": 2,
        "-16": 1,
        "-25": 2,
        "-36": 3,
        "-37": 2,
        "-4": 7,
        "-40": 1,
        "-44": 2,
        "-45": 2,
        "-48": 2,
        "-5": 3,
        "-52": 1,
        "-56": 1,
        "-64": 2,
    }
    assert Counter(
        record["radical_quotient"]["even_unimodular_E10_extension_count"]
        for record in payload["class_records"]
    ) == {0: 16, 1: 51, 2: 3}
    assert denominator["classes_without_even_unimodular_E10_extension"] == 16
    assert denominator["classes_rejected_by_nonsquare_quotient_determinant"] == 16
    assert denominator["square_determinant_classes_without_E10_extension"] == 0
    assert denominator["integral_E10_extension_branch_count"] == 57


def test_all_seventy_root_cones_are_hyperideal() -> None:
    payload = _archive()
    records = payload["class_records"]
    assert payload["root_cone_volume_ledger"] == {
        "infinite_volume_hyperideal": 70
    }
    assert Counter(
        record["root_cone"]["hyperideal_vertex_count"] for record in records
    ) == {
        10: 2,
        12: 3,
        13: 1,
        14: 3,
        15: 2,
        16: 1,
        17: 2,
        18: 1,
        19: 1,
        20: 5,
        21: 1,
        22: 2,
        24: 2,
        25: 3,
        26: 3,
        27: 1,
        29: 3,
        30: 1,
        31: 2,
        32: 1,
        33: 1,
        34: 3,
        35: 3,
        37: 2,
        38: 3,
        39: 2,
        40: 3,
        41: 3,
        42: 2,
        44: 6,
        45: 1,
        49: 1,
    }
    assert Counter(
        record["root_positive_strict_interior"]["finite_ray_multiplier"]
        for record in records
    ) == {0: 68, 1: 1, 2: 1}
    assert all(
        record["root_cone_classification"] == "infinite_volume_hyperideal"
        and record["root_cone"]["all_input_walls_irredundant"]
        and record["root_cone"]["strict_interior_from_ray_sum"]
        and record["root_positive_strict_interior"]["square"] > 0
        and record["root_positive_strict_interior"][
            "all_wall_slacks_strictly_positive"
        ]
        and record["root_cone"]["all_causal_rays_same_time_cone"]
        and record["root_literal_vinberg"]["stopped_at_hyperideal_endpoint"]
        for record in records
    )
    assert len(
        {
            record["root_cone_ray_hashes"]["hyperideal_ray_records_sha256"]
            for record in records
        }
    ) == 70


def test_all_fifty_seven_e10_branches_are_hyperideal() -> None:
    payload = _archive()
    branches = [
        branch
        for record in payload["class_records"]
        for branch in record["extension_branches"]
    ]
    assert len(branches) == 57
    assert payload["E10_extension_volume_ledger"] == {
        "infinite_volume_hyperideal": 57
    }
    assert Counter(branch["classification"] for branch in branches) == {
        "infinite_volume_hyperideal": 57
    }
    assert all(
        branch["cone"]["hyperideal_vertex_count"] > 0
        and branch["positive_strict_interior"]["square"] > 0
        and branch["positive_strict_interior"][
            "all_wall_slacks_strictly_positive"
        ]
        and branch["literal_vinberg"]["stopped_at_hyperideal_endpoint"]
        and branch["full_cusp_lift"] is None
        for branch in branches
    )
    assert all(
        (
            branch["cone"]["finite_vertex_count"],
            branch["cone"]["ideal_vertex_count"],
            branch["cone"]["hyperideal_vertex_count"],
        )
        == (
            record["root_cone"]["finite_vertex_count"],
            record["root_cone"]["ideal_vertex_count"],
            record["root_cone"]["hyperideal_vertex_count"],
        )
        for record in payload["class_records"]
        for branch in record["extension_branches"]
    )
    assert payload["full_D_1_split_branch_count"] == 0
    assert payload["full_D_1_split_isomorphism_class_count"] == 0
    assert payload["full_survivor_ledger"] == []


def test_archive_replays_byte_for_byte() -> None:
    assert build() == _archive()


def test_scope_and_lineage_are_not_overstated() -> None:
    payload = _archive()
    claims = payload["claims"]
    assert claims["strict_L4_for_full_D_1_split_nodot"] is False
    assert claims["new_abstract_hyperbolic_polytope_claimed"] is False
    assert claims["every_Lorentzian_root_cone_has_a_hyperideal_ray"] is True
    assert claims["no_finite_volume_survivor_in_the_scoped_lane"] is True
    assert "exactly-two-pairing-two" in claims["reason"]
    assert all(
        record["known_root_diagram_matches"] == []
        for record in payload["class_records"]
    )
