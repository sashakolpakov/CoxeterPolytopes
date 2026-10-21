from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json

from scripts.certify_coble_nodot_15edge_three_pair2_lane import OUTPUT, build


ARCHIVE_SHA256 = "581480510c340cb77ad21f017aa48ace1e072d0f07ba35c6f54eb895a5b6be2c"


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_complete_exact_unordered_edge_triple_denominator() -> None:
    assert sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert payload["scope"] == {
        "all_Lorentzian_corank_two_determinants_retained": True,
        "all_even_unimodular_E10_extensions_retained": True,
        "all_other_nonzero_pairings_equal_one": True,
        "all_unordered_marked_edge_triples_of_every_unlabelled_support_tested": True,
        "connected_ordinary_support": True,
        "exact_rank_on_every_marked_trial_without_float_prefilter": True,
        "not_the_full_D_1_split_nodot_denominator": True,
        "ordinary_support_edge_count": 15,
        "ordinary_support_maximum_degree": 3,
        "pairing_alphabet": [0, 1, 2],
        "pairing_two_edge_count": 3,
        "vertices": 12,
        "weighted_trials_are_edge_triple_markings_not_vertex_labellings": True,
    }
    denominator = payload["denominator"]
    assert denominator["underlying_graph_count"] == 5_741
    assert denominator["weighted_trial_count"] == 2_612_155
    assert denominator["exact_rank_histogram"] == {
        "9": 3,
        "10": 208,
        "11": 27_509,
        "12": 2_584_435,
    }
    assert denominator["rank_ten_exact_inertia_histogram"] == {
        "8,2,2": 84,
        "9,1,2": 124,
    }
    assert denominator["exact_Lorentzian_corank_two_marked_trial_count"] == 124
    assert denominator["weighted_diagram_isometry_class_count"] == 111
    assert Counter(
        record["marked_trial_member_count"] for record in payload["class_records"]
    ) == {1: 100, 2: 10, 4: 1}


def test_determinant_and_all_e10_saturations() -> None:
    payload = _archive()
    denominator = payload["denominator"]
    assert denominator["radical_quotient_determinant_histogram"] == {
        "-1": 59,
        "-13": 8,
        "-16": 9,
        "-28": 1,
        "-32": 1,
        "-4": 21,
        "-40": 1,
        "-52": 4,
        "-56": 1,
        "-64": 7,
        "-72": 3,
        "-76": 2,
        "-8": 3,
        "-80": 1,
        "-9": 3,
    }
    assert Counter(
        record["radical_quotient"]["even_unimodular_E10_extension_count"]
        for record in payload["class_records"]
    ) == {0: 19, 1: 85, 2: 7}
    assert denominator["classes_without_even_unimodular_E10_extension"] == 19
    assert denominator["classes_rejected_by_nonsquare_quotient_determinant"] == 19
    assert denominator["square_determinant_classes_without_E10_extension"] == 0
    assert denominator["integral_E10_extension_branch_count"] == 99


def test_all_111_root_cones_are_hyperideal_with_timelike_witnesses() -> None:
    payload = _archive()
    records = payload["class_records"]
    assert payload["root_cone_volume_ledger"] == {
        "infinite_volume_hyperideal": 111
    }
    assert Counter(
        record["root_cone"]["hyperideal_vertex_count"] for record in records
    ) == {
        11: 4,
        12: 3,
        13: 2,
        14: 1,
        16: 1,
        17: 2,
        18: 5,
        19: 4,
        20: 3,
        21: 3,
        22: 5,
        23: 2,
        24: 7,
        25: 1,
        26: 2,
        28: 3,
        30: 3,
        31: 3,
        32: 11,
        33: 2,
        34: 3,
        35: 3,
        36: 5,
        38: 4,
        39: 6,
        40: 1,
        41: 1,
        42: 5,
        43: 1,
        44: 5,
        45: 3,
        46: 1,
        47: 1,
        48: 2,
        49: 1,
        58: 1,
        61: 1,
    }
    assert Counter(
        (
            len(
                record["root_positive_strict_interior"][
                    "timelike_seed_ray_indices_zero_based"
                ]
            ),
            record["root_positive_strict_interior"]["timelike_seed_multiplier"],
        )
        for record in records
    ) == {
        (1, 0): 76,
        (1, 1): 3,
        (1, 4): 2,
        (1, 8): 1,
        (2, 0): 23,
        (2, 1): 1,
        (2, 2): 2,
        (2, 4): 1,
        (2, 16): 1,
        (2, 64): 1,
    }
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
    ) == 111


def test_all_99_e10_branches_are_hyperideal() -> None:
    payload = _archive()
    branches = [
        branch
        for record in payload["class_records"]
        for branch in record["extension_branches"]
    ]
    assert len(branches) == 99
    assert payload["E10_extension_volume_ledger"] == {
        "infinite_volume_hyperideal": 99
    }
    assert Counter(branch["classification"] for branch in branches) == {
        "infinite_volume_hyperideal": 99
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
    assert "exactly-three-pairing-two" in claims["reason"]
    assert all(
        record["known_root_diagram_matches"] == []
        for record in payload["class_records"]
    )
