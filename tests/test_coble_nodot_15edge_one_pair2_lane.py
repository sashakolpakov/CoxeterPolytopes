from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json

from scripts.certify_coble_nodot_15edge_one_pair2_lane import OUTPUT, build


ARCHIVE_SHA256 = "1e3e20026c418c6857ceb7a834693fd9217479e05506ddcab93d16304c38b084"


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_complete_exact_marked_edge_denominator() -> None:
    assert sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert payload["scope"] == {
        "all_Lorentzian_corank_two_determinants_retained": True,
        "all_even_unimodular_E10_extensions_retained": True,
        "all_marked_edges_of_every_unlabelled_support_tested": True,
        "all_other_nonzero_pairings_equal_one": True,
        "connected_ordinary_support": True,
        "exact_rank_on_every_marked_trial_without_float_prefilter": True,
        "not_the_full_D_1_split_nodot_denominator": True,
        "ordinary_support_edge_count": 15,
        "ordinary_support_maximum_degree": 3,
        "pairing_alphabet": [0, 1, 2],
        "pairing_two_edge_count": 1,
        "vertices": 12,
        "weighted_trials_are_edge_markings_not_vertex_labellings": True,
    }
    denominator = payload["denominator"]
    assert denominator["underlying_graph_count"] == 5_741
    assert denominator["weighted_trial_count"] == 86_115
    assert denominator["exact_rank_histogram"] == {
        "10": 29,
        "11": 1_462,
        "12": 84_624,
    }
    assert denominator["rank_ten_exact_inertia_histogram"] == {"9,1,2": 29}
    assert denominator["exact_Lorentzian_corank_two_marked_trial_count"] == 29
    assert denominator["weighted_diagram_isometry_class_count"] == 28
    assert sum(
        record["marked_trial_member_count"] for record in payload["class_records"]
    ) == 29


def test_every_even_unimodular_saturation_is_retained() -> None:
    payload = _archive()
    denominator = payload["denominator"]
    assert denominator["radical_quotient_determinant_histogram"] == {
        "-1": 9,
        "-12": 3,
        "-16": 2,
        "-36": 1,
        "-4": 1,
        "-40": 1,
        "-48": 2,
        "-49": 1,
        "-56": 1,
        "-61": 1,
        "-64": 3,
        "-72": 1,
        "-76": 1,
        "-80": 1,
        "-9": 1,
    }
    extension_counts = Counter(
        record["radical_quotient"]["even_unimodular_E10_extension_count"]
        for record in payload["class_records"]
    )
    assert extension_counts == {0: 10, 1: 15, 2: 3}
    assert denominator["classes_without_even_unimodular_E10_extension"] == 10
    assert denominator["classes_rejected_by_nonsquare_quotient_determinant"] == 10
    assert denominator["square_determinant_classes_without_E10_extension"] == 0
    assert denominator["integral_E10_extension_branch_count"] == 21


def test_all_twenty_one_branches_are_exactly_infinite_volume() -> None:
    payload = _archive()
    root_records = payload["class_records"]
    assert payload["root_cone_volume_ledger"] == {
        "infinite_volume_hyperideal": 28
    }
    assert Counter(
        record["root_cone"]["hyperideal_vertex_count"] for record in root_records
    ) == {
        7: 1,
        11: 1,
        13: 2,
        15: 1,
        16: 2,
        17: 1,
        18: 1,
        19: 1,
        20: 2,
        22: 4,
        24: 1,
        25: 1,
        26: 1,
        27: 3,
        28: 1,
        30: 2,
        31: 1,
        34: 2,
    }
    assert all(
        record["root_cone_classification"] == "infinite_volume_hyperideal"
        and record["root_cone"]["all_input_walls_irredundant"]
        and record["root_cone"]["strict_interior_from_ray_sum"]
        and record["root_cone"]["strict_interior_square"] > 0
        and record["root_cone"]["all_causal_rays_same_time_cone"]
        and record["root_cone"]["hyperideal_vertex_count"] > 0
        and record["root_literal_vinberg"]["stopped_at_hyperideal_endpoint"]
        and len(record["root_cone_ray_hashes"]["extreme_ray_records_sha256"])
        == 64
        and len(record["root_cone_ray_hashes"]["hyperideal_ray_records_sha256"])
        == 64
        for record in root_records
    )
    branches = [
        branch
        for record in payload["class_records"]
        for branch in record["extension_branches"]
    ]
    assert len(branches) == 21
    assert Counter(branch["classification"] for branch in branches) == {
        "infinite_volume_hyperideal": 21
    }
    assert Counter(branch["cone"]["hyperideal_vertex_count"] for branch in branches) == {
        7: 2,
        13: 1,
        15: 2,
        16: 3,
        17: 1,
        18: 1,
        22: 1,
        25: 1,
        26: 1,
        27: 2,
        28: 1,
        30: 2,
        31: 1,
        34: 2,
    }
    assert all(
        branch["cone"]["all_input_walls_irredundant"]
        and branch["cone"]["strict_interior_from_ray_sum"]
        and branch["cone"]["strict_interior_square"] > 0
        and branch["cone"]["all_causal_rays_same_time_cone"]
        and branch["cone"]["all_ideal_supports_parabolic_rank_n_minus_one"]
        and not branch["cone"][
            "complete_cdd_V_representation_has_no_hyperideal_ray"
        ]
        and branch["literal_vinberg"]["stopped_at_hyperideal_endpoint"]
        and branch["full_cusp_lift"] is None
        for branch in branches
    )
    assert payload["E10_extension_volume_ledger"] == {
        "infinite_volume_hyperideal": 21
    }
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
        for record in root_records
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
    assert claims["hyperideal_branches_are_infinite_volume"] is True
    assert claims["every_Lorentzian_root_cone_has_a_hyperideal_ray"] is True
    assert claims["no_finite_volume_survivor_in_the_scoped_lane"] is True
    assert "exactly-one-pairing-two" in claims["reason"]
    assert all(
        record["known_root_diagram_matches"] == []
        for record in payload["class_records"]
    )
