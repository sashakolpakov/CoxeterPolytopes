from __future__ import annotations

from collections import Counter
from functools import lru_cache
import json

from scripts.certify_h9_14_t2_surface_unfolding_pilot import (
    OUTPUT,
    build_certificate,
)


@lru_cache(maxsize=1)
def _payload() -> dict:
    return build_certificate()


def test_archive_is_an_exact_replay() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    source = payload["surface_free_logical_input"]
    assert source["construction_id"] == "A18_4_0_C01"
    assert source["facet_gram_string_canonical_json_sha256"] == (
        "6979b16de7583c20277077a142df863e9f9a7209c443e9252505234b15d66c56"
    )
    assert source["rank_and_inertia"] == [10, 1, 9, 4]
    assert source["normalized_coxeter_relation_counts"] == {
        "m2": 75,
        "m3": 8,
        "m4": 8,
    }


def test_complete_scaling_and_fixed_overlattice_census() -> None:
    compiler = _payload()["finite_surface_free_compiler"]
    assert compiler["scaling_ray_count"] == 18
    assert compiler["diagram_automorphism_group_order"] == 48
    assert compiler["diagram_automorphism_orbit_sizes"] == [1, 1, 2, 4, 4, 6]
    assert compiler["two_sheet_compatible_scaling_count"] == 2
    assert compiler["norm_eight_obstructed_scaling_count"] == 16
    lanes = {row["label"]: row for row in compiler["lanes"]}
    assert lanes["archived_orientation"]["scaling_number"] == 1
    assert lanes["dual_orientation"]["scaling_number"] == 9
    assert len(lanes["archived_orientation"]["candidates"]) == 16
    assert len(lanes["dual_orientation"]["candidates"]) == 5
    assert Counter(
        row["glue_dimension"]
        for row in lanes["archived_orientation"]["candidates"]
    ) == {0: 1, 1: 7, 2: 7, 3: 1}
    assert Counter(
        row["glue_dimension"]
        for row in lanes["dual_orientation"]["candidates"]
    ) == {0: 1, 1: 3, 2: 1}


def test_obstructions_are_separate_and_leave_one_candidate() -> None:
    compiler = _payload()["finite_surface_free_compiler"]
    lanes = {row["label"]: row for row in compiler["lanes"]}
    archived_root_span = lanes["archived_orientation"]["candidates"][0]
    assert archived_root_span["long_class_span_dimension"] == 7
    assert archived_root_span["long_class_radical_dimension"] == 3
    assert archived_root_span["E8_class_obstruction"] == (
        "radical_dimension_exceeds_8_minus_dimension"
    )
    dual_root_span = lanes["dual_orientation"]["candidates"][0]
    assert dual_root_span["gates"]["embeddable_as_nonsingular_E8_classes"]
    assert dual_root_span["minimal_rank18_discriminant_length"] == 8
    assert dual_root_span["gates"]["rank4_K3_complement_length"] is False
    survivors = [
        row
        for lane in lanes.values()
        for row in lane["candidates"]
        if row["survives_minimal_two_sheet_weak_K3_gate"]
    ]
    assert len(survivors) == 1
    assert survivors[0]["candidate_number"] == 8
    assert survivors[0]["glue_basis_little_endian_hex"] == ["0x28a"]


def test_unique_u_e8_minus2_unfolding_and_weak_k3_embedding() -> None:
    survivor = _payload()["unique_survivor"]
    fixed = survivor["fixed_lattice"]
    assert fixed["isometry_class"] == "U+E8(-2)"
    assert fixed["two_elementary_invariants"] == [10, 8, 0]
    assert fixed["determinant"] == -256
    assert fixed["smith_diagonal_absolute"] == [1, 1] + [2] * 8
    e8 = survivor["E8_two_sheet_data"]
    assert e8["target_basis_classes_little_endian_hex"] == [
        "01", "02", "04", "10", "20", "40"
    ]
    assert e8["all_target_classes_little_endian_hex"] == [
        "01", "02", "04", "10", "20", "40", "45", "54", "15", "51"
    ]
    assert e8["absolute_anti_pairings_equal_folded_long_pairings"] is True
    rank18 = survivor["rank18_unfolding"]
    assert rank18["rank_signature_determinant"] == [18, 1, 17, -16]
    assert rank18["smith_diagonal_absolute"] == [1] * 14 + [2] * 4
    assert rank18["two_elementary_invariants"] == [18, 4, 0]
    assert rank18["fixed_and_coinvariant_embeddings_primitive"] is True
    roots = survivor["source_roots_and_fold"]
    assert roots["source_root_count"] == 24
    assert roots["one_sheet_orbit_count"] == 4
    assert roots["two_sheet_orbit_count"] == 10
    assert roots["off_diagonal_pairing_histogram"] == {"0": 244, "1": 32}
    assert roots["source_roots_span_N_with_index_one"] is True
    k3 = survivor["weak_K3_embedding_and_involution_extension"]
    assert k3["k3_lattice_rank_signature_determinant"] == [22, 3, 19, -1]
    assert k3["N_and_T_embeddings_primitive"] is True
    assert k3["extended_involution_characteristic_polynomial"] == (
        "(lambda-1)^14*(lambda+1)^8"
    )


def test_exact_marked_bridge_and_claim_boundary() -> None:
    payload = _payload()
    bridge = payload["comparison_only_marked_bridge"]
    assert bridge["active_source_wall_orbits_one_based"] == [
        [1], [3, 9], [4, 8], [5, 7], [6], [11, 17], [12, 16],
        [13, 15], [14], [18, 24], [19, 25], [22, 32], [23, 31], [30],
    ]
    assert bridge["sheet_assignment_count_checked"] == 1024
    assert bridge["integral_marked_isometry_count"] == 4
    witness = bridge["canonical_witness"]
    assert witness["sheet_swap_mask_little_endian_hex"] == "0x0"
    assert witness["determinant"] == -1
    assert witness["matrix_canonical_json_sha256"] == (
        "15f683beded0f7366fe8be371298e5c17bf42238e36c2bec7833398acf6d6e85"
    )
    assert all(witness[key] is True for key in (
        "maps_all_24_source_roots_exactly",
        "intertwines_involutions_exactly",
        "is_integral_unimodular_isometry",
    ))
    geometry = payload["known_geometry_from_archives"]
    assert geometry["classification"] == "finite_volume_noncompact"
    assert geometry["literal_vinberg"]["finite_endpoint_count"] == 90
    assert geometry["literal_vinberg"]["ideal_endpoint_count"] == 20
    assert geometry["literal_vinberg"]["bad_edge_count"] == 0
    assert geometry["literal_vinberg"]["bad_parabolic_neighborhood_count"] == 0
    boundary = payload["claim_boundary"]
    assert boundary["new_hyperbolic_polytope_claim"] is False
    assert boundary["class_wide_surface_correspondence_claim"] is False
    assert boundary["class_wide_IC2_claim"] is False
    assert boundary["class_wide_IC4_claim"] is False
