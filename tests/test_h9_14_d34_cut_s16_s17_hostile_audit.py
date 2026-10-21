from functools import lru_cache
import json
from pathlib import Path

from scripts.audit_h9_14_d34_cut_s16_s17_hostile import (
    EXPECTED_SURVIVOR_GRAM_HASHES,
    EXPECTED_T1_PACKAGE_HASHES,
    EXPECTED_T2_PACKAGE_HASHES,
    OUTPUT,
    build_certificate,
)


@lru_cache(maxsize=1)
def payload() -> dict:
    return build_certificate()


def supports() -> dict[str, dict]:
    return {row["support_id"]: row for row in payload()["support_audits"]}


def test_archive_and_clean_room_boundary() -> None:
    result = payload()
    assert result == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert result["all_checks_pass"] is True
    assert all(result["checks"].values())
    method = result["method_independence"]
    assert method["production_denominator_module_imported"] is False
    assert method["old_D34_corpus_or_compiler_read"] is False
    assert method["T2_loader_module_imported"] is False
    assert method["S16_S17_synthesis_module_imported"] is False
    source = Path(
        "scripts/audit_h9_14_d34_cut_s16_s17_hostile.py"
    ).read_text(encoding="utf-8")
    forbidden_imports = (
        "import scripts.certify_h9_14_d34_cut_denominator",
        "import scripts.certify_h9_14_t1_t2_d34_unfolding_census",
        "import scripts.certify_h9_14_t2_surface_unfolding_pilot",
        "import scripts.certify_h9_14_s16_strong_transfer_corollary",
    )
    assert not any(statement in source for statement in forbidden_imports)


def test_custom_graph_groups_and_double_orbit_census() -> None:
    rows = supports()
    assert rows["S16"]["graph"]["automorphism_order"] == 48
    assert rows["S17"]["graph"]["automorphism_order"] == 8
    assert all(row["graph"]["group_closure_checked"] for row in rows.values())
    assert all(row["graph"]["edge_action_faithful"] for row in rows.values())
    assert rows["S16"]["all_eight_m4_labellings"]["labelled_count"] == 12870
    assert rows["S17"]["all_eight_m4_labellings"]["labelled_count"] == 24310
    assert rows["S16"]["all_eight_m4_labellings"]["explicit_orbit_count"] == 376
    assert rows["S17"]["all_eight_m4_labellings"]["explicit_orbit_count"] == 3234
    for row in rows.values():
        raw = row["all_eight_m4_labellings"]
        assert raw["explicit_orbit_count"] == raw["burnside_orbit_count"]


def test_cut_iff_normalization_all_lanes_inertia_and_selectors() -> None:
    rows = supports()
    s16 = rows["S16"]["primitive_norm_2_4_normalization"]
    s17 = rows["S17"]["primitive_norm_2_4_normalization"]
    assert (s16["labelled_compatible_count"], s16["explicit_orbit_count"]) == (3270, 136)
    assert (s17["labelled_compatible_count"], s17["explicit_orbit_count"]) == (1519, 238)
    assert s16["iff_cut_delta_X"] and s17["iff_cut_delta_X"]
    s16_lanes = rows["S16"]["oriented_integral_gram_lanes"]
    s17_lanes = rows["S17"]["oriented_integral_gram_lanes"]
    assert (s16_lanes["count"], s17_lanes["count"]) == (272, 476)
    assert s16_lanes["exact_rank_inertia_histogram"] == {
        "rank10/inertia_1_9_4": 2,
        "rank12/inertia_2_10_2": 8,
        "rank13/inertia_2_11_1": 46,
        "rank14/inertia_2_12_0": 216,
    }
    assert s17_lanes["exact_rank_inertia_histogram"] == {
        "rank10/inertia_1_9_4": 2,
        "rank12/inertia_2_10_2": 6,
        "rank13/inertia_2_11_1": 36,
        "rank13/inertia_3_10_1": 4,
        "rank14/inertia_2_12_0": 326,
        "rank14/inertia_3_11_0": 102,
    }
    assert s16_lanes["short4_long10_selector_lane_count_before_inertia"] == 25
    assert s17_lanes["short4_long10_selector_lane_count_before_inertia"] == 40
    assert s16_lanes["short4_long10_selector_lorentzian_lane_count"] == 1
    assert s17_lanes["short4_long10_selector_lorentzian_lane_count"] == 1


def test_survivor_grams_orientations_and_public_lineage() -> None:
    for support_id, row in supports().items():
        survivor = row["unique_lorentzian_diagram_orbit"]
        assert survivor["fixed_by_full_support_automorphism_group"] is True
        lanes = survivor["scaled_dual_lanes"]
        assert [lane["gram_sha256"] for lane in lanes] == list(
            EXPECTED_SURVIVOR_GRAM_HASHES[support_id]
        )
        assert lanes[0]["normal_square_histogram"] == {"2": 4, "4": 10}
        assert lanes[1]["normal_square_histogram"] == {"2": 10, "4": 4}
        assert [lane["short4_long10_selector"] for lane in lanes] == [True, False]
    public = payload()["public_matrix_and_synthesis_cross_checks"]
    assert public["S16_generated_Gram_equals_public_T2_entry_for_entry"] is True
    assert public["S17_generated_Coxeter_relation_matrix_equals_public_T1_label_for_label"] is True


def test_independent_active_set_rays_and_full_vinberg_gate() -> None:
    expected = {
        "S16": (110, 90, 20, 683, {"0,2": 97, "1,1": 362, "2,0": 224}),
        "S17": (127, 112, 15, 707, {"0,2": 53, "1,1": 300, "2,0": 354}),
    }
    for support_id, row in supports().items():
        total, finite, ideal, edges, endpoint_histogram = expected[support_id]
        fresh = row["fresh_cone_and_vinberg"]
        cone = fresh["independent_active_set_cone"]
        vinberg = fresh["generic_exact_vinberg"]
        assert fresh["active_set_endpoints_equal_vinberg_endpoints"] is True
        assert (cone["extreme_ray_count"], cone["finite_ray_count"], cone["ideal_ray_count"]) == (
            total,
            finite,
            ideal,
        )
        assert cone["hyperideal_ray_count"] == 0
        assert cone["all_14_facets_irredundant"] is True
        assert vinberg["elliptic_rank8_edge_count"] == edges
        assert vinberg["edge_endpoint_type_histogram"] == endpoint_histogram
        assert vinberg["bad_edge_count"] == 0
        assert vinberg["bad_full_parabolic_neighborhood_count"] == 0
        assert vinberg["decision"] == "finite_volume_noncompact"


def test_synthesis_hashes_and_honest_transported_boundary() -> None:
    public = payload()["public_matrix_and_synthesis_cross_checks"]
    assert public["T2_package_hashes_recomputed_directly_from_public_matrix_arrays"] == (
        EXPECTED_T2_PACKAGE_HASHES
    )
    bridge = public["T2_marked_bridge_matrix"]
    assert bridge["published_string_matrix_sha256_recomputed"] == bridge[
        "published_matrix_sha256"
    ]
    assert bridge["determinant"] == 1 and bridge["inverse_integral"]
    assert public["T1_package_hash_manifest"] == EXPECTED_T1_PACKAGE_HASHES
    boundary = payload()["scope_boundary"]
    assert "16/5 stable-overclosure" in boundary["transported_not_reproved"]
    assert "does not pretend" in boundary["T1_manifest_limit"]
    assert boundary["new_abstract_polytope"] is False
