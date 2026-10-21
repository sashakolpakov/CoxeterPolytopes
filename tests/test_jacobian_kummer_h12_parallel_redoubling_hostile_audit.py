from __future__ import annotations

from functools import lru_cache
import json

from scripts.audit_jacobian_kummer_h12_parallel_redoubling_hostile import (
    OUTPUT,
    build,
)


@lru_cache(maxsize=1)
def _payload():
    return build()


def test_hostile_archive_rebuilds_exactly() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["independence"]["production_module_imported"] is False
    assert payload["independence"]["other_project_module_imported"] is False


def test_complete_single_double_replay_matches_production() -> None:
    replay = _payload()["single_doubling_hostile_replay"]
    assert replay["doubling_wall_count"] == 124
    assert replay["non_doubling_wall_count"] == 48
    assert {
        count: len(walls)
        for count, walls in replay["single_double_facet_count_classes"].items()
    } == {"296": 12, "317": 16, "326": 96}
    assert [row["selected_wall_one_based"] for row in replay["representatives"]] == [
        31,
        61,
        77,
    ]
    assert [row["facet_count"] for row in replay["representatives"]] == [
        296,
        317,
        326,
    ]
    assert replay["all_production_representatives_match_exactly"] is True
    assert replay["seed_roots_excluding_W36_smith_diagonal"] == [1] * 13
    assert all(
        row["production_boundary_ordered_coordinates_exact_match"]
        and row["boundary_root_smith_diagonal"] == [1] * 13
        and row["integral_involutive_isometry_reflection_count"]
        == row["facet_count"]
        and row["strict_witness_count"] == row["facet_count"]
        and row["all_exact_checks_pass"]
        for row in replay["representatives"]
    )


def test_parallel_pair_and_all_integer_collision_classification() -> None:
    pair = _payload()["parallel_redoubling_pair"]
    collision = pair["all_integer_collision_audit"]
    assert pair["squares_and_pairing"] == [-4, -4, 4]
    assert pair["primitive_integral_common_null_direction"] == [
        -3,
        -3,
        -3,
        -3,
        -2,
        -2,
        5,
        6,
        5,
        6,
        5,
        6,
        8,
    ]
    assert pair["translation_characteristic_polynomial"] == "(lambda-1)^13"
    assert pair["rank_T_minus_I"] == 2
    assert pair["rank_T_minus_I_squared"] == 1
    assert pair["T_minus_I_cubed_is_zero"] is True
    assert pair["generated_reflection_group"] == "D_infinity"
    assert pair["side_facet_category_counts"] == {
        "orthogonal_to_W31_only": 32,
        "orthogonal_to_W36_only": 32,
        "orthogonal_to_both": 14,
        "orthogonal_to_neither": 92,
    }
    assert collision["bounded_exponent_scan_used"] is False
    assert collision["finite_solution_count"] == 220
    assert collision["finite_solution_histogram"] == {
        "T^z*S36;z=-1;sign=1": 32,
        "T^z*S36;z=0;sign=1": 32,
        "T^z;z=0;sign=1": 156,
    }
    assert collision["all_z_identity_count"] == 28
    assert collision["all_z_identity_histogram"] == {
        "T^z": 14,
        "T^z*S36": 14,
    }
    assert collision["negative_root_collision_count"] == 0
    assert collision["cross_seed_facet_collision_count"] == 0
    assert collision["endpoint_to_side_orbit_collision_count"] == 0
    assert collision["outer_endpoint_equations"] == {
        "T^z(W36)=plus_or_minus_W31": [],
        "T^z*S36(W31)=plus_or_minus_W31": [{"sign": -1, "z": -1}],
    }


def test_universal_facet_formula_and_full_lattice() -> None:
    family = _payload()["infinite_path_family"]
    assert family["facet_formula"] == "f(N)=124*N+48 for every integer N>=1"
    assert family["full_fixed_Picard_lattice_for_every_N"] is True
    assert family["coxeter_for_every_N"] is True
    assert family["finite_volume"] is True
    assert family["compact"] is False
    assert family["pairwise_nonisometric"] is True
    assert family["direct_exact_reconstruction_checks_N_1_through_12"] == [
        {"chamber_count_N": n, "facet_count": 124 * n + 48}
        for n in range(1, 13)
    ]
    decision = _payload()["decision"]
    assert decision["facet_formula_for_every_N_at_least_1"] == "124*N+48"
    assert decision["full_nef_or_surface_automorphism_claimed"] is False
    assert decision["worldwide_priority"] == "open"


def test_bounded_coxiter_replay_has_no_literal_family_header() -> None:
    corpus = _payload()["bounded_coxiter_corpus_replay"]
    assert corpus["zip_sha256"] == (
        "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
    )
    assert corpus["bundled_index_two_tar_gz_sha256"] == (
        "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
    )
    assert corpus["ordinary"]["coxiter_file_count"] == 1214
    assert corpus["ordinary"]["h12_facet_histogram"] == {
        "13": 6,
        "14": 3,
        "16": 1,
        "18": 1,
    }
    assert corpus["index_two"]["coxiter_file_count"] == 2959
    assert corpus["index_two"]["h12_facet_histogram"] == {
        "18": 1,
        "20": 1,
        "25": 1,
        "35": 1,
    }
    assert corpus["literal_seed_or_descendant_match_absent_by_header"] is True
    assert corpus["does_not_claim_worldwide_priority"] is True
