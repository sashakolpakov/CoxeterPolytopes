from __future__ import annotations

from functools import lru_cache
import hashlib
import json

from scripts.certify_enriques_type_i_explicit_lift_package import (
    OSCAR_ARCHIVE_SHA256,
    OSCAR_DATA,
    OSCAR_UPSTREAM_SHA256,
    OUTPUT,
    analyze,
)


ARCHIVE_SHA256 = "d9313e37687e355c9ca8f86d1ee11d74f0810a8c6906a0443f5377d2481fece7"


@lru_cache(maxsize=1)
def _payload() -> dict:
    return analyze()


def test_archive_replays_from_the_pinned_oscar_and_shimada_data() -> None:
    raw = OSCAR_DATA.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == OSCAR_ARCHIVE_SHA256
    assert hashlib.sha256(raw[:-1]).hexdigest() == OSCAR_UPSTREAM_SHA256
    assert raw.endswith(b"\n") and not raw.endswith(b"\n\n")
    assert _payload() == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert hashlib.sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    assert _payload()["all_checks_pass"] is True


def test_global_enriques_involution_and_rootless_anti_lattice() -> None:
    payload = _payload()
    surface = payload["surface_and_lattices"]
    assert surface["K3_cover_Picard_rank"] == 19
    assert surface["K3_cover_Picard_determinant"] == 4
    assert surface["period_domain_complex_dimension"] == 1
    assert surface["invariant_plus_anti_index_in_SX"] == 512
    anti = surface["anti_invariant_algebraic_lattice"]
    assert anti["isometry_type"] == "E8(-2)+<-4>"
    assert anti["rank"] == 9
    assert anti["determinant_absolute"] == 1024
    assert anti["smith_diagonal"] == [2] * 8 + [4]
    assert anti["primitive_in_SX_smith_diagonal"] == [1] * 9
    assert anti["negative_definite"] is True
    assert anti["all_norms_divisible_by"] == 4
    assert anti["contains_minus_two_root"] is False
    involution = payload["global_involution"]
    assert involution["tau_squared_is_identity"] is True
    assert involution["tau_is_Picard_isometry"] is True
    assert involution["projection_formula"] == "r=2a+cD"


def test_exact_twenty_four_lifts_and_saturated_difference_span() -> None:
    compiler = _payload()["finite_lift_compiler_calibration"]
    assert compiler["fiber_root_system"] == "E8+A1"
    assert compiler["fiber_root_count"] == 242
    assert compiler["fiber_enumeration"]["signed_vector_count"] == 242
    assert compiler["each_downstairs_wall_has_exactly_two_lifts"] is True
    assert len(compiler["lift_records"]) == 12
    assert all(record["both_squares"] == [-2, -2] for record in compiler["lift_records"])
    assert all(record["mutual_intersection"] == 0 for record in compiler["lift_records"])
    assert compiler["all_24_lifts_distinct"] is True
    assert compiler["all_24_pairwise_intersections_nonnegative"] is True
    assert compiler["all_24_unordered_pair_intersection_histogram"] == {
        "0": 246,
        "1": 28,
        "2": 2,
    }
    assert compiler["all_4096_sheet_selections_pass_pairwise_nonnegativity"] is True
    assert compiler["frozen_selected_c_smith_diagonal"] == [1] * 9
    assert compiler["frozen_differences_generate_full_saturated_anti_lattice"] is True


def test_small_ample_orientation_and_complete_L4_enumeration() -> None:
    gate = _payload()["ample_orientation_and_L4"]
    assert gate["downstairs_H_coordinates"] == [
        38,
        23,
        46,
        70,
        57,
        45,
        34,
        24,
        15,
        7,
    ]
    assert gate["downstairs_H_square"] == 70
    assert gate["K3_pullback_H_square"] == 140
    assert gate["wall_degrees"] == [1] * 10 + [8, 8]
    assert gate["root_regular"] is True
    enumeration = gate["bounded_enumeration"]
    assert enumeration["bound"] == 378
    assert enumeration["signed_vector_count"] == 472
    assert enumeration["sign_representative_count"] == 236
    assert gate["positive_root_degree_histogram"] == {
        "1": 20,
        "2": 20,
        "3": 24,
        "4": 28,
        "5": 36,
        "6": 44,
        "7": 64,
    }
    assert gate["degree_zero_root_count"] == 0
    assert len(gate["L4_records"]) == 24
    assert all(record["L4_pass"] for record in gate["L4_records"])
    assert all(record["negative_pairing_witness_count"] == 0 for record in gate["L4_records"])
    assert gate["all_24_lifts_irreducible_by_L4"] is True


def test_joined_downstairs_vinberg_gate_and_scope_are_exact() -> None:
    payload = _payload()
    gate = payload["downstairs_polytope_gate"]
    assert gate == {
        "compact": False,
        "complete_full_nef_wall_set": True,
        "facet_count": 12,
        "finite_vertex_count": 14,
        "finite_volume": True,
        "hyperbolic_dimension": 9,
        "hyperideal_vertex_count": 0,
        "ideal_vertex_count": 9,
        "literal_vinberg_affine_audit_count": 12,
        "literal_vinberg_elliptic_edge_count": 114,
        "shimada_record": "12A",
    }
    assert all(payload["good_lift_package_checks"].values())
    scope = payload["scope_and_lineage"]
    assert scope["new_abstract_coxeter_polytope_claimed"] is False
    assert scope["new_algebraic_surface_claimed"] is False
    assert "first exact worked positive package" in scope["project_gain"]
