from __future__ import annotations

from hashlib import sha256
import json

from scripts.certify_act_real_cubic_h10_congruence_obstruction import (
    OUTPUT,
    build_certificate,
)
from scripts.independent_act_real_cubic_h10_congruence_replay import replay


ARCHIVE_SHA256 = "c48878ca13d547e951b3796d1a6c3b07bdc24e43f661b06509a8cd29adea7d60"


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_frozen_archive_and_exact_live_replay() -> None:
    assert sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    assert build_certificate() == _archive()


def test_extension_gate_is_the_six_dimensional_ternary_code() -> None:
    gate = _archive()["exact_extension_congruence"]
    assert gate["accepted_residue_count"] == 3**6 == 729
    assert gate["ambient_residue_count"] == 3**11 == 177_147
    assert gate["direct_sum_decomposition"] == (
        "F3*e0 + T4 + T4 + <e9+t>"
    )
    assert gate["tetracode_weight_histogram"] == {"0": 1, "3": 8}
    assert gate["quadratic_identity_on_every_codeword"] == "Q(r)=x0^2 mod 3"
    assert gate["bilinear_identity_on_every_pair"] == (
        "B(r,s)=x0(r)*x0(s) mod 3"
    )


def test_global_root_and_finite_volume_obstruction() -> None:
    payload = _archive()
    root_gate = payload["primitive_reflection_root_gate"]
    assert root_gate["only_possible_positive_norms"] == [1, 2]
    assert root_gate["norm_two_excluded_by_code_identity"] is True
    assert root_gate["every_extendable_root_has_norm"] == 1
    assert root_gate["orthogonal_pair_count_in_the_complete_root_set"] == 0

    theorem = payload["global_finite_volume_theorem"]
    assert theorem["maximum_elliptic_rank"] == 1
    assert theorem["maximum_parabolic_Coxeter_rank"] == 1
    assert theorem["Vinberg_first_finite_volume_gate_passes"] is False
    assert theorem["compact_chamber_exists"] is False
    assert theorem["finite_volume_noncompact_chamber_exists"] is False
    assert theorem["new_hyperbolic_example"] is False
    assert theorem["Vinberg_search_status"] == (
        "TERMINATED_BY_GLOBAL_CONGRUENCE_OBSTRUCTION"
    )


def test_symmetry_compression_is_complete_through_height_three() -> None:
    payload = _archive()
    symmetry = payload["safe_symmetry_compression"]
    assert symmetry["tetracode_monomial_automorphism_order"] == 48
    assert symmetry["basepoint_and_oriented_initial_wall_stabilizer_order"] == 4608
    assert symmetry["wall_orbit_count_through_height_three"] == 10
    assert symmetry["wall_orbits_by_height"] == {
        "0": 1,
        "1": 1,
        "2": 2,
        "3": 6,
    }
    assert symmetry["every_computed_orbit_is_complete"] is True
    assert sum(
        record["orbit_size"] for record in symmetry["orbit_representatives"]
    ) == 580
    bounded = payload["bounded_vinberg_corroboration"]
    assert bounded["streamed_Gram_summary"]["zero_offdiagonal_pairs"] == 0
    assert bounded["bounded_exact_cdd_observation_2026_08_09"] == {
        "wall_count": 580,
        "extreme_rays": 1680,
        "finite_timelike": 0,
        "ideal_lightlike": 0,
        "hyperideal_spacelike": 1680,
        "scope": "partial chamber only; corroboration, not the proof",
    }


def test_scope_does_not_promote_the_lattice_obstruction_to_geometry() -> None:
    scope = _archive()["scope_and_nonclaims"]
    assert scope["proved_for_all_heights"] is True
    assert scope["proved_for_this_explicit_antiunitary_class_only"] is True
    assert scope["classification_of_all_ACT_antiunitary_classes"] is False
    assert scope["real_cubic_threefold_component_realization_proved"] is False
    assert scope["algebro_geometric_leverage_delta"] == 0
    assert scope["large_cdd_needed"] is False
    assert scope["height_four_cdd_permitted"] is False


def test_no_project_import_independent_arithmetic_replay() -> None:
    assert replay() == {
        "fixed_Gram": "2*diag(1^10,-1)",
        "literal_residue_count": 729,
        "linear_residue_count": 729,
        "literal_equals_five_equations": True,
        "Q_equals_x0_squared_on_code": True,
        "norm_two_roots": 0,
        "orthogonal_norm_one_root_pairs": 0,
        "Vinberg_first_finite_volume_gate": False,
        "status": "PASS",
    }
