from __future__ import annotations

from functools import lru_cache
import hashlib
import json

from scripts.audit_enriques_type_ii_downstairs_anti_gram_census_hostile import (
    OUTPUT,
    analyze,
)


ARCHIVE_SHA256 = "a292b28be9352e2d21841d0ae97f1586308c42b29c53fbc08abf72acf2df712f"


@lru_cache(maxsize=1)
def _payload() -> dict:
    return analyze()


def test_hostile_archive_replays_and_freezes_type_ii_compiler() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert hashlib.sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    assert payload["verdict"] == "PASS/FREEZE"
    assert payload["all_checks_pass"] is True
    assert "32,768 anti-Grams" in payload["frozen_conclusion"]


def test_independent_parse_switching_gauge_and_sturm_inertia() -> None:
    parsed = _payload()["independent_downstairs_parse"]
    assert parsed["parser"].endswith("no project GAP parser")
    assert parsed["wall_count"] == 12
    assert parsed["m1_edge_count"] == 15
    assert parsed["wall_rank_mod_2"] == 9
    assert parsed["binary_left_kernel_size"] == 8

    quotient = _payload()["independent_switching_census"]
    assert quotient["different_from_production_negative_BFS_gauge"] is True
    assert quotient["raw_assignment_count"] == 32768
    assert quotient["orbit_count"] == 16
    assert quotient["orbit_size_histogram"] == {"2048": 16}
    assert quotient["canonical_key_set_sha256"] == (
        "fa590a8af74202c4b03297b7e78362428cdea02a9437022a50aadb60c02bc6d3"
    )

    inertia = _payload()["independent_exact_inertia"]
    assert inertia["engine"].endswith("Sturm root counts")
    assert inertia["raw_assignment_histogram"] == {
        "p0_n9_z3": 4096,
        "p1_n10_z1": 12288,
        "p1_n9_z2": 2048,
        "p2_n10_z0": 14336,
    }
    assert inertia["switch_orbit_histogram"] == {
        "p0_n9_z3": 2,
        "p1_n10_z1": 6,
        "p1_n9_z2": 1,
        "p2_n10_z0": 7,
    }


def test_one_stage_full_half_sum_lattices_are_D9_and_cubic() -> None:
    lattices = _payload()["full_half_sum_lattices"]
    survivor = lattices["survivor"]
    assert survivor["different_from_production_two_stage_hnf"] is True
    assert survivor["reverse_lex_delta_basis_indices_one_based"] == [
        12,
        11,
        10,
        9,
        8,
        7,
        5,
        4,
        3,
    ]
    assert survivor["absolute_determinant"] == 2048
    assert survivor["smith_diagonal"] == [2] * 8 + [8]
    classification = survivor["classification"]
    assert classification["isometry_class"] == "D9(-2)"
    assert classification["rootless_through_norm_two"] is True
    assert classification["signed_norm_four_vector_count"] == 144
    assert classification["simple_root_basis_determinant_absolute"] == 1
    assert classification["dynkin_arm_lengths"] == [1, 1, 6]
    assert classification["minimum_enumeration"]["engine"].startswith(
        "independent exact rational LDL"
    )

    eliminated = lattices["eliminated"]
    assert eliminated["fractional_coset_count_mod_Z9"] == 2
    assert eliminated["absolute_determinant"] == 512
    assert eliminated["smith_diagonal"] == [2] * 9
    classification = eliminated["classification"]
    assert classification["isometry_class"] == "<-2>^9"
    assert classification["signed_minus_two_vector_count"] == 18
    assert classification["minus_two_vectors_modulo_sign_count"] == 9
    assert classification["unimodular_root_basis_determinant_absolute"] == 1


def test_all_36_witnesses_and_all_sheets_replay_with_binary_gauge_warning() -> None:
    ternary = _payload()["ternary_root_audit"]
    assert ternary["survivor_witness_count_modulo_sign"] == 0
    assert ternary["eliminated_witness_count_modulo_sign"] == 36
    assert ternary["minimum_support_size"] == 9
    assert ternary["distinct_lattice_root_pairs"] == 9
    assert ternary["coefficient_presentations_per_root_pair"] == {"4": 9}
    assert ternary["production_gauge_witness_hash_replayed"] == (
        "6c7a35b3bafc46c5569742bcbe83c915681e7ef641175e0112a350654e357415"
    )
    binary = ternary["binary_shortcut_gauge_dependence"]
    assert binary["production_negative_BFS_survivor_count"] == 0
    assert binary["production_negative_BFS_eliminated_count"] == 0
    assert binary["hostile_positive_DFS_survivor_count"] == 0
    assert binary["hostile_positive_DFS_eliminated_count"] == 3
    assert "depend on sheet labels" in binary["conclusion"]

    sheets = _payload()["sheet_invariance"]
    assert sheets["sheet_switches_checked"] == 2048
    assert sheets["witnesses_per_switch"] == 36
    assert sheets["direct_parity_and_norm_checks"] == 73728
    assert sheets[
        "survivor_rootlessness_invariant_by_the_same_integral_lattice_isometry"
    ] is True


def test_primary_provenance_is_pinned_and_scope_is_not_overclaimed() -> None:
    provenance = _payload()["geometric_provenance"]
    brandhorst = provenance["brandhorst_shimada"]
    assert brandhorst["hostile_local_text_extraction_checked"] is True
    assert "D9(-2)" in brandhorst["checked_statement"]
    assert "No.184" in brandhorst["checked_statement"]
    assert "12B" in brandhorst["checked_statement"]
    martin = provenance["martin"]
    assert martin["source_tex_checked_in_memory"] is True
    assert martin["no_lift_gram_imported"] is True
    assert "beta != 0,-64" in martin["checked_statement"]

    scope = _payload()["scope"]
    assert scope["new_surface_claimed"] is False
    assert scope["new_abstract_hyperbolic_Coxeter_polytope_claimed"] is False
    assert scope["period_or_ample_cone_reconstructed_from_downstairs_alone"] is False
