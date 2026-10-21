"""Independent hostile tests for the compact (2,4,6) Property-4 no-go."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.audit_compact_246_property4_counterexample_independent_hostile import (
    AUDIT_ARCHIVE,
    AUDIT_REPORT,
    OWNER_HASHES,
    anisotropy_descent,
    build_certificate,
    digest,
    encoded_certificate,
    exact_minimal_trace_obstruction,
    owner_hash_audit,
    rational_full_nef_obstruction,
    reconstruct_raw_packet,
    reflection_area_maximality,
)


AUDITOR = (
    ROOT
    / "scripts/audit_compact_246_property4_counterexample_independent_hostile.py"
)


def test_owner_bytes_are_frozen_and_replay_is_independent() -> None:
    assert owner_hash_audit() == OWNER_HASHES
    source = AUDITOR.read_text(encoding="utf-8")
    assert "certify_compact_246_property4_counterexample import" not in source
    assert "import sympy" not in source


def test_raw_gram_angles_reflections_and_field_reconstruct() -> None:
    raw = reconstruct_raw_packet()
    assert raw["wall_gram"] == [[-12, 6, 0], [6, -6, 3], [0, 3, -2]]
    assert raw["determinant"] == 36
    assert raw["LDL_pivots"] == [-12, -3, 1]
    assert raw["signature"] == [1, 2]
    assert raw["normalized_cosine_squares"] == ["1/2", "3/4", "0"]
    assert raw["edge_labels"] == [4, 6, 2]
    assert raw["principal_pair_determinants"] == [36, 3, 24]
    assert raw["all_vertices_finite"] is True
    assert raw["area_over_pi"] == "1/12"
    assert raw["cartan_coefficients"] == {
        "C_12": "-1",
        "C_21": "-2",
        "C_23": "-1",
        "C_32": "-3",
    }
    assert raw["cycle_products"] == [2, 3]
    assert raw["minimal_vinberg_field"] == "Q"
    assert raw["reflection_matrices"] == [
        [[-1, 1, 0], [0, 1, 0], [0, 0, 1]],
        [[1, 0, 0], [2, -1, 1], [0, 0, 1]],
        [[1, 0, 0], [0, 1, 0], [0, 3, -1]],
    ]
    assert raw["integral_reflection_host"] is True


def test_mod_three_descent_proves_rational_anisotropy() -> None:
    descent = anisotropy_descent()
    assert descent["descent_prime"] == 3
    assert descent["q_mod_3_zero_count"] == 9
    assert descent["q_mod_3_forces_z_zero"] is True
    assert descent["residual_binary_form"] == "2x^2+xy+y^2"
    assert descent["residual_binary_zeros_mod_3"] == [[0, 0]]
    assert descent["primitive_integer_zero_would_be_divisible_by_3"] is True
    assert descent["rational_form_anisotropic"] is True


def test_absolute_maximality_has_one_candidate_and_no_target_double() -> None:
    maximality = reflection_area_maximality()
    assert maximality["source_area_over_pi"] == "1/12"
    assert maximality["proper_overgroup_area_upper_bound_over_pi"] == "1/24"
    assert maximality["compact_polygon_s_ge_4_area_floor_over_pi"] == "1/6"
    assert maximality["overgroup_chamber_must_be_triangle"] is True
    assert maximality["compact_triangle_area_floor_over_pi"] == "1/42"
    assert maximality["derived_index_range"] == [2, 3]
    assert maximality["solutions_by_index"] == {
        "2": [[2, 3, 8]],
        "3": [],
    }
    assert [case["union"] for case in maximality["doubling_cases"]] == [
        "triangle (3,3,4)",
        "triangle angles (pi/8,pi/8,2pi/3)",
        "quadrilateral",
    ]
    assert not any(case["equals_target"] for case in maximality["doubling_cases"])
    assert maximality["proper_discrete_reflection_overgroup_exists"] is False
    assert maximality["absolute_reflection_maximality"] is True


def test_both_rational_exact_full_nef_branches_are_empty_twice_over() -> None:
    obstruction = rational_full_nef_obstruction(anisotropy_descent())
    assert obstruction["target_hyperbolic_dimension"] == 2
    assert obstruction["literal_full_nef_forces_picard_rank"] == 3
    assert obstruction["rank_three_rational_NS"] == "I_(1,2)"
    assert obstruction["I_1_2_isotropic_vector"] == [1, 1, 0]
    assert obstruction["I_1_2_isotropic_vector_square"] == 0
    assert obstruction["target_trace_form_anisotropic"] is True
    assert obstruction["full_rank_similarity_preserves_isotropy"] is True
    assert obstruction["rational_HOST_exists"] is False
    assert obstruction[
        "every_rational_surface_rho_ge_2_has_nef_isotropic_ray"
    ] is True
    assert obstruction["compact_target_has_ideal_ray"] is False
    assert obstruction["direct_rational_exact_full_nef"] == "EMPTY"
    assert (
        obstruction["certified_weak_log_rational_exact_full_nef"] == "EMPTY"
    )


def test_exact_minimal_trace_branch_has_only_labels_two_and_three() -> None:
    obstruction = exact_minimal_trace_obstruction()
    assert obstruction["trace_rank"] == 3
    assert obstruction["literal_full_nef_forces_picard_rank"] == 3
    assert obstruction["nontrivial_contraction_endpoint_rank_at_most"] == 2
    assert obstruction["rank_three_trace_cannot_descend_nontrivially"] is True
    assert obstruction["exact_carrier_is_already_minimal"] is True
    assert obstruction["trace_space_equals_NS_over_Q"] is True
    assert (
        obstruction["canonical_orthogonality_forces_K_numerically_zero"]
        is True
    )
    assert obstruction["facet_curves_are_smooth_rational_minus_two"] is True
    assert obstruction["allowed_finite_angle_rows"] == [
        {
            "intersection": 0,
            "cosine_square": "0",
            "finite_Coxeter_label": 2,
        },
        {
            "intersection": 1,
            "cosine_square": "1/4",
            "finite_Coxeter_label": 3,
        },
    ]
    assert obstruction["allowed_finite_Coxeter_labels"] == [2, 3]
    assert obstruction["forbidden_target_labels"] == [4, 6]
    assert obstruction["exact_return_minimalizable_full_trace"] == "EMPTY"
    assert obstruction["bare_carrier_relative_trace_injection_excluded"] is False


def test_four_property_scope_refutes_only_current_outer_surjectivity() -> None:
    payload = build_certificate()
    scope = payload["four_property_scope"]
    assert scope["current_three_species_all_empty"] is True
    assert scope["conditional_AG_dimension_theorem_refuted"] is False
    assert scope["all_conceivable_AG_carriers_excluded"] is False
    assert scope["proper_invariant_or_eigen_nef_slice_for_k_Q_excluded"] is False
    assert scope["bare_carrier_relative_third_row_excluded"] is False
    assert scope["projected_larger_Picard_return_excluded"] is False
    assert scope["higher_dimensional_variety_return_excluded"] is False
    assert scope["positive_characteristic_return_excluded"] is False
    assert payload["verdict"] == "PASS"


def test_report_states_exact_return_scope_and_all_escape_clauses() -> None:
    text = AUDIT_REPORT.read_text(encoding="utf-8")
    required = (
        "PASS at the owner's exact scope",
        "Absolute reflection maximality",
        "\\operatorname{area}(P)",
        "N=2,3",
        "(2,3,8)",
        "Both rational exact full-nef species are empty",
        "The exact-return minimalizable full-trace species is empty",
        "word **exact-return** does essential work",
        "Exact four-property verdict",
        "conditional AG dimension bounds",
        "proper invariant/eigen-nef slice",
        "bare carrier-relative third row",
        "did not edit any owner byte",
    )
    for needle in required:
        assert needle in text
    assert text.count("\\[") == text.count("\\]")
    assert "../scripts/audit_compact_246" in text
    assert "../data/enumerations/compact_246" in text
    assert "../tests/test_compact_246" in text


def test_archive_is_exact_replay() -> None:
    archived = json.loads(AUDIT_ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert AUDIT_ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_cli_check() -> None:
    subprocess.run(
        [sys.executable, str(AUDITOR), "--check-output"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
