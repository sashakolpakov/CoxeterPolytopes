"""Independent hostile tests for the revised four-property completion ledger."""

from __future__ import annotations

from scripts.audit_canonical_chamber_four_property_completion_independent_hostile import (
    PINS,
    ROOT,
    anticanonical_record,
    bianticanonical_record,
    bounded_spanning_root_record,
    build_audit,
    compact_field_record,
    compact_full_picard_compiler_record,
    conditional_dimension_record,
    corrected_genus_one_record,
    fixed_jacobian_index_record,
    fixed_data_outer_routing_record,
    noncompact_cusp_host_record,
    noncompact_field_record,
    outer_gate_record,
    physical_star_host_record,
    projected_delta246_record,
    projected_scaling_record,
    same_host_composition_record,
    sha256_file,
)


AUDIT_REPORT = (
    ROOT
    / "results/canonical_chamber_four_property_completion_independent_hostile_audit.md"
)
AUDIT_REPORT_SHA256 = (
    "8ee64106a10655c9126990365e78e1546f7d50806a31ebab7b5936ee6baf80a0"
)


def test_all_seventy_evidence_and_synthesis_artifacts_are_byte_pinned() -> None:
    assert len(PINS) == 70
    for relative, expected in PINS.values():
        assert sha256_file(ROOT / relative) == expected
    assert sha256_file(AUDIT_REPORT) == AUDIT_REPORT_SHA256


def test_dimension_transfer_stays_conditional() -> None:
    assert conditional_dimension_record() == {
        "direct_rational_full_NS_n_max": 11,
        "weak_log_rational_n_max": 25,
        "full_trace_complex_K3": "d(n+1)<=20",
        "nonrational_K3_n_max": 9,
        "conditional_global_n_max": 25,
        "exhaustive_property_4_entry_required": True,
    }


def test_primitive_anticanonical_reflection_gate_is_exact() -> None:
    record = anticanonical_record()
    assert record["unimodular_divisibility"] == 1
    assert record["dimension_bound"] == 11
    assert record["P11_sharp"] is True
    assert [
        (row["dimension"], row["Picard_rank"], row["F_square"], row["arithmetic_genus"])
        for row in record["allowed_negative_cases"]
    ] == [(10, 11, -1, 1), (11, 12, -2, 1)]


def test_bianticanonical_same_frame_has_opposite_section_algebras() -> None:
    record = bianticanonical_record()
    assert (
        record["K_square"],
        record["F_square"],
        record["B_square"],
        record["K_dot_B"],
        record["B_genus"],
    ) == (-1, -1, -4, 2, 0)
    assert record["cubic_matrix_shape"] == [10, 10]
    assert record["sextic_double_jet_matrix_shape"] == [30, 28]
    assert record["elliptic_branch"] == {
        "h0_minus_K": 1,
        "h0_minus_2K": 1,
        "unique_divisor": "nonreduced 2F",
    }
    assert record["coble_branch"] == {
        "h0_minus_K": 0,
        "h0_minus_2K": 1,
        "unique_divisor": "geometrically integral rational (-4)",
    }
    assert record["same_full_nef_cone_claimed"] is False
    assert record["prime_sections_alone_imply_exact_nef"] is False
    assert record["full_handoff_hypotheses_also_required"] is True
    assert record["universal_prime_locus_nonempty"] is False


def test_corrected_genus_one_packet_is_rank_blind_but_abel_cutoff_is_not() -> None:
    records = [corrected_genus_one_record(n) for n in range(12, 65)]
    assert {
        (row["C_square"], row["K_dot_C"], row["C_genus"], row["C_dot_Z"])
        for row in records
    } == {(-2, 2, 1, 0)}
    assert records[0]["Z_square"] == -1
    assert records[-1]["Z_square"] == 11 - 64
    boundary = [corrected_genus_one_record(n) for n in range(12, 20)]
    assert [row["full_nef_carrier"] for row in boundary] == [
        True,
        True,
        True,
        True,
        True,
        False,
        False,
        False,
    ]
    assert [row["quartic_exceptional_count"] for row in boundary[-3:]] == [17, 18, 18]
    assert [row["A_square"] for row in boundary[-3:]] == [-1, -2, -2]
    assert {row["A_genus"] for row in boundary[-3:]} == {3}
    assert {row["C_dot_A"] for row in boundary[-3:]} == {1}


def test_noncompact_full_cusp_packet_forces_the_rational_field() -> None:
    assert noncompact_field_record() == {
        "hyperbolic_dimension_min": 2,
        "full_incident_cusp_packet_required": True,
        "bounded_horosphere_section": True,
        "common_fixed_space_dimension": 1,
        "k_rational_null_line": True,
        "definite_nonidentity_place_contradiction": True,
        "coefficient_field": "Q",
        "K3_needed_for_field_collapse": False,
        "full_trace_K3_n_max": 19,
    }


def test_compact_narrow_facet_degree_and_height_table() -> None:
    record = compact_field_record()
    assert record["m_min"] == 3
    assert record["narrow_pairing_bound_strict"] == 14
    assert record["star_spans"] is True
    assert record["lorentz_cycle_field_equals_coefficient_field"] is True
    assert record["separate_minimal_field_hypothesis"] is False
    assert record["cycles_integral_by_reflection_word_trace"] is True
    assert [row["degree_cap"] for row in record["rows"]] == [6, 5, 4, 3, 2, 2, 2, 2]
    assert [row["discriminant_exponent_cap"] for row in record["rows"]] == [
        30,
        20,
        12,
        6,
        2,
        2,
        2,
        2,
    ]
    assert record["degree_cap"] == 6
    assert record["m_cap_for_nonrational_fields"] == 10


def test_compact_primitive_grid_and_discriminant_constant() -> None:
    record = compact_field_record()
    assert record["minimal_generator_count_cap"] == 2
    assert record["primitive_grid_q_max"] == 15
    assert record["primitive_grid_size"] == 16
    assert record["uniform_cycle_height"] == 289_254_654_976
    assert record["primitive_conjugate_bound"] == 4_628_074_479_616
    assert record["conjugate_difference_bound"] == 9_256_148_959_232
    assert record["discriminant_base"] == 32 * 14**10
    assert record["discriminant_exponent"] == 30
    assert record["discriminant_bound"] == (32 * 14**10) ** 30


def test_physical_star_finitizes_only_the_native_same_span_trace_host() -> None:
    record = physical_star_host_record()
    assert record["formal_star_hypothesis"] == "every minimal-aperture-star normal"
    assert record["same_vector_in_trace_module"] is True
    assert record["physical_curve_square"] == -2
    assert record["placewise_energy_sum"] == 2
    assert record["distinguished_off_diagonal_bound_strict"] == 14
    assert record["definite_place_off_diagonal_bound"] == 2
    assert record["field_gram_entries_in_order_codifferent"] is True
    assert record["root_frame_sandwich"] == "F subset M subset F#"
    assert record["finite_same_span_trace_hosts"] is True
    assert record["finite_ambient_Picard_extensions"] is False
    assert record["projected_wall_line_suffices"] is False
    assert record["property4_proved"] is False


def test_same_host_composition_retains_exact_return_and_native_maximality() -> None:
    record = same_host_composition_record()
    assert record["complete_oriented_native_host_required"] is True
    assert record["exact_native_return_required"] is True
    assert record["same_host_k_reflection_maximality_required"] is True
    assert record["facet_group_equals_complete_native_reflection_group"] is True
    assert record["one_unmarked_normalized_chamber_orbit_per_host"] is True
    assert record["finite_index_hypothesis_required"] is False
    assert record["ordinary_rank_one_K3_root_maximality_suffices_for_d_gt_1"] is False
    assert record["complement_period_extendable_maximality_suffices_without_more_data"] is False
    assert record["counts_surfaces_or_periods"] is False


def test_projected_physical_roots_leave_the_compact_trace_scale_unbounded() -> None:
    rows = [projected_scaling_record(scale) for scale in range(2, 40)]
    assert all(row["visible_root_half_norm_cases_in_checked_window"] == [1] for row in rows)
    assert all(row["visible_root_trace_square"] == -2 for row in rows)
    assert all(row["trace_contains_minus_two_vector"] is False for row in rows)
    assert all(row["physical_facet_roots_lie_outside_trace"] is True for row in rows)
    assert all(row["fixed_field_order_cycles_and_chamber"] is True for row in rows)
    assert all(row["infinitely_many_abstract_chambers_claimed"] is False for row in rows)
    assert [row["glue_index"] for row in rows] == [scale**3 for scale in range(2, 40)]
    assert [row["ambient_absolute_discriminant"] for row in rows] == [
        64 * scale**3 * (scale - 1) ** 3 for scale in range(2, 40)
    ]


def test_compact_full_picard_compiler_is_finite_iff_empty_not_property_four() -> None:
    record = compact_full_picard_compiler_record()
    assert record["one_pinned_exact_compact_nonrational_packet"] is True
    assert record["field_degree_gt_one"] is True
    assert record["hyperbolic_dimension_min"] == 2
    assert record["exact_algebraic_wall_lines_required"] is True
    assert record["finite_marked_package_list"] is True
    assert record["may_return_EMPTY"] is True
    assert record["gates"] == ["PHYS", "EMB", "ROOT"]
    assert record["full_Picard_means_NS_equals_trace_lattice"] is True
    assert record["ample_eigenframe_required"] is True
    assert record["p_in_Amp_intersect_E"] is True
    assert record["projective_p_in_polytope_interior"] is True
    assert record["returned_chamber_is_exact_eigen_nef_section"] is True
    assert record["returned_chamber_is_literal_full_nef_in_general"] is False
    assert record["FV13_applicable_in_general"] is False
    assert record["raw_maximality_replaces_ROOT"] is False
    assert record["zero_projection_and_facet_witness_tests_remain"] is True
    assert record["universal_property4_nonemptiness"] is False


def test_fv13_closes_literal_full_nef_hosts_without_closing_property_four() -> None:
    rows = [bounded_spanning_root_record(rank) for rank in range(3, 21)]
    assert [row["Picard_rank"] for row in rows] == list(range(3, 21))
    assert all(row["finite_sided"] and row["finite_volume"] for row in rows)
    assert all(row["literal_full_Picard_full_nef"] for row in rows)
    assert all(row["selected_actual_facet_roots"] == row["Picard_rank"] for row in rows)
    assert {row["root_square"] for row in rows} == {-2}
    assert {row["off_diagonal_intersection_bound"] for row in rows} == {13}
    assert {row["root_frame_sandwich"] for row in rows} == {
        "R subset NS subset R*"
    }
    assert all(row["finite_Picard_hosts"] for row in rows)
    assert all(row["one_complete_root_Weyl_orbit_per_oriented_host"] for row in rows)
    assert all(row["reflection_maximality_required"] is False for row in rows)
    assert all(row["proper_eigen_nef_section_covered"] is False for row in rows)
    assert all(row["property4_proved"] is False for row in rows)
    rank_18 = bounded_spanning_root_record(18)
    assert rank_18["row_norm_squared"] == 2877
    assert rank_18["hadamard_determinant_bound_squared"] == 2877**18


def test_noncompact_affine_packet_is_local_not_a_counterexample_to_fv13() -> None:
    rows = [noncompact_cusp_host_record(d) for d in range(2, 40)]
    assert {tuple(map(tuple, row["affine_A1_gram"])) for row in rows} == {
        ((-2, 2), (2, -2))
    }
    assert [row["primitive_fibre_divisibility"] for row in rows] == list(range(2, 40))
    assert [row["absolute_discriminant"] for row in rows] == [
        2 * d**2 for d in range(2, 40)
    ]
    assert all(row["transverse_root_square"] == -2 for row in rows)
    assert [row["fibre_transverse_pairing"] for row in rows] == list(range(2, 40))
    assert all(row["finite_volume_claimed"] is False for row in rows)
    assert all(row["reflection_maximality_claimed"] is False for row in rows)
    assert all(row["family_proves_finite_host_theorem"] is False for row in rows)
    assert all(row["global_fv13_full_nef_theorem_proved"] is True for row in rows)
    assert all(row["bounded_spanning_root_lemma_proved"] is True for row in rows)
    assert all(
        row["bounded_spanning_root_lemma_applies_without_finite_volume"] is False
        for row in rows
    )


def test_fixed_jacobian_local_data_are_index_blind_with_correct_fv13_corollary() -> None:
    rows = [fixed_jacobian_index_record(index) for index in range(2, 40)]
    assert all(row["projective_K3_torsor"] for row in rows)
    assert {row["Picard_rank"] for row in rows} == {18}
    assert [row["absolute_discriminant"] for row in rows] == [
        index**2 for index in range(2, 40)
    ]
    assert [row["fibre_divisibility"] for row in rows] == list(range(2, 40))
    assert {row["Kodaira_fibres"] for row in rows} == {"2II*+4I1"}
    assert {row["Euler_number"] for row in rows} == {24}
    assert {row["Mordell_Weil_rank"] for row in rows} == {0}
    assert {row["discriminant_length"] for row in rows} == {2}
    assert {row["physical_affine_E8_packets"] for row in rows} == {2}
    assert all(row["finite_volume_complete_chamber_claimed"] is False for row in rows)
    assert all(row["cusp_local_data_bound_index"] is False for row in rows)
    assert {
        row["conditional_finite_volume_index_squared_bound"] for row in rows
    } == {2877**9}
    assert all(row["q_le_390_requires_cusp_adapted_transverse_wall"] for row in rows)
    assert all(row["q_le_390_automatic_from_FV13"] is False for row in rows)


def test_outer_gate_distinguishes_false_three_way_from_open_fourth_category() -> None:
    outer = outer_gate_record()
    assert outer["weak_log_species_has_finite_iff_EMPTY_compiler"] is True
    assert outer["weak_log_EMPTY_implies_all_species_empty"] is False
    assert outer["compact_256_same_vector_K3_PHYS_empty"] is True
    assert outer["compact_256_excludes_projected_K3"] is False
    assert outer["compact_246_absolutely_reflection_maximal"] is True
    assert outer["old_exact_three_branch_property4_universal"] is False
    assert outer["bare_carrier_relative_third_row_disproved"] is False
    assert outer["projected_K3_fourth_category_universal"] is False
    assert outer["all_surface_or_variety_carriers_excluded"] is False
    assert outer["P17_P19_declared_hosts_exhausted"] is True
    assert outer["P17_P19_absolute_arithmetic_maximality_certified"] is False


def test_projected_delta246_has_exact_revised_property_map_and_scope() -> None:
    record = projected_delta246_record()
    assert record["ambient_determinant"] == -393336
    assert record["ambient_signature"] == [1, 5]
    assert record["projector_idempotent"] is True
    assert record["projector_self_adjoint"] is True
    assert record["ellipsoid_reduces_to_zero_and_signed_basis"] is True
    assert record["P1_full_raw_module_primitive_in_NS"] is True
    assert record["old_literal_rational_P2"] is False
    assert record["revised_exact_projected_eigen_nef_P2"] is True
    assert record["P3_minimal_K3_K_numerically_zero"] is True
    assert record["P4_this_input"] is True
    assert record["P4_universal_projected_category"] is False
    assert record["raw_reflections_integral_on_NS"] is False
    assert record["raw_reflections_rational_on_NS_Q"] is True
    assert record["ambient_root_reflections_are_raw_reflections"] is False
    assert record["finite_projected_host_control"] is False
    assert record["K3_projected_category_n_max"] == 19
    assert record["nonrational_members_require_full_trace_admissibility"] is True


def test_fixed_data_closure_and_outer_routing_quantifiers_are_exact() -> None:
    record = fixed_data_outer_routing_record()
    assert len(record["documents_pinned_individually"]) == 6
    assert all(
        all(checks.values()) for checks in record["links_by_document"].values()
    )
    finite_pair = record["finite_pair"]
    assert finite_pair["finite_native_pairs_plus_reflection_maximality"] is True
    assert finite_pair["maximal_discrete_or_arithmetic_scope"] is True
    assert finite_pair["fixed_field_rank_similarity_classification_imported"] is True
    assert finite_pair["maximality_forces_W_equals_R"] is True
    assert finite_pair["one_complete_chamber_orbit_per_oriented_similarity_host"] is True
    assert finite_pair["global_AG_only_ABSW"] is False
    x3 = record["x3"]
    assert x3["one_counterexample_to_all_finite_volume_universality"] is True
    assert x3["obstructs_every_finite_volume_input"] is False
    assert x3["absolute_reflection_maximality_proved"] is False
    assert x3["maximal_only_no_go_proved"] is False
    assert x3["mixed_divisor_K3_positive"] is True
    h21 = record["h21"]
    assert h21["absolute_reflection_maximality"] is True
    assert h21["current_complex_stack_empty"] is True
    assert h21["characteristic_two_mixed_curve_Cremona_positive"] is True
    assert h21["all_AG_carriers_empty"] is False
    assert (
        record["projected_ambient_host_finiteness_needed_for_abstract_count"]
        is False
    )
    assert record["conditional_ABSW_recovery"] is True
    assert record["global_AG_only_ABSW_recovered"] is False


def test_final_scope_adds_conditional_fixed_data_closure_without_global_claim() -> None:
    audit = build_audit()
    assert audit["verdict"] == (
        "PASS_FIXED_PAIR_CLASSIFICATION_IMPORT_OUTER_ROUTING_NONEXHAUSTIVE"
    )
    assert audit["conditional_ABSW_recovery"] is True
    assert audit["global_AG_only_ABSW_recovered"] is False
    assert audit["global_ABSW_recovered"] is False
    assert audit["four_properties"]["raw_exhaustive_return"] == (
        "false for the current exact three-branch disjunction"
    )
    assert audit["four_properties"][
        "universal_exact_module_projected_root_complex_K3"
    ] == "false over the all-finite-volume proposal"
    assert audit["four_properties"][
        "mixed_divisor_action_or_other_revised_category"
    ] == (
        "open class-wide"
    )
    assert len(audit["narrow_scope_repairs"]) == 22
    fields = audit["field"]
    assert fields["all_finite_volume_exact_admissible_full_trace_K3_n_ge_2_closed"] is True
    assert fields["dimension_one_exception"] is True
    assert fields["mixed_height_needed_for_exact_finite_volume_K3_field_set"] is False
    assert fields["partial_packet_fallback_theorems_retained"] is True
    assert fields["carrier_construction_proved_by_field_theorems"] is False
    assert fields["property4_proved_by_field_theorems"] is False
    assert fields["finite_host_or_chamber_set_proved_by_field_theorems"] is False
    compact_hosts = audit["compact_host_and_chamber"]
    assert compact_hosts["finite_trace_hosts_under_physical_star"] is True
    assert compact_hosts[
        "finite_chambers_also_require_exact_return_and_same_host_maximality"
    ] is True
    assert compact_hosts["finite_ambient_Picard_hosts_claimed"] is False
    compiler = audit["compact_pinned_compiler"]
    assert compiler["gates"] == ["PHYS", "EMB", "ROOT"]
    assert compiler["may_return_EMPTY"] is True
    assert compiler["ample_eigenframe_required"] is True
    assert compiler["FV13_applicable_in_general"] is False
    assert compiler["universal_property4_nonemptiness"] is False
    noncompact_hosts = audit["noncompact_host"]
    assert noncompact_hosts["finite_host_theorem_proved_in_literal_FV13_scope"] is True
    assert noncompact_hosts["universal_raw_return_proved"] is False
    fv13 = noncompact_hosts["literal_full_Picard_full_nef_FV13"]
    assert fv13["ranks_checked"] == list(range(3, 21))
    assert fv13["uniform_intersection_bound"] == 13
    assert fv13["finite_Picard_hosts"] is True
    assert fv13["one_Weyl_orbit_per_oriented_host"] is True
    assert fv13["reflection_maximality_required"] is False
    assert fv13["proper_eigenslices_covered"] is False
    assert fv13["property4_proved"] is False
    torsors = noncompact_hosts["fixed_Jacobian_torsor_calibration"]
    assert torsors["fixed_local_data"] is True
    assert torsors["local_data_bound_index"] is False
    assert torsors["family_claims_finite_volume"] is False
    assert torsors["conditional_q_squared_bound"] == 2877**9
    assert torsors["q_le_390_automatic"] is False
