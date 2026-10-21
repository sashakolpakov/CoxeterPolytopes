from __future__ import annotations

import json

from scripts.certify_polytope_first_coble_h10_hostile import OUTPUT


def _payload() -> dict[str, object]:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def _presentation(chamber_id: str) -> dict[str, object]:
    return next(
        item
        for item in _payload()["presentations"]
        if item["chamber_id"] == chamber_id
    )


def test_frozen_denominator_and_decision() -> None:
    payload = _payload()
    assert payload["all_checks_pass"] is True
    assert payload["summary"] == {
        "Coble_incidence_pass_count": 2,
        "Kondo_catalogue_obstruction_count": 1,
        "Kondo_nonempty_complex_carrier_count": 1,
        "chamber_type_count": 6,
        "characteristic_norm_minus_one_candidate_count": 154,
        "complete_P_package_regeneration_count": 2,
        "decision": (
            "D16-D5-H10-30 is the unique positive in this frozen local "
            "denominator; P10/12 is the only other arithmetic survivor and "
            "is obstructed by its ten-root CM chamber"
        ),
        "integral_presentation_count": 24,
    }


def test_d16_is_reconstructed_from_the_polytope() -> None:
    presentation = _presentation("D16-D5-H10-30")
    assert presentation["characteristic_norm_minus_one_facet_count"] == 1
    candidate = presentation["candidate_audits"][0]
    assert candidate["boundary_facet_one_based"] == 26
    assert len(candidate["root_facet_indices_one_based"]) == 12
    assert len(candidate["exceptional_facet_indices_one_based"]) == 17
    assert candidate["b_perp"]["is_E10"] is True
    assert candidate["cm_cone"]["kind_histogram"] == {
        "finite": 38,
        "ideal": 17,
    }
    assert candidate["cusp_regeneration"]["sets_equal"] is True
    assert candidate["kondo_type_II_root_gram_match"] is True
    assert candidate["surface_classification_gate"] == {
        "first_surface_input": (
            "Kondo's existence/classification of the complex one-boundary "
            "Type-II Coble surface"
        ),
        "kondo_type": "II",
        "status": "nonempty_complex_Coble_carrier",
    }


def test_integral_saturation_U_gate_is_closed() -> None:
    candidate = _presentation("D16-D5-H10-30")["candidate_audits"][0]
    saturation = candidate["integral_saturation_uniqueness"]
    assert saturation["root_span_index_in_E10"] == 2
    assert saturation["root_lattice_determinant"] == -4
    assert saturation["root_lattice_smith_diagonal"] == [1] * 8 + [2, 2]
    assert saturation["discriminant_group"] == "(Z/2)^2"
    assert saturation["nonzero_isotropic_index_two_extensions"] == 2
    assert saturation["weighted_diagram_automorphism_order"] == 48
    assert saturation[
        "diagram_group_transitive_on_even_unimodular_extensions"
    ] is True
    assert saturation["unique_E10_saturation_up_to_diagram_automorphism"] is True


def test_p10_failure_is_only_a_scoped_coble_obstruction() -> None:
    presentation = _presentation("R12-E8A1x2-C2-H10-12")
    candidate = presentation["candidate_audits"][0]
    assert candidate["coble_incidence_pass"] is True
    assert candidate["complete_package_regenerated_from_P"] is True
    assert len(candidate["root_facet_indices_one_based"]) == 10
    assert candidate["cm_cone"]["kind_histogram"] == {
        "finite": 9,
        "ideal": 1,
    }
    assert candidate["surface_classification_gate"]["status"] == (
        "Coble_carrier_obstructed"
    )
    assert "obstruction only" in _payload()["scope"]["failure_means"]
    assert _payload()["logic_layers"]["no_surface_first_wall_list_imported"] is True
    assert _payload()["logic_layers"]["no_project_script_imports"] is True

