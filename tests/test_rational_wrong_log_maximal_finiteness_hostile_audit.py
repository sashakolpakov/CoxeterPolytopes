import json
from pathlib import Path

import pytest

from scripts.audit_rational_wrong_log_maximal_finiteness_hostile import (
    ARCHIVE,
    OWNER_PINS,
    adapted_wall_audit,
    ag_rank_elimination,
    binary_subspace_counts,
    build_audit,
    chamber_orbit_audit,
    field_audit,
    is_totally_isotropic,
    literal_binary_subspaces,
    rendered,
    sha256,
)


ROOT = Path(__file__).resolve().parents[1]
OWNER_REPORT = (
    ROOT / "results/rational_wrong_log_reflection_maximal_finiteness.md"
)
AUDIT_REPORT = (
    ROOT / "results/rational_wrong_log_maximal_finiteness_hostile_audit.md"
)


def test_owner_artifact_pins_and_archive_replay():
    payload = build_audit()
    assert payload["owner_pins_match"]
    for relative, expected in OWNER_PINS.items():
        assert sha256(ROOT / relative) == expected
    assert ARCHIVE.read_text(encoding="utf-8") == rendered(payload)
    assert json.loads(ARCHIVE.read_text(encoding="utf-8"))["all_checks_pass"]


def test_symbolic_rank_elimination_keeps_only_the_exact_26_endpoint():
    result = ag_rank_elimination(5000)
    assert result["symbolic_elimination"] == (
        "7(rho-10)<=4rho+8, hence 3rho<=78"
    )
    assert result["maximum_feasible_rank"] == 26
    assert result["rank_26_boundary_counts"] == [16]
    assert result["rank_27_boundary_counts"] == []
    assert result["no_feasible_rank_27_through_horizon"]
    assert not result["endpoint_is_eliminated_as_nonexistent"]


def test_all_adapted_wall_reflections_and_missing_W2_control():
    walls = adapted_wall_audit()
    assert walls["boundary_root_norm"] == "-2"
    assert walls["minus_one_root_norm"] == "-2"
    assert walls["minus_two_root_norm"] == "-4"
    assert walls["boundary_reflection_coefficients_integral"]
    assert walls["minus_one_reflection_coefficients_integral"]
    assert walls["minus_two_lattice_coefficients_integral"]
    assert walls["minus_two_half_boundary_coefficient"] == "0"
    assert walls["missing_W2_control_coefficient"] == "1/2"
    assert walls["missing_W2_control_is_nonintegral"]
    assert not walls["W3_used_for_host_integrality"]
    assert walls["W3_used_for_AG_nefness"]


def test_host_count_recurrence_against_literal_subspace_enumeration():
    recurrence = binary_subspace_counts(26)
    assert recurrence[1:6] == [2, 5, 16, 67, 374]
    literal = [
        len(literal_binary_subspaces(rank)) for rank in range(1, 6)
    ]
    isotropic = [
        sum(
            is_totally_isotropic(space)
            for space in literal_binary_subspaces(rank)
        )
        for rank in range(1, 6)
    ]
    assert literal == recurrence[1:6]
    assert isotropic == [1, 2, 4, 11, 31]
    assert 2 * recurrence[2] + sum(recurrence[3:27]) == (
        5515561289142524504797715176506813442521589215285286
    )


def test_hostile_helper_input_rejection():
    for bad in (-1, True):
        with pytest.raises(ValueError):
            literal_binary_subspaces(bad)
    with pytest.raises(ValueError):
        literal_binary_subspaces(7)
    with pytest.raises(ValueError):
        ag_rank_elimination(26)


def test_all_three_maximality_modes_and_negative_controls():
    modes = build_audit()["maximality_modes"]
    assert modes["absolute_discrete_reflection_maximality_saturates"]
    assert modes["host_relative_reflection_maximality_saturates"]
    assert modes["arithmetic_reflection_maximality_saturates"]
    assert "finite index" in modes["arithmetic_reason"]
    assert not modes["lattice_maximality_alone_saturates"]
    assert not modes["known_effective_walls_alone_saturate"]


def test_complete_chamber_orbit_and_rational_field_controls():
    orbit = chamber_orbit_audit()
    assert orbit == {
        "A2_complete_chamber_count": 6,
        "A2_complete_orbit_size": 6,
        "A2_incomplete_orbit_size": 2,
        "complete_arrangement_transitive": True,
        "selected_subarrangement_need_not_be_transitive": True,
    }
    field = field_audit()
    assert field["all_entries_and_cycles_rational"]
    assert field["minimal_Vinberg_field"] == "Q"
    assert not field["surface_field_of_definition_claimed"]


def test_scope_and_repaired_reports_are_clean():
    payload = build_audit()
    scope = payload["scope"]
    assert scope["carrier_relative"]
    assert scope["literal_complete_full_nef_required"]
    assert scope["finite_sided_finite_volume_required"]
    assert not scope["raw_carrier_surjectivity_claimed"]
    assert not scope["hostile_incidence_bound_claimed"]
    assert not scope["terminating_reflectivity_census_claimed"]
    assert not scope["surface_classification_avoided"]
    assert scope["reflective_lattice_classification_avoided"]

    for path in (OWNER_REPORT, AUDIT_REPORT):
        raw = path.read_bytes()
        assert b"\r" not in raw
        assert not any(byte < 32 and byte not in (9, 10) for byte in raw)

    owner = OWNER_REPORT.read_text(encoding="utf-8")
    for stale in (
        "preserves (M)",
        "Let (P=",
        "and (W(P))",
        "arrangement of (R(M))",
        "Regard (L(2))",
        "If (H\\leq",
        "(k(P)=\\mathbf Q)",
    ):
        assert stale not in owner
    for needle in (
        "reflective-classification-free",
        "three senses of Section 3",
        r"\boxed{W(P)=R(M)}",
        r"\boxed{n\leq25}",
        "raw reflection-maximal arithmetic chamber is not yet proved",
    ):
        assert needle in owner

    audit = AUDIT_REPORT.read_text(encoding="utf-8")
    for needle in (
        r"(\rho,b)=(26,16)",
        "it does not prove that rank \\(26\\) itself is empty",
        "all three physical wall types",
        "Saturation under each maximality mode",
        "Exactly one complete chamber orbit",
        r"k(P)=\mathbf Q",
        "carrier-relative",
        "essential surjectivity",
    ):
        assert needle in audit
