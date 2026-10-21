"""Independent hostile tests for fixed-Jacobian K3 cusp-index blindness."""

from __future__ import annotations

import ast
import json
from pathlib import Path
import re
import subprocess
import sys

import pytest

from scripts.audit_noncompact_k3_fixed_jacobian_index_blindness_independent_hostile import (
    PINS,
    ROOT,
    archive_certificate,
    bounded_spanning_interface_certificate,
    build_audit,
    family_lattice_certificate,
    fixed_jacobian_certificate,
    marked_binary_certificate,
    sha256,
)


AUDIT = ROOT / "results/noncompact_k3_fixed_jacobian_index_blindness_independent_hostile_audit.md"
SCRIPT = ROOT / "scripts/audit_noncompact_k3_fixed_jacobian_index_blindness_independent_hostile.py"


def test_all_owner_and_global_input_bytes_are_frozen() -> None:
    for relative, expected in PINS.values():
        assert sha256(ROOT / relative) == expected


def test_owner_archive_self_hash_scope_and_replay() -> None:
    record = archive_certificate()
    assert record["owner_archive_self_hash_valid"] is True
    assert record["negative_scope_flags_valid"] is True
    assert record["sample_d_values"] == [2, 3, 4, 5, 7, 11, 19, 31]
    result = subprocess.run(
        [
            sys.executable,
            str(ROOT / PINS["owner_certifier"][0]),
            "--check",
            str(ROOT / PINS["owner_archive"][0]),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "PASS archive replay" in result.stdout


def test_fixed_jacobian_has_two_II_star_and_four_I1_fibres() -> None:
    record = fixed_jacobian_certificate()
    assert record["orders_at_zero"] == [4, 5, 10]
    assert record["orders_at_infinity"] == [4, 5, 10]
    assert record["residual_discriminant_degree"] == 4
    assert record["kodaira_fibres"] == ["II*", "II*", "I1", "I1", "I1", "I1"]
    assert record["euler_sum"] == 24
    assert record["trivial_lattice"] == "U+E8(-1)^2"
    assert record["very_general_picard_rank"] == 18
    assert record["very_general_transcendental_lattice"] == "U^2"
    assert record["shioda_tate_MW_rank"] == 0


@pytest.mark.parametrize("d", range(2, 51))
def test_exact_brauer_kernel_primitive_complement_and_marked_packet(d: int) -> None:
    record = family_lattice_certificate(d)
    assert record["brauer_character_order"] == d
    assert record["brauer_d_torsion_rank"] == 4
    assert record["transcendental_lattice"] == "U(d)+U"
    assert abs(record["transcendental_primitive_unit_minor"]) == 1
    assert record["host_lattice"] == "U(d)+E8(-1)^2"
    assert record["host_primitive_unit_minor"] == 1
    assert record["host_is_exact_orthogonal_complement"] is True
    assert record["absolute_host_discriminant"] == d**2
    assert record["discriminant_group_invariant_factors"] == [d, d]
    assert record["discriminant_length"] == 2
    assert record["fibre_divisibility"] == record["multisection_index"] == d
    assert record["affine_component_count"] == 18
    assert record["affine_vector_span_rank"] == 17
    assert record["affine_gram_rank"] == record["finite_ADE_rank"] == 16
    assert record["common_orthogonal_dimension"] == 1
    assert record["formal_transverse_root_square"] == -2
    assert record["formal_transverse_root_degree"] == d
    assert record["formal_transverse_root_nonnegative_on_packet"] is True


@pytest.mark.parametrize("d", range(2, 51))
def test_full_discriminant_form_forces_the_actual_marked_Ud(d: int) -> None:
    record = marked_binary_certificate(d)
    assert record["smith_group_forces_d_divides_2a"] is True
    assert record["basis_shift_changes_a_by_multiple_of_d"] is True
    assert record["standard_Ud_has_isotropic_generating_pair"] is True
    assert record["full_discriminant_form_forces_a_zero_mod_d"] is True
    if d % 2:
        assert record["candidate_a_mod_d"] == [0]
        assert record["exceptional_even_residue_possible_before_q_check"] is False
    else:
        assert record["candidate_a_mod_d"] == [0, d // 2]
        assert record["exceptional_every_isotropic_first_coordinate_even"] is True
        assert record["exceptional_has_isotropic_generating_pair"] is False


def test_corrected_B13_interface_is_hadamard_not_automatic_390() -> None:
    record = bounded_spanning_interface_certificate()
    assert record["rank"] == 18
    assert record["physical_pairing_bound"] == 13
    assert record["euclidean_row_norm_squared_bound"] == 2877
    assert record["hadamard_formula"] == "abs(det G)<=2877^9"
    assert record["automatic_finite_volume_consequence"] == "d^2<=2877^9"
    assert record["selected_II_star_mark_sum"] == 30
    assert record["cusp_adapted_conditional_consequence"] == "d<=30*13=390"
    assert record["arbitrary_spanning_packet_contains_selected_II_star_packet"] is False
    assert record["d_391_not_excluded_by_hadamard_bound"] is True
    assert record["d_le_390_automatic_from_B13"] is False


def test_full_audit_scope_forbids_a_global_counterfamily_claim() -> None:
    audit = build_audit()
    assert audit["verdict"] == "PASS_LOCAL_INDEX_BLINDNESS_EXACT_SCOPE"
    scope = audit["scope"]
    for positive in (
        "same_fixed_relative_jacobian",
        "projective_torsors_of_every_order",
        "same_two_II_star_and_four_I1_fibres",
        "same_full_local_parabolic_packet",
        "unbounded_fibre_divisibility_without_global_volume",
        "conditional_d_squared_le_2877_to_the_9th_under_finite_volume",
    ):
        assert scope[positive] is True
    for negative in (
        "finite_volume_claimed_for_family",
        "reflection_maximality_claimed_for_family",
        "formal_transverse_root_claimed_physical",
        "automatic_d_le_390_claimed",
        "Property_4_proved",
    ):
        assert scope[negative] is False
    assert audit["source_theorems"]["machine_proved_here"] is False


def test_independent_script_imports_no_owner_implementation() -> None:
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    imported_modules: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported_modules.append(node.module)
    assert not any(
        module.endswith("certify_noncompact_k3_fixed_jacobian_index_blindness")
        for module in imported_modules
    )


def test_independent_command_line_rechecks_pins() -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--check-pins"],
        check=True,
        capture_output=True,
        text=True,
    )
    payload = json.loads(result.stdout)
    assert payload["verdict"] == "PASS_LOCAL_INDEX_BLINDNESS_EXACT_SCOPE"
    assert set(payload["pins"]) == set(PINS)


def test_report_contains_every_hostile_gate_and_primary_source_boundary() -> None:
    text = " ".join(AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "PASS, exactly in the asserted cusp-local",
        "projective genus-one K3 torsors of every index",
        "not a finite-volume family",
        r"d^2\leq2877^9",
        r"not \(d\leq390\)",
        "surjectivity in Theorem 4.13",
        r"\operatorname {ind}(X_d)=\operatorname {ord}(\alpha_d)=d",
        r"T(X_d)\simeq U(d)\oplus U",
        "one primitive-embedding orbit",
        "wrong discriminant form",
        r"J^k(X)_t\simeq X_t",
        "strongest possible local parabolic rank",
        "not call it an irreducible curve",
        r"4+17\cdot13^2=2877",
        "arbitrary spanning packet",
        "Property 4",
    ):
        assert phrase in text


def test_transport_and_control_bytes_are_clean() -> None:
    payload = AUDIT.read_bytes() + SCRIPT.read_bytes() + Path(__file__).read_bytes()
    forbidden = set(range(0, 9)) | {11, 12} | set(range(14, 32)) | {127}
    assert not any(byte in forbidden for byte in payload)
    text = payload.decode("utf-8")
    assert "," + "quad" not in text
    assert "," + "qquad" not in text
    assert "|" + "le" not in text
    assert "TO" + "DO" not in text
    assert "FIX" + "ME" not in text


def test_no_inline_tex_delimiters_were_stripped() -> None:
    source = AUDIT.read_text(encoding="utf-8")
    outside_display: list[str] = []
    in_display = False
    for line in source.splitlines():
        if line.strip() == r"\[":
            in_display = True
            continue
        if line.strip() == r"\]":
            in_display = False
            continue
        if not in_display:
            outside_display.append(line)
    prose = "\n".join(outside_display)
    prose = re.sub(r"\\\(.*?\\\)", "", prose)
    prose = re.sub(r"\x60[^\x60]*\x60", "", prose)
    for command in (
        r"\operatorname",
        r"\mathbf",
        r"\simeq",
        r"\perp",
        r"\leq",
        r"\geq",
        r"\alpha",
        r"\theta",
    ):
        assert command not in prose


def test_invalid_index_is_rejected() -> None:
    with pytest.raises(ValueError):
        family_lattice_certificate(1)
    with pytest.raises(ValueError):
        marked_binary_certificate(1)
