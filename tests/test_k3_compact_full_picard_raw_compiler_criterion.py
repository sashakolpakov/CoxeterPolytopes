"""Regression tests for the compact full-Picard K3 iff compiler."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_k3_compact_full_picard_raw_compiler_criterion import (
    ARCHIVE,
    REPORT,
    build_certificate,
    digest,
    encoded_certificate,
    golden_energy_shell,
    hidden_root_ellipsoid_toy,
    nonrational_rank_window,
    quadratic_sign,
)


CERTIFIER = (
    ROOT / "scripts/certify_k3_compact_full_picard_raw_compiler_criterion.py"
)


def test_quadratic_sign_is_exact_on_both_sides() -> None:
    assert quadratic_sign(3, -1) == 1
    assert quadratic_sign(2, -1) == -1
    assert quadratic_sign(-3, 1) == -1
    assert quadratic_sign(-2, 1) == 1
    assert quadratic_sign(0, 4) == 1
    assert quadratic_sign(0, -4) == -1


def test_golden_trace_minus_two_energy_shell_is_finite() -> None:
    shell = golden_energy_shell()
    assert shell["required_trace"] == -2
    assert shell["integer_u_lower"] == -1
    assert shell["integer_u_upper"] == 3
    assert shell["coordinates"] == [
        [-1, -2],
        [0, -2],
        [1, -2],
        [2, -2],
        [3, -2],
    ]
    assert shell["cardinality"] == 5
    assert shell["lower_excluded_u"] == -2
    assert shell["upper_excluded_u"] == 4


def test_hidden_root_ellipsoid_is_positive_and_exhaustive_in_toy() -> None:
    toy = hidden_root_ellipsoid_toy()
    gram_qp = sp.Matrix(toy["positive_definite_Qp_gram"])
    assert gram_qp.is_positive_definite
    assert toy["ellipsoid_bound_2cosh_2R"] == 14
    assert toy["signed_root_count"] == 12
    assert toy["unsigned_mirror_count"] == 6
    assert all(
        2 * (a * a + b * b + c * c) <= 14
        for a, b, c in toy["signed_roots"]
    )
    assert all(abs(a) <= 1 for a, _, _ in toy["signed_roots"])
    assert toy["first_excluded_Qp"] == 18
    assert toy["first_excluded_signed_root_count"] == 16


def test_nonrational_k3_rank_window() -> None:
    rows = nonrational_rank_window()
    by_rank = {row["module_rank"]: row for row in rows}
    assert min(by_rank) == 3
    assert max(by_rank) == 10
    assert [by_rank[m]["field_degree_ceiling"] for m in range(3, 11)] == [
        6,
        5,
        4,
        3,
        2,
        2,
        2,
        2,
    ]
    assert all(row["nonrational_possible"] for row in rows)


def test_compiler_gate_order_and_empty_scope() -> None:
    payload = build_certificate()
    assert payload["physical_normalization_gate"]["place_energies_sum"] == 2
    assert payload["module_gate"]["finite_interval"] == (
        "F0 subset M subset F0#"
    )
    assert payload["module_gate"]["multiplier_order"] == (
        "O_M={a in k:aM subset M}"
    )
    assert payload["module_gate"]["projectivity_over_O0_required"] is False
    assert payload["embedding_gate"]["necessary_and_sufficient"] is True
    assert payload["hidden_root_gate"]["Qp_positive_definite"] is True
    assert payload["hidden_root_gate"]["exact_lattice_enumeration_finite"] is True
    assert payload["construction"]["selected_p_lies_in_E_and_is_ample"] is True
    assert (
        payload["construction"][
            "return_is_closure_of_projectivized_Amp_intersection"
        ]
        is True
    )
    verdict = payload["verdict"]
    assert verdict["necessary_and_sufficient_for_pinned_realization"] is True
    assert verdict["compiler_may_return_EMPTY"] is True
    assert verdict["universal_nonemptiness_for_every_raw_maximal_P"] is False
    assert verdict["Property_4_essential_surjectivity_proved"] is False
    assert verdict["first_nonautomatic_gate"] == "PHYS"
    assert verdict["later_exact_gates"] == ["EMB", "ROOT"]


def test_maximality_species_and_cusped_comparison_are_separated() -> None:
    payload = build_certificate()
    maximality = payload["maximality_scope"]
    assert maximality["ordinary_root_Q_minus_dimension"] == 1
    assert maximality["native_root_Q_minus_dimension"] == "[k:Q]"
    assert (
        maximality["compact_hidden_root_test_removed_by_raw_maximality_alone"]
        is False
    )
    assert (
        maximality[
            "same_host_native_closure_plus_relative_maximality_is_shortcut"
        ]
        is True
    )
    cusped = payload["cusped_comparison"]
    assert cusped["exact_finite_volume_full_trace_forces_k_Q"] is True
    assert cusped["bounded_closed_star_physical_basis_constant"] == 13
    assert cusped["finite_rational_Picard_host_census"] is True
    assert cusped["separate_from_compact_nonrational_compiler"] is True


def test_report_pins_the_iff_and_nonclaim() -> None:
    text = REPORT.read_text(encoding="utf-8")
    for needle in (
        "finite decision theorem for a pinned exact",
        "\\mathscr E_{-2}(\\mathcal O_0)",
        "\\frac{u_i/h_i}{u_1/h_1}\\in k^{\\times2}",
        "F_0\\subseteq M\\subseteq F_0^\\#",
        "q_T\\simeq-q_N",
        "Q_p(\\delta)\\le2\\cosh(2R)",
        "zero distinguished projection",
        "p\\in E\\cap\\operatorname {Amp}(S)",
        "\\overline{\\mathbf P\\bigl(\\operatorname {Amp}(S)\\cap E\\bigr)}",
        "is part of the return contract",
        "finite iff/EMPTY compiler",
        "The compiler can return `EMPTY`",
        "not Property 4 by itself",
        "ordinary K3-root maximality",
        "essential-surjectivity/nonemptiness",
    ):
        assert needle in text


def test_archive_is_exact_replay() -> None:
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived["payload"] == build_certificate()
    assert archived["payload_sha256"] == digest(archived["payload"])
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()


def test_cli_check() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check-output"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
