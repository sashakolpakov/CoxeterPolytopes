from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

from scripts.audit_coble_type_i_two_boundary_plane_inverse_hostile import (
    OUTPUT,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def test_hostile_certificate_is_an_exact_replay() -> None:
    fresh = build_certificate()
    frozen = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert fresh == frozen
    assert frozen["all_checks_pass"] is True
    assert frozen["verdict"] == "PASS_AFTER_H11_SPECIFIC_RIEMANN_ROCH_REPAIR"


def test_independent_rank_locus_and_whole_bad_component() -> None:
    certificate = build_certificate()["determinantal_and_plane_geometry"]
    matrix = certificate["matrix"]
    assert matrix["shape"] == [23, 21]
    assert matrix["block_sizes"] == [3, 2, 9, 9]
    assert matrix["witness_minor"] == "-64*(a - b)*(a + b)"
    assert matrix["ranks_plus_minus_origin_off"] == [20, 20, 20, 21]
    assert matrix["exact_rank_drop_locus_over_Q"] == "a^2-b^2=0"
    assert matrix["both_components_have_constant_rank_20"] is True
    assert matrix["plus_component_cover_gcd"] == "32"
    assert matrix["minus_component_cover_gcd"] == "16"
    bad = certificate["bad_component"]
    assert bad["coefficientwise_kernel_identity"] is True
    assert bad["every_parameter_is_nonreduced"] is True


def test_quintic_A5_boundary_and_strong_coble_geometry() -> None:
    geometry = build_certificate()["determinantal_and_plane_geometry"]
    good = geometry["good_component"]
    assert good["absolute_audit"]["absolutely_irreducible"] is True
    assert good["absolute_audit"]["projective_singular_scheme_length"] == 10
    assert good["absolute_audit"]["singular_radical_degree"] == 2
    assert [record["multiplicity_sequence"] for record in good["local_resolutions"]] == [
        [2, 2, 2],
        [2, 2, 2],
    ]
    assert [record["terminal_discriminant"] for record in good["local_resolutions"]] == [
        "32",
        "-32",
    ]
    assert good["delta_sum"] == good["arithmetic_genus"] == 6
    assert geometry["contacts"]["boundary_line_restriction"] == "16*x**3*(x - 1)**2"
    boundary = geometry["boundary"]
    assert boundary["line_square"] == boundary["quintic_square"] == -4
    assert boundary["line_quintic"] == 0
    assert boundary["sum_equals_minus_2K"] is True
    assert boundary["strong_Coble"] is True


def test_missing_lifts_are_effective_irreducible_before_exact_return() -> None:
    transfer = build_certificate()["transfer_and_return"]
    gap = transfer["generic_D_CobVin_handoff_audit"]
    assert gap["verdict"] == "FAIL_FOR_ROOT_SEED_ONLY_INVOCATION"
    assert "not challenged" in gap["full_incidence_compiler_theorem"]
    branches = transfer["four_plane_branches"]
    assert len(branches) == 4
    expected_degrees = [2, 2, 2, 2, 3, 3, 4, 4, 5, 5]
    for branch in branches:
        assert branch["basis_determinant_absolute"] == 1
        assert branch["missing_lift_count"] == 10
        assert branch["missing_degree_multiset"] == expected_degrees
        assert branch["all_ten_effective_by_P2_RR"] is True
        assert all(
            lift["square"] == -1
            and lift["K_dot"] == -1
            and lift["chi_O_E"] == 1
            and lift["degree_K_minus_E"] < 0
            and lift["effective_by_Riemann_Roch"]
            for lift in branch["missing_lifts"]
        )
    repair = transfer["H11_specific_RR_repair"]
    assert repair["Riemann_Roch_effectivity"] is True
    assert repair["all_distinct_raw_wall_pairings_nonnegative"] is True
    assert repair["all_26_raw_walls_effective_irreducible"] is True
    assert repair["no_hidden_negative_walls"] is True
    assert repair["literal_same_marked_P"] is True


def test_bounded_cli_frozen_check() -> None:
    run = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/audit_coble_type_i_two_boundary_plane_inverse_hostile.py"),
            "--check-output",
            str(OUTPUT),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr
