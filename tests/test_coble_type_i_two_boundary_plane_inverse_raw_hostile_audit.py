from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

from scripts.independent_coble_type_i_two_boundary_plane_inverse_raw_hostile import (
    OUTPUT,
    build_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def test_raw_hostile_archive_replays_and_states_both_verdicts() -> None:
    fresh = json.loads(json.dumps(build_certificate()))
    frozen = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert fresh == frozen
    assert frozen["all_checks_pass"] is True
    assert frozen["verdict"]["original_primary_certificate"].startswith("FAIL_AS_WRITTEN")
    assert frozen["verdict"]["repaired_theorem"].startswith("PASS")


def test_raw_endpoint_and_all_twenty_six_class_transport() -> None:
    certificate = json.loads(OUTPUT.read_text(encoding="utf-8"))
    raw = certificate["raw_endpoint_replay"]
    assert [
        raw["all_sequential_P2_endpoint_count"],
        raw["endpoint_degree_signature_count"],
        raw["selected_line_quintic_endpoint_count"],
    ] == [52, 9, 4]
    branches = certificate["branches"]
    assert len(branches) == 4
    assert len({branch["canonical_wall_set_sha256"] for branch in branches}) == 1
    assert all(
        abs(branch["plane_basis_determinant"]) == 1
        and branch["all_26_classes_reconstruct_raw_vectors"]
        for branch in branches
    )


def test_explicit_missing_curves_and_exact_cones() -> None:
    certificate = json.loads(OUTPUT.read_text(encoding="utf-8"))
    repair = certificate["explicit_missing_curve_repair"]
    assert repair["degree_multiset"] == [2, 2, 2, 2, 3, 3, 4, 4, 5, 5]
    assert len(repair["curves"]) == 10
    assert all(
        curve["linear_system_dimension_projective"] == 0
        and curve["absolute_factor_audit"]["absolutely_irreducible"]
        and curve["arithmetic_genus_after_blowup"] == 0
        and curve["strict_transform_smooth_geometrically_integral_rational"]
        for curve in repair["curves"]
    )
    cm = certificate["CM_root_and_cusp_replay"]
    assert cm["CM_finite_ideal_hyperideal"] == [14, 9, 0]
    assert cm["split_and_together_cusps"] == [6, 3]
    assert cm["all_and_only_integral_cusp_lifts_equal_twelve_raw_walls"] is True
    full = certificate["full_cone_replay"]
    assert full["extreme_ray_count"] == 550
    assert full["finite_ideal_hyperideal"] == [508, 42, 0]
    assert full["all_26_facets_irredundant"] is True


def test_raw_hostile_cli_under_external_guard_contract() -> None:
    run = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/independent_coble_type_i_two_boundary_plane_inverse_raw_hostile.py"),
            "--check-output",
            str(OUTPUT),
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr
