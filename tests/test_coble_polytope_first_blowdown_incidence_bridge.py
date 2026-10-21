from __future__ import annotations

import json
import shutil

import pytest

from scripts.certify_coble_polytope_first_blowdown_incidence_bridge import (
    OUTPUT,
    build_certificate,
)


pytestmark = pytest.mark.skipif(
    shutil.which("cddexec_gmp") is None,
    reason="cddexec_gmp is required for the exact skeleton-cone audits",
)


def _archive() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_frozen_blowdown_incidence_certificate_replays() -> None:
    assert build_certificate() == _archive()
    assert _archive()["all_checks_pass"]


def test_raw_mii_compiler_is_exhaustive_and_forces_characteristic_three() -> None:
    mii = _archive()["type_MII_raw_compiler"]
    search = mii["boundary_frame_search"]
    assert (
        search["all_eight_subsets"],
        search["pairwise_orthogonal"],
        search["q_over_two_integral_on_P_generated_module"],
    ) == (735_471, 21_211, 1)
    assert search[
        "odd_pair_coverage_is_equivalent_to_integrality_on_wall_generators"
    ]
    assert mii["reconstructed_Picard"]["identification"] == "odd unimodular I_(1,17)"
    assert mii["P_compiled_blowdown"]["endpoint"] == "P1xP1"
    assert mii["P_compiled_blowdown"]["grid_edge_count"] == 16
    assert mii["P_compiled_blowdown"]["all_24_permutations_occur"]
    assert mii["incidence_nonemptiness_iff"]["allowed_characteristics_not_two"] == [3]
    assert mii["incidence_nonemptiness_iff"]["iff_proved_without_Kondo_template"]
    assert all(mii["direct_strong_Coble"]["checks"].values())


def test_type_vi_has_five_plane_blowdowns_and_both_skeletons_fail_volume() -> None:
    payload = _archive()
    vi = payload["type_VI_plane_compiler"]
    blowdowns = vi["boundary_preserving_P2_blowdowns"]
    assert blowdowns["count"] == 5
    assert {
        row["conic_boundary_index_zero_based"] for row in blowdowns["blowdowns"]
    } == set(range(5))
    assert all(sorted(row["boundary_degrees"]) == [1, 1, 1, 1, 2] for row in blowdowns["blowdowns"])
    assert vi["residual_root_incidence"]["allowed_characteristics_not_two"] == [3]
    assert vi["residual_root_incidence"]["first_two_already_force"] == "3=0 and lambda=2=-1"
    assert vi["boundary_skeleton_cone"]["hyperideal"] == 161
    assert payload["type_MII_raw_compiler"]["boundary_skeleton_cone"]["hyperideal"] == 1_464
    assert vi["full_cone"]["hyperideal"] == 0
    assert payload["logic_classification"]["level"].endswith("not L4")
