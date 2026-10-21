from __future__ import annotations

import json
import shutil

import pytest

from scripts.audit_coble_type_mii_h17_hostile import (
    ALEXEEV_MATRIX,
    ALEXEEV_TEX,
    KONDO_TEX,
    OUTPUT,
    build_audit,
)


pytestmark = pytest.mark.skipif(
    shutil.which("cddexec_gmp") is None
    or not KONDO_TEX.exists()
    or not ALEXEEV_MATRIX.exists()
    or not ALEXEEV_TEX.exists(),
    reason="exact cddlib and pinned primary-source extractions are required",
)


@pytest.fixture(scope="module")
def replay() -> dict:
    return build_audit()


def test_hostile_audit_replays_without_the_type_mii_certificate(replay: dict) -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert replay == archived
    assert replay["independence"] == {
        "imports_Type_MII_certificate": False,
        "imports_Alexeev_reconstruction_archive": False,
        "uses_raw_pinned_primary_sources": True,
        "reconstructs_wall_permutation": True,
        "performs_full_48_wall_double_description": True,
    }


def test_literal_full_cone_excludes_infinite_volume(replay: dict) -> None:
    cone = replay["literal_full_cone"]
    assert (
        cone["extreme_ray_count"],
        cone["finite_vertex_count"],
        cone["ideal_vertex_count"],
        cone["hyperideal_vertex_count"],
    ) == (196_052, 190_808, 5_244, 0)
    assert cone["checks"]["finite_volume"]
    assert cone["checks"]["noncompact"]
    assert replay["checks"]["ideal_count_matches_published_parabolics"]


def test_hidden_single_boundary_minus_one_curves_are_exhausted(replay: dict) -> None:
    cm = replay["CM_cusp_audit"]
    assert (
        cm["finite_vertex_count"],
        cm["ideal_vertex_count"],
        cm["hyperideal_vertex_count"],
    ) == (4_752, 202, 0)
    assert cm["checks"]["entire_linear_cone_is_one_causal_sheet"]
    assert cm["checks"]["no_single_boundary_integral_lift"]
    assert replay["negative_curve_exhaustion"]["exact_cusp_boundary_lift_tests"] == 1_616
    assert replay["negative_curve_exhaustion"]["surviving_lifts"] == 0


def test_raw_alexeev_isometry_was_rediscovered(replay: dict) -> None:
    lineage = replay["raw_Alexeev_lineage"]
    assert lineage["weighted_Gram_isomorphism_count"] == 1_152
    assert lineage["checks"]["integral_unimodular_marking"]
    assert lineage["checks"]["all_48_marked_roots"]
    assert replay["decision"]["abstract_lineage"] == "known Alexeev/Keum--Kondo (18,4,0)"
    assert not replay["decision"]["new_abstract_polytope"]
