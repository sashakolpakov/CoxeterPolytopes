from __future__ import annotations

import json
from pathlib import Path
import shutil

import pytest

from scripts.certify_coble_type_mii_full_ns_h17 import (
    OUTPUT,
    build_certificate,
)


pytestmark = pytest.mark.skipif(
    shutil.which("cddexec_gmp") is None,
    reason="cddexec_gmp is required for the exact CM cusp enumeration",
)


def test_type_mii_full_ns_h17_archive_replays_exactly() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert build_certificate() == archived
    assert archived["all_checks_pass"]


def test_type_mii_has_complete_48_wall_full_nef_chamber() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert archived["AG_native_stratum"]["hyperbolic_dimension"] == 17
    assert archived["curve_exhaustion"]["complete_negative_curve_count"] == 48
    assert archived["curve_exhaustion"]["integral_single_boundary_lifts"] == 0
    assert archived["curve_exhaustion"]["total_single_boundary_lift_tests"] == 1_616
    assert archived["facets"]["norm_histogram"] == {
        "-4": 8,
        "-2": 24,
        "-1": 16,
    }
    assert archived["facets"]["all_facets_irredundant"]
    assert archived["finite_volume"]["decision"] == "finite_volume_noncompact"


def test_type_mii_cm_cusps_and_known_lineage_are_frozen() -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    cone = archived["Coble_Mukai_chamber"]["exact_double_description"]
    assert (
        cone["finite_vertex_count"],
        cone["ideal_vertex_count"],
        cone["hyperideal_vertex_count"],
    ) == (4_752, 202, 0)
    assert sorted(
        row["size"]
        for row in archived["Coble_Mukai_chamber"][
            "cusp_orbits_under_C2_semidirect_S4xS4"
        ]
    ) == [16, 18, 72, 96]
    lineage = archived["lineage"]
    assert lineage["exact_identification"] == "Alexeev/Keum--Kondo (18,4,0), H17/48"
    assert lineage["checks"]["all_2304_weighted_Gram_entries_equal"]
    assert lineage["checks"]["integral_marking_unimodular"]
    assert not lineage["priority"]["new_abstract_polytope"]
    assert not archived["inverse_program_status"]["inverse_theorem_credit"]
