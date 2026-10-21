from __future__ import annotations

import json

from coxeter_polytopes.combinatorics.h9_14_record12_residual import (
    EXPECTED_COLOURED_BRANCH_ORBIT_COUNT,
    EXPECTED_FULL_BLOCKER_PASS_COUNT,
    EXPECTED_FULL_ORBIT_COUNT,
    EXPECTED_SCHEDULE_ORBIT_COUNTS,
    add_y_to_branch_blocker,
    minimal_transversals,
    verify_archive,
)
from scripts.certify_h9_14_record12_residual import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT as CERTIFICATE,
    RESIDUAL_ARCHIVE,
    verify as verify_certificate,
)


def test_record12_corrected_exact_371_archive() -> None:
    document = json.loads(RESIDUAL_ARCHIVE.read_text(encoding="utf-8"))
    summary = verify_archive(document)
    assert summary["coloured_branch_orbits"] == EXPECTED_COLOURED_BRANCH_ORBIT_COUNT
    assert summary["relative_embedding_multiplicity_sum"] == (
        EXPECTED_FULL_BLOCKER_PASS_COUNT
    )
    assert summary["residual_orbits"] == EXPECTED_FULL_ORBIT_COUNT
    assert summary["schedule_orbit_counts"] == EXPECTED_SCHEDULE_ORBIT_COUNTS
    assert len(document["residual_orbits"]) == 371
    assert len(
        {row["selected_complement_sha256"] for row in document["residual_orbits"]}
    ) == 371
    joined_nonclaims = " ".join(document["nonclaims"])
    assert "not certified PL spheres" in joined_nonclaims
    assert "ordinary convex-polytopal realizability is not claimed" in joined_nonclaims
    assert "No compact or finite-volume hyperbolic Coxeter polytope" in joined_nonclaims


def test_record12_empty_clutter_and_y_blocker_edge_cases() -> None:
    # The blocker of the empty clutter is {empty}; an empty list was the bug in
    # an exploratory branch stream and must never silently return here again.
    assert minimal_transversals(()) == (0,)
    assert add_y_to_branch_blocker((0,), 0) == (0,)
    y = (1 << 2) | (1 << 7)
    assert add_y_to_branch_blocker((0,), y) == (1 << 10, 1 << 11, y)


def test_record12_corrected_371_frozen_certificate() -> None:
    manifest = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    summary = verify_certificate(manifest)
    assert summary["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert summary["theorem_status"] == "certified_exact_371_orbit_abstract_residual"
    assert summary["complete_blocker_passes"] == 42_000
    assert summary["residual_orbits"] == 371
    assert manifest["cut_normalization_audit"][
        "both_q_and_complement_q_enumerated"
    ] is True
    assert manifest["correction_history"]["new_relative_transpose_orbits"] == 27
