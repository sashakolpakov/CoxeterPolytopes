from __future__ import annotations

import json

from coxeter_polytopes.combinatorics.h9_14_record0 import (
    TARGET_F_VECTOR,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    audit_record0_polytopality,
    target_facet_complements,
)
from scripts.certify_h9_14_record0 import OUTPUT, verify
from scripts.independent_h9_14_record0_audit import replay


def test_record0_exact_rational_construction() -> None:
    audit = audit_record0_polytopality()
    assert audit.verified
    assert audit.stellar_facet_counts == (20, 28, 34)
    assert audit.base_face_vector == (10, 41, 84, 85, 34)
    assert audit.blocker_order_histogram == ((5, 168),)
    assert len(audit.expanded_minimal_nonfaces) == (
        len(TARGET_MISSING_EDGES)
        + len(TARGET_MISSING_FOURS)
        + len(TARGET_MISSING_FIVES)
    )
    assert len(target_facet_complements()) == 168


def test_record0_frozen_manifest() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify(manifest)
    assert summary["facet_count"] == 168
    assert tuple(summary["face_vector"]) == TARGET_F_VECTOR
    assert manifest["theorem_status"] == "certified_explicit_ordinary_polytopal"
    assert manifest["target"]["missing_four_schedule"] == 6
    joined_nonclaims = " ".join(manifest["nonclaims"])
    assert "No new hyperbolic Coxeter polytope" in joined_nonclaims
    assert "one completion" in joined_nonclaims


def test_record0_hostile_replay() -> None:
    summary = replay()
    assert summary["ordinary_rational_polytopal"] is True
    assert summary["blocker_order_histogram"] == {"5": 168}
    assert summary["ridge_degree_histogram"] == {"2": 756}
    assert tuple(summary["target_face_vector"]) == TARGET_F_VECTOR
    assert summary["all_links_mod2_homology_spheres"] is True

