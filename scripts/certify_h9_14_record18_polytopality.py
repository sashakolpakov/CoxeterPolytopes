#!/usr/bin/env python3
"""Freeze the exact rational ordinary-polytopality proof for record 18."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record18_polytopality import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_RECORD_INDEX,
    audit_record18_polytopality,
)


ROOT = Path(__file__).resolve().parents[1]
PROOF_SOURCE = (
    ROOT
    / "src/coxeter_polytopes/combinatorics/h9_14_record18_polytopality.py"
)
OUTPUT = ROOT / "data/enumerations/h9_14_record18_polytopality.json"
EXPECTED_PROOF_SOURCE_SHA256 = (
    "9de308c9a4dc54b84daaa2a8e4c38ca062da4dee28f352a41a2beceabaff165e"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "6dae56855a8380ee8cffd8b3e2572df90164f98e483e717f9323de8289d4cb90"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _certificate_digest(manifest: Mapping[str, object]) -> str:
    payload = dict(manifest)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def build() -> dict[str, object]:
    source_hash = _file_sha256(PROOF_SOURCE)
    if source_hash != EXPECTED_PROOF_SOURCE_SHA256:
        raise ValueError("record-18 polytopality proof source is not frozen")
    audit = audit_record18_polytopality()
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record18_exact_rational_polytopality",
        "claim": (
            "The selected H9/14 record-18 nerve is the boundary of a rational "
            "simplicial convex 9-polytope; dually it is the nerve of a rational "
            "simple convex 9-polytope with fourteen facets."
        ),
        "theorem_status": "certified_explicit_ordinary_polytopal",
        "target": {
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
            "missing_edges": [list(face) for face in TARGET_MISSING_EDGES],
            "missing_fives": [list(face) for face in TARGET_MISSING_FIVES],
        },
        "proof_source": {
            "path": (
                "src/coxeter_polytopes/combinatorics/"
                "h9_14_record18_polytopality.py"
            ),
            "sha256": source_hash,
        },
        "proof_backend": (
            "standard-library exact rational arithmetic: Gale evenness, two "
            "abstract stellar subdivisions, exhaustive supporting-hyperplane "
            "determinants, radial-beyond coordinate identities, and inverse-"
            "multiwedge minimal-nonface equality"
        ),
        "audit": audit.to_record(),
        "construction": (
            "Take C(10,7) on moment-curve labels "
            "[1,0,10,5,9,2,3,4,7,11]. Stellar-subdivide {2,4,5,11} by "
            "vertex 12 and then {0,5,7} by vertex 13. Replace carriers 5 "
            "and 7 by clone pairs {5,6} and {7,8}. The two simplicial "
            "wedges raise dimension seven to nine and reproduce the target "
            "clutter exactly."
        ),
        "proof_lemmas": [
            (
                "A radial-beyond placement y+epsilon(y-p), with y in the "
                "relative interior of the stellar face and p in the polytope "
                "interior, realizes its boundary stellar subdivision for "
                "sufficiently small positive rational epsilon; epsilon=1/32 "
                "is verified exactly at both steps."
            ),
            (
                "A simplicial wedge preserves rational convex polytopality; "
                "iterating realizes every positive integral multiwedge."
            ),
        ],
        "algebro_geometric_corollary": {
            "claim": (
                "A projective Q-factorial toric 9-fold has toric boundary-"
                "stratum dual complex equal to this nerve."
            ),
            "caveat": "No Delzant or smoothness assertion is made.",
        },
        "nonclaims": [
            "No Coxeter labels are supplied.",
            "No Lorentzian Coxeter Gram matrix is supplied.",
            (
                "Ordinary polytopality does not imply compactness or finite "
                "volume in hyperbolic space."
            ),
            "No new hyperbolic Coxeter polytope is claimed.",
        ],
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-18 polytopality certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-18 polytopality certificate digest mismatch")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-18 polytopality certificate is not frozen")
    return {
        "base_facet_count": rebuilt["audit"]["coordinate_support_facet_count"],
        "base_face_vector": rebuilt["audit"]["base_face_vector"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    del argv
    manifest = build()
    OUTPUT.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"output": str(OUTPUT), **verify(manifest)}, indent=2))


if __name__ == "__main__":
    main()
