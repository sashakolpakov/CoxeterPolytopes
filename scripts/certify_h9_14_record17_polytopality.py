#!/usr/bin/env python3
"""Freeze the exact rational ordinary-polytopality proof for record 17."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record17_polytopality import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_RECORD_INDEX,
    audit_record17_polytopality,
)


ROOT = Path(__file__).resolve().parents[1]
PROOF_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record17_polytopality.py"
OUTPUT = ROOT / "data/enumerations/h9_14_record17_polytopality.json"
EXPECTED_PROOF_SOURCE_SHA256 = "00722328898ff4b0fe7cdf1efade2264300680cf8abccdbbc55388215efd9aae"
EXPECTED_CERTIFICATE_SHA256 = "d17332be82e6bddb0ebfe320d678238cc0fa182f2a1b03d81931217826c72e31"


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
        raise ValueError("record-17 polytopality proof source is not frozen")
    audit = audit_record17_polytopality()
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record17_exact_rational_polytopality",
        "claim": (
            "The selected H9/14 record-17 nerve is the boundary of a rational "
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
            "path": "src/coxeter_polytopes/combinatorics/h9_14_record17_polytopality.py",
            "sha256": source_hash,
        },
        "proof_backend": (
            "standard-library exact rational arithmetic: Gale evenness, two "
            "abstract stellar subdivisions, exhaustive supporting-hyperplane "
            "determinants, and inverse-multiwedge minimal-nonface equality"
        ),
        "audit": audit.to_record(),
        "construction": (
            "Take C(8,5) on moment-curve labels [1,0,11,15,17,14,16,2]. "
            "Stellar-subdivide {1,14,15} by vertex 13 and then {2,14,15} by "
            "vertex 12. Replace carriers 14,15,16,17 by clone pairs "
            "{3,4},{5,6},{7,8},{9,10}. The four simplicial wedges raise "
            "dimension five to nine and reproduce the target clutter exactly."
        ),
        "proof_lemmas": [
            (
                "The boundary stellar subdivision at a proper face of a rational "
                "polytope is rational-polytopal; here the supplied coordinates "
                "are independently checked against every supporting hyperplane."
            ),
            (
                "A simplicial wedge preserves rational convex polytopality; "
                "iterating realizes every positive integral multiwedge."
            ),
        ],
        "algebro_geometric_corollary": {
            "claim": (
                "A projective Q-factorial toric 9-fold has toric boundary-stratum "
                "dual complex equal to this nerve."
            ),
            "caveat": "No Delzant or smoothness assertion is made.",
        },
        "nonclaims": [
            "No Coxeter labels are supplied.",
            "No Lorentzian Coxeter Gram matrix is supplied.",
            "Ordinary polytopality does not imply compactness or finite volume in hyperbolic space.",
            "No new hyperbolic Coxeter polytope is claimed.",
        ],
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-17 polytopality certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-17 polytopality certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-17 polytopality certificate is not frozen")
    return {
        "base_facet_count": rebuilt["audit"]["coordinate_support_facet_count"],
        "base_face_vector": rebuilt["audit"]["base_face_vector"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    del argv
    manifest = build()
    OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), **verify(manifest)}, indent=2))


if __name__ == "__main__":
    main()
