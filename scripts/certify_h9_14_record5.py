#!/usr/bin/env python3
"""Freeze the exact ordinary-polytopal completion of H9/14 record 5."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_record5 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    TARGET_RECORD_INDEX,
    audit_record5_polytopality,
    target_facet_complements,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
PROOF_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record5.py"
SUPPORT_SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record5_polytopality.json"
EXPECTED_PROOF_SOURCE_SHA256 = (
    "37a9c0b6139bd03f1dddb0f2e6ba0298ef4a0c6566ecb4718d62940bb551e598"
)
EXPECTED_SUPPORT_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "6cdd1b7951093bb6765a3d77376e3bbad50d9f470d1752d0abc80eb061cdb015"
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


def _support_record() -> tuple[Mapping[str, object], str]:
    source_hash = _file_sha256(SUPPORT_SOURCE)
    if source_hash != EXPECTED_SUPPORT_SOURCE_SHA256:
        raise ValueError("record-5 support source is not frozen")
    source = json.loads(SUPPORT_SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("record-5 survivor index or graph6 changed")
    return record, source_hash


def build() -> dict[str, object]:
    proof_hash = _file_sha256(PROOF_SOURCE)
    if proof_hash != EXPECTED_PROOF_SOURCE_SHA256:
        raise ValueError("record-5 proof source is not frozen")
    support, support_hash = _support_record()
    polytopality = audit_record5_polytopality()
    facets = target_facet_complements()
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        support["ridge_support_core_masks"],
        support["dehn_sommerville_schedule"],
        facets,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("record-5 selected nerve failed its global audit")
    if global_audit.face_vector != TARGET_F_VECTOR:
        raise AssertionError("record-5 selected face vector changed")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record5_exact_rational_polytopality",
        "claim": (
            "The selected H9/14 record-5 nerve is the boundary of a rational "
            "simplicial convex 9-polytope; dually it is the nerve of a rational "
            "simple convex 9-polytope with fourteen facets."
        ),
        "theorem_status": "certified_explicit_ordinary_polytopal",
        "target": {
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
            "missing_edges": [list(face) for face in TARGET_MISSING_EDGES],
            "missing_fours": [list(face) for face in TARGET_MISSING_FOURS],
            "missing_fives": [list(face) for face in TARGET_MISSING_FIVES],
            "facet_complement_masks": list(facets),
            "face_vector": list(TARGET_F_VECTOR),
            "missing_four_schedule": 3,
        },
        "sources": {
            "proof": {
                "path": "src/coxeter_polytopes/combinatorics/h9_14_record5.py",
                "sha256": proof_hash,
            },
            "ridge_support": {
                "path": "data/enumerations/h9_14_nerve_support.json",
                "sha256": support_hash,
            },
        },
        "proof_backend": (
            "standard-library exact arithmetic: Gale evenness, two abstract "
            "edge stellar subdivisions, exhaustive rational supporting-"
            "hyperplane determinants, radial-beyond identities, four "
            "simplicial wedges with explicit relabelling, minimal-nonface "
            "equality, complete blocker enumeration, and the generic exact "
            "global face/link audit"
        ),
        "construction": (
            "Take C(8,5) in moment-curve order [0,1,2,3,4,5,6,7]. "
            "Stellar-subdivide edge {0,2} by vertex 8 and then edge {1,3} "
            "by vertex 9. Replace carriers 0,2,4,8 by clone pairs "
            "{0,10},{2,11},{4,12},{8,13}; then relabel construction vertices "
            "by {0:3,1:0,2:4,3:12,4:5,5:1,6:2,7:6,8:7,9:13,10:8,"
            "11:9,12:10,13:11}."
        ),
        "polytopality_audit": polytopality.to_record(),
        "global_incidence_audit": global_audit.to_record(),
        "proof_lemmas": [
            (
                "A sufficiently small rational radial-beyond placement "
                "realizes boundary stellar subdivision; epsilons 1/16 and "
                "1/32 are verified by exhaustive exact determinants."
            ),
            (
                "A simplicial wedge preserves rational convex polytopality, "
                "and replaces every minimal nonface containing its carrier "
                "by the same nonface containing both clones."
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
        raise ValueError("record-5 certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-5 certificate digest mismatch")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-5 certificate is not frozen")
    return {
        "facet_count": rebuilt["global_incidence_audit"]["facet_count"],
        "face_vector": rebuilt["global_incidence_audit"]["face_vector"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    del argv
    manifest = build()
    OUTPUT.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"output": str(OUTPUT), **verify(manifest)}, indent=2))


if __name__ == "__main__":
    main()
