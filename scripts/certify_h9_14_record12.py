#!/usr/bin/env python3
"""Freeze the exact rational ordinary-polytopal completion of H9/14 record 12."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_record12 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    TARGET_RECORD_INDEX,
    audit_record12_polytopality,
    target_facet_complements,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
PROOF_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record12.py"
HOSTILE_SOURCE = ROOT / "scripts/independent_h9_14_record12_audit.py"
SUPPORT_SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record12_polytopality.json"
EXPECTED_PROOF_SOURCE_SHA256 = (
    "d119c3cdc011d16142856dc5a25ce2d713ee480fa1289788d954e5928e548b22"
)
EXPECTED_HOSTILE_SOURCE_SHA256 = (
    "4e23de4d43843c6d9da446876ce97be50f60946e78bf333c054aee08b36e63cc"
)
EXPECTED_SUPPORT_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "cd98e07ccb15d04654c951e3271239be88d842b21c7bfa5f3c7f4c48d5875ab1"
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


def _require_frozen(path: Path, expected: str, description: str) -> str:
    observed = _file_sha256(path)
    if expected != "TO_BE_FROZEN" and observed != expected:
        raise ValueError(f"record-12 {description} source is not frozen")
    return observed


def _support_record() -> tuple[Mapping[str, object], str]:
    source_hash = _require_frozen(
        SUPPORT_SOURCE, EXPECTED_SUPPORT_SOURCE_SHA256, "support"
    )
    source = json.loads(SUPPORT_SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    survivors = [record for record in source["records"] if record["survives_second_layer"]]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("record-12 survivor index or graph6 changed")
    return record, source_hash


def build() -> dict[str, object]:
    proof_hash = _require_frozen(PROOF_SOURCE, EXPECTED_PROOF_SOURCE_SHA256, "proof")
    hostile_hash = _require_frozen(
        HOSTILE_SOURCE, EXPECTED_HOSTILE_SOURCE_SHA256, "hostile replay"
    )
    support, support_hash = _support_record()
    polytopality = audit_record12_polytopality()
    facets = target_facet_complements()
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        support["ridge_support_core_masks"],
        support["dehn_sommerville_schedule"],
        facets,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("record-12 selected nerve failed its global audit")
    if global_audit.face_vector != TARGET_F_VECTOR:
        raise AssertionError("record-12 selected face vector changed")
    if global_audit.missing_four_count != 2:
        raise AssertionError("record-12 selected schedule changed")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record12_exact_rational_polytopality",
        "claim": (
            "The selected H9/14 record-12 nerve is the boundary of a rational "
            "simplicial convex 9-polytope; dually it is the nerve of a rational "
            "simple convex 9-polytope with fourteen facets."
        ),
        "theorem_status": "certified_explicit_ordinary_polytopal",
        "target": {
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
            "missing_edges": [list(value) for value in TARGET_MISSING_EDGES],
            "missing_fours": [list(value) for value in TARGET_MISSING_FOURS],
            "missing_fives": [list(value) for value in TARGET_MISSING_FIVES],
            "facet_complement_masks": list(facets),
            "face_vector": list(TARGET_F_VECTOR),
            "missing_four_schedule": 2,
        },
        "sources": {
            "proof": {
                "path": "src/coxeter_polytopes/combinatorics/h9_14_record12.py",
                "sha256": proof_hash,
            },
            "hostile_replay": {
                "path": "scripts/independent_h9_14_record12_audit.py",
                "sha256": hostile_hash,
            },
            "ridge_support": {
                "path": "data/enumerations/h9_14_nerve_support.json",
                "sha256": support_hash,
            },
        },
        "proof_backend": (
            "standard-library exact arithmetic: Gale evenness, two abstract "
            "stellar subdivisions, exhaustive rational supporting-hyperplane "
            "determinants, exact radial-beyond identities, four simplicial-"
            "wedge minimal-nonface expansions, complete blocker enumeration, "
            "and the generic exact global face/link audit"
        ),
        "construction": (
            "Take C(8,5) in moment-curve order [0,1,2,3,4,5,6,7]. "
            "Stellar-subdivide triangle {0,2,4} by vertex 8 with epsilon 1/8, "
            "then triangle {0,1,8} by vertex 9 with epsilon 1/32. Apply one "
            "simplicial wedge at each of vertices 1,2,3,4 and use the frozen "
            "target relabeling."
        ),
        "polytopality_audit": polytopality.to_record(),
        "global_incidence_audit": global_audit.to_record(),
        "proof_lemmas": [
            (
                "The displayed rational radial-beyond placements realize the "
                "two boundary stellar subdivisions; every supporting "
                "hyperplane is recomputed with exact determinants."
            ),
            (
                "A simplicial wedge preserves rational convex polytopality and "
                "expands each minimal nonface containing its carrier across "
                "the two-vertex clone block."
            ),
        ],
        "algebro_geometric_corollary": {
            "claim": (
                "The normal fan of the dual rational simple 9-polytope gives "
                "a projective Q-factorial toric 9-fold whose toric boundary-"
                "stratum dual complex is this nerve."
            ),
            "caveat": "No Delzant or smoothness assertion is made.",
        },
        "nonclaims": [
            "No Coxeter labels or Lorentzian Coxeter Gram matrix are supplied.",
            (
                "Ordinary polytopality does not imply compactness or finite "
                "volume in hyperbolic space."
            ),
            "No new hyperbolic Coxeter polytope is claimed.",
            (
                "The positive witness settles only one completion of graph "
                "shadow record 12, not every completion of that shadow."
            ),
        ],
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-12 certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-12 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-12 certificate is not frozen")
    return {
        "facet_count": rebuilt["global_incidence_audit"]["facet_count"],
        "face_vector": rebuilt["global_incidence_audit"]["face_vector"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    if args.verify_only:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    else:
        manifest = build()
        OUTPUT.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps({"output": str(OUTPUT), **verify(manifest)}, indent=2))


if __name__ == "__main__":
    main()
