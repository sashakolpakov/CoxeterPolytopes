#!/usr/bin/env python3
"""Freeze the exact rational ordinary-polytopal completion of H9/14 record 4."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import audit_global_facet_selection
from coxeter_polytopes.combinatorics.h9_14_record4 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    TARGET_RECORD_INDEX,
    audit_record4_polytopality,
    target_facet_complements,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
PROOF_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record4.py"
HOSTILE_SOURCE = ROOT / "scripts/independent_h9_14_record4_audit.py"
SUPPORT_SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record4_polytopality.json"
EXPECTED_PROOF_SOURCE_SHA256 = (
    "ed6aba49d3bed3c0149c3ddcaaa00c3c6c8b41eefa11da88ddb2b93c83fde43b"
)
EXPECTED_HOSTILE_SOURCE_SHA256 = (
    "329c4123c311824fee39f489e7551960af7898a1d4bb585c7657cc395f3b6ced"
)
EXPECTED_SUPPORT_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "0625af6609333638803df22868791179f3a659c47819bd762a07f1bd069cb1fc"
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
        raise ValueError(f"record-4 {description} source is not frozen")
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
        raise ValueError("record-4 survivor index or graph6 changed")
    return record, source_hash


def build() -> dict[str, object]:
    proof_hash = _require_frozen(PROOF_SOURCE, EXPECTED_PROOF_SOURCE_SHA256, "proof")
    hostile_hash = _require_frozen(
        HOSTILE_SOURCE, EXPECTED_HOSTILE_SOURCE_SHA256, "hostile replay"
    )
    support, support_hash = _support_record()
    polytopality = audit_record4_polytopality()
    facets = target_facet_complements()
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        support["ridge_support_core_masks"],
        support["dehn_sommerville_schedule"],
        facets,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("record-4 selected nerve failed its global audit")
    if global_audit.face_vector != TARGET_F_VECTOR or global_audit.missing_four_count != 2:
        raise AssertionError("record-4 selected schedule changed")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record4_exact_rational_polytopality",
        "claim": (
            "The selected H9/14 record-4 nerve is the boundary of a rational "
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
                "path": "src/coxeter_polytopes/combinatorics/h9_14_record4.py",
                "sha256": proof_hash,
            },
            "hostile_replay": {
                "path": "scripts/independent_h9_14_record4_audit.py",
                "sha256": hostile_hash,
            },
            "ridge_support": {
                "path": "data/enumerations/h9_14_nerve_support.json",
                "sha256": support_hash,
            },
        },
        "proof_backend": (
            "standard-library exact arithmetic: Gale evenness, two abstract "
            "edge stellar subdivisions, exhaustive rational supporting-"
            "hyperplane determinants, exact radial-beyond identities, four "
            "simplicial-wedge minimal-nonface expansions, complete blocker "
            "enumeration, and the generic exact global face/link audit"
        ),
        "construction": (
            "Take C(8,5) in moment-curve order [0,1,2,3,4,5,6,7]. "
            "Stellar-subdivide edge {0,1} by vertex 8 and edge {0,7} by "
            "vertex 9, using epsilon 1/16 at both steps. Apply one simplicial "
            "wedge at each of vertices 3,4,5,6 and use the frozen relabeling."
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
                "shadow record 4, not every completion of that shadow."
            ),
        ],
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-4 certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-4 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-4 certificate is not frozen")
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
