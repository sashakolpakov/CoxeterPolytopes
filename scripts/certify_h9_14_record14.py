#!/usr/bin/env python3
"""Freeze and verify one ordinary-polytopal H9/14 record-14 completion."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_record14 import (
    TARGET_F_VECTOR,
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
    TARGET_MISSING_FIVES,
    TARGET_MISSING_FOURS,
    TARGET_RECORD_INDEX,
    audit_record14_polytopality,
    mask,
    record14_facet_complements,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest
from scripts.independent_h9_14_record14_audit import (
    EXPECTED_SUMMARY_SHA256,
    audit as independent_audit,
)


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record14_polytopality.json"
PROOF_SOURCES = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record14.py",
    ROOT / "src/coxeter_polytopes/combinatorics/facet_completion.py",
    ROOT / "src/coxeter_polytopes/combinatorics/local_link.py",
    ROOT / "src/coxeter_polytopes/combinatorics/missing_edges.py",
    ROOT / "src/coxeter_polytopes/combinatorics/nerve_support.py",
    ROOT / "src/coxeter_polytopes/combinatorics/simplicial_homology.py",
    ROOT / "scripts/independent_h9_14_record14_audit.py",
)

EXPECTED_SUPPORT_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_PROOF_SOURCE_SHA256 = {
    "src/coxeter_polytopes/combinatorics/h9_14_record14.py": (
        "9c83d807df43f364455164e7a120c2600558fb965400b2aa61d166de2a7c1891"
    ),
    "src/coxeter_polytopes/combinatorics/facet_completion.py": (
        "5f2d078cfe400b39f3e194352f7a0c6d3bd8cb6bb8e7cee37be159e8766d9400"
    ),
    "src/coxeter_polytopes/combinatorics/local_link.py": (
        "8f8169edab8b87aff682785667c2a8eab5c52b92a7c68cef130aeec4d3ae807d"
    ),
    "src/coxeter_polytopes/combinatorics/missing_edges.py": (
        "231e85fa98e27e864f997cd4ebe59136707f6aa8dd2d8136b80cff543dfbfbb8"
    ),
    "src/coxeter_polytopes/combinatorics/nerve_support.py": (
        "aad8fb25c2f9a56e45121fa1d1cbbdedfdff8ba1f67e379be74eed73940ad2a2"
    ),
    "src/coxeter_polytopes/combinatorics/simplicial_homology.py": (
        "547f51b58c1146aac80c6e08910754bf93ace5fa808dae5a342da7e042958110"
    ),
    "scripts/independent_h9_14_record14_audit.py": (
        "14cdcf1236eeda0f92fc857b5ef98c53cf513911b6d4188af58a23558dbef0d4"
    ),
}
EXPECTED_CERTIFICATE_SHA256 = (
    "0fba773d1449780048625c0243eb4740efc0ccd23d87b266e5f3ca06deb7edb1"
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def certificate_digest(payload: Mapping[str, object]) -> str:
    normalized = dict(payload)
    normalized.pop("certificate_sha256", None)
    return sha256(
        json.dumps(normalized, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def support_record() -> tuple[dict[str, object], str]:
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    verify_support_manifest(support)
    survivors = [
        record for record in support["records"] if record["survives_second_layer"]
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("second-layer survivor record 14 changed")
    return record, file_sha256(SUPPORT)


def proof_source_records() -> list[dict[str, str]]:
    return [
        {"path": str(path.relative_to(ROOT)), "sha256": file_sha256(path)}
        for path in PROOF_SOURCES
    ]


def build() -> dict[str, object]:
    source, support_sha = support_record()
    facets = record14_facet_complements()
    polytopality = audit_record14_polytopality()
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        source["ridge_support_core_masks"],
        source["dehn_sommerville_schedule"],
        facets,
    )
    if not global_audit.passes_connected_homology_audit:
        raise AssertionError("record-14 completion failed the global/link audit")
    if global_audit.face_vector != TARGET_F_VECTOR:
        raise AssertionError("record-14 completion face vector changed")
    expected_nonfaces = {
        2: tuple(sorted(mask(face) for face in TARGET_MISSING_EDGES)),
        4: tuple(sorted(mask(face) for face in TARGET_MISSING_FOURS)),
        5: tuple(sorted(mask(face) for face in TARGET_MISSING_FIVES)),
    }
    if global_audit.minimal_nonfaces_by_size != expected_nonfaces:
        raise AssertionError("record-14 completion minimal nonfaces changed")
    independent = independent_audit()
    if independent["summary_sha256"] != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("record-14 hostile replay digest changed")
    facet_digest = sha256(
        ("\n".join(map(str, facets)) + "\n").encode()
    ).hexdigest()
    if facet_digest != independent["facet_complement_sha256"]:
        raise AssertionError("record-14 independent facet transcript differs")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record14_ordinary_polytopal_completion",
        "claim": (
            "The H9/14 graph shadow at second-layer record 14 has one "
            "explicit completion that is the boundary of a convex simplicial "
            "9-polytope, obtained as a five-fold simplicial multiwedge of an "
            "explicitly Gale-realized simplicial 4-polytope."
        ),
        "status": "certified_ordinary_polytopal_completion",
        "scope": (
            "One explicitly displayed completion of record 14 only; this does "
            "not classify every completion of the graph shadow."
        ),
        "nonclaims": [
            "No Coxeter labels are supplied or classified.",
            "No Lorentzian Gram matrix or hyperbolic realization is supplied.",
            "No compactness or finite-volume hyperbolic conclusion is made.",
            "No whole-shadow uniqueness or classification is claimed.",
            "The archived bounded Z3 timeout contributes no logical step.",
        ],
        "source": {
            "support_path": "data/enumerations/h9_14_nerve_support.json",
            "support_sha256": support_sha,
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
        },
        "proof_sources": proof_source_records(),
        "minimal_nonfaces": {
            "2": [list(face) for face in TARGET_MISSING_EDGES],
            "4": [list(face) for face in TARGET_MISSING_FOURS],
            "5": [list(face) for face in TARGET_MISSING_FIVES],
        },
        "selected_facet_complement_masks": list(facets),
        "facet_complement_sha256": facet_digest,
        "global_selection_audit": global_audit.to_record(),
        "ordinary_polytopality_audit": polytopality.to_record(),
        "independent_hostile_replay": independent,
        "construction": {
            "base_dimension": 4,
            "base_vertex_count": 9,
            "simplicial_wedge_count": 5,
            "target_dimension": 9,
            "target_vertex_count": 14,
            "exact_minimal_nonface_match": True,
        },
        "vinberg_scope": {
            "coxeter_diagram_present": False,
            "elliptic_subdiagram_test_run": False,
            "parabolic_subdiagram_test_run": False,
            "finite_volume_criterion_checked": False,
            "reason": (
                "Vinberg's compactness/finite-volume criteria apply only after "
                "Coxeter labels and a Lorentzian Gram matrix are supplied."
            ),
        },
    }
    manifest = json.loads(json.dumps(manifest, sort_keys=True))
    manifest["certificate_sha256"] = certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_record14_ordinary_polytopal_completion":
        raise ValueError("wrong record-14 certificate kind")
    if manifest.get("status") != "certified_ordinary_polytopal_completion":
        raise ValueError("record-14 certificate status changed")
    observed_support = file_sha256(SUPPORT)
    if (
        EXPECTED_SUPPORT_SHA256 != "TO_BE_FROZEN"
        and observed_support != EXPECTED_SUPPORT_SHA256
    ):
        raise ValueError("record-14 support source changed")
    for path in PROOF_SOURCES:
        relative = str(path.relative_to(ROOT))
        observed = file_sha256(path)
        expected = EXPECTED_PROOF_SOURCE_SHA256[relative]
        if expected != "TO_BE_FROZEN" and observed != expected:
            raise ValueError(f"record-14 frozen proof source changed: {relative}")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-14 certificate did not reproduce")
    digest = certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-14 certificate digest mismatch")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-14 certificate is not frozen")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "status": "certified_ordinary_polytopal_completion",
        "facet_count": len(rebuilt["selected_facet_complement_masks"]),
        "face_vector": rebuilt["global_selection_audit"]["face_vector"],
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
    print(json.dumps(verify(manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
