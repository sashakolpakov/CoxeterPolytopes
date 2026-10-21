"""Build and verify the solver-independent H9/14 record-74 exclusion."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record74 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_record74,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
BACKEND_SOURCE = ROOT / "data/enumerations/h9_14_facet_completion.json"
AUDIT_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record74.py"
OUTPUT = ROOT / "data/enumerations/h9_14_record74_certificate.json"
EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_BACKEND_SOURCE_SHA256 = (
    "2d329439e5baeffb574932ff7bbab2dd78daf706b2917a0674bdc9f01a2f5123"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "0a8812c4dbade45cabe3c544c968d265581e4d3f0fdd6b48673ff7847148ceb4"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    return sha256(
        json.dumps(critical, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _source_record(source: Mapping[str, object]) -> Mapping[str, object]:
    records = source.get("records")
    if not isinstance(records, list):
        raise ValueError("nerve-support source lacks records")
    survivors = [
        record
        for record in records
        if isinstance(record, Mapping) and record.get("survives_second_layer")
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("survivor index 74 has changed")
    return record


def _backend_corroboration() -> dict[str, object]:
    backend = json.loads(BACKEND_SOURCE.read_text(encoding="utf-8"))
    matching = [
        record
        for record in backend["records"]
        if record.get("index") == TARGET_RECORD_INDEX
    ]
    if len(matching) != 1 or matching[0]["solver"]["status"] != "unknown":
        raise ValueError("record-74 bounded backend disposition changed")
    solver = matching[0]["solver"]
    return {
        "path": "data/enumerations/h9_14_facet_completion.json",
        "sha256": _file_sha256(BACKEND_SOURCE),
        "status": solver["status"],
        "reason_unknown": solver["reason_unknown"],
        "backend": solver["backend"],
        "backend_version": solver["backend_version"],
        "assertion_count": solver["assertion_count"],
        "elapsed_seconds": solver["elapsed_seconds"],
        "trusted_for_proof": False,
    }


def build() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    audit = audit_record74(
        TARGET_GRAPH6, record["dehn_sommerville_schedule"]
    ).to_record()
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record74_solver_independent_exclusion",
        "claim": (
            "H9/14 survivor record 74 cannot be the nerve of a compact "
            "hyperbolic Coxeter 9-polytope."
        ),
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_unsat",
        "proof_backend": "exact local g-vectors and blocker divisibility",
        "source": {
            "path": "data/enumerations/h9_14_nerve_support.json",
            "sha256": _file_sha256(SOURCE),
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
        },
        "audit_source": {
            "path": "src/coxeter_polytopes/combinatorics/h9_14_record74.py",
            "sha256": _file_sha256(AUDIT_SOURCE),
        },
        "proof_dependencies": [
            "every compact Coxeter missing face is Lanner of order at most five",
            "strict generalized lower-bound inequalities in nonsimplex links",
            "Dehn--Sommerville face-vector schedules",
            (
                "maximal faces are complementary to inclusion-minimal "
                "transversals of the complete minimal-nonface clutter"
            ),
        ],
        "proof": audit,
        "backend_corroboration": _backend_corroboration(),
        "updated_theorem_partition": {
            "independently_certified_impossible": 118,
            "certified_explicit_witness": 2,
            "theorem_undecided": 22,
            "total_graph_shadows": 142,
        },
        "nonclaims": [
            "This excludes one graph shadow, not all compact H9/14 nerves.",
            "The archived bounded Z3 timeout contributes no logical step.",
            "No enumeration of A5 families or SAT/SMT backend is required.",
            "This is not a general non-existence theorem in dimension 9.",
        ],
    }
    manifest["certificate_sha256"] = _canonical_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("schema_version") != "1.0.0" or manifest.get("kind") != (
        "h9_14_record74_solver_independent_exclusion"
    ):
        raise ValueError("unexpected record-74 certificate metadata")
    if manifest.get("record_index") != TARGET_RECORD_INDEX or manifest.get(
        "graph6"
    ) != TARGET_GRAPH6:
        raise ValueError("record-74 target mismatch")
    if manifest.get("theorem_status") != "certified_unsat":
        raise ValueError("record-74 theorem status was not promoted")
    if _file_sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise ValueError("nerve-support source is not frozen")
    if _file_sha256(BACKEND_SOURCE) != EXPECTED_BACKEND_SOURCE_SHA256:
        raise ValueError("backend source is not frozen")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    expected_proof = audit_record74(
        TARGET_GRAPH6, record["dehn_sommerville_schedule"]
    ).to_record()
    if manifest.get("proof") != expected_proof:
        raise ValueError("record-74 structural proof does not verify")
    if manifest.get("audit_source") != {
        "path": "src/coxeter_polytopes/combinatorics/h9_14_record74.py",
        "sha256": _file_sha256(AUDIT_SOURCE),
    }:
        raise ValueError("record-74 audit-source provenance mismatch")
    if manifest.get("backend_corroboration") != _backend_corroboration():
        raise ValueError("record-74 backend corroboration mismatch")
    backend = manifest["backend_corroboration"]
    if not isinstance(backend, Mapping) or backend.get("trusted_for_proof") is not False:
        raise ValueError("bounded backend timeout was incorrectly trusted")
    digest = _canonical_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-74 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-74 certificate is not frozen")
    return expected_proof


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
    proof = verify(manifest)
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "certificate_sha256": manifest["certificate_sha256"],
                "verified": proof["verified"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
