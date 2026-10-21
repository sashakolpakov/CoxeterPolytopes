"""Build and verify the solver-independent H9/14 record-54 exclusion."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record54 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_record54,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
BACKEND_SOURCE = ROOT / "data/enumerations/h9_14_facet_completion.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record54_certificate.json"
EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_BACKEND_SOURCE_SHA256 = (
    "2d329439e5baeffb574932ff7bbab2dd78daf706b2917a0674bdc9f01a2f5123"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "911955924b1ef679e5a12c0c22a0c11184f0b1ed35ca9b6d2459502ef0b2f295"
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
        json.dumps(critical, sort_keys=True, separators=(",", ":")).encode("utf-8")
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
    if TARGET_RECORD_INDEX >= len(survivors):
        raise ValueError("nerve-support source lacks survivor record 54")
    record = survivors[TARGET_RECORD_INDEX]
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("survivor index 54 has changed")
    return record


def _backend_corroboration() -> dict[str, object]:
    backend = json.loads(BACKEND_SOURCE.read_text(encoding="utf-8"))
    records = backend.get("records")
    if not isinstance(records, list):
        raise ValueError("facet-completion archive lacks records")
    matching = [
        record
        for record in records
        if isinstance(record, Mapping) and record.get("index") == TARGET_RECORD_INDEX
    ]
    if len(matching) != 1:
        raise ValueError("facet-completion archive does not isolate record 54")
    solver = matching[0].get("solver")
    if not isinstance(solver, Mapping) or solver.get("status") != "unsat":
        raise ValueError("frozen record-54 backend disposition changed")
    return {
        "path": "data/enumerations/h9_14_facet_completion.json",
        "sha256": _file_sha256(BACKEND_SOURCE),
        "status": solver.get("status"),
        "backend": solver.get("backend"),
        "backend_version": solver.get("backend_version"),
        "assertion_count": solver.get("assertion_count"),
        "trusted_for_proof": False,
    }


def build() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    schedules = record.get("dehn_sommerville_schedule")
    if not isinstance(schedules, list):
        raise ValueError("record-54 source lacks schedules")
    proof = audit_record54(TARGET_GRAPH6, schedules).to_record()
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record54_solver_independent_exclusion",
        "claim": (
            "H9/14 survivor record 54 cannot be the nerve of a compact "
            "hyperbolic Coxeter 9-polytope."
        ),
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_unsat",
        "proof_backend": "finite integer enumeration; no SAT/SMT backend",
        "source": {
            "path": "data/enumerations/h9_14_nerve_support.json",
            "sha256": _file_sha256(SOURCE),
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
        },
        "proof_dependencies": [
            (
                "every missing face of a compact Coxeter polytope is Lanner "
                "and has order at most five"
            ),
            (
                "the generalized lower-bound equality theorem makes g4 "
                "strictly positive in every nonsimplex vertex link"
            ),
            "the g-theorem Macaulay upper inequalities for vertex links",
            "Dehn--Sommerville face-vector schedules",
            (
                "maximal faces are complementary to inclusion-minimal "
                "transversals of the complete minimal-nonface clutter"
            ),
        ],
        "proof": proof,
        "backend_corroboration": _backend_corroboration(),
        "updated_theorem_partition": {
            "independently_certified_impossible": 116,
            "certified_explicit_witness": 2,
            "theorem_undecided": 24,
            "total_graph_shadows": 142,
        },
        "nonclaims": [
            "This excludes one graph shadow, not all compact H9/14 nerves.",
            (
                "The archived Z3 UNSAT result corroborates the theorem but "
                "is not used by the proof or verifier."
            ),
            (
                "No Coxeter labeling or Lorentzian Gram analysis is needed "
                "because this shadow already fails necessary nerve conditions."
            ),
            "This is not a general non-existence theorem in dimension 9.",
        ],
    }
    manifest["certificate_sha256"] = _canonical_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("schema_version") != "1.0.0" or manifest.get("kind") != (
        "h9_14_record54_solver_independent_exclusion"
    ):
        raise ValueError("unexpected record-54 certificate metadata")
    if manifest.get("record_index") != TARGET_RECORD_INDEX or manifest.get(
        "graph6"
    ) != TARGET_GRAPH6:
        raise ValueError("record-54 certificate target mismatch")
    if manifest.get("theorem_status") != "certified_unsat":
        raise ValueError("record-54 theorem status was not promoted")

    source_sha = _file_sha256(SOURCE)
    if source_sha != EXPECTED_SOURCE_SHA256:
        raise ValueError("nerve-support source is not frozen")
    source_metadata = manifest.get("source")
    if not isinstance(source_metadata, Mapping) or source_metadata.get(
        "sha256"
    ) != source_sha:
        raise ValueError("record-54 source provenance mismatch")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    expected_proof = audit_record54(
        TARGET_GRAPH6, record["dehn_sommerville_schedule"]
    ).to_record()
    if manifest.get("proof") != expected_proof:
        raise ValueError("record-54 arithmetic proof does not verify")

    if _file_sha256(BACKEND_SOURCE) != EXPECTED_BACKEND_SOURCE_SHA256:
        raise ValueError("facet-completion corroboration source is not frozen")
    if manifest.get("backend_corroboration") != _backend_corroboration():
        raise ValueError("record-54 backend corroboration mismatch")
    backend = manifest["backend_corroboration"]
    if not isinstance(backend, Mapping) or backend.get("trusted_for_proof") is not False:
        raise ValueError("backend evidence was incorrectly promoted")

    digest = _canonical_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-54 certificate digest mismatch")
    if digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-54 certificate is not frozen")
    return dict(expected_proof)


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
