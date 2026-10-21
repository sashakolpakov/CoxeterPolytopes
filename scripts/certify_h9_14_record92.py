"""Build and verify the solver-independent H9/14 record-92 exclusion."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_record92,
)
from coxeter_polytopes.combinatorics.nerve_support import (
    verify_support_manifest,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
BACKEND_SOURCE = ROOT / "data" / "enumerations" / "h9_14_schedule_split.json"
OUTPUT = ROOT / "data" / "enumerations" / "h9_14_record92_certificate.json"
EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_BACKEND_SOURCE_SHA256 = (
    "f9e3d149f4da9bc05220c84c1ffdde5ba8cb8ed9d388df94992d20d5f200b6f5"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "c0e0903caa31469d3bd95d1e10cd42a62704a1de1859cfe370749065f05a5c25"
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
    encoded = json.dumps(
        critical, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return sha256(encoded).hexdigest()


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
        raise ValueError("nerve-support source lacks survivor record 92")
    record = survivors[TARGET_RECORD_INDEX]
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("survivor index 92 has changed")
    return record


def _backend_corroboration() -> dict[str, object]:
    backend = json.loads(BACKEND_SOURCE.read_text(encoding="utf-8"))
    records = backend.get("records")
    if not isinstance(records, list):
        raise ValueError("schedule-split archive lacks records")
    matching = [
        record
        for record in records
        if isinstance(record, Mapping)
        and record.get("index") == TARGET_RECORD_INDEX
    ]
    if len(matching) != 1:
        raise ValueError("schedule-split archive does not isolate record 92")
    record = matching[0]
    branches = record.get("branches")
    if not isinstance(branches, list):
        raise ValueError("record-92 backend entry lacks branches")
    statuses = {
        int(branch["missing_four_count"]): branch["solver"]["status"]
        for branch in branches
    }
    if statuses != {0: "unsat", 1: "unsat", 2: "unsat"}:
        raise ValueError("frozen record-92 backend branches have changed")
    return {
        "path": "data/enumerations/h9_14_schedule_split.json",
        "sha256": _file_sha256(BACKEND_SOURCE),
        "branch_statuses": {
            str(a): statuses[a] for a in sorted(statuses)
        },
        "aggregate_backend_status": record.get("aggregate_backend_status"),
        "trusted_for_proof": False,
    }


def build() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    graph6 = record.get("graph6")
    schedules = record.get("dehn_sommerville_schedule")
    if not isinstance(graph6, str) or not isinstance(schedules, list):
        raise ValueError("record-92 source fields are malformed")
    proof = audit_record92(graph6, schedules).to_record()
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record92_solver_independent_exclusion",
        "claim": (
            "H9/14 survivor record 92 cannot be the nerve of a compact "
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
                "strictly positive for every nonsimplex vertex link"
            ),
            "Dehn--Sommerville face-vector schedules",
            (
                "maximal faces are complementary to inclusion-minimal "
                "transversals of the minimal-nonface clutter"
            ),
        ],
        "proof": proof,
        "backend_corroboration": _backend_corroboration(),
        "updated_theorem_partition": {
            "independently_certified_impossible": 114,
            "certified_explicit_witness": 2,
            "theorem_undecided": 26,
            "total_graph_shadows": 142,
        },
        "nonclaims": [
            "This excludes one graph shadow, not all compact H9/14 nerves.",
            (
                "The archived Z3 UNSAT run corroborates the theorem but is "
                "not used by the proof or verifier."
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
        "h9_14_record92_solver_independent_exclusion"
    ):
        raise ValueError("unexpected record-92 certificate metadata")
    if manifest.get("record_index") != TARGET_RECORD_INDEX or manifest.get(
        "graph6"
    ) != TARGET_GRAPH6:
        raise ValueError("record-92 certificate target mismatch")
    if manifest.get("theorem_status") != "certified_unsat":
        raise ValueError("record-92 theorem status was not promoted")

    source_sha = _file_sha256(SOURCE)
    if source_sha != EXPECTED_SOURCE_SHA256:
        raise ValueError("nerve-support source is not frozen")
    source_metadata = manifest.get("source")
    if not isinstance(source_metadata, Mapping) or source_metadata.get(
        "sha256"
    ) != source_sha:
        raise ValueError("record-92 source provenance mismatch")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    expected_proof = audit_record92(
        str(record["graph6"]), record["dehn_sommerville_schedule"]
    ).to_record()
    if manifest.get("proof") != expected_proof:
        raise ValueError("record-92 arithmetic proof does not verify")

    if _file_sha256(BACKEND_SOURCE) != EXPECTED_BACKEND_SOURCE_SHA256:
        raise ValueError("schedule-split corroboration source is not frozen")
    if manifest.get("backend_corroboration") != _backend_corroboration():
        raise ValueError("record-92 backend corroboration mismatch")
    backend = manifest["backend_corroboration"]
    if not isinstance(backend, Mapping) or backend.get("trusted_for_proof") is not False:
        raise ValueError("backend evidence was incorrectly promoted")

    digest = _canonical_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-92 certificate digest mismatch")
    if digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-92 certificate is not frozen")
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
