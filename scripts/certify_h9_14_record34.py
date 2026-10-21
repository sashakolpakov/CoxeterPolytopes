#!/usr/bin/env python3
"""Build and verify the solver-independent H9/14 record-34 exclusion."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record34 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_record34_positive_branches,
)
from coxeter_polytopes.combinatorics.h9_14_record34_a0 import verify_manifest
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
BACKEND_SOURCE = ROOT / "data/enumerations/h9_14_facet_completion.json"
A0_ARCHIVE = ROOT / "data/enumerations/h9_14_record34_a0_census.json"
POSITIVE_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record34.py"
)
A0_REPLAY_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record34_a0.py"
)
A0_ENUMERATOR_SOURCE = ROOT / "scripts/certify_h9_14_record34_a0.cpp"
OUTPUT = ROOT / "data/enumerations/h9_14_record34_certificate.json"

EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_BACKEND_SOURCE_SHA256 = (
    "2d329439e5baeffb574932ff7bbab2dd78daf706b2917a0674bdc9f01a2f5123"
)
EXPECTED_A0_ARCHIVE_SHA256 = (
    "a364009ba514668d21ad3664cd3a798ad820e9d3d46d8e21256fadf14a16bb4a"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "79d4de8b4c9d463bc0aca44bfd78a3c8eb1ee6bfdc211f90b47ad661e5837db2"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _certificate_digest(manifest: Mapping[str, object]) -> str:
    critical = dict(manifest)
    critical.pop("certificate_sha256", None)
    return _canonical_sha256(critical)


def _source_record(source: Mapping[str, object]) -> Mapping[str, object]:
    survivors = [
        record
        for record in source["records"]
        if isinstance(record, Mapping) and record.get("survives_second_layer")
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("second-layer survivor record 34 changed")
    return record


def _backend_corroboration() -> dict[str, object]:
    backend = json.loads(BACKEND_SOURCE.read_text(encoding="utf-8"))
    records = [
        record
        for record in backend["records"]
        if record.get("index") == TARGET_RECORD_INDEX
    ]
    if len(records) != 1 or records[0]["solver"]["status"] != "unknown":
        raise ValueError("record-34 bounded-backend disposition changed")
    solver = records[0]["solver"]
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


def _a0_summary(document: Mapping[str, object]) -> dict[str, object]:
    totals = document["totals"]
    if (
        totals["genbg_representatives"] != 1_288_153
        or totals["active_colourings"] != 26_598_204
        or totals["distinct_row_states"] != 25_501_701
        or totals["no_small_transversal_states"] != 3_101_794
        or totals["target_114_states"] != 18
        or totals["pure_target_states"] != 18
        or totals["ridge_valid_states"] != 0
    ):
        raise AssertionError("record-34 a=0 aggregate census changed")
    return {
        "path": "data/enumerations/h9_14_record34_a0_census.json",
        "sha256": _file_sha256(A0_ARCHIVE),
        "certificate_sha256": document["certificate_sha256"],
        "generator_provenance": document["generator_provenance"],
        "strata": document["strata"],
        "totals": totals,
        "target_profile_counts": document["target_profile_counts"],
        "canonical_orbit_count": document["canonical_orbit_count"],
        "canonical_orbit_multiplicities": document[
            "canonical_orbit_multiplicities"
        ],
        "conclusion": (
            "all 18 pure target-count blockers have a ridge of degree "
            "different from two; no a=0 state is a pseudomanifold"
        ),
    }


def build() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    positive = audit_record34_positive_branches(
        TARGET_GRAPH6, record["dehn_sommerville_schedule"]
    )
    a0_document = verify_manifest(A0_ARCHIVE)
    a0 = _a0_summary(a0_document)
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record34_solver_independent_exclusion",
        "claim": (
            "H9/14 survivor record 34 cannot be the nerve of a compact "
            "hyperbolic Coxeter 9-polytope."
        ),
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_unsat",
        "proof_backend": (
            "exact local g-vectors, forced minimal-nonface shadows, a "
            "labelled a=1 blocker census, nauty genbg bicoloured-graph "
            "generation, and an independent 2^14 blocker/ridge replay; "
            "no SAT/SMT disposition is trusted"
        ),
        "source": {
            "path": "data/enumerations/h9_14_nerve_support.json",
            "sha256": _file_sha256(SOURCE),
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
        },
        "proof_sources": [
            {
                "path": "src/coxeter_polytopes/combinatorics/h9_14_record34.py",
                "sha256": _file_sha256(POSITIVE_SOURCE),
            },
            {
                "path": (
                    "src/coxeter_polytopes/combinatorics/"
                    "h9_14_record34_a0.py"
                ),
                "sha256": _file_sha256(A0_REPLAY_SOURCE),
            },
            {
                "path": "scripts/certify_h9_14_record34_a0.cpp",
                "sha256": _file_sha256(A0_ENUMERATOR_SOURCE),
            },
        ],
        "proof_dependencies": [
            "Burcroff's no-missing-triples theorem for compact H9/14 nerves",
            "every compact Coxeter missing face is Lanner of order at most five",
            "strict generalized lower-bound and local Macaulay inequalities",
            "Dehn--Sommerville face-vector schedules",
            (
                "maximal faces are complementary to inclusion-minimal "
                "transversals of the complete minimal-nonface clutter"
            ),
            "every ridge of a simplicial sphere has degree exactly two",
        ],
        "positive_a_exclusion": positive.to_record(),
        "a0_exhaustive_census": a0,
        "backend_corroboration": _backend_corroboration(),
        "updated_theorem_partition": {
            "independently_certified_impossible": 121,
            "certified_explicit_ordinary_polytopal": 6,
            "theorem_undecided": 15,
            "total_graph_shadows": 142,
        },
        "nonclaims": [
            "This excludes one graph shadow, not all compact H9/14 nerves.",
            "The archived bounded Z3 timeout contributes no logical step.",
            "No Coxeter labels or Lorentzian Gram matrices are classified.",
            "This is not a general non-existence theorem in dimension 9.",
        ],
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_record34_solver_independent_exclusion":
        raise ValueError("wrong record-34 certificate kind")
    if manifest.get("record_index") != TARGET_RECORD_INDEX or manifest.get(
        "graph6"
    ) != TARGET_GRAPH6:
        raise ValueError("record-34 certificate target changed")
    if _file_sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise ValueError("record-34 support source is not frozen")
    if _file_sha256(BACKEND_SOURCE) != EXPECTED_BACKEND_SOURCE_SHA256:
        raise ValueError("record-34 backend source is not frozen")
    if EXPECTED_A0_ARCHIVE_SHA256 != "TO_BE_FROZEN" and _file_sha256(
        A0_ARCHIVE
    ) != EXPECTED_A0_ARCHIVE_SHA256:
        raise ValueError("record-34 a=0 archive is not frozen")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-34 certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-34 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-34 certificate is not frozen")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_unsat",
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
