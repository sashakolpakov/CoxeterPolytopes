#!/usr/bin/env python3
"""Build and verify the analytic H9/14 record-16 ``a=0`` certificate."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record16 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_record16_a0,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.nerve_support import (
    verify_support_manifest,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
BACKEND_SOURCE = ROOT / "data/enumerations/h9_14_facet_completion.json"
PROOF_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record16.py"
)
OUTPUT = ROOT / "data/enumerations/h9_14_record16_a0_certificate.json"

EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_BACKEND_SOURCE_SHA256 = (
    "2d329439e5baeffb574932ff7bbab2dd78daf706b2917a0674bdc9f01a2f5123"
)
EXPECTED_PROOF_SOURCE_SHA256 = (
    "fd2d02c223d683d9e72ad2e6e72edd9d928c6b0a4fe2525dc254b3b5bab8e564"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "883e2bff979d1c6681950c4e2027b1c5a5fdb5856b09542af56a166a63975965"
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
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
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
    record = survivors[TARGET_RECORD_INDEX]
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("second-layer survivor record 16 changed")
    return record


def _backend_corroboration() -> dict[str, object]:
    backend = json.loads(BACKEND_SOURCE.read_text(encoding="utf-8"))
    records = [
        record
        for record in backend["records"]
        if record.get("index") == TARGET_RECORD_INDEX
    ]
    if len(records) != 1 or records[0]["solver"]["status"] != "unknown":
        raise ValueError("record-16 bounded-backend disposition changed")
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


def build() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    schedules = record["dehn_sommerville_schedule"]
    a0 = audit_record16_a0(TARGET_GRAPH6, schedules)
    schedule_table = [
        {
            "missing_four_count": schedule["missing_four_count"],
            "missing_five_count_offset": schedule[
                "missing_five_count_offset"
            ],
            "facet_count": schedule["f_vector"][-1],
            "summed_vertex_link_g4": summed_vertex_link_g4(
                schedule["f_vector"]
            ),
        }
        for schedule in schedules
    ]
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record16_a0_analytic_exclusion",
        "claim": (
            "The a=0 Dehn--Sommerville branch of H9/14 survivor record "
            "16 cannot be the nerve of a compact hyperbolic Coxeter "
            "9-polytope."
        ),
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_branch_unsat",
        "excluded_missing_four_counts": [0, 5, 6],
        "remaining_missing_four_counts": [1, 2, 3, 4],
        "proof_backend": (
            "exact blocker duality, complete symbolic ridge incidence, "
            "Boolean truth tables, and a Cauchy--Schwarz bound; no "
            "SAT/SMT or exhaustive missing-five census"
        ),
        "source": {
            "path": "data/enumerations/h9_14_nerve_support.json",
            "sha256": _file_sha256(SOURCE),
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
        },
        "proof_sources": [
            {
                "path": (
                    "src/coxeter_polytopes/combinatorics/"
                    "h9_14_record16.py"
                ),
                "sha256": _file_sha256(PROOF_SOURCE),
            }
        ],
        "proof_dependencies": [
            "Burcroff's no-missing-triples theorem for compact H9/14 nerves",
            "every compact Coxeter missing face is Lanner of order at most five",
            "strict generalized lower-bound inequalities in every vertex link",
            "Dehn--Sommerville face-vector schedules",
            (
                "facets are complementary to the inclusion-minimal "
                "transversals of the complete minimal-nonface clutter"
            ),
            "every ridge of a simplicial sphere has degree exactly two",
            "Cauchy--Schwarz for the degree sequence of a finite graph",
        ],
        "schedule_table": schedule_table,
        "strict_link_g4_exclusion": {
            "required_sum": 14,
            "excluded": {
                "5": summed_vertex_link_g4(schedules[5]["f_vector"]),
                "6": summed_vertex_link_g4(schedules[6]["f_vector"]),
            },
        },
        "a0_analytic_audit": a0.to_record(),
        "backend_corroboration": _backend_corroboration(),
        "frontier_effect": {
            "record16_schedule_count_before": 5,
            "record16_schedule_count_after": 4,
            "graph_shadow_partition_changed": False,
            "current_global_partition": {
                "independently_certified_impossible": 121,
                "certified_explicit_ordinary_polytopal": 7,
                "theorem_undecided": 14,
                "total_graph_shadows": 142,
            },
        },
        "nonclaims": [
            "Record 16 as a whole is not excluded; branches a=1,2,3,4 remain.",
            "The archived bounded Z3 timeout contributes no logical step.",
            "No Coxeter labels or Lorentzian Gram matrices are classified.",
            "This is not a general compact H9/14 non-existence theorem.",
        ],
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_record16_a0_analytic_exclusion":
        raise ValueError("wrong record-16 a=0 certificate kind")
    if (
        manifest.get("record_index") != TARGET_RECORD_INDEX
        or manifest.get("graph6") != TARGET_GRAPH6
    ):
        raise ValueError("record-16 a=0 certificate target changed")
    if manifest.get("theorem_status") != "certified_branch_unsat":
        raise ValueError("record-16 a=0 branch status changed")
    if _file_sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise ValueError("record-16 support source is not frozen")
    if _file_sha256(BACKEND_SOURCE) != EXPECTED_BACKEND_SOURCE_SHA256:
        raise ValueError("record-16 backend source is not frozen")
    proof_digest = _file_sha256(PROOF_SOURCE)
    if (
        EXPECTED_PROOF_SOURCE_SHA256 != "TO_BE_FROZEN"
        and proof_digest != EXPECTED_PROOF_SOURCE_SHA256
    ):
        raise ValueError("record-16 proof source is not frozen")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-16 a=0 certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-16 a=0 certificate digest mismatch")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-16 a=0 certificate is not frozen")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_branch_unsat",
        "excluded_missing_four_counts": [0, 5, 6],
        "remaining_missing_four_counts": [1, 2, 3, 4],
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
