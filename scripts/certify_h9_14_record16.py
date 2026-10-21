#!/usr/bin/env python3
"""Build and verify the complete solver-independent record-16 exclusion."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record16 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
)
from coxeter_polytopes.combinatorics.h9_14_record16_positive import (
    audit_record16_positive_branches,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.nerve_support import (
    verify_support_manifest,
)
from scripts.certify_h9_14_record16_a0 import (
    OUTPUT as A0_CERTIFICATE,
    verify as verify_a0,
)
from scripts.independent_h9_14_record16_audit import (
    EXPECTED_SUMMARY_SHA256,
    audit as independent_audit,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
BACKEND_SOURCE = ROOT / "data/enumerations/h9_14_facet_completion.json"
X_ARCHIVE = ROOT / "data/enumerations/h9_14_record16_x_graphs.json"
POSITIVE_SOURCE = (
    ROOT
    / "src/coxeter_polytopes/combinatorics/h9_14_record16_positive.py"
)
A0_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record16.py"
LOCAL_LINK_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/local_link.py"
MISSING_EDGES_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/missing_edges.py"
)
RECORD92_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record92.py"
NERVE_SUPPORT_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/nerve_support.py"
)
A0_CERTIFIER_SOURCE = ROOT / "scripts/certify_h9_14_record16_a0.py"
GENERATOR_SOURCE = ROOT / "scripts/enumerate_h9_14_record16_x_graphs.py"
INDEPENDENT_SOURCE = ROOT / "scripts/independent_h9_14_record16_audit.py"
OUTPUT = ROOT / "data/enumerations/h9_14_record16_certificate.json"

EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_BACKEND_SOURCE_SHA256 = (
    "2d329439e5baeffb574932ff7bbab2dd78daf706b2917a0674bdc9f01a2f5123"
)
EXPECTED_X_ARCHIVE_SHA256 = (
    "1d6915efe97c2ec5fc81e0022d415c1073e17bbed0e012038a1171132630a1ce"
)
EXPECTED_A0_SOURCE_SHA256 = (
    "fd2d02c223d683d9e72ad2e6e72edd9d928c6b0a4fe2525dc254b3b5bab8e564"
)
EXPECTED_POSITIVE_SOURCE_SHA256 = (
    "6fd9d9842db6b5e907a9e5f65a4f5c04dd9b0485b4d71a9682167075993abe19"
)
EXPECTED_LOCAL_LINK_SOURCE_SHA256 = (
    "8f8169edab8b87aff682785667c2a8eab5c52b92a7c68cef130aeec4d3ae807d"
)
EXPECTED_MISSING_EDGES_SOURCE_SHA256 = (
    "231e85fa98e27e864f997cd4ebe59136707f6aa8dd2d8136b80cff543dfbfbb8"
)
EXPECTED_RECORD92_SOURCE_SHA256 = (
    "6b165e2a9e3ed120949c10784e42786364a7197014bc57da6b2b0a1cb30e1e29"
)
EXPECTED_NERVE_SUPPORT_SOURCE_SHA256 = (
    "aad8fb25c2f9a56e45121fa1d1cbbdedfdff8ba1f67e379be74eed73940ad2a2"
)
EXPECTED_A0_CERTIFIER_SOURCE_SHA256 = (
    "9b032b4afce9a76bf67bd06ad52374b2fcaf8c28dcc74f1333bae7c23c92f9c2"
)
EXPECTED_GENERATOR_SOURCE_SHA256 = (
    "d9b6457e483e24e0fe78ab89e7f10a5d8f802d9536dcfd5731d4bb83f8f7171f"
)
EXPECTED_INDEPENDENT_SOURCE_SHA256 = (
    "2961fe9223689d744e2a2ab8a55e92459bb3520f21aef5859a0ac0442c1f3139"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "ac7f386a36ffe8a0ddbeab74c0fdd60913eb891365735588ee7a111fc8b7c8d5"
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


def _source_record(source: Mapping[str, object]) -> Mapping[str, object]:
    records = source.get("records")
    if not isinstance(records, list):
        raise ValueError("record-16 support source lacks records")
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
    matching = [
        record
        for record in backend["records"]
        if record.get("index") == TARGET_RECORD_INDEX
    ]
    if len(matching) != 1 or matching[0]["solver"]["status"] != "unknown":
        raise ValueError("record-16 bounded backend disposition changed")
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


def _proof_sources() -> list[dict[str, str]]:
    return [
        {
            "path": str(path.relative_to(ROOT)),
            "sha256": _file_sha256(path),
        }
        for path in (
            A0_SOURCE,
            POSITIVE_SOURCE,
            LOCAL_LINK_SOURCE,
            MISSING_EDGES_SOURCE,
            RECORD92_SOURCE,
            NERVE_SUPPORT_SOURCE,
            A0_CERTIFIER_SOURCE,
            GENERATOR_SOURCE,
            INDEPENDENT_SOURCE,
        )
    ]


def build() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    schedules = record["dehn_sommerville_schedule"]

    a0_manifest = json.loads(A0_CERTIFICATE.read_text(encoding="utf-8"))
    a0_summary = verify_a0(a0_manifest)
    positive = audit_record16_positive_branches(schedules, X_ARCHIVE)
    independent = independent_audit()
    if independent["summary_sha256"] != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("record-16 hostile replay digest changed")
    x_document = json.loads(X_ARCHIVE.read_text(encoding="utf-8"))
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
        "kind": "h9_14_record16_solver_independent_exclusion",
        "claim": (
            "H9/14 second-layer survivor record 16 has no simplicial-sphere "
            "completion satisfying the compact Coxeter missing-face "
            "reductions, and hence cannot be the nerve of a compact "
            "hyperbolic Coxeter 9-polytope."
        ),
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_unsat",
        "proof_backend": (
            "blocker duality, complete ridge incidence, an analytic graph "
            "degree bound, exact nauty triangle-free graph streams through "
            "order eight with independent property replay, exhaustive A5 "
            "and full-blocker dualization, and a separate hostile replay; "
            "no SAT/SMT result is trusted"
        ),
        "source": {
            "path": "data/enumerations/h9_14_nerve_support.json",
            "sha256": _file_sha256(SOURCE),
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
        },
        "proof_sources": _proof_sources(),
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
            "Cauchy--Schwarz for finite graph degree sequences",
            "nauty geng exhaustive generation of unlabeled triangle-free graphs",
        ],
        "schedule_table": schedule_table,
        "strict_link_g4_exclusion": {
            "required_sum": 14,
            "excluded_missing_four_counts": [5, 6],
            "observed_sums": {
                "5": summed_vertex_link_g4(schedules[5]["f_vector"]),
                "6": summed_vertex_link_g4(schedules[6]["f_vector"]),
            },
        },
        "a0_subcertificate": {
            "path": str(A0_CERTIFICATE.relative_to(ROOT)),
            "file_sha256": _file_sha256(A0_CERTIFICATE),
            "summary": a0_summary,
        },
        "positive_branch_audit": positive.to_record(),
        "x_graph_archive": {
            "path": str(X_ARCHIVE.relative_to(ROOT)),
            "file_sha256": _file_sha256(X_ARCHIVE),
            "certificate_sha256": x_document["certificate_sha256"],
            "generator": x_document["generator"],
            "streams": {
                order: {
                    "command": stream["command"],
                    "graph_count": stream["graph_count"],
                    "raw_stream_sha256": stream["raw_stream_sha256"],
                }
                for order, stream in x_document["streams"].items()
            },
        },
        "independent_hostile_replay": {
            "source": str(INDEPENDENT_SOURCE.relative_to(ROOT)),
            "source_sha256": _file_sha256(INDEPENDENT_SOURCE),
            "summary_sha256": independent["summary_sha256"],
            "blocker_and_ridges": independent["blocker_and_ridges"],
            "z_cut": independent["z_cut"],
            "schedules": independent["schedules"],
            "x_archive": independent["x_archive"],
            "final_case": independent["final_case"],
        },
        "backend_corroboration": _backend_corroboration(),
        "updated_theorem_partition": {
            "independently_certified_impossible": 122,
            "certified_explicit_ordinary_polytopal": 7,
            "theorem_undecided": 13,
            "total_graph_shadows": 142,
        },
        "nonclaims": [
            "This excludes one graph shadow, not all compact H9/14 nerves.",
            "The archived bounded Z3 timeout contributes no logical step.",
            "No Coxeter labels or Lorentzian Gram matrices are classified.",
            "This is not a general compact non-existence theorem in dimension 9.",
        ],
    }
    # Normalize tuples returned by the independent replay before comparing a
    # rebuilt certificate with its JSON round trip.
    manifest = json.loads(json.dumps(manifest, sort_keys=True))
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_record16_solver_independent_exclusion":
        raise ValueError("wrong record-16 certificate kind")
    if (
        manifest.get("record_index") != TARGET_RECORD_INDEX
        or manifest.get("graph6") != TARGET_GRAPH6
    ):
        raise ValueError("record-16 certificate target changed")
    if manifest.get("theorem_status") != "certified_unsat":
        raise ValueError("record-16 theorem status changed")
    expected_files = {
        SOURCE: EXPECTED_SOURCE_SHA256,
        BACKEND_SOURCE: EXPECTED_BACKEND_SOURCE_SHA256,
        X_ARCHIVE: EXPECTED_X_ARCHIVE_SHA256,
        A0_SOURCE: EXPECTED_A0_SOURCE_SHA256,
        POSITIVE_SOURCE: EXPECTED_POSITIVE_SOURCE_SHA256,
        LOCAL_LINK_SOURCE: EXPECTED_LOCAL_LINK_SOURCE_SHA256,
        MISSING_EDGES_SOURCE: EXPECTED_MISSING_EDGES_SOURCE_SHA256,
        RECORD92_SOURCE: EXPECTED_RECORD92_SOURCE_SHA256,
        NERVE_SUPPORT_SOURCE: EXPECTED_NERVE_SUPPORT_SOURCE_SHA256,
        A0_CERTIFIER_SOURCE: EXPECTED_A0_CERTIFIER_SOURCE_SHA256,
        GENERATOR_SOURCE: EXPECTED_GENERATOR_SOURCE_SHA256,
        INDEPENDENT_SOURCE: EXPECTED_INDEPENDENT_SOURCE_SHA256,
    }
    for path, expected in expected_files.items():
        observed = _file_sha256(path)
        if expected != "TO_BE_FROZEN" and observed != expected:
            raise ValueError(f"record-16 frozen source changed: {path}")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-16 certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-16 certificate digest mismatch")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-16 certificate is not frozen")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_unsat",
        "certificate_sha256": digest,
        "updated_theorem_partition": manifest["updated_theorem_partition"],
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
