#!/usr/bin/env python3
"""Freeze the exact incidence no-go for H9/14 record 6."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record6 import (
    EXPECTED_NUMERICAL_BRANCHES,
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_link_census,
    audit_ten_branches,
    derive_numerical_branches,
    e_link_orbit,
    h_link_orbit,
    verify_target_graph,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
PROOF_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record6.py"
SUPPORT_SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record6_incidence_no_go.json"
EXPECTED_PROOF_SOURCE_SHA256 = "e9d10285b73e92ca323a65a69126b9766cbda2945a67e956f6c84364ba29a367"
EXPECTED_SUPPORT_SOURCE_SHA256 = "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
EXPECTED_GENG_STDOUT_SHA256 = "162cd507b93e74307ed9e88052da7bfb6f68e6cfeae6ccd2254c949c84ce0214"
EXPECTED_CERTIFICATE_SHA256 = "cc0c3186b16cafc3c3a9572ea15f6d24880ad7bb3b18da40cf8c635d064f78ff"


def file_hash(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def certificate_digest(manifest: Mapping[str, object]) -> str:
    payload = dict(manifest)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _geng_words() -> tuple[tuple[str, ...], str, str]:
    executable = shutil.which("geng")
    if executable is None:
        raise RuntimeError("nauty geng is required for the frozen link census")
    completed = subprocess.run(
        [executable, "-q", "-t", "10"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    digest = sha256(completed.stdout).hexdigest()
    if digest != EXPECTED_GENG_STDOUT_SHA256:
        raise ValueError("record-6 triangle-free graph transcript changed")
    words = tuple(line.decode("ascii") for line in completed.stdout.splitlines())
    return words, digest, executable


def _support_record() -> tuple[Mapping[str, object], str]:
    digest = file_hash(SUPPORT_SOURCE)
    if digest != EXPECTED_SUPPORT_SOURCE_SHA256:
        raise ValueError("record-6 support source is not frozen")
    document = json.loads(SUPPORT_SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(document)
    record = [row for row in document["records"] if row["survives_second_layer"]][
        TARGET_RECORD_INDEX
    ]
    if record["graph6"] != TARGET_GRAPH6 or not verify_target_graph():
        raise ValueError("record-6 survivor index or graph changed")
    return record, digest


def build() -> dict[str, object]:
    proof_hash = file_hash(PROOF_SOURCE)
    if proof_hash != EXPECTED_PROOF_SOURCE_SHA256:
        raise ValueError("record-6 proof source is not frozen")
    support, support_hash = _support_record()
    words, transcript_hash, executable = _geng_words()
    link_census = audit_link_census(words)
    if not link_census.verified:
        raise AssertionError("record-6 special-link census changed")
    numerical = derive_numerical_branches()
    if numerical != EXPECTED_NUMERICAL_BRANCHES:
        raise AssertionError("record-6 numerical branch derivation changed")
    branches = audit_ten_branches()
    if len(branches) != 10 or not all(branch.excluded for branch in branches):
        raise AssertionError("record-6 incidence branch survived")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record6_exact_incidence_no_go",
        "claim": (
            "No simplicial 9-polytope boundary satisfying the frozen Coxeter-"
            "nerve incidence necessities has missing-edge graph K3+11I; in "
            "particular H9/14 record 6 cannot be a compact hyperbolic Coxeter "
            "polytope nerve."
        ),
        "theorem_status": "certified_incidence_impossible",
        "target": {
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
            "graph6": TARGET_GRAPH6,
            "missing_edges": [[0, 12], [0, 13], [12, 13]],
            "support_core_size": len(support["ridge_support_core_masks"]),
        },
        "sources": {
            "proof": {
                "path": "src/coxeter_polytopes/combinatorics/h9_14_record6.py",
                "sha256": proof_hash,
            },
            "ridge_support": {
                "path": "data/enumerations/h9_14_nerve_support.json",
                "sha256": support_hash,
            },
        },
        "generator": {
            "command": "geng -q -t 10",
            "resolved_executable": executable,
            "graph6_stdout_sha256": transcript_hash,
        },
        "link_census": link_census.to_record(),
        "source_link_orbits": {
            "E_3_plus_4_plus_4_labelled_count": len(e_link_orbit()),
            "H_3_plus_3_plus_5_labelled_count": len(h_link_orbit()),
            "completeness_check": (
                "Every fixed-base E survivor belongs to the full labelled E "
                "partition orbit, and every fixed-base H survivor belongs to "
                "the full labelled H partition orbit. F/G/A occur only as "
                "derived targets and are accepted by the complete link predicate."
            ),
        },
        "numerical_branches": [
            {
                "link_signatures": [list(signature) for signature in signatures],
                "isolate_missing_four_count": a0,
                "z_size": z_size,
            }
            for signatures, a0, z_size in numerical
        ],
        "branch_audits": [branch.to_record() for branch in branches],
        "proof_lemmas": [
            (
                "Every special-vertex link is a two-graph delta(X); its 0/2 "
                "ridge rule is equivalent to X being triangle-free and induced-"
                "2K2-free. The exact census leaves 42 fixed-base types."
            ),
            (
                "For link signature (y,q,r), Dehn--Sommerville forces z = "
                "150-2*(sum(q)+a0)-sum(y). Exhausting all admissible signatures "
                "derives exactly ten branches."
            ),
            (
                "After one common isolate relabelling fixes an E or H link, all "
                "5,775 or 4,620 labelled copies of the second link and every "
                "labelled Z support are enumerated. No branch satisfies both "
                "the no-missing-triples theorem and its forced a0 value."
            ),
        ],
        "nonclaims": [
            "This supplies no new hyperbolic Coxeter example.",
            "The obstruction is incidence-theoretic; no Gram-matrix search is needed.",
        ],
    }
    manifest["certificate_sha256"] = certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-6 certificate did not reproduce")
    digest = certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-6 certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("record-6 certificate is not frozen")
    return {
        "fixed_base_link_type_count": rebuilt["link_census"]["audited_fixed_base_link_types"],
        "numerical_branch_count": len(rebuilt["numerical_branches"]),
        "excluded_branch_count": sum(row["excluded"] for row in rebuilt["branch_audits"]),
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    del argv
    manifest = build()
    OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), **verify(manifest)}, indent=2))


if __name__ == "__main__":
    main()
