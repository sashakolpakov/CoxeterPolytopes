#!/usr/bin/env python3
"""Freeze the corrected exact 371-orbit record-12 abstract census."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record12_residual import (
    EXPECTED_COLOURED_BRANCH_ORBIT_COUNT,
    EXPECTED_FULL_BLOCKER_PASS_COUNT,
    EXPECTED_FULL_ORBIT_COUNT,
    EXPECTED_LOCAL_PAIR_FEASIBLE_STATE_COUNT,
    EXPECTED_NORMALIZED_CUT_STATE_COUNT,
    EXPECTED_PROFILE_FACET_TYPE_PAIR_COUNT,
    EXPECTED_RELATIVE_EMBEDDING_COUNT,
    EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT,
    EXPECTED_SCHEDULE_ORBIT_COUNTS,
    EXPECTED_X_GRAPH_COUNT,
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    verify_archive,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import summed_vertex_link_g4
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
X_ARCHIVE = ROOT / "data/enumerations/h9_14_record12_x_graphs.json"
RESIDUAL_ARCHIVE = ROOT / "data/enumerations/h9_14_record12_residual.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record12_residual_certificate.json"

PROOF_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record12_residual.py"
HOMOLOGY_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/simplicial_homology.py"
ARCHIVE_BUILDER_SOURCE = ROOT / "scripts/build_h9_14_record12_residual_archive.py"
X_GENERATOR_SOURCE = ROOT / "scripts/enumerate_h9_14_record12_x_graphs.py"
BRANCH_GENERATOR_SOURCE = ROOT / "scripts/explore_h9_14_record12_branch.cpp"
COLOURED_QUOTIENT_SOURCE = ROOT / "scripts/explore_h9_14_record12_unfiltered_orbits.py"
PAIR_GENERATOR_SOURCE = ROOT / "scripts/explore_h9_14_record12_safe_pairs.cpp"
ORBIT_GENERATOR_SOURCE = ROOT / "scripts/explore_h9_14_record12_pair_orbits.cpp"
HOSTILE_BRANCH_SOURCE = ROOT / "scripts/independent_h9_14_record12_371_branches.py"
HOSTILE_PAIR_SOURCE = ROOT / "scripts/independent_h9_14_record12_371_pairs.cpp"
HOSTILE_QUOTIENT_SOURCE = ROOT / "scripts/independent_h9_14_record12_371_quotient.cpp"
HOSTILE_TOPOLOGY_SOURCE = ROOT / "scripts/independent_h9_14_record12_371_orbits.py"
WITHDRAWAL_REPORT = ROOT / "results/h9_14_record12_residual_hostile_audit.md"
HOSTILE_REPORT = ROOT / "results/h9_14_record12_corrected_371_hostile_audit.md"

EXPECTED_FILES = {
    SOURCE: "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d",
    X_ARCHIVE: "a055814261d48c85532baf38e8d51af7fca4970f390892baabe66cccc61cb80e",
    RESIDUAL_ARCHIVE: "faca367738bba58afc0ec78a43e9d4aa4b02721802c5d9e77b7d155c8e7faa05",
    PROOF_SOURCE: "8d71264321ede7378193e6be8c086700b8235b1d0c1216fdc9ed7ef495db86f4",
    HOMOLOGY_SOURCE: "547f51b58c1146aac80c6e08910754bf93ace5fa808dae5a342da7e042958110",
    ARCHIVE_BUILDER_SOURCE: "07d5a2fecb1ddf78d4881c36a765ededc05845ef364b0ae4e55f9789bcb79fd7",
    X_GENERATOR_SOURCE: "ac4d08dfe557e13d130d791b4b4bdbfec8f3d9fcdeb5f86810ca104635f6df0c",
    BRANCH_GENERATOR_SOURCE: "93fae6156dd8e4ea2febfbfe614bfaa37a55cae76b5888ca268b03b523bd5032",
    COLOURED_QUOTIENT_SOURCE: "f67cfcf614140df7531ca7df1e0a475794290e73f7ad541ca17767a407f7a06f",
    PAIR_GENERATOR_SOURCE: "2d1eff38d2a431fad3d5119b99ac458b3cf4436dc7deb430de0b3013d5965c7b",
    ORBIT_GENERATOR_SOURCE: "1b351273bf34eb76e2bf2357ca70175596e36a874a31449221c38b13f56d2cd7",
    HOSTILE_BRANCH_SOURCE: "edd501dabf990f8da8aac20580e57f928d8b44507b64f413623a0427147a8ff5",
    HOSTILE_PAIR_SOURCE: "f60d924e95f364377187165595b9d343d805144e382b48fc3ccaff3ff232976b",
    HOSTILE_QUOTIENT_SOURCE: "5da00e6ecff30164fde90afece4ac2d5d7e9de0728f6da660172f92ee7cdd197",
    HOSTILE_TOPOLOGY_SOURCE: "993b0cfc9634eae2c2422f409843984cc5b0b895bfe04a70aa8473424a11e5b4",
    WITHDRAWAL_REPORT: "2d0b484e26622a970f9f8b7823d4d00e7055f582e9f16ce4975ed8119274b2fd",
    HOSTILE_REPORT: "20284682d8956d27e5c36a45d75b7e2fba2b92c0804935dad5edf59455a51970",
}

EXPECTED_BRANCH_SUMMARY_SHA256 = (
    "de2ce0116ee97e2b9f958f6e3b616803744b0d96817f221db341931b2dfb3e3c"
)
EXPECTED_PAIR_SET_SHA256 = (
    "a189faa3d1faeb1f4db21630d2fbb816f8d3d4f1421d04eb251112d37945f336"
)
EXPECTED_TOPOLOGY_SUMMARY_SHA256 = (
    "6a58b02d39324916d24de5d43a8b1091393f83cadc7963de8581023266c66c84"
)
EXPECTED_TOPOLOGY_OUTPUT_SHA256 = (
    "3119558a2e170ef2b2419b7049409e800d4ad3177fa2538db65ba8fad95145fc"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "c759980c4a18cbc9c5c5cdf892b7fe68a142bec472f7b8be5f7e9abf9b6c9a98"
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _record(source: Mapping[str, object]) -> Mapping[str, object]:
    survivors = [
        row
        for row in source["records"]
        if isinstance(row, Mapping) and row.get("survives_second_layer")
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("record-12 support row changed")
    return record


def _proof_sources() -> list[dict[str, str]]:
    paths = (
        PROOF_SOURCE,
        HOMOLOGY_SOURCE,
        ARCHIVE_BUILDER_SOURCE,
        X_GENERATOR_SOURCE,
        BRANCH_GENERATOR_SOURCE,
        COLOURED_QUOTIENT_SOURCE,
        PAIR_GENERATOR_SOURCE,
        ORBIT_GENERATOR_SOURCE,
        HOSTILE_BRANCH_SOURCE,
        HOSTILE_PAIR_SOURCE,
        HOSTILE_QUOTIENT_SOURCE,
        HOSTILE_TOPOLOGY_SOURCE,
    )
    return [
        {"path": str(path.relative_to(ROOT)), "sha256": _file_sha256(path)}
        for path in paths
    ]


def build() -> dict[str, object]:
    support_document = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(support_document)
    support = _record(support_document)
    x_document = json.loads(X_ARCHIVE.read_text(encoding="utf-8"))
    residual_document = json.loads(RESIDUAL_ARCHIVE.read_text(encoding="utf-8"))
    residual_summary = verify_archive(residual_document)

    schedule_table = [
        {
            "missing_four_count": row["missing_four_count"],
            "missing_five_count": row["missing_five_count_offset"],
            "facet_count": row["f_vector"][-1],
            "summed_vertex_link_g4": summed_vertex_link_g4(row["f_vector"]),
            "residual_orbit_count": EXPECTED_SCHEDULE_ORBIT_COUNTS.get(
                int(row["missing_four_count"]), 0
            ),
        }
        for row in support["dehn_sommerville_schedule"]
    ]

    manifest: dict[str, object] = {
        "schema_version": "2.0.0",
        "kind": "h9_14_record12_exact_371_orbit_whole_shadow_certificate",
        "claim": (
            "Under the published compact-H9/14 missing-face reductions, the "
            "entire record-12 graph shadow has exactly 371 abstract residual "
            "incidence-complex orbits, with a=0,...,4 histogram "
            "48,95,144,63,21. Every representative passes the archived "
            "incidence, local g/Macaulay, Euler-link, and mod-two homology-link "
            "audits. PL, ordinary-polytopal, and Coxeter realizability are not "
            "asserted."
        ),
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_exact_371_orbit_abstract_residual",
        "proof_backend": (
            "blocker duality; all local zero-or-two ridge equations; exhaustive "
            "nauty order-ten triangle-free graph stream; all normalized cuts "
            "and y supports; only the proved special-meeting branch-blocker "
            "filter; exact coloured branch quotient; both direct and "
            "complementary q-cut orientations; complete coupled blocker before "
            "filtering isolate-only rows; full D8 x S10 quotient; direct face "
            "tables and all-face GF(2) link homology; no SAT/SMT result"
        ),
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": _file_sha256(SOURCE),
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
        },
        "x_graph_archive": {
            "path": str(X_ARCHIVE.relative_to(ROOT)),
            "file_sha256": _file_sha256(X_ARCHIVE),
            "certificate_sha256": x_document["certificate_sha256"],
            "raw_stream_sha256": x_document["raw_stream_sha256"],
            "graph_count": x_document["graph_count"],
        },
        "residual_archive": {
            "path": str(RESIDUAL_ARCHIVE.relative_to(ROOT)),
            "file_sha256": _file_sha256(RESIDUAL_ARCHIVE),
            **residual_summary,
        },
        "proof_sources": _proof_sources(),
        "proof_dependencies": [
            "Burcroff's no-missing-triples theorem for compact H9/14 nerves",
            "compact Coxeter missing faces are Lanner diagrams of order at most five",
            "strict generalized lower-bound and Macaulay inequalities in all links",
            "Dehn--Sommerville schedules for simplicial 9-polytopes",
            "facet complements are the blocker of the complete minimal-nonface clutter",
            "every ridge of a simplicial sphere has degree two",
            "nauty geng supplies one representative of each unlabeled triangle-free graph",
        ],
        "safe_special_filter_lemma": residual_document[
            "safe_special_filter_lemma"
        ],
        "generation_streams": residual_document["generation_streams"],
        "census": residual_document["census"],
        "cut_normalization_audit": residual_document[
            "cut_normalization_audit"
        ],
        "pairing_audit": residual_document["pairing_audit"],
        "automorphism_quotient": residual_document["automorphism_quotient"],
        "schedule_table": schedule_table,
        "schedule_orbit_counts": residual_document["schedule_orbit_counts"],
        "excluded_schedule_branches": {
            "a=5,6,7": (
                "no complete coupled blocker survives after exact branch "
                "pairing, both q-cut orientations, and the global quotient"
            ),
            "a=8,9": (
                "strict summed vertex-link g4 requires at least 14, while the "
                "two sums are 12 and 6"
            ),
        },
        "independent_hostile_replay": {
            "report": {
                "path": str(HOSTILE_REPORT.relative_to(ROOT)),
                "sha256": _file_sha256(HOSTILE_REPORT),
            },
            "branch_summary_sha256": EXPECTED_BRANCH_SUMMARY_SHA256,
            "semantic_pair_set_sha256": EXPECTED_PAIR_SET_SHA256,
            "topology_summary_sha256": EXPECTED_TOPOLOGY_SUMMARY_SHA256,
            "topology_output_file_sha256": EXPECTED_TOPOLOGY_OUTPUT_SHA256,
            "x_graphs": EXPECTED_X_GRAPH_COUNT,
            "locally_feasible_rows": EXPECTED_LOCAL_PAIR_FEASIBLE_STATE_COUNT,
            "safe_special_rows": EXPECTED_SAFE_SPECIAL_BRANCH_STATE_COUNT,
            "coloured_branch_orbits": EXPECTED_COLOURED_BRANCH_ORBIT_COUNT,
            "compatible_ordered_type_pairs": EXPECTED_PROFILE_FACET_TYPE_PAIR_COUNT,
            "relative_embeddings": EXPECTED_RELATIVE_EMBEDDING_COUNT,
            "complete_blocker_passes": EXPECTED_FULL_BLOCKER_PASS_COUNT,
            "full_symmetry_orbits": EXPECTED_FULL_ORBIT_COUNT,
            "all_full_blockers_verified": True,
            "all_ridge_degrees_two": True,
            "all_dual_graphs_connected": True,
            "all_strict_local_g_and_macaulay_audits": True,
            "all_face_links_euler_spheres": True,
            "all_face_links_mod2_homology_spheres": True,
        },
        "correction_history": residual_document["correction_history"],
        "relation_to_known_example": residual_document[
            "relation_to_known_example"
        ],
        "nonclaims": residual_document["nonclaims"],
    }
    manifest = json.loads(json.dumps(manifest, sort_keys=True))
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if (
        manifest.get("kind")
        != "h9_14_record12_exact_371_orbit_whole_shadow_certificate"
    ):
        raise ValueError("wrong record-12 residual certificate kind")
    if (
        manifest.get("record_index") != TARGET_RECORD_INDEX
        or manifest.get("graph6") != TARGET_GRAPH6
        or manifest.get("theorem_status")
        != "certified_exact_371_orbit_abstract_residual"
    ):
        raise ValueError("record-12 residual certificate target changed")
    for path, expected in EXPECTED_FILES.items():
        if _file_sha256(path) != expected:
            raise ValueError(f"record-12 residual frozen source changed: {path}")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-12 residual certificate digest mismatch")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-12 residual certificate is not frozen")
    if manifest != build():
        raise ValueError("record-12 residual certificate did not reproduce")
    census = manifest["census"]
    if not isinstance(census, Mapping):
        raise ValueError("record-12 residual census is malformed")
    if (
        census.get("x_graph_count") != EXPECTED_X_GRAPH_COUNT
        or census.get("normalized_cut_states")
        != EXPECTED_NORMALIZED_CUT_STATE_COUNT
        or census.get("full_graph_shadow_orbits") != EXPECTED_FULL_ORBIT_COUNT
    ):
        raise ValueError("record-12 residual headline census changed")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": manifest["theorem_status"],
        "coloured_branch_orbits": census["coloured_branch_orbits"],
        "complete_blocker_passes": census[
            "full_blocker_passing_relative_embeddings"
        ],
        "residual_orbits": census["full_graph_shadow_orbits"],
        "schedule_orbit_counts": manifest["schedule_orbit_counts"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args(argv)
    document = (
        json.loads(OUTPUT.read_text(encoding="utf-8"))
        if args.verify_only
        else build()
    )
    if not args.verify_only:
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
