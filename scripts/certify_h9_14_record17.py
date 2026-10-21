#!/usr/bin/env python3
"""Build and verify the exact whole-shadow residual certificate for record 17."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_record17 import (
    TARGET_GRAPH6,
    TARGET_RECORD_INDEX,
    audit_record17,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    summed_vertex_link_g4,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest
from scripts.independent_h9_14_record17_audit import (
    EXPECTED_SUMMARY_SHA256,
    audit as independent_audit,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
X_ARCHIVE = ROOT / "data/enumerations/h9_14_record17_x_graphs.json"
PROOF_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record17.py"
FACET_AUDIT_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/facet_completion.py"
HOMOLOGY_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/simplicial_homology.py"
LOCAL_LINK_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/local_link.py"
MISSING_EDGES_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/missing_edges.py"
RECORD92_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record92.py"
NERVE_SUPPORT_SOURCE = ROOT / "src/coxeter_polytopes/combinatorics/nerve_support.py"
GENERATOR_SOURCE = ROOT / "scripts/enumerate_h9_14_record17_x_graphs.py"
INDEPENDENT_SOURCE = ROOT / "scripts/independent_h9_14_record17_audit.py"
OUTPUT = ROOT / "data/enumerations/h9_14_record17_certificate.json"

EXPECTED_SOURCE_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
EXPECTED_X_ARCHIVE_SHA256 = (
    "60350a9c2e25043845810a58ef61d2c0fe3ea51bb5cf8c0c93e4086635c82395"
)
EXPECTED_PROOF_SOURCE_SHA256 = (
    "e01895b4fc4c5767154516c0b189754e72bd945bd145b40436ec8357382fbede"
)
EXPECTED_FACET_AUDIT_SOURCE_SHA256 = (
    "5f2d078cfe400b39f3e194352f7a0c6d3bd8cb6bb8e7cee37be159e8766d9400"
)
EXPECTED_HOMOLOGY_SOURCE_SHA256 = (
    "547f51b58c1146aac80c6e08910754bf93ace5fa808dae5a342da7e042958110"
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
EXPECTED_GENERATOR_SOURCE_SHA256 = (
    "1bc78a18dd6a653f962c29838284b743d1b25e5654f7a20c6c6198cae0c4aea9"
)
EXPECTED_INDEPENDENT_SOURCE_SHA256 = (
    "1d30692f545c6044080eb5da218de76cf13f284837db98d0c53ae3afcd53f093"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "b98e58b390cbf8d9b40e08da16da29f78a2f04805bb0b13d7f2efb623b9afb51"
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
        raise ValueError("record-17 support source lacks records")
    survivors = [
        record
        for record in records
        if isinstance(record, Mapping) and record.get("survives_second_layer")
    ]
    record = survivors[TARGET_RECORD_INDEX]
    if record.get("graph6") != TARGET_GRAPH6:
        raise ValueError("second-layer survivor record 17 changed")
    return record


def _proof_sources() -> list[dict[str, str]]:
    return [
        {"path": str(path.relative_to(ROOT)), "sha256": _file_sha256(path)}
        for path in (
            PROOF_SOURCE,
            FACET_AUDIT_SOURCE,
            HOMOLOGY_SOURCE,
            LOCAL_LINK_SOURCE,
            MISSING_EDGES_SOURCE,
            RECORD92_SOURCE,
            NERVE_SUPPORT_SOURCE,
            GENERATOR_SOURCE,
            INDEPENDENT_SOURCE,
        )
    ]


def _topology_audits(
    record: Mapping[str, object], proof_audit
) -> list[dict[str, object]]:
    support = record["ridge_support_core_masks"]
    schedules = record["dehn_sommerville_schedule"]
    result = []
    for orbit in proof_audit.enumeration.residual_orbits:
        representative = orbit.representative
        audit = audit_global_facet_selection(
            TARGET_GRAPH6,
            support,
            schedules,
            representative.selected_facet_complement_masks,
        )
        if not audit.passes_connected_homology_audit:
            raise AssertionError("a record-17 residual failed the full link audit")
        result.append(
            {
                "missing_four_count": representative.missing_four_count,
                "x_graph6": representative.x_graph6,
                "p_local_vertices": [
                    vertex
                    for vertex in range(9)
                    if representative.p_local_mask & (1 << vertex)
                ],
                "representation_count": orbit.representation_count,
                "facet_count": audit.facet_count,
                "face_vector": list(audit.face_vector),
                "ridge_degree_histogram": {
                    str(degree): count
                    for degree, count in audit.ridge_degree_histogram.items()
                },
                "dual_component_count": audit.dual_component_count,
                "minimal_nonface_histogram": {
                    str(size): len(values)
                    for size, values in audit.minimal_nonfaces_by_size.items()
                },
                "matching_schedule_index": audit.matching_schedule_index,
                "passes_selection_model": audit.passes_selection_model,
                "all_links_have_spherical_euler_characteristic": (
                    audit.all_links_have_spherical_euler_characteristic
                ),
                "all_links_are_mod2_homology_spheres": (
                    audit.all_links_are_mod2_homology_spheres
                ),
                "dual_graph_connected": audit.dual_graph_connected,
                "passes_connected_homology_audit": (
                    audit.passes_connected_homology_audit
                ),
            }
        )
    return result


def build() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    record = _source_record(source)
    schedules = record["dehn_sommerville_schedule"]
    proof = audit_record17(schedules, X_ARCHIVE)
    hostile = independent_audit()
    if hostile["summary_sha256"] != EXPECTED_SUMMARY_SHA256:
        if EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN":
            raise AssertionError("record-17 hostile replay digest changed")
    topology = _topology_audits(record, proof)
    schedule_table = [
        {
            "missing_four_count": row["missing_four_count"],
            "missing_five_count": row["missing_five_count_offset"],
            "facet_count": row["f_vector"][-1],
            "summed_vertex_link_g4": summed_vertex_link_g4(row["f_vector"]),
        }
        for row in schedules
    ]
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record17_exact_whole_shadow_residual",
        "claim": (
            "Under the compact-H9/14 missing-face reductions, second-layer "
            "record 17 has exactly eight residual incidence-complex orbits: "
            "one in a=0, four in a=1, and three in a=2. Branches a=3,4 are "
            "excluded by full blocker dualization and a=5,6 by strict "
            "vertex-link g4. Every residual is a connected closed "
            "pseudomanifold whose every face link is a mod-two homology sphere."
        ),
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_exact_eight_orbit_residual",
        "proof_backend": (
            "blocker duality, all 3003 ridge-complement equations, exact "
            "nauty order-nine triangle-free graph quotient, exhaustive cut-pair "
            "enumeration, direct full-blocker dualization, exact automorphism "
            "quotient, and independent GF(2) link homology; no SAT/SMT result"
        ),
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": _file_sha256(SOURCE),
            "second_layer_survivor_index": TARGET_RECORD_INDEX,
        },
        "x_graph_archive": {
            "path": str(X_ARCHIVE.relative_to(ROOT)),
            "file_sha256": _file_sha256(X_ARCHIVE),
            **proof.x_archive.to_record(),
            "trusted_scope": (
                "exhaustive unlabeled generation only; decoding, graph "
                "properties, and residual constraints are independently replayed"
            ),
        },
        "proof_sources": _proof_sources(),
        "proof_dependencies": [
            "Burcroff's no-missing-triples theorem for compact H9/14 nerves",
            "every compact Coxeter missing face is Lanner of order at most five",
            "strict generalized lower-bound inequalities in every vertex link",
            "Dehn--Sommerville schedules for simplicial 9-polytopes",
            "facet complements are the blocker of the full minimal-nonface clutter",
            "every simplicial-sphere ridge has degree exactly two",
            "a graph with zero-or-two edges on every triple is a complete cut",
            "blocker involution for finite clutters",
            "nauty geng exhaustive isomorphism quotient",
        ],
        "schedule_table": schedule_table,
        "strict_vertex_link_g4_exclusion": {
            "required_sum": 14,
            "excluded_missing_four_counts": [5, 6],
            "observed_sums": {
                "5": summed_vertex_link_g4(schedules[5]["f_vector"]),
                "6": summed_vertex_link_g4(schedules[6]["f_vector"]),
            },
        },
        "whole_shadow_audit": proof.to_record(),
        "automorphism_quotient": {
            "missing_edge_graph_automorphism_group": "S9 x C2",
            "reason_distinct_x_words_are_distinct_orbits": (
                "the S9 factor relabels isolates, the C2 swaps the two branches "
                "while fixing vertex 0, and X_ij=y_{0ij} is invariant"
            ),
            "labelled_canonical_x_representations": len(
                proof.enumeration.witnesses
            ),
            "incidence_complex_orbits": len(
                proof.enumeration.residual_orbits
            ),
            "branch_histogram": {"0": 1, "1": 4, "2": 3},
        },
        "full_topology_audits": topology,
        "independent_hostile_replay": {
            "source": str(INDEPENDENT_SOURCE.relative_to(ROOT)),
            "source_sha256": _file_sha256(INDEPENDENT_SOURCE),
            "summary_sha256": hostile["summary_sha256"],
            "blocker_and_ridges": hostile["blocker_and_ridges"],
            "x_archive": hostile["x_archive"],
            "census": hostile["census"],
            "residual": hostile["residual"],
        },
        "correction_note": (
            "An exploratory intermediate report incorrectly called the a=0 "
            "state unique by treating the schedule's missing-five count as "
            "the number of inclusion-minimal A5 rows. For a>0 the schedule "
            "counts all missing five-sets, including the A4 shadow. The direct "
            "blocker and face-table replay retains 17 a=1 and 3 a=2 labelled "
            "states, giving seven additional automorphism orbits."
        ),
        "nonclaims": [
            "A mod-two homology sphere need not be a PL sphere.",
            "No convex realization is supplied by this certificate.",
            "No Coxeter labels, Lorentzian Gram matrix, or compact hyperbolic polytope is supplied.",
            "Finite-volume conclusions are outside this compact-nerve certificate.",
            "This resolves one graph shadow only, not all compact H9/14 nerves.",
        ],
    }
    manifest = json.loads(json.dumps(manifest, sort_keys=True))
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_record17_exact_whole_shadow_residual":
        raise ValueError("wrong record-17 certificate kind")
    if (
        manifest.get("record_index") != TARGET_RECORD_INDEX
        or manifest.get("graph6") != TARGET_GRAPH6
    ):
        raise ValueError("record-17 certificate target changed")
    if manifest.get("theorem_status") != "certified_exact_eight_orbit_residual":
        raise ValueError("record-17 theorem status changed")
    expected_files = {
        SOURCE: EXPECTED_SOURCE_SHA256,
        X_ARCHIVE: EXPECTED_X_ARCHIVE_SHA256,
        PROOF_SOURCE: EXPECTED_PROOF_SOURCE_SHA256,
        FACET_AUDIT_SOURCE: EXPECTED_FACET_AUDIT_SOURCE_SHA256,
        HOMOLOGY_SOURCE: EXPECTED_HOMOLOGY_SOURCE_SHA256,
        LOCAL_LINK_SOURCE: EXPECTED_LOCAL_LINK_SOURCE_SHA256,
        MISSING_EDGES_SOURCE: EXPECTED_MISSING_EDGES_SOURCE_SHA256,
        RECORD92_SOURCE: EXPECTED_RECORD92_SOURCE_SHA256,
        NERVE_SUPPORT_SOURCE: EXPECTED_NERVE_SUPPORT_SOURCE_SHA256,
        GENERATOR_SOURCE: EXPECTED_GENERATOR_SOURCE_SHA256,
        INDEPENDENT_SOURCE: EXPECTED_INDEPENDENT_SOURCE_SHA256,
    }
    for path, expected in expected_files.items():
        observed = _file_sha256(path)
        if expected != "TO_BE_FROZEN" and observed != expected:
            raise ValueError(f"record-17 frozen source changed: {path}")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("record-17 certificate digest mismatch")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and digest != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-17 certificate is not frozen")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("record-17 certificate did not reproduce")
    return {
        "record_index": TARGET_RECORD_INDEX,
        "graph6": TARGET_GRAPH6,
        "theorem_status": manifest["theorem_status"],
        "labelled_representations": manifest["automorphism_quotient"][
            "labelled_canonical_x_representations"
        ],
        "incidence_complex_orbits": manifest["automorphism_quotient"][
            "incidence_complex_orbits"
        ],
        "branch_histogram": manifest["automorphism_quotient"][
            "branch_histogram"
        ],
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
