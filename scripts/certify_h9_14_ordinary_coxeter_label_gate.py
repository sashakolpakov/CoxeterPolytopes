#!/usr/bin/env python3
"""Replay the unified five-record H9/14 Coxeter-label theorem."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from scripts.certify_h9_14_product_coxeter_labels import (
    OUTPUT as PRODUCT_OUTPUT,
    verify as verify_product,
)
from scripts.certify_h9_14_sparse_coxeter_labels import (
    OUTPUT as MULTIWEDGE_OUTPUT,
    verify as verify_multiwedge,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_ordinary_coxeter_label_gate.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "9c3943ea538733020bc09f34c0150cdf6a7d66af743d5643a32ac6033057f4ea"
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


def build() -> dict[str, object]:
    product = json.loads(PRODUCT_OUTPUT.read_text(encoding="utf-8"))
    multiwedge = json.loads(MULTIWEDGE_OUTPUT.read_text(encoding="utf-8"))
    product_summary = verify_product(product)
    multiwedge_summary = verify_multiwedge(multiwedge)
    product_by_index = {record["record_index"]: record for record in product["records"]}
    multiwedge_by_index = {
        record["record_index"]: record for record in multiwedge["records"]
    }

    records = []
    for index in (30, 42):
        source = product_by_index[index]
        family = source["orthogonal_join_subfamily"]
        records.append(
            {
                "record_index": index,
                "rank_nine_vertices_checked": source["rank_nine_vertex_count"],
                "sparse_labelled_solutions": family["raw_labelled_placements"],
                "full_nerve_symmetry_orbits": family["isomorphism_orbits"],
                "local_discrete_status": "sat",
                "sparse_global_gram_status": family["global_gram_obstruction"][
                    "status"
                ],
                "unrestricted_minimum_cross_component_nonright_edges": source[
                    "unrestricted_forced_cross_edge_gate"
                ]["minimum_nonright_meeting_edges_between_components"],
                "forced_next_clause": (
                    "at least "
                    f"{source['unrestricted_forced_cross_edge_gate']['minimum_nonright_meeting_edges_between_components']} "
                    "cross-component non-right meeting edges"
                ),
            }
        )
    for index in (43, 44, 47):
        source = multiwedge_by_index[index]
        records.append(
            {
                "record_index": index,
                "rank_nine_vertices_checked": source["vertex_count"],
                "sparse_labelled_solutions": source["full_labelled_assignment_count"],
                "full_nerve_symmetry_orbits": source["symmetry_orbit_count"],
                "local_discrete_status": "sat",
                "sparse_global_gram_status": "unsat",
                "unrestricted_minimum_cross_component_nonright_edges": source[
                    "unrestricted_forced_cross_edge_gate"
                ]["minimum_nonright_meeting_edges_between_components"],
                "forced_next_clause": (
                    "at least one non-right edge between the isolated "
                    "simplex Lanner factor and the residual factor"
                ),
            }
        )
    if [record["record_index"] for record in records] != [30, 42, 43, 44, 47]:
        raise AssertionError("the five-record label frontier changed")
    if not all(
        record["local_discrete_status"] == "sat"
        and record["sparse_global_gram_status"] == "unsat"
        for record in records
    ):
        raise AssertionError("a five-record local/global status changed")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_five_ordinary_nerve_coxeter_label_gate",
        "claim": (
            "All five known ordinary-polytopal H9/14 nerves admit exact "
            "full-domain factor-supported Coxeter label families satisfying "
            "every finite-vertex and Lanner-minimal-nonface constraint.  Every "
            "member of all five families is globally non-Lorentzian, so a "
            "genuine candidate requires cross-factor non-right coupling.  The "
            "unrestricted minimum counts are 4,2,1,1,1 for records "
            "30,42,43,44,47."
        ),
        "finite_label_domain": {
            "labels": [2, 3, 4, 5],
            "larger_label_exclusion": (
                "Burcroff Lemma 10.10 excludes multi-multiple edges m>=6 "
                "for a compact 9-polytope with 14 facets."
            ),
            "golden_only": False,
        },
        "counts": {
            "ordinary_polytopal_nerves": 5,
            "locally_sat_sparse_families": 5,
            "globally_unsat_sparse_families": 5,
            "sparse_labelled_solutions": sum(
                record["sparse_labelled_solutions"] for record in records
            ),
            "full_nerve_symmetry_orbits": sum(
                record["full_nerve_symmetry_orbits"] for record in records
            ),
        },
        "records": records,
        "detailed_certificates": [
            {
                "path": str(PRODUCT_OUTPUT.relative_to(ROOT)),
                "file_sha256": _file_sha256(PRODUCT_OUTPUT),
                "embedded_certificate_sha256": product_summary[
                    "certificate_sha256"
                ],
                "records": product_summary["record_indices"],
            },
            {
                "path": str(MULTIWEDGE_OUTPUT.relative_to(ROOT)),
                "file_sha256": _file_sha256(MULTIWEDGE_OUTPUT),
                "embedded_certificate_sha256": multiwedge_summary[
                    "certificate_sha256"
                ],
                "records": multiwedge_summary["record_indices"],
            },
        ],
        "proof_summary": {
            "local": (
                "Complete finite Coxeter classification, exhaustive order-4/5 "
                "Lanner generation over {2,3,4,5}, exact placement-table joins, "
                "and all 498 rank-nine vertex replays."
            ),
            "symmetry": (
                "Full minimal-nonface-hypergraph automorphism groups and exact "
                "orbit traversal; complete only within each stated sparse ansatz."
            ),
            "global": (
                "Orthogonal Lanner-bearing blocks force negative index at least "
                "two; record 30 has inertia exactly (9,5,0).  More generally, "
                "minimal-nonface components give unrestricted cross-edge lower "
                "bounds 4,2,1,1,1."
            ),
        },
        "scope": {
            "unrestricted_cross_coupled_label_csp": "open",
            "dotted_gram_variables": "open",
            "rank_10_inertia_9_1_4": "open",
            "new_hyperbolic_coxeter_polytope": False,
        },
    }
    if manifest["counts"] != {
        "ordinary_polytopal_nerves": 5,
        "locally_sat_sparse_families": 5,
        "globally_unsat_sparse_families": 5,
        "sparse_labelled_solutions": 46_540_944,
        "full_nerve_symmetry_orbits": 6_210,
    }:
        raise AssertionError("the unified five-record counts changed")
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_five_ordinary_nerve_coxeter_label_gate":
        raise ValueError("unexpected unified label-gate certificate kind")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("unified label-gate certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("unified label-gate certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("unified label-gate certificate is not frozen")
    return {
        "record_indices": [record["record_index"] for record in rebuilt["records"]],
        "sparse_labelled_solutions": rebuilt["counts"]["sparse_labelled_solutions"],
        "full_nerve_symmetry_orbits": rebuilt["counts"][
            "full_nerve_symmetry_orbits"
        ],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--print-certificate", action="store_true")
    args = parser.parse_args(argv)
    if args.print_certificate:
        print(json.dumps(build(), indent=2, sort_keys=True))
        return
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps(verify(manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
