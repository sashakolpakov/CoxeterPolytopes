#!/usr/bin/env python3
"""Freeze the sparse full-domain Coxeter-label censuses for records 43/44/47."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_sparse_coxeter_labels import (
    audit_all_sparse_h9_14_labels,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_sparse_coxeter_labels.json"
SOURCES = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_sparse_coxeter_labels.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record43.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record44.py",
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record47.py",
)
EXPECTED_CERTIFICATE_SHA256 = (
    "4362c36bc471c9dfd9b414fb4ba0439880a66b5fec0310466e7d1119dddc2c63"
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
    audits = audit_all_sparse_h9_14_labels()
    if [audit.record_index for audit in audits] != [43, 44, 47] or not all(
        audit.verified for audit in audits
    ):
        raise AssertionError("the three sparse label censuses did not verify")
    records = [audit.to_record() for audit in audits]
    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_sparse_full_domain_coxeter_label_census",
        "claim": (
            "For each of the ordinary-polytopal H9/14 multiwedge records "
            "43, 44, and 47, the sparse ansatz of Lanner labels only on the "
            "four missing five-faces is locally satisfiable.  Its complete "
            "full-domain census has 12600, 60480, and 7560 labelled solutions, "
            "or 15, 25, and 5 full-nerve symmetry orbits.  Every such solution "
            "fails the global Lorentzian gate by negative index at least two."
        ),
        "finite_label_domain": {
            "labels": [2, 3, 4, 5],
            "scope": "the complete Burcroff domain, not the golden subalphabet",
            "exclusion_of_larger_labels": (
                "Burcroff Lemma 10.10 excludes multi-multiple edges, hence "
                "finite Coxeter labels m>=6, in the compact H9/14 problem."
            ),
            "reference": (
                "A. Burcroff, Compact hyperbolic Coxeter d-polytopes with "
                "d+4 facets and related dimension bounds, arXiv:2201.03437, "
                "Lemma 10.10"
            ),
        },
        "ansatz": {
            "definition": (
                "Put an arbitrary order-five Lanner placement on each of the "
                "four missing five-faces, require equality on overlapping "
                "pairs, and label every other meeting pair by 2."
            ),
            "complete_generation": (
                "Join all 420-row Lanner tables on the three overlapping "
                "supports without heuristic pruning; test all 18 residual "
                "vertex diagrams; independently multiply by all 420 placements "
                "on the isolated boundary-of-Delta4 factor.  The archived 90 "
                "vertices are replayed as the exact Cartesian product 5x18."
            ),
            "symmetry_quotient": (
                "Explicit generators are checked on the complete vertex list "
                "and minimal-nonface hypergraph.  Their generated groups have "
                "orders 1440, 5760, and 4320; complete orbit traversal exports "
                "one canonical label map per orbit."
            ),
        },
        "counts": {
            "records": 3,
            "full_labelled_solutions": sum(
                audit.full_labelled_assignment_count for audit in audits
            ),
            "symmetry_orbits": sum(audit.symmetry_orbit_count for audit in audits),
            "all_90_vertex_representatives_pass": True,
            "all_sparse_families_globally_obstructed": True,
        },
        "records": records,
        "sources": [
            {"path": str(path.relative_to(ROOT)), "sha256": _file_sha256(path)}
            for path in SOURCES
        ],
        "scope": {
            "discrete_sparse_label_status": "sat",
            "global_sparse_gram_status": "unsat",
            "full_unrestricted_discrete_label_census": "not claimed",
            "dotted_values": "not assigned",
            "new_hyperbolic_coxeter_polytope": False,
            "next_gate": (
                "At least one non-right edge joining the isolated simplex "
                "Lanner factor to the residual factor is necessary."
            ),
        },
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_sparse_full_domain_coxeter_label_census":
        raise ValueError("unexpected sparse Coxeter-label certificate kind")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("sparse Coxeter-label certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("sparse Coxeter-label certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != (
        EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("sparse Coxeter-label certificate is not frozen")
    return {
        "record_indices": [record["record_index"] for record in rebuilt["records"]],
        "full_labelled_solutions": rebuilt["counts"]["full_labelled_solutions"],
        "symmetry_orbits": rebuilt["counts"]["symmetry_orbits"],
        "certificate_sha256": digest,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--print-certificate",
        action="store_true",
        help="print a freshly rebuilt certificate instead of replaying the archive",
    )
    args = parser.parse_args(argv)
    if args.print_certificate:
        print(json.dumps(build(), indent=2, sort_keys=True))
        return
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps(verify(manifest), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
