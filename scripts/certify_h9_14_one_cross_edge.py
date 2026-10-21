#!/usr/bin/env python3
"""Freeze the exact one-cross-edge shell around five H9/14 label families."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from coxeter_polytopes.combinatorics.h9_14_one_cross_edge import (
    audit_all_one_cross_edges,
)


ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "src/coxeter_polytopes/combinatorics/h9_14_one_cross_edge.py"
PRODUCT_BASELINE = ROOT / "data/enumerations/h9_14_product_coxeter_labels.json"
MULTIWEDGE_BASELINE = ROOT / "data/enumerations/h9_14_sparse_coxeter_labels.json"
OUTPUT = ROOT / "data/enumerations/h9_14_one_cross_edge.json"
EXPECTED_CERTIFICATE_SHA256 = (
    "b6b440c8649725bbc49a43662036dc8765423b252e7c8092d85f8f7436e3e4f7"
)


EXPECTED = {
    30: (2_659_392, 1_218, 575_424, 180, 4),
    42: (9_057_484_800, 603_840, 78_520_320, 6_120, 2),
    43: (1_701_000, 1_440, 19_440, 16, 2),
    44: (8_164_800, 2_220, 95_040, 20, 2),
    47: (1_020_600, 360, 16_200, 6, 2),
}


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
    records = list(audit_all_one_cross_edges())
    if [record["record_index"] for record in records] != [30, 42, 43, 44, 47]:
        raise AssertionError("the one-edge record order changed")
    for record in records:
        index = record["record_index"]
        shell = record["one_edge_shell"]
        metric = record["metric_principal_gate"]
        observed = (
            shell["full_labelled_trial_count"],
            shell["full_nerve_symmetry_orbit_count"],
            shell["locally_admissible_full_labelled_count"],
            shell["locally_admissible_full_nerve_symmetry_orbit_count"],
            metric.get("record_negative_index_lower_bound", 2),
        )
        if observed != EXPECTED[index]:
            raise AssertionError(f"record-{index} one-edge census changed: {observed}")
        if metric["not_excluded_orbit_count"] != 0:
            raise AssertionError(f"record-{index} retained a metric survivor")

    total_trials = sum(
        record["one_edge_shell"]["full_labelled_trial_count"] for record in records
    )
    total_orbits = sum(
        record["one_edge_shell"]["full_nerve_symmetry_orbit_count"]
        for record in records
    )
    total_local_labelled = sum(
        record["one_edge_shell"]["locally_admissible_full_labelled_count"]
        for record in records
    )
    total_local_orbits = sum(
        record["one_edge_shell"][
            "locally_admissible_full_nerve_symmetry_orbit_count"
        ]
        for record in records
    )
    if (total_trials, total_orbits, total_local_labelled, total_local_orbits) != (
        9_071_030_592,
        609_078,
        79_226_424,
        6_342,
    ):
        raise AssertionError("the five-record one-edge totals changed")

    schur_methods = {
        str(record["record_index"]): dict(
            sorted(
                Counter(
                    survivor["metric_principal_test"]["method"]
                    for survivor in record["survivor_orbit_representatives"]
                ).items()
            )
        )
        for record in records
    }
    if schur_methods["43"] != {
        "exact_two_lanner_schur_factor": 1,
        "untouched_orthogonal_lanner_pair": 15,
    } or schur_methods["44"] != {
        "exact_two_lanner_schur_factor": 1,
        "untouched_orthogonal_lanner_pair": 19,
    } or schur_methods["47"] != {"untouched_orthogonal_lanner_pair": 6}:
        raise AssertionError("the multiwedge metric-witness partition changed")

    manifest: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_exact_one_cross_edge_shell",
        "claim": (
            "Starting from every full-nerve symmetry orbit in the five frozen "
            "factor-orthogonal/sparse Coxeter-label families, change exactly "
            "one right-angled cross-factor pair to m=3,4,or 5.  The exact "
            "labelled-stabilizer quotient has 609078 extension orbits.  Of "
            "these, 6342 satisfy every rank-nine finite-vertex and higher "
            "Lanner-minimal-nonface test, but every one has a principal Gram "
            "submatrix of negative index at least two.  Therefore no member "
            "of this one-cross-edge shell has Lorentzian negative index one."
        ),
        "theorem": {
            "statement": (
                "Relative to each frozen baseline family, any Lorentzian Gram "
                "candidate must change at least two cross-factor right angles "
                "to non-right finite Coxeter edges."
            ),
            "scope_warning": (
                "This is a complete one-edge shell, not an unrestricted label "
                "census.  In particular, for records 43,44,47 it does not "
                "cover simultaneous changes to residual meeting labels outside "
                "the four Lanner supports."
            ),
            "not_a_new_polytope": True,
        },
        "method": {
            "discrete_domain": [2, 3, 4, 5],
            "new_cross_edge_labels": [3, 4, 5],
            "quotient": (
                "For each frozen baseline orbit representative, generate its "
                "full labelled stabilizer inside the full minimal-nonface "
                "automorphism group and orbit the cross pairs under that "
                "stabilizer. Orbit-stabilizer weights independently recover "
                "the full labelled baseline and extension counts."
            ),
            "local_gate": (
                "All unchanged vertices reuse the frozen exact audit. Every "
                "rank-nine vertex containing the changed pair is retested with "
                "the complete finite Coxeter classifier. A cross-factor pair "
                "lies in no higher minimal nonface, so every Lanner support is "
                "unchanged and replayed from the baseline."
            ),
            "metric_gate": (
                "Use Cauchy interlacing on explicit principal blocks. Product "
                "records retain enough orthogonal negative factor blocks. "
                "Multiwedge survivors use either two untouched orthogonal "
                "Lanner blocks or the exact one-edge Schur determinant factor "
                "for two Lanner blocks over Q(sqrt(2),sqrt(5))."
            ),
            "multiwedge_metric_witness_method_counts": schur_methods,
            "nontrivial_exact_schur_factor": (
                "(50 - 22*sqrt(5))/(16*(2 - sqrt(5))**2) > 0; "
                "positivity follows from 50^2 > (22*sqrt(5))^2."
            ),
        },
        "totals": {
            "full_labelled_trials": total_trials,
            "full_nerve_symmetry_orbits": total_orbits,
            "locally_admissible_full_labelled_extensions": total_local_labelled,
            "locally_admissible_full_nerve_symmetry_orbits": total_local_orbits,
            "metric_survivor_orbits": 0,
        },
        "records": records,
        "source": {
            "module_path": str(MODULE.relative_to(ROOT)),
            "module_sha256": _file_sha256(MODULE),
            "product_baseline_path": str(PRODUCT_BASELINE.relative_to(ROOT)),
            "product_baseline_sha256": _file_sha256(PRODUCT_BASELINE),
            "multiwedge_baseline_path": str(MULTIWEDGE_BASELINE.relative_to(ROOT)),
            "multiwedge_baseline_sha256": _file_sha256(MULTIWEDGE_BASELINE),
        },
    }
    manifest["certificate_sha256"] = _certificate_digest(manifest)
    return manifest


def verify(manifest: Mapping[str, object]) -> dict[str, object]:
    if manifest.get("kind") != "h9_14_exact_one_cross_edge_shell":
        raise ValueError("unexpected one-cross-edge certificate kind")
    rebuilt = build()
    if manifest != rebuilt:
        raise ValueError("one-cross-edge certificate did not reproduce")
    digest = _certificate_digest(manifest)
    if manifest.get("certificate_sha256") != digest:
        raise ValueError("one-cross-edge certificate digest mismatch")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("one-cross-edge certificate is not frozen")
    return {
        "record_indices": [record["record_index"] for record in rebuilt["records"]],
        "extension_orbits": rebuilt["totals"]["full_nerve_symmetry_orbits"],
        "local_survivor_orbits": rebuilt["totals"][
            "locally_admissible_full_nerve_symmetry_orbits"
        ],
        "metric_survivor_orbits": rebuilt["totals"]["metric_survivor_orbits"],
        "certificate_sha256": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
        print(json.dumps(verify(manifest), indent=2, sort_keys=True))
        return
    manifest = build()
    if args.write:
        OUTPUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(OUTPUT.relative_to(ROOT)),
                "written": args.write,
                "totals": manifest["totals"],
                "certificate_sha256": manifest["certificate_sha256"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
