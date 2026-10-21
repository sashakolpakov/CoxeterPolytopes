#!/usr/bin/env python3
"""Build or verify the all-orbit compact-Coxeter no-go for record 17."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record17_all_coxeter_labels import (
    audit_record17_all_coxeter_labels,
)
from scripts.independent_h9_14_record17_all_coxeter_audit import (
    replay as independent_replay,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/enumerations/h9_14_record17_certificate.json"
OUTPUT = ROOT / "data/enumerations/h9_14_record17_all_coxeter_labels.json"
SOURCE_PATHS = (
    "src/coxeter_polytopes/combinatorics/h9_14_record17_all_coxeter_labels.py",
    "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    "scripts/independent_h9_14_record17_all_coxeter_audit.py",
)
INPUT_PATHS = ("data/enumerations/h9_14_record17_certificate.json",)
EXPECTED_HOSTILE_SUMMARY_SHA256 = (
    "e19c45bb7134690b3649b9476c69829072076a777335c124216aeee5b0e3b609"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "59a0ebdd63a28914412cda34d60289e06a539dcf3ae4c8d65debc9e7b2f7d820"
)


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _path_records(paths: Sequence[str]) -> list[dict[str, str]]:
    return [
        {"path": relative, "sha256": _file_sha256(ROOT / relative)}
        for relative in paths
    ]


def certificate_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


def generate() -> dict[str, object]:
    audit = audit_record17_all_coxeter_labels(INPUT)
    if not audit.verified or not audit.whole_shadow_compact_coxeter_no_go:
        raise AssertionError("record-17 whole-shadow compact no-go failed")
    hostile = independent_replay(INPUT)
    if (
        not hostile["all_eight_orbits_excluded"]
        or hostile["summary_sha256"] != EXPECTED_HOSTILE_SUMMARY_SHA256
    ):
        raise AssertionError("record-17 independent all-orbit replay changed")

    document: dict[str, object] = {
        "schema_version": 1,
        "kind": "h9_14_record17_all_orbit_compact_coxeter_no_go",
        "claim": (
            "None of the exactly eight residual incidence-complex orbits for "
            "H9/14 second-layer record 17 is the nerve of a compact "
            "hyperbolic Coxeter polytope; hence the whole graph shadow is "
            "compact-Coxeter impossible."
        ),
        "theorem": audit.to_record(),
        "proof": {
            "compact_inputs": [
                "Every order-four/five minimal nonface of a compact Coxeter nerve supports a classical Lanner diagram.",
                "Every actual rank-nine vertex diagram of a compact Coxeter polytope is elliptic.",
                "The complete classical Lanner tables have 108 and 420 labelled placements in orders four and five.",
                "The frozen whole-shadow theorem leaves exactly the eight audited incidence-complex orbits and no others.",
            ],
            "finite_search": (
                "For each orbit, join the complete Lanner placement tables "
                "on the listed decisive supports and replay every restriction "
                "induced by every actual rank-nine face. Five orbits have "
                "zero elliptic states after three supports; the other three "
                "have 18, 8, and 18 states after three supports and zero "
                "after the fourth."
            ),
            "perron_frobenius_step": (
                "On any actual-vertex restriction write the normalized "
                "Coxeter Gram matrix as G=I-A, where A is symmetric and "
                "entrywise nonnegative. G is positive definite exactly when "
                "lambda_max(A)=rho(A)<1. Replacing an unforced right angle "
                "by any finite non-right angle gives A'>=A entrywise, hence "
                "rho(A')>=rho(A). A sparse restriction that is already "
                "non-elliptic therefore cannot become elliptic. This proves "
                "that incremental sparse-shadow pruning is lossless."
            ),
            "why_no_global_gram_search_is_needed": (
                "Every possible finite Coxeter labeling already fails a "
                "necessary positive-definite principal vertex condition."
            ),
        },
        "independent_hostile_replay": hostile,
        "scope": {
            "record_index": 17,
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "selected_nerve_only": False,
            "all_eight_residual_incidence_orbits_covered": True,
            "compact_only": True,
            "finite_volume_noncompact_chamber_excluded": False,
            "pl_sphericity_of_abstract_residuals_claimed": False,
            "convex_polytopality_of_abstract_residuals_claimed": False,
            "global_lorentzian_gram_values_enumerated": False,
        },
        "sources": _path_records(SOURCE_PATHS),
        "inputs": _path_records(INPUT_PATHS),
    }
    document["certificate_sha256"] = certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record17_all_orbit_compact_coxeter_no_go":
        raise ValueError("unexpected record-17 all-orbit certificate kind")
    if document.get("sources") != _path_records(SOURCE_PATHS):
        raise ValueError("record-17 all-orbit proof source closure changed")
    if document.get("inputs") != _path_records(INPUT_PATHS):
        raise ValueError("record-17 all-orbit proof input closure changed")
    if document.get("certificate_sha256") != certificate_digest(document):
        raise ValueError("record-17 all-orbit certificate digest changed")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and document.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-17 all-orbit certificate is not frozen")
    rebuilt = generate()
    if document != rebuilt:
        raise ValueError("record-17 all-orbit certificate does not replay")
    theorem = document["theorem"]
    if (
        not theorem["verified"]
        or not theorem["whole_shadow_compact_coxeter_no_go"]
        or theorem["orbits_excluded"] != 8
        or not document["independent_hostile_replay"][
            "all_eight_orbits_excluded"
        ]
    ):
        raise ValueError("record-17 all-orbit theorem disappeared")
    return {
        "record_index": document["scope"]["record_index"],
        "residual_orbits": len(theorem["orbit_gates"]),
        "final_survivors": [
            orbit["incremental_lanner_and_actual_vertex_gate"][-1][
                "elliptic_actual_vertex_survivor_count"
            ]
            for orbit in theorem["orbit_gates"]
        ],
        "whole_graph_shadow_compact_coxeter_excluded": document["scope"][
            "whole_graph_shadow_compact_coxeter_excluded"
        ],
        "hostile_summary_sha256": document["independent_hostile_replay"][
            "summary_sha256"
        ],
        "certificate_sha256": document["certificate_sha256"],
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    document = (
        json.loads(OUTPUT.read_text(encoding="utf-8"))
        if args.verify
        else generate()
    )
    if args.write:
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
