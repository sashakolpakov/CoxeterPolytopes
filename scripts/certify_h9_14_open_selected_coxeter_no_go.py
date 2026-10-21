#!/usr/bin/env python3
"""Build or verify seven selected H9/14 compact-Coxeter no-gos."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_open_selected_coxeter_no_go import (
    TARGET_RECORDS,
    audit_all_open_selected_nerves,
)
from scripts.independent_h9_14_open_selected_coxeter_audit import replay


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_open_selected_coxeter_no_go.json"
INPUT_PATHS = tuple(
    f"data/enumerations/h9_14_record{record}_polytopality.json"
    for record in TARGET_RECORDS
)
SOURCE_PATHS = (
    "src/coxeter_polytopes/combinatorics/h9_14_open_selected_coxeter_no_go.py",
    "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    "scripts/independent_h9_14_open_selected_coxeter_audit.py",
    "scripts/independent_h9_14_record17_all_coxeter_audit.py",
)
EXPECTED_HOSTILE_SUMMARY_SHA256 = (
    "29c5bc8f79024d92c33b2d2d9fb4f561600658322475e91e0d710ec9d9ea4b26"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "c0f693a5ec92abb0f146d7919d7c2a22910fa07c0576d402a9700b18ef142ed7"
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


def _compare_replays(
    production: Mapping[str, object], hostile: Mapping[str, object]
) -> None:
    if hostile["record_indices"] != production["record_indices"]:
        raise AssertionError("production/hostile record coverage differs")
    if hostile["lanner_table_counts"] != {"4": 108, "5": 420}:
        raise AssertionError("hostile Lanner census changed")
    for first, second in zip(
        production["record_gates"], hostile["record_gates"], strict=True
    ):
        if (
            first["record_index"] != second["record_index"]
            or first["rank_nine_vertex_count"] != second["rank_nine_vertex_count"]
            or first["minimal_nonface_histogram"]
            != second["minimal_nonface_histogram"]
            or first["decisive_support_indices"]
            != second["decisive_support_indices"]
            or first["incremental_lanner_and_actual_vertex_gate"]
            != second["stages"]
            or not second["selected_nerve_compact_coxeter_no_go"]
        ):
            raise AssertionError(
                f"record {first['record_index']}: hostile replay disagrees"
            )


def generate() -> dict[str, object]:
    audit = audit_all_open_selected_nerves()
    if not audit.verified or not audit.all_seven_selected_nerves_excluded:
        raise AssertionError("the production selected-nerve theorem failed")
    theorem = audit.to_record()
    hostile = replay()
    if (
        hostile["summary_sha256"] != EXPECTED_HOSTILE_SUMMARY_SHA256
        or not hostile["all_seven_selected_nerves_excluded"]
    ):
        raise AssertionError("the hostile selected-nerve replay changed")
    _compare_replays(theorem, hostile)

    document: dict[str, object] = {
        "schema_version": 1,
        "kind": "h9_14_seven_selected_nerve_compact_coxeter_no_go",
        "claim": (
            "The certified ordinary-polytopal selected nerves for H9/14 "
            "records 0, 1, 3, 4, 5, 7, and 9 admit no compact hyperbolic "
            "Coxeter labeling."
        ),
        "theorem": theorem,
        "proof": {
            "compact_inputs": [
                "Every inclusion-minimal order-four/five nonface in a compact Coxeter nerve supports a classical Lanner diagram.",
                "Every actual rank-nine vertex diagram of a compact hyperbolic Coxeter 9-polytope is elliptic.",
                "For compact H9/14 the complete finite label alphabet is 2, 3, 4, 5.",
                "The complete labelled classical Lanner tables contain 108 order-four and 420 order-five placements.",
            ],
            "finite_search": (
                "For each selected nerve, reconstruct all rank-nine vertices "
                "from its exact minimal-nonface clutter. Join complete Lanner "
                "placement tables in the frozen overlap order and, after each "
                "join, require every restriction induced by every actual "
                "rank-nine vertex to remain elliptic. The final state counts "
                "for all seven records are zero."
            ),
            "perron_frobenius_step": (
                "On an actual-vertex restriction write the normalized Coxeter "
                "Gram matrix as G=I-A, with A symmetric and entrywise "
                "nonnegative. Ellipticity is rho(A)<1. Replacing an unforced "
                "right angle by any finite non-right angle gives A'>=A "
                "entrywise and rho(A')>=rho(A). Therefore a sparse restriction "
                "already rejected as non-elliptic cannot be repaired by labels "
                "outside the joined supports."
            ),
            "why_no_gram_stage_is_needed": (
                "Every finite Coxeter labeling fails a necessary positive-"
                "definite principal vertex condition before global Lorentzian "
                "rank or inertia can be tested."
            ),
        },
        "independent_hostile_replay": hostile,
        "scope": {
            "record_indices": list(TARGET_RECORDS),
            "selected_nerves_only": True,
            "whole_graph_shadows_compact_coxeter_excluded": False,
            "alternative_completions_excluded": False,
            "compact_only": True,
            "finite_volume_noncompact_chambers_excluded": False,
            "global_lorentzian_gram_values_enumerated": False,
        },
        "sources": _path_records(SOURCE_PATHS),
        "inputs": _path_records(INPUT_PATHS),
    }
    document["certificate_sha256"] = certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_seven_selected_nerve_compact_coxeter_no_go":
        raise ValueError("unexpected selected-nerve certificate kind")
    if document.get("sources") != _path_records(SOURCE_PATHS):
        raise ValueError("selected-nerve proof source closure changed")
    if document.get("inputs") != _path_records(INPUT_PATHS):
        raise ValueError("selected-nerve proof input closure changed")
    if document.get("certificate_sha256") != certificate_digest(document):
        raise ValueError("selected-nerve certificate digest changed")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and document.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("selected-nerve certificate is not frozen")
    rebuilt = generate()
    if document != rebuilt:
        raise ValueError("selected-nerve certificate does not replay")
    theorem = document["theorem"]
    if (
        not theorem["verified"]
        or not theorem["all_seven_selected_nerves_excluded"]
        or theorem["selected_nerves_excluded"] != 7
        or not document["independent_hostile_replay"][
            "all_seven_selected_nerves_excluded"
        ]
    ):
        raise ValueError("the sevenfold selected-nerve theorem disappeared")
    return {
        "record_indices": document["scope"]["record_indices"],
        "final_survivors": [
            row["incremental_lanner_and_actual_vertex_gate"][-1][
                "elliptic_actual_vertex_survivor_count"
            ]
            for row in theorem["record_gates"]
        ],
        "selected_nerves_only": document["scope"]["selected_nerves_only"],
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
