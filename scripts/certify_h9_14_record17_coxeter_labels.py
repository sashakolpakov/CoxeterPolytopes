#!/usr/bin/env python3
"""Build or verify the selected record-17 compact-Coxeter no-go."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record17_coxeter_labels import (
    audit_record17_coxeter_labels,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record17_coxeter_labels.json"
SOURCE_PATHS = (
    "src/coxeter_polytopes/combinatorics/h9_14_record17_coxeter_labels.py",
    "src/coxeter_polytopes/combinatorics/h9_14_record17_polytopality.py",
    "src/coxeter_polytopes/combinatorics/h9_14_record35_coxeter_labels.py",
    "src/coxeter_polytopes/combinatorics/h9_14_coxeter_labels.py",
    "src/coxeter_polytopes/combinatorics/h9_14_record35.py",
    "scripts/independent_h9_14_record17_coxeter_labels_audit.py",
)
INPUT_PATHS = (
    "data/enumerations/h9_14_record17_polytopality.json",
)
EXPECTED_CERTIFICATE_SHA256 = (
    "5ae3eb4bb3dabeb0438fc56e8c92a228d88220a9e946040fc5f55a650e933437"
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
    audit = audit_record17_coxeter_labels()
    if not audit.verified or not audit.unrestricted_compact_coxeter_no_go:
        raise AssertionError("record-17 compact-Coxeter no-go failed")
    document: dict[str, object] = {
        "schema_version": 1,
        "kind": "h9_14_record17_selected_nerve_compact_coxeter_no_go",
        "claim": (
            "The exact ordinary-polytopal selected H9/14 record-17 nerve is "
            "not the nerve of a compact hyperbolic Coxeter polytope."
        ),
        "theorem": audit.to_record(),
        "proof": {
            "compact_inputs": [
                "Every order-five missing face of a compact Coxeter polytope supports a Lanner diagram.",
                "Every actual rank-nine vertex diagram is elliptic.",
                "The complete classical order-five Lanner table has 420 labelled placements.",
            ],
            "exhaustive_gate": (
                "Join all placements on the three supports with common core "
                "{3,4,5,6}; the compatible-state counts are 420, 1332, and "
                "4620. Exact finite-Coxeter classification leaves 80 states "
                "after rank-five restrictions and zero after the four "
                "rank-six restrictions arising in actual vertices."
            ),
            "unrestricted_step": (
                "The decisive seven-facet union has no dotted pair. In "
                "normalized Gram form G=I-A, all unforced entries are zero "
                "in the sparse shadow, and every finite non-right label "
                "increases A entrywise. Perron-Frobenius monotonicity cannot "
                "restore ellipticity."
            ),
        },
        "scope": {
            "record_index": 17,
            "selected_nerve_only": True,
            "compact_coxeter_realization_exists": False,
            "whole_graph_shadow_excluded": False,
            "finite_volume_noncompact_chamber_excluded": False,
            "global_dotted_gram_values_enumerated": False,
            "why_no_global_gram_search_is_needed": (
                "The contradiction is already a finite actual-vertex "
                "obstruction before dotted Gram entries occur."
            ),
        },
        "sources": _path_records(SOURCE_PATHS),
        "inputs": _path_records(INPUT_PATHS),
    }
    document["certificate_sha256"] = certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record17_selected_nerve_compact_coxeter_no_go":
        raise ValueError("unexpected record-17 Coxeter-label certificate kind")
    if document.get("sources") != _path_records(SOURCE_PATHS):
        raise ValueError("record-17 Coxeter-label source closure changed")
    if document.get("inputs") != _path_records(INPUT_PATHS):
        raise ValueError("record-17 Coxeter-label input closure changed")
    if document.get("certificate_sha256") != certificate_digest(document):
        raise ValueError("record-17 Coxeter-label certificate digest changed")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and document.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("record-17 Coxeter-label certificate is not frozen")
    rebuilt = generate()
    if document != rebuilt:
        raise ValueError("record-17 Coxeter-label certificate does not replay")
    theorem = document["theorem"]
    no_go = theorem["unrestricted_compact_coxeter_no_go"]
    gate = theorem["decisive_common_four_core_gate"]
    if not theorem["verified"] or not no_go["proved"]:
        raise ValueError("record-17 compact-Coxeter theorem disappeared")
    return {
        "record_index": document["scope"]["record_index"],
        "rank_nine_vertices": theorem["rank_nine_vertex_count"],
        "joined_lanner_states": gate["join_counts"][-1],
        "rank_five_survivors": gate["cumulative_restriction_stages"][1][
            "cumulative_passing_state_count"
        ],
        "rank_six_survivors": gate["cumulative_restriction_stages"][2][
            "cumulative_passing_state_count"
        ],
        "compact_coxeter_realization_exists": document["scope"][
            "compact_coxeter_realization_exists"
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
