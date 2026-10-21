#!/usr/bin/env python3
"""Build or verify strict link-g4 exclusions in seven open H9/14 shadows."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.h9_14_record92 import (
    summed_vertex_link_g4,
)
from scripts.independent_h9_14_open_high_schedule_g4 import (
    RECORD_GRAPHS,
    replay,
)


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = ROOT / "data/enumerations/h9_14_open_high_schedule_g4_no_go.json"
SOURCE_PATHS = (
    "src/coxeter_polytopes/combinatorics/h9_14_record92.py",
    "scripts/independent_h9_14_open_high_schedule_g4.py",
)
EXPECTED_EXCLUSIONS = {
    0: (17, 18),
    1: (14, 15),
    3: (14, 15),
    4: (11, 12),
    5: (9, 10),
    7: (8, 9),
    9: (8, 9),
}
EXPECTED_HOSTILE_SUMMARY_SHA256 = (
    "787e3680c788620538c108fdec894662b674315b8df2dd30d3725b64364c5b22"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "a7f9d6c4ca31c91d6d456555703a273a67a70bf16ad74ab54b26afeb735156f5"
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


def _production_rows() -> list[dict[str, object]]:
    document = json.loads(INPUT.read_text(encoding="utf-8"))
    by_graph = {row["graph6"]: row for row in document["records"]}
    rows = []
    for record, graph6 in RECORD_GRAPHS.items():
        source = by_graph[graph6]
        schedules = []
        for schedule in source["dehn_sommerville_schedule"]:
            f_vector = tuple(map(int, schedule["f_vector"]))
            total = summed_vertex_link_g4(f_vector)
            schedules.append(
                {
                    "missing_four_count": int(schedule["missing_four_count"]),
                    "missing_five_count": int(
                        schedule["missing_five_count_offset"]
                    ),
                    "facet_count": f_vector[-1],
                    "summed_vertex_link_g4": total,
                    "strict_sum_feasible": total >= 14,
                }
            )
        exclusions = tuple(
            row["missing_four_count"]
            for row in schedules
            if not row["strict_sum_feasible"]
        )
        if exclusions != EXPECTED_EXCLUSIONS[record]:
            raise AssertionError(f"record {record}: strict link-g4 cutoff changed")
        rows.append(
            {
                "record_index": record,
                "graph6": graph6,
                "schedules": schedules,
                "strict_link_g4_excluded_schedules": list(exclusions),
            }
        )
    return rows


def generate() -> dict[str, object]:
    production = _production_rows()
    hostile = replay()
    if (
        hostile["summary_sha256"] != EXPECTED_HOSTILE_SUMMARY_SHA256
        or hostile["records"] != production
        or hostile["total_schedule_branches_excluded"] != 14
    ):
        raise AssertionError("independent strict link-g4 replay changed")
    input_document = json.loads(INPUT.read_text(encoding="utf-8"))
    document: dict[str, object] = {
        "schema_version": 1,
        "kind": "h9_14_seven_open_shadows_high_schedule_strict_g4_no_go",
        "claim": (
            "The top two Dehn--Sommerville schedules of each of H9/14 "
            "records 0, 1, 3, 4, 5, 7, and 9 are impossible for a compact "
            "hyperbolic Coxeter polytope, excluding fourteen schedule "
            "branches in total."
        ),
        "theorem": {
            "required_strict_summed_vertex_link_g4": 14,
            "records": production,
            "total_schedule_branches_excluded": 14,
            "verified": True,
        },
        "proof": {
            "strict_local_bound": (
                "For every vertex v, lk(v) is the boundary of a nonsimplex "
                "simplicial 8-polytope and inherits the Coxeter--Lanner "
                "missing-face order bound five. If g4(lk(v))=0, the "
                "generalized lower-bound equality theorem gives a nontrivial "
                "3-stacked triangulation and hence a missing link face of "
                "order at least six. Therefore g4(lk(v))>=1 at all fourteen "
                "vertices."
            ),
            "double_count": (
                "For an 8-polytope boundary, g4=f3-6f2+21f1-56f0+126. "
                "Summing over vertex links and double-counting incidences "
                "gives sum_v g4(lk(v)) = 5 f4 - 24 f3 + 63 f2 - "
                "112 f1 + 14*126. Hence every compact schedule must have "
                "summed value at least fourteen."
            ),
            "independent_replay": (
                "A no-project-import verifier derives h3 and h4 directly "
                "from the link h-vector definition, separately checks the "
                "expanded polynomial, and reproduces every schedule row."
            ),
        },
        "independent_replay": hostile,
        "scope": {
            "compact_only": True,
            "record_indices": list(RECORD_GRAPHS),
            "listed_schedule_branches_excluded": True,
            "whole_graph_shadows_excluded": False,
            "remaining_schedules_excluded": False,
            "finite_volume_noncompact_chambers_excluded": False,
            "new_hyperbolic_coxeter_polytope": False,
        },
        "input": {
            "path": str(INPUT.relative_to(ROOT)),
            "file_sha256": _file_sha256(INPUT),
            "certificate_sha256": input_document["certificate_sha256"],
        },
        "sources": _path_records(SOURCE_PATHS),
    }
    document["certificate_sha256"] = certificate_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_seven_open_shadows_high_schedule_strict_g4_no_go":
        raise ValueError("unexpected strict link-g4 certificate kind")
    if document.get("sources") != _path_records(SOURCE_PATHS):
        raise ValueError("strict link-g4 source closure changed")
    input_document = json.loads(INPUT.read_text(encoding="utf-8"))
    if document.get("input") != {
        "path": str(INPUT.relative_to(ROOT)),
        "file_sha256": _file_sha256(INPUT),
        "certificate_sha256": input_document["certificate_sha256"],
    }:
        raise ValueError("strict link-g4 input closure changed")
    if document.get("certificate_sha256") != certificate_digest(document):
        raise ValueError("strict link-g4 certificate digest changed")
    if (
        EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN"
        and document.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256
    ):
        raise ValueError("strict link-g4 certificate is not frozen")
    if document != generate():
        raise ValueError("strict link-g4 certificate does not replay")
    return {
        "record_exclusions": {
            str(row["record_index"]): row["strict_link_g4_excluded_schedules"]
            for row in document["theorem"]["records"]
        },
        "total_schedule_branches_excluded": document["theorem"][
            "total_schedule_branches_excluded"
        ],
        "hostile_summary_sha256": document["independent_replay"][
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
