#!/usr/bin/env python3
"""No-project-import sparse Coxeter-label census for all 77 Record-30 a=2 targets."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import argparse
import json
from pathlib import Path
import sys
from typing import Mapping


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.independent_h9_14_record30_a2_orbit14_labels as engine


EXPECTED_SURVIVORS = ((1, 0), (7, 1), (10, 0))


def _canonical_sha(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if mask >> vertex & 1)


def replay(document: Mapping[str, object]) -> dict[str, object]:
    raw_rows = document.get("orbit_audits")
    if not isinstance(raw_rows, list) or len(raw_rows) != 21:
        raise ValueError("expected all 21 positive-A4 orbit audits")
    audits = []
    coordinates = []
    global_index = 0
    for orbit_index, raw_row in enumerate(raw_rows):
        if not isinstance(raw_row, Mapping):
            raise ValueError(f"orbit row {orbit_index} is malformed")
        if int(raw_row.get("positive_A4_orbit_index", -1)) != orbit_index:
            raise ValueError("positive-A4 orbit order changed")
        fixed_masks = raw_row.get("fixed_A4_masks")
        families = raw_row.get("pure_target_representatives")
        if (
            not isinstance(fixed_masks, list)
            or len(fixed_masks) != 2
            or not isinstance(families, list)
        ):
            raise ValueError(f"orbit row {orbit_index} target data is malformed")
        fours = tuple(_vertices(int(value)) for value in fixed_masks)
        for target_index, raw_family in enumerate(families):
            if not isinstance(raw_family, list):
                raise ValueError("target family is malformed")
            family = tuple(sorted(int(face) for face in raw_family))
            audit = engine._label_family(
                global_index,
                family,
                fours=fours,
                require_no_go=False,
            )
            audit["positive_A4_orbit_index"] = orbit_index
            audit["orbit_target_index"] = target_index
            audits.append(audit)
            coordinates.append(
                {
                    "positive_A4_orbit_index": orbit_index,
                    "orbit_target_index": target_index,
                    "fixed_A4_masks": [int(value) for value in fixed_masks],
                    "additional_A5_masks": list(family),
                }
            )
            global_index += 1
    if global_index != 77:
        raise AssertionError("complete a=2 target count changed")
    survivors = tuple(
        (int(row["positive_A4_orbit_index"]), int(row["orbit_target_index"]))
        for row in audits
        if int(row["terminal_sparse_state_count"]) > 0
    )
    if survivors != EXPECTED_SURVIVORS:
        raise AssertionError("sparse label survivor list changed")
    if any(int(row["actual_rank_nine_vertex_count"]) != 128 for row in audits):
        raise AssertionError("a target lost the scheduled 128 actual vertices")
    transcript = [
        {
            "positive_A4_orbit_index": row["positive_A4_orbit_index"],
            "orbit_target_index": row["orbit_target_index"],
            "stages": row["stages"],
            "terminal_sparse_state_count": row["terminal_sparse_state_count"],
        }
        for row in audits
    ]
    terminal_histogram = Counter(
        int(row["terminal_sparse_state_count"]) for row in audits
    )
    result = {
        "kind": "h9_14_record30_a2_all_target_sparse_label_replay",
        "compact_only": True,
        "imports_coxeter_polytopes": False,
        "complete_order_four_lanner_placement_count": len(
            engine._lanner_codes(4)
        ),
        "complete_order_five_lanner_placement_count": len(
            engine._lanner_codes(5)
        ),
        "target_orbit_count": len(audits),
        "target_coordinate_sha256": _canonical_sha(coordinates),
        "target_audits": audits,
        "label_stage_transcript_sha256": _canonical_sha(transcript),
        "terminal_sparse_state_count_histogram": {
            str(key): value for key, value in sorted(terminal_histogram.items())
        },
        "sparse_label_no_go_target_count": sum(
            bool(row["unrestricted_compact_label_no_go"]) for row in audits
        ),
        "sparse_label_survivor_count": len(survivors),
        "sparse_label_survivors": [
            {
                "positive_A4_orbit_index": orbit_index,
                "orbit_target_index": target_index,
                "terminal_sparse_state_count": next(
                    int(row["terminal_sparse_state_count"])
                    for row in audits
                    if row["positive_A4_orbit_index"] == orbit_index
                    and row["orbit_target_index"] == target_index
                ),
            }
            for orbit_index, target_index in survivors
        ],
        "scope": (
            "Complete 108/420 Lanner-placement and actual-rank-nine-vertex "
            "sparse join for all 77 compact Record-30 a=2 incidence target "
            "orbits. Seventy-four targets die by Perron--Frobenius "
            "monotonicity; the three named survivors require separate "
            "whole-nerve metric certificates. No finite-volume inference."
        ),
    }
    result["summary_sha256"] = _canonical_sha(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    args = parser.parse_args()
    document = json.loads(args.input.read_text(encoding="utf-8"))
    print(json.dumps(replay(document), sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
