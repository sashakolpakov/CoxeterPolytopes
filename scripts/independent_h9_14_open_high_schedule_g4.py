#!/usr/bin/env python3
"""No-project-import replay of strict link-g4 cutoffs in seven open shadows."""

from __future__ import annotations

from hashlib import sha256
from math import comb
import json
from pathlib import Path
from typing import Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/enumerations/h9_14_nerve_support.json"
VERTEX_COUNT = 14
LINK_POLYTOPE_DIMENSION = 8
RECORD_GRAPHS = {
    0: "M?????????????o??",
    1: "M?????????????w??",
    3: "M???????????_?O??",
    4: "M???????????_?o??",
    5: "M???????????_?W??",
    7: "M???????????_?w??",
    9: "M???????????_?o?_",
}


def _summary_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("summary_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


def _summed_link_g4_from_h(f_vector: Sequence[int]) -> int:
    summed_link_f = {0: VERTEX_COUNT}
    for order in range(1, 5):
        summed_link_f[order] = (order + 1) * int(f_vector[order])

    def h(index: int) -> int:
        return sum(
            (-1) ** (index - order)
            * comb(
                LINK_POLYTOPE_DIMENSION - order,
                LINK_POLYTOPE_DIMENSION - index,
            )
            * summed_link_f[order]
            for order in range(index + 1)
        )

    return h(4) - h(3)


def _summed_link_g4_polynomial(f_vector: Sequence[int]) -> int:
    return (
        5 * int(f_vector[4])
        - 24 * int(f_vector[3])
        + 63 * int(f_vector[2])
        - 112 * int(f_vector[1])
        + 126 * VERTEX_COUNT
    )


def replay() -> dict[str, object]:
    raw = INPUT.read_bytes()
    document = json.loads(raw)
    by_graph = {row["graph6"]: row for row in document["records"]}
    rows = []
    for record, graph6 in RECORD_GRAPHS.items():
        source = by_graph[graph6]
        schedules = []
        for schedule in source["dehn_sommerville_schedule"]:
            f_vector = tuple(map(int, schedule["f_vector"]))
            direct = _summed_link_g4_from_h(f_vector)
            polynomial = _summed_link_g4_polynomial(f_vector)
            if direct != polynomial:
                raise AssertionError("independent link-g4 formulas disagree")
            schedules.append(
                {
                    "missing_four_count": int(schedule["missing_four_count"]),
                    "missing_five_count": int(
                        schedule["missing_five_count_offset"]
                    ),
                    "facet_count": f_vector[-1],
                    "summed_vertex_link_g4": direct,
                    "strict_sum_feasible": direct >= VERTEX_COUNT,
                }
            )
        excluded = [
            row["missing_four_count"]
            for row in schedules
            if not row["strict_sum_feasible"]
        ]
        rows.append(
            {
                "record_index": record,
                "graph6": graph6,
                "schedules": schedules,
                "strict_link_g4_excluded_schedules": excluded,
            }
        )
    result: dict[str, object] = {
        "implementation": (
            "Reads the frozen support archive directly, derives h3 and h4 "
            "from the eight-dimensional link h-vector formula, and separately "
            "checks the expanded global-f-vector polynomial."
        ),
        "input_file_sha256": sha256(raw).hexdigest(),
        "required_strict_sum": VERTEX_COUNT,
        "records": rows,
        "total_schedule_branches_excluded": sum(
            len(row["strict_link_g4_excluded_schedules"]) for row in rows
        ),
    }
    result["summary_sha256"] = _summary_digest(result)
    return result


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
