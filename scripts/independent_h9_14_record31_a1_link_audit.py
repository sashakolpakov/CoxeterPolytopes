#!/usr/bin/env python3
"""Independent exact strict face-link audit for Record 31 compact a=1.

This script deliberately imports no project module.  It reconstructs every
vertex-, edge-, and triangle-link prefix directly from the graph shadow and
the six one-A4 representatives frozen in the positive-A5 archive.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
from typing import Iterable, Mapping


N = 14
FULL = (1 << N) - 1
MISSING_EDGES = ((0, 11), (1, 12), (0, 13), (1, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in MISSING_EDGES)
EXPECTED_F_VECTOR = [14, 87, 319, 770, 1281, 1477, 1142, 540, 120]
EXPECTED_MISSING_FIVES = 15
EXPECTED_FACETS = 120

EXPECTED_FAILURES = (1,)
EXPECTED_WITNESS_COUNTS = (0, 4, 0, 0, 0, 0)


def mask(vertices: Iterable[int]) -> int:
    return sum(1 << int(vertex) for vertex in vertices)


def independent(face: int) -> bool:
    return all(face & edge != edge for edge in EDGE_MASKS)


def macaulay_raise(value: int, degree: int) -> int:
    """Return the exact degree-``degree`` Macaulay successor of ``value``."""

    remaining = value
    next_top: int | None = None
    answer = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (
            (next_top is None or top + 1 < next_top)
            and comb(top + 1, lower) <= remaining
        ):
            top += 1
        remaining -= comb(top, lower)
        answer += comb(top + 1, lower + 1)
        next_top = top
    if remaining:
        raise AssertionError("inexact Macaulay expansion")
    return answer


def link_prefix(linked: int, fours: tuple[int, ...]) -> tuple[list[int], list[int], list[int]]:
    linked_size = linked.bit_count()
    dimension = 9 - linked_size
    top = 5 - linked_size
    f = [0] * (top + 1)
    f[0] = 1
    rest = FULL ^ linked
    extension = rest
    while True:
        size = extension.bit_count()
        if 1 <= size <= top:
            face = linked | extension
            if independent(face) and all(face & four != four for four in fours):
                f[size] += 1
        if not extension:
            break
        extension = (extension - 1) & rest
    h = [
        sum(
            (-1 if (index - order) & 1 else 1)
            * comb(dimension - order, dimension - index)
            * f[order]
            for order in range(index + 1)
        )
        for index in range(top + 1)
    ]
    g = [h[0], *(h[index] - h[index - 1] for index in range(1, top + 1))]
    return f, h, g


def _digest(rows: object) -> str:
    return sha256(
        json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def audit(document: Mapping[str, object]) -> dict[str, object]:
    production = document.get("production")
    if not isinstance(production, Mapping):
        raise ValueError("positive-A4 archive changed")

    schedules = production.get("schedules")
    if not isinstance(schedules, list):
        raise ValueError("Record 31 schedule table changed")
    schedule_matches = [row for row in schedules if row.get("a") == 1]
    if len(schedule_matches) != 1:
        raise ValueError("Record 31 a=1 schedule multiplicity changed")
    schedule = schedule_matches[0]
    if (
        schedule.get("f_vector") != EXPECTED_F_VECTOR
        or schedule.get("missing_independent_fives") != EXPECTED_MISSING_FIVES
        or EXPECTED_F_VECTOR[-1] != EXPECTED_FACETS
    ):
        raise ValueError("Record 31 a=1 Dehn--Sommerville schedule changed")

    branches = production.get("branches")
    if not isinstance(branches, list):
        raise ValueError("positive-A4 branch table changed")
    matches = [row for row in branches if row.get("a") == 1]
    if len(matches) != 1 or len(matches[0].get("orbits", [])) != 6:
        raise ValueError("Record 31 a=1 positive-A4 orbit table changed")

    rows = []
    failures = []
    for orbit_index, raw in enumerate(matches[0]["orbits"]):
        fours = tuple(mask(face) for face in raw["representative"])
        if len(fours) != 1 or fours[0].bit_count() != 4:
            raise ValueError("Record 31 a=1 A4 family changed")
        witnesses = []
        constraints = []
        intervals: Counter[tuple[int, int]] = Counter()
        size_histogram: Counter[int] = Counter()
        for size in (1, 2, 3):
            for vertices in combinations(range(N), size):
                linked = mask(vertices)
                if not independent(linked):
                    continue
                f, h, g = link_prefix(linked, fours)
                top = 5 - size
                previous, before = g[top - 1], g[top]
                lower = max(0, before - macaulay_raise(previous, top - 1))
                upper = before - 1
                interval = (lower, upper)
                intervals[interval] += 1
                size_histogram[size] += 1
                constraint = {
                    "subset": list(vertices),
                    "link_dimension": 9 - size,
                    "f_prefix": f,
                    "h_prefix": h,
                    "g_prefix": g,
                    "additional_A5_incidence_interval": [lower, upper],
                }
                constraints.append(constraint)
                if lower > upper:
                    witnesses.append(constraint)
        if len(constraints) != 420 or size_histogram != {1: 14, 2: 87, 3: 319}:
            raise AssertionError("independent local constraint ground changed")
        if witnesses:
            failures.append(orbit_index)
        rows.append(
            {
                "orbit_index": orbit_index,
                "positive_A4_masks": list(fours),
                "constraint_count": len(constraints),
                "constraint_size_histogram": {
                    str(size): size_histogram[size] for size in (1, 2, 3)
                },
                "constraint_table_sha256": _digest(constraints),
                "additional_A5_incidence_interval_histogram": {
                    f"{lower},{upper}": count
                    for (lower, upper), count in sorted(intervals.items())
                },
                "strict_face_link_failure": bool(witnesses),
                "strict_face_link_witness_count": len(witnesses),
                "witnesses": witnesses,
            }
        )
    if tuple(failures) != EXPECTED_FAILURES:
        raise AssertionError("independent strict-link failure partition changed")
    if tuple(row["strict_face_link_witness_count"] for row in rows) != EXPECTED_WITNESS_COUNTS:
        raise AssertionError("independent strict-link witness multiplicities changed")
    return {
        "kind": "independent_h9_14_record31_a1_strict_face_link_audit",
        "compact_only": True,
        "schedule": {
            "a": 1,
            "missing_five_count": EXPECTED_MISSING_FIVES,
            "facet_count": EXPECTED_FACETS,
            "f_vector": EXPECTED_F_VECTOR,
        },
        "positive_A4_orbit_count": len(rows),
        "constraints_per_orbit": 420,
        "total_constraints": sum(row["constraint_count"] for row in rows),
        "strict_face_link_excluded_orbit_count": len(failures),
        "strict_face_link_excluded_orbit_indices": failures,
        "residual_orbit_indices": [index for index in range(6) if index not in failures],
        "orbits": rows,
        "scope": {
            "strict_top_link_g_positivity_used": True,
            "macaulay_transition_used": True,
            "finite_volume_noncompact_inference": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("positive_archive", type=Path)
    args = parser.parse_args()
    document = json.loads(args.positive_archive.read_text(encoding="utf-8"))
    print(json.dumps(audit(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
