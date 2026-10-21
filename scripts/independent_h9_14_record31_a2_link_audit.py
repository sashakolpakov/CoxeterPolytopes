#!/usr/bin/env python3
"""Independent exact strict face-link audit for Record 31 compact a=2."""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from math import comb
from pathlib import Path
from typing import Iterable, Mapping


N = 14
FULL = (1 << N) - 1
MISSING_EDGES = ((0, 11), (1, 12), (0, 13), (1, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in MISSING_EDGES)
EXPECTED_FAILURES = (
    4, 7, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31, 32, 36, 38, 39, 41, 43, 45, 46, 48,
)


def mask(vertices: Iterable[int]) -> int:
    return sum(1 << int(vertex) for vertex in vertices)


def independent(face: int) -> bool:
    return all(face & edge != edge for edge in EDGE_MASKS)


def macaulay_raise(value: int, degree: int) -> int:
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


def link_prefix(linked: int, fours: tuple[int, int]) -> tuple[list[int], list[int], list[int]]:
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


def audit(document: Mapping[str, object]) -> dict[str, object]:
    production = document.get("production")
    if not isinstance(production, Mapping):
        raise ValueError("positive-A4 archive changed")
    branches = production.get("branches")
    if not isinstance(branches, list):
        raise ValueError("positive-A4 branch table changed")
    matches = [row for row in branches if row.get("a") == 2]
    if len(matches) != 1 or len(matches[0].get("orbits", [])) != 49:
        raise ValueError("Record 31 a=2 positive-A4 orbit table changed")

    rows = []
    failures = []
    for orbit_index, raw in enumerate(matches[0]["orbits"]):
        fours = tuple(mask(face) for face in raw["representative"])
        if len(fours) != 2:
            raise ValueError("Record 31 a=2 A4 family changed")
        witnesses = []
        constraint_count = 0
        for size in (1, 2, 3):
            for vertices in combinations(range(N), size):
                linked = mask(vertices)
                if not independent(linked):
                    continue
                constraint_count += 1
                f, h, g = link_prefix(linked, fours)
                top = 5 - size
                previous, before = g[top - 1], g[top]
                lower = max(0, before - macaulay_raise(previous, top - 1))
                upper = before - 1
                if lower > upper:
                    witnesses.append(
                        {
                            "subset": list(vertices),
                            "link_dimension": 9 - size,
                            "f_prefix": f,
                            "h_prefix": h,
                            "g_prefix": g,
                            "additional_A5_incidence_interval": [lower, upper],
                        }
                    )
        if constraint_count != 420:
            raise AssertionError("independent local constraint ground changed")
        if witnesses:
            failures.append(orbit_index)
        rows.append(
            {
                "orbit_index": orbit_index,
                "positive_A4_masks": list(fours),
                "strict_face_link_failure": bool(witnesses),
                "witnesses": witnesses,
            }
        )
    if tuple(failures) != EXPECTED_FAILURES:
        raise AssertionError("independent strict-link failure partition changed")
    return {
        "kind": "independent_h9_14_record31_a2_strict_face_link_audit",
        "compact_only": True,
        "positive_A4_orbit_count": len(rows),
        "strict_face_link_excluded_orbit_count": len(failures),
        "strict_face_link_excluded_orbit_indices": failures,
        "residual_orbit_indices": [index for index in range(49) if index not in failures],
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
