#!/usr/bin/env python3
"""Compare producer and independent record-32 a=0 targets by true symmetry.

The two streams intentionally choose different literal representatives of an
``S8`` isolate orbit, so their raw representative counts need not agree.  This
no-project-import comparator canonically exhausts the genuine
``Aut(P4 disjoint union K2) x S8 x S6`` action and compares exact orbit-key
sets.  The independent C++ target file contains one comma-separated six-mask
family per line.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
from typing import Iterable, Sequence


ACTIVE = (0, 1, 2, 11, 12, 13)
ISOLATES = tuple(range(3, 11))
ACTIVE_MASK = sum(1 << vertex for vertex in ACTIVE)
EDGES = {
    frozenset(edge)
    for edge in ((0, 11), (1, 12), (0, 13), (2, 13))
}
EXPECTED_ORBIT_COUNT = 292
EXPECTED_ORBIT_KEY_SHA256 = (
    "97f630bd464ec01d77087019c5dd02b549f0cbaebd466d5e621ae065797dc848"
)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(14) if value & (1 << vertex))


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


@lru_cache(maxsize=1)
def active_automorphisms() -> tuple[dict[int, int], ...]:
    answers = []
    for image in permutations(ACTIVE):
        action = dict(zip(ACTIVE, image, strict=True))
        if {
            frozenset((action[left], action[right]))
            for left, right in map(tuple, EDGES)
        } == EDGES:
            answers.append(action)
    if len(answers) != 4:
        raise AssertionError("record-32 active automorphism group changed")
    return tuple(answers)


def orbit_key(family: Sequence[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    if len(family) != 6:
        raise ValueError("record-32 a=0 family does not have six rows")
    answer = None
    for action in active_automorphisms():
        transformed = tuple(
            mask(action.get(vertex, vertex) for vertex in vertices(face))
            for face in family
        )
        for row_order in permutations(range(6)):
            active_parts = tuple(
                transformed[index] & ACTIVE_MASK for index in row_order
            )
            isolate_columns = tuple(
                sorted(
                    sum(
                        bool(transformed[row_order[row]] & (1 << isolate)) << row
                        for row in range(6)
                    )
                    for isolate in ISOLATES
                )
            )
            key = (active_parts, isolate_columns)
            if answer is None or key < answer:
                answer = key
    if answer is None:
        raise AssertionError("record-32 orbit key escaped its action")
    return answer


def key_digest(keys: Iterable[tuple[tuple[int, ...], tuple[int, ...]]]) -> str:
    digest = sha256()
    for active_parts, isolate_columns in sorted(keys):
        for value in active_parts:
            digest.update(value.to_bytes(2, "little"))
        digest.update(bytes(isolate_columns))
    return digest.hexdigest()


def compare(production_path: Path, independent_path: Path) -> dict[str, object]:
    document = json.loads(production_path.read_text(encoding="utf-8"))
    production = tuple(
        tuple(map(int, family)) for family in document["target_family_masks"]
    )
    independent = tuple(
        tuple(map(int, line.split(",")))
        for raw in independent_path.read_text(encoding="utf-8").splitlines()
        if (line := raw.strip())
    )
    production_keys = {orbit_key(family) for family in production}
    independent_keys = {orbit_key(family) for family in independent}
    digest = key_digest(production_keys)
    result = {
        "production_literal_representative_count": len(production),
        "independent_literal_representative_count": len(independent),
        "production_orbit_count": len(production_keys),
        "independent_orbit_count": len(independent_keys),
        "production_only_orbit_count": len(production_keys - independent_keys),
        "independent_only_orbit_count": len(independent_keys - production_keys),
        "orbit_key_sha256": digest,
        "target_orbit_sets_equal": production_keys == independent_keys,
    }
    if (
        result["production_orbit_count"] != EXPECTED_ORBIT_COUNT
        or result["independent_orbit_count"] != EXPECTED_ORBIT_COUNT
        or not result["target_orbit_sets_equal"]
        or EXPECTED_ORBIT_KEY_SHA256 not in ("TO_BE_FROZEN", digest)
    ):
        raise AssertionError("record-32 independent a=0 target orbit replay changed")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("production", type=Path)
    parser.add_argument("independent", type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(args.production, args.independent), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
