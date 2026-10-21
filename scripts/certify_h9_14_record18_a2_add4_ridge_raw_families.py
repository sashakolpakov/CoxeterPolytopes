#!/usr/bin/env python3
"""Freeze and verify all twelve literal record-18 a=2/add4 target sets."""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/explore_h9_14_record18_ridge_aware_additions.cpp"
INDEPENDENT_SOURCE = (
    ROOT / "scripts/independent_h9_14_record18_a2_add5_audit.cpp"
)
INDEPENDENT = (
    ROOT / "data/enumerations/h9_14_record18_high_addition_independent_audit.json"
)
POSITIVE_A4 = (
    ROOT / "data/enumerations/h9_14_record18_positive_A4_scout.json"
)
CERTIFIER_SOURCE = Path(__file__).resolve()
DEFAULT_OUTPUT = (
    ROOT / "data/enumerations/h9_14_record18_a2_add4_ridge_raw_families.json"
)
REPRESENTATIVES = (
    (53, 197),
    (53, 86),
    (53, 198),
    (23, 101),
    (53, 58),
    (53, 90),
    (53, 202),
    (23, 57),
    (23, 105),
    (15, 53),
    (23, 43),
    (23, 45),
)
EXPECTED_RAW_COUNTS = {(53, 198): 6, (53, 202): 2}


Family = tuple[int, ...]
EDGE_MASKS = (
    (1 << 0) | (1 << 12),
    (1 << 1) | (1 << 12),
    (1 << 2) | (1 << 13),
    (1 << 3) | (1 << 13),
    (1 << 12) | (1 << 13),
)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_files() -> list[dict[str, str]]:
    return [
        {
            "path": str(path.relative_to(ROOT)),
            "sha256": _file_sha256(path),
        }
        for path in (
            SOURCE,
            INDEPENDENT_SOURCE,
            INDEPENDENT,
            POSITIVE_A4,
            CERTIFIER_SOURCE,
        )
    ]


def _positive_add4_representatives() -> tuple[tuple[int, int], ...]:
    document = json.loads(POSITIVE_A4.read_text(encoding="utf-8"))
    branches = [branch for branch in document["branches"] if branch["a"] == 2]
    if len(branches) != 1:
        raise AssertionError("positive archive lost its unique a=2 branch")
    representatives = []
    for orbit in branches[0]["orbits"]:
        fours = tuple(
            sorted(_mask(row) for row in orbit["representative"])
        )
        shadow_count = sum(
            any(value & four == four for four in fours)
            for value in range(1 << 14)
            if value.bit_count() == 5
            and not any(value & edge == edge for edge in EDGE_MASKS)
        )
        if shadow_count == 16:
            representatives.append(fours)
        elif shadow_count != 15:
            raise AssertionError("positive a=2 orbit has an unknown shadow size")
    return tuple(representatives)


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << int(vertex) for vertex in vertices)


def _identity() -> tuple[int, ...]:
    return tuple(range(14))


def _compose(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(14))


def _swap(first: int, second: int) -> tuple[int, ...]:
    move = list(range(14))
    move[first], move[second] = move[second], move[first]
    return tuple(move)


def _exchange_leaf_pairs() -> tuple[int, ...]:
    move = list(range(14))
    for first, second in ((0, 2), (1, 3), (12, 13)):
        move[first], move[second] = move[second], move[first]
    return tuple(move)


def _generated(generators: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    identity = _identity()
    found = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            for composite in (
                _compose(current, generator),
                _compose(generator, current),
            ):
                if composite not in found:
                    found.add(composite)
                    frontier.append(composite)
    return tuple(sorted(found))


ACTIVE = _generated((_swap(0, 1), _swap(2, 3), _exchange_leaf_pairs()))


def _move_mask(value: int, move: Sequence[int]) -> int:
    answer = 0
    for source in range(14):
        if value >> source & 1:
            answer |= 1 << move[source]
    return answer


def _small_stabilizer(fours: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    used = sorted(
        vertex
        for vertex in range(4, 12)
        if any(four >> vertex & 1 for four in fours)
    )
    if used and used != list(range(4, max(used) + 1)):
        raise AssertionError("A4 representative has noncontiguous used isolates")
    target = frozenset(fours)
    actions = set()
    for images in permutations(used):
        for active in ACTIVE:
            move = list(active)
            for source, image in zip(used, images, strict=True):
                move[source] = image
            if frozenset(_move_mask(four, move) for four in fours) == target:
                actions.add(tuple(move))
    if not actions:
        raise AssertionError("empty A4-pair stabilizer")
    return tuple(sorted(actions))


def _normalize_raw(raw: Sequence[Sequence[Sequence[int]]]) -> tuple[Family, ...]:
    answer = []
    for family in raw:
        rows = tuple(sorted(_mask(row) for row in family))
        if len(rows) != 4 or len(set(rows)) != 4:
            raise AssertionError("a2/add4 target is not four distinct rows")
        if any(row.bit_count() != 5 for row in rows):
            raise AssertionError("a2/add4 target contains a non-five-set")
        answer.append(rows)
    if len(answer) != len(set(answer)):
        raise AssertionError("duplicate literal a2/add4 family")
    return tuple(sorted(answer))


def _expand(
    fours: Sequence[int], canonical: Sequence[Sequence[Sequence[int]]]
) -> tuple[Family, ...]:
    used = [
        vertex
        for vertex in range(4, 12)
        if any(four >> vertex & 1 for four in fours)
    ]
    residual_start = max(used, default=3) + 1
    actions = _small_stabilizer(fours)
    answer: set[Family] = set()
    for residual_images in permutations(range(residual_start, 12)):
        for fixed in actions:
            move = list(fixed)
            for source, image in zip(
                range(residual_start, 12), residual_images, strict=True
            ):
                move[source] = image
            for family in canonical:
                rows = tuple(sorted(_move_mask(_mask(row), move) for row in family))
                answer.add(rows)
    return tuple(sorted(answer))


def _independent_sets() -> dict[tuple[int, int], tuple[Family, ...]]:
    document = json.loads(INDEPENDENT.read_text(encoding="utf-8"))
    cases = {
        tuple(case["representative_A4_masks"]): case
        for case in document["cases"]
    }
    answer = {}
    for representative in REPRESENTATIVES:
        case = cases[representative]
        expanded = _expand(representative, case["canonical_target_additions"])
        expected = EXPECTED_RAW_COUNTS.get(representative, 0)
        if len(expanded) != expected:
            raise AssertionError(
                f"independent raw weight changed at {representative}: {len(expanded)}"
            )
        answer[representative] = expanded
    return answer


def _run_production() -> tuple[list[dict[str, object]], dict[tuple[int, int], tuple[Family, ...]]]:
    with tempfile.TemporaryDirectory(prefix="record18-a2-add4-raw-") as directory:
        binary = Path(directory) / "ridge_audit"
        subprocess.run(
            [
                "c++",
                "-std=c++20",
                "-O3",
                "-Wall",
                "-Wextra",
                "-pedantic",
                str(SOURCE),
                "-o",
                str(binary),
            ],
            cwd=ROOT,
            check=True,
        )
        cases = []
        raw_sets = {}
        for representative in REPRESENTATIVES:
            completed = subprocess.run(
                [
                    str(binary),
                    "2",
                    *map(str, representative),
                    "--symmetry-expand",
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            case = json.loads(completed.stdout)
            raw_sets[representative] = _normalize_raw(
                case.pop("expanded_literal_target_additions")
            )
            cases.append(case)
    return cases, raw_sets


def _all_family_digest(
    raw_sets: dict[tuple[int, int], tuple[Family, ...]],
) -> str:
    payload = [
        [list(representative), [list(family) for family in raw_sets[representative]]]
        for representative in REPRESENTATIVES
    ]
    return sha256(json.dumps(payload, separators=(",", ":")).encode("ascii")).hexdigest()


def fresh_record() -> dict[str, object]:
    if _positive_add4_representatives() != REPRESENTATIVES:
        raise AssertionError("the twelve add4 representatives are not exhaustive")
    cases, production = _run_production()
    independent = _independent_sets()
    for representative in REPRESENTATIVES:
        if production[representative] != independent[representative]:
            missing = set(independent[representative]) - set(production[representative])
            extra = set(production[representative]) - set(independent[representative])
            raise AssertionError(
                f"a2/add4 set mismatch at {representative}: "
                f"missing={len(missing)} extra={len(extra)}"
            )
    return {
        "kind": "h9_14_record18_a2_add4_literal_raw_family_audit",
        "backend": (
            "collision-free byte-matrix production quotient, full literal "
            "orbit expansion, and independent 128-bit quotient comparison"
        ),
        "production_cases": cases,
        "normalized_raw_target_sets": [
            {
                "representative_A4_masks": list(representative),
                "target_count": len(production[representative]),
                "target_families": [
                    list(family) for family in production[representative]
                ],
            }
            for representative in REPRESENTATIVES
        ],
        "normalized_all_target_sha256": _all_family_digest(production),
        "source_files": _source_files(),
        "production_independent_literal_set_equality": True,
        "representative_count": 12,
        "zero_target_representative_count": 10,
        "target_representatives": [[53, 198], [53, 202]],
        "literal_target_count": 8,
        "compact_only": True,
        "finite_volume_noncompact_checked": False,
        "scope": (
            "Exact literal compact-incidence/ridge family audit for all twelve "
            "a=2/add4 representatives; no Coxeter or ideal-vertex conclusion."
        ),
    }


def verify_record(document: dict[str, object]) -> None:
    if document.get("kind") != "h9_14_record18_a2_add4_literal_raw_family_audit":
        raise AssertionError("wrong a2/add4 raw-family certificate kind")
    raw_sets = {}
    for entry in document["normalized_raw_target_sets"]:  # type: ignore[index]
        representative = tuple(int(value) for value in entry["representative_A4_masks"])
        families = tuple(tuple(int(row) for row in family) for family in entry["target_families"])
        if tuple(sorted(set(families))) != families:
            raise AssertionError(f"unnormalized frozen set at {representative}")
        raw_sets[representative] = families
    if tuple(raw_sets) != REPRESENTATIVES:
        raise AssertionError("frozen a2/add4 representative order changed")
    if _all_family_digest(raw_sets) != document["normalized_all_target_sha256"]:
        raise AssertionError("frozen a2/add4 literal family payload changed")
    if raw_sets != _independent_sets():
        raise AssertionError("frozen a2/add4 sets differ from independent expansion")
    if sum(map(len, raw_sets.values())) != 8:
        raise AssertionError("frozen a2/add4 literal target total is not eight")
    if _positive_add4_representatives() != REPRESENTATIVES:
        raise AssertionError("positive archive no longer gives the twelve add4 reps")
    if document.get("source_files") != _source_files():
        raise AssertionError("a2/add4 source-file provenance changed")
    if document.get("production_independent_literal_set_equality") is not True:
        raise AssertionError("a2/add4 literal set equality flag is not true")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    if args.fresh:
        document = fresh_record()
    else:
        document = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        verify_record(document)
    rendered = json.dumps(document, indent=2, sort_keys=True) + "\n"
    if args.write is not None:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
