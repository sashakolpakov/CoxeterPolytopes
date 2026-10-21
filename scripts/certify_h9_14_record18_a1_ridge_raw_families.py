#!/usr/bin/env python3
"""Freeze and verify the literal record-18 a=1/add6 target families.

The production source performs its own exact incidence/ridge search, quotients
only for speed with a collision-free byte-matrix key, expands every terminal
under the full graph/A4 stabilizer, and blocker-dualizes every expanded family.
This wrapper independently expands the canonical families in the standalone
128-bit audit and requires literal equality of the two normalized raw sets.
"""

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
    ROOT / "data/enumerations/h9_14_record18_a1_ridge_raw_families.json"
)
REPRESENTATIVES = (23, 15, 53)


Family = tuple[int, ...]


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


def _positive_a1_representatives() -> set[int]:
    document = json.loads(POSITIVE_A4.read_text(encoding="utf-8"))
    branches = [branch for branch in document["branches"] if branch["a"] == 1]
    if len(branches) != 1:
        raise AssertionError("positive archive lost its unique a=1 branch")
    representatives = {
        _mask(orbit["representative"][0]) for orbit in branches[0]["orbits"]
    }
    if any(len(orbit["representative"]) != 1 for orbit in branches[0]["orbits"]):
        raise AssertionError("positive a=1 orbit has the wrong row count")
    return representatives


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << int(vertex) for vertex in vertices)


def _normalized_families(raw: Sequence[Sequence[Sequence[int]]]) -> tuple[Family, ...]:
    answer = []
    for family in raw:
        masks = tuple(sorted(_mask(row) for row in family))
        if len(masks) != 6 or len(set(masks)) != 6:
            raise AssertionError("a raw a1 target is not six distinct rows")
        if any(mask.bit_count() != 5 for mask in masks):
            raise AssertionError("a raw a1 target contains a non-five-set")
        answer.append(masks)
    if len(answer) != len(set(answer)):
        raise AssertionError("raw a1 target list contains a duplicate family")
    return tuple(sorted(answer))


def _identity() -> tuple[int, ...]:
    return tuple(range(14))


def _compose(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(14))


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


def _swap(first: int, second: int) -> tuple[int, ...]:
    move = list(range(14))
    move[first], move[second] = move[second], move[first]
    return tuple(move)


def _exchange_leaf_pairs() -> tuple[int, ...]:
    move = list(range(14))
    for first, second in ((0, 2), (1, 3), (12, 13)):
        move[first], move[second] = move[second], move[first]
    return tuple(move)


def _move_mask(value: int, move: Sequence[int]) -> int:
    answer = 0
    for source in range(14):
        if value >> source & 1:
            answer |= 1 << move[source]
    return answer


def _expand_independent_a1(
    canonical: Sequence[Sequence[Sequence[int]]],
) -> tuple[Family, ...]:
    fixed_actions = _generated((_exchange_leaf_pairs(), _swap(4, 5)))
    if len(fixed_actions) != 4:
        raise AssertionError("rep53 small stabilizer is not order four")
    answer: set[Family] = set()
    for residual_images in permutations(range(6, 12)):
        for fixed in fixed_actions:
            move = list(fixed)
            for source, image in zip(range(6, 12), residual_images, strict=True):
                move[source] = image
            for family in canonical:
                rows = tuple(sorted(_move_mask(_mask(row), move) for row in family))
                answer.add(rows)
    return tuple(sorted(answer))


def _independent_raw_set() -> tuple[Family, ...]:
    document = json.loads(INDEPENDENT.read_text(encoding="utf-8"))
    cases = {
        tuple(case["representative_A4_masks"]): case
        for case in document["cases"]
    }
    case = cases[(53,)]
    if (
        case["canonical_target_count"] != 8
        or case["literal_target_orbit_weight_sum"] != 1080
    ):
        raise AssertionError("independent a1 target summary changed")
    expanded = _expand_independent_a1(case["canonical_target_additions"])
    if len(expanded) != 1080:
        raise AssertionError("independent a1 canonical expansion is not 1080")
    return expanded


def _run_production() -> tuple[list[dict[str, object]], tuple[Family, ...]]:
    with tempfile.TemporaryDirectory(prefix="record18-a1-raw-") as directory:
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
        raw53: tuple[Family, ...] = ()
        for representative in REPRESENTATIVES:
            completed = subprocess.run(
                [
                    str(binary),
                    "1",
                    str(representative),
                    "--symmetry-expand",
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            case = json.loads(completed.stdout)
            raw = _normalized_families(case.pop("expanded_literal_target_additions"))
            if representative == 53:
                raw53 = raw
            elif raw:
                raise AssertionError("an expected-unsat a1 representative has targets")
            cases.append(case)
    return cases, raw53


def _family_digest(families: Sequence[Family]) -> str:
    payload = json.dumps(families, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def fresh_record() -> dict[str, object]:
    if _positive_a1_representatives() != set(REPRESENTATIVES):
        raise AssertionError("the three a1 representatives are not exhaustive")
    cases, production = _run_production()
    independent = _independent_raw_set()
    if production != independent:
        missing = sorted(set(independent) - set(production))
        extra = sorted(set(production) - set(independent))
        raise AssertionError(
            f"literal a1 set mismatch: missing={len(missing)} extra={len(extra)}"
        )
    by_rep = {tuple(case["representative_A4_masks"]): case for case in cases}
    if (
        by_rep[(23,)]["canonical_literal_pure_target_count"] != 0
        or by_rep[(15,)]["canonical_literal_pure_target_count"] != 0
        or by_rep[(53,)]["canonical_literal_pure_target_count"] != 8
        or by_rep[(53,)]["expanded_literal_pure_target_count"] != 1080
    ):
        raise AssertionError("production a1 target transcript changed")
    return {
        "kind": "h9_14_record18_a1_literal_raw_family_audit",
        "backend": (
            "collision-free byte-matrix production quotient, full literal "
            "orbit expansion, and independent 128-bit quotient comparison"
        ),
        "production_cases": cases,
        "normalized_raw_target_count": len(production),
        "normalized_raw_target_families": [list(family) for family in production],
        "normalized_raw_target_sha256": _family_digest(production),
        "source_files": _source_files(),
        "production_independent_literal_set_equality": True,
        "a1_unsat_representatives": [23, 15],
        "a1_target_representative": 53,
        "a1_target_orbit_count": 8,
        "a1_target_literal_count": 1080,
        "compact_only": True,
        "finite_volume_noncompact_checked": False,
        "scope": (
            "Exact literal compact-incidence/ridge family audit for all three "
            "a=1 representatives; no Coxeter or ideal-vertex conclusion."
        ),
    }


def verify_record(document: dict[str, object]) -> None:
    if document.get("kind") != "h9_14_record18_a1_literal_raw_family_audit":
        raise AssertionError("wrong a1 raw-family certificate kind")
    families = tuple(
        tuple(int(value) for value in family)
        for family in document["normalized_raw_target_families"]  # type: ignore[index]
    )
    if len(families) != 1080 or tuple(sorted(set(families))) != families:
        raise AssertionError("frozen a1 literal family set is not normalized")
    if _family_digest(families) != document["normalized_raw_target_sha256"]:
        raise AssertionError("frozen a1 literal family payload changed")
    if families != _independent_raw_set():
        raise AssertionError("frozen a1 literal set differs from independent expansion")
    if document.get("production_independent_literal_set_equality") is not True:
        raise AssertionError("a1 literal set equality flag is not true")
    if _positive_a1_representatives() != set(REPRESENTATIVES):
        raise AssertionError("positive archive no longer gives the three a1 reps")
    if document.get("source_files") != _source_files():
        raise AssertionError("a1 source-file provenance changed")


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
