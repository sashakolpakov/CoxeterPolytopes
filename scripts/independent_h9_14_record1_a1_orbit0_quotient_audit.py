#!/usr/bin/env python3
"""Independent bounded audit of the Record-1 source-orbit-0 quotient engine.

This replay does not use the production C++ canonicalizer to derive its
orbits.  It reconstructs the candidate ground from the graph and local-link
Macaulay intervals, derives the exact S2 x S3 x S7 action as explicit ground
permutations, and compares the resulting equivalence relations with nauty.

The replay is intentionally bounded:

* all 1,278 singleton candidate rows are partitioned under the full group;
* all 1,277 depth-two families containing seed 31 are partitioned under the
  seed stabilizer; and
* a deterministic multirow pool is canonicalized by streaming all 60,480
  source-group elements.

It also performs a 64-node production cutoff run and verifies that neither
the class terminal nor the run terminal can be mistaken for completeness.
It never launches the full 21-class incidence census.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, islice, permutations
import json
from math import comb, factorial
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
CPP = ROOT / "scripts/explore_h9_14_record1_a1_orbit0_quotient.cpp"
ARCHIVE = (
    ROOT
    / "data/enumerations/h9_14_record1_a1_orbit0_quotient_bounded_audit.json"
)
N = 14
FULL = (1 << N) - 1
EDGES = ((0, 13), (1, 13), (2, 13))
EDGE_MASKS = tuple(sum(1 << value for value in edge) for edge in EDGES)
FOUR = sum(1 << value for value in (0, 3, 4, 5))
SOURCE_BLOCKS = ((0,), (13,), (1, 2), (3, 4, 5), tuple(range(6, 13)))
EXPECTED_GROUP_ORDER = factorial(2) * factorial(3) * factorial(7)

TARGET0 = (
    (1, 3, 4, 7, 8),
    (2, 3, 4, 7, 8),
    (3, 4, 5, 7, 8),
    (1, 3, 4, 8, 9),
    (0, 4, 5, 9, 10),
    (1, 4, 5, 9, 10),
    (0, 5, 6, 9, 10),
    (1, 2, 7, 8, 11),
    (1, 2, 8, 9, 11),
    (1, 4, 8, 9, 11),
    (1, 4, 9, 10, 11),
    (1, 6, 9, 10, 11),
    (0, 3, 5, 7, 12),
    (0, 2, 6, 7, 12),
    (2, 3, 7, 8, 12),
    (0, 5, 6, 10, 12),
    (2, 6, 7, 11, 12),
    (2, 6, 10, 11, 12),
    (6, 9, 10, 11, 12),
)


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def independent(value: int) -> bool:
    return all(value & edge != edge for edge in EDGE_MASKS)


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
        raise ArithmeticError("inexact independent Macaulay expansion")
    return answer


def local_interval(linked: int) -> tuple[int, int]:
    linked_size = linked.bit_count()
    dimension = 9 - linked_size
    top = 5 - linked_size
    rest = FULL ^ linked
    f = [0] * (top + 1)
    f[0] = 1
    extension = rest
    while True:
        size = extension.bit_count()
        if 1 <= size <= top:
            face = linked | extension
            f[size] += independent(face) and face & FOUR != FOUR
        if extension == 0:
            break
        extension = (extension - 1) & rest
    h = [
        sum(
            (-1) ** (index - order)
            * comb(dimension - order, dimension - index)
            * f[order]
            for order in range(index + 1)
        )
        for index in range(top + 1)
    ]
    g = [h[0], *(h[index] - h[index - 1] for index in range(1, top + 1))]
    previous = g[top - 1]
    before = g[top]
    lower = max(0, before - macaulay_raise(previous, top - 1))
    return lower, before - 1


def independent_candidates() -> tuple[int, ...]:
    intervals = {}
    for size in (1, 2, 3):
        for values in combinations(range(N), size):
            subset = mask(values)
            if independent(subset):
                intervals[subset] = local_interval(subset)
    if len(intervals) != 433 or any(left > right for left, right in intervals.values()):
        raise ArithmeticError("independent local-interval reconstruction changed")

    candidates = []
    shadow = 0
    independent_fives = 0
    for values in combinations(range(N), 5):
        five = mask(values)
        if not independent(five):
            continue
        independent_fives += 1
        if five & FOUR == FOUR:
            shadow += 1
            continue
        allowed = all(
            intervals[mask(subset)][1] >= 1
            for size in (1, 2, 3)
            for subset in combinations(values, size)
        )
        if allowed:
            candidates.append(five)
    if (independent_fives, shadow, len(candidates)) != (1497, 9, 1278):
        raise ArithmeticError("independent candidate reconstruction changed")
    return tuple(candidates)


def signature(value: int, blocks: Sequence[Sequence[int]]) -> tuple[int, ...]:
    return tuple(sum((value >> vertex) & 1 for vertex in block) for block in blocks)


def source_permutations() -> tuple[tuple[int, ...], ...]:
    answer = []
    for order_a in permutations((1, 2)):
        for order_b in permutations((3, 4, 5)):
            for order_c in permutations(range(6, 13)):
                action = list(range(N))
                for source, target in zip((1, 2), order_a, strict=True):
                    action[source] = target
                for source, target in zip((3, 4, 5), order_b, strict=True):
                    action[source] = target
                for source, target in zip(range(6, 13), order_c, strict=True):
                    action[source] = target
                answer.append(tuple(action))
    if len(answer) != EXPECTED_GROUP_ORDER or len(set(answer)) != EXPECTED_GROUP_ORDER:
        raise ArithmeticError("direct source group changed")
    return tuple(answer)


def moved(value: int, action: Sequence[int]) -> int:
    answer = 0
    for vertex in range(N):
        if value >> vertex & 1:
            answer |= 1 << action[vertex]
    return answer


def direct_family_keys(
    families: Sequence[tuple[int, ...]], actions: Sequence[Sequence[int]]
) -> tuple[tuple[int, ...], ...]:
    distinct_masks = sorted({value for family in families for value in family})
    image_table = {
        value: tuple(moved(value, action) for action in actions)
        for value in distinct_masks
    }
    answers = []
    for family in families:
        answers.append(
            min(
                tuple(sorted(image_table[value][index] for value in family))
                for index in range(len(actions))
            )
        )
    return tuple(answers)


def nauty_paths() -> tuple[Path, Path]:
    candidates = []
    configured = os.environ.get("NAUTY_PREFIX")
    if configured:
        candidates.append(Path(configured))
    candidates.extend((Path("/opt/homebrew"), Path("/usr/local")))
    for prefix in candidates:
        header = prefix / "include/nauty/nauty.h"
        library = prefix / "lib/libnauty.a"
        if header.is_file() and library.is_file():
            return header.parent, library
    raise RuntimeError("nauty headers/static library not found")


def compile_engine(binary: Path) -> None:
    compiler = shutil.which("clang++") or shutil.which("c++")
    if compiler is None:
        raise RuntimeError("a C++20 compiler is required")
    include, library = nauty_paths()
    subprocess.run(
        [
            compiler,
            "-std=c++20",
            "-O2",
            "-DNDEBUG",
            "-DRECORD1_VERIFY_PACKED_KEY",
            "-Wall",
            "-Wextra",
            "-Werror",
            f"-I{include}",
            str(CPP),
            str(library),
            "-o",
            str(binary),
        ],
        cwd=ROOT,
        check=True,
    )


def parse_jsonl(text: str) -> list[dict[str, object]]:
    return [json.loads(line) for line in text.splitlines() if line.strip()]


def run_engine(
    binary: Path, arguments: Sequence[str], stdin: str = ""
) -> list[dict[str, object]]:
    completed = subprocess.run(
        [str(binary), *arguments],
        input=stdin,
        text=True,
        capture_output=True,
        cwd=ROOT,
        check=True,
        timeout=180,
    )
    return parse_jsonl(completed.stdout)


def equivalence_bijection(
    direct: Sequence[object], production: Sequence[object]
) -> bool:
    if len(direct) != len(production):
        return False
    direct_to_production: dict[object, set[object]] = defaultdict(set)
    production_to_direct: dict[object, set[object]] = defaultdict(set)
    for left, right in zip(direct, production, strict=True):
        direct_to_production[left].add(right)
        production_to_direct[right].add(left)
    return all(len(values) == 1 for values in direct_to_production.values()) and all(
        len(values) == 1 for values in production_to_direct.values()
    )


def canonicalize(
    binary: Path,
    families: Sequence[tuple[int, ...]],
    seed_class: int | None = None,
) -> tuple[str, ...]:
    stdin = "".join(
        f"{len(family)} {' '.join(map(str, family))}\n" for family in families
    )
    arguments = (
        ["--canonicalize-full-stdin"]
        if seed_class is None
        else ["--canonicalize-seed-stdin", str(seed_class)]
    )
    events = run_engine(binary, arguments, stdin)
    terminal = events[-1]
    rows = events[:-1]
    if terminal != {
        "event": "canonical_terminal",
        "family_count": len(families),
        "terminal_status": "complete_eof",
        "proof_complete": True,
    }:
        raise ArithmeticError("production canonical stream did not close")
    if [row.get("row_index") for row in rows] != list(range(len(families))):
        raise ArithmeticError("production canonical stream changed order")
    return tuple(str(row["canonical_key"]) for row in rows)


def build_audit() -> dict[str, object]:
    candidates = independent_candidates()
    direct_singleton_groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for value in candidates:
        direct_singleton_groups[signature(value, SOURCE_BLOCKS)].append(value)
    direct_inventory = {
        min(values): len(values) for values in direct_singleton_groups.values()
    }
    if len(direct_inventory) != 21:
        raise ArithmeticError("direct singleton quotient changed")

    with tempfile.TemporaryDirectory(prefix="record1-a1-quotient-audit-") as raw:
        binary = Path(raw) / "engine"
        compile_engine(binary)
        inventory_events = run_engine(binary, ["--inventory"])
        header, inventory, terminal = inventory_events
        if not (
            header.get("protocol")
            == "h9_14_record1_a1_orbit0_incidence_quotient_jsonl_v1"
            and inventory.get("event") == "class_inventory"
            and terminal.get("terminal_status") == "inventory_complete"
            and terminal.get("whole_census_complete") is False
        ):
            raise ArithmeticError("production inventory protocol changed")
        production_inventory = {
            int(row["seed_mask"]): int(row["candidate_orbit_size"])
            for row in inventory["classes"]
        }
        if production_inventory != direct_inventory:
            raise ArithmeticError("direct and nauty candidate quotients differ")
        if sum(production_inventory.values()) != len(candidates):
            raise ArithmeticError("production candidate quotient is not a partition")

        seed = 31
        seed_row = next(row for row in inventory["classes"] if row["seed_mask"] == seed)
        seed_class = int(seed_row["class_index"])
        refined_seed_blocks = (
            (0,),
            (13,),
            (1, 2),
            (3, 4),
            (5,),
            tuple(range(6, 13)),
        )
        depth_two_families = tuple(
            (seed, value) for value in candidates if value != seed
        )
        direct_depth_two = tuple(
            signature(value, refined_seed_blocks)
            for value in candidates
            if value != seed
        )
        production_depth_two = canonicalize(
            binary, depth_two_families, seed_class=seed_class
        )
        if not equivalence_bijection(direct_depth_two, production_depth_two):
            raise ArithmeticError("depth-two seed-stabilizer quotients differ")

        representatives = tuple(sorted(direct_inventory))
        deterministic_pairs = tuple(
            islice(
                (
                    (candidates[left], candidates[right])
                    for left in range(0, len(candidates), 37)
                    for right in range(left + 1, len(candidates), 53)
                ),
                24,
            )
        )
        deterministic_triples = tuple(
            (
                candidates[index],
                candidates[(index + 137) % len(candidates)],
                candidates[(index + 503) % len(candidates)],
            )
            for index in range(0, 192, 24)
        )
        target_masks = tuple(mask(row) for row in TARGET0)
        target_prefixes = tuple(
            tuple(sorted(target_masks[:length])) for length in (2, 3, 5, 9, 19)
        )
        multirow_families = (
            tuple((value,) for value in representatives)
            + deterministic_pairs
            + deterministic_triples
            + target_prefixes
        )
        actions = source_permutations()
        direct_multirow = direct_family_keys(multirow_families, actions)
        production_multirow = canonicalize(binary, multirow_families)
        if not equivalence_bijection(direct_multirow, production_multirow):
            raise ArithmeticError("direct and nauty multirow quotients differ")

        cutoff_events = run_engine(
            binary,
            [
                "--class",
                str(seed_class),
                "--max-nodes",
                "64",
                "--checkpoint-nodes",
                "16",
            ],
        )
        checkpoints = [row for row in cutoff_events if row["event"] == "checkpoint"]
        class_terminal = next(
            row for row in cutoff_events if row["event"] == "class_terminal"
        )
        run_terminal = cutoff_events[-1]
        if not (
            [row["nodes"] for row in checkpoints] == [16, 32, 48, 64]
            and class_terminal["terminal_status"] == "cutoff_max_nodes"
            and class_terminal["class_proof_complete"] is False
            and run_terminal["terminal_status"] == "cutoff_max_nodes"
            and run_terminal["enumeration_complete"] is False
            and run_terminal["whole_census_complete"] is False
        ):
            raise ArithmeticError("cutoff protocol could be misread as complete")

    transcript_payload = json.dumps(
        {
            "inventory": inventory,
            "depth_two_keys": production_depth_two,
            "multirow_keys": production_multirow,
            "cutoff": cutoff_events,
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return {
        "schema": "h9_14_record1_a1_orbit0_quotient_bounded_audit_v1",
        "claim": (
            "An independently derived S2xS3xS7 action agrees exactly with "
            "the production nauty quotient on all candidate singletons, all "
            "depth-two families over seed 31, and a deterministic multirow "
            "pool; a bounded production cutoff cannot emit a completeness "
            "terminal."
        ),
        "source_group": {
            "abstract_structure": "S2 x S3 x S7",
            "explicit_permutation_count": EXPECTED_GROUP_ORDER,
            "candidate_count": len(candidates),
            "candidate_orbit_count": len(direct_inventory),
            "candidate_orbit_size_histogram": {
                str(size): count
                for size, count in sorted(
                    Counter(direct_inventory.values()).items()
                )
            },
        },
        "seed_stabilizer_depth_two_replay": {
            "seed_mask": seed,
            "seed_class_index": seed_class,
            "seed_orbit_size": production_inventory[seed],
            "seed_stabilizer_order": EXPECTED_GROUP_ORDER // production_inventory[seed],
            "literal_family_count": len(depth_two_families),
            "direct_quotient_count": len(set(direct_depth_two)),
            "production_quotient_count": len(set(production_depth_two)),
            "equivalence_relations_identical": True,
        },
        "full_group_multirow_replay": {
            "family_count": len(multirow_families),
            "direct_quotient_count": len(set(direct_multirow)),
            "production_quotient_count": len(set(production_multirow)),
            "equivalence_relations_identical": True,
            "includes_complete_known_19_row_target": True,
        },
        "packed_key_dense_adjacency_replay_enabled": True,
        "cutoff_protocol": {
            "node_limit": 64,
            "checkpoint_nodes": [16, 32, 48, 64],
            "class_terminal": "cutoff_max_nodes",
            "class_proof_complete": False,
            "run_terminal": "cutoff_max_nodes",
            "enumeration_complete": False,
            "whole_census_complete": False,
        },
        "engine_sha256": sha256(CPP.read_bytes()).hexdigest(),
        "bounded_transcript_sha256": sha256(transcript_payload).hexdigest(),
        "scope": (
            "Exact bounded quotient/canonicalization and terminal-protocol "
            "audit only.  No class was exhausted and no Record-1 source "
            "closure is claimed."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build_audit()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not ARCHIVE.is_file() or ARCHIVE.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {ARCHIVE}")
    else:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
