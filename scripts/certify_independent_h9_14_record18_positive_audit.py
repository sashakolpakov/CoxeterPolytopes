#!/usr/bin/env python3
"""Freeze and replay the no-project record-18 positive-branch hostile audit.

Both C++ sources used here derive their data from the five missing edges.  The
production archive is read only after the independent positive census has run,
and only to compare complete canonical key/size sets.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from math import comb
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
POSITIVE_SOURCE = (
    ROOT / "scripts/independent_h9_14_record18_positive_orbit_audit.cpp"
)
EXCEPTION_SOURCE = (
    ROOT / "scripts/independent_h9_14_record18_a4_exception_audit.cpp"
)
PRODUCTION_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_record18_positive_A4_scout.json"
)
OUTPUT = (
    ROOT / "data/enumerations/h9_14_record18_positive_independent_audit.json"
)

N = 14
DIMENSION = 9
EDGES = (
    frozenset((0, 12)),
    frozenset((1, 12)),
    frozenset((2, 13)),
    frozenset((3, 13)),
    frozenset((12, 13)),
)
EXPECTED_BRANCHES = (
    (1, 495, 495, 145, 3),
    (2, 122_265, 120_969, 9_144, 14),
    (3, 20_092_215, 19_449_199, 188_944, 40),
    (4, 2_471_342_445, 2_310_730_275, 14_812, 13),
)
EXPECTED_SCHEDULES = (
    (0, 8, 110, 42),
    (1, 14, 108, 36),
    (2, 20, 106, 30),
    (3, 26, 104, 24),
    (4, 32, 102, 18),
    (5, 38, 100, 12),
    (6, 44, 98, 6),
)
EXPECTED_PAIR_CLASSES = {
    (252, 31_626, 3_501, 0): 3,
    (246, 30_135, 3_475, 0): 2,
    (126, 7_875, 1_375, 0): 5,
    (56, 1_540, 540, 0): 2,
}
EXPECTED_PAIR_CLASS_RECORD = {
    f"{candidate}/{raw}/{cap}/{no_small}": count
    for (candidate, raw, cap, no_small), count in sorted(EXPECTED_PAIR_CLASSES.items())
}
EXPECTED_EXCEPTION_HISTOGRAM = {
    "82": 160,
    "83": 1_920,
    "85": 240,
    "87": 120,
    "88": 900,
    "92": 260,
    "94": 360,
    "102": 0,
}


def _mask(vertices: Iterable[int]) -> int:
    value = 0
    for vertex in vertices:
        value |= 1 << vertex
    return value


EDGE_MASKS = tuple(_mask(edge) for edge in EDGES)


def _independent(value: int) -> bool:
    return all(value & edge != edge for edge in EDGE_MASKS)


def _independent_masks(order: int) -> tuple[int, ...]:
    return tuple(
        value
        for face in combinations(range(N), order)
        for value in (_mask(face),)
        if _independent(value)
    )


def _macaulay_raise(value: int, degree: int) -> int:
    remainder = value
    next_top: int | None = None
    raised = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (
            (next_top is None or top + 1 < next_top)
            and comb(top + 1, lower) <= remainder
        ):
            top += 1
        remainder -= comb(top, lower)
        raised += comb(top + 1, lower + 1)
        next_top = top
    if remainder:
        raise AssertionError("independent Macaulay expansion failed")
    return raised


def _h_prefix(a: int) -> tuple[int, ...]:
    partial_f = (1, 14, 86, 310, 735 - a)
    return tuple(
        sum(
            (-1) ** (index - order)
            * comb(DIMENSION - order, DIMENSION - index)
            * partial_f[order]
            for order in range(index + 1)
        )
        for index in range(5)
    )


def _summed_link_g4(f_vector: Sequence[int]) -> int:
    values = {0: N}
    for order in range(1, 5):
        values[order] = (order + 1) * int(f_vector[order])

    def h(index: int) -> int:
        return sum(
            (-1) ** (index - order)
            * comb(8 - order, 8 - index)
            * values[order]
            for order in range(index + 1)
        )

    return h(4) - h(3)


def _schedules() -> tuple[tuple[int, int, int, int], ...]:
    independent_fives = len(_independent_masks(5))
    rows = []
    for a in range(7):
        prefix = _h_prefix(a)
        if prefix[4] - prefix[3] < 1:
            break
        h_vector = prefix + tuple(reversed(prefix))
        f_vector = tuple(
            sum(
                comb(DIMENSION - order, DIMENSION - index) * h_vector[order]
                for order in range(index + 1)
            )
            for index in range(1, DIMENSION + 1)
        )
        rows.append(
            (
                a,
                independent_fives - f_vector[4],
                f_vector[-1],
                _summed_link_g4(f_vector),
            )
        )
    if tuple(rows) != EXPECTED_SCHEDULES:
        raise AssertionError("independent record-18 schedule derivation changed")
    return tuple(rows)


def _compile_and_run(source: Path) -> Mapping[str, object]:
    compiler = shutil.which("clang++") or shutil.which("g++")
    if compiler is None:
        raise RuntimeError("no C++20 compiler is available")
    with tempfile.TemporaryDirectory(prefix="record18-independent-") as directory:
        executable = Path(directory) / "audit"
        subprocess.run(
            [
                compiler,
                "-std=c++20",
                "-O3",
                "-Wall",
                "-Wextra",
                "-Wpedantic",
                str(source),
                "-o",
                str(executable),
            ],
            cwd=ROOT,
            check=True,
        )
        completed = subprocess.run(
            [str(executable)],
            cwd=ROOT,
            check=True,
            text=True,
            capture_output=True,
        )
    value = json.loads(completed.stdout)
    if not isinstance(value, Mapping):
        raise ValueError(f"independent executable {source.name} returned malformed JSON")
    return value


def _branch_rows(document: Mapping[str, object]) -> tuple[Mapping[str, object], ...]:
    raw = document.get("branches")
    if not isinstance(raw, list) or not all(isinstance(row, Mapping) for row in raw):
        raise ValueError("independent positive branches are malformed")
    return tuple(raw)  # type: ignore[return-value]


def _key_size_set(row: Mapping[str, object], field: str) -> set[tuple[tuple[int, ...], int]]:
    raw = row.get(field)
    if not isinstance(raw, list):
        raise ValueError(f"orbit field {field} is malformed")
    answer = set()
    for value in raw:
        if not isinstance(value, Mapping):
            raise ValueError(f"orbit field {field} contains a malformed row")
        size_field = "literal_size" if "literal_size" in value else "size"
        answer.add((tuple(map(int, value["key"])), int(value[size_field])))
    return answer


def _compare_positive_keys(
    independent: Mapping[str, object], production: Mapping[str, object]
) -> list[dict[str, object]]:
    replay_rows = _branch_rows(independent)
    production_rows = _branch_rows(production)
    summary = tuple(
        (
            int(row["a"]),
            int(row["raw_literal_families"]),
            int(row["interval_literal_families"]),
            int(row["viable_literal_families"]),
            int(row["viable_orbits"]),
        )
        for row in replay_rows
    )
    if summary != EXPECTED_BRANCHES:
        raise AssertionError("independent record-18 positive branch counts changed")
    answers = []
    for replay, archived in zip(replay_rows, production_rows, strict=True):
        if int(replay["a"]) != int(archived["a"]):
            raise AssertionError("positive branch ordering differs")
        replay_set = _key_size_set(replay, "viable_orbit_keys")
        archived_set = _key_size_set(archived, "orbits")
        if replay_set != archived_set:
            raise AssertionError(f"a={replay['a']} complete orbit key set differs")
        answers.append(
            {
                "a": int(replay["a"]),
                "orbit_count": len(replay_set),
                "producer_only": 0,
                "independent_only": 0,
                "complete_key_and_literal_size_set_equal": True,
            }
        )
    return answers


def _reconstruct_family(key: Sequence[int], a: int) -> tuple[int, ...]:
    active = tuple(map(int, key[:a]))
    columns = tuple(map(int, key[a:]))
    if len(columns) != 8:
        raise ValueError("canonical positive key does not have eight isolate columns")
    rows = []
    for row in range(a):
        value = active[row]
        for isolate, pattern in enumerate(columns, start=4):
            if (pattern >> row) & 1:
                value |= 1 << isolate
        if value.bit_count() != 4:
            raise AssertionError("canonical key did not reconstruct an A4")
        rows.append(value)
    if len(set(rows)) != a:
        raise AssertionError("canonical key reconstructed repeated A4 rows")
    return tuple(rows)


def _base_link_g() -> tuple[tuple[int, ...], tuple[int, ...]]:
    independent = {order: _independent_masks(order) for order in range(2, 6)}
    g3, g4 = [], []
    for vertex in range(N):
        counts = {
            order: sum(bool(face & (1 << vertex)) for face in faces)
            for order, faces in independent.items()
        }
        f0, f1, f2, f3 = (counts[order] for order in range(2, 6))
        g3.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    expected_g3 = (5,) * 4 + (6,) * 8 + (1,) * 2
    expected_g4 = (6,) * 4 + (7,) * 8 + (1,) * 2
    if tuple(g3) != expected_g3 or tuple(g4) != expected_g4:
        raise AssertionError("record-18 base link-g values changed")
    return tuple(g3), tuple(g4)


def _addition_intervals(
    family: Sequence[int], shadow: Sequence[int]
) -> tuple[tuple[int, int], ...]:
    base_g3, base_g4 = _base_link_g()
    answer = []
    for vertex in range(N):
        q = sum(bool(face & (1 << vertex)) for face in family)
        s = sum(bool(face & (1 << vertex)) for face in shadow)
        previous = base_g3[vertex] - q
        before = base_g4[vertex] + 6 * q - s
        answer.append(
            (max(0, before - _macaulay_raise(previous, 3)), before - 1)
        )
    return tuple(answer)


def _small_transversals(clutter: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        value
        for order in range(5)
        for vertices in combinations(range(N), order)
        for value in (_mask(vertices),)
        if all(value & edge for edge in clutter)
    )


def _pair_addition_audit(a4_row: Mapping[str, object]) -> dict[str, object]:
    key = tuple(map(int, a4_row["key"]))
    family = _reconstruct_family(key, 4)
    independent_fives = _independent_masks(5)
    shadow = tuple(
        five for five in independent_fives if any(five & four == four for four in family)
    )
    if len(shadow) != 30:
        raise AssertionError("ordinary a=4 orbit does not have shadow 30")
    intervals = _addition_intervals(family, shadow)
    shadow_set = set(shadow)
    candidates = tuple(
        five
        for five in independent_fives
        if five not in shadow_set
        and all(int(bool(five & (1 << vertex))) <= intervals[vertex][1] for vertex in range(N))
    )
    clutter = (*EDGE_MASKS, *family, *shadow)
    small = _small_transversals(clutter)
    cap_feasible = 0
    no_small = 0
    for left, right in combinations(candidates, 2):
        if not all(
            lower
            <= int(bool(left & (1 << vertex))) + int(bool(right & (1 << vertex)))
            <= upper
            for vertex, (lower, upper) in enumerate(intervals)
        ):
            continue
        cap_feasible += 1
        if not any((transversal & left) and (transversal & right) for transversal in small):
            no_small += 1
    return {
        "key": list(key),
        "positive_orbit_literal_size": int(a4_row["literal_size"]),
        "forced_shadow_count": len(shadow),
        "candidate_count": len(candidates),
        "raw_pair_count": comb(len(candidates), 2),
        "cap_feasible_pair_count": cap_feasible,
        "base_small_transversal_count": len(small),
        "no_small_transversal_pair_count": no_small,
    }


def _ordinary_a4_additions(
    independent: Mapping[str, object]
) -> tuple[list[dict[str, object]], dict[str, int]]:
    a4 = _branch_rows(independent)[3]
    raw_orbits = a4.get("viable_orbit_keys")
    if not isinstance(raw_orbits, list):
        raise ValueError("independent a=4 orbit list is malformed")
    ordinary = [
        orbit
        for orbit in raw_orbits
        if isinstance(orbit, Mapping) and int(orbit["additional_fives"]) == 2
    ]
    rows = [_pair_addition_audit(orbit) for orbit in ordinary]
    classes = Counter(
        (
            int(row["candidate_count"]),
            int(row["raw_pair_count"]),
            int(row["cap_feasible_pair_count"]),
            int(row["no_small_transversal_pair_count"]),
        )
        for row in rows
    )
    if len(rows) != 12 or dict(classes) != EXPECTED_PAIR_CLASSES:
        raise AssertionError("independent ordinary a=4 pair transcript changed")
    if any(row["no_small_transversal_pair_count"] != 0 for row in rows):
        raise AssertionError("ordinary a=4 orbit survived independent pair audit")
    return rows, {
        f"{candidate}/{raw}/{cap}/{no_small}": count
        for (candidate, raw, cap, no_small), count in sorted(classes.items())
    }


def build() -> dict[str, object]:
    schedules = _schedules()
    independent = _compile_and_run(POSITIVE_SOURCE)
    production = json.loads(PRODUCTION_ARCHIVE.read_text(encoding="utf-8"))
    if not isinstance(production, Mapping):
        raise ValueError("production positive archive is malformed")
    comparisons = _compare_positive_keys(independent, production)
    pair_rows, pair_classes = _ordinary_a4_additions(independent)
    exception = dict(_compile_and_run(EXCEPTION_SOURCE))
    if exception.get("no_small_transversal_addition_count") != 3_960 or (
        exception.get("facet_count_histogram") != EXPECTED_EXCEPTION_HISTOGRAM
        or exception.get("scheduled_102_facet_target_count") != 0
    ):
        raise AssertionError("independent exceptional a=4 transcript changed")

    branches = _branch_rows(independent)
    return {
        "schema_version": "1.0.0",
        "kind": "independent_h9_14_record18_positive_hostile_audit",
        "claim": (
            "A no-project incidence-matrix quotient independently reproduces "
            "all positive-A4 branch counts and every canonical orbit key/size; "
            "independent addition searches exclude the complete a=4 branch at "
            "a necessary compact-incidence gate."
        ),
        "compact_only": True,
        "theorem_certificate": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "remaining_unfinished_schedules": [0, 1, 2, 3],
        "schedules": [list(row) for row in schedules],
        "strict_compact_excluded_a": [5, 6],
        "strict_compact_surviving_a": [0, 1, 2, 3, 4],
        "positive_replay": {
            "backend": (
                "standalone C++20 incidence-matrix column-multiplicity quotient; "
                "no project headers, modules, or production tables"
            ),
            "active_automorphism_order": independent["active_automorphism_order"],
            "local_bound_count": independent["local_bound_count"],
            "independent_A4_count": independent["independent_A4_count"],
            "interval_A4_count": independent["interval_A4_count"],
            "branches": [
                {
                    "a": int(row["a"]),
                    "raw": int(row["raw_literal_families"]),
                    "interval_feasible": int(row["interval_literal_families"]),
                    "strict_aggregate_viable": int(row["viable_literal_families"]),
                    "orbit_count": int(row["viable_orbits"]),
                    "raw_orbit_count": int(row["raw_orbits"]),
                    "interval_orbit_count": int(row["interval_orbits"]),
                    "orbit_size_histogram": row["viable_orbit_size_histogram"],
                    "orbits": row["viable_orbit_keys"],
                }
                for row in branches
            ],
            "complete_production_key_set_comparison": comparisons,
        },
        "a4_ordinary_pair_addition_audit": {
            "backend": "standard-library Python; direct reverse-order pair replay",
            "orbit_count": len(pair_rows),
            "class_histogram": pair_classes,
            "orbits": pair_rows,
            "all_no_small_transversal_pair_counts_zero": True,
        },
        "a4_exception_audit": exception,
        "a4_schedule_incidence_empty": True,
        "scope": (
            "Necessary compact-incidence result only. The a=0,1,2,3 schedules "
            "remain open here; no whole-shadow Coxeter theorem and no finite-"
            "volume noncompact inference follows."
        ),
    }


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "independent_h9_14_record18_positive_hostile_audit":
        raise ValueError("wrong independent record-18 audit kind")
    if (
        document.get("compact_only") is not True
        or document.get("theorem_certificate") is not False
        or document.get("whole_shadow_compact_coxeter_theorem") is not False
        or document.get("finite_volume_noncompact_chamber_excluded") is not False
        or document.get("remaining_unfinished_schedules") != [0, 1, 2, 3]
    ):
        raise ValueError("independent record-18 audit scope was over-promoted")
    if tuple(map(tuple, document.get("schedules", ()))) != EXPECTED_SCHEDULES:
        raise ValueError("independent record-18 schedules changed")
    positive = document.get("positive_replay")
    if not isinstance(positive, Mapping):
        raise ValueError("independent record-18 positive replay is malformed")
    summary = tuple(
        (
            int(row["a"]),
            int(row["raw"]),
            int(row["interval_feasible"]),
            int(row["strict_aggregate_viable"]),
            int(row["orbit_count"]),
        )
        for row in positive["branches"]
    )
    if summary != EXPECTED_BRANCHES:
        raise ValueError("independent record-18 positive summary changed")
    if any(
        row.get("complete_key_and_literal_size_set_equal") is not True
        or row.get("producer_only") != 0
        or row.get("independent_only") != 0
        for row in positive["complete_production_key_set_comparison"]
    ):
        raise ValueError("independent record-18 orbit linkage changed")
    pairs = document.get("a4_ordinary_pair_addition_audit")
    if not isinstance(pairs, Mapping) or pairs.get("orbit_count") != 12 or (
        pairs.get("all_no_small_transversal_pair_counts_zero") is not True
        or pairs.get("class_histogram") != EXPECTED_PAIR_CLASS_RECORD
    ):
        raise ValueError("independent record-18 ordinary a=4 audit changed")
    exception = document.get("a4_exception_audit")
    if not isinstance(exception, Mapping) or (
        exception.get("no_small_transversal_addition_count") != 3_960
        or exception.get("facet_count_histogram") != EXPECTED_EXCEPTION_HISTOGRAM
        or exception.get("scheduled_102_facet_target_count") != 0
    ):
        raise ValueError("independent record-18 exceptional a=4 audit changed")
    if document.get("a4_schedule_incidence_empty") is not True:
        raise ValueError("independent record-18 a=4 conclusion changed")
    return {
        "positive_branch_counts_reproduced": True,
        "complete_orbit_key_sets_reproduced": True,
        "positive_orbit_count": 70,
        "a4_orbit_count": 13,
        "a4_schedule_incidence_empty": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "finite_volume_noncompact_chamber_excluded": False,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--fresh", action="store_true")
    args = parser.parse_args(argv)
    if args.write:
        document = build()
        OUTPUT.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    else:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    result = verify(document)
    if args.fresh and not args.write:
        rebuilt = build()
        if rebuilt != document:
            raise AssertionError("fresh independent record-18 audit differs from archive")
        result["fresh_no_project_replay_matches_archive"] = True
    print(json.dumps({"output": str(OUTPUT), **result}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
