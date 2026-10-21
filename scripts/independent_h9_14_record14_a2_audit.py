#!/usr/bin/env python3
"""Independent CP-SAT replay of the complete Record-14 ``a=2`` census.

This file deliberately imports no production Record-14 completion engine.
It rebuilds the graph, local arithmetic, source representatives, Boolean
incidence model, ridge equations, and blocker verification directly.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable, Sequence

from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[1]
POSITIVE = ROOT / "data/enumerations/h9_14_record14_positive_A4_scout.json"
N = 14
FULL = (1 << N) - 1
EDGE_ROWS = (
    (1 << 0) | (1 << 12),
    (1 << 1) | (1 << 12),
    (1 << 0) | (1 << 13),
    (1 << 1) | (1 << 13),
    (1 << 2) | (1 << 13),
)
LINK_G3 = (4, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 1)
LINK_G4 = (5, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 2, 1)


def _mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def _vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value >> vertex & 1)


def _independent(value: int) -> bool:
    return not any(value & edge == edge for edge in EDGE_ROWS)


def _hits(value: int, rows: Sequence[int]) -> bool:
    return all(value & row for row in rows)


def _choose(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    k = min(k, n - k)
    answer = 1
    for index in range(1, k + 1):
        answer = answer * (n - k + index) // index
    return answer


def _raise_macaulay(value: int, degree: int) -> int:
    rest = value
    ceiling: int | None = None
    result = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (ceiling is None or top + 1 < ceiling) and _choose(top + 1, lower) <= rest:
            top += 1
        rest -= _choose(top, lower)
        result += _choose(top + 1, lower + 1)
        ceiling = top
    if rest:
        raise AssertionError("independent Macaulay expansion failed")
    return result


def _minimal_blocker(rows: Sequence[int]) -> tuple[int, ...]:
    answers = []
    for value in range(1, FULL + 1):
        if not _hits(value, rows):
            continue
        minimal = True
        for vertex in _vertices(value):
            if _hits(value ^ (1 << vertex), rows):
                minimal = False
                break
        if minimal:
            answers.append(value)
    return tuple(answers)


def _source_representative(index: int) -> tuple[int, int]:
    archive = json.loads(POSITIVE.read_text(encoding="utf-8"))
    branch = next(row for row in archive["branches"] if row["a"] == 2)
    faces = branch["orbits"][index]["representative"]
    if len(faces) != 2:
        raise AssertionError("independent a=2 source changed")
    return tuple(_mask(face) for face in faces)  # type: ignore[return-value]


class _Collector(cp_model.CpSolverSolutionCallback):
    def __init__(self, variables: Sequence[cp_model.IntVar], candidates: Sequence[int]):
        super().__init__()
        self._variables = tuple(variables)
        self._candidates = tuple(candidates)
        self.families: set[tuple[int, ...]] = set()

    def on_solution_callback(self) -> None:
        family = tuple(
            sorted(
                face
                for variable, face in zip(self._variables, self._candidates, strict=True)
                if self.value(variable)
            )
        )
        self.families.add(family)


def replay_orbit(index: int) -> dict[str, object]:
    fours = _source_representative(index)
    five_ground = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and _independent(value)
    )
    forced = tuple(
        five for five in five_ground if any(five & four == four for four in fours)
    )
    additions = 20 - len(forced)
    lower = []
    upper = []
    for vertex in range(N):
        four_degree = sum(bool(face & (1 << vertex)) for face in fours)
        forced_degree = sum(bool(face & (1 << vertex)) for face in forced)
        previous = LINK_G3[vertex] - four_degree
        before = LINK_G4[vertex] + 6 * four_degree - forced_degree
        lower.append(max(0, before - _raise_macaulay(previous, 3)))
        upper.append(before - 1)
    candidates = tuple(
        face
        for face in five_ground
        if face not in forced
        and all(not (face & (1 << vertex)) or upper[vertex] > 0 for vertex in range(N))
    )
    base = (*EDGE_ROWS, *fours, *forced)
    transversals = tuple(value for value in range(1, FULL + 1) if _hits(value, base))
    small = tuple(value for value in transversals if value.bit_count() <= 4)
    facets = tuple(value for value in transversals if value.bit_count() == 5)

    model = cp_model.CpModel()
    chosen = [model.new_bool_var(f"choose_{position}") for position in range(len(candidates))]
    live = [model.new_bool_var(f"live_{position}") for position in range(len(facets))]
    model.add(sum(chosen) == additions)
    for vertex in range(N):
        incidence = sum(
            chosen[position]
            for position, face in enumerate(candidates)
            if face & (1 << vertex)
        )
        model.add(incidence >= lower[vertex])
        model.add(incidence <= upper[vertex])
    for transversal in small:
        destroyers = [
            chosen[position]
            for position, face in enumerate(candidates)
            if not face & transversal
        ]
        if not destroyers:
            model.add(0 == 1)
        else:
            model.add_bool_or(destroyers)
    for facet_index, facet in enumerate(facets):
        destroyers = [
            chosen[position]
            for position, face in enumerate(candidates)
            if not face & facet
        ]
        if destroyers:
            for destroyer in destroyers:
                model.add(destroyer == 0).only_enforce_if(live[facet_index])
            model.add_bool_or([live[facet_index], *destroyers])
        else:
            model.add(live[facet_index] == 1)
    model.add(sum(live) == 106)
    six_count = 0
    for six in range(FULL + 1):
        if six.bit_count() != 6:
            continue
        members = [
            live[position]
            for position, facet in enumerate(facets)
            if facet & six == facet
        ]
        if not members:
            continue
        degree = model.new_int_var(0, len(members), f"ridge_{six}")
        model.add(degree == sum(members))
        model.add_allowed_assignments([degree], [(0,), (2,)])
        six_count += 1

    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.num_search_workers = 1
    collector = _Collector(chosen, candidates)
    status = solver.solve(model, collector)
    if status not in (cp_model.OPTIMAL, cp_model.INFEASIBLE):
        raise AssertionError(f"independent CP-SAT replay incomplete: {solver.status_name(status)}")
    pure = []
    hist: Counter[str] = Counter()
    for family in sorted(collector.families):
        blocker = _minimal_blocker((*base, *family))
        sizes = Counter(value.bit_count() for value in blocker)
        key = ",".join(f"{size}:{count}" for size, count in sorted(sizes.items()))
        hist[key] += 1
        if sizes == {5: 106}:
            pure.append(family)
    if len(pure) != len(collector.families):
        raise AssertionError("independent model escaped literal blocker purity")
    return {
        "source_positive_A4_orbit_index": index,
        "representative_A4s": [list(_vertices(face)) for face in fours],
        "forced_A5_shadow_count": len(forced),
        "additional_A5_count": additions,
        "candidate_count": len(candidates),
        "base_small_transversal_count": len(small),
        "base_order_five_transversal_count": len(facets),
        "six_set_constraint_count": six_count,
        "solver_status": solver.status_name(status),
        "literal_model_count": len(collector.families),
        "blocker_size_histogram": dict(sorted(hist.items())),
        "literal_pure_targets": [
            [list(_vertices(face)) for face in family] for family in pure
        ],
        "enumeration_complete": True,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbit_index", type=int)
    args = parser.parse_args(argv)
    if not 0 <= args.orbit_index < 36:
        raise SystemExit("orbit index must lie in 0..35")
    print(json.dumps(replay_orbit(args.orbit_index), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
