#!/usr/bin/env python3
"""Independent CP-SAT replay of Record 14's compact ``a=0`` schedule."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from math import factorial
from typing import Sequence

from ortools.sat.python import cp_model

from scripts.independent_h9_14_record14_a1_audit import (
    _add_involution_lex_leader,
)
from scripts.independent_h9_14_record14_a2_audit import (
    EDGE_ROWS,
    FULL,
    LINK_G3,
    LINK_G4,
    N,
    _Collector,
    _hits,
    _independent,
    _minimal_blocker,
    _raise_macaulay,
    _vertices,
)


def _build_a0_data() -> dict[str, object]:
    five_ground = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and _independent(value)
    )
    lower = tuple(
        max(0, LINK_G4[vertex] - _raise_macaulay(LINK_G3[vertex], 3))
        for vertex in range(N)
    )
    upper = tuple(LINK_G4[vertex] - 1 for vertex in range(N))
    candidates = tuple(
        face
        for face in five_ground
        if all(
            not (face & (1 << vertex)) or upper[vertex] > 0
            for vertex in range(N)
        )
    )
    base = tuple(EDGE_ROWS)
    transversals = tuple(
        value for value in range(1, FULL + 1) if _hits(value, base)
    )
    small = tuple(value for value in transversals if value.bit_count() <= 4)
    facets = tuple(value for value in transversals if value.bit_count() == 5)

    return {
        "lower": lower,
        "upper": upper,
        "candidates": candidates,
        "base": base,
        "small": small,
        "facets": facets,
    }


def _move(face: int, action: Sequence[int]) -> int:
    return sum(1 << action[vertex] for vertex in _vertices(face))


def _swap(left: int, right: int) -> tuple[int, ...]:
    action = list(range(N))
    action[left], action[right] = action[right], action[left]
    return tuple(action)


def _source_actions() -> tuple[tuple[int, ...], ...]:
    return tuple([_swap(left, left + 1) for left in range(3, 11)] + [_swap(0, 1)])


def replay_a0() -> dict[str, object]:
    data = _build_a0_data()
    lower = tuple(data["lower"])
    upper = tuple(data["upper"])
    candidates = tuple(data["candidates"])
    base = tuple(data["base"])
    small = tuple(data["small"])
    facets = tuple(data["facets"])
    model = cp_model.CpModel()
    chosen = [
        model.new_bool_var(f"choose_{position}")
        for position in range(len(candidates))
    ]
    live = [
        model.new_bool_var(f"live_{position}")
        for position in range(len(facets))
    ]
    model.add(sum(chosen) == 8)
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
        if destroyers:
            model.add_bool_or(destroyers)
        else:
            model.add(0 == 1)
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
    model.add(sum(live) == 110)
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

    actions = _source_actions()
    candidate_index = {face: position for position, face in enumerate(candidates)}
    moved_pair_counts = []
    for action_index, action in enumerate(actions):
        permutation = tuple(candidate_index[_move(face, action)] for face in candidates)
        moved_pair_counts.append(
            _add_involution_lex_leader(
                model, chosen, permutation, f"lex_{action_index}"
            )
        )

    solver = cp_model.CpSolver()
    solver.parameters.enumerate_all_solutions = True
    solver.parameters.num_search_workers = 1
    collector = _Collector(chosen, candidates)
    status = solver.solve(model, collector)
    if status not in (cp_model.OPTIMAL, cp_model.INFEASIBLE):
        raise AssertionError(
            f"independent a=0 CP-SAT replay incomplete: {solver.status_name(status)}"
        )
    histogram: Counter[str] = Counter()
    pure = []
    for family in sorted(collector.families):
        blocker = _minimal_blocker((*base, *family))
        sizes = Counter(value.bit_count() for value in blocker)
        key = ",".join(f"{size}:{count}" for size, count in sorted(sizes.items()))
        histogram[key] += 1
        if sizes == {5: 110}:
            pure.append(family)
    return {
        "kind": "independent_h9_14_record14_a0_full_source_generator_lex_slice_cp_sat_census",
        "schedule_a": 0,
        "source_branch_count": 1,
        "additional_A5_count": 8,
        "scheduled_facet_count": 110,
        "candidate_count": len(candidates),
        "base_small_transversal_count": len(small),
        "base_order_five_transversal_count": len(facets),
        "six_set_constraint_count": six_count,
        "stabilizer_generator_count": len(actions),
        "source_stabilizer_order": 2 * factorial(9),
        "generator_moved_candidate_pair_counts": moved_pair_counts,
        "lex_leader_statement": (
            "x is lexicographically no larger than every independently "
            "derived S9 x S2 generator image; every finite group orbit has "
            "at least one such element"
        ),
        "solver_status": solver.status_name(status),
        "generator_lex_slice_model_count": len(collector.families),
        "generator_lex_slice_pure_target_count": len(pure),
        "blocker_size_histogram": dict(sorted(histogram.items())),
        "generator_lex_slice_pure_targets": [
            [list(_vertices(face)) for face in family] for family in pure
        ],
        "enumeration_complete": True,
    }


def main(argv: Sequence[str] | None = None) -> None:
    argparse.ArgumentParser(description=__doc__).parse_args(argv)
    print(json.dumps(replay_a0(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
