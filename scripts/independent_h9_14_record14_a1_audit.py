#!/usr/bin/env python3
"""Independent CP-SAT replay of Record 14's five ``a=1`` source orbits."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from math import factorial
from typing import Sequence

from ortools.sat.python import cp_model

from scripts.independent_h9_14_record14_a2_audit import (
    EDGE_ROWS,
    FULL,
    LINK_G3,
    LINK_G4,
    N,
    POSITIVE,
    _Collector,
    _hits,
    _independent,
    _mask,
    _minimal_blocker,
    _raise_macaulay,
    _vertices,
)


def _add_involution_lex_leader(
    model: cp_model.CpModel,
    chosen: Sequence[cp_model.IntVar],
    permutation: Sequence[int],
    tag: str,
) -> int:
    if sorted(permutation) != list(range(len(chosen))):
        raise AssertionError("independent candidate action is not a permutation")
    if any(permutation[permutation[position]] != position for position in range(len(chosen))):
        raise AssertionError("independent candidate action is not an involution")
    representatives = tuple(
        position
        for position, partner in enumerate(permutation)
        if position < partner
    )
    prefix = model.new_bool_var(f"{tag}_prefix_0")
    model.add(prefix == 1)
    for step, position in enumerate(representatives):
        partner = permutation[position]
        model.add(chosen[position] <= chosen[partner]).only_enforce_if(prefix)
        equal = model.new_bool_var(f"{tag}_equal_{step}")
        model.add(chosen[position] == chosen[partner]).only_enforce_if(equal)
        model.add(chosen[position] != chosen[partner]).only_enforce_if(equal.Not())
        following = model.new_bool_var(f"{tag}_prefix_{step + 1}")
        model.add(following <= prefix)
        model.add(following <= equal)
        model.add(following >= prefix + equal - 1)
        prefix = following
    return len(representatives)


def _move(face: int, action: Sequence[int]) -> int:
    return sum(1 << action[vertex] for vertex in _vertices(face))


def _swap(left: int, right: int) -> tuple[int, ...]:
    action = list(range(N))
    action[left], action[right] = action[right], action[left]
    return tuple(action)


def _source_actions(four: int) -> tuple[tuple[int, ...], ...]:
    inside = [vertex for vertex in range(3, 12) if four & (1 << vertex)]
    outside = [vertex for vertex in range(3, 12) if not four & (1 << vertex)]
    actions = [
        _swap(block[position], block[position + 1])
        for block in (inside, outside)
        for position in range(len(block) - 1)
    ]
    active_swap = _swap(0, 1)
    if _move(four, active_swap) == four:
        actions.append(active_swap)
    return tuple(actions)


def _build_a1_data(index: int) -> dict[str, object]:
    archive = json.loads(POSITIVE.read_text(encoding="utf-8"))
    branch = next(row for row in archive["branches"] if row["a"] == 1)
    source = branch["orbits"][index]["representative"]
    if len(source) != 1:
        raise AssertionError("independent a=1 source changed")
    fours = (_mask(source[0]),)
    five_ground = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and _independent(value)
    )
    forced = tuple(five for five in five_ground if five & fours[0] == fours[0])
    addition_count = 14 - len(forced)
    lower = []
    upper = []
    for vertex in range(N):
        four_degree = int(bool(fours[0] & (1 << vertex)))
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

    return {
        "source": source,
        "fours": fours,
        "forced": forced,
        "addition_count": addition_count,
        "lower": tuple(lower),
        "upper": tuple(upper),
        "candidates": candidates,
        "base": base,
        "small": small,
        "facets": facets,
    }


def replay_orbit(index: int) -> dict[str, object]:
    data = _build_a1_data(index)
    source = data["source"]
    fours = tuple(data["fours"])
    forced = tuple(data["forced"])
    addition_count = int(data["addition_count"])
    lower = tuple(data["lower"])
    upper = tuple(data["upper"])
    candidates = tuple(data["candidates"])
    base = tuple(data["base"])
    small = tuple(data["small"])
    facets = tuple(data["facets"])

    model = cp_model.CpModel()
    chosen = [model.new_bool_var(f"choose_{position}") for position in range(len(candidates))]
    live = [model.new_bool_var(f"live_{position}") for position in range(len(facets))]
    model.add(sum(chosen) == addition_count)
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
    model.add(sum(live) == 108)
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

    inside = [vertex for vertex in range(3, 12) if fours[0] & (1 << vertex)]
    outside = [vertex for vertex in range(3, 12) if not fours[0] & (1 << vertex)]
    actions = _source_actions(fours[0])
    active_swap = _swap(0, 1)
    source_stabilizer_order = (
        factorial(len(inside))
        * factorial(len(outside))
        * (2 if _move(fours[0], active_swap) == fours[0] else 1)
    )
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
        raise AssertionError(f"independent CP-SAT replay incomplete: {solver.status_name(status)}")
    histogram: Counter[str] = Counter()
    pure = []
    for family in sorted(collector.families):
        blocker = _minimal_blocker((*base, *family))
        sizes = Counter(value.bit_count() for value in blocker)
        key = ",".join(f"{size}:{count}" for size, count in sorted(sizes.items()))
        histogram[key] += 1
        if sizes == {5: 108}:
            pure.append(family)
    return {
        "kind": "independent_h9_14_record14_a1_stabilizer_generator_lex_slice_cp_sat_census",
        "schedule_a": 1,
        "source_positive_A4_orbit_index": index,
        "representative_A4s": source,
        "forced_A5_shadow_count": len(forced),
        "additional_A5_count": addition_count,
        "scheduled_facet_count": 108,
        "candidate_count": len(candidates),
        "base_small_transversal_count": len(small),
        "base_order_five_transversal_count": len(facets),
        "six_set_constraint_count": six_count,
        "stabilizer_generator_count": len(actions),
        "source_stabilizer_order": source_stabilizer_order,
        "generator_moved_candidate_pair_counts": moved_pair_counts,
        "lex_leader_statement": (
            "x is lexicographically no larger than every independently "
            "derived stabilizer-generator image; every finite group orbit "
            "has at least one such element"
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbit_index", type=int)
    args = parser.parse_args(argv)
    if not 0 <= args.orbit_index < 5:
        raise SystemExit("orbit index must lie in 0..4")
    print(json.dumps(replay_orbit(args.orbit_index), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
