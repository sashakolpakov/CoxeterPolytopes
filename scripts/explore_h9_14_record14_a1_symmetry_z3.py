#!/usr/bin/env python3
"""Lossless stabilizer-generator lex-slice census for Record 14 ``a=1``."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from math import factorial
import sys
from typing import Sequence

import z3

from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    FULL,
    N,
    build_problem,
    minimal_transversals,
    vertices,
)


def _move(value: int, action: Sequence[int]) -> int:
    return sum(1 << action[vertex] for vertex in vertices(value))


def _swap_action(left: int, right: int) -> tuple[int, ...]:
    action = list(range(N))
    action[left], action[right] = action[right], action[left]
    return tuple(action)


def stabilizer_generators(four: int) -> tuple[tuple[int, ...], ...]:
    inside = [vertex for vertex in range(3, 12) if four & (1 << vertex)]
    outside = [vertex for vertex in range(3, 12) if not four & (1 << vertex)]
    generators = [
        _swap_action(block[index], block[index + 1])
        for block in (inside, outside)
        for index in range(len(block) - 1)
    ]
    active_swap = _swap_action(0, 1)
    if _move(four, active_swap) == four:
        generators.append(active_swap)
    if not generators:
        raise AssertionError("Record-14 a=1 stabilizer has no generators")
    return tuple(generators)


def _add_lex_leader(
    solver: z3.Solver,
    selected: Sequence[z3.BoolRef],
    involution: Sequence[int],
    tag: str,
    context: z3.Context,
) -> int:
    """Impose ``selected <=lex involution(selected)`` losslessly.

    Every generator used here is a vertex transposition, so its induced
    candidate permutation is an involution.  A moved candidate pair
    ``i < j`` produces the same disagreement twice, first at ``i`` and later
    at ``j``.  Consequently only the smaller endpoint of each pair can be the
    first lexicographic disagreement.  Omitting fixed coordinates and larger
    endpoints is an exact simplification, not an additional assumption.
    """

    if sorted(involution) != list(range(len(selected))):
        raise AssertionError("candidate action is not a permutation")
    if any(involution[involution[index]] != index for index in range(len(selected))):
        raise AssertionError("candidate action is not an involution")
    representatives = tuple(
        index for index, partner in enumerate(involution) if index < partner
    )
    prefix = z3.Bool(f"{tag}_prefix_0", ctx=context)
    solver.add(prefix)
    for step, index in enumerate(representatives):
        left = selected[index]
        right = selected[involution[index]]
        solver.add(z3.Implies(prefix, z3.Or(z3.Not(left), right)))
        following = z3.Bool(f"{tag}_prefix_{step + 1}", ctx=context)
        solver.add(following == z3.And(prefix, left == right))
        prefix = following
    return len(representatives)


def enumerate_orbit(
    orbit_index: int,
    max_models: int | None = None,
    *,
    problem_override: dict[str, object] | None = None,
    generator_override: Sequence[Sequence[int]] | None = None,
    result_kind: str = "h9_14_record14_a1_stabilizer_generator_lex_slice_z3_census",
) -> dict[str, object]:
    problem = (
        build_problem(1, orbit_index)
        if problem_override is None
        else problem_override
    )
    fours = tuple(problem["fours"])
    if problem_override is None and len(fours) != 1:
        raise AssertionError("Record-14 a=1 source changed")
    candidates = tuple(problem["candidates"])
    small = tuple(problem["small"])
    facets = tuple(problem["fives"])
    lower = tuple(problem["lower"])
    upper = tuple(problem["upper"])
    scheduled_facets = int(problem["scheduled_facets"])
    context = z3.Context()
    selected = [z3.Bool(f"x_{index}", ctx=context) for index in range(len(candidates))]
    live = [z3.Bool(f"y_{index}", ctx=context) for index in range(len(facets))]
    solver = z3.Solver(ctx=context)
    solver.add(z3.PbEq([(variable, 1) for variable in selected], int(problem["addition_count"])))
    for vertex in range(N):
        terms = [
            (selected[index], 1)
            for index, face in enumerate(candidates)
            if face & (1 << vertex)
        ]
        if terms:
            solver.add(z3.PbGe(terms, lower[vertex]))
            solver.add(z3.PbLe(terms, upper[vertex]))
        else:
            solver.add(z3.BoolVal(lower[vertex] <= 0 <= upper[vertex], ctx=context))
    for transversal in small:
        destroyers = [
            selected[index]
            for index, face in enumerate(candidates)
            if not face & transversal
        ]
        solver.add(z3.Or(destroyers) if destroyers else z3.BoolVal(False, ctx=context))
    for facet_index, facet in enumerate(facets):
        destroyers = [
            selected[index]
            for index, face in enumerate(candidates)
            if not face & facet
        ]
        solver.add(
            live[facet_index]
            == (
                z3.And([z3.Not(value) for value in destroyers])
                if destroyers
                else z3.BoolVal(True, ctx=context)
            )
        )
    solver.add(z3.PbEq([(variable, 1) for variable in live], scheduled_facets))
    six_count = 0
    for six in range(FULL + 1):
        if six.bit_count() != 6:
            continue
        members = [
            (live[index], 1)
            for index, facet in enumerate(facets)
            if facet & six == facet
        ]
        if members:
            solver.add(z3.Or(z3.PbEq(members, 0), z3.PbEq(members, 2)))
            six_count += 1

    candidate_index = {face: index for index, face in enumerate(candidates)}
    generators = (
        stabilizer_generators(fours[0])
        if generator_override is None
        else tuple(tuple(action) for action in generator_override)
    )
    isolate_inside = sum(
        bool(fours[0] & (1 << vertex)) for vertex in range(3, 12)
    ) if fours else 0
    active_factor = (
        2
        if not fours
        or bool(fours[0] & (1 << 0)) == bool(fours[0] & (1 << 1))
        else 1
    )
    source_stabilizer_order = (
        factorial(isolate_inside)
        * factorial(9 - isolate_inside)
        * active_factor
    )
    generator_permutations = []
    generator_moved_pair_counts = []
    for generator_index, action in enumerate(generators):
        permutation = tuple(candidate_index[_move(face, action)] for face in candidates)
        generator_moved_pair_counts.append(
            _add_lex_leader(
                solver, selected, permutation, f"g{generator_index}", context
            )
        )
        generator_permutations.append(permutation)

    models = 0
    targets: set[tuple[int, ...]] = set()
    histogram: Counter[str] = Counter()
    complete = True
    terminal_status = solver.check()
    while terminal_status == z3.sat:
        model = solver.model()
        indices = tuple(
            index for index, variable in enumerate(selected) if z3.is_true(model.eval(variable))
        )
        family = tuple(sorted(candidates[index] for index in indices))
        blocker = minimal_transversals((*problem["base"], *family))
        sizes = Counter(value.bit_count() for value in blocker)
        key = ",".join(f"{size}:{count}" for size, count in sorted(sizes.items()))
        histogram[key] += 1
        if sizes == {5: scheduled_facets}:
            targets.add(family)
        models += 1
        solver.add(z3.Or([z3.Not(selected[index]) for index in indices]))
        if models % 100 == 0:
            print(f"models {models} targets {len(targets)}", file=sys.stderr, flush=True)
        if max_models is not None and models >= max_models:
            complete = False
            break
        terminal_status = solver.check()
    if complete and terminal_status != z3.unsat:
        raise RuntimeError(
            "Z3 generator-lex-slice census incomplete: "
            f"{terminal_status}; reason={solver.reason_unknown()}"
        )
    return {
        "kind": result_kind,
        "schedule_a": problem["a"],
        "source_positive_A4_orbit_index": orbit_index,
        "representative_A4s": [list(vertices(face)) for face in fours],
        "forced_A5_shadow_count": len(problem["shadow"]),
        "additional_A5_count": problem["addition_count"],
        "candidate_count": len(candidates),
        "base_small_transversal_count": len(small),
        "base_order_five_transversal_count": len(facets),
        "scheduled_facet_count": scheduled_facets,
        "six_set_constraint_count": six_count,
        "stabilizer_generator_count": len(generators),
        "source_stabilizer_order": source_stabilizer_order,
        "generator_moved_candidate_pair_counts": generator_moved_pair_counts,
        "generator_candidate_permutations_verified": all(
            sorted(permutation) == list(range(len(candidates)))
            for permutation in generator_permutations
        ),
        "lex_leader_statement": (
            "x is lexicographically no larger than every stabilizer-generator "
            "image; every finite group orbit has at least one such element"
        ),
        "generator_lex_slice_model_count": models,
        "generator_lex_slice_pure_target_count": len(targets),
        "generator_lex_slice_pure_targets": [
            [list(vertices(face)) for face in family] for family in sorted(targets)
        ],
        "blocker_size_histogram": dict(sorted(histogram.items())),
        "enumeration_complete": complete,
        "terminal_solver_status": (
            "cutoff"
            if not complete
            else str(terminal_status)
        ),
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbit_index", type=int)
    parser.add_argument("--max-models", type=int)
    args = parser.parse_args(argv)
    if not 0 <= args.orbit_index < 5:
        raise SystemExit("orbit index must lie in 0..4")
    print(json.dumps(enumerate_orbit(args.orbit_index, args.max_models), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
