#!/usr/bin/env python3
"""Exact ridge-aware Boolean scout for record-18 ``a=1/add6`` orbits.

This is a discovery/cross-check backend, not the terminal solver-independent
certificate.  It selects six additional missing five-sets, imposes every
local interval and order-at-most-four transversal constraint, tracks the
surviving order-five transversals, and requires the scheduled 108 facets and
literal degree-zero-or-two ridge equations.  Every Boolean model is finally
dualized directly.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from typing import Sequence

import z3

from coxeter_polytopes.combinatorics.h9_14_record18 import (
    EDGE_MASKS,
    TARGET_BASE_LINK_G3,
    TARGET_BASE_LINK_G4,
)
from scripts.explore_h9_14_record18_a2_add5_z3 import (
    FULL,
    N,
    hits_all,
    independent,
    macaulay_raise,
    mask,
    minimal_transversals,
    vertices,
)


ADDITION_COUNT = 6
SCHEDULED_FACETS = 108


def derive(a4_mask: int) -> dict[str, object]:
    if a4_mask.bit_count() != 4 or not independent(a4_mask):
        raise ValueError("invalid independent A4 mask")
    fives = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and independent(value)
    )
    shadow = tuple(value for value in fives if value & a4_mask == a4_mask)
    if len(shadow) != 8:
        raise AssertionError("this backend is only for shadow-8/add6 cases")
    lower = []
    upper = []
    for vertex in range(N):
        a4_incidence = int(bool(a4_mask >> vertex & 1))
        forced = sum(value >> vertex & 1 for value in shadow)
        previous = TARGET_BASE_LINK_G3[vertex] - a4_incidence
        before = TARGET_BASE_LINK_G4[vertex] + 6 * a4_incidence - forced
        lower.append(max(0, before - macaulay_raise(previous, 3)))
        upper.append(before - 1)
    candidates = tuple(
        value
        for value in fives
        if value not in shadow
        and all(not (value >> vertex & 1) or upper[vertex] >= 1 for vertex in range(N))
    )
    base = (*EDGE_MASKS, a4_mask, *shadow)
    small = tuple(
        value
        for value in range(1, FULL + 1)
        if value.bit_count() <= 4 and hits_all(value, base)
    )
    order_five = tuple(
        value
        for value in range(1, FULL + 1)
        if value.bit_count() == 5 and hits_all(value, base)
    )
    return {
        "shadow": shadow,
        "lower": tuple(lower),
        "upper": tuple(upper),
        "candidates": candidates,
        "base": base,
        "small": small,
        "order_five": order_five,
    }


def solve(a4_mask: int) -> dict[str, object]:
    derived = derive(a4_mask)
    candidates: tuple[int, ...] = derived["candidates"]  # type: ignore[assignment]
    small: tuple[int, ...] = derived["small"]  # type: ignore[assignment]
    order_five: tuple[int, ...] = derived["order_five"]  # type: ignore[assignment]
    lower: tuple[int, ...] = derived["lower"]  # type: ignore[assignment]
    upper: tuple[int, ...] = derived["upper"]  # type: ignore[assignment]
    base: tuple[int, ...] = derived["base"]  # type: ignore[assignment]

    solver = z3.Solver()
    chosen = tuple(z3.Bool(f"x_{index}") for index in range(len(candidates)))
    solver.add(z3.PbEq([(variable, 1) for variable in chosen], ADDITION_COUNT))
    for vertex in range(N):
        incident = [
            (chosen[index], 1)
            for index, value in enumerate(candidates)
            if value >> vertex & 1
        ]
        if incident:
            solver.add(z3.PbGe(incident, lower[vertex]))
            solver.add(z3.PbLe(incident, upper[vertex]))
        elif lower[vertex] > 0 or upper[vertex] < 0:
            solver.add(z3.BoolVal(False))
    for value in small:
        coverers = [
            chosen[index]
            for index, row in enumerate(candidates)
            if not row & value
        ]
        if not coverers:
            raise AssertionError("an order-at-most-four transversal is uncoverable")
        solver.add(z3.Or(coverers))

    survives = tuple(z3.Bool(f"y_{index}") for index in range(len(order_five)))
    for index, value in enumerate(order_five):
        killers = [
            chosen[row]
            for row, candidate in enumerate(candidates)
            if not candidate & value
        ]
        solver.add(survives[index] == z3.Not(z3.Or(killers)))
    solver.add(z3.PbEq([(variable, 1) for variable in survives], SCHEDULED_FACETS))

    ridge_constraint_count = 0
    for six_vertices in combinations(range(N), 6):
        six = mask(six_vertices)
        group = [
            survives[index]
            for index, value in enumerate(order_five)
            if value & six == value
        ]
        if not group:
            continue
        degree = z3.Sum([z3.If(variable, 1, 0) for variable in group])
        solver.add(z3.Or(degree == 0, degree == 2))
        ridge_constraint_count += 1

    model_count = 0
    pure_targets: list[list[list[int]]] = []
    blocker_histogram: Counter[str] = Counter()
    status = solver.check()
    while status == z3.sat:
        model = solver.model()
        indices = tuple(
            index for index, variable in enumerate(chosen) if z3.is_true(model[variable])
        )
        additions = tuple(candidates[index] for index in indices)
        blocker = minimal_transversals((*base, *additions))
        histogram = Counter(value.bit_count() for value in blocker)
        key = ",".join(f"{size}:{count}" for size, count in sorted(histogram.items()))
        blocker_histogram[key] += 1
        if histogram == {5: SCHEDULED_FACETS}:
            pure_targets.append([list(vertices(value)) for value in additions])
        model_count += 1
        solver.add(z3.Or([z3.Not(chosen[index]) for index in indices]))
        status = solver.check()
    if status != z3.unsat:
        raise RuntimeError(f"Z3 did not terminate exactly: {status}")

    return {
        "kind": "h9_14_record18_a1_add6_ridge_aware_z3_scout",
        "theorem_certificate": False,
        "independent_replay_complete": False,
        "compact_only": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "z3_version": z3.get_version_string(),
        "representative_A4_mask": a4_mask,
        "representative_A4": list(vertices(a4_mask)),
        "forced_A5_shadow_count": len(derived["shadow"]),
        "additional_A5_count": ADDITION_COUNT,
        "addition_lower": list(lower),
        "addition_upper": list(upper),
        "candidate_count": len(candidates),
        "base_small_transversal_count": len(small),
        "base_order_five_transversal_count": len(order_five),
        "ridge_constraint_count": ridge_constraint_count,
        "ridge_and_scheduled_model_count": model_count,
        "literal_blocker_histogram": dict(sorted(blocker_histogram.items())),
        "scheduled_pure_target_count": len(pure_targets),
        "scheduled_target_additions": pure_targets,
        "solver_terminal_status": "unsat",
        "scope": (
            "Exact ridge-aware compact incidence scout and literal blocker replay on "
            "every model. Solver-independent replay, topology, labels, Gram data, "
            "and finite-volume noncompact conclusions are outside this artifact."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("a4_mask", type=lambda value: int(value, 0))
    args = parser.parse_args()
    print(json.dumps(solve(args.a4_mask), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
