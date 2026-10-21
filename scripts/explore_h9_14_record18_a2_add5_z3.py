#!/usr/bin/env python3
"""Exact ridge-aware Boolean scout for record-18 ``a=2/add5`` orbits.

This is an independent formulation of the two slow cover-DFS cases.  It
selects five additional minimal A5s, imposes the literal local degree caps,
eliminates every base transversal of order at most four, tracks every
surviving base transversal of order five, requires the scheduled count 106,
and imposes degree zero or two on every six-subset ridge group.  Every model
is then blocker-dualized directly; model blocking continues until UNSAT.

The result is a compact-incidence scout only.  It does not address ideal
vertices or finite-volume noncompact Coxeter chambers.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from typing import Iterable, Sequence

import z3

from coxeter_polytopes.combinatorics.h9_14_record18 import (
    EDGE_MASKS,
    TARGET_BASE_LINK_G3,
    TARGET_BASE_LINK_G4,
)


N = 14
FULL = (1 << N) - 1
SCHEDULED_FACETS = 106


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value >> vertex & 1)


def independent(value: int) -> bool:
    return not any(value & edge == edge for edge in EDGE_MASKS)


def hits_all(value: int, clutter: Sequence[int]) -> bool:
    return all(value & row for row in clutter)


def binomial(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    answer = 1
    for index in range(1, min(k, n - k) + 1):
        answer = answer * (n - min(k, n - k) + index) // index
    return answer


def macaulay_raise(value: int, degree: int) -> int:
    remaining = value
    next_top: int | None = None
    raised = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (
            (next_top is None or top + 1 < next_top)
            and binomial(top + 1, lower) <= remaining
        ):
            top += 1
        remaining -= binomial(top, lower)
        raised += binomial(top + 1, lower + 1)
        next_top = top
    if remaining:
        raise AssertionError("inexact Macaulay expansion")
    return raised


def minimal_transversals(clutter: Sequence[int]) -> tuple[int, ...]:
    answers = []
    for value in range(1, FULL + 1):
        if not hits_all(value, clutter):
            continue
        if any(
            hits_all(value ^ (1 << vertex), clutter)
            for vertex in vertices(value)
        ):
            continue
        answers.append(value)
    return tuple(answers)


def derive(a4_masks: Sequence[int]) -> dict[str, object]:
    if len(a4_masks) != 2 or len(set(a4_masks)) != 2:
        raise ValueError("exactly two distinct A4 masks are required")
    if any(value.bit_count() != 4 or not independent(value) for value in a4_masks):
        raise ValueError("invalid independent A4 mask")
    fives = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and independent(value)
    )
    shadow = tuple(
        value
        for value in fives
        if any(value & four == four for four in a4_masks)
    )
    if len(shadow) != 15:
        raise AssertionError("this backend is only for the shadow-15/add5 cases")
    lower = []
    upper = []
    for vertex in range(N):
        a4_incidence = sum(value >> vertex & 1 for value in a4_masks)
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
    base = (*EDGE_MASKS, *a4_masks, *shadow)
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


def solve(a4_masks: Sequence[int]) -> dict[str, object]:
    derived = derive(a4_masks)
    candidates: tuple[int, ...] = derived["candidates"]  # type: ignore[assignment]
    small: tuple[int, ...] = derived["small"]  # type: ignore[assignment]
    order_five: tuple[int, ...] = derived["order_five"]  # type: ignore[assignment]
    lower: tuple[int, ...] = derived["lower"]  # type: ignore[assignment]
    upper: tuple[int, ...] = derived["upper"]  # type: ignore[assignment]
    base: tuple[int, ...] = derived["base"]  # type: ignore[assignment]

    solver = z3.Solver()
    chosen = tuple(z3.Bool(f"x_{index}") for index in range(len(candidates)))
    solver.add(z3.PbEq([(variable, 1) for variable in chosen], 5))
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
        if len(indices) != 5:
            raise AssertionError("model lost the five-row cardinality")
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
        "kind": "h9_14_record18_a2_add5_ridge_aware_z3_scout",
        "theorem_certificate": False,
        "independent_replay_complete": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "compact_only": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "z3_version": z3.get_version_string(),
        "representative_A4_masks": list(a4_masks),
        "representative_A4s": [list(vertices(value)) for value in a4_masks],
        "forced_A5_shadow_count": len(derived["shadow"]),
        "additional_A5_count": 5,
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
            "Exact ridge-aware compact incidence scout; direct blocker audit "
            "on every Boolean model. No topology, Coxeter-label, Gram, or "
            "finite-volume noncompact conclusion."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("a4_mask", nargs=2, type=lambda value: int(value, 0))
    args = parser.parse_args()
    print(json.dumps(solve(tuple(args.a4_mask)), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
