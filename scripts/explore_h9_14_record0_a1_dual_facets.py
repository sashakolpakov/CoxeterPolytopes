#!/usr/bin/env python3
"""Exact dual-facet scout for Record 0, schedule ``a=1``.

The older primal model chooses 22--24 additional five-nonfaces among roughly
1,600 candidates.  This equivalent dual model instead chooses the scheduled
178 facet complements among roughly 750 base transversals.  It enforces:

* exactly the fixed two missing edges and one fixed four-nonface below rank 5;
* the exact Dehn--Sommerville number of missing independent five-sets;
* the exact facet count; and
* degree zero or two over every six-set (the ridge law).

A stabilizer-generator lex slice is lossless: every orbit of the finite fixed
source action contains its global lexicographic minimum, which satisfies every
generator comparison.  This remains an exploratory compact-incidence scout
until an independent model and frozen certificate are supplied.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from math import factorial
from typing import Iterable, Sequence

from ortools.sat.python import cp_model

from scripts.explore_h9_14_open_positive_completion import (
    audit_compact_labels,
    build_problem,
)
from scripts.explore_h9_14_record14_a1_symmetry_z3 import (
    _move,
    _swap_action,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    minimal_transversals,
    vertices,
)
from scripts.independent_h9_14_record14_a1_audit import (
    _add_involution_lex_leader,
)


N = 14
FULL = (1 << N) - 1
MISSING_EDGES = ((1 << 0) | (1 << 13), (1 << 1) | (1 << 13))


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def contains_base(value: int, fours: Sequence[int]) -> bool:
    return any(value & row == row for row in (*MISSING_EDGES, *fours))


def independent(value: int) -> bool:
    return all(value & edge != edge for edge in MISSING_EDGES)


def source_generators(four: int) -> tuple[tuple[int, ...], ...]:
    isolates = tuple(range(2, 13))
    inside = tuple(vertex for vertex in isolates if four & (1 << vertex))
    outside = tuple(vertex for vertex in isolates if not four & (1 << vertex))
    generators = [
        _swap_action(block[index], block[index + 1])
        for block in (inside, outside)
        for index in range(len(block) - 1)
    ]
    active_swap = _swap_action(0, 1)
    if _move(four, active_swap) == four:
        generators.append(active_swap)
    return tuple(generators)


def build_dual_model(
    orbit_index: int, *, symmetry_breaking: bool = True
) -> dict[str, object]:
    problem = build_problem(0, 1, orbit_index)
    fours = tuple(problem["fours"])
    if len(fours) != 1:
        raise AssertionError("Record-0 a=1 source changed")
    facets = tuple(problem["fives"])
    scheduled_facets = int(problem["scheduled_facets"])
    total_missing_fives = len(problem["shadow"]) + int(
        problem["addition_count"]
    )

    model = cp_model.CpModel()
    live = [
        model.new_bool_var(f"facet_{index}")
        for index in range(len(facets))
    ]
    model.add(sum(live) == scheduled_facets)

    face_witness_constraint_count = 0
    for order in range(5):
        for chosen in combinations(range(N), order):
            value = mask(chosen)
            if contains_base(value, fours):
                continue
            witnesses = [
                live[index]
                for index, facet in enumerate(facets)
                if not facet & value
            ]
            if witnesses:
                model.add_bool_or(witnesses)
            else:
                model.add(0 == 1)
            face_witness_constraint_count += 1

    five_ground = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and independent(value)
    )
    missing_five = []
    missing_five_by_mask = {}
    for position, five in enumerate(five_ground):
        variable = model.new_bool_var(f"missing5_{position}")
        disjoint_facets = [
            live[index]
            for index, facet in enumerate(facets)
            if not facet & five
        ]
        if disjoint_facets:
            model.add(sum(disjoint_facets) == 0).only_enforce_if(variable)
            model.add_bool_or(disjoint_facets).only_enforce_if(variable.Not())
        else:
            model.add(variable == 1)
        missing_five.append(variable)
        missing_five_by_mask[five] = variable
    model.add(sum(missing_five) == total_missing_fives)

    # Exact Lanner-order closure.  A minimal nonface of a compact
    # hyperbolic Coxeter polytope has order at most five.  Equivalently, any
    # set not already containing a fixed edge/four must either lie in a
    # facet or contain one of the missing independent five-sets.  Enforcing
    # this for every order >=6 rules out *all* higher minimal blockers, not
    # merely the first one found after a solve.
    lanner_order_constraint_count = 0
    for value in range(FULL + 1):
        if value.bit_count() < 6 or contains_base(value, fours):
            continue
        contained_missing_fives = [
            missing_five_by_mask[mask(chosen)]
            for chosen in combinations(vertices(value), 5)
        ]
        facet_witnesses = [
            live[index]
            for index, facet in enumerate(facets)
            if not facet & value
        ]
        model.add_bool_or([*contained_missing_fives, *facet_witnesses])
        lanner_order_constraint_count += 1

    ridge_constraint_count = 0
    for six in range(FULL + 1):
        if six.bit_count() != 6:
            continue
        members = [
            live[index]
            for index, facet in enumerate(facets)
            if facet & six == facet
        ]
        if not members:
            continue
        degree = model.new_int_var(0, len(members), f"ridge_{six}")
        model.add(degree == sum(members))
        model.add_allowed_assignments([degree], [(0,), (2,)])
        ridge_constraint_count += 1

    generators = source_generators(fours[0]) if symmetry_breaking else ()
    facet_index = {facet: index for index, facet in enumerate(facets)}
    moved_pair_counts = []
    for generator_index, action in enumerate(generators):
        permutation = tuple(
            facet_index[_move(facet, action)] for facet in facets
        )
        moved_pair_counts.append(
            _add_involution_lex_leader(
                model, live, permutation, f"lex_{generator_index}"
            )
        )

    isolates_inside = sum(
        bool(fours[0] & (1 << vertex)) for vertex in range(2, 13)
    )
    active_swap = _swap_action(0, 1)
    stabilizer_order = (
        factorial(isolates_inside)
        * factorial(11 - isolates_inside)
        * (2 if _move(fours[0], active_swap) == fours[0] else 1)
    )
    return {
        "problem": problem,
        "model": model,
        "live": live,
        "facets": facets,
        "five_ground": five_ground,
        "missing_five": missing_five,
        "face_witness_constraint_count": face_witness_constraint_count,
        "ridge_constraint_count": ridge_constraint_count,
        "lanner_order_constraint_count": lanner_order_constraint_count,
        "stabilizer_generator_count": len(generators),
        "stabilizer_order": stabilizer_order,
        "generator_moved_pair_counts": moved_pair_counts,
    }


def solve_one(
    orbit_index: int,
    time_limit_seconds: float | None,
    workers: int,
) -> dict[str, object]:
    data = build_dual_model(orbit_index)
    model = data["model"]
    live = tuple(data["live"])
    facets = tuple(data["facets"])
    problem = data["problem"]
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    if time_limit_seconds is not None:
        solver.parameters.max_time_in_seconds = time_limit_seconds
    status = solver.solve(model)
    status_name = solver.status_name(status)
    target: dict[str, object] | None = None
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        selected_facets = tuple(
            facet
            for index, facet in enumerate(facets)
            if solver.value(live[index])
        )
        blocker = minimal_transversals(selected_facets)
        sizes = Counter(value.bit_count() for value in blocker)
        fours = tuple(problem["fours"])
        additions = tuple(
            value
            for value in blocker
            if value.bit_count() == 5
        )
        label_audit = (
            audit_compact_labels(
                problem["edge_rows"], fours, additions
            )
            if sizes == {2: 2, 4: 1, 5: int(problem["addition_count"])}
            else None
        )
        target = {
            "selected_facet_complements": [
                list(vertices(value)) for value in selected_facets
            ],
            "blocker_size_histogram": {
                str(size): count for size, count in sorted(sizes.items())
            },
            "minimal_nonfaces": [list(vertices(value)) for value in blocker],
            "pure_expected_blocker": sizes
            == {2: 2, 4: 1, 5: int(problem["addition_count"])},
            "compact_label_audit": label_audit,
        }
    return {
        "kind": "h9_14_record0_a1_dual_facet_exploratory_scout",
        "theorem_certificate": False,
        "record_index": 0,
        "graph6": problem["graph6"],
        "schedule_a": 1,
        "source_orbit_index": orbit_index,
        "representative_A4s": [
            list(vertices(value)) for value in problem["fours"]
        ],
        "scheduled_facet_count": problem["scheduled_facets"],
        "total_missing_independent_five_count": len(problem["shadow"])
        + int(problem["addition_count"]),
        "candidate_facet_complement_count": len(facets),
        "face_witness_constraint_count": data[
            "face_witness_constraint_count"
        ],
        "ridge_constraint_count": data["ridge_constraint_count"],
        "lanner_order_constraint_count": data[
            "lanner_order_constraint_count"
        ],
        "stabilizer_generator_count": data["stabilizer_generator_count"],
        "source_stabilizer_order": data["stabilizer_order"],
        "generator_moved_facet_pair_counts": data[
            "generator_moved_pair_counts"
        ],
        "solver_status": status_name,
        "solver_wall_time_seconds": solver.wall_time,
        "solver_branch_count": solver.num_branches,
        "solver_conflict_count": solver.num_conflicts,
        "target": target,
        "scope": (
            "Exact compact necessary-incidence solve for one fixed source. "
            "UNKNOWN under a time limit is not an exclusion; no finite-volume "
            "noncompact inference is made."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbit_index", type=int, choices=range(4))
    parser.add_argument("--seconds", type=float)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    print(
        json.dumps(
            solve_one(args.orbit_index, args.seconds, args.workers),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
