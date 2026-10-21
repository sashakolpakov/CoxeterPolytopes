#!/usr/bin/env python3
"""Independent Z3 replay of the Record-0, ``a=1`` dual nerve model.

This file deliberately does not import the production CP-SAT builder.  A
Boolean variable selects each possible complement of a rank-nine facet of
the nerve.  The model imposes the scheduled 178 facets, every required face
through rank four, exactly 32 independent missing five-sets, the simplicial
ridge degree law, and closure of every nonface above order five.  Optional
generator lex constraints are lossless because the lexicographically least
member of every fixed-source orbit satisfies all of them.

An ``unsat`` result is an incidence obstruction for the fixed source in the
compact Coxeter problem.  It does not address finite-volume noncompact
polytopes.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from time import perf_counter
from typing import Iterable, Sequence

import z3


N = 14
FULL = (1 << N) - 1
MISSING_EDGES = (
    (1 << 0) | (1 << 13),
    (1 << 1) | (1 << 13),
)
SOURCE_FOURS = (
    (2, 3, 4, 5),
    (0, 2, 3, 4),
    (0, 1, 2, 3),
    (2, 3, 4, 13),
)
SCHEDULED_FACETS = 178
SCHEDULED_MISSING_FIVES = 32


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(index for index in range(N) if value & (1 << index))


def contains_base(value: int, four: int) -> bool:
    return any(value & row == row for row in (*MISSING_EDGES, four))


def hits_base(value: int, four: int) -> bool:
    return all(value & row for row in (*MISSING_EDGES, four))


def independent(value: int) -> bool:
    return all(value & edge != edge for edge in MISSING_EDGES)


def move(value: int, action: Sequence[int]) -> int:
    return mask(action[vertex] for vertex in vertices(value))


def swap(left: int, right: int) -> tuple[int, ...]:
    action = list(range(N))
    action[left], action[right] = action[right], action[left]
    return tuple(action)


def source_generators(four: int) -> tuple[tuple[int, ...], ...]:
    isolates = tuple(range(2, 13))
    inside = tuple(vertex for vertex in isolates if four & (1 << vertex))
    outside = tuple(vertex for vertex in isolates if not four & (1 << vertex))
    answer = [
        swap(block[index], block[index + 1])
        for block in (inside, outside)
        for index in range(len(block) - 1)
    ]
    active = swap(0, 1)
    if move(four, active) == four:
        answer.append(active)
    return tuple(answer)


def add_involution_lex_leader(
    solver: z3.Solver,
    variables: Sequence[z3.BoolRef],
    permutation: Sequence[int],
) -> int:
    """Impose ``variables <=lex permuted variables`` independently."""

    if sorted(permutation) != list(range(len(variables))):
        raise AssertionError("candidate action is not a permutation")
    if any(permutation[permutation[i]] != i for i in range(len(variables))):
        raise AssertionError("generator is not an involution")
    representatives = tuple(
        index for index, partner in enumerate(permutation) if index < partner
    )
    prefix: z3.BoolRef = z3.BoolVal(True)
    for index in representatives:
        partner = permutation[index]
        # False < True in the characteristic-vector order.
        solver.add(
            z3.Implies(
                prefix,
                z3.Or(z3.Not(variables[index]), variables[partner]),
            )
        )
        prefix = z3.And(prefix, variables[index] == variables[partner])
    return len(representatives)


def pb_eq(variables: Sequence[z3.BoolRef], value: int) -> z3.BoolRef:
    return z3.PbEq([(variable, 1) for variable in variables], value)


def build_model(
    source_orbit: int, *, symmetry_breaking: bool
) -> tuple[z3.Solver, dict[str, object]]:
    four = mask(SOURCE_FOURS[source_orbit])
    facets = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and hits_base(value, four)
    )
    facet_index = {facet: index for index, facet in enumerate(facets)}
    selected = tuple(z3.Bool(f"facet_{index}") for index in range(len(facets)))
    solver = z3.Solver()
    solver.add(pb_eq(selected, SCHEDULED_FACETS))

    face_witness_count = 0
    for order in range(5):
        for chosen in combinations(range(N), order):
            face = mask(chosen)
            if contains_base(face, four):
                continue
            witnesses = tuple(
                selected[index]
                for index, facet in enumerate(facets)
                if not facet & face
            )
            solver.add(z3.Or(*witnesses))
            face_witness_count += 1

    five_ground = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and independent(value)
    )
    missing = tuple(
        z3.Bool(f"missing5_{index}") for index in range(len(five_ground))
    )
    missing_by_mask = dict(zip(five_ground, missing, strict=True))
    for five, variable in zip(five_ground, missing, strict=True):
        witnesses = tuple(
            selected[index]
            for index, facet in enumerate(facets)
            if not facet & five
        )
        solver.add(variable == z3.Not(z3.Or(*witnesses)))
    solver.add(pb_eq(missing, SCHEDULED_MISSING_FIVES))

    lanner_order_count = 0
    for value in range(FULL + 1):
        if value.bit_count() < 6 or contains_base(value, four):
            continue
        contained_missing = tuple(
            missing_by_mask[mask(chosen)]
            for chosen in combinations(vertices(value), 5)
        )
        witnesses = tuple(
            selected[index]
            for index, facet in enumerate(facets)
            if not facet & value
        )
        solver.add(z3.Or(*contained_missing, *witnesses))
        lanner_order_count += 1

    ridge_count = 0
    for chosen in combinations(range(N), 6):
        six = mask(chosen)
        members = tuple(
            selected[index]
            for index, facet in enumerate(facets)
            if facet & six == facet
        )
        if not members:
            continue
        solver.add(z3.Or(pb_eq(members, 0), pb_eq(members, 2)))
        ridge_count += 1

    moved_pair_counts = []
    generators = source_generators(four) if symmetry_breaking else ()
    for action in generators:
        permutation = tuple(
            facet_index[move(facet, action)] for facet in facets
        )
        moved_pair_counts.append(
            add_involution_lex_leader(solver, selected, permutation)
        )

    return solver, {
        "source_orbit_index": source_orbit,
        "representative_A4": list(SOURCE_FOURS[source_orbit]),
        "candidate_facet_complement_count": len(facets),
        "independent_five_count": len(five_ground),
        "face_witness_constraint_count": face_witness_count,
        "lanner_order_constraint_count": lanner_order_count,
        "ridge_constraint_count": ridge_count,
        "symmetry_breaking": symmetry_breaking,
        "stabilizer_generator_count": len(generators),
        "generator_moved_pair_counts": moved_pair_counts,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_orbit", type=int, choices=range(4))
    parser.add_argument("--timeout-seconds", type=float)
    parser.add_argument("--no-symmetry", action="store_true")
    args = parser.parse_args()
    solver, result = build_model(
        args.source_orbit, symmetry_breaking=not args.no_symmetry
    )
    if args.timeout_seconds is not None:
        solver.set(timeout=round(1000 * args.timeout_seconds))
    before = perf_counter()
    status = solver.check()
    result.update(
        {
            "solver": "Z3",
            "solver_version": z3.get_version_string(),
            "solver_status": str(status),
            "wall_time_seconds": perf_counter() - before,
            "reason_unknown": solver.reason_unknown()
            if status == z3.unknown
            else None,
            "compact_incidence_no_go": status == z3.unsat,
            "finite_volume_noncompact_checked": False,
        }
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
