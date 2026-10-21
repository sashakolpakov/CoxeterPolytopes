#!/usr/bin/env python3
"""Independent Z3 replay for selected Record-0/a=2 source orbits.

This intentionally does not import the CP-SAT production builder.  It
reconstructs the complete dual-facet model and the exact setwise source
stabilizer from three literal pairs of missing four-sets.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from math import factorial
from time import perf_counter
from typing import Iterable, Sequence

import z3


N = 14
FULL = (1 << N) - 1
MISSING_EDGES = ((1 << 0) | (1 << 13), (1 << 1) | (1 << 13))
SCHEDULED_FACETS = 176
SCHEDULED_MISSING_FIVES = 38
SOURCES = {
    0: ((2, 3, 4, 5), (2, 3, 4, 6)),
    1: ((2, 3, 4, 5), (2, 3, 6, 7)),
    8: ((0, 1, 2, 3), (2, 3, 4, 5)),
    9: ((0, 1, 2, 3), (2, 4, 5, 6)),
    10: ((0, 1, 2, 3), (4, 5, 6, 7)),
    12: ((2, 3, 4, 5), (2, 3, 6, 13)),
    23: ((0, 1, 2, 3), (0, 2, 4, 5)),
    24: ((0, 1, 2, 3), (0, 4, 5, 6)),
    29: ((0, 1, 2, 3), (0, 1, 2, 4)),
    30: ((0, 1, 2, 3), (0, 1, 4, 5)),
    31: ((0, 1, 2, 3), (2, 3, 4, 13)),
    32: ((0, 1, 2, 3), (2, 4, 5, 13)),
    34: ((2, 3, 4, 13), (2, 5, 6, 13)),
}


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(v for v in range(N) if value & (1 << v))


def move(value: int, action: Sequence[int]) -> int:
    return mask(action[v] for v in vertices(value))


def contains_base(value: int, fours: Sequence[int]) -> bool:
    return any(value & row == row for row in (*MISSING_EDGES, *fours))


def independent(value: int) -> bool:
    return all(value & edge != edge for edge in MISSING_EDGES)


def column(fours: Sequence[int], vertex: int) -> tuple[bool, bool]:
    return tuple(bool(four & (1 << vertex)) for four in fours)  # type: ignore[return-value]


def source_generators(
    fours: Sequence[int],
) -> tuple[tuple[tuple[int, ...], ...], int]:
    isolates = tuple(range(2, 13))
    columns = ((False, False), (False, True), (True, False), (True, True))
    classes = {
        code: tuple(v for v in isolates if column(fours, v) == code)
        for code in columns
    }
    actions = []
    for block in classes.values():
        for index in range(len(block) - 1):
            action = list(range(N))
            action[block[index]], action[block[index + 1]] = (
                action[block[index + 1]], action[block[index]]
            )
            actions.append(tuple(action))
    quotient_order = 0
    for active_swap in (False, True):
        for row_swap in (False, True):
            def transform(code: tuple[bool, bool]) -> tuple[bool, bool]:
                return (code[1], code[0]) if row_swap else code

            if any(
                column(fours, 1 - v if active_swap else v)
                != transform(column(fours, v))
                for v in (0, 1)
            ):
                continue
            if column(fours, 13) != transform(column(fours, 13)):
                continue
            if any(
                len(classes[code]) != len(classes[transform(code)])
                for code in columns
            ):
                continue
            quotient_order += 1
            if not active_swap and not row_swap:
                continue
            action = list(range(N))
            if active_swap:
                action[0], action[1] = 1, 0
            for code, block in classes.items():
                for source, target in zip(
                    block, classes[transform(code)], strict=True
                ):
                    action[source] = target
            actions.append(tuple(action))
    actions = list(dict.fromkeys(actions))
    order = quotient_order
    for block in classes.values():
        order *= factorial(len(block))
    for action in actions:
        if any(action[action[v]] != v for v in range(N)):
            raise AssertionError("stabilizer generator is not involutory")
        if {move(four, action) for four in fours} != set(fours):
            raise AssertionError("stabilizer generator moved source")
    return tuple(actions), order


def pb_eq(variables: Sequence[z3.BoolRef], value: int) -> z3.BoolRef:
    return z3.PbEq([(variable, 1) for variable in variables], value)


def add_lex(
    solver: z3.Solver,
    variables: Sequence[z3.BoolRef],
    permutation: Sequence[int],
) -> int:
    representatives = tuple(
        i for i, partner in enumerate(permutation) if i < partner
    )
    prefix: z3.BoolRef = z3.BoolVal(True)
    for i in representatives:
        partner = permutation[i]
        solver.add(
            z3.Implies(prefix, z3.Or(z3.Not(variables[i]), variables[partner]))
        )
        prefix = z3.And(prefix, variables[i] == variables[partner])
    return len(representatives)


def build_model(source_orbit: int) -> tuple[z3.Solver, dict[str, object]]:
    fours = tuple(mask(face) for face in SOURCES[source_orbit])
    shadow = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5
        and independent(value)
        and any(value & four == four for four in fours)
    )
    base = (*MISSING_EDGES, *fours, *shadow)
    facets = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and all(value & row for row in base)
    )
    selected = tuple(z3.Bool(f"facet_{i}") for i in range(len(facets)))
    solver = z3.Solver()
    solver.add(pb_eq(selected, SCHEDULED_FACETS))

    face_count = 0
    for order in range(5):
        for chosen in combinations(range(N), order):
            face = mask(chosen)
            if contains_base(face, fours):
                continue
            witnesses = tuple(
                selected[i]
                for i, facet in enumerate(facets)
                if not facet & face
            )
            solver.add(z3.Or(*witnesses))
            face_count += 1

    five_ground = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and independent(value)
    )
    missing = tuple(z3.Bool(f"missing5_{i}") for i in range(len(five_ground)))
    missing_by_mask = dict(zip(five_ground, missing, strict=True))
    for five, variable in zip(five_ground, missing, strict=True):
        witnesses = tuple(
            selected[i]
            for i, facet in enumerate(facets)
            if not facet & five
        )
        solver.add(variable == z3.Not(z3.Or(*witnesses)))
    solver.add(pb_eq(missing, SCHEDULED_MISSING_FIVES))

    closure_count = 0
    for value in range(FULL + 1):
        if value.bit_count() < 6 or contains_base(value, fours):
            continue
        contained = tuple(
            missing_by_mask[mask(chosen)]
            for chosen in combinations(vertices(value), 5)
        )
        witnesses = tuple(
            selected[i]
            for i, facet in enumerate(facets)
            if not facet & value
        )
        solver.add(z3.Or(*contained, *witnesses))
        closure_count += 1

    ridge_count = 0
    for chosen in combinations(range(N), 6):
        six = mask(chosen)
        members = tuple(
            selected[i]
            for i, facet in enumerate(facets)
            if facet & six == facet
        )
        if members:
            solver.add(z3.Or(pb_eq(members, 0), pb_eq(members, 2)))
            ridge_count += 1

    generators, group_order = source_generators(fours)
    facet_index = {facet: i for i, facet in enumerate(facets)}
    moved_counts = []
    for action in generators:
        permutation = tuple(facet_index[move(facet, action)] for facet in facets)
        moved_counts.append(add_lex(solver, selected, permutation))
    return solver, {
        "source_orbit_index": source_orbit,
        "representative_A4s": [list(face) for face in SOURCES[source_orbit]],
        "forced_A5_shadow_count": len(shadow),
        "additional_A5_count": SCHEDULED_MISSING_FIVES - len(shadow),
        "candidate_facet_complement_count": len(facets),
        "face_witness_constraint_count": face_count,
        "lanner_order_constraint_count": closure_count,
        "ridge_constraint_count": ridge_count,
        "source_stabilizer_order": group_order,
        "stabilizer_generator_count": len(generators),
        "generator_moved_facet_pair_counts": moved_counts,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_orbit", type=int, choices=tuple(SOURCES))
    parser.add_argument("--timeout-seconds", type=float)
    args = parser.parse_args()
    solver, record = build_model(args.source_orbit)
    if args.timeout_seconds is not None:
        solver.set(timeout=round(1000 * args.timeout_seconds))
    before = perf_counter()
    status = solver.check()
    record.update(
        {
            "solver": "Z3",
            "solver_version": z3.get_version_string(),
            "solver_status": str(status),
            "wall_time_seconds": perf_counter() - before,
            "reason_unknown": solver.reason_unknown()
            if status == z3.unknown
            else None,
        }
    )
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
