#!/usr/bin/env python3
"""Exact compact-incidence scout for Record 0, schedule ``a=2``.

The Boolean variables are the possible complements of rank-nine facets.
For one of the 36 exact source orbits of two missing four-sets, the model
imposes the scheduled 176 facets, every required face through rank four,
exactly 38 independent missing five-sets, closure of every larger nonface by
a missing set of order at most five, and the simplicial ridge law.  This is
an exploratory scout: a terminal infeasible solve is exact for the fixed
source, but is not promoted here to a frozen theorem certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from math import factorial
from typing import Iterable, Sequence

from ortools.sat.python import cp_model

from scripts.explore_h9_14_open_positive_completion import audit_compact_labels
from scripts.explore_h9_14_record0_a1_dual_facets import _move
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    minimal_transversals,
    pair,
    placement_options,
    restriction_finite,
    vertices,
)
from scripts.independent_h9_14_record14_a1_audit import (
    _add_involution_lex_leader,
)


N = 14
FULL = (1 << N) - 1
MISSING_EDGES = ((1 << 0) | (1 << 13), (1 << 1) | (1 << 13))
SCHEDULED_FACETS = 176
SCHEDULED_MISSING_FIVES = 38
SOURCE_A4_FAMILIES = (
    ((2, 3, 4, 5), (2, 3, 4, 6)),
    ((2, 3, 4, 5), (2, 3, 6, 7)),
    ((2, 3, 4, 5), (2, 6, 7, 8)),
    ((2, 3, 4, 5), (6, 7, 8, 9)),
    ((0, 2, 3, 4), (2, 3, 4, 5)),
    ((0, 2, 3, 4), (2, 3, 5, 6)),
    ((0, 2, 3, 4), (2, 5, 6, 7)),
    ((0, 2, 3, 4), (5, 6, 7, 8)),
    ((0, 1, 2, 3), (2, 3, 4, 5)),
    ((0, 1, 2, 3), (2, 4, 5, 6)),
    ((0, 1, 2, 3), (4, 5, 6, 7)),
    ((2, 3, 4, 5), (2, 3, 4, 13)),
    ((2, 3, 4, 5), (2, 3, 6, 13)),
    ((2, 3, 4, 5), (2, 6, 7, 13)),
    ((2, 3, 4, 5), (6, 7, 8, 13)),
    ((0, 2, 3, 4), (0, 2, 3, 5)),
    ((0, 2, 3, 4), (0, 2, 5, 6)),
    ((0, 2, 3, 4), (0, 5, 6, 7)),
    ((0, 2, 3, 4), (1, 2, 3, 4)),
    ((0, 2, 3, 4), (1, 2, 3, 5)),
    ((0, 2, 3, 4), (1, 2, 5, 6)),
    ((0, 2, 3, 4), (1, 5, 6, 7)),
    ((0, 1, 2, 3), (0, 2, 3, 4)),
    ((0, 1, 2, 3), (0, 2, 4, 5)),
    ((0, 1, 2, 3), (0, 4, 5, 6)),
    ((0, 2, 3, 4), (2, 3, 4, 13)),
    ((0, 2, 3, 4), (2, 3, 5, 13)),
    ((0, 2, 3, 4), (2, 5, 6, 13)),
    ((0, 2, 3, 4), (5, 6, 7, 13)),
    ((0, 1, 2, 3), (0, 1, 2, 4)),
    ((0, 1, 2, 3), (0, 1, 4, 5)),
    ((0, 1, 2, 3), (2, 3, 4, 13)),
    ((0, 1, 2, 3), (2, 4, 5, 13)),
    ((0, 1, 2, 3), (4, 5, 6, 13)),
    ((2, 3, 4, 13), (2, 5, 6, 13)),
    ((2, 3, 4, 13), (5, 6, 7, 13)),
)


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def independent(value: int) -> bool:
    return all(value & edge != edge for edge in MISSING_EDGES)


def contains_base(value: int, fours: Sequence[int]) -> bool:
    return any(value & row == row for row in (*MISSING_EDGES, *fours))


def _column(fours: Sequence[int], vertex: int) -> tuple[bool, bool]:
    if len(fours) != 2:
        raise AssertionError("a=2 needs exactly two A4 rows")
    return (
        bool(fours[0] & (1 << vertex)),
        bool(fours[1] & (1 << vertex)),
    )


def source_stabilizer_generators(
    fours: Sequence[int],
) -> tuple[tuple[tuple[int, ...], ...], int]:
    """Return involutory generators and exact setwise stabilizer order."""

    isolates = tuple(range(2, 13))
    columns = ((False, False), (False, True), (True, False), (True, True))
    classes = {
        column: tuple(v for v in isolates if _column(fours, v) == column)
        for column in columns
    }
    actions: list[tuple[int, ...]] = []
    for block in classes.values():
        for index in range(len(block) - 1):
            action = list(range(N))
            action[block[index]], action[block[index + 1]] = (
                action[block[index + 1]],
                action[block[index]],
            )
            actions.append(tuple(action))

    quotient_order = 0
    for active_swap in (False, True):
        for row_swap in (False, True):
            def transformed(column: tuple[bool, bool]) -> tuple[bool, bool]:
                return (column[1], column[0]) if row_swap else column

            if any(
                _column(fours, 1 - vertex if active_swap else vertex)
                != transformed(_column(fours, vertex))
                for vertex in (0, 1)
            ):
                continue
            if _column(fours, 13) != transformed(_column(fours, 13)):
                continue
            if any(
                len(classes[column]) != len(classes[transformed(column)])
                for column in columns
            ):
                continue
            quotient_order += 1
            if not active_swap and not row_swap:
                continue
            action = list(range(N))
            if active_swap:
                action[0], action[1] = 1, 0
            for column, block in classes.items():
                target = classes[transformed(column)]
                for source, image in zip(block, target, strict=True):
                    action[source] = image
            action_tuple = tuple(action)
            if any(action_tuple[action_tuple[v]] != v for v in range(N)):
                raise AssertionError("canonical source-stabilizer lift is not involutory")
            actions.append(action_tuple)

    actions = list(dict.fromkeys(actions))
    stabilizer_order = quotient_order
    for block in classes.values():
        stabilizer_order *= factorial(len(block))
    source = set(fours)
    for action in actions:
        if {_move(four, action) for four in fours} != source:
            raise AssertionError("source-stabilizer generator moved the A4 family")
    return tuple(actions), stabilizer_order


def source_representatives() -> tuple[tuple[int, int], ...]:
    representatives = tuple(
        tuple(mask(face) for face in family)
        for family in SOURCE_A4_FAMILIES
    )
    if len(representatives) != 36:
        raise AssertionError("Record-0/a=2 source-orbit count changed")
    return representatives


def audit_source_local_labels(fours: Sequence[int]) -> dict[str, object]:
    """Test the complete two-Lanner overlap against every forced small face."""

    supports = tuple(vertices(four) for four in fours)
    cluster = tuple(sorted(set().union(*map(set, supports))))
    restrictions = tuple(
        chosen
        for order in range(1, min(4, len(cluster)) + 1)
        for chosen in combinations(cluster, order)
        if not contains_base(mask(chosen), fours)
    )
    states: set[tuple[tuple[tuple[int, int], int], ...]] = {()}
    stage_counts = []
    for support in supports:
        passing: set[tuple[tuple[tuple[int, int], int], ...]] = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in placement_options(support):
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(
                    restriction_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = passing
        stage_counts.append(
            {
                "support": list(support),
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
            }
        )
    return {
        "union": list(cluster),
        "union_order": len(cluster),
        "forced_small_face_restriction_count": len(restrictions),
        "stage_counts": stage_counts,
        "terminal_sparse_label_state_count": len(states),
        "local_compact_label_no_go": not states,
        "scope": (
            "Complete order-four Lanner placement tables and every forced "
            "face of order at most four. Unassigned meeting labels are set "
            "to 2; entrywise Gram monotonicity makes a nonfinite restriction "
            "at 2 terminal. Surviving states are not Gram realizations."
        ),
    }


def build_model(
    fours: Sequence[int], *, include_lanner_closure: bool = True,
    symmetry_breaking: bool = True,
) -> dict[str, object]:
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
    model = cp_model.CpModel()
    selected = tuple(
        model.new_bool_var(f"facet_{index}")
        for index in range(len(facets))
    )
    model.add(sum(selected) == SCHEDULED_FACETS)

    face_witness_count = 0
    for order in range(5):
        for chosen in combinations(range(N), order):
            face = mask(chosen)
            if contains_base(face, fours):
                continue
            witnesses = tuple(
                selected[index]
                for index, facet in enumerate(facets)
                if not facet & face
            )
            if witnesses:
                model.add_bool_or(witnesses)
            else:
                model.add(0 == 1)
            face_witness_count += 1

    five_ground = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and independent(value)
    )
    missing = []
    missing_by_mask = {}
    for index, five in enumerate(five_ground):
        variable = model.new_bool_var(f"missing5_{index}")
        witnesses = tuple(
            selected[position]
            for position, facet in enumerate(facets)
            if not facet & five
        )
        if witnesses:
            model.add(sum(witnesses) == 0).only_enforce_if(variable)
            model.add_bool_or(witnesses).only_enforce_if(variable.Not())
        else:
            model.add(variable == 1)
        missing.append(variable)
        missing_by_mask[five] = variable
    model.add(sum(missing) == SCHEDULED_MISSING_FIVES)

    lanner_order_count = 0
    if include_lanner_closure:
        for value in range(FULL + 1):
            if value.bit_count() < 6 or contains_base(value, fours):
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
            model.add_bool_or((*contained_missing, *witnesses))
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
        degree = model.new_int_var(0, len(members), f"ridge_{six}")
        model.add(degree == sum(members))
        model.add_allowed_assignments([degree], ((0,), (2,)))
        ridge_count += 1

    generators, stabilizer_order = source_stabilizer_generators(fours)
    if not symmetry_breaking:
        generators = ()
    facet_index = {facet: index for index, facet in enumerate(facets)}
    moved_pair_counts = []
    for index, action in enumerate(generators):
        permutation = tuple(
            facet_index[_move(facet, action)] for facet in facets
        )
        moved_pair_counts.append(
            _add_involution_lex_leader(
                model, selected, permutation, f"lex_{index}"
            )
        )

    return {
        "model": model,
        "selected": selected,
        "facets": facets,
        "five_ground": five_ground,
        "missing": tuple(missing),
        "shadow": shadow,
        "face_witness_constraint_count": face_witness_count,
        "lanner_order_constraint_count": lanner_order_count,
        "ridge_constraint_count": ridge_count,
        "stabilizer_order": stabilizer_order,
        "stabilizer_generator_count": len(generators),
        "generator_moved_pair_counts": moved_pair_counts,
    }


def solve_source(
    orbit_index: int,
    fours: Sequence[int],
    seconds: float | None,
    workers: int,
    include_lanner_closure: bool,
    symmetry_breaking: bool,
    audit_labels: bool,
) -> dict[str, object]:
    data = build_model(
        fours,
        include_lanner_closure=include_lanner_closure,
        symmetry_breaking=symmetry_breaking,
    )
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    if seconds is not None:
        solver.parameters.max_time_in_seconds = seconds
    status = solver.solve(data["model"])
    target = None
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        facets = tuple(
            facet
            for facet, variable in zip(
                data["facets"], data["selected"], strict=True
            )
            if solver.value(variable)
        )
        blocker = minimal_transversals(facets)
        additions = tuple(face for face in blocker if face.bit_count() == 5)
        target = {
            "selected_facet_complements": [list(vertices(x)) for x in facets],
            "minimal_nonfaces": [list(vertices(x)) for x in blocker],
            "blocker_size_histogram": {
                str(k): v
                for k, v in sorted(Counter(x.bit_count() for x in blocker).items())
            },
            "compact_label_audit": (
                audit_compact_labels(MISSING_EDGES, fours, additions)
                if audit_labels
                else None
            ),
        }
    return {
        "source_orbit_index": orbit_index,
        "representative_A4s": [list(vertices(x)) for x in fours],
        "forced_A5_shadow_count": len(data["shadow"]),
        "additional_A5_count": SCHEDULED_MISSING_FIVES - len(data["shadow"]),
        "candidate_facet_complement_count": len(data["facets"]),
        "face_witness_constraint_count": data["face_witness_constraint_count"],
        "lanner_order_constraint_count": data["lanner_order_constraint_count"],
        "ridge_constraint_count": data["ridge_constraint_count"],
        "lanner_order_closure_installed": include_lanner_closure,
        "source_stabilizer_order": data["stabilizer_order"],
        "stabilizer_generator_count": data["stabilizer_generator_count"],
        "generator_moved_facet_pair_counts": data["generator_moved_pair_counts"],
        "solver_status": solver.status_name(status),
        "solver_wall_time_seconds": solver.wall_time,
        "solver_branch_count": solver.num_branches,
        "solver_conflict_count": solver.num_conflicts,
        "target": target,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbits", nargs="*", type=int)
    parser.add_argument("--seconds", type=float)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--core-only", action="store_true")
    parser.add_argument("--no-symmetry", action="store_true")
    parser.add_argument("--local-labels-only", action="store_true")
    parser.add_argument("--labels", action="store_true")
    args = parser.parse_args()
    representatives = source_representatives()
    indices = args.orbits or list(range(len(representatives)))
    if any(index not in range(len(representatives)) for index in indices):
        parser.error("source orbit must lie in 0,...,35")
    rows = []
    for index in indices:
        if args.local_labels_only:
            audit = audit_source_local_labels(representatives[index])
            row = {
                "source_orbit_index": index,
                "representative_A4s": [
                    list(vertices(x)) for x in representatives[index]
                ],
                "local_label_audit": audit,
            }
            rows.append(row)
            print(
                f"source {index}: local-label states="
                f"{audit['terminal_sparse_label_state_count']}",
                flush=True,
            )
            continue
        row = solve_source(
            index,
            representatives[index],
            args.seconds,
            args.workers,
            not args.core_only,
            not args.no_symmetry,
            args.labels,
        )
        rows.append(row)
        print(
            f"source {index}: {row['solver_status']} "
            f"({row['solver_wall_time_seconds']:.3f}s)",
            flush=True,
        )
    print(
        json.dumps(
            {
                "kind": "h9_14_record0_a2_dual_facet_exploratory_scout",
                "theorem_certificate": False,
                "record_index": 0,
                "graph6": "M?????????????o??",
                "schedule_a": 2,
                "source_orbit_count": 36,
                "rows": rows,
                "scope": (
                    "Exact necessary compact-incidence solve for each listed "
                    "source. UNKNOWN is not an exclusion. No finite-volume "
                    "noncompact inference is made."
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
