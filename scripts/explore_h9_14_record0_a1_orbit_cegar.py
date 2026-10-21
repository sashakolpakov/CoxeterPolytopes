#!/usr/bin/env python3
"""Full-source orbit CEGAR for the four Record-0 ``a=1`` sources.

Each iteration asks the exact dual-facet model for one pure incidence target,
runs the complete unrestricted compact label gate, and—when labels die—adds
no-goods for every image under the full fixed-source group.  Thus a later
target is never merely a relabelling of an already excluded one.

Terminal INFEASIBLE after the orbit cuts closes the fixed source.  A timeout
is reported as UNKNOWN, and a surviving label state is handed off to the
exact Gram stage; neither is silently converted into a theorem.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import permutations
import json
from math import factorial
from typing import Iterable, Iterator, Sequence

from ortools.sat.python import cp_model

from scripts.explore_h9_14_open_positive_completion import audit_compact_labels
from scripts.explore_h9_14_record0_a1_dual_facets import (
    build_dual_model,
    contains_base,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    minimal_transversals,
    vertices,
)


N = 14


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _move(value: int, action: Sequence[int]) -> int:
    return sum(1 << action[vertex] for vertex in vertices(value))


def source_actions(four: int) -> Iterator[tuple[int, ...]]:
    isolates = tuple(range(2, 13))
    inside = tuple(vertex for vertex in isolates if four & (1 << vertex))
    outside = tuple(vertex for vertex in isolates if not four & (1 << vertex))
    active_swaps = (
        (False, True)
        if bool(four & 1) == bool(four & 2)
        else (False,)
    )
    for inside_image in permutations(inside):
        for outside_image in permutations(outside):
            base = list(range(N))
            for source, image in zip(inside, inside_image, strict=True):
                base[source] = image
            for source, image in zip(outside, outside_image, strict=True):
                base[source] = image
            for swap_active in active_swaps:
                action = base.copy()
                if swap_active:
                    action[0], action[1] = 1, 0
                yield tuple(action)


def source_group_order(four: int) -> int:
    inside = sum(bool(four & (1 << vertex)) for vertex in range(2, 13))
    active = 2 if bool(four & 1) == bool(four & 2) else 1
    return factorial(inside) * factorial(11 - inside) * active


def orbit_images(
    family: Sequence[int], four: int
) -> set[tuple[int, ...]]:
    return {
        tuple(sorted(_move(face, action) for face in family))
        for action in source_actions(four)
    }


def run_cegar(
    orbit_index: int,
    maximum_orbits: int | None,
    seconds_per_solve: float | None,
) -> dict[str, object]:
    data = build_dual_model(orbit_index)
    problem = data["problem"]
    model = data["model"]
    live = tuple(data["live"])
    facets = tuple(data["facets"])
    missing_five = tuple(data["missing_five"])
    five_ground = tuple(data["five_ground"])
    five_index = {face: index for index, face in enumerate(five_ground)}
    fours = tuple(problem["fours"])
    four = fours[0]
    group_order = source_group_order(four)
    iterations = []
    terminal_status = "not_started"
    stopped_reason = "not_started"

    while maximum_orbits is None or len(iterations) < maximum_orbits:
        solver = cp_model.CpSolver()
        solver.parameters.num_search_workers = 1
        if seconds_per_solve is not None:
            solver.parameters.max_time_in_seconds = seconds_per_solve
        status = solver.solve(model)
        terminal_status = solver.status_name(status)
        if status == cp_model.INFEASIBLE:
            stopped_reason = "terminal_infeasible_after_full_orbit_cuts"
            break
        if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
            stopped_reason = "solver_unknown_or_timeout"
            break

        selected_facets = tuple(
            facet
            for index, facet in enumerate(facets)
            if solver.value(live[index])
        )
        additions = tuple(
            face
            for face, variable in zip(
                five_ground, missing_five, strict=True
            )
            if solver.value(variable) and not contains_base(face, fours)
        )
        blocker = minimal_transversals(selected_facets)
        expected_blocker = {
            *problem["edge_rows"],
            four,
            *additions,
        }
        if set(blocker) != expected_blocker:
            raise AssertionError("dual target failed exact blocker involution")
        label_audit = audit_compact_labels(
            problem["edge_rows"], fours, additions
        )
        images = orbit_images(additions, four)
        if group_order % len(images):
            raise AssertionError("target orbit size does not divide source group")
        iteration = {
            "target_orbit_index": len(iterations),
            "representative_additional_A5s": [
                list(vertices(face)) for face in additions
            ],
            "representative_digest": _digest(additions),
            "selected_facet_complement_count": len(selected_facets),
            "source_group_order": group_order,
            "target_stabilizer_order": group_order // len(images),
            "full_target_orbit_size": len(images),
            "solver_status": terminal_status,
            "solver_wall_time_seconds": solver.wall_time,
            "solver_branch_count": solver.num_branches,
            "solver_conflict_count": solver.num_conflicts,
            "compact_label_audit": label_audit,
        }
        iterations.append(iteration)
        print(
            f"orbit {len(iterations) - 1}: size={len(images)} "
            f"labels={label_audit['terminal_sparse_label_state_count']}",
            flush=True,
        )
        if not label_audit["unrestricted_compact_label_no_go"]:
            stopped_reason = "compact_label_survivor_requires_exact_Gram"
            break
        # At an explicit exploratory cutoff, retain the exact orbit size and
        # label audit above but do not spend time/memory installing clauses
        # that cannot be used by a subsequent solve in this run.
        if maximum_orbits is not None and len(iterations) >= maximum_orbits:
            stopped_reason = "orbit_cutoff"
            break
        for image in images:
            model.add_bool_or(
                [missing_five[five_index[face]].Not() for face in image]
            )
        print(
            f"blocked {len(images)} full-source images",
            flush=True,
        )
    all_killed = all(
        row["compact_label_audit"]["unrestricted_compact_label_no_go"]
        for row in iterations
    )
    source_closed = (
        stopped_reason == "terminal_infeasible_after_full_orbit_cuts"
        and all_killed
    )
    return {
        "kind": "h9_14_record0_a1_full_source_orbit_cegar_exploratory",
        "theorem_certificate": False,
        "record_index": 0,
        "graph6": problem["graph6"],
        "schedule_a": 1,
        "source_orbit_index": orbit_index,
        "representative_A4s": [list(vertices(four))],
        "source_group_order": group_order,
        "dual_model_counts": {
            "candidate_facet_complements": len(facets),
            "face_witness_constraints": data[
                "face_witness_constraint_count"
            ],
            "lanner_order_constraints": data[
                "lanner_order_constraint_count"
            ],
            "ridge_constraints": data["ridge_constraint_count"],
        },
        "excluded_full_target_orbit_count": len(iterations),
        "excluded_labelled_target_count": sum(
            row["full_target_orbit_size"] for row in iterations
        ),
        "iterations": iterations,
        "terminal_solver_status": terminal_status,
        "stopped_reason": stopped_reason,
        "fixed_source_compact_coxeter_no_go": source_closed,
        "finite_volume_noncompact_checked": False,
        "scope": (
            "Exact compact incidence plus unrestricted compact labels for one "
            "fixed a=1 source. Only terminal INFEASIBLE after all full-orbit "
            "cuts closes the source; no noncompact inference is made."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbit_index", type=int, choices=range(4))
    parser.add_argument("--max-orbits", type=int)
    parser.add_argument("--seconds-per-solve", type=float)
    args = parser.parse_args()
    print(
        json.dumps(
            run_cegar(
                args.orbit_index,
                args.max_orbits,
                args.seconds_per_solve,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
