#!/usr/bin/env python3
"""Local compact-label-core CEGAR for Record 0, schedule ``a=1``.

For each exact pure dual-nerve target, the unrestricted compact label audit
either finds a survivor or supplies a small contradiction involving only a
few minimal Lanner supports and actual vertices.  Every full-source image of
that local contradiction is installed as a no-good, usually excluding many
complete A5 families at once.

Only terminal ``INFEASIBLE`` after these lossless cuts closes a fixed source.
A label survivor is sent to exact Gram realization.  Timeout/cutoff remains
explicitly inconclusive, and finite-volume noncompact polytopes are outside
the scope of this compact calculation.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Sequence

from ortools.sat.python import cp_model

from scripts.explore_h9_14_open_positive_completion import audit_compact_labels
from scripts.explore_h9_14_record0_a1_dual_facets import (
    build_dual_model,
    contains_base,
)
from scripts.explore_h9_14_record0_a1_label_core import extract_label_core, mask
from scripts.explore_h9_14_record0_a1_orbit_cegar import (
    _move,
    source_actions,
    source_group_order,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    minimal_transversals,
    vertices,
)


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def core_orbit_images(
    additions: Sequence[int], restrictions: Sequence[int], four: int
) -> set[tuple[tuple[int, ...], tuple[int, ...]]]:
    return {
        (
            tuple(sorted(_move(value, action) for value in additions)),
            tuple(sorted(_move(value, action) for value in restrictions)),
        )
        for action in source_actions(four)
    }


def run_core_cegar(
    source_orbit: int,
    maximum_cores: int | None,
    seconds_per_solve: float | None,
    seeded_iterations: Sequence[dict[str, object]] = (),
) -> dict[str, object]:
    data = build_dual_model(source_orbit)
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
    iterations = [dict(row) for row in seeded_iterations]
    terminal_status = "not_started"
    stopped_reason = "not_started"
    face_variables: dict[int, cp_model.IntVar] = {}

    def face_variable(restriction: int) -> cp_model.IntVar:
        """Whether the restriction extends to a selected rank-nine vertex."""

        known = face_variables.get(restriction)
        if known is not None:
            return known
        variable = model.new_bool_var(f"core_face_{restriction}")
        witnesses = tuple(
            live[index]
            for index, facet in enumerate(facets)
            if not facet & restriction
        )
        if witnesses:
            model.add(sum(witnesses) >= 1).only_enforce_if(variable)
            model.add(sum(witnesses) == 0).only_enforce_if(variable.Not())
        else:
            model.add(variable == 0)
        face_variables[restriction] = variable
        return variable

    def restriction_masks(core: dict[str, object]) -> tuple[int, ...]:
        stored = core.get("core_face_restrictions")
        if stored is not None:
            return tuple(sorted(mask(face) for face in stored))
        union = 0
        for support in core["core_supports"]:
            union |= mask(support)
        restrictions = tuple(
            sorted(mask(actual) & union for actual in core["core_actual_vertices"])
        )
        core["core_face_restrictions"] = [
            list(vertices(value)) for value in restrictions
        ]
        return restrictions

    for row in iterations:
        core = row.get("local_core")
        if not isinstance(core, dict):
            raise ValueError("a resumed row lacks a local label core")
        core_additions = tuple(
            mask(face) for face in core["core_additional_A5s"]
        )
        core_restrictions = restriction_masks(core)
        images = core_orbit_images(core_additions, core_restrictions, four)
        legacy_size = row.get("full_core_orbit_size")
        if "full_face_core_orbit_size" not in row:
            row["legacy_full_vertex_core_orbit_size"] = legacy_size
        row["full_core_orbit_size"] = len(images)
        row["full_face_core_orbit_size"] = len(images)
        row["core_stabilizer_order"] = group_order // len(images)
        for image_additions, image_restrictions in images:
            model.add_bool_or(
                [
                    *(missing_five[five_index[face]].Not() for face in image_additions),
                    *(face_variable(face).Not() for face in image_restrictions),
                ]
            )
        print(
            f"resumed core {row['iteration']}: installed {len(images)} images",
            flush=True,
        )

    while maximum_cores is None or len(iterations) < maximum_cores:
        solver = cp_model.CpSolver()
        solver.parameters.num_search_workers = 1
        if seconds_per_solve is not None:
            solver.parameters.max_time_in_seconds = seconds_per_solve
        status = solver.solve(model)
        terminal_status = solver.status_name(status)
        if status == cp_model.INFEASIBLE:
            stopped_reason = "terminal_infeasible_after_full_core_orbit_cuts"
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
            for face, variable in zip(five_ground, missing_five, strict=True)
            if solver.value(variable) and not contains_base(face, fours)
        )
        blocker = minimal_transversals(selected_facets)
        expected_blocker = {*problem["edge_rows"], four, *additions}
        if set(blocker) != expected_blocker:
            raise AssertionError("dual target failed exact blocker involution")
        label_audit = audit_compact_labels(
            problem["edge_rows"], fours, additions
        )
        if not label_audit["unrestricted_compact_label_no_go"]:
            iterations.append(
                {
                    "iteration": len(iterations),
                    "representative_digest": _digest(additions),
                    "compact_label_audit": label_audit,
                    "solver_wall_time_seconds": solver.wall_time,
                    "local_core": None,
                }
            )
            stopped_reason = "compact_label_survivor_requires_exact_Gram"
            break

        core = extract_label_core(four, additions, selected_facets)
        core_additions = tuple(
            mask(row) for row in core["core_additional_A5s"]
        )
        core_restrictions = restriction_masks(core)
        images = core_orbit_images(core_additions, core_restrictions, four)
        if group_order % len(images):
            raise AssertionError("core orbit size does not divide source group")
        row = {
            "iteration": len(iterations),
            "representative_digest": _digest(additions),
            "selected_facet_complement_count": len(selected_facets),
            "rank_nine_actual_vertex_count": label_audit[
                "rank_nine_actual_vertex_count"
            ],
            "label_stage_state_counts": [
                stage["passing_state_count"] for stage in label_audit["stages"]
            ],
            "solver_wall_time_seconds": solver.wall_time,
            "solver_branch_count": solver.num_branches,
            "solver_conflict_count": solver.num_conflicts,
            "local_core": core,
            "full_core_orbit_size": len(images),
            "full_face_core_orbit_size": len(images),
            "core_stabilizer_order": group_order // len(images),
        }
        iterations.append(row)
        print(
            f"core {len(iterations) - 1}: "
            f"{len(core_additions)} A5 + {len(core_restrictions)} faces; "
            f"orbit={len(images)}",
            flush=True,
        )
        if maximum_cores is not None and len(iterations) >= maximum_cores:
            stopped_reason = "core_cutoff"
            break
        for image_additions, image_restrictions in images:
            model.add_bool_or(
                [
                    *(missing_five[five_index[face]].Not() for face in image_additions),
                    *(face_variable(face).Not() for face in image_restrictions),
                ]
            )
        print(f"installed {len(images)} full-source core images", flush=True)

    source_closed = (
        stopped_reason == "terminal_infeasible_after_full_core_orbit_cuts"
        and all(row["local_core"] is not None for row in iterations)
    )
    return {
        "kind": "h9_14_record0_a1_local_label_core_cegar_exploratory",
        "theorem_certificate": False,
        "record_index": 0,
        "graph6": problem["graph6"],
        "schedule_a": 1,
        "source_orbit_index": source_orbit,
        "representative_A4s": [list(vertices(four))],
        "source_group_order": group_order,
        "dual_model_counts": {
            "candidate_facet_complements": len(facets),
            "face_witness_constraints": data["face_witness_constraint_count"],
            "lanner_order_constraints": data["lanner_order_constraint_count"],
            "ridge_constraints": data["ridge_constraint_count"],
        },
        "excluded_local_core_orbit_count": len(
            [row for row in iterations if row["local_core"] is not None]
        ),
        "iterations": iterations,
        "terminal_solver_status": terminal_status,
        "stopped_reason": stopped_reason,
        "fixed_source_compact_coxeter_no_go": source_closed,
        "finite_volume_noncompact_checked": False,
        "scope": (
            "Exact compact incidence plus complete local Lanner/elliptic "
            "label cores for one fixed a=1 source. Only terminal INFEASIBLE "
            "after all full-source core cuts closes the source."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_orbit", type=int, choices=range(4))
    parser.add_argument("--max-cores", type=int)
    parser.add_argument("--seconds-per-solve", type=float)
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--resume", type=Path)
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    seeded_iterations: Sequence[dict[str, object]] = ()
    if args.resume is not None:
        seed = json.loads(args.resume.read_text(encoding="utf-8"))
        if seed.get("source_orbit_index") != args.source_orbit:
            raise ValueError("resume checkpoint belongs to another source orbit")
        seeded_iterations = seed["iterations"]
    result = run_core_cegar(
        args.source_orbit,
        args.max_cores,
        args.seconds_per_solve,
        seeded_iterations,
    )
    if args.write is not None:
        args.write.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.summary_only:
        result = {
            "source_orbit_index": result["source_orbit_index"],
            "excluded_local_core_orbit_count": result[
                "excluded_local_core_orbit_count"
            ],
            "core_summaries": [
                {
                    "iteration": row["iteration"],
                    "representative_digest": row["representative_digest"],
                    "label_stage_state_counts": row.get(
                        "label_stage_state_counts"
                    ),
                    "core_A5_count": len(
                        row["local_core"]["core_additional_A5s"]
                    )
                    if row["local_core"] is not None
                    else None,
                    "core_actual_vertex_count": len(
                        row["local_core"]["core_actual_vertices"]
                    )
                    if row["local_core"] is not None
                    else None,
                    "full_core_orbit_size": row.get("full_core_orbit_size"),
                    "solver_wall_time_seconds": row[
                        "solver_wall_time_seconds"
                    ],
                }
                for row in result["iterations"]
            ],
            "terminal_solver_status": result["terminal_solver_status"],
            "stopped_reason": result["stopped_reason"],
            "fixed_source_compact_coxeter_no_go": result[
                "fixed_source_compact_coxeter_no_go"
            ],
            "finite_volume_noncompact_checked": False,
        }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
