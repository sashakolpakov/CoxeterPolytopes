#!/usr/bin/env python3
"""Staged missing-five shadow design for the four open Record-0/a=2 sources.

This deliberately omits facet/ridge/purity conditions and asks first whether
the additional A5 family can realize the scheduled six-, seven-, and
eight-face shadow sizes together with every strict vertex/edge/triple link
g/Macaulay inequality through triple g3.  UNSAT is therefore already a valid
necessary obstruction; SAT supplies an exact A5 seed for the full model.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json

from ortools.sat.python import cp_model

from scripts.explore_h9_14_record0_a1_dual_facets import _move
from scripts.explore_h9_14_record0_a2_dual_facets import (
    MISSING_EDGES,
    mask,
    source_stabilizer_generators,
)
from scripts.explore_h9_14_record0_a2_open_primal_strong import (
    N,
    OPEN_SOURCES,
    SCHEDULED_F_VECTOR,
    _and_not_equivalence,
    _install_strict_link_constraints,
    source_problem,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import vertices
from scripts.independent_h9_14_record14_a1_audit import (
    _add_involution_lex_leader,
)


# Frozen in ``certify_h9_14_record0_a2_open_strong_audit.py``.  It realizes
# the exact order-6/7/8 shadow counts but fails 48 strict edge/triple links.
# Supplying it as a CP-SAT hint is a search aid only, never a certificate.
SOURCE29_SHADOW_WITNESS = (
    158,
    3107,
    2087,
    5157,
    244,
    541,
    188,
    4141,
    1549,
    1174,
    1543,
    6179,
    1306,
    12836,
    14432,
    2466,
    6592,
    9496,
    11584,
    1558,
    4649,
)


def build_model(
    source: int,
    *,
    face_orders: tuple[int, ...],
    strict_links: bool,
    symmetry_breaking: bool,
    source29_hint: bool = False,
) -> dict[str, object]:
    problem = source_problem(source)
    candidates = tuple(problem["candidates"])
    candidate_index = {five: index for index, five in enumerate(candidates)}
    model = cp_model.CpModel()
    selected = tuple(
        model.new_bool_var(f"missing5_{index}")
        for index in range(len(candidates))
    )
    if source29_hint:
        if source != 29:
            raise ValueError("the archived shadow hint is defined only for source 29")
        chosen = set(SOURCE29_SHADOW_WITNESS)
        if not chosen <= set(candidates):
            raise AssertionError("the archived source-29 hint left the candidate ground")
        for five, variable in zip(candidates, selected, strict=True):
            model.add_hint(variable, int(five in chosen))
    model.add(sum(selected) == int(problem["addition_count"]))
    for vertex in range(N):
        members = tuple(
            selected[index]
            for index, five in enumerate(candidates)
            if five & (1 << vertex)
        )
        model.add(sum(members) >= problem["lower"][vertex])
        model.add(sum(members) <= problem["upper"][vertex])

    base_low = (*MISSING_EDGES, *problem["fours"])
    face_values_by_order: dict[int, dict[int, cp_model.IntVar]] = {}
    face_counts = {}
    for order in face_orders:
        order_values: dict[int, cp_model.IntVar] = {}
        for chosen in combinations(range(N), order):
            face = mask(chosen)
            if any(face & row == row for row in base_low):
                continue
            contained = tuple(
                selected[candidate_index[five]]
                for five in (
                    mask(subset) for subset in combinations(chosen, 5)
                )
                if five in candidate_index
            )
            value = model.new_bool_var(f"face{order}_{face}")
            _and_not_equivalence(model, value, contained)
            order_values[face] = value
        model.add(
            sum(order_values.values()) == SCHEDULED_F_VECTOR[order - 1]
        )
        face_values_by_order[order] = order_values
        face_counts[str(order)] = len(order_values)

    if strict_links and 6 not in face_values_by_order:
        raise ValueError("strict links require face order 6")
    link_audit = (
        _install_strict_link_constraints(
            model,
            problem,
            candidates,
            selected,
            face_values_by_order[6],
        )
        if strict_links
        else {
            "edge_and_triple_incidence_constraint_counts": {},
            "triple_g2_g3_macaulay_table_count": 0,
        }
    )

    generators, stabilizer_order = source_stabilizer_generators(problem["fours"])
    if not symmetry_breaking:
        generators = ()
    moved_counts = []
    for index, action in enumerate(generators):
        permutation = tuple(
            candidate_index[_move(five, action)] for five in candidates
        )
        moved_counts.append(
            _add_involution_lex_leader(
                model, selected, permutation, f"lex_{index}"
            )
        )
    return {
        "model": model,
        "problem": problem,
        "candidates": candidates,
        "selected": selected,
        "face_counts": face_counts,
        "strict_links": strict_links,
        "link_audit": link_audit,
        "stabilizer_order": stabilizer_order,
        "generator_count": len(generators),
        "moved_counts": moved_counts,
        "source29_hint": source29_hint,
    }


def solve_source(
    source: int,
    *,
    seconds: float | None,
    workers: int,
    face_orders: tuple[int, ...],
    strict_links: bool,
    symmetry_breaking: bool,
    presolve: bool,
    probing_level: int,
    internal_symmetry_level: int,
    source29_hint: bool = False,
) -> dict[str, object]:
    data = build_model(
        source,
        face_orders=face_orders,
        strict_links=strict_links,
        symmetry_breaking=symmetry_breaking,
        source29_hint=source29_hint,
    )
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.cp_model_presolve = presolve
    solver.parameters.cp_model_probing_level = probing_level
    solver.parameters.symmetry_level = internal_symmetry_level
    if seconds is not None:
        solver.parameters.max_time_in_seconds = seconds
    status = solver.solve(data["model"])
    target = None
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        additions = tuple(
            five
            for five, variable in zip(
                data["candidates"], data["selected"], strict=True
            )
            if solver.value(variable)
        )
        target = {
            "additional_A5_masks": list(additions),
            "additional_A5s": [list(vertices(value)) for value in additions],
        }
    problem = data["problem"]
    return {
        "kind": "h9_14_record0_a2_open_shadow_design_exploratory",
        "theorem_certificate": False,
        "source_orbit_index": source,
        "representative_A4s": [
            list(vertices(value)) for value in problem["fours"]
        ],
        "additional_A5_count": problem["addition_count"],
        "candidate_A5_count": len(data["candidates"]),
        "scheduled_face_orders": list(face_orders),
        "scheduled_face_candidate_counts": data["face_counts"],
        "strict_links": data["strict_links"],
        "strict_link_audit": data["link_audit"],
        "source_stabilizer_order": data["stabilizer_order"],
        "stabilizer_generator_count": data["generator_count"],
        "generator_moved_candidate_pair_counts": data["moved_counts"],
        "cp_model_presolve": presolve,
        "cp_model_probing_level": probing_level,
        "internal_symmetry_level": internal_symmetry_level,
        "archived_source29_shadow_witness_used_as_hint": data["source29_hint"],
        "solver_status": solver.status_name(status),
        "solver_wall_time_seconds": solver.wall_time,
        "solver_branch_count": solver.num_branches,
        "solver_conflict_count": solver.num_conflicts,
        "target": target,
        "scope": (
            "Exact necessary A5-shadow/link design only. SAT is a seed, not "
            "a compact polytope incidence; UNKNOWN is not exclusion."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=int, choices=OPEN_SOURCES)
    parser.add_argument("--seconds", type=float)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--face-orders", type=int, nargs="*", default=[6, 7, 8])
    parser.add_argument("--no-strict-links", action="store_true")
    parser.add_argument("--no-symmetry", action="store_true")
    parser.add_argument("--no-presolve", action="store_true")
    parser.add_argument("--probing-level", type=int, default=2)
    parser.add_argument("--internal-symmetry-level", type=int, default=0)
    parser.add_argument("--source29-hint", action="store_true")
    args = parser.parse_args()
    if any(order not in (6, 7, 8) for order in args.face_orders):
        parser.error("face orders must lie in 6,7,8")
    print(
        json.dumps(
            solve_source(
                args.source,
                seconds=args.seconds,
                workers=args.workers,
                face_orders=tuple(args.face_orders),
                strict_links=not args.no_strict_links,
                symmetry_breaking=not args.no_symmetry,
                presolve=not args.no_presolve,
                probing_level=args.probing_level,
                internal_symmetry_level=args.internal_symmetry_level,
                source29_hint=args.source29_hint,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
