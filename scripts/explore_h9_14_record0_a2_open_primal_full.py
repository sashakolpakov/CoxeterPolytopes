#!/usr/bin/env python3
"""Specialized full compact-incidence model for open Record-0/a=2 sources.

The variables are the additional independent missing five-sets.  Relative to
the lighter primal scout, this model also requires every set avoiding the
chosen missing faces to extend to a scheduled rank-nine facet, fixes the full
Dehn--Sommerville face vector, and installs the strict g/Macaulay conditions
for every vertex, edge, and triple link.  It is a CP-SAT translation of the
fixed-source consequences in the production facet-selection builder.

SAT remains incidence feasibility only.  UNSAT is an exact disposition of
this necessary model, but must be replayed independently before promotion.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from math import comb
from typing import Sequence

from ortools.sat.python import cp_model

from coxeter_polytopes.combinatorics.missing_edges import macaulay_raise
from scripts.explore_h9_14_record0_a1_dual_facets import _move
from scripts.explore_h9_14_record0_a2_dual_facets import (
    MISSING_EDGES,
    SOURCE_A4_FAMILIES,
    mask,
    source_stabilizer_generators,
)
from scripts.explore_h9_14_record0_a2_open_primal_strong import (
    FULL,
    N,
    OPEN_SOURCES,
    SCHEDULED_F_VECTOR,
    _and_not_equivalence,
    independent,
    source_problem,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    minimal_transversals,
    vertices,
)
from scripts.independent_h9_14_record14_a1_audit import (
    _add_involution_lex_leader,
)


def _face_literal(
    model: cp_model.CpModel,
    face_mask: int,
    candidates: Sequence[int],
    candidate_index: dict[int, int],
    selected: Sequence[cp_model.IntVar],
) -> cp_model.IntVar | int:
    """Return the exact indicator that ``face_mask`` avoids selected A5s."""

    contained = tuple(
        selected[candidate_index[five]]
        for five in (
            mask(chosen) for chosen in combinations(vertices(face_mask), 5)
        )
        if five in candidate_index
    )
    if not contained:
        return 1
    if face_mask.bit_count() == 5 and len(contained) == 1:
        return contained[0].Not()
    variable = model.new_bool_var(f"face_{face_mask}")
    _and_not_equivalence(model, variable, contained)
    return variable


def _enforce_true(
    model: cp_model.CpModel,
    clause: Sequence[cp_model.IntVar],
    condition: cp_model.IntVar | int,
) -> None:
    """Enforce one disjunction whenever a Boolean expression is true."""

    if isinstance(condition, int) and condition == 0:
        return
    if not clause:
        if isinstance(condition, int):
            model.add(0 == 1)
        else:
            model.add(0 == 1).only_enforce_if(condition)
        return
    constraint = model.add_bool_or(clause)
    if not isinstance(condition, int):
        constraint.only_enforce_if(condition)


def _strict_link_constraints(
    model: cp_model.CpModel,
    faces: dict[int, cp_model.IntVar | int],
) -> dict[str, object]:
    """Install the production strict local g/Macaulay inequalities."""

    strict_orders = {1: (2, 3, 4), 2: (2, 3), 3: (2,)}
    figure_counts: dict[str, int] = {}
    macaulay_table_count = 0
    for face_size, strict in strict_orders.items():
        figure_count = 0
        link_dimension = 9 - face_size
        for chosen in combinations(range(N), face_size):
            face_mask = mask(chosen)
            if not independent(face_mask):
                continue
            figure_count += 1
            h_expressions = []
            for index in range(link_dimension // 2 + 1):
                terms: list[object] = []
                for face_order in range(index + 1):
                    coefficient = (
                        (-1) ** (index - face_order)
                        * comb(
                            link_dimension - face_order,
                            link_dimension - index,
                        )
                    )
                    if face_order == 0:
                        terms.append(coefficient)
                        continue
                    supersets = []
                    outside = tuple(
                        vertex
                        for vertex in range(N)
                        if not face_mask & (1 << vertex)
                    )
                    for extension in combinations(outside, face_order):
                        superset = face_mask | mask(extension)
                        value = faces.get(superset, 0)
                        if not (isinstance(value, int) and value == 0):
                            supersets.append(value)
                    terms.append(coefficient * sum(supersets))
                h_expressions.append(sum(terms))

            maxima = [4]
            for degree in range(1, link_dimension // 2):
                maxima.append(macaulay_raise(maxima[-1], degree))
            g_variables = [
                model.new_int_var(0, maxima[index - 1], f"g_{face_mask}_{index}")
                for index in range(1, link_dimension // 2 + 1)
            ]
            for index, variable in enumerate(g_variables, start=1):
                model.add(
                    h_expressions[index]
                    == 1 + sum(g_variables[:index])
                )
                if index == 1 or index in strict:
                    model.add(variable >= 1)
            for degree in range(1, len(g_variables)):
                previous = g_variables[degree - 1]
                following = g_variables[degree]
                allowed = [
                    (value, next_value)
                    for value in range(maxima[degree - 1] + 1)
                    for next_value in range(macaulay_raise(value, degree) + 1)
                ]
                model.add_allowed_assignments([previous, following], allowed)
                macaulay_table_count += 1
        figure_counts[str(face_size)] = figure_count
    return {
        "figure_counts_by_face_size": figure_counts,
        "macaulay_table_constraint_count": macaulay_table_count,
    }


def build_model(
    source: int,
    *,
    symmetry_breaking: bool,
    strict_links: bool,
    purity_orders: Sequence[int],
) -> dict[str, object]:
    problem = source_problem(source)
    candidates = tuple(problem["candidates"])
    facets = tuple(problem["facet_complements"])
    candidate_index = {five: index for index, five in enumerate(candidates)}
    facet_index = {facet: index for index, facet in enumerate(facets)}

    model = cp_model.CpModel()
    selected = tuple(
        model.new_bool_var(f"missing5_{index}")
        for index in range(len(candidates))
    )
    model.add(sum(selected) == int(problem["addition_count"]))
    for vertex in range(N):
        members = tuple(
            selected[index]
            for index, five in enumerate(candidates)
            if five & (1 << vertex)
        )
        model.add(sum(members) >= problem["lower"][vertex])
        model.add(sum(members) <= problem["upper"][vertex])

    survives = tuple(
        model.new_bool_var(f"facet_{index}") for index in range(len(facets))
    )
    for target, facet in zip(survives, facets, strict=True):
        destroyers = tuple(
            selected[index]
            for index, five in enumerate(candidates)
            if five & facet == 0
        )
        _and_not_equivalence(model, target, destroyers)

    base_low = (*MISSING_EDGES, *problem["fours"])
    faces: dict[int, cp_model.IntVar | int] = {0: 1}
    face_counts: dict[str, int] = {}
    for order in range(1, 10):
        order_values = []
        for chosen in combinations(range(N), order):
            face_mask = mask(chosen)
            if any(face_mask & row == row for row in base_low):
                continue
            if order <= 4:
                value: cp_model.IntVar | int = 1
            elif order < 9:
                value = _face_literal(
                    model,
                    face_mask,
                    candidates,
                    candidate_index,
                    selected,
                )
            else:
                complement = FULL ^ face_mask
                value = survives[facet_index[complement]]
            faces[face_mask] = value
            order_values.append(value)
        model.add(sum(order_values) == SCHEDULED_F_VECTOR[order - 1])
        face_counts[str(order)] = len(order_values)

    # Every combinatorially allowed face in the requested range must be a
    # subset of at least one selected rank-nine facet.  The reverse direction
    # already follows from the exact facet-survival equivalences.
    extension_counts: dict[str, int] = {}
    for order in purity_orders:
        count = 0
        for face_mask, value in faces.items():
            if face_mask.bit_count() != order:
                continue
            witnesses = tuple(
                survives[index]
                for index, facet in enumerate(facets)
                if facet & face_mask == 0
            )
            _enforce_true(model, witnesses, value)
            count += 1
        extension_counts[str(order)] = count

    # Complements of order at most four are precisely the possible sets of
    # more than nine vertices.  Destroying every such base transversal is the
    # remaining Lanner-order closure branch.
    large_closure_count = 0
    for transversal in problem["small_transversals"]:
        destroyers = tuple(
            selected[index]
            for index, five in enumerate(candidates)
            if five & transversal == 0
        )
        if destroyers:
            model.add_bool_or(destroyers)
        else:
            model.add(0 == 1)
        large_closure_count += 1

    ridge_count = 0
    for six in (mask(chosen) for chosen in combinations(range(N), 6)):
        members = tuple(
            survives[index]
            for index, facet in enumerate(facets)
            if facet & six == facet
        )
        if not members:
            continue
        total = sum(members)
        model.add(total <= 2)
        for variable in members:
            model.add(total >= 2 * variable)
        ridge_count += 1

    link_audit = (
        _strict_link_constraints(model, faces)
        if strict_links
        else {
            "figure_counts_by_face_size": {},
            "macaulay_table_constraint_count": 0,
        }
    )

    generators, stabilizer_order = source_stabilizer_generators(problem["fours"])
    if not symmetry_breaking:
        generators = ()
    moved_counts = []
    for generator_index, action in enumerate(generators):
        permutation = tuple(
            candidate_index[_move(five, action)] for five in candidates
        )
        moved_counts.append(
            _add_involution_lex_leader(
                model, selected, permutation, f"lex_{generator_index}"
            )
        )

    return {
        "problem": problem,
        "model": model,
        "selected": selected,
        "survives": survives,
        "candidates": candidates,
        "facets": facets,
        "face_counts": face_counts,
        "extension_counts": extension_counts,
        "large_closure_count": large_closure_count,
        "ridge_count": ridge_count,
        "link_audit": link_audit,
        "stabilizer_order": stabilizer_order,
        "generator_count": len(generators),
        "moved_counts": moved_counts,
    }


def solve_source(
    source: int,
    *,
    seconds: float | None,
    workers: int,
    symmetry_breaking: bool,
    strict_links: bool,
    purity_orders: Sequence[int],
    cp_model_presolve: bool = True,
    internal_symmetry_level: int = 2,
    log_search_progress: bool = False,
) -> dict[str, object]:
    data = build_model(
        source,
        symmetry_breaking=symmetry_breaking,
        strict_links=strict_links,
        purity_orders=purity_orders,
    )
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.cp_model_presolve = cp_model_presolve
    solver.parameters.symmetry_level = internal_symmetry_level
    solver.parameters.log_search_progress = log_search_progress
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
        facets = tuple(
            facet
            for facet, variable in zip(
                data["facets"], data["survives"], strict=True
            )
            if solver.value(variable)
        )
        blocker = minimal_transversals((*data["problem"]["base"], *additions))
        target = {
            "additional_A5s": [list(vertices(value)) for value in additions],
            "facet_complements": [list(vertices(value)) for value in facets],
            "blocker_order_histogram": dict(
                sorted(Counter(value.bit_count() for value in blocker).items())
            ),
        }
    problem = data["problem"]
    return {
        "kind": "h9_14_record0_a2_open_primal_full_exploratory",
        "theorem_certificate": False,
        "record_index": 0,
        "schedule_a": 2,
        "source_orbit_index": source,
        "representative_A4s": [
            list(vertices(value)) for value in problem["fours"]
        ],
        "forced_A5_shadow_count": len(problem["shadow"]),
        "additional_A5_count": problem["addition_count"],
        "candidate_A5_count": len(data["candidates"]),
        "candidate_facet_complement_count": len(data["facets"]),
        "scheduled_face_candidate_counts": data["face_counts"],
        "purity_orders": list(purity_orders),
        "face_extension_constraint_counts": data["extension_counts"],
        "large_set_closure_constraint_count": data["large_closure_count"],
        "ridge_constraint_count": data["ridge_count"],
        "strict_links": strict_links,
        "cp_model_presolve": cp_model_presolve,
        "internal_symmetry_level": internal_symmetry_level,
        "strict_link_audit": data["link_audit"],
        "source_stabilizer_order": data["stabilizer_order"],
        "stabilizer_generator_count": data["generator_count"],
        "generator_moved_candidate_pair_counts": data["moved_counts"],
        "solver_status": solver.status_name(status),
        "solver_wall_time_seconds": solver.wall_time,
        "solver_branch_count": solver.num_branches,
        "solver_conflict_count": solver.num_conflicts,
        "target": target,
        "scope": (
            "Exact fixed-source necessary compact-incidence/full-link model. "
            "SAT is incidence feasibility only; UNKNOWN is not exclusion; "
            "UNSAT needs independent replay before theorem promotion."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=int, choices=OPEN_SOURCES)
    parser.add_argument("--seconds", type=float)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--no-symmetry", action="store_true")
    parser.add_argument("--no-strict-links", action="store_true")
    parser.add_argument("--no-presolve", action="store_true")
    parser.add_argument("--internal-symmetry-level", type=int, default=2)
    parser.add_argument("--log-search-progress", action="store_true")
    parser.add_argument(
        "--purity-orders", type=int, nargs="*", default=list(range(1, 9))
    )
    args = parser.parse_args()
    if any(order not in range(1, 9) for order in args.purity_orders):
        parser.error("purity orders must lie in 1,...,8")
    print(
        json.dumps(
            solve_source(
                args.source,
                seconds=args.seconds,
                workers=args.workers,
                symmetry_breaking=not args.no_symmetry,
                strict_links=not args.no_strict_links,
                purity_orders=tuple(args.purity_orders),
                cp_model_presolve=not args.no_presolve,
                internal_symmetry_level=args.internal_symmetry_level,
                log_search_progress=args.log_search_progress,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
