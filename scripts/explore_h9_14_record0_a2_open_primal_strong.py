#!/usr/bin/env python3
"""Strong primal compact-incidence scout for four open Record-0/a=2 sources.

Variables select the additional minimal five-nonfaces, rather than dual
facets.  Exact facet-survival, linear ridge, scheduled intermediate-face
counts, and a source-stabilizer lex slice give an independent formulation of
the compact incidence overcover.  SAT is incidence feasibility only.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from math import comb
from typing import Iterable, Sequence

from ortools.sat.python import cp_model

from coxeter_polytopes.combinatorics.h9_14_record46 import _base_link_g3_g4
from coxeter_polytopes.combinatorics.missing_edges import macaulay_raise
from scripts.explore_h9_14_record0_a2_dual_facets import (
    MISSING_EDGES,
    SOURCE_A4_FAMILIES,
    mask,
    source_stabilizer_generators,
)
from scripts.explore_h9_14_record0_a1_dual_facets import _move
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    minimal_transversals,
    vertices,
)
from scripts.independent_h9_14_record14_a1_audit import (
    _add_involution_lex_leader,
)


N = 14
FULL = (1 << N) - 1
GRAPH6 = "M?????????????o??"
OPEN_SOURCES = (24, 29, 30, 32)
SCHEDULED_F_VECTOR = (14, 89, 341, 878, 1579, 1967, 1618, 792, 176)
SCHEDULED_MISSING_FIVES = 38


def independent(value: int) -> bool:
    return all(value & edge != edge for edge in MISSING_EDGES)


def hits_all(value: int, rows: Sequence[int]) -> bool:
    return all(value & row for row in rows)


def source_problem(source: int) -> dict[str, object]:
    if source not in OPEN_SOURCES:
        raise ValueError(f"source must lie in {OPEN_SOURCES}")
    fours = tuple(mask(face) for face in SOURCE_A4_FAMILIES[source])
    all_fives = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and independent(value)
    )
    shadow = tuple(
        five
        for five in all_fives
        if any(five & four == four for four in fours)
    )
    addition_count = SCHEDULED_MISSING_FIVES - len(shadow)
    base_g3, base_g4 = _base_link_g3_g4(GRAPH6)
    lower = []
    upper = []
    for vertex in range(N):
        q = sum(bool(four & (1 << vertex)) for four in fours)
        forced = sum(bool(five & (1 << vertex)) for five in shadow)
        previous = base_g3[vertex] - q
        before = base_g4[vertex] + 6 * q - forced
        lower.append(max(0, before - macaulay_raise(previous, 3)))
        upper.append(before - 1)
    candidates = tuple(
        five
        for five in all_fives
        if five not in shadow
        and all(
            not five & (1 << vertex) or upper[vertex] >= 1
            for vertex in range(N)
        )
    )
    base = (*MISSING_EDGES, *fours, *shadow)
    transversals = tuple(
        value for value in range(1, FULL + 1) if hits_all(value, base)
    )
    return {
        "source": source,
        "fours": fours,
        "shadow": shadow,
        "addition_count": addition_count,
        "lower": tuple(lower),
        "upper": tuple(upper),
        "candidates": candidates,
        "base": base,
        "small_transversals": tuple(
            value for value in transversals if value.bit_count() <= 4
        ),
        "facet_complements": tuple(
            value for value in transversals if value.bit_count() == 5
        ),
    }


def _and_not_equivalence(
    model: cp_model.CpModel,
    target: cp_model.IntVar,
    destroyers: Sequence[cp_model.IntVar],
) -> None:
    """Install ``target <=> AND(NOT destroyers)`` exactly."""

    if not destroyers:
        model.add(target == 1)
        return
    for variable in destroyers:
        model.add_implication(target, variable.Not())
    model.add_bool_or([target, *destroyers])


def _base_link_g_vector(
    problem: dict[str, object], face: int
) -> tuple[int, ...]:
    """Return the link g-vector before selecting additional A5s."""

    face_size = face.bit_count()
    link_dimension = 9 - face_size
    shadow = set(problem["shadow"])
    outside = tuple(vertex for vertex in range(N) if not face & (1 << vertex))
    h_entries = []
    for index in range(link_dimension // 2 + 1):
        value = 0
        for face_order in range(index + 1):
            coefficient = (
                (-1) ** (index - face_order)
                * comb(
                    link_dimension - face_order,
                    link_dimension - index,
                )
            )
            if face_order == 0:
                count = 1
            else:
                count = 0
                for extension in combinations(outside, face_order):
                    superset = face | mask(extension)
                    if not independent(superset):
                        continue
                    if any(
                        superset & four == four for four in problem["fours"]
                    ):
                        continue
                    if superset.bit_count() == 5 and superset in shadow:
                        continue
                    count += 1
            value += coefficient * count
        h_entries.append(value)
    return tuple(
        h_entries[index] - h_entries[index - 1]
        for index in range(1, len(h_entries))
    )


def _install_strict_link_constraints(
    model: cp_model.CpModel,
    problem: dict[str, object],
    candidates: Sequence[int],
    selected: Sequence[cp_model.IntVar],
    face_six: dict[int, cp_model.IntVar],
) -> dict[str, object]:
    """Install edge/triple strict-link consequences in compact form.

    Vertex-link g3/g4 constraints are already the local incidence intervals
    used to form ``problem``.  For every edge, g3 is the pre-addition value
    minus the selected-A5 incidence.  For every triple, the analogous g2 is
    linear, while g3 is computed from the exact six-face indicators and is
    bounded by Macaulay's inequality.
    """

    incidence_counts = {"2": 0, "3": 0}
    triple_table_count = 0
    for face_size in (2, 3):
        for chosen in combinations(range(N), face_size):
            face = mask(chosen)
            if not independent(face):
                continue
            g_vector = _base_link_g_vector(problem, face)
            incidence_degree = 3 if face_size == 2 else 2
            before = g_vector[incidence_degree - 1]
            previous = g_vector[incidence_degree - 2]
            lower = max(
                0,
                before
                - macaulay_raise(previous, incidence_degree - 1),
            )
            upper = before - 1
            members = tuple(
                selected[index]
                for index, five in enumerate(candidates)
                if five & face == face
            )
            incidence = sum(members)
            model.add(incidence >= lower)
            model.add(incidence <= upper)
            incidence_counts[str(face_size)] += 1

            if face_size != 3:
                continue
            # The last triple-link entry g3 depends on global six-faces.
            # Recompute h2,h3 from their defining face-incidence formula;
            # g2 is simultaneously tied to the selected-A5 incidence above.
            link_dimension = 6
            outside = tuple(
                vertex for vertex in range(N) if not face & (1 << vertex)
            )
            n1 = 0
            for extension in combinations(outside, 1):
                superset = face | mask(extension)
                if independent(superset) and not any(
                    superset & four == four for four in problem["fours"]
                ):
                    n1 += 1
            # Before additions, every valid five-superset contributes; each
            # selected candidate containing the triple removes one.
            base_n2 = 0
            shadow = set(problem["shadow"])
            for extension in combinations(outside, 2):
                superset = face | mask(extension)
                if not independent(superset):
                    continue
                if any(
                    superset & four == four for four in problem["fours"]
                ):
                    continue
                if superset in shadow:
                    continue
                base_n2 += 1
            n2 = base_n2 - incidence
            n3 = sum(
                face_six.get(face | mask(extension), 0)
                for extension in combinations(outside, 3)
            )
            h1 = -link_dimension + n1
            h2 = comb(link_dimension, 2) - (link_dimension - 1) * n1 + n2
            h3 = (
                -comb(link_dimension, 3)
                + comb(link_dimension - 1, 2) * n1
                - (link_dimension - 2) * n2
                + n3
            )
            g2_upper = macaulay_raise(g_vector[0], 1)
            g2 = model.new_int_var(1, g2_upper, f"triple_g2_{face}")
            g3_upper = macaulay_raise(g2_upper, 2)
            g3 = model.new_int_var(0, g3_upper, f"triple_g3_{face}")
            model.add(g2 == h2 - h1)
            model.add(g3 == h3 - h2)
            allowed = [
                (value, following)
                for value in range(1, g2_upper + 1)
                for following in range(macaulay_raise(value, 2) + 1)
            ]
            model.add_allowed_assignments([g2, g3], allowed)
            triple_table_count += 1
    return {
        "edge_and_triple_incidence_constraint_counts": incidence_counts,
        "triple_g2_g3_macaulay_table_count": triple_table_count,
    }


def build_model(
    source: int,
    *,
    face_orders: Sequence[int] = (6, 7),
    symmetry_breaking: bool = True,
    strict_links: bool = False,
) -> dict[str, object]:
    problem = source_problem(source)
    candidates = tuple(problem["candidates"])
    facets = tuple(problem["facet_complements"])
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

    small_constraint_count = 0
    for transversal in problem["small_transversals"]:
        destroyers = tuple(
            selected[index]
            for index, five in enumerate(candidates)
            if not five & transversal
        )
        if destroyers:
            model.add_bool_or(destroyers)
        else:
            model.add(0 == 1)
        small_constraint_count += 1

    survives = tuple(
        model.new_bool_var(f"facet_{index}")
        for index in range(len(facets))
    )
    for target, facet in zip(survives, facets, strict=True):
        destroyers = tuple(
            selected[index]
            for index, five in enumerate(candidates)
            if not five & facet
        )
        _and_not_equivalence(model, target, destroyers)
    model.add(sum(survives) == SCHEDULED_F_VECTOR[-1])

    ridge_constraint_count = 0
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
        ridge_constraint_count += 1

    face_count_constraints = {}
    face_values_by_order: dict[int, dict[int, cp_model.IntVar]] = {}
    base_low = (*MISSING_EDGES, *problem["fours"])
    candidate_index = {five: index for index, five in enumerate(candidates)}
    for order in face_orders:
        variables = []
        order_values = {}
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
            variable = model.new_bool_var(f"face{order}_{face}")
            _and_not_equivalence(model, variable, contained)
            variables.append(variable)
            order_values[face] = variable
        model.add(sum(variables) == SCHEDULED_F_VECTOR[order - 1])
        face_count_constraints[str(order)] = len(variables)
        face_values_by_order[order] = order_values

    if strict_links and 6 not in face_values_by_order:
        raise ValueError("strict link constraints require scheduled order 6")
    strict_link_audit = (
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

    generators, stabilizer_order = source_stabilizer_generators(
        problem["fours"]
    )
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
        "candidate_fives": candidates,
        "facet_complements": facets,
        "small_transversal_constraint_count": small_constraint_count,
        "ridge_constraint_count": ridge_constraint_count,
        "face_count_constraints": face_count_constraints,
        "strict_links": strict_links,
        "strict_link_audit": strict_link_audit,
        "source_stabilizer_order": stabilizer_order,
        "stabilizer_generator_count": len(generators),
        "generator_moved_candidate_pair_counts": moved_counts,
    }


def solve(
    source: int,
    *,
    seconds: float | None,
    workers: int,
    face_orders: Sequence[int],
    symmetry_breaking: bool,
    strict_links: bool = False,
) -> dict[str, object]:
    data = build_model(
        source,
        face_orders=face_orders,
        symmetry_breaking=symmetry_breaking,
        strict_links=strict_links,
    )
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    if seconds is not None:
        solver.parameters.max_time_in_seconds = seconds
    status = solver.solve(data["model"])
    target = None
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        additions = tuple(
            five
            for five, variable in zip(
                data["candidate_fives"], data["selected"], strict=True
            )
            if solver.value(variable)
        )
        facets = tuple(
            facet
            for facet, variable in zip(
                data["facet_complements"], data["survives"], strict=True
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
            "pure_scheduled_facets": (
                len(facets) == SCHEDULED_F_VECTOR[-1]
                and Counter(value.bit_count() for value in blocker)
                == {5: SCHEDULED_F_VECTOR[-1]}
            ),
        }
    problem = data["problem"]
    return {
        "kind": "h9_14_record0_a2_open_primal_strong_exploratory",
        "theorem_certificate": False,
        "record_index": 0,
        "schedule_a": 2,
        "source_orbit_index": source,
        "representative_A4s": [
            list(vertices(value)) for value in problem["fours"]
        ],
        "forced_A5_shadow_count": len(problem["shadow"]),
        "additional_A5_count": problem["addition_count"],
        "candidate_A5_count": len(problem["candidates"]),
        "candidate_facet_complement_count": len(problem["facet_complements"]),
        "small_transversal_constraint_count": data[
            "small_transversal_constraint_count"
        ],
        "ridge_constraint_count": data["ridge_constraint_count"],
        "scheduled_face_orders_installed": list(face_orders),
        "face_count_variable_counts": data["face_count_constraints"],
        "strict_links": data["strict_links"],
        "strict_link_audit": data["strict_link_audit"],
        "source_stabilizer_order": data["source_stabilizer_order"],
        "stabilizer_generator_count": data["stabilizer_generator_count"],
        "generator_moved_candidate_pair_counts": data[
            "generator_moved_candidate_pair_counts"
        ],
        "solver_status": solver.status_name(status),
        "solver_wall_time_seconds": solver.wall_time,
        "solver_branch_count": solver.num_branches,
        "solver_conflict_count": solver.num_conflicts,
        "target": target,
        "scope": (
            "Exact compact-incidence overcover for one fixed source. SAT is "
            "incidence feasibility only; UNKNOWN is not an exclusion; no "
            "finite-volume noncompact inference is made."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=int, choices=OPEN_SOURCES)
    parser.add_argument("--seconds", type=float)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--face-orders", type=int, nargs="*", default=[6, 7])
    parser.add_argument("--no-symmetry", action="store_true")
    parser.add_argument("--strict-links", action="store_true")
    args = parser.parse_args()
    if any(order not in (6, 7, 8) for order in args.face_orders):
        parser.error("face orders must lie in 6,7,8")
    print(
        json.dumps(
            solve(
                args.source,
                seconds=args.seconds,
                workers=args.workers,
                face_orders=tuple(args.face_orders),
                symmetry_breaking=not args.no_symmetry,
                strict_links=args.strict_links,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
