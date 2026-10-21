#!/usr/bin/env python3
"""Search for a literal source-29 exact-shadow seed passing all strict links.

SAT/FEASIBLE produces a replayable incidence seed.  INFEASIBLE is reported
only for the explicitly fixed branch.  UNKNOWN is always exploratory and is
never promoted to mathematical evidence.  The model intentionally stops
before blocker purity, facets/ridges, Coxeter labels, or Gram geometry.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
import resource
import sys

from ortools.sat.python import cp_model

from scripts.certify_h9_14_record0_a2_open_strong_audit import (
    _shadow_counts,
    _strict_failures,
)
from scripts.explore_h9_14_record0_a2_dual_facets import mask
from scripts.explore_h9_14_record0_a2_open_primal_strong import (
    N,
    _install_strict_link_constraints,
    independent,
    source_problem,
)


RSS_LIMIT_BYTES = 3 * (1 << 28)
TARGETS = {6: 147, 7: 420, 8: 629}
T012 = mask((0, 1, 2))
FREE = tuple(range(5, 13))
F567 = (5, 6, 7)
PATTERNS = {
    "zero": (),
    "227": (227,),
    "229": (229,),
    "227-229": (227, 229),
    "229-230": (229, 230),
    "227-229-230": (227, 229, 230),
}


def _rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def _rss_guard() -> None:
    if _rss_bytes() > RSS_LIMIT_BYTES:
        raise MemoryError("strict-shadow model exceeded 0.75 GiB")


def build(branch: str) -> tuple[cp_model.CpModel, dict[str, object]]:
    problem = source_problem(29)
    candidates = tuple(problem["candidates"])
    cindex = {value: index for index, value in enumerate(candidates)}
    model = cp_model.CpModel()
    selected = tuple(model.new_bool_var(f"a5_{value}") for value in candidates)
    model.add(sum(selected) == int(problem["addition_count"]))
    for vertex in range(N):
        degree = sum(
            selected[index]
            for index, value in enumerate(candidates)
            if value >> vertex & 1
        )
        model.add(degree >= int(problem["lower"][vertex]))
        model.add(degree <= int(problem["upper"][vertex]))

    # These variables have the same polarity as the primal model: 1 means
    # the face survives.  The strict triple-link formula sums surviving
    # six-faces, so passing missing-face indicators here would be a fatal
    # polarity error.
    face_variables: dict[int, dict[int, cp_model.IntVar]] = {}
    base_low = (*problem["fours"],)
    for order, target in TARGETS.items():
        rows: dict[int, cp_model.IntVar] = {}
        for chosen in combinations(range(N), order):
            face = mask(chosen)
            if not independent(face) or any(face & row == row for row in base_low):
                continue
            destroyers = tuple(
                selected[index]
                for index, five in enumerate(candidates)
                if face & five == five
            )
            survives = model.new_bool_var(f"survives_{order}_{face}")
            if destroyers:
                for variable in destroyers:
                    model.add_implication(survives, variable.Not())
                model.add_bool_or([survives, *destroyers])
            else:
                model.add(survives == 1)
            rows[face] = survives
        model.add(sum(rows.values()) == len(rows) - target)
        face_variables[order] = rows
        _rss_guard()

    strict = _install_strict_link_constraints(
        model,
        problem,
        candidates,
        selected,
        face_variables[6],
    )

    # Analytic triple-012 cut: no selected A5 may contain 012, and the 168
    # noncontaining killers are either all zero or supported exactly on 567.
    containing012 = tuple(value for value in candidates if value & T012 == T012)
    killer_fibres = {}
    for free_triple in combinations(FREE, 3):
        six = T012 | mask(free_triple)
        killer_fibres[free_triple] = tuple(
            value
            for value in candidates
            if value & T012 != T012 and value & six == value
        )
    for value in containing012:
        model.add(selected[cindex[value]] == 0)
    pattern = set(PATTERNS[branch])
    for free_triple, values in killer_fibres.items():
        for value in values:
            model.add(
                selected[cindex[value]]
                == int(free_triple == F567 and value in pattern)
            )

    metadata = {
        "problem": problem,
        "candidates": candidates,
        "selected": selected,
        "face_variable_counts": {
            str(order): len(rows) for order, rows in face_variables.items()
        },
        "strict_constraint_counts": strict,
        "branch": branch,
        "fixed_pattern": sorted(pattern),
    }
    _rss_guard()
    return model, metadata


def solve(branch: str, seconds: float) -> dict[str, object]:
    model, metadata = build(branch)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 290123
    solver.parameters.max_memory_in_mb = 600
    status = solver.solve(model)
    status_name = solver.status_name(status)
    result: dict[str, object] = {
        "branch": branch,
        "status": status_name,
        "terminal": status in (cp_model.OPTIMAL, cp_model.INFEASIBLE),
        "literal_strict_seed_found": status in (cp_model.OPTIMAL, cp_model.FEASIBLE),
        "wall_time_seconds": solver.wall_time,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "face_variable_counts": metadata["face_variable_counts"],
        "strict_constraint_counts": metadata["strict_constraint_counts"],
        "peak_rss_bytes": _rss_bytes(),
        "scope": "exact shadows plus vertex/edge/triple strict links only",
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        candidates = metadata["candidates"]
        variables = metadata["selected"]
        additions = tuple(
            value
            for value, variable in zip(candidates, variables, strict=True)
            if solver.boolean_value(variable)
        )
        problem = metadata["problem"]
        shadows, missing_sixes = _shadow_counts(problem, additions)
        failures = _strict_failures(problem, additions, missing_sixes)
        if shadows != {"6": 147, "7": 420, "8": 629} or failures:
            raise AssertionError("solver seed failed independent project replay")
        result.update(
            {
                "additional_A5_masks": list(additions),
                "scheduled_missing_face_counts": shadows,
                "strict_failure_count": len(failures),
            }
        )
    _rss_guard()
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--branch", choices=tuple(PATTERNS), default="zero")
    parser.add_argument("--seconds", type=float, default=60.0)
    args = parser.parse_args()
    print(json.dumps(solve(args.branch, args.seconds), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
