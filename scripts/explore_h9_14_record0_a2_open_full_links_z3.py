#!/usr/bin/env python3
"""Strict-link/full-face Z3 audit for the four open Record-0/a=2 sources.

This script composes the production facet-selection constraint system with
an exact choice of the two missing independent four-sets.  Thus it includes
the complete Dehn--Sommerville face vector, ridge degree 0/2, all low-face
witnesses, Lanner-order closure, and the strict local g/Macaulay conditions
for every vertex, edge, and triple link.  Generator lex constraints are a
lossless source-stabilizer slice and can be disabled for diagnosis.

An UNSAT answer is an exact disposition of this necessary incidence model;
it is not by itself the independent certificate required for a theorem.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
from typing import Sequence

import z3

from coxeter_polytopes.combinatorics.facet_completion import (
    build_facet_selection_constraint_system,
)
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest
from scripts.explore_h9_14_record0_a1_dual_facets import _move
from scripts.explore_h9_14_record0_a2_dual_facets import (
    SOURCE_A4_FAMILIES,
    mask,
    source_stabilizer_generators,
)


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = ROOT / "data" / "enumerations" / "h9_14_nerve_support.json"
GRAPH6 = "M?????????????o??"
OPEN_SOURCES = (24, 29, 30, 32)
SCHEDULE_A = 2
FULL = (1 << 14) - 1


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _record0() -> dict[str, object]:
    document = json.loads(SUPPORT_PATH.read_text(encoding="utf-8"))
    verify_support_manifest(document)
    survivors = [
        row for row in document["records"] if row["survives_second_layer"]
    ]
    row = survivors[0]
    if row["graph6"] != GRAPH6:
        raise AssertionError("Record-0 support graph changed")
    return row


def _add_involution_lex_leader(
    solver: z3.Solver,
    variables: Sequence[z3.BoolRef],
    permutation: Sequence[int],
) -> int:
    """Impose ``variables <=lex permutation(variables)`` exactly."""

    if sorted(permutation) != list(range(len(variables))):
        raise AssertionError("facet action is not a permutation")
    if any(permutation[permutation[index]] != index for index in range(len(variables))):
        raise AssertionError("source-stabilizer generator is not involutory")
    representatives = tuple(
        index for index, partner in enumerate(permutation) if index < partner
    )
    prefix: z3.BoolRef = z3.BoolVal(True)
    for index in representatives:
        partner = permutation[index]
        solver.add(
            z3.Implies(
                prefix,
                z3.Or(z3.Not(variables[index]), variables[partner]),
            )
        )
        prefix = z3.And(prefix, variables[index] == variables[partner])
    return len(representatives)


def build_source(
    source: int,
    *,
    timeout_ms: int | None,
    symmetry_breaking: bool,
) -> tuple[z3.Solver, dict[str, object]]:
    if source not in OPEN_SOURCES:
        raise ValueError(f"source must lie in {OPEN_SOURCES}")
    support = _record0()
    schedule = next(
        row
        for row in support["dehn_sommerville_schedule"]
        if row["missing_four_count"] == SCHEDULE_A
    )
    system = build_facet_selection_constraint_system(
        GRAPH6,
        support["ridge_support_core_masks"],
        [schedule],
        timeout_ms=timeout_ms,
    )
    solver = system.solver
    fours = tuple(mask(face) for face in SOURCE_A4_FAMILIES[source])

    # A four-set is a face exactly when some selected facet complement is
    # disjoint from it.  Fix both source representatives to be nonfaces.
    killed_literals = []
    for four in fours:
        witnesses = [
            literal
            for complement, literal in system.facet_literals.items()
            if complement & four == 0
        ]
        if not witnesses:
            raise AssertionError("source four-set has no support witness")
        solver.add(z3.Not(z3.Or(witnesses)))
        killed_literals.append(len(witnesses))

    generator_count = 0
    moved_pair_counts: list[int] = []
    stabilizer_order = 1
    if symmetry_breaking:
        generators, stabilizer_order = source_stabilizer_generators(fours)
        core_index = {
            complement: index for index, complement in enumerate(system.core)
        }
        for action in generators:
            permutation = tuple(
                core_index[_move(complement, action)]
                for complement in system.core
            )
            moved_pair_counts.append(
                _add_involution_lex_leader(
                    solver, system.variables, permutation
                )
            )
        generator_count = len(generators)
    else:
        _, stabilizer_order = source_stabilizer_generators(fours)

    return solver, {
        "source_orbit_index": source,
        "representative_A4s": [list(face) for face in SOURCE_A4_FAMILIES[source]],
        "support_core_count": len(system.core),
        "source_four_witness_counts": killed_literals,
        "assertion_count": len(solver.assertions()),
        "source_stabilizer_order": stabilizer_order,
        "symmetry_breaking": symmetry_breaking,
        "stabilizer_generator_count": generator_count,
        "generator_moved_facet_pair_counts": moved_pair_counts,
        "variables": system.variables,
        "core": system.core,
    }


def solve_source(
    source: int,
    *,
    timeout_ms: int | None,
    symmetry_breaking: bool,
) -> dict[str, object]:
    start = perf_counter()
    solver, data = build_source(
        source,
        timeout_ms=timeout_ms,
        symmetry_breaking=symmetry_breaking,
    )
    build_seconds = perf_counter() - start
    check_start = perf_counter()
    status = solver.check()
    check_seconds = perf_counter() - check_start
    selected: list[int] | None = None
    if status == z3.sat:
        model = solver.model()
        selected = [
            complement
            for complement, variable in zip(
                data["core"], data["variables"], strict=True
            )
            if z3.is_true(model.eval(variable, model_completion=True))
        ]
    reason_unknown = solver.reason_unknown() if status == z3.unknown else None
    del data["variables"]
    del data["core"]
    return {
        **data,
        "solver_status": str(status),
        "reason_unknown": reason_unknown,
        "build_seconds": build_seconds,
        "check_seconds": check_seconds,
        "selected_facet_complement_masks": selected,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", nargs="*", type=int)
    parser.add_argument("--timeout-ms", type=int, default=120_000)
    parser.add_argument("--no-timeout", action="store_true")
    parser.add_argument("--no-symmetry", action="store_true")
    args = parser.parse_args()
    sources = args.sources or list(OPEN_SOURCES)
    if any(source not in OPEN_SOURCES for source in sources):
        parser.error(f"sources must lie in {OPEN_SOURCES}")
    timeout_ms = None if args.no_timeout else args.timeout_ms
    rows = []
    for source in sources:
        row = solve_source(
            source,
            timeout_ms=timeout_ms,
            symmetry_breaking=not args.no_symmetry,
        )
        rows.append(row)
        print(
            f"source {source}: {row['solver_status']} "
            f"build={row['build_seconds']:.3f}s "
            f"check={row['check_seconds']:.3f}s",
            flush=True,
        )
    print(
        json.dumps(
            {
                "kind": "h9_14_record0_a2_open_full_links_z3_exploration",
                "theorem_certificate": False,
                "record_index": 0,
                "graph6": GRAPH6,
                "schedule_a": SCHEDULE_A,
                "source": {
                    "path": str(SUPPORT_PATH.relative_to(ROOT)),
                    "sha256": _sha256(SUPPORT_PATH),
                },
                "rows": rows,
                "scope": (
                    "Exact necessary compact-incidence plus strict-link model. "
                    "UNKNOWN is not an exclusion; UNSAT still requires an "
                    "independent reproducible certificate before promotion."
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
