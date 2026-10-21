#!/usr/bin/env python3
"""Exploratory exact a=0 model stream for H9/14 record 30.

Eight A5 rows are ordered lexicographically and the eight isolate incidence
columns are ordered lexicographically.  This double-lex condition is a
lossless (though potentially colliding) quotient by row and isolate
permutations: a lexicographically least matrix in every finite orbit obeys
both orders.  Z3 is only a search backend; final theorem use requires a
solver-independent replay.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import subprocess
import tempfile
from time import monotonic

import z3


N = 14
ROWS = 8
ACTIVE = frozenset((0, 1, 2, 11, 12, 13))
ISOLATES = tuple(range(3, 11))
EDGES = ((0, 11), (1, 12), (2, 13))
EDGE_MASKS = tuple((1 << a) | (1 << b) for a, b in EDGES)
ROOT = Path(__file__).resolve().parents[1]
ORBIT_SOURCE = ROOT / "scripts/h9_14_record30_a0_orbit_forms.cpp"


def mask(values) -> int:
    return sum(1 << value for value in values)


def independent(value: int) -> bool:
    return all(value & edge != edge for edge in EDGE_MASKS)


def local_cap(chosen: tuple[int, ...]) -> int:
    active = sum(vertex in ACTIVE for vertex in chosen)
    if len(chosen) == 1:
        return {0: 7, 1: 3}[active]
    if len(chosen) == 2:
        return {0: 7, 1: 3, 2: 1}[active]
    return {0: 6, 1: 3, 2: 1, 3: 0}[active]


def build() -> tuple[z3.Solver, list[list[z3.BoolRef]]]:
    selected = [[z3.Bool(f"r{row}v{vertex}") for vertex in range(N)] for row in range(ROWS)]
    solver = z3.Solver()
    solver.set(timeout=300_000)
    for row in range(ROWS):
        solver.add(z3.PbEq([(selected[row][v], 1) for v in range(N)], 5))
        for a, b in EDGES:
            solver.add(z3.Or(z3.Not(selected[row][a]), z3.Not(selected[row][b])))

    # Every vertex/edge/triple local A5 incidence interval has lower end zero.
    for size in (1, 2, 3):
        for chosen in combinations(range(N), size):
            value = mask(chosen)
            if not independent(value):
                continue
            incidences = [z3.And(*(selected[row][v] for v in chosen)) for row in range(ROWS)]
            solver.add(z3.PbLe([(item, 1) for item in incidences], local_cap(chosen)))

    # Destroy every order-three/four transversal of the missing matching.
    for size in (3, 4):
        for chosen in combinations(range(N), size):
            value = mask(chosen)
            if not all(value & edge for edge in EDGE_MASKS):
                continue
            solver.add(z3.Or(*(
                z3.And(*(z3.Not(selected[row][v]) for v in chosen))
                for row in range(ROWS)
            )))

    five_covers = tuple(
        chosen for chosen in combinations(range(N), 5)
        if all(mask(chosen) & edge for edge in EDGE_MASKS)
    )
    facets = []
    for index, chosen in enumerate(five_covers):
        facet = z3.Bool(f"f{index}")
        solver.add(facet == z3.And(*(
            z3.Or(*(selected[row][v] for v in chosen))
            for row in range(ROWS)
        )))
        facets.append(facet)
    solver.add(z3.PbEq([(facet, 1) for facet in facets], 132))

    # A ridge is indexed by a six-set containing its five-set complements.
    cover_index = {mask(chosen): index for index, chosen in enumerate(five_covers)}
    for chosen in combinations(range(N), 6):
        six = mask(chosen)
        indices = [cover_index[six ^ (1 << v)] for v in chosen if six ^ (1 << v) in cover_index]
        if indices:
            degree = z3.Sum(*(z3.If(facets[index], 1, 0) for index in indices))
            solver.add(z3.Or(degree == 0, degree == 2))

    # Lossless row/isolate double-lex quotient.
    row_values = [z3.Sum(*(z3.If(selected[row][v], 1 << v, 0) for v in range(N))) for row in range(ROWS)]
    for row in range(ROWS - 1):
        solver.add(row_values[row] < row_values[row + 1])
    column_values = [z3.Sum(*(z3.If(selected[row][v], 1 << row, 0) for row in range(ROWS))) for v in ISOLATES]
    for left, right in zip(column_values, column_values[1:]):
        solver.add(left <= right)
    return solver, selected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-orbits", type=int, default=10000)
    args = parser.parse_args()
    solver, selected = build()
    started = monotonic()
    orbits = []
    status = "max_orbits"
    with tempfile.TemporaryDirectory(prefix="record30_a0_orbits_") as temporary:
        executable = Path(temporary) / "orbit_forms"
        subprocess.run(
            ["clang++", "-std=c++20", "-O3", "-Wall", "-Wextra", "-Werror",
             str(ORBIT_SOURCE), "-o", str(executable)], check=True, cwd=ROOT
        )
        while len(orbits) < args.max_orbits:
            result = solver.check()
            if result == z3.unsat:
                status = "exhausted"
                break
            if result != z3.sat:
                status = str(result)
                break
            model = solver.model()
            family = tuple(
                sum(z3.is_true(model.eval(selected[row][v], model_completion=True)) << v for v in range(N))
                for row in range(ROWS)
            )
            orbit = json.loads(subprocess.run(
                [str(executable), *map(str, family)], check=True, text=True,
                capture_output=True, cwd=ROOT
            ).stdout)
            forms = tuple(tuple(map(int, row)) for row in orbit["forms"])
            if not forms or family not in forms:
                raise AssertionError("orbit helper lost the Z3 double-lex model")
            representative = min(forms)
            if any(item["representative"] == representative for item in orbits):
                raise AssertionError("an already blocked full orbit returned")
            orbits.append({"representative": representative, "double_lex_form_count": len(forms)})
            for form in forms:
                solver.add(z3.Or(*(
                    selected[row][v] != bool(form[row] & (1 << v))
                    for row in range(ROWS) for v in range(N)
                )))
            print(json.dumps({"orbit": len(orbits), "representative": representative,
                              "double_lex_forms": len(forms),
                              "elapsed": round(monotonic() - started, 3)}), flush=True)
    print(json.dumps({
        "kind": "exploratory_h9_14_record30_a0_double_lex_models",
        "status": status,
        "orbit_count": len(orbits),
        "orbits": orbits,
        "elapsed_seconds": round(monotonic() - started, 6),
        "solver_independent_certificate": False,
    }), flush=True)


if __name__ == "__main__":
    main()
