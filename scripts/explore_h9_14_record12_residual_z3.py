"""Exploratory reduced residual solver for H9/14 record 12."""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

import z3

from coxeter_polytopes.combinatorics.facet_completion import (
    build_facet_selection_constraint_system,
)


ROOT = Path(__file__).resolve().parents[1]
A = {0, 1}
B = {12, 13}
I = tuple(range(2, 12))


def parity(values):
    result = z3.BoolVal(False)
    for value in values:
        result = z3.Xor(result, value)
    return result


def solve(sizes: tuple[int, int, int, int], phase: int, limit: int) -> None:
    source = json.loads(
        (ROOT / "data/enumerations/h9_14_nerve_support.json").read_text()
    )
    record = [item for item in source["records"] if item["survives_second_layer"]][12]
    system = build_facet_selection_constraint_system(
        record["graph6"],
        record["ridge_support_core_masks"],
        record["dehn_sommerville_schedule"][:8],
    )
    solver = system.solver
    cells = ((0, 0), (0, 1), (1, 0), (1, 1))
    p = {}
    q = {}
    y = {}
    cursor = 0
    for cell, size in zip(cells, sizes):
        for vertex in I[cursor : cursor + size]:
            p[vertex], q[vertex] = cell
            y[vertex] = cell[0] ^ cell[1] ^ phase
        cursor += size
    if cursor != 10:
        raise ValueError("cell sizes must sum to ten")
    x = {
        pair: z3.Bool(f"rx_{pair[0]}_{pair[1]}") for pair in combinations(I, 2)
    }
    w = {
        pair: z3.Bool(f"rz_{pair[0]}_{pair[1]}") for pair in combinations(I, 2)
    }

    def edge(table, first, second):
        return table[tuple(sorted((first, second)))]

    for complement, literal in system.facet_literals.items():
        vertices = {vertex for vertex in range(14) if complement >> vertex & 1}
        if A <= vertices and B <= vertices:
            isolate = next(iter(vertices - A - B))
            expression = z3.BoolVal(bool(y[isolate]))
        elif A <= vertices:
            triple = vertices - A
            special = triple & B
            isolate = sorted(triple - B)
            if special == {13}:
                expression = edge(x, *isolate)
            elif special == {12}:
                expression = z3.Xor(
                    edge(x, *isolate),
                    z3.BoolVal(bool(p[isolate[0]] ^ p[isolate[1]])),
                )
            elif not special:
                expression = parity(edge(x, *pair) for pair in combinations(isolate, 2))
            else:
                raise AssertionError((vertices, triple))
        elif B <= vertices:
            triple = vertices - B
            special = triple & A
            isolate = sorted(triple - A)
            if special == {1}:
                expression = edge(w, *isolate)
            elif special == {0}:
                expression = z3.Xor(
                    edge(w, *isolate),
                    z3.BoolVal(bool(q[isolate[0]] ^ q[isolate[1]])),
                )
            elif not special:
                expression = parity(edge(w, *pair) for pair in combinations(isolate, 2))
            else:
                raise AssertionError((vertices, triple))
        else:
            raise AssertionError("support core contains a non-cover")
        solver.add(literal == expression)

    count = 0
    while count < limit and solver.check() == z3.sat:
        model = solver.model()
        selected = tuple(
            complement
            for complement, literal in system.facet_literals.items()
            if z3.is_true(model.eval(literal, model_completion=True))
        )
        print(
            json.dumps(
                {
                    "index": count,
                    "facet_count": len(selected),
                    "selected": selected,
                    "x_edges": [list(pair) for pair, value in x.items() if z3.is_true(model.eval(value, model_completion=True))],
                    "z_edges": [list(pair) for pair, value in w.items() if z3.is_true(model.eval(value, model_completion=True))],
                }
            ),
            flush=True,
        )
        assignment = [
            value != z3.is_true(model.eval(value, model_completion=True))
            for value in (*x.values(), *w.values())
        ]
        solver.add(z3.Or(assignment))
        count += 1
    print(json.dumps({"count": count, "final": str(solver.check())}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", default="5,0,0,5")
    parser.add_argument("--phase", type=int, default=0)
    parser.add_argument("--limit", type=int, default=10)
    args = parser.parse_args()
    solve(tuple(map(int, args.sizes.split(","))), args.phase, args.limit)
