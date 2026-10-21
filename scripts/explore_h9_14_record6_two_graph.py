#!/usr/bin/env python3
"""Solve H9/14 record 6 through its exact three-two-graph reduction.

For the missing graph ``K3 + 11I``, every possible five-complement is either
two special vertices plus three isolates, or all three specials plus two
isolates.  The ridge equations imply that each of the three 3-uniform
families is a two-graph.  After fixing one isolate as switching base it is the
triangle-coboundary of an ordinary graph.  The remaining two-complements are
the xor of those three graphs and a cut.

This script adds that *equivalent redundant formulation* to the complete Z3
facet-selection system.  SAT witnesses receive the independent exact global
audit.  Backend UNSAT/UNKNOWN is exploratory and is not a final theorem
certificate.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
    build_facet_selection_constraint_system,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
RECORD_INDEX = 6


def _mask(vertices: tuple[int, ...]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _xor_all(*values):
    """Return parity of every operand; ``z3.Xor`` itself is binary."""

    import z3  # type: ignore[import-not-found]

    if not values:
        return z3.BoolVal(False)
    result = values[0]
    for value in values[1:]:
        result = z3.Xor(result, value)
    return result


def solve(timeout_ms: int) -> dict[str, object]:
    import z3  # type: ignore[import-not-found]

    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    record = [item for item in source["records"] if item["survives_second_layer"]][
        RECORD_INDEX
    ]
    graph6 = record["graph6"]
    _, edges = decode_graph6(graph6)
    specials = tuple(sorted({vertex for edge in edges for vertex in edge}))
    if len(specials) != 3 or set(edges) != set(combinations(specials, 2)):
        raise AssertionError("record 6 is no longer K3 plus isolates")
    isolates = tuple(vertex for vertex in range(14) if vertex not in specials)
    if len(isolates) != 11:
        raise AssertionError("record 6 isolate count changed")

    system = build_facet_selection_constraint_system(
        graph6,
        record["ridge_support_core_masks"],
        record["dehn_sommerville_schedule"],
        timeout_ms=timeout_ms,
    )
    solver = system.solver
    facet = system.facet_literals

    # Canonical switching representatives: x_c(0,i)=0 and
    # y_c(0,i,j)=x_c(i,j).  Local isolate labels are 0,...,10.
    x: dict[tuple[int, int, int], object] = {}
    for colour in range(3):
        for first, second in combinations(range(1, 11), 2):
            x[colour, first, second] = z3.Bool(
                f"twograph_x_{colour}_{first}_{second}"
            )

    def x_edge(colour: int, first: int, second: int):
        if first > second:
            first, second = second, first
        if first == 0:
            return z3.BoolVal(False)
        return x[colour, first, second]

    def y_expr(colour: int, triple: tuple[int, int, int]):
        first, second, third = triple
        return _xor_all(
            x_edge(colour, first, second),
            x_edge(colour, first, third),
            x_edge(colour, second, third),
        )

    y_literals: dict[tuple[int, tuple[int, int, int]], object] = {}
    for colour, omitted_special in enumerate(specials):
        special_pair = tuple(
            vertex for vertex in specials if vertex != omitted_special
        )
        for local_triple in combinations(range(11), 3):
            global_triple = tuple(isolates[index] for index in local_triple)
            complement = _mask(tuple(sorted((*special_pair, *global_triple))))
            literal = facet.get(complement)
            if literal is None:
                raise AssertionError("record-6 pair/triple cover left the core")
            expression = y_expr(colour, local_triple)
            solver.add(literal == expression)
            y_literals[colour, local_triple] = literal

    # The xor of Z with the three graph representatives has zero triangle
    # coboundary, hence is a cut.  Fix the cut gauge p_0=0.
    cut = {0: z3.BoolVal(False)}
    cut.update({index: z3.Bool(f"cut_{index}") for index in range(1, 11)})
    z_literals: dict[tuple[int, int], object] = {}
    for first, second in combinations(range(11), 2):
        global_pair = (isolates[first], isolates[second])
        complement = _mask(tuple(sorted((*specials, *global_pair))))
        literal = facet.get(complement)
        if literal is None:
            raise AssertionError("record-6 triple/pair cover left the core")
        solver.add(
            literal
            == _xor_all(
                x_edge(0, first, second),
                x_edge(1, first, second),
                x_edge(2, first, second),
                cut[first],
                cut[second],
            )
        )
        z_literals[first, second] = literal

    assertion_count = len(solver.assertions())
    disposition = solver.check()
    result: dict[str, object] = {
        "record_index": RECORD_INDEX,
        "graph6": graph6,
        "backend": "z3_with_exact_two_graph_extension",
        "backend_version": z3.get_version_string(),
        "assertion_count": assertion_count,
        "status": str(disposition),
    }
    if disposition == z3.sat:
        model = solver.model()
        selected = tuple(
            complement
            for complement in system.core
            if z3.is_true(model.eval(facet[complement], model_completion=True))
        )
        audit = audit_global_facet_selection(
            graph6,
            record["ridge_support_core_masks"],
            record["dehn_sommerville_schedule"],
            selected,
        )
        result["selected_facet_complement_masks"] = list(selected)
        result["audit"] = audit.to_record()
    elif disposition == z3.unknown:
        result["reason_unknown"] = solver.reason_unknown()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout-ms", type=int, default=300_000)
    args = parser.parse_args()
    print(json.dumps(solve(args.timeout_ms), sort_keys=True))


if __name__ == "__main__":
    main()
