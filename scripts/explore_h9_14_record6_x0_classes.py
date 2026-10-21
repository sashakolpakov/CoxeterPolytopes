#!/usr/bin/env python3
"""Canonical-X0 decomposition of the corrected record-6 two-graph model.

This remains exploratory until every canonical branch is decided and the
result receives an independent replay.  It fixes one minimum-cardinality
colour graph to each unlabeled triangle-free induced-2K2-free graph on ten
vertices.  Common isolate relabelling and colour permutation make the split
complete.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from pathlib import Path
import shutil
import subprocess
import sys
from time import monotonic

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
    build_facet_selection_constraint_system,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
RECORD_INDEX = 6


def _mask(vertices) -> int:
    return sum(1 << vertex for vertex in vertices)


def _adjacency(order: int, edges: tuple[tuple[int, int], ...]) -> tuple[int, ...]:
    result = [0] * order
    for first, second in edges:
        result[first] |= 1 << second
        result[second] |= 1 << first
    return tuple(result)


def _has_induced_2k2(adjacency: tuple[int, ...]) -> bool:
    for chosen in combinations(range(len(adjacency)), 4):
        degrees = [
            sum(bool(adjacency[vertex] & (1 << other)) for other in chosen)
            for vertex in chosen
        ]
        if sorted(degrees) == [1, 1, 1, 1]:
            return True
    return False


def _y_support(adjacency: tuple[int, ...]) -> frozenset[tuple[int, int, int]]:
    # Local vertex 0 is the switching base; graph vertices are local 1..10.
    result = set()
    for triple in combinations(range(11), 3):
        parity = 0
        for first, second in combinations(triple, 2):
            if first and second:
                parity ^= bool(adjacency[first - 1] & (1 << (second - 1)))
        if parity:
            result.add(triple)
    return frozenset(result)


def _special_a4_count(
    y_support: frozenset[tuple[int, int, int]],
) -> int:
    return sum(
        all(set(transversal) & set(triple) for triple in y_support)
        for transversal in combinations(range(11), 3)
    )


def _canonical_classes() -> list[dict[str, object]]:
    executable = shutil.which("geng")
    if executable is None:
        raise RuntimeError("nauty geng is required for the exploratory split")
    completed = subprocess.run(
        [executable, "-q", "-t", "10"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    retained = []
    for word in completed.stdout.splitlines():
        order, edges = decode_graph6(word)
        if order != 10:
            raise AssertionError("geng mixed graph orders")
        adjacency = _adjacency(order, edges)
        if _has_induced_2k2(adjacency):
            continue
        y_support = _y_support(adjacency)
        special_a4 = _special_a4_count(y_support)
        retained.append(
            {
                "graph6": word,
                "edges": edges,
                "adjacency": adjacency,
                "y_count": len(y_support),
                "special_a4_count": special_a4,
            }
        )
    if len(retained) != 308:
        raise AssertionError("record-6 canonical X census changed")
    # X0 is chosen to have minimum Y cardinality among the three colours, so
    # |Y0| <= 150/3.  The special-link interval independently gives A4<=3.
    return [
        row
        for row in retained
        if row["special_a4_count"] <= 3 and row["y_count"] <= 50
    ]


def solve(
    per_class_timeout_ms: int,
    start: int,
    stop: int | None,
    *,
    audit_sat: bool = False,
) -> dict:
    import z3  # type: ignore[import-not-found]

    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    record = [item for item in source["records"] if item["survives_second_layer"]][
        RECORD_INDEX
    ]
    _, raw_edges = decode_graph6(record["graph6"])
    specials = tuple(sorted({vertex for edge in raw_edges for vertex in edge}))
    if set(raw_edges) != set(combinations(specials, 2)):
        raise AssertionError("record 6 is no longer K3 plus isolates")
    isolates = tuple(vertex for vertex in range(14) if vertex not in specials)

    system = build_facet_selection_constraint_system(
        record["graph6"],
        record["ridge_support_core_masks"],
        record["dehn_sommerville_schedule"],
    )
    solver = system.solver
    solver.set(timeout=per_class_timeout_ms)
    facet = system.facet_literals

    x = {
        (colour, first, second): z3.Bool(f"record6_x_{colour}_{first}_{second}")
        for colour in range(3)
        for first, second in combinations(range(1, 11), 2)
    }

    def x_edge(colour: int, first: int, second: int):
        if first > second:
            first, second = second, first
        return z3.BoolVal(False) if first == 0 else x[colour, first, second]

    def xor_all(*values):
        result = values[0]
        for value in values[1:]:
            result = z3.Xor(result, value)
        return result

    y = {}
    for colour, omitted in enumerate(specials):
        pair = tuple(vertex for vertex in specials if vertex != omitted)
        for triple in combinations(range(11), 3):
            first, second, third = triple
            expression = xor_all(
                x_edge(colour, first, second),
                x_edge(colour, first, third),
                x_edge(colour, second, third),
            )
            complement = _mask(
                (*pair, *(isolates[index] for index in triple))
            )
            solver.add(facet[complement] == expression)
            y[colour, triple] = expression
    cut = {0: z3.BoolVal(False)}
    cut.update(
        {index: z3.Bool(f"record6_cut_{index}") for index in range(1, 11)}
    )
    for first, second in combinations(range(11), 2):
        complement = _mask(
            (*specials, isolates[first], isolates[second])
        )
        solver.add(
            facet[complement]
            == xor_all(
                x_edge(0, first, second),
                x_edge(1, first, second),
                x_edge(2, first, second),
                cut[first],
                cut[second],
            )
        )

    y_counts = [
        z3.Sum([z3.If(y[colour, triple], 1, 0) for triple in combinations(range(11), 3)])
        for colour in range(3)
    ]
    solver.add(y_counts[0] <= y_counts[1], y_counts[0] <= y_counts[2])

    classes = _canonical_classes()
    upper = len(classes) if stop is None else min(stop, len(classes))
    classes = classes[start:upper]
    outcomes = []
    started = monotonic()
    for offset, row in enumerate(classes, start=start):
        adjacency = row["adjacency"]
        assumptions = []
        for first, second in combinations(range(1, 11), 2):
            value = bool(adjacency[first - 1] & (1 << (second - 1)))
            variable = x[0, first, second]
            assumptions.append(variable if value else z3.Not(variable))
        # ``push()`` on this large mixed Boolean/LIA instance itself triggers
        # an unbounded LRA propagation pass. Assumption solving is equivalent
        # for the fixed X0 branch and keeps the per-branch timeout effective.
        disposition = solver.check(*assumptions)
        outcome = {
            "class_index": offset,
            "graph6": row["graph6"],
            "y_count": row["y_count"],
            "special_a4_count": row["special_a4_count"],
            "status": str(disposition),
        }
        if disposition == z3.sat:
            model = solver.model()
            selected = tuple(
                complement
                for complement in system.core
                if z3.is_true(model.eval(facet[complement], model_completion=True))
            )
            outcome["selected_facet_complement_masks"] = list(selected)
            # Preserve the raw model before invoking the intentionally much
            # stronger all-links homology audit: a malformed SAT witness can
            # make that audit far slower than the bounded solver branch.
            print(
                json.dumps(
                    {
                        "raw_sat_witness": True,
                        "class_index": offset,
                        "graph6": row["graph6"],
                        "selected_facet_complement_masks": list(selected),
                    },
                    sort_keys=True,
                ),
                file=sys.stderr,
                flush=True,
            )
            if audit_sat:
                audit = audit_global_facet_selection(
                    record["graph6"],
                    record["ridge_support_core_masks"],
                    record["dehn_sommerville_schedule"],
                    selected,
                )
                outcome["audit"] = audit.to_record()
            outcomes.append(outcome)
            break
        if disposition == z3.unknown:
            outcome["reason_unknown"] = solver.reason_unknown()
        outcomes.append(outcome)
        if (offset - start + 1) % 10 == 0:
            print(
                json.dumps(
                    {
                        "progress": offset + 1,
                        "branch_count": upper,
                        "last_status": str(disposition),
                    },
                    sort_keys=True,
                ),
                file=sys.stderr,
                flush=True,
            )
    return {
        "record_index": RECORD_INDEX,
        "graph6": record["graph6"],
        "correct_nested_xor": True,
        "raw_triangle_free_class_count": 12_172,
        "triangle_free_induced_2k2_free_class_count": 308,
        "admissible_minimum_colour_class_count": len(_canonical_classes()),
        "range": [start, upper],
        "per_class_timeout_ms": per_class_timeout_ms,
        "elapsed_seconds": monotonic() - started,
        "outcomes": outcomes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout-ms", type=int, default=5_000)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--audit-sat", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            solve(
                args.timeout_ms,
                args.start,
                args.stop,
                audit_sat=args.audit_sat,
            ),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
