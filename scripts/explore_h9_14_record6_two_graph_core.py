#!/usr/bin/env python3
"""Probe which record-6 constraints already contradict the two-graph core.

This omits the expensive full face/link system.  It retains exactly the
ridge equations in their canonical two-graph form, the allowed facet-count
schedules, and optionally the independent-triple coverage constraints.  Its
purpose is conflict localization; neither SAT nor UNSAT is a frontier claim.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json


def _xor_all(*values):
    """Return parity of every operand; ``z3.Xor`` itself is binary."""

    import z3  # type: ignore[import-not-found]

    if not values:
        return z3.BoolVal(False)
    result = values[0]
    for value in values[1:]:
        result = z3.Xor(result, value)
    return result


def solve(*, coverage: bool, timeout_ms: int) -> dict[str, object]:
    import z3  # type: ignore[import-not-found]

    x = {
        (colour, first, second): z3.Bool(f"x_{colour}_{first}_{second}")
        for colour in range(3)
        for first, second in combinations(range(1, 11), 2)
    }

    def x_edge(colour: int, first: int, second: int):
        if first > second:
            first, second = second, first
        return z3.BoolVal(False) if first == 0 else x[colour, first, second]

    def y(colour: int, triple: tuple[int, int, int]):
        first, second, third = triple
        return _xor_all(
            x_edge(colour, first, second),
            x_edge(colour, first, third),
            x_edge(colour, second, third),
        )

    cut = {0: z3.BoolVal(False)}
    cut.update({index: z3.Bool(f"cut_{index}") for index in range(1, 11)})

    def z(first: int, second: int):
        return _xor_all(
            x_edge(0, first, second),
            x_edge(1, first, second),
            x_edge(2, first, second),
            cut[first],
            cut[second],
        )

    solver = z3.Solver()
    solver.set(timeout=timeout_ms)
    triples = tuple(combinations(range(11), 3))
    pairs = tuple(combinations(range(11), 2))
    y_value = {
        (colour, triple): y(colour, triple)
        for colour in range(3)
        for triple in triples
    }
    z_value = {pair: z(*pair) for pair in pairs}

    # Parity is built into the representation.  These upper bounds therefore
    # sharpen the even ridge sums to exactly zero or two.
    for colour in range(3):
        for four_set in combinations(range(11), 4):
            solver.add(
                z3.Sum(
                    [
                        z3.If(y_value[colour, triple], 1, 0)
                        for triple in combinations(four_set, 3)
                    ]
                )
                <= 2
            )
    for triple in triples:
        first, second, third = triple
        solver.add(
            z3.Sum(
                [z3.If(y_value[colour, triple], 1, 0) for colour in range(3)]
                + [
                    z3.If(z_value[first, second], 1, 0),
                    z3.If(z_value[first, third], 1, 0),
                    z3.If(z_value[second, third], 1, 0),
                ]
            )
            <= 2
        )

    facet_count = z3.Sum(
        [
            z3.If(value, 1, 0)
            for value in (*y_value.values(), *z_value.values())
        ]
    )
    solver.add(z3.Or(*(facet_count == value for value in range(130, 151, 2))))

    if coverage:
        for colour in range(3):
            for pair in pairs:
                pair_set = set(pair)
                solver.add(
                    z3.Or(
                        *(
                            y_value[colour, triple]
                            for triple in triples
                            if pair_set.isdisjoint(triple)
                        )
                    )
                )
        for missing_triple in triples:
            missing_set = set(missing_triple)
            solver.add(
                z3.Or(
                    *(
                        value
                        for (_, triple), value in y_value.items()
                        if missing_set.isdisjoint(triple)
                    ),
                    *(
                        value
                        for pair, value in z_value.items()
                        if missing_set.isdisjoint(pair)
                    ),
                )
            )

    assertion_count = len(solver.assertions())
    disposition = solver.check()
    result: dict[str, object] = {
        "status": str(disposition),
        "coverage": coverage,
        "assertion_count": assertion_count,
    }
    if disposition == z3.sat:
        model = solver.model()
        result["facet_count"] = model.eval(facet_count).as_long()
        result["y_counts"] = [
            sum(
                z3.is_true(model.eval(y_value[colour, triple]))
                for triple in triples
            )
            for colour in range(3)
        ]
        result["z_count"] = sum(
            z3.is_true(model.eval(value)) for value in z_value.values()
        )
    elif disposition == z3.unknown:
        result["reason_unknown"] = solver.reason_unknown()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage", action="store_true")
    parser.add_argument("--timeout-ms", type=int, default=300_000)
    args = parser.parse_args()
    print(
        json.dumps(
            solve(coverage=args.coverage, timeout_ms=args.timeout_ms),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
