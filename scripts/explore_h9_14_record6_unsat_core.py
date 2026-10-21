#!/usr/bin/env python3
"""Extract a diagnostic Z3 core for the exact record-6 two-graph model.

This is deliberately exploratory.  An SMT unsat core is useful for locating
the finite combinatorial contradiction, but it is not the solver-independent
certificate required by the theorem ledger.

The production builder does not expose its auxiliary literals.  We capture
the fully constructed solver without changing that builder, replay each
assertion under a tracking literal, and archive only the assertion indices and
text needed to translate a small core into a direct checker or hand lemma.
"""

from __future__ import annotations

import argparse
import json
from typing import Any


def captured_assertions(timeout_ms: int) -> tuple[Any, ...]:
    """Run the existing constructor once and return its exact assertions."""

    import z3  # type: ignore[import-not-found]

    from scripts import explore_h9_14_record6_two_graph as model

    original_solver = z3.Solver
    captured: list[Any] = []

    def factory(*args: Any, **kwargs: Any):
        solver = original_solver(*args, **kwargs)
        captured.append(solver)
        return solver

    z3.Solver = factory  # type: ignore[misc]
    try:
        initial = model.solve(timeout_ms)
    finally:
        z3.Solver = original_solver  # type: ignore[misc]
    if not captured:
        raise AssertionError("record-6 constructor did not create a solver")
    if len(captured) != 1:
        raise AssertionError(f"expected one captured solver, found {len(captured)}")
    if initial["status"] != "unsat":
        raise RuntimeError(
            "tracked replay requires the untracked model to be UNSAT; "
            f"observed {initial['status']}"
        )
    return tuple(captured[0].assertions())


def tracked_core(timeout_ms: int) -> dict[str, object]:
    import z3  # type: ignore[import-not-found]

    assertions = captured_assertions(timeout_ms)
    solver = z3.Solver()
    solver.set(timeout=timeout_ms, unsat_core=True)
    labels = []
    for index, assertion in enumerate(assertions):
        label = z3.Bool(f"track_{index}")
        labels.append(label)
        solver.assert_and_track(assertion, label)
    disposition = solver.check()
    result: dict[str, object] = {
        "status": str(disposition),
        "assertion_count": len(assertions),
    }
    if disposition == z3.unsat:
        core_indices = sorted(
            int(str(label).removeprefix("track_"))
            for label in solver.unsat_core()
        )
        result["core_size"] = len(core_indices)
        result["core_indices"] = core_indices
        result["core_assertions"] = [str(assertions[index]) for index in core_indices]
    elif disposition == z3.unknown:
        result["reason_unknown"] = solver.reason_unknown()
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout-ms", type=int, default=600_000)
    args = parser.parse_args()
    print(json.dumps(tracked_core(args.timeout_ms), sort_keys=True))


if __name__ == "__main__":
    main()
