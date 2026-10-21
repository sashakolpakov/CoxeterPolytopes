#!/usr/bin/env python3
"""Build, run, and freeze the no-SAT record-18 a2/add5 replay."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/independent_h9_14_record18_a2_add5_audit.cpp"
DEFAULT_OUTPUT = (
    ROOT / "data/enumerations/h9_14_record18_a2_add5_independent_audit.json"
)
REPRESENTATIVES = ((53, 54), (23, 27))


def run() -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="record18-a2-add5-") as directory:
        binary = Path(directory) / "audit"
        subprocess.run(
            ["c++", "-std=c++20", "-O3", str(SOURCE), "-o", str(binary)],
            cwd=ROOT,
            check=True,
        )
        cases = []
        for representative in REPRESENTATIVES:
            completed = subprocess.run(
                [str(binary), *map(str, representative)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            cases.append(json.loads(completed.stdout))
    expected = {
        (53, 54): {
            "incidence_symmetry_order": 2880,
            "distinct_state_counts": [1, 6, 81, 997, 17947, 215909],
            "search_node_count": 413064,
        },
        (23, 27): {
            "incidence_symmetry_order": 20160,
            "distinct_state_counts": [1, 4, 41, 421, 6641, 40659],
            "search_node_count": 87819,
        },
    }
    for case in cases:
        key = tuple(case["representative_A4_masks"])
        if key not in expected:
            raise AssertionError("unexpected a2/add5 representative")
        for field, value in expected[key].items():
            if case[field] != value:
                raise AssertionError(f"a2/add5 transcript changed at {key}/{field}")
        if (
            case["candidate_count"],
            case["base_small_transversal_count"],
            case["base_order_five_transversal_count"],
            case["nonempty_ridge_constraint_count"],
            case["scheduled_ridge_terminal_count"],
            case["literal_pure_target_count"],
            case["compact_incidence_unsat"],
        ) != (777, 56, 233, 630, 0, 0, True):
            raise AssertionError(f"a2/add5 conclusion changed at {key}")
    return {
        "kind": "independent_h9_14_record18_a2_add5_combined_audit",
        "backend": "standalone C++ incidence-matrix quotient; no SAT/SMT",
        "representative_count": len(cases),
        "cases": cases,
        "all_a2_add5_compact_incidence_branches_unsat": True,
        "compact_only": True,
        "finite_volume_noncompact_checked": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "scope": (
            "Complete exact audit of the two a=2/add5 compact incidence "
            "branches. Other record-18 branches and ideal vertices are not "
            "disposed of by this artifact."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    record = run()
    rendered = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if args.write is not None:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
