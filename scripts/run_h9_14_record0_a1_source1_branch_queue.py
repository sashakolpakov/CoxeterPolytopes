#!/usr/bin/env python3
"""Run the lossless source-1 labelled-A4 branch queue fail-fast."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
WORKER = ROOT / "scripts/explore_h9_14_record0_a1_joint_label_cegar.py"
DEFAULT_SEED = (
    ROOT
    / "data/enumerations/h9_14_record0_a1_source1_joint_label_cegar_checkpoint.json"
)


def branch_path(index: int) -> Path:
    return (
        ROOT
        / "data/enumerations"
        / f"h9_14_record0_a1_source1_branch{index:02d}_joint_checkpoint.json"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--end", type=int, default=18)
    parser.add_argument("--seed", type=Path, default=DEFAULT_SEED)
    parser.add_argument("--seconds-per-solve", type=float, default=600.0)
    parser.add_argument("--max-rounds", type=int, default=5)
    parser.add_argument("--rerun", action="store_true")
    args = parser.parse_args()
    if not 0 <= args.start <= args.end <= 18:
        raise ValueError("source-1 branch range must lie in 0..18")

    environment = os.environ.copy()
    environment["PYTHONPATH"] = "src:."
    ledger = []
    for index in range(args.start, args.end + 1):
        output = branch_path(index)
        if output.exists() and not args.rerun:
            document = json.loads(output.read_text(encoding="utf-8"))
            if document.get("a4_branch_compact_coxeter_no_go") is True:
                print(f"branch {index}: already INFEASIBLE", flush=True)
                ledger.append(
                    {
                        "branch": index,
                        "status": document["terminal_solver_status"],
                        "wall_time_seconds": document.get(
                            "terminal_solver_wall_time_seconds"
                        ),
                        "placement_orbit_size": document["a4_branch"][
                            "placement_orbit_size"
                        ],
                    }
                )
                continue
        command = [
            sys.executable,
            str(WORKER),
            "1",
            "--a4-orbit-index",
            str(index),
            "--seed-cuts",
            str(args.seed),
            "--max-rounds",
            str(args.max_rounds),
            "--seconds-per-solve",
            str(args.seconds_per_solve),
            "--write",
            str(output),
            "--summary-only",
        ]
        completed = subprocess.run(
            command,
            cwd=ROOT,
            env=environment,
            check=False,
        )
        if completed.returncode:
            raise SystemExit(completed.returncode)
        document = json.loads(output.read_text(encoding="utf-8"))
        ledger.append(
            {
                "branch": index,
                "status": document["terminal_solver_status"],
                "wall_time_seconds": document.get(
                    "terminal_solver_wall_time_seconds"
                ),
                "placement_orbit_size": document["a4_branch"][
                    "placement_orbit_size"
                ],
            }
        )
        print(
            f"branch {index}: {document['terminal_solver_status']} "
            f"wall={document.get('terminal_solver_wall_time_seconds')}",
            flush=True,
        )
        if document.get("a4_branch_compact_coxeter_no_go") is not True:
            print("queue stopped at first non-closed branch", flush=True)
            break
    print(json.dumps(ledger, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
