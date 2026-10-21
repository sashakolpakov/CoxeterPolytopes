#!/usr/bin/env python3
"""Exploratory parallel full-link replay for record-12 residual orbits."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
from pathlib import Path

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from scripts.explore_h9_14_record12_pairs import (
    PAIR_INDEX,
    selected_complements,
)


ROOT = Path(__file__).resolve().parents[1]
ORBIT_PATH = Path("/tmp/record12_full_orbits.txt")
SUPPORT_PATH = ROOT / "data/enumerations/h9_14_nerve_support.json"
GRAPH6 = "M???????????o?o??"


def _pair_mask(matrix: int) -> int:
    return sum(
        1 << PAIR_INDEX[first, second + 5]
        for first in range(5)
        for second in range(5)
        if matrix & (1 << (5 * first + second))
    )


def _support() -> tuple[tuple[int, ...], tuple[dict[str, object], ...]]:
    document = json.loads(SUPPORT_PATH.read_text(encoding="utf-8"))
    record = [row for row in document["records"] if row["survives_second_layer"]][12]
    if record["graph6"] != GRAPH6:
        raise AssertionError("record-12 support row changed")
    return (
        tuple(record["ridge_support_core_masks"]),
        tuple(record["dehn_sommerville_schedule"]),
    )


SUPPORT, SCHEDULES = _support()


def _audit(task: tuple[int, int, int]) -> dict[str, object]:
    x_matrix, z_matrix, multiplicity = task
    complements = selected_complements(
        _pair_mask(x_matrix), _pair_mask(z_matrix), 992, 992, 0
    )
    audit = audit_global_facet_selection(
        GRAPH6, SUPPORT, SCHEDULES, complements
    )
    return {
        "x_matrix": x_matrix,
        "z_matrix": z_matrix,
        "source_double_coset_multiplicity": multiplicity,
        "missing_four_count": audit.missing_four_count,
        "facet_count": audit.facet_count,
        "passes_selection_model": audit.passes_selection_model,
        "passes_connected_euler_audit": audit.passes_connected_euler_audit,
        "passes_connected_homology_audit": audit.passes_connected_homology_audit,
        "first_strict_local_g_failure": audit.first_strict_local_g_failure,
        "first_local_macaulay_failure": audit.first_local_macaulay_failure,
        "first_local_euler_failure": audit.first_local_euler_failure,
        "first_mod2_homology_failure": audit.first_mod2_homology_failure,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("orbit_path", nargs="?", type=Path, default=ORBIT_PATH)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    tasks = []
    for line in args.orbit_path.read_text(encoding="utf-8").splitlines():
        marker, x_value, z_value, multiplicity = line.split()
        if marker != "ORBIT":
            raise ValueError("malformed exploratory orbit stream")
        tasks.append((int(x_value), int(z_value), int(multiplicity)))
    if args.workers == 1:
        records = list(map(_audit, tasks))
    else:
        context = multiprocessing.get_context("fork")
        with ProcessPoolExecutor(max_workers=args.workers, mp_context=context) as executor:
            records = list(executor.map(_audit, tasks, chunksize=1))
    summary = {
        "orbit_count": len(records),
        "schedule_orbit_counts": dict(
            sorted(Counter(row["missing_four_count"] for row in records).items())
        ),
        "selection_model_passing": sum(row["passes_selection_model"] for row in records),
        "connected_euler_passing": sum(
            row["passes_connected_euler_audit"] for row in records
        ),
        "connected_homology_passing": sum(
            row["passes_connected_homology_audit"] for row in records
        ),
        "failure_histogram": dict(
            Counter(
                str(
                    row["first_strict_local_g_failure"]
                    or row["first_local_macaulay_failure"]
                    or row["first_local_euler_failure"]
                    or row["first_mod2_homology_failure"]
                )
                for row in records
                if not row["passes_connected_homology_audit"]
            )
        ),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
