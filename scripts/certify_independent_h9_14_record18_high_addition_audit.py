#!/usr/bin/env python3
"""Build, run, and freeze the no-SAT record-18 high-addition replay.

This covers all fourteen ``a=2`` representatives (the twelve add4 cases and
the two add5 cases) as well as all three ``a=1`` representatives.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/independent_h9_14_record18_a2_add5_audit.cpp"
DEFAULT_OUTPUT = (
    ROOT
    / "data/enumerations/h9_14_record18_high_addition_independent_audit.json"
)
A2_ADD4_CASES = (
    (53, 197),
    (53, 86),
    (53, 198),
    (23, 101),
    (53, 58),
    (53, 90),
    (53, 202),
    (23, 57),
    (23, 105),
    (15, 53),
    (23, 43),
    (23, 45),
)
A2_ADD5_CASES = ((53, 54), (23, 27))
A1_ADD6_CASES = ((23,), (15,), (53,))
CASES = (*A2_ADD4_CASES, *A2_ADD5_CASES, *A1_ADD6_CASES)

EXPECTED = {
    (53, 197): {
        "search_node_count": 393014,
        "distinct_state_counts": [1, 31, 1428, 20503, 219241],
        "target_count": 0,
    },
    (53, 86): {
        "search_node_count": 281903,
        "distinct_state_counts": [1, 16, 498, 8126, 136913],
        "target_count": 0,
    },
    (53, 198): {
        "search_node_count": 1118698,
        "distinct_state_counts": [1, 31, 1428, 32871, 647128],
        "target_count": 1,
    },
    (23, 101): {
        "search_node_count": 134461,
        "distinct_state_counts": [1, 16, 498, 5984, 65967],
        "target_count": 0,
    },
    (53, 58): {
        "search_node_count": 54909,
        "distinct_state_counts": [1, 8, 165, 1836, 22953],
        "target_count": 0,
    },
    (53, 90): {
        "search_node_count": 243493,
        "distinct_state_counts": [1, 16, 498, 8101, 122216],
        "target_count": 0,
    },
    (53, 202): {
        "search_node_count": 949885,
        "distinct_state_counts": [1, 31, 1428, 32546, 572810],
        "target_count": 2,
    },
    (23, 57): {
        "search_node_count": 31744,
        "distinct_state_counts": [1, 8, 165, 1340, 12840],
        "target_count": 0,
    },
    (23, 105): {
        "search_node_count": 112491,
        "distinct_state_counts": [1, 16, 498, 5785, 54858],
        "target_count": 0,
    },
    (15, 53): {
        "search_node_count": 40013,
        "distinct_state_counts": [1, 8, 165, 1595, 16823],
        "target_count": 0,
    },
    (23, 43): {
        "search_node_count": 13439,
        "distinct_state_counts": [1, 8, 128, 898, 4837],
        "target_count": 0,
    },
    (23, 45): {
        "search_node_count": 40352,
        "distinct_state_counts": [1, 8, 165, 1595, 16906],
        "target_count": 0,
    },
    (53, 54): {
        "search_node_count": 413064,
        "distinct_state_counts": [1, 6, 81, 997, 17947, 215909],
        "target_count": 0,
    },
    (23, 27): {
        "search_node_count": 87819,
        "distinct_state_counts": [1, 4, 41, 421, 6641, 40659],
        "target_count": 0,
    },
    (23,): {
        "search_node_count": 3232195,
        "distinct_state_counts": [1, 4, 41, 430, 8342, 214722, 2041840],
        "target_count": 0,
    },
    (15,): {
        "search_node_count": 361673,
        "distinct_state_counts": [1, 2, 16, 121, 1930, 41280, 190029],
        "target_count": 0,
    },
    (53,): {
        "search_node_count": 17415998,
        "distinct_state_counts": [1, 6, 109, 1409, 33929, 900755, 11680696],
        "target_count": 8,
    },
}


def run() -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="record18-high-addition-") as directory:
        binary = Path(directory) / "audit"
        subprocess.run(
            ["c++", "-std=c++20", "-O3", str(SOURCE), "-o", str(binary)],
            cwd=ROOT,
            check=True,
        )
        records = []
        for arguments in CASES:
            completed = subprocess.run(
                [str(binary), *map(str, arguments)],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            record = json.loads(completed.stdout)
            expected = EXPECTED[arguments]
            for field in ("search_node_count", "distinct_state_counts"):
                if record[field] != expected[field]:
                    raise AssertionError(f"transcript changed at {arguments}/{field}")
            target_count = expected["target_count"]
            if (
                record["scheduled_ridge_terminal_count"],
                record["literal_pure_target_count"],
                record["canonical_target_count"],
            ) != (target_count, target_count, target_count):
                raise AssertionError(f"target census changed at {arguments}")
            records.append(record)

    by_key = {
        tuple(record["representative_A4_masks"]): record for record in records
    }
    for key in ((53, 54), (23, 27)):
        record = by_key[key]
        if (
            record["candidate_count"],
            record["base_small_transversal_count"],
            record["base_order_five_transversal_count"],
            record["nonempty_ridge_constraint_count"],
            record["compact_incidence_unsat"],
        ) != (777, 56, 233, 630, True):
            raise AssertionError(f"a2/add5 conclusion changed at {key}")
    for key in A2_ADD4_CASES:
        record = by_key[key]
        if (
            record["candidate_count"] != 776
            or record["forced_A5_shadow_count"] != 16
            or record["additional_A5_count"] != 4
            or record["addition_lower"] != [0] * 14
        ):
            raise AssertionError(f"a2/add4 base schedule changed at {key}")
    nonempty_add4 = {
        key: (
            by_key[key]["canonical_target_count"],
            by_key[key]["canonical_target_orbit_sizes"],
            by_key[key]["literal_target_orbit_weight_sum"],
        )
        for key in A2_ADD4_CASES
        if by_key[key]["canonical_target_count"]
    }
    if nonempty_add4 != {
        (53, 198): (1, [6], 6),
        (53, 202): (2, [1, 1], 2),
    }:
        raise AssertionError(f"a2/add4 target census changed: {nonempty_add4}")
    for key in ((23,), (15,), (53,)):
        record = by_key[key]
        if (
            record["candidate_count"],
            record["base_small_transversal_count"],
            record["base_order_five_transversal_count"],
            record["nonempty_ridge_constraint_count"],
        ) != (784, 64, 254, 665):
            raise AssertionError(f"a1/add6 base census changed at {key}")
    target = by_key[(53,)]
    if (
        sorted(target["canonical_target_orbit_sizes"])
        != [20, 40, 60, 120, 180, 180, 240, 240]
        or target["literal_target_orbit_weight_sum"] != 1080
        or target["compact_incidence_unsat"] is not False
    ):
        raise AssertionError("a1/add6 target orbit census changed")
    return {
        "kind": "independent_h9_14_record18_high_addition_combined_audit",
        "backend": (
            "standalone collision-free 128-bit incidence-matrix quotient; "
            "no SAT/SMT"
        ),
        "cases": records,
        "a2_add4_representative_count": len(A2_ADD4_CASES),
        "a2_add4_zero_target_representative_count": 10,
        "a2_add4_target_representatives": [[53, 198], [53, 202]],
        "a2_add4_target_literal_weight": 8,
        "a2_add4_target_subgroup_class_count": 3,
        "a2_add5_all_compact_incidence_unsat": True,
        "a1_add6_unsat_representatives": [[15], [23]],
        "a1_add6_target_representative": [53],
        "a1_add6_target_literal_weight": 1080,
        "a1_add6_target_orbit_count": 8,
        "compact_only": True,
        "finite_volume_noncompact_checked": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "scope": (
            "Complete exact incidence/ridge audit of record-18 a=1/add6 and "
            "all a=2 cases. Topology and Coxeter labels are separate; ideal "
            "vertices are not tested."
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
