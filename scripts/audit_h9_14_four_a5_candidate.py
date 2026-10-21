#!/usr/bin/env python3
"""Audit an explicit H9/14 ``a=0`` blocker candidate.

This checker is independent of the compressed discovery census.  It scans all
``2^14`` subsets to reconstruct the complete blocker, then passes the resulting
five-transversals to the generic global incidence, Euler, and mod-two link
homology auditor.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
    generated_face_table,
)
from coxeter_polytopes.combinatorics.gale_chirotope import (
    audit_uniform_rank_four_gale,
    solve_uniform_rank_four_gale,
)
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6
from coxeter_polytopes.combinatorics.nerve_support import verify_support_manifest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
VERTEX_COUNT = 14
FULL_MASK = (1 << VERTEX_COUNT) - 1


def _vertices(mask: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if mask & (1 << vertex))


def _minimal_transversals(hyperedges: tuple[int, ...]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, FULL_MASK + 1):
        if not all(candidate & edge for edge in hyperedges):
            continue
        remaining = candidate
        minimal = True
        while remaining:
            bit = remaining & -remaining
            if all((candidate ^ bit) & edge for edge in hyperedges):
                minimal = False
                break
            remaining ^= bit
        if minimal:
            answers.append(candidate)
    return tuple(answers)


def _packed_face_table_sha256(facets: Iterable[int]) -> str:
    table = generated_face_table(facets)
    packed = bytes(
        sum(bool(table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(table), 8)
    )
    return sha256(packed).hexdigest()


def audit(index: int, missing_five_masks: Sequence[int]) -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    verify_support_manifest(source)
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    record = survivors[index]
    graph6 = record["graph6"]
    order, edges = decode_graph6(graph6)
    if order != VERTEX_COUNT:
        raise ValueError("candidate graph does not have fourteen vertices")
    fives = tuple(sorted(set(int(mask) for mask in missing_five_masks)))
    a0_schedule = next(
        schedule
        for schedule in record["dehn_sommerville_schedule"]
        if schedule["missing_four_count"] == 0
    )
    expected_a5_count = a0_schedule["missing_five_count_offset"]
    expected_facet_count = a0_schedule["f_vector"][8]
    if len(fives) != expected_a5_count or any(
        mask <= 0 or mask & ~FULL_MASK or mask.bit_count() != 5 for mask in fives
    ):
        raise ValueError(
            f"candidate must contain {expected_a5_count} distinct five-set masks"
        )
    edge_masks = tuple(sum(1 << vertex for vertex in edge) for edge in edges)
    if any(any(face & edge == edge for edge in edge_masks) for face in fives):
        raise ValueError("candidate A5 contains a missing edge")

    blocker = _minimal_transversals(edge_masks + fives)
    histogram = Counter(mask.bit_count() for mask in blocker)
    five_blocker = tuple(mask for mask in blocker if mask.bit_count() == 5)
    global_audit = audit_global_facet_selection(
        graph6,
        record["ridge_support_core_masks"],
        record["dehn_sommerville_schedule"],
        five_blocker,
    )
    return {
        "record_index": index,
        "graph6": graph6,
        "missing_edges": [list(edge) for edge in edges],
        "missing_five_masks": list(fives),
        "missing_five_vertices": [list(_vertices(mask)) for mask in fives],
        "minimal_transversal_size_histogram": {
            str(size): count for size, count in sorted(histogram.items())
        },
        "expected_missing_five_count": expected_a5_count,
        "expected_facet_count": expected_facet_count,
        "blocker_is_pure_order_five": histogram == {5: expected_facet_count},
        "facet_complement_masks": list(five_blocker),
        "facet_complement_sha256": sha256(
            b"".join(mask.to_bytes(2, "little") for mask in five_blocker)
        ).hexdigest(),
        "packed_face_table_sha256": _packed_face_table_sha256(five_blocker),
        "global_audit": global_audit.to_record(),
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("index", type=int)
    parser.add_argument("missing_five_masks", nargs="+", type=int)
    parser.add_argument(
        "--gale-timeout-ms",
        type=int,
        help="also search for and independently audit a rank-four Gale sign word",
    )
    args = parser.parse_args(argv)
    result = audit(args.index, args.missing_five_masks)
    if args.gale_timeout_ms is not None:
        solve = solve_uniform_rank_four_gale(
            result["facet_complement_masks"],
            timeout_ms=args.gale_timeout_ms,
        )
        gale_record: dict[str, object] = {
            "status": solve.status,
            "backend": solve.backend,
            "backend_version": solve.backend_version,
            "elapsed_seconds": solve.elapsed_seconds,
            "assertion_count": solve.assertion_count,
            "reason_unknown": solve.reason_unknown,
            "signs": list(solve.signs),
        }
        if solve.signs:
            independent = audit_uniform_rank_four_gale(
                result["facet_complement_masks"], solve.signs
            )
            gale_record["independent_audit"] = {
                "passes": independent.passes,
                "sign_count": independent.sign_count,
                "first_grassmann_pluecker_failure": (
                    list(independent.first_grassmann_pluecker_failure)
                    if independent.first_grassmann_pluecker_failure is not None
                    else None
                ),
                "positive_circuit_count": len(
                    independent.positive_circuit_masks
                ),
                "positive_circuits_match_facets": (
                    independent.positive_circuits_match_facets
                ),
            }
        result["gale_search"] = gale_record
    print(
        json.dumps(
            result,
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
