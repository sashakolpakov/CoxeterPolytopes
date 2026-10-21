#!/usr/bin/env python3
"""Exact exploratory A5 completion for an open H9/14 positive-A4 source.

This generalizes the already frozen Record-14 fixed-source Z3 census to the
seven still-open support shadows.  A zero-target complete row is a compact
incidence no-go for that source; a surviving target is passed to the complete
Lanner/actual-vertex compact-label audit.  This script is exploratory until a
separate independent replay and frozen certifier are supplied.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
from typing import Sequence

from coxeter_polytopes.combinatorics.h9_14_record46 import _base_link_g3_g4
from coxeter_polytopes.combinatorics.h9_14_record92 import independent_masks
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6
from scripts.explore_h9_14_open_positive_a4 import (
    RECORD_GRAPHS,
    audit_record,
    mask,
    vertices,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    enumerate_models,
    finite_code,
    hits_all,
    macaulay_raise,
    pair,
    placement_options,
    rank_nine_vertices,
    restriction_finite,
)


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
N = 14
FULL = (1 << N) - 1


def source_representative(record: int, a: int, orbit_index: int) -> tuple[int, ...]:
    scout = audit_record(record, a)
    branch = next(row for row in scout["branches"] if row["a"] == a)
    representatives = branch["representatives"]
    if not 0 <= orbit_index < len(representatives):
        raise ValueError(
            f"orbit index must lie in 0..{len(representatives) - 1} "
            f"for record {record}, a={a}"
        )
    return tuple(mask(face) for face in representatives[orbit_index])


def build_problem(record: int, a: int, orbit_index: int) -> dict[str, object]:
    graph6 = RECORD_GRAPHS[record]
    _, edge_pairs = decode_graph6(graph6)
    edge_rows = tuple(mask(edge) for edge in edge_pairs)
    fours = source_representative(record, a, orbit_index)
    all_fives = independent_masks(graph6, 5)
    shadow = tuple(
        five for five in all_fives if any(five & four == four for four in fours)
    )
    support = json.loads(SUPPORT.read_text(encoding="utf-8"))
    source = next(row for row in support["records"] if row["graph6"] == graph6)
    schedule = next(
        row for row in source["dehn_sommerville_schedule"]
        if int(row["missing_four_count"]) == a
    )
    total_fives = int(schedule["missing_five_count_offset"])
    addition_count = total_fives - len(shadow)
    scheduled_facets = int(schedule["f_vector"][-1])
    base_g3, base_g4 = _base_link_g3_g4(graph6)
    lower = []
    upper = []
    for vertex in range(N):
        q = sum(bool(four & (1 << vertex)) for four in fours)
        forced = sum(bool(five & (1 << vertex)) for five in shadow)
        previous = base_g3[vertex] - q
        before = base_g4[vertex] + 6 * q - forced
        lower.append(max(0, before - macaulay_raise(previous, 3)))
        upper.append(before - 1)
    if addition_count < 0 or any(
        left > right for left, right in zip(lower, upper, strict=True)
    ):
        raise AssertionError("archived positive source has an empty A5 interval")
    candidates = tuple(
        five for five in all_fives
        if five not in shadow
        and all(
            not five & (1 << vertex) or upper[vertex] >= 1
            for vertex in range(N)
        )
    )
    base = (*edge_rows, *fours, *shadow)
    base_transversals = tuple(
        value for value in range(1, FULL + 1) if hits_all(value, base)
    )
    small = tuple(value for value in base_transversals if value.bit_count() <= 4)
    fives = tuple(value for value in base_transversals if value.bit_count() == 5)
    return {
        "record": record,
        "graph6": graph6,
        "a": a,
        "orbit_index": orbit_index,
        "edge_rows": edge_rows,
        "fours": fours,
        "shadow": shadow,
        "addition_count": addition_count,
        "scheduled_facets": scheduled_facets,
        "lower": tuple(lower),
        "upper": tuple(upper),
        "candidates": candidates,
        "base": base,
        "small": small,
        "fives": fives,
    }


def audit_compact_labels(
    edge_rows: Sequence[int], fours: Sequence[int], additions: Sequence[int]
) -> dict[str, object]:
    supports = tuple(vertices(value) for value in (*fours, *additions))
    actual_vertices = rank_nine_vertices((*edge_rows, *fours, *additions))
    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[tuple[tuple[tuple[int, int], int], ...], ...] = ((),)
    stages = []
    while states and remaining:
        if not union:
            support_index = min(remaining)
        else:
            support_index = max(
                remaining,
                key=lambda index: (
                    len(set(supports[index]) & union),
                    -len(set(supports[index]) - union),
                    -index,
                ),
            )
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {
                    tuple(vertex for vertex in actual if vertex in union)
                    for actual in actual_vertices
                },
                key=lambda row: (len(row), row),
            )
        )
        options = placement_options(support)
        passing: set[tuple[tuple[tuple[int, int], int], ...]] = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in options:
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(
                    restriction_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "complete_lanner_placement_count": len(options),
                "restriction_size_histogram": dict(
                    sorted(Counter(map(len, restrictions)).items())
                ),
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
            }
        )
        remaining.remove(support_index)
    return {
        "rank_nine_actual_vertex_count": len(actual_vertices),
        "minimal_lanner_support_count": len(supports),
        "stages": stages,
        "terminal_sparse_label_state_count": len(states),
        "unrestricted_compact_label_no_go": not states,
        "monotonicity_scope": (
            "unassigned meeting labels are 2; increasing any to 3, 4, or 5 "
            "cannot repair a non-elliptic actual-vertex restriction"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=int, choices=tuple(RECORD_GRAPHS))
    parser.add_argument("a", type=int, choices=(1, 2))
    parser.add_argument("orbit_index", type=int)
    parser.add_argument("--describe-only", action="store_true")
    parser.add_argument("--max-models", type=int)
    parser.add_argument("--labels", action="store_true")
    parser.add_argument("--label-target", type=int)
    args = parser.parse_args()
    problem = build_problem(args.record, args.a, args.orbit_index)
    description = {
        "record_index": args.record,
        "graph6": problem["graph6"],
        "a": args.a,
        "source_orbit_index": args.orbit_index,
        "representative_A4s": [list(vertices(value)) for value in problem["fours"]],
        "forced_A5_shadow_count": len(problem["shadow"]),
        "additional_A5_count": problem["addition_count"],
        "scheduled_facet_count": problem["scheduled_facets"],
        "addition_lower": list(problem["lower"]),
        "addition_upper": list(problem["upper"]),
        "candidate_count": len(problem["candidates"]),
        "base_small_transversal_count": len(problem["small"]),
        "base_order_five_transversal_count": len(problem["fives"]),
    }
    if args.describe_only:
        print(json.dumps(description, indent=2, sort_keys=True))
        return
    census = enumerate_models(problem, args.max_models)
    target_label_audits = []
    if args.labels and census["literal_pure_targets"]:
        targets = list(enumerate(census["literal_pure_targets"]))
        if args.label_target is not None:
            targets = [row for row in targets if row[0] == args.label_target]
        for target_index, family in targets:
            audit = audit_compact_labels(
                problem["edge_rows"],
                problem["fours"],
                tuple(mask(face) for face in family),
            )
            target_label_audits.append({"target_index": target_index, **audit})
    output = {
        "kind": "h9_14_open_positive_completion_exploratory",
        "theorem_certificate": False,
        **description,
        **census,
        "target_label_audits": target_label_audits,
        "scope": (
            "Exact necessary compact incidence/ridge census for one source; "
            "exploratory until independent replay and frozen certification."
        ),
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
