#!/usr/bin/env python3
"""Exact local-link triage for the unresolved H9/14 record 37.

This is initially a theorem-search utility.  It enumerates all locally
admissible missing-four families in the only positive schedules ``a=1,2``
and applies the forced five-shadow, strict vertex-link ``g3,g4``, and
aggregate Macaulay incidence gates without a SAT/SMT backend.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json
from math import comb
from pathlib import Path

from coxeter_polytopes.combinatorics.h9_14_record46 import (
    _a4_interval_feasible,
    _base_link_g3_g4,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import independent_masks
from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.missing_edges import macaulay_raise


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"
TARGET_INDEX = 37
TARGET_GRAPH6 = "M?????????O?O?o@?"


def _mask(vertices: tuple[int, ...]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _family_interval_feasible(family: tuple[int, ...], local) -> bool:
    for figure in local.vertex_figures:
        incidence = sum(bool(face & (1 << figure.vertex)) for face in family)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.edge_figures:
        subset = _mask(figure.edge)
        incidence = sum(face & subset == subset for face in family)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.triple_figures:
        subset = _mask(figure.triple)
        incidence = sum(face & subset == subset for face in family)
        if incidence > figure.missing_four_incidence_max:
            return False
    return True


def _forced_shadow(family: tuple[int, ...], fives: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        five
        for five in fives
        if any(five & four == four for four in family)
    )


def _strict_shadow_gate(
    family: tuple[int, ...],
    shadow: tuple[int, ...],
    base_g3: tuple[int, ...],
    base_g4: tuple[int, ...],
    total_missing_fives: int,
) -> tuple[
    bool,
    bool,
    int,
    tuple[int, ...],
    tuple[int, ...],
    tuple[int, ...],
]:
    additional = total_missing_fives - len(shadow)
    if additional < 0:
        return False, False, additional, (), (), ()
    lower = []
    upper = []
    nonpositive = []
    strict = True
    for vertex in range(14):
        a4_incidence = sum(bool(face & (1 << vertex)) for face in family)
        previous = base_g3[vertex] - a4_incidence
        forced_a5_incidence = sum(bool(face & (1 << vertex)) for face in shadow)
        before = base_g4[vertex] + 6 * a4_incidence - forced_a5_incidence
        if previous < 1 or before < 1:
            strict = False
            nonpositive.append(vertex)
        lower.append(
            max(0, before - macaulay_raise(previous, 3)) if previous >= 0 else 0
        )
        upper.append(before - 1)
    aggregate = (
        strict
        and sum(lower) <= 5 * additional <= sum(upper)
        and all(value >= 0 for value in upper)
    )
    return (
        strict,
        aggregate,
        additional,
        tuple(lower),
        tuple(upper),
        tuple(nonpositive),
    )


def enumerate_positive_branches() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    record = survivors[TARGET_INDEX]
    if record["graph6"] != TARGET_GRAPH6:
        raise ValueError("record-37 graph6 string changed")
    schedules = {
        schedule["missing_four_count"]: schedule
        for schedule in record["dehn_sommerville_schedule"]
    }
    local = analyze_local_link_graph(TARGET_GRAPH6)
    fours = tuple(
        face
        for face in independent_masks(TARGET_GRAPH6, 4)
        if _a4_interval_feasible(face, local)
    )
    fives = independent_masks(TARGET_GRAPH6, 5)
    base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)
    output: dict[str, object] = {
        "record_index": TARGET_INDEX,
        "graph6": TARGET_GRAPH6,
        "individually_interval_feasible_a4": len(fours),
        "branches": {},
    }
    for a in (1, 2):
        schedule = schedules[a]
        total_missing_fives = schedule["missing_five_count_offset"]
        locally_feasible = 0
        shadow_strict = 0
        aggregate_feasible = 0
        forced_count_overrun = 0
        shadow_sizes: Counter[int] = Counter()
        additional_counts: Counter[int] = Counter()
        aggregate_margin: Counter[int] = Counter()
        nonpositive_histogram: Counter[tuple[int, ...]] = Counter()
        survivors_out = []
        for family in combinations(fours, a):
            if not _family_interval_feasible(family, local):
                continue
            locally_feasible += 1
            shadow = _forced_shadow(family, fives)
            strict, aggregate, additional, lower, upper, nonpositive = _strict_shadow_gate(
                family,
                shadow,
                base_g3,
                base_g4,
                total_missing_fives,
            )
            shadow_sizes[len(shadow)] += 1
            additional_counts[additional] += 1
            if additional < 0:
                forced_count_overrun += 1
            else:
                nonpositive_histogram[nonpositive] += 1
            if strict:
                shadow_strict += 1
            if aggregate:
                aggregate_feasible += 1
                margin = min(5 * additional - sum(lower), sum(upper) - 5 * additional)
                aggregate_margin[margin] += 1
                if len(survivors_out) < 20:
                    survivors_out.append(
                        {
                            "a4_masks": list(family),
                            "forced_shadow_size": len(shadow),
                            "additional_minimal_a5": additional,
                            "lower_incidence_sum": sum(lower),
                            "upper_incidence_sum": sum(upper),
                        }
                    )
        output["branches"][str(a)] = {
            "total_a4_families": comb(len(fours), a),
            "locally_feasible_a4_families": locally_feasible,
            "forced_shadow_size_histogram": dict(sorted(shadow_sizes.items())),
            "additional_minimal_a5_histogram": dict(
                sorted(additional_counts.items())
            ),
            "shadow_strict_families": shadow_strict,
            "forced_shadow_count_overrun_families": forced_count_overrun,
            "nonpositive_vertex_set_histogram": {
                ",".join(map(str, key)): value
                for key, value in sorted(nonpositive_histogram.items())
            },
            "aggregate_macaulay_feasible_families": aggregate_feasible,
            "aggregate_margin_histogram": dict(sorted(aggregate_margin.items())),
            "first_aggregate_survivors": survivors_out,
        }
    return output


def main() -> None:
    output = enumerate_positive_branches()
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
