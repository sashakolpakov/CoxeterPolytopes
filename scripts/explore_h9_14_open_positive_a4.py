#!/usr/bin/env python3
"""Exact a=1,2 positive-A4 source-orbit scout for seven open H9/14 shadows.

This is an exploratory necessary-incidence reduction, not a theorem archive.
It applies the complete local vertex/edge/triple intervals, forced A5 shadow,
strict link, and aggregate Macaulay-capacity gates, then quotients surviving
A4 families by the full missing-graph automorphism group (including the
symmetric isolate block through binary column multisets).
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable, Sequence

from coxeter_polytopes.combinatorics.h9_14_record46 import (
    _a4_interval_feasible,
    _base_link_g3_g4,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    independent_masks,
)
from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.missing_edges import (
    decode_graph6,
    macaulay_raise,
)


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
RECORD_GRAPHS = {
    0: "M?????????????o??",
    1: "M?????????????w??",
    3: "M???????????_?O??",
    4: "M???????????_?o??",
    5: "M???????????_?W??",
    7: "M???????????_?w??",
    9: "M???????????_?o?_",
}
N = 14


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value >> vertex & 1)


def family_interval_feasible(family: Sequence[int], local: object) -> bool:
    for figure in local.vertex_figures:
        incidence = sum(bool(face & (1 << figure.vertex)) for face in family)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.edge_figures:
        subset = mask(figure.edge)
        incidence = sum(face & subset == subset for face in family)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.triple_figures:
        subset = mask(figure.triple)
        incidence = sum(face & subset == subset for face in family)
        if incidence > figure.missing_four_incidence_max:
            return False
    return True


def forbidden_A4_pairs(
    fours: Sequence[int], local: object
) -> frozenset[tuple[int, int]]:
    """Compile exactly the capacity-one constraints relevant when a=2.

    Every member of ``fours`` already passes every interval on its own.  Thus
    adding a second member can fail only when both hit a constraint whose
    upper endpoint is one; constraints of capacity at least two cannot bind.
    """

    constraints: list[tuple[int, str]] = []
    for figure in local.vertex_figures:
        interval = figure.missing_four_incidence_interval
        if interval is not None:
            assert interval[0] == 0
        if interval is not None and interval[1] == 1:
            constraints.append((1 << figure.vertex, "meet"))
    for figure in local.edge_figures:
        interval = figure.missing_four_incidence_interval
        if interval is not None:
            assert interval[0] == 0
        if interval is not None and interval[1] == 1:
            constraints.append((mask(figure.edge), "contain"))
    for figure in local.triple_figures:
        if figure.missing_four_incidence_max == 1:
            constraints.append((mask(figure.triple), "contain"))

    answer: set[tuple[int, int]] = set()
    for subset, relation in constraints:
        hits = (
            [face for face in fours if face & subset]
            if relation == "meet"
            else [face for face in fours if face & subset == subset]
        )
        answer.update(combinations(hits, 2))
    return frozenset(answer)


def shadow_gate(
    family: Sequence[int],
    shadow_bits_by_four: dict[int, int],
    five_vertex_bits: Sequence[int],
    base_g3: Sequence[int],
    base_g4: Sequence[int],
    total_fives: int,
) -> tuple[bool, bool, int, int]:
    shadow_bits = 0
    for four in family:
        shadow_bits |= shadow_bits_by_four[four]
    forced_count = shadow_bits.bit_count()
    additional = total_fives - forced_count
    if additional < 0:
        return False, False, forced_count, additional
    lower = []
    upper = []
    strict = True
    for vertex in range(N):
        a4_incidence = sum(bool(face & (1 << vertex)) for face in family)
        previous = base_g3[vertex] - a4_incidence
        forced_incidence = (shadow_bits & five_vertex_bits[vertex]).bit_count()
        before = base_g4[vertex] + 6 * a4_incidence - forced_incidence
        strict &= previous >= 1 and before >= 1
        lower.append(max(0, before - macaulay_raise(previous, 3)) if previous >= 0 else 0)
        upper.append(before - 1)
    aggregate = (
        strict
        and all(value >= 0 for value in upper)
        and sum(lower) <= 5 * additional <= sum(upper)
    )
    return strict, aggregate, forced_count, additional


def active_automorphisms(
    edges: Sequence[tuple[int, int]], active: Sequence[int]
) -> tuple[dict[int, int], ...]:
    edge_set = {frozenset(edge) for edge in edges}
    answer = []
    for image_order in permutations(active):
        action = dict(zip(active, image_order, strict=True))
        moved = {frozenset((action[left], action[right])) for left, right in edges}
        if moved == edge_set:
            answer.append(action)
    return tuple(answer)


def family_orbit_key(
    family: Sequence[int],
    active: Sequence[int],
    isolates: Sequence[int],
    automorphisms: Sequence[dict[int, int]],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    active_mask = mask(active)
    keys = []
    for action in automorphisms:
        transformed = tuple(
            mask(action.get(vertex, vertex) for vertex in vertices(face))
            for face in family
        )
        for order in permutations(range(len(family))):
            active_rows = tuple(transformed[index] & active_mask for index in order)
            columns = tuple(
                sorted(
                    sum(
                        bool(transformed[order[row]] & (1 << isolate)) << row
                        for row in range(len(family))
                    )
                    for isolate in isolates
                )
            )
            keys.append((active_rows, columns))
    return min(keys)


def audit_record(record: int, maximum_a: int) -> dict[str, object]:
    graph6 = RECORD_GRAPHS[record]
    _, edges = decode_graph6(graph6)
    active = tuple(sorted({vertex for edge in edges for vertex in edge}))
    isolates = tuple(vertex for vertex in range(N) if vertex not in active)
    automorphisms = active_automorphisms(edges, active)
    local = analyze_local_link_graph(graph6)
    fours = tuple(
        face
        for face in independent_masks(graph6, 4)
        if _a4_interval_feasible(face, local)
    )
    fives = independent_masks(graph6, 5)
    shadow_bits_by_four = {
        four: sum(
            1 << index
            for index, five in enumerate(fives)
            if five & four == four
        )
        for four in fours
    }
    five_vertex_bits = tuple(
        sum(
            1 << index
            for index, five in enumerate(fives)
            if five & (1 << vertex)
        )
        for vertex in range(N)
    )
    base_g3, base_g4 = _base_link_g3_g4(graph6)
    assert all(family_interval_feasible((four,), local) for four in fours)
    forbidden_pairs = forbidden_A4_pairs(fours, local)
    support_document = json.loads(SUPPORT.read_text(encoding="utf-8"))
    source = next(row for row in support_document["records"] if row["graph6"] == graph6)
    schedule_fives = {
        int(row["missing_four_count"]): int(row["missing_five_count_offset"])
        for row in source["dehn_sommerville_schedule"]
    }

    branches = []
    for a in range(1, maximum_a + 1):
        raw = 0
        interval = 0
        gate_histogram: Counter[tuple[bool, bool, int, int]] = Counter()
        orbits: dict[tuple[tuple[int, ...], tuple[int, ...]], list[tuple[int, ...]]] = defaultdict(list)
        for family in combinations(fours, a):
            raw += 1
            if a == 1:
                feasible = True
            else:
                feasible = family not in forbidden_pairs
            if not feasible:
                continue
            interval += 1
            gate = shadow_gate(
                family,
                shadow_bits_by_four,
                five_vertex_bits,
                base_g3,
                base_g4,
                schedule_fives[a],
            )
            gate_histogram[gate] += 1
            if gate[0] and gate[1]:
                key = family_orbit_key(family, active, isolates, automorphisms)
                orbits[key].append(family)
        branches.append(
            {
                "a": a,
                "raw_family_count": raw,
                "interval_feasible_family_count": interval,
                "strict_aggregate_viable_labelled_count": sum(map(len, orbits.values())),
                "source_orbit_count": len(orbits),
                "source_orbit_size_histogram": {
                    str(size): count
                    for size, count in sorted(Counter(map(len, orbits.values())).items())
                },
                "representatives": [
                    [list(vertices(face)) for face in min(families)]
                    for _, families in sorted(orbits.items())
                ],
                "gate_histogram": [
                    {
                        "strict": gate[0],
                        "aggregate": gate[1],
                        "forced_shadow": gate[2],
                        "additional": gate[3],
                        "count": count,
                    }
                    for gate, count in sorted(gate_histogram.items())
                ],
            }
        )
    return {
        "record_index": record,
        "graph6": graph6,
        "missing_edges": [list(edge) for edge in edges],
        "active_vertices": list(active),
        "isolate_count": len(isolates),
        "active_automorphism_order": len(automorphisms),
        "interval_feasible_A4_count": len(fours),
        "branches": branches,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("records", nargs="*", type=int, default=[7])
    parser.add_argument("--max-a", type=int, default=2, choices=(1, 2))
    parser.add_argument(
        "--summary",
        action="store_true",
        help="omit representatives and detailed gate histograms",
    )
    args = parser.parse_args()
    if any(record not in RECORD_GRAPHS for record in args.records):
        parser.error(f"records must lie in {tuple(RECORD_GRAPHS)}")
    result = {
        "kind": "h9_14_open_positive_A4_exploratory_scout",
        "theorem_certificate": False,
        "compact_only": True,
        "records": [audit_record(record, args.max_a) for record in args.records],
    }
    if args.summary:
        for record in result["records"]:
            for branch in record["branches"]:
                branch.pop("representatives")
                branch.pop("gate_histogram")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
