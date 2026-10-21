#!/usr/bin/env python3
"""Independent compact-label replay for the first Record-1 a=1 target.

No production completion, Coxeter-label, or open-shadow module is imported.
The classical Lanner tables and finite Coxeter classifier come from the
separately written hostile Record-17 replay.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json

from scripts.independent_h9_14_record17_all_coxeter_audit import (
    LabelState,
    _finite,
    _lanner_codes,
    _options,
    _pair,
    _state_digest,
)


N = 14
DIMENSION = 9
EDGES = ((0, 13), (1, 13), (2, 13))
FOURS = ((0, 3, 4, 5),)
ADDITIONS = (
    (1, 3, 4, 7, 8),
    (2, 3, 4, 7, 8),
    (3, 4, 5, 7, 8),
    (1, 3, 4, 8, 9),
    (0, 4, 5, 9, 10),
    (1, 4, 5, 9, 10),
    (0, 5, 6, 9, 10),
    (1, 2, 7, 8, 11),
    (1, 2, 8, 9, 11),
    (1, 4, 8, 9, 11),
    (1, 4, 9, 10, 11),
    (1, 6, 9, 10, 11),
    (0, 3, 5, 7, 12),
    (0, 2, 6, 7, 12),
    (2, 3, 7, 8, 12),
    (0, 5, 6, 10, 12),
    (2, 6, 7, 11, 12),
    (2, 6, 10, 11, 12),
    (6, 9, 10, 11, 12),
)


def replay() -> dict[str, object]:
    minimal_nonfaces = tuple(map(frozenset, (*EDGES, *FOURS, *ADDITIONS)))
    actual_vertices = tuple(
        face for face in combinations(range(N), DIMENSION)
        if not any(nonface <= set(face) for nonface in minimal_nonfaces)
    )
    if len(actual_vertices) != 158:
        raise AssertionError("independent rank-nine vertex count changed")
    supports = (*FOURS, *ADDITIONS)
    tables = {4: _lanner_codes(4), 5: _lanner_codes(5)}
    remaining = set(range(len(supports)))
    cluster: tuple[int, ...] = ()
    states: tuple[LabelState, ...] = ((),)
    stages = []
    while states and remaining:
        if not cluster:
            support_index = min(remaining)
        else:
            union = set(cluster)
            support_index = max(
                remaining,
                key=lambda index: (
                    len(set(supports[index]) & union),
                    -len(set(supports[index]) - union),
                    -index,
                ),
            )
        support = supports[support_index]
        cluster = tuple(sorted(set(cluster) | set(support)))
        cluster_set = set(cluster)
        restrictions = tuple(
            sorted(
                {
                    tuple(vertex for vertex in face if vertex in cluster_set)
                    for face in actual_vertices
                },
                key=lambda value: (len(value), value),
            )
        )
        options = _options(support, tables)
        compatible = 0
        passing: set[LabelState] = set()
        for state in states:
            old = dict(state)
            for option in options:
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    _pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(_finite(restriction, labels) for restriction in restrictions):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "complete_lanner_placement_count": len(options),
                "compatible_extension_count": compatible,
                "restriction_size_histogram": {
                    str(size): count
                    for size, count in sorted(Counter(map(len, restrictions)).items())
                },
                "passing_state_count": len(states),
                "passing_state_sha256": _state_digest(states),
            }
        )
        remaining.remove(support_index)
    return {
        "kind": "h9_14_record1_a1_target0_independent_compact_label_audit",
        "rank_nine_actual_vertex_count": len(actual_vertices),
        "lanner_table_counts": {str(order): len(table) for order, table in tables.items()},
        "stages": stages,
        "terminal_state_count": len(states),
        "unrestricted_compact_label_no_go": bool(stages) and not states,
        "scope": "This excludes one literal compact target, not its source orbit or shadow.",
    }


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
