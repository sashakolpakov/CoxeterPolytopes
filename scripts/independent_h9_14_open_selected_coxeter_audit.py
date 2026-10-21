#!/usr/bin/env python3
"""Independent replay of seven selected H9/14 compact-label no-gos.

This script does not import the production selected-nerve gate or the project
Coxeter/Lanner classifier.  It reads the seven frozen polytopality archives,
reconstructs their rank-nine faces from the archived minimal nonfaces, and
uses the hostile classical-table generator and finite Coxeter classifier that
were written independently for the record-17 replay.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from scripts.independent_h9_14_record17_all_coxeter_audit import (
    LabelState,
    _finite,
    _lanner_codes,
    _options,
    _pair,
    _state_digest,
)


ROOT = Path(__file__).resolve().parents[1]
RECORDS = (0, 1, 3, 4, 5, 7, 9)
FACET_COUNT = 14
DIMENSION = 9


def _summary_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("summary_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


def _next_support_index(
    supports: Sequence[tuple[int, ...]],
    remaining: Iterable[int],
    cluster_union: Iterable[int],
) -> int:
    choices = tuple(sorted(remaining))
    if not choices:
        raise ValueError("no support remains")
    union = set(cluster_union)
    if not union:
        return choices[0]
    return max(
        choices,
        key=lambda index: (
            len(set(supports[index]) & union),
            -len(set(supports[index]) - union),
            -index,
        ),
    )


def _rank_nine_vertices(
    minimal_nonfaces: Sequence[tuple[int, ...]],
) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(map(frozenset, minimal_nonfaces))
    return tuple(
        face
        for face in combinations(range(FACET_COUNT), DIMENSION)
        if not any(nonface <= set(face) for nonface in forbidden)
    )


def _load_nerve(record: int) -> tuple[
    tuple[tuple[int, int], ...],
    tuple[tuple[int, ...], ...],
    tuple[tuple[int, ...], ...],
    str,
]:
    path = ROOT / f"data/enumerations/h9_14_record{record}_polytopality.json"
    raw = path.read_bytes()
    document = json.loads(raw)
    audit = document["polytopality_audit"]
    if not audit["ordinary_rational_polytopal"]:
        raise AssertionError(f"record {record}: archived ordinary polytope vanished")
    nonfaces = tuple(tuple(map(int, face)) for face in audit["target_minimal_nonfaces"])
    if any(len(face) not in (2, 4, 5) for face in nonfaces):
        raise AssertionError(f"record {record}: unexpected minimal-nonface order")
    edges = tuple(face for face in nonfaces if len(face) == 2)
    supports = tuple(face for face in nonfaces if len(face) in (4, 5))
    vertices = _rank_nine_vertices(nonfaces)
    target_vector = audit.get("target_face_vector", document["target"]["face_vector"])
    if len(vertices) != int(target_vector[-1]):
        raise AssertionError(f"record {record}: facet count changed")
    return edges, supports, vertices, sha256(raw).hexdigest()


def _audit_record(
    record: int, tables: Mapping[int, Sequence[Sequence[int]]]
) -> dict[str, object]:
    edges, supports, vertices, input_sha = _load_nerve(record)
    remaining = set(range(len(supports)))
    states: tuple[LabelState, ...] = ((),)
    cluster_union: tuple[int, ...] = ()
    stages: list[dict[str, object]] = []
    while states and remaining:
        support_index = _next_support_index(supports, remaining, cluster_union)
        support = supports[support_index]
        cluster_union = tuple(sorted(set(cluster_union) | set(support)))
        cluster_set = set(cluster_union)
        restrictions = tuple(
            sorted(
                {
                    tuple(vertex for vertex in face if vertex in cluster_set)
                    for face in vertices
                },
                key=lambda value: (len(value), value),
            )
        )
        if any(
            set(dotted) <= set(restriction)
            for dotted in edges
            for restriction in restrictions
        ):
            raise AssertionError("an actual vertex contains a dotted pair")

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
                    for left, right in combinations(cluster_union, 2)
                }
                labels.update(merged)
                if all(_finite(restriction, labels) for restriction in restrictions):
                    passing.add(tuple(sorted(merged.items())))
        passing_states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "support_order": len(support),
                "cumulative_cluster_union": list(cluster_union),
                "actual_vertex_restriction_size_histogram": {
                    str(size): count
                    for size, count in sorted(Counter(map(len, restrictions)).items())
                },
                "complete_lanner_placement_count": len(options),
                "raw_extension_count": len(states) * len(options),
                "compatible_lanner_extension_count": compatible,
                "elliptic_actual_vertex_survivor_count": len(passing_states),
                "passing_state_sha256": _state_digest(passing_states),
            }
        )
        states = passing_states
        remaining.remove(support_index)

    return {
        "record_index": record,
        "input_polytopality_file_sha256": input_sha,
        "rank_nine_vertex_count": len(vertices),
        "minimal_nonface_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(len, (*edges, *supports))).items())
        },
        "decisive_support_indices": [stage["support_index"] for stage in stages],
        "stages": stages,
        "selected_nerve_compact_coxeter_no_go": bool(stages)
        and stages[-1]["elliptic_actual_vertex_survivor_count"] == 0,
    }


def replay() -> dict[str, object]:
    tables = {4: _lanner_codes(4), 5: _lanner_codes(5)}
    records = tuple(_audit_record(record, tables) for record in RECORDS)
    result: dict[str, object] = {
        "implementation": (
            "No production selected-nerve gate or project Coxeter/Lanner "
            "classifier is imported. Classical Lanner representatives are "
            "permuted independently, finite Coxeter diagrams are classified "
            "independently, and rank-nine faces are rebuilt from the seven "
            "polytopality archives."
        ),
        "record_indices": list(RECORDS),
        "lanner_table_counts": {str(order): len(table) for order, table in tables.items()},
        "record_gates": list(records),
        "all_seven_selected_nerves_excluded": all(
            row["selected_nerve_compact_coxeter_no_go"] for row in records
        ),
    }
    result["summary_sha256"] = _summary_digest(result)
    return result


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
