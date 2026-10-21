#!/usr/bin/env python3
"""Hostile quotient/label replay for the Record-14 ``a=0`` theorem.

This replay imports neither production completion code, the production
symmetry canonicalizer, nor the project Coxeter/Lanner classifier.  It reads
the frozen literal target families, quotients them by colored bipartite graph
isomorphism, regenerates the classical order-five Lanner table, and applies an
independent finite-Coxeter classifier to every actual rank-nine vertex.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Mapping, Sequence

import networkx as nx

from scripts.independent_h9_14_record17_all_coxeter_audit import (
    LabelState,
    _finite,
    _lanner_codes,
    _options,
    _pair,
    _state_digest,
)


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record14_a0_whole_shadow_no_go.json"
OUTPUT = (
    ROOT
    / "data/enumerations/h9_14_record14_a0_downstream_independent_audit.json"
)
EDGES = ((0, 12), (1, 12), (0, 13), (1, 13), (2, 13))
EXPECTED_GROUPS = (
    (0, 1, 7, 16, 19, 23, 26, 30, 32, 34, 38, 44, 45, 51),
    (2, 5, 6, 11, 12, 15, 18, 21, 25, 29, 33, 35, 36, 40, 41, 42, 47, 49),
    (3, 4, 8, 9, 10, 13, 14, 17, 20, 22, 24, 27, 28, 31, 37, 39, 43, 46, 48, 50),
)
EXPECTED_SUMMARY_SHA256 = "e871f807999f07e3e89a213669d34fab8ff717abec74c16fe9a96a61337966c3"


def _summary_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("summary_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _target_graph(family: Sequence[Sequence[int]]) -> nx.Graph:
    graph = nx.Graph()
    for row in range(len(family)):
        graph.add_node(("row", row), color="row")
    for vertex in range(14):
        if vertex in (0, 1):
            color = "active_twin"
        elif vertex in range(3, 12):
            color = "isolate"
        else:
            color = f"fixed_{vertex}"
        graph.add_node(("column", vertex), color=color)
    for row, face in enumerate(family):
        for vertex in face:
            graph.add_edge(("row", row), ("column", vertex))
    return graph


def _isomorphic(first: nx.Graph, second: nx.Graph) -> bool:
    return nx.is_isomorphic(
        first,
        second,
        node_match=nx.algorithms.isomorphism.categorical_node_match("color", None),
    )


def _independent_groups(
    targets: Sequence[Sequence[Sequence[int]]],
) -> tuple[tuple[int, ...], ...]:
    graphs = tuple(_target_graph(family) for family in targets)
    representatives: list[int] = []
    groups: list[list[int]] = []
    for index, graph in enumerate(graphs):
        match = next(
            (
                position
                for position, representative in enumerate(representatives)
                if _isomorphic(graph, graphs[representative])
            ),
            None,
        )
        if match is None:
            representatives.append(index)
            groups.append([index])
        else:
            groups[match].append(index)
    return tuple(sorted(tuple(group) for group in groups))


def _actual_vertices(
    supports: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(map(frozenset, (*EDGES, *supports)))
    return tuple(
        chosen
        for chosen in combinations(range(14), 9)
        if not any(row <= set(chosen) for row in forbidden)
    )


def _label_audit(
    supports: Sequence[tuple[int, ...]], tables: Mapping[int, Sequence[Sequence[int]]]
) -> dict[str, object]:
    vertices = _actual_vertices(supports)
    remaining = set(range(len(supports)))
    states: tuple[LabelState, ...] = ((),)
    union: set[int] = set()
    stages = []
    while states and remaining:
        support_index = (
            min(remaining)
            if not union
            else max(
                remaining,
                key=lambda index: (
                    len(set(supports[index]) & union),
                    -len(set(supports[index]) - union),
                    -index,
                ),
            )
        )
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {
                    tuple(vertex for vertex in actual if vertex in union)
                    for actual in vertices
                },
                key=lambda value: (len(value), value),
            )
        )
        if any(
            set(edge) <= set(restriction)
            for edge in EDGES
            for restriction in restrictions
        ):
            raise AssertionError("hostile actual vertex contains a missing edge")
        options = _options(support, tables)
        passing: set[LabelState] = set()
        compatible = 0
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
                "actual_vertex_restriction_size_histogram": {
                    str(size): count
                    for size, count in sorted(Counter(map(len, restrictions)).items())
                },
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
                "passing_state_sha256": _state_digest(states),
            }
        )
        remaining.remove(support_index)
    return {
        "rank_nine_actual_vertex_count": len(vertices),
        "stages": stages,
        "terminal_sparse_label_state_count": len(states),
        "unrestricted_compact_label_no_go": not states,
    }


def replay() -> dict[str, object]:
    raw = ARCHIVE.read_bytes()
    archive = json.loads(raw)
    targets = tuple(
        tuple(tuple(map(int, face)) for face in family)
        for family in archive["evidence"]["generator_lex_slice_pure_targets"]
    )
    if len(targets) != 52:
        raise AssertionError("hostile Record-14 target count changed")
    groups = _independent_groups(targets)
    if groups != EXPECTED_GROUPS:
        raise AssertionError("hostile Record-14 quotient changed")
    tables = {5: _lanner_codes(5)}
    audits = tuple(_label_audit(targets[group[0]], tables) for group in groups)
    result: dict[str, object] = {
        "kind": "independent_h9_14_record14_a0_downstream_audit",
        "input_archive_sha256": sha256(raw).hexdigest(),
        "implementation": (
            "No production completion model, symmetry canonicalizer, or "
            "project Coxeter classifier is imported. The quotient is colored "
            "bipartite graph isomorphism; the Lanner table and finite Coxeter "
            "classifier come from the hostile classical-table replay."
        ),
        "generator_lex_target_count": len(targets),
        "full_target_orbit_groups": [list(group) for group in groups],
        "full_target_orbit_count": len(groups),
        "lanner_table_count": len(tables[5]),
        "label_audits": list(audits),
        "all_target_orbits_unrestricted_compact_label_excluded": all(
            audit["unrestricted_compact_label_no_go"] for audit in audits
        ),
    }
    result["summary_sha256"] = _summary_digest(result)
    if EXPECTED_SUMMARY_SHA256 not in ("TO_BE_FROZEN", result["summary_sha256"]):
        raise AssertionError("hostile Record-14 summary digest changed")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", type=Path, nargs="?", const=OUTPUT)
    args = parser.parse_args()
    rendered = json.dumps(replay(), indent=2, sort_keys=True) + "\n"
    if args.write is not None:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
