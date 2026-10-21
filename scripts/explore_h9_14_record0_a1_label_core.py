#!/usr/bin/env python3
"""Extract transferable compact-label no-good cores for Record 0, ``a=1``.

The ordinary label audit can kill one complete incidence target after seeing
only a handful of its minimal Lanner supports.  This module shrinks that
failure to:

* a small family of missing four/five supports that must be Lanner; and
* a small family of selected rank-nine vertices that must be elliptic.

That conjunction is itself impossible, independently of every other face of
the target.  It can therefore be fed back to the incidence solver as a much
stronger no-good than blocking the complete 22--24 member A5 family.

This is an exploratory accelerator.  Every reported core is replayed by the
same complete finite/Lanner classifiers used by the unrestricted label audit;
an independent frozen replay is still required before a whole-source theorem.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from typing import Iterable, Sequence

from scripts.explore_h9_14_open_positive_completion import audit_compact_labels
from scripts.explore_h9_14_record0_a1_dual_facets import MISSING_EDGES, solve_one
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    placement_options,
    restriction_finite,
    vertices,
)


N = 14
FULL = (1 << N) - 1


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def _support_order(supports: Sequence[int]) -> tuple[int, ...]:
    remaining = set(range(len(supports)))
    union: set[int] = set()
    answer = []
    while remaining:
        if not union:
            index = min(remaining)
        else:
            index = max(
                remaining,
                key=lambda candidate: (
                    len(set(vertices(supports[candidate])) & union),
                    -len(set(vertices(supports[candidate])) - union),
                    -candidate,
                ),
            )
        answer.append(index)
        union.update(vertices(supports[index]))
        remaining.remove(index)
    return tuple(answer)


def label_state_count(
    supports: Sequence[int], actual_vertices: Sequence[Sequence[int]]
) -> tuple[int, tuple[dict[str, int], ...]]:
    """Count sparse label states satisfying these local requirements.

    Unassigned meeting pairs are set to label two.  If such a restriction is
    non-elliptic, increasing an unassigned label cannot repair it, so pruning
    is complete for the alphabet ``{2,3,4,5}``.
    """

    order = _support_order(supports)
    states: tuple[tuple[tuple[tuple[int, int], int], ...], ...] = ((),)
    union: set[int] = set()
    stages = []
    for source_index in order:
        support = vertices(supports[source_index])
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
        passing: set[tuple[tuple[tuple[int, int], int], ...]] = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in placement_options(support):
                if not all(
                    old.get(edge, label) == label for edge, label in option
                ):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    (left, right): 2
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
                "support_index": source_index,
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
                "restriction_count": len(restrictions),
            }
        )
        if not states:
            break
    return len(states), tuple(stages)


def _deduplicate_dominated_vertices(
    actual_vertices: Sequence[Sequence[int]], union: set[int]
) -> tuple[tuple[int, ...], ...]:
    by_restriction: dict[frozenset[int], tuple[int, ...]] = {}
    for actual in actual_vertices:
        actual_tuple = tuple(actual)
        restriction = frozenset(vertex for vertex in actual if vertex in union)
        by_restriction.setdefault(restriction, actual_tuple)
    restrictions = tuple(by_restriction)
    maximal = tuple(
        restriction
        for restriction in restrictions
        if not any(restriction < other for other in restrictions)
    )
    return tuple(by_restriction[restriction] for restriction in maximal)


def extract_label_core(
    four: int,
    additions: Sequence[int],
    selected_facet_complements: Sequence[int],
) -> dict[str, object]:
    supports = (four, *additions)
    audit = audit_compact_labels(MISSING_EDGES, (four,), additions)
    if not audit["unrestricted_compact_label_no_go"]:
        raise ValueError("target has a compact-label survivor, not a no-go core")
    death_indices = tuple(
        int(stage["support_index"]) for stage in audit["stages"]
    )
    core_supports = tuple(supports[index] for index in death_indices)
    actual_vertices = tuple(
        tuple(vertex for vertex in range(N) if not facet & (1 << vertex))
        for facet in selected_facet_complements
    )
    union = set().union(*(set(vertices(support)) for support in core_supports))
    core_vertices = list(
        _deduplicate_dominated_vertices(actual_vertices, union)
    )
    initial_vertex_count = len(core_vertices)
    state_count, _ = label_state_count(core_supports, core_vertices)
    if state_count:
        raise AssertionError("deduplicated vertex family lost the label no-go")

    # Monotonic deletion: once removal makes the system satisfiable, deleting
    # still more constraints cannot make that particular vertex redundant.
    index = 0
    while index < len(core_vertices):
        trial = (*core_vertices[:index], *core_vertices[index + 1 :])
        trial_count, _ = label_state_count(core_supports, trial)
        if not trial_count:
            core_vertices.pop(index)
        else:
            index += 1

    # Then remove any Lanner support not needed by the local contradiction,
    # retaining the fixed A4 as the geometric source marker.
    core_support_list = list(core_supports)
    index = 0
    while index < len(core_support_list):
        if core_support_list[index] == four:
            index += 1
            continue
        trial = (*core_support_list[:index], *core_support_list[index + 1 :])
        trial_count, _ = label_state_count(trial, core_vertices)
        if not trial_count:
            core_support_list.pop(index)
        else:
            index += 1
    core_supports = tuple(core_support_list)
    final_count, final_stages = label_state_count(core_supports, core_vertices)
    if final_count:
        raise AssertionError("minimized local label core is satisfiable")

    required_facets = tuple(
        sorted(FULL ^ mask(actual) for actual in core_vertices)
    )
    union_mask = mask(union)
    face_restrictions = tuple(
        sorted(mask(actual) & union_mask for actual in core_vertices)
    )
    selected_set = set(selected_facet_complements)
    if any(facet not in selected_set for facet in required_facets):
        raise AssertionError("core vertex does not belong to the target")
    return {
        "initial_death_support_count": len(death_indices),
        "initial_maximal_vertex_restriction_count": initial_vertex_count,
        "core_supports": [list(vertices(value)) for value in core_supports],
        "core_additional_A5s": [
            list(vertices(value)) for value in core_supports if value != four
        ],
        "core_actual_vertices": [list(value) for value in core_vertices],
        "core_face_restrictions": [
            list(vertices(value)) for value in face_restrictions
        ],
        "core_facet_complements": [
            list(vertices(value)) for value in required_facets
        ],
        "core_support_count": len(core_supports),
        "core_actual_vertex_count": len(core_vertices),
        "terminal_sparse_label_state_count": final_count,
        "replay_stages": list(final_stages),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_orbit", type=int, choices=range(4))
    parser.add_argument("--seconds", type=float, default=180.0)
    args = parser.parse_args()
    result = solve_one(args.source_orbit, args.seconds, 1)
    target = result["target"]
    if target is None or not target["pure_expected_blocker"]:
        print(json.dumps({"solve": result, "core": None}, indent=2))
        return
    problem_four = mask(result["representative_A4s"][0])
    blockers = tuple(mask(row) for row in target["minimal_nonfaces"])
    additions = tuple(value for value in blockers if value.bit_count() == 5)
    selected = tuple(mask(row) for row in target["selected_facet_complements"])
    core = extract_label_core(problem_four, additions, selected)
    print(
        json.dumps(
            {
                "source_orbit_index": args.source_orbit,
                "solver_status": result["solver_status"],
                "core": core,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
