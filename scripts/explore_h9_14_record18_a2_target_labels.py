#!/usr/bin/env python3
"""Audit the residual record-18 ``a=2`` compact-incidence targets.

The fixed-addition search for positive orbit ``[53, 198]`` leaves six pure
order-five blockers.  This script independently dualizes those six clutters,
checks their global face/link topology, quotients them by the complete
automorphism group of the missing-edge graph, and joins every labelled Lanner
placement on every inclusion-minimal missing four- or five-face.

At each join stage, all restrictions induced by actual rank-nine vertices are
required to be elliptic.  Unassigned finite pairs are set to label two.  Thus
an uncapped empty terminal is an unrestricted finite-label no-go: increasing
an unassigned cosine adjacency entry cannot restore positive definiteness.

This is compact-only.  It says nothing about finite-volume noncompact
chambers, whose ideal vertices require parabolic rather than elliptic tests.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_record18 import (
    TARGET_GRAPH6,
    TARGET_MISSING_EDGES,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT / "data/enumerations/h9_14_record18_a2_idx3_additions.json"
)
VERTEX_COUNT = 14
DIMENSION = 9
FULL_MASK = (1 << VERTEX_COUNT) - 1
EXPECTED_F_VECTOR = (14, 86, 310, 733, 1192, 1344, 1020, 477, 106)
ACTIVE = (0, 1, 2, 3, 12, 13)
ISOLATES = tuple(range(4, 12))

LabelState = tuple[tuple[Pair, int], ...]
Support = tuple[int, ...]


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1)


def minimal_transversals(supports: Sequence[Support]) -> tuple[int, ...]:
    hyperedges = tuple(mask(support) for support in supports)
    answers = []
    for candidate in range(1, FULL_MASK + 1):
        if not all(candidate & edge for edge in hyperedges):
            continue
        if any(
            all((candidate ^ (1 << vertex)) & edge for edge in hyperedges)
            for vertex in vertices(candidate)
        ):
            continue
        answers.append(candidate)
    return tuple(answers)


def rank_nine_vertices(supports: Sequence[Support]) -> tuple[Support, ...]:
    forbidden = tuple(map(frozenset, (*TARGET_MISSING_EDGES, *supports)))
    return tuple(
        chosen
        for chosen in combinations(range(VERTEX_COUNT), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


@lru_cache(maxsize=None)
def placement_options(support: Support) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(support))
    )


@lru_cache(maxsize=524_288)
def finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return is_finite_coxeter_labeling(range(order), labels)


def restriction_finite(restriction: Sequence[int], labels: Mapping[Pair, int]) -> bool:
    return finite_code(
        len(restriction),
        tuple(
            labels[pair(left, right)]
            for left, right in combinations(restriction, 2)
        ),
    )


def next_support(
    supports: Sequence[Support], remaining: set[int], union: set[int]
) -> int:
    if not union:
        return min(remaining)
    return max(
        remaining,
        key=lambda index: (
            len(set(supports[index]) & union),
            -len(set(supports[index]) - union),
            -index,
        ),
    )


def sparse_label_gate(supports: Sequence[Support], state_cap: int) -> dict[str, object]:
    actual_vertices = rank_nine_vertices(supports)
    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[LabelState, ...] = ((),)
    stages = []
    capped = False
    while states and remaining:
        support_index = next_support(supports, remaining, union)
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {
                    tuple(facet for facet in vertex if facet in union)
                    for vertex in actual_vertices
                },
                key=lambda row: (len(row), row),
            )
        )
        options = placement_options(support)
        passing: set[LabelState] = set()
        compatible = 0
        raw = len(states) * len(options)
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
                    if len(passing) > state_cap:
                        capped = True
                        break
            if capped:
                break
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "raw_extension_count": raw,
                "compatible_extension_count_before_cap": compatible,
                "passing_state_count": len(states),
                "restriction_size_histogram": {
                    str(size): count
                    for size, count in sorted(Counter(map(len, restrictions)).items())
                },
                "state_cap_exceeded": capped,
            }
        )
        remaining.remove(support_index)
        if capped:
            break
    return {
        "actual_rank_nine_vertex_count": len(actual_vertices),
        "stages": stages,
        "terminal_sparse_state_count": len(states),
        "state_cap_exceeded": capped,
        "unrestricted_compact_label_no_go": not capped and not states,
    }


@lru_cache(maxsize=1)
def active_automorphisms() -> tuple[tuple[int, ...], ...]:
    edges = {tuple(sorted(edge)) for edge in TARGET_MISSING_EDGES}
    answers = []
    for image in permutations(ACTIVE):
        permutation = dict(zip(ACTIVE, image, strict=True))
        mapped = {
            tuple(sorted((permutation[left], permutation[right])))
            for left, right in edges
        }
        if mapped == edges:
            answers.append(image)
    if len(answers) != 8:
        raise AssertionError("record-18 active automorphism order changed")
    return tuple(answers)


def canonical_supports(supports: Sequence[Support]) -> tuple[Support, ...]:
    best: tuple[Support, ...] | None = None
    for active_image in active_automorphisms():
        active_map = dict(zip(ACTIVE, active_image, strict=True))
        for isolate_image in permutations(ISOLATES):
            permutation = active_map | dict(zip(ISOLATES, isolate_image, strict=True))
            code = tuple(
                sorted(
                    tuple(sorted(permutation[vertex] for vertex in support))
                    for support in supports
                )
            )
            if best is None or code < best:
                best = code
    if best is None:
        raise AssertionError("empty canonicalization orbit")
    return best


def audit_target(fours: Sequence[Support], additions: Sequence[Support]) -> dict[str, object]:
    supports = tuple(fours) + tuple(additions)
    clutter = tuple(TARGET_MISSING_EDGES) + supports
    blocker = minimal_transversals(clutter)
    blocker_histogram = Counter(value.bit_count() for value in blocker)
    global_audit = audit_global_facet_selection(
        TARGET_GRAPH6,
        blocker,
        [{"missing_four_count": 2, "f_vector": EXPECTED_F_VECTOR}],
        blocker,
    )
    return {
        "supports": [list(support) for support in supports],
        "blocker_size_histogram": {
            str(size): count for size, count in sorted(blocker_histogram.items())
        },
        "global_audit": global_audit.to_record(),
        "canonical_supports": [list(support) for support in canonical_supports(supports)],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, nargs="?", default=DEFAULT_INPUT)
    parser.add_argument("--state-cap", type=int, default=250_000)
    args = parser.parse_args()

    document = json.loads(args.input.read_text(encoding="utf-8"))
    if document.get("a") != 2 or document.get("representative_A4_masks") != [53, 198]:
        raise AssertionError("expected record-18 a2 orbit [53,198]")
    fours = tuple(tuple(map(int, row)) for row in document["representative_A4s"])
    additions = tuple(
        tuple(tuple(map(int, row)) for row in family)
        for family in document["scheduled_target_additions"]
    )
    if document.get("scheduled_pure_target_count") != 6 or len(additions) != 6:
        raise AssertionError("the six residual addition families changed")

    audits = tuple(audit_target(fours, family) for family in additions)
    orbit_groups: dict[tuple[Support, ...], list[int]] = {}
    for index, audit in enumerate(audits):
        key = tuple(tuple(row) for row in audit["canonical_supports"])
        orbit_groups.setdefault(key, []).append(index)
    if len(orbit_groups) != 1:
        raise AssertionError("the six residual targets no longer form one orbit")
    canonical, members = next(iter(orbit_groups.items()))
    label_gate = sparse_label_gate(canonical, args.state_cap)

    output = {
        "kind": "h9_14_record18_a2_idx3_residual_compact_label_scout",
        "theorem_certificate": False,
        "independent_replay_complete": False,
        "compact_only": True,
        "finite_volume_noncompact_chamber_excluded": False,
        "source": str(args.input),
        "target_count": len(audits),
        "target_audits": audits,
        "full_graph_automorphism_order": 8 * 40320,
        "target_orbit_count": 1,
        "target_orbits": [
            {
                "members": members,
                "canonical_supports": [list(support) for support in canonical],
                "sparse_label_gate": label_gate,
            }
        ],
        "all_targets_pass_connected_homology_audit": all(
            audit["global_audit"]["passes_connected_homology_audit"]
            for audit in audits
        ),
        "all_target_orbits_have_unrestricted_compact_label_no_go": (
            label_gate["unrestricted_compact_label_no_go"]
        ),
        "scope": (
            "Exact compact-only incidence/topology and finite-label Lanner/"
            "actual-vertex scout. An uncapped zero terminal is an unrestricted "
            "compact Coxeter-label no-go by entrywise Gram monotonicity. No "
            "finite-volume noncompact conclusion and no theorem promotion "
            "without independent replay."
        ),
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
