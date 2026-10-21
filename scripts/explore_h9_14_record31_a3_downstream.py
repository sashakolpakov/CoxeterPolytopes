#!/usr/bin/env python3
"""Exact symmetry, topology, and compact-label gates for Record 31 ``a=3``.

The input is the JSON emitted by
``explore_h9_14_record31_a3_additions.cpp``.  Its literal pure blocker
targets are quotiented by the complete missing-edge graph action
``C2 x S9``.  Every quotient representative is dualized again and receives
the exact ridge gate.  Ridge survivors receive the complete global
facet-selection audit (including every-link Euler and mod-two homology), then
an unrestricted sparse compact Coxeter-label gate on all order-four and
order-five Lanner supports.

An uncapped zero sparse state count is unrestricted in every as-yet
unassigned finite label: in normalized Gram form ``G=I-A``, replacing a
right angle by any finite non-right angle only increases the nonnegative
matrix ``A`` and cannot restore ellipticity of a failed principal
restriction.  This script makes no finite-volume or cusp claim.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
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


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_ARCHIVE = ROOT / "data/enumerations/h9_14_nerve_support.json"
GRAPH6 = "M?????????O?O?o??"
N = 14
DIMENSION = 9
FULL_MASK = (1 << N) - 1
MISSING_EDGES = ((0, 11), (1, 12), (0, 13), (1, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in MISSING_EDGES)
ACTIVE = (0, 1, 11, 12, 13)
ISOLATED = tuple(range(2, 11))
ACTIVE_FLIP = {0: 1, 1: 0, 11: 12, 12: 11, 13: 13}

LabelState = tuple[tuple[Pair, int], ...]


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value >> vertex & 1)


def minimal_transversals(hyperedges: Sequence[int]) -> tuple[int, ...]:
    """Build the complete literal blocker on fourteen vertices."""

    result = []
    for candidate in range(1, FULL_MASK + 1):
        if not all(candidate & face for face in hyperedges):
            continue
        if all(
            any(candidate & face == bit for face in hyperedges)
            for bit in (1 << vertex for vertex in vertices(candidate))
        ):
            result.append(candidate)
    return tuple(result)


def ridge_degree_histogram(
    facet_complements: Sequence[int],
) -> dict[str, int]:
    ridges: Counter[int] = Counter()
    for complement in facet_complements:
        outside = FULL_MASK ^ complement
        while outside:
            bit = outside & -outside
            ridges[complement | bit] += 1
            outside ^= bit
    return {
        str(degree): multiplicity
        for degree, multiplicity in sorted(Counter(ridges.values()).items())
    }


def _flip_active(face: int) -> int:
    active_mask = mask(ACTIVE)
    return (face & ~active_mask) | mask(
        ACTIVE_FLIP[vertex] for vertex in ACTIVE if face >> vertex & 1
    )


def canonical_family_key(family: Sequence[int]) -> tuple[int, ...]:
    """Complete invariant for the ``C2 x S9`` graph action.

    Rows of equal order are unordered.  After either active ``P5``
    automorphism is applied, the nine isolated vertices are the unordered
    binary columns of the row-incidence matrix.
    """

    groups = {
        order: tuple(face for face in family if face.bit_count() == order)
        for order in (4, 5)
    }
    if len(groups[4]) != 3 or not 1 <= len(groups[5]) <= 4:
        raise ValueError("Record 31 a=3 target family shape changed")
    best: tuple[int, ...] | None = None
    for flip in (False, True):
        moved = {
            order: tuple(_flip_active(face) if flip else face for face in groups[order])
            for order in (4, 5)
        }
        for ordered_fours in permutations(moved[4]):
            for ordered_fives in permutations(moved[5]):
                rows = ordered_fours + ordered_fives
                active_rows = tuple(
                    sum(((face >> vertex) & 1) << index
                        for index, vertex in enumerate(ACTIVE))
                    for face in rows
                )
                isolate_columns = tuple(
                    sorted(
                        sum(((face >> vertex) & 1) << index
                            for index, face in enumerate(rows))
                        for vertex in ISOLATED
                    )
                )
                candidate = active_rows + isolate_columns
                if best is None or candidate < best:
                    best = candidate
    if best is None:
        raise AssertionError("canonicalization produced no image")
    return best


def production_literal_families(
    document: Mapping[str, object],
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    rows = document.get("representative_audits")
    if not isinstance(rows, list) or len(rows) != 143:
        raise ValueError("Record 31 a=3 production row coverage changed")
    result = []
    for expected_index, raw in enumerate(rows):
        if not isinstance(raw, Mapping):
            raise ValueError("Record 31 production row is malformed")
        index = int(raw["orbit_representative_index"])
        if index != expected_index:
            raise ValueError("Record 31 production indices changed")
        raw_families = raw.get("pure_target_minimal_nonfaces")
        if not isinstance(raw_families, list):
            raise ValueError("Record 31 production target list is malformed")
        if int(raw["pure_116_facet_target_count"]) != len(raw_families):
            raise ValueError("Record 31 production target count changed")
        for raw_family in raw_families:
            if not isinstance(raw_family, list):
                raise ValueError("Record 31 target family is malformed")
            family = tuple(
                sorted(mask(map(int, face)) for face in raw_family)
            )
            result.append((index, family))
    if (
        int(document.get("pure_116_facet_target_count_across_representatives", -1))
        != len(result)
    ):
        raise ValueError("Record 31 production total changed")
    return tuple(result)


def literal_family_sha256(
    families: Sequence[tuple[int, tuple[int, ...]]],
) -> str:
    payload = json.dumps(
        [[index, list(family)] for index, family in families],
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return sha256(payload).hexdigest()


def support_record() -> Mapping[str, object]:
    document = json.loads(SUPPORT_ARCHIVE.read_text(encoding="utf-8"))
    rows = document.get("records")
    if not isinstance(rows, list):
        raise ValueError("nerve-support archive changed")
    matches = [row for row in rows if row.get("graph6") == GRAPH6]
    if len(matches) != 1:
        raise ValueError("Record 31 support row is not unique")
    return matches[0]


def rank_nine_vertices(minimal_faces: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(frozenset(vertices(face)) for face in (*EDGE_MASKS, *minimal_faces))
    return tuple(
        chosen
        for chosen in combinations(range(N), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


@lru_cache(maxsize=None)
def placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
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


def restriction_finite(
    restriction: Sequence[int], labels: Mapping[Pair, int]
) -> bool:
    return finite_code(
        len(restriction),
        tuple(
            labels[pair(left, right)]
            for left, right in combinations(restriction, 2)
        ),
    )


def _next_support(
    supports: Sequence[tuple[int, ...]], remaining: set[int], union: set[int]
) -> int:
    # Start with a five-support; thereafter maximize overlap.  This is merely
    # a join order and does not alter the exhaustive state set.
    if not union:
        return max(remaining, key=lambda index: (len(supports[index]), -index))
    return max(
        remaining,
        key=lambda index: (
            len(set(supports[index]) & union),
            len(supports[index]),
            -len(set(supports[index]) - union),
            -index,
        ),
    )


def sparse_label_gate(minimal_faces: Sequence[int]) -> dict[str, object]:
    supports = tuple(vertices(face) for face in minimal_faces)
    actual_vertices = rank_nine_vertices(minimal_faces)
    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[LabelState, ...] = ((),)
    stages = []
    while states and remaining:
        support_index = _next_support(supports, remaining, union)
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
        raw_extension_count = len(states) * len(options)
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
                "raw_extension_count": raw_extension_count,
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
                "restriction_size_histogram": {
                    str(size): multiplicity
                    for size, multiplicity in sorted(
                        Counter(map(len, restrictions)).items()
                    )
                },
            }
        )
        remaining.remove(support_index)
    return {
        "rank_nine_vertex_count": len(actual_vertices),
        "stages": stages,
        "terminal_sparse_state_count": len(states),
        "state_cap_used": False,
        "unrestricted_compact_label_no_go": not states,
    }


def audit(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record31_a3_pivoted_addition_scout":
        raise ValueError("wrong Record 31 a=3 production kind")
    families = production_literal_families(document)
    grouped: dict[tuple[int, ...], list[tuple[int, tuple[int, ...]]]] = defaultdict(list)
    for row in families:
        grouped[canonical_family_key(row[1])].append(row)

    support = support_record()
    support_masks = support.get("ridge_support_core_masks")
    schedules = support.get("dehn_sommerville_schedule")
    if not isinstance(support_masks, list) or not isinstance(schedules, list):
        raise ValueError("Record 31 support metadata changed")

    orbit_rows = []
    for orbit_index, key in enumerate(sorted(grouped)):
        members = grouped[key]
        representative_index, family = min(members)
        blocker = minimal_transversals((*EDGE_MASKS, *family))
        blocker_histogram = Counter(value.bit_count() for value in blocker)
        ridge_histogram = ridge_degree_histogram(blocker)
        ridge_valid = blocker_histogram == {5: 116} and ridge_histogram == {"2": 522}
        row: dict[str, object] = {
            "target_orbit_index": orbit_index,
            "canonical_key": list(key),
            "literal_multiplicity": len(members),
            "source_representative_indices": sorted(
                Counter(index for index, _family in members).items()
            ),
            "representative_source_index": representative_index,
            "minimal_nonface_masks": list(family),
            "minimal_nonfaces": [list(vertices(face)) for face in family],
            "minimal_transversal_size_histogram": {
                str(order): multiplicity
                for order, multiplicity in sorted(blocker_histogram.items())
            },
            "ridge_degree_histogram": ridge_histogram,
            "ridge_pseudomanifold": ridge_valid,
        }
        if ridge_valid:
            topology = audit_global_facet_selection(
                GRAPH6, support_masks, schedules, blocker
            )
            row["global_topology_audit"] = topology.to_record()
            if topology.passes_connected_homology_audit:
                row["sparse_compact_label_gate"] = sparse_label_gate(family)
        orbit_rows.append(row)

    ridge_rows = [row for row in orbit_rows if row["ridge_pseudomanifold"]]
    topology_rows = [
        row
        for row in ridge_rows
        if row.get("global_topology_audit", {}).get(
            "passes_connected_homology_audit", False
        )
    ]
    label_no_go_rows = [
        row
        for row in topology_rows
        if row.get("sparse_compact_label_gate", {}).get(
            "unrestricted_compact_label_no_go", False
        )
    ]
    return {
        "kind": "h9_14_record31_a3_downstream_compact_audit",
        "compact_only": True,
        "literal_target_count": len(families),
        "literal_family_sha256": literal_family_sha256(families),
        "full_graph_action": {
            "active_automorphism_order": 2,
            "isolated_vertex_symmetric_group_degree": 9,
            "description": "C2 x S9",
        },
        "target_orbit_count": len(orbit_rows),
        "ridge_pseudomanifold_orbit_count": len(ridge_rows),
        "connected_mod2_homology_sphere_orbit_count": len(topology_rows),
        "unrestricted_compact_label_no_go_orbit_count": len(label_no_go_rows),
        "compact_label_survivor_orbit_count": len(topology_rows) - len(label_no_go_rows),
        "target_orbits": orbit_rows,
        "scope": {
            "compact_only": True,
            "ridge_failure_is_only_a_pseudomanifold_rejection": True,
            "mod2_homology_sphere_is_not_claimed_to_prove_polytopality": True,
            "zero_uncapped_sparse_state_is_an_unrestricted_compact_label_no_go": True,
            "finite_volume_noncompact_inference": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("production_json", type=Path)
    args = parser.parse_args()
    document = json.loads(args.production_json.read_text(encoding="utf-8"))
    print(json.dumps(audit(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
