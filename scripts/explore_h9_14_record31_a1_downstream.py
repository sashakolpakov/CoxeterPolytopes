#!/usr/bin/env python3
"""Exact topology and compact-label gates for Record 31 compact a=1 targets."""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import permutations
import json
from pathlib import Path
from typing import Mapping

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from scripts.explore_h9_14_record31_a3_downstream import (
    EDGE_MASKS,
    GRAPH6,
    mask,
    minimal_transversals,
    ridge_degree_histogram,
    sparse_label_gate,
    support_record,
    vertices,
)


EXPECTED_F_VECTOR = (14, 87, 319, 770, 1281, 1477, 1142, 540, 120)
ACTIVE = (0, 1, 11, 12, 13)
ISOLATES = tuple(range(2, 11))
ACTIVE_FLIP = {0: 1, 1: 0, 11: 12, 12: 11, 13: 13}


def _flip_active(face: int) -> int:
    active_mask = mask(ACTIVE)
    return (face & ~active_mask) | mask(
        ACTIVE_FLIP[vertex] for vertex in ACTIVE if face >> vertex & 1
    )


def canonical_family_key(family: tuple[int, ...]) -> tuple[int, ...]:
    """Complete invariant for one A4 and five-to-seven A5s under C2 x S9."""

    fours = tuple(face for face in family if face.bit_count() == 4)
    fives = tuple(face for face in family if face.bit_count() == 5)
    if len(fours) != 1 or not 5 <= len(fives) <= 7:
        raise ValueError("Record 31 a=1 target family shape changed")
    best: tuple[int, ...] | None = None
    for flip in (False, True):
        moved_four = _flip_active(fours[0]) if flip else fours[0]
        moved_fives = tuple(_flip_active(face) if flip else face for face in fives)
        for ordered_fives in permutations(moved_fives):
            rows = (moved_four, *ordered_fives)
            active_rows = tuple(
                sum(
                    ((face >> vertex) & 1) << index
                    for index, vertex in enumerate(ACTIVE)
                )
                for face in rows
            )
            isolate_columns = tuple(
                sorted(
                    sum(
                        ((face >> vertex) & 1) << index
                        for index, face in enumerate(rows)
                    )
                    for vertex in ISOLATES
                )
            )
            candidate = active_rows + isolate_columns
            if best is None or candidate < best:
                best = candidate
    if best is None:
        raise AssertionError("Record 31 a=1 canonicalization failed")
    return best


def audit(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record31_a1_pivoted_addition_scout":
        raise ValueError("wrong Record 31 a=1 incidence kind")
    raw_rows = document.get("orbits")
    if not isinstance(raw_rows, list) or not raw_rows:
        raise ValueError("Record 31 a=1 incidence rows changed")

    support = support_record()
    support_masks = support.get("ridge_support_core_masks")
    schedules = support.get("dehn_sommerville_schedule")
    if not isinstance(support_masks, list) or not isinstance(schedules, list):
        raise ValueError("Record 31 support archive changed")

    rows = []
    target_index = 0
    for raw_source in raw_rows:
        if not isinstance(raw_source, Mapping):
            raise ValueError("Record 31 a=1 source row is malformed")
        source_index = int(raw_source["orbit_index"])
        raw_targets = raw_source.get("ridge_valid_target_minimal_nonfaces")
        if not isinstance(raw_targets, list):
            raise ValueError("Record 31 a=1 target list changed")
        if int(raw_source["ridge_valid_target_orbit_count"]) != len(raw_targets):
            raise ValueError("Record 31 a=1 target count changed")
        for raw_family in raw_targets:
            if not isinstance(raw_family, list):
                raise ValueError("Record 31 a=1 family is malformed")
            family = tuple(mask(map(int, face)) for face in raw_family)
            if len([face for face in family if face.bit_count() == 4]) != 1:
                raise ValueError("Record 31 a=1 A4 count changed")
            if not 5 <= len(
                [face for face in family if face.bit_count() == 5]
            ) <= 7:
                raise ValueError("Record 31 a=1 added-A5 count changed")

            blocker = minimal_transversals((*EDGE_MASKS, *family))
            blocker_histogram = Counter(face.bit_count() for face in blocker)
            ridge_histogram = ridge_degree_histogram(blocker)
            if blocker_histogram != {5: 120} or ridge_histogram != {"2": 540}:
                raise ValueError("Record 31 a=1 ridge target failed literal replay")
            topology = audit_global_facet_selection(
                GRAPH6, support_masks, schedules, blocker
            )
            row: dict[str, object] = {
                "target_orbit_index": target_index,
                "source_positive_A4_orbit_index": source_index,
                "minimal_nonface_masks": list(family),
                "minimal_nonfaces": [list(vertices(face)) for face in family],
                "canonical_family_key": list(canonical_family_key(family)),
                "minimal_transversal_size_histogram": {"5": 120},
                "ridge_degree_histogram": {"2": 540},
                "global_topology_audit": topology.to_record(),
            }
            if topology.passes_connected_homology_audit:
                row["sparse_compact_label_gate"] = sparse_label_gate(family)
            rows.append(row)
            target_index += 1

    topology_rows = [
        row
        for row in rows
        if row["global_topology_audit"]["passes_connected_homology_audit"]
    ]
    label_no_go_rows = [
        row
        for row in topology_rows
        if row.get("sparse_compact_label_gate", {}).get(
            "unrestricted_compact_label_no_go", False
        )
    ]
    return {
        "kind": "h9_14_record31_a1_downstream_compact_audit",
        "compact_only": True,
        "schedule": {
            "a": 1,
            "f_vector": list(EXPECTED_F_VECTOR),
            "facet_count": 120,
            "ridge_count": 540,
        },
        "input_target_orbit_count": len(rows),
        "connected_mod2_homology_sphere_orbit_count": len(topology_rows),
        "unrestricted_compact_label_no_go_orbit_count": len(label_no_go_rows),
        "compact_label_survivor_orbit_count": len(topology_rows) - len(label_no_go_rows),
        "target_orbits": rows,
        "scope": {
            "compact_only": True,
            "mod2_homology_sphere_is_not_claimed_to_prove_polytopality": True,
            "zero_uncapped_sparse_state_is_an_unrestricted_compact_label_no_go": True,
            "finite_volume_noncompact_inference": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("incidence", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    document = json.loads(args.incidence.read_text(encoding="utf-8"))
    payload = json.dumps(audit(document), indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
