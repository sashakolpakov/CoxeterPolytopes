#!/usr/bin/env python3
"""No-project-import blocker and compact-label replay for Record 31 a=2."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from typing import Mapping

from scripts import independent_h9_14_record32_downstream_audit as hostile
from scripts.independent_h9_14_record31_a3_downstream_audit import (
    fixed_order_sparse_join,
)


N = 14
MISSING_EDGES = ((0, 11), (1, 12), (0, 13), (1, 13))
EDGE_MASKS = tuple(hostile.mask(edge) for edge in MISSING_EDGES)


def alternate_support_order(family: tuple[int, ...]) -> tuple[int, ...]:
    """Overlap-greedy order with opposite start/tie rules from production."""

    supports = tuple(set(hostile.vertices(face)) for face in family)
    remaining = set(range(len(family)))
    first = max(remaining, key=lambda index: (len(supports[index]), index))
    order = [first]
    remaining.remove(first)
    union = set(supports[first])
    while remaining:
        chosen = max(
            remaining,
            key=lambda index: (
                len(supports[index] & union),
                len(supports[index]),
                -len(supports[index] - union),
                index,
            ),
        )
        order.append(chosen)
        remaining.remove(chosen)
        union.update(supports[chosen])
    return tuple(order)


def audit(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record31_a2_working_checkpoint":
        raise ValueError("wrong independent Record 31 a=2 input kind")
    if document.get("remaining_residual_orbit_indices") != []:
        raise ValueError("independent Record 31 a=2 input is incomplete")
    raw_targets = document.get("ridge_valid_target_orbits")
    if not isinstance(raw_targets, list) or len(raw_targets) != 12:
        raise ValueError("independent Record 31 a=2 target list changed")

    hostile.MISSING_EDGES = MISSING_EDGES
    hostile.EDGE_MASKS = EDGE_MASKS
    if (len(hostile.lanner_codes(4)), len(hostile.lanner_codes(5))) != (108, 420):
        raise AssertionError("independent classical Lanner tables changed")

    label_rows = []
    for target_index, raw_family in enumerate(raw_targets):
        if not isinstance(raw_family, list):
            raise ValueError("independent target family is malformed")
        family = tuple(hostile.mask(map(int, face)) for face in raw_family)
        blocker = hostile.minimal_transversals((*EDGE_MASKS, *family))
        if Counter(face.bit_count() for face in blocker) != {5: 118}:
            raise AssertionError("independent a=2 pure blocker replay failed")
        if hostile.ridge_histogram(blocker) != {"2": 531}:
            raise AssertionError("independent a=2 ridge replay failed")
        # Production starts at the first A5 and breaks ties toward lower
        # indices.  This replay starts at the last A5 and breaks ties toward
        # higher indices, then freezes that independently chosen order.
        order = alternate_support_order(family)
        row = fixed_order_sparse_join(family, order)
        row["target_orbit_index"] = target_index
        label_rows.append(row)

    if any(row["terminal_sparse_state_count"] for row in label_rows):
        raise AssertionError("an independent a=2 sparse label state survived")
    return {
        "kind": "independent_h9_14_record31_a2_downstream_audit",
        "compact_only": True,
        "classical_lanner_placement_counts": {"4": 108, "5": 420},
        "literal_target_orbit_count": len(label_rows),
        "ridge_survivor_label_audits": label_rows,
        "all_targets_have_unrestricted_compact_label_no_go": True,
        "compact_label_survivor_count": 0,
        "scope": {
            "no_coxeter_polytopes_project_import": True,
            "fixed_join_order_differs_from_production": True,
            "compact_only": True,
            "finite_volume_noncompact_inference": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    args = parser.parse_args()
    document = json.loads(args.checkpoint.read_text(encoding="utf-8"))
    print(json.dumps(audit(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
