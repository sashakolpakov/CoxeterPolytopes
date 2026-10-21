#!/usr/bin/env python3
"""Exact topology and unrestricted compact-label gates for Record 31 a=2.

The input is the crash-safe production checkpoint emitted by the independent
incidence work.  Every stored full-action ridge target is dualized literally,
checked against the scheduled blocker and ridge counts, sent through the
complete global face/link homology audit, and finally through the uncapped
Lanner/actual-rank-nine-vertex join.  A zero sparse join is unrestricted by
Perron--Frobenius monotonicity.  This script makes no noncompact finite-volume
or Vinberg endpoint claim.
"""

from __future__ import annotations

import argparse
from collections import Counter
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


EXPECTED_F_VECTOR = (14, 87, 319, 769, 1275, 1463, 1126, 531, 118)


def audit(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record31_a2_working_checkpoint":
        raise ValueError("wrong Record 31 a=2 checkpoint kind")
    if document.get("remaining_residual_orbit_indices") != []:
        raise ValueError("Record 31 a=2 incidence census is not complete")
    if len(document.get("strict_face_link_excluded_orbit_indices", [])) != 31:
        raise ValueError("Record 31 a=2 strict-link partition changed")
    completed = document.get("completed_residual_orbit_indices")
    if not isinstance(completed, list) or len(completed) != 18:
        raise ValueError("Record 31 a=2 residual coverage changed")
    if sorted(
        map(int, document["strict_face_link_excluded_orbit_indices"] + completed)
    ) != list(range(49)):
        raise ValueError("Record 31 a=2 positive-A4 partition is not exhaustive")

    schedule = document.get("schedule")
    if not isinstance(schedule, Mapping) or tuple(schedule.get("f_vector", ())) != EXPECTED_F_VECTOR:
        raise ValueError("Record 31 a=2 schedule changed")
    raw_targets = document.get("ridge_valid_target_orbits")
    if not isinstance(raw_targets, list) or len(raw_targets) != 12:
        raise ValueError("Record 31 a=2 ridge target list changed")

    support = support_record()
    support_masks = support.get("ridge_support_core_masks")
    schedules = support.get("dehn_sommerville_schedule")
    if not isinstance(support_masks, list) or not isinstance(schedules, list):
        raise ValueError("Record 31 support archive changed")

    rows = []
    for target_index, raw_family in enumerate(raw_targets):
        if not isinstance(raw_family, list):
            raise ValueError("Record 31 a=2 family is malformed")
        family = tuple(mask(map(int, face)) for face in raw_family)
        if len([face for face in family if face.bit_count() == 4]) != 2:
            raise ValueError("Record 31 a=2 A4 count changed")
        if not 3 <= len([face for face in family if face.bit_count() == 5]) <= 6:
            raise ValueError("Record 31 a=2 added-A5 count changed")
        blocker = minimal_transversals((*EDGE_MASKS, *family))
        blocker_histogram = Counter(face.bit_count() for face in blocker)
        ridge_histogram = ridge_degree_histogram(blocker)
        if blocker_histogram != {5: 118} or ridge_histogram != {"2": 531}:
            raise ValueError("archived Record 31 a=2 ridge target failed replay")
        topology = audit_global_facet_selection(
            GRAPH6, support_masks, schedules, blocker
        )
        row: dict[str, object] = {
            "target_orbit_index": target_index,
            "minimal_nonface_masks": list(family),
            "minimal_nonfaces": [list(vertices(face)) for face in family],
            "minimal_transversal_size_histogram": {"5": 118},
            "ridge_degree_histogram": {"2": 531},
            "global_topology_audit": topology.to_record(),
        }
        if topology.passes_connected_homology_audit:
            row["sparse_compact_label_gate"] = sparse_label_gate(family)
        rows.append(row)

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
        "kind": "h9_14_record31_a2_downstream_compact_audit",
        "compact_only": True,
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
    parser.add_argument("checkpoint", type=Path)
    args = parser.parse_args()
    document = json.loads(args.checkpoint.read_text(encoding="utf-8"))
    print(json.dumps(audit(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
