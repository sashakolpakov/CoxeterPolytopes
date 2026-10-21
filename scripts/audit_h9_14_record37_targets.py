#!/usr/bin/env python3
"""Audit every target-count state from the exact record-37 A5 census."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json

from coxeter_polytopes.combinatorics.facet_completion import (
    audit_global_facet_selection,
)
from scripts.audit_h9_14_four_a5_candidate import (
    _minimal_transversals,
    _packed_face_table_sha256,
)
from scripts.enumerate_h9_14_record37_a0 import (
    SOURCE,
    TARGET_INDEX,
    enumerate_a0,
)


def audit_target_families(census: dict[str, object]) -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [
        record for record in source["records"] if record["survives_second_layer"]
    ]
    record = survivors[TARGET_INDEX]
    edge_masks = tuple(
        sum(1 << vertex for vertex in edge) for edge in census["missing_edges"]
    )
    unique_families = sorted(
        {
            tuple(sorted(state["missing_five_masks"]))
            for state in census["target-states"]
        }
    )
    family_index = {family: index for index, family in enumerate(unique_families)}
    raw_to_unique = [
        {
            "family_index": family_index[
                tuple(sorted(state["missing_five_masks"]))
            ],
            "row_labels": list(state["row_labels"]),
            "isolate_pattern_counts": list(state["isolate_pattern_counts"]),
        }
        for state in census["target-states"]
    ]
    results = []
    for family in unique_families:
        blocker = _minimal_transversals(edge_masks + family)
        histogram = Counter(mask.bit_count() for mask in blocker)
        five_blocker = tuple(mask for mask in blocker if mask.bit_count() == 5)
        audit = audit_global_facet_selection(
            record["graph6"],
            record["ridge_support_core_masks"],
            record["dehn_sommerville_schedule"],
            five_blocker,
        )
        results.append(
            {
                "missing_five_masks": list(family),
                "minimal_transversal_size_histogram": dict(sorted(histogram.items())),
                "facet_complement_masks": list(five_blocker),
                "facet_complement_sha256": sha256(
                    b"".join(
                        mask.to_bytes(2, "little") for mask in five_blocker
                    )
                ).hexdigest(),
                "packed_face_table_sha256": _packed_face_table_sha256(
                    five_blocker
                ),
                "passes_selection_model": audit.passes_selection_model,
                "all_ridge_degrees_two": audit.all_ridge_degrees_two,
                "ridge_degree_histogram": dict(audit.ridge_degree_histogram),
                "dual_component_count": audit.dual_component_count,
                "first_strict_local_g_failure": (
                    list(audit.first_strict_local_g_failure)
                    if audit.first_strict_local_g_failure is not None
                    else None
                ),
                "first_local_euler_failure": (
                    list(audit.first_local_euler_failure)
                    if audit.first_local_euler_failure is not None
                    else None
                ),
                "first_mod2_homology_failure": (
                    [
                        audit.first_mod2_homology_failure[0],
                        audit.first_mod2_homology_failure[1],
                        list(audit.first_mod2_homology_failure[2]),
                    ]
                    if audit.first_mod2_homology_failure is not None
                    else None
                ),
                "global_audit": audit.to_record(),
            }
        )
    return {
        "record_index": TARGET_INDEX,
        "raw_target_state_count": census["target-state-count"],
        "unique_target_family_count": len(unique_families),
        "raw_to_unique_family_map": raw_to_unique,
        "passing_family_count": sum(
            item["passes_selection_model"] for item in results
        ),
        "families": results,
    }


def main() -> None:
    census = enumerate_a0()
    output = audit_target_families(census)
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
