#!/usr/bin/env python3
"""Independent literal downstream replay for Record 31 compact ``a=3``.

This audit consumes the output of the independent least-destroyer C++
enumerator.  It deliberately does *not* quotient the 373 targets: every
literal family is dualized again and receives the ridge gate.  For the three
literal ridge survivors it uses the separately implemented classical Lanner
tables and finite-diagram classifier from the no-project-import Record-32
audit, with different fixed support orders from production, and exhausts the
actual-vertex ellipticity join.

Only uncapped empty sparse joins are promoted to unrestricted compact-label
no-gos.  No finite-volume noncompact conclusion is made.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
from typing import Mapping, Sequence

from scripts import independent_h9_14_record32_downstream_audit as hostile


N = 14
DIMENSION = 9
FULL_MASK = (1 << N) - 1
MISSING_EDGES = ((0, 11), (1, 12), (0, 13), (1, 13))
EDGE_MASKS = tuple(hostile.mask(edge) for edge in MISSING_EDGES)

# The production join starts on support 1 in all three cases.  These fixed
# replay orders start on the last support and give the same empty terminal by
# a distinct state-join transcript.
EXPECTED_RIDGE_FAMILIES_AND_ORDERS = {
    (60, 1984, 6212, 8248, 10112): (4, 1, 2, 3, 0),
    (60, 1984, 6336, 9984): (3, 1, 2, 0),
    (60, 1984, 7712, 8220, 8640): (4, 1, 2, 3, 0),
}


def literal_families(
    document: Mapping[str, object],
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    if document.get("kind") != "h9_14_record31_a3_independent_partition_audit":
        raise ValueError("wrong independent Record 31 a=3 kind")
    rows = document.get("representative_audits")
    if not isinstance(rows, list) or len(rows) != 143:
        raise ValueError("independent representative coverage changed")
    result = []
    for expected_index, raw in enumerate(rows):
        if not isinstance(raw, Mapping):
            raise ValueError("independent row is malformed")
        index = int(raw["orbit_representative_index"])
        if index != expected_index:
            raise ValueError("independent row order changed")
        raw_families = raw.get("pure_target_minimal_nonface_masks")
        if not isinstance(raw_families, list):
            raise ValueError("independent target list is malformed")
        if int(raw["pure_116_facet_target_count"]) != len(raw_families):
            raise ValueError("independent target count changed")
        for raw_family in raw_families:
            family = tuple(sorted(map(int, raw_family)))
            if len(family) not in (4, 5, 6, 7):
                raise ValueError("independent target family shape changed")
            result.append((index, family))
    if (
        int(document.get("pure_116_facet_target_count_across_representatives", -1))
        != len(result)
    ):
        raise ValueError("independent target total changed")
    return tuple(result)


def fixed_order_sparse_join(
    family: Sequence[int], support_order: Sequence[int]
) -> dict[str, object]:
    supports = tuple(hostile.vertices(face) for face in family)
    if sorted(support_order) != list(range(len(supports))):
        raise ValueError("fixed support order is not a permutation")
    actual_vertices = hostile.rank_nine_vertices(family)
    states: tuple[hostile.LabelState, ...] = ((),)
    union: set[int] = set()
    stages = []
    for support_index in support_order:
        if not states:
            break
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = hostile.cluster_restrictions(actual_vertices, cluster)
        passing: set[hostile.LabelState] = set()
        compatible = 0
        raw_extensions = len(states) * len(hostile.placement_options(support))
        for state in states:
            old = dict(state)
            for option in hostile.placement_options(support):
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    hostile.pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(
                    hostile.restriction_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "raw_extension_count": raw_extensions,
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
    return {
        "minimal_nonface_masks": list(family),
        "fixed_support_order": list(support_order),
        "rank_nine_vertex_count": len(actual_vertices),
        "stages": stages,
        "terminal_sparse_state_count": len(states),
        "state_cap_used": False,
        "unrestricted_compact_label_no_go": not states,
    }


def audit(document: Mapping[str, object]) -> dict[str, object]:
    # Redirect the reusable no-project-import primitives to Record 31 before
    # any cached actual-vertex or Lanner state is requested.
    hostile.MISSING_EDGES = MISSING_EDGES
    hostile.EDGE_MASKS = EDGE_MASKS
    if (len(hostile.lanner_codes(4)), len(hostile.lanner_codes(5))) != (108, 420):
        raise AssertionError("classical Lanner placement counts changed")

    families = literal_families(document)
    ridge_survivors = []
    blocker_histograms: Counter[tuple[tuple[int, int], ...]] = Counter()
    ridge_failure_histograms: Counter[str] = Counter()
    for source_index, family in families:
        blocker = hostile.minimal_transversals((*EDGE_MASKS, *family))
        blocker_histogram = Counter(face.bit_count() for face in blocker)
        blocker_histograms[tuple(sorted(blocker_histogram.items()))] += 1
        ridge = hostile.ridge_histogram(blocker)
        if blocker_histogram == {5: 116} and ridge == {"2": 522}:
            ridge_survivors.append((source_index, family))
        else:
            ridge_failure_histograms[json.dumps(ridge, sort_keys=True)] += 1

    expected_families = set(EXPECTED_RIDGE_FAMILIES_AND_ORDERS)
    observed_families = {family for _source, family in ridge_survivors}
    if observed_families != expected_families or len(ridge_survivors) != 3:
        raise AssertionError("independent Record 31 ridge survivor set changed")
    label_rows = [
        {
            "source_representative_index": source_index,
            **fixed_order_sparse_join(
                family, EXPECTED_RIDGE_FAMILIES_AND_ORDERS[family]
            ),
        }
        for source_index, family in sorted(ridge_survivors)
    ]
    if any(row["terminal_sparse_state_count"] for row in label_rows):
        raise AssertionError("a Record 31 ridge survivor retained a sparse state")

    return {
        "kind": "independent_h9_14_record31_a3_literal_downstream_audit",
        "compact_only": True,
        "classical_lanner_placement_counts": {"4": 108, "5": 420},
        "literal_target_count": len(families),
        "minimal_transversal_histogram_census": {
            json.dumps(dict(histogram), sort_keys=True): multiplicity
            for histogram, multiplicity in sorted(blocker_histograms.items())
        },
        "ridge_pseudomanifold_literal_count": len(ridge_survivors),
        "ridge_failure_histogram_census": dict(sorted(ridge_failure_histograms.items())),
        "ridge_survivor_label_audits": label_rows,
        "all_ridge_survivors_have_unrestricted_compact_label_no_go": True,
        "compact_label_survivor_count": 0,
        "scope": {
            "every_literal_target_replayed_without_symmetry_quotient": True,
            "compact_only": True,
            "ridge_failure_is_only_a_pseudomanifold_rejection": True,
            "finite_volume_noncompact_inference": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("independent_json", type=Path)
    args = parser.parse_args()
    document = json.loads(args.independent_json.read_text(encoding="utf-8"))
    print(json.dumps(audit(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
