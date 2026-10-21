#!/usr/bin/env python3
"""Close the sole strict-shadow orbit in the record-32 ``a=3`` branch.

The representative is the unique orbit returned by the exact positive-branch
scout.  This file exhausts its two additional minimal A5s, applies the exact
order-at-most-four transversal and scheduled-facet gates, and sends any
scheduled literal target through blocker, ridge, and unrestricted sparse
Lanner/actual-vertex checks.

This remains a scout: the orbit reduction and any label transcript require an
independent no-project-import replay before they can support a theorem.  No
finite-volume or noncompact conclusion is made.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json
from math import comb

from coxeter_polytopes.combinatorics.h9_14_record32 import TARGET_GRAPH6
from coxeter_polytopes.combinatorics.h9_14_record35_all_completions import (
    _forced_five_shadow,
    _shadow_gate,
)
from coxeter_polytopes.combinatorics.h9_14_record46 import _base_link_g3_g4
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
)
from scripts.explore_h9_14_record32_a0_target_gate import (
    EDGE_MASKS,
    minimal_transversals,
    ridge_degree_histogram,
    sparse_label_gate,
)
from scripts.explore_h9_14_record32_a2_additions import (
    at_least_masks,
    hit_word,
    mask,
    small_base_transversals,
    vertices,
)


VERTEX_COUNT = 14
TOTAL_FIVES = 24
SCHEDULED_FACETS = 108
A3_REPRESENTATIVE = (
    (1, 2, 3, 4),
    (2, 3, 4, 11),
    (3, 4, 11, 12),
)


def audit_a3_representative(state_cap: int = 250_000) -> dict[str, object]:
    family = tuple(mask(face) for face in A3_REPRESENTATIVE)
    fives = independent_masks(TARGET_GRAPH6, 5)
    base_g3, base_g4 = _base_link_g3_g4(TARGET_GRAPH6)
    shadow = _forced_five_shadow(family, fives)
    strict, aggregate, additional, lower, upper = _shadow_gate(
        family, shadow, base_g3, base_g4, TOTAL_FIVES
    )
    if not strict or not aggregate or additional != 2:
        raise AssertionError("record-32 a=3 representative gate changed")

    candidates = tuple(
        five
        for five in fives
        if five not in shadow
        and not any(five & four == four for four in family)
        and all(
            not (five & (1 << vertex)) or upper[vertex] >= 1
            for vertex in range(VERTEX_COUNT)
        )
    )
    lower_at_least = tuple(
        mask(vertex for vertex in range(VERTEX_COUNT) if lower[vertex] >= threshold)
        for threshold in range(1, additional + 1)
    )
    upper_below = tuple(
        mask(vertex for vertex in range(VERTEX_COUNT) if upper[vertex] < threshold)
        for threshold in range(1, additional + 1)
    )

    small = small_base_transversals((*EDGE_MASKS, *family, *shadow))
    small_words = {five: hit_word(five, small) for five in candidates}
    five_covers = graph_vertex_cover_masks(TARGET_GRAPH6, 5)
    five_words = {face: hit_word(face, five_covers) for face in (*family, *candidates)}
    fixed_five_word = (1 << len(five_covers)) - 1
    for four in family:
        fixed_five_word &= five_words[four]

    cap_feasible = 0
    no_small = 0
    facet_histogram: Counter[int] = Counter()
    targets: list[tuple[int, ...]] = []
    for additions in combinations(candidates, additional):
        present = at_least_masks(additions)
        if any(
            required & actual != required or forbidden & actual
            for required, forbidden, actual in zip(
                lower_at_least, upper_below, present, strict=True
            )
        ):
            continue
        cap_feasible += 1
        small_word = (1 << len(small)) - 1
        for face in additions:
            small_word &= small_words[face]
        if small_word:
            continue
        no_small += 1
        word = fixed_five_word
        for face in additions:
            word &= five_words[face]
        facet_count = word.bit_count()
        facet_histogram[facet_count] += 1
        if facet_count == SCHEDULED_FACETS:
            targets.append(tuple(additions))

    target_audits = []
    for additions in targets:
        minimal_faces = (*family, *additions)
        transversals = minimal_transversals((*EDGE_MASKS, *minimal_faces))
        transversal_histogram = tuple(
            sorted(Counter(value.bit_count() for value in transversals).items())
        )
        five_transversals = tuple(
            value for value in transversals if value.bit_count() == 5
        )
        ridge_histogram = ridge_degree_histogram(five_transversals)
        pure_ridge = (
            transversal_histogram == ((5, SCHEDULED_FACETS),)
            and bool(ridge_histogram)
            and set(dict(ridge_histogram)) == {2}
        )
        target_audits.append(
            {
                "additional_A5s": [list(vertices(face)) for face in additions],
                "minimal_transversal_size_histogram": dict(transversal_histogram),
                "ridge_degree_histogram": dict(ridge_histogram),
                "pure_ridge_pseudomanifold": pure_ridge,
                "sparse_label_gate": (
                    sparse_label_gate(minimal_faces, state_cap) if pure_ridge else None
                ),
            }
        )

    return {
        "missing_fours": [list(face) for face in A3_REPRESENTATIVE],
        "forced_five_shadow_count": len(shadow),
        "additional_minimal_A5_count": additional,
        "individual_candidate_count": len(candidates),
        "raw_unordered_addition_count": comb(len(candidates), additional),
        "vertex_cap_feasible_count": cap_feasible,
        "no_transversal_of_order_at_most_four_count": no_small,
        "facet_count_histogram": dict(sorted(facet_histogram.items())),
        "scheduled_target_count": len(targets),
        "target_audits": target_audits,
    }


def main() -> None:
    audit = audit_a3_representative()
    print(
        json.dumps(
            {
                "kind": "h9_14_record32_a3_additional_A5_scout",
                "theorem_certificate": False,
                "strict_shadow_orbit_count": 1,
                "representative_audit": audit,
                "scope": (
                    "Exact on the sole strict-shadow orbit representative; the "
                    "orbit reduction and label search are not independently "
                    "replayed. Compact-only, with no finite-volume conclusion."
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
