#!/usr/bin/env python3
"""Independent topology and compact-Coxeter replay for record-18 a=1."""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
import importlib.util
import json
from pathlib import Path
from typing import Mapping, Sequence

ROOT = Path(__file__).resolve().parents[1]
BASE_SOURCE = ROOT / "scripts/independent_h9_14_record18_a0_downstream_audit.py"
BASE_SPEC = importlib.util.spec_from_file_location(
    "record18_a0_independent_primitives", BASE_SOURCE
)
if BASE_SPEC is None or BASE_SPEC.loader is None:
    raise RuntimeError("could not load independent topology primitives")
base = importlib.util.module_from_spec(BASE_SPEC)
BASE_SPEC.loader.exec_module(base)

DEFAULT_INPUT = (
    ROOT
    / "data/enumerations/h9_14_record18_high_addition_independent_audit.json"
)
DEFAULT_OUTPUT = (
    ROOT / "data/enumerations/h9_14_record18_a1_downstream_independent_audit.json"
)
A4_MASK = 53
EXPECTED_F_VECTOR = (14, 86, 310, 734, 1198, 1358, 1036, 486, 108)

Pair = tuple[int, int]
LabelState = tuple[tuple[Pair, int], ...]


LANNER_REPRESENTATIVES = {
    4: (
        (((0, 1), 3), ((0, 2), 3), ((0, 3), 5)),
        (((0, 1), 5), ((0, 2), 3), ((1, 3), 3)),
        (((0, 1), 3), ((0, 2), 4), ((1, 3), 5)),
        (((0, 1), 3), ((0, 2), 5), ((1, 3), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 3), 4)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 3), 5)),
        (((0, 1), 3), ((0, 2), 4), ((1, 3), 4), ((2, 3), 3)),
        (((0, 1), 3), ((0, 2), 4), ((1, 3), 5), ((2, 3), 3)),
        (((0, 1), 3), ((0, 2), 5), ((1, 3), 5), ((2, 3), 3)),
    ),
    5: base.LANNER_FIVE_REPRESENTATIVES,
}


@lru_cache(maxsize=None)
def lanner_codes(order: int) -> tuple[tuple[int, ...], ...]:
    edges = tuple(combinations(range(order), 2))
    codes = set()
    for representative in LANNER_REPRESENTATIVES[order]:
        labels = {edge: 2 for edge in edges}
        labels.update(representative)
        for move in permutations(range(order)):
            moved = {
                base.pair(move[left], move[right]): label
                for (left, right), label in labels.items()
            }
            codes.add(tuple(moved[edge] for edge in edges))
    expected = {4: 108, 5: 420}[order]
    if len(codes) != expected:
        raise AssertionError(f"rank-{order} Lanner placement count changed")
    return tuple(sorted(codes))


@lru_cache(maxsize=None)
def placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_codes(len(support))
    )


def next_support(
    supports: Sequence[tuple[int, ...]], remaining: set[int], union: set[int]
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


def restriction_finite(
    restriction: Sequence[int], labels: Mapping[Pair, int]
) -> bool:
    return base.finite_code(
        len(restriction),
        tuple(
            labels[base.pair(left, right)]
            for left, right in combinations(restriction, 2)
        ),
    )


def label_gate(
    supports: Sequence[tuple[int, ...]], complements: Sequence[int]
) -> dict[str, object]:
    actual_vertices = tuple(
        tuple(
            vertex
            for vertex in range(base.VERTEX_COUNT)
            if not complement >> vertex & 1
        )
        for complement in complements
    )
    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[LabelState, ...] = ((),)
    stages = []
    while states and remaining:
        support_index = next_support(supports, remaining, union)
        support = supports[support_index]
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
        options = placement_options(support)
        raw = len(states) * len(options)
        compatible = 0
        passing: set[LabelState] = set()
        for state in states:
            old = dict(state)
            for option in options:
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    base.pair(left, right): 2
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
                "lanner_placement_count": len(options),
                "raw_extension_count": raw,
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
                "restriction_size_histogram": dict(
                    sorted(Counter(map(len, restrictions)).items())
                ),
            }
        )
        remaining.remove(support_index)
    return {
        "actual_rank_nine_vertex_count": len(actual_vertices),
        "stages": stages,
        "processed_support_count": len(stages),
        "terminal_state_count": len(states),
        "unrestricted_compact_coxeter_no_go": not states,
        "monotonicity_basis": (
            "Unassigned labels are 2. For G=I-A, the symmetric matrix A is "
            "entrywise nonnegative; increasing a Coxeter label cannot "
            "decrease its Perron root, so a failed elliptic restriction "
            "cannot be repaired by later labels."
        ),
    }


def topology_audit(
    additions: Sequence[int],
) -> tuple[dict[str, object], tuple[int, ...]]:
    clutter = tuple(sorted((*base.EDGE_MASKS, A4_MASK, *additions)))
    complements = base.minimal_transversals(clutter)
    blocker_sizes = Counter(value.bit_count() for value in complements)
    table = base.face_table(complements)
    vector = base.face_vector(table)
    recovered = base.minimal_nonfaces(table)
    ridges = base.ridge_histogram(complements)
    components = base.dual_component_count(complements)
    triples = all(table[value] for value in base.independent_masks(3))
    strict_failure = base.first_strict_local_g_failure(table)
    macaulay_failure = base.first_local_macaulay_failure(table)
    euler_failure = base.first_link_euler_failure(table)
    homology_failure = base.first_mod2_homology_failure(table)
    valid = (
        blocker_sizes == {5: 108}
        and vector == EXPECTED_F_VECTOR
        and set(recovered) == set(clutter)
        and ridges == {2: 486}
        and components == 1
        and triples
        and strict_failure is None
        and macaulay_failure is None
        and euler_failure is None
        and homology_failure is None
    )
    if not valid:
        raise AssertionError(
            "independent a1 topology failed: "
            f"blocker={blocker_sizes}, f={vector}, ridges={ridges}, "
            f"components={components}, strict={strict_failure}, "
            f"Macaulay={macaulay_failure}, Euler={euler_failure}, "
            f"homology={homology_failure}"
        )
    return (
        {
            "facet_complement_count": len(complements),
            "blocker_size_histogram": dict(sorted(blocker_sizes.items())),
            "face_vector": list(vector),
            "minimal_nonface_size_histogram": dict(
                sorted(Counter(value.bit_count() for value in recovered).items())
            ),
            "minimal_nonfaces_match_input_clutter": True,
            "ridge_degree_histogram": ridges,
            "dual_component_count": components,
            "covers_all_independent_triples": triples,
            "strict_local_g_constraints_hold": True,
            "local_macaulay_constraints_hold": True,
            "all_links_have_spherical_euler_characteristic": True,
            "all_links_are_mod2_homology_spheres": True,
            "passes_connected_homology_audit": True,
        },
        complements,
    )


def addition_mask(row: Sequence[int]) -> int:
    return base.mask(map(int, row))


def build_record(source: Path) -> dict[str, object]:
    document = json.loads(source.read_text(encoding="utf-8"))
    target_cases = [
        case
        for case in document["cases"]
        if case["representative_A4_masks"] == [A4_MASK]
    ]
    if len(target_cases) != 1:
        raise AssertionError("unique a1 target case is missing")
    target = target_cases[0]
    families = target["canonical_target_additions"]
    orbit_sizes = target["canonical_target_orbit_sizes"]
    if len(families) != len(orbit_sizes) or len(families) != 8:
        raise AssertionError("a1 target orbit archive changed")
    audits = []
    for additions_raw, orbit_size in zip(families, orbit_sizes, strict=True):
        additions = tuple(addition_mask(row) for row in additions_raw)
        topology, complements = topology_audit(additions)
        supports = (base.vertices(A4_MASK), *(base.vertices(row) for row in additions))
        labels = label_gate(supports, complements)
        if not labels["unrestricted_compact_coxeter_no_go"]:
            raise AssertionError("an a1 target class survived compact labels")
        audits.append(
            {
                "literal_orbit_size": orbit_size,
                "representative_addition_masks": list(additions),
                "representative_additions": [list(base.vertices(row)) for row in additions],
                "topology": topology,
                "coxeter_labels": labels,
            }
        )
    observed = sorted(
        (
            audit["literal_orbit_size"],
            tuple(
                stage["passing_state_count"]
                for stage in audit["coxeter_labels"]["stages"]
            ),
        )
        for audit in audits
    )
    expected = sorted(
        (
            (180, (108, 84, 8, 0)),
            (180, (108, 84, 24, 16, 72, 0)),
            (240, (108, 84, 24, 16, 0)),
            (240, (108, 84, 8, 0)),
            (40, (108, 0)),
            (120, (108, 84, 8, 0)),
            (60, (108, 84, 8, 0)),
            (20, (108, 0)),
        )
    )
    if observed != expected:
        raise AssertionError(f"a1 compact label transcripts changed: {observed}")
    return {
        "kind": "independent_h9_14_record18_a1_downstream_audit",
        "source": str(source),
        "raw_target_weight": sum(orbit_sizes),
        "canonical_class_count": len(audits),
        "class_orbit_sizes": sorted(orbit_sizes),
        "all_classes_pass_connected_mod2_homology_audit": True,
        "all_classes_have_unrestricted_compact_coxeter_no_go": True,
        "classes": audits,
        "a1_compact_branch_excluded": True,
        "compact_only": True,
        "finite_volume_noncompact_checked": False,
        "whole_shadow_compact_coxeter_theorem": False,
        "scope": (
            "Independent downstream audit of all eight a=1 target classes. "
            "No ideal-vertex or finite-volume noncompact conclusion."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, nargs="?", default=DEFAULT_INPUT)
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    record = build_record(args.input)
    rendered = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if args.write is not None:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
