#!/usr/bin/env python3
"""Exact symmetry quotient and compact-label gate for Record 14 ``a=0``.

The incidence source is supplied by the independent CP-SAT full-source
generator-lex census.  This script quotients its finite target list by the
full ``S9 x S2`` missing-graph action and applies the complete sparse
Lanner/actual-vertex compact-label gate to one representative per orbit.
It is exploratory until the Z3 literal target list is compared and a frozen
certificate is issued.
"""

from __future__ import annotations

import argparse
from itertools import permutations, product
import json
from typing import Iterable, Sequence

from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    audit_compact_labels,
    mask,
    vertices,
)
from scripts.independent_h9_14_record14_a0_audit import replay_a0


N = 14
ISOLATES = tuple(range(3, 12))
FIXED = (0, 1, 2, 12, 13)


def _move_active(value: int, swap: bool) -> int:
    if not swap:
        return value
    zero = bool(value & 1)
    one = bool(value & 2)
    value &= ~3
    if zero:
        value |= 2
    if one:
        value |= 1
    return value


def _fixed_code(value: int) -> int:
    return sum(bool(value & (1 << vertex)) << index for index, vertex in enumerate(FIXED))


def _block_orders(blocks: Sequence[Sequence[int]]) -> Iterable[tuple[int, ...]]:
    for choices in product(*(permutations(block) for block in blocks)):
        yield tuple(index for choice in choices for index in choice)


def full_source_orbit_key(family: Sequence[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Canonical unordered incidence-matrix key under ``S9 x S2``.

    After either active twin action, row permutations are needed only within
    equal fixed-column codes.  For each such order, sorting the nine isolate
    column words is the exact quotient by ``S9``.
    """

    best: tuple[tuple[int, ...], tuple[int, ...]] | None = None
    for swap in (False, True):
        moved = tuple(_move_active(value, swap) for value in family)
        grouped: dict[int, list[int]] = {}
        for index, value in enumerate(moved):
            grouped.setdefault(_fixed_code(value), []).append(index)
        codes = tuple(sorted(grouped))
        blocks = tuple(grouped[code] for code in codes)
        for order in _block_orders(blocks):
            fixed_rows = tuple(_fixed_code(moved[index]) for index in order)
            isolate_columns = tuple(
                sorted(
                    sum(
                        bool(moved[index] & (1 << isolate)) << position
                        for position, index in enumerate(order)
                    )
                    for isolate in ISOLATES
                )
            )
            key = (fixed_rows, isolate_columns)
            if best is None or key < best:
                best = key
    if best is None:
        raise AssertionError("empty Record-14 a=0 target family")
    return best


def downstream() -> dict[str, object]:
    incidence = replay_a0()
    targets = tuple(
        tuple(sorted(mask(face) for face in family))
        for family in incidence["generator_lex_slice_pure_targets"]
    )
    groups: dict[tuple[tuple[int, ...], tuple[int, ...]], list[int]] = {}
    for index, family in enumerate(targets):
        groups.setdefault(full_source_orbit_key(family), []).append(index)
    ordered_groups = tuple(tuple(indices) for _key, indices in sorted(groups.items()))
    audits = []
    for orbit_index, indices in enumerate(ordered_groups):
        representative = targets[indices[0]]
        audit = audit_compact_labels((), representative)
        audits.append(
            {
                "target_orbit_index": orbit_index,
                "generator_lex_target_indices_covered": list(indices),
                "representative_additional_A5s": [
                    list(vertices(face)) for face in representative
                ],
                **audit,
            }
        )
    return {
        "kind": "h9_14_record14_a0_downstream_exploratory",
        "incidence_backend": incidence["kind"],
        "incidence_solver_status": incidence["solver_status"],
        "incidence_enumeration_complete": incidence["enumeration_complete"],
        "generator_lex_slice_pure_target_count": len(targets),
        "full_source_symmetry_order": 2 * 362880,
        "full_target_orbit_count": len(ordered_groups),
        "target_orbit_groups": [list(group) for group in ordered_groups],
        "compact_label_audits": audits,
        "all_target_orbits_unrestricted_compact_label_excluded": all(
            audit["unrestricted_compact_label_no_go"] for audit in audits
        ),
        "terminal_sparse_label_state_count": sum(
            audit["terminal_sparse_label_state_count"] for audit in audits
        ),
        "scope": (
            "Complete CP-SAT incidence slice and exact full-source symmetry "
            "quotient. Compact-only label gate; no noncompact finite-volume "
            "inference. A frozen theorem additionally requires literal Z3 "
            "target equality."
        ),
    }


def main(argv: Sequence[str] | None = None) -> None:
    argparse.ArgumentParser(description=__doc__).parse_args(argv)
    print(json.dumps(downstream(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
