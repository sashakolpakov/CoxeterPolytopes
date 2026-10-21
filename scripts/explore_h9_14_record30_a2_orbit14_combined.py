#!/usr/bin/env python3
"""Conflict-driven compact scout for Record 30, ``a=2`` orbit 14.

The older literal searches choose five additional missing ``A5`` supports
from 966 candidates before applying the Coxeter-label gate.  Their raw search
space has more than 6.9 trillion five-subsets.  This scout instead combines
the exact incidence equations with the complete sparse Lanner join:

* Z3 enumerates one symmetry orbit of scheduled incidence targets at a time;
* the full blocker is replayed literally on fourteen vertices;
* every order-four/five Lanner placement is joined on the minimal nonfaces;
* every actual rank-nine vertex restriction is required to be elliptic; and
* an entire stabilizer orbit of a rejected incidence target is blocked.

The sparse join is a sound compact no-go: all otherwise unassigned meeting
pairs are set to label two.  Replacing any of them by 3, 4, or 5 only increases
the nonnegative cosine adjacency matrix entrywise and cannot repair a
non-elliptic principal vertex diagram.  A surviving sparse state is only an
open candidate, because its remaining labels and global Gram equations have
not been solved.

This is an exploratory enumerator, not yet a frozen solver-independent
certificate.  It makes no finite-volume noncompact inference.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
import argparse
import json
from typing import Iterable, Mapping, Sequence

from z3 import And, Bool, If, Implies, Not, Or, Solver, Sum, is_true, sat, unsat

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)


N = 14
DIMENSION = 9
SCHEDULED_FIVES = 20
SCHEDULED_FACETS = 128
FIXED_MISSING_EDGES = ((0, 11), (1, 12), (2, 13))
FIXED_MISSING_FOURS = ((0, 1, 3, 4), (0, 2, 3, 4))
BASE_LINK_G3 = (4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 4, 4, 4)
BASE_LINK_G4 = BASE_LINK_G3

Mask = int
LabelState = tuple[tuple[Pair, int], ...]


def mask(vertices: Iterable[int]) -> Mask:
    value = 0
    for vertex in vertices:
        value |= 1 << vertex
    return value


def vertices(value: Mask) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value >> vertex & 1)


MISSING_EDGE_MASKS = tuple(mask(edge) for edge in FIXED_MISSING_EDGES)
MISSING_FOUR_MASKS = tuple(mask(face) for face in FIXED_MISSING_FOURS)


def _binomial(n: int, k: int) -> int:
    if not 0 <= k <= n:
        return 0
    import math

    return math.comb(n, k)


def _macaulay_raise(value: int, degree: int) -> int:
    remaining = value
    next_top = -1
    answer = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (
            (next_top < 0 or top + 1 < next_top)
            and _binomial(top + 1, lower) <= remaining
        ):
            top += 1
        remaining -= _binomial(top, lower)
        answer += _binomial(top + 1, lower + 1)
        next_top = top
    if remaining:
        raise AssertionError("Macaulay expansion changed")
    return answer


def _hits_all(transversal: Mask, clutter: Sequence[Mask]) -> bool:
    return all(transversal & face for face in clutter)


def _independent_fives() -> tuple[Mask, ...]:
    answer = tuple(
        mask(face)
        for face in combinations(range(N), 5)
        if all(mask(face) & edge != edge for edge in MISSING_EDGE_MASKS)
    )
    if len(answer) != 1372:
        raise AssertionError("independent-five ground changed")
    return answer


def _five_vertex_covers() -> tuple[Mask, ...]:
    return tuple(
        mask(face)
        for face in combinations(range(N), 5)
        if all(mask(face) & edge for edge in MISSING_EDGE_MASKS)
    )


def _minimal_transversals(clutter: Sequence[Mask]) -> tuple[Mask, ...]:
    answer = []
    for candidate in range(1, 1 << N):
        if not _hits_all(candidate, clutter):
            continue
        if all(
            not _hits_all(candidate ^ (1 << vertex), clutter)
            for vertex in vertices(candidate)
        ):
            answer.append(candidate)
    return tuple(answer)


@lru_cache(maxsize=None)
def _placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    ordered = tuple(sorted(support))
    edges = tuple(combinations(ordered, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(ordered))
    )


@lru_cache(maxsize=1_048_576)
def _finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return is_finite_coxeter_labeling(range(order), labels)


def _restriction_finite(
    restriction: Sequence[int], labels: Mapping[Pair, int]
) -> bool:
    return _finite_code(
        len(restriction),
        tuple(
            labels[pair(left, right)]
            for left, right in combinations(restriction, 2)
        ),
    )


def _next_support(
    supports: Sequence[tuple[int, ...]], remaining: set[int], union: set[int]
) -> int:
    if not union:
        return 0
    return max(
        remaining,
        key=lambda index: (
            len(set(supports[index]) & union),
            len(supports[index]),
            -len(set(supports[index]) - union),
            -index,
        ),
    )


def _sparse_label_gate(additions: Sequence[Mask]) -> dict[str, object]:
    supports = tuple(FIXED_MISSING_FOURS) + tuple(vertices(face) for face in additions)
    forbidden = tuple(map(frozenset, (*FIXED_MISSING_EDGES, *supports)))
    actual_vertices = tuple(
        chosen
        for chosen in combinations(range(N), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )
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
                    tuple(vertex for vertex in actual if vertex in union)
                    for actual in actual_vertices
                },
                key=lambda row: (len(row), row),
            )
        )
        options = _placement_options(tuple(sorted(support)))
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
                    pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(
                    _restriction_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "raw_extension_count": len(options)
                * (stages[-1]["passing_state_count"] if stages else 1),
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
                "restriction_size_histogram": dict(
                    sorted(Counter(map(len, restrictions)).items())
                ),
            }
        )
        remaining.remove(support_index)
    return {
        "rank_nine_vertex_count": len(actual_vertices),
        "stages": stages,
        "terminal_sparse_state_count": len(states),
        "compact_label_no_go": not states,
    }


def _stabilizer() -> tuple[tuple[int, ...], ...]:
    """The full 2 x 2 x S6 stabilizer of the orbit-14 four-support pair."""

    answer = []
    for swap_12 in (False, True):
        for swap_34 in (False, True):
            for isolate_tail in permutations(range(5, 11)):
                permutation = list(range(N))
                if swap_12:
                    permutation[1], permutation[2] = 2, 1
                    permutation[12], permutation[13] = 13, 12
                if swap_34:
                    permutation[3], permutation[4] = 4, 3
                for source, target in zip(range(5, 11), isolate_tail, strict=True):
                    permutation[source] = target
                answer.append(tuple(permutation))
    if len(answer) != 2880 or len(set(answer)) != 2880:
        raise AssertionError("orbit-14 stabilizer changed")
    return tuple(answer)


def _stabilizer_generators() -> tuple[tuple[int, ...], ...]:
    generators = []
    for swaps in (
        ((1, 2), (12, 13)),
        ((3, 4),),
        ((5, 6),),
        ((6, 7),),
        ((7, 8),),
        ((8, 9),),
        ((9, 10),),
    ):
        permutation = list(range(N))
        for left, right in swaps:
            permutation[left], permutation[right] = right, left
        generators.append(tuple(permutation))
    return tuple(generators)


def _permuted_mask(value: Mask, permutation: Sequence[int]) -> Mask:
    return mask(permutation[vertex] for vertex in vertices(value))


def _family_orbit(
    family: Sequence[Mask], stabilizer: Sequence[Sequence[int]]
) -> tuple[tuple[Mask, ...], ...]:
    return tuple(
        sorted(
            {
                tuple(sorted(_permuted_mask(face, permutation) for face in family))
                for permutation in stabilizer
            }
        )
    )


def _build_incidence_solver() -> tuple[
    Solver,
    tuple[object, ...],
    tuple[Mask, ...],
    dict[str, object],
]:
    independent_fives = _independent_fives()
    shadow = tuple(
        face
        for face in independent_fives
        if any(face & four == four for four in MISSING_FOUR_MASKS)
    )
    additions = SCHEDULED_FIVES - len(shadow)
    if len(shadow) != 15 or additions != 5:
        raise AssertionError("orbit-14 addition schedule changed")

    upper = []
    for vertex in range(N):
        four_incidence = sum(four >> vertex & 1 for four in MISSING_FOUR_MASKS)
        shadow_incidence = sum(face >> vertex & 1 for face in shadow)
        previous = BASE_LINK_G3[vertex] - four_incidence
        before = (
            BASE_LINK_G4[vertex]
            + 6 * four_incidence
            - shadow_incidence
        )
        lower = max(0, before - _macaulay_raise(previous, 3))
        if lower:
            raise AssertionError("orbit-14 lower incidence cap changed")
        upper.append(before - 1)

    candidates = tuple(
        face
        for face in independent_fives
        if not any(face & four == four for four in MISSING_FOUR_MASKS)
        and all(not (face >> vertex & 1) or upper[vertex] >= 1 for vertex in range(N))
    )
    if len(candidates) != 966:
        raise AssertionError("orbit-14 candidate count changed")

    base = (*MISSING_EDGE_MASKS, *MISSING_FOUR_MASKS, *shadow)
    # Cardinality-first ordering is deterministic and materially improves the
    # SAT branching: the strongest small-cover clauses enter first.
    small = tuple(
        mask(support)
        for order in range(1, 5)
        for support in combinations(range(N), order)
        if _hits_all(mask(support), base)
    )
    if len(small) != 61:
        raise AssertionError("orbit-14 small-transversal count changed")

    covers = tuple(
        cover
        for cover in _five_vertex_covers()
        if all(cover & four for four in MISSING_FOUR_MASKS)
    )
    if len(covers) != 269:
        raise AssertionError("orbit-14 initial five-cover count changed")

    solver = Solver()
    selected = tuple(Bool(f"x_{index}") for index in range(len(candidates)))
    surviving = tuple(Bool(f"y_{index}") for index in range(len(covers)))
    solver.add(Sum([If(value, 1, 0) for value in selected]) == additions)
    for vertex, cap in enumerate(upper):
        solver.add(
            Sum(
                [
                    If(selected[index], 1, 0)
                    for index, face in enumerate(candidates)
                    if face >> vertex & 1
                ]
            )
            <= cap
        )
    for transversal in small:
        solver.add(
            Or(
                [
                    selected[index]
                    for index, face in enumerate(candidates)
                    if not face & transversal
                ]
            )
        )
    for cover_index, cover in enumerate(covers):
        destroyers = [
            selected[index]
            for index, face in enumerate(candidates)
            if not face & cover
        ]
        solver.add(surviving[cover_index] == Not(Or(destroyers)))
    solver.add(Sum([If(value, 1, 0) for value in surviving]) == SCHEDULED_FACETS)

    # Every stabilizer orbit has a lexicographically least characteristic
    # vector.  Requiring the current vector to precede each generator image is
    # therefore lossless and sharply reduces repeated symmetric branching.
    candidate_index = {face: index for index, face in enumerate(candidates)}
    for permutation in _stabilizer_generators():
        image_index = tuple(
            candidate_index[_permuted_mask(face, permutation)]
            for face in candidates
        )
        prefix_equal = True
        for index, image in enumerate(image_index):
            solver.add(
                Implies(
                    prefix_equal,
                    Or(Not(selected[index]), selected[image]),
                )
            )
            prefix_equal = And(
                prefix_equal,
                selected[index] == selected[image],
            )
    metadata = {
        "forced_shadow_count": len(shadow),
        "additional_A5_count": additions,
        "candidate_count": len(candidates),
        "raw_unordered_addition_count": _binomial(len(candidates), additions),
        "base_small_transversal_count": len(small),
        "initial_order_five_cover_count": len(covers),
        "local_upper_caps": upper,
        "lossless_lex_generator_count": len(_stabilizer_generators()),
    }
    return solver, selected, candidates, metadata


def run(max_orbits: int | None, timeout_ms: int | None) -> dict[str, object]:
    solver, selected, candidates, metadata = _build_incidence_solver()
    if timeout_ms is not None:
        solver.set(timeout=timeout_ms)
    candidate_index = {face: index for index, face in enumerate(candidates)}
    stabilizer = _stabilizer()
    rows = []
    blocked_families: set[tuple[Mask, ...]] = set()
    pure_orbit_count = 0
    label_no_go_orbit_count = 0
    while max_orbits is None or len(rows) < max_orbits:
        status = solver.check()
        if status == unsat:
            conclusion = "exhausted"
            break
        if status != sat:
            conclusion = "solver_unknown"
            break
        model = solver.model()
        family = tuple(
            candidates[index]
            for index, variable in enumerate(selected)
            if is_true(model.eval(variable))
        )
        if len(family) != 5:
            raise AssertionError("model did not choose five additions")
        clutter = (*MISSING_EDGE_MASKS, *MISSING_FOUR_MASKS, *family)
        blocker = _minimal_transversals(clutter)
        blocker_histogram = dict(sorted(Counter(face.bit_count() for face in blocker).items()))
        pure = len(blocker) == SCHEDULED_FACETS and blocker_histogram == {5: 128}
        label_gate = _sparse_label_gate(family) if pure else None
        orbit = _family_orbit(family, stabilizer)
        new_orbit = tuple(item for item in orbit if item not in blocked_families)
        if not new_orbit:
            raise AssertionError("a blocked family returned from the solver")
        for transformed in new_orbit:
            solver.add(
                Or(
                    [
                        Not(selected[candidate_index[face]])
                        for face in transformed
                    ]
                )
            )
            blocked_families.add(transformed)
        if pure:
            pure_orbit_count += 1
            if label_gate is not None and label_gate["compact_label_no_go"]:
                label_no_go_orbit_count += 1
        row = {
            "orbit_index": len(rows),
            "representative_additional_A5s": [list(vertices(face)) for face in sorted(family)],
            "stabilizer_orbit_size": len(orbit),
            "minimal_transversal_size_histogram": blocker_histogram,
            "pure_128_facet_target": pure,
            "sparse_label_gate": label_gate,
        }
        rows.append(row)
        print(json.dumps(row, sort_keys=True), flush=True)
        if pure and label_gate is not None and not label_gate["compact_label_no_go"]:
            conclusion = "sparse_label_survivor"
            break
    else:
        conclusion = "orbit_limit"
    return {
        "kind": "h9_14_record30_a2_orbit14_combined_scout",
        "compact_only": True,
        "incidence_model": metadata,
        "stabilizer_order": len(stabilizer),
        "enumerated_incidence_orbit_count": len(rows),
        "enumerated_pure_target_orbit_count": pure_orbit_count,
        "label_no_go_pure_target_orbit_count": label_no_go_orbit_count,
        "blocked_literal_family_count": len(blocked_families),
        "conclusion": conclusion,
        "whole_orbit14_compact_coxeter_excluded": (
            conclusion == "exhausted"
            and pure_orbit_count == label_no_go_orbit_count
        ),
        "scope": (
            "Conflict-driven Z3 incidence scout plus exact sparse Lanner and "
            "actual-vertex gate for Record-30 a=2 orbit 14. A theorem claim "
            "still requires a frozen solver-independent replay. No "
            "finite-volume inference."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-orbits", type=int)
    parser.add_argument("--timeout-ms", type=int)
    args = parser.parse_args()
    print(
        json.dumps(
            run(max_orbits=args.max_orbits, timeout_ms=args.timeout_ms),
            sort_keys=True,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
