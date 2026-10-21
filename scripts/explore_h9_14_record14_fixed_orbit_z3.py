#!/usr/bin/env python3
"""Exact fixed-positive-A4 incidence census for H9/14 Record 14.

For one archived positive-A4 orbit, this program enumerates every literal
family of additional minimal A5 rows satisfying the complete necessary
compact incidence model.  The Boolean model is lossless:

* the number of additions and all graph-derived local incidence intervals;
* destruction of every base transversal of order at most four;
* an iff variable for every surviving order-five transversal;
* the exact Dehn--Sommerville facet count; and
* degree zero or two for every six-set (the ridge law).

Every model is then blocker-dualized literally.  The program is an exact
incidence enumerator, not a polytopality or hyperbolic-volume theorem.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import sys
from typing import Iterable, Sequence

import z3

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    Pair,
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)


ROOT = Path(__file__).resolve().parents[1]
POSITIVE = ROOT / "data/enumerations/h9_14_record14_positive_A4_scout.json"
N = 14
FULL = (1 << N) - 1
EDGES = (
    (1 << 0) | (1 << 12),
    (1 << 1) | (1 << 12),
    (1 << 0) | (1 << 13),
    (1 << 1) | (1 << 13),
    (1 << 2) | (1 << 13),
)
BASE_G3 = (4, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 1)
BASE_G4 = (5, 5, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 2, 1)


def mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value & (1 << vertex))


def independent(value: int) -> bool:
    return all(value & edge != edge for edge in EDGES)


def hits_all(value: int, clutter: Sequence[int]) -> bool:
    return all(value & row for row in clutter)


def macaulay_raise(value: int, degree: int) -> int:
    remaining = value
    next_top: int | None = None
    raised = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (
            (next_top is None or top + 1 < next_top)
            and _comb(top + 1, lower) <= remaining
        ):
            top += 1
        remaining -= _comb(top, lower)
        raised += _comb(top + 1, lower + 1)
        next_top = top
    if remaining:
        raise AssertionError("inexact Macaulay expansion")
    return raised


def _comb(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    answer = 1
    for index in range(1, min(k, n - k) + 1):
        answer = answer * (n - min(k, n - k) + index) // index
    return answer


def minimal_transversals(clutter: Sequence[int]) -> tuple[int, ...]:
    answers = []
    for value in range(1, FULL + 1):
        if not hits_all(value, clutter):
            continue
        if all(
            not hits_all(value ^ (1 << vertex), clutter)
            for vertex in vertices(value)
        ):
            answers.append(value)
    return tuple(answers)


def rank_nine_vertices(minimal_nonfaces: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        chosen
        for chosen in combinations(range(N), 9)
        if all(mask(chosen) & nonface != nonface for nonface in minimal_nonfaces)
    )


@lru_cache(maxsize=None)
def placement_options(support: tuple[int, ...]) -> tuple[tuple[tuple[Pair, int], ...], ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True))
        for code in lanner_placement_codes(len(support))
    )


@lru_cache(maxsize=262_144)
def finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return is_finite_coxeter_labeling(range(order), labels)


def restriction_finite(restriction: Sequence[int], labels: dict[Pair, int]) -> bool:
    return finite_code(
        len(restriction),
        tuple(labels[pair(left, right)] for left, right in combinations(restriction, 2)),
    )


def audit_compact_labels(fours: Sequence[int], additions: Sequence[int]) -> dict[str, object]:
    supports = tuple(vertices(value) for value in (*fours, *additions))
    actual_vertices = rank_nine_vertices((*EDGES, *fours, *additions))
    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[tuple[tuple[Pair, int], ...], ...] = ((),)
    stages = []
    while states and remaining:
        if not union:
            support_index = min(remaining)
        else:
            support_index = max(
                remaining,
                key=lambda index: (
                    len(set(supports[index]) & union),
                    -len(set(supports[index]) - union),
                    -index,
                ),
            )
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
        passing: set[tuple[tuple[Pair, int], ...]] = set()
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
                if all(restriction_finite(restriction, labels) for restriction in restrictions):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "complete_lanner_placement_count": len(options),
                "restriction_size_histogram": dict(sorted(Counter(map(len, restrictions)).items())),
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
            }
        )
        remaining.remove(support_index)
    return {
        "rank_nine_actual_vertex_count": len(actual_vertices),
        "minimal_lanner_support_count": len(supports),
        "stages": stages,
        "terminal_sparse_label_state_count": len(states),
        "unrestricted_compact_label_no_go": not states,
        "monotonicity_scope": (
            "all unassigned meeting labels are 2; increasing any to 3, 4, or 5 "
            "increases the nonnegative cosine adjacency entrywise and cannot "
            "repair a non-elliptic actual-vertex restriction"
        ),
    }


def source_orbit(a: int, orbit_index: int) -> tuple[int, ...]:
    document = json.loads(POSITIVE.read_text(encoding="utf-8"))
    branch = next(row for row in document["branches"] if row["a"] == a)
    orbit = branch["orbits"][orbit_index]
    return tuple(mask(face) for face in orbit["representative"])


def build_problem(a: int, orbit_index: int) -> dict[str, object]:
    fours = source_orbit(a, orbit_index)
    if len(fours) != a or len(set(fours)) != a:
        raise AssertionError("positive-A4 source orbit changed")
    all_fives = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and independent(value)
    )
    shadow = tuple(
        five for five in all_fives if any(five & four == four for four in fours)
    )
    addition_count = 8 + 6 * a - len(shadow)
    scheduled_facets = 110 - 2 * a
    lower = []
    upper = []
    for vertex in range(N):
        q = sum(bool(four & (1 << vertex)) for four in fours)
        forced = sum(bool(five & (1 << vertex)) for five in shadow)
        previous = BASE_G3[vertex] - q
        before = BASE_G4[vertex] + 6 * q - forced
        lower.append(max(0, before - macaulay_raise(previous, 3)))
        upper.append(before - 1)
    if any(left > right for left, right in zip(lower, upper, strict=True)):
        raise AssertionError("archived source orbit has an empty addition interval")
    candidates = tuple(
        five
        for five in all_fives
        if five not in shadow
        and all(not (five & (1 << vertex)) or upper[vertex] >= 1 for vertex in range(N))
    )
    base = (*EDGES, *fours, *shadow)
    base_transversals = tuple(value for value in range(1, FULL + 1) if hits_all(value, base))
    small = tuple(value for value in base_transversals if value.bit_count() <= 4)
    fives = tuple(value for value in base_transversals if value.bit_count() == 5)
    return {
        "a": a,
        "orbit_index": orbit_index,
        "fours": fours,
        "shadow": shadow,
        "addition_count": addition_count,
        "scheduled_facets": scheduled_facets,
        "lower": tuple(lower),
        "upper": tuple(upper),
        "candidates": candidates,
        "base": base,
        "small": small,
        "fives": fives,
    }


def enumerate_models(problem: dict[str, object], max_models: int | None = None) -> dict[str, object]:
    candidates = tuple(problem["candidates"])
    small = tuple(problem["small"])
    fives = tuple(problem["fives"])
    addition_count = int(problem["addition_count"])
    scheduled_facets = int(problem["scheduled_facets"])
    lower = tuple(problem["lower"])
    upper = tuple(problem["upper"])
    context = z3.Context()
    selected = [z3.Bool(f"x_{index}", ctx=context) for index in range(len(candidates))]
    survives = [z3.Bool(f"y_{index}", ctx=context) for index in range(len(fives))]
    solver = z3.Solver(ctx=context)
    solver.add(z3.PbEq([(variable, 1) for variable in selected], addition_count))
    for vertex in range(N):
        terms = [
            (selected[index], 1)
            for index, face in enumerate(candidates)
            if face & (1 << vertex)
        ]
        if terms:
            solver.add(z3.PbGe(terms, lower[vertex]))
            solver.add(z3.PbLe(terms, upper[vertex]))
        else:
            solver.add(z3.BoolVal(lower[vertex] <= 0 <= upper[vertex], ctx=context))
    for transversal in small:
        destroyers = [
            selected[index]
            for index, face in enumerate(candidates)
            if not face & transversal
        ]
        solver.add(z3.Or(destroyers) if destroyers else z3.BoolVal(False, ctx=context))
    for index, transversal in enumerate(fives):
        destroyers = [
            selected[candidate]
            for candidate, face in enumerate(candidates)
            if not face & transversal
        ]
        solver.add(
            survives[index]
            == (z3.And([z3.Not(value) for value in destroyers]) if destroyers else z3.BoolVal(True, ctx=context))
        )
    solver.add(z3.PbEq([(variable, 1) for variable in survives], scheduled_facets))
    six_constraint_count = 0
    for six in range(FULL + 1):
        if six.bit_count() != 6:
            continue
        members = [
            (survives[index], 1)
            for index, facet in enumerate(fives)
            if facet & six == facet
        ]
        if not members:
            continue
        six_constraint_count += 1
        solver.add(z3.Or(z3.PbEq(members, 0), z3.PbEq(members, 2)))

    model_count = 0
    pure_targets: set[tuple[int, ...]] = set()
    blocker_histogram: Counter[str] = Counter()
    complete = True
    while solver.check() == z3.sat:
        model = solver.model()
        indices = tuple(
            index for index, variable in enumerate(selected) if z3.is_true(model.eval(variable))
        )
        if len(indices) != addition_count:
            raise AssertionError("Z3 model has the wrong addition cardinality")
        family = tuple(sorted(candidates[index] for index in indices))
        clutter = (*problem["base"], *family)
        blocker = minimal_transversals(clutter)
        sizes = Counter(value.bit_count() for value in blocker)
        key = ",".join(f"{size}:{count}" for size, count in sorted(sizes.items()))
        blocker_histogram[key] += 1
        if sizes == {5: scheduled_facets}:
            pure_targets.add(family)
        model_count += 1
        solver.add(z3.Or([z3.Not(selected[index]) for index in indices]))
        if model_count % 1000 == 0:
            print(f"models {model_count} pure {len(pure_targets)}", file=sys.stderr, flush=True)
        if max_models is not None and model_count >= max_models:
            complete = False
            break
    terminal = "cutoff" if not complete else "unsat_after_blocking_all_models"
    return {
        "six_set_constraint_count": six_constraint_count,
        "z3_model_count": model_count,
        "blocker_size_histogram": dict(sorted(blocker_histogram.items())),
        "literal_pure_target_count": len(pure_targets),
        "literal_pure_targets": [
            [list(vertices(face)) for face in family] for family in sorted(pure_targets)
        ],
        "enumeration_complete": complete,
        "terminal_solver_status": terminal,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("a", type=int)
    parser.add_argument("orbit_index", type=int)
    parser.add_argument("--max-models", type=int)
    parser.add_argument("--labels", action="store_true")
    parser.add_argument("--label-target", type=int)
    args = parser.parse_args(argv)
    if args.a not in (1, 2, 3, 4):
        raise SystemExit("a must lie in 1..4")
    problem = build_problem(args.a, args.orbit_index)
    census = enumerate_models(problem, args.max_models)
    if args.labels and census["enumeration_complete"]:
        target_label_audits = []
        indexed_targets = list(enumerate(census["literal_pure_targets"]))
        if args.label_target is not None:
            indexed_targets = [
                row for row in indexed_targets if row[0] == args.label_target
            ]
            if not indexed_targets:
                raise SystemExit("label target index is outside the literal target list")
        for target_index, family in indexed_targets:
            additions = tuple(mask(face) for face in family)
            label_audit = audit_compact_labels(problem["fours"], additions)
            target_label_audits.append(
                {"target_index": target_index, **label_audit}
            )
            print(
                f"target {target_index} terminal labels "
                f"{label_audit['terminal_sparse_label_state_count']}",
                file=sys.stderr,
                flush=True,
            )
    else:
        target_label_audits = []
    output = {
        "kind": "h9_14_record14_fixed_positive_A4_z3_incidence_census",
        "a": args.a,
        "source_positive_A4_orbit_index": args.orbit_index,
        "representative_A4s": [list(vertices(value)) for value in problem["fours"]],
        "forced_A5_shadow_count": len(problem["shadow"]),
        "additional_A5_count": problem["addition_count"],
        "scheduled_facet_count": problem["scheduled_facets"],
        "addition_lower": list(problem["lower"]),
        "addition_upper": list(problem["upper"]),
        "candidate_count": len(problem["candidates"]),
        "base_small_transversal_count": len(problem["small"]),
        "base_order_five_transversal_count": len(problem["fives"]),
        **census,
        "target_label_audits": target_label_audits,
        "all_literal_targets_compact_label_excluded": bool(
            census["enumeration_complete"]
            and len(target_label_audits) == census["literal_pure_target_count"]
            and target_label_audits
            and all(row["unrestricted_compact_label_no_go"] for row in target_label_audits)
        ),
        "scope": (
            "Exact necessary compact incidence/ridge census. A complete zero-target "
            "row closes this fixed source orbit before topology or Coxeter labels; "
            "no finite-volume noncompact inference."
        ),
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
