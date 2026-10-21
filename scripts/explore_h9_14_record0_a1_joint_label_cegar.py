#!/usr/bin/env python3
"""Complete joint incidence/compact-label search for Record 0, ``a=1``.

Unlike the target-by-target local-core queue, this model existentially chooses
the dual nerve and one global Coxeter labelling at the same time.  The fixed
``A4`` is constrained to one of all 108 labelled order-four Lanner
placements; every genuinely minimal missing independent five-set is
conditionally constrained to one of all 420 order-five placements.

Ellipticity of faces is completed by exact cutting planes.  Whenever a
selected rank-nine vertex is non-elliptic, an inclusion-minimal non-elliptic
subset is extracted.  Perron--Frobenius monotonicity makes the cut global:
that labelled subset (or any componentwise larger labels) cannot be a face of
any compact Coxeter polytope.  Face variables are defined iff a selected
rank-nine vertex extends the subset.

There are finitely many incidence and label assignments.  Hence unrestricted
execution has a sharp terminal outcome: ``INFEASIBLE`` closes the fixed
compact source, while a model whose every selected rank-nine vertex is
elliptic is an exact-label survivor for the Gram stage.  Timeout or an
explicit round cutoff is reported as inconclusive.  Finite-volume noncompact
polytopes are not addressed.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable, Sequence

from ortools.sat.python import cp_model

from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)
from scripts.explore_h9_14_record0_a1_dual_facets import (
    build_dual_model,
    contains_base,
    source_generators,
)
from scripts.explore_h9_14_record14_a1_symmetry_z3 import _move, _swap_action
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    minimal_transversals,
    vertices,
)


N = 14
FULL = (1 << N) - 1
LabelCode = tuple[int, ...]
AbstractPattern = tuple[int, LabelCode]


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _minimal_nonfinite(
    actual: Sequence[int], labels: dict[tuple[int, int], int]
) -> tuple[int, ...]:
    subset = list(actual)
    index = 0
    while index < len(subset):
        trial = (*subset[:index], *subset[index + 1 :])
        if trial and not is_finite_coxeter_labeling(trial, labels):
            subset.pop(index)
        else:
            index += 1
    answer = tuple(subset)
    if is_finite_coxeter_labeling(answer, labels):
        raise AssertionError("minimal nonfinite extraction became finite")
    if any(
        not is_finite_coxeter_labeling(
            (*answer[:index], *answer[index + 1 :]), labels
        )
        for index in range(len(answer))
    ):
        raise AssertionError("nonfinite subset is not deletion-minimal")
    return answer


def _permuted_code(code: LabelCode, permutation: Sequence[int]) -> LabelCode:
    order = len(permutation)
    edges = tuple(combinations(range(order), 2))
    source = dict(zip(edges, code, strict=True))
    target = {
        pair(permutation[left], permutation[right]): label
        for (left, right), label in source.items()
    }
    return tuple(target[edge] for edge in edges)


@lru_cache(maxsize=None)
def _code_orbit(code: LabelCode) -> tuple[LabelCode, ...]:
    order = 1
    while order * (order - 1) // 2 < len(code):
        order += 1
    if order * (order - 1) // 2 != len(code):
        raise ValueError("label code has no integral order")
    return tuple(
        sorted(
            {
                _permuted_code(code, permutation)
                for permutation in permutations(range(order))
            }
        )
    )


def _canonical_pattern(
    subset: Sequence[int], labels: dict[tuple[int, int], int]
) -> AbstractPattern:
    subset = tuple(sorted(subset))
    code = tuple(labels[pair(left, right)] for left, right in combinations(subset, 2))
    return len(subset), min(_code_orbit(code))


def _pattern_from_cut(cut: dict[str, object]) -> AbstractPattern:
    subset = tuple(int(vertex) for vertex in cut["restriction"])
    labels = {
        pair(left, right): 2 for left, right in combinations(subset, 2)
    }
    for left, right, value in cut["positive_labels"]:
        labels[pair(int(left), int(right))] = int(value)
    return _canonical_pattern(subset, labels)


def _pattern_document(pattern: AbstractPattern) -> dict[str, object]:
    order, code = pattern
    return {
        "order": order,
        "canonical_code": list(code),
        "positive_labels": [
            [left, right, value]
            for (left, right), value in zip(
                combinations(range(order), 2), code, strict=True
            )
            if value > 2
        ],
        "labelled_embedding_orbit_size": len(_code_orbit(code)),
    }


def _pattern_from_document(document: dict[str, object]) -> AbstractPattern:
    pattern = (int(document["order"]), tuple(document["canonical_code"]))
    if pattern[1] != min(_code_orbit(pattern[1])):
        raise ValueError("stored nonfinite pattern is not canonical")
    return pattern


def _move_code_on_support(
    code: LabelCode, support: Sequence[int], action: Sequence[int]
) -> LabelCode:
    support = tuple(sorted(support))
    source = dict(zip(combinations(support, 2), code, strict=True))
    target = {
        pair(action[left], action[right]): value
        for (left, right), value in source.items()
    }
    return tuple(target[edge] for edge in combinations(support, 2))


def a4_label_orbits(source_orbit: int) -> tuple[tuple[LabelCode, ...], ...]:
    """Exact placement orbits under the source action induced on the A4."""

    probe = build_dual_model(source_orbit, symmetry_breaking=False)
    four = tuple(probe["problem"]["fours"])[0]
    support = vertices(four)
    inside = tuple(vertex for vertex in range(2, 13) if four & (1 << vertex))
    active_swaps = (
        (False, True)
        if bool(four & 1) == bool(four & 2)
        else (False,)
    )
    actions = []
    for image in permutations(inside):
        base = list(range(N))
        for source, target in zip(inside, image, strict=True):
            base[source] = target
        for active_swap in active_swaps:
            action = base.copy()
            if active_swap:
                action[0], action[1] = 1, 0
            actions.append(tuple(action))
    remaining = set(lanner_placement_codes(4))
    answer = []
    while remaining:
        representative = min(remaining)
        orbit = tuple(
            sorted(
                {
                    _move_code_on_support(representative, support, action)
                    for action in actions
                }
            )
        )
        if any(code not in remaining for code in orbit):
            raise AssertionError("A4 placement action left its partition block")
        remaining.difference_update(orbit)
        answer.append(orbit)
    return tuple(answer)


def _add_general_lex_leader(
    model: cp_model.CpModel,
    chosen: Sequence[cp_model.IntVar],
    permutation: Sequence[int],
    tag: str,
) -> int:
    if sorted(permutation) != list(range(len(chosen))):
        raise AssertionError("facet action is not a permutation")
    prefix = model.new_bool_var(f"{tag}_prefix_0")
    model.add(prefix == 1)
    steps = 0
    for position, partner in enumerate(permutation):
        if position == partner:
            continue
        model.add(chosen[position] <= chosen[partner]).only_enforce_if(prefix)
        equal = model.new_bool_var(f"{tag}_equal_{steps}")
        model.add(chosen[position] == chosen[partner]).only_enforce_if(equal)
        model.add(chosen[position] != chosen[partner]).only_enforce_if(equal.Not())
        following = model.new_bool_var(f"{tag}_prefix_{steps + 1}")
        model.add(following <= prefix)
        model.add(following <= equal)
        model.add(following >= prefix + equal - 1)
        prefix = following
        steps += 1
    return steps


def _branch_stabilizer_actions(four: int, code: LabelCode) -> tuple[tuple[int, ...], ...]:
    """Small generating superset for the exact labelled-A4 stabilizer."""

    isolates = tuple(range(2, 13))
    inside = tuple(vertex for vertex in isolates if four & (1 << vertex))
    outside = tuple(vertex for vertex in isolates if not four & (1 << vertex))
    actions = [
        _swap_action(outside[index], outside[index + 1])
        for index in range(len(outside) - 1)
    ]
    active_swaps = (
        (False, True)
        if bool(four & 1) == bool(four & 2)
        else (False,)
    )
    support = vertices(four)
    for image in permutations(inside):
        base = list(range(N))
        for source, target in zip(inside, image, strict=True):
            base[source] = target
        for active_swap in active_swaps:
            action = base.copy()
            if active_swap:
                action[0], action[1] = 1, 0
            action_tuple = tuple(action)
            if action_tuple == tuple(range(N)):
                continue
            if _move_code_on_support(code, support, action_tuple) == code:
                actions.append(action_tuple)
    return tuple(dict.fromkeys(actions))


def run_joint_cegar(
    source_orbit: int,
    maximum_rounds: int | None,
    seconds_per_solve: float | None,
    *,
    a4_orbit_index: int | None = None,
    seeded_rounds: Sequence[dict[str, object]] = (),
    seed_rounds_are_history: bool = False,
    maximum_seed_pattern_order: int | None = None,
    checkpoint_path: Path | None = None,
) -> dict[str, object]:
    data = build_dual_model(
        source_orbit, symmetry_breaking=a4_orbit_index is None
    )
    problem = data["problem"]
    model = data["model"]
    live = tuple(data["live"])
    facets = tuple(data["facets"])
    missing_five = tuple(data["missing_five"])
    five_ground = tuple(data["five_ground"])
    fours = tuple(problem["fours"])
    four = fours[0]
    a4_orbits = a4_label_orbits(source_orbit)
    branch_code: LabelCode | None = None
    branch_lex_counts: tuple[int, ...] = ()
    if a4_orbit_index is not None:
        if not 0 <= a4_orbit_index < len(a4_orbits):
            raise ValueError(
                f"A4 orbit index must lie in 0..{len(a4_orbits) - 1}"
            )
        branch_code = a4_orbits[a4_orbit_index][0]
        facet_index = {facet: index for index, facet in enumerate(facets)}
        branch_counts = []
        for action_index, action in enumerate(
            _branch_stabilizer_actions(four, branch_code)
        ):
            permutation = tuple(
                facet_index[_move(facet, action)] for facet in facets
            )
            branch_counts.append(
                _add_general_lex_leader(
                    model,
                    live,
                    permutation,
                    f"branch_lex_{action_index}",
                )
            )
        branch_lex_counts = tuple(branch_counts)

    meeting_pairs = tuple(
        pair(left, right)
        for left, right in combinations(range(N), 2)
        if mask((left, right)) not in problem["edge_rows"]
    )
    label_variables = {
        edge: model.new_int_var(2, 5, f"label_{edge[0]}_{edge[1]}")
        for edge in meeting_pairs
    }
    nonright: dict[tuple[int, int], cp_model.IntVar] = {}
    label_four: dict[tuple[int, int], cp_model.IntVar] = {}
    label_five: dict[tuple[int, int], cp_model.IntVar] = {}
    for edge in meeting_pairs:
        nonright[edge] = model.new_bool_var(
            f"label_{edge[0]}_{edge[1]}_nonright"
        )
        label_four[edge] = model.new_bool_var(
            f"label_{edge[0]}_{edge[1]}_is4"
        )
        label_five[edge] = model.new_bool_var(
            f"label_{edge[0]}_{edge[1]}_is5"
        )
        model.add_allowed_assignments(
            [
                label_variables[edge],
                nonright[edge],
                label_four[edge],
                label_five[edge],
            ],
            (
                (2, 0, 0, 0),
                (3, 1, 0, 0),
                (4, 1, 1, 0),
                (5, 1, 0, 1),
            ),
        )

    four_vertices = vertices(four)
    four_pairs = tuple(combinations(four_vertices, 2))
    model.add_allowed_assignments(
        [label_variables[edge] for edge in four_pairs],
        (branch_code,) if branch_code is not None else lanner_placement_codes(4),
    )

    conditional_A5_encoding_count = 0
    star_degrees = tuple(sorted(set(permutations((3, 2, 1, 1, 1)))))
    path_degrees = tuple(sorted(set(permutations((2, 2, 2, 1, 1)))))
    cycle_degrees = ((2, 2, 2, 2, 2),)
    endpoint_degree_pairs = ((1, 1), (1, 2), (2, 1))
    star_five_degree_pairs = ((1, 2), (2, 1))
    for five, missing in zip(five_ground, missing_five, strict=True):
        if contains_base(five, fours):
            continue
        support = vertices(five)
        edges = tuple(combinations(support, 2))
        degree = tuple(
            model.new_int_var(0, 4, f"lanner5_{five}_degree_{vertex}")
            for vertex in support
        )
        for local, vertex in enumerate(support):
            model.add(
                degree[local]
                == sum(nonright[pair(vertex, other)] for other in support if other != vertex)
            )
        types = tuple(
            model.new_bool_var(f"lanner5_{five}_type_{index}")
            for index in range(5)
        )
        model.add(sum(types) == missing)
        model.add_allowed_assignments(degree, star_degrees).only_enforce_if(types[0])
        for type_variable in types[1:4]:
            model.add_allowed_assignments(degree, path_degrees).only_enforce_if(
                type_variable
            )
        model.add_allowed_assignments(degree, cycle_degrees).only_enforce_if(types[4])

        label_counts = (
            (0, 1),  # branched tree: three 3-edges and a terminal 5-edge
            (0, 1),  # path 3-3-3-5
            (1, 1),  # path 4-3-3-5
            (0, 2),  # path 5-3-3-5
            (1, 0),  # five-cycle with one 4-edge
        )
        for type_variable, (four_count, five_count) in zip(
            types, label_counts, strict=True
        ):
            model.add(
                sum(label_four[edge] for edge in edges) == four_count
            ).only_enforce_if(type_variable)
            model.add(
                sum(label_five[edge] for edge in edges) == five_count
            ).only_enforce_if(type_variable)

        tree = model.new_bool_var(f"lanner5_{five}_tree")
        model.add(tree == sum(types[:4]))
        for left, middle, right in combinations(support, 3):
            model.add_bool_or(
                [
                    tree.Not(),
                    nonright[pair(left, middle)].Not(),
                    nonright[pair(left, right)].Not(),
                    nonright[pair(middle, right)].Not(),
                ]
            )

        local_index = {vertex: index for index, vertex in enumerate(support)}
        for edge in edges:
            endpoint_degrees = [
                degree[local_index[edge[0]]], degree[local_index[edge[1]]]
            ]
            model.add_allowed_assignments(
                endpoint_degrees, star_five_degree_pairs
            ).only_enforce_if([types[0], label_five[edge]])
            for type_index in (1, 2, 3):
                model.add_allowed_assignments(
                    endpoint_degrees, endpoint_degree_pairs
                ).only_enforce_if([types[type_index], label_five[edge]])
            model.add_allowed_assignments(
                endpoint_degrees, endpoint_degree_pairs
            ).only_enforce_if([types[2], label_four[edge]])
        conditional_A5_encoding_count += 1

    face_variables: dict[int, cp_model.IntVar] = {}
    lower_literals: dict[tuple[tuple[int, int], int], cp_model.IntVar] = {}

    def face_variable(restriction: int) -> cp_model.IntVar:
        known = face_variables.get(restriction)
        if known is not None:
            return known
        variable = model.new_bool_var(f"finite_cut_face_{restriction}")
        witnesses = tuple(
            live[index]
            for index, facet in enumerate(facets)
            if not facet & restriction
        )
        if witnesses:
            model.add(sum(witnesses) >= 1).only_enforce_if(variable)
            model.add(sum(witnesses) == 0).only_enforce_if(variable.Not())
        else:
            model.add(variable == 0)
        face_variables[restriction] = variable
        return variable

    def lower_literal(edge: tuple[int, int], value: int) -> cp_model.IntVar:
        key = (edge, value)
        known = lower_literals.get(key)
        if known is not None:
            return known
        variable = model.new_bool_var(
            f"label_{edge[0]}_{edge[1]}_below_{value}"
        )
        model.add(label_variables[edge] <= value - 1).only_enforce_if(variable)
        model.add(label_variables[edge] >= value).only_enforce_if(variable.Not())
        lower_literals[key] = variable
        return variable

    rounds: list[dict[str, object]] = []
    installed_patterns: set[AbstractPattern] = set()
    installed_expanded_cut_count = 0

    def install_patterns(
        patterns: Sequence[AbstractPattern],
    ) -> tuple[int, str]:
        nonlocal installed_expanded_cut_count
        digest = sha256()
        count = 0
        for pattern in sorted(patterns):
            if pattern in installed_patterns:
                continue
            order, code = pattern
            for target in combinations(range(N), order):
                restriction = mask(target)
                if contains_base(restriction, fours):
                    continue
                target_pairs = tuple(combinations(target, 2))
                for placement in _code_orbit(code):
                    positive = tuple(
                        (edge, value)
                        for edge, value in zip(
                            target_pairs, placement, strict=True
                        )
                        if value > 2
                    )
                    model.add_bool_or(
                        [
                            face_variable(restriction).Not(),
                            *(
                                lower_literal(edge, value)
                                for edge, value in positive
                            ),
                        ]
                    )
                    digest.update(
                        json.dumps(
                            [
                                restriction,
                                [
                                    [edge[0], edge[1], value]
                                    for edge, value in positive
                                ],
                            ],
                            separators=(",", ":"),
                        ).encode()
                    )
                    digest.update(b"\n")
                    count += 1
            installed_patterns.add(pattern)
        installed_expanded_cut_count += count
        return count, digest.hexdigest()

    seed_expansion_transcript = []
    for seed_row in seeded_rounds:
        row = dict(seed_row)
        raw_patterns = row.get("new_abstract_nonfinite_patterns")
        if raw_patterns is None:
            raw_patterns = [
                _pattern_document(_pattern_from_cut(cut))
                for cut in row["new_global_monotone_face_cuts"]
            ]
        patterns = tuple(
            sorted(
                {
                    _pattern_from_document(document)
                    for document in raw_patterns
                    if maximum_seed_pattern_order is None
                    or int(document.get("order", 0))
                    <= maximum_seed_pattern_order
                }
                - installed_patterns
            )
        )
        expanded_count, expanded_digest = install_patterns(patterns)
        seed_expansion_transcript.append(
            {
                "seed_round": seed_row.get("round"),
                "new_abstract_pattern_count": len(patterns),
                "expanded_global_face_cut_count": expanded_count,
                "expanded_global_face_cut_digest": expanded_digest,
            }
        )
        if "expanded_global_face_cut_count" in row and row[
            "expanded_global_face_cut_count"
        ] != expanded_count:
            raise ValueError("resumed expanded cut count changed")
        if "expanded_global_face_cut_digest" in row and row[
            "expanded_global_face_cut_digest"
        ] != expanded_digest:
            raise ValueError("resumed expanded cut digest changed")
        row.setdefault(
            "legacy_literal_face_cut_count",
            row.get("new_global_monotone_face_cut_count"),
        )
        row["new_abstract_nonfinite_patterns"] = [
            _pattern_document(pattern) for pattern in patterns
        ]
        row["new_abstract_nonfinite_pattern_count"] = len(patterns)
        row["abstract_nonfinite_pattern_digest"] = _digest(
            row["new_abstract_nonfinite_patterns"]
        )
        row["expanded_global_face_cut_count"] = expanded_count
        row["expanded_global_face_cut_digest"] = expanded_digest
        row["new_global_monotone_face_cut_count"] = expanded_count
        if seed_rounds_are_history:
            rounds.append(row)
    seeded_abstract_pattern_count = len(installed_patterns)
    seeded_expanded_cut_count = installed_expanded_cut_count
    seeded_pattern_documents = [
        _pattern_document(pattern) for pattern in sorted(installed_patterns)
    ]
    stopped_reason = "not_started"
    terminal_status = "not_started"
    terminal_solver_wall_time_seconds: float | None = None
    terminal_solver_branch_count: int | None = None
    terminal_solver_conflict_count: int | None = None
    survivor: dict[str, object] | None = None

    def result_document() -> dict[str, object]:
        branch_closed = (
            stopped_reason == "terminal_joint_incidence_label_infeasible"
        )
        source_closed = branch_closed and a4_orbit_index is None
        return {
            "kind": "h9_14_record0_a1_joint_incidence_label_cegar_exploratory",
            "theorem_certificate": False,
            "record_index": 0,
            "graph6": problem["graph6"],
            "schedule_a": 1,
            "source_orbit_index": source_orbit,
            "representative_A4s": [list(four_vertices)],
            "a4_placement_partition": {
                "literal_placement_count": sum(map(len, a4_orbits)),
                "orbit_count": len(a4_orbits),
                "orbit_sizes": [len(orbit) for orbit in a4_orbits],
                "orbit_representatives": [list(orbit[0]) for orbit in a4_orbits],
                "partition_digest": _digest(a4_orbits),
                "coverage_is_disjoint_and_complete": (
                    sum(map(len, a4_orbits)) == 108
                    and len({code for orbit in a4_orbits for code in orbit}) == 108
                ),
            },
            "a4_branch": {
                "orbit_index": a4_orbit_index,
                "representative_code": list(branch_code)
                if branch_code is not None
                else None,
                "placement_orbit_size": len(a4_orbits[a4_orbit_index])
                if a4_orbit_index is not None
                else 108,
                "incidence_symmetry_slice": (
                    "exact labelled-A4 stabilizer"
                    if a4_orbit_index is not None
                    else "full unlabelled source stabilizer"
                ),
                "lex_generator_step_counts": list(branch_lex_counts),
                "lossless_reason": (
                    "Move the A4 placement to its unique placement-orbit "
                    "representative, then choose a lexicographically least "
                    "incidence vector under the exact-code stabilizer."
                    if a4_orbit_index is not None
                    else "Every full source orbit has a lexicographically least incidence vector."
                ),
            },
            "dual_model_counts": {
                "candidate_facet_complements": len(facets),
                "independent_five_sets": len(five_ground),
                "conditional_minimal_A5_compact_lanner_encodings": (
                    conditional_A5_encoding_count
                ),
                "global_meeting_label_variables": len(label_variables),
                "face_witness_constraints": data["face_witness_constraint_count"],
                "lanner_order_constraints": data["lanner_order_constraint_count"],
                "ridge_constraints": data["ridge_constraint_count"],
            },
            "completed_cut_round_count": len(rounds),
            "seeded_abstract_nonfinite_pattern_count": (
                seeded_abstract_pattern_count
            ),
            "maximum_seed_pattern_order": maximum_seed_pattern_order,
            "seeded_abstract_nonfinite_patterns": seeded_pattern_documents,
            "seeded_abstract_nonfinite_pattern_digest": _digest(
                seeded_pattern_documents
            ),
            "seeded_global_monotone_face_cut_count": seeded_expanded_cut_count,
            "seed_expansion_transcript": seed_expansion_transcript,
            "seed_expansion_transcript_digest": _digest(
                seed_expansion_transcript
            ),
            "installed_abstract_nonfinite_pattern_count": len(
                installed_patterns
            ),
            "installed_global_monotone_face_cut_count": (
                installed_expanded_cut_count
            ),
            "rounds": rounds,
            "terminal_solver_status": terminal_status,
            "terminal_solver_wall_time_seconds": terminal_solver_wall_time_seconds,
            "terminal_solver_branch_count": terminal_solver_branch_count,
            "terminal_solver_conflict_count": terminal_solver_conflict_count,
            "stopped_reason": stopped_reason,
            "compact_label_survivor": survivor,
            "a4_branch_compact_coxeter_no_go": branch_closed,
            "fixed_source_compact_coxeter_no_go": source_closed,
            "finite_volume_noncompact_checked": False,
            "scope": (
                "Complete joint compact incidence and finite-label formulation "
                "for one fixed a=1 source/A4-label branch. Only terminal "
                "INFEASIBLE closes the branch; all placement branches must "
                "close to infer a fixed-source theorem. An all-elliptic label "
                "model still requires exact Gram."
            ),
        }

    def write_checkpoint() -> None:
        if checkpoint_path is not None:
            checkpoint_path.write_text(
                json.dumps(result_document(), indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )

    while maximum_rounds is None or len(rounds) < maximum_rounds:
        solver = cp_model.CpSolver()
        solver.parameters.num_search_workers = 1
        if seconds_per_solve is not None:
            solver.parameters.max_time_in_seconds = seconds_per_solve
        status = solver.solve(model)
        terminal_status = solver.status_name(status)
        terminal_solver_wall_time_seconds = solver.wall_time
        terminal_solver_branch_count = solver.num_branches
        terminal_solver_conflict_count = solver.num_conflicts
        if status == cp_model.INFEASIBLE:
            stopped_reason = "terminal_joint_incidence_label_infeasible"
            break
        if status not in (cp_model.FEASIBLE, cp_model.OPTIMAL):
            stopped_reason = "solver_unknown_or_timeout"
            break

        selected_facets = tuple(
            facet
            for index, facet in enumerate(facets)
            if solver.value(live[index])
        )
        additions = tuple(
            five
            for five, variable in zip(five_ground, missing_five, strict=True)
            if solver.value(variable) and not contains_base(five, fours)
        )
        blocker = minimal_transversals(selected_facets)
        expected_blocker = {*problem["edge_rows"], four, *additions}
        if set(blocker) != expected_blocker:
            raise AssertionError("joint target failed exact blocker involution")
        labels = {
            edge: solver.value(variable)
            for edge, variable in label_variables.items()
        }
        actual_vertices = tuple(
            tuple(vertex for vertex in range(N) if not facet & (1 << vertex))
            for facet in selected_facets
        )
        nonfinite = tuple(
            actual
            for actual in actual_vertices
            if not is_finite_coxeter_labeling(actual, labels)
        )
        if not nonfinite:
            survivor = {
                "representative_additional_A5s": [
                    list(vertices(value)) for value in additions
                ],
                "selected_facet_complements": [
                    list(vertices(value)) for value in selected_facets
                ],
                "meeting_labels": [
                    [left, right, labels[(left, right)]]
                    for left, right in meeting_pairs
                ],
                "target_digest": _digest(additions),
                "rank_nine_actual_vertex_count": len(actual_vertices),
            }
            stopped_reason = "all_actual_vertices_finite_requires_exact_Gram"
            break

        found_patterns: set[AbstractPattern] = set()
        for actual in nonfinite:
            subset = _minimal_nonfinite(actual, labels)
            found_patterns.add(_canonical_pattern(subset, labels))
        new_patterns = tuple(sorted(found_patterns - installed_patterns))
        if not new_patterns:
            raise AssertionError(
                "nonfinite model produced no new abstract finite cut"
            )
        expanded_count, expanded_digest = install_patterns(new_patterns)
        pattern_documents = [
            _pattern_document(pattern) for pattern in new_patterns
        ]
        rounds.append(
            {
                "round": len(rounds),
                "representative_target_digest": _digest(additions),
                "nonfinite_actual_vertex_count": len(nonfinite),
                "new_abstract_nonfinite_pattern_count": len(new_patterns),
                "new_abstract_nonfinite_patterns": pattern_documents,
                "abstract_nonfinite_pattern_digest": _digest(
                    pattern_documents
                ),
                "expanded_global_face_cut_count": expanded_count,
                "expanded_global_face_cut_digest": expanded_digest,
                "new_global_monotone_face_cut_count": expanded_count,
                "solver_wall_time_seconds": solver.wall_time,
                "solver_branch_count": solver.num_branches,
                "solver_conflict_count": solver.num_conflicts,
            }
        )
        print(
            f"round {len(rounds) - 1}: nonfinite={len(nonfinite)} "
            f"new_patterns={len(new_patterns)} expanded_cuts={expanded_count}",
            flush=True,
        )
        write_checkpoint()

    if (
        maximum_rounds is not None
        and len(rounds) >= maximum_rounds
        and stopped_reason == "not_started"
    ):
        stopped_reason = "round_cutoff"
    write_checkpoint()
    return result_document()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_orbit", type=int, choices=range(4))
    parser.add_argument("--max-rounds", type=int)
    parser.add_argument("--seconds-per-solve", type=float)
    parser.add_argument("--a4-orbit-index", type=int)
    parser.add_argument("--seed-cuts", type=Path)
    parser.add_argument("--resume", type=Path)
    parser.add_argument("--max-seed-pattern-order", type=int)
    parser.add_argument("--write", type=Path)
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()
    if args.seed_cuts is not None and args.resume is not None:
        raise ValueError("choose either --seed-cuts or --resume")
    seeded_rounds: Sequence[dict[str, object]] = ()
    seed_rounds_are_history = False
    source_path = args.resume or args.seed_cuts
    if source_path is not None:
        seed = json.loads(source_path.read_text(encoding="utf-8"))
        if seed.get("source_orbit_index") != args.source_orbit:
            raise ValueError("cut checkpoint belongs to another source orbit")
        seeded_rounds = seed["rounds"]
        if args.resume is not None:
            branch = seed.get("a4_branch", {})
            if branch.get("orbit_index") != args.a4_orbit_index:
                raise ValueError("resume checkpoint belongs to another A4 branch")
            seed_rounds_are_history = True
    result = run_joint_cegar(
        args.source_orbit,
        args.max_rounds,
        args.seconds_per_solve,
        a4_orbit_index=args.a4_orbit_index,
        seeded_rounds=seeded_rounds,
        seed_rounds_are_history=seed_rounds_are_history,
        maximum_seed_pattern_order=args.max_seed_pattern_order,
        checkpoint_path=args.write,
    )
    if args.write is not None:
        args.write.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    if args.summary_only:
        result = {
            "source_orbit_index": result["source_orbit_index"],
            "completed_cut_round_count": result["completed_cut_round_count"],
            "installed_global_monotone_face_cut_count": result[
                "installed_global_monotone_face_cut_count"
            ],
            "installed_abstract_nonfinite_pattern_count": result[
                "installed_abstract_nonfinite_pattern_count"
            ],
            "round_summaries": [
                {
                    "round": row["round"],
                    "nonfinite_actual_vertex_count": row[
                        "nonfinite_actual_vertex_count"
                    ],
                    "new_global_monotone_face_cut_count": row[
                        "new_global_monotone_face_cut_count"
                    ],
                    "new_abstract_nonfinite_pattern_count": row[
                        "new_abstract_nonfinite_pattern_count"
                    ],
                    "solver_wall_time_seconds": row["solver_wall_time_seconds"],
                }
                for row in result["rounds"]
            ],
            "terminal_solver_status": result["terminal_solver_status"],
            "stopped_reason": result["stopped_reason"],
            "compact_label_survivor_found": result["compact_label_survivor"]
            is not None,
            "a4_branch_compact_coxeter_no_go": result[
                "a4_branch_compact_coxeter_no_go"
            ],
            "fixed_source_compact_coxeter_no_go": result[
                "fixed_source_compact_coxeter_no_go"
            ],
            "finite_volume_noncompact_checked": False,
        }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
