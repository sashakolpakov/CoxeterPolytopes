#!/usr/bin/env python3
"""Exact whole-nerve replay for the two exceptional Record-30 a=2 nerves.

This replays positive-A4 orbit 7 target 1 and orbit 10 target 0 without any
``coxeter_polytopes`` import.  A standalone Lanner/finite-diagram engine
generates the 648 sparse states.  The full nerve automorphism action reduces
them to exact state orbits; every locally finite connected shell of additional
non-right edges is then enumerated.

For every quotient shell assignment an eleven-node dotted-free principal
Gram block is materialized.  Its determinant is computed in the integral
basis ``1,sqrt(2),sqrt(5),sqrt(10)`` by a division-free subset recurrence.
The exact four integer coefficients, exact integer field norm, and a rigorous
rational-interval sign are recorded.  A second-prime residue is only a replay
check, never the proof of nonsingularity.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from fractions import Fraction
from math import isqrt
import json
from pathlib import Path
import sys
from typing import Iterable, Mapping, Sequence

import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import scripts.independent_h9_14_record30_a2_orbit14_labels as labels_engine


N = 14
DIMENSION = 9
DOTTED = ((0, 11), (1, 12), (2, 13))
FINITE_NONRIGHT = (3, 4, 5)
JOIN_ORDER = (3, 1, 2, 0)
FIXED_ORDER = tuple(range(3, 14))
MODULUS = 1129
SQRT2_RESIDUE = 406
SQRT5_RESIDUE = 474

Pair = tuple[int, int]
State = tuple[tuple[Pair, int], ...]
Assignment = tuple[tuple[int, int], ...]
AlgebraicInteger = tuple[int, int, int, int]

CASES: dict[str, dict[str, object]] = {
    "orbit7_target1": {
        "positive_A4_orbit_index": 7,
        "target_index": 1,
        "fixed_A4_masks": [57, 2104],
        "additional_A5_masks": [248, 1984],
        "supports": (
            (0, 3, 4, 5),
            (3, 4, 5, 11),
            (3, 4, 5, 6, 7),
            (6, 7, 8, 9, 10),
        ),
    },
    "orbit10_target0": {
        "positive_A4_orbit_index": 10,
        "target_index": 0,
        "fixed_A4_masks": [57, 2496],
        "additional_A5_masks": [1592, 1984],
        "supports": (
            (0, 3, 4, 5),
            (6, 7, 8, 11),
            (3, 4, 5, 9, 10),
            (6, 7, 8, 9, 10),
        ),
    },
}


def _pair(left: int, right: int) -> Pair:
    return (left, right) if left < right else (right, left)


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _canonical_sha(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _actual_vertices(supports: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(map(frozenset, DOTTED + tuple(supports)))
    return tuple(
        chosen
        for chosen in combinations(range(N), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


def _restriction_finite(
    restriction: Sequence[int], finite_labels: Mapping[Pair, int]
) -> bool:
    return labels_engine._finite_code(
        len(restriction),
        tuple(
            finite_labels[_pair(left, right)]
            for left, right in combinations(restriction, 2)
        ),
    )


def _sparse_states(
    supports: Sequence[Sequence[int]],
    vertices: Sequence[Sequence[int]],
) -> tuple[State, ...]:
    states: tuple[State, ...] = ((),)
    union: set[int] = set()
    for support_index in JOIN_ORDER:
        support = tuple(supports[support_index])
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {
                    tuple(facet for facet in vertex if facet in union)
                    for vertex in vertices
                },
                key=lambda row: (len(row), row),
            )
        )
        passing: set[State] = set()
        for state in states:
            old = dict(state)
            for option in labels_engine._placement_options(support):
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                merged = old | dict(option)
                finite_labels = {
                    _pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                finite_labels.update(merged)
                if all(
                    _restriction_finite(restriction, finite_labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
    return states


def _automorphisms(supports: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    graph = nx.Graph()
    for vertex in range(N):
        graph.add_node(("v", vertex), color="facet")
    for index, face in enumerate(DOTTED + tuple(supports)):
        graph.add_node(("f", index), color=f"missing_{len(face)}")
        for vertex in face:
            graph.add_edge(("v", vertex), ("f", index))
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph,
        graph,
        node_match=lambda left, right: left["color"] == right["color"],
    )
    result = {
        tuple(move[("v", vertex)][1] for vertex in range(N))
        for move in matcher.isomorphisms_iter()
    }
    if len(result) != 1152:
        raise AssertionError("exceptional nerve automorphism order changed")
    return tuple(sorted(result))


def _move_state(state: State, move: Sequence[int]) -> State:
    return tuple(
        sorted(
            (_pair(move[edge[0]], move[edge[1]]), label)
            for edge, label in state
        )
    )


def _state_orbits(
    states: tuple[State, ...], moves: Sequence[Sequence[int]]
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    lookup = {state: index for index, state in enumerate(states)}
    remaining = set(range(len(states)))
    result = []
    while remaining:
        representative = min(remaining)
        orbit = {
            lookup[_move_state(states[representative], move)] for move in moves
        }
        result.append((representative, tuple(sorted(orbit))))
        remaining.difference_update(orbit)
    return tuple(result)


def _components(
    supports: Sequence[Sequence[int]],
) -> tuple[tuple[int, ...], ...]:
    components: list[set[int]] = []
    for face in DOTTED + tuple(supports):
        support = set(face)
        touching = [
            index
            for index, component in enumerate(components)
            if component & support
        ]
        for index in reversed(touching):
            support.update(components.pop(index))
        components.append(support)
    result = tuple(sorted(tuple(sorted(component)) for component in components))
    if sorted(map(len, result)) != [2, 2, 10]:
        raise AssertionError("exceptional minimal-nonface components changed")
    return result


def _cross_edges(components: Sequence[Sequence[int]]) -> tuple[Pair, ...]:
    return tuple(
        sorted(
            _pair(left, right)
            for first, second in combinations(range(len(components)), 2)
            for left in components[first]
            for right in components[second]
        )
    )


def _residual_internal_pairs(
    supports: Sequence[Sequence[int]], components: Sequence[Sequence[int]]
) -> tuple[Pair, ...]:
    fixed = {
        _pair(left, right)
        for support in supports
        for left, right in combinations(support, 2)
    }
    dotted = {_pair(*edge) for edge in DOTTED}
    return tuple(
        sorted(
            edge
            for component in components
            for edge in combinations(component, 2)
            if edge not in fixed and edge not in dotted
        )
    )


def _base_labels(state: State) -> dict[Pair, int]:
    dotted = {_pair(*edge) for edge in DOTTED}
    result = {
        _pair(left, right): 2
        for left, right in combinations(range(N), 2)
        if _pair(left, right) not in dotted
    }
    result.update(state)
    return result


def _all_vertices_finite(
    state: State,
    assignment: Assignment,
    edges: Sequence[Pair],
    vertices: Sequence[Sequence[int]],
) -> bool:
    finite_labels = _base_labels(state)
    for edge_index, label in assignment:
        finite_labels[edges[edge_index]] = label
    return all(
        _restriction_finite(vertex, finite_labels) for vertex in vertices
    )


def _single_choices(
    state: State,
    edges: Sequence[Pair],
    vertices: Sequence[Sequence[int]],
) -> frozenset[tuple[int, int]]:
    return frozenset(
        (index, label)
        for index in range(len(edges))
        for label in FINITE_NONRIGHT
        if _all_vertices_finite(state, ((index, label),), edges, vertices)
    )


def _category(edge: Pair, components: Sequence[Sequence[int]]) -> tuple[int, int]:
    membership = {
        vertex: index
        for index, component in enumerate(components)
        for vertex in component
    }
    return tuple(sorted((membership[edge[0]], membership[edge[1]])))


def _connected_two_seeds(
    state: State,
    choices: frozenset[tuple[int, int]],
    edges: Sequence[Pair],
    vertices: Sequence[Sequence[int]],
    components: Sequence[Sequence[int]],
) -> set[Assignment]:
    result = set()
    for first, second in combinations(sorted(choices), 2):
        if first[0] == second[0]:
            continue
        if _category(edges[first[0]], components) == _category(
            edges[second[0]], components
        ):
            continue
        assignment = tuple(sorted((first, second)))
        if _all_vertices_finite(state, assignment, edges, vertices):
            result.add(assignment)
    return result


def _extend(
    state: State,
    current: set[Assignment],
    choices: frozenset[tuple[int, int]],
    edges: Sequence[Pair],
    vertices: Sequence[Sequence[int]],
) -> set[Assignment]:
    result = set()
    tested = set()
    for assignment in current:
        occupied = {index for index, _label in assignment}
        for choice in choices:
            if choice[0] in occupied:
                continue
            candidate = tuple(sorted((*assignment, choice)))
            if candidate in tested:
                continue
            tested.add(candidate)
            if _all_vertices_finite(state, candidate, edges, vertices):
                result.add(candidate)
    return result


def _aadd(left: AlgebraicInteger, right: AlgebraicInteger) -> AlgebraicInteger:
    return tuple(a + b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def _aneg(value: AlgebraicInteger) -> AlgebraicInteger:
    return tuple(-entry for entry in value)  # type: ignore[return-value]


def _amul(left: AlgebraicInteger, right: AlgebraicInteger) -> AlgebraicInteger:
    a, b, c, d = left
    e, f, g, h = right
    return (
        a * e + 2 * b * f + 5 * c * g + 10 * d * h,
        a * f + b * e + 5 * (c * h + d * g),
        a * g + c * e + 2 * (b * h + d * f),
        a * h + d * e + b * g + c * f,
    )


def _gram_entry(label: int) -> AlgebraicInteger:
    return {
        2: (0, 0, 0, 0),
        3: (-2, 0, 0, 0),
        4: (0, -2, 0, 0),
        5: (-1, 0, -1, 0),
    }[label]


@lru_cache(maxsize=None)
def _determinant_from_code(code: tuple[int, ...]) -> AlgebraicInteger:
    order = len(FIXED_ORDER)
    pairs = tuple(combinations(range(order), 2))
    finite_labels = dict(zip(pairs, code, strict=True))
    matrix = []
    for row in range(order):
        values = []
        for column in range(order):
            if row == column:
                values.append((4, 0, 0, 0))
            else:
                values.append(_gram_entry(finite_labels[_pair(row, column)]))
        matrix.append(values)
    zero: AlgebraicInteger = (0, 0, 0, 0)
    one: AlgebraicInteger = (1, 0, 0, 0)
    dynamic = [zero] * (1 << order)
    dynamic[0] = one
    for chosen in range(1 << order):
        row = chosen.bit_count()
        if row == order or dynamic[chosen] == zero:
            continue
        for column in range(order):
            if chosen >> column & 1 or matrix[row][column] == zero:
                continue
            term = _amul(dynamic[chosen], matrix[row][column])
            if (chosen >> (column + 1)).bit_count() & 1:
                term = _aneg(term)
            target = chosen | (1 << column)
            dynamic[target] = _aadd(dynamic[target], term)
    return dynamic[-1]


def _field_norm(value: AlgebraicInteger) -> int:
    a, b, c, d = value
    # First take the norm over sqrt(5), then over sqrt(2).
    u = a * a + 2 * b * b - 5 * (c * c + 2 * d * d)
    v = 2 * a * b - 10 * c * d
    return u * u - 2 * v * v


def _radical_bounds(radicand: int, denominator: int = 10**12) -> tuple[Fraction, Fraction]:
    floor_value = isqrt(radicand * denominator * denominator)
    lower = Fraction(floor_value, denominator)
    upper = Fraction(floor_value + 1, denominator)
    if not lower * lower < radicand < upper * upper:
        raise AssertionError("rational radical enclosure failed")
    return lower, upper


def _exact_sign(value: AlgebraicInteger) -> int:
    intervals = (
        (Fraction(1), Fraction(1)),
        _radical_bounds(2),
        _radical_bounds(5),
        _radical_bounds(10),
    )
    lower = Fraction(0)
    upper = Fraction(0)
    for coefficient, (term_lower, term_upper) in zip(value, intervals, strict=True):
        if coefficient >= 0:
            lower += coefficient * term_lower
            upper += coefficient * term_upper
        else:
            lower += coefficient * term_upper
            upper += coefficient * term_lower
    if lower > 0:
        return 1
    if upper < 0:
        return -1
    raise AssertionError("rational interval did not decide determinant sign")


def _residue(value: AlgebraicInteger) -> int:
    a, b, c, d = value
    return (
        a
        + b * SQRT2_RESIDUE
        + c * SQRT5_RESIDUE
        + d * SQRT2_RESIDUE * SQRT5_RESIDUE
    ) % MODULUS


def _fixed_code(finite_labels: Mapping[Pair, int]) -> tuple[int, ...]:
    return tuple(
        finite_labels[_pair(left, right)]
        for left, right in combinations(FIXED_ORDER, 2)
    )


def _case_replay(name: str, case: Mapping[str, object]) -> dict[str, object]:
    supports = tuple(tuple(face) for face in case["supports"])  # type: ignore[arg-type]
    vertices = _actual_vertices(supports)
    states = _sparse_states(supports, vertices)
    moves = _automorphisms(supports)
    state_orbits = _state_orbits(states, moves)
    components = _components(supports)
    edges = _cross_edges(components)
    internal = _residual_internal_pairs(supports, components)

    internal_passing = []
    for state_index, state in enumerate(states):
        finite_labels = _base_labels(state)
        for edge in internal:
            for label in FINITE_NONRIGHT:
                finite_labels[edge] = label
                if all(
                    _restriction_finite(vertex, finite_labels)
                    for vertex in vertices
                ):
                    internal_passing.append([state_index, list(edge), label])
                finite_labels[edge] = 2

    choice_sets = tuple(_single_choices(state, edges, vertices) for state in states)
    if len(set(choice_sets)) != 1:
        raise AssertionError(f"{name}: single-choice sets vary by sparse state")
    choices = choice_sets[0]
    common_shells: dict[int, set[Assignment]] | None = None
    quotient_rows = []
    determinant_rows = []
    for representative, orbit in state_orbits:
        state = states[representative]
        current = _connected_two_seeds(
            state, choices, edges, vertices, components
        )
        shells = {}
        size = 2
        while current:
            shells[size] = current
            current = _extend(state, current, choices, edges, vertices)
            size += 1
        if common_shells is None:
            common_shells = shells
        elif shells != common_shells:
            raise AssertionError(f"{name}: quotient shell sets disagree")
        quotient_rows.append(
            {
                "representative_state_index": representative,
                "orbit_size": len(orbit),
                "shell_counts": [len(shells[key]) for key in sorted(shells)],
                "first_empty_shell_size": size,
            }
        )
        for shell_size, shell in sorted(shells.items()):
            for assignment in sorted(shell):
                finite_labels = _base_labels(state)
                for edge_index, label in assignment:
                    finite_labels[edges[edge_index]] = label
                code = _fixed_code(finite_labels)
                determinant = _determinant_from_code(code)
                norm = _field_norm(determinant)
                sign = _exact_sign(determinant)
                residue = _residue(determinant)
                if norm == 0 or residue == 0:
                    raise AssertionError(f"{name}: an exact fixed minor vanished")
                determinant_rows.append(
                    {
                        "representative_state_index": representative,
                        "state_orbit_size": len(orbit),
                        "shell_size": shell_size,
                        "assignment": [
                            [list(edges[edge_index]), label]
                            for edge_index, label in assignment
                        ],
                        "determinant_coefficients_1_sqrt2_sqrt5_sqrt10": list(
                            determinant
                        ),
                        "determinant_field_norm": norm,
                        "determinant_sign": sign,
                        "second_prime_residue": residue,
                    }
                )
    if common_shells is None:
        raise AssertionError(f"{name}: no quotient shell was produced")
    shell_counts = {size: len(shell) for size, shell in sorted(common_shells.items())}
    if (
        len(vertices) != 128
        or len(states) != 648
        or len(moves) != 1152
        or internal_passing
        or shell_counts != {2: 56, 3: 128, 4: 38}
        or sum(len(orbit) for _representative, orbit in state_orbits) != 648
    ):
        raise AssertionError(f"{name}: exceptional whole-nerve census changed")
    if len(determinant_rows) != len(state_orbits) * sum(shell_counts.values()):
        raise AssertionError(f"{name}: exact determinant coverage changed")
    determinant_payload = [
        {
            key: row[key]
            for key in (
                "representative_state_index",
                "shell_size",
                "assignment",
                "determinant_coefficients_1_sqrt2_sqrt5_sqrt10",
                "determinant_field_norm",
                "determinant_sign",
            )
        }
        for row in determinant_rows
    ]
    return {
        "case": name,
        "positive_A4_orbit_index": case["positive_A4_orbit_index"],
        "target_index": case["target_index"],
        "fixed_A4_masks": case["fixed_A4_masks"],
        "additional_A5_masks": case["additional_A5_masks"],
        "minimal_nonfaces": [list(face) for face in DOTTED + supports],
        "actual_rank_nine_vertex_count": len(vertices),
        "actual_vertices_sha256": _canonical_sha(vertices),
        "terminal_sparse_state_count": len(states),
        "terminal_states_sha256": _canonical_sha(states),
        "nerve_automorphism_order": len(moves),
        "sparse_state_orbit_count": len(state_orbits),
        "sparse_state_orbit_sizes": [len(orbit) for _rep, orbit in state_orbits],
        "quotient_rows": quotient_rows,
        "minimal_nonface_component_vertex_sets": [list(row) for row in components],
        "residual_internal_pair_count": len(internal),
        "residual_internal_trial_count": len(states) * len(internal) * 3,
        "residual_internal_passing_count": len(internal_passing),
        "cross_meeting_pair_count": len(edges),
        "common_single_choice_count": len(choices),
        "common_single_choices_sha256": _canonical_sha(
            [[list(edges[index]), label] for index, label in sorted(choices)]
        ),
        "connected_shell_counts_per_sparse_state": {
            str(size): count for size, count in shell_counts.items()
        },
        "first_empty_connected_shell_size": 5,
        "exact_fixed_principal_order": list(FIXED_ORDER),
        "exact_quotient_principal_block_count": len(determinant_rows),
        "exact_determinant_rows": determinant_rows,
        "exact_determinant_transcript_sha256": _canonical_sha(determinant_payload),
        "determinant_sign_histogram": {
            str(sign): count
            for sign, count in sorted(
                Counter(row["determinant_sign"] for row in determinant_rows).items()
            )
        },
        "zero_exact_field_norm_count": sum(
            row["determinant_field_norm"] == 0 for row in determinant_rows
        ),
        "zero_second_prime_residue_count": sum(
            row["second_prime_residue"] == 0 for row in determinant_rows
        ),
        "all_literal_assignments_have_rank_at_least_eleven_by_symmetry_transport": True,
        "whole_nerve_compact_coxeter_impossible": True,
    }


def replay() -> dict[str, object]:
    if SQRT2_RESIDUE * SQRT2_RESIDUE % MODULUS != 2:
        raise AssertionError("sqrt(2) second-prime residue changed")
    if SQRT5_RESIDUE * SQRT5_RESIDUE % MODULUS != 5:
        raise AssertionError("sqrt(5) second-prime residue changed")
    cases = [_case_replay(name, case) for name, case in CASES.items()]
    result = {
        "kind": "h9_14_record30_a2_exceptional_whole_nerves_exact_replay",
        "compact_only": True,
        "imports_coxeter_polytopes": False,
        "complete_order_four_lanner_placement_count": len(
            labels_engine._lanner_codes(4)
        ),
        "complete_order_five_lanner_placement_count": len(
            labels_engine._lanner_codes(5)
        ),
        "exact_coefficient_basis": ["1", "sqrt(2)", "sqrt(5)", "sqrt(10)"],
        "exact_determinant_method": (
            "division-free signed subset recurrence over the integral "
            "biquadratic basis; nonzero exact integer field norm and rigorous "
            "rational radical enclosure determine every determinant sign"
        ),
        "second_modulus_check": {
            "modulus": MODULUS,
            "sqrt2_residue": SQRT2_RESIDUE,
            "sqrt5_residue": SQRT5_RESIDUE,
            "proof_role": "redundant replay check only",
        },
        "cases": cases,
        "theorem": {
            "orbit7_target1_whole_nerve_compact_coxeter_impossible": True,
            "orbit10_target0_whole_nerve_compact_coxeter_impossible": True,
            "finite_volume_noncompact_chamber_excluded": False,
        },
        "scope": (
            "All finite labels 2,3,4,5 on the two named compact nerves, "
            "including every locally finite connected cross-component shell. "
            "Exact nonzero eleven-node principal determinants force Gram rank "
            "at least eleven, impossible for a chamber in H^9. No ideal-vertex "
            "or finite-volume noncompact conclusion is made."
        ),
    }
    result["summary_sha256"] = _canonical_sha(result)
    return result


def main() -> None:
    print(json.dumps(replay(), sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
