#!/usr/bin/env python3
"""Hostile replay of the Record 30 orbit-10 exactly-two-cross theorem.

This implementation does not import the production enumerator or the Coxeter
label module.  It regenerates the classical Lanner placements from the nine
order-four and five order-five representatives, recognizes finite Coxeter
forests directly, uses the opposite support-join order, and certifies the
decisive fixed principal minor over a second finite field.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations, product
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record30_a2_orbit10_two_cross.json"

N = 14
DIMENSION = 9
DOTTED = ((0, 11), (1, 12), (2, 13))
SUPPORTS = (
    (0, 3, 4, 5),
    (6, 7, 8, 11),
    (3, 4, 5, 9, 10),
    (6, 7, 8, 9, 10),
)
COMPONENTS = (
    (0, 3, 4, 5, 6, 7, 8, 9, 10, 11),
    (1, 12),
    (2, 13),
)
JOIN_ORDER = (3, 1, 2, 0)
LABELS = (3, 4, 5)
Pair = tuple[int, int]
State = tuple[tuple[Pair, int], ...]

MODULUS = 1129
SQRT2 = 406
SQRT5 = 474


LANNER_REPRESENTATIVES: dict[int, tuple[tuple[tuple[Pair, int], ...], ...]] = {
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
    5: (
        (((0, 1), 3), ((0, 2), 3), ((0, 3), 3), ((1, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 4), ((2, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 5), ((2, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 3), ((3, 4), 4)),
    ),
}


def pair(left: int, right: int) -> Pair:
    return (left, right) if left < right else (right, left)


def canonical_sha(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def placement_codes(order: int) -> tuple[tuple[int, ...], ...]:
    source_pairs = tuple(combinations(range(order), 2))
    result = set()
    for representative in LANNER_REPRESENTATIVES[order]:
        base = {edge: 2 for edge in source_pairs}
        base.update(representative)
        for move in permutations(range(order)):
            moved = {
                pair(move[left], move[right]): label
                for (left, right), label in base.items()
            }
            result.add(tuple(moved[edge] for edge in source_pairs))
    expected = {4: 108, 5: 420}[order]
    if len(result) != expected:
        raise AssertionError("independent Lanner placement count changed")
    return tuple(sorted(result))


def placements(support: Sequence[int]) -> tuple[State, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True)) for code in placement_codes(len(support))
    )


def path_labels(component: set[int], adjacency: Mapping[int, list[tuple[int, int]]]) -> tuple[int, ...]:
    current = min(vertex for vertex in component if len(adjacency[vertex]) == 1)
    previous = -1
    result = []
    while True:
        onward = [(vertex, label) for vertex, label in adjacency[current] if vertex != previous]
        if not onward:
            break
        if len(onward) != 1:
            raise AssertionError("independent path unexpectedly branched")
        following, label = onward[0]
        result.append(label)
        previous, current = current, following
    forward = tuple(result)
    return min(forward, tuple(reversed(forward)))


def branch_arms(component: set[int], adjacency: Mapping[int, list[tuple[int, int]]]) -> tuple[int, ...]:
    branch = next(vertex for vertex in component if len(adjacency[vertex]) == 3)
    result = []
    for neighbor, _ in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(vertex for vertex, _ in adjacency[current] if vertex != previous)
            previous, current = current, following
            length += 1
        result.append(length)
    return tuple(sorted(result))


@lru_cache(maxsize=1_048_576)
def finite_code(order: int, code: tuple[int, ...]) -> bool:
    edges = tuple(combinations(range(order), 2))
    labels = dict(zip(edges, code, strict=True))
    adjacency = {
        vertex: [
            (other, labels[pair(vertex, other)])
            for other in range(order)
            if other != vertex and labels[pair(vertex, other)] > 2
        ]
        for vertex in range(order)
    }
    unseen = set(range(order))
    components = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        component = {seed}
        stack = [seed]
        while stack:
            vertex = stack.pop()
            new = {other for other, _ in adjacency[vertex] if other in unseen}
            unseen.difference_update(new)
            component.update(new)
            stack.extend(new)
        components.append(component)

    for component in components:
        size = len(component)
        if size == 1:
            continue
        if sum(len(adjacency[vertex]) for vertex in component) != 2 * (size - 1):
            return False
        degrees = Counter(len(adjacency[vertex]) for vertex in component)
        maximum = max(degrees)
        if maximum <= 2:
            sequence = path_labels(component, adjacency)
            if size == 2:
                continue
            if all(label == 3 for label in sequence):
                continue
            if sequence == (3,) * (size - 2) + (4,):
                continue
            if size == 3 and sequence == (3, 5):
                continue
            if size == 4 and sequence in ((3, 3, 5), (3, 4, 3)):
                continue
            return False
        if degrees[3] != 1 or maximum > 3:
            return False
        if any(label != 3 for vertex in component for _, label in adjacency[vertex]):
            return False
        arms = branch_arms(component, adjacency)
        if arms != (1, 1, size - 3) and arms not in ((1, 2, 2), (1, 2, 3), (1, 2, 4)):
            return False
    return True


def restriction_finite(vertices: Sequence[int], labels: Mapping[Pair, int]) -> bool:
    vertices = tuple(vertices)
    return finite_code(
        len(vertices),
        tuple(labels[pair(a, b)] for a, b in combinations(vertices, 2)),
    )


def actual_vertices() -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(set(face) for face in DOTTED + SUPPORTS)
    return tuple(
        chosen
        for chosen in combinations(range(N), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


def sparse_states(vertices: tuple[tuple[int, ...], ...]) -> tuple[State, ...]:
    states: tuple[State, ...] = ((),)
    union: set[int] = set()
    for support_index in JOIN_ORDER:
        support = SUPPORTS[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {tuple(facet for facet in vertex if facet in union) for vertex in vertices},
                key=lambda row: (len(row), row),
            )
        )
        passing = set()
        for state in states:
            old = dict(state)
            for option in placements(support):
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                merged = old | dict(option)
                labels = {pair(a, b): 2 for a, b in combinations(cluster, 2)}
                labels.update(merged)
                if all(restriction_finite(row, labels) for row in restrictions):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
    return states


def component_categories() -> tuple[tuple[Pair, ...], ...]:
    return tuple(
        tuple(pair(a, b) for a in COMPONENTS[left] for b in COMPONENTS[right])
        for left, right in ((0, 1), (0, 2), (1, 2))
    )


def edge_placements() -> tuple[tuple[Pair, Pair], ...]:
    categories = component_categories()
    return tuple(
        (first, second)
        for left, right in combinations(range(3), 2)
        for first in categories[left]
        for second in categories[right]
    )


def base_labels(state: State) -> dict[Pair, int]:
    dotted = set(DOTTED)
    labels = {
        pair(a, b): 2
        for a, b in combinations(range(N), 2)
        if pair(a, b) not in dotted
    }
    labels.update(state)
    return labels


def gram_entry(label: int) -> int:
    return {2: 0, 3: -2, 4: -2 * SQRT2, 5: -(1 + SQRT5)}[label]


def determinant_mod(matrix: list[list[int]]) -> int:
    p = MODULUS
    matrix = [[entry % p for entry in row] for row in matrix]
    result = 1
    for column in range(len(matrix)):
        pivot = next((row for row in range(column, len(matrix)) if matrix[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            result = -result
        value = matrix[column][column]
        result = result * value % p
        inverse = pow(value, p - 2, p)
        for row in range(column + 1, len(matrix)):
            multiplier = matrix[row][column] * inverse % p
            for entry in range(column, len(matrix)):
                matrix[row][entry] = (
                    matrix[row][entry] - multiplier * matrix[column][entry]
                ) % p
    return result % p


def fixed_residue(labels: Mapping[Pair, int]) -> int:
    order = tuple(range(3, 14))
    matrix = [
        [4 if a == b else gram_entry(labels[pair(a, b)]) for b in order]
        for a in order
    ]
    return determinant_mod(matrix)


def replay() -> dict[str, object]:
    if any(MODULUS % divisor == 0 for divisor in range(2, int(MODULUS**0.5) + 1)):
        raise AssertionError("independent modulus is composite")
    if SQRT2 * SQRT2 % MODULUS != 2 or SQRT5 * SQRT5 % MODULUS != 5:
        raise AssertionError("independent radical residues changed")

    vertices = actual_vertices()
    states = sparse_states(vertices)
    placements_two = edge_placements()
    vertex_sets = tuple(frozenset(vertex) for vertex in vertices)
    affected = {
        edge: tuple(index for index, vertex in enumerate(vertex_sets) if set(edge) <= vertex)
        for category in component_categories()
        for edge in category
    }

    single_pass = 0
    residual_trials = 0
    survivors = []
    residues: Counter[int] = Counter()
    per_state = Counter()
    label_pairs = Counter()
    zero_residue = []
    for state_index, state in enumerate(states):
        labels = base_labels(state)
        valid_single = {}
        for category in component_categories():
            for edge in category:
                for label in LABELS:
                    labels[edge] = label
                    valid = all(restriction_finite(vertices[index], labels) for index in affected[edge])
                    labels[edge] = 2
                    valid_single[(edge, label)] = valid
                    single_pass += valid
        state_count = 0
        for first, second in placements_two:
            common = tuple(sorted(set(affected[first]) & set(affected[second])))
            for label_first, label_second in product(LABELS, repeat=2):
                if not valid_single[(first, label_first)] or not valid_single[(second, label_second)]:
                    continue
                residual_trials += 1
                labels[first] = label_first
                labels[second] = label_second
                if all(restriction_finite(vertices[index], labels) for index in common):
                    payload = [state_index, list(first), label_first, list(second), label_second]
                    survivors.append(payload)
                    state_count += 1
                    label_pairs[tuple(sorted((label_first, label_second)))] += 1
                    residue = fixed_residue(labels)
                    residues[residue] += 1
                    if not residue:
                        zero_residue.append(payload)
                labels[first] = 2
                labels[second] = 2
        per_state[state_count] += 1

    survivors.sort()
    result = {
        "kind": "h9_14_record30_a2_orbit10_two_cross_independent_replay",
        "implementation": {
            "imports_production_enumerator": False,
            "imports_coxeter_label_module": False,
            "support_join_order": list(JOIN_ORDER),
            "finite_recognizer": "standalone ADEBFH forest classifier",
            "gram_modulus": MODULUS,
            "sqrt2_residue": SQRT2,
            "sqrt5_residue": SQRT5,
        },
        "actual_rank_nine_vertex_count": len(vertices),
        "terminal_sparse_state_count": len(states),
        "terminal_states_sha256": canonical_sha(states),
        "edge_pair_placement_count": len(placements_two),
        "literal_labelled_trial_count": len(states) * len(placements_two) * 9,
        "single_edge_passing_count": single_pass,
        "two_edge_residual_trial_count": residual_trials,
        "local_survivor_count": len(survivors),
        "local_survivor_sha256": canonical_sha(survivors),
        "local_survivor_histogram_per_sparse_state": {
            str(value): count for value, count in sorted(per_state.items())
        },
        "local_survivor_label_pair_histogram": {
            f"{key[0]},{key[1]}": value for key, value in sorted(label_pairs.items())
        },
        "fixed_order": list(range(3, 14)),
        "fixed_determinant_residue_histogram": {
            str(residue): count for residue, count in sorted(residues.items())
        },
        "zero_fixed_determinant_residue_count": len(zero_residue),
        "all_local_survivors_rank_at_least_eleven": not zero_residue,
    }
    if ARCHIVE.exists():
        source = json.loads(ARCHIVE.read_text(encoding="utf-8"))
        shell = source["two_cross_shell"]
        sparse = source["sparse_gate"]
        gate = source["fixed_minor_gate"]
        result["agrees_with_frozen_production"] = (
            result["actual_rank_nine_vertex_count"]
            == source["nerve"]["actual_rank_nine_vertex_count"]
            and result["terminal_sparse_state_count"] == sparse["terminal_state_count"]
            and result["terminal_states_sha256"] == sparse["terminal_states_sha256"]
            and result["literal_labelled_trial_count"] == shell["literal_labelled_trial_count"]
            and result["single_edge_passing_count"] == shell["single_edge_passing_count"]
            and result["two_edge_residual_trial_count"] == shell["two_edge_residual_trial_count"]
            and result["local_survivor_count"] == shell["local_survivor_count"]
            and result["local_survivor_sha256"] == shell["local_survivor_sha256"]
            and gate["no_fixed_minor_witness_count"] == 0
            and result["all_local_survivors_rank_at_least_eleven"]
        )
    return result


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
