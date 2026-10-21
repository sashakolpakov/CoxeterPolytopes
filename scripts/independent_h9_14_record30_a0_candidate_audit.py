#!/usr/bin/env python3
"""Standalone audit of the explicit record-30 a=0 completion.

Only the Python standard library is used.  The script reconstructs the
minimal-transversal facets, the full face table, local g/Macaulay and
mod-two homology-sphere-link tests, the two-class multiwedge quotient, and a
reverse-order order-five-Lanner join proving an unrestricted compact Coxeter
label obstruction.  It deliberately makes no convex-polytopality or
finite-volume noncompact claim.
"""

from __future__ import annotations

from collections import Counter, deque
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import comb
from typing import Iterable, Mapping, Sequence


N = 14
DIMENSION = 9
FULL = (1 << N) - 1
EDGES = ((0, 11), (1, 12), (2, 13))
A5S = (
    (0, 1, 4, 6, 7),
    (0, 4, 6, 7, 13),
    (1, 4, 5, 6, 7),
    (2, 3, 9, 10, 12),
    (3, 5, 8, 9, 10),
    (3, 8, 9, 10, 11),
    (3, 9, 10, 11, 12),
    (4, 5, 6, 7, 8),
)
EDGE_MASKS = tuple(sum(1 << v for v in edge) for edge in EDGES)
A5_MASKS = tuple(sum(1 << v for v in face) for face in A5S)
CLUTTER = (*EDGE_MASKS, *A5_MASKS)

LANNER_REPRESENTATIVES = (
    (((0, 1), 3), ((0, 2), 3), ((0, 3), 3), ((1, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 4), ((2, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 5), ((2, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 3), ((3, 4), 4)),
)

Pair = tuple[int, int]
State = tuple[tuple[Pair, int], ...]


def pair(a: int, b: int) -> Pair:
    return (a, b) if a < b else (b, a)


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int, n: int = N) -> tuple[int, ...]:
    return tuple(v for v in range(n) if value & (1 << v))


def canonical_sha(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def minimal_transversals(hyperedges: Sequence[int], n: int = N) -> tuple[int, ...]:
    result = []
    for candidate in range(1, 1 << n):
        if not all(candidate & edge for edge in hyperedges):
            continue
        if any(
            all((candidate ^ (1 << v)) & edge for edge in hyperedges)
            for v in vertices(candidate, n)
        ):
            continue
        result.append(candidate)
    return tuple(result)


def face_table(complements: Sequence[int], n: int = N) -> tuple[bool, ...]:
    full = (1 << n) - 1
    table = [False] * (1 << n)
    for complement in complements:
        facet = full ^ complement
        subset = facet
        while True:
            table[subset] = True
            if subset == 0:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def minimal_nonfaces(table: Sequence[bool], n: int = N) -> tuple[int, ...]:
    return tuple(
        candidate
        for candidate in range(1, 1 << n)
        if not table[candidate]
        and all(table[candidate ^ (1 << v)] for v in vertices(candidate, n))
    )


def ridge_groups(complements: Sequence[int]) -> dict[int, tuple[int, ...]]:
    groups: dict[int, list[int]] = {}
    for complement in complements:
        outside = FULL ^ complement
        while outside:
            bit = outside & -outside
            groups.setdefault(complement | bit, []).append(complement)
            outside ^= bit
    return {key: tuple(sorted(value)) for key, value in groups.items()}


def dual_components(complements: Sequence[int]) -> int:
    adjacency = {value: set() for value in complements}
    for group in ridge_groups(complements).values():
        if len(group) == 2:
            a, b = group
            adjacency[a].add(b)
            adjacency[b].add(a)
    unseen = set(complements)
    count = 0
    while unseen:
        count += 1
        root = min(unseen)
        unseen.remove(root)
        queue = deque([root])
        while queue:
            current = queue.popleft()
            for other in adjacency[current] & unseen:
                unseen.remove(other)
                queue.append(other)
    return count


def macaulay_raise(value: int, degree: int) -> int:
    remaining = value
    next_top: int | None = None
    raised = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (
            (next_top is None or top + 1 < next_top)
            and comb(top + 1, lower) <= remaining
        ):
            top += 1
        remaining -= comb(top, lower)
        raised += comb(top + 1, lower + 1)
        next_top = top
    if remaining:
        raise AssertionError("inexact Macaulay expansion")
    return raised


def link_g(table: Sequence[bool], face: int) -> tuple[int, tuple[int, ...]]:
    face_size = face.bit_count()
    polytope_dimension = DIMENSION - face_size
    rest = vertices(FULL ^ face)
    link_f = {-1: 1}
    for dimension in range(polytope_dimension):
        link_f[dimension] = sum(
            table[face | mask(extension)]
            for extension in combinations(rest, dimension + 1)
        )
    h = tuple(
        sum(
            (-1) ** (index - order)
            * comb(polytope_dimension - order, polytope_dimension - index)
            * link_f[order - 1]
            for order in range(index + 1)
        )
        for index in range(polytope_dimension // 2 + 1)
    )
    return link_f[0], tuple(h[i] - h[i - 1] for i in range(1, len(h)))


def local_g_macaulay_audit(table: Sequence[bool]) -> dict[str, object]:
    strict_orders = {1: (2, 3, 4), 2: (2, 3), 3: (2,)}
    tested = 0
    transitions = 0
    for size in (1, 2, 3):
        for chosen in combinations(range(N), size):
            face = mask(chosen)
            if any(face & edge == edge for edge in EDGE_MASKS):
                continue
            if not table[face]:
                raise AssertionError("an independent local base is absent")
            link_vertices, g = link_g(table, face)
            if link_vertices < DIMENSION - size + 2:
                raise AssertionError("a forbidden simplex link appeared")
            for order, value in enumerate(g, start=1):
                if value < (1 if order in strict_orders[size] else 0):
                    raise AssertionError("strict local g inequality failed")
            for index in range(1, len(g)):
                transitions += 1
                if g[index] > macaulay_raise(g[index - 1], index):
                    raise AssertionError("local Macaulay inequality failed")
            tested += 1
    return {"face_count": tested, "Macaulay_transition_count": transitions}


def euler_link_audit(table: Sequence[bool]) -> int:
    tested = 0
    for face, present in enumerate(table):
        if not present or face.bit_count() > 8:
            continue
        rest = FULL ^ face
        observed = 0
        subset = rest
        while subset:
            if table[face | subset]:
                observed += 1 if subset.bit_count() & 1 else -1
            subset = (subset - 1) & rest
        expected = 1 + (-1) ** (8 - face.bit_count())
        if observed != expected:
            raise AssertionError("face-link Euler characteristic failed")
        tested += 1
    return tested


def gf2_rank(columns: Iterable[int]) -> int:
    pivots: dict[int, int] = {}
    for raw in columns:
        column = raw
        while column:
            pivot = column.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = column
                break
            column ^= pivots[pivot]
    return len(pivots)


def link_faces(table: Sequence[bool], face: int) -> tuple[tuple[int, ...], ...]:
    groups: list[list[int]] = []
    subset = FULL ^ face
    while subset:
        if table[face | subset]:
            dimension = subset.bit_count() - 1
            while len(groups) <= dimension:
                groups.append([])
            groups[dimension].append(subset)
        subset = (subset - 1) & (FULL ^ face)
    return tuple(tuple(sorted(group)) for group in groups)


def reduced_betti(groups: Sequence[Sequence[int]]) -> tuple[int, ...]:
    if not groups:
        return ()
    ranks = [0] * (len(groups) + 1)
    for dimension in range(1, len(groups)):
        lower = {face: i for i, face in enumerate(groups[dimension - 1])}
        columns = []
        for face in groups[dimension]:
            column = 0
            bits = face
            while bits:
                bit = bits & -bits
                column |= 1 << lower[face ^ bit]
                bits ^= bit
            columns.append(column)
        ranks[dimension] = gf2_rank(columns)
    answer = [
        len(groups[d]) - ranks[d] - ranks[d + 1]
        for d in range(len(groups))
    ]
    answer[0] -= 1
    return tuple(answer)


def homology_link_audit(table: Sequence[bool]) -> int:
    tested = 0
    for face, present in enumerate(table):
        if not present or face.bit_count() > 9:
            continue
        expected_dimension = 8 - face.bit_count()
        groups = link_faces(table, face)
        if expected_dimension == -1:
            if groups:
                raise AssertionError("facet link is not void")
        else:
            observed = reduced_betti(groups)
            length = expected_dimension + 1
            padded = observed + (0,) * max(0, length - len(observed))
            expected = tuple(int(d == expected_dimension) for d in range(length))
            if padded[:length] != expected or any(padded[length:]):
                raise AssertionError("mod-two face-link homology failed")
        tested += 1
    return tested


def path_labels(
    component: frozenset[int], adjacency: Mapping[int, tuple[tuple[int, int], ...]]
) -> tuple[int, ...]:
    current = min(v for v in component if len(adjacency[v]) == 1)
    previous: int | None = None
    labels = []
    while True:
        following = [item for item in adjacency[current] if item[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("finite path branched")
        other, label = following[0]
        labels.append(label)
        previous, current = current, other
    result = tuple(labels)
    return min(result, tuple(reversed(result)))


def finite_component(component: frozenset[int], labels: Mapping[Pair, int]) -> bool:
    if len(component) == 1:
        return True
    adjacency = {
        v: tuple((w, labels[pair(v, w)]) for w in component if w != v and labels[pair(v, w)] > 2)
        for v in component
    }
    if sum(map(len, adjacency.values())) // 2 != len(component) - 1:
        return False
    degrees = Counter(map(len, adjacency.values()))
    if max(degrees) <= 2:
        sequence = path_labels(component, adjacency)
        order = len(component)
        return (
            (order == 2 and sequence[0] in (3, 4, 5))
            or all(label == 3 for label in sequence)
            or sequence == (3,) * (order - 2) + (4,)
            or (order == 3 and sequence == (3, 5))
            or (order == 4 and sequence in ((3, 3, 5), (3, 4, 3)))
        )
    if degrees[3] != 1 or any(degree > 3 for degree in degrees):
        return False
    if any(labels[pair(a, b)] != 3 for a, b in combinations(component, 2) if labels[pair(a, b)] > 2):
        return False
    branch = next(v for v in component if len(adjacency[v]) == 3)
    arms = []
    for neighbor, _ in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(v for v, _ in adjacency[current] if v != previous)
            previous, current = current, following
            length += 1
        arms.append(length)
    shape = tuple(sorted(arms))
    return shape == (1, 1, len(component) - 3) or shape in ((1, 2, 2), (1, 2, 3), (1, 2, 4))


def finite(facets: Iterable[int], labels: Mapping[Pair, int]) -> bool:
    unseen = set(facets)
    components = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        frontier = [start]
        while frontier:
            current = frontier.pop()
            neighbors = {v for v in unseen if labels[pair(current, v)] > 2}
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        components.append(frozenset(component))
    return all(finite_component(component, labels) for component in components)


def representative_code(representative: Sequence[tuple[Pair, int]]) -> tuple[int, ...]:
    edges = tuple(combinations(range(5), 2))
    labels = {edge: 2 for edge in edges}
    labels.update(representative)
    return tuple(labels[edge] for edge in edges)


def permuted_code(code: Sequence[int], permutation: Sequence[int]) -> tuple[int, ...]:
    edges = tuple(combinations(range(5), 2))
    source = dict(zip(edges, code, strict=True))
    moved = {pair(permutation[a], permutation[b]): label for (a, b), label in source.items()}
    return tuple(moved[edge] for edge in edges)


@lru_cache(maxsize=None)
def lanner_codes() -> tuple[tuple[int, ...], ...]:
    result = tuple(sorted({
        permuted_code(representative_code(rep), permutation)
        for rep in LANNER_REPRESENTATIVES
        for permutation in permutations(range(5))
    }))
    canonical = {min(permuted_code(code, p) for p in permutations(range(5))) for code in result}
    if (len(result), len(canonical)) != (420, 5):
        raise AssertionError("order-five Lanner table changed")
    return result


def options(support: Sequence[int]) -> tuple[State, ...]:
    edges = tuple(combinations(sorted(support), 2))
    return tuple(tuple(zip(edges, code, strict=True)) for code in lanner_codes())


def label_obstruction(actual_vertices: Sequence[Sequence[int]]) -> dict[str, object]:
    def run(order: tuple[int, ...]) -> tuple[list[dict[str, object]], int]:
        states: tuple[State, ...] = ((),)
        union: set[int] = set()
        stages = []
        for support_index in order:
            support = A5S[support_index]
            union.update(support)
            cluster = tuple(sorted(union))
            restrictions = tuple(sorted({
                tuple(v for v in actual if v in union) for actual in actual_vertices
            }, key=lambda row: (len(row), row)))
            passing: set[State] = set()
            compatible = 0
            for state in states:
                old = dict(state)
                for option in options(support):
                    if any(edge in old and old[edge] != label for edge, label in option):
                        continue
                    compatible += 1
                    merged = old | dict(option)
                    labels = {pair(a, b): 2 for a, b in combinations(cluster, 2)}
                    labels.update(merged)
                    if all(finite(restriction, labels) for restriction in restrictions):
                        passing.add(tuple(sorted(merged.items())))
            states = tuple(sorted(passing))
            stages.append({
                "support_index": support_index,
                "compatible_join_count": compatible,
                "locally_elliptic_state_count": len(states),
                "restriction_size_histogram": {
                    str(k): v for k, v in sorted(Counter(map(len, restrictions)).items())
                },
            })
        return stages, len(states)

    production_stages, production_terminal = run((0, 1, 2, 7))
    reverse_stages, reverse_terminal = run((7, 2, 1))
    if [row["locally_elliptic_state_count"] for row in production_stages] != [420, 96, 6, 0]:
        raise AssertionError("production-order Lanner obstruction changed")
    if [row["locally_elliptic_state_count"] for row in reverse_stages] != [420, 96, 0]:
        raise AssertionError("reverse-order Lanner obstruction changed")
    return {
        "production_order": {"support_order": [0, 1, 2, 7], "stages": production_stages,
                             "terminal_state_count": production_terminal},
        "independent_reverse_order": {"support_order": [7, 2, 1], "stages": reverse_stages,
                                      "terminal_state_count": reverse_terminal},
        "terminal_state_count": 0,
        "minimal_right_angle_baseline_fails": True,
        "Perron_Frobenius_monotonicity_promotes_to_unrestricted_labels_2_through_5": True,
        "unrestricted_compact_coxeter_no_go_for_this_completion": True,
    }


def multiwedge_audit(minimal_faces: Sequence[int]) -> dict[str, object]:
    clone_classes = ((3, 9, 10), (4, 6, 7))
    singles = (0, 1, 2, 5, 8, 11, 12, 13)
    blocks = tuple((v,) for v in singles) + clone_classes
    collapse = {v: index for index, block in enumerate(blocks) for v in block}
    base_nonfaces = tuple(sorted({mask({collapse[v] for v in vertices(face)}) for face in minimal_faces}))
    expanded = tuple(sorted(mask(v for index in vertices(face, 10) for v in blocks[index]) for face in base_nonfaces))
    if expanded != tuple(sorted(minimal_faces)):
        raise AssertionError("multiwedge expansion does not recover the clutter")
    base_table = tuple(
        not any(value & nonface == nonface for nonface in base_nonfaces)
        for value in range(1 << 10)
    )
    base_facets = tuple(value for value in range(1 << 10) if base_table[value] and not any(
        base_table[value | (1 << v)] for v in range(10) if not value & (1 << v)
    ))
    base_vector = tuple(sum(base_table[value] and value.bit_count() == size for value in range(1 << 10)) for size in range(1, 6))
    if Counter(map(int.bit_count, base_facets)) != {5: 36} or base_vector != (10, 42, 88, 90, 36):
        raise AssertionError("base multiwedge sphere census changed")
    return {
        "clone_classes": [list(block) for block in clone_classes],
        "base_vertex_blocks": [list(block) for block in blocks],
        "base_minimal_nonfaces": [list(vertices(face, 10)) for face in base_nonfaces],
        "base_f_vector": list(base_vector),
        "base_facet_count": len(base_facets),
        "exact_multiwedge_identity": True,
        "convex_polytopality_proved": False,
    }


def audit() -> dict[str, object]:
    blocker = minimal_transversals(CLUTTER)
    if Counter(map(int.bit_count, blocker)) != {5: 132}:
        raise AssertionError("candidate blocker is not pure 5^132")
    table = face_table(blocker)
    vector = tuple(sum(present and value.bit_count() == size for value, present in enumerate(table)) for size in range(1, 10))
    if vector != (14, 88, 328, 806, 1364, 1596, 1248, 594, 132):
        raise AssertionError("candidate face vector changed")
    recovered = minimal_nonfaces(table)
    if recovered != tuple(sorted(CLUTTER)):
        raise AssertionError("candidate minimal nonfaces changed")
    ridges = ridge_groups(blocker)
    if Counter(map(len, ridges.values())) != {2: 594} or dual_components(blocker) != 1:
        raise AssertionError("candidate pseudomanifold gate failed")
    actual_vertices = tuple(vertices(FULL ^ complement) for complement in blocker)
    local = local_g_macaulay_audit(table)
    euler_count = euler_link_audit(table)
    homology_count = homology_link_audit(table)
    payload = b"".join(value.to_bytes(2, "little") for value in blocker)
    result = {
        "kind": "independent_h9_14_record30_a0_candidate_audit",
        "compact_only": True,
        "missing_edges": [list(edge) for edge in EDGES],
        "missing_A5s": [list(face) for face in A5S],
        "minimal_transversal_size_histogram": {"5": 132},
        "facet_complement_sha256": sha256(payload).hexdigest(),
        "face_vector": list(vector),
        "ridge_degree_histogram": {"2": 594},
        "dual_component_count": 1,
        "minimal_nonface_size_histogram": {"2": 3, "5": 8},
        "local_g_Macaulay_audit": local,
        "Euler_link_count": euler_count,
        "mod2_homology_sphere_link_count": homology_count,
        "abstract_mod2_homology_sphere_gate_passes": True,
        "multiwedge": multiwedge_audit(recovered),
        "coxeter_label_obstruction": label_obstruction(actual_vertices),
        "convex_polytopality_proved": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "scope": (
            "One explicit abstract compact-nerve completion only. The exact "
            "unrestricted Coxeter-label no-go is compact-only; convex "
            "polytopality and finite-volume noncompact chambers are not decided."
        ),
    }
    result["certificate_sha256"] = canonical_sha(result)
    return result


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
