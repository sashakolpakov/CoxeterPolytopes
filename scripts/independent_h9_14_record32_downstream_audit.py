#!/usr/bin/env python3
"""No-project-import downstream hostile replay for H9/14 record 32.

This script independently dualizes all 1,025 a=0 target families, applies
the literal blocker and ridge gates, rebuilds the classical order-four and
order-five Lanner placement tables, and runs actual-vertex ellipticity gates.
It also reconstructs both ridge-valid a=1 targets and exhausts their twenty
single-m=3 cross-edge tests.  Only the Python standard library is imported.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
TARGETS = ROOT / "data/enumerations/h9_14_record32_a0_all_targets_scout.json"
SELECTED_A2 = ROOT / "data/enumerations/h9_14_record32_coxeter_labels.json"
TARGETS_SHA256 = "ae84f9c8dde0f66cd4d8fc8c0b9dee2afab985c7d484f2e052709ad3c3734938"
SELECTED_A2_SHA256 = "728d3bd5640b69e7584a061766db102c11e74cb117fd1a852ebdba8433094632"

N = 14
DIMENSION = 9
FULL_MASK = (1 << N) - 1
MISSING_EDGES = ((0, 11), (1, 12), (0, 13), (2, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in MISSING_EDGES)

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
    5: (
        (((0, 1), 3), ((0, 2), 3), ((0, 3), 3), ((1, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 4), ((2, 4), 5)),
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 5), ((2, 4), 5)),
        (
            ((0, 1), 3),
            ((0, 2), 3),
            ((1, 3), 3),
            ((2, 4), 3),
            ((3, 4), 4),
        ),
    ),
}

Pair = tuple[int, int]
LabelState = tuple[tuple[Pair, int], ...]


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pair(left: int, right: int) -> Pair:
    return (left, right) if left < right else (right, left)


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value & (1 << vertex))


def path_labels(
    component: frozenset[int],
    adjacency: Mapping[int, tuple[tuple[int, int], ...]],
) -> tuple[int, ...]:
    endpoint = min(vertex for vertex in component if len(adjacency[vertex]) == 1)
    labels: list[int] = []
    previous: int | None = None
    current = endpoint
    while True:
        following = [item for item in adjacency[current] if item[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("finite-path traversal branched")
        next_vertex, label = following[0]
        labels.append(label)
        previous, current = current, next_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def arm_lengths(
    component: frozenset[int],
    adjacency: Mapping[int, tuple[tuple[int, int], ...]],
) -> tuple[int, ...]:
    branch = next(vertex for vertex in component if len(adjacency[vertex]) == 3)
    lengths = []
    for neighbor, _ in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(
                vertex for vertex, _ in adjacency[current] if vertex != previous
            )
            previous, current = current, following
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def finite_component(
    component: frozenset[int], labels: Mapping[Pair, int]
) -> bool:
    order = len(component)
    if order == 1:
        return True
    adjacency = {
        vertex: tuple(
            (other, labels[pair(vertex, other)])
            for other in component
            if other != vertex and labels[pair(vertex, other)] > 2
        )
        for vertex in component
    }
    if sum(map(len, adjacency.values())) // 2 != order - 1:
        return False
    degrees = Counter(len(adjacency[vertex]) for vertex in component)
    if max(degrees) <= 2:
        sequence = path_labels(component, adjacency)
        if order == 2:
            return sequence[0] >= 3
        if all(label == 3 for label in sequence):
            return True
        if sequence == (3,) * (order - 2) + (4,):
            return True
        if order == 3 and sequence == (3, 5):
            return True
        if order == 4 and sequence in ((3, 3, 5), (3, 4, 3)):
            return True
        return False
    if degrees[3] != 1 or any(degree > 3 for degree in degrees):
        return False
    if any(
        labels[pair(left, right)] != 3
        for left, right in combinations(component, 2)
        if labels[pair(left, right)] > 2
    ):
        return False
    arms = arm_lengths(component, adjacency)
    return arms == (1, 1, order - 3) or arms in {
        (1, 2, 2),
        (1, 2, 3),
        (1, 2, 4),
    }


def finite(facets: Iterable[int], labels: Mapping[Pair, int]) -> bool:
    values = tuple(sorted(set(facets)))
    if any(pair(left, right) not in labels for left, right in combinations(values, 2)):
        return False
    unseen = set(values)
    components = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        frontier = [start]
        while frontier:
            vertex = frontier.pop()
            neighbors = {
                other
                for other in unseen
                if labels[pair(vertex, other)] > 2
            }
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        components.append(frozenset(component))
    return all(finite_component(component, labels) for component in components)


def representative_code(
    order: int, representative: Sequence[tuple[Pair, int]]
) -> tuple[int, ...]:
    edges = tuple(combinations(range(order), 2))
    labels = {edge: 2 for edge in edges}
    labels.update(representative)
    return tuple(labels[edge] for edge in edges)


def permuted_code(
    order: int, code: Sequence[int], permutation: Sequence[int]
) -> tuple[int, ...]:
    edges = tuple(combinations(range(order), 2))
    source = dict(zip(edges, code, strict=True))
    moved = {
        pair(permutation[left], permutation[right]): label
        for (left, right), label in source.items()
    }
    return tuple(moved[edge] for edge in edges)


@lru_cache(maxsize=None)
def lanner_codes(order: int) -> tuple[tuple[int, ...], ...]:
    codes = tuple(
        sorted(
            {
                permuted_code(
                    order,
                    representative_code(order, representative),
                    permutation,
                )
                for representative in LANNER_REPRESENTATIVES[order]
                for permutation in permutations(range(order))
            }
        )
    )
    expected = {4: (108, 9), 5: (420, 5)}[order]
    canonical = {
        min(permuted_code(order, code, p) for p in permutations(range(order)))
        for code in codes
    }
    if (len(codes), len(canonical)) != expected:
        raise AssertionError("independent Lanner table census changed")
    local_edges = tuple(combinations(range(order), 2))
    for code in codes:
        labels = dict(zip(local_edges, code, strict=True))
        if finite(range(order), labels):
            raise AssertionError("Lanner row classified finite")
        for deleted in range(order):
            if not finite((v for v in range(order) if v != deleted), labels):
                raise AssertionError("Lanner row is not minimal")
    return codes


@lru_cache(maxsize=None)
def placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True)) for code in lanner_codes(len(support))
    )


@lru_cache(maxsize=524_288)
def finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return finite(range(order), labels)


def restriction_finite(restriction: Sequence[int], labels: Mapping[Pair, int]) -> bool:
    return finite_code(
        len(restriction),
        tuple(labels[pair(left, right)] for left, right in combinations(restriction, 2)),
    )


def minimal_transversals(hyperedges: Sequence[int]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, FULL_MASK + 1):
        if not all(candidate & edge for edge in hyperedges):
            continue
        if any(
            all((candidate ^ (1 << vertex)) & edge for edge in hyperedges)
            for vertex in vertices(candidate)
        ):
            continue
        answers.append(candidate)
    return tuple(answers)


def ridge_histogram(facet_complements: Sequence[int]) -> dict[str, int]:
    ridges: Counter[int] = Counter()
    for complement in facet_complements:
        facet = FULL_MASK ^ complement
        for removed in vertices(facet):
            ridges[facet ^ (1 << removed)] += 1
    return {str(degree): count for degree, count in sorted(Counter(ridges.values()).items())}


def rank_nine_vertices(minimal_faces: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(frozenset(vertices(face)) for face in (*EDGE_MASKS, *minimal_faces))
    return tuple(
        chosen
        for chosen in combinations(range(N), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


def cluster_restrictions(
    actual_vertices: Sequence[Sequence[int]], cluster: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    cluster_set = set(cluster)
    return tuple(
        sorted(
            {
                tuple(vertex for vertex in actual if vertex in cluster_set)
                for actual in actual_vertices
            },
            key=lambda row: (len(row), row),
        )
    )


def hostile_next_support(
    supports: Sequence[tuple[int, ...]], remaining: set[int], union: set[int]
) -> int:
    if not union:
        return min(remaining)
    # Reverse the producer's final index tie-break deliberately.
    return max(
        remaining,
        key=lambda index: (
            len(set(supports[index]) & union),
            -len(set(supports[index]) - union),
            index,
        ),
    )


def sparse_states(
    minimal_faces: Sequence[int],
) -> tuple[tuple[LabelState, ...], tuple[tuple[int, ...], ...], list[list[int]]]:
    supports = tuple(vertices(face) for face in minimal_faces)
    actual_vertices = rank_nine_vertices(minimal_faces)
    remaining = set(range(len(supports)))
    union: set[int] = set()
    states: tuple[LabelState, ...] = ((),)
    stages: list[list[int]] = []
    while states and remaining:
        support_index = hostile_next_support(supports, remaining, union)
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = cluster_restrictions(actual_vertices, cluster)
        passing: set[LabelState] = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in placement_options(support):
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(restriction_finite(row, labels) for row in restrictions):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append([support_index, compatible, len(states)])
        remaining.remove(support_index)
    return states, actual_vertices, stages


def audit_a0() -> dict[str, object]:
    if file_sha256(TARGETS) != TARGETS_SHA256:
        raise AssertionError("a0 target archive digest changed")
    document = json.loads(TARGETS.read_text(encoding="utf-8"))
    families = tuple(tuple(map(int, row)) for row in document["target_family_masks"])
    if len(families) != 1025:
        raise AssertionError("a0 target-family count changed")

    pure = 0
    impure: Counter[str] = Counter()
    ridge_indices = []
    label_rows = []
    for index, family in enumerate(families):
        blocker = minimal_transversals((*EDGE_MASKS, *family))
        histogram = Counter(value.bit_count() for value in blocker)
        if histogram == {5: 114}:
            pure += 1
        else:
            if histogram.get(5) != 114 or set(histogram) != {5, 6}:
                raise AssertionError(f"unexpected blocker at a0 family {index}")
            impure[f"6:{histogram[6]}"] += 1
        fives = tuple(value for value in blocker if value.bit_count() == 5)
        ridge = ridge_histogram(fives)
        if histogram == {5: 114} and ridge == {"2": 513}:
            ridge_indices.append(index)
            states, actual_vertices, stages = sparse_states(family)
            label_rows.append(
                {
                    "family_index": index,
                    "actual_vertex_count": len(actual_vertices),
                    "stage_support_compatible_passing": stages,
                    "terminal_state_count": len(states),
                    "unrestricted_compact_label_no_go": not states,
                }
            )
    expected_indices = [6, 27, 43, 57, 132, 139, 436, 439, 481, 574, 626, 629]
    if pure != 941 or ridge_indices != expected_indices:
        raise AssertionError("a0 blocker/ridge census changed")
    if Counter(impure) != Counter({"6:1": 3, "6:2": 30, "6:3": 8,
                                   "6:4": 20, "6:8": 19, "6:16": 4}):
        raise AssertionError("a0 impure-blocker histogram changed")
    if any(row["terminal_state_count"] for row in label_rows):
        raise AssertionError("an a0 ridge family retained a sparse label state")
    return {
        "target_family_count": len(families),
        "pure_blocker_count": pure,
        "impure_blocker_histogram": dict(sorted(impure.items())),
        "ridge_valid_indices": ridge_indices,
        "ridge_histogram": {"2": 513},
        "label_audits": label_rows,
    }


A1_TARGETS = (
    (
        (3, 4, 5, 6),
        (2, 7, 8, 9, 10),
        (7, 8, 9, 10, 11),
    ),
    (
        (2, 3, 4, 5),
        (6, 7, 8, 9, 10),
        (3, 4, 5, 6, 11),
        (7, 8, 9, 10, 13),
    ),
)


def audit_a1_target(faces: Sequence[Sequence[int]]) -> dict[str, object]:
    minimal_faces = tuple(mask(face) for face in faces)
    blocker = minimal_transversals((*EDGE_MASKS, *minimal_faces))
    blocker_hist = Counter(value.bit_count() for value in blocker)
    ridge = ridge_histogram(blocker)
    if blocker_hist != {5: 112} or ridge != {"2": 504}:
        raise AssertionError("a1 target failed independent blocker/ridge replay")

    states, actual_vertices, stages = sparse_states(minimal_faces)
    left, right = tuple(faces[0]), tuple(faces[1])
    if set(left) & set(right):
        raise AssertionError("a1 cross supports are not disjoint")
    cross_pairs = tuple(pair(a, b) for a in left for b in right)
    if len(cross_pairs) != 20 or any(mask(edge) in EDGE_MASKS for edge in cross_pairs):
        raise AssertionError("a1 cross-pair scope changed")
    cluster = tuple(sorted(set().union(*map(set, faces))))
    restrictions = cluster_restrictions(actual_vertices, cluster)
    restrictions_by_edge = {
        edge: tuple(row for row in restrictions if set(edge) <= set(row))
        for edge in cross_pairs
    }
    template = {
        pair(left_vertex, right_vertex): 2
        for left_vertex, right_vertex in combinations(cluster, 2)
    }
    passing = 0
    bad_sizes: Counter[int] = Counter()
    for state in states:
        base = dict(template)
        base.update(dict(state))
        for edge in cross_pairs:
            labels = dict(base)
            labels[edge] = 3
            bad = next(
                (
                    row
                    for row in restrictions_by_edge[edge]
                    if not restriction_finite(row, labels)
                ),
                None,
            )
            if bad is None:
                passing += 1
            else:
                bad_sizes[len(bad)] += 1
    return {
        "minimal_faces": [list(face) for face in faces],
        "minimal_transversal_size_histogram": {"5": 112},
        "ridge_histogram": ridge,
        "actual_vertex_count": len(actual_vertices),
        "sparse_join_stages": stages,
        "sparse_state_count": len(states),
        "cross_pair_count": len(cross_pairs),
        "single_m3_test_count": len(states) * len(cross_pairs),
        "single_m3_passing_count": passing,
        "first_bad_restriction_size_histogram": {
            str(size): count for size, count in sorted(bad_sizes.items())
        },
        "unrestricted_compact_label_no_go": passing == 0,
    }


def audit_a1() -> dict[str, object]:
    rows = [audit_a1_target(faces) for faces in A1_TARGETS]
    expected = (
        (10368, 207360, {"7": 206208, "8": 1152}),
        (864, 17280, {"7": 17280}),
    )
    for row, (states, tests, histogram) in zip(rows, expected):
        if (
            row["sparse_state_count"] != states
            or row["single_m3_test_count"] != tests
            or row["single_m3_passing_count"] != 0
            or row["first_bad_restriction_size_histogram"] != histogram
        ):
            raise AssertionError("a1 hostile cross-gate transcript changed")
    return {"target_audits": rows}


def audit_selected_a2_dependency() -> dict[str, object]:
    if file_sha256(SELECTED_A2) != SELECTED_A2_SHA256:
        raise AssertionError("selected-a2 compact-label certificate digest changed")
    document = json.loads(SELECTED_A2.read_text(encoding="utf-8"))
    theorem = document["theorem"]
    faces = theorem["minimal_nonfaces"]
    expected_fours = [[2, 3, 4, 5], [3, 4, 5, 11]]
    expected_five = [6, 7, 8, 9, 10]
    if (
        faces["order_four_lanner_supports"] != expected_fours
        or faces["order_five_lanner_support"] != expected_five
        or not theorem["unrestricted_compact_coxeter_no_go"]
        or not theorem["verified"]
    ):
        raise AssertionError("selected-a2 theorem identity changed")
    return {
        "certificate_path": str(SELECTED_A2.relative_to(ROOT)),
        "file_sha256": SELECTED_A2_SHA256,
        "certificate_sha256": document["certificate_sha256"],
        "order_four_minimal_faces": expected_fours,
        "order_five_minimal_face": expected_five,
        "unrestricted_compact_coxeter_no_go": True,
    }


def main() -> None:
    # Build and self-check both classical tables even if a later gate exits early.
    if (len(lanner_codes(4)), len(lanner_codes(5))) != (108, 420):
        raise AssertionError("classical Lanner placement counts changed")
    output = {
        "kind": "independent_h9_14_record32_downstream_hostile_audit",
        "classical_lanner_placement_counts": {"4": 108, "5": 420},
        "a0": audit_a0(),
        "a1": audit_a1(),
        "selected_a2_dependency": audit_selected_a2_dependency(),
        "scope": {
            "compact_only": True,
            "ridge_rejections_are_only_pseudomanifold_rejections": True,
            "abstract_target_polytopality_claimed": False,
            "abstract_target_PL_sphere_status_claimed": False,
            "finite_volume_noncompact_inference": False,
        },
        "status": "pass",
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
