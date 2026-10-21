#!/usr/bin/env python3
"""No-project-import topology and Coxeter audit of record-18 ``a=0``.

The input is the exact-family archive produced by the independently written
balanced-cut census.  This program uses only the Python standard library.  It
reconstructs the facet complex, checks every face link over GF(2), regenerates
the rank-five Lanner placement table from the classical five representatives,
and joins the placements while testing every actual rank-nine vertex against
the finite Coxeter classification.

Only compact chambers are tested: their vertices must be elliptic.  No claim
about finite-volume noncompact chambers is made here.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from functools import lru_cache
from itertools import combinations, permutations
import json
from math import comb
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    ROOT / "data/enumerations/h9_14_record18_a0_independent_audit.json"
)
DEFAULT_OUTPUT = (
    ROOT / "data/enumerations/h9_14_record18_a0_downstream_independent_audit.json"
)
VERTEX_COUNT = 14
DIMENSION = 9
FULL_MASK = (1 << VERTEX_COUNT) - 1
EDGES = ((0, 12), (1, 12), (2, 13), (3, 13), (12, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in EDGES)
EXPECTED_F_VECTOR = (14, 86, 310, 735, 1204, 1372, 1052, 495, 110)

Pair = tuple[int, int]
LabelState = tuple[tuple[Pair, int], ...]


def mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(
        vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1
    )


def pair(left: int, right: int) -> Pair:
    return (left, right) if left < right else (right, left)


def hits_all(value: int, clutter: Sequence[int]) -> bool:
    return all(value & row for row in clutter)


def minimal_transversals(clutter: Sequence[int]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, FULL_MASK + 1):
        if not hits_all(candidate, clutter):
            continue
        if any(
            hits_all(candidate ^ (1 << vertex), clutter)
            for vertex in vertices(candidate)
        ):
            continue
        answers.append(candidate)
    return tuple(answers)


def face_table(complements: Sequence[int]) -> tuple[bool, ...]:
    table = [False] * (1 << VERTEX_COUNT)
    table[0] = True
    for complement in sorted(set(complements)):
        if complement.bit_count() != 5:
            raise AssertionError("facet complement is not a five-set")
        facet = FULL_MASK ^ complement
        subset = facet
        while True:
            table[subset] = True
            if not subset:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def face_vector(table: Sequence[bool]) -> tuple[int, ...]:
    return tuple(
        sum(
            table[mask(chosen)]
            for chosen in combinations(range(VERTEX_COUNT), order)
        )
        for order in range(1, DIMENSION + 1)
    )


def minimal_nonfaces(table: Sequence[bool]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, 1 << VERTEX_COUNT):
        if table[candidate]:
            continue
        if all(
            table[candidate ^ (1 << vertex)]
            for vertex in vertices(candidate)
        ):
            answers.append(candidate)
    return tuple(answers)


def ridge_groups(complements: Sequence[int]) -> dict[int, tuple[int, ...]]:
    groups: dict[int, list[int]] = {}
    for complement in sorted(set(complements)):
        facet = FULL_MASK ^ complement
        remaining = facet
        while remaining:
            bit = remaining & -remaining
            groups.setdefault(complement | bit, []).append(complement)
            remaining ^= bit
    return {
        ridge: tuple(sorted(group)) for ridge, group in groups.items()
    }


def ridge_histogram(complements: Sequence[int]) -> dict[int, int]:
    return dict(
        sorted(Counter(map(len, ridge_groups(complements).values())).items())
    )


def dual_component_count(complements: Sequence[int]) -> int:
    selected = tuple(sorted(set(complements)))
    adjacency = {complement: set() for complement in selected}
    for group in ridge_groups(selected).values():
        if len(group) != 2:
            continue
        first, second = group
        adjacency[first].add(second)
        adjacency[second].add(first)
    unseen = set(selected)
    components = 0
    while unseen:
        components += 1
        root = min(unseen)
        unseen.remove(root)
        queue = deque((root,))
        while queue:
            current = queue.popleft()
            for neighbor in adjacency[current] & unseen:
                unseen.remove(neighbor)
                queue.append(neighbor)
    return components


@lru_cache(maxsize=None)
def independent_masks(order: int) -> tuple[int, ...]:
    return tuple(
        candidate
        for chosen in combinations(range(VERTEX_COUNT), order)
        for candidate in (mask(chosen),)
        if all(candidate & edge != edge for edge in EDGE_MASKS)
    )


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
        raise AssertionError("Macaulay expansion failed")
    return raised


def link_g_entries(
    table: Sequence[bool], face: int
) -> tuple[int, tuple[int, ...]]:
    face_size = face.bit_count()
    polytope_dimension = DIMENSION - face_size
    remaining = vertices(FULL_MASK ^ face)
    link_f = {-1: 1}
    for link_dimension in range(polytope_dimension):
        order = link_dimension + 1
        link_f[link_dimension] = sum(
            table[face | mask(extension)]
            for extension in combinations(remaining, order)
        )
    h = tuple(
        sum(
            (-1) ** (index - face_order)
            * comb(
                polytope_dimension - face_order,
                polytope_dimension - index,
            )
            * link_f[face_order - 1]
            for face_order in range(index + 1)
        )
        for index in range(polytope_dimension // 2 + 1)
    )
    return link_f[0], tuple(
        h[index] - h[index - 1] for index in range(1, len(h))
    )


def first_strict_local_g_failure(
    table: Sequence[bool],
) -> tuple[int, int, int, int] | None:
    strict_orders = {1: (2, 3, 4), 2: (2, 3), 3: (2,)}
    for face_size, required in strict_orders.items():
        polytope_dimension = DIMENSION - face_size
        for face in independent_masks(face_size):
            if not table[face]:
                return face, 0, -1, 0
            link_vertices, entries = link_g_entries(table, face)
            if link_vertices < polytope_dimension + 2:
                return face, 0, link_vertices, polytope_dimension + 2
            for order, observed in enumerate(entries, start=1):
                minimum = int(order in required)
                if observed < minimum:
                    return face, order, observed, minimum
    return None


def first_local_macaulay_failure(
    table: Sequence[bool],
) -> tuple[int, int, int, int] | None:
    for face_size in (1, 2, 3):
        for face in independent_masks(face_size):
            _link_vertices, entries = link_g_entries(table, face)
            for previous_order in range(1, len(entries)):
                previous = entries[previous_order - 1]
                observed = entries[previous_order]
                if previous < 0:
                    return face, previous_order + 1, observed, -1
                upper = macaulay_raise(previous, previous_order)
                if observed > upper:
                    return face, previous_order + 1, observed, upper
    return None


def first_link_euler_failure(
    table: Sequence[bool],
) -> tuple[int, int, int] | None:
    for face, is_face in enumerate(table):
        if not is_face or face.bit_count() > DIMENSION - 1:
            continue
        remaining = FULL_MASK ^ face
        observed = 0
        subset = remaining
        while subset:
            if table[face | subset]:
                observed += 1 if subset.bit_count() % 2 else -1
            subset = (subset - 1) & remaining
        link_dimension = DIMENSION - 1 - face.bit_count()
        expected = 1 + (-1) ** link_dimension
        if observed != expected:
            return face, observed, expected
    return None


def gf2_rank(columns: Iterable[int]) -> int:
    pivots: dict[int, int] = {}
    for raw in columns:
        column = raw
        while column:
            pivot = column.bit_length() - 1
            previous = pivots.get(pivot)
            if previous is None:
                pivots[pivot] = column
                break
            column ^= previous
    return len(pivots)


def link_face_masks(
    table: Sequence[bool], face: int
) -> tuple[tuple[int, ...], ...]:
    remaining = FULL_MASK ^ face
    grouped: list[list[int]] = []
    subset = remaining
    while subset:
        if table[face | subset]:
            dimension = subset.bit_count() - 1
            while len(grouped) <= dimension:
                grouped.append([])
            grouped[dimension].append(subset)
        subset = (subset - 1) & remaining
    return tuple(tuple(sorted(group)) for group in grouped)


def reduced_betti_mod2(
    faces_by_dimension: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    if not faces_by_dimension:
        return ()
    normalized = tuple(tuple(sorted(set(group))) for group in faces_by_dimension)
    boundary_ranks = [0] * (len(normalized) + 1)
    for dimension in range(1, len(normalized)):
        lower_index = {
            face: index for index, face in enumerate(normalized[dimension - 1])
        }
        columns = []
        for face in normalized[dimension]:
            column = 0
            bits = face
            while bits:
                bit = bits & -bits
                boundary = face ^ bit
                if boundary not in lower_index:
                    raise AssertionError("link is not downward closed")
                column |= 1 << lower_index[boundary]
                bits ^= bit
            columns.append(column)
        boundary_ranks[dimension] = gf2_rank(columns)
    ordinary = tuple(
        len(normalized[dimension])
        - boundary_ranks[dimension]
        - boundary_ranks[dimension + 1]
        for dimension in range(len(normalized))
    )
    if ordinary[0] <= 0:
        raise AssertionError("nonempty link has no H_0")
    return (ordinary[0] - 1, *ordinary[1:])


def first_mod2_homology_failure(
    table: Sequence[bool],
) -> tuple[int, int, tuple[int, ...]] | None:
    for face, is_face in enumerate(table):
        if not is_face or face.bit_count() > DIMENSION:
            continue
        expected_dimension = DIMENSION - 1 - face.bit_count()
        link_faces = link_face_masks(table, face)
        if expected_dimension == -1:
            if link_faces:
                return face, expected_dimension, ()
            continue
        observed = reduced_betti_mod2(link_faces)
        expected_length = expected_dimension + 1
        padded = observed + (0,) * max(0, expected_length - len(observed))
        expected = tuple(
            int(degree == expected_dimension)
            for degree in range(expected_length)
        )
        if padded[:expected_length] != expected or any(padded[expected_length:]):
            return face, expected_dimension, observed
    return None


def topology_audit(fives: Sequence[int]) -> tuple[dict[str, object], tuple[int, ...]]:
    clutter = tuple(sorted((*EDGE_MASKS, *fives)))
    complements = minimal_transversals(clutter)
    blocker_sizes = Counter(value.bit_count() for value in complements)
    table = face_table(complements)
    vector = face_vector(table)
    recovered = minimal_nonfaces(table)
    ridges = ridge_histogram(complements)
    components = dual_component_count(complements)
    triples = all(table[value] for value in independent_masks(3))
    strict_failure = first_strict_local_g_failure(table)
    macaulay_failure = first_local_macaulay_failure(table)
    euler_failure = first_link_euler_failure(table)
    homology_failure = first_mod2_homology_failure(table)
    valid = (
        blocker_sizes == {5: 110}
        and vector == EXPECTED_F_VECTOR
        and set(recovered) == set(clutter)
        and ridges == {2: 495}
        and components == 1
        and triples
        and strict_failure is None
        and macaulay_failure is None
        and euler_failure is None
        and homology_failure is None
    )
    if not valid:
        raise AssertionError(
            "independent a0 topology gate failed: "
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
            "face_count_including_empty": sum(table),
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


LANNER_FIVE_REPRESENTATIVES = (
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
)


@lru_cache(maxsize=1)
def lanner_five_codes() -> tuple[tuple[int, ...], ...]:
    edges = tuple(combinations(range(5), 2))
    codes = set()
    for representative in LANNER_FIVE_REPRESENTATIVES:
        base = {edge: 2 for edge in edges}
        base.update(representative)
        for move in permutations(range(5)):
            labels = {
                pair(move[left], move[right]): label
                for (left, right), label in base.items()
            }
            codes.add(tuple(labels[edge] for edge in edges))
    if len(codes) != 420:
        raise AssertionError("rank-five Lanner placement count changed")
    return tuple(sorted(codes))


def path_labels(component: set[int], adjacency) -> tuple[int, ...]:
    endpoint = min(vertex for vertex in component if len(adjacency[vertex]) == 1)
    labels = []
    previous = None
    current = endpoint
    while True:
        following = [item for item in adjacency[current] if item[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("finite path unexpectedly branched")
        next_vertex, label = following[0]
        labels.append(label)
        previous, current = current, next_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def branch_arms(component: set[int], adjacency) -> tuple[int, ...]:
    branch = next(vertex for vertex in component if len(adjacency[vertex]) == 3)
    lengths = []
    for neighbor, _label in adjacency[branch]:
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


def finite_component(component: set[int], labels: Mapping[Pair, int]) -> bool:
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
    degree_counts = Counter(map(len, adjacency.values()))
    if max(degree_counts) <= 2:
        sequence = path_labels(component, adjacency)
        if order == 2:
            return sequence[0] in (3, 4, 5)
        if all(label == 3 for label in sequence):
            return True
        if sequence == (3,) * (order - 2) + (4,):
            return True
        return (order, sequence) in {
            (3, (3, 5)),
            (4, (3, 3, 5)),
            (4, (3, 4, 3)),
        }
    if degree_counts[3] != 1 or any(degree > 3 for degree in degree_counts):
        return False
    if any(
        labels[pair(left, right)] != 3
        for left in component
        for right in component
        if left < right and labels[pair(left, right)] > 2
    ):
        return False
    arms = branch_arms(component, adjacency)
    return arms == (1, 1, order - 3) or arms in {
        (1, 2, 2),
        (1, 2, 3),
        (1, 2, 4),
    }


@lru_cache(maxsize=262_144)
def finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    unseen = set(range(order))
    components = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        frontier = [start]
        while frontier:
            vertex = frontier.pop()
            neighbors = {
                other for other in unseen if labels[pair(vertex, other)] > 2
            }
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        components.append(component)
    return all(finite_component(component, labels) for component in components)


@lru_cache(maxsize=None)
def placement_options(support: tuple[int, ...]) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True)) for code in lanner_five_codes()
    )


def restriction_finite(
    restriction: Sequence[int], labels: Mapping[Pair, int]
) -> bool:
    return finite_code(
        len(restriction),
        tuple(
            labels[pair(left, right)]
            for left, right in combinations(restriction, 2)
        ),
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


def label_gate(
    fives: Sequence[int], complements: Sequence[int]
) -> dict[str, object]:
    supports = tuple(vertices(value) for value in fives)
    actual_vertices = tuple(
        tuple(vertex for vertex in range(VERTEX_COUNT) if not complement >> vertex & 1)
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
        passing: set[LabelState] = set()
        compatible = 0
        raw = len(states) * len(options)
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
                    restriction_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
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
        "rank_five_lanner_placement_count": len(lanner_five_codes()),
        "stages": stages,
        "terminal_state_count": len(states),
        "processed_support_count": len(stages),
        "unrestricted_compact_coxeter_no_go": not states,
        "monotonicity_basis": (
            "Unassigned labels are set to 2. In a Coxeter Gram matrix G=I-A, "
            "A is symmetric and entrywise nonnegative; increasing any label "
            "increases A and cannot decrease its Perron root. A failed finite "
            "restriction therefore cannot be repaired by later labels."
        ),
    }


def build_record(source: Path) -> dict[str, object]:
    document = json.loads(source.read_text(encoding="utf-8"))
    if not document.get("literal_family_sets_equal"):
        raise AssertionError("exact independent/production family equality is absent")
    representatives = document.get("canonical_representatives")
    if not isinstance(representatives, list) or len(representatives) != 3:
        raise AssertionError("expected three corrected a0 canonical classes")
    audits = []
    for representative in representatives:
        fives = tuple(map(int, representative["representative_missing_five_masks"]))
        topology, complements = topology_audit(fives)
        labels = label_gate(fives, complements)
        if not labels["unrestricted_compact_coxeter_no_go"]:
            raise AssertionError("a corrected a0 class survived the label gate")
        audits.append(
            {
                "canonical_family_key": representative["canonical_family_key"],
                "raw_family_multiplicity": representative["multiplicity"],
                "representative_missing_five_masks": list(fives),
                "representative_missing_fives": [
                    list(vertices(value)) for value in fives
                ],
                "topology": topology,
                "coxeter_labels": labels,
            }
        )
    return {
        "kind": "independent_h9_14_record18_a0_downstream_audit",
        "source": str(source),
        "compact_only": True,
        "finite_volume_noncompact_checked": False,
        "canonical_class_count": len(audits),
        "class_multiplicities": sorted(
            audit["raw_family_multiplicity"] for audit in audits
        ),
        "all_classes_pass_connected_mod2_homology_audit": all(
            audit["topology"]["passes_connected_homology_audit"]
            for audit in audits
        ),
        "all_classes_have_unrestricted_compact_coxeter_no_go": all(
            audit["coxeter_labels"]["unrestricted_compact_coxeter_no_go"]
            for audit in audits
        ),
        "classes": audits,
        "a0_compact_branch_excluded": True,
        "whole_shadow_compact_coxeter_theorem": False,
        "scope": (
            "Independent compact-only downstream audit of all corrected a=0 "
            "classes. This is not a finite-volume noncompact audit and does "
            "not by itself close the other record-18 branches."
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
