#!/usr/bin/env python3
"""No-import compact-Coxeter label probe for corrected record-12 orbits.

This script is deliberately independent of all project combinatorics,
Coxeter, Lanner, and reconstruction modules.  It reads only the frozen
371-orbit archive, reconstructs each nerve from its two 5x5 matrices, closes
the classical rank-4/rank-5 Lanner tables under all permutations, and tests
actual rank-nine faces with a separate finite-Coxeter classifier.

The default frozen claim concerns all 371 corrected incidence orbits.  The
``--a0-only`` option provides a faster, separately frozen 48-orbit replay.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/enumerations/h9_14_record12_residual.json"
INPUT_FILE_SHA256 = "faca367738bba58afc0ec78a43e9d4aa4b02721802c5d9e77b7d155c8e7faa05"
INPUT_CERTIFICATE_SHA256 = "50ca012334e7a4973a4c29bde4226d4460bc8ade2d644a3b9dbf56b1bf92e8f4"
EXPECTED_A0_SUMMARY_SHA256 = "efb4d61db80f84505cbf56e93b6f9ed47c2358061680b6c4ae17e9c50f3dd22a"
EXPECTED_ALL_SUMMARY_SHA256 = "06a0d502d2f739190762f5b7e2e6442a2cadb35a4c4a96577ee3885d9d822258"

FACET_COUNT = 14
DIMENSION = 9
FULL14 = (1 << FACET_COUNT) - 1
LOCAL_ORDER = 10
BALANCED_CUT = sum(1 << vertex for vertex in range(5, 10))
PAIR_INDEX = {
    pair: index for index, pair in enumerate(combinations(range(LOCAL_ORDER), 2))
}
GRAPH_EDGES = {
    (1 << 0) | (1 << 12),
    (1 << 1) | (1 << 12),
    (1 << 0) | (1 << 13),
    (1 << 1) | (1 << 13),
}

# Unlisted pairs have label 2.  These are the complete classical rank-4 and
# rank-5 Lanner representatives used only as seeds for permutation closure.
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
State = tuple[tuple[Pair, int], ...]


def pair(left: int, right: int) -> Pair:
    return (left, right) if left < right else (right, left)


def canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(FACET_COUNT) if value & (1 << vertex))


def pair_mask_from_matrix(matrix: int) -> int:
    return sum(
        1 << PAIR_INDEX[row, column + 5]
        for row in range(5)
        for column in range(5)
        if matrix & (1 << (5 * row + column))
    )


def graph_edge(value: int, first: int, second: int) -> bool:
    return bool(value & (1 << PAIR_INDEX[pair(first, second)]))


def cut_edge(first: int, second: int) -> bool:
    return bool(BALANCED_CUT & (1 << first)) != bool(
        BALANCED_CUT & (1 << second)
    )


def branch_triples(value: int) -> tuple[int, ...]:
    result = set()
    for first, second in combinations(range(LOCAL_ORDER), 2):
        edge = graph_edge(value, first, second)
        delta = cut_edge(first, second)
        if edge:
            result.add((1 << 11) | (1 << first) | (1 << second))
        if edge != delta:
            result.add((1 << 10) | (1 << first) | (1 << second))
    for triple in combinations(range(LOCAL_ORDER), 3):
        if sum(graph_edge(value, *edge) for edge in combinations(triple, 2)) % 2:
            result.add(mask(triple))
    return tuple(sorted(result))


def selected_complements(x_matrix: int, z_matrix: int) -> tuple[int, ...]:
    result = set()
    for fixed, extra, matrix in (
        ((0, 1), (12, 13), x_matrix),
        ((12, 13), (0, 1), z_matrix),
    ):
        for triple in branch_triples(pair_mask_from_matrix(matrix)):
            moved = 0
            for local in range(LOCAL_ORDER):
                if triple & (1 << local):
                    moved |= 1 << (local + 2)
            if triple & (1 << 10):
                moved |= 1 << extra[0]
            if triple & (1 << 11):
                moved |= 1 << extra[1]
            result.add(mask(fixed) | moved)
    return tuple(sorted(result))


def face_table(complements: Sequence[int]) -> tuple[bool, ...]:
    table = [False] * (1 << FACET_COUNT)
    for complement in complements:
        facet = FULL14 ^ complement
        subset = facet
        while True:
            table[subset] = True
            if not subset:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def minimal_nonfaces(table: Sequence[bool]) -> dict[int, tuple[int, ...]]:
    result: dict[int, list[int]] = {}
    for value in range(1, 1 << FACET_COUNT):
        if table[value]:
            continue
        if all(
            table[value ^ (1 << vertex)]
            for vertex in range(FACET_COUNT)
            if value & (1 << vertex)
        ):
            result.setdefault(value.bit_count(), []).append(value)
    return {order: tuple(rows) for order, rows in sorted(result.items())}


def rank_nine_faces(
    complements: Sequence[int], nonfaces: Mapping[int, Sequence[int]]
) -> tuple[tuple[int, ...], ...]:
    direct = tuple(sorted(vertices(FULL14 ^ complement) for complement in complements))
    forbidden = tuple(
        frozenset(vertices(value)) for rows in nonfaces.values() for value in rows
    )
    enumerated = tuple(
        chosen
        for chosen in combinations(range(FACET_COUNT), DIMENSION)
        if not any(row <= set(chosen) for row in forbidden)
    )
    if direct != enumerated:
        raise AssertionError("actual rank-nine faces disagree in two reconstructions")
    return direct


def path_labels(adjacency: Mapping[int, tuple[tuple[int, int], ...]]) -> tuple[int, ...]:
    endpoint = min(vertex for vertex, rows in adjacency.items() if len(rows) == 1)
    labels = []
    previous = None
    current = endpoint
    while True:
        following = [row for row in adjacency[current] if row[0] != previous]
        if not following:
            break
        if len(following) != 1:
            return ()
        following_vertex, label = following[0]
        labels.append(label)
        previous, current = current, following_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def arm_lengths(adjacency: Mapping[int, tuple[tuple[int, int], ...]]) -> tuple[int, ...]:
    branch = next(vertex for vertex, rows in adjacency.items() if len(rows) == 3)
    lengths = []
    for neighbor, _label in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(v for v, _ in adjacency[current] if v != previous)
            previous, current = current, following
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


@lru_cache(maxsize=None)
def finite_code(order: int, code: tuple[int, ...]) -> bool:
    edges = tuple(combinations(range(order), 2))
    if len(code) != len(edges):
        raise ValueError("finite code has the wrong length")
    labels = dict(zip(edges, code, strict=True))
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
                other
                for other in unseen
                if labels[pair(vertex, other)] > 2
            }
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        components.append(frozenset(component))

    for component in components:
        size = len(component)
        if size == 1:
            continue
        adjacency = {
            vertex: tuple(
                (other, labels[pair(vertex, other)])
                for other in component
                if other != vertex and labels[pair(vertex, other)] > 2
            )
            for vertex in component
        }
        if sum(map(len, adjacency.values())) // 2 != size - 1:
            return False
        degree_histogram = Counter(map(len, adjacency.values()))
        if max(degree_histogram) <= 2:
            sequence = path_labels(adjacency)
            if size == 2:
                if sequence[0] not in (3, 4, 5):
                    return False
            elif all(label == 3 for label in sequence):
                pass
            elif sequence == (3,) * (size - 2) + (4,):
                pass
            elif size == 3 and sequence == (3, 5):
                pass
            elif size == 4 and sequence in ((3, 3, 5), (3, 4, 3)):
                pass
            else:
                return False
            continue
        if degree_histogram[3] != 1 or any(degree > 3 for degree in degree_histogram):
            return False
        if any(
            labels[pair(left, right)] != 3
            for left, right in combinations(component, 2)
            if labels[pair(left, right)] > 2
        ):
            return False
        arms = arm_lengths(adjacency)
        if not (
            arms == (1, 1, size - 3)
            or arms in {(1, 2, 2), (1, 2, 3), (1, 2, 4)}
        ):
            return False
    return True


def representative_code(
    order: int, representative: Sequence[tuple[Pair, int]]
) -> tuple[int, ...]:
    labels = {edge: 2 for edge in combinations(range(order), 2)}
    labels.update(representative)
    return tuple(labels[edge] for edge in combinations(range(order), 2))


def permuted_code(order: int, code: Sequence[int], permutation: Sequence[int]) -> tuple[int, ...]:
    edges = tuple(combinations(range(order), 2))
    source = dict(zip(edges, code, strict=True))
    moved = {
        pair(permutation[left], permutation[right]): label
        for (left, right), label in source.items()
    }
    return tuple(moved[edge] for edge in edges)


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
        min(permuted_code(order, code, permutation) for permutation in permutations(range(order)))
        for code in codes
    }
    if (len(codes), len(canonical)) != expected:
        raise AssertionError("independent Lanner table count changed")
    for code in codes:
        if finite_code(order, code):
            raise AssertionError("Lanner placement is finite")
        labels = dict(zip(combinations(range(order), 2), code, strict=True))
        for deleted in range(order):
            retained = tuple(vertex for vertex in range(order) if vertex != deleted)
            subcode = tuple(labels[pair(*edge)] for edge in combinations(retained, 2))
            if not finite_code(order - 1, subcode):
                raise AssertionError("Lanner placement is not minimally nonfinite")
    return codes


def options(support: tuple[int, ...], tables: Mapping[int, Sequence[Sequence[int]]]) -> tuple[State, ...]:
    edges = tuple(combinations(support, 2))
    return tuple(tuple(zip(edges, code, strict=True)) for code in tables[len(support)])


def state_digest(states: Sequence[State]) -> str:
    payload = [
        [[left, right, label] for (left, right), label in state]
        for state in sorted(states)
    ]
    return sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()


def support_order(supports: Sequence[tuple[int, ...]]) -> tuple[int, ...]:
    """Deterministic touching-first order; it never removes a support.

    Start with the lexicographically first shortest support.  Thereafter take
    the first remaining support meeting the accumulated vertex union,
    deferring a disjoint component until necessary.  If the state set reaches
    zero, the unused supports need not be imposed for a no-go certificate.
    """

    remaining = set(range(len(supports)))
    first = min(remaining, key=lambda index: (len(supports[index]), supports[index]))
    order = [first]
    remaining.remove(first)
    union = set(supports[first])
    while remaining:
        touching = [index for index in remaining if union & set(supports[index])]
        chosen = min(
            touching or remaining,
            key=lambda index: (len(supports[index]), supports[index]),
        )
        order.append(chosen)
        remaining.remove(chosen)
        union.update(supports[chosen])
    return tuple(order)


def extend(
    states: Sequence[State],
    assigned_pairs: frozenset[Pair],
    union: frozenset[int],
    support_index: int,
    support: tuple[int, ...],
    actual_faces: Sequence[tuple[int, ...]],
    dotted_pairs: Sequence[tuple[int, ...]],
    tables: Mapping[int, Sequence[Sequence[int]]],
) -> tuple[tuple[State, ...], frozenset[Pair], frozenset[int], dict[str, object]]:
    new_union = union | frozenset(support)
    restrictions = tuple(
        sorted(
            {
                tuple(vertex for vertex in face if vertex in new_union)
                for face in actual_faces
            },
            key=lambda row: (len(row), row),
        )
    )
    if any(
        set(dotted) <= set(restriction)
        for dotted in dotted_pairs
        for restriction in restrictions
    ):
        raise AssertionError("an actual rank-nine restriction contains a dotted pair")

    placements = options(support, tables)
    support_pairs = frozenset(combinations(support, 2))
    overlap = tuple(sorted(assigned_pairs & support_pairs))
    buckets: dict[tuple[int, ...], list[State]] = defaultdict(list)
    for placement in placements:
        row = dict(placement)
        buckets[tuple(row[edge] for edge in overlap)].append(placement)

    passing = set()
    compatible = 0
    for state in states:
        old = dict(state)
        key = tuple(old[edge] for edge in overlap)
        for placement in buckets.get(key, ()):
            compatible += 1
            merged = old | dict(placement)
            finite = True
            for restriction in restrictions:
                code = tuple(
                    merged.get(pair(left, right), 2)
                    for left, right in combinations(restriction, 2)
                )
                if not finite_code(len(restriction), code):
                    finite = False
                    break
            if finite:
                passing.add(tuple(sorted(merged.items())))
    passing_states = tuple(sorted(passing))
    record = {
        "support_index": support_index,
        "support": list(support),
        "support_order": len(support),
        "complete_lanner_placement_count": len(placements),
        "overlap_pair_count": len(overlap),
        "cumulative_cluster_union": sorted(new_union),
        "actual_vertex_restriction_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(len, restrictions)).items())
        },
        "raw_extension_count": len(states) * len(placements),
        "compatible_lanner_extension_count": compatible,
        "elliptic_actual_vertex_survivor_count": len(passing_states),
        "passing_state_sha256": state_digest(passing_states),
    }
    return passing_states, assigned_pairs | support_pairs, new_union, record


def audit_orbit(
    index: int,
    orbit: Mapping[str, object],
    tables: Mapping[int, Sequence[Sequence[int]]],
) -> dict[str, object]:
    complements = selected_complements(int(orbit["x_matrix"]), int(orbit["z_matrix"]))
    digest = sha256(
        b"".join(value.to_bytes(2, "little") for value in complements)
    ).hexdigest()
    if digest != orbit["selected_complement_sha256"]:
        raise AssertionError("selected-complement digest changed")
    table = face_table(complements)
    nonfaces = minimal_nonfaces(table)
    observed_histogram = {str(order): len(rows) for order, rows in nonfaces.items()}
    if observed_histogram != orbit["minimal_nonface_histogram"]:
        raise AssertionError("minimal-nonface histogram changed")
    if set(nonfaces.get(2, ())) != GRAPH_EDGES:
        raise AssertionError("rank-two minimal nonfaces changed")
    actual_faces = rank_nine_faces(complements, nonfaces)
    if len(actual_faces) != int(orbit["facet_count"]):
        raise AssertionError("actual rank-nine face count changed")
    higher = tuple(
        sorted(
            (
                vertices(value)
                for order in (4, 5)
                for value in nonfaces.get(order, ())
            ),
            key=lambda row: (len(row), row),
        )
    )
    dotted = tuple(vertices(value) for value in nonfaces[2])
    order = support_order(higher)
    states: tuple[State, ...] = ((),)
    assigned_pairs: frozenset[Pair] = frozenset()
    union: frozenset[int] = frozenset()
    stages = []
    for support_index in order:
        states, assigned_pairs, union, stage = extend(
            states,
            assigned_pairs,
            union,
            support_index,
            higher[support_index],
            actual_faces,
            dotted,
            tables,
        )
        stages.append(stage)
        if not states:
            break
    return {
        "orbit_index": index,
        "x_matrix": int(orbit["x_matrix"]),
        "z_matrix": int(orbit["z_matrix"]),
        "missing_four_count": int(orbit["missing_four_count"]),
        "selected_complement_sha256": digest,
        "minimal_nonface_supports": [list(row) for row in higher],
        "actual_rank_nine_face_count": len(actual_faces),
        "deterministic_support_order": list(order),
        "stages": stages,
        "compact_coxeter_no_go": not states,
        "all_supports_processed_if_states_remain": not states or len(stages) == len(higher),
    }


@lru_cache(maxsize=2)
def replay(all_schedules: bool = True) -> dict[str, object]:
    raw = INPUT.read_bytes()
    if sha256(raw).hexdigest() != INPUT_FILE_SHA256:
        raise ValueError("corrected record-12 input file hash changed")
    document = json.loads(raw)
    if (
        document.get("certificate_sha256") != INPUT_CERTIFICATE_SHA256
        or canonical_digest(document) != INPUT_CERTIFICATE_SHA256
        or document.get("schema_version") != "2.0.0"
        or len(document.get("residual_orbits", ())) != 371
        or document.get("schedule_orbit_counts")
        != {"0": 48, "1": 95, "2": 144, "3": 63, "4": 21}
    ):
        raise ValueError("corrected 371-orbit certificate closure failed")
    tables = {4: lanner_codes(4), 5: lanner_codes(5)}
    selected = [
        (index, orbit)
        for index, orbit in enumerate(document["residual_orbits"])
        if all_schedules or int(orbit["missing_four_count"]) == 0
    ]
    orbits = [audit_orbit(index, orbit, tables) for index, orbit in selected]
    result: dict[str, object] = {
        "kind": (
            "h9_14_record12_all_schedule_compact_coxeter_probe"
            if all_schedules
            else "h9_14_record12_a0_compact_coxeter_no_go"
        ),
        "scope": {
            "compact_only": True,
            "burcroff_lemma_10_10_rank_four_and_five_lanners": True,
            "perron_frobenius_monotone_sparse_ellipticity_pruning": True,
            "finite_volume_noncompact_excluded": False,
            "global_lorentzian_gram_search_performed": False,
        },
        "method": (
            "standard-library only; corrected matrix reconstruction, two-way "
            "actual rank-nine face equality, independent 9+5 Lanner closure, "
            "independent finite Coxeter classifier, no label-state symmetry "
            "quotient, and deterministic touching-first support joins"
        ),
        "input_file_sha256": INPUT_FILE_SHA256,
        "input_certificate_sha256": INPUT_CERTIFICATE_SHA256,
        "lanner_tables": {
            "rank_four_placements": len(tables[4]),
            "rank_four_isomorphism_types": len(LANNER_REPRESENTATIVES[4]),
            "rank_five_placements": len(tables[5]),
            "rank_five_isomorphism_types": len(LANNER_REPRESENTATIVES[5]),
        },
        "orbit_count": len(orbits),
        "schedule_histogram": {
            str(value): count
            for value, count in sorted(
                Counter(row["missing_four_count"] for row in orbits).items()
            )
        },
        "orbits": orbits,
        "all_selected_orbits_excluded": all(row["compact_coxeter_no_go"] for row in orbits),
    }
    payload = json.dumps(result, sort_keys=True, separators=(",", ":")).encode()
    result["summary_sha256"] = sha256(payload).hexdigest()
    expected = (
        EXPECTED_ALL_SUMMARY_SHA256
        if all_schedules
        else EXPECTED_A0_SUMMARY_SHA256
    )
    if result["summary_sha256"] != expected:
        scope = "all-371" if all_schedules else "a=0"
        raise AssertionError(f"frozen {scope} hostile summary changed")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--all",
        action="store_true",
        help="replay all 371 schedules (the default)",
    )
    selection.add_argument(
        "--a0-only",
        action="store_true",
        help="replay only the frozen 48-orbit a=0 schedule",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = replay(not args.a0_only)
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
