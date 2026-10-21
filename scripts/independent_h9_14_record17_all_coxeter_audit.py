#!/usr/bin/env python3
"""Hostile no-import replay of the all-orbit record-17 Coxeter no-go.

This script deliberately does not import the production all-orbit gate, the
project's Coxeter classifier, or its Lanner tables.  It reconstructs the
order-four and order-five Lanner placement tables from the classical 9+5
representatives, implements a separate finite Coxeter classifier, rebuilds
the eight nerves from the frozen whole-shadow certificate, and repeats every
incremental actual-vertex restriction.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/enumerations/h9_14_record17_certificate.json"
ARCHIVE = ROOT / "data/enumerations/h9_14_record17_all_coxeter_labels.json"
INPUT_FILE_SHA256 = (
    "db613bd1b3f6606f7b140f4183f5fe2897db6ef8292d72670e38523f7c397d94"
)
INPUT_CERTIFICATE_SHA256 = (
    "b98e58b390cbf8d9b40e08da16da29f78a2f04805bb0b13d7f2efb623b9afb51"
)
FACET_COUNT = 14
DIMENSION = 9

SUPPORT_SEQUENCES = (
    (1, 2, 0),
    (2, 3, 4, 0),
    (2, 3, 0),
    (2, 3, 4),
    (0, 1, 2),
    (3, 4, 5, 0),
    (4, 5, 6, 0),
    (1, 2, 0),
)

ORBIT_SIGNATURES = (
    (1, "H???FB_", (0, 1, 2, 3, 4), 6, 108),
    (1, "H???FB_", (3, 4, 5, 7, 8), 4, 108),
    (2, "H???F?^", (3, 4, 5, 6, 7), 1, 106),
    (2, "H???Fbo", (4, 5, 6, 7, 8), 1, 106),
    (0, "H??EEBo", (4, 5, 6, 7, 8), 1, 110),
    (1, "H??EEBw", (0, 1, 2, 3, 4), 1, 108),
    (1, "H??FFB_", (0, 1, 2, 3, 4), 6, 108),
    (2, "H??FF?^", (3, 4, 5, 6, 7), 1, 106),
)

# Unlisted pairs have label two.  These are the classical nine order-four
# and five order-five Lanner isomorphism representatives.
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


def _pair(left: int, right: int) -> Pair:
    if left == right:
        raise ValueError("a pair needs distinct vertices")
    return (left, right) if left < right else (right, left)


def _canonical_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


def _summary_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("summary_sha256", None)
    return sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode("utf-8")
    ).hexdigest()


def _state_digest(states: Sequence[LabelState]) -> str:
    payload = [
        [[left, right, label] for (left, right), label in state]
        for state in sorted(states)
    ]
    return sha256(
        json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode(
            "utf-8"
        )
    ).hexdigest()


def _path_labels(
    component: frozenset[int],
    adjacency: Mapping[int, tuple[tuple[int, int], ...]],
) -> tuple[int, ...]:
    endpoint = min(
        vertex for vertex in component if len(adjacency[vertex]) == 1
    )
    labels: list[int] = []
    previous: int | None = None
    current = endpoint
    while True:
        following = [item for item in adjacency[current] if item[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("independent finite path traversal branched")
        next_vertex, label = following[0]
        labels.append(label)
        previous, current = current, next_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def _arm_lengths(
    component: frozenset[int],
    adjacency: Mapping[int, tuple[tuple[int, int], ...]],
) -> tuple[int, ...]:
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


def _finite_component(
    component: frozenset[int], labels: Mapping[Pair, int]
) -> bool:
    order = len(component)
    if order == 1:
        return True
    adjacency = {
        vertex: tuple(
            (other, labels[_pair(vertex, other)])
            for other in component
            if other != vertex and labels[_pair(vertex, other)] > 2
        )
        for vertex in component
    }
    if sum(map(len, adjacency.values())) // 2 != order - 1:
        return False
    degrees = Counter(len(adjacency[vertex]) for vertex in component)
    if max(degrees) <= 2:
        sequence = _path_labels(component, adjacency)
        if order == 2:
            return sequence[0] in (3, 4, 5)
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
        labels[_pair(left, right)] != 3
        for left, right in combinations(component, 2)
        if labels[_pair(left, right)] > 2
    ):
        return False
    arms = _arm_lengths(component, adjacency)
    return arms == (1, 1, order - 3) or arms in {
        (1, 2, 2),
        (1, 2, 3),
        (1, 2, 4),
    }


def _finite(facets: Iterable[int], labels: Mapping[Pair, int]) -> bool:
    vertices = tuple(sorted(set(facets)))
    if any(
        _pair(left, right) not in labels
        for left, right in combinations(vertices, 2)
    ):
        return False
    unseen = set(vertices)
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
                if labels[_pair(vertex, other)] > 2
            }
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        components.append(frozenset(component))
    return all(_finite_component(component, labels) for component in components)


def _representative_code(
    order: int, representative: Sequence[tuple[Pair, int]]
) -> tuple[int, ...]:
    edges = tuple(combinations(range(order), 2))
    labels = {edge: 2 for edge in edges}
    labels.update(representative)
    return tuple(labels[edge] for edge in edges)


def _permuted_code(
    order: int, code: Sequence[int], permutation: Sequence[int]
) -> tuple[int, ...]:
    edges = tuple(combinations(range(order), 2))
    source = dict(zip(edges, code, strict=True))
    moved = {
        _pair(permutation[left], permutation[right]): label
        for (left, right), label in source.items()
    }
    return tuple(moved[edge] for edge in edges)


def _lanner_codes(order: int) -> tuple[tuple[int, ...], ...]:
    codes = tuple(
        sorted(
            {
                _permuted_code(
                    order,
                    _representative_code(order, representative),
                    permutation,
                )
                for representative in LANNER_REPRESENTATIVES[order]
                for permutation in permutations(range(order))
            }
        )
    )
    expected = {4: (108, 9), 5: (420, 5)}[order]
    canonical = {
        min(
            _permuted_code(order, code, permutation)
            for permutation in permutations(range(order))
        )
        for code in codes
    }
    if (len(codes), len(canonical)) != expected:
        raise AssertionError("independent Lanner table census changed")
    local_edges = tuple(combinations(range(order), 2))
    for code in codes:
        labels = dict(zip(local_edges, code, strict=True))
        if _finite(range(order), labels):
            raise AssertionError("an independent Lanner row is finite")
        for deleted in range(order):
            if not _finite(
                (vertex for vertex in range(order) if vertex != deleted),
                labels,
            ):
                raise AssertionError("an independent Lanner row is not minimal")
    return codes


def _options(
    support: tuple[int, ...], tables: Mapping[int, Sequence[Sequence[int]]]
) -> tuple[LabelState, ...]:
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True)) for code in tables[len(support)]
    )


def _rank_nine_vertices(
    minimal_nonfaces: Sequence[tuple[int, ...]],
) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(frozenset(face) for face in minimal_nonfaces)
    return tuple(
        facets
        for facets in combinations(range(FACET_COUNT), DIMENSION)
        if not any(face <= set(facets) for face in forbidden)
    )


def _extract_orbits(document: Mapping[str, object]) -> tuple[dict[str, object], ...]:
    enumeration = document["whole_shadow_audit"]["enumeration"]
    if enumeration["residual_orbit_count"] != 8:
        raise AssertionError("independent residual orbit count changed")
    result = []
    for index, orbit in enumerate(enumeration["residual_orbits"]):
        signature3 = (
            int(orbit["missing_four_count"]),
            str(orbit["x_graph6"]),
            tuple(int(value) for value in orbit["p_local_vertices"]),
        )
        matches = [
            witness
            for witness in enumeration["witnesses"]
            if (
                int(witness["missing_four_count"]),
                str(witness["x_graph6"]),
                tuple(int(value) for value in witness["p_local_vertices"]),
            )
            == signature3
        ]
        if len(matches) != 1:
            raise AssertionError("independent residual representative is ambiguous")
        witness = matches[0]
        if witness["p_local_vertices"] != witness["q_local_vertices"]:
            raise AssertionError("independent residual cuts are unequal")
        by_size = {
            int(size): tuple(tuple(face) for face in faces)
            for size, faces in witness[
                "minimal_nonface_vertices_by_size"
            ].items()
        }
        all_nonfaces = tuple(face for faces in by_size.values() for face in faces)
        vertices = _rank_nine_vertices(all_nonfaces)
        signature = (
            *signature3,
            int(orbit["representation_count"]),
            int(orbit["facet_count"]),
        )
        if signature != ORBIT_SIGNATURES[index] or len(vertices) != signature[-1]:
            raise AssertionError("independent residual orbit signature changed")
        if {size: len(faces) for size, faces in by_size.items()} != {
            int(size): int(count)
            for size, count in orbit["minimal_nonface_histogram"].items()
        }:
            raise AssertionError("independent minimal-nonface histogram changed")
        result.append(
            {
                "signature": signature,
                "by_size": by_size,
                "higher": by_size.get(4, ()) + by_size.get(5, ()),
                "dotted": by_size[2],
                "vertices": vertices,
            }
        )
    return tuple(result)


def _extend(
    states: Sequence[LabelState],
    previous_union: Iterable[int],
    support_index: int,
    support: tuple[int, ...],
    vertices: Sequence[tuple[int, ...]],
    dotted_pairs: Sequence[tuple[int, ...]],
    tables: Mapping[int, Sequence[Sequence[int]]],
) -> tuple[tuple[LabelState, ...], dict[str, object]]:
    union = tuple(sorted(set(previous_union) | set(support)))
    union_set = set(union)
    restrictions = tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in union_set)
                for vertex in vertices
            },
            key=lambda value: (len(value), value),
        )
    )
    if any(
        set(dotted) <= set(restriction)
        for dotted in dotted_pairs
        for restriction in restrictions
    ):
        raise AssertionError("independent actual restriction contains a dotted pair")
    options = _options(support, tables)
    passing = set()
    compatible = 0
    for state in states:
        old = dict(state)
        for option in options:
            if not all(old.get(edge, label) == label for edge, label in option):
                continue
            compatible += 1
            merged = old | dict(option)
            labels = {
                _pair(left, right): 2 for left, right in combinations(union, 2)
            }
            labels.update(merged)
            if all(_finite(restriction, labels) for restriction in restrictions):
                passing.add(tuple(sorted(merged.items())))
    passing_states = tuple(sorted(passing))
    return passing_states, {
        "support_index": support_index,
        "support": list(support),
        "support_order": len(support),
        "complete_lanner_placement_count": len(options),
        "cumulative_cluster_union": list(union),
        "actual_vertex_restriction_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(map(len, restrictions)).items())
        },
        "raw_extension_count": len(states) * len(options),
        "compatible_lanner_extension_count": compatible,
        "elliptic_actual_vertex_survivor_count": len(passing_states),
        "passing_state_sha256": _state_digest(passing_states),
    }


def replay(input_path: Path = INPUT) -> dict[str, object]:
    raw = input_path.read_bytes()
    if sha256(raw).hexdigest() != INPUT_FILE_SHA256:
        raise ValueError("independent whole-shadow input file hash changed")
    document = json.loads(raw)
    if (
        document.get("certificate_sha256") != INPUT_CERTIFICATE_SHA256
        or _canonical_digest(document) != INPUT_CERTIFICATE_SHA256
    ):
        raise ValueError("independent whole-shadow certificate replay failed")

    tables = {4: _lanner_codes(4), 5: _lanner_codes(5)}
    orbit_records = []
    for index, orbit in enumerate(_extract_orbits(document)):
        states: tuple[LabelState, ...] = ((),)
        union: tuple[int, ...] = ()
        stages = []
        for support_index in SUPPORT_SEQUENCES[index]:
            support = orbit["higher"][support_index]
            states, stage = _extend(
                states,
                union,
                support_index,
                support,
                orbit["vertices"],
                orbit["dotted"],
                tables,
            )
            stages.append(stage)
            union = tuple(stage["cumulative_cluster_union"])
            if not states:
                break
        orbit_records.append(
            {
                "orbit_index": index,
                "orbit_signature": [
                    orbit["signature"][0],
                    orbit["signature"][1],
                    list(orbit["signature"][2]),
                    orbit["signature"][3],
                    orbit["signature"][4],
                ],
                "rank_nine_vertex_count": len(orbit["vertices"]),
                "decisive_support_indices": list(SUPPORT_SEQUENCES[index]),
                "stages": stages,
                "compact_coxeter_no_go": not states,
            }
        )
    result: dict[str, object] = {
        "method": (
            "no production imports; independent 9+5 Lanner representative "
            "closure, finite Coxeter classifier, nerve reconstruction, and "
            "incremental actual-vertex replay"
        ),
        "input_file_sha256": INPUT_FILE_SHA256,
        "input_certificate_sha256": INPUT_CERTIFICATE_SHA256,
        "lanner_tables": {
            "order_four_placements": len(tables[4]),
            "order_four_isomorphism_types": len(LANNER_REPRESENTATIVES[4]),
            "order_five_placements": len(tables[5]),
            "order_five_isomorphism_types": len(LANNER_REPRESENTATIVES[5]),
        },
        "orbit_gates": orbit_records,
        "all_eight_orbits_excluded": len(orbit_records) == 8
        and all(orbit["compact_coxeter_no_go"] for orbit in orbit_records),
    }
    result["summary_sha256"] = _summary_digest(result)
    return result


def verify_archive(archive_path: Path = ARCHIVE) -> dict[str, object]:
    replayed = replay()
    archived = json.loads(archive_path.read_text(encoding="utf-8"))
    if archived.get("kind") != "h9_14_record17_all_orbit_compact_coxeter_no_go":
        raise ValueError("unexpected all-orbit archive kind")
    hostile = archived["independent_hostile_replay"]
    if hostile != replayed:
        raise ValueError("archived hostile replay differs from no-import replay")
    theorem_orbits = archived["theorem"]["orbit_gates"]
    if len(theorem_orbits) != len(replayed["orbit_gates"]):
        raise ValueError("production/hostile orbit counts differ")
    for production, independent in zip(
        theorem_orbits, replayed["orbit_gates"], strict=True
    ):
        if (
            production["orbit_index"] != independent["orbit_index"]
            or production["rank_nine_vertex_count"]
            != independent["rank_nine_vertex_count"]
            or production["decisive_support_indices"]
            != independent["decisive_support_indices"]
            or production["incremental_lanner_and_actual_vertex_gate"]
            != independent["stages"]
            or not production["compact_coxeter_no_go"]
        ):
            raise ValueError("production/hostile orbit stage replay differs")
    return {
        "all_eight_orbits_excluded": replayed["all_eight_orbits_excluded"],
        "orbit_final_survivors": [
            orbit["stages"][-1]["elliptic_actual_vertex_survivor_count"]
            for orbit in replayed["orbit_gates"]
        ],
        "summary_sha256": replayed["summary_sha256"],
    }


def main() -> None:
    result = verify_archive() if ARCHIVE.exists() else replay()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
