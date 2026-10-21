#!/usr/bin/env python3
"""Hostile no-import replay of the all-371 record-12 Coxeter label gate.

No production module, project Coxeter classifier, Lanner table, or residual
reconstructor is imported.  The script independently rebuilds the 371 nerves
from the frozen matrices, regenerates the classical labelled Lanner tables,
classifies finite diagrams, and replays every deterministic support stage.
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
INPUT = ROOT / "data/enumerations/h9_14_record12_residual.json"
ARCHIVE = ROOT / "data/enumerations/h9_14_record12_all_coxeter_labels.json"
INPUT_FILE_SHA256 = (
    "faca367738bba58afc0ec78a43e9d4aa4b02721802c5d9e77b7d155c8e7faa05"
)
INPUT_CERTIFICATE_SHA256 = (
    "50ca012334e7a4973a4c29bde4226d4460bc8ade2d644a3b9dbf56b1bf92e8f4"
)

FACET_COUNT = 14
DIMENSION = 9
FULL_MASK = (1 << FACET_COUNT) - 1
ACTIVE_A = (0, 1)
ACTIVE_B = (12, 13)
LOCAL_ISOLATES = tuple(range(10))
BALANCED_CUT_MASK = sum(1 << vertex for vertex in range(5, 10))
PAIR_INDEX = {
    edge: index for index, edge in enumerate(combinations(LOCAL_ISOLATES, 2))
}
EXPECTED_SCHEDULE = {0: 48, 1: 95, 2: 144, 3: 63, 4: 21}
EXPECTED_DOTTED_MASKS = {
    _mask_pair
    for _mask_pair in (
        (1 << 0) | (1 << 12),
        (1 << 1) | (1 << 12),
        (1 << 0) | (1 << 13),
        (1 << 1) | (1 << 13),
    )
}


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


def _mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def _vertices(value: int, order: int = FACET_COUNT) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(order) if value & (1 << vertex))


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


def _matrix_to_pair_mask(matrix: int) -> int:
    return sum(
        1 << PAIR_INDEX[first, second + 5]
        for first in range(5)
        for second in range(5)
        if matrix & (1 << (5 * first + second))
    )


def _pair_edge(value: int, first: int, second: int) -> bool:
    return bool(value & (1 << PAIR_INDEX[_pair(first, second)]))


def _cut_edge(cut_mask: int, first: int, second: int) -> bool:
    return bool(cut_mask & (1 << first)) != bool(cut_mask & (1 << second))


def _branch_triples(x_mask: int, cut_mask: int) -> tuple[int, ...]:
    result: set[int] = set()
    for first, second in combinations(LOCAL_ISOLATES, 2):
        x = _pair_edge(x_mask, first, second)
        delta = _cut_edge(cut_mask, first, second)
        if x:
            result.add(_mask((11, first, second)))
        if x != delta:
            result.add(_mask((10, first, second)))
    for triple in combinations(LOCAL_ISOLATES, 3):
        if sum(
            _pair_edge(x_mask, first, second)
            for first, second in combinations(triple, 2)
        ) % 2:
            result.add(_mask(triple))
    return tuple(sorted(result))


def _selected_complements(x_matrix: int, z_matrix: int) -> tuple[int, ...]:
    result: set[int] = set()
    for triple in _branch_triples(
        _matrix_to_pair_mask(x_matrix), BALANCED_CUT_MASK
    ):
        global_triple = 0
        for local in LOCAL_ISOLATES:
            if triple & (1 << local):
                global_triple |= 1 << (local + 2)
        if triple & (1 << 10):
            global_triple |= 1 << 12
        if triple & (1 << 11):
            global_triple |= 1 << 13
        result.add(_mask(ACTIVE_A) | global_triple)
    for triple in _branch_triples(
        _matrix_to_pair_mask(z_matrix), BALANCED_CUT_MASK
    ):
        global_triple = 0
        for local in LOCAL_ISOLATES:
            if triple & (1 << local):
                global_triple |= 1 << (local + 2)
        if triple & (1 << 10):
            global_triple |= 1 << 0
        if triple & (1 << 11):
            global_triple |= 1 << 1
        result.add(_mask(ACTIVE_B) | global_triple)
    return tuple(sorted(result))


def _face_table(complements: Sequence[int]) -> tuple[bool, ...]:
    table = [False] * (1 << FACET_COUNT)
    for complement in complements:
        facet = FULL_MASK ^ complement
        subset = facet
        while True:
            table[subset] = True
            if subset == 0:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def _minimal_nonfaces(table: Sequence[bool]) -> dict[int, tuple[int, ...]]:
    result: dict[int, list[int]] = {}
    for value in range(1, 1 << FACET_COUNT):
        if table[value]:
            continue
        bits = value
        while bits:
            bit = bits & -bits
            if not table[value ^ bit]:
                break
            bits ^= bit
        else:
            result.setdefault(value.bit_count(), []).append(value)
    return {size: tuple(rows) for size, rows in sorted(result.items())}


def _path_labels(
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
    values = tuple(sorted(set(facets)))
    if any(
        _pair(left, right) not in labels
        for left, right in combinations(values, 2)
    ):
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
                (vertex for vertex in range(order) if vertex != deleted), labels
            ):
                raise AssertionError("an independent Lanner row is not minimal")
    return codes


def _options(
    support: tuple[int, ...], tables: Mapping[int, Sequence[Sequence[int]]]
) -> tuple[LabelState, ...]:
    edges = tuple(combinations(tuple(sorted(support)), 2))
    return tuple(
        tuple(zip(edges, code, strict=True)) for code in tables[len(support)]
    )


@lru_cache(maxsize=262_144)
def _finite_code(order: int, code: tuple[int, ...]) -> bool:
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return _finite(range(order), labels)


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
    rows = document.get("residual_orbits")
    if not isinstance(rows, list) or len(rows) != 371:
        raise AssertionError("independent record-12 orbit count changed")
    result = []
    schedule: Counter[int] = Counter()
    for index, row in enumerate(rows):
        x_matrix = int(row["x_matrix"])
        z_matrix = int(row["z_matrix"])
        complements = _selected_complements(x_matrix, z_matrix)
        nonface_masks = _minimal_nonfaces(_face_table(complements))
        by_size = {
            size: tuple(_vertices(value) for value in masks)
            for size, masks in nonface_masks.items()
        }
        all_nonfaces = tuple(face for faces in by_size.values() for face in faces)
        reconstructed = _rank_nine_vertices(all_nonfaces)
        direct = tuple(
            sorted(_vertices(FULL_MASK ^ complement) for complement in complements)
        )
        histogram = {size: len(faces) for size, faces in by_size.items()}
        archived_histogram = {
            int(size): int(count)
            for size, count in row["minimal_nonface_histogram"].items()
        }
        complement_digest = sha256(
            b"".join(value.to_bytes(2, "little") for value in complements)
        ).hexdigest()
        if (
            reconstructed != direct
            or histogram != archived_histogram
            or set(histogram) - {2, 4, 5}
            or set(nonface_masks.get(2, ())) != EXPECTED_DOTTED_MASKS
            or len(direct) != int(row["facet_count"])
            or complement_digest != row["selected_complement_sha256"]
        ):
            raise AssertionError(f"independent record-12 orbit {index} changed")
        missing_four = int(row["missing_four_count"])
        if histogram.get(4, 0) != missing_four:
            raise AssertionError("independent record-12 schedule changed")
        schedule[missing_four] += 1
        result.append(
            {
                "index": index,
                "signature": {
                    "x_matrix": x_matrix,
                    "z_matrix": z_matrix,
                    "missing_four_count": missing_four,
                    "relative_embedding_multiplicity": int(
                        row["relative_embedding_multiplicity"]
                    ),
                    "facet_count": int(row["facet_count"]),
                    "selected_complement_sha256": complement_digest,
                },
                "by_size": by_size,
                "higher": by_size.get(4, ()) + by_size.get(5, ()),
                "dotted": by_size[2],
                "vertices": direct,
            }
        )
    if dict(sorted(schedule.items())) != EXPECTED_SCHEDULE:
        raise AssertionError("independent record-12 schedule histogram changed")
    return tuple(result)


def _next_support_index(
    supports: Sequence[tuple[int, ...]],
    remaining: Iterable[int],
    cluster_union: Iterable[int],
) -> int:
    choices = tuple(sorted(remaining))
    union = set(cluster_union)
    if not union:
        return choices[0]
    return max(
        choices,
        key=lambda index: (
            len(set(supports[index]) & union),
            -len(set(supports[index]) - union),
            -index,
        ),
    )


def _extend(
    states: Sequence[LabelState],
    previous_union: Iterable[int],
    support_index: int,
    support: tuple[int, ...],
    actual_vertices: Sequence[tuple[int, ...]],
    dotted_pairs: Sequence[tuple[int, ...]],
    tables: Mapping[int, Sequence[Sequence[int]]],
) -> tuple[tuple[LabelState, ...], dict[str, object]]:
    union = tuple(sorted(set(previous_union) | set(support)))
    union_set = set(union)
    restrictions = tuple(
        sorted(
            {
                tuple(facet for facet in vertex if facet in union_set)
                for vertex in actual_vertices
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
    passing: set[LabelState] = set()
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
            finite = True
            for restriction in restrictions:
                code = tuple(
                    labels[_pair(left, right)]
                    for left, right in combinations(restriction, 2)
                )
                if not _finite_code(len(restriction), code):
                    finite = False
                    break
            if finite:
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


@lru_cache(maxsize=1)
def replay(input_path: Path = INPUT) -> dict[str, object]:
    raw = input_path.read_bytes()
    if sha256(raw).hexdigest() != INPUT_FILE_SHA256:
        raise ValueError("independent record-12 input file hash changed")
    document = json.loads(raw)
    if (
        document.get("certificate_sha256") != INPUT_CERTIFICATE_SHA256
        or _canonical_digest(document) != INPUT_CERTIFICATE_SHA256
    ):
        raise ValueError("independent record-12 certificate replay failed")

    tables = {4: _lanner_codes(4), 5: _lanner_codes(5)}
    orbit_records = []
    support_histogram: Counter[int] = Counter()
    schedule_histogram: Counter[int] = Counter()
    for orbit in _extract_orbits(document):
        supports = orbit["higher"]
        remaining = set(range(len(supports)))
        states: tuple[LabelState, ...] = ((),)
        union: tuple[int, ...] = ()
        stages = []
        while states and remaining:
            support_index = _next_support_index(supports, remaining, union)
            states, stage = _extend(
                states,
                union,
                support_index,
                supports[support_index],
                orbit["vertices"],
                orbit["dotted"],
                tables,
            )
            stages.append(stage)
            remaining.remove(support_index)
            union = tuple(stage["cumulative_cluster_union"])
        support_histogram[len(stages)] += 1
        schedule_histogram[orbit["signature"]["missing_four_count"]] += 1
        orbit_records.append(
            {
                "orbit_index": orbit["index"],
                "orbit_signature": orbit["signature"],
                "rank_nine_vertex_count": len(orbit["vertices"]),
                "minimal_nonface_histogram": {
                    str(size): len(faces)
                    for size, faces in orbit["by_size"].items()
                },
                "decisive_support_indices": [
                    stage["support_index"] for stage in stages
                ],
                "stages": stages,
                "compact_coxeter_no_go": not states,
            }
        )
    result: dict[str, object] = {
        "method": (
            "no project imports; independent matrix-to-nerve reconstruction, "
            "9+5 classical Lanner closure, finite Coxeter classifier, "
            "deterministic overlap-first support choice, and every actual "
            "rank-nine vertex restriction"
        ),
        "input_file_sha256": INPUT_FILE_SHA256,
        "input_certificate_sha256": INPUT_CERTIFICATE_SHA256,
        "lanner_tables": {
            "order_four_placements": len(tables[4]),
            "order_four_isomorphism_types": len(LANNER_REPRESENTATIVES[4]),
            "order_five_placements": len(tables[5]),
            "order_five_isomorphism_types": len(LANNER_REPRESENTATIVES[5]),
        },
        "schedule_exclusion_histogram": {
            str(key): value for key, value in sorted(schedule_histogram.items())
        },
        "support_count_histogram": {
            str(key): value for key, value in sorted(support_histogram.items())
        },
        "orbit_gates": orbit_records,
        "all_371_orbits_excluded": len(orbit_records) == 371
        and all(orbit["compact_coxeter_no_go"] for orbit in orbit_records),
    }
    result["summary_sha256"] = _summary_digest(result)
    return result


def compact_summary(replayed: Mapping[str, object]) -> dict[str, object]:
    return {
        "method": replayed["method"],
        "input_file_sha256": replayed["input_file_sha256"],
        "input_certificate_sha256": replayed["input_certificate_sha256"],
        "lanner_tables": replayed["lanner_tables"],
        "schedule_exclusion_histogram": replayed[
            "schedule_exclusion_histogram"
        ],
        "support_count_histogram": replayed["support_count_histogram"],
        "residual_orbit_count": len(replayed["orbit_gates"]),
        "all_371_orbits_excluded": replayed["all_371_orbits_excluded"],
        "full_replay_summary_sha256": replayed["summary_sha256"],
    }


def verify_archive(archive_path: Path = ARCHIVE) -> dict[str, object]:
    replayed = replay()
    archived = json.loads(archive_path.read_text(encoding="utf-8"))
    if archived.get("kind") != "h9_14_record12_all_orbit_compact_coxeter_no_go":
        raise ValueError("unexpected record-12 all-orbit archive kind")
    if archived["independent_hostile_replay"] != compact_summary(replayed):
        raise ValueError("archived hostile summary differs from no-import replay")
    theorem_orbits = archived["theorem"]["orbit_gates"]
    independent_orbits = replayed["orbit_gates"]
    if len(theorem_orbits) != len(independent_orbits):
        raise ValueError("production/hostile record-12 orbit counts differ")
    for production, independent in zip(
        theorem_orbits, independent_orbits, strict=True
    ):
        if (
            production["orbit_index"] != independent["orbit_index"]
            or production["orbit_signature"] != independent["orbit_signature"]
            or production["rank_nine_vertex_count"]
            != independent["rank_nine_vertex_count"]
            or production["minimal_nonface_histogram"]
            != independent["minimal_nonface_histogram"]
            or production["decisive_support_indices"]
            != independent["decisive_support_indices"]
            or production["incremental_lanner_and_actual_vertex_gate"]
            != independent["stages"]
            or not production["compact_coxeter_no_go"]
        ):
            raise ValueError("production/hostile record-12 stage replay differs")
    return {
        "residual_orbits": len(independent_orbits),
        "all_371_orbits_excluded": replayed["all_371_orbits_excluded"],
        "support_count_histogram": replayed["support_count_histogram"],
        "summary_sha256": replayed["summary_sha256"],
    }


def main() -> None:
    result = verify_archive() if ARCHIVE.exists() else compact_summary(replay())
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
