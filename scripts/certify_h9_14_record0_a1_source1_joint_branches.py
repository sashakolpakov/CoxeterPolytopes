#!/usr/bin/env python3
"""Certify the 19-way labelled-A4 exhaustion for Record 0, a=1, source 1.

This is an offline aggregate verifier.  It never invokes CP-SAT or performs a
branch search.  Instead it independently reconstructs the labelled A4 orbit
partition, the exact stabilizer generators and model-shape counts, and the
abstract nonfinite-pattern/cut expansion used by every terminal checkpoint.
Only 19 audited terminal ``INFEASIBLE`` branches imply the fixed-source
compact no-go.  Finite-volume noncompact chambers are outside the scope.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import factorial, isfinite
from pathlib import Path
import re
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
ENUMERATIONS = ROOT / "data/enumerations"
SEED = ENUMERATIONS / "h9_14_record0_a1_source1_joint_label_cegar_checkpoint.json"
OUTPUT = (
    ENUMERATIONS
    / "h9_14_record0_a1_source1_joint_branch_certificate.json"
)
CHECKPOINT_GLOB = "h9_14_record0_a1_source1_branch*_joint_checkpoint.json"
CHECKPOINT_RE = re.compile(
    r"h9_14_record0_a1_source1_branch(?P<branch>[0-9]{2})_joint_checkpoint\.json"
)

N = 14
EXPECTED_BRANCH_INDICES = tuple(range(19))
EXPECTED_PARTITION_DIGEST = (
    "df9a46605bdf4899a5994ef4e8ce9a8002f9c1cc0d832341c547eefba433f938"
)
EXPECTED_PATTERN_DIGEST = (
    "f9767cc4cd2bd873de238b315c0bac6c5d5509bc73549967da55c7f101bb0971"
)
EXPECTED_TRANSCRIPT_DIGEST = (
    "f3758234e404cb557a89d9caebf8e66e3a2c2065cd4f831a2a6e656076d45010"
)
EXPECTED_EXPANDED_CUT_COUNT = 1_769_801
EXPECTED_TRANSCRIPT = [
    {
        "seed_round": 0,
        "new_abstract_pattern_count": 16,
        "expanded_global_face_cut_count": 1_707_843,
        "expanded_global_face_cut_digest": (
            "6ff340050b5485917e88e150e1f67432a42aa529616cfea85848943ede442408"
        ),
    },
    {
        "seed_round": 1,
        "new_abstract_pattern_count": 6,
        "expanded_global_face_cut_count": 51_410,
        "expanded_global_face_cut_digest": (
            "fe00fefeaa6b4b0bfd41a7aa25ae3b3b917b147c2e6abc3f10ddfe10de666169"
        ),
    },
    {
        "seed_round": 2,
        "new_abstract_pattern_count": 1,
        "expanded_global_face_cut_count": 10_548,
        "expanded_global_face_cut_digest": (
            "6b97c1c26c882d15fb236143061924145141f1ab4f05490875349fad6dee8d90"
        ),
    },
]
EXPECTED_REPRESENTATIVES = (
    (2, 2, 3, 2, 3, 5),
    (2, 2, 3, 3, 2, 5),
    (2, 2, 4, 5, 2, 3),
    (2, 2, 5, 2, 3, 3),
    (2, 2, 5, 4, 2, 3),
    (2, 2, 5, 5, 2, 3),
    (2, 3, 3, 3, 4, 2),
    (2, 3, 3, 3, 5, 2),
    (2, 3, 4, 3, 3, 2),
    (2, 3, 4, 4, 3, 2),
    (2, 3, 4, 5, 2, 2),
    (2, 3, 4, 5, 3, 2),
    (2, 3, 5, 2, 3, 2),
    (2, 3, 5, 3, 3, 2),
    (2, 3, 5, 4, 2, 2),
    (2, 3, 5, 4, 3, 2),
    (2, 3, 5, 5, 2, 2),
    (2, 3, 5, 5, 3, 2),
    (3, 3, 5, 2, 2, 2),
)
EXPECTED_ORBIT_SIZES = (6, 6, 6, 3, *([6] * 14), 3)
EXPECTED_MODEL_COUNTS = {
    "candidate_facet_complements": 754,
    "conditional_minimal_A5_compact_lanner_encodings": 1608,
    "face_witness_constraints": 1324,
    "global_meeting_label_variables": 89,
    "independent_five_sets": 1617,
    "lanner_order_constraints": 6796,
    "ridge_constraints": 1491,
}
EXPECTED_CERTIFICATE_SHA256 = (
    "da685ee5fcfbe0ee6ba1e8951b5c115dd73cb6ca6e3bccebf2f4487c24f05777"
)

Pair = tuple[int, int]
Code = tuple[int, ...]
Pattern = tuple[int, Code]
Action = tuple[int, ...]


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def _pair(left: int, right: int) -> Pair:
    if left == right:
        raise ValueError("a pair needs distinct vertices")
    return (left, right) if left < right else (right, left)


MISSING_EDGES = (_mask((0, 13)), _mask((1, 13)))
FIXED_A4 = _mask((0, 2, 3, 4))


def _permuted_code(order: int, code: Code, action: Sequence[int]) -> Code:
    pairs = tuple(combinations(range(order), 2))
    if len(code) != len(pairs) or sorted(action) != list(range(order)):
        raise ValueError("malformed labelled-diagram action")
    source = dict(zip(pairs, code, strict=True))
    target = {
        _pair(action[left], action[right]): label
        for (left, right), label in source.items()
    }
    return tuple(target[edge] for edge in pairs)


@lru_cache(maxsize=None)
def _code_orbit(order: int, code: Code) -> tuple[Code, ...]:
    return tuple(
        sorted(
            {
                _permuted_code(order, code, action)
                for action in permutations(range(order))
            }
        )
    )


# The classical nine order-four Lanner diagrams, stated independently of the
# production label module.  Omitted edges have label 2.
LANNER4_TYPES: tuple[tuple[tuple[Pair, int], ...], ...] = (
    (((0, 1), 3), ((0, 2), 3), ((0, 3), 5)),
    (((0, 1), 5), ((0, 2), 3), ((1, 3), 3)),
    (((0, 1), 3), ((0, 2), 4), ((1, 3), 5)),
    (((0, 1), 3), ((0, 2), 5), ((1, 3), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 3), 4)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 3), 5)),
    (((0, 1), 3), ((0, 2), 4), ((1, 3), 4), ((2, 3), 3)),
    (((0, 1), 3), ((0, 2), 4), ((1, 3), 5), ((2, 3), 3)),
    (((0, 1), 3), ((0, 2), 5), ((1, 3), 5), ((2, 3), 3)),
)


def _lanner4_placements() -> tuple[Code, ...]:
    pairs = tuple(combinations(range(4), 2))
    placements: set[Code] = set()
    for nonright in LANNER4_TYPES:
        labels = {edge: 2 for edge in pairs}
        labels.update(nonright)
        code = tuple(labels[edge] for edge in pairs)
        placements.update(_code_orbit(4, code))
    if len(placements) != 108:
        raise AssertionError("independent order-four Lanner census changed")
    return tuple(sorted(placements))


@lru_cache(maxsize=None)
def placement_orbits() -> tuple[tuple[Code, ...], ...]:
    # Local vertex 0 is the endpoint 0 of the missing-edge path and is fixed;
    # the source action on the other three A4 vertices is the full S3.
    actions = tuple((0, *image) for image in permutations((1, 2, 3)))
    remaining = set(_lanner4_placements())
    answer: list[tuple[Code, ...]] = []
    while remaining:
        representative = min(remaining)
        orbit = tuple(
            sorted(
                {
                    _permuted_code(4, representative, action)
                    for action in actions
                }
            )
        )
        if not set(orbit) <= remaining:
            raise AssertionError("placement orbit blocks overlap")
        remaining.difference_update(orbit)
        answer.append(orbit)
    result = tuple(answer)
    if tuple(orbit[0] for orbit in result) != EXPECTED_REPRESENTATIVES:
        raise AssertionError("source-1 placement representatives changed")
    if tuple(map(len, result)) != EXPECTED_ORBIT_SIZES:
        raise AssertionError("source-1 placement orbit sizes changed")
    if sum(map(len, result)) != 108 or _digest(result) != EXPECTED_PARTITION_DIGEST:
        raise AssertionError("source-1 placement partition changed")
    return result


def partition_record() -> dict[str, object]:
    orbits = placement_orbits()
    return {
        "literal_placement_count": sum(map(len, orbits)),
        "orbit_count": len(orbits),
        "orbit_sizes": list(map(len, orbits)),
        "orbit_representatives": [list(orbit[0]) for orbit in orbits],
        "partition_digest": _digest(orbits),
        "coverage_is_disjoint_and_complete": (
            len({code for orbit in orbits for code in orbit}) == 108
            and {code for orbit in orbits for code in orbit}
            == set(_lanner4_placements())
        ),
    }


def _contains_base(value: int) -> bool:
    return any(value & row == row for row in (*MISSING_EDGES, FIXED_A4))


@lru_cache(maxsize=None)
def _candidate_facets() -> tuple[int, ...]:
    independent_fives = tuple(
        _mask(chosen)
        for chosen in combinations(range(N), 5)
        if all(_mask(chosen) & edge != edge for edge in MISSING_EDGES)
    )
    shadow = tuple(value for value in independent_fives if value & FIXED_A4 == FIXED_A4)
    base = (*MISSING_EDGES, FIXED_A4, *shadow)
    return tuple(
        value
        for value in (_mask(chosen) for chosen in combinations(range(N), 5))
        if all(value & row for row in base)
    )


@lru_cache(maxsize=None)
def independent_model_counts() -> dict[str, int]:
    independent_fives = tuple(
        value
        for value in (_mask(chosen) for chosen in combinations(range(N), 5))
        if all(value & edge != edge for edge in MISSING_EDGES)
    )
    facets = _candidate_facets()
    counts = {
        "candidate_facet_complements": len(facets),
        "conditional_minimal_A5_compact_lanner_encodings": sum(
            not _contains_base(value) for value in independent_fives
        ),
        "face_witness_constraints": sum(
            not _contains_base(_mask(chosen))
            for order in range(5)
            for chosen in combinations(range(N), order)
        ),
        "global_meeting_label_variables": sum(
            _mask(edge) not in MISSING_EDGES
            for edge in combinations(range(N), 2)
        ),
        "independent_five_sets": len(independent_fives),
        "lanner_order_constraints": sum(
            not _contains_base(_mask(chosen))
            for order in range(6, N + 1)
            for chosen in combinations(range(N), order)
        ),
        "ridge_constraints": sum(
            any(facet & six == facet for facet in facets)
            for six in (_mask(chosen) for chosen in combinations(range(N), 6))
        ),
    }
    if counts != EXPECTED_MODEL_COUNTS:
        raise AssertionError("independent joint-model counts changed")
    return counts


def _swap(left: int, right: int) -> Action:
    action = list(range(N))
    action[left], action[right] = action[right], action[left]
    return tuple(action)


def _move_mask(value: int, action: Sequence[int]) -> int:
    return sum(1 << action[vertex] for vertex in range(N) if value & (1 << vertex))


def _move_code_on_A4(code: Code, action: Sequence[int]) -> Code:
    support = (0, 2, 3, 4)
    source = dict(zip(combinations(support, 2), code, strict=True))
    target = {
        _pair(action[left], action[right]): label
        for (left, right), label in source.items()
    }
    return tuple(target[edge] for edge in combinations(support, 2))


def stabilizer_record(branch: int) -> dict[str, object]:
    orbit = placement_orbits()[branch]
    code = orbit[0]
    outside_generators = tuple(_swap(left, left + 1) for left in range(5, 12))
    internal_stabilizer: list[Action] = []
    for image in permutations((2, 3, 4)):
        action = list(range(N))
        for source, target in zip((2, 3, 4), image, strict=True):
            action[source] = target
        action_tuple = tuple(action)
        if _move_code_on_A4(code, action_tuple) == code:
            internal_stabilizer.append(action_tuple)
    identity = tuple(range(N))
    internal_generators = tuple(
        action for action in internal_stabilizer if action != identity
    )
    # Here the exact internal stabilizer has order one or two.  Consequently
    # its sole nonidentity element (when present), together with adjacent
    # transpositions generating S8, generates the full labelled-A4 stabilizer.
    if len(internal_stabilizer) not in (1, 2):
        raise AssertionError("unexpected internal code stabilizer")
    actions = (*outside_generators, *internal_generators)
    facets = set(_candidate_facets())
    for action in actions:
        if _move_mask(FIXED_A4, action) != FIXED_A4:
            raise AssertionError("stabilizer generator moved the A4")
        if _move_code_on_A4(code, action) != code:
            raise AssertionError("stabilizer generator moved the A4 code")
        if {_move_mask(facet, action) for facet in facets} != facets:
            raise AssertionError("stabilizer generator moved the facet ground")
    steps = [
        sum(_move_mask(facet, action) != facet for facet in facets)
        for action in actions
    ]
    if len(orbit) * len(internal_stabilizer) != factorial(3):
        raise AssertionError("S3 orbit-stabilizer identity failed")
    return {
        "source_action_structure": "S3({2,3,4}) x S8({5,...,12})",
        "source_action_order": factorial(3) * factorial(8),
        "exact_internal_code_stabilizer_order": len(internal_stabilizer),
        "exact_labelled_A4_stabilizer_order": (
            len(internal_stabilizer) * factorial(8)
        ),
        "generator_actions": [list(action) for action in actions],
        "generator_actions_digest": _digest(actions),
        "lex_generator_step_counts": steps,
        "generators_generate_exact_stabilizer": True,
        "all_generators_preserve_fixed_source_A4_code_and_facet_ground": True,
        "orbit_stabilizer_identity": (
            f"{len(orbit)} * {len(internal_stabilizer)} = 6"
        ),
    }


def _path_labels(
    component: frozenset[int], adjacency: Mapping[int, tuple[tuple[int, int], ...]]
) -> tuple[int, ...]:
    endpoint = min(vertex for vertex in component if len(adjacency[vertex]) == 1)
    labels: list[int] = []
    previous: int | None = None
    current = endpoint
    while True:
        following = [row for row in adjacency[current] if row[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("path traversal branched")
        following_vertex, label = following[0]
        labels.append(label)
        previous, current = current, following_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def _finite_component(component: frozenset[int], labels: Mapping[Pair, int]) -> bool:
    if len(component) == 1:
        return True
    adjacency = {
        vertex: tuple(
            (other, labels[_pair(vertex, other)])
            for other in component
            if other != vertex and labels[_pair(vertex, other)] > 2
        )
        for vertex in component
    }
    if sum(map(len, adjacency.values())) // 2 != len(component) - 1:
        return False
    degree_histogram = Counter(map(len, adjacency.values()))
    if max(degree_histogram) <= 2:
        sequence = _path_labels(component, adjacency)
        order = len(component)
        return (
            (order == 2 and sequence[0] in (3, 4, 5))
            or all(label == 3 for label in sequence)
            or sequence == (3,) * (order - 2) + (4,)
            or (order == 3 and sequence == (3, 5))
            or (order == 4 and sequence in ((3, 3, 5), (3, 4, 3)))
        )
    if degree_histogram[3] != 1 or any(degree > 3 for degree in degree_histogram):
        return False
    if any(
        labels[_pair(left, right)] != 3
        for left, right in combinations(component, 2)
        if labels[_pair(left, right)] > 2
    ):
        return False
    branch = next(vertex for vertex in component if len(adjacency[vertex]) == 3)
    lengths: list[int] = []
    for neighbor, _label in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(v for v, _ in adjacency[current] if v != previous)
            previous, current = current, following
            length += 1
        lengths.append(length)
    arms = tuple(sorted(lengths))
    order = len(component)
    return arms == (1, 1, order - 3) or arms in ((1, 2, 2), (1, 2, 3), (1, 2, 4))


def _finite(vertices: Iterable[int], labels: Mapping[Pair, int]) -> bool:
    vertices = tuple(sorted(vertices))
    unseen = set(vertices)
    components: list[frozenset[int]] = []
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


def _pattern_from_cut(cut: Mapping[str, object]) -> Pattern:
    raw_restriction = cut.get("restriction")
    raw_positive = cut.get("positive_labels")
    if not isinstance(raw_restriction, list) or not isinstance(raw_positive, list):
        raise ValueError("malformed legacy finite cut")
    restriction = tuple(int(vertex) for vertex in raw_restriction)
    if restriction != tuple(sorted(set(restriction))) or not all(
        0 <= vertex < N for vertex in restriction
    ):
        raise ValueError("malformed cut restriction")
    actual_labels = {
        _pair(left, right): 2 for left, right in combinations(restriction, 2)
    }
    seen: set[Pair] = set()
    for raw_row in raw_positive:
        if not isinstance(raw_row, list) or len(raw_row) != 3:
            raise ValueError("malformed positive-label row")
        left, right, label = map(int, raw_row)
        edge = _pair(left, right)
        if edge not in actual_labels or edge in seen or label not in (3, 4, 5):
            raise ValueError("invalid positive label in legacy cut")
        seen.add(edge)
        actual_labels[edge] = label
    actual_pairs = tuple(combinations(restriction, 2))
    code = tuple(actual_labels[edge] for edge in actual_pairs)
    order = len(restriction)
    return order, min(_code_orbit(order, code))


def _pattern_document(pattern: Pattern) -> dict[str, object]:
    order, code = pattern
    return {
        "order": order,
        "canonical_code": list(code),
        "positive_labels": [
            [left, right, label]
            for (left, right), label in zip(
                combinations(range(order), 2), code, strict=True
            )
            if label > 2
        ],
        "labelled_embedding_orbit_size": len(_code_orbit(order, code)),
    }


def _check_minimal_nonfinite(pattern: Pattern) -> None:
    order, code = pattern
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    if _finite(range(order), labels):
        raise ValueError("stored obstruction pattern is finite")
    for deleted in range(order):
        if not _finite((vertex for vertex in range(order) if vertex != deleted), labels):
            raise ValueError("stored obstruction pattern is not deletion-minimal")


def _expand_pattern_batch(
    patterns: Sequence[Pattern], aggregate: object
) -> tuple[int, str]:
    digest = sha256()
    count = 0
    for order, code in sorted(patterns):
        for target in combinations(range(N), order):
            restriction = _mask(target)
            if _contains_base(restriction):
                continue
            target_pairs = tuple(combinations(target, 2))
            for placement in _code_orbit(order, code):
                positive = [
                    [left, right, label]
                    for (left, right), label in zip(
                        target_pairs, placement, strict=True
                    )
                    if label > 2
                ]
                encoded = (
                    json.dumps([restriction, positive], separators=(",", ":")).encode()
                    + b"\n"
                )
                digest.update(encoded)
                aggregate.update(encoded)
                count += 1
    return count, digest.hexdigest()


@lru_cache(maxsize=None)
def replay_seed_cuts() -> dict[str, object]:
    seed = json.loads(SEED.read_text(encoding="utf-8"))
    if (
        seed.get("record_index") != 0
        or seed.get("graph6") != "M?????????????o??"
        or seed.get("schedule_a") != 1
        or seed.get("source_orbit_index") != 1
        or seed.get("representative_A4s") != [[0, 2, 3, 4]]
    ):
        raise ValueError("seed archive belongs to another fixed source")
    rounds = seed.get("rounds")
    if not isinstance(rounds, list) or len(rounds) != 3:
        raise ValueError("seed archive does not contain the three cut rounds")
    installed: set[Pattern] = set()
    transcript: list[dict[str, object]] = []
    aggregate = sha256()
    legacy_cut_count = 0
    round_pattern_digests: list[str] = []
    for expected_round, row in enumerate(rounds):
        if not isinstance(row, Mapping) or row.get("round") != expected_round:
            raise ValueError("seed round sequence changed")
        cuts = row.get("new_global_monotone_face_cuts")
        if not isinstance(cuts, list) or row.get(
            "new_global_monotone_face_cut_count"
        ) != len(cuts):
            raise ValueError("legacy seed-cut count changed")
        legacy_cut_count += len(cuts)
        patterns = tuple(sorted({_pattern_from_cut(cut) for cut in cuts} - installed))
        for pattern in patterns:
            _check_minimal_nonfinite(pattern)
        installed.update(patterns)
        documents = [_pattern_document(pattern) for pattern in patterns]
        round_pattern_digests.append(_digest(documents))
        count, cut_digest = _expand_pattern_batch(patterns, aggregate)
        transcript.append(
            {
                "seed_round": expected_round,
                "new_abstract_pattern_count": len(patterns),
                "expanded_global_face_cut_count": count,
                "expanded_global_face_cut_digest": cut_digest,
            }
        )
    documents = [_pattern_document(pattern) for pattern in sorted(installed)]
    pattern_digest = _digest(documents)
    transcript_digest = _digest(transcript)
    if len(documents) != 23 or pattern_digest != EXPECTED_PATTERN_DIGEST:
        raise AssertionError("independent abstract-pattern replay changed")
    if transcript != EXPECTED_TRANSCRIPT or transcript_digest != EXPECTED_TRANSCRIPT_DIGEST:
        raise AssertionError("independent expanded-cut transcript changed")
    if sum(row["expanded_global_face_cut_count"] for row in transcript) != EXPECTED_EXPANDED_CUT_COUNT:
        raise AssertionError("expanded cut total changed")
    return {
        "seed_path": str(SEED.relative_to(ROOT)),
        "seed_file_sha256": _file_sha256(SEED),
        "legacy_concrete_seed_cut_count": legacy_cut_count,
        "abstract_nonfinite_pattern_count": len(documents),
        "abstract_nonfinite_pattern_order_histogram": {
            str(order): count
            for order, count in sorted(
                Counter(document["order"] for document in documents).items()
            )
        },
        "all_abstract_patterns_deletion_minimal_nonfinite": True,
        "abstract_nonfinite_patterns": documents,
        "abstract_nonfinite_pattern_digest": pattern_digest,
        "round_abstract_pattern_digests": round_pattern_digests,
        "expansion_transcript": transcript,
        "expansion_transcript_digest": transcript_digest,
        "expanded_global_monotone_face_cut_count": EXPECTED_EXPANDED_CUT_COUNT,
        "aggregate_expanded_cut_stream_digest": aggregate.hexdigest(),
        "independent_streaming_replay": True,
    }


REQUIRED_CHECKPOINT_KEYS = {
    "a4_branch",
    "a4_branch_compact_coxeter_no_go",
    "a4_placement_partition",
    "compact_label_survivor",
    "completed_cut_round_count",
    "dual_model_counts",
    "finite_volume_noncompact_checked",
    "fixed_source_compact_coxeter_no_go",
    "graph6",
    "installed_abstract_nonfinite_pattern_count",
    "installed_global_monotone_face_cut_count",
    "kind",
    "record_index",
    "representative_A4s",
    "rounds",
    "schedule_a",
    "scope",
    "seed_expansion_transcript",
    "seed_expansion_transcript_digest",
    "seeded_abstract_nonfinite_pattern_count",
    "seeded_abstract_nonfinite_pattern_digest",
    "seeded_abstract_nonfinite_patterns",
    "seeded_global_monotone_face_cut_count",
    "source_orbit_index",
    "stopped_reason",
    "terminal_solver_status",
    "theorem_certificate",
}
OPTIONAL_CHECKPOINT_KEYS = {
    "maximum_seed_pattern_order",
    "terminal_solver_branch_count",
    "terminal_solver_conflict_count",
    "terminal_solver_wall_time_seconds",
}
BRANCH_SCOPE = (
    "Complete joint compact incidence and finite-label formulation for one "
    "fixed a=1 source/A4-label branch. Only terminal INFEASIBLE closes the "
    "branch; all placement branches must close to infer a fixed-source theorem. "
    "An all-elliptic label model still requires exact Gram."
)
LOSSLESS_REASON = (
    "Move the A4 placement to its unique placement-orbit representative, then "
    "choose a lexicographically least incidence vector under the exact-code "
    "stabilizer."
)


def audit_checkpoint(path: Path, cut_audit: Mapping[str, object]) -> dict[str, object]:
    match = CHECKPOINT_RE.fullmatch(path.name)
    if match is None:
        raise ValueError(f"malformed branch checkpoint name: {path.name}")
    branch = int(match.group("branch"))
    if branch not in EXPECTED_BRANCH_INDICES:
        raise ValueError(f"out-of-range branch checkpoint: {path.name}")
    document = json.loads(path.read_text(encoding="utf-8"))
    keys = set(document)
    if not REQUIRED_CHECKPOINT_KEYS <= keys or not keys <= (
        REQUIRED_CHECKPOINT_KEYS | OPTIONAL_CHECKPOINT_KEYS
    ):
        raise ValueError(f"branch {branch}: checkpoint schema changed")
    common = {
        "kind": "h9_14_record0_a1_joint_incidence_label_cegar_exploratory",
        "theorem_certificate": False,
        "record_index": 0,
        "graph6": "M?????????????o??",
        "schedule_a": 1,
        "source_orbit_index": 1,
        "representative_A4s": [[0, 2, 3, 4]],
        "dual_model_counts": independent_model_counts(),
        "completed_cut_round_count": 0,
        "seeded_abstract_nonfinite_pattern_count": 23,
        "seeded_abstract_nonfinite_patterns": cut_audit[
            "abstract_nonfinite_patterns"
        ],
        "seeded_abstract_nonfinite_pattern_digest": EXPECTED_PATTERN_DIGEST,
        "seeded_global_monotone_face_cut_count": EXPECTED_EXPANDED_CUT_COUNT,
        "seed_expansion_transcript": EXPECTED_TRANSCRIPT,
        "seed_expansion_transcript_digest": EXPECTED_TRANSCRIPT_DIGEST,
        "installed_abstract_nonfinite_pattern_count": 23,
        "installed_global_monotone_face_cut_count": EXPECTED_EXPANDED_CUT_COUNT,
        "rounds": [],
        "terminal_solver_status": "INFEASIBLE",
        "stopped_reason": "terminal_joint_incidence_label_infeasible",
        "compact_label_survivor": None,
        "a4_branch_compact_coxeter_no_go": True,
        "fixed_source_compact_coxeter_no_go": False,
        "finite_volume_noncompact_checked": False,
        "scope": BRANCH_SCOPE,
    }
    for key, expected in common.items():
        if document.get(key) != expected:
            raise ValueError(f"branch {branch}: {key} changed")
    if "maximum_seed_pattern_order" in document and document[
        "maximum_seed_pattern_order"
    ] is not None:
        raise ValueError(f"branch {branch}: seed-pattern order filter changed")
    partition = partition_record()
    if document.get("a4_placement_partition") != partition:
        raise ValueError(f"branch {branch}: placement partition changed")
    stabilizer = stabilizer_record(branch)
    expected_branch = {
        "orbit_index": branch,
        "representative_code": list(EXPECTED_REPRESENTATIVES[branch]),
        "placement_orbit_size": EXPECTED_ORBIT_SIZES[branch],
        "incidence_symmetry_slice": "exact labelled-A4 stabilizer",
        "lex_generator_step_counts": stabilizer["lex_generator_step_counts"],
        "lossless_reason": LOSSLESS_REASON,
    }
    if document.get("a4_branch") != expected_branch:
        raise ValueError(f"branch {branch}: branch/stabilizer metadata changed")
    statistics_keys = (
        "terminal_solver_wall_time_seconds",
        "terminal_solver_branch_count",
        "terminal_solver_conflict_count",
    )
    present = [key in document for key in statistics_keys]
    if any(present) and not all(present):
        raise ValueError(f"branch {branch}: incomplete terminal statistics")
    if not all(present) and branch != 0:
        raise ValueError(f"branch {branch}: terminal statistics are missing")
    statistics: dict[str, object] | None = None
    if all(present):
        wall = document[statistics_keys[0]]
        branches = document[statistics_keys[1]]
        conflicts = document[statistics_keys[2]]
        if (
            isinstance(wall, bool)
            or not isinstance(wall, (int, float))
            or not isfinite(wall)
            or wall < 0
            or isinstance(branches, bool)
            or not isinstance(branches, int)
            or branches < 0
            or isinstance(conflicts, bool)
            or not isinstance(conflicts, int)
            or conflicts < 0
        ):
            raise ValueError(f"branch {branch}: invalid terminal statistics")
        statistics = {
            "wall_time_seconds": wall,
            "solver_branch_count": branches,
            "solver_conflict_count": conflicts,
        }
    try:
        checkpoint_path = str(path.relative_to(ROOT))
    except ValueError:  # Enables hostile/legacy fixture checks outside the tree.
        checkpoint_path = str(path)
    return {
        "branch_index": branch,
        "checkpoint_path": checkpoint_path,
        "checkpoint_file_sha256": _file_sha256(path),
        "representative_code": list(EXPECTED_REPRESENTATIVES[branch]),
        "placement_orbit_size": EXPECTED_ORBIT_SIZES[branch],
        "terminal_solver_status": "INFEASIBLE",
        "compact_label_survivor": None,
        "branch_compact_coxeter_no_go": True,
        "checkpoint_fixed_source_compact_coxeter_no_go": False,
        "terminal_statistics_available": statistics is not None,
        "terminal_statistics": statistics,
        "exact_stabilizer": stabilizer,
    }


def aggregate_coverage(branches: Sequence[Mapping[str, object]]) -> dict[str, object]:
    indices = [int(row["branch_index"]) for row in branches]
    if len(indices) != len(set(indices)):
        raise ValueError("duplicate branch index in aggregate")
    present = sorted(indices)
    missing = sorted(set(EXPECTED_BRANCH_INDICES) - set(present))
    unexpected = sorted(set(present) - set(EXPECTED_BRANCH_INDICES))
    closed = sorted(
        int(row["branch_index"])
        for row in branches
        if row.get("branch_compact_coxeter_no_go") is True
        and row.get("terminal_solver_status") == "INFEASIBLE"
        and row.get("compact_label_survivor") is None
    )
    orbit_size_sum = sum(int(row["placement_orbit_size"]) for row in branches)
    complete = (
        not missing
        and not unexpected
        and present == list(EXPECTED_BRANCH_INDICES)
        and closed == list(EXPECTED_BRANCH_INDICES)
        and orbit_size_sum == 108
    )
    return {
        "expected_branch_indices": list(EXPECTED_BRANCH_INDICES),
        "present_branch_indices": present,
        "missing_branch_indices": missing,
        "unexpected_branch_indices": unexpected,
        "terminally_closed_branch_indices": closed,
        "individual_branch_closure_count": len(closed),
        "present_placement_orbit_size_sum": orbit_size_sum,
        "required_placement_orbit_size_sum": 108,
        "all_19_terminal_branches_present": not missing and not unexpected,
        "all_19_terminal_branches_closed": complete,
        "fixed_source_compact_coxeter_no_go": complete,
        "aggregation_rule": (
            "A checkpoint closes only its labelled-A4 placement orbit. The "
            "fixed source closes iff the exact 19-orbit partition has size "
            "sum 108 and every branch is terminal INFEASIBLE with no survivor."
        ),
    }


def _checkpoint_paths() -> tuple[Path, ...]:
    return tuple(sorted(ENUMERATIONS.glob(CHECKPOINT_GLOB)))


def fresh_record(
    *,
    require_complete: bool = True,
    checkpoint_paths: Sequence[Path] | None = None,
) -> dict[str, object]:
    cut_audit = replay_seed_cuts()
    paths = tuple(checkpoint_paths) if checkpoint_paths is not None else _checkpoint_paths()
    branches = sorted(
        (audit_checkpoint(path, cut_audit) for path in paths),
        key=lambda row: int(row["branch_index"]),
    )
    coverage = aggregate_coverage(branches)
    if require_complete and coverage["fixed_source_compact_coxeter_no_go"] is not True:
        raise ValueError(
            "fixed-source closure requires all 19 terminal branches; missing "
            f"{coverage['missing_branch_indices']}"
        )
    fixed_closed = bool(coverage["fixed_source_compact_coxeter_no_go"])
    document: dict[str, object] = {
        "schema_version": "1.0.0",
        "kind": "h9_14_record0_a1_source1_joint_branch_aggregate_certificate",
        "claim": (
            "For Record 0, schedule a=1, fixed source orbit 1, no compact "
            "joint incidence/Coxeter labelling survives any of the 19 exact "
            "labelled-A4 placement branches."
            if fixed_closed
            else "Incomplete branch ledger; no fixed-source conclusion."
        ),
        "verifier": str(Path(__file__).resolve().relative_to(ROOT)),
        "proof_scope": (
            "One fixed compact Record-0 a=1 source orbit, jointly choosing "
            "incidence and finite Coxeter labels. Terminal solver statuses are "
            "audited artifacts; the expensive branch searches are not rerun."
        ),
        "fixed_source": {
            "record_index": 0,
            "graph6": "M?????????????o??",
            "schedule_a": 1,
            "source_orbit_index": 1,
            "representative_A4": [0, 2, 3, 4],
            "source_action_structure": "S3({2,3,4}) x S8({5,...,12})",
            "source_action_order": factorial(3) * factorial(8),
        },
        "placement_partition": partition_record(),
        "independent_model_counts": independent_model_counts(),
        "cut_library_audit": cut_audit,
        "branch_evidence": branches,
        "coverage": coverage,
        "theorem": {
            "closed_labelled_A4_branch_indices": coverage[
                "terminally_closed_branch_indices"
            ],
            "fixed_source_compact_coxeter_no_go": fixed_closed,
            "whole_record0_a1_schedule_excluded": False,
            "whole_record0_shadow_excluded": False,
            "finite_volume_noncompact_checked": False,
            "finite_volume_noncompact_chamber_excluded": False,
        },
        "nonclaims": [
            "A terminal branch checkpoint is not by itself a fixed-source theorem.",
            "This certificate does not close the other Record-0 a=1 source orbits.",
            "This certificate does not close the whole Record-0 graph shadow.",
            "No finite-volume noncompact chamber, ideal vertex, affine cusp, or maximal parabolic endpoint is checked.",
        ],
    }
    document["certificate_sha256"] = _digest(document)
    return document


def verify_record(
    document: Mapping[str, object], *, require_complete: bool = True
) -> dict[str, object]:
    if document.get("kind") != (
        "h9_14_record0_a1_source1_joint_branch_aggregate_certificate"
    ):
        raise ValueError("wrong source-1 joint-branch certificate kind")
    observed = _digest({key: value for key, value in document.items() if key != "certificate_sha256"})
    if document.get("certificate_sha256") != observed:
        raise ValueError("source-1 aggregate internal digest changed")
    if EXPECTED_CERTIFICATE_SHA256 != "TO_BE_FROZEN" and observed != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("source-1 aggregate certificate is not frozen")
    current = fresh_record(require_complete=require_complete)
    if dict(document) != current:
        raise ValueError("source-1 aggregate evidence or replay changed")
    theorem = document.get("theorem")
    if not isinstance(theorem, Mapping):
        raise ValueError("source-1 aggregate theorem record is malformed")
    if require_complete and theorem.get("fixed_source_compact_coxeter_no_go") is not True:
        raise ValueError("source-1 aggregate does not close the fixed source")
    if (
        theorem.get("whole_record0_a1_schedule_excluded") is not False
        or theorem.get("whole_record0_shadow_excluded") is not False
        or theorem.get("finite_volume_noncompact_checked") is not False
        or theorem.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise ValueError("source-1 aggregate scope changed")
    return {
        "record_index": 0,
        "schedule_a": 1,
        "source_orbit_index": 1,
        "closed_branch_count": len(document["branch_evidence"]),
        "placement_orbit_size_sum": document["coverage"][
            "present_placement_orbit_size_sum"
        ],
        "terminal_solver_status": "INFEASIBLE (all 19 branches)",
        "compact_label_survivor_count": 0,
        "fixed_source_compact_coxeter_no_go": theorem[
            "fixed_source_compact_coxeter_no_go"
        ],
        "whole_record0_a1_schedule_excluded": False,
        "finite_volume_noncompact_checked": False,
        "certificate_sha256": observed,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--allow-incomplete", action="store_true")
    args = parser.parse_args(argv)
    require_complete = not args.allow_incomplete
    if args.build or args.write:
        document = fresh_record(require_complete=require_complete)
    else:
        document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    summary = verify_record(document, require_complete=require_complete)
    if args.write:
        OUTPUT.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps({"output": str(OUTPUT), **summary}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
