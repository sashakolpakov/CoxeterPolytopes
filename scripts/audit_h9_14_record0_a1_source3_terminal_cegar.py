#!/usr/bin/env python3
"""No-trust replay of the terminal Record-0/a=1/source-3 compact CEGAR no-go.

The input checkpoint contains 41 local compact-label contradictions.  This
audit does not import the exploratory CEGAR driver, its dual-model builder,
or its local-core extractor.  It independently:

* rebuilds the fixed-source dual CP-SAT model from graph6 and the pinned
  Dehn--Sommerville schedule;
* enumerates the complete order-four/order-five Lanner placement tables and
  replays every local elliptic/Lanner contradiction;
* generates every image of every core under the full fixed-source
  ``S3 x S8 x C2`` action by generator-orbit BFS;
* installs the exact lossless clauses and freezes the terminal model proto;
* solves that terminal integer model with one worker and no time limit.

The quantifier is intentionally narrow: this proves a compact Coxeter no-go
for one fixed source orbit only.  It says nothing about the other Record-0
sources, the rest of H9/14, or finite-volume noncompact polytopes.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from copy import deepcopy
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
import math
from pathlib import Path
import time
import tempfile
from typing import Any, Iterable, Mapping, Sequence

import ortools
from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source3_label_core_cegar_checkpoint_44_more.json"
)
SUPPORT = ROOT / "data/enumerations/h9_14_nerve_support.json"
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "h9_14_record0_a1_source3_terminal_cegar_audit.json"
)

SCHEMA = "coxeter.h9_14.record0_a1_source3.terminal_cegar_audit.v1"
STATUS_DATE = "2026-08-08"
CHECKPOINT_SHA256 = (
    "f0b38d563b0b532f39af73203e40050cfe65314f3db50894cac84e4a7f27ef7a"
)
SUPPORT_SHA256 = (
    "8121f869c182b1780f3d812d34823f1579a5c37b085a4a8ecc313ea94937c47d"
)
PINNED_ORTOOLS_VERSION = "9.15.6755"

N = 14
FULL = (1 << N) - 1
GRAPH6 = "M?????????????o??"
FIXED_FOUR_VERTICES = (2, 3, 4, 13)
FIXED_FOUR = sum(1 << vertex for vertex in FIXED_FOUR_VERTICES)
EXPECTED_EDGE_ROWS = ((1 << 0) | (1 << 13), (1 << 1) | (1 << 13))
EXPECTED_SOURCE_GROUP_ORDER = math.factorial(3) * math.factorial(8) * 2
EXPECTED_BASE_COUNTS = {
    "candidate_facet_complements": 824,
    "face_witness_constraints": 1324,
    "lanner_order_constraints": 7051,
    "ridge_constraints": 1547,
}
EXPECTED_CORE_COUNT = 41
EXPECTED_ORBIT_IMAGE_CLAUSE_COUNT = 427_140

Pair = tuple[int, int]
CoreSignature = tuple[tuple[int, ...], tuple[int, ...]]


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _digest(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _mask(values: Iterable[int]) -> int:
    return sum(1 << int(value) for value in values)


def _vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value & (1 << vertex))


def _pair(left: int, right: int) -> Pair:
    if left == right:
        raise ValueError("a Coxeter pair needs distinct vertices")
    return (left, right) if left < right else (right, left)


def _decode_graph6(text: str) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Decode the small-order graph6 format without a project helper."""

    if not text or ord(text[0]) - 63 >= 63:
        raise ValueError("this audit expects a small-order graph6 record")
    order = ord(text[0]) - 63
    bits = []
    for character in text[1:]:
        value = ord(character) - 63
        if not 0 <= value <= 63:
            raise ValueError("invalid graph6 character")
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))
    edges = []
    position = 0
    for right in range(1, order):
        for left in range(right):
            if bits[position]:
                edges.append((left, right))
            position += 1
    return order, tuple(edges)


def _contains_base(value: int) -> bool:
    return any(value & row == row for row in (*EXPECTED_EDGE_ROWS, FIXED_FOUR))


def _independent(value: int) -> bool:
    return all(value & edge != edge for edge in EXPECTED_EDGE_ROWS)


def _hits_all(value: int, clutter: Sequence[int]) -> bool:
    return all(value & row for row in clutter)


def _move(value: int, action: Sequence[int]) -> int:
    answer = 0
    remaining = int(value)
    while remaining:
        bit = remaining & -remaining
        vertex = bit.bit_length() - 1
        answer |= 1 << action[vertex]
        remaining ^= bit
    return answer


def _swap_action(left: int, right: int) -> tuple[int, ...]:
    action = list(range(N))
    action[left], action[right] = action[right], action[left]
    return tuple(action)


def _source_generators() -> tuple[tuple[int, ...], ...]:
    isolates = tuple(range(2, 13))
    inside = tuple(vertex for vertex in isolates if FIXED_FOUR & (1 << vertex))
    outside = tuple(vertex for vertex in isolates if not FIXED_FOUR & (1 << vertex))
    generators = [
        _swap_action(block[index], block[index + 1])
        for block in (inside, outside)
        for index in range(len(block) - 1)
    ]
    generators.append(_swap_action(0, 1))
    if len(generators) != 10:
        raise ArithmeticError("fixed-source generator count changed")
    return tuple(generators)


# Classical Lanner representatives.  Missing local pairs have label two.
_LANNER_REPRESENTATIVES: dict[
    int, tuple[tuple[tuple[tuple[int, int], int], ...], ...]
] = {
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


def _representative_code(
    order: int, representative: Sequence[tuple[tuple[int, int], int]]
) -> tuple[int, ...]:
    labels = {
        _pair(left, right): 2 for left, right in combinations(range(order), 2)
    }
    labels.update(representative)
    return tuple(labels[_pair(left, right)] for left, right in combinations(range(order), 2))


def _permuted_code(
    order: int, code: tuple[int, ...], permutation: tuple[int, ...]
) -> tuple[int, ...]:
    pairs = tuple(combinations(range(order), 2))
    source = dict(zip(pairs, code, strict=True))
    target = {
        _pair(permutation[left], permutation[right]): label
        for (left, right), label in source.items()
    }
    return tuple(target[_pair(left, right)] for left, right in pairs)


@lru_cache(maxsize=None)
def _lanner_codes(order: int) -> tuple[tuple[int, ...], ...]:
    if order not in (4, 5):
        raise ValueError("only order-four/order-five Lanner supports occur")
    codes = {
        _permuted_code(order, _representative_code(order, representative), permutation)
        for representative in _LANNER_REPRESENTATIVES[order]
        for permutation in permutations(range(order))
    }
    expected = {4: 108, 5: 420}[order]
    if len(codes) != expected:
        raise ArithmeticError(f"order-{order} Lanner placement count changed")
    return tuple(sorted(codes))


def _path_labels(
    component: frozenset[int],
    adjacency: Mapping[int, tuple[tuple[int, int], ...]],
) -> tuple[int, ...]:
    endpoint = min(vertex for vertex in component if len(adjacency[vertex]) == 1)
    labels = []
    previous: int | None = None
    current = endpoint
    while True:
        following = [
            (vertex, label)
            for vertex, label in adjacency[current]
            if vertex != previous
        ]
        if not following:
            break
        if len(following) != 1:
            raise ArithmeticError("finite path traversal branched")
        next_vertex, label = following[0]
        labels.append(label)
        previous, current = current, next_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def _branch_arms(
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


def _finite_connected(
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
    edge_count = sum(map(len, adjacency.values())) // 2
    if edge_count != order - 1:
        return False
    degree_histogram = Counter(len(adjacency[vertex]) for vertex in component)
    if max(degree_histogram) <= 2:
        sequence = _path_labels(component, adjacency)
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
    if degree_histogram[3] != 1 or any(degree > 3 for degree in degree_histogram):
        return False
    if any(
        labels[_pair(left, right)] != 3
        for left in component
        for right in component
        if left < right and labels[_pair(left, right)] > 2
    ):
        return False
    arms = _branch_arms(component, adjacency)
    return arms == (1, 1, order - 3) or arms in {
        (1, 2, 2),
        (1, 2, 3),
        (1, 2, 4),
    }


@lru_cache(maxsize=262_144)
def _finite_code(order: int, code: tuple[int, ...]) -> bool:
    if len(code) != order * (order - 1) // 2 or any(
        label not in (2, 3, 4, 5) for label in code
    ):
        return False
    pairs = tuple(combinations(range(order), 2))
    labels = dict(zip(pairs, code, strict=True))
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
                if labels[_pair(vertex, other)] > 2
            }
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        components.append(frozenset(component))
    return all(_finite_connected(component, labels) for component in components)


def _restriction_finite(
    restriction: Sequence[int], labels: Mapping[Pair, int]
) -> bool:
    ordered = tuple(sorted(restriction))
    code = tuple(
        labels[_pair(left, right)] for left, right in combinations(ordered, 2)
    )
    return _finite_code(len(ordered), code)


@lru_cache(maxsize=None)
def _placement_options(
    support: int,
) -> tuple[tuple[tuple[Pair, int], ...], ...]:
    support_vertices = _vertices(support)
    local_pairs = tuple(combinations(range(len(support_vertices)), 2))
    return tuple(
        tuple(
            (
                _pair(support_vertices[left], support_vertices[right]),
                label,
            )
            for (left, right), label in zip(local_pairs, code, strict=True)
        )
        for code in _lanner_codes(len(support_vertices))
    )


def _support_order(supports: Sequence[int]) -> tuple[int, ...]:
    remaining = set(range(len(supports)))
    union = 0
    answer = []
    while remaining:
        if not union:
            index = min(remaining)
        else:
            index = max(
                remaining,
                key=lambda candidate: (
                    (supports[candidate] & union).bit_count(),
                    -(supports[candidate] & ~union).bit_count(),
                    -candidate,
                ),
            )
        answer.append(index)
        union |= supports[index]
        remaining.remove(index)
    return tuple(answer)


def _label_state_count(
    supports: Sequence[int],
    restriction_sources: Sequence[int],
) -> tuple[int, tuple[dict[str, int], ...]]:
    """Exhaust all compatible Lanner placements for one local cut.

    ``restriction_sources`` may be full actual vertices or the stored face
    restrictions.  At each stage only their intersection with the current
    support union matters.  Every unassigned meeting pair is set to label two;
    failure there is monotone under increasing a label to 3, 4, or 5.
    """

    order = _support_order(supports)
    states: tuple[tuple[tuple[Pair, int], ...], ...] = ((),)
    union = 0
    stages = []
    for source_index in order:
        support = supports[source_index]
        union |= support
        cluster = _vertices(union)
        restrictions = tuple(
            sorted(
                {
                    tuple(vertex for vertex in _vertices(actual) if union & (1 << vertex))
                    for actual in restriction_sources
                },
                key=lambda row: (len(row), row),
            )
        )
        passing: set[tuple[tuple[Pair, int], ...]] = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in _placement_options(support):
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    _pair(left, right): 2
                    for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(
                    _restriction_finite(restriction, labels)
                    for restriction in restrictions
                ):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": source_index,
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
                "restriction_count": len(restrictions),
            }
        )
        if not states:
            break
    return len(states), tuple(stages)


def _audit_lanner_tables() -> dict[str, Any]:
    checks: dict[str, bool] = {}
    summaries = {}
    for order in (4, 5):
        codes = _lanner_codes(order)
        proper_checks = []
        full_checks = []
        pairs = tuple(combinations(range(order), 2))
        for code in codes:
            labels = dict(zip(pairs, code, strict=True))
            full_checks.append(not _finite_code(order, code))
            proper_checks.extend(
                _restriction_finite(subset, labels)
                for size in range(1, order)
                for subset in combinations(range(order), size)
            )
        summaries[str(order)] = {
            "classical_isomorphism_type_count": len(_LANNER_REPRESENTATIVES[order]),
            "labelled_placement_count": len(codes),
            "all_full_diagrams_nonelliptic": all(full_checks),
            "all_proper_subdiagrams_elliptic": all(proper_checks),
            "placement_code_sha256": _digest(codes),
        }
        checks[f"order_{order}_placement_count"] = len(codes) == {4: 108, 5: 420}[order]
        checks[f"order_{order}_isomorphism_count"] = len(
            _LANNER_REPRESENTATIVES[order]
        ) == {4: 9, 5: 5}[order]
        checks[f"order_{order}_minimal_nonelliptic"] = all(full_checks) and all(
            proper_checks
        )
    if not all(checks.values()):
        raise ArithmeticError(f"independent Lanner table audit failed: {checks}")
    return {"orders": summaries, "checks": checks}


def _load_inputs() -> tuple[dict[str, Any], dict[str, Any]]:
    if _sha256(CHECKPOINT) != CHECKPOINT_SHA256:
        raise ArithmeticError("terminal CEGAR checkpoint hash changed")
    if _sha256(SUPPORT) != SUPPORT_SHA256:
        raise ArithmeticError("H9/14 support schedule hash changed")
    if ortools.__version__ != PINNED_ORTOOLS_VERSION:
        raise ArithmeticError(
            f"OR-Tools version changed: {ortools.__version__}"
        )
    return (
        json.loads(CHECKPOINT.read_text(encoding="utf-8")),
        json.loads(SUPPORT.read_text(encoding="utf-8")),
    )


def _add_involution_lex_leader(
    model: cp_model.CpModel,
    chosen: Sequence[cp_model.IntVar],
    permutation: Sequence[int],
    tag: str,
) -> int:
    """Impose x <=lex g(x) for one involutive source generator."""

    if sorted(permutation) != list(range(len(chosen))):
        raise ArithmeticError("source generator does not permute facet variables")
    if any(
        permutation[permutation[position]] != position
        for position in range(len(chosen))
    ):
        raise ArithmeticError("source generator is not an involution")
    representatives = tuple(
        position
        for position, partner in enumerate(permutation)
        if position < partner
    )
    prefix = model.new_bool_var(f"{tag}_prefix_0")
    model.add(prefix == 1)
    for step, position in enumerate(representatives):
        partner = permutation[position]
        model.add(chosen[position] <= chosen[partner]).only_enforce_if(prefix)
        equal = model.new_bool_var(f"{tag}_equal_{step}")
        model.add(chosen[position] == chosen[partner]).only_enforce_if(equal)
        model.add(chosen[position] != chosen[partner]).only_enforce_if(equal.Not())
        following = model.new_bool_var(f"{tag}_prefix_{step + 1}")
        model.add(following <= prefix)
        model.add(following <= equal)
        model.add(following >= prefix + equal - 1)
        prefix = following
    return len(representatives)


def _build_terminal_base_model(
    checkpoint: dict[str, Any], support_document: dict[str, Any]
) -> dict[str, Any]:
    order, edges = _decode_graph6(GRAPH6)
    edge_rows = tuple(_mask(edge) for edge in edges)
    record = next(
        row for row in support_document["records"] if row["graph6"] == GRAPH6
    )
    schedule = next(
        row
        for row in record["dehn_sommerville_schedule"]
        if int(row["missing_four_count"]) == 1
    )
    scheduled_facets = int(schedule["f_vector"][-1])
    total_missing_fives = int(schedule["missing_five_count_offset"])

    five_ground = tuple(
        value
        for value in range(FULL + 1)
        if value.bit_count() == 5 and _independent(value)
    )
    shadow = tuple(
        five for five in five_ground if five & FIXED_FOUR == FIXED_FOUR
    )
    base = (*edge_rows, FIXED_FOUR, *shadow)
    facets = tuple(
        value
        for value in range(1, FULL + 1)
        if value.bit_count() == 5 and _hits_all(value, base)
    )

    model = cp_model.CpModel()
    live = [
        model.new_bool_var(f"facet_{index}") for index in range(len(facets))
    ]
    model.add(sum(live) == scheduled_facets)

    face_witness_constraint_count = 0
    for size in range(5):
        for chosen in combinations(range(N), size):
            value = _mask(chosen)
            if _contains_base(value):
                continue
            witnesses = [
                live[index]
                for index, facet in enumerate(facets)
                if not facet & value
            ]
            if witnesses:
                model.add_bool_or(witnesses)
            else:
                model.add(0 == 1)
            face_witness_constraint_count += 1

    missing_five = []
    missing_five_by_mask: dict[int, cp_model.IntVar] = {}
    for position, five in enumerate(five_ground):
        variable = model.new_bool_var(f"missing5_{position}")
        disjoint_facets = [
            live[index]
            for index, facet in enumerate(facets)
            if not facet & five
        ]
        if disjoint_facets:
            model.add(sum(disjoint_facets) == 0).only_enforce_if(variable)
            model.add_bool_or(disjoint_facets).only_enforce_if(variable.Not())
        else:
            model.add(variable == 1)
        missing_five.append(variable)
        missing_five_by_mask[five] = variable
    model.add(sum(missing_five) == total_missing_fives)

    lanner_order_constraint_count = 0
    for value in range(FULL + 1):
        if value.bit_count() < 6 or _contains_base(value):
            continue
        contained_missing = [
            missing_five_by_mask[_mask(chosen)]
            for chosen in combinations(_vertices(value), 5)
        ]
        facet_witnesses = [
            live[index]
            for index, facet in enumerate(facets)
            if not facet & value
        ]
        model.add_bool_or([*contained_missing, *facet_witnesses])
        lanner_order_constraint_count += 1

    ridge_constraint_count = 0
    for six in range(FULL + 1):
        if six.bit_count() != 6:
            continue
        members = [
            live[index]
            for index, facet in enumerate(facets)
            if facet & six == facet
        ]
        if not members:
            continue
        degree = model.new_int_var(0, len(members), f"ridge_{six}")
        model.add(degree == sum(members))
        model.add_allowed_assignments([degree], [(0,), (2,)])
        ridge_constraint_count += 1

    generators = _source_generators()
    facet_index = {facet: index for index, facet in enumerate(facets)}
    lex_moved_pairs = []
    generator_checks = []
    five_ground_set = set(five_ground)
    facet_set = set(facets)
    for generator_index, action in enumerate(generators):
        permutation = tuple(
            facet_index[_move(facet, action)] for facet in facets
        )
        lex_moved_pairs.append(
            _add_involution_lex_leader(
                model, live, permutation, f"lex_{generator_index}"
            )
        )
        generator_checks.append(
            _move(FIXED_FOUR, action) == FIXED_FOUR
            and {_move(edge, action) for edge in edge_rows} == set(edge_rows)
            and {_move(five, action) for five in five_ground} == five_ground_set
            and {_move(facet, action) for facet in facets} == facet_set
        )

    counts = {
        "candidate_facet_complements": len(facets),
        "face_witness_constraints": face_witness_constraint_count,
        "lanner_order_constraints": lanner_order_constraint_count,
        "ridge_constraints": ridge_constraint_count,
    }
    checks = {
        "graph6_order_is_14": order == N,
        "graph6_edges_are_exactly_0_13_and_1_13": edge_rows == EXPECTED_EDGE_ROWS,
        "checkpoint_identity_matches_fixed_quantifier": (
            checkpoint["record_index"] == 0
            and checkpoint["graph6"] == GRAPH6
            and checkpoint["schedule_a"] == 1
            and checkpoint["source_orbit_index"] == 3
            and checkpoint["representative_A4s"] == [list(FIXED_FOUR_VERTICES)]
        ),
        "schedule_has_178_vertices_and_32_missing_fives": (
            scheduled_facets == 178 and total_missing_fives == 32
        ),
        "fixed_A4_shadow_has_8_fives": len(shadow) == 8,
        "remaining_additional_A5_count_is_24": total_missing_fives - len(shadow)
        == 24,
        "base_model_counts_match_checkpoint": counts == checkpoint["dual_model_counts"]
        == EXPECTED_BASE_COUNTS,
        "ten_generators_preserve_the_full_fixed_source_model": (
            len(generators) == 10 and all(generator_checks)
        ),
        "generator_blocks_give_full_source_group_order_483840": (
            EXPECTED_SOURCE_GROUP_ORDER == checkpoint["source_group_order"] == 483_840
        ),
        "lex_slice_is_lossless": (
            len(generators) == 10
            and all(generator_checks)
            and EXPECTED_SOURCE_GROUP_ORDER == 483_840
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"independent terminal base model failed: {checks}")
    return {
        "model": model,
        "live": tuple(live),
        "facets": facets,
        "five_ground": five_ground,
        "missing_five": tuple(missing_five),
        "missing_five_by_mask": missing_five_by_mask,
        "generators": generators,
        "scheduled_facets": scheduled_facets,
        "total_missing_fives": total_missing_fives,
        "shadow": shadow,
        "counts": counts,
        "lex_moved_pair_counts": lex_moved_pairs,
        "checks": checks,
    }


def _move_signature(
    signature: CoreSignature, action: Sequence[int]
) -> CoreSignature:
    additions, restrictions = signature
    return (
        tuple(sorted(_move(value, action) for value in additions)),
        tuple(sorted(_move(value, action) for value in restrictions)),
    )


def _core_orbit(
    signature: CoreSignature,
    generators: Sequence[Sequence[int]],
) -> set[CoreSignature]:
    """Enumerate the complete generator orbit, not a sampled action list."""

    orbit = {signature}
    queue = deque([signature])
    while queue:
        current = queue.popleft()
        for generator in generators:
            image = _move_signature(current, generator)
            if image not in orbit:
                orbit.add(image)
                queue.append(image)
    return orbit


def _checkpoint_stage_projection(stages: Sequence[dict[str, Any]]) -> list[dict[str, int]]:
    return [
        {
            "support_index": int(stage["support_index"]),
            "compatible_extension_count": int(stage["compatible_extension_count"]),
            "passing_state_count": int(stage["passing_state_count"]),
            "restriction_count": int(stage["restriction_count"]),
        }
        for stage in stages
    ]


def _audit_one_core(
    row: dict[str, Any],
    base: dict[str, Any],
) -> tuple[dict[str, Any], set[CoreSignature]]:
    core = row.get("local_core")
    if not isinstance(core, dict):
        raise ArithmeticError("terminal checkpoint contains a non-core iteration")
    supports = tuple(_mask(face) for face in core["core_supports"])
    additions = tuple(sorted(_mask(face) for face in core["core_additional_A5s"]))
    actual_vertices = tuple(_mask(face) for face in core["core_actual_vertices"])
    restrictions = tuple(sorted(_mask(face) for face in core["core_face_restrictions"]))
    facet_complements = tuple(
        sorted(_mask(face) for face in core["core_facet_complements"])
    )
    support_union = 0
    for support in supports:
        support_union |= support

    terminal_from_actual, actual_stages = _label_state_count(
        supports, actual_vertices
    )
    terminal_from_cut, cut_stages = _label_state_count(
        supports, restrictions
    )
    projected_actual = tuple(
        sorted({actual & support_union for actual in actual_vertices})
    )
    projected_cut = tuple(
        sorted({restriction & support_union for restriction in restrictions})
    )

    expected_replay = _checkpoint_stage_projection(core["replay_stages"])
    actual_replay = _checkpoint_stage_projection(actual_stages)
    cut_replay = _checkpoint_stage_projection(cut_stages)
    five_ground = set(base["five_ground"])
    shadow = set(base["shadow"])
    local_checks = {
        "iteration_number_is_sequential": int(row["iteration"]) >= 0,
        "fixed_A4_is_present_once": supports.count(FIXED_FOUR) == 1,
        "all_other_supports_are_exactly_stored_additions": (
            tuple(sorted(value for value in supports if value != FIXED_FOUR))
            == additions
        ),
        "all_additions_are_independent_order_five": all(
            value.bit_count() == 5 and value in five_ground for value in additions
        ),
        "no_core_addition_is_forced_shadow": not (set(additions) & shadow),
        "all_actual_vertices_have_order_nine": all(
            value.bit_count() == 9 for value in actual_vertices
        ),
        "stored_facet_complements_match_actual_vertices": facet_complements
        == tuple(sorted(FULL ^ value for value in actual_vertices)),
        "stored_cut_restrictions_project_to_actual_restrictions": (
            projected_cut == projected_actual
        ),
        "actual_vertex_replay_matches_checkpoint_stages": actual_replay
        == expected_replay,
        "cut_semantics_replay_matches_checkpoint_stages": cut_replay
        == expected_replay,
        "actual_vertex_label_state_count_is_zero": terminal_from_actual == 0,
        "installed_cut_label_state_count_is_zero": terminal_from_cut == 0,
        "stored_terminal_state_count_is_zero": (
            int(core["terminal_sparse_label_state_count"]) == 0
        ),
        "stored_core_counts_match": (
            int(core["core_support_count"]) == len(supports)
            and int(core["core_actual_vertex_count"]) == len(actual_vertices)
        ),
    }
    if not all(local_checks.values()):
        raise ArithmeticError(
            f"local core {row['iteration']} failed independent replay: {local_checks}"
        )

    signature = (additions, restrictions)
    orbit_started = time.monotonic()
    orbit = _core_orbit(signature, base["generators"])
    orbit_seconds = time.monotonic() - orbit_started
    orbit_size = len(orbit)
    orbit_checks = {
        "orbit_size_matches_checkpoint": orbit_size
        == int(row["full_face_core_orbit_size"])
        == int(row["full_core_orbit_size"]),
        "orbit_size_divides_full_source_group": (
            EXPECTED_SOURCE_GROUP_ORDER % orbit_size == 0
        ),
        "stabilizer_order_matches_checkpoint": (
            EXPECTED_SOURCE_GROUP_ORDER // orbit_size
            == int(row["core_stabilizer_order"])
        ),
        "every_image_addition_remains_in_missing_five_ground": all(
            all(value in five_ground for value in image_additions)
            for image_additions, _ in orbit
        ),
        "every_image_clause_has_a_literal": all(
            bool(image_additions or image_restrictions)
            for image_additions, image_restrictions in orbit
        ),
    }
    if not all(orbit_checks.values()):
        raise ArithmeticError(
            f"core orbit {row['iteration']} failed: {orbit_checks}"
        )

    sorted_orbit = sorted(orbit)
    summary = {
        "iteration": int(row["iteration"]),
        "support_count": len(supports),
        "additional_A5_count": len(additions),
        "face_restriction_count": len(restrictions),
        "terminal_sparse_label_state_count": terminal_from_cut,
        "replay_stages": list(cut_stages),
        "full_source_core_orbit_size": orbit_size,
        "core_stabilizer_order": EXPECTED_SOURCE_GROUP_ORDER // orbit_size,
        "full_source_orbit_sha256": _digest(sorted_orbit),
        "orbit_rebuild_wall_time_seconds": orbit_seconds,
        "local_checks": local_checks,
        "orbit_checks": orbit_checks,
        "all_checks_pass": all(local_checks.values()) and all(orbit_checks.values()),
    }
    return summary, orbit


def _constraint_type_histogram(model: cp_model.CpModel) -> dict[str, int]:
    histogram: Counter[str] = Counter()
    kinds = (
        "all_diff", "at_most_one", "automaton", "bool_and", "bool_or",
        "bool_xor", "circuit", "cumulative", "dummy_constraint", "element",
        "exactly_one", "int_div", "int_mod", "int_prod", "interval",
        "inverse", "lin_max", "linear", "no_overlap", "no_overlap_2d",
        "reservoir", "routes", "table",
    )
    for constraint in model.proto.constraints:
        present = [
            kind for kind in kinds if getattr(constraint, f"has_{kind}")()
        ]
        if len(present) != 1:
            raise ArithmeticError(
                f"terminal constraint has {len(present)} active kinds"
            )
        histogram[present[0]] += 1
    return dict(sorted(histogram.items()))


def _model_proto_digest(model: cp_model.CpModel) -> tuple[str, int]:
    """Export the native deterministic protobuf and hash it without wrappers."""

    with tempfile.TemporaryDirectory(prefix="h914_terminal_proto_") as directory:
        path = Path(directory) / "terminal_model.pb"
        if not model.export_to_file(str(path)):
            raise ArithmeticError("OR-Tools failed to export the terminal model")
        return _sha256(path), path.stat().st_size


def _install_and_audit_cores(
    checkpoint: dict[str, Any], base: dict[str, Any]
) -> dict[str, Any]:
    model: cp_model.CpModel = base["model"]
    live: tuple[cp_model.IntVar, ...] = base["live"]
    facets: tuple[int, ...] = base["facets"]
    missing_by_mask: dict[int, cp_model.IntVar] = base["missing_five_by_mask"]
    iterations = checkpoint["iterations"]
    face_variables: dict[int, cp_model.IntVar] = {}

    def face_variable(restriction: int) -> cp_model.IntVar:
        known = face_variables.get(restriction)
        if known is not None:
            return known
        variable = model.new_bool_var(f"core_face_{restriction}")
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

    if len(iterations) != EXPECTED_CORE_COUNT:
        raise ArithmeticError("checkpoint does not contain exactly 41 cores")
    core_summaries = []
    installed_clause_count = 0
    cut_arity_histogram: Counter[int] = Counter()
    combined_orbit_digests = []
    for expected_iteration, row in enumerate(iterations):
        if int(row["iteration"]) != expected_iteration:
            raise ArithmeticError("checkpoint core iteration order changed")
        summary, orbit = _audit_one_core(row, base)
        core_summaries.append(summary)
        combined_orbit_digests.append(summary["full_source_orbit_sha256"])
        for additions, restrictions in sorted(orbit):
            literals = [
                *(missing_by_mask[value].Not() for value in additions),
                *(face_variable(value).Not() for value in restrictions),
            ]
            model.add_bool_or(literals)
            installed_clause_count += 1
            cut_arity_histogram[len(literals)] += 1

    checks = {
        "exactly_41_local_cores_replayed": len(core_summaries)
        == EXPECTED_CORE_COUNT,
        "every_local_core_is_a_complete_label_contradiction": all(
            row["terminal_sparse_label_state_count"] == 0
            and all(row["local_checks"].values())
            for row in core_summaries
        ),
        "every_full_source_orbit_rebuilt_exactly": all(
            all(row["orbit_checks"].values()) for row in core_summaries
        ),
        "installed_orbit_image_clause_count_is_427140": installed_clause_count
        == EXPECTED_ORBIT_IMAGE_CLAUSE_COUNT,
        "checkpoint_declares_41_excluded_core_orbits": (
            checkpoint["excluded_local_core_orbit_count"] == EXPECTED_CORE_COUNT
        ),
        "every_cut_is_lossless_for_compact_labels": (
            all(
                row["terminal_sparse_label_state_count"] == 0
                and all(row["local_checks"].values())
                and all(row["orbit_checks"].values())
                for row in core_summaries
            )
            and all(base["checks"].values())
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"terminal core installation failed: {checks}")
    proto_sha256, proto_size = _model_proto_digest(model)
    return {
        "model": model,
        "core_summaries": core_summaries,
        "face_variable_count": len(face_variables),
        "installed_clause_count": installed_clause_count,
        "cut_arity_histogram": {
            str(key): value for key, value in sorted(cut_arity_histogram.items())
        },
        "combined_core_orbit_digest_sha256": _digest(combined_orbit_digests),
        "terminal_model_proto_sha256": proto_sha256,
        "terminal_model_serialized_byte_count": proto_size,
        "terminal_model_variable_count": len(model.proto.variables),
        "terminal_model_constraint_count": len(model.proto.constraints),
        "terminal_model_constraint_type_histogram": _constraint_type_histogram(model),
        "checks": checks,
    }


def _solve_terminal(model: cp_model.CpModel) -> dict[str, Any]:
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = 1
    solver.parameters.random_seed = 0
    started = time.monotonic()
    status = solver.solve(model)
    external_seconds = time.monotonic() - started
    status_name = solver.status_name(status)
    response = solver.response_proto
    return {
        "status": status_name,
        "infeasible": status == cp_model.INFEASIBLE,
        "num_search_workers": 1,
        "time_limit_seconds": None,
        "random_seed": 0,
        "solver_wall_time_seconds": solver.wall_time,
        "external_wall_time_seconds": external_seconds,
        "deterministic_time": response.deterministic_time,
        "branch_count": solver.num_branches,
        "conflict_count": solver.num_conflicts,
        "integer_propagation_count": response.num_integer_propagations,
        "boolean_propagation_count": response.num_binary_propagations,
    }


def build_certificate(*, solve_terminal: bool = True) -> dict[str, Any]:
    total_started = time.monotonic()
    checkpoint, support = _load_inputs()
    lanner_audit = _audit_lanner_tables()
    base = _build_terminal_base_model(checkpoint, support)
    installed = _install_and_audit_cores(checkpoint, base)
    terminal = (
        _solve_terminal(installed["model"])
        if solve_terminal
        else {
            "status": "NOT_RUN",
            "infeasible": False,
            "num_search_workers": 1,
            "time_limit_seconds": None,
            "random_seed": 0,
        }
    )
    checks = {
        "pinned_checkpoint_hash_verified": _sha256(CHECKPOINT)
        == CHECKPOINT_SHA256,
        "pinned_support_hash_verified": _sha256(SUPPORT) == SUPPORT_SHA256,
        "pinned_ortools_version_verified": ortools.__version__
        == PINNED_ORTOOLS_VERSION,
        "independent_lanner_and_elliptic_tables_green": all(
            lanner_audit["checks"].values()
        ),
        "independent_base_model_green": all(base["checks"].values()),
        "all_local_cores_and_full_orbits_green": all(
            installed["checks"].values()
        ),
        "terminal_cp_sat_status_is_INFEASIBLE": (
            terminal["infeasible"] if solve_terminal else True
        ),
        "checkpoint_terminal_status_agrees": (
            checkpoint["terminal_solver_status"] == "INFEASIBLE"
            and checkpoint["stopped_reason"]
            == "terminal_infeasible_after_full_core_orbit_cuts"
        ),
        "finite_volume_noncompact_not_checked": (
            checkpoint["finite_volume_noncompact_checked"] is False
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"terminal CEGAR audit failed: {checks}")

    core_summaries = installed["core_summaries"]
    stable_terminal_model = {
        key: value
        for key, value in installed.items()
        if key not in {"model", "core_summaries"}
    }
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": all(checks.values()),
        "theorem": {
            "statement": (
                "There is no compact H9 Coxeter polytope with 14 facets whose "
                "missing-edge graph is Record 0, whose Dehn--Sommerville schedule "
                "has a=1, and whose unique missing A4 belongs to fixed source orbit "
                "3 represented by {2,3,4,13}."
            ),
            "quantifier": "one fixed Record-0/a=1 source orbit only",
            "record_index": 0,
            "graph6": GRAPH6,
            "schedule_a": 1,
            "source_orbit_index": 3,
            "representative_A4": list(FIXED_FOUR_VERTICES),
            "compact_only": True,
            "finite_volume_noncompact_checked": False,
            "not_global_H9_14": True,
        },
        "pinned_inputs": {
            "terminal_checkpoint": {
                "path": CHECKPOINT.relative_to(ROOT).as_posix(),
                "sha256": CHECKPOINT_SHA256,
            },
            "nerve_support_schedule": {
                "path": SUPPORT.relative_to(ROOT).as_posix(),
                "sha256": SUPPORT_SHA256,
            },
            "ortools_version": PINNED_ORTOOLS_VERSION,
        },
        "independence": {
            "exploratory_CEGAR_driver_imported": False,
            "exploratory_dual_builder_imported": False,
            "exploratory_core_extractor_imported": False,
            "method": (
                "direct graph6/schedule reconstruction, independent classical "
                "finite/Lanner classifier, generator-orbit BFS, fresh CP-SAT model"
            ),
        },
        "base_model": {
            "scheduled_rank_nine_vertex_count": base["scheduled_facets"],
            "total_missing_independent_five_count": base["total_missing_fives"],
            "forced_A5_shadow_count": len(base["shadow"]),
            "additional_A5_count": base["total_missing_fives"]
            - len(base["shadow"]),
            "counts": base["counts"],
            "source_group": {
                "structure": "S3 x S8 x C2",
                "order": EXPECTED_SOURCE_GROUP_ORDER,
                "generator_count": len(base["generators"]),
                "lex_generator_moved_pair_counts": base["lex_moved_pair_counts"],
                "lossless_lex_reason": (
                    "the global lexicographic minimum of every full source orbit "
                    "satisfies x <= g(x) for each installed generator"
                ),
            },
            "checks": base["checks"],
        },
        "local_label_logic": {
            "finite_label_alphabet": [2, 3, 4, 5],
            "unassigned_pair_rule": (
                "set to 2; increasing a cosine adjacency entry cannot turn a "
                "nonelliptic restriction elliptic"
            ),
            "lanner_tables": lanner_audit,
            "core_count": len(core_summaries),
            "core_summaries": core_summaries,
        },
        "terminal_model": stable_terminal_model,
        "terminal_solver": terminal,
        "total_replay_wall_time_seconds": time.monotonic() - total_started,
        "decision": {
            "fixed_source_compact_coxeter_no_go": bool(
                terminal["infeasible"] if solve_terminal else False
            ),
            "terminal_solver_status": terminal["status"],
            "excluded_local_core_orbit_count": len(core_summaries),
            "installed_full_source_core_image_count": installed[
                "installed_clause_count"
            ],
            "scope": (
                "Record 0, a=1, source orbit 3 only; compact incidence/labels. "
                "Finite-volume noncompact cases and all other H9/14 sources remain open."
            ),
        },
        "checks": checks,
    }


def _stable_for_comparison(payload: dict[str, Any]) -> dict[str, Any]:
    stable = deepcopy(payload)
    stable.pop("total_replay_wall_time_seconds", None)
    solver = stable.get("terminal_solver", {})
    for key in (
        "solver_wall_time_seconds",
        "external_wall_time_seconds",
        "deterministic_time",
        "branch_count",
        "conflict_count",
        "integer_propagation_count",
        "boolean_propagation_count",
    ):
        solver.pop(key, None)
    for row in stable.get("local_label_logic", {}).get("core_summaries", []):
        row.pop("orbit_rebuild_wall_time_seconds", None)
    return stable


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument(
        "--skip-solve",
        action="store_true",
        help="rebuild all proofs/model but do not claim the theorem",
    )
    args = parser.parse_args()
    payload = build_certificate(solve_terminal=not args.skip_solve)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if _stable_for_comparison(archived) != _stable_for_comparison(payload):
            raise ArithmeticError("archived terminal CEGAR audit changed")
        print(f"verified {args.output}")
    elif args.skip_solve:
        print(json.dumps(payload["terminal_model"], indent=2, sort_keys=True))
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    print(json.dumps(payload["decision"], sort_keys=True))
    print(f"payload_sha256={sha256(rendered.encode('utf-8')).hexdigest()}")


if __name__ == "__main__":
    main()
