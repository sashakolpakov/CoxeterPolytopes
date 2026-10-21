#!/usr/bin/env python3
"""No-project-import replay of the whole-shadow record-43 compact no-go.

The incidence part of this verifier is deliberately rebuilt from the six
hard-coded missing edges.  It does not import ``coxeter_polytopes`` and it
does not trust a SAT/SMT answer, a graph package, or a symmetry package.

The replay proves the following bounded statement.

* Dehn--Sommerville and the strict global g-vector condition leave the four
  schedules ``a=0,1,2,3``.
* The sum of the fourteen strict vertex-link ``g4`` entries excludes
  ``a=2,3``.
* All 690 possible missing four-faces reduce to 378 under independently
  reconstructed local-link intervals.  Twenty-eight survive the forced
  shadow gate in one exact S8 orbit.  For a representative, all 540
  vertex-feasible pairs of additional minimal five-faces have a transversal
  of order four, so ``a=1`` is empty.
* A lossless row/column census of the ``a=0`` branch has 1,073 enumeration
  paths after duplicate rows are removed.  Exact-mask deduplication leaves
  887 clutters, and canonicalization under row permutations and all eight
  isolate-column permutations leaves 394 states.  Exactly two states at
  every level have the scheduled 90 facet complements.

The first state is matched exactly to the frozen, principal-block-re-audited
selected-nerve theorem.  For
the second, this script independently verifies the S4 x S4 nerve
automorphisms, reconstructs all 20-by-20 distinguished-support Lanner
placement orbits, runs the finite-label and monotonicity gates, and checks
the exact index-two witness formulas without importing project code.

Nothing here decides whether a noncompact chamber has finite or infinite
volume.  In particular, no parabolic cusp or Vinberg endpoint-extension
audit is performed.
"""

from __future__ import annotations

from collections import Counter, deque
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations
import json
from math import comb
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import networkx as nx
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SELECTION_ARCHIVE = ROOT / "data/enumerations/h9_14_record43_selection.json"
COXETER_ARCHIVE = (
    ROOT / "data/enumerations/h9_14_selected_nerve_coxeter_no_go.json"
)
COXETER_RECORD43_SOURCE = (
    ROOT / "src/coxeter_polytopes/combinatorics/h9_14_record43.py"
)

VERTEX_COUNT = 14
POLYTOPE_DIMENSION = 9
FULL_MASK = (1 << VERTEX_COUNT) - 1
GRAPH6 = "M?????????O?o?w??"
EDGES = (
    (0, 11),
    (0, 12),
    (1, 12),
    (0, 13),
    (1, 13),
    (2, 13),
)
FIVES = (
    (3, 4, 5, 6, 7),
    (1, 2, 8, 9, 10),
    (2, 8, 9, 10, 11),
    (8, 9, 10, 11, 12),
)
SECOND_FIVES = (
    (1, 3, 4, 5, 6),
    (2, 7, 8, 9, 10),
    (3, 4, 5, 6, 11),
    (7, 8, 9, 10, 12),
)
ACTIVE = (0, 1, 2, 11, 12, 13)
ISOLATES = (3, 4, 5, 6, 7, 8, 9, 10)

EXPECTED_SELECTION_FILE_SHA256 = (
    "2c6463598ab4ae2073b54fc8a9345ab099acba0492c786ba163d0ad5d29d3c66"
)
EXPECTED_SELECTION_CERTIFICATE_SHA256 = (
    "af521ed4e7476c95945df9b770f36a343dcf3fc79dabed557aaab0fd97a671a3"
)
EXPECTED_COXETER_FILE_SHA256 = (
    "028bbc8e266145e59fa1ee9f9796a6e2c241f4a4f09637bd73f98bced9908629"
)
EXPECTED_COXETER_CERTIFICATE_SHA256 = (
    "c2dfcbfc03dc5dbae34938ab32f3a1c1111fd729d7844fd76d04dfc6cfa44de8"
)
EXPECTED_COXETER_RECORD43_PROJECTION_SHA256 = (
    "6437c7a8881def45ba1b5abb620b1c0574f1e831f34806b303336a6dd6a6f66b"
)
EXPECTED_COXETER_KIND = (
    "h9_14_seven_selected_nerve_compact_coxeter_no_go_reaudited"
)
EXPECTED_SUPPORT_SYMMETRY_ORDERS = [120]
EXPECTED_SCOPE_REAUDIT_STATUS = (
    "hostile_false_alarm_rejected_by_principal_block_proof"
)
EXPECTED_COXETER_RECORD43_SOURCE_SHA256 = (
    "49623375ad2f9a3aa19865e740719daceb93d399db1bd9e38ec21b1f9e4d2ed5"
)
EXPECTED_FACET_COMPLEMENT_SHA256 = (
    "e93bd0a0353dc03b18ba499454c08f7f267637a18ccb0e953d22d8faec66effd"
)
EXPECTED_PACKED_FACE_TABLE_SHA256 = (
    "32cb562b346841aa4c4a850a911a04da96d907cf6503ed6947f4960bd0212b51"
)
EXPECTED_SECOND_FACET_COMPLEMENT_SHA256 = (
    "7afc91cca45c5d5796350082aa06d63e0df593734b2fb3fa200c615b6a20faa9"
)
EXPECTED_SECOND_PACKED_FACE_TABLE_SHA256 = (
    "13c95733594fcb5f91ef74a7c261cbc75683f42283884696e574616107dd8b2b"
)

EXPECTED_SUMMARY_SHA256 = "3193e0570170f2e4896b27f5404bee1d23f2240c47dd6ed62f42830254223d61"


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1)


EDGE_MASKS = tuple(_mask(edge) for edge in EDGES)
FIVE_MASKS = tuple(_mask(face) for face in FIVES)


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        ).encode("utf-8")
    ).hexdigest()


def _is_independent(value: int) -> bool:
    return all(value & edge != edge for edge in EDGE_MASKS)


@lru_cache(maxsize=None)
def _independent_masks(size: int) -> tuple[int, ...]:
    return tuple(
        value
        for chosen in combinations(range(VERTEX_COUNT), size)
        for value in (_mask(chosen),)
        if _is_independent(value)
    )


def _macaulay_raise(value: int, degree: int) -> int:
    if value < 0 or degree < 1:
        raise ValueError("invalid Macaulay raising input")
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
        raise AssertionError("Macaulay expansion did not terminate")
    return raised


def _h_prefix(missing_four_count: int) -> tuple[int, ...]:
    """Return h0,...,h4 from f[-1],f0,...,f3."""

    partial_f = (1, 14, 85, 300, 690 - missing_four_count)
    return tuple(
        sum(
            (-1) ** (index - order)
            * comb(POLYTOPE_DIMENSION - order, POLYTOPE_DIMENSION - index)
            * partial_f[order]
            for order in range(index + 1)
        )
        for index in range(5)
    )


def _schedule(missing_four_count: int) -> dict[str, object]:
    prefix = _h_prefix(missing_four_count)
    h_vector = prefix + tuple(reversed(prefix))
    f_vector = tuple(
        sum(
            comb(POLYTOPE_DIMENSION - order, POLYTOPE_DIMENSION - index)
            * h_vector[order]
            for order in range(index + 1)
        )
        for index in range(1, POLYTOPE_DIMENSION + 1)
    )
    independent_fives = len(_independent_masks(5))
    return {
        "missing_four_count": missing_four_count,
        "global_g4": prefix[4] - prefix[3],
        "h_vector": list(h_vector),
        "f_vector": list(f_vector),
        "missing_independent_five_count": independent_fives - f_vector[4],
        "summed_vertex_link_g4": _summed_vertex_link_g4(f_vector),
    }


def _summed_vertex_link_g4(f_vector: Sequence[int]) -> int:
    """Double-count the sum of g4 over the fourteen vertex links."""

    summed_link_f = {0: VERTEX_COUNT}
    for order in range(1, 5):
        summed_link_f[order] = (order + 1) * int(f_vector[order])

    def h_entry(index: int) -> int:
        return sum(
            (-1) ** (index - order)
            * comb(8 - order, 8 - index)
            * summed_link_f[order]
            for order in range(index + 1)
        )

    return h_entry(4) - h_entry(3)


def _enumerate_schedules() -> tuple[dict[str, object], ...]:
    schedules = []
    for missing_four_count in range(len(_independent_masks(4)) + 1):
        if _h_prefix(missing_four_count)[4] - _h_prefix(missing_four_count)[3] < 1:
            break
        schedules.append(_schedule(missing_four_count))
    if tuple(row["missing_four_count"] for row in schedules) != (0, 1, 2, 3):
        raise AssertionError("record-43 raw Dehn--Sommerville schedules changed")
    if tuple(row["summed_vertex_link_g4"] for row in schedules) != (24, 18, 12, 6):
        raise AssertionError("record-43 strict link-g4 sums changed")
    return tuple(schedules)


def _strict_interval(
    previous_g: int,
    next_g_without_missing_four: int,
    candidate_count: int,
    macaulay_degree: int,
) -> tuple[int, int] | None:
    lower = max(
        0,
        next_g_without_missing_four
        - _macaulay_raise(previous_g, macaulay_degree),
    )
    upper = min(candidate_count, next_g_without_missing_four - 1)
    return (lower, upper) if lower <= upper else None


def _meeting(first: int, second: int) -> bool:
    pair = (1 << first) | (1 << second)
    return all(pair != edge for edge in EDGE_MASKS)


def _local_link_constraints() -> tuple[
    tuple[tuple[int, tuple[int, int]], ...],
    tuple[tuple[int, tuple[int, int]], ...],
    tuple[tuple[int, int], ...],
]:
    """Reconstruct the vertex, edge, and triple missing-A4 intervals."""

    vertex_constraints: list[tuple[int, tuple[int, int]]] = []
    for vertex in range(VERTEX_COUNT):
        link = tuple(
            other
            for other in range(VERTEX_COUNT)
            if other != vertex and _meeting(vertex, other)
        )
        inherited_edges = sum(
            not _meeting(first, second)
            for first, second in combinations(link, 2)
        )
        f1 = comb(len(link), 2) - inherited_edges
        candidate_f2 = sum(
            _is_independent(_mask(triple))
            for triple in combinations(link, 3)
        )
        g1 = len(link) - 9
        g2 = f1 - 8 * len(link) + 36
        beta = candidate_f2 - 7 * f1 + 28 * len(link) - 84
        if not (g1 > 0 and 1 <= g2 <= _macaulay_raise(g1, 1)):
            raise AssertionError("record-43 vertex link fails its graph-only gate")
        interval = _strict_interval(g2, beta, candidate_f2, 2)
        if interval is None:
            raise AssertionError("record-43 vertex link has no A4 interval")
        vertex_constraints.append((1 << vertex, interval))

    edge_constraints: list[tuple[int, tuple[int, int]]] = []
    for first, second in combinations(range(VERTEX_COUNT), 2):
        if not _meeting(first, second):
            continue
        link = tuple(
            other
            for other in range(VERTEX_COUNT)
            if other not in (first, second)
            and _meeting(other, first)
            and _meeting(other, second)
        )
        candidate_f1 = sum(
            _meeting(one, two) for one, two in combinations(link, 2)
        )
        g1 = len(link) - 8
        beta = candidate_f1 - 7 * len(link) + 28
        interval = (
            _strict_interval(g1, beta, candidate_f1, 1) if g1 > 0 else None
        )
        if interval is None:
            raise AssertionError("record-43 edge link has no A4 interval")
        edge_constraints.append(((1 << first) | (1 << second), interval))

    triple_constraints: list[tuple[int, int]] = []
    for triple in combinations(range(VERTEX_COUNT), 3):
        triple_mask = _mask(triple)
        if not _is_independent(triple_mask):
            continue
        candidate_vertices = sum(
            all(_meeting(other, endpoint) for endpoint in triple)
            for other in range(VERTEX_COUNT)
            if other not in triple
        )
        maximum = candidate_vertices - 8
        if maximum < 0:
            raise AssertionError("record-43 triple link is already impossible")
        triple_constraints.append((triple_mask, maximum))

    return (
        tuple(vertex_constraints),
        tuple(edge_constraints),
        tuple(triple_constraints),
    )


def _a4_local_feasible(
    face: int,
    constraints: tuple[
        tuple[tuple[int, tuple[int, int]], ...],
        tuple[tuple[int, tuple[int, int]], ...],
        tuple[tuple[int, int], ...],
    ],
) -> bool:
    vertex_constraints, edge_constraints, triple_constraints = constraints
    for subset, interval in vertex_constraints:
        incidence = int(bool(face & subset))
        if not interval[0] <= incidence <= interval[1]:
            return False
    for subset, interval in edge_constraints:
        incidence = int(face & subset == subset)
        if not interval[0] <= incidence <= interval[1]:
            return False
    return all(
        int(face & subset == subset) <= maximum
        for subset, maximum in triple_constraints
    )


def _base_vertex_link_g3_g4() -> tuple[tuple[int, ...], tuple[int, ...]]:
    independent = {size: _independent_masks(size) for size in range(2, 6)}
    g3_values = []
    g4_values = []
    for vertex in range(VERTEX_COUNT):
        counts = {
            size: sum(bool(face & (1 << vertex)) for face in independent[size])
            for size in independent
        }
        f0, f1, f2, f3 = (counts[size] for size in range(2, 6))
        g3_values.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4_values.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    return tuple(g3_values), tuple(g4_values)


def _a1_shadow_gate(
    missing_four: int,
    base_g3: Sequence[int],
    base_g4: Sequence[int],
) -> tuple[bool, bool, int, int, int | None]:
    shadow = tuple(
        face for face in _independent_masks(5) if face & missing_four == missing_four
    )
    additional = 10 - len(shadow)
    if additional < 0:
        return False, False, len(shadow), additional, None
    lower_sum = 0
    upper_sum = 0
    first_strict_failure: int | None = None
    for vertex in range(VERTEX_COUNT):
        incidence = int(bool(missing_four & (1 << vertex)))
        previous = base_g3[vertex] - incidence
        shadow_incidence = sum(bool(face & (1 << vertex)) for face in shadow)
        before_additional_fives = (
            base_g4[vertex] + 6 * incidence - shadow_incidence
        )
        if previous < 1 or before_additional_fives < 1:
            if first_strict_failure is None:
                first_strict_failure = vertex
        if previous >= 0:
            lower_sum += max(
                0,
                before_additional_fives - _macaulay_raise(previous, 3),
            )
        upper_sum += before_additional_fives - 1
    strict = first_strict_failure is None
    aggregate = strict and lower_sum <= 5 * additional <= upper_sum
    return strict, aggregate, len(shadow), additional, first_strict_failure


def _audit_a1() -> dict[str, object]:
    independent_fours = _independent_masks(4)
    if len(independent_fours) != 690:
        raise AssertionError("record-43 independent A4 count changed")
    constraints = _local_link_constraints()
    local = tuple(face for face in independent_fours if _a4_local_feasible(face, constraints))
    if len(local) != 378:
        raise AssertionError("record-43 local A4 count changed")
    base_g3, base_g4 = _base_vertex_link_g3_g4()
    expected_base = (1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 2, 1)
    if base_g3 != expected_base or base_g4 != expected_base:
        raise AssertionError("record-43 base vertex-link g-vectors changed")

    rows = tuple(_a1_shadow_gate(face, base_g3, base_g4) for face in local)
    strict_count = sum(row[0] for row in rows)
    aggregate_count = sum(row[1] for row in rows)
    viable = tuple(
        face for face, row in zip(local, rows, strict=True) if row[1]
    )
    expected_viable = tuple(
        _mask((2, 11, first, second))
        for first, second in combinations(ISOLATES, 2)
    )
    if strict_count != 28 or aggregate_count != 28 or viable != expected_viable:
        raise AssertionError("record-43 viable a=1 orbit changed")

    representative = _mask((2, 3, 4, 11))
    shadow = tuple(
        face
        for face in _independent_masks(5)
        if face & representative == representative
    )
    candidates = tuple(
        face
        for face in _independent_masks(5)
        if face & representative != representative
    )
    lower = []
    upper = []
    for vertex in range(VERTEX_COUNT):
        incidence = int(bool(representative & (1 << vertex)))
        previous = base_g3[vertex] - incidence
        shadow_incidence = sum(
            bool(face & (1 << vertex)) for face in shadow
        )
        before_additional = (
            base_g4[vertex] + 6 * incidence - shadow_incidence
        )
        lower.append(max(0, before_additional - _macaulay_raise(previous, 3)))
        upper.append(before_additional - 1)
    feasible_pairs = tuple(
        (first, second)
        for first, second in combinations(candidates, 2)
        if all(
            lower[vertex]
            <= int(bool(first & (1 << vertex)))
            + int(bool(second & (1 << vertex)))
            <= upper[vertex]
            for vertex in range(VERTEX_COUNT)
        )
    )
    transversal_sizes = []
    for first, second in feasible_pairs:
        clutter = (*EDGE_MASKS, representative, first, second)
        size = next(
            order
            for order in range(1, 5)
            if any(
                all(_mask(chosen) & hyperedge for hyperedge in clutter)
                for chosen in combinations(range(VERTEX_COUNT), order)
            )
        )
        transversal_sizes.append(size)
    if (
        len(shadow) != 8
        or len(candidates) != 1084
        or comb(len(candidates), 2) != 586_986
        or len(feasible_pairs) != 540
        or Counter(transversal_sizes) != {4: 540}
    ):
        raise AssertionError("record-43 a=1 addition-pair census changed")
    return {
        "independent_missing_four_count": len(independent_fours),
        "local_interval_passing_count": len(local),
        "strict_shadow_passing_count": strict_count,
        "aggregate_macaulay_passing_count": aggregate_count,
        "shadow_size_histogram": {
            str(key): value for key, value in sorted(Counter(row[2] for row in rows).items())
        },
        "additional_minimal_five_count_histogram": {
            str(key): value for key, value in sorted(Counter(row[3] for row in rows).items())
        },
        "first_strict_failure_vertex_histogram": {
            str(key): value
            for key, value in sorted(
                Counter(row[4] for row in rows).items(),
                key=lambda item: str(item[0]),
            )
        },
        "viable_missing_four_orbit": [
            list(_vertices(face)) for face in viable
        ],
        "verified_s8_orbit_size": len(viable),
        "representative_missing_four": list(_vertices(representative)),
        "representative_forced_shadow_count": len(shadow),
        "representative_additional_a5_candidate_count": len(candidates),
        "representative_raw_addition_pair_count": 586_986,
        "representative_vertex_feasible_pair_count": len(feasible_pairs),
        "representative_transversal_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(transversal_sizes).items())
        },
        "surviving_addition_pair_count": 0,
    }


def _minimal_graph_covers() -> tuple[int, ...]:
    covers = []
    for candidate in range(1 << VERTEX_COUNT):
        if not all(candidate & edge for edge in EDGE_MASKS):
            continue
        if all(
            not all((candidate ^ bit) & edge for edge in EDGE_MASKS)
            for bit in (
                1 << vertex
                for vertex in range(VERTEX_COUNT)
                if candidate & (1 << vertex)
            )
        ):
            covers.append(candidate)
    result = tuple(covers)
    expected = (
        _mask((0, 1, 2)),
        _mask((0, 1, 13)),
        _mask((0, 12, 13)),
        _mask((11, 12, 13)),
    )
    if result != expected:
        raise AssertionError("record-43 minimal graph covers changed")
    return result


@lru_cache(maxsize=None)
def _column_distributions(
    allowed: tuple[int, ...], row_sums: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    counts = [0] * len(allowed)
    answers: list[tuple[int, ...]] = []

    def recurse(position: int, columns_left: int, sums: tuple[int, ...]) -> None:
        if position == len(allowed):
            if columns_left == 0 and not any(sums):
                answers.append(tuple(counts))
            return
        pattern = allowed[position]
        maximum = columns_left
        for row in range(4):
            if pattern & (1 << row):
                maximum = min(maximum, sums[row])
        for multiplicity in range(maximum + 1):
            next_left = columns_left - multiplicity
            next_sums = tuple(
                sums[row] - multiplicity * bool(pattern & (1 << row))
                for row in range(4)
            )
            if all(0 <= value <= next_left for value in next_sums):
                counts[position] = multiplicity
                recurse(position + 1, next_left, next_sums)
        counts[position] = 0

    recurse(0, len(ISOLATES), row_sums)
    return tuple(answers)


def _has_small_transversal(hyperedges: Sequence[int]) -> bool:
    return any(
        all(candidate & edge for edge in hyperedges)
        for size in range(1, 5)
        for chosen in combinations(range(VERTEX_COUNT), size)
        for candidate in (_mask(chosen),)
    )


def _canonical_a5_state(faces: Sequence[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Canonicalize four A5 rows under ``S4`` rows and ``S8`` isolates.

    For a fixed row order, an isolate is completely described by its four-bit
    incidence column.  Sorting the eight column words is exactly the quotient
    by arbitrary isolate permutations.  Taking the least result over all 24
    row orders then gives the full ``S4 x S8`` canonical key, without relying
    on a graph-isomorphism or group-action package.
    """

    active_mask = _mask(ACTIVE)
    candidates = []
    for ordered in permutations(tuple(faces)):
        active_rows = tuple(face & active_mask for face in ordered)
        columns = tuple(
            sorted(
                sum(bool(face & (1 << isolate)) << row for row, face in enumerate(ordered))
                for isolate in ISOLATES
            )
        )
        candidates.append((active_rows, columns))
    return min(candidates)


def _active_graph_automorphisms() -> tuple[tuple[int, ...], ...]:
    edge_set = {frozenset(edge) for edge in EDGES}
    answers = []
    for image_order in permutations(ACTIVE):
        action = dict(zip(ACTIVE, image_order, strict=True))
        if {
            frozenset((action[left], action[right])) for left, right in EDGES
        } == edge_set:
            answers.append(image_order)
    if len(answers) != 2:
        raise AssertionError("record-43 active graph automorphism count changed")
    return tuple(answers)


def _full_graph_target_key(faces: Sequence[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    keys = []
    for image_order in _active_graph_automorphisms():
        action = dict(zip(ACTIVE, image_order, strict=True))
        transformed = tuple(
            _mask(action.get(vertex, vertex) for vertex in _vertices(face))
            for face in faces
        )
        keys.append(_canonical_a5_state(transformed))
    return min(keys)


def _a0_census() -> tuple[dict[str, object], tuple[int, ...]]:
    covers = _minimal_graph_covers()
    _, base_g4 = _base_vertex_link_g3_g4()
    caps = tuple(value - 1 for value in base_g4)
    active_mask = _mask(ACTIVE)
    active_edges = tuple(edge for edge in EDGE_MASKS if edge & active_mask == edge)
    independent_active = tuple(
        value
        for value in range(1 << VERTEX_COUNT)
        if not value & ~active_mask
        and all(value & edge != edge for edge in active_edges)
        and all(bool(value & (1 << vertex)) <= caps[vertex] for vertex in ACTIVE)
    )
    if len(independent_active) != 12:
        raise AssertionError("record-43 active row alphabet changed")

    active_configs: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for config in combinations_with_replacement(independent_active, 4):
        if any(
            sum(bool(row & (1 << vertex)) for row in config) > caps[vertex]
            for vertex in ACTIVE
        ):
            continue
        missed_rows = tuple(
            sum(
                1 << row_index
                for row_index, row in enumerate(config)
                if not row & cover
            )
            for cover in covers
        )
        if any(not rows for rows in missed_rows):
            continue
        if any(
            all(
                not (rows & (1 << row_index))
                or config[row_index] & (1 << vertex)
                for row_index in range(4)
            )
            for rows in missed_rows
            for vertex in ACTIVE
        ):
            continue
        active_configs.append((config, missed_rows))
    if len(active_configs) != 121:
        raise AssertionError("record-43 active row configurations changed")

    graph_cover_fives = tuple(
        value
        for chosen in combinations(range(VERTEX_COUNT), 5)
        for value in (_mask(chosen),)
        if all(value & edge for edge in EDGE_MASKS)
    )
    path_histogram: Counter[int] = Counter()
    exact_states: dict[tuple[int, ...], int] = {}
    canonical_states: dict[tuple[tuple[int, ...], tuple[int, ...]], int] = {}
    duplicate_row_states = 0
    enumeration_paths = 0
    small_transversal_paths = 0
    target_paths: list[tuple[int, ...]] = []
    for config, missed_rows in active_configs:
        allowed = tuple(
            column
            for column in range(1 << 4)
            if column.bit_count() <= 3
            if all(column & rows != rows for rows in missed_rows)
        )
        row_sums = tuple(5 - row.bit_count() for row in config)
        for multiplicities in _column_distributions(allowed, row_sums):
            columns = tuple(
                column
                for column, multiplicity in zip(allowed, multiplicities, strict=True)
                for _ in range(multiplicity)
            )
            if len(columns) != len(ISOLATES):
                raise AssertionError("record-43 isolate column total changed")
            faces = list(config)
            for isolate, column in zip(ISOLATES, columns, strict=True):
                for row in range(4):
                    if column & (1 << row):
                        faces[row] |= 1 << isolate
            face_tuple = tuple(faces)
            if any(face.bit_count() != 5 for face in face_tuple):
                raise AssertionError("record-43 A5 row size changed")
            if len(set(face_tuple)) != 4:
                duplicate_row_states += 1
                continue
            enumeration_paths += 1
            clutter = (*EDGE_MASKS, *face_tuple)
            if _has_small_transversal(clutter):
                small_transversal_paths += 1
                continue
            facet_count = sum(
                all(transversal & face for face in face_tuple)
                for transversal in graph_cover_fives
            )
            path_histogram[facet_count] += 1
            exact_key = tuple(sorted(face_tuple))
            previous = exact_states.setdefault(exact_key, facet_count)
            if previous != facet_count:
                raise AssertionError("one exact record-43 clutter has two facet counts")
            canonical_key = _canonical_a5_state(exact_key)
            previous = canonical_states.setdefault(canonical_key, facet_count)
            if previous != facet_count:
                raise AssertionError("one record-43 quotient state has two facet counts")
            if facet_count == 90:
                target_paths.append(exact_key)

    expected_histogram = {
        58: 1,
        60: 27,
        61: 26,
        62: 20,
        63: 8,
        64: 155,
        65: 24,
        66: 64,
        67: 76,
        68: 206,
        69: 20,
        70: 161,
        71: 30,
        72: 59,
        73: 73,
        74: 26,
        75: 16,
        76: 7,
        77: 30,
        79: 4,
        80: 21,
        81: 2,
        82: 2,
        83: 12,
        86: 1,
        90: 2,
    }
    expected_exact_histogram = {
        58: 1,
        60: 27,
        61: 26,
        62: 20,
        63: 6,
        64: 126,
        65: 24,
        66: 64,
        67: 52,
        68: 163,
        69: 19,
        70: 125,
        71: 28,
        72: 41,
        73: 57,
        74: 19,
        75: 16,
        76: 5,
        77: 30,
        79: 4,
        80: 15,
        81: 2,
        82: 2,
        83: 12,
        86: 1,
        90: 2,
    }
    expected_canonical_histogram = {
        58: 1,
        60: 6,
        61: 8,
        62: 7,
        63: 4,
        64: 32,
        65: 10,
        66: 22,
        67: 24,
        68: 41,
        69: 9,
        70: 55,
        71: 16,
        72: 18,
        73: 43,
        74: 13,
        75: 14,
        76: 5,
        77: 30,
        79: 2,
        80: 15,
        81: 2,
        82: 2,
        83: 12,
        86: 1,
        90: 2,
    }
    exact_histogram = Counter(exact_states.values())
    canonical_histogram = Counter(canonical_states.values())
    if enumeration_paths != 1073 or small_transversal_paths != 0:
        raise AssertionError("record-43 a=0 enumeration-path census changed")
    if len(exact_states) != 887 or len(canonical_states) != 394:
        raise AssertionError("record-43 a=0 deduplicated census changed")
    if dict(sorted(path_histogram.items())) != expected_histogram:
        raise AssertionError("record-43 path facet-count histogram changed")
    if dict(sorted(exact_histogram.items())) != expected_exact_histogram:
        raise AssertionError("record-43 exact-mask facet-count histogram changed")
    if dict(sorted(canonical_histogram.items())) != expected_canonical_histogram:
        raise AssertionError("record-43 canonical facet-count histogram changed")
    if len(target_paths) != 2:
        raise AssertionError("record-43 scheduled target count changed")
    if sum(count == 90 for count in exact_states.values()) != 2:
        raise AssertionError("record-43 exact-mask target count changed")
    if sum(count == 90 for count in canonical_states.values()) != 2:
        raise AssertionError("record-43 canonical target count changed")
    targets = tuple(dict.fromkeys(target_paths))
    expected_targets = (
        tuple(sorted(FIVE_MASKS)),
        tuple(sorted(_mask(face) for face in SECOND_FIVES)),
    )
    if targets != expected_targets:
        raise AssertionError("record-43 target nerves changed")
    full_target_orbit_count = len({_full_graph_target_key(target) for target in targets})
    if full_target_orbit_count != 2:
        raise AssertionError("record-43 target full-graph orbit count changed")
    return (
        {
            "active_row_alphabet_count": len(independent_active),
            "active_configuration_count": len(active_configs),
            "column_signature_count": _column_distributions.cache_info().currsize,
            "duplicate_row_enumeration_path_count": duplicate_row_states,
            "enumeration_path_count": enumeration_paths,
            "small_transversal_enumeration_path_count": small_transversal_paths,
            "exact_mask_clutter_count": len(exact_states),
            "s4_s8_canonical_state_count": len(canonical_states),
            "active_graph_automorphism_count": len(_active_graph_automorphisms()),
            "scheduled_target_full_graph_orbit_count": full_target_orbit_count,
            "active_graph_automorphism_quotient_applied": False,
            "enumeration_path_facet_count_histogram": {
                str(key): value for key, value in sorted(path_histogram.items())
            },
            "exact_mask_facet_count_histogram": {
                str(key): value for key, value in sorted(exact_histogram.items())
            },
            "s4_s8_canonical_facet_count_histogram": {
                str(key): value for key, value in sorted(canonical_histogram.items())
            },
            "scheduled_target_enumeration_path_count": len(target_paths),
            "scheduled_target_exact_mask_count": sum(
                count == 90 for count in exact_states.values()
            ),
            "scheduled_target_s4_s8_canonical_count": sum(
                count == 90 for count in canonical_states.values()
            ),
            "target_missing_five_states": [
                [list(_vertices(face)) for face in target]
                for target in targets
            ],
        },
        targets,
    )


def _minimal_transversals(hyperedges: Sequence[int]) -> tuple[int, ...]:
    result = []
    for candidate in range(1, 1 << VERTEX_COUNT):
        if not all(candidate & edge for edge in hyperedges):
            continue
        if all(
            not all((candidate ^ bit) & edge for edge in hyperedges)
            for bit in (
                1 << vertex
                for vertex in range(VERTEX_COUNT)
                if candidate & (1 << vertex)
            )
        ):
            result.append(candidate)
    return tuple(result)


def _face_table(facet_complements: Sequence[int]) -> tuple[bool, ...]:
    table = [False] * (1 << VERTEX_COUNT)
    table[0] = True
    for complement in facet_complements:
        facet = FULL_MASK ^ complement
        subset = facet
        while True:
            table[subset] = True
            if not subset:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def _ridge_groups(facet_complements: Sequence[int]) -> dict[int, tuple[int, ...]]:
    groups: dict[int, list[int]] = {}
    for complement in facet_complements:
        facet = FULL_MASK ^ complement
        remaining = facet
        while remaining:
            bit = remaining & -remaining
            groups.setdefault(complement | bit, []).append(complement)
            remaining ^= bit
    return {key: tuple(sorted(value)) for key, value in groups.items()}


def _dual_component_count(facet_complements: Sequence[int]) -> int:
    adjacency = {complement: set() for complement in facet_complements}
    for group in _ridge_groups(facet_complements).values():
        if len(group) == 2:
            first, second = group
            adjacency[first].add(second)
            adjacency[second].add(first)
    unseen = set(facet_complements)
    components = 0
    while unseen:
        components += 1
        root = min(unseen)
        unseen.remove(root)
        queue = deque([root])
        while queue:
            current = queue.popleft()
            for neighbor in adjacency[current] & unseen:
                unseen.remove(neighbor)
                queue.append(neighbor)
    return components


def _minimal_nonfaces(face_table: Sequence[bool]) -> tuple[int, ...]:
    result = []
    for candidate in range(1, 1 << VERTEX_COUNT):
        if face_table[candidate]:
            continue
        if all(
            face_table[candidate ^ (1 << vertex)]
            for vertex in range(VERTEX_COUNT)
            if candidate & (1 << vertex)
        ):
            result.append(candidate)
    return tuple(result)


def _audit_target(
    target_fives: Sequence[int],
    expected_facet_digest: str,
    expected_packed_digest: str,
) -> tuple[dict[str, object], tuple[int, ...]]:
    clutter = (*EDGE_MASKS, *target_fives)
    blocker = _minimal_transversals(clutter)
    blocker_histogram = Counter(value.bit_count() for value in blocker)
    if blocker_histogram != {5: 90}:
        raise AssertionError("record-43 target blocker is not pure 5^90")
    facet_payload = b"".join(value.to_bytes(2, "little") for value in blocker)
    facet_digest = sha256(facet_payload).hexdigest()
    if facet_digest != expected_facet_digest:
        raise AssertionError("record-43 facet-complement digest changed")

    face_table = _face_table(blocker)
    packed = bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )
    packed_digest = sha256(packed).hexdigest()
    if packed_digest != expected_packed_digest:
        raise AssertionError("record-43 packed face-table digest changed")
    face_vector = tuple(
        sum(face_table[value] for value in range(1 << VERTEX_COUNT) if value.bit_count() == size)
        for size in range(1, POLYTOPE_DIMENSION + 1)
    )
    expected_vector = (14, 85, 300, 690, 1088, 1190, 880, 405, 90)
    if face_vector != expected_vector:
        raise AssertionError("record-43 face vector changed")
    minimal_nonfaces = _minimal_nonfaces(face_table)
    if set(minimal_nonfaces) != set(clutter):
        raise AssertionError("record-43 generated complex has different minimal nonfaces")
    ridge_histogram = Counter(len(group) for group in _ridge_groups(blocker).values())
    if ridge_histogram != {2: 405}:
        raise AssertionError("record-43 ridge pseudomanifold audit changed")
    dual_components = _dual_component_count(blocker)
    if dual_components != 1:
        raise AssertionError("record-43 dual graph disconnected")
    return (
        {
            "minimal_transversal_size_histogram": {
                str(key): value for key, value in sorted(blocker_histogram.items())
            },
            "facet_complement_sha256": facet_digest,
            "packed_face_table_sha256": packed_digest,
            "face_table_true_count": sum(face_table),
            "face_vector": list(face_vector),
            "minimal_nonface_size_histogram": {
                str(key): value
                for key, value in sorted(Counter(value.bit_count() for value in minimal_nonfaces).items())
            },
            "ridge_degree_histogram": {
                str(key): value for key, value in sorted(ridge_histogram.items())
            },
            "dual_component_count": dual_components,
        },
        blocker,
    )


def _verify_selection_archive(
    blocker: Sequence[int], target_fives: Sequence[int]
) -> dict[str, object]:
    file_digest = _file_sha256(SELECTION_ARCHIVE)
    if file_digest != EXPECTED_SELECTION_FILE_SHA256:
        raise AssertionError("record-43 selection archive file changed")
    archive = json.loads(SELECTION_ARCHIVE.read_text(encoding="utf-8"))
    critical_keys = (
        "kind",
        "claim",
        "source",
        "incidence",
        "global_selection_audit",
        "gale_chirotope_signs",
        "gale_audit",
    )
    critical_digest = _canonical_sha256({key: archive[key] for key in critical_keys})
    if critical_digest != EXPECTED_SELECTION_CERTIFICATE_SHA256:
        raise AssertionError("record-43 selection critical digest changed")
    if archive.get("certificate_sha256") != critical_digest:
        raise AssertionError("record-43 selection internal digest is invalid")
    incidence = archive["incidence"]
    if tuple(map(tuple, incidence["missing_edges"])) != EDGES:
        raise AssertionError("record-43 archived edges differ")
    if tuple(sorted(_mask(face) for face in incidence["missing_fives"])) != tuple(target_fives):
        raise AssertionError("record-43 archived A5s differ")
    if tuple(incidence["facet_complement_masks"]) != tuple(blocker):
        raise AssertionError("record-43 archived facets differ")
    if incidence["facet_complement_sha256"] != EXPECTED_FACET_COMPLEMENT_SHA256:
        raise AssertionError("record-43 archived facet digest differs")
    if incidence["packed_face_table_sha256"] != EXPECTED_PACKED_FACE_TABLE_SHA256:
        raise AssertionError("record-43 archived face-table digest differs")

    global_audit = archive["global_selection_audit"]
    expected_nonfaces = {
        "2": list(sorted(EDGE_MASKS)),
        "5": list(sorted(target_fives)),
    }
    if global_audit["minimal_nonfaces_by_size"] != expected_nonfaces:
        raise AssertionError("record-43 archived minimal nonfaces differ")
    if (
        global_audit["graph6"] != GRAPH6
        or global_audit["facet_count"] != 90
        or global_audit["face_vector"] != [14, 85, 300, 690, 1088, 1190, 880, 405, 90]
        or global_audit["ridge_degree_histogram"] != {"2": 405}
        or global_audit["dual_component_count"] != 1
        or not global_audit["passes_selection_model"]
    ):
        raise AssertionError("record-43 archived global incidence projection differs")
    return {
        "path": str(SELECTION_ARCHIVE.relative_to(ROOT)),
        "file_sha256": file_digest,
        "certificate_sha256": critical_digest,
        "exact_incidence_equality": True,
    }


def _verify_coxeter_projection() -> dict[str, object]:
    file_digest = _file_sha256(COXETER_ARCHIVE)
    if file_digest != EXPECTED_COXETER_FILE_SHA256:
        raise AssertionError("selected-nerve Coxeter archive file changed")
    archive = json.loads(COXETER_ARCHIVE.read_text(encoding="utf-8"))
    payload = dict(archive)
    internal = payload.pop("certificate_sha256", None)
    certificate_digest = _canonical_sha256(payload)
    if certificate_digest != EXPECTED_COXETER_CERTIFICATE_SHA256 or internal != certificate_digest:
        raise AssertionError("selected-nerve Coxeter certificate digest changed")
    scope_reaudit = archive.get("scope_reaudit")
    if (
        archive.get("kind") != EXPECTED_COXETER_KIND
        or not isinstance(scope_reaudit, dict)
        or scope_reaudit.get("status") != EXPECTED_SCOPE_REAUDIT_STATUS
        or 43 not in scope_reaudit.get("reaudited_indices", [])
        or scope_reaudit.get("normalizing_support_symmetry_orders", {}).get("43")
        != EXPECTED_SUPPORT_SYMMETRY_ORDERS
        or "principal submatrix" not in scope_reaudit.get(
            "principal_block_lemma", ""
        )
        or "inertia" not in scope_reaudit.get("principal_block_lemma", "")
    ):
        raise AssertionError("selected-nerve principal-block scope re-audit changed")
    matches = [
        row for row in archive["six_nerve_census"] if row.get("record_index") == 43
    ]
    if len(matches) != 1:
        raise AssertionError("selected-nerve Coxeter archive lacks unique record 43")
    projection = matches[0]
    projection_digest = _canonical_sha256(projection)
    if projection_digest != EXPECTED_COXETER_RECORD43_PROJECTION_SHA256:
        raise AssertionError("record-43 Coxeter projection digest changed")
    expected_right_supports = {
        (1, 2, 8, 9, 10),
        (2, 8, 9, 10, 11),
        (8, 9, 10, 11, 12),
    }
    if {
        tuple(gate["right_support"]) for gate in projection["support_gates"]
    } != expected_right_supports:
        raise AssertionError("record-43 Coxeter support gates changed")
    if any(tuple(gate["left_support"]) != FIVES[0] for gate in projection["support_gates"]):
        raise AssertionError("record-43 isolated Lanner support changed")
    if (
        projection["rank_nine_vertex_count"] != 90
        or not projection["selected_nerve_only"]
        or not projection["unrestricted_compact_coxeter_no_go"]
        or projection.get("normalizing_support_symmetry_orders")
        != EXPECTED_SUPPORT_SYMMETRY_ORDERS
        or projection["two_edge_passing_count"] != 0
        or projection["monotonicity_residual_two_edge_trial_count"] != 0
        or projection["stronger_single_edge_survivors"] != {"4": 0, "5": 0}
        or any(witness["inertia"] != [8, 2, 0] for witness in projection["exact_index_two_witnesses"])
    ):
        raise AssertionError("record-43 compact-Coxeter theorem projection changed")

    source_digest = _file_sha256(COXETER_RECORD43_SOURCE)
    if source_digest != EXPECTED_COXETER_RECORD43_SOURCE_SHA256:
        raise AssertionError("record-43 selected-nerve source file changed")
    source_rows = [
        row
        for row in archive["sources"]["code"]
        if row["path"] == "src/coxeter_polytopes/combinatorics/h9_14_record43.py"
    ]
    if source_rows != [
        {
            "path": "src/coxeter_polytopes/combinatorics/h9_14_record43.py",
            "sha256": source_digest,
        }
    ]:
        raise AssertionError("record-43 Coxeter source provenance changed")
    return {
        "path": str(COXETER_ARCHIVE.relative_to(ROOT)),
        "file_sha256": file_digest,
        "certificate_sha256": certificate_digest,
        "record43_projection_sha256": projection_digest,
        "record43_source_sha256": source_digest,
        "aggregate_kind": archive["kind"],
        "scope_reaudit_status": scope_reaudit["status"],
        "principal_block_lemma": scope_reaudit["principal_block_lemma"],
        "normalizing_support_symmetry_orders": projection[
            "normalizing_support_symmetry_orders"
        ],
        "rank_nine_vertex_count": projection["rank_nine_vertex_count"],
        "unrestricted_compact_coxeter_no_go": True,
        "finite_volume_noncompact_checked": False,
    }


# Five order-five Lanner isomorphism types in lexicographic edge order.
_A5_LANNER_CODES = (
    (2, 2, 2, 3, 2, 2, 3, 5, 2, 3),
    (2, 2, 2, 3, 2, 3, 3, 5, 2, 2),
    (2, 2, 2, 4, 2, 3, 3, 5, 2, 2),
    (2, 2, 2, 5, 2, 3, 3, 5, 2, 2),
    (2, 2, 3, 3, 3, 2, 3, 4, 2, 2),
)


def _pair(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left < right else (right, left)


def _all_a5_states(support: Sequence[int]) -> tuple[dict[tuple[int, int], int], ...]:
    ordered = tuple(sorted(support))
    abstract_edges = tuple(combinations(range(5), 2))
    actual_edges = tuple(combinations(ordered, 2))
    states = set()
    for code in _A5_LANNER_CODES:
        abstract = dict(zip(abstract_edges, code, strict=True))
        for image_order in permutations(range(5)):
            inverse = {image: source for source, image in enumerate(image_order)}
            states.add(
                tuple(
                    abstract[_pair(inverse[left], inverse[right])]
                    for left, right in abstract_edges
                )
            )
    if len(states) != 420:
        raise AssertionError("independent A5 placement count changed")
    return tuple(
        dict(zip(actual_edges, code, strict=True)) for code in sorted(states)
    )


def _distinguished_a5_states(
    support: tuple[int, ...],
) -> tuple[tuple[dict[tuple[int, int], int], ...], dict[int, int]]:
    distinguished = support[0]
    interchangeable = support[1:]
    edges = tuple(combinations(sorted(support), 2))
    orbit_counts: Counter[tuple[int, ...]] = Counter()
    for state in _all_a5_states(support):
        codes = []
        for image_order in permutations(interchangeable):
            action = {
                distinguished: distinguished,
                **dict(zip(interchangeable, image_order, strict=True)),
            }
            inverse = {image: source for source, image in action.items()}
            codes.append(
                tuple(state[_pair(inverse[a], inverse[b])] for a, b in edges)
            )
        orbit_counts[min(codes)] += 1
    histogram = dict(sorted(Counter(orbit_counts.values()).items()))
    if len(orbit_counts) != 20 or histogram != {12: 5, 24: 15}:
        raise AssertionError("distinguished A5 orbit partition changed")
    return (
        tuple(
            dict(zip(edges, code, strict=True)) for code in sorted(orbit_counts)
        ),
        histogram,
    )


def _finite_component(
    component: frozenset[int], labels: Mapping[tuple[int, int], int]
) -> bool:
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
    degrees = Counter(map(len, adjacency.values()))
    if max(degrees) <= 2:
        endpoint = min(vertex for vertex in component if len(adjacency[vertex]) == 1)
        sequence = []
        previous = None
        current = endpoint
        while True:
            following = [row for row in adjacency[current] if row[0] != previous]
            if not following:
                break
            if len(following) != 1:
                return False
            next_vertex, label = following[0]
            sequence.append(label)
            previous, current = current, next_vertex
        path = min(tuple(sequence), tuple(reversed(sequence)))
        order = len(component)
        return (
            (order == 2 and path[0] in (3, 4, 5))
            or all(label == 3 for label in path)
            or path == (3,) * (order - 2) + (4,)
            or (order == 3 and path == (3, 5))
            or (order == 4 and path in ((3, 3, 5), (3, 4, 3)))
        )
    if degrees[3] != 1 or any(degree > 3 for degree in degrees):
        return False
    if any(
        labels[_pair(a, b)] != 3
        for a, b in combinations(component, 2)
        if labels[_pair(a, b)] > 2
    ):
        return False
    branch = next(vertex for vertex in component if len(adjacency[vertex]) == 3)
    arms = []
    for neighbor, _label in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(v for v, _ in adjacency[current] if v != previous)
            previous, current = current, following
            length += 1
        arms.append(length)
    arm_tuple = tuple(sorted(arms))
    return arm_tuple == (1, 1, len(component) - 3) or arm_tuple in {
        (1, 2, 2), (1, 2, 3), (1, 2, 4)
    }


def _finite_labeling(
    facets: Iterable[int], labels: Mapping[tuple[int, int], int]
) -> bool:
    vertices = tuple(sorted(set(facets)))
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
                for other in tuple(unseen)
                if labels[_pair(vertex, other)] > 2
            }
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        components.append(frozenset(component))
    return all(_finite_component(component, labels) for component in components)


def _block_graph(
    left: Sequence[int],
    right: Sequence[int],
    left_state: Mapping[tuple[int, int], int],
    right_state: Mapping[tuple[int, int], int],
    edge: tuple[int, int],
) -> nx.Graph:
    union = tuple((*left, *right))
    labels = {_pair(a, b): 2 for a, b in combinations(union, 2)}
    labels.update(left_state)
    labels.update(right_state)
    labels[edge] = 3
    graph = nx.Graph()
    graph.add_nodes_from(union)
    for a, b in combinations(union, 2):
        graph.add_edge(a, b, label=labels[_pair(a, b)])
    return graph


def _exact_sign(value: sp.Expr) -> str:
    bounds = {
        sp.sqrt(2): (sp.Rational(7, 5), sp.Rational(3, 2)),
        sp.sqrt(5): (sp.Rational(11, 5), sp.Rational(9, 4)),
        sp.sqrt(10): (sp.Rational(79, 25), sp.Rational(19, 6)),
    }
    remainder = sp.expand(value)
    lower = upper = sp.Integer(0)
    for radical, (lo, hi) in bounds.items():
        coefficient = remainder.coeff(radical)
        remainder = sp.expand(remainder - coefficient * radical)
        if coefficient >= 0:
            lower += coefficient * lo
            upper += coefficient * hi
        else:
            lower += coefficient * hi
            upper += coefficient * lo
    lower += remainder
    upper += remainder
    if lower > 0:
        return "+"
    if upper < 0:
        return "-"
    raise AssertionError("exact radical enclosure did not decide a sign")


def _exact_witness(
    left: tuple[int, ...], right: tuple[int, ...], kind: int
) -> dict[str, object]:
    canonical_left = tuple(
        dict(zip(combinations(sorted(left), 2), code, strict=True))
        for code in _A5_LANNER_CODES
    )
    canonical_right = tuple(
        dict(zip(combinations(sorted(right), 2), code, strict=True))
        for code in _A5_LANNER_CODES
    )
    if kind == 0:
        order = tuple((*left, *right))
        edge = _pair(left[0], right[0])
        left_state, right_state = canonical_left[1], canonical_right[1]
        expected = (
            2, 4, 8, 6 - 2 * sp.sqrt(5), 4 - 2 * sp.sqrt(5),
            (9 - 5 * sp.sqrt(5)) / 2, 9 - 5 * sp.sqrt(5),
            18 - 10 * sp.sqrt(5), 26 - 12 * sp.sqrt(5),
            (25 - 11 * sp.sqrt(5)) / 2,
        )
        name = "A5_A5_type_1_1"
    else:
        order = tuple(reversed(tuple((*left, *right))))
        edge = _pair(left[4], right[0])
        left_state, right_state = canonical_left[4], canonical_right[1]
        expected = (
            2, 4, 5 - sp.sqrt(5), (7 - 3 * sp.sqrt(5)) / 2,
            4 - 2 * sp.sqrt(5), (9 - 5 * sp.sqrt(5)) / 2,
            9 - 5 * sp.sqrt(5), 9 - 5 * sp.sqrt(5), 1 - sp.sqrt(5),
            (-23 - 16 * sp.sqrt(2) + 11 * sp.sqrt(5) + 8 * sp.sqrt(10)) / 2,
        )
        name = "A5_A5_type_4_1"
    labels = {_pair(a, b): 2 for a, b in combinations((*left, *right), 2)}
    labels.update(left_state)
    labels.update(right_state)
    labels[edge] = 3
    gram_value = {2: 0, 3: -1, 4: -sp.sqrt(2), 5: -(1 + sp.sqrt(5)) / 2}
    matrix = sp.Matrix([
        [2 if a == b else gram_value[labels[_pair(a, b)]] for b in order]
        for a in order
    ])
    observed = tuple(sp.factor(matrix[:size, :size].det()) for size in range(1, 11))
    if any(sp.simplify(a - b) != 0 for a, b in zip(observed, expected, strict=True)):
        raise AssertionError("independent exact exceptional-block determinant changed")
    signs = "".join(_exact_sign(value) for value in expected)
    negative = sum(a != b for a, b in zip(("+" + signs)[:-1], ("+" + signs)[1:], strict=True))
    if negative != 2:
        raise AssertionError("independent exceptional block lost index two")
    return {"name": name, "signs": signs, "inertia": [8, 2, 0]}


def _audit_second_label_gate(second_blocker: Sequence[int]) -> dict[str, object]:
    left, right = SECOND_FIVES[:2]
    left_states, left_orbits = _distinguished_a5_states(left)
    right_states, right_orbits = _distinguished_a5_states(right)
    clutter = {frozenset(face) for face in (*EDGES, *SECOND_FIVES)}
    automorphisms = 0
    for left_order in permutations(left[1:]):
        for right_order in permutations(right[1:]):
            action = dict(zip(left[1:], left_order, strict=True))
            action.update(dict(zip(right[1:], right_order, strict=True)))
            image = {
                frozenset(action.get(vertex, vertex) for vertex in face)
                for face in clutter
            }
            if image != clutter:
                raise AssertionError("claimed S4 x S4 action is not a nerve automorphism")
            automorphisms += 1
    vertices = tuple(
        tuple(vertex for vertex in range(VERTEX_COUNT) if not complement >> vertex & 1)
        for complement in second_blocker
    )
    union = tuple(sorted((*left, *right)))
    restrictions = tuple(sorted({
        tuple(facet for facet in vertex if facet in union) for vertex in vertices
    }, key=lambda row: (len(row), row)))
    if Counter(map(len, restrictions)) != {6: 16, 7: 24, 8: 25}:
        raise AssertionError("second-nerve actual restriction census changed")
    cross_pairs = tuple(_pair(a, b) for a in left for b in right)
    by_edge = {
        edge: tuple(r for r in restrictions if set(edge) <= set(r))
        for edge in cross_pairs
    }
    passing = []
    template = {_pair(a, b): 2 for a, b in combinations(union, 2)}
    for left_index, left_state in enumerate(left_states):
        for right_index, right_state in enumerate(right_states):
            sparse = dict(template)
            sparse.update(left_state)
            sparse.update(right_state)
            for edge in cross_pairs:
                labels = dict(sparse)
                labels[edge] = 3
                if all(_finite_labeling(r, labels) for r in by_edge[edge]):
                    passing.append((left_index, right_index, edge))
    if len(passing) != 55:
        raise AssertionError("second-nerve label-three survivor count changed")
    stronger = {}
    for label in (4, 5):
        stronger[label] = sum(
            all(
                _finite_labeling(
                    restriction,
                    {
                        **template,
                        **left_states[left_index],
                        **right_states[right_index],
                        edge: label,
                    },
                )
                for restriction in by_edge[edge]
            )
            for left_index, right_index, edge in passing
        )
    by_state: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for left_index, right_index, edge in passing:
        by_state.setdefault((left_index, right_index), []).append(edge)
    residual_pairs = sum(
        comb(len(set(edges)), 2) for edges in by_state.values()
    )
    if stronger != {4: 0, 5: 0} or residual_pairs != 0:
        raise AssertionError("second-nerve monotonicity gate changed")

    canonical_left = tuple(
        dict(zip(combinations(sorted(left), 2), code, strict=True))
        for code in _A5_LANNER_CODES
    )
    canonical_right = tuple(
        dict(zip(combinations(sorted(right), 2), code, strict=True))
        for code in _A5_LANNER_CODES
    )
    references = (
        _block_graph(left, right, canonical_left[1], canonical_right[1], _pair(left[0], right[0])),
        _block_graph(left, right, canonical_left[4], canonical_right[1], _pair(left[4], right[0])),
    )
    classes = Counter()
    names = ("A5_A5_type_1_1", "A5_A5_type_4_1")
    edge_match = nx.algorithms.isomorphism.categorical_edge_match("label", None)
    for left_index, right_index, edge in passing:
        graph = _block_graph(
            left, right, left_states[left_index], right_states[right_index], edge
        )
        matches = [
            index
            for index, reference in enumerate(references)
            if nx.is_isomorphic(graph, reference, edge_match=edge_match)
        ]
        if len(matches) != 1:
            raise AssertionError("independent exceptional-block classification changed")
        classes[names[matches[0]]] += 1
    witnesses = (_exact_witness(left, right, 0), _exact_witness(left, right, 1))
    if classes != {"A5_A5_type_1_1": 25, "A5_A5_type_4_1": 30}:
        raise AssertionError("independent exceptional-class histogram changed")
    return {
        "verified_nerve_automorphism_count": automorphisms,
        "full_placement_count_per_support": 420,
        "representative_count_per_support": 20,
        "left_orbit_size_histogram": {str(k): v for k, v in left_orbits.items()},
        "right_orbit_size_histogram": {str(k): v for k, v in right_orbits.items()},
        "label_three_test_count": 20 * 20 * 25,
        "label_three_passing_orbit_count": len(passing),
        "stronger_single_edge_survivors": {str(k): v for k, v in stronger.items()},
        "raw_two_edge_orbit_trial_count": 20 * 20 * comb(25, 2),
        "monotonicity_residual_two_edge_trial_count": residual_pairs,
        "two_edge_passing_count": 0,
        "exceptional_class_histogram": dict(sorted(classes.items())),
        "exact_index_two_witnesses": list(witnesses),
        "unrestricted_compact_coxeter_no_go": True,
    }


def replay() -> dict[str, object]:
    """Run and return the complete deterministic hostile-replay summary."""

    if tuple(len(_independent_masks(size)) for size in range(2, 6)) != (
        85,
        300,
        690,
        1092,
    ):
        raise AssertionError("record-43 independent-set counts changed")
    schedules = _enumerate_schedules()
    strict_survivors = tuple(
        int(row["missing_four_count"])
        for row in schedules
        if int(row["summed_vertex_link_g4"]) >= VERTEX_COUNT
    )
    if strict_survivors != (0, 1):
        raise AssertionError("record-43 strict schedule reduction changed")
    a1 = _audit_a1()
    a0, targets = _a0_census()
    first_audit, first_blocker = _audit_target(
        targets[0],
        EXPECTED_FACET_COMPLEMENT_SHA256,
        EXPECTED_PACKED_FACE_TABLE_SHA256,
    )
    second_audit, second_blocker = _audit_target(
        targets[1],
        EXPECTED_SECOND_FACET_COMPLEMENT_SHA256,
        EXPECTED_SECOND_PACKED_FACE_TABLE_SHA256,
    )
    selection = _verify_selection_archive(first_blocker, targets[0])
    coxeter = _verify_coxeter_projection()
    second_labels = _audit_second_label_gate(second_blocker)

    summary: dict[str, object] = {
        "kind": "h9_14_record43_all_completions_no_import_hostile_replay",
        "graph6": GRAPH6,
        "missing_edges": [list(edge) for edge in EDGES],
        "independent_set_counts": {
            str(size): len(_independent_masks(size)) for size in range(2, 6)
        },
        "raw_schedules": list(schedules),
        "strict_schedule_survivors": list(strict_survivors),
        "strict_schedule_exclusions": [2, 3],
        "a1_census": a1,
        "a0_census": a0,
        "first_target_audit": first_audit,
        "second_target_audit": second_audit,
        "selection_archive": selection,
        "first_target_compact_coxeter_archive": coxeter,
        "second_target_independent_label_gate": second_labels,
        "theorem": {
            "whole_record43_graph_shadow_compact_coxeter_impossible": True,
            "incidence_completion_orbit_count": 2,
            "both_incidence_targets_compact_label_impossible": True,
        },
        "scope": {
            "compact_case_excluded": True,
            "finite_volume_noncompact_checked": False,
            "maximal_parabolic_cusp_audit_run": False,
            "elliptic_edge_endpoint_audit_run": False,
        },
    }
    digest = _canonical_sha256(summary)
    if EXPECTED_SUMMARY_SHA256 != "TO_BE_FROZEN" and digest != EXPECTED_SUMMARY_SHA256:
        raise AssertionError("record-43 hostile replay summary digest changed")
    summary["summary_sha256"] = digest
    return summary


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2, sort_keys=True))
