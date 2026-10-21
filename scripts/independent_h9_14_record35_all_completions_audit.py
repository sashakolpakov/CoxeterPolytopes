#!/usr/bin/env python3
"""No-project-import replay of the record-35 whole-shadow compact no-go.

The verifier uses only the Python standard library.  It reconstructs every
Dehn--Sommerville branch, permits every independent active subset allowed by
the exact incidence caps in an A5 row, repeats the complete row/column blocker
census, regenerates the
classical order-four/five Lanner placement tables, and applies every actual
rank-nine vertex test.  The residual ``a=1`` exceptional blocks are checked
in exact ``Q(sqrt(2),sqrt(5))`` arithmetic.

No ideal-face or finite-volume noncompact claim is made.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations
import json
from math import comb
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
VERTEX_COUNT = 14
DIMENSION = 9
FULL_MASK = (1 << VERTEX_COUNT) - 1
GRAPH6 = "M?????????O?O?w??"
EDGES = ((0, 11), (1, 12), (0, 13), (1, 13), (2, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in EDGES)
ACTIVE = (0, 1, 2, 11, 12, 13)
ISOLATES = tuple(range(3, 11))
ROW_COUNT = 5

A1_Q = (2, 3, 4, 5)
A1_R = (3, 4, 5, 11, 12)
A1_B = (6, 7, 8, 9, 10)
A1_LEFT = tuple(sorted(set(A1_Q) | set(A1_R)))
A1_UNION = tuple(sorted(set(A1_LEFT) | set(A1_B)))
A1_CROSS = tuple(
    tuple(sorted((left, right))) for left in A1_LEFT for right in A1_B
)

EXPECTED_TARGET_FAMILY_SHA256 = (
    "a3678a721ea0e69b1fa82381eaa4625e996c826c06f52c168d94f20a962671d1"
)
EXPECTED_A0_LABEL_TRANSCRIPT_SHA256 = (
    "edc677f407f3a011579b42fff7e0990473d6bd38191cc3d9e852705ae9f6d14a"
)
EXPECTED_A1_LABEL_TRANSCRIPT_SHA256 = (
    "45b32fb880ba66a77202682b5db5532a0ba055017274ceb44e99a3e35cc2804e"
)
EXPECTED_SUMMARY_SHA256 = (
    "b01b632d01569c847194a8fb592eca4e91b0cc9865e00d296d029147b3e25332"
)

EXPECTED_ACTIVE_MASKS = (
    0,
    1,
    2,
    4,
    2_048,
    4_096,
    3,
    5,
    4_097,
    6,
    2_050,
    2_052,
    4_100,
    6_144,
    7,
    4_101,
    2_054,
    6_148,
)
EXPECTED_ACTIVE_MASK_SIZE_HISTOGRAM = {0: 1, 1: 5, 2: 8, 3: 4}
EXPECTED_ACTIVE_TRIPLES = (
    (0, 1, 2),
    (0, 2, 12),
    (1, 2, 11),
    (2, 11, 12),
)


def _mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def _vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1)


def _pair(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left < right else (right, left)


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode()
    ).hexdigest()


def _is_independent(value: int) -> bool:
    return all(value & edge != edge for edge in EDGE_MASKS)


@lru_cache(maxsize=None)
def _independent_masks(order: int) -> tuple[int, ...]:
    return tuple(
        value
        for chosen in combinations(range(VERTEX_COUNT), order)
        for value in (_mask(chosen),)
        if _is_independent(value)
    )


def _macaulay_raise(value: int, degree: int) -> int:
    remaining = value
    next_top = None
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
        raise AssertionError("hostile Macaulay expansion failed")
    return raised


def _h_prefix(a: int) -> tuple[int, ...]:
    partial_f = (1, 14, 86, 309, 726 - a)
    return tuple(
        sum(
            (-1) ** (index - order)
            * comb(DIMENSION - order, DIMENSION - index)
            * partial_f[order]
            for order in range(index + 1)
        )
        for index in range(5)
    )


def _summed_link_g4(f_vector: Sequence[int]) -> int:
    values = {0: VERTEX_COUNT}
    for order in range(1, 5):
        values[order] = (order + 1) * int(f_vector[order])

    def h(index: int) -> int:
        return sum(
            (-1) ** (index - order)
            * comb(8 - order, 8 - index)
            * values[order]
            for order in range(index + 1)
        )

    return h(4) - h(3)


def _schedules() -> tuple[tuple[int, int, int, int], ...]:
    answers = []
    independent_fives = len(_independent_masks(5))
    for a in range(6):
        prefix = _h_prefix(a)
        if prefix[4] - prefix[3] < 1:
            break
        h_vector = prefix + tuple(reversed(prefix))
        f_vector = tuple(
            sum(
                comb(DIMENSION - order, DIMENSION - index) * h_vector[order]
                for order in range(index + 1)
            )
            for index in range(1, DIMENSION + 1)
        )
        answers.append(
            (
                a,
                independent_fives - f_vector[4],
                f_vector[-1],
                _summed_link_g4(f_vector),
            )
        )
    expected = (
        (0, 5, 102, 30),
        (1, 11, 100, 24),
        (2, 17, 98, 18),
        (3, 23, 96, 12),
        (4, 29, 94, 6),
    )
    if tuple(answers) != expected:
        raise AssertionError("record-35 hostile schedules changed")
    return tuple(answers)


def _strict_interval(
    previous: int, baseline: int, candidates: int, degree: int
) -> tuple[int, int] | None:
    lower = max(0, baseline - _macaulay_raise(previous, degree))
    upper = min(candidates, baseline - 1)
    return (lower, upper) if lower <= upper else None


@lru_cache(maxsize=1)
def _local_constraints():
    vertex_constraints = []
    for vertex in range(VERTEX_COUNT):
        link = tuple(
            other
            for other in range(VERTEX_COUNT)
            if other != vertex and _is_independent(_mask((vertex, other)))
        )
        n = len(link)
        f1 = sum(_is_independent(_mask(edge)) for edge in combinations(link, 2))
        f2 = sum(_is_independent(_mask(triple)) for triple in combinations(link, 3))
        g1 = n - 9
        g2 = f1 - 8 * n + 36
        beta = f2 - 7 * f1 + 28 * n - 84
        if not (g1 > 0 and 1 <= g2 <= _macaulay_raise(g1, 1)):
            raise AssertionError("hostile vertex-link baseline failed")
        vertex_constraints.append(
            (vertex, _strict_interval(g2, beta, f2, 2))
        )
    edge_constraints = []
    for edge in combinations(range(VERTEX_COUNT), 2):
        if not _is_independent(_mask(edge)):
            continue
        link = tuple(
            other
            for other in range(VERTEX_COUNT)
            if other not in edge and _is_independent(_mask((*edge, other)))
        )
        n = len(link)
        f1 = sum(_is_independent(_mask(pair)) for pair in combinations(link, 2))
        g1 = n - 8
        beta = f1 - 7 * n + 28
        edge_constraints.append((edge, _strict_interval(g1, beta, f1, 1)))
    triple_constraints = []
    for triple in combinations(range(VERTEX_COUNT), 3):
        if not _is_independent(_mask(triple)):
            continue
        candidate = sum(
            _is_independent(_mask((*triple, other)))
            for other in range(VERTEX_COUNT)
            if other not in triple
        )
        triple_constraints.append((triple, candidate - 8))
    return tuple(vertex_constraints), tuple(edge_constraints), tuple(triple_constraints)


def _a4_interval_feasible(value: int, constraints) -> bool:
    vertex_constraints, edge_constraints, triple_constraints = constraints
    for vertex, interval in vertex_constraints:
        incidence = int(bool(value & (1 << vertex)))
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for edge, interval in edge_constraints:
        incidence = int(value & _mask(edge) == _mask(edge))
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    return all(
        int(value & _mask(triple) == _mask(triple)) <= maximum
        for triple, maximum in triple_constraints
    )


def _family_interval_feasible(family: Sequence[int], constraints) -> bool:
    vertex_constraints, edge_constraints, triple_constraints = constraints
    for vertex, interval in vertex_constraints:
        incidence = sum(bool(face & (1 << vertex)) for face in family)
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for edge, interval in edge_constraints:
        subset = _mask(edge)
        incidence = sum(face & subset == subset for face in family)
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    return all(
        sum(face & _mask(triple) == _mask(triple) for face in family) <= maximum
        for triple, maximum in triple_constraints
    )


def _base_link_g3_g4() -> tuple[tuple[int, ...], tuple[int, ...]]:
    independent = {order: _independent_masks(order) for order in range(2, 6)}
    g3, g4 = [], []
    for vertex in range(VERTEX_COUNT):
        counts = {
            order: sum(bool(value & (1 << vertex)) for value in values)
            for order, values in independent.items()
        }
        f0, f1, f2, f3 = (counts[order] for order in range(2, 6))
        g3.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    return tuple(g3), tuple(g4)


def _shadow_gate(family, shadow, base_g3, base_g4, total):
    additional = total - len(shadow)
    if additional < 0:
        return False, False, additional, (), ()
    lower, upper = [], []
    strict = True
    for vertex in range(VERTEX_COUNT):
        a4 = sum(bool(face & (1 << vertex)) for face in family)
        previous = base_g3[vertex] - a4
        forced = sum(bool(face & (1 << vertex)) for face in shadow)
        before = base_g4[vertex] + 6 * a4 - forced
        strict &= previous >= 1 and before >= 1
        lower.append(max(0, before - _macaulay_raise(previous, 3)) if previous >= 0 else 0)
        upper.append(before - 1)
    aggregate = strict and all(x >= 0 for x in upper) and sum(lower) <= 5 * additional <= sum(upper)
    return strict, aggregate, additional, tuple(lower), tuple(upper)


def _active_automorphisms() -> tuple[dict[int, int], ...]:
    edge_set = {frozenset(edge) for edge in EDGES}
    answers = []
    for image_order in permutations(ACTIVE):
        action = dict(zip(ACTIVE, image_order, strict=True))
        if {frozenset((action[a], action[b])) for a, b in EDGES} == edge_set:
            answers.append(action)
    if len(answers) != 2:
        raise AssertionError("hostile active graph automorphism count changed")
    return tuple(answers)


def _orbit_key(family: Sequence[int]):
    active_mask = _mask(ACTIVE)
    keys = []
    for action in _active_automorphisms():
        transformed = tuple(
            _mask(action.get(vertex, vertex) for vertex in _vertices(face))
            for face in family
        )
        for row_order in permutations(range(len(family))):
            active_rows = tuple(transformed[index] & active_mask for index in row_order)
            columns = tuple(
                sorted(
                    sum(
                        bool(transformed[row_order[row]] & (1 << isolate)) << row
                        for row in range(len(family))
                    )
                    for isolate in ISOLATES
                )
            )
            keys.append((active_rows, columns))
    return min(keys)


def _minimal_graph_covers() -> tuple[int, ...]:
    answers = []
    for value in range(1, 1 << VERTEX_COUNT):
        if not all(value & edge for edge in EDGE_MASKS):
            continue
        if all(
            not all((value ^ (1 << vertex)) & edge for edge in EDGE_MASKS)
            for vertex in _vertices(value)
        ):
            answers.append(value)
    expected = tuple(
        _mask(face)
        for face in (
            (0, 1, 2),
            (0, 1, 13),
            (1, 11, 13),
            (0, 12, 13),
            (11, 12, 13),
        )
    )
    if tuple(answers) != expected:
        raise AssertionError("hostile minimal graph covers changed")
    return tuple(answers)


def _small_transversals(hyperedges: Sequence[int]) -> tuple[int, ...]:
    return tuple(
        candidate
        for order in range(1, 5)
        for chosen in combinations(range(VERTEX_COUNT), order)
        for candidate in (_mask(chosen),)
        if all(candidate & edge for edge in hyperedges)
    )


def _hit_word(face: int, covers: Sequence[int]) -> int:
    return sum(bool(face & cover) << index for index, cover in enumerate(covers))


def _addition_census(family, total_fives, scheduled_facets):
    fives = _independent_masks(5)
    base_g3, base_g4 = _base_link_g3_g4()
    shadow = tuple(
        five for five in fives if any(five & four == four for four in family)
    )
    strict, aggregate, additional, lower, upper = _shadow_gate(
        family, shadow, base_g3, base_g4, total_fives
    )
    if not (strict and aggregate):
        raise AssertionError("hostile addition census received bad shadow")
    candidates = tuple(
        five
        for five in fives
        if five not in shadow
        and not any(five & four == four for four in family)
        and all(not (five & (1 << vertex)) or upper[vertex] >= 1 for vertex in range(VERTEX_COUNT))
    )
    small = _small_transversals((*EDGE_MASKS, *family, *shadow))
    covers5 = tuple(
        candidate
        for chosen in combinations(range(VERTEX_COUNT), 5)
        for candidate in (_mask(chosen),)
        if all(candidate & edge for edge in EDGE_MASKS)
    )
    hits = {}
    fixed = (1 << len(covers5)) - 1
    for face in (*family, *shadow):
        hits.setdefault(face, _hit_word(face, covers5))
        fixed &= hits[face]
    lower_at_least = tuple(
        _mask(vertex for vertex in range(VERTEX_COUNT) if lower[vertex] >= threshold)
        for threshold in range(1, additional + 1)
    )
    upper_below = tuple(
        _mask(vertex for vertex in range(VERTEX_COUNT) if upper[vertex] < threshold)
        for threshold in range(1, additional + 1)
    )
    cap = no_small = 0
    histogram = Counter()
    targets = []
    for additions in combinations(candidates, additional):
        if additional == 2:
            first, second = additions
            at_least = (first | second, first & second)
        elif additional == 3:
            first, second, third = additions
            at_least = (
                first | second | third,
                (first & second) | (first & third) | (second & third),
                first & second & third,
            )
        else:
            raise AssertionError("unexpected hostile addition count")
        if any(
            required & present != required or forbidden & present
            for required, forbidden, present in zip(lower_at_least, upper_below, at_least, strict=True)
        ):
            continue
        cap += 1
        if any(all(transversal & face for face in additions) for transversal in small):
            continue
        no_small += 1
        word = fixed
        for face in additions:
            hits.setdefault(face, _hit_word(face, covers5))
            word &= hits[face]
        facets = word.bit_count()
        histogram[facets] += 1
        if facets == scheduled_facets:
            targets.append(tuple(additions))
    return {
        "representative_fours": [list(_vertices(face)) for face in family],
        "shadow": len(shadow),
        "additional": additional,
        "candidates": len(candidates),
        "raw": comb(len(candidates), additional),
        "cap": cap,
        "no_small": no_small,
        "facet_histogram": dict(sorted(histogram.items())),
        "targets": [[list(_vertices(face)) for face in additions] for additions in targets],
    }


def _positive_branches() -> dict[str, object]:
    constraints = _local_constraints()
    vertex_constraints, edge_constraints, _triple_constraints = constraints
    if any(
        interval is None or interval[0] != 0
        for _support, interval in (*vertex_constraints, *edge_constraints)
    ):
        raise AssertionError("hostile A4 lower-bound audit changed")
    fours = tuple(
        face for face in _independent_masks(4) if _a4_interval_feasible(face, constraints)
    )
    fives = _independent_masks(5)
    base_g3, base_g4 = _base_link_g3_g4()
    gates = []
    viable_a1 = []
    for four in fours:
        shadow = tuple(five for five in fives if five & four == four)
        strict, aggregate, additional, _lower, _upper = _shadow_gate(
            (four,), shadow, base_g3, base_g4, 11
        )
        gates.append((strict, aggregate, len(shadow), additional))
        if strict and aggregate:
            viable_a1.append(four)
    a1_orbits = defaultdict(list)
    for four in viable_a1:
        a1_orbits[_orbit_key((four,))].append(four)
    representatives = tuple(sorted(min(values) for values in a1_orbits.values()))
    a1_censuses = tuple(_addition_census((value,), 11, 100) for value in representatives)

    local_pairs = 0
    shadow_histogram = Counter()
    viable_a2 = []
    for family in combinations(fours, 2):
        if not _family_interval_feasible(family, constraints):
            continue
        local_pairs += 1
        shadow = tuple(
            five for five in fives if any(five & four == four for four in family)
        )
        strict, aggregate, _additional, _lower, _upper = _shadow_gate(
            family, shadow, base_g3, base_g4, 17
        )
        shadow_histogram[len(shadow)] += 1
        if strict and aggregate:
            viable_a2.append(family)
    a2_orbits = defaultdict(list)
    for family in viable_a2:
        a2_orbits[_orbit_key(family)].append(family)
    if len(a2_orbits) != 1:
        raise AssertionError("hostile a2 orbit count changed")
    a2_representative = min(next(iter(a2_orbits.values())))
    a2_census = _addition_census(a2_representative, 17, 98)
    expected_gates = {
        (False, False, 8, 3): 196,
        (False, False, 9, 2): 112,
        (False, False, 10, 1): 70,
        (True, True, 8, 3): 56,
        (True, True, 9, 2): 56,
    }
    if (
        len(fours) != 490
        or dict(Counter(gates)) != expected_gates
        or len(viable_a1) != 112
        or sorted(map(len, a1_orbits.values())) != [56, 56]
        or a1_censuses[0]["raw"] != 31_626
        or a1_censuses[0]["cap"] != 5_626
        or a1_censuses[0]["no_small"] != 1
        or len(a1_censuses[0]["targets"]) != 1
        or a1_censuses[1]["raw"] != 16_328_620
        or a1_censuses[1]["cap"] != 1_794_605
        or a1_censuses[1]["no_small"] != 840
        or a1_censuses[1]["targets"]
        or comb(len(fours), 2) != 119_805
        or local_pairs != 111_699
        or len(viable_a2) != 28
        or len(next(iter(a2_orbits.values()))) != 28
        or a2_census["raw"] != 1_540
        or a2_census["cap"] != 540
        or a2_census["no_small"] != 0
    ):
        raise AssertionError("hostile positive branch census changed")
    return {
        "a4_count": len(fours),
        "a1_gate_histogram": [
            {"gate": list(gate), "count": count}
            for gate, count in sorted(Counter(gates).items())
        ],
        "a1_orbit_sizes": sorted(map(len, a1_orbits.values())),
        "a1_censuses": list(a1_censuses),
        "a2_pair_count": comb(len(fours), 2),
        "a2_local_count": local_pairs,
        "a2_shadow_histogram": dict(sorted(shadow_histogram.items())),
        "a2_viable_count": len(viable_a2),
        "a2_orbit_sizes": sorted(map(len, a2_orbits.values())),
        "a2_census": a2_census,
    }


def _distribution_iterator(allowed, row_sums):
    counts = [0] * len(allowed)

    def visit(position, left, sums):
        if position == len(allowed):
            if left == 0 and not any(sums):
                yield tuple(counts)
            return
        column = allowed[position]
        maximum = left
        for row in range(ROW_COUNT):
            if column >> row & 1:
                maximum = min(maximum, sums[row])
        for multiplicity in range(maximum + 1):
            next_left = left - multiplicity
            next_sums = tuple(
                sums[row] - multiplicity * bool(column >> row & 1)
                for row in range(ROW_COUNT)
            )
            if all(0 <= value <= next_left for value in next_sums):
                counts[position] = multiplicity
                yield from visit(position + 1, next_left, next_sums)
        counts[position] = 0

    yield from visit(0, len(ISOLATES), row_sums)


@lru_cache(maxsize=None)
def _column_distributions(allowed, row_sums):
    return tuple(_distribution_iterator(allowed, row_sums))


def _faces(active_rows, allowed, counts):
    columns = tuple(
        column
        for column, multiplicity in zip(allowed, counts, strict=True)
        for _ in range(multiplicity)
    )
    rows = list(active_rows)
    for isolate, column in zip(ISOLATES, columns, strict=True):
        for row in range(ROW_COUNT):
            if column >> row & 1:
                rows[row] |= 1 << isolate
    return tuple(sorted(rows))


def _has_small_transversal(fives, graph_covers):
    for cover in graph_covers:
        missed = tuple(face for face in fives if not face & cover)
        if not missed:
            return True
        common = missed[0]
        for face in missed[1:]:
            common &= face
        if common:
            return True
    return False


def _minimal_transversals(hyperedges):
    answers = []
    for candidate in range(1, FULL_MASK + 1):
        if not all(candidate & edge for edge in hyperedges):
            continue
        if any(
            all((candidate ^ (1 << vertex)) & edge for edge in hyperedges)
            for vertex in _vertices(candidate)
        ):
            continue
        answers.append(candidate)
    return tuple(answers)


def _ridge_histogram(complements):
    facets = tuple(FULL_MASK ^ complement for complement in complements)
    ridges = Counter()
    for facet in facets:
        remaining = facet
        while remaining:
            bit = remaining & -remaining
            ridges[facet ^ bit] += 1
            remaining ^= bit
    return dict(sorted(Counter(ridges.values()).items()))


def _face_table(complements: Sequence[int]) -> tuple[bool, ...]:
    """Generate the complete downward-closed complex from rank-nine facets."""

    table = [False] * (1 << VERTEX_COUNT)
    table[0] = True
    for complement in sorted(set(complements)):
        if complement < 0 or complement & ~FULL_MASK or complement.bit_count() != 5:
            raise ValueError("hostile facet complement is not a five-set")
        facet = FULL_MASK ^ complement
        subset = facet
        while True:
            table[subset] = True
            if not subset:
                break
            subset = (subset - 1) & facet
    return tuple(table)


def _face_vector(face_table: Sequence[bool]) -> tuple[int, ...]:
    return tuple(
        sum(
            face_table[_mask(chosen)]
            for chosen in combinations(range(VERTEX_COUNT), size)
        )
        for size in range(1, DIMENSION + 1)
    )


def _minimal_nonfaces_from_table(face_table: Sequence[bool]) -> tuple[int, ...]:
    answers = []
    for candidate in range(1, 1 << VERTEX_COUNT):
        if face_table[candidate]:
            continue
        if all(
            face_table[candidate ^ (1 << vertex)]
            for vertex in _vertices(candidate)
        ):
            answers.append(candidate)
    return tuple(answers)


def _ridge_groups(complements: Sequence[int]) -> dict[int, tuple[int, ...]]:
    groups: dict[int, list[int]] = {}
    for complement in sorted(set(complements)):
        facet = FULL_MASK ^ complement
        remaining = facet
        while remaining:
            bit = remaining & -remaining
            groups.setdefault(complement | bit, []).append(complement)
            remaining ^= bit
    return {
        ridge_complement: tuple(sorted(group))
        for ridge_complement, group in groups.items()
    }


def _dual_component_count(complements: Sequence[int]) -> int:
    selected = tuple(sorted(set(complements)))
    if not selected:
        return 0
    adjacency = {complement: set() for complement in selected}
    for group in _ridge_groups(selected).values():
        if len(group) == 2:
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


def _first_link_euler_failure(
    face_table: Sequence[bool],
) -> tuple[int, int, int] | None:
    """Check the Euler characteristic of the complex and every proper link."""

    for face, is_face in enumerate(face_table):
        if not is_face or face.bit_count() > DIMENSION - 1:
            continue
        remaining = FULL_MASK ^ face
        observed = 0
        subset = remaining
        while subset:
            if face_table[face | subset]:
                observed += 1 if subset.bit_count() % 2 else -1
            subset = (subset - 1) & remaining
        link_dimension = DIMENSION - 1 - face.bit_count()
        expected = 1 + (-1) ** link_dimension
        if observed != expected:
            return face, observed, expected
    return None


def _gf2_rank(columns: Iterable[int]) -> int:
    pivots: dict[int, int] = {}
    for raw_column in columns:
        if raw_column < 0:
            raise ValueError("hostile GF(2) vector is negative")
        column = raw_column
        while column:
            pivot = column.bit_length() - 1
            previous = pivots.get(pivot)
            if previous is None:
                pivots[pivot] = column
                break
            column ^= previous
    return len(pivots)


def _link_face_masks(
    face_table: Sequence[bool], face: int
) -> tuple[tuple[int, ...], ...]:
    remaining = FULL_MASK ^ face
    grouped: list[list[int]] = []
    subset = remaining
    while subset:
        if face_table[face | subset]:
            dimension = subset.bit_count() - 1
            while len(grouped) <= dimension:
                grouped.append([])
            grouped[dimension].append(subset)
        subset = (subset - 1) & remaining
    return tuple(tuple(sorted(group)) for group in grouped)


def _reduced_betti_mod2(
    faces_by_dimension: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    if not faces_by_dimension:
        return ()
    normalized = tuple(tuple(sorted(set(group))) for group in faces_by_dimension)
    if any(
        face <= 0 or face.bit_count() != dimension + 1
        for dimension, group in enumerate(normalized)
        for face in group
    ):
        raise ValueError("hostile link faces have the wrong dimensions")
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
                    raise ValueError("hostile link face collection is not downward closed")
                column |= 1 << lower_index[boundary]
                bits ^= bit
            columns.append(column)
        boundary_ranks[dimension] = _gf2_rank(columns)
    ordinary = tuple(
        len(normalized[dimension])
        - boundary_ranks[dimension]
        - boundary_ranks[dimension + 1]
        for dimension in range(len(normalized))
    )
    if ordinary[0] <= 0:
        raise ValueError("hostile nonempty link has zero-dimensional H_0")
    return (ordinary[0] - 1, *ordinary[1:])


def _first_mod2_homology_link_failure(
    face_table: Sequence[bool],
) -> tuple[int, int, tuple[int, ...]] | None:
    """Check the complex and every face link as an exact GF(2) homology sphere."""

    for face, is_face in enumerate(face_table):
        if not is_face or face.bit_count() > DIMENSION:
            continue
        expected_dimension = DIMENSION - 1 - face.bit_count()
        link_faces = _link_face_masks(face_table, face)
        if expected_dimension == -1:
            if link_faces:
                return face, expected_dimension, ()
            continue
        observed = _reduced_betti_mod2(link_faces)
        expected_length = expected_dimension + 1
        padded = observed + (0,) * max(0, expected_length - len(observed))
        expected = tuple(
            int(degree == expected_dimension)
            for degree in range(expected_length)
        )
        if padded[:expected_length] != expected or any(padded[expected_length:]):
            return face, expected_dimension, observed
    return None


def _link_g_entries(
    face_table: Sequence[bool], face: int
) -> tuple[int, tuple[int, ...]]:
    face_size = face.bit_count()
    polytope_dimension = DIMENSION - face_size
    remaining = _vertices(FULL_MASK ^ face)
    link_f = {-1: 1}
    for dimension in range(polytope_dimension):
        order = dimension + 1
        link_f[dimension] = sum(
            face_table[face | _mask(extension)]
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
    return link_f[0], tuple(h[index] - h[index - 1] for index in range(1, len(h)))


def _first_strict_local_g_failure(
    face_table: Sequence[bool],
) -> tuple[int, int, int, int] | None:
    strict_orders = {1: (2, 3, 4), 2: (2, 3), 3: (2,)}
    for face_size, required_orders in strict_orders.items():
        polytope_dimension = DIMENSION - face_size
        for face in _independent_masks(face_size):
            if not face_table[face]:
                return face, 0, -1, 0
            link_vertices, g_entries = _link_g_entries(face_table, face)
            if link_vertices < polytope_dimension + 2:
                return face, 0, link_vertices, polytope_dimension + 2
            for order, observed in enumerate(g_entries, start=1):
                minimum = int(order in required_orders)
                if observed < minimum:
                    return face, order, observed, minimum
    return None


def _first_local_macaulay_failure(
    face_table: Sequence[bool],
) -> tuple[int, int, int, int] | None:
    for face_size in (1, 2, 3):
        for face in _independent_masks(face_size):
            _link_vertices, g_entries = _link_g_entries(face_table, face)
            for previous_order in range(1, len(g_entries)):
                previous = g_entries[previous_order - 1]
                observed = g_entries[previous_order]
                if previous < 0:
                    return face, previous_order + 1, observed, -1
                upper = _macaulay_raise(previous, previous_order)
                if observed > upper:
                    return face, previous_order + 1, observed, upper
    return None


def _topology_audit(
    complements: Sequence[int],
    clutter: Sequence[int],
    missing_fours: Sequence[int],
    expected_vector: Sequence[int],
) -> dict[str, object]:
    """Full standard-library necessary topology audit of one completion."""

    selected = tuple(sorted(set(complements)))
    face_table = _face_table(selected)
    vector = _face_vector(face_table)
    minimal_nonfaces = _minimal_nonfaces_from_table(face_table)
    ridge_histogram = _ridge_histogram(selected)
    dual_components = _dual_component_count(selected)
    covers_triples = all(
        face_table[triple] for triple in _independent_masks(3)
    )
    interval_gate = _family_interval_feasible(
        tuple(missing_fours), _local_constraints()
    )
    strict_failure = _first_strict_local_g_failure(face_table)
    macaulay_failure = _first_local_macaulay_failure(face_table)
    euler_failure = _first_link_euler_failure(face_table)
    homology_failure = _first_mod2_homology_link_failure(face_table)
    packed = bytes(
        sum(bool(face_table[index + offset]) << offset for offset in range(8))
        for index in range(0, len(face_table), 8)
    )
    topology_valid = (
        vector == tuple(expected_vector)
        and set(minimal_nonfaces) == set(clutter)
        and ridge_histogram == {2: int(expected_vector[-2])}
        and dual_components == 1
        and covers_triples
        and interval_gate
        and strict_failure is None
        and macaulay_failure is None
        and euler_failure is None
        and homology_failure is None
    )
    if not topology_valid:
        raise AssertionError(
            "hostile full topology audit failed: "
            f"vector={vector}, ridge={ridge_histogram}, dual={dual_components}, "
            f"interval={interval_gate}, strict={strict_failure}, "
            f"Macaulay={macaulay_failure}, Euler={euler_failure}, "
            f"homology={homology_failure}"
        )
    return {
        "facet_count": len(selected),
        "face_vector": list(vector),
        "face_table_true_count": sum(face_table),
        "packed_face_table_sha256": sha256(packed).hexdigest(),
        "minimal_nonface_masks": list(minimal_nonfaces),
        "minimal_nonface_size_histogram": dict(
            sorted(Counter(face.bit_count() for face in minimal_nonfaces).items())
        ),
        "minimal_nonfaces_match_clutter": True,
        "largest_minimal_nonface_order": max(
            (face.bit_count() for face in minimal_nonfaces), default=0
        ),
        "lanner_order_bound_holds": all(
            face.bit_count() <= 5 for face in minimal_nonfaces
        ),
        "matches_dehn_sommerville_schedule": True,
        "ridge_degree_histogram": ridge_histogram,
        "dual_component_count": dual_components,
        "dual_graph_connected": True,
        "covers_all_independent_triples": covers_triples,
        "missing_four_interval_gate": interval_gate,
        "strict_local_g_constraints_hold": True,
        "local_macaulay_constraints_hold": True,
        "all_links_have_spherical_euler_characteristic": True,
        "all_links_are_mod2_homology_spheres": True,
        "passes_connected_homology_audit": True,
        "topology_valid": True,
    }


def _a0_census():
    _column_distributions.cache_clear()
    graph_covers = _minimal_graph_covers()
    five_covers = tuple(
        candidate
        for chosen in combinations(range(VERTEX_COUNT), 5)
        for candidate in (_mask(chosen),)
        if all(candidate & edge for edge in EDGE_MASKS)
    )
    all_word = (1 << len(five_covers)) - 1
    hit_words = {}
    _g3, base_g4 = _base_link_g3_g4()
    caps = tuple(value - 1 for value in base_g4)
    active_edges = tuple(
        edge for edge in EDGE_MASKS if edge & _mask(ACTIVE) == edge
    )
    active_masks = tuple(
        value
        for order in range(len(ACTIVE) + 1)
        for chosen in combinations(ACTIVE, order)
        for value in (_mask(chosen),)
        if all(value & edge != edge for edge in active_edges)
        and all(
            int(bool(value & (1 << vertex))) <= caps[vertex]
            for vertex in ACTIVE
        )
    )
    active_mask_histogram = dict(
        sorted(Counter(value.bit_count() for value in active_masks).items())
    )
    active_triples = tuple(
        _vertices(value) for value in active_masks if value.bit_count() == 3
    )
    if (
        active_masks != EXPECTED_ACTIVE_MASKS
        or active_mask_histogram != EXPECTED_ACTIVE_MASK_SIZE_HISTOGRAM
        or active_triples != EXPECTED_ACTIVE_TRIPLES
    ):
        raise AssertionError("hostile exact active-row alphabet changed")
    configurations = []
    for rows in combinations_with_replacement(active_masks, ROW_COUNT):
        if any(
            sum(bool(row & (1 << vertex)) for row in rows) > caps[vertex]
            for vertex in ACTIVE
        ):
            continue
        missed = tuple(
            sum(1 << row for row, active in enumerate(rows) if not active & cover)
            for cover in graph_covers
        )
        if any(
            all(
                not (pattern & (1 << row)) or rows[row] & (1 << vertex)
                for row in range(ROW_COUNT)
            )
            for pattern in missed
            for vertex in ACTIVE
        ):
            continue
        configurations.append((rows, missed))

    total = repeated = small = valid = 0
    histogram = Counter()
    targets = []
    path_hasher = sha256()
    for rows, missed in configurations:
        allowed = tuple(
            column
            for column in range(1 << ROW_COUNT)
            if column.bit_count() <= ROW_COUNT - 1
            and all(column & pattern != pattern for pattern in missed)
        )
        row_sums = tuple(5 - row.bit_count() for row in rows)
        for counts in _column_distributions(allowed, row_sums):
            total += 1
            family = _faces(rows, allowed, counts)
            path_hasher.update(
                b"".join(value.to_bytes(2, "little") for value in family)
            )
            if len(set(family)) != ROW_COUNT:
                repeated += 1
                path_hasher.update(b"R")
                continue
            if _has_small_transversal(family, graph_covers):
                small += 1
                path_hasher.update(b"S")
                continue
            valid += 1
            word = all_word
            for face in family:
                if face not in hit_words:
                    hit_words[face] = _hit_word(face, five_covers)
                word &= hit_words[face]
            facets = word.bit_count()
            path_hasher.update(b"V" + facets.to_bytes(2, "little"))
            histogram[facets] += 1
            if facets == 102:
                targets.append(family)
    families = tuple(sorted(set(targets)))
    family_digest = _canonical_sha256(
        [[list(_vertices(face)) for face in family] for family in families]
    )
    blockers = []
    ridge_histograms = []
    for family in families:
        transversals = _minimal_transversals((*EDGE_MASKS, *family))
        if Counter(value.bit_count() for value in transversals) != {5: 102}:
            raise AssertionError("hostile a0 target blocker is not pure")
        blockers.append(transversals)
        ridge_histograms.append(_ridge_histogram(transversals))
    orbit_keys = {_orbit_key(family) for family in families}
    ridge_degree_two_indices = tuple(
        index for index, ridge in enumerate(ridge_histograms) if ridge == {2: 459}
    )
    topology_audits = tuple(
        (
            index,
            _topology_audit(
                blockers[index],
                (*EDGE_MASKS, *families[index]),
                (),
                (14, 86, 309, 726, 1_171, 1_309, 986, 459, 102),
            ),
        )
        for index in ridge_degree_two_indices
    )
    topology_valid_indices = tuple(
        index for index, audit in topology_audits if audit["topology_valid"]
    )
    topology_valid_keys = {
        _orbit_key(families[index]) for index in topology_valid_indices
    }
    if (
        len(active_masks) != 18
        or len(configurations) != 1_221
        or _column_distributions.cache_info().currsize != 235
        or total != 375_582
        or repeated != 7_239
        or small != 0
        or valid != 368_343
        or len(targets) != 20
        or len(families) != 18
        or ridge_degree_two_indices != (0, 1, 2, 3, 4)
        or topology_valid_indices != (0, 1, 2, 3, 4)
        or len(orbit_keys) != 11
        or len(topology_valid_keys) != 3
        or EXPECTED_TARGET_FAMILY_SHA256 not in ("TO_BE_FROZEN", family_digest)
    ):
        raise AssertionError("hostile a0 census changed")
    return (
        {
            "independent_active_mask_count": len(active_masks),
            "independent_active_masks": list(active_masks),
            "independent_active_mask_vertices": [
                list(_vertices(value)) for value in active_masks
            ],
            "independent_active_mask_size_histogram": active_mask_histogram,
            "independent_active_triples": [list(face) for face in active_triples],
            "active_row_configuration_count": len(configurations),
            "complete_active_row_scope": (
                "all independent active subsets allowed by exact incidence "
                "caps; no singleton-row normalization"
            ),
            "column_distribution_problem_count": _column_distributions.cache_info().currsize,
            "total_distribution_path_count": total,
            "repeated_row_path_count": repeated,
            "small_transversal_path_count": small,
            "blocker_valid_path_count": valid,
            "facet_count_histogram": dict(sorted(histogram.items())),
            "path_transcript_sha256": path_hasher.hexdigest(),
            "target_path_count": len(targets),
            "exact_target_family_count": len(families),
            "target_family_sha256": family_digest,
            "ridge_histograms": ridge_histograms,
            "full_topology_audits": [
                {"family_index": index, **audit}
                for index, audit in topology_audits
            ],
            "topology_valid_family_indices": list(topology_valid_indices),
            "all_target_orbit_count": len(orbit_keys),
            "topology_valid_orbit_count": len(topology_valid_keys),
        },
        families,
    )


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
        (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 3), ((3, 4), 4)),
    ),
}


@lru_cache(maxsize=None)
def _lanner_codes(order: int) -> tuple[tuple[int, ...], ...]:
    edges = tuple(combinations(range(order), 2))
    codes = set()
    for representative in LANNER_REPRESENTATIVES[order]:
        base = {edge: 2 for edge in edges}
        base.update(representative)
        for permutation in permutations(range(order)):
            labels = {
                _pair(permutation[left], permutation[right]): label
                for (left, right), label in base.items()
            }
            codes.add(tuple(labels[edge] for edge in edges))
    expected = {4: 108, 5: 420}[order]
    if len(codes) != expected:
        raise AssertionError("hostile Lanner placement table changed")
    return tuple(sorted(codes))


def _path_labels(component, adjacency):
    endpoint = min(vertex for vertex in component if len(adjacency[vertex]) == 1)
    labels = []
    previous = None
    current = endpoint
    while True:
        following = [item for item in adjacency[current] if item[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("hostile finite path branched")
        next_vertex, label = following[0]
        labels.append(label)
        previous, current = current, next_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def _branch_arms(component, adjacency):
    branch = next(vertex for vertex in component if len(adjacency[vertex]) == 3)
    lengths = []
    for neighbor, _label in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(vertex for vertex, _ in adjacency[current] if vertex != previous)
            previous, current = current, following
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def _finite_component(component, labels):
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
    degrees = Counter(map(len, adjacency.values()))
    if max(degrees) <= 2:
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
    if degrees[3] != 1 or any(degree > 3 for degree in degrees):
        return False
    if any(
        labels[_pair(left, right)] != 3
        for left in component
        for right in component
        if left < right and labels[_pair(left, right)] > 2
    ):
        return False
    arms = _branch_arms(component, adjacency)
    return arms == (1, 1, order - 3) or arms in {(1, 2, 2), (1, 2, 3), (1, 2, 4)}


@lru_cache(maxsize=262_144)
def _finite_code(order: int, code: tuple[int, ...]) -> bool:
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
                other for other in unseen if labels[_pair(vertex, other)] > 2
            }
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        components.append(component)
    return all(_finite_component(component, labels) for component in components)


@lru_cache(maxsize=None)
def _placement_options(support: tuple[int, ...]):
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(
        tuple(zip(edges, code, strict=True)) for code in _lanner_codes(len(support))
    )


def _restriction_finite(restriction, labels):
    return _finite_code(
        len(restriction),
        tuple(labels[_pair(left, right)] for left, right in combinations(restriction, 2)),
    )


def _rank_nine_vertices(nonfaces):
    forbidden = tuple(frozenset(face) for face in nonfaces)
    return tuple(
        chosen
        for chosen in combinations(range(VERTEX_COUNT), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


def _state_digest(states):
    return _canonical_sha256(
        [
            [[left, right, label] for (left, right), label in state]
            for state in sorted(states)
        ]
    )


def _label_family(index, family):
    supports = tuple(_vertices(face) for face in family)
    actual = _rank_nine_vertices((*EDGES, *supports))
    remaining = set(range(len(supports)))
    union = set()
    states = ((),)
    stages = []
    while states and remaining:
        support_index = min(remaining) if not union else max(
            remaining,
            key=lambda value: (
                len(set(supports[value]) & union),
                -len(set(supports[value]) - union),
                -value,
            ),
        )
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {tuple(facet for facet in vertex if facet in union) for vertex in actual},
                key=lambda row: (len(row), row),
            )
        )
        passing = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in _placement_options(support):
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {
                    _pair(left, right): 2 for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                if all(_restriction_finite(restriction, labels) for restriction in restrictions):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "support": list(support),
                "raw": (stages[-1]["passing"] if stages else 1) * len(_placement_options(support)),
                "compatible": compatible,
                "passing": len(states),
                "state_sha256": _state_digest(states),
            }
        )
        remaining.remove(support_index)
    if not stages or stages[-1]["passing"] != 0:
        raise AssertionError(f"hostile a0 family {index} survived labels")
    return {
        "family_index": index,
        "rank_nine_vertex_count": len(actual),
        "stages": stages,
        "compact_no_go": True,
    }


def _exact_witness_state(support, type_index):
    """Fixed rooted representatives used by the exact two-block witnesses.

    The complete placement table above is generated from an independent
    hand-coded list of unrooted Lanner diagrams.  Its representative order is
    immaterial to that enumeration, but the exact leading-minor formulas need
    a particular rooted representative.  Spell those two rooted states out
    literally instead of silently identifying the two unrelated orderings.
    """

    support = tuple(sorted(support))
    patterns = {
        1: (
            ((0, 4), 3),
            ((1, 3), 3),
            ((1, 4), 3),
            ((2, 3), 5),
        ),
        4: (
            ((0, 3), 3),
            ((0, 4), 3),
            ((1, 2), 3),
            ((1, 4), 3),
            ((2, 3), 4),
        ),
    }
    labels = {
        _pair(support[left], support[right]): 2
        for left, right in combinations(range(len(support)), 2)
    }
    labels.update(
        {
            _pair(support[left], support[right]): label
            for (left, right), label in patterns[type_index]
        }
    )
    code = tuple(labels[edge] for edge in combinations(support, 2))
    if code not in _lanner_codes(5):
        raise AssertionError("hostile rooted exact witness is not Lanner")
    return labels


def _rooted_code(support, state, root):
    support = tuple(sorted(support))
    others = tuple(vertex for vertex in support if vertex != root)
    return min(
        tuple(
            state[_pair(order[left], order[right])]
            for left, right in combinations(range(5), 2)
        )
        for tail in permutations(others)
        for order in ((root, *tail),)
    )


def _bridge_key(left_support, right_support, left_state, right_state, edge):
    left_endpoint = next(vertex for vertex in edge if vertex in left_support)
    right_endpoint = next(vertex for vertex in edge if vertex in right_support)
    return tuple(
        sorted(
            (
                _rooted_code(left_support, left_state, left_endpoint),
                _rooted_code(right_support, right_state, right_endpoint),
            )
        )
    )


Surd = tuple[Fraction, Fraction, Fraction, Fraction]


def _s(a=0, b=0, c=0, d=0) -> Surd:
    return tuple(map(Fraction, (a, b, c, d)))  # type: ignore[return-value]


def _sadd(left: Surd, right: Surd) -> Surd:
    return tuple(a + b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def _sneg(value: Surd) -> Surd:
    return tuple(-entry for entry in value)  # type: ignore[return-value]


def _ssub(left: Surd, right: Surd) -> Surd:
    return _sadd(left, _sneg(right))


def _smul(left: Surd, right: Surd) -> Surd:
    answer = [Fraction(0) for _ in range(4)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            factor = 1
            if i & 1 and j & 1:
                factor *= 2
            if i & 2 and j & 2:
                factor *= 5
            answer[i ^ j] += factor * a * b
    return tuple(answer)  # type: ignore[return-value]


def _sconjugate(value: Surd, flip2: bool, flip5: bool) -> Surd:
    return tuple(
        coefficient
        * (-1 if flip2 and index & 1 else 1)
        * (-1 if flip5 and index & 2 else 1)
        for index, coefficient in enumerate(value)
    )  # type: ignore[return-value]


def _sinverse(value: Surd) -> Surd:
    numerator = _smul(
        _smul(_sconjugate(value, True, False), _sconjugate(value, False, True)),
        _sconjugate(value, True, True),
    )
    norm = _smul(value, numerator)
    if norm[1:] != (0, 0, 0) or not norm[0]:
        raise ArithmeticError("hostile biquadratic inverse failed")
    return tuple(coefficient / norm[0] for coefficient in numerator)  # type: ignore[return-value]


def _sdiv(left: Surd, right: Surd) -> Surd:
    return _smul(left, _sinverse(right))


def _determinant(matrix: Sequence[Sequence[Surd]]) -> Surd:
    rows = [list(row) for row in matrix]
    order = len(rows)
    determinant = _s(1)
    for column in range(order):
        pivot = next((row for row in range(column, order) if rows[row][column] != _s()), None)
        if pivot is None:
            return _s()
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            determinant = _sneg(determinant)
        pivot_value = rows[column][column]
        determinant = _smul(determinant, pivot_value)
        for row in range(column + 1, order):
            factor = _sdiv(rows[row][column], pivot_value)
            for later in range(column, order):
                rows[row][later] = _ssub(rows[row][later], _smul(factor, rows[column][later]))
    return determinant


def _twice_gram(order, labels):
    values = {
        2: _s(),
        3: _s(-1),
        4: _s(0, -1),
        5: _s(Fraction(-1, 2), 0, Fraction(-1, 2)),
    }
    return tuple(
        tuple(
            _s(2) if left == right else values[labels[_pair(left, right)]]
            for right in order
        )
        for left in order
    )


def _leading_minors(order, labels):
    matrix = _twice_gram(order, labels)
    return tuple(
        _determinant(tuple(row[:size] for row in matrix[:size]))
        for size in range(1, len(order) + 1)
    )


def _linear_sign(value: Surd) -> str:
    bounds = (
        (Fraction(1), Fraction(1)),
        (Fraction(7, 5), Fraction(3, 2)),
        (Fraction(11, 5), Fraction(9, 4)),
        (Fraction(79, 25), Fraction(19, 6)),
    )
    lower = upper = Fraction(0)
    for coefficient, (lo, hi) in zip(value, bounds, strict=True):
        if coefficient >= 0:
            lower += coefficient * lo
            upper += coefficient * hi
        else:
            lower += coefficient * hi
            upper += coefficient * lo
    if lower > 0:
        return "+"
    if upper < 0:
        return "-"
    raise AssertionError(f"hostile radical bounds do not decide {value}")


def _block_labels(left, right, left_state, right_state, bridge):
    labels = {
        _pair(a, b): 2 for a, b in combinations(tuple((*left, *right)), 2)
    }
    labels.update(left_state)
    labels.update(right_state)
    labels[bridge] = 3
    return labels


def _exact_index_two_witnesses():
    canonical_r = {
        index: _exact_witness_state(A1_R, index) for index in (1, 4)
    }
    canonical_b = {1: _exact_witness_state(A1_B, 1)}
    labels_a = _block_labels(A1_R, A1_B, canonical_r[1], canonical_b[1], _pair(A1_R[0], A1_B[0]))
    expected_a = (
        _s(2), _s(4), _s(8), _s(6, 0, -2), _s(4, 0, -2),
        _s(Fraction(9, 2), 0, Fraction(-5, 2)), _s(9, 0, -5),
        _s(18, 0, -10), _s(26, 0, -12),
        _s(Fraction(25, 2), 0, Fraction(-11, 2)),
    )
    labels_b = _block_labels(A1_R, A1_B, canonical_r[4], canonical_b[1], _pair(A1_R[4], A1_B[0]))
    expected_b = (
        _s(2), _s(4), _s(5, 0, -1),
        _s(Fraction(7, 2), 0, Fraction(-3, 2)), _s(4, 0, -2),
        _s(Fraction(9, 2), 0, Fraction(-5, 2)), _s(9, 0, -5),
        _s(9, 0, -5), _s(1, 0, -1),
        _s(Fraction(-23, 2), -8, Fraction(11, 2), 4),
    )
    rows = (
        ("A5_A5_type_1_1", tuple((*A1_R, *A1_B)), labels_a, expected_a),
        ("A5_A5_type_4_1", tuple(reversed(tuple((*A1_R, *A1_B)))), labels_b, expected_b),
    )
    answer = []
    for name, order, labels, expected in rows:
        observed = _leading_minors(order, labels)
        if observed != expected:
            raise AssertionError(f"hostile exact minors changed for {name}: {observed}")
        signs = "".join(_linear_sign(value) for value in observed)
        augmented = "+" + signs
        negative = sum(
            left != right
            for left, right in zip(augmented[:-1], augmented[1:], strict=True)
        )
        if negative != 2:
            raise AssertionError("hostile exact inertia is no longer index two")
        answer.append({"name": name, "signs": signs, "inertia": [8, 2, 0]})
    return tuple(answer)


def _merge(left, right):
    answer = dict(left)
    for edge, label in right.items():
        if edge in answer and answer[edge] != label:
            return None
        answer[edge] = label
    return answer


def _placement_maps(support):
    return tuple(dict(option) for option in _placement_options(tuple(support)))


def _cluster_restrictions(actual, cluster):
    support = set(cluster)
    return tuple(
        sorted(
            {tuple(facet for facet in vertex if facet in support) for vertex in actual},
            key=lambda row: (len(row), row),
        )
    )


def _a1_label_audit():
    actual = _rank_nine_vertices((*EDGES, A1_Q, A1_R, A1_B))
    q_states = _placement_maps(A1_Q)
    r_states = _placement_maps(A1_R)
    b_states = _placement_maps(A1_B)
    left_restrictions = _cluster_restrictions(actual, A1_LEFT)
    compatible = 0
    left_passing = set()
    for q_state in q_states:
        for r_state in r_states:
            joined = _merge(q_state, r_state)
            if joined is None:
                continue
            compatible += 1
            labels = {
                _pair(left, right): 2 for left, right in combinations(A1_LEFT, 2)
            }
            labels.update(joined)
            if all(_restriction_finite(row, labels) for row in left_restrictions):
                left_passing.add(tuple(sorted(joined.items())))
    left_states = tuple(dict(state) for state in sorted(left_passing))
    restrictions = _cluster_restrictions(actual, A1_UNION)
    by_edge = {
        edge: tuple(row for row in restrictions if set(edge) <= set(row))
        for edge in A1_CROSS
    }
    template = {
        _pair(left, right): 2 for left, right in combinations(A1_UNION, 2)
    }
    sparse_failures = 0
    passing = []
    passing_by_state = defaultdict(list)
    bad_sizes = Counter()
    for left_index, left_state in enumerate(left_states):
        for right_index, right_state in enumerate(b_states):
            base = dict(template)
            base.update(left_state)
            base.update(right_state)
            sparse_failures += not all(_restriction_finite(row, base) for row in restrictions)
            for edge in A1_CROSS:
                labels = dict(base)
                labels[edge] = 3
                bad = next((row for row in by_edge[edge] if not _restriction_finite(row, labels)), None)
                if bad is None:
                    passing.append((left_index, right_index, edge))
                    passing_by_state[(left_index, right_index)].append(edge)
                else:
                    bad_sizes[len(bad)] += 1
    stronger = []
    for label in (4, 5):
        count = 0
        for left_index, right_index, edge in passing:
            labels = dict(template)
            labels.update(left_states[left_index])
            labels.update(b_states[right_index])
            labels[edge] = label
            count += all(_restriction_finite(row, labels) for row in by_edge[edge])
        stronger.append((label, count))

    canonical_r = {
        index: _exact_witness_state(A1_R, index) for index in (1, 4)
    }
    canonical_b = {1: _exact_witness_state(A1_B, 1)}
    references = (
        _bridge_key(A1_R, A1_B, canonical_r[1], canonical_b[1], _pair(A1_R[0], A1_B[0])),
        _bridge_key(A1_R, A1_B, canonical_r[4], canonical_b[1], _pair(A1_R[4], A1_B[0])),
    )
    names = ("A5_A5_type_1_1", "A5_A5_type_4_1")
    classes = Counter()
    endpoints = Counter()
    transcript = []
    for left_index, right_index, edge in passing:
        r_endpoint = tuple(vertex for vertex in edge if vertex in A1_R)
        b_endpoint = tuple(vertex for vertex in edge if vertex in A1_B)
        if len(r_endpoint) != 1 or len(b_endpoint) != 1:
            raise AssertionError("hostile passing edge escaped R x B")
        key = _bridge_key(
            A1_R, A1_B, left_states[left_index], b_states[right_index], edge
        )
        matches = tuple(index for index, reference in enumerate(references) if key == reference)
        if len(matches) != 1:
            raise AssertionError("hostile exceptional class is not unique")
        name = names[matches[0]]
        classes[name] += 1
        endpoints[r_endpoint[0]] += 1
        transcript.append((left_index, right_index, edge, name))
    witnesses = _exact_index_two_witnesses()
    transcript_digest = _canonical_sha256(transcript)
    if (
        len(actual) != 100
        or len(q_states) != 108
        or (len(r_states), len(b_states)) != (420, 420)
        or compatible != 1_542
        or len(left_states) != 84
        or Counter(map(len, restrictions)) != {7: 20, 8: 55, 9: 15}
        or sparse_failures != 0
        or len(passing) != 5_040
        or Counter(map(len, passing_by_state.values())) != {1: 5_040}
        or bad_sizes != {7: 705_600, 8: 346_320, 9: 1_440}
        or tuple(stronger) != ((4, 0), (5, 0))
        or endpoints != {11: 2_520, 12: 2_520}
        or classes != {"A5_A5_type_1_1": 4_320, "A5_A5_type_4_1": 720}
        or EXPECTED_A1_LABEL_TRANSCRIPT_SHA256 not in ("TO_BE_FROZEN", transcript_digest)
    ):
        raise AssertionError("hostile a1 label gate changed")
    return {
        "rank_nine_vertex_count": len(actual),
        "placement_counts": [len(q_states), len(r_states), len(b_states)],
        "compatible_overlapping_join_count": compatible,
        "left_cluster_passing_count": len(left_states),
        "restriction_histogram": dict(sorted(Counter(map(len, restrictions)).items())),
        "sparse_state_count": len(left_states) * len(b_states),
        "sparse_failure_count": sparse_failures,
        "label_three_test_count": len(left_states) * len(b_states) * len(A1_CROSS),
        "label_three_passing_count": len(passing),
        "passing_edges_per_state_histogram": dict(sorted(Counter(map(len, passing_by_state.values())).items())),
        "bad_restriction_size_histogram": dict(sorted(bad_sizes.items())),
        "stronger_survivors": dict(stronger),
        "R_endpoint_histogram": dict(sorted(endpoints.items())),
        "exceptional_class_histogram": dict(sorted(classes.items())),
        "exact_witnesses": list(witnesses),
        "transcript_sha256": transcript_digest,
        "compact_no_go": True,
    }


def _a1_incidence_audit():
    clutter = (*EDGE_MASKS, _mask(A1_Q), _mask(A1_R), _mask(A1_B))
    transversals = _minimal_transversals(clutter)
    sizes = Counter(value.bit_count() for value in transversals)
    ridges = _ridge_histogram(transversals)
    if sizes != {5: 100} or ridges != {2: 450}:
        raise AssertionError("hostile a1 incidence audit changed")
    topology = _topology_audit(
        transversals,
        clutter,
        (_mask(A1_Q),),
        (14, 86, 309, 725, 1_165, 1_295, 970, 450, 100),
    )
    return {
        "minimal_transversal_size_histogram": dict(sorted(sizes.items())),
        "ridge_degree_histogram": ridges,
        "facet_complement_sha256": sha256(
            b"".join(value.to_bytes(2, "little") for value in transversals)
        ).hexdigest(),
        "full_topology_audit": topology,
    }


@lru_cache(maxsize=1)
def replay() -> dict[str, object]:
    schedules = _schedules()
    positive = _positive_branches()
    a0, families = _a0_census()
    a0_labels = tuple(_label_family(index, family) for index, family in enumerate(families))
    a0_label_digest = _canonical_sha256(a0_labels)
    if EXPECTED_A0_LABEL_TRANSCRIPT_SHA256 not in ("TO_BE_FROZEN", a0_label_digest):
        raise AssertionError("hostile a0 label transcript changed")
    a1_incidence = _a1_incidence_audit()
    a1_labels = _a1_label_audit()
    result = {
        "schedules": [list(row) for row in schedules],
        "strict_schedule_survivors": [0, 1, 2],
        "positive_branches": positive,
        "a0_census": a0,
        "a0_all_target_label_audits": list(a0_labels),
        "a0_label_transcript_sha256": a0_label_digest,
        "a1_incidence_audit": a1_incidence,
        "a1_label_audit": a1_labels,
        "theorem": {
            "exact_topology_valid_residual_incidence_orbit_count": (
                a0["topology_valid_orbit_count"]
                + int(
                    a1_incidence["full_topology_audit"]["topology_valid"]
                )
            ),
            "whole_graph_shadow_compact_coxeter_excluded": True,
            "finite_volume_noncompact_chamber_excluded": False,
        },
    }
    if (
        result["theorem"][
            "exact_topology_valid_residual_incidence_orbit_count"
        ]
        != 4
    ):
        raise AssertionError("hostile topology-valid residual orbit count changed")
    result["summary_sha256"] = _canonical_sha256(result)
    if EXPECTED_SUMMARY_SHA256 not in ("TO_BE_FROZEN", result["summary_sha256"]):
        raise AssertionError("hostile record-35 summary changed")
    return result


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2, sort_keys=True))
