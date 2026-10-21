#!/usr/bin/env python3
"""No-project-import replay of the record-42 whole-shadow compact no-go.

This verifier uses only the Python standard library.  It independently
reconstructs the five schedules, the ``a=1,2`` shadow gates, the cover-pruned
``a=0`` blocker census, the five ridge-valid target families and their exact
``D5 x S9 x S5`` row/column orbit.  Most importantly, it regenerates all 420
labelled order-five Lanner diagrams and obtains zero compact label states for
each of the fifteen exact ``a=0`` target families.

The selected ``a=1`` product completion is linked to the re-audited
principal-block theorem by file hash, internal certificate digest, exact
semantic projection, and the explicit full support-symmetry orders.  No
finite-volume noncompact claim is made.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations
import json
from math import comb
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
COXETER_ARCHIVE = ROOT / "data/enumerations/h9_14_selected_nerve_coxeter_no_go.json"

VERTEX_COUNT = 14
DIMENSION = 9
FULL_MASK = (1 << VERTEX_COUNT) - 1
GRAPH6 = "M?????????O?o?O@?"
EDGES = ((0, 11), (0, 12), (1, 12), (1, 13), (11, 13))
ACTIVE = (0, 1, 11, 12, 13)
ACTIVE_CYCLE = (0, 11, 13, 1, 12)
ISOLATES = tuple(range(2, 11))
ROW_COUNT = 5
PRODUCT_FOUR = (2, 3, 4, 5)
PRODUCT_FIVE = (6, 7, 8, 9, 10)
REPRESENTATIVE_FIVES = (
    (0, 3, 4, 5, 6),
    (1, 3, 4, 5, 6),
    (2, 3, 4, 5, 6),
    (2, 7, 8, 9, 10),
    (7, 8, 9, 10, 12),
)

EXPECTED_COXETER_FILE_SHA256 = (
    "028bbc8e266145e59fa1ee9f9796a6e2c241f4a4f09637bd73f98bced9908629"
)
EXPECTED_COXETER_CERTIFICATE_SHA256 = (
    "c2dfcbfc03dc5dbae34938ab32f3a1c1111fd729d7844fd76d04dfc6cfa44de8"
)
EXPECTED_RECORD42_PROJECTION_SHA256 = (
    "9238c596b14426d120c25d8041c770b6dcce1dbf2526c2e2ea43c04894f66598"
)
EXPECTED_COXETER_KIND = (
    "h9_14_seven_selected_nerve_compact_coxeter_no_go_reaudited"
)
EXPECTED_SUPPORT_SYMMETRY_ORDERS = [24, 120]
EXPECTED_SCOPE_REAUDIT_STATUS = (
    "hostile_false_alarm_rejected_by_principal_block_proof"
)
EXPECTED_TARGET_FAMILY_SHA256 = (
    "27412b42ae9cbe2c09e915c70b2986f3d5ca4b90a11b758d0df674727534d9d5"
)
EXPECTED_ORBIT_KEY_SHA256 = (
    "05c1d904e9aab3257c140300f94b7a3d3d13029736f70ab730db0b286cd0b9e2"
)
EXPECTED_SUMMARY_SHA256 = (
    "55b0206407df3aa6c3d9c14a76070db8a03f165a06798b670b148eb9ce5e30e6"
)


def _mask(values: Iterable[int]) -> int:
    return sum(1 << value for value in values)


def _vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(VERTEX_COUNT) if value >> vertex & 1)


EDGE_MASKS = tuple(_mask(edge) for edge in EDGES)


def _canonical_sha256(value: object) -> str:
    return sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode()
    ).hexdigest()


def _require_frozen_digest(expected: str, observed: str, message: str) -> None:
    if expected not in ("TO_BE_FROZEN", observed):
        raise ValueError(message)


def _is_independent(value: int) -> bool:
    return all(value & edge != edge for edge in EDGE_MASKS)


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


@lru_cache(maxsize=None)
def _independent_masks(order: int) -> tuple[int, ...]:
    return tuple(
        value
        for chosen in combinations(range(VERTEX_COUNT), order)
        for value in (_mask(chosen),)
        if _is_independent(value)
    )


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
            (a, independent_fives - f_vector[4], f_vector[-1], _summed_link_g4(f_vector))
        )
    expected = (
        (0, 5, 102, 30),
        (1, 11, 100, 24),
        (2, 17, 98, 18),
        (3, 23, 96, 12),
        (4, 29, 94, 6),
    )
    if tuple(answers) != expected:
        raise AssertionError("record-42 hostile schedule census changed")
    return tuple(answers)


def _base_link_g3_g4() -> tuple[tuple[int, ...], tuple[int, ...]]:
    independent = {order: _independent_masks(order) for order in range(2, 6)}
    g3 = []
    g4 = []
    for vertex in range(VERTEX_COUNT):
        counts = {
            order: sum(bool(value & (1 << vertex)) for value in rows)
            for order, rows in independent.items()
        }
        f0, f1, f2, f3 = (counts[order] for order in range(2, 6))
        g3.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    return tuple(g3), tuple(g4)


def _positive_branches() -> dict[str, object]:
    # The independently reconstructed local intervals say exactly: an A4
    # contains at most one active C5 vertex.  In a family, each active vertex
    # can occur in at most one A4.  This gives 126+5*84=546 individual A4s
    # and removes 5*C(84,2)=17430 of their unordered pairs.
    fours = tuple(
        value
        for value in _independent_masks(4)
        if (value & _mask(ACTIVE)).bit_count() <= 1
    )
    fives = _independent_masks(5)
    base_g3, base_g4 = _base_link_g3_g4()

    def shadow(family: Sequence[int]) -> tuple[int, ...]:
        return tuple(
            five for five in fives if any(five & four == four for four in family)
        )

    def strict(family: Sequence[int], forced: Sequence[int]) -> tuple[bool, tuple[int, ...]]:
        upper = []
        okay = True
        for vertex in range(VERTEX_COUNT):
            a4 = sum(bool(face & (1 << vertex)) for face in family)
            previous = base_g3[vertex] - a4
            a5 = sum(bool(face & (1 << vertex)) for face in forced)
            before = base_g4[vertex] + 6 * a4 - a5
            okay &= previous >= 1 and before >= 1
            upper.append(before - 1)
        return okay, tuple(upper)

    a1_hist: Counter[int] = Counter()
    a1_survivors = []
    for four in fours:
        forced = shadow((four,))
        a1_hist[len(forced)] += 1
        okay, upper = strict((four,), forced)
        if not okay:
            continue
        candidates = tuple(
            five
            for five in fives
            if five not in forced
            and all(
                not (five & (1 << vertex)) or upper[vertex] >= 1
                for vertex in range(VERTEX_COUNT)
            )
        )
        if len(candidates) != 1:
            raise AssertionError("record-42 hostile a1 completion is not forced")
        a1_survivors.append((four, candidates[0]))
    expected_a1 = {
        (_mask(four), _mask(set(ISOLATES) - set(four)))
        for four in combinations(ISOLATES, 4)
    }
    if set(a1_survivors) != expected_a1:
        raise AssertionError("record-42 hostile a1 product orbit changed")

    pair_count = 0
    local_count = 0
    overrun = 0
    strict_count = 0
    a2_hist: Counter[int] = Counter()
    for first, second in combinations(fours, 2):
        pair_count += 1
        if any(
            first & (1 << active) and second & (1 << active)
            for active in ACTIVE
        ):
            continue
        local_count += 1
        forced = shadow((first, second))
        a2_hist[len(forced)] += 1
        overrun += len(forced) > 17
        strict_count += strict((first, second), forced)[0]
    expected_a2_hist = {
        15: 420, 16: 70_140, 17: 2_520,
        18: 50_400, 19: 1_260, 20: 6_615,
    }
    if (
        len(fours) != 546
        or pair_count != 148_785
        or local_count != 131_355
        or dict(sorted(a2_hist.items())) != expected_a2_hist
        or overrun != 58_275
        or strict_count != 0
    ):
        raise AssertionError("record-42 hostile a2 census changed")
    return {
        "individual_a4_count": len(fours),
        "a1_shadow_histogram": dict(sorted(a1_hist.items())),
        "a1_strict_count": len(a1_survivors),
        "a1_product_orbit_count": 1,
        "a1_representative": {
            "missing_four": list(PRODUCT_FOUR),
            "added_minimal_five": list(PRODUCT_FIVE),
        },
        "a2_pair_count": pair_count,
        "a2_family_interval_count": local_count,
        "a2_shadow_histogram": dict(sorted(a2_hist.items())),
        "a2_overrun_count": overrun,
        "a2_strict_count": strict_count,
    }


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
    result = tuple(answers)
    if result != (2051, 6146, 8195, 12289, 14336):
        raise AssertionError("record-42 hostile graph covers changed")
    return result


def _local_a5_interval(subset: tuple[int, ...]) -> tuple[int, int]:
    size = len(subset)
    link_dimension = DIMENSION - size
    high_order = 5 - size
    rest = tuple(vertex for vertex in range(VERTEX_COUNT) if vertex not in subset)
    link_f = {-1: 1}
    for order in range(1, high_order + 1):
        link_f[order - 1] = sum(
            _is_independent(_mask(subset + extension))
            for extension in combinations(rest, order)
        )
    h = tuple(
        sum(
            (-1) ** (index - order)
            * comb(link_dimension - order, link_dimension - index)
            * link_f[order - 1]
            for order in range(index + 1)
        )
        for index in range(high_order + 1)
    )
    g = tuple(h[index] - h[index - 1] for index in range(1, len(h)))
    previous, baseline = g[-2], g[-1]
    return (
        max(0, baseline - _macaulay_raise(previous, high_order - 1)),
        baseline - 1,
    )


def _active_row_configurations() -> tuple[tuple[int, ...], ...]:
    active_masks = tuple(
        value
        for order in range(3)
        for chosen in combinations(ACTIVE, order)
        for value in (_mask(chosen),)
        if _is_independent(value)
    )
    configurations = tuple(
        rows
        for rows in combinations_with_replacement(active_masks, ROW_COUNT)
        if all(not (left & right) for left, right in combinations(rows, 2))
    )
    if len(active_masks) != 11 or len(configurations) != 82:
        raise AssertionError("hostile record-42 active-row quotient changed")
    if tuple(_local_a5_interval((vertex,)) for vertex in ACTIVE) != ((0, 1),) * 5:
        raise AssertionError("hostile active singleton interval changed")
    if tuple(_local_a5_interval((vertex,)) for vertex in ISOLATES) != ((0, 4),) * 9:
        raise AssertionError("hostile isolate singleton interval changed")
    return configurations


@lru_cache(maxsize=None)
def _column_distributions(
    allowed: tuple[int, ...], row_sums: tuple[int, ...]
) -> tuple[tuple[int, ...], ...]:
    counts = [0] * len(allowed)
    answers = []

    def visit(position: int, left: int, sums: tuple[int, ...]) -> None:
        if position == len(allowed):
            if left == 0 and not any(sums):
                answers.append(tuple(counts))
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
                visit(position + 1, next_left, next_sums)
        counts[position] = 0

    visit(0, len(ISOLATES), row_sums)
    return tuple(answers)


def _faces(
    active_rows: Sequence[int], allowed: Sequence[int], counts: Sequence[int]
) -> tuple[int, ...]:
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


def _minimal_transversals(hyperedges: Sequence[int]) -> tuple[int, ...]:
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


def _has_small_transversal(
    fives: Sequence[int], graph_covers: Sequence[int]
) -> bool:
    """Detect a minimal transversal of order below five.

    Every transversal of the five graph edges contains a minimal graph cover,
    so testing the five fixed covers against the candidate A5 rows is exact.
    """

    return any(all(cover & five for five in fives) for cover in graph_covers)


def _ridge_histogram(complements: Sequence[int]) -> dict[int, int]:
    facets = tuple(FULL_MASK ^ complement for complement in complements)
    ridges: Counter[int] = Counter()
    for facet in facets:
        remaining = facet
        while remaining:
            bit = remaining & -remaining
            ridges[facet ^ bit] += 1
            remaining ^= bit
    return dict(sorted(Counter(ridges.values()).items()))


def _cycle_automorphisms() -> tuple[tuple[int, ...], ...]:
    answers = []
    for direction in (1, -1):
        for shift in range(5):
            image = list(range(VERTEX_COUNT))
            for index, source in enumerate(ACTIVE_CYCLE):
                image[source] = ACTIVE_CYCLE[(shift + direction * index) % 5]
            answers.append(tuple(image))
    return tuple(answers)


def _orbit_key(family: Sequence[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    active_mask = _mask(ACTIVE)
    keys = []
    for automorphism in _cycle_automorphisms():
        transformed = tuple(
            sum(1 << automorphism[vertex] for vertex in _vertices(face))
            for face in family
        )
        for order in permutations(range(ROW_COUNT)):
            active_rows = tuple(transformed[index] & active_mask for index in order)
            columns = tuple(
                sorted(
                    sum(
                        bool(transformed[order[row]] & (1 << isolate)) << row
                        for row in range(ROW_COUNT)
                    )
                    for isolate in ISOLATES
                )
            )
            keys.append((active_rows, columns))
    return min(keys)


def _a0_census() -> tuple[dict[str, object], tuple[tuple[int, ...], ...]]:
    graph_covers = _minimal_graph_covers()
    five_covers = tuple(
        value
        for chosen in combinations(range(VERTEX_COUNT), 5)
        for value in (_mask(chosen),)
        if all(value & edge for edge in EDGE_MASKS)
    )
    all_five = (1 << len(five_covers)) - 1
    hit_words: dict[int, int] = {}
    paths = repeated = small_transversal_failures = blocker = 0
    histogram: Counter[int] = Counter()
    targets = []
    active_configurations = _active_row_configurations()
    for active_rows in active_configurations:
        missed_rows = tuple(
            sum(1 << row for row, active in enumerate(active_rows) if not active & cover)
            for cover in graph_covers
        )
        allowed = tuple(
            column
            for column in range(1 << 5)
            if column.bit_count() <= 4
            and all(column & missed != missed for missed in missed_rows)
        )
        row_sums = tuple(5 - row.bit_count() for row in active_rows)
        for counts in _column_distributions(allowed, row_sums):
            paths += 1
            family = _faces(active_rows, allowed, counts)
            if len(set(family)) != 5:
                repeated += 1
                continue
            if _has_small_transversal(family, graph_covers):
                small_transversal_failures += 1
                continue
            blocker += 1
            word = all_five
            for face in family:
                if face not in hit_words:
                    hit_words[face] = sum(
                        bool(face & cover) << index
                        for index, cover in enumerate(five_covers)
                    )
                word &= hit_words[face]
            facet_count = word.bit_count()
            histogram[facet_count] += 1
            if facet_count == 102:
                targets.append(family)
    families = tuple(sorted(set(targets)))
    family_digest = _canonical_sha256(
        [[list(_vertices(face)) for face in family] for family in families]
    )
    if (
        paths != 1_098_455
        or repeated != 28_805
        or small_transversal_failures != 0
        or blocker != 1_069_650
        or len(targets) != 35
        or len(families) != 15
        or family_digest != EXPECTED_TARGET_FAMILY_SHA256
    ):
        raise AssertionError("record-42 hostile a0 blocker census changed")
    ridge_valid = []
    ridge_histograms = []
    for index, family in enumerate(families):
        transversals = _minimal_transversals(EDGE_MASKS + family)
        size_histogram = Counter(value.bit_count() for value in transversals)
        if size_histogram != {5: 102}:
            raise AssertionError("record-42 target blocker lost purity")
        ridge = _ridge_histogram(transversals)
        ridge_histograms.append(ridge)
        if ridge == {2: 459}:
            ridge_valid.append(index)
        elif ridge != {1: 8, 2: 443, 3: 8}:
            raise AssertionError("record-42 target ridge histogram changed")
    orbit_keys = {_orbit_key(families[index]) for index in ridge_valid}
    orbit_digest = _canonical_sha256(tuple(sorted(orbit_keys)))
    if tuple(ridge_valid) != (0, 1, 2, 3, 6) or len(orbit_keys) != 1:
        raise AssertionError("record-42 hostile topology orbit changed")
    if orbit_digest != EXPECTED_ORBIT_KEY_SHA256:
        raise AssertionError("record-42 hostile topology orbit digest changed")
    if tuple(sorted(_mask(face) for face in REPRESENTATIVE_FIVES)) != families[0]:
        raise AssertionError("record-42 named representative changed")
    return (
        {
            "active_row_configuration_count": len(active_configurations),
            "cover_prefilter_paths": paths,
            "repeated_row_paths": repeated,
            "small_transversal_failure_count": small_transversal_failures,
            "blocker_valid_paths": blocker,
            "facet_count_histogram": dict(sorted(histogram.items())),
            "target_path_count": len(targets),
            "exact_target_family_count": len(families),
            "target_family_sha256": family_digest,
            "ridge_histograms": ridge_histograms,
            "ridge_valid_family_indices": ridge_valid,
            "ridge_valid_orbit_count": len(orbit_keys),
            "orbit_key_sha256": orbit_digest,
        },
        families,
    )


# Five classical order-five Lanner representatives.  Unlisted pairs are 2.
LANNER5_REPRESENTATIVES = (
    (((0, 1), 3), ((0, 2), 3), ((0, 3), 3), ((1, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 4), ((2, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 5), ((2, 4), 5)),
    (((0, 1), 3), ((0, 2), 3), ((1, 3), 3), ((2, 4), 3), ((3, 4), 4)),
)


def _pair(left: int, right: int) -> tuple[int, int]:
    return (left, right) if left < right else (right, left)


@lru_cache(maxsize=1)
def _lanner5_codes() -> tuple[tuple[int, ...], ...]:
    codes = set()
    edges = tuple(combinations(range(5), 2))
    for representative in LANNER5_REPRESENTATIVES:
        base = {edge: 2 for edge in edges}
        base.update(representative)
        for permutation in permutations(range(5)):
            labels = {
                _pair(permutation[left], permutation[right]): label
                for (left, right), label in base.items()
            }
            codes.add(tuple(labels[edge] for edge in edges))
    if len(codes) != 420:
        raise AssertionError("hostile order-five Lanner table changed")
    return tuple(sorted(codes))


def _path_labels(component: set[int], adjacency: Mapping[int, tuple[tuple[int, int], ...]]) -> tuple[int, ...]:
    endpoint = min(vertex for vertex in component if len(adjacency[vertex]) == 1)
    labels = []
    previous = None
    current = endpoint
    while True:
        following = [row for row in adjacency[current] if row[0] != previous]
        if not following:
            break
        if len(following) != 1:
            raise AssertionError("hostile finite path branched")
        next_vertex, label = following[0]
        labels.append(label)
        previous, current = current, next_vertex
    forward = tuple(labels)
    return min(forward, tuple(reversed(forward)))


def _arms(component: set[int], adjacency: Mapping[int, tuple[tuple[int, int], ...]]) -> tuple[int, ...]:
    branch = next(vertex for vertex in component if len(adjacency[vertex]) == 3)
    lengths = []
    for neighbor, _ in adjacency[branch]:
        length = 1
        previous, current = branch, neighbor
        while len(adjacency[current]) == 2:
            following = next(vertex for vertex, _ in adjacency[current] if vertex != previous)
            previous, current = current, following
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def _finite_component(component: set[int], labels: Mapping[tuple[int, int], int]) -> bool:
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
        sequence = _path_labels(component, adjacency)
        order = len(component)
        if order == 2:
            return sequence[0] in (3, 4, 5)
        return (
            all(label == 3 for label in sequence)
            or sequence == (3,) * (order - 2) + (4,)
            or (order == 3 and sequence == (3, 5))
            or (order == 4 and sequence in ((3, 3, 5), (3, 4, 3)))
        )
    if degrees[3] != 1 or any(degree > 3 for degree in degrees):
        return False
    if any(
        labels[_pair(left, right)] != 3
        for left in component
        for right in component
        if left < right and labels[_pair(left, right)] > 2
    ):
        return False
    arms = _arms(component, adjacency)
    return arms == (1, 1, len(component) - 3) or arms in {
        (1, 2, 2), (1, 2, 3), (1, 2, 4)
    }


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


def _rank_nine(nonfaces: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    forbidden = tuple(set(face) for face in nonfaces)
    return tuple(
        chosen
        for chosen in combinations(range(VERTEX_COUNT), DIMENSION)
        if not any(face <= set(chosen) for face in forbidden)
    )


def _state_digest(states: Sequence[tuple[tuple[tuple[int, int], int], ...]]) -> str:
    return _canonical_sha256(
        [
            [[left, right, label] for (left, right), label in state]
            for state in sorted(states)
        ]
    )


def _label_family(index: int, family: Sequence[int]) -> dict[str, object]:
    supports = tuple(_vertices(face) for face in family)
    rank_nine = _rank_nine(EDGES + supports)
    remaining = set(range(5))
    union: set[int] = set()
    states: tuple[tuple[tuple[tuple[int, int], int], ...], ...] = ((),)
    stages = []
    while states and remaining:
        support_index = (
            min(remaining)
            if not union
            else max(
                remaining,
                key=lambda item: (
                    len(set(supports[item]) & union),
                    -len(set(supports[item]) - union),
                    -item,
                ),
            )
        )
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {tuple(vertex for vertex in face if vertex in union) for face in rank_nine},
                key=lambda row: (len(row), row),
            )
        )
        support_edges = tuple(combinations(support, 2))
        options = tuple(
            tuple(zip(support_edges, code, strict=True)) for code in _lanner5_codes()
        )
        old_count = len(states)
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
                    _pair(left, right): 2 for left, right in combinations(cluster, 2)
                }
                labels.update(merged)
                good = True
                for restriction in restrictions:
                    code = tuple(
                        labels[_pair(left, right)]
                        for left, right in combinations(restriction, 2)
                    )
                    if not _finite_code(len(restriction), code):
                        good = False
                        break
                if good:
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        stages.append(
            {
                "support_index": support_index,
                "raw_extension_count": old_count * 420,
                "compatible_extension_count": compatible,
                "passing_state_count": len(states),
                "passing_state_sha256": _state_digest(states),
            }
        )
        remaining.remove(support_index)
    if not stages or stages[-1]["passing_state_count"] != 0:
        raise AssertionError(f"hostile family {index} survived the full label gate")
    return {
        "family_index": index,
        "rank_nine_vertex_count": len(rank_nine),
        "stages": stages,
        "unrestricted_compact_coxeter_no_go": True,
    }


def _selected_product_chain() -> dict[str, object]:
    raw = COXETER_ARCHIVE.read_bytes()
    file_digest = sha256(raw).hexdigest()
    _require_frozen_digest(
        EXPECTED_COXETER_FILE_SHA256,
        file_digest,
        "selected-nerve archive file hash changed",
    )
    document = json.loads(raw)
    payload = dict(document)
    embedded = payload.pop("certificate_sha256")
    digest = _canonical_sha256(payload)
    if embedded != digest:
        raise ValueError("selected-nerve archive digest changed")
    _require_frozen_digest(
        EXPECTED_COXETER_CERTIFICATE_SHA256,
        digest,
        "selected-nerve archive certificate changed",
    )
    scope_reaudit = document.get("scope_reaudit")
    if (
        document.get("kind") != EXPECTED_COXETER_KIND
        or not isinstance(scope_reaudit, dict)
        or scope_reaudit.get("status") != EXPECTED_SCOPE_REAUDIT_STATUS
        or 42 not in scope_reaudit.get("reaudited_indices", [])
        or scope_reaudit.get("normalizing_support_symmetry_orders", {}).get("42")
        != EXPECTED_SUPPORT_SYMMETRY_ORDERS
        or "principal submatrix" not in scope_reaudit.get(
            "principal_block_lemma", ""
        )
        or "inertia" not in scope_reaudit.get("principal_block_lemma", "")
    ):
        raise ValueError("selected-nerve principal-block scope re-audit changed")
    record = next(
        row for row in document["six_nerve_census"] if row["record_index"] == 42
    )
    projection_digest = _canonical_sha256(record)
    gate = record["support_gates"][0]
    if (
        gate["left_support"] != list(PRODUCT_FOUR)
        or gate["right_support"] != list(PRODUCT_FIVE)
        or gate["single_edge_test_count"] != 900
        or gate["single_edge_passing_count"] != 2
        or record.get("normalizing_support_symmetry_orders")
        != EXPECTED_SUPPORT_SYMMETRY_ORDERS
        or record["raw_two_edge_trial_count"] != 8_550
        or record["monotonicity_residual_two_edge_trial_count"] != 1
        or record["two_edge_passing_count"] != 0
        or record["stronger_single_edge_survivors"] != {"4": 0, "5": 0}
        or record["exact_index_two_witnesses"][0]["inertia"] != [7, 2, 0]
        or record["unrestricted_compact_coxeter_no_go"] is not True
    ):
        raise ValueError("selected record-42 unrestricted theorem changed")
    _require_frozen_digest(
        EXPECTED_RECORD42_PROJECTION_SHA256,
        projection_digest,
        "selected record-42 projection changed",
    )
    return {
        "file_sha256": file_digest,
        "certificate_sha256": digest,
        "record42_projection_sha256": projection_digest,
        "aggregate_kind": document["kind"],
        "scope_reaudit_status": scope_reaudit["status"],
        "principal_block_lemma": scope_reaudit["principal_block_lemma"],
        "normalizing_support_symmetry_orders": (
            record["normalizing_support_symmetry_orders"]
        ),
        "single_edge_test_count": gate["single_edge_test_count"],
        "raw_two_edge_trial_count": record["raw_two_edge_trial_count"],
        "two_edge_passing_count": record["two_edge_passing_count"],
        "exceptional_inertia": record["exact_index_two_witnesses"][0]["inertia"],
        "unrestricted_compact_coxeter_no_go": True,
    }


def replay() -> dict[str, object]:
    schedules = _schedules()
    positive = _positive_branches()
    a0, families = _a0_census()
    labels = tuple(_label_family(index, family) for index, family in enumerate(families))
    label_transcript = _canonical_sha256(labels)
    product = _selected_product_chain()
    result: dict[str, object] = {
        "graph6": GRAPH6,
        "missing_edges": [list(edge) for edge in EDGES],
        "schedules": [list(row) for row in schedules],
        "strict_schedule_survivors": [0, 1, 2],
        "positive_branches": positive,
        "a0_census": a0,
        "a0_all_target_label_audits": list(labels),
        "a0_label_stage_transcript_sha256": label_transcript,
        "a1_product_label_chain": product,
        "theorem": {
            "whole_record42_graph_shadow_compact_coxeter_impossible": True,
            "a0_exact_target_families_excluded": len(labels),
            "a1_product_orbit_excluded": True,
            "finite_volume_noncompact_excluded": False,
        },
    }
    result["summary_sha256"] = _canonical_sha256(result)
    if EXPECTED_SUMMARY_SHA256 not in ("TO_BE_FROZEN", result["summary_sha256"]):
        raise ValueError("record-42 hostile summary is not frozen")
    return result


def main() -> None:
    print(json.dumps(replay(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
