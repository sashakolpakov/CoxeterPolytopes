#!/usr/bin/env python3
"""Standard-library-only a=1,2,3 completion replay for H9/14 record 32.

The a=1 branch uses the full fixed-A4 stabilizer quotient: every independent
active row part and every isolate-column multiset in the two genuine isolate
classes is enumerated.  The a=2 and a=3 branches exhaust all additional A5
choices.  Literal small-transversal, five-cover, blocker, and ridge tests are
implemented locally; no project module or producer table is imported.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations
import json
from math import comb


N = 14
FULL_MASK = (1 << N) - 1
MISSING_EDGES = ((0, 11), (1, 12), (0, 13), (2, 13))
EDGE_MASKS = tuple(sum(1 << vertex for vertex in edge) for edge in MISSING_EDGES)
ACTIVE_VERTICES = (0, 1, 2, 11, 12, 13)
ISOLATES = tuple(range(3, 11))
ACTIVE_MASK = sum(1 << vertex for vertex in ACTIVE_VERTICES)

A1_REPRESENTATIVES = (
    (3, 4, 5, 6),
    (2, 3, 4, 5),
    (1, 2, 3, 4),
    (2, 3, 4, 11),
)
A2_REPRESENTATIVES = (
    ((2, 3, 4, 11), (5, 6, 7, 8)),
    ((2, 3, 4, 5), (3, 4, 5, 11)),
    ((1, 2, 3, 4), (2, 3, 4, 11)),
    ((1, 2, 3, 4), (3, 4, 11, 12)),
    ((1, 2, 3, 4), (3, 5, 11, 12)),
    ((1, 2, 3, 4), (5, 6, 11, 12)),
)
A3_REPRESENTATIVE = (
    (1, 2, 3, 4),
    (2, 3, 4, 11),
    (3, 4, 11, 12),
)


def mask(values) -> int:
    return sum(1 << value for value in values)


def vertices(value: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(N) if value & (1 << vertex))


@lru_cache(maxsize=None)
def independent_masks(size: int) -> tuple[int, ...]:
    return tuple(
        value
        for chosen in combinations(range(N), size)
        for value in (mask(chosen),)
        if all(value & edge != edge for edge in EDGE_MASKS)
    )


@lru_cache(maxsize=None)
def graph_cover_masks(size: int) -> tuple[int, ...]:
    return tuple(
        value
        for chosen in combinations(range(N), size)
        for value in (mask(chosen),)
        if all(value & edge for edge in EDGE_MASKS)
    )


def binomial(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def macaulay_raise(value: int, degree: int) -> int:
    remaining = value
    next_top = -1
    raised = 0
    for lower in range(degree, 0, -1):
        top = lower - 1
        while (next_top < 0 or top + 1 < next_top) and binomial(top + 1, lower) <= remaining:
            top += 1
        remaining -= binomial(top, lower)
        raised += binomial(top + 1, lower + 1)
        next_top = top
    if remaining:
        raise AssertionError("inexact Macaulay expansion")
    return raised


@lru_cache(maxsize=1)
def base_link_g3_g4() -> tuple[tuple[int, ...], tuple[int, ...]]:
    independent = {size: independent_masks(size) for size in range(2, 6)}
    g3 = []
    g4 = []
    for vertex in range(N):
        counts = {
            size: sum(bool(face & (1 << vertex)) for face in independent[size])
            for size in independent
        }
        f0, f1, f2, f3 = (counts[size] for size in range(2, 6))
        g3.append(f2 - 7 * f1 + 28 * f0 - 84)
        g4.append(f3 - 6 * f2 + 21 * f1 - 56 * f0 + 126)
    return tuple(g3), tuple(g4)


def forced_shadow(family: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        five
        for five in independent_masks(5)
        if any(five & four == four for four in family)
    )


def shadow_gate(
    family: tuple[int, ...], shadow: tuple[int, ...], total_fives: int
) -> tuple[bool, bool, int, tuple[int, ...], tuple[int, ...]]:
    base_g3, base_g4 = base_link_g3_g4()
    additional = total_fives - len(shadow)
    if additional < 0:
        return False, False, additional, (), ()
    lower = []
    upper = []
    strict = True
    for vertex in range(N):
        incidence = sum(bool(face & (1 << vertex)) for face in family)
        previous = base_g3[vertex] - incidence
        forced = sum(bool(face & (1 << vertex)) for face in shadow)
        before = base_g4[vertex] + 6 * incidence - forced
        strict &= previous >= 1 and before >= 1
        lower.append(max(0, before - macaulay_raise(previous, 3)) if previous >= 0 else 0)
        upper.append(before - 1)
    aggregate = (
        strict
        and all(value >= 0 for value in upper)
        and sum(lower) <= 5 * additional <= sum(upper)
    )
    return strict, aggregate, additional, tuple(lower), tuple(upper)


def small_transversals(hyperedges: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        candidate
        for order in range(1, 5)
        for chosen in combinations(range(N), order)
        for candidate in (mask(chosen),)
        if all(candidate & edge for edge in hyperedges)
    )


def minimal_transversals(hyperedges: tuple[int, ...]) -> tuple[int, ...]:
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


def ridge_histogram(facets: tuple[int, ...]) -> dict[str, int]:
    degrees: Counter[int] = Counter()
    for complement in facets:
        facet = FULL_MASK ^ complement
        for removed in vertices(facet):
            degrees[facet ^ (1 << removed)] += 1
    return {str(key): value for key, value in sorted(Counter(degrees.values()).items())}


def hit_word(face: int, candidates: tuple[int, ...]) -> int:
    return sum(bool(face & candidate) << index for index, candidate in enumerate(candidates))


def independent_active_masks() -> tuple[int, ...]:
    return tuple(
        value
        for size in range(len(ACTIVE_VERTICES) + 1)
        for chosen in combinations(ACTIVE_VERTICES, size)
        for value in (mask(chosen),)
        if all(value & edge != edge for edge in EDGE_MASKS)
    )


@lru_cache(maxsize=1)
def active_automorphisms() -> tuple[dict[int, int], ...]:
    edge_set = {frozenset(edge) for edge in MISSING_EDGES}
    answers = []
    for image in permutations(ACTIVE_VERTICES):
        action = dict(zip(ACTIVE_VERTICES, image, strict=True))
        if {
            frozenset((action[left], action[right])) for left, right in MISSING_EDGES
        } == edge_set:
            answers.append(action)
    if len(answers) != 4:
        raise AssertionError("active automorphism count changed")
    return tuple(answers)


@lru_cache(maxsize=None)
def column_multisets(
    row_count: int, column_count: int, column_cap: int
) -> dict[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    allowed = tuple(
        column for column in range(1 << row_count) if column.bit_count() <= column_cap
    )
    grouped: defaultdict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
    for columns in combinations_with_replacement(allowed, column_count):
        sums = tuple(
            sum(bool(column & (1 << row)) for column in columns)
            for row in range(row_count)
        )
        grouped[sums].append(columns)
    return {key: tuple(value) for key, value in grouped.items()}


def rows_from_columns(
    active_rows: tuple[int, ...],
    groups: tuple[tuple[tuple[int, ...], tuple[int, ...]], ...],
) -> tuple[int, ...]:
    rows = list(active_rows)
    for group_vertices, columns in groups:
        for vertex, column in zip(group_vertices, columns, strict=True):
            for row in range(len(rows)):
                if column & (1 << row):
                    rows[row] |= 1 << vertex
    return tuple(rows)


def active_stabilizer(four: int) -> tuple[dict[int, int], ...]:
    result = []
    for action in active_automorphisms():
        moved = mask(action.get(vertex, vertex) for vertex in vertices(four))
        if moved == four:
            result.append(action)
    return tuple(result)


def full_family_orbit_key(family: tuple[int, ...]) -> list[int]:
    """Match the independent C++ full ``Aut(P4+K2) x S8 x S_a`` key."""

    keys = []
    for action in active_automorphisms():
        moved = tuple(
            mask(action.get(vertex, vertex) for vertex in vertices(face))
            for face in family
        )
        for order in permutations(range(len(family))):
            active_rows = tuple(moved[row] & ACTIVE_MASK for row in order)
            columns = tuple(sorted(
                sum(bool(moved[order[row]] & (1 << isolate)) << row
                    for row in range(len(order)))
                for isolate in ISOLATES
            ))
            keys.append((*active_rows, *columns))
    return list(min(keys))


def full_family_orbit_size(family: tuple[int, ...]) -> int:
    """Compute the literal orbit size by orbit--stabilizer under the full action."""

    target = tuple(sorted(family))
    face_vertices = tuple(vertices(face) for face in family)
    stabilizer = 0
    for action in active_automorphisms():
        for isolate_image in permutations(ISOLATES):
            moved = []
            for face in face_vertices:
                value = 0
                for vertex in face:
                    destination = (
                        action[vertex]
                        if vertex in action
                        else isolate_image[vertex - ISOLATES[0]]
                    )
                    value |= 1 << destination
                moved.append(value)
            stabilizer += tuple(sorted(moved)) == target
    group_order = len(active_automorphisms()) * 40320
    if stabilizer == 0 or group_order % stabilizer:
        raise AssertionError("invalid full-action stabilizer")
    return group_order // stabilizer


def family_orbit_key(four: int, additions: tuple[int, ...]) -> object:
    inside = tuple(vertex for vertex in ISOLATES if four & (1 << vertex))
    outside = tuple(vertex for vertex in ISOLATES if not four & (1 << vertex))
    keys = []
    for action in active_stabilizer(four):
        transformed = tuple(
            mask(action.get(vertex, vertex) for vertex in vertices(face))
            for face in additions
        )
        for ordered in permutations(transformed):
            active_parts = tuple(face & ACTIVE_MASK for face in ordered)
            inside_columns = tuple(sorted(
                sum(bool(ordered[row] & (1 << vertex)) << row for row in range(len(ordered)))
                for vertex in inside
            ))
            outside_columns = tuple(sorted(
                sum(bool(ordered[row] & (1 << vertex)) << row for row in range(len(ordered)))
                for vertex in outside
            ))
            keys.append((active_parts, inside_columns, outside_columns))
    return min(keys)


def audit_a1_representative(index: int, values: tuple[int, ...]) -> dict[str, object]:
    four = mask(values)
    family = (four,)
    shadow = forced_shadow(family)
    strict, aggregate, additional, lower, upper = shadow_gate(family, shadow, 12)
    if not strict or not aggregate or any(lower):
        raise AssertionError("a1 shadow gate changed")
    candidates = {
        five
        for five in independent_masks(5)
        if five not in shadow
        and five & four != four
        and all(not (five & (1 << vertex)) or upper[vertex] >= 1 for vertex in range(N))
    }
    inside = tuple(vertex for vertex in ISOLATES if four & (1 << vertex))
    outside = tuple(vertex for vertex in ISOLATES if not four & (1 << vertex))
    inside_cap = min((upper[vertex] for vertex in inside), default=additional)
    outside_cap = min((upper[vertex] for vertex in outside), default=additional)
    inside_patterns = column_multisets(additional, len(inside), inside_cap)
    outside_patterns = column_multisets(additional, len(outside), outside_cap)
    allowed_active = tuple(
        row
        for row in independent_active_masks()
        if all(not (row & (1 << vertex)) or upper[vertex] >= 1 for vertex in ACTIVE_VERTICES)
    )

    small = small_transversals((*EDGE_MASKS, four, *shadow))
    full_small_word = (1 << len(small)) - 1
    small_words = {five: hit_word(five, small) for five in candidates}
    covers = graph_cover_masks(5)
    fixed_word = hit_word(four, covers)
    five_words = {five: hit_word(five, covers) for five in candidates}

    active_count = 0
    column_paths = 0
    distinct_paths = 0
    no_small_paths = 0
    facet_hist: Counter[int] = Counter()
    target_families: set[tuple[int, ...]] = set()
    for active_rows in combinations_with_replacement(allowed_active, additional):
        if any(
            sum(bool(row & (1 << vertex)) for row in active_rows) > upper[vertex]
            for vertex in ACTIVE_VERTICES
        ):
            continue
        active_count += 1
        needed = tuple(5 - row.bit_count() for row in active_rows)
        for inside_sums, inside_rows in inside_patterns.items():
            outside_sums = tuple(needed[row] - inside_sums[row] for row in range(additional))
            if min(outside_sums, default=0) < 0:
                continue
            outside_rows = outside_patterns.get(outside_sums, ())
            for inside_columns in inside_rows:
                for outside_columns in outside_rows:
                    column_paths += 1
                    rows = rows_from_columns(
                        active_rows,
                        ((inside, inside_columns), (outside, outside_columns)),
                    )
                    if len(set(rows)) != additional or any(row not in candidates for row in rows):
                        continue
                    distinct_paths += 1
                    word = full_small_word
                    for row in rows:
                        word &= small_words[row]
                    if word:
                        continue
                    no_small_paths += 1
                    facet_word = fixed_word
                    for row in rows:
                        facet_word &= five_words[row]
                    facets = facet_word.bit_count()
                    facet_hist[facets] += 1
                    if facets == 112:
                        target_families.add(tuple(sorted(rows)))

    orbits: defaultdict[object, list[tuple[int, ...]]] = defaultdict(list)
    for additions in sorted(target_families):
        orbits[family_orbit_key(four, additions)].append(additions)
    representatives = tuple(min(rows) for rows in orbits.values())
    target_audits = []
    for additions in representatives:
        blocker = minimal_transversals((*EDGE_MASKS, four, *additions))
        histogram = Counter(value.bit_count() for value in blocker)
        ridge = ridge_histogram(tuple(value for value in blocker if value.bit_count() == 5))
        target_audits.append({
            "additions": [list(vertices(face)) for face in additions],
            "minimal_transversal_size_histogram": {str(k): v for k, v in sorted(histogram.items())},
            "ridge_histogram": ridge,
            "pure_ridge": histogram == {5: 112} and ridge == {"2": 504},
        })
    return {
        "representative_index": index,
        "A4": list(values),
        "shadow": len(shadow),
        "additional": additional,
        "candidates": len(candidates),
        "raw": comb(len(candidates), additional),
        "active_configurations": active_count,
        "column_paths": column_paths,
        "distinct_row_paths": distinct_paths,
        "no_small_paths": no_small_paths,
        "facet_histogram": {str(k): v for k, v in sorted(facet_hist.items())},
        "literal_targets": len(target_families),
        "target_orbits": len(orbits),
        "target_audits": target_audits,
    }


def at_least_masks(additions: tuple[int, ...]) -> tuple[int, ...]:
    if not additions:
        return ()
    if len(additions) == 1:
        return additions
    if len(additions) == 2:
        first, second = additions
        return first | second, first & second
    if len(additions) == 3:
        first, second, third = additions
        return (
            first | second | third,
            (first & second) | (first & third) | (second & third),
            first & second & third,
        )
    raise AssertionError("record-32 addition count escaped 0..3")


def audit_addition_representative(
    representative: tuple[tuple[int, ...], ...], total_fives: int, scheduled: int
) -> dict[str, object]:
    family = tuple(mask(face) for face in representative)
    shadow = forced_shadow(family)
    strict, aggregate, additional, lower, upper = shadow_gate(family, shadow, total_fives)
    if not strict or not aggregate or additional not in range(4):
        raise AssertionError("a2/a3 shadow gate changed")
    candidates = tuple(
        five
        for five in independent_masks(5)
        if five not in shadow
        and not any(five & four == four for four in family)
        and all(not (five & (1 << vertex)) or upper[vertex] >= 1 for vertex in range(N))
    )
    lower_at_least = tuple(
        mask(vertex for vertex in range(N) if lower[vertex] >= threshold)
        for threshold in range(1, additional + 1)
    )
    upper_below = tuple(
        mask(vertex for vertex in range(N) if upper[vertex] < threshold)
        for threshold in range(1, additional + 1)
    )
    small = small_transversals((*EDGE_MASKS, *family, *shadow))
    small_words = {five: hit_word(five, small) for five in candidates}
    covers = graph_cover_masks(5)
    all_words = {face: hit_word(face, covers) for face in (*family, *candidates)}
    fixed_word = (1 << len(covers)) - 1
    for four in family:
        fixed_word &= all_words[four]

    cap_feasible = 0
    no_small = 0
    facet_hist: Counter[int] = Counter()
    targets = []
    choices = ((),) if additional == 0 else combinations(candidates, additional)
    for additions in choices:
        additions = tuple(additions)
        present = at_least_masks(additions)
        if any(
            required & actual != required or forbidden & actual
            for required, forbidden, actual in zip(lower_at_least, upper_below, present, strict=True)
        ):
            continue
        cap_feasible += 1
        word = (1 << len(small)) - 1
        for face in additions:
            word &= small_words[face]
        if word:
            continue
        no_small += 1
        facet_word = fixed_word
        for face in additions:
            facet_word &= all_words[face]
        facets = facet_word.bit_count()
        facet_hist[facets] += 1
        if facets == scheduled:
            targets.append(additions)
    return {
        "missing_fours": [list(face) for face in representative],
        "shadow": len(shadow),
        "additional": additional,
        "candidates": len(candidates),
        "raw": comb(len(candidates), additional),
        "cap_feasible": cap_feasible,
        "no_small": no_small,
        "facet_histogram": {str(k): v for k, v in sorted(facet_hist.items())},
        "targets": [[list(vertices(face)) for face in row] for row in targets],
    }


def main() -> None:
    a1 = [audit_a1_representative(index, representative)
          for index, representative in enumerate(A1_REPRESENTATIVES)]
    expected_a1 = (
        (58, 1161, 1102, 9, 1, 1),
        (70, 43229, 42865, 271, 3, 3),
        (155, 2323448, 2278508, 12873, 0, 0),
        (410, 6068459, 5934396, 19836, 8, 6),
    )
    for row, expected in zip(a1, expected_a1):
        actual = tuple(row[key] for key in (
            "active_configurations", "column_paths", "distinct_row_paths",
            "no_small_paths", "literal_targets", "target_orbits"
        ))
        if actual != expected:
            raise AssertionError(f"a1 fixed-A4 transcript changed: {actual} != {expected}")
    ridge_targets = [
        {"A4": row["A4"], "additions": audit["additions"]}
        for row in a1
        for audit in row["target_audits"]
        if audit["pure_ridge"]
    ]
    expected_ridge_targets = [
        {"A4": [3, 4, 5, 6], "additions": [[2, 7, 8, 9, 10], [7, 8, 9, 10, 11]]},
        {"A4": [2, 3, 4, 5], "additions": [[6, 7, 8, 9, 10], [3, 4, 5, 6, 11], [7, 8, 9, 10, 13]]},
    ]
    if ridge_targets != expected_ridge_targets:
        raise AssertionError(
            f"derived a1 ridge-target identities changed: {ridge_targets!r}"
        )

    a2 = [audit_addition_representative(row, 18, 110) for row in A2_REPRESENTATIVES]
    expected_a2 = (
        (18, 0, 0, 1, 1, 0, 0),
        (17, 1, 1, 1, 1, 1, 1),
        (15, 3, 462, 16328620, 1794605, 840, 0),
        (16, 2, 56, 1540, 540, 0, 0),
        (16, 2, 56, 1540, 945, 0, 0),
        (16, 2, 56, 1540, 1540, 0, 0),
    )
    for row, expected in zip(a2, expected_a2):
        actual = (
            row["shadow"], row["additional"], row["candidates"], row["raw"],
            row["cap_feasible"], row["no_small"], len(row["targets"]),
        )
        if actual != expected:
            raise AssertionError(f"a2 transcript changed: {actual} != {expected}")
    selected_a2 = {
        "order_four_minimal_faces": a2[1]["missing_fours"],
        "order_five_minimal_face": a2[1]["targets"][0][0],
    }
    if selected_a2 != {
        "order_four_minimal_faces": [[2, 3, 4, 5], [3, 4, 5, 11]],
        "order_five_minimal_face": [6, 7, 8, 9, 10],
    }:
        raise AssertionError("selected a2 target identity changed")

    a3 = audit_addition_representative(A3_REPRESENTATIVE, 24, 108)
    actual_a3 = (
        a3["shadow"], a3["additional"], a3["candidates"], a3["raw"],
        a3["cap_feasible"], a3["no_small"], len(a3["targets"]),
    )
    if actual_a3 != (22, 2, 56, 1540, 540, 0, 0):
        raise AssertionError("a3 addition transcript changed")

    representative_families = {
        "a1": [tuple((mask(face),)) for face in A1_REPRESENTATIVES],
        "a2": [tuple(mask(face) for face in family) for family in A2_REPRESENTATIVES],
        "a3": [tuple(mask(face) for face in A3_REPRESENTATIVE)],
    }
    representative_keys = {
        branch: [full_family_orbit_key(family) for family in families]
        for branch, families in representative_families.items()
    }
    representative_sizes = {
        branch: [full_family_orbit_size(family) for family in families]
        for branch, families in representative_families.items()
    }
    expected_size_histograms = {
        "a1": Counter({28: 1, 70: 1, 112: 2}),
        "a2": Counter({56: 2, 112: 1, 420: 1, 672: 1, 840: 1}),
        "a3": Counter({56: 1}),
    }
    expected_viable_counts = {"a1": 322, "a2": 2156, "a3": 56}
    for branch in ("a1", "a2", "a3"):
        if len({tuple(key) for key in representative_keys[branch]}) != len(
            representative_keys[branch]
        ):
            raise AssertionError(f"{branch} representatives are not inequivalent")
        if Counter(representative_sizes[branch]) != expected_size_histograms[branch]:
            raise AssertionError(f"{branch} representative orbit sizes changed")
        if sum(representative_sizes[branch]) != expected_viable_counts[branch]:
            raise AssertionError(f"{branch} representatives do not exhaust viable families")

    print(json.dumps({
        "kind": "independent_h9_14_record32_additions_hostile_audit",
        "positive_representative_orbit_keys": representative_keys,
        "positive_representative_orbit_sizes": representative_sizes,
        "a1": {"representatives": a1, "ridge_targets": ridge_targets},
        "a2": {"representatives": a2, "selected_target": selected_a2},
        "a3": a3,
        "scope": {
            "lossless_fixed_A4_column_quotient": True,
            "lossless_a2_exhaustion": True,
            "lossless_a3_cap_zero_one_two_gate": True,
            "incidence_only": True,
            "compact_only": True,
            "finite_volume_noncompact_inference": False,
        },
        "status": "pass",
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
