#!/usr/bin/env python3
"""Read-only exact scout for the record-35 whole-shadow completion problem.

This file is intentionally exploratory.  It prints exact counts that are later
frozen by the production certificate and is not itself a theorem artifact.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations
from math import comb
import json
import sys

from coxeter_polytopes.combinatorics.facet_completion import audit_global_facet_selection
from coxeter_polytopes.combinatorics.h9_14_coxeter_labels import (
    is_finite_coxeter_labeling,
    lanner_placement_codes,
    pair,
)
from coxeter_polytopes.combinatorics.h9_14_record46 import (
    _a4_interval_feasible,
    _base_link_g3_g4,
)
from coxeter_polytopes.combinatorics.h9_14_record92 import (
    graph_vertex_cover_masks,
    independent_masks,
)
from coxeter_polytopes.combinatorics.local_link import analyze_local_link_graph
from coxeter_polytopes.combinatorics.missing_edges import decode_graph6, macaulay_raise


GRAPH6 = "M?????????O?O?w??"
EDGES = ((0, 11), (1, 12), (0, 13), (1, 13), (2, 13))
EDGE_MASKS = tuple(sum(1 << v for v in edge) for edge in EDGES)
ACTIVE = (0, 1, 2, 11, 12, 13)
ISOLATES = tuple(range(3, 11))
FULL = (1 << 14) - 1


def mask(values):
    return sum(1 << value for value in values)


def vertices(value):
    return tuple(v for v in range(14) if value >> v & 1)


def minimal_graph_covers():
    answer = []
    for value in range(1, 1 << 14):
        if not all(value & edge for edge in EDGE_MASKS):
            continue
        if all(
            not all((value ^ (1 << v)) & edge for edge in EDGE_MASKS)
            for v in vertices(value)
        ):
            answer.append(value)
    return tuple(answer)


def forced_shadow(family, fives):
    return tuple(five for five in fives if any(five & four == four for four in family))


def family_interval_feasible(family, local):
    for figure in local.vertex_figures:
        incidence = sum(bool(face & (1 << figure.vertex)) for face in family)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.edge_figures:
        subset = mask(figure.edge)
        incidence = sum(face & subset == subset for face in family)
        interval = figure.missing_four_incidence_interval
        if interval is None or not interval[0] <= incidence <= interval[1]:
            return False
    for figure in local.triple_figures:
        subset = mask(figure.triple)
        incidence = sum(face & subset == subset for face in family)
        if incidence > figure.missing_four_incidence_max:
            return False
    return True


def shadow_gate(family, shadow, base_g3, base_g4, total_fives):
    additional = total_fives - len(shadow)
    if additional < 0:
        return False, False, additional, (), ()
    lower, upper = [], []
    strict = True
    for vertex in range(14):
        a4_incidence = sum(bool(face & (1 << vertex)) for face in family)
        previous = base_g3[vertex] - a4_incidence
        forced_incidence = sum(bool(face & (1 << vertex)) for face in shadow)
        before = base_g4[vertex] + 6 * a4_incidence - forced_incidence
        strict &= previous >= 1 and before >= 1
        lower.append(max(0, before - macaulay_raise(previous, 3)) if previous >= 0 else 0)
        upper.append(before - 1)
    aggregate = strict and all(x >= 0 for x in upper) and sum(lower) <= 5 * additional <= sum(upper)
    return strict, aggregate, additional, tuple(lower), tuple(upper)


def has_small_transversal(hyperedges, maximum=4):
    return next(
        (
            order
            for order in range(1, maximum + 1)
            if any(
                all(mask(chosen) & edge for edge in hyperedges)
                for chosen in combinations(range(14), order)
            )
        ),
        None,
    )


def has_small_transversal_from_graph_covers(fives, graph_covers):
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


def active_automorphisms():
    edge_set = {frozenset(edge) for edge in EDGES}
    answer = []
    for image_order in permutations(ACTIVE):
        action = dict(zip(ACTIVE, image_order, strict=True))
        if {frozenset((action[a], action[b])) for a, b in EDGES} == edge_set:
            answer.append(action)
    return tuple(answer)


AUTS = active_automorphisms()


def family_orbit_key(family):
    active_mask = mask(ACTIVE)
    keys = []
    for action in AUTS:
        transformed = tuple(
            mask(action.get(v, v) for v in vertices(face)) for face in family
        )
        for order in permutations(range(len(family))):
            active_rows = tuple(transformed[index] & active_mask for index in order)
            columns = tuple(
                sorted(
                    sum(bool(transformed[order[row]] & (1 << isolate)) << row for row in range(len(family)))
                    for isolate in ISOLATES
                )
            )
            keys.append((active_rows, columns))
    return min(keys)


def enumerate_additions(family, total_fives):
    fives = independent_masks(GRAPH6, 5)
    base_g3, base_g4 = _base_link_g3_g4(GRAPH6)
    shadow = forced_shadow(family, fives)
    strict, aggregate, additional, lower, upper = shadow_gate(
        family, shadow, base_g3, base_g4, total_fives
    )
    candidates = tuple(
        five
        for five in fives
        if five not in shadow
        and not any(five & four == four for four in family)
        and all(not (five & (1 << vertex)) or upper[vertex] >= 1 for vertex in range(14))
    )
    raw = comb(len(candidates), additional)
    cap = []
    no_small = []
    target = []
    covers5 = graph_vertex_cover_masks(GRAPH6, 5)
    base_clutter = (*EDGE_MASKS, *family, *shadow)
    small_base_transversals = tuple(
        mask(chosen)
        for order in range(1, 5)
        for chosen in combinations(range(14), order)
        if all(mask(chosen) & edge for edge in base_clutter)
    )
    for additions in combinations(candidates, additional):
        incidence = tuple(sum(bool(face & (1 << v)) for face in additions) for v in range(14))
        if not all(lower[v] <= incidence[v] <= upper[v] for v in range(14)):
            continue
        cap.append(additions)
        if any(all(transversal & face for face in additions) for transversal in small_base_transversals):
            continue
        no_small.append(additions)
        facet_count = sum(all(cover & face for face in (*family, *shadow, *additions)) for cover in covers5)
        if facet_count == 100 if len(family) == 1 else facet_count == 98:
            target.append(additions)
    return {
        "strict": strict,
        "aggregate": aggregate,
        "shadow": shadow,
        "additional": additional,
        "candidate_count": len(candidates),
        "raw": raw,
        "cap": tuple(cap),
        "no_small": tuple(no_small),
        "target": tuple(target),
    }


def positive_branches():
    local = analyze_local_link_graph(GRAPH6)
    fours = tuple(face for face in independent_masks(GRAPH6, 4) if _a4_interval_feasible(face, local))
    fives = independent_masks(GRAPH6, 5)
    base_g3, base_g4 = _base_link_g3_g4(GRAPH6)
    gates1 = []
    viable1 = []
    for four in fours:
        shadow = forced_shadow((four,), fives)
        gate = shadow_gate((four,), shadow, base_g3, base_g4, 11)
        gates1.append((gate[0], gate[1], len(shadow), gate[2]))
        if gate[0] and gate[1]:
            viable1.append(four)
    orbits1 = defaultdict(list)
    for four in viable1:
        orbits1[family_orbit_key((four,))].append(four)

    viable2 = []
    hist2 = Counter()
    local2 = 0
    for family in combinations(fours, 2):
        if not family_interval_feasible(family, local):
            continue
        local2 += 1
        shadow = forced_shadow(family, fives)
        gate = shadow_gate(family, shadow, base_g3, base_g4, 17)
        hist2[len(shadow)] += 1
        if gate[0] and gate[1]:
            viable2.append(family)
    orbits2 = defaultdict(list)
    for family in viable2:
        orbits2[family_orbit_key(family)].append(family)

    print("AUT", len(AUTS), "FOURS", len(fours))
    print("A1_GATES", tuple(sorted(Counter(gates1).items())))
    print("A1_VIABLE", len(viable1), "ORBITS", [len(v) for v in orbits1.values()])
    for orbit in orbits1.values():
        representative = min(orbit)
        audit = enumerate_additions((representative,), 11)
        print("A1_REP", vertices(representative), "SHADOW", len(audit["shadow"]), "ADD", audit["additional"], "CAND", audit["candidate_count"], "RAW", audit["raw"], "CAP", len(audit["cap"]), "NOSMALL", len(audit["no_small"]), "TARGET", len(audit["target"]))
        for additions in audit["target"]:
            print("A1_TARGET", [vertices(value) for value in additions])
    print("A2_RAW", comb(len(fours), 2), "LOCAL", local2, "HIST", tuple(sorted(hist2.items())), "VIABLE", len(viable2), "ORBITS", [len(v) for v in orbits2.values()])
    for orbit in orbits2.values():
        representative = min(orbit)
        audit = enumerate_additions(representative, 17)
        print("A2_REP", [vertices(value) for value in representative], "SHADOW", len(audit["shadow"]), "ADD", audit["additional"], "CAND", audit["candidate_count"], "RAW", audit["raw"], "CAP", len(audit["cap"]), "NOSMALL", len(audit["no_small"]), "TARGET", len(audit["target"]))


def distribution_iterator(allowed, row_sums):
    counts = [0] * len(allowed)
    def visit(position, left, sums):
        if position == len(allowed):
            if left == 0 and not any(sums):
                yield tuple(counts)
            return
        column = allowed[position]
        maximum = left
        for row in range(5):
            if column & (1 << row):
                maximum = min(maximum, sums[row])
        for multiplicity in range(maximum + 1):
            next_left = left - multiplicity
            next_sums = tuple(sums[row] - multiplicity * bool(column & (1 << row)) for row in range(5))
            if all(0 <= value <= next_left for value in next_sums):
                counts[position] = multiplicity
                yield from visit(position + 1, next_left, next_sums)
        counts[position] = 0
    yield from visit(0, len(ISOLATES), row_sums)


@lru_cache(maxsize=None)
def distributions(allowed, row_sums):
    return tuple(distribution_iterator(allowed, row_sums))


def a0_census():
    covers = minimal_graph_covers()
    covers5 = graph_vertex_cover_masks(GRAPH6, 5)
    active_edges = tuple(edge for edge in EDGE_MASKS if edge & mask(ACTIVE) == edge)
    _, base_g4 = _base_link_g3_g4(GRAPH6)
    caps = tuple(value - 1 for value in base_g4)
    independent_active = tuple(
        value
        for order in range(len(ACTIVE) + 1)
        for chosen in combinations(ACTIVE, order)
        for value in (mask(chosen),)
        if all(value & edge != edge for edge in active_edges)
        and all(int(bool(value & (1 << v))) <= caps[v] for v in ACTIVE)
    )
    configs = []
    for config in combinations_with_replacement(independent_active, 5):
        if any(sum(bool(row & (1 << v)) for row in config) > caps[v] for v in ACTIVE):
            continue
        missed = tuple(sum(1 << i for i, row in enumerate(config) if not row & cover) for cover in covers)
        if any(
            all(not (pattern & (1 << i)) or config[i] & (1 << v) for i in range(5))
            for pattern in missed for v in ACTIVE
        ):
            continue
        configs.append((config, missed))
    states = []
    repeated = small = 0
    targets = []
    for config, missed in configs:
        allowed = tuple(column for column in range(1 << 5) if column.bit_count() <= 4 and all(column & rows != rows for rows in missed))
        sums = tuple(5 - row.bit_count() for row in config)
        for mult in distributions(allowed, sums):
            columns = tuple(column for column, count in zip(allowed, mult, strict=True) for _ in range(count))
            fives = list(config)
            for isolate, column in zip(ISOLATES, columns, strict=True):
                for row in range(5):
                    if column & (1 << row):
                        fives[row] |= 1 << isolate
            family = tuple(sorted(fives))
            if len(set(family)) != 5:
                repeated += 1
                continue
            if has_small_transversal_from_graph_covers(family, covers):
                small += 1
                continue
            facet_count = sum(all(cover & face for face in family) for cover in covers5)
            states.append((family, facet_count))
            if facet_count == 102:
                targets.append(family)
    exact_targets = tuple(sorted(set(targets)))
    orbit_groups = defaultdict(list)
    for family in exact_targets:
        orbit_groups[family_orbit_key(family)].append(family)
    print("A0 ACTIVE", len(independent_active), "CONFIG", len(configs), "CACHE", distributions.cache_info().currsize)
    print("A0 PATH", len(states), "REPEAT", repeated, "SMALL", small, "HIST", tuple(sorted(Counter(count for _, count in states).items())))
    print("A0 TARGET PATH", len(targets), "EXACT", len(exact_targets), "ORBITS", [len(v) for v in orbit_groups.values()])
    for index, family in enumerate(exact_targets):
        print("A0_TARGET", index, [vertices(value) for value in family])
    return exact_targets


def minimal_transversals(hyperedges):
    answer = []
    for candidate in range(1, FULL + 1):
        if not all(candidate & edge for edge in hyperedges):
            continue
        if any(
            all((candidate ^ (1 << vertex)) & edge for edge in hyperedges)
            for vertex in vertices(candidate)
        ):
            continue
        answer.append(candidate)
    return tuple(answer)


@lru_cache(maxsize=None)
def placement_options(support):
    support = tuple(sorted(support))
    edges = tuple(combinations(support, 2))
    return tuple(tuple(zip(edges, code, strict=True)) for code in lanner_placement_codes(len(support)))


def rank_nine_vertices(nonfaces):
    frozen = tuple(frozenset(face) for face in nonfaces)
    return tuple(
        chosen
        for chosen in combinations(range(14), 9)
        if not any(face <= set(chosen) for face in frozen)
    )


@lru_cache(maxsize=262_144)
def finite_code(order, code):
    labels = dict(zip(combinations(range(order), 2), code, strict=True))
    return is_finite_coxeter_labeling(range(order), labels)


def restriction_finite(row, labels):
    return finite_code(
        len(row), tuple(labels[pair(a, b)] for a, b in combinations(row, 2))
    )


def audit_family_labels(index, family):
    supports = tuple(vertices(face) for face in family)
    actual = rank_nine_vertices((*EDGES, *supports))
    remaining = set(range(5))
    union = set()
    states = ((),)
    transcript = []
    while states and remaining:
        support_index = min(remaining) if not union else max(
            remaining,
            key=lambda i: (len(set(supports[i]) & union), -len(set(supports[i]) - union), -i),
        )
        support = supports[support_index]
        union.update(support)
        cluster = tuple(sorted(union))
        restrictions = tuple(
            sorted(
                {tuple(f for f in vertex if f in union) for vertex in actual},
                key=lambda row: (len(row), row),
            )
        )
        passing = set()
        compatible = 0
        for state in states:
            old = dict(state)
            for option in placement_options(support):
                if not all(old.get(edge, label) == label for edge, label in option):
                    continue
                compatible += 1
                merged = old | dict(option)
                labels = {pair(a, b): 2 for a, b in combinations(cluster, 2)}
                labels.update(merged)
                if all(restriction_finite(row, labels) for row in restrictions):
                    passing.add(tuple(sorted(merged.items())))
        states = tuple(sorted(passing))
        transcript.append((support_index, len(states), compatible, tuple(sorted(Counter(map(len, restrictions)).items()))))
        remaining.remove(support_index)
    print("A0_LABEL", index, "VERT", len(actual), "STAGES", transcript)
    return not states


def audit_a0_targets(targets):
    document = json.load(open("data/enumerations/h9_14_nerve_support.json", encoding="utf-8"))
    record = [row for row in document["records"] if row["survives_second_layer"]][35]
    schedules = record["dehn_sommerville_schedule"]
    support = record["ridge_support_core_masks"]
    topology = []
    for index, family in enumerate(targets):
        transversals = minimal_transversals((*EDGE_MASKS, *family))
        histogram = tuple(sorted(Counter(value.bit_count() for value in transversals).items()))
        audit = audit_global_facet_selection(GRAPH6, support, schedules, tuple(value for value in transversals if value.bit_count() == 5))
        row = audit.to_record()
        print("A0_GLOBAL", index, histogram, row["ridge_degree_histogram"], row["passes_connected_homology_audit"])
        if row["passes_connected_homology_audit"]:
            topology.append(index)
    groups = defaultdict(list)
    for index in topology:
        groups[family_orbit_key(targets[index])].append(index)
    print("A0_TOPOLOGY", topology, "ORBITS", list(groups.values()))
    for index, family in enumerate(targets):
        assert audit_family_labels(index, family)


if __name__ == "__main__":
    order, decoded = decode_graph6(GRAPH6)
    assert order == 14 and decoded == EDGES
    if "--a0-only" not in sys.argv:
        positive_branches()
    audit_a0_targets(a0_census())
