#!/usr/bin/env python3
"""OOM-safe residual-support quotient for source 29 after the 012/123 cuts.

The first stage enumerates, for every normalized T=123 pair configuration,
the orbits of the still-new T-killer six-fibre supports under the literal
configuration stabilizer.  It is an exact finite reduction, not a survivor
or UNSAT claim.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from math import comb, factorial
from pathlib import Path
import resource
import sys
from typing import Iterable, Mapping, Sequence

from scripts.independent_h9_14_record0_a2_source29_triple123_hostile_audit import (
    EDGES,
    EDGE_INDEX,
    FREE_FIBRE,
    Q,
    REMAINDER,
    T,
    VERTICES,
    _allowed,
    _block_group,
)


ROOT = Path(__file__).resolve().parents[1]
RSS_LIMIT_BYTES = 3 * (1 << 28)  # 0.75 GiB
FIBRES = tuple(combinations(Q, 3))
FIBRE_INDEX = {frozenset(row): index for index, row in enumerate(FIBRES)}
FIXED_FIBRE_INDEX = FIBRE_INDEX[FREE_FIBRE]
EXPECTED_MATHEMATICAL_DIGEST = "1273598d67613ffbf02afd367576e002ee51f8f3b3742a37ea8556ee9bdbcb8a"


def _mask(vertices: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in vertices)


def _decode(value: int) -> frozenset[int]:
    return frozenset(vertex for vertex in VERTICES if value >> vertex & 1)


def _base_g(face: frozenset[int]) -> tuple[int, ...]:
    link_dimension = 9 - len(face)
    outside = tuple(vertex for vertex in VERTICES if vertex not in face)
    h = []
    for index in range(link_dimension // 2 + 1):
        entry = 0
        for order in range(index + 1):
            count = 1 if order == 0 else sum(
                _allowed(face | frozenset(extension))
                for extension in combinations(outside, order)
            )
            entry += (
                (-1) ** (index - order)
                * comb(link_dimension - order, link_dimension - index)
                * count
            )
        h.append(entry)
    return tuple(h[index] - h[index - 1] for index in range(1, len(h)))


def _degree_caps() -> dict[frozenset[int], int]:
    result = {}
    for size, entry in ((1, 3), (2, 2), (3, 1)):
        for row in combinations(VERTICES, size):
            face = frozenset(row)
            if _allowed(face):
                result[face] = _base_g(face)[entry] - 1
    return result


DEGREE_CAPS = _degree_caps()


def _rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def _rss_guard() -> None:
    if _rss_bytes() > RSS_LIMIT_BYTES:
        raise MemoryError("triple-123 support quotient exceeded 0.75 GiB")


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _group(blocks: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    """Return the full group by listing the image of every vertex in Q."""

    return tuple(
        tuple(action.get(vertex, vertex) for vertex in Q)
        for action in _block_group(*(tuple(block) for block in blocks))
    )


def _adjacent_generators(blocks: Sequence[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    generators = []
    for block in blocks:
        for index in range(len(block) - 1):
            left, right = block[index], block[index + 1]
            generators.append(
                tuple(
                    right if vertex == left else left if vertex == right else vertex
                    for vertex in Q
                )
            )
    return tuple(generators)


def _induced_actions(
    vertex_group: Sequence[Sequence[int]], rows: Sequence[Sequence[int]]
) -> tuple[tuple[int, ...], ...]:
    index = {frozenset(row): i for i, row in enumerate(rows)}
    return tuple(
        tuple(
            index[frozenset(action[Q.index(vertex)] for vertex in row)]
            for row in rows
        )
        for action in vertex_group
    )


def _configuration_orbits(
    degree: int, edge_actions: Sequence[Sequence[int]]
) -> tuple[tuple[int, ...], ...]:
    if degree == 0:
        return ((),)
    remaining = set(combinations(range(len(EDGES)), degree))
    representatives: list[tuple[int, ...]] = []
    while remaining:
        representative = min(remaining)
        orbit = {
            tuple(sorted(action[index] for index in representative))
            for action in edge_actions
        }
        remaining.difference_update(orbit)
        representatives.append(representative)
    return tuple(representatives)


def _stabilizer_indices(
    configuration: tuple[int, ...], edge_actions: Sequence[Sequence[int]]
) -> tuple[int, ...]:
    return tuple(
        index
        for index, action in enumerate(edge_actions)
        if tuple(sorted(action[item] for item in configuration)) == configuration
    )


def _support_orbits(
    universe: tuple[int, ...],
    fixed: tuple[int, ...],
    maximum_size: int,
    fibre_actions: Sequence[Sequence[int]],
    generator_actions: Sequence[Sequence[int]] = (),
) -> tuple[tuple[tuple[int, ...], ...], ...]:
    if any(item not in universe for item in fixed):
        raise AssertionError("fixed support left the new-fibre universe")
    fixed_set = set(fixed)
    if any(
        {action[item] for item in fixed} != fixed_set
        for action in fibre_actions
    ):
        raise AssertionError("configuration stabilizer does not preserve fixed support")
    available = tuple(item for item in universe if item not in fixed_set)
    result = []
    for size in range(len(fixed), maximum_size + 1):
        choose = size - len(fixed)
        remaining = {
            tuple(sorted((*fixed, *row)))
            for row in combinations(available, choose)
        }
        _rss_guard()
        representatives = []
        while remaining:
            seed = next(iter(remaining))
            if generator_actions:
                orbit = {seed}
                stack = [seed]
                while stack:
                    current = stack.pop()
                    for action in generator_actions:
                        moved = tuple(sorted(action[item] for item in current))
                        if moved not in orbit:
                            orbit.add(moved)
                            stack.append(moved)
            else:
                orbit = {
                    tuple(sorted(action[item] for item in seed))
                    for action in fibre_actions
                }
            remaining.difference_update(orbit)
            representatives.append(min(orbit))
        result.append(tuple(sorted(representatives)))
        del remaining
        _rss_guard()
    return tuple(result)


def _burnside_support_counts(
    universe: tuple[int, ...],
    fixed: tuple[int, ...],
    maximum_size: int,
    fibre_actions: Sequence[Sequence[int]],
) -> tuple[int, ...]:
    """Independent cycle-index count for the same bounded support orbits."""

    fixed_set = set(fixed)
    available = set(universe) - fixed_set
    maximum_extra = maximum_size - len(fixed)
    sums = [0] * (maximum_extra + 1)
    for action in fibre_actions:
        if {action[item] for item in available} != available:
            raise AssertionError("support universe is not stabilizer-invariant")
        seen: set[int] = set()
        cycle_lengths = []
        for seed in available:
            if seed in seen:
                continue
            current = seed
            length = 0
            while current not in seen:
                seen.add(current)
                length += 1
                current = action[current]
            cycle_lengths.append(length)
        coefficients = [1] + [0] * maximum_extra
        for length in cycle_lengths:
            for degree in range(maximum_extra, length - 1, -1):
                coefficients[degree] += coefficients[degree - length]
        sums = [left + right for left, right in zip(sums, coefficients, strict=True)]
    if any(value % len(fibre_actions) for value in sums):
        raise AssertionError("support Burnside average was not integral")
    return tuple(value // len(fibre_actions) for value in sums)


def _minimal_killer_gate(
    pattern: tuple[int, ...],
    configuration: tuple[int, ...],
    support: tuple[int, ...],
) -> bool:
    """Necessary monotone degree/shadow gate using one killer per new fibre.

    If no choice of one killer on each not-already-killed support fibre obeys
    every vertex/edge/triple incidence cap, adding further killers cannot
    repair the violation.  The fixed minimal shadow must also not already
    exceed any global target.  Passing is deliberately not called survival.
    """

    fixed = set(pattern)
    fixed.update(_mask(T | frozenset(EDGES[index])) for index in configuration)
    groups = []
    for fibre_index in support:
        six = T | frozenset(FIBRES[fibre_index])
        if any(_decode(value) < six for value in fixed):
            continue
        groups.append(
            tuple(
                _mask((T - {vertex}) | frozenset(FIBRES[fibre_index]))
                for vertex in sorted(T)
            )
        )

    counts: Counter[frozenset[int]] = Counter()

    def install(value: int, delta: int) -> bool:
        vertices = _decode(value)
        okay = True
        for size in (1, 2, 3):
            for row in combinations(vertices, size):
                face = frozenset(row)
                counts[face] += delta
                if delta > 0 and counts[face] > DEGREE_CAPS[face]:
                    okay = False
        return okay

    if any(not install(value, 1) for value in sorted(fixed)):
        return False

    def shadow_upper_ok(selected: set[int]) -> bool:
        for order, target in ((6, 147), (7, 420), (8, 629)):
            missing = set()
            for value in selected:
                vertices = _decode(value)
                outside = tuple(vertex for vertex in VERTICES if vertex not in vertices)
                for extension in combinations(outside, order - 5):
                    face = vertices | frozenset(extension)
                    if _allowed(face):
                        missing.add(face)
            if len(missing) > target:
                return False
        return True

    selected = set(fixed)

    def search(index: int) -> bool:
        if index == len(groups):
            return shadow_upper_ok(selected)
        for value in groups[index]:
            okay = install(value, 1)
            if okay:
                selected.add(value)
                if search(index + 1):
                    return True
                selected.remove(value)
            install(value, -1)
        return False

    return search(0)


def _row(
    family: str,
    configuration: tuple[int, ...],
    edge_actions: Sequence[Sequence[int]],
    all_fibre_actions: Sequence[Sequence[int]],
    base_generator_fibre_actions: Sequence[Sequence[int]],
    *,
    fixed_230: bool,
) -> dict[str, object] | None:
    degree = len(configuration)
    pairs = tuple(frozenset(EDGES[index]) for index in configuration)
    covered = {
        index
        for index, fibre in enumerate(FIBRES)
        if any(pair <= frozenset(fibre) for pair in pairs)
    }
    coverage = len(covered)
    expected_coverage = (0, 7, None)[degree]
    if degree == 2:
        expected_coverage = 13 if pairs[0] & pairs[1] else 14
    if coverage != expected_coverage:
        raise AssertionError("pair base coverage changed")
    maximum = {0: 4, 1: 2, 2: 1 if coverage == 13 else 0}[degree]
    universe = tuple(index for index in range(len(FIBRES)) if index not in covered)
    fixed = (
        (FIXED_FIBRE_INDEX,)
        if fixed_230 and FIXED_FIBRE_INDEX not in covered
        else ()
    )
    if len(fixed) > maximum:
        return None
    stabilizer = _stabilizer_indices(configuration, edge_actions)
    fibre_actions = tuple(all_fibre_actions[index] for index in stabilizer)
    # The only million-row raw layer is d=0, whose stabilizer is the complete
    # base group generated by the adjacent block transpositions.  Generator
    # BFS visits each raw support a bounded number of times.  All other rows
    # have maximum support size at most two and can use literal full images.
    generators = base_generator_fibre_actions if degree == 0 else ()
    support_orbits = _support_orbits(
        universe,
        fixed,
        maximum,
        fibre_actions,
        generators,
    )
    burnside_counts = _burnside_support_counts(
        universe, fixed, maximum, fibre_actions
    )
    if burnside_counts != tuple(map(len, support_orbits)):
        raise AssertionError("literal support orbits disagree with Burnside")
    size_counts = {
        str(len(fixed) + offset): len(rows)
        for offset, rows in enumerate(support_orbits)
    }
    transcript = [
        [list(FIBRES[index]) for index in support]
        for rows in support_orbits
        for support in rows
    ]
    if family == "zero":
        patterns = ((),)
    elif fixed_230:
        patterns = ((229, 230), (227, 229, 230))
    else:
        patterns = ((227,), (229,), (227, 229))
    minimal_gate_rows = []
    flat_supports = tuple(support for rows in support_orbits for support in rows)
    for pattern in patterns:
        passing = sum(
            _minimal_killer_gate(pattern, configuration, support)
            for support in flat_supports
        )
        minimal_gate_rows.append(
            {
                "triple012_pattern": list(pattern),
                "necessary_minimal_killer_gate_pass": passing,
                "necessary_minimal_killer_gate_reject": len(flat_supports) - passing,
            }
        )
    result: dict[str, object] = {
        "family": family,
        "configuration_pairs": [list(EDGES[index]) for index in configuration],
        "degree": degree,
        "pair_intersection_size": (
            len(pairs[0] & pairs[1]) if degree == 2 else None
        ),
        "base_covered_fibre_count": coverage,
        "outside_fibre_budget": maximum,
        "fixed_230": fixed_230,
        "fixed_new_fibre_567": bool(fixed),
        "configuration_stabilizer_order": len(stabilizer),
        "support_orbit_count_by_total_size": size_counts,
        "support_orbit_count": sum(size_counts.values()),
        "support_representative_digest": _digest(transcript),
        "minimal_killer_necessary_gate_by_pattern": minimal_gate_rows,
    }
    return result


def quotient() -> dict[str, object]:
    rows: list[dict[str, object]] = []
    group_specs = (
        ("zero", (tuple(range(5, 13)),), (False,)),
        ("active", (tuple(FREE_FIBRE), REMAINDER), (False, True)),
    )
    for family, blocks, fixed_options in group_specs:
        vertex_group = _group(blocks)
        if len(vertex_group) != (factorial(8) if family == "zero" else factorial(3) * factorial(5)):
            raise AssertionError("residual group order changed")
        edge_actions = _induced_actions(vertex_group, EDGES)
        fibre_actions = _induced_actions(vertex_group, FIBRES)
        generator_fibre_actions = _induced_actions(
            _adjacent_generators(blocks), FIBRES
        )
        expected = {"zero": (1, 2, 5), "active": (1, 5, 24)}[family]
        for degree in range(3):
            configurations = _configuration_orbits(degree, edge_actions)
            if len(configurations) != expected[degree]:
                raise AssertionError("pair-configuration orbit count changed")
            for configuration in configurations:
                for fixed_230 in fixed_options:
                    row = _row(
                        family,
                        configuration,
                        edge_actions,
                        fibre_actions,
                        generator_fibre_actions,
                        fixed_230=fixed_230,
                    )
                    if row is not None:
                        rows.append(row)
        del vertex_group, edge_actions, fibre_actions, generator_fibre_actions
        _rss_guard()

    summary = Counter()
    for row in rows:
        summary[(str(row["family"]), bool(row["fixed_230"]), int(row["degree"]))] += int(
            row["support_orbit_count"]
        )
    result = {
        "kind": "source29_triple123_residual_support_quotient",
        "claim_scope": (
            "Exact quotient of bounded new T=123 killer-fibre supports only; "
            "rows are not completed strict-link survivors."
        ),
        "rows": rows,
        "summary": [
            {
                "family": key[0],
                "fixed_230": key[1],
                "degree": key[2],
                "support_orbit_count": value,
            }
            for key, value in sorted(summary.items())
        ],
        "row_count": len(rows),
        "support_orbit_count": sum(int(row["support_orbit_count"]) for row in rows),
        "configuration_branch_instance_count": sum(
            1 if row["family"] == "zero" else 2 if row["fixed_230"] else 3
            for row in rows
        ),
        "support_branch_instance_count": sum(
            int(row["support_orbit_count"])
            * (1 if row["family"] == "zero" else 2 if row["fixed_230"] else 3)
            for row in rows
        ),
        "minimal_killer_necessary_gate_pass_count": sum(
            int(gate["necessary_minimal_killer_gate_pass"])
            for row in rows
            for gate in row["minimal_killer_necessary_gate_by_pattern"]
        ),
        "minimal_killer_necessary_gate_reject_count": sum(
            int(gate["necessary_minimal_killer_gate_reject"])
            for row in rows
            for gate in row["minimal_killer_necessary_gate_by_pattern"]
        ),
        "source29_status": "OPEN",
        "peak_rss_bytes": _rss_bytes(),
    }
    result["mathematical_digest"] = _digest(
        {key: value for key, value in result.items() if key != "peak_rss_bytes"}
    )
    if result["mathematical_digest"] != EXPECTED_MATHEMATICAL_DIGEST:
        raise AssertionError("residual-support quotient transcript changed")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()
    document = quotient()
    if args.summary_only:
        document = {key: value for key, value in document.items() if key != "rows"}
    print(json.dumps(document, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
