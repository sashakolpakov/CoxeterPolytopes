#!/usr/bin/env python3
"""Certify the source-29 triple-123 fibre-budget refinement."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Mapping, Sequence

from coxeter_polytopes.combinatorics.missing_edges import macaulay_raise
from scripts.certify_h9_14_record0_a2_source29_triple012_cut import (
    NEAR_SEED,
    _analytic_cut as triple012_cut,
)
from scripts.explore_h9_14_record0_a2_dual_facets import mask
from scripts.explore_h9_14_record0_a2_open_primal_strong import (
    _base_link_g_vector,
    source_problem,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import vertices


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record0_a2_source29_triple123_cut.json"
SOURCE = 29
T_VERTICES = (1, 2, 3)
T = mask(T_VERTICES)
OUTSIDE = tuple(range(4, 13))
FREE_FIBRE = (5, 6, 7)
FIXED_FIBRE_MASKS = (227, 229, 230)
EXPECTED_MATHEMATICAL_DIGEST = "c181454ab3136546b4727196d618f988c9c93c929e59f6327e59de6b142fa3dc"


def _digest(value: object) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _move_pair(pair: tuple[int, int], action: Mapping[int, int]) -> tuple[int, int]:
    return tuple(sorted(action.get(vertex, vertex) for vertex in pair))  # type: ignore[return-value]


def _configuration_orbits(
    degree: int, generators: Sequence[Mapping[int, int]]
) -> tuple[tuple[tuple[tuple[int, int], ...], ...], ...]:
    pairs = tuple(combinations(OUTSIDE, 2))
    configurations = tuple(combinations(pairs, degree))
    remaining = set(configurations)
    orbits = []
    while remaining:
        representative = min(remaining)
        orbit = {representative}
        stack = [representative]
        while stack:
            row = stack.pop()
            for action in generators:
                moved = tuple(sorted(_move_pair(pair, action) for pair in row))
                if moved not in orbit:
                    orbit.add(moved)
                    stack.append(moved)
        remaining -= orbit
        orbits.append(tuple(sorted(orbit)))
    return tuple(orbits)


def _generators(blocks: Sequence[Sequence[int]]) -> tuple[dict[int, int], ...]:
    result = []
    for block in blocks:
        for index in range(len(block) - 1):
            result.append({block[index]: block[index + 1], block[index + 1]: block[index]})
    return tuple(result)


def _strict_table(base_g: tuple[int, ...]) -> list[dict[str, int | str]]:
    rows = []
    for degree in range(base_g[1]):
        g2 = base_g[1] - degree
        constant = base_g[2] + 5 * degree
        lower = constant - macaulay_raise(g2, 2)
        upper = constant
        if degree == 0:
            base_coverage = 0
            kind = "none"
        elif degree == 1:
            base_coverage = 7
            kind = "one outside pair"
        else:
            base_coverage = -1
            kind = "two outside pairs: 13 if intersecting, 14 if disjoint"
        rows.append(
            {
                "selected_T_containing_A5_count": degree,
                "final_g2": g2,
                "allowed_missing_six_min": lower,
                "allowed_missing_six_max": upper,
                "T_block_base_coverage": base_coverage,
                "configuration": kind,
            }
        )
    return rows


def _exact_structure() -> dict[str, object]:
    problem = source_problem(SOURCE)
    candidates = tuple(problem["candidates"])
    candidate_set = set(candidates)
    base_g = _base_link_g_vector(problem, T)
    if base_g != (2, 3, 4):
        raise AssertionError("triple-123 base link changed")

    sixes = tuple(T | mask(outside) for outside in combinations(OUTSIDE, 3))
    if len(sixes) != 84:
        raise AssertionError("triple-123 six-fibre count changed")
    containing = tuple(value for value in candidates if value & T == T)
    if len(containing) != 36:
        raise AssertionError("triple-123 containing-A5 count changed")
    pair_by_block = {
        value: tuple(vertex for vertex in vertices(value) if vertex not in T_VERTICES)
        for value in containing
    }
    if set(pair_by_block.values()) != set(combinations(OUTSIDE, 2)):
        raise AssertionError("triple-123 block/pair bijection changed")

    killers: dict[tuple[int, int, int], tuple[int, ...]] = {}
    for outside in combinations(OUTSIDE, 3):
        six = T | mask(outside)
        rows = tuple(
            value for value in candidates
            if value & T != T and six & value == value
        )
        killers[outside] = rows
    if any(len(rows) != 3 for rows in killers.values()):
        raise AssertionError("triple-123 noncontaining killer fibres changed")
    killer_union = set().union(*map(set, killers.values()))
    if len(killer_union) != 252:
        raise AssertionError("triple-123 noncontaining killer count changed")
    if any(value not in candidate_set for value in killer_union):
        raise AssertionError("triple-123 killer left candidate ground")

    coverage_histogram = Counter()
    pair_pair_histogram = Counter()
    for degree in (0, 1, 2):
        for selected in combinations(containing, degree):
            pairs = tuple(pair_by_block[value] for value in selected)
            covered = {
                outside for outside in killers
                if any(set(pair) <= set(outside) for pair in pairs)
            }
            coverage_histogram[(degree, len(covered))] += 1
            if degree == 2:
                pair_pair_histogram[
                    "intersecting" if set(pairs[0]) & set(pairs[1]) else "disjoint"
                ] += 1
    if coverage_histogram != {(0, 0): 1, (1, 7): 36, (2, 13): 252, (2, 14): 378}:
        raise AssertionError("triple-123 base-coverage census changed")
    if pair_pair_histogram != {"intersecting": 252, "disjoint": 378}:
        raise AssertionError("triple-123 pair-pair census changed")

    zero_generators = _generators((tuple(range(5, 13)),))
    active_generators = _generators(((5, 6, 7), tuple(range(8, 13))))
    zero_orbits = {degree: _configuration_orbits(degree, zero_generators) for degree in (1, 2)}
    active_orbits = {degree: _configuration_orbits(degree, active_generators) for degree in (1, 2)}
    if (len(zero_orbits[1]), len(zero_orbits[2])) != (2, 5):
        raise AssertionError("zero-support T-configuration orbit count changed")
    if (len(active_orbits[1]), len(active_orbits[2])) != (5, 24):
        raise AssertionError("active-support T-configuration orbit count changed")

    active_two_profiles = Counter()
    active_two_rows = []
    for orbit in active_orbits[2]:
        representative = min(orbit)
        intersection = len(set(representative[0]) & set(representative[1]))
        covers_fixed = any(set(pair) <= set(FREE_FIBRE) for pair in representative)
        active_two_profiles[(intersection, covers_fixed)] += 1
        active_two_rows.append(
            {
                "representative_outside_pairs": [list(pair) for pair in representative],
                "orbit_size": len(orbit),
                "pair_intersection_size": intersection,
                "base_covers_fixed_567_fibre": covers_fixed,
                "outside_fibre_budget": 1 if intersection else 0,
            }
        )
    if active_two_profiles != {(1, False): 10, (1, True): 3, (0, True): 4, (0, False): 7}:
        raise AssertionError("active d=2 profile census changed")

    triple012 = triple012_cut()
    patterns = tuple(
        tuple(row) for row in triple012["normalized_active_subset_orbit_representatives"]
    )
    if patterns != ((227,), (229,), (227, 229), (229, 230), (227, 229, 230)):
        raise AssertionError("triple-012 active branch patterns changed")
    branch_rows = [
        {
            "triple012_branch": "zero_support",
            "T123_configuration_orbits_by_d": {"0": 1, "1": 2, "2": 5},
            "refined_configuration_branch_count": 8,
        }
    ]
    for pattern in patterns:
        fixed_567 = 230 in pattern
        d2_count = 17 if fixed_567 else 24
        branch_rows.append(
            {
                "triple012_branch": list(pattern),
                "fixed_123_killer_fibre_567": fixed_567,
                "T123_configuration_orbits_by_d": {"0": 1, "1": 5, "2": d2_count},
                "excluded_disjoint_d2_orbits": 7 if fixed_567 else 0,
                "refined_configuration_branch_count": 1 + 5 + d2_count,
            }
        )
    if sum(int(row["refined_configuration_branch_count"]) for row in branch_rows) != 144:
        raise AssertionError("combined 012/123 branch count changed")

    return {
        "triple": list(T_VERTICES),
        "base_link_g": list(base_g),
        "outside_vertex_set": list(OUTSIDE),
        "allowed_six_fibre_count": len(sixes),
        "T_containing_candidate_A5_count": len(containing),
        "noncontaining_killer_candidate_A5_count": len(killer_union),
        "noncontaining_killers_per_fibre": 3,
        "strict_table": _strict_table(base_g),
        "T_block_base_coverage_census": [
            {"d": degree, "base_coverage": coverage, "configuration_count": count}
            for (degree, coverage), count in sorted(coverage_histogram.items())
        ],
        "outside_fibre_budget": {
            "d0": 4,
            "d1": 2,
            "d2_intersecting_pairs": 1,
            "d2_disjoint_pairs": 0,
        },
        "zero_support_configuration_orbits": {"d0": 1, "d1": 2, "d2": 5, "total": 8},
        "active_support_configuration_orbits": {"d0": 1, "d1": 5, "d2": 24, "total": 30},
        "active_d2_orbits": active_two_rows,
        "active_d2_profile_histogram": [
            {
                "pair_intersection_size": key[0],
                "base_covers_fixed_567_fibre": key[1],
                "orbit_count": count,
            }
            for key, count in sorted(active_two_profiles.items())
        ],
        "combined_triple012_triple123_branch_rows": branch_rows,
        "combined_refined_configuration_branch_count": 144,
    }


def _near_seed_localization() -> dict[str, object]:
    problem = source_problem(SOURCE)
    containing = tuple(value for value in NEAR_SEED if value & T == T)
    pairs = tuple(
        tuple(vertex for vertex in vertices(value) if vertex not in T_VERTICES)
        for value in containing
    )
    base_fibres = {
        outside for outside in combinations(OUTSIDE, 3)
        if any(set(pair) <= set(outside) for pair in pairs)
    }
    outside_killers = []
    for value in NEAR_SEED:
        if value & T == T:
            continue
        support = set(vertices(value))
        for outside in combinations(OUTSIDE, 3):
            if support <= set(T_VERTICES) | set(outside) and outside not in base_fibres:
                outside_killers.append((value, outside))
    if containing != (2094, 302) or pairs != ((5, 11), (5, 8)):
        raise AssertionError("near-seed T-containing blocks changed")
    if sorted(outside_killers) != [(342, (4, 6, 8)), (1578, (5, 9, 10))]:
        raise AssertionError("near-seed outside killer fibres changed")
    if len(base_fibres) != 13:
        raise AssertionError("near-seed T-base coverage changed")
    return {
        "T_containing_masks": list(containing),
        "outside_pair_labels": [list(pair) for pair in pairs],
        "pair_intersection": [5],
        "base_missing_six_fibres": len(base_fibres),
        "selected_outside_killer_fibres": [
            {"mask": value, "outside_triple": list(outside)}
            for value, outside in sorted(outside_killers)
        ],
        "actual_missing_six_fibres": len(base_fibres) + len(outside_killers),
        "strict_max_missing_six_fibres": 14,
        "strict_defect": 1,
        "repair_conclusion_holding_T_blocks": "At least one of masks 342 and 1578 must be removed or moved onto a base-covered fibre.",
    }


def _mathematical_digest(document: Mapping[str, object]) -> str:
    return _digest({key: value for key, value in document.items() if key != "mathematical_digest"})


def fresh_document() -> dict[str, object]:
    document: dict[str, object] = {
        "schema": "coxeter.h9_14.record0_a2.source29_triple123_cut.v1",
        "kind": "h9_14_record0_a2_source29_triple123_exact_refinement",
        "claim": "The strict triple-123 link refines all six triple-012 branches into 144 exact pair-configuration branches with outside-fibre budgets 4,2,1,0; source 29 remains OPEN.",
        "record_index": 0,
        "schedule_a": 2,
        "source_orbit_index": SOURCE,
        "triple123_structure": _exact_structure(),
        "near_seed_localization": _near_seed_localization(),
        "frontier": {"status": "OPEN", "frontier_changed": False},
        "theorem": {
            "source29_excluded": False,
            "strict_link_seed_found": False,
            "compact_incidence_witness_found": False,
            "hyperbolic_coxeter_polytope_found": False,
            "exact_branch_refinement_proved": True,
        },
        "next_finite_subproblem": [
            "Within each of the 144 pair-configuration branches, quotient the at-most 4/2/1/0 outside killer-fibre supports by its residual stabilizer and impose the exact order-6/7/8 shadows.",
            "If a strict seed survives, only then impose blocker purity, 176 facets, and ridge degree 0/2.",
        ],
        "nonclaims": [
            "No search timeout, UNKNOWN status, or failed heuristic repair is mathematical evidence.",
            "The 144 rows refine the six branches; they do not exclude them.",
            "No compact incidence, Coxeter-label, Gram, or hyperbolic witness is produced.",
        ],
    }
    document["mathematical_digest"] = _mathematical_digest(document)
    return document


def verify_document(document: Mapping[str, object]) -> dict[str, object]:
    if document != fresh_document():
        raise AssertionError("frozen triple-123 refinement differs from exact reconstruction")
    if document.get("mathematical_digest") != EXPECTED_MATHEMATICAL_DIGEST:
        raise AssertionError("triple-123 mathematical digest changed")
    return {
        "status": "PASS",
        "source29_status": "OPEN",
        "combined_refined_configuration_branch_count": 144,
        "strict_link_seed_found": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--print-fresh", action="store_true")
    parser.add_argument("--write-fresh", action="store_true")
    args = parser.parse_args()
    if args.print_fresh:
        print(json.dumps(fresh_document(), indent=2, sort_keys=True))
    elif args.write_fresh:
        OUTPUT.write_text(json.dumps(fresh_document(), indent=2, sort_keys=True) + "\n")
        print(OUTPUT)
    else:
        print(json.dumps(verify_document(json.loads(OUTPUT.read_text())), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
