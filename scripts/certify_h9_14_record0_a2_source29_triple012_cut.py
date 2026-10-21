#!/usr/bin/env python3
"""Certify the source-29 triple-012 cut and a radius-two localization.

The default audit is small.  ``--replay-radius-two`` additionally exhausts
the canonical one- and two-replacement neighbourhood of the frozen near
seed.  The replay uses integer bitsets and checks an explicit 1 GiB RSS
ceiling at every outer stage; it never promotes a timeout to evidence.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import reduce
from hashlib import sha256
from itertools import combinations
import json
from operator import or_
from pathlib import Path
import resource
import sys
from typing import Mapping

from coxeter_polytopes.combinatorics.missing_edges import macaulay_raise
from scripts.certify_h9_14_record0_a2_open_strong_audit import (
    _shadow_counts,
    _strict_failures,
)
from scripts.explore_h9_14_record0_a2_dual_facets import (
    MISSING_EDGES,
    SOURCE_A4_FAMILIES,
    mask,
    source_stabilizer_generators,
)
from scripts.explore_h9_14_record0_a2_open_primal_strong import (
    N,
    _base_link_g_vector,
    independent,
    source_problem,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    minimal_transversals,
    vertices,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record0_a2_source29_triple012_cut.json"
RSS_LIMIT_BYTES = 1 << 30
SOURCE = 29
TRIPLE = mask((0, 1, 2))
FREE = tuple(range(5, 13))
NEAR_SEED = (
    3296, 7328, 6337, 1578, 2094, 4944, 1872, 6277, 6157, 302,
    15488, 595, 6188, 7312, 789, 3114, 12332, 9040, 856, 342, 3176,
)
EXPECTED_RADIUS_TWO = {
    "all_ground_one_replacement_shadow_exact": 599,
    "triple012_forbidden_one_replacement_shadow_exact": 3,
    "post_cut_one_replacement_shadow_exact": 596,
    "post_cut_one_replacement_strict_pass": 0,
    "post_cut_two_replacement_shadow_exact": 515354,
    "post_cut_two_replacement_degree_reject": 167656,
    "post_cut_two_replacement_triple_g3_reject": 347698,
    "post_cut_two_replacement_strict_pass": 0,
    "minimum_triple_penalty": 1,
    "minimum_failure_count": 1,
    "minimum_trial_masks": [
        7328, 6337, 1578, 2094, 4944, 1872, 6277, 6157, 302, 15488,
        595, 6188, 7312, 789, 3114, 9040, 856, 342, 3176, 8364, 13444,
    ],
}
EXPECTED_MATHEMATICAL_DIGEST = "584e9901f4091e0ddca8b294cda77d70c6498847fe3e92f46b6b26164630617c"


def _digest(value: object) -> str:
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _rss_bytes() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value * 1024)


def _rss_guard() -> None:
    if _rss_bytes() > RSS_LIMIT_BYTES:
        raise MemoryError("source-29 replay exceeded its explicit 1 GiB RSS ceiling")


def _signature(value: int) -> tuple[int, int, int, int, int]:
    return (
        sum(value >> i & 1 for i in (0, 1)),
        value >> 2 & 1,
        sum(value >> i & 1 for i in (3, 4)),
        sum(value >> i & 1 for i in FREE),
        value >> 13 & 1,
    )


def _analytic_cut() -> dict[str, object]:
    problem = source_problem(SOURCE)
    candidates = tuple(problem["candidates"])
    generators, stabilizer_order = source_stabilizer_generators(problem["fours"])
    census = Counter(_signature(value) for value in candidates)
    if len(census) != 22 or sum(census.values()) != 1600:
        raise AssertionError("source-29 candidate orbit census changed")
    base_g = _base_link_g_vector(problem, TRIPLE)
    if base_g != (1, 1, 1):
        raise AssertionError("triple-012 base link changed")

    forbidden = tuple(value for value in candidates if value & TRIPLE == TRIPLE)
    killers: dict[tuple[int, int, int], tuple[int, ...]] = {}
    for free_triple in combinations(FREE, 3):
        six = TRIPLE | mask(free_triple)
        rows = tuple(
            value for value in candidates
            if value & TRIPLE != TRIPLE and six & value == value
        )
        killers[free_triple] = rows
    if len(forbidden) != 28 or {_signature(v) for v in forbidden} != {(2, 1, 0, 2, 0)}:
        raise AssertionError("triple-012 forbidden orbit changed")
    if len(killers) != 56 or any(len(rows) != 3 for rows in killers.values()):
        raise AssertionError("triple-012 six-shadow fibres changed")
    killer_union = set().union(*map(set, killers.values()))
    if len(killer_union) != 168 or Counter(_signature(v) for v in killer_union) != {
        (2, 0, 0, 3, 0): 56,
        (1, 1, 0, 3, 0): 112,
    }:
        raise AssertionError("triple-012 killer orbits changed")

    canonical = tuple(sorted(killers[(5, 6, 7)]))
    if canonical != (227, 229, 230):
        raise AssertionError("canonical killer fibre changed")
    subset_representatives = ((227,), (229,), (227, 229), (229, 230), canonical)
    return {
        "source_stabilizer": {
            "factors": ["S2_{0,1}", "S2_{3,4}", "S8_{5,...,12}"],
            "order": stabilizer_order,
            "involutory_generator_count": len(generators),
        },
        "candidate_A5_orbit_count": len(census),
        "candidate_A5_count": len(candidates),
        "candidate_A5_orbits": [
            {"signature": list(signature), "size": census[signature]}
            for signature in sorted(census)
        ],
        "triple": [0, 1, 2],
        "base_link_g": list(base_g),
        "selected_incidence_conclusion": "g2=1-d>=1 forces d=0",
        "forbidden_containing_triple_orbit": {
            "signature": [2, 1, 0, 2, 0],
            "size": len(forbidden),
        },
        "allowed_sixes_containing_triple": len(killers),
        "surviving_killer_orbits": [
            {"signature": [2, 0, 0, 3, 0], "size": 56},
            {"signature": [1, 1, 0, 3, 0], "size": 112},
        ],
        "strict_g3_conclusion": "g3=1-m>=0 forces at most one free-triple fibre",
        "canonical_free_triple": [5, 6, 7],
        "canonical_killer_masks": list(canonical),
        "canonical_killers": [list(vertices(value)) for value in canonical],
        "normalized_active_subset_orbit_representatives": [list(row) for row in subset_representatives],
        "normalized_active_subset_orbit_count": len(subset_representatives),
        "normalized_support_stabilizer_order": 2 * 2 * 6 * 120,
        "exact_branch_count": 6,
        "branches": {
            "zero_support": {
                "undecided_candidates": 1600 - 28 - 168,
                "remaining_selection_count": 21,
            },
            "five_normalized_active_support_branches": {
                "undecided_candidates_after_fixing_three_fibre_bits": 1600 - 28 - 168,
                "pre_fix_live_candidates": 1600 - 28 - (168 - 3),
                "fixed_selected_counts": [1, 1, 2, 2, 3],
                "remaining_selection_counts": [20, 20, 19, 19, 18],
            },
        },
    }


def _near_seed() -> dict[str, object]:
    problem = source_problem(SOURCE)
    additions = tuple(NEAR_SEED)
    if len(additions) != 21 or len(set(additions)) != 21 or not set(additions) <= set(problem["candidates"]):
        raise AssertionError("near seed left the source-29 candidate ground")
    shadow_counts, missing_sixes = _shadow_counts(problem, additions)
    failures = _strict_failures(problem, additions, missing_sixes)
    if shadow_counts != {"6": 147, "7": 420, "8": 629}:
        raise AssertionError("near-seed shadow schedule changed")
    if failures != [{"face": [1, 2, 3], "g": [2, 1, -1], "selected_A5_incidence": 2}]:
        raise AssertionError("near-seed strict transcript changed")
    blocker = minimal_transversals((*problem["base"], *additions))
    histogram = {str(k): v for k, v in sorted(Counter(x.bit_count() for x in blocker).items())}
    if histogram != {"4": 9, "5": 99, "6": 5}:
        raise AssertionError("near-seed blocker histogram changed")
    return {
        "additional_A5_masks": list(additions),
        "additional_A5s": [list(vertices(value)) for value in additions],
        "scheduled_missing_face_counts": shadow_counts,
        "vertex_strict_links_pass": True,
        "edge_strict_links_pass": True,
        "strict_failure_count": 1,
        "strict_failure": failures[0],
        "strict_penalty": 1,
        "minimal_transversal_order_histogram": histogram,
        "compact_incidence_witness": False,
    }


def _radius_two_replay() -> dict[str, object]:
    """Exhaust the canonical distance-one/two neighbourhood of NEAR_SEED."""

    problem = source_problem(SOURCE)
    candidates = tuple(problem["candidates"])
    cindex = {value: i for i, value in enumerate(candidates)}
    seed = tuple(cindex[value] for value in NEAR_SEED)
    selected = set(seed)
    base_low = (*MISSING_EDGES, *problem["fours"])
    grounds = {
        order: tuple(
            mask(chosen) for chosen in combinations(range(N), order)
            if not any(mask(chosen) & row == row for row in base_low)
        )
        for order in (6, 7, 8)
    }
    ground_index = {order: {face: i for i, face in enumerate(rows)} for order, rows in grounds.items()}
    covers = {
        order: tuple(
            tuple(
                ground_index[order][face]
                for extension in combinations(tuple(v for v in range(N) if not five >> v & 1), order - 5)
                if (face := five | mask(extension)) in ground_index[order]
            )
            for five in candidates
        )
        for order in (6, 7, 8)
    }
    cover_masks = {order: tuple(sum(1 << i for i in row) for row in covers[order]) for order in covers}
    faces = {
        size: tuple(mask(chosen) for chosen in combinations(range(N), size) if independent(mask(chosen)))
        for size in (1, 2, 3)
    }
    incidence = {
        size: tuple(tuple(i for i, face in enumerate(faces[size]) if five & face == face) for five in candidates)
        for size in faces
    }
    base_g = {size: tuple(_base_link_g_vector(problem, face) for face in faces[size]) for size in faces}
    caps = {
        1: tuple(g[3] - 1 for g in base_g[1]),
        2: tuple(g[2] - 1 for g in base_g[2]),
        3: tuple(g[1] - 1 for g in base_g[3]),
    }
    six_triples = tuple(
        tuple(faces[3].index(mask(chosen)) for chosen in combinations(vertices(six), 3))
        for six in grounds[6]
    )
    six_by_triple = tuple(
        sum(1 << i for i, rows in enumerate(six_triples) if j in rows)
        for j in range(len(faces[3]))
    )
    _rss_guard()

    one_shadow = one_forbidden_shadow = one_strict = 0
    for slot, old in enumerate(seed):
        remaining = tuple(seed[i] for i in range(21) if i != slot)
        union_base = tuple(reduce(or_, (cover_masks[o][i] for i in remaining), 0) for o in (6, 7, 8))
        remaining_degrees = {
            size: [sum(j in incidence[size][candidate] for candidate in remaining) for j in range(len(faces[size]))]
            for size in faces
        }
        for new in range(len(candidates)):
            if new in selected:
                continue
            unions = tuple(union_base[j] | cover_masks[o][new] for j, o in enumerate((6, 7, 8)))
            if tuple(value.bit_count() for value in unions) != (147, 420, 629):
                continue
            if candidates[new] & TRIPLE == TRIPLE:
                one_forbidden_shadow += 1
                continue
            one_shadow += 1
            degree_bad = any(
                remaining_degrees[size][j] + 1 > caps[size][j]
                for size in faces for j in incidence[size][new]
            )
            if degree_bad:
                continue
            penalty = 0
            new_triples = set(incidence[3][new])
            for j, g in enumerate(base_g[3]):
                d = remaining_degrees[3][j] + (j in new_triples)
                missing = (unions[0] & six_by_triple[j]).bit_count()
                g2 = g[1] - d
                g3 = g[2] + 5 * d - missing
                penalty += max(0, -g3) + max(0, g3 - macaulay_raise(g2, 2))
            if not penalty:
                one_strict += 1
    if one_strict:
        raise AssertionError("strict source-29 seed found at radius one")
    _rss_guard()

    two_shadow = degree_reject = triple_reject = strict_pass = 0
    minimum: tuple[int, int] | None = None
    minimum_trial: tuple[int, ...] | None = None
    targets = (147, 420, 629)
    orders = (6, 7, 8)
    for removal_number, (left, right) in enumerate(combinations(range(21), 2), 1):
        remaining = tuple(seed[i] for i in range(21) if i not in (left, right))
        base = tuple(reduce(or_, (cover_masks[o][i] for i in remaining), 0) for o in orders)
        deficits = tuple(target - value.bit_count() for target, value in zip(targets, base, strict=True))
        remaining_degrees = {
            size: [sum(j in incidence[size][candidate] for candidate in remaining) for j in range(len(faces[size]))]
            for size in faces
        }
        eligible = []
        for candidate in range(len(candidates)):
            if candidate in selected or candidates[candidate] & TRIPLE == TRIPLE:
                continue
            increments = tuple((cover_masks[o][candidate] & ~base[j]).bit_count() for j, o in enumerate(orders))
            if all(x <= d for x, d in zip(increments, deficits, strict=True)):
                eligible.append((candidate, increments))
        for ai, (a, ia) in enumerate(eligible):
            ia_sets = {size: set(incidence[size][a]) for size in faces}
            for b, ib in eligible[ai + 1:]:
                if any(x + y < d for x, y, d in zip(ia, ib, deficits, strict=True)):
                    continue
                unions = tuple(base[j] | cover_masks[o][a] | cover_masks[o][b] for j, o in enumerate(orders))
                if any(unions[j].bit_count() != targets[j] for j in range(3)):
                    continue
                two_shadow += 1
                degree_bad = False
                for size in faces:
                    affected = ia_sets[size] | set(incidence[size][b])
                    if any(
                        remaining_degrees[size][j] + (j in ia_sets[size]) + (j in incidence[size][b]) > caps[size][j]
                        for j in affected
                    ):
                        degree_bad = True
                        break
                if degree_bad:
                    degree_reject += 1
                    continue
                penalty = failure_count = 0
                for j, g in enumerate(base_g[3]):
                    d = remaining_degrees[3][j] + (j in ia_sets[3]) + (j in incidence[3][b])
                    missing = (unions[0] & six_by_triple[j]).bit_count()
                    g2 = g[1] - d
                    g3 = g[2] + 5 * d - missing
                    bad = max(0, -g3) + max(0, g3 - macaulay_raise(g2, 2))
                    penalty += bad
                    failure_count += bool(bad)
                if not penalty:
                    strict_pass += 1
                else:
                    triple_reject += 1
                if minimum is None or (penalty, failure_count) < minimum:
                    minimum = (penalty, failure_count)
                    minimum_trial = (*remaining, a, b)
        if removal_number % 10 == 0:
            _rss_guard()
    result = {
        "all_ground_one_replacement_shadow_exact": one_shadow + one_forbidden_shadow,
        "triple012_forbidden_one_replacement_shadow_exact": one_forbidden_shadow,
        "post_cut_one_replacement_shadow_exact": one_shadow,
        "post_cut_one_replacement_strict_pass": one_strict,
        "post_cut_two_replacement_shadow_exact": two_shadow,
        "post_cut_two_replacement_degree_reject": degree_reject,
        "post_cut_two_replacement_triple_g3_reject": triple_reject,
        "post_cut_two_replacement_strict_pass": strict_pass,
        "minimum_triple_penalty": minimum[0] if minimum else None,
        "minimum_failure_count": minimum[1] if minimum else None,
        "minimum_trial_masks": [candidates[i] for i in minimum_trial] if minimum_trial else [],
    }
    if result != EXPECTED_RADIUS_TWO:
        raise AssertionError(f"radius-two transcript changed: {result}")
    return result


def _mathematical_digest(document: Mapping[str, object]) -> str:
    return _digest({key: value for key, value in document.items() if key != "mathematical_digest"})


def fresh_document() -> dict[str, object]:
    document: dict[str, object] = {
        "schema": "coxeter.h9_14.record0_a2.source29_triple012_cut.v2",
        "kind": "h9_14_record0_a2_source29_exact_localization",
        "claim": "Source 29 remains OPEN, but the strict triple-012 link gives six exhaustive symmetry-normalized branches and a frozen exact-shadow near seed has no strict repair within two replacements.",
        "record_index": 0,
        "schedule_a": 2,
        "source_orbit_index": SOURCE,
        "representative_A4s": [list(face) for face in SOURCE_A4_FAMILIES[SOURCE]],
        "analytic_triple012_cut": _analytic_cut(),
        "near_seed": _near_seed(),
        "radius_two_localization": {
            **EXPECTED_RADIUS_TWO,
            "canonical_distance_definition": "After applying the proved triple-012 cut, replace exactly one or exactly two selected A5s by distinct candidates outside the seed and inside the post-cut ground; distance zero is the seed itself.",
            "all_ground_radius_one_cross_check": "Exactly three of the 599 all-ground shadow-exact one-replacements use the analytically forbidden d_012>0 orbit; 599=596+3. They are outside the post-cut domain and cannot be strict.",
            "conclusion": "Every strict-link scheduled-shadow seed lies in the post-triple-012-cut ground and has at least three A5 replacements from the frozen near seed.",
            "optional_full_replay": "python3 scripts/certify_h9_14_record0_a2_source29_triple012_cut.py --replay-radius-two",
            "rss_ceiling_bytes": RSS_LIMIT_BYTES,
        },
        "frontier": {"status": "OPEN", "frontier_changed": False},
        "theorem": {
            "source29_excluded": False,
            "compact_incidence_witness_found": False,
            "hyperbolic_coxeter_polytope_found": False,
            "exact_branch_reduction_proved": True,
            "radius_two_strict_repair_excluded": True,
        },
        "next_finite_subproblem": [
            "Solve the zero-support branch and the five normalized active-fibre branches with the exact order-6/7/8 shadows and strict links.",
            "Equivalently search repairs at three or more A5 replacements from the frozen near seed; then impose blocker purity, 176 facets, and ridge degree 0/2.",
        ],
        "nonclaims": [
            "This is not an exclusion of source 29 or of Record-0/a=2.",
            "The near seed is not a compact-incidence, Coxeter-label, Gram, or hyperbolic witness.",
            "No timeout or bounded search disposition is used as mathematical evidence.",
        ],
    }
    document["mathematical_digest"] = _mathematical_digest(document)
    return document


def verify_document(document: Mapping[str, object], *, replay_radius_two: bool = False) -> dict[str, object]:
    if document.get("kind") != "h9_14_record0_a2_source29_exact_localization":
        raise ValueError("wrong source-29 localization kind")
    fresh = fresh_document()
    if document != fresh:
        raise AssertionError("frozen source-29 localization differs from exact reconstruction")
    if document.get("mathematical_digest") != EXPECTED_MATHEMATICAL_DIGEST:
        raise AssertionError("source-29 mathematical digest changed")
    replay = _radius_two_replay() if replay_radius_two else None
    _rss_guard()
    return {
        "status": "PASS",
        "source29_status": "OPEN",
        "exact_branch_count": 6,
        "radius_two_replayed": replay is not None,
        "peak_rss_bytes": _rss_bytes(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay-radius-two", action="store_true")
    parser.add_argument("--print-fresh", action="store_true")
    parser.add_argument("--write-fresh", action="store_true")
    args = parser.parse_args()
    if args.print_fresh:
        print(json.dumps(fresh_document(), indent=2, sort_keys=True))
        return
    if args.write_fresh:
        OUTPUT.write_text(json.dumps(fresh_document(), indent=2, sort_keys=True) + "\n")
        print(OUTPUT)
        return
    document = json.loads(OUTPUT.read_text())
    print(json.dumps(verify_document(document, replay_radius_two=args.replay_radius_two), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
