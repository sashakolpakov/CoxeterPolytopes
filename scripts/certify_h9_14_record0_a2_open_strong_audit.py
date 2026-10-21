#!/usr/bin/env python3
"""Certify the stronger exact audit of the four open Record-0/a=2 sources.

This certificate freezes structural reductions and one positive shadow-design
witness.  It deliberately does not promote solver timeouts to exclusions.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Mapping

from coxeter_polytopes.combinatorics.missing_edges import macaulay_raise
from scripts.explore_h9_14_record0_a2_dual_facets import (
    MISSING_EDGES,
    SOURCE_A4_FAMILIES,
    mask,
)
from scripts.explore_h9_14_record0_a2_open_primal_strong import (
    FULL,
    N,
    OPEN_SOURCES,
    SCHEDULED_F_VECTOR,
    _base_link_g_vector,
    independent,
    source_problem,
)
from scripts.explore_h9_14_record14_fixed_orbit_z3 import (
    minimal_transversals,
    vertices,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "h9_14_record0_a2_open_strong_audit.json"
)
EXPECTED_MATHEMATICAL_DIGEST = (
    "c45702ef85e17f990ec8216c85292b09459059738ceba98a8acaffe8cf066475"
)
SOURCE29_SHADOW_WITNESS = (
    158,
    3107,
    2087,
    5157,
    244,
    541,
    188,
    4141,
    1549,
    1174,
    1543,
    6179,
    1306,
    12836,
    14432,
    2466,
    6592,
    9496,
    11584,
    1558,
    4649,
)
EXPECTED_SOURCE29_STRICT_FAILURE_DIGEST = (
    "c87d6941843fca0a43d64039e87e2149431b8bbbf0caa4bf9c8c09e99d618102"
)
SOURCE_FILES = (
    ROOT / "scripts" / "explore_h9_14_record0_a2_dual_facets.py",
    ROOT / "scripts" / "explore_h9_14_record0_a2_open_primal_strong.py",
    ROOT / "scripts" / "explore_h9_14_record0_a2_open_primal_full.py",
    ROOT / "scripts" / "explore_h9_14_record0_a2_open_full_links_z3.py",
    ROOT / "scripts" / "explore_h9_14_record0_a2_open_shadow_design.py",
    ROOT / "src" / "coxeter_polytopes" / "combinatorics" / "facet_completion.py",
    ROOT / "src" / "coxeter_polytopes" / "combinatorics" / "missing_edges.py",
)


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _mathematical_digest(document: Mapping[str, object]) -> str:
    return _digest(
        {
            key: value
            for key, value in document.items()
            if key not in {
                "mathematical_digest",
                "run_metadata",
                "source_files",
            }
        }
    )


def _source_counts(source: int) -> dict[str, object]:
    problem = source_problem(source)
    base_low = (*MISSING_EDGES, *problem["fours"])
    face_candidates = {}
    for order in range(1, 10):
        face_candidates[str(order)] = sum(
            1
            for chosen in combinations(range(N), order)
            if not any(mask(chosen) & row == row for row in base_low)
        )
    facets = tuple(problem["facet_complements"])
    ridge_count = sum(
        any(facet & mask(chosen) == facet for facet in facets)
        for chosen in combinations(range(N), 6)
    )
    if face_candidates["9"] != len(facets):
        raise AssertionError("order-nine face ground differs from facet ground")
    return {
        "forced_A5_shadow_count": len(problem["shadow"]),
        "additional_A5_count": problem["addition_count"],
        "candidate_A5_count": len(problem["candidates"]),
        "candidate_facet_complement_count": len(facets),
        "scheduled_face_candidate_counts": face_candidates,
        "face_extension_constraint_counts": {
            str(order): face_candidates[str(order)] for order in range(1, 9)
        },
        "large_set_closure_constraint_count": len(
            problem["small_transversals"]
        ),
        "ridge_constraint_count": ridge_count,
        "strict_link_figure_counts_by_face_size": {
            "1": 14,
            "2": 89,
            "3": 341,
        },
        "strict_link_macaulay_table_constraint_count": 902,
    }


def _shadow_counts(
    problem: Mapping[str, object], additions: tuple[int, ...]
) -> tuple[dict[str, int], tuple[int, ...]]:
    base_low = (*MISSING_EDGES, *problem["fours"])
    counts = {}
    missing_sixes: tuple[int, ...] = ()
    for order in (6, 7, 8):
        ground = tuple(
            mask(chosen)
            for chosen in combinations(range(N), order)
            if not any(mask(chosen) & row == row for row in base_low)
        )
        missing = tuple(
            face
            for face in ground
            if any(face & addition == addition for addition in additions)
        )
        counts[str(order)] = len(missing)
        if order == 6:
            missing_sixes = missing
    return counts, missing_sixes


def _strict_failures(
    problem: dict[str, object],
    additions: tuple[int, ...],
    missing_sixes: tuple[int, ...],
) -> list[dict[str, object]]:
    failures = []
    for face_size in (1, 2, 3):
        for chosen in combinations(range(N), face_size):
            face = mask(chosen)
            if not independent(face):
                continue
            g_vector = list(_base_link_g_vector(problem, face))
            incidence = sum(
                addition & face == face for addition in additions
            )
            if face_size == 1:
                g_vector[3] -= incidence
                passes = (
                    1
                    <= g_vector[3]
                    <= macaulay_raise(g_vector[2], 3)
                )
            elif face_size == 2:
                g_vector[2] -= incidence
                passes = (
                    1
                    <= g_vector[2]
                    <= macaulay_raise(g_vector[1], 2)
                )
            else:
                g_vector[1] -= incidence
                g_vector[2] += 5 * incidence - sum(
                    six & face == face for six in missing_sixes
                )
                passes = (
                    1
                    <= g_vector[1]
                    <= macaulay_raise(g_vector[0], 1)
                    and 0
                    <= g_vector[2]
                    <= macaulay_raise(g_vector[1], 2)
                )
            if not passes:
                failures.append(
                    {
                        "face": list(chosen),
                        "g": g_vector,
                        "selected_A5_incidence": incidence,
                    }
                )
    return failures


def _source29_witness() -> dict[str, object]:
    problem = source_problem(29)
    additions = tuple(SOURCE29_SHADOW_WITNESS)
    if len(additions) != problem["addition_count"] or len(set(additions)) != len(
        additions
    ):
        raise AssertionError("source-29 witness has wrong cardinality")
    if not set(additions) <= set(problem["candidates"]):
        raise AssertionError("source-29 witness leaves the exact candidate ground")
    vertex_incidence = tuple(
        sum(addition & (1 << vertex) != 0 for addition in additions)
        for vertex in range(N)
    )
    if any(
        not problem["lower"][vertex]
        <= vertex_incidence[vertex]
        <= problem["upper"][vertex]
        for vertex in range(N)
    ):
        raise AssertionError("source-29 witness fails a vertex-link interval")
    shadow_counts, missing_sixes = _shadow_counts(problem, additions)
    expected_shadow_counts = {
        "6": 2114 - SCHEDULED_F_VECTOR[5],
        "7": 2038 - SCHEDULED_F_VECTOR[6],
        "8": 1421 - SCHEDULED_F_VECTOR[7],
    }
    if shadow_counts != expected_shadow_counts:
        raise AssertionError("source-29 shadow counts changed")
    failures = _strict_failures(problem, additions, missing_sixes)
    failure_digest = _digest(failures)
    if failure_digest != EXPECTED_SOURCE29_STRICT_FAILURE_DIGEST:
        raise AssertionError("source-29 strict-failure transcript changed")
    blocker = minimal_transversals((*problem["base"], *additions))
    blocker_histogram = {
        str(order): count
        for order, count in sorted(
            Counter(value.bit_count() for value in blocker).items()
        )
    }
    if blocker_histogram != {"4": 10, "5": 90, "6": 27}:
        raise AssertionError("source-29 blocker audit changed")
    return {
        "source_orbit_index": 29,
        "representative_A4s": [list(face) for face in SOURCE_A4_FAMILIES[29]],
        "additional_A5_masks": list(additions),
        "additional_A5s": [list(vertices(value)) for value in additions],
        "vertex_A5_incidence": list(vertex_incidence),
        "scheduled_missing_face_counts": shadow_counts,
        "scheduled_shadow_counts_pass": True,
        "vertex_link_intervals_pass": True,
        "strict_link_failure_count": len(failures),
        "strict_link_failure_histogram_by_face_size": {
            str(size): sum(len(row["face"]) == size for row in failures)
            for size in (1, 2, 3)
        },
        "strict_link_failure_digest": failure_digest,
        "minimal_transversal_order_histogram": blocker_histogram,
        "full_compact_incidence_model_pass": False,
        "scope": (
            "Exact A5 family realizing the scheduled six/seven/eight shadow "
            "counts and all vertex intervals. It is not a facet/ridge, "
            "strict-link, Coxeter-label, Gram, or hyperbolic witness."
        ),
    }


def fresh_document() -> dict[str, object]:
    rows = []
    for source in OPEN_SOURCES:
        rows.append(
            {
                "source_orbit_index": source,
                "representative_A4s": [
                    list(face) for face in SOURCE_A4_FAMILIES[source]
                ],
                "model_counts": _source_counts(source),
                "strong_exact_disposition": "open",
                "compact_incidence_sat_witness_found": False,
                "compact_incidence_unsat_certificate_found": False,
            }
        )
    document: dict[str, object] = {
        "schema": "coxeter.h9_14.record0_a2.open_strong_audit.v1",
        "kind": "h9_14_record0_a2_open_sources_strong_exact_audit",
        "claim": (
            "The four residual sources survive a strictly stronger exact "
            "audit only as OPEN; source 29 has an exact scheduled A5-shadow "
            "witness, but that witness fails strict links and purity."
        ),
        "record_index": 0,
        "graph6": "M?????????????o??",
        "schedule_a": 2,
        "strong_model": {
            "additional_A5_selection": True,
            "full_dehn_sommerville_face_vector": list(SCHEDULED_F_VECTOR),
            "every_allowed_face_order_1_through_8_extends_to_a_facet": True,
            "larger_set_lanner_order_closure": True,
            "ridge_degree_zero_or_two": True,
            "strict_vertex_edge_triple_link_g_and_macaulay": True,
            "strict_link_figure_count": 444,
            "strict_link_macaulay_table_count": 902,
            "lossless_source_stabilizer_lex_optional": True,
        },
        "rows": rows,
        "positive_shadow_witnesses": [_source29_witness()],
        "frontier": {
            "two_solver_unsat_source_orbits": [0, 8, 9, 10, 12, 23, 31, 34],
            "incidence_sat_source_orbits": [
                1, 2, 3, 4, 5, 6, 7, 11, 13, 14, 15, 16, 17, 18,
                19, 20, 21, 22, 25, 26, 27, 28, 33, 35,
            ],
            "strong_audit_open_source_orbits": list(OPEN_SOURCES),
            "partition_identity": "8 + 24 + 4 = 36",
            "frontier_changed": False,
        },
        "theorem": {
            "new_compact_incidence_exclusion": False,
            "new_compact_incidence_sat_witness": False,
            "new_hyperbolic_coxeter_polytope": False,
            "finite_volume_noncompact_checked": False,
        },
        "nonclaims": [
            "UNKNOWN and timeout outcomes are not exclusions.",
            "The source-29 shadow witness is not a compact-incidence witness.",
            "No Coxeter labels or Lorentzian Gram matrix are produced.",
            "No finite-volume noncompact chamber or cusp is tested.",
        ],
        "run_metadata": {
            "canonical": False,
            "excluded_from_mathematical_digest": True,
            "full_primal_no_lex": [
                {
                    "source_orbit_index": 24,
                    "requested_limit_seconds": 60,
                    "status": "UNKNOWN",
                    "wall_time_seconds": 60.910150,
                    "branches": 0,
                    "conflicts": 0,
                },
                {
                    "source_orbit_index": 29,
                    "requested_limit_seconds": 60,
                    "status": "UNKNOWN",
                    "wall_time_seconds": 52.111674,
                    "branches": 0,
                    "conflicts": 0,
                },
                {
                    "source_orbit_index": 30,
                    "requested_limit_seconds": 60,
                    "status": "UNKNOWN",
                    "wall_time_seconds": 51.889377,
                    "branches": 0,
                    "conflicts": 0,
                },
                {
                    "source_orbit_index": 32,
                    "requested_limit_seconds": 60,
                    "status": "UNKNOWN",
                    "wall_time_seconds": 48.512746,
                    "branches": 0,
                    "conflicts": 0,
                },
            ],
            "source24_full_primal_with_lex": {
                "requested_limit_seconds": 60,
                "status": "UNKNOWN",
                "wall_time_seconds": 60.090808,
                "branches": 13978,
                "conflicts": 0,
            },
            "source24_production_z3_composition": {
                "requested_limit_seconds": 30,
                "status": "unknown",
                "reason_unknown": "timeout",
                "assertion_count": 22687,
                "build_seconds": 98.529874,
                "check_seconds": 31.933578,
            },
            "source29_shadow_staging": {
                "cp_sat_full_shadow_and_strict_links": {
                    "requested_limit_seconds": 60,
                    "status": "UNKNOWN",
                    "wall_time_seconds": 60.056414,
                    "branches": 172156,
                    "conflicts": 85,
                },
                "highs_full_shadow_and_strict_links": {
                    "requested_limit_seconds": 300,
                    "status": "time_limit",
                    "nodes": 0,
                    "primal_solution_found": False,
                },
            },
        },
        "source_files": {
            str(path.relative_to(ROOT)): _file_sha256(path)
            for path in SOURCE_FILES
        },
    }
    document["mathematical_digest"] = _mathematical_digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record0_a2_open_sources_strong_exact_audit":
        raise ValueError("wrong certificate kind")
    sources = document.get("source_files")
    if not isinstance(sources, Mapping):
        raise ValueError("source hash map is malformed")
    for path in SOURCE_FILES:
        key = str(path.relative_to(ROOT))
        if sources.get(key) != _file_sha256(path):
            raise ValueError(f"source file changed: {key}")
    digest = document.get("mathematical_digest")
    if digest != _mathematical_digest(document):
        raise ValueError("mathematical digest mismatch")
    if EXPECTED_MATHEMATICAL_DIGEST and digest != EXPECTED_MATHEMATICAL_DIGEST:
        raise ValueError("certificate is not the frozen audit")
    fresh = fresh_document()
    for key, value in fresh.items():
        if key in {"run_metadata", "source_files"}:
            continue
        if document.get(key) != value:
            raise ValueError(f"canonical field changed: {key}")
    frontier = document.get("frontier")
    theorem = document.get("theorem")
    if not isinstance(frontier, Mapping) or not isinstance(theorem, Mapping):
        raise ValueError("frontier/theorem ledger is malformed")
    if frontier.get("partition_identity") != "8 + 24 + 4 = 36":
        raise ValueError("frontier partition changed")
    if frontier.get("frontier_changed") is not False:
        raise ValueError("timeout was promoted to a frontier change")
    if any(
        theorem.get(key) is not False
        for key in (
            "new_compact_incidence_exclusion",
            "new_compact_incidence_sat_witness",
            "new_hyperbolic_coxeter_polytope",
            "finite_volume_noncompact_checked",
        )
    ):
        raise ValueError("nonclaim was promoted")
    return {
        "mathematical_digest": digest,
        "strong_audit_open_source_orbits": list(OPEN_SOURCES),
        "frontier_changed": False,
        "source29_scheduled_shadow_witness": True,
        "source29_strict_link_failure_count": 48,
        "new_hyperbolic_coxeter_polytope": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write:
        OUTPUT.write_text(
            json.dumps(fresh_document(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    if args.check or not args.write:
        print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
