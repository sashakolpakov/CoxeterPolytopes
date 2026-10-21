#!/usr/bin/env python3
"""Freeze the compact incidence no-go for Record-0 ``a=1`` source 2."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Mapping

from ortools.sat.python import cp_model
import z3

from scripts.explore_h9_14_record0_a1_dual_facets import build_dual_model
from scripts.explore_h9_14_record0_a1_orbit_cegar import source_group_order
from scripts.independent_h9_14_record0_a1_dual_z3 import (
    SOURCE_FOURS,
    build_model as build_independent_model,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record0_a1_source2_incidence_no_go.json"
SOURCE_PATHS = (
    ROOT / "scripts/explore_h9_14_record0_a1_dual_facets.py",
    ROOT / "scripts/independent_h9_14_record0_a1_dual_z3.py",
    ROOT / "scripts/explore_h9_14_open_positive_completion.py",
    ROOT / "scripts/explore_h9_14_open_positive_a4.py",
    ROOT / "data/enumerations/h9_14_nerve_support.json",
    ROOT / "src/coxeter_polytopes/combinatorics/missing_edges.py",
)
EXPECTED_COUNTS = {
    "candidate_facet_complements": 754,
    "independent_five_sets": 1617,
    "face_witness_constraints": 1324,
    "all_order_lanner_closure_constraints": 6796,
    "ridge_zero_or_two_constraints": 1491,
    "stabilizer_generators": 10,
    "generator_moved_facet_pair_counts": [
        201,
        166,
        166,
        166,
        166,
        166,
        166,
        166,
        166,
        165,
    ],
}
EXPECTED_CERTIFICATE_DIGEST = "9c30e26f9466338adac2acf7aa1fc033cdeee201520002baa350bdfb1b362ab5"


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _production_counts(data: Mapping[str, object]) -> dict[str, object]:
    return {
        "candidate_facet_complements": len(data["facets"]),
        "independent_five_sets": len(data["five_ground"]),
        "face_witness_constraints": data["face_witness_constraint_count"],
        "all_order_lanner_closure_constraints": data[
            "lanner_order_constraint_count"
        ],
        "ridge_zero_or_two_constraints": data["ridge_constraint_count"],
        "stabilizer_generators": data["stabilizer_generator_count"],
        "generator_moved_facet_pair_counts": data[
            "generator_moved_pair_counts"
        ],
    }


def _independent_counts(data: Mapping[str, object]) -> dict[str, object]:
    return {
        "candidate_facet_complements": data[
            "candidate_facet_complement_count"
        ],
        "independent_five_sets": data["independent_five_count"],
        "face_witness_constraints": data["face_witness_constraint_count"],
        "all_order_lanner_closure_constraints": data[
            "lanner_order_constraint_count"
        ],
        "ridge_zero_or_two_constraints": data["ridge_constraint_count"],
        "stabilizer_generators": data["stabilizer_generator_count"],
        "generator_moved_facet_pair_counts": data[
            "generator_moved_pair_counts"
        ],
    }


def fresh_record() -> dict[str, object]:
    production = build_dual_model(2)
    production_counts = _production_counts(production)
    if production_counts != EXPECTED_COUNTS:
        raise AssertionError("production source-2 model counts changed")
    cp_solver = cp_model.CpSolver()
    cp_solver.parameters.num_search_workers = 1
    cp_status = cp_solver.solve(production["model"])
    if cp_status != cp_model.INFEASIBLE:
        raise AssertionError(
            "production source-2 model did not terminate INFEASIBLE"
        )

    independent_solver, independent = build_independent_model(
        2, symmetry_breaking=True
    )
    independent_counts = _independent_counts(independent)
    if independent_counts != EXPECTED_COUNTS:
        raise AssertionError("independent source-2 model counts changed")
    z3_status = independent_solver.check()
    if z3_status != z3.unsat:
        raise AssertionError("independent source-2 model did not terminate unsat")

    four_mask = sum(1 << vertex for vertex in SOURCE_FOURS[2])
    document: dict[str, object] = {
        "kind": "h9_14_record0_a1_source2_compact_incidence_no_go",
        "record_index": 0,
        "graph6": "M?????????????o??",
        "schedule": {
            "a": 1,
            "missing_independent_five_sets": 32,
            "scheduled_facet_count": 178,
        },
        "fixed_source": {
            "source_orbit_index": 2,
            "representative_A4": list(SOURCE_FOURS[2]),
            "source_stabilizer": "S2({2,3}) x S9({4,...,12}) x S2({0,1})",
            "source_stabilizer_order": source_group_order(four_mask),
            "a1_source_orbit_count": 4,
            "other_source_orbits_not_closed_by_this_certificate": [0, 1, 3],
        },
        "source_files": {
            str(path.relative_to(ROOT)): _sha(path) for path in SOURCE_PATHS
        },
        "proof_scope": (
            "Two independently written exact Boolean models select the 178 "
            "rank-nine facet complements for the fixed A4 source. They "
            "impose all required faces through rank four, exactly 32 missing "
            "independent five-sets, ridge degree zero or two, and closure of "
            "every nonface of order at least six by a missing set of order "
            "at most five. Both use the lossless generator-lex slice of the "
            "fixed-source stabilizer and terminate UNSAT. Thus this one "
            "source cannot be the nerve of a compact Coxeter polytope."
        ),
        "evidence": {
            "production_solver": "OR-Tools CP-SAT, one worker",
            "production_terminal_status": cp_solver.status_name(cp_status),
            "independent_solver": "Z3",
            "independent_terminal_status": str(z3_status),
            "production_counts": production_counts,
            "independent_counts": independent_counts,
            "independent_reconstruction_exact_count_agreement": True,
            "symmetry_slice_lossless_reason": (
                "Every orbit under the finite fixed-source stabilizer has a "
                "global lexicographic minimum; that minimum is no larger "
                "than its image under each involutory generator."
            ),
        },
        "theorem": {
            "fixed_source_compact_incidence_no_go": True,
            "excluded_a1_source_orbits": [2],
            "remaining_a1_source_orbits": [0, 1, 3],
            "whole_record0_a1_schedule_excluded": False,
            "whole_shadow_compact_coxeter_no_go": False,
            "finite_volume_noncompact_chamber_excluded": False,
        },
    }
    document["certificate_digest"] = _digest(document)
    return document


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != (
        "h9_14_record0_a1_source2_compact_incidence_no_go"
    ):
        raise ValueError("wrong Record-0 source-2 certificate kind")
    sources = document.get("source_files")
    if not isinstance(sources, Mapping):
        raise ValueError("source-file map is malformed")
    for path in SOURCE_PATHS:
        if sources.get(str(path.relative_to(ROOT))) != _sha(path):
            raise ValueError(f"source file changed: {path.name}")
    fixed = document.get("fixed_source")
    evidence = document.get("evidence")
    theorem = document.get("theorem")
    if not isinstance(fixed, Mapping) or not isinstance(evidence, Mapping):
        raise ValueError("fixed-source evidence is malformed")
    if not isinstance(theorem, Mapping):
        raise ValueError("theorem scope is malformed")
    if (
        fixed.get("source_orbit_index") != 2
        or fixed.get("representative_A4") != [0, 1, 2, 3]
        or fixed.get("source_stabilizer_order") != 1451520
        or fixed.get("other_source_orbits_not_closed_by_this_certificate")
        != [0, 1, 3]
    ):
        raise ValueError("fixed source changed")
    if (
        evidence.get("production_terminal_status") != "INFEASIBLE"
        or evidence.get("independent_terminal_status") != "unsat"
        or evidence.get("production_counts") != EXPECTED_COUNTS
        or evidence.get("independent_counts") != EXPECTED_COUNTS
        or evidence.get("independent_reconstruction_exact_count_agreement")
        is not True
    ):
        raise ValueError("two-solver incidence evidence changed")
    if (
        theorem.get("fixed_source_compact_incidence_no_go") is not True
        or theorem.get("remaining_a1_source_orbits") != [0, 1, 3]
        or theorem.get("whole_record0_a1_schedule_excluded") is not False
        or theorem.get("whole_shadow_compact_coxeter_no_go") is not False
        or theorem.get("finite_volume_noncompact_chamber_excluded") is not False
    ):
        raise ValueError("fixed-source theorem scope changed")
    digest = document.get("certificate_digest")
    undigested = dict(document)
    undigested.pop("certificate_digest", None)
    if digest != _digest(undigested):
        raise ValueError("certificate digest does not match document")
    if digest != EXPECTED_CERTIFICATE_DIGEST:
        raise ValueError("certificate digest changed")
    return {
        "record_index": 0,
        "schedule_a": 1,
        "source_orbit_index": 2,
        "representative_A4": [0, 1, 2, 3],
        "production_terminal_status": "INFEASIBLE",
        "independent_terminal_status": "unsat",
        "fixed_source_compact_incidence_no_go": True,
        "remaining_a1_source_orbits": [0, 1, 3],
        "whole_record0_a1_schedule_excluded": False,
        "whole_shadow_compact_coxeter_no_go": False,
        "finite_volume_noncompact_chamber_excluded": False,
        "certificate_digest": digest,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--write", type=Path)
    args = parser.parse_args()
    if args.fresh:
        document = fresh_record()
        if args.write is not None:
            args.write.write_text(
                json.dumps(document, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        print(json.dumps(document, indent=2, sort_keys=True))
        return
    document = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
