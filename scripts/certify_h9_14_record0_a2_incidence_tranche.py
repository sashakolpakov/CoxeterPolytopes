#!/usr/bin/env python3
"""Freeze a two-solver compact-incidence no-go tranche for Record 0/a=2."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import isfinite
from pathlib import Path
from typing import Mapping, Sequence

from ortools.sat.python import cp_model
import z3

from scripts.explore_h9_14_open_positive_a4 import audit_record
from scripts.explore_h9_14_record0_a2_dual_facets import (
    SOURCE_A4_FAMILIES,
    build_model as build_cp_model,
    mask,
    source_representatives,
)
from scripts.independent_h9_14_record0_a2_tranche_z3 import (
    SOURCES as Z3_SOURCES,
    build_model as build_z3_model,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/h9_14_record0_a2_incidence_tranche.json"
EXCLUDED_SOURCE_ORBITS = (0, 8, 9, 10, 12, 23, 31, 34)
INCIDENCE_SAT_SOURCE_ORBITS = (
    1, 2, 3, 4, 5, 6, 7, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
    25, 26, 27, 28, 33, 35,
)
SOLVER_OPEN_SOURCE_ORBITS = (24, 29, 30, 32)
EXPECTED_SOURCE_REDUCTION = {
    "raw_family_count": 386_760,
    "interval_feasible_family_count": 384_780,
    "strict_aggregate_viable_labelled_count": 384_780,
    "source_orbit_count": 36,
    "source_orbit_size_histogram": {
        "165": 1,
        "330": 1,
        "495": 2,
        "990": 2,
        "1320": 1,
        "1980": 1,
        "2640": 1,
        "3960": 3,
        "4620": 3,
        "5775": 1,
        "6930": 2,
        "7920": 2,
        "9240": 4,
        "11550": 1,
        "13860": 3,
        "18480": 1,
        "20790": 1,
        "23100": 2,
        "27720": 3,
        "55440": 1,
    },
}
EXPECTED_SOURCE_REPRESENTATIVE_DIGEST = (
    "ec1e72a5b58761142864fb30ddabe65a1cc272587f942c0918a710f3e54c7b4f"
)
EXPECTED_MATHEMATICAL_DIGEST = (
    "f8691cb17d07cd47ca66519fdd4737641b8c5fdb69fcb68d957a0ad5dd878288"
)
SOURCE_PATHS = (
    ROOT / "scripts/explore_h9_14_open_positive_a4.py",
    ROOT / "scripts/explore_h9_14_record0_a2_dual_facets.py",
    ROOT / "scripts/independent_h9_14_record0_a2_tranche_z3.py",
    ROOT / "data/enumerations/h9_14_nerve_support.json",
    ROOT / "src/coxeter_polytopes/combinatorics/local_link.py",
    ROOT / "src/coxeter_polytopes/combinatorics/missing_edges.py",
)


def _digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _mathematical_digest(document: Mapping[str, object]) -> str:
    """Digest theorem data, excluding provenance and nondeterministic runs."""

    payload = {
        key: value
        for key, value in document.items()
        if key not in {
            "certificate_digest",
            "mathematical_digest",
            "run_metadata",
            "source_files",
        }
    }
    return _digest(payload)


def normalize_document(document: Mapping[str, object]) -> dict[str, object]:
    """Move solver telemetry outside the canonical mathematical payload."""

    normalized = json.loads(json.dumps(document))
    evidence = normalized.get("evidence")
    if not isinstance(evidence, list):
        raise ValueError("evidence is malformed")
    existing_metadata = normalized.get("run_metadata")
    existing_rows: dict[int, dict[str, object]] = {}
    if isinstance(existing_metadata, Mapping):
        rows = existing_metadata.get("per_source")
        if isinstance(rows, list):
            existing_rows = {
                int(row["source_orbit_index"]): dict(row)
                for row in rows
                if isinstance(row, Mapping) and "source_orbit_index" in row
            }
    run_rows = []
    for row in evidence:
        if not isinstance(row, dict):
            raise ValueError("evidence row is malformed")
        source = int(row["source_orbit_index"])
        metadata = existing_rows.get(
            source, {"source_orbit_index": source}
        )
        for key in (
            "production_wall_time_seconds",
            "production_branch_count",
            "production_conflict_count",
            "independent_solver_version",
        ):
            if key in row:
                metadata[key] = row.pop(key)
        run_rows.append(metadata)
    normalized["run_metadata"] = {
        "canonical": False,
        "excluded_from_mathematical_digest": True,
        "per_source": run_rows,
    }
    normalized["scout_status_ledger"] = {
        "two_solver_unsat_source_orbits": list(EXCLUDED_SOURCE_ORBITS),
        "explicit_incidence_sat_source_orbits": list(
            INCIDENCE_SAT_SOURCE_ORBITS
        ),
        "solver_open_source_orbits": list(SOLVER_OPEN_SOURCE_ORBITS),
        "partition_identity": "8 + 24 + 4 = 36",
        "sat_scope": (
            "An explicit solution of the compact-incidence overcover was "
            "found; this is not a Coxeter-label or Gram survivor."
        ),
        "open_scope": "Recorded timeouts are not exclusions.",
    }
    coverage = normalized.get("coverage")
    if not isinstance(coverage, dict):
        raise ValueError("coverage is malformed")
    coverage["incidence_sat_source_orbit_count"] = len(
        INCIDENCE_SAT_SOURCE_ORBITS
    )
    coverage["solver_open_source_orbit_count"] = len(
        SOLVER_OPEN_SOURCE_ORBITS
    )
    normalized.pop("certificate_digest", None)
    normalized.pop("mathematical_digest", None)
    normalized["mathematical_digest"] = _mathematical_digest(normalized)
    return normalized


def _source_reduction() -> dict[str, object]:
    audit = audit_record(0, 2)
    branch = next(row for row in audit["branches"] if row["a"] == 2)
    observed = {
        key: branch[key] for key in EXPECTED_SOURCE_REDUCTION
    }
    if observed != EXPECTED_SOURCE_REDUCTION:
        raise AssertionError("Record-0/a=2 source reduction changed")
    representatives = tuple(
        tuple(tuple(face) for face in family)
        for family in branch["representatives"]
    )
    if representatives != SOURCE_A4_FAMILIES:
        raise AssertionError("frozen source representatives changed")
    digest = _digest(representatives)
    if digest != EXPECTED_SOURCE_REPRESENTATIVE_DIGEST:
        raise AssertionError(
            f"source representative digest changed: {digest}"
        )
    if sum(
        int(size) * count
        for size, count in branch["source_orbit_size_histogram"].items()
    ) != branch["strict_aggregate_viable_labelled_count"]:
        raise AssertionError("source orbit sizes do not cover the viable ground")
    return {
        **observed,
        "representative_digest": digest,
        "orbit_sizes_cover_every_viable_labelled_source": True,
        "full_missing_graph_automorphism_quotient": True,
    }


def _counts(data: Mapping[str, object]) -> dict[str, object]:
    return {
        "forced_A5_shadow_count": len(data["shadow"]),
        "additional_A5_count": 38 - len(data["shadow"]),
        "candidate_facet_complement_count": len(data["facets"]),
        "face_witness_constraint_count": data[
            "face_witness_constraint_count"
        ],
        "lanner_order_constraint_count": data[
            "lanner_order_constraint_count"
        ],
        "ridge_constraint_count": data["ridge_constraint_count"],
        "source_stabilizer_order": data["stabilizer_order"],
        "stabilizer_generator_count": data["stabilizer_generator_count"],
        "generator_moved_facet_pair_counts": data[
            "generator_moved_pair_counts"
        ],
    }


def fresh_record() -> dict[str, object]:
    reduction = _source_reduction()
    representatives = source_representatives()
    evidence = []
    for source in EXCLUDED_SOURCE_ORBITS:
        fours = representatives[source]
        if tuple(tuple(face) for face in SOURCE_A4_FAMILIES[source]) != Z3_SOURCES[source]:
            raise AssertionError("independent literal source changed")

        cp_data = build_cp_model(
            fours, include_lanner_closure=True, symmetry_breaking=True
        )
        cp_counts = _counts(cp_data)
        cp_solver = cp_model.CpSolver()
        cp_solver.parameters.num_search_workers = 1
        cp_status = cp_solver.solve(cp_data["model"])
        if cp_status != cp_model.INFEASIBLE:
            raise AssertionError(f"source {source}: CP-SAT is not INFEASIBLE")

        z3_solver, z3_counts = build_z3_model(source)
        if z3_counts != {
            "source_orbit_index": source,
            "representative_A4s": [
                list(face) for face in SOURCE_A4_FAMILIES[source]
            ],
            **cp_counts,
        }:
            raise AssertionError(f"source {source}: independent counts differ")
        z3_status = z3_solver.check()
        if z3_status != z3.unsat:
            raise AssertionError(f"source {source}: Z3 is not unsat")
        evidence.append(
            {
                "source_orbit_index": source,
                "representative_A4s": [
                    list(face) for face in SOURCE_A4_FAMILIES[source]
                ],
                "model_counts": cp_counts,
                "production_solver": "OR-Tools CP-SAT; one worker; no time limit",
                "production_terminal_status": cp_solver.status_name(cp_status),
                "production_wall_time_seconds": cp_solver.wall_time,
                "production_branch_count": cp_solver.num_branches,
                "production_conflict_count": cp_solver.num_conflicts,
                "independent_solver": "Z3",
                "independent_solver_version": z3.get_version_string(),
                "independent_terminal_status": str(z3_status),
                "two_solver_exact_count_agreement": True,
                "fixed_source_compact_incidence_no_go": True,
            }
        )

    document: dict[str, object] = {
        "schema": "coxeter.h9_14.record0_a2.incidence_tranche.v1",
        "kind": "h9_14_record0_a2_compact_incidence_no_go_tranche",
        "claim": (
            "Eight of the 36 exact Record-0/a=2 source orbits admit no "
            "compact Coxeter nerve incidence completion."
        ),
        "record_index": 0,
        "graph6": "M?????????????o??",
        "schedule": {
            "a": 2,
            "scheduled_facet_count": 176,
            "scheduled_missing_independent_five_count": 38,
        },
        "source_reduction": reduction,
        "excluded_source_orbit_indices": list(EXCLUDED_SOURCE_ORBITS),
        "remaining_source_orbit_indices": [
            source for source in range(36)
            if source not in EXCLUDED_SOURCE_ORBITS
        ],
        "evidence": evidence,
        "coverage": {
            "exact_source_orbit_count": 36,
            "two_solver_closed_source_orbit_count": len(evidence),
            "remaining_source_orbit_count": 36 - len(evidence),
            "whole_record0_a2_schedule_excluded": False,
            "whole_record0_compact_shadow_excluded": False,
        },
        "proof": {
            "incidence_model": (
                "Select 176 complements of rank-nine facets; witness every "
                "required face through order four; enforce exactly 38 "
                "independent missing five-sets, closure of every larger "
                "nonface by a missing set of order at most five, and ridge "
                "degree zero or two."
            ),
            "symmetry_slice": (
                "Adjacent swaps generate each isolate-column symmetric "
                "kernel; every valid active-swap/A4-row-swap quotient lift "
                "is included. The global lexicographic minimum of each "
                "finite source orbit satisfies every generator inequality."
            ),
            "compact_consequence": (
                "Every compact hyperbolic Coxeter polytope has a simplicial "
                "nerve whose minimal nonfaces are classical Lanner diagrams "
                "of order at most five, so infeasibility of this overcover "
                "excludes the fixed source."
            ),
        },
        "theorem": {
            "excluded_compact_source_orbits": list(EXCLUDED_SOURCE_ORBITS),
            "whole_record0_a2_schedule_excluded": False,
            "whole_record0_compact_shadow_excluded": False,
            "finite_volume_noncompact_checked": False,
            "compact_label_or_gram_survivor_produced": False,
        },
        "nonclaims": [
            "The other 28 Record-0/a=2 source orbits are not excluded by this certificate.",
            "No other Record-0 compact schedule is affected.",
            "No finite-volume noncompact chamber, ideal vertex, or Vinberg cusp is tested.",
            "No hyperbolic Coxeter polytope is constructed.",
        ],
        "source_files": {
            str(path.relative_to(ROOT)): _file_sha256(path)
            for path in SOURCE_PATHS
        },
    }
    return normalize_document(document)


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("kind") != "h9_14_record0_a2_compact_incidence_no_go_tranche":
        raise ValueError("wrong certificate kind")
    sources = document.get("source_files")
    if not isinstance(sources, Mapping):
        raise ValueError("source hash map is malformed")
    for path in SOURCE_PATHS:
        if sources.get(str(path.relative_to(ROOT))) != _file_sha256(path):
            raise ValueError(f"source file changed: {path}")
    digest = document.get("mathematical_digest")
    observed = _mathematical_digest(document)
    if digest != observed:
        raise ValueError("mathematical digest mismatch")
    if EXPECTED_MATHEMATICAL_DIGEST and digest != EXPECTED_MATHEMATICAL_DIGEST:
        raise ValueError("certificate is not the frozen production record")
    reduction = document.get("source_reduction")
    coverage = document.get("coverage")
    theorem = document.get("theorem")
    evidence = document.get("evidence")
    if not isinstance(reduction, Mapping) or not isinstance(coverage, Mapping):
        raise ValueError("coverage is malformed")
    if not isinstance(theorem, Mapping) or not isinstance(evidence, list):
        raise ValueError("theorem evidence is malformed")
    if any(reduction.get(key) != value for key, value in EXPECTED_SOURCE_REDUCTION.items()):
        raise ValueError("source reduction changed")
    if reduction.get("representative_digest") != EXPECTED_SOURCE_REPRESENTATIVE_DIGEST:
        raise ValueError("source representative digest changed")
    if document.get("excluded_source_orbit_indices") != list(EXCLUDED_SOURCE_ORBITS):
        raise ValueError("excluded tranche changed")
    if len(evidence) != len(EXCLUDED_SOURCE_ORBITS):
        raise ValueError("evidence coverage changed")
    for expected_source, row in zip(EXCLUDED_SOURCE_ORBITS, evidence, strict=True):
        if (
            row.get("source_orbit_index") != expected_source
            or row.get("production_terminal_status") != "INFEASIBLE"
            or row.get("independent_terminal_status") != "unsat"
            or row.get("two_solver_exact_count_agreement") is not True
            or row.get("fixed_source_compact_incidence_no_go") is not True
        ):
            raise ValueError(f"source {expected_source}: evidence changed")
    metadata = document.get("run_metadata")
    if not isinstance(metadata, Mapping) or metadata.get("canonical") is not False:
        raise ValueError("noncanonical run metadata is malformed")
    metadata_rows = metadata.get("per_source")
    if not isinstance(metadata_rows, list) or len(metadata_rows) != len(evidence):
        raise ValueError("run metadata coverage changed")
    for expected_source, row in zip(EXCLUDED_SOURCE_ORBITS, metadata_rows, strict=True):
        if not isinstance(row, Mapping) or row.get("source_orbit_index") != expected_source:
            raise ValueError("run metadata source order changed")
        for key in (
            "production_wall_time_seconds",
            "production_branch_count",
            "production_conflict_count",
        ):
            value = row.get(key)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"source {expected_source}: invalid statistics")
            if not isfinite(value) or value < 0:
                raise ValueError(f"source {expected_source}: invalid statistics")
    ledger = document.get("scout_status_ledger")
    if not isinstance(ledger, Mapping) or (
        ledger.get("two_solver_unsat_source_orbits") != list(EXCLUDED_SOURCE_ORBITS)
        or ledger.get("explicit_incidence_sat_source_orbits")
        != list(INCIDENCE_SAT_SOURCE_ORBITS)
        or ledger.get("solver_open_source_orbits")
        != list(SOLVER_OPEN_SOURCE_ORBITS)
    ):
        raise ValueError("scout status partition changed")
    if (
        coverage.get("two_solver_closed_source_orbit_count") != 8
        or coverage.get("remaining_source_orbit_count") != 28
        or coverage.get("whole_record0_a2_schedule_excluded") is not False
        or theorem.get("whole_record0_a2_schedule_excluded") is not False
        or theorem.get("finite_volume_noncompact_checked") is not False
    ):
        raise ValueError("theorem scope changed")
    return {
        "mathematical_digest": digest,
        "exact_source_orbit_count": 36,
        "excluded_source_orbit_indices": list(EXCLUDED_SOURCE_ORBITS),
        "two_solver_closed_source_orbit_count": 8,
        "remaining_source_orbit_count": 28,
        "whole_record0_a2_schedule_excluded": False,
        "whole_record0_compact_shadow_excluded": False,
        "finite_volume_noncompact_checked": False,
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresh", action="store_true")
    parser.add_argument("--write", type=Path, default=OUTPUT)
    parser.add_argument("--normalize-archive", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args(argv)
    if sum((args.fresh, args.normalize_archive, args.check_archive)) > 1:
        parser.error("choose at most one archive mode")
    if args.fresh:
        document = fresh_record()
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        document = json.loads(args.write.read_text(encoding="utf-8"))
        if args.normalize_archive:
            document = normalize_document(document)
            args.write.write_text(
                json.dumps(document, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        elif args.check_archive and normalize_document(document) != document:
            raise ValueError("archive is not in normalized canonical form")
    print(json.dumps(verify(document), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
