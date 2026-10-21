#!/usr/bin/env python3
"""Run Vinberg's literal finite-volume gate on the dual H11/29 sections."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any

import sympy as sp

from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    _parse_output,
)


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT / "data/enumerations/alexeev_118_dual_h11_29_construction.json"
)
CPP_SOURCE = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
OUTPUT = ROOT / "data/enumerations/alexeev_118_dual_h11_29_vinberg.json"
SCHEMA = "coxeter.alexeev_118_dual_h11_29.vinberg.v1"


EXPECTED = {
    identifier: {
        "dimension": 11,
        "facets": 29,
        "finite_special_counts": [
            29,
            375,
            2744,
            12668,
            39131,
            83290,
            122888,
            123297,
            79509,
            28798,
            3876,
        ],
        "affine_components": 668,
        "cusps": 253,
        "endpoint_histogram": {"0,2": 1816, "1,1": 11328, "2,0": 15654},
    }
    for identifier in ("A16_6_1_C01", "A16_12_1_C01")
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _input_stream(case: dict[str, Any]) -> str:
    dimension = int(case["case"]["hyperbolic_dimension"])
    gram = case["restricted_mirrors"]["positive_wall_gram"]
    rows = [f"{dimension} {len(gram)}"]
    for row in gram:
        integers = []
        for value in row:
            rational = sp.Rational(value)
            if rational.q != 1:
                raise ArithmeticError("Vinberg engine expected an integral wall Gram")
            integers.append(str(int(rational)))
        rows.append(" ".join(integers))
    return "\n".join(rows) + "\n"


def _audit_case(
    executable: Path, case: dict[str, Any], expected: dict[str, Any]
) -> dict[str, Any]:
    completed = subprocess.run(
        [str(executable)],
        input=_input_stream(case),
        text=True,
        capture_output=True,
        check=True,
    )
    record = _parse_output(completed.stdout)
    special = record["finite_special_subset_counts_by_rank"]
    endpoints = record["edge_endpoint_type_histogram"]
    elliptic_edges = special[-2]
    finite_vertices = special[-1]
    checks = {
        "dimension_and_facets": (
            record["dimension"] == expected["dimension"]
            and record["facet_count"] == expected["facets"]
        ),
        "all_finite_special_subsets_enumerated": (
            special == expected["finite_special_counts"]
            and record["all_finite_subset_search_node_count"] == 1 + sum(special)
        ),
        "all_connected_affine_components_enumerated": (
            record["connected_affine_component_count"]
            == expected["affine_components"]
        ),
        "all_affine_zero_neighborhoods_are_full_rank_parabolic": (
            record["bad_parabolic_neighborhood_count"] == 0
            and record["ideal_vertex_count"] == expected["cusps"]
        ),
        "every_elliptic_edge_has_exactly_two_endpoints": (
            record["bad_edge_count"] == 0
            and endpoints == expected["endpoint_histogram"]
            and sum(endpoints.values()) == elliptic_edges
            and sum(
                (int(kind.split(",")[0]) + int(kind.split(",")[1])) * count
                for kind, count in endpoints.items()
            )
            == 2 * elliptic_edges
        ),
        "finite_vertex_edge_incidence": (
            2 * endpoints.get("2,0", 0) + endpoints.get("1,1", 0)
            == record["dimension"] * finite_vertices
        ),
        "noncompact": record["ideal_vertex_count"] > 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Vinberg replay changed for {case['id']}: {checks}")
    return record | {
        "id": case["id"],
        "elliptic_edge_count": elliptic_edges,
        "finite_vertex_count": finite_vertices,
        "criterion": (
            "Vinberg edge-endpoint criterion plus every irreducible affine "
            "component's full zero-neighborhood"
        ),
        "finite_volume": True,
        "compact": False,
        "classification": "finite_volume_noncompact_Coxeter_polytope",
        "checks": checks,
        "all_checks_pass": True,
    }


def build_certificate(
    construction_path: Path = CONSTRUCTION,
    cpp_source: Path = CPP_SOURCE,
    compiler: str | None = None,
) -> dict[str, Any]:
    construction = json.loads(construction_path.read_text(encoding="utf-8"))
    cases = {case["id"]: case for case in construction["cases"]}
    if set(cases) != set(EXPECTED) or construction["all_checks_pass"] is not True:
        raise ArithmeticError("the pinned dual H11/29 construction changed")
    selected = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected is None:
        raise RuntimeError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="alexeev-h11-29-vinberg-") as raw:
        executable = Path(raw) / "vinberg"
        subprocess.run(
            [
                selected,
                "-O3",
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(cpp_source),
                "-o",
                str(executable),
            ],
            check=True,
        )
        audits = [
            _audit_case(executable, cases[identifier], EXPECTED[identifier])
            for identifier in EXPECTED
        ]
    invariant_signatures = {
        (
            tuple(row["finite_special_subset_counts_by_rank"]),
            row["connected_affine_component_count"],
            row["ideal_vertex_count"],
            tuple(sorted(row["edge_endpoint_type_histogram"].items())),
        )
        for row in audits
    }
    if len(invariant_signatures) != 1:
        raise ArithmeticError("the two dual volume invariants no longer agree")
    return {
        "schema": SCHEMA,
        "source": {
            "construction_path": str(construction_path.relative_to(ROOT)),
            "construction_sha256": _sha256(construction_path),
            "cpp_source_path": str(cpp_source.relative_to(ROOT)),
            "cpp_source_sha256": _sha256(cpp_source),
            "compiler": Path(selected).name,
            "arithmetic": (
                "integer finite-Coxeter recognition and overflow-guarded "
                "signed 128-bit Bareiss determinants"
            ),
        },
        "method": {
            "finite_enumeration": (
                "all finite Coxeter special subsets, retaining every rank-10 "
                "edge and rank-11 finite endpoint"
            ),
            "affine_completeness": (
                "all irreducible affine subsets plus parallel affine A1 pairs; "
                "each full zero-neighborhood checked at parabolic rank 10"
            ),
            "volume_criterion": (
                "every elliptic edge has exactly two finite or ideal endpoints"
            ),
        },
        "audits": audits,
        "summary": {
            "audit_count": 2,
            "finite_volume_count": 2,
            "compact_count": 0,
            "noncompact_count": 2,
            "shared_elliptic_edge_count": 28798,
            "shared_finite_vertex_count": 3876,
            "shared_cusp_count": 253,
            "shared_bad_edge_count": 0,
            "shared_bad_parabolic_neighborhood_count": 0,
            "all_invariant_vinberg_counts_agree": True,
        },
        "all_checks_pass": all(row["all_checks_pass"] for row in audits),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--construction", type=Path, default=CONSTRUCTION)
    parser.add_argument("--cpp-source", type=Path, default=CPP_SOURCE)
    parser.add_argument("--compiler")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.construction, args.cpp_source, args.compiler)
    if args.check_archive:
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt H11/29 Vinberg archive changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
