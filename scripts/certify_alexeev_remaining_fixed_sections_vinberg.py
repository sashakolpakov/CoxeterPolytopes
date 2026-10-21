#!/usr/bin/env python3
"""Run Vinberg's literal finite-volume gate on six Alexeev sections.

The companion construction archive supplies exact positive-normal Gram
matrices.  The shared C++ engine enumerates every elliptic rank ``n-1`` edge,
every elliptic rank ``n`` finite endpoint, every irreducible affine special
subdiagram, and each affine component's full zero-neighborhood.  Thus no
volume conclusion is inferred from facet counts or from a partial face list.
"""

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
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_remaining_high_fixed_section_construction.json"
)
CPP_SOURCE = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_remaining_high_fixed_section_vinberg.json"
)
SCHEMA = "coxeter.alexeev_remaining_high_fixed_sections.vinberg.v1"


EXPECTED = {
    "A18_4_0_C03": {
        "dimension": 11,
        "facets": 28,
        "finite_special_counts": [
            28, 292, 1698, 6405, 16794, 31498, 42316, 39764, 24660, 8851, 1252
        ],
        "affine_components": 186,
        "cusps": 64,
        "endpoint_histogram": {"0,2": 405, "1,1": 3120, "2,0": 5326},
    },
    "A18_4_0_C04": {
        "dimension": 10,
        "facets": 28,
        "finite_special_counts": [
            28, 295, 1687, 5993, 14091, 22529, 24428, 17212, 6956, 1087
        ],
        "affine_components": 204,
        "cusps": 74,
        "endpoint_histogram": {"0,2": 282, "1,1": 2478, "2,0": 4196},
    },
    "A18_4_0_C06": {
        "dimension": 13,
        "facets": 31,
        "finite_special_counts": [
            31,
            399,
            2925,
            13964,
            46543,
            112606,
            201943,
            270518,
            269179,
            194317,
            96280,
            28766,
            3433,
        ],
        "affine_components": 529,
        "cusps": 238,
        "endpoint_histogram": {"0,2": 1681, "1,1": 9541, "2,0": 17544},
    },
    "A18_4_0_C18": {
        "dimension": 10,
        "facets": 21,
        "finite_special_counts": [
            21, 173, 806, 2443, 5125, 7581, 7792, 5284, 2088, 329
        ],
        "affine_components": 57,
        "cusps": 22,
        "endpoint_histogram": {"0,2": 76, "1,1": 734, "2,0": 1278},
    },
    "A15_11_1_C02": {
        "dimension": 10,
        "facets": 21,
        "finite_special_counts": [
            21, 194, 1042, 3617, 8503, 13728, 15006, 10535, 4159, 607
        ],
        "affine_components": 132,
        "cusps": 50,
        "endpoint_histogram": {"0,2": 276, "1,1": 1696, "2,0": 2187},
    },
    "A15_7_1_C01": {
        "dimension": 10,
        "facets": 27,
        "finite_special_counts": [
            27, 264, 1435, 5034, 12061, 19934, 22234, 15784, 6221, 895
        ],
        "affine_components": 198,
        "cusps": 67,
        "endpoint_histogram": {"0,2": 432, "1,1": 2628, "2,0": 3161},
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _input_stream(case: dict[str, Any]) -> str:
    dimension = int(case["case"]["hyperbolic_dimension"])
    gram = case["restricted_mirrors"]["positive_wall_gram"]
    facets = len(gram)
    rows = [f"{dimension} {facets}"]
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
    edges = special[-2]
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
            and sum(endpoints.values()) == edges
            and sum(
                (int(kind.split(",")[0]) + int(kind.split(",")[1])) * count
                for kind, count in endpoints.items()
            )
            == 2 * edges
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
        "elliptic_edge_count": edges,
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
        raise ArithmeticError("the pinned remaining-stock construction changed")
    selected_compiler = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected_compiler is None:
        raise RuntimeError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="alexeev-remaining-vinberg-") as raw:
        executable = Path(raw) / "vinberg"
        subprocess.run(
            [
                selected_compiler,
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
    return {
        "schema": SCHEMA,
        "source": {
            "construction_path": str(construction_path.relative_to(ROOT)),
            "construction_sha256": _sha256(construction_path),
            "cpp_source_path": str(cpp_source.relative_to(ROOT)),
            "cpp_source_sha256": _sha256(cpp_source),
            "compiler": Path(selected_compiler).name,
            "arithmetic": (
                "integer finite-Coxeter recognition and overflow-guarded "
                "signed 128-bit Bareiss determinants"
            ),
        },
        "method": {
            "finite_enumeration": (
                "all finite Coxeter special subsets, retaining every rank "
                "n-1 edge and rank n finite endpoint"
            ),
            "affine_completeness": (
                "every irreducible affine set reached at its first nonfinite "
                "extension; parallel affine A1 pairs scanned separately"
            ),
            "volume_gate": (
                "all full affine zero-neighborhoods have parabolic rank n-1 "
                "and every elliptic edge has two finite/ideal endpoints"
            ),
            "ordering_warning": (
                "ordering_dependent_edge_search_prefix_counts are search "
                "diagnostics, not face numbers; finite_special_subset_counts "
                "are the invariant special-subset census"
            ),
        },
        "cases": audits,
        "summary": {
            "case_count": len(audits),
            "finite_volume_noncompact_count": sum(
                row["classification"]
                == "finite_volume_noncompact_Coxeter_polytope"
                for row in audits
            ),
            "bad_edge_count": sum(row["bad_edge_count"] for row in audits),
            "bad_parabolic_neighborhood_count": sum(
                row["bad_parabolic_neighborhood_count"] for row in audits
            ),
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
            raise ArithmeticError("rebuilt remaining-stock Vinberg certificate changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
