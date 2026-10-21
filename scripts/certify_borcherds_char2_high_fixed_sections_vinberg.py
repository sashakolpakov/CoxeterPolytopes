#!/usr/bin/env python3
"""Run the literal Vinberg edge/cusp audit for H13/50 and H14/98."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT
    / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_construction.json"
)
CPP_SOURCE = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
OUTPUT = (
    ROOT
    / "data/enumerations/borcherds_char2_2a_h13_50_2c_h14_98_vinberg.json"
)
SCHEMA = "coxeter.borcherds_char2_2a_h13_50_2c_h14_98.vinberg.v1"


EXPECTED = {
    "SS2-2A-H13-50": {
        "dimension": 13,
        "facets": 50,
        "finite_special_counts": [
            50,
            949,
            9528,
            59446,
            250786,
            749748,
            1627172,
            2581459,
            2965280,
            2393834,
            1275684,
            391614,
            45954,
        ],
        "affine_components": 6210,
        "cusps": 2598,
        "endpoint_histogram": {"0,2": 23340, "1,1": 139146, "2,0": 229128},
    },
    "SS2-2C-H14-98": {
        "dimension": 14,
        "facets": 98,
        "finite_special_counts": [
            98,
            1666,
            14616,
            83559,
            342888,
            1056573,
            2490674,
            4501714,
            6177500,
            6299076,
            4598006,
            2249968,
            646786,
            75898,
        ],
        "affine_components": 5838,
        "cusps": 2114,
        "endpoint_histogram": {"0,2": 19628, "1,1": 191744, "2,0": 435414},
    },
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_output(raw: str) -> dict[str, Any]:
    lines = [line.strip() for line in raw.splitlines() if line.strip()]
    if len(lines) != 5:
        raise ArithmeticError(f"unexpected C++ output: {raw}")
    header = re.fullmatch(
        r"dimension (\d+) facets (\d+) finite_nodes (\d+) affine_nodes (\d+)",
        lines[0],
    )
    prefix = re.fullmatch(r"edge_search_prefix_counts((?: \d+)+)", lines[1])
    special = re.fullmatch(r"finite_special_counts((?: \d+)+)", lines[2])
    affine = re.fullmatch(
        r"affine_components (\d+) cusps (\d+) bad_neighborhoods (\d+) "
        r"bad_edges (\d+)",
        lines[3],
    )
    endpoint = re.fullmatch(r"endpoint_hist((?: f\di\d \d+)+)", lines[4])
    if not all((header, prefix, special, affine, endpoint)):
        raise ArithmeticError(f"could not parse C++ output: {raw}")
    endpoint_histogram = {
        f"{finite},{ideal}": int(count)
        for finite, ideal, count in re.findall(r"f(\d)i(\d) (\d+)", lines[4])
    }
    return {
        "dimension": int(header.group(1)),
        "facet_count": int(header.group(2)),
        "edge_search_node_count": int(header.group(3)),
        "all_finite_subset_search_node_count": int(header.group(4)),
        "ordering_dependent_edge_search_prefix_counts": [
            int(value) for value in prefix.group(1).split()
        ],
        "finite_special_subset_counts_by_rank": [
            int(value) for value in special.group(1).split()
        ],
        "connected_affine_component_count": int(affine.group(1)),
        "ideal_vertex_count": int(affine.group(2)),
        "bad_parabolic_neighborhood_count": int(affine.group(3)),
        "bad_edge_count": int(affine.group(4)),
        "edge_endpoint_type_histogram": endpoint_histogram,
    }


def _input_stream(case: dict[str, Any]) -> str:
    dimension = int(case["decision"]["hyperbolic_dimension"])
    negative_gram = case["restricted_chamber"]["facet_gram"]
    facets = len(negative_gram)
    rows = [f"{dimension} {facets}"]
    rows.extend(
        " ".join(str(-int(value)) for value in row) for row in negative_gram
    )
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
    edge_count = special[-2]
    finite_vertices = special[-1]
    checks = {
        "dimension_and_facet_count": (
            record["dimension"] == expected["dimension"]
            and record["facet_count"] == expected["facets"]
        ),
        "all_finite_special_subsets_enumerated": (
            special == expected["finite_special_counts"]
            and record["all_finite_subset_search_node_count"]
            == 1 + sum(special)
        ),
        "elliptic_edges_and_finite_vertices": (
            edge_count == expected["finite_special_counts"][-2]
            and finite_vertices == expected["finite_special_counts"][-1]
        ),
        "all_connected_affine_components_enumerated": (
            record["connected_affine_component_count"]
            == expected["affine_components"]
        ),
        "all_full_zero_neighborhoods_parabolic": (
            record["bad_parabolic_neighborhood_count"] == 0
            and record["ideal_vertex_count"] == expected["cusps"]
        ),
        "every_edge_has_exactly_two_endpoints": (
            record["bad_edge_count"] == 0
            and endpoints == expected["endpoint_histogram"]
            and sum(endpoints.values()) == edge_count
            and sum(
                (int(kind.split(",")[0]) + int(kind.split(",")[1])) * count
                for kind, count in endpoints.items()
            )
            == 2 * edge_count
        ),
        "finite_vertex_edge_incidence": (
            2 * endpoints.get("2,0", 0) + endpoints.get("1,1", 0)
            == record["dimension"] * finite_vertices
        ),
        "cusps_exist": record["ideal_vertex_count"] > 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Vinberg census changed for {case['id']}: {checks}")
    record |= {
        "id": case["id"],
        "elliptic_edge_count": edge_count,
        "finite_vertex_count": finite_vertices,
        "criterion": (
            "Vinberg edge-endpoint criterion plus every irreducible affine "
            "component's full zero-neighborhood"
        ),
        "finite_volume": True,
        "compact": False,
        "classification": "finite_volume_noncompact_Coxeter_polytope",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    return record


def build_certificate(
    construction_path: Path = CONSTRUCTION,
    cpp_source: Path = CPP_SOURCE,
    compiler: str | None = None,
) -> dict[str, Any]:
    construction = json.loads(construction_path.read_text(encoding="utf-8"))
    cases = {case["id"]: case for case in construction["cases"]}
    if set(cases) != set(EXPECTED) or not construction["all_checks_pass"]:
        raise ArithmeticError("the pinned H13/H14 construction archive changed")
    selected_compiler = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected_compiler is None:
        raise RuntimeError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="borcherds-char2-vinberg-") as raw:
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
                "integer Coxeter classification and overflow-guarded signed "
                "128-bit Bareiss determinants"
            ),
        },
        "method": {
            "finite_enumeration": (
                "all finite Coxeter special subsets, with every rank n-1 "
                "edge and rank n finite endpoint retained"
            ),
            "affine_completeness": (
                "every irreducible affine set is reached as the first "
                "nonfinite extension of an elliptic proper subset; parallel "
                "affine A1 pairs are scanned separately"
            ),
            "volume_gate": (
                "all full affine zero-neighborhoods have parabolic rank n-1 "
                "and every elliptic edge has exactly two finite/ideal endpoints"
            ),
        },
        "cases": audits,
        "all_checks_pass": all(record["all_checks_pass"] for record in audits),
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
            raise ArithmeticError("rebuilt H13/H14 Vinberg certificate changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                record["id"]: {
                    "edges": record["elliptic_edge_count"],
                    "finite_vertices": record["finite_vertex_count"],
                    "cusps": record["ideal_vertex_count"],
                    "bad_edges": record["bad_edge_count"],
                }
                for record in payload["cases"]
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
