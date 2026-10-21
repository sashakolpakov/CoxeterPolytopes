#!/usr/bin/env python3
"""Certify two unexpanded rank-ten Roulleau full-nef K3 chambers.

Roulleau's ``U+E8(-1)`` and ``U+D8(-1)`` rows each have ten complete
``(-2)`` walls.  Their published curve graphs are respectively the classical
JKRT ``T9`` and ``Q9`` finite-volume Coxeter simplices.  This script checks
the integral wall lattices, ample orientation, every projective vertex ray,
literal graph identifications, and Vinberg's edge/cusp criterion.

The completeness of the ten curve lists is a cited theorem of Roulleau, not
something inferred from the ten-by-ten matrices below.
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
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    _parse_output,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CPP_SOURCE = (
    ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
)
OUTPUT = (
    ROOT
    / "data/enumerations/roulleau_rank10_two_full_nef_simplices.json"
)
SCHEMA = "coxeter.roulleau_rank10_two_full_nef_simplices.v1"

# Hashes of the pinned primary-source snapshots used for the transcription.
ROULLEAU_ARXIV = "https://arxiv.org/abs/2003.08985"
ROULLEAU_SOURCE_TAR_SHA256 = (
    "3e445ed01446b47ded6c007b395308de34e566b736f0ba982470bc497adeb1e9"
)
ROULLEAU_SOURCE_TEX_SHA256 = (
    "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
)
JKRT_DOI = "https://doi.org/10.1007/BF01238563"


CASES: dict[str, dict[str, Any]] = {
    "R10-UE8-NEF-H9-10": {
        "lattice": "U+E8(-1)",
        "absolute_discriminant": 1,
        "edges": (
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (6, 8),
            (8, 9),
        ),
        "ample_coordinates": (30, 61, 93, 126, 160, 195, 231, 115, 153, 76),
        "ample_square": 1240,
        "fiber_coordinates": (0, 1, 2, 3, 4, 5, 6, 3, 4, 2),
        "section_wall": 1,
        "branch_model": (
            "double cover of F4 branched over the negative section and "
            "B in |3s+12f| with an E8 singularity"
        ),
        "reference": {
            "name": "JKRT T9",
            "coxiter_file": "9-jkrt_351-01-T9.coxiter",
            "snapshot_sha256": (
                "ae7ceb0ef167f0f59370a559cd36338ed8206e700631c8c26e5be93ba0dac583"
            ),
            "edges": (
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 6),
                (6, 7),
                (6, 8),
                (8, 9),
            ),
            "candidate_to_reference_permutation": tuple(range(10)),
        },
        "expected": {
            "determinant": -1,
            "smith": (1,) * 10,
            "vertex_norms": (0, 2, 6, 12, 20, 30, 42, 10, 18, 4),
            "finite_vertices": 9,
            "ideal_vertices": 1,
            "affine_components": 1,
            "endpoint_histogram": {"1,1": 9, "2,0": 36},
            "finite_search_nodes": 384,
            "all_subset_search_nodes": 1022,
        },
    },
    "R10-UD8-NEF-H9-10": {
        "lattice": "U+D8(-1)",
        "absolute_discriminant": 4,
        "edges": (
            (0, 2),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 7),
            (6, 8),
            (8, 9),
        ),
        "ample_coordinates": (15, 15, 31, 33, 36, 40, 45, 22, 29, 14),
        "ample_square": 280,
        "fiber_coordinates": (1, 1, 2, 2, 2, 2, 2, 1, 1, 0),
        "section_wall": 10,
        "branch_model": (
            "double cover of F4 branched over the negative section and "
            "B in |3s+12f| with a D8 singularity"
        ),
        "reference": {
            "name": "JKRT Q9",
            "coxiter_file": "9-jkrt_351-03-Q9.coxiter",
            "snapshot_sha256": (
                "0940d126d99987c236ee31a06c12f37c04ab786c0ab8b089ab2a9b941988b1bc"
            ),
            "edges": (
                (0, 1),
                (1, 2),
                (1, 3),
                (3, 4),
                (4, 5),
                (5, 6),
                (6, 7),
                (6, 8),
                (8, 9),
            ),
            "candidate_to_reference_permutation": (0, 2, 1, 3, 4, 5, 6, 7, 8, 9),
        },
        "expected": {
            "determinant": -4,
            "smith": (1, 1, 1, 1, 1, 1, 1, 1, 2, 2),
            "vertex_norms": (0, 0, 2, 3, 4, 5, 6, 1, 2, 0),
            "finite_vertices": 7,
            "ideal_vertices": 3,
            "affine_components": 3,
            "endpoint_histogram": {"0,2": 3, "1,1": 21, "2,0": 21},
            "finite_search_nodes": 382,
            "all_subset_search_nodes": 1020,
        },
    },
}


def _wall_gram(edges: tuple[tuple[int, int], ...]) -> sp.Matrix:
    gram = -2 * sp.eye(10)
    for left, right in edges:
        gram[left, right] = gram[right, left] = 1
    return gram


def _payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _canonical_sha256(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _input_stream(gram: sp.Matrix) -> str:
    rows = ["9 10"]
    rows.extend(
        " ".join(str(-int(gram[row, column])) for column in range(10))
        for row in range(10)
    )
    return "\n".join(rows) + "\n"


def _mapped_edges(
    edges: tuple[tuple[int, int], ...], permutation: tuple[int, ...]
) -> set[tuple[int, int]]:
    return {
        tuple(sorted((permutation[left], permutation[right])))
        for left, right in edges
    }


def _matrix_certificate(case: dict[str, Any]) -> dict[str, Any]:
    gram = _wall_gram(case["edges"])
    expected = case["expected"]
    ample = sp.Matrix(case["ample_coordinates"])
    ones = sp.ones(10, 1)
    inverse = gram.inv()
    vertex_norms = tuple(int(inverse[index, index]) for index in range(10))

    # The i-th column of G^{-1} is the ray opposite wall i.  The sum of all
    # other dual rays is a strict relative-interior witness for facet i.
    dual_rays = [inverse[:, index] for index in range(10)]
    facet_witnesses = []
    for omitted in range(10):
        witness = sum(
            (dual_rays[index] for index in range(10) if index != omitted),
            sp.zeros(10, 1),
        )
        pairings = tuple(int(value) for value in gram * witness)
        facet_witnesses.append(
            {
                "wall": omitted + 1,
                "square": int((witness.T * gram * witness)[0]),
                "pairings": pairings,
                "same_sheet_pairing_with_ample": int(
                    (witness.T * gram * ample)[0]
                ),
            }
        )

    smith = tuple(
        abs(int(value))
        for value in smith_normal_form(gram, domain=ZZ).diagonal()
    )
    mapped = _mapped_edges(
        case["edges"], case["reference"]["candidate_to_reference_permutation"]
    )
    reference_edges = {tuple(sorted(edge)) for edge in case["reference"]["edges"]}
    checks = {
        "ten_primitive_minus_two_walls": all(gram[index, index] == -2 for index in range(10)),
        "coxeter_angles_are_only_pi_over_2_or_pi_over_3": all(
            gram[row, column] in (0, 1)
            for row in range(10)
            for column in range(row + 1, 10)
        ),
        "rank_and_inertia": gram.rank() == 10 and exact_inertia(gram) == (1, 9, 0),
        "determinant": int(gram.det()) == expected["determinant"],
        "absolute_discriminant": abs(int(gram.det())) == case["absolute_discriminant"],
        "smith_normal_form": smith == expected["smith"],
        "published_ample_pairings": gram * ample == ones,
        "published_ample_square": int((ample.T * gram * ample)[0]) == case["ample_square"],
        "all_vertex_rays_finite_or_ideal": all(value >= 0 for value in vertex_norms),
        "vertex_norms": vertex_norms == expected["vertex_norms"],
        "finite_vertex_count": sum(value > 0 for value in vertex_norms)
        == expected["finite_vertices"],
        "ideal_vertex_count": sum(value == 0 for value in vertex_norms)
        == expected["ideal_vertices"],
        "every_wall_has_strict_facet_witness": all(
            witness["square"] > 0
            and witness["same_sheet_pairing_with_ample"] > 0
            and witness["pairings"][witness["wall"] - 1] == 0
            and all(
                pairing == 1
                for index, pairing in enumerate(witness["pairings"])
                if index != witness["wall"] - 1
            )
            for witness in facet_witnesses
        ),
        "literal_jkrt_graph_isomorphism": mapped == reference_edges,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ArithmeticError(f"matrix checks failed: {failed}")
    payload = {
        "wall_gram": _payload(gram),
        "wall_gram_sha256": _canonical_sha256(_payload(gram)),
        "determinant": int(gram.det()),
        "smith_invariants": smith,
        "inertia": exact_inertia(gram),
        "ample_coordinates_in_wall_basis": case["ample_coordinates"],
        "ample_wall_pairings": tuple(int(value) for value in gram * ample),
        "ample_square": int((ample.T * gram * ample)[0]),
        "opposite_vertex_norms": vertex_norms,
        "facet_witnesses": facet_witnesses,
        "reference": case["reference"],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    # Keep the in-memory replay byte-for-byte comparable with the JSON
    # archive: tuples in the source constants become JSON arrays here.
    return json.loads(json.dumps(payload))


def _vinberg_certificate(
    executable: Path, case: dict[str, Any], gram: sp.Matrix
) -> dict[str, Any]:
    raw_input = _input_stream(gram)
    completed = subprocess.run(
        [str(executable)],
        input=raw_input,
        text=True,
        capture_output=True,
        check=True,
    )
    record = _parse_output(completed.stdout)
    expected = case["expected"]
    special = record["finite_special_subset_counts_by_rank"]
    endpoint_histogram = record["edge_endpoint_type_histogram"]
    checks = {
        "dimension_and_facets": record["dimension"] == 9 and record["facet_count"] == 10,
        "all_finite_subset_search": record["all_finite_subset_search_node_count"]
        == expected["all_subset_search_nodes"],
        "edge_search": record["edge_search_node_count"]
        == expected["finite_search_nodes"],
        "all_45_elliptic_edges": special[-2] == 45,
        "finite_vertices": special[-1] == expected["finite_vertices"],
        "ideal_vertices": record["ideal_vertex_count"] == expected["ideal_vertices"],
        "affine_components": record["connected_affine_component_count"]
        == expected["affine_components"],
        "endpoint_histogram": endpoint_histogram == expected["endpoint_histogram"],
        "every_edge_has_two_endpoints": record["bad_edge_count"] == 0,
        "every_full_parabolic_neighborhood_passes": record[
            "bad_parabolic_neighborhood_count"
        ]
        == 0,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ArithmeticError(f"Vinberg checks failed: {failed}")
    payload = {
        **record,
        "elliptic_edge_count": special[-2],
        "finite_vertex_count": special[-1],
        "input_sha256": hashlib.sha256(raw_input.encode()).hexdigest(),
        "stdout_sha256": hashlib.sha256(completed.stdout.encode()).hexdigest(),
        "criterion": (
            "every elliptic rank-8 edge has exactly two finite/ideal endpoints "
            "and every connected affine component has a full parabolic "
            "zero-neighborhood of rank 8"
        ),
        "decision": "finite_volume_noncompact",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    return json.loads(json.dumps(payload))


def build_certificate(
    cpp_source: Path = CPP_SOURCE, compiler: str | None = None
) -> dict[str, Any]:
    selected_compiler = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected_compiler is None:
        raise RuntimeError("a C++17 compiler is required")
    records = []
    with tempfile.TemporaryDirectory(prefix="roulleau-r10-simplices-") as raw:
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
        for identifier, case in CASES.items():
            gram = _wall_gram(case["edges"])
            matrix = _matrix_certificate(case)
            vinberg = _vinberg_certificate(executable, case, gram)
            records.append(
                {
                    "id": identifier,
                    "picard_rank": 10,
                    "hyperbolic_dimension": 9,
                    "facet_count": 10,
                    "lattice": case["lattice"],
                    "complete_curve_source": (
                        "Roulleau, rank-10 subsection and final curve-count table"
                    ),
                    "complete_full_nef_reason": (
                        "the published ten curves are the complete irreducible "
                        "(-2)-curve list; the K3 ample-cone theorem therefore "
                        "identifies their chamber with the entire nef cone"
                    ),
                    "branch_model": case["branch_model"],
                    "elliptic_fiber_coordinates": case["fiber_coordinates"],
                    "section_wall": case["section_wall"],
                    "matrix_certificate": matrix,
                    "literal_vinberg": vinberg,
                    "classification": "finite_volume_noncompact_full_nef_K3_Coxeter_simplex",
                    "novelty": (
                        "known JKRT abstract simplex and published Roulleau K3 "
                        "surface/curve configuration; new explicit project-ledger "
                        "alignment and literal certificate only"
                    ),
                    "all_checks_pass": matrix["all_checks_pass"]
                    and vinberg["all_checks_pass"],
                }
            )

    checks = {
        "two_newly_expanded_rank10_rows": len(records) == 2,
        "both_complete_full_nef": all(
            row["complete_full_nef_reason"] and row["all_checks_pass"]
            for row in records
        ),
        "both_literal_finite_volume_noncompact": all(
            row["literal_vinberg"]["decision"] == "finite_volume_noncompact"
            for row in records
        ),
        "distinct_weighted_diagrams": records[0]["matrix_certificate"][
            "wall_gram_sha256"
        ]
        != records[1]["matrix_certificate"]["wall_gram_sha256"],
        "rank10_denominator_is_two_of_nine": len(records) == 2 and 9 - len(records) == 7,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global checks failed: {checks}")
    payload = {
        "schema": SCHEMA,
        "artifact_id": "roulleau_rank10_two_full_nef_simplices",
        "status_date": "2026-08-09",
        "source": {
            "roulleau_arxiv": ROULLEAU_ARXIV,
            "roulleau_source_tar_sha256": ROULLEAU_SOURCE_TAR_SHA256,
            "roulleau_source_tex_sha256": ROULLEAU_SOURCE_TEX_SHA256,
            "roulleau_tex_lines": {
                "U+E8(-1)": "7105-7170",
                "U+D8(-1)": "7173-7234",
                "rank10_table_rows": "10864-10879",
            },
            "jkrt_doi": JKRT_DOI,
            "cpp_verifier": str(cpp_source.relative_to(ROOT)),
            "cpp_verifier_sha256": hashlib.sha256(cpp_source.read_bytes()).hexdigest(),
            "compiler": Path(selected_compiler).name,
        },
        "ledger_comparison": {
            "global_root2fv_types": 118,
            "published_rank10_types": 9,
            "rows_expanded_here": 2,
            "rank10_rows_left_for_the_next_finite_batch": 7,
            "warning": (
                "The global 118 theorem previously covered these two types only "
                "collectively; T9 and Q9 had no dedicated positive-ledger rows. "
                "This is not a global alias census of every unrelated project construction."
            ),
        },
        "cases": records,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    return json.loads(json.dumps(payload))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler")
    parser.add_argument("--cpp-source", type=Path, default=CPP_SOURCE)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build_certificate(arguments.cpp_source, arguments.compiler)
    if arguments.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if payload != archived:
            raise ArithmeticError("frozen rank-ten simplex archive changed")
    if arguments.write:
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
