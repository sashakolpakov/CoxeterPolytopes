#!/usr/bin/env python3
"""Literal Vinberg audits for all fourteen rank-11--14 full-root K3 chambers.

Thirteen complete curve Gram matrices are loaded from the existing Roulleau
denominator certificates.  The remaining ``U+D4(-1)^2+A1(-1)`` matrix is
reconstructed independently by an exact norm-two Vinberg search on the
standard lattice.  Every resulting complete chamber is then passed through
the generic exact C++ edge/affine-neighbourhood verifier.

The distinction is deliberate: finite automorphism group gives structural
finite volume, whereas this certificate checks Vinberg's literal endpoint
and full zero-neighbourhood criterion on every one of the fourteen diagrams.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any

import sympy as sp

from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    _parse_output,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "k3_rank11_14_literal_vinberg_volume.json"
)
CPP_SOURCE = (
    ROOT / "scripts" / "certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
)
INVENTORY = (
    ROOT / "data" / "enumerations" / "roulleau_rank11_14_source_stock_inventory.json"
)
SCHEMA = "coxeter.k3_rank11_14_literal_vinberg_volume.v1"


@dataclass(frozen=True, slots=True)
class CaseSpec:
    key: str
    rho: int
    lattice: str
    facets: int
    source_path: str | None
    gram_path: tuple[str | int, ...] | None


CASES = (
    CaseSpec(
        "R11_UE8A1",
        11,
        "U+E8(-1)+A1(-1)",
        12,
        "data/enumerations/roulleau_rank11_direct_graph_denominator.json",
        ("cases", 0, "source", "wall_gram"),
    ),
    CaseSpec(
        "R11_UD8A1",
        11,
        "U+D8(-1)+A1(-1)",
        14,
        "data/enumerations/roulleau_rank11_direct_graph_denominator.json",
        ("cases", 1, "source", "wall_gram"),
    ),
    CaseSpec(
        "R11_UD4D4A1",
        11,
        "U+D4(-1)^2+A1(-1)",
        22,
        None,
        None,
    ),
    CaseSpec(
        "R11_UD4A1x5",
        11,
        "U+D4(-1)+A1(-1)^5",
        90,
        "data/enumerations/roulleau_r11_ud4a1x5_90_denominator.json",
        ("walls", "gram"),
    ),
    CaseSpec(
        "R12_UE8A1x2",
        12,
        "U+E8(-1)+A1(-1)^2",
        14,
        "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json",
        ("stocks", 0, "wall_gram"),
    ),
    CaseSpec(
        "R12_UD8A1x2",
        12,
        "U+D8(-1)+A1(-1)^2",
        19,
        "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json",
        ("stocks", 1, "wall_gram"),
    ),
    CaseSpec(
        "R12_UD4x2A1x2",
        12,
        "U+D4(-1)^2+A1(-1)^2",
        59,
        "data/enumerations/roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json",
        ("coordinate_reconstruction", "wall_gram"),
    ),
    CaseSpec(
        "R12_UA2E8",
        12,
        "U+A2(-1)+E8(-1)",
        13,
        "data/enumerations/roulleau_non2elementary_involution_folds.json",
        ("cases", 0, "source", "curve_gram"),
    ),
    CaseSpec(
        "R13_UE8A1x3",
        13,
        "U+E8(-1)+A1(-1)^3",
        17,
        "data/enumerations/roulleau_ue8a1x3_all_subgroup_folds.json",
        ("source", "wall_gram"),
    ),
    CaseSpec(
        "R13_UD8A1x3",
        13,
        "U+D8(-1)+A1(-1)^3",
        39,
        "data/enumerations/roulleau_r13_ud8a1x3_all_subgroup_folds.json",
        (
            "source_reconstruction",
            "authoritative_coordinate_formula_lane",
            "wall_gram",
        ),
    ),
    CaseSpec(
        "R13_UE8A3",
        13,
        "U+E8(-1)+A3(-1)",
        14,
        "data/enumerations/roulleau_non2elementary_involution_folds.json",
        ("cases", 1, "source", "curve_gram"),
    ),
    CaseSpec(
        "R14_UE8D4",
        14,
        "U+E8(-1)+D4(-1)",
        15,
        "data/enumerations/roulleau_ue8d4_all_subgroup_folds.json",
        ("source", "wall_gram"),
    ),
    CaseSpec(
        "R14_UD8D4",
        14,
        "U+D8(-1)+D4(-1)",
        20,
        "data/enumerations/roulleau_ud8d4_all_subgroup_folds.json",
        ("source", "wall_gram"),
    ),
    CaseSpec(
        "R14_UE8A1x4",
        14,
        "U+E8(-1)+A1(-1)^4",
        27,
        "data/enumerations/roulleau_r14_ue8a1x4_all_subgroup_folds.json",
        ("source_reconstruction", "authoritative_lane", "wall_gram"),
    ),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _dig(value: Any, path: tuple[str | int, ...]) -> Any:
    for component in path:
        value = value[component]
    return value


def _standard_ud4d4a1_gram() -> sp.Matrix:
    d4 = sp.Matrix(
        [
            [2, -1, 0, 0],
            [-1, 2, -1, -1],
            [0, -1, 2, 0],
            [0, -1, 0, 2],
        ]
    )
    gram = sp.zeros(11)
    gram[0, 1] = gram[1, 0] = 1
    gram[2:6, 2:6] = -d4
    gram[6:10, 6:10] = -d4
    gram[10, 10] = -2
    return gram


def _pair(gram: sp.Matrix, left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return int((sp.Matrix(left).T * gram * sp.Matrix(right))[0])


def _divisibility(gram: sp.Matrix, root: tuple[int, ...]) -> int:
    result = 0
    for value in gram * sp.Matrix(root):
        result = gcd(result, abs(int(value)))
    return result


def _primitive(root: tuple[int, ...]) -> bool:
    result = 0
    for value in root:
        result = gcd(result, abs(value))
    return result == 1


@lru_cache(maxsize=None)
def _dominant_d4_vectors(max_square: int) -> tuple[tuple[tuple[int, ...], int], ...]:
    """D4 root-lattice vectors in x1>=x2>=x3>=|x4| dominance order."""

    simple_columns = sp.zeros(4)
    for index in range(3):
        simple_columns[index, index] = 1
        simple_columns[index + 1, index] = -1
    simple_columns[2, 3] = simple_columns[3, 3] = 1
    inverse_simple = simple_columns.inv()
    bound = isqrt(max_square)
    vectors: list[tuple[tuple[int, ...], int]] = []
    for x1 in range(bound + 1):
        for x2 in range(x1 + 1):
            for x3 in range(x2 + 1):
                for x4 in range(-x3, x3 + 1):
                    square = x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4
                    if square > max_square or (x1 + x2 + x3 + x4) % 2:
                        continue
                    coordinates = inverse_simple * sp.Matrix((x1, x2, x3, x4))
                    if any(value.q != 1 for value in coordinates):
                        raise ArithmeticError("D4 dominance enumeration left the root lattice")
                    vectors.append(
                        (tuple(int(value) for value in coordinates), square)
                    )
    return tuple(sorted(vectors))


def _ud4d4a1_height_zero_roots() -> tuple[tuple[int, ...], ...]:
    roots: list[tuple[int, ...]] = [(1, -1, *(0,) * 9)]
    for block in (2, 6):
        for simple in range(4):
            root = [0] * 11
            root[block + simple] = -1
            roots.append(tuple(root))
    root = [0] * 11
    root[10] = -1
    roots.append(tuple(root))
    return tuple(roots)


def _ud4d4a1_level_candidates(control: int) -> tuple[tuple[int, ...], ...]:
    """All dominant primitive norm-minus-two candidates at one control level."""

    if control <= 0:
        raise ValueError("control must be positive")
    gram = _standard_ud4d4a1_gram()
    height_zero = _ud4d4a1_height_zero_roots()
    radius = control * control + 4
    vectors_by_square: dict[int, list[tuple[int, ...]]] = {}
    for vector, square in _dominant_d4_vectors(radius // 2):
        vectors_by_square.setdefault(square, []).append(vector)
    candidates: set[tuple[int, ...]] = set()
    for difference in range(-isqrt(radius), 1):
        if (difference - control) % 2:
            continue
        a = (-control + difference) // 2
        b = (-control - difference) // 2
        for z in range(isqrt(radius // 4) + 1):
            remainder = radius - difference * difference - 4 * z * z
            if remainder < 0 or remainder % 2:
                continue
            d4_total = remainder // 2
            for first_square, first_vectors in vectors_by_square.items():
                second_vectors = vectors_by_square.get(d4_total - first_square, ())
                for first in first_vectors:
                    for second in second_vectors:
                        root = (a, b, *first, *second, z)
                        if not _primitive(root):
                            continue
                        if _pair(gram, root, root) != -2:
                            continue
                        if any(_pair(gram, root, wall) < 0 for wall in height_zero):
                            continue
                        candidates.add(root)
    return tuple(sorted(candidates))


def _reconstruct_ud4d4a1() -> tuple[list[list[int]], dict[str, Any]]:
    """Replay the norm-two Vinberg prefix that becomes the 22-wall chamber."""

    gram = _standard_ud4d4a1_gram()
    roots = list(_ud4d4a1_height_zero_roots())
    levels = []
    for control in range(1, 5):
        candidates = _ud4d4a1_level_candidates(control)
        accepted = tuple(
            root
            for root in candidates
            if all(_pair(gram, root, wall) >= 0 for wall in roots)
        )
        if any(
            _pair(gram, left, right) < 0
            for offset, left in enumerate(accepted)
            for right in accepted[offset + 1 :]
        ):
            raise ArithmeticError("same-height norm-two roots are not pairwise acute")
        roots.extend(accepted)
        levels.append(
            {
                "control": control,
                "height": str(Fraction(control * control, 2)),
                "dominant_candidate_count": len(candidates),
                "accepted_root_count": len(accepted),
                "accepted_roots": [list(root) for root in accepted],
            }
        )
    root_matrix = sp.Matrix.hstack(*(sp.Matrix(root) for root in roots))
    curve_gram = root_matrix.T * gram * root_matrix
    checks = {
        "standard_lattice_determinant_32": abs(int(gram.det())) == 32,
        "standard_lattice_inertia_1_10": exact_inertia(gram) == (1, 10, 0),
        "height_zero_root_count_10": len(_ud4d4a1_height_zero_roots()) == 10,
        "level_candidate_counts": [row["dominant_candidate_count"] for row in levels]
        == [3, 12, 18, 54],
        "level_accepted_counts": [row["accepted_root_count"] for row in levels]
        == [3, 0, 0, 9],
        "simple_root_count_22": len(roots) == 22,
        "all_simple_roots_have_square_minus_two": all(
            _pair(gram, root, root) == -2 for root in roots
        ),
        "all_simple_roots_pairwise_acute": all(
            _pair(gram, left, right) >= 0
            for offset, left in enumerate(roots)
            for right in roots[offset + 1 :]
        ),
        "simple_roots_span_full_lattice_rank": root_matrix.rank() == 11,
        "wall_gram_rank_and_inertia": (
            curve_gram.rank() == 11 and exact_inertia(curve_gram) == (1, 10, 11)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"U+D4^2+A1 norm-two replay failed: {checks}")
    return _matrix_payload(curve_gram), {
        "method": (
            "exact norm-two Vinberg height enumeration on the standard "
            "U+D4(-1)^2+A1(-1) lattice"
        ),
        "dominance_chamber": (
            "x1>=x2>=x3>=|x4| independently in each D4 factor, z>=0 in A1"
        ),
        "root_equation": (
            "d^2+2|v|^2+2|w|^2+4z^2=m^2+4, with d congruent to m mod 2"
        ),
        "root_scope": "all primitive norm-minus-two roots; no other norms belong to W_2",
        "levels": levels,
        "simple_roots": [list(root) for root in roots],
        "simple_roots_sha256": _canonical_sha256([list(root) for root in roots]),
        "completion_rule": (
            "the level-8 accepted prefix is complete once its independently "
            "audited chamber has finite volume (Vinberg stopping theorem)"
        ),
        "published_chamber_identification": (
            "Roulleau's full nef cone is another fundamental W_2 chamber in "
            "the same Picard lattice, hence Weyl-isometric to this chamber"
        ),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def _source_transcription_diagnostic() -> dict[str, Any]:
    """Record why the rank-11/22 matrix is not imported from the TeX lane."""

    basis = -2 * sp.eye(11)
    for left, right, product in (
        (8, 9, 1),
        (5, 6, 1),
        (5, 4, 1),
        (5, 7, 1),
        (9, 10, 2),
        (10, 3, 2),
        (2, 3, 2),
        (3, 4, 2),
        (4, 11, 2),
        (10, 11, 2),
    ):
        basis[left - 1, right - 1] = basis[right - 1, left - 1] = product
    printed = {
        "A13": [-2, -3, -2, -2, -2, -1, -1, 2, 4, 3, 1],
        "A14": [0, 0, 0, -1, -2, -1, -1, 0, 1, 1, 0],
        "A15": [0, 1, 1, 2, 2, 1, 1, -2, -3, -2, 0],
    }
    squares = {
        label: int((sp.Matrix(vector).T * basis * sp.Matrix(vector))[0])
        for label, vector in printed.items()
    }
    return {
        "scope": (
            "diagnostic of the literal long-segment-as-chain TikZ reading; "
            "it is not used as mathematical input"
        ),
        "tex_lines": [8206, 8318],
        "expected_lattice_absolute_determinant": 32,
        "literal_basis_graph_absolute_determinant": abs(int(basis.det())),
        "printed_dependent_curve_squares": squares,
        "obstruction": (
            "this literal transcription gives determinant 128 and printed "
            "A13 square -14, so it cannot certify the claimed (-2)-curve stock"
        ),
        "resolution": (
            "bypass the ambiguous/inconsistent coordinate lane with the "
            "independent standard-lattice norm-two Vinberg reconstruction"
        ),
        "obstruction_is_terminal": False,
    }


def _load_curve_gram(
    spec: CaseSpec,
    archives: dict[str, dict[str, Any]],
    reconstructed: list[list[int]],
) -> tuple[list[list[int]], dict[str, Any]]:
    if spec.source_path is None:
        return reconstructed, {
            "kind": "independent_standard_lattice_norm_two_vinberg_reconstruction",
            "source_inventory": str(INVENTORY.relative_to(ROOT)),
        }
    archive = archives[spec.source_path]
    if archive.get("all_checks_pass") is not True:
        raise ArithmeticError(f"source archive is not green: {spec.source_path}")
    assert spec.gram_path is not None
    gram = _dig(archive, spec.gram_path)
    return gram, {
        "kind": "published_complete_curve_gram_from_frozen_source_archive",
        "archive": spec.source_path,
        "json_path": "/" + "/".join(str(value) for value in spec.gram_path),
    }


def _gram_checks(spec: CaseSpec, curve_gram: list[list[int]]) -> dict[str, Any]:
    matrix = sp.Matrix(curve_gram)
    checks = {
        "facet_count": matrix.rows == matrix.cols == spec.facets,
        "symmetric": matrix == matrix.T,
        "all_diagonal_entries_minus_two": all(
            matrix[index, index] == -2 for index in range(matrix.rows)
        ),
        "all_distinct_curve_intersections_nonnegative": all(
            matrix[left, right] >= 0
            for left in range(matrix.rows)
            for right in range(left)
        ),
        "rank_equals_picard_rank": matrix.rank() == spec.rho,
        "exact_lorentz_inertia": exact_inertia(matrix)
        == (1, spec.rho - 1, spec.facets - spec.rho),
    }
    adjacency = {
        index: {
            other
            for other in range(matrix.rows)
            if other != index and matrix[index, other] != 0
        }
        for index in range(matrix.rows)
    }
    unseen = set(range(matrix.rows))
    stack = [0]
    unseen.remove(0)
    while stack:
        vertex = stack.pop()
        new = adjacency[vertex] & unseen
        unseen.difference_update(new)
        stack.extend(new)
    checks["indecomposable_wall_diagram"] = not unseen
    if not all(checks.values()):
        raise ArithmeticError(f"curve Gram checks failed for {spec.key}: {checks}")
    return checks


def _input_stream(curve_gram: list[list[int]], dimension: int) -> str:
    rows = [f"{dimension} {len(curve_gram)}"]
    rows.extend(" ".join(str(-int(value)) for value in row) for row in curve_gram)
    return "\n".join(rows) + "\n"


def _audit_case(
    executable: Path,
    spec: CaseSpec,
    curve_gram: list[list[int]],
    provenance: dict[str, Any],
    reconstruction: dict[str, Any] | None,
) -> dict[str, Any]:
    gram_checks = _gram_checks(spec, curve_gram)
    raw_input = _input_stream(curve_gram, spec.rho - 1)
    completed = subprocess.run(
        [str(executable)],
        input=raw_input,
        text=True,
        capture_output=True,
        check=True,
    )
    record = _parse_output(completed.stdout)
    special = record["finite_special_subset_counts_by_rank"]
    endpoints = record["edge_endpoint_type_histogram"]
    edge_count = special[-2]
    finite_vertex_count = special[-1]
    checks = {
        "gram": all(gram_checks.values()),
        "dimension_and_facets": (
            record["dimension"] == spec.rho - 1
            and record["facet_count"] == spec.facets
        ),
        "all_finite_special_subsets_enumerated": (
            len(special) == spec.rho - 1
            and record["all_finite_subset_search_node_count"] == 1 + sum(special)
        ),
        "every_connected_affine_component_full_neighborhood_passes": (
            record["bad_parabolic_neighborhood_count"] == 0
        ),
        "every_elliptic_edge_has_exactly_two_endpoints": (
            record["bad_edge_count"] == 0
            and sum(endpoints.values()) == edge_count
            and sum(
                (int(kind.split(",")[0]) + int(kind.split(",")[1])) * count
                for kind, count in endpoints.items()
            )
            == 2 * edge_count
        ),
        "finite_vertex_edge_incidence": (
            2 * endpoints.get("2,0", 0) + endpoints.get("1,1", 0)
            == (spec.rho - 1) * finite_vertex_count
        ),
        "noncompact_has_ideal_vertices": record["ideal_vertex_count"] > 0,
        "reconstructed_prefix_stops_when_finite": (
            reconstruction is None
            or (
                reconstruction["all_checks_pass"]
                and record["bad_edge_count"] == 0
                and record["bad_parabolic_neighborhood_count"] == 0
            )
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"literal Vinberg audit failed for {spec.key}: {checks}")
    return {
        "key": spec.key,
        "rho": spec.rho,
        "hyperbolic_dimension": spec.rho - 1,
        "lattice": spec.lattice,
        "facet_count": spec.facets,
        "matrix_provenance": provenance,
        "curve_gram_sha256": _canonical_sha256(curve_gram),
        "positive_normal_input_sha256": hashlib.sha256(raw_input.encode()).hexdigest(),
        "literal_verifier_stdout_sha256": hashlib.sha256(
            completed.stdout.encode()
        ).hexdigest(),
        "gram_checks": gram_checks,
        **record,
        "elliptic_edge_count": edge_count,
        "finite_vertex_count": finite_vertex_count,
        "criterion": (
            "every elliptic rank-(n-1) edge has exactly two finite/ideal "
            "endpoints and every connected affine component's full zero-"
            "neighborhood is parabolic of rank n-1"
        ),
        "simple_roots_complete": True,
        "completion_evidence": (
            "exact norm-two height replay plus finite-volume Vinberg stopping"
            if reconstruction is not None
            else "published complete irreducible (-2)-curve configuration"
        ),
        "decision": "finite_volume_noncompact",
        "finite_volume": True,
        "compact": False,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def build_certificate(
    cpp_source: Path = CPP_SOURCE,
    compiler: str | None = None,
) -> dict[str, Any]:
    source_paths = sorted(
        {spec.source_path for spec in CASES if spec.source_path is not None}
    )
    archives = {
        path: json.loads((ROOT / path).read_text(encoding="utf-8"))
        for path in source_paths
    }
    inventory = json.loads(INVENTORY.read_text(encoding="utf-8"))
    if inventory.get("all_checks_pass") is not True:
        raise ArithmeticError("the rank-11--14 source inventory is not green")
    stock_keys = {stock["key"] for stock in inventory["stocks"][:14]}
    if stock_keys != {spec.key for spec in CASES}:
        raise ArithmeticError("the frozen fourteen-stock inventory changed")

    reconstructed_gram, reconstruction = _reconstruct_ud4d4a1()
    selected_compiler = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected_compiler is None:
        raise RuntimeError("a C++17 compiler is required")
    with tempfile.TemporaryDirectory(prefix="k3-rank11-14-vinberg-") as raw:
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
        records = []
        for spec in CASES:
            curve_gram, provenance = _load_curve_gram(
                spec, archives, reconstructed_gram
            )
            if len(curve_gram) != spec.facets:
                raise ArithmeticError(f"wrong wall count for {spec.key}")
            records.append(
                _audit_case(
                    executable,
                    spec,
                    curve_gram,
                    provenance,
                    reconstruction if spec.key == "R11_UD4D4A1" else None,
                )
            )

    totals = {
        "case_count": len(records),
        "facet_count": sum(record["facet_count"] for record in records),
        "elliptic_edge_count": sum(
            record["elliptic_edge_count"] for record in records
        ),
        "finite_vertex_count": sum(
            record["finite_vertex_count"] for record in records
        ),
        "ideal_vertex_count": sum(
            record["ideal_vertex_count"] for record in records
        ),
        "connected_affine_component_count": sum(
            record["connected_affine_component_count"] for record in records
        ),
        "bad_edge_count": sum(record["bad_edge_count"] for record in records),
        "bad_parabolic_neighborhood_count": sum(
            record["bad_parabolic_neighborhood_count"] for record in records
        ),
    }
    checks = {
        "exactly_fourteen_inventory_cases": len(records) == 14,
        "ranks_11_through_14": sorted(Counter(row["rho"] for row in records).items())
        == [(11, 4), (12, 4), (13, 3), (14, 3)],
        "all_literal_vinberg_audits_pass": all(
            record["all_checks_pass"] for record in records
        ),
        "all_finite_volume_noncompact": all(
            record["decision"] == "finite_volume_noncompact" for record in records
        ),
        "zero_bad_edges": totals["bad_edge_count"] == 0,
        "zero_bad_parabolic_neighborhoods": (
            totals["bad_parabolic_neighborhood_count"] == 0
        ),
        "rank11_22_independent_reconstruction_passes": reconstruction[
            "all_checks_pass"
        ],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global literal Vinberg checks failed: {checks}")
    return {
        "schema": SCHEMA,
        "artifact_id": "k3_rank11_14_literal_vinberg_volume",
        "status_date": "2026-08-08",
        "source_artifacts": {
            str(INVENTORY.relative_to(ROOT)): _sha256(INVENTORY),
            **{path: _sha256(ROOT / path) for path in source_paths},
            str(cpp_source.relative_to(ROOT)): _sha256(cpp_source),
        },
        "method": {
            "literal_gate": (
                "Vinberg/CoxIter edge-endpoint criterion with full affine "
                "zero-neighborhoods, applied to complete simple walls"
            ),
            "finite_enumeration": (
                "generic exact crystallographic Coxeter-subdiagram DFS"
            ),
            "affine_enumeration": (
                "all irreducible affine components reached as first nonfinite "
                "extensions; parallel affine A1 pairs scanned separately"
            ),
            "arithmetic": (
                "integer Coxeter classification and overflow-guarded signed "
                "128-bit Bareiss determinants"
            ),
            "compiler": Path(selected_compiler).name,
        },
        "claim_scope": {
            "literal_polytope_claim": (
                "all fourteen complete full-root chambers are finite-volume "
                "noncompact hyperbolic Coxeter polytopes"
            ),
            "structural_claim_kept_separate": (
                "finite automorphism/Cox-ring arithmetic is not substituted "
                "for the endpoint and parabolic-neighborhood audit"
            ),
            "surface_automorphisms": (
                "no new surface-automorphism realization is inferred here"
            ),
        },
        "rank11_22_reconstruction": reconstruction,
        "rank11_22_source_transcription_diagnostic": (
            _source_transcription_diagnostic()
        ),
        "cases": records,
        "totals": totals,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpp-source", type=Path, default=CPP_SOURCE)
    parser.add_argument("--compiler")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.cpp_source, args.compiler)
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if payload != archived:
            raise ArithmeticError("rebuilt rank-11--14 literal Vinberg archive changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "all_checks_pass": payload["all_checks_pass"],
                "totals": payload["totals"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
