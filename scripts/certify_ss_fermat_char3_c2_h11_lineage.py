#!/usr/bin/env python3
"""Certify the exact lineage of the 206-wall char-3 C2 fixed section.

The saturated fixed lattice is identified by an explicit determinant-one
change of basis with ``U + E8(-2) + A1(-2)^2``.  A second explicit change of
basis embeds Alexeev's published 14-root Coxeter Gram for the reflective
2-elementary lattice ``(r,a,delta)=(12,10,1)`` in the same fixed lattice.
Finally all 206 facet mirrors are reconstructed from the complete published
5,944-wall stock and checked to be primitive integral reflections.

No claim that the 206-wall subgroup chamber was previously unknown is made:
the bounded corpus checks only exclude a literal occurrence in the pinned
CoxIter snapshots.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import tarfile
from typing import Any, Sequence

import numpy as np
import sympy as sp

from scripts.certify_higher_dimensional_rational_lineage import _compare
from scripts.certify_jacobian_kummer_h10_reflection_index import _relation_matrix
from scripts.certify_ss_fermat_char3_survivor_fixed_lattices import (
    _primitive_integer_column,
)
from scripts.screen_ss_fermat_quartic_char3_c2_survivors import wall_cast


ROOT = Path(__file__).resolve().parents[1]
FIXED = (
    ROOT
    / "data/enumerations/ss_fermat_quartic_char3_survivor_fixed_lattices.json"
)
CYCLIC = (
    ROOT
    / "data/enumerations/ss_fermat_quartic_char3_c2_element_class_screen.json"
)
REDUCTION = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_survivor_screen.json"
)
COXETER_SCREEN = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_coxeter_screen.json"
)
COMPONENT_AUDIT = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_component_audit.json"
)
VINBERG = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_vinberg_certificate.json"
)
OUTPUT = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_lineage.json"
)

DEFAULT_ALEXEEV_SOURCE = Path("/private/tmp/alexeev-2209.09110v4.tar")
DEFAULT_COXITER_CORPUS = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")

ALEXEEV_PAPER_URL = "https://arxiv.org/abs/2209.09110"
ALEXEEV_SOURCE_SHA256 = (
    "10e4002f6518f8bf6d669946775bbe42731f69277a52847c5602fa60ea33c2a6"
)
ALEXEEV_MATRIX_FILE_SHA256 = (
    "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
)
ALEXEEV_TEX_SHA256 = (
    "006e2a1eac2aa3c4916ef66aceb27e8e80a35145162e6785a4874843f249f1df"
)
ALEXEEV_MATRIX_SHA256 = (
    "a12618b59ad911429dc4dfe1a18a1d7b6b64a7f5dae53fcf0c677f28c43eb9e5"
)
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
DOUBLING_ARCHIVE_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)


# Columns are the standard basis of U + E8(-2) + A1(-2)^2 in the archived
# short fixed basis.  This certifies P^T G_fixed P = G_standard.
C2_H11_TO_STANDARD = sp.Matrix(
    [
        [0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        [0, -1, -1, 0, 0, -1, 1, -1, 1, -1, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0],
        [0, 0, -1, 1, 0, 0, 0, 0, -1, 1, 0, 0],
        [0, 0, 2, -1, 0, 0, 0, 1, -1, 1, 0, 0],
        [0, 0, 0, -1, 0, 1, 0, 0, -1, 1, 1, 0],
        [-1, -1, 1, 1, 1, -3, 2, -2, 0, 0, 0, -1],
        [0, 0, -1, 1, 0, -1, 1, -1, 1, -2, -1, 0],
        [0, 0, 1, -1, 1, 0, 0, 0, -1, 1, 0, 0],
        [0, 0, 0, -1, 0, 1, -1, 1, 0, 1, 1, 0],
        [0, 0, 0, -1, 0, 0, 0, 0, 1, -1, -1, 0],
        [0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 0],
    ]
)


# Alexeev, ``Reflexive hyperbolic 2-elementary lattices, K3 surfaces and
# hyperkahler manifolds'', distributed Sage data, coxeter_matrix[(12,10,1)].
# The sign convention is signature (1,11), with negative root squares.
ALEXEEV_COXETER_GRAM = sp.Matrix(
    [
        [-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0],
        [0, -4, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, -4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2],
        [0, 0, 2, -4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 2, 0, 2, -4, 2, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 2, -4, 2, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 2, -4, 2, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 2, -4, 2, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 2, -4, 0, 0, 0, 2, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 2, 0, 0, 2],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -4, 2, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0],
        [2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -4, 0],
        [0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -4],
    ]
)


# Columns give the first twelve Alexeev roots (which form an integral basis)
# in the archived C2-H11 fixed basis.  The final two roots are recovered from
# the full Gram rather than trusted as additional constants.
ALEXEEV_BASIS_ROOTS_IN_FIXED_BASIS = sp.Matrix(
    [
        [1, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, 1],
        [1, 1, -1, -1, 1, -1, 1, -1, 0, 0, 0, 0],
        [0, -1, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0],
        [0, -2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, -1, 0, 0, 0, 1, 0, 1, -1, 0],
        [0, -1, 0, 1, 1, -2, 2, -3, 1, 0, -1, 2],
        [0, 1, -2, 1, 0, -1, 1, -1, 0, -1, 1, 0],
        [0, -1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 1, -1, 0, 1, -1, 1, 0, 1, -1, 0],
        [0, 0, -1, -1, 1, 0, 0, 0, 0, -1, 1, 0],
        [0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0],
    ]
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.is_Integer is not True for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return [[int(value) for value in row] for row in matrix.tolist()]


def _matrix_sha256(matrix: sp.Matrix) -> str:
    payload = ";".join(
        ",".join(str(int(matrix[row, column])) for column in range(matrix.cols))
        for row in range(matrix.rows)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _e8_cartan() -> sp.Matrix:
    matrix = 2 * sp.eye(8)
    for left, right in (
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (4, 6),
        (6, 7),
    ):
        matrix[left, right] = matrix[right, left] = -1
    return matrix


def _standard_gram() -> sp.Matrix:
    return sp.diag(
        sp.Matrix([[0, 1], [1, 0]]),
        -2 * _e8_cartan(),
        sp.diag(-2, -2),
    )


def _root_type(gram: sp.Matrix, root: sp.Matrix) -> tuple[int, int]:
    square = int((root.T * gram * root)[0])
    divisibility = math.gcd(*(abs(int(value)) for value in gram * root))
    return square, divisibility


def _reflection_is_integral(gram: sp.Matrix, root: sp.Matrix) -> bool:
    square = int((root.T * gram * root)[0])
    reflection = sp.eye(gram.rows) - sp.Rational(2, square) * root * (
        root.T * gram
    )
    return (
        all(value.is_Integer is True for value in reflection)
        and reflection.T * gram * reflection == gram
        and reflection * reflection == sp.eye(gram.rows)
        and reflection * root == -root
    )


def _alexeev_roots(gram: sp.Matrix) -> tuple[sp.Matrix, list[list[int]]]:
    basis_gram = ALEXEEV_COXETER_GRAM[:12, :12]
    basis_roots = ALEXEEV_BASIS_ROOTS_IN_FIXED_BASIS
    if not (
        abs(int(basis_roots.det())) == 1
        and basis_roots.T * gram * basis_roots == basis_gram
    ):
        raise ArithmeticError("Alexeev basis-root isometry failed")
    coordinates: list[sp.Matrix] = []
    for root_index in range(14):
        coordinate = basis_gram.inv() * ALEXEEV_COXETER_GRAM[:12, root_index]
        if any(value.is_Integer is not True for value in coordinate):
            raise ArithmeticError("an Alexeev root left the basis-root lattice")
        if (coordinate.T * basis_gram * coordinate)[0] != ALEXEEV_COXETER_GRAM[
            root_index, root_index
        ]:
            raise ArithmeticError("an Alexeev root square failed")
        coordinates.append(coordinate)
    roots = basis_roots * sp.Matrix.hstack(*coordinates)
    if roots.T * gram * roots != ALEXEEV_COXETER_GRAM:
        raise ArithmeticError("full Alexeev Coxeter Gram failed")
    return roots, [[int(value) for value in coordinate] for coordinate in coordinates]


def _reconstruct_206_mirrors(
    gram: sp.Matrix, fixed_basis: sp.Matrix
) -> tuple[list[dict[str, Any]], Counter[tuple[int, int]]]:
    cyclic = json.loads(CYCLIC.read_text(encoding="utf-8"))
    survivor = next(
        record
        for record in cyclic["high_fixed_rank_casts"]
        if record["bad_pair_count"] == 0
        and record["hyperbolic_dimension"] == 11
    )
    cast = wall_cast(np.asarray(survivor["representative_matrix"], dtype=np.int64))
    reduction = json.loads(REDUCTION.read_text(encoding="utf-8"))
    reduced = next(
        row for row in reduction["casts"] if row["hyperbolic_dimension"] == 11
    )
    records: list[dict[str, Any]] = []
    histogram: Counter[tuple[int, int]] = Counter()
    standard_inverse = C2_H11_TO_STANDARD.inv()
    for restricted_ray_index in reduced["facet_restricted_ray_indices"]:
        ambient_normal = sp.Matrix(cast["rays"][int(restricted_ray_index) - 1])
        fixed_coordinates = _primitive_integer_column(
            gram.inv() * fixed_basis * ambient_normal
        )
        standard_coordinates = standard_inverse * fixed_coordinates
        if any(value.is_Integer is not True for value in standard_coordinates):
            raise ArithmeticError("a C2-H11 mirror left the standard lattice")
        square, divisibility = _root_type(gram, fixed_coordinates)
        reflection_valid = _reflection_is_integral(gram, fixed_coordinates)
        histogram[(square, divisibility)] += 1
        records.append(
            {
                "restricted_ray_index": int(restricted_ray_index),
                "fixed_basis_coordinates": [
                    int(value) for value in fixed_coordinates
                ],
                "standard_lattice_coordinates": [
                    int(value) for value in standard_coordinates
                ],
                "square": square,
                "divisibility": divisibility,
                "reflection_integral_involutive_isometry": reflection_valid,
            }
        )
    return records, histogram


def _first_header(lines: Sequence[str]) -> tuple[int, int]:
    content = next(
        stripped
        for line in lines
        if (stripped := line.split("#", 1)[0].strip())
    )
    facets, dimension = map(int, content.split()[:2])
    return facets, dimension


def _verify_external_sources(
    alexeev_source: Path, coxiter_corpus: Path
) -> None:
    if _sha256(alexeev_source) != ALEXEEV_SOURCE_SHA256:
        raise ArithmeticError("Alexeev source archive hash mismatch")
    with tarfile.open(alexeev_source) as archive:
        matrix_raw = archive.extractfile("all_cases_in_2el_paper.txt")
        tex_raw = archive.extractfile("g0.tex")
        if matrix_raw is None or tex_raw is None:
            raise ArithmeticError("Alexeev source member missing")
        matrix_payload = matrix_raw.read()
        tex_payload = tex_raw.read()
    if hashlib.sha256(matrix_payload).hexdigest() != ALEXEEV_MATRIX_FILE_SHA256:
        raise ArithmeticError("Alexeev matrix-file hash mismatch")
    if hashlib.sha256(tex_payload).hexdigest() != ALEXEEV_TEX_SHA256:
        raise ArithmeticError("Alexeev TeX hash mismatch")
    if _sha256(COXITER_ARCHIVE) != COXITER_ARCHIVE_SHA256:
        raise ArithmeticError("CoxIter archive hash mismatch")

    ordinary: Counter[tuple[int, int]] = Counter()
    ordinary_total = 0
    for path in coxiter_corpus.rglob("*.coxiter"):
        facets, dimension = _first_header(
            path.read_text(encoding="utf-8").splitlines()
        )
        ordinary[(dimension, facets)] += 1
        ordinary_total += 1
    expected_ordinary_h11 = Counter({(11, 12): 14, (11, 13): 6, (11, 15): 2})
    if ordinary_total != 1_214 or Counter(
        {key: value for key, value in ordinary.items() if key[0] == 11}
    ) != expected_ordinary_h11:
        raise ArithmeticError("ordinary CoxIter H11 census changed")

    doubling_path = coxiter_corpus / "graphs/testing_index2.tar.gz"
    if _sha256(doubling_path) != DOUBLING_ARCHIVE_SHA256:
        raise ArithmeticError("CoxIter doubling-archive hash mismatch")
    doubling: Counter[tuple[int, int]] = Counter()
    doubling_total = 0
    with tarfile.open(doubling_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError("failed to read a doubling record")
            lines = [raw.decode("utf-8") for raw in handle]
            facets, dimension = _first_header(lines)
            doubling[(dimension, facets)] += 1
            doubling_total += 1
    expected_doubling_h11 = Counter(
        {
            (11, 15): 2,
            (11, 16): 2,
            (11, 18): 2,
            (11, 22): 2,
            (11, 30): 2,
            (11, 46): 2,
        }
    )
    if doubling_total != 2_959 or Counter(
        {key: value for key, value in doubling.items() if key[0] == 11}
    ) != expected_doubling_h11:
        raise ArithmeticError("doubling CoxIter H11 census changed")


def build_certificate() -> dict[str, Any]:
    fixed = json.loads(FIXED.read_text(encoding="utf-8"))
    case = next(row for row in fixed["cases"] if row["id"] == "char3-C2-H11")
    gram = sp.Matrix(case["gram"])
    fixed_basis = sp.Matrix(case["basis_rows_in_published_NS_basis"])
    standard = _standard_gram()
    standard_isometry = C2_H11_TO_STANDARD
    standard_isometry_valid = (
        abs(int(standard_isometry.det())) == 1
        and standard_isometry.T * gram * standard_isometry == standard
    )

    alexeev_roots, alexeev_root_basis_coordinates = _alexeev_roots(gram)
    alexeev_types = Counter(
        _root_type(gram, alexeev_roots[:, index]) for index in range(14)
    )
    alexeev_reflections_integral = all(
        _reflection_is_integral(gram, alexeev_roots[:, index])
        for index in range(14)
    )
    alexeev_relations = _relation_matrix(
        [
            [str(int(-ALEXEEV_COXETER_GRAM[row, column])) for column in range(14)]
            for row in range(14)
        ]
    )

    mirrors, mirror_types = _reconstruct_206_mirrors(gram, fixed_basis)
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    components = json.loads(COMPONENT_AUDIT.read_text(encoding="utf-8"))
    coxeter_screen = json.loads(COXETER_SCREEN.read_text(encoding="utf-8"))
    reduction_payload = json.loads(REDUCTION.read_text(encoding="utf-8"))
    h11_reduction = next(
        row
        for row in reduction_payload["casts"]
        if row["hyperbolic_dimension"] == 11
    )
    rational_lineage = _compare(
        gram, sp.diag(1, *([-1] * 11)), Fraction(1)
    )

    checks = {
        "upstream_fixed_lattice_is_saturated_even_lorentzian": (
            case["smith_kernel"]["saturated_smith_kernel"]
            and case["signature"] == [1, 11]
            and case["even"]
            and case["determinant_absolute"] == 1_024
            and case["smith_diagonal_absolute"] == [1, 1, *([2] * 10)]
        ),
        "explicit_standard_lattice_isometry": standard_isometry_valid,
        "alexeev_matrix_literal_hash": (
            _matrix_sha256(ALEXEEV_COXETER_GRAM) == ALEXEEV_MATRIX_SHA256
        ),
        "alexeev_first_twelve_roots_are_an_integral_basis": (
            abs(int(ALEXEEV_COXETER_GRAM[:12, :12].det())) == 1_024
            and abs(int(ALEXEEV_BASIS_ROOTS_IN_FIXED_BASIS.det())) == 1
        ),
        "full_alexeev_14_root_gram_matches_entrywise": (
            alexeev_roots.T * gram * alexeev_roots == ALEXEEV_COXETER_GRAM
        ),
        "all_alexeev_roots_are_integral_reflecting_roots": (
            alexeev_reflections_integral
            and alexeev_types == Counter({(-4, 2): 11, (-2, 1): 2, (-2, 2): 1})
        ),
        "exactly_206_distinct_primitive_facet_mirrors": (
            len(mirrors) == 206
            and len(
                {tuple(record["fixed_basis_coordinates"]) for record in mirrors}
            )
            == 206
        ),
        "all_206_facet_reflections_are_integral": (
            all(
                record["reflection_integral_involutive_isometry"]
                for record in mirrors
            )
            and mirror_types
            == Counter({(-2, 2): 164, (-4, 2): 30, (-2, 1): 12})
        ),
        "complete_206_wall_coxeter_gate": (
            h11_reduction["exact_complete"]
            and h11_reduction["facet_count"] == 206
            and h11_reduction["unresolved_count"] == 0
            and coxeter_screen["facet_count"] == 206
            and coxeter_screen["status"]["all_facet_pairs_are_coxeter"]
        ),
        "literal_vinberg_finite_volume_noncompact": (
            vinberg["finite_volume"]
            and not vinberg["compact"]
            and vinberg["bad_edge_count"] == 0
            and vinberg["unrecognized_endpoint_edge_count"] == 0
            and components["invalid_full_zero_neighborhood_count"] == 0
        ),
        "rational_space_is_odd_unimodular_H11_space": rational_lineage[
            "rationally_similar"
        ],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"C2-H11 lineage checks failed: {checks}")

    return {
        "schema": "ss_fermat_quartic_char3_c2_h11_lineage_v1",
        "candidate": {
            "id": "F3-C2-H11-206",
            "dimension": 11,
            "facet_count": 206,
            "finite_vertex_count": int(vinberg["finite_vertex_count"]),
            "ideal_vertex_count": int(vinberg["ideal_vertex_count"]),
            "volume_decision": vinberg["vinberg_decision"],
        },
        "fixed_lattice_identification": {
            "target_lattice": "U + E8(-2) + A1(-2)^2",
            "nikulin_invariants": {"r": 12, "a": 10, "delta": 1},
            "target_gram": _matrix_payload(standard),
            "identity": "P^T * G_fixed * P = G_target",
            "explicit_isometry_P": _matrix_payload(standard_isometry),
            "explicit_isometry_determinant": int(standard_isometry.det()),
            "explicit_isometry_sha256": _matrix_sha256(standard_isometry),
            "rational_comparison_with_I_1_11": rational_lineage,
        },
        "known_maximal_reflection_chamber": {
            "source": {
                "author": "Valery Alexeev",
                "paper": "Reflexive hyperbolic 2-elementary lattices, K3 surfaces and hyperkahler manifolds",
                "url": ALEXEEV_PAPER_URL,
                "source_archive_sha256": ALEXEEV_SOURCE_SHA256,
                "distributed_matrix_file_sha256": ALEXEEV_MATRIX_FILE_SHA256,
                "paper_tex_sha256": ALEXEEV_TEX_SHA256,
                "source_key": "coxeter_matrix[(12,10,1)]",
                "paper_statements": {
                    "reflective": True,
                    "facet_root_count": 14,
                    "maximal_parabolic_subdiagram_count": 5,
                    "diagram_automorphism_group": "trivial",
                },
            },
            "coxeter_gram": _matrix_payload(ALEXEEV_COXETER_GRAM),
            "coxeter_gram_sha256": _matrix_sha256(ALEXEEV_COXETER_GRAM),
            "normalized_relation_matrix": alexeev_relations,
            "first_twelve_root_isometry_to_fixed_lattice": _matrix_payload(
                ALEXEEV_BASIS_ROOTS_IN_FIXED_BASIS
            ),
            "all_root_coordinates_in_first_twelve_root_basis": (
                alexeev_root_basis_coordinates
            ),
            "all_root_coordinates_in_fixed_basis": _matrix_payload(
                alexeev_roots.T
            ),
            "primitive_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(alexeev_types.items())
            },
            "all_reflections_integral_involutive_isometries": (
                alexeev_reflections_integral
            ),
            "maximality_reason": (
                "Alexeev's reflective classification identifies this 14-root "
                "matrix as the full Gamma_r chamber of the unique even "
                "2-elementary lattice (12,10,1); the explicit isometry above "
                "places that chamber in the candidate fixed lattice."
            ),
        },
        "induced_206_wall_reflection_subgroup": {
            "facet_count": len(mirrors),
            "primitive_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(mirror_types.items())
            },
            "all_reflections_integral_involutive_isometries": True,
            "mirrors": mirrors,
            "subgroup_conclusion": (
                "W_206 is a literal reflection subgroup of the known "
                "14-wall full reflection group W(12,10,1)."
            ),
            "finite_index_reason": (
                "Both the 206-wall chamber and Alexeev's full reflection "
                "chamber have finite hyperbolic covolume."
            ),
            "numeric_index_status": (
                "not computed: in H11 both rational Euler characteristics "
                "vanish, so a Solomon-Euler ratio is 0/0; use explicit "
                "cosets/tiling or an odd-dimensional hyperbolic-volume "
                "comparison"
            ),
        },
        "bounded_weighted_diagram_audit": {
            "coxiter_snapshot": {
                "repository": "https://github.com/rgugliel/CoxIter",
                "commit": COXITER_COMMIT,
                "archive_sha256": COXITER_ARCHIVE_SHA256,
                "ordinary_file_count": 1_214,
                "ordinary_H11_facet_histogram": {"12": 14, "13": 6, "15": 2},
                "ordinary_global_maximum_facets": 50,
                "doubling_archive_sha256": DOUBLING_ARCHIVE_SHA256,
                "doubling_file_count": 2_959,
                "doubling_H11_facet_histogram": {
                    "15": 2,
                    "16": 2,
                    "18": 2,
                    "22": 2,
                    "30": 2,
                    "46": 2,
                },
                "doubling_global_maximum_facets": 264,
            },
            "literal_206_facet_record_present": False,
            "scope_warning": (
                "Facet-count exclusion proves absence only from these two "
                "pinned corpora; it is not a worldwide chamber-priority proof."
            ),
        },
        "lineage_decision": {
            "new_lattice": False,
            "new_maximal_reflection_group": False,
            "new_arithmetic_commensurability_class": False,
            "positive_algebraic_construction": True,
            "subgroup_chamber_priority": "open after bounded corpus exclusion",
        },
        "upstream_artifacts": {
            str(path.relative_to(ROOT)): _sha256(path)
            for path in (
                FIXED,
                CYCLIC,
                REDUCTION,
                COXETER_SCREEN,
                COMPONENT_AUDIT,
                VINBERG,
            )
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--verify-external", action="store_true")
    parser.add_argument("--alexeev-source", type=Path, default=DEFAULT_ALEXEEV_SOURCE)
    parser.add_argument("--coxiter-corpus", type=Path, default=DEFAULT_COXITER_CORPUS)
    arguments = parser.parse_args()
    if arguments.verify_external:
        _verify_external_sources(arguments.alexeev_source, arguments.coxiter_corpus)
    result = build_certificate()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if (
            not arguments.output.is_file()
            or arguments.output.read_text(encoding="utf-8") != rendered
        ):
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(arguments.output)


if __name__ == "__main__":
    main()
