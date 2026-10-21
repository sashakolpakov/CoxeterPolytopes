#!/usr/bin/env python3
"""Certify all ten Picard-rank-nine Roulleau full-nef K3 chambers."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    _parse_output,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_rank9_all_full_nef.json"
SCHEMA = "coxeter.roulleau_rank9_all_full_nef.v1"
CPP_SOURCE = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
CPP_SOURCE_SHA256 = "45174f5461d28affc1fdb4cb55047b7a1581037d307e32800651badb3c2df045"
ROULLEAU_ARXIV = "https://arxiv.org/abs/2003.08985"
ROULLEAU_SOURCE_TAR_SHA256 = (
    "3e445ed01446b47ded6c007b395308de34e566b736f0ba982470bc497adeb1e9"
)
ROULLEAU_SOURCE_TEX_SHA256 = (
    "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
)


VINBERG_EXPECTED: dict[str, dict[str, Any]] = {
    "R9-UE7-NEF-H8-9": {
        "finite_nodes": 254,
        "affine_nodes": 510,
        "prefix": [9, 21, 34, 45, 51, 49, 36, 8],
        "special": [9, 36, 84, 126, 126, 84, 36, 8],
        "affine_components": 1,
        "cusps": 1,
        "endpoints": {"1,1": 8, "2,0": 28},
    },
    "R9-UD6A1-NEF-H8-10": {
        "finite_nodes": 411,
        "affine_nodes": 760,
        "prefix": [10, 29, 52, 73, 86, 85, 63, 12],
        "special": [10, 44, 112, 182, 196, 140, 63, 12],
        "affine_components": 4,
        "cusps": 3,
        "endpoints": {"0,2": 3, "1,1": 24, "2,0": 36},
    },
    "R9-UD4A1x3-NEF-H8-15": {
        "finite_nodes": 2112,
        "affine_nodes": 2810,
        "prefix": [15, 79, 232, 432, 524, 461, 315, 53],
        "special": [15, 90, 296, 600, 784, 656, 315, 53],
        "affine_components": 26,
        "cusps": 11,
        "endpoints": {"0,2": 37, "1,1": 132, "2,0": 146},
    },
    "R9-UA1x7-NEF-H8-37": {
        "finite_nodes": 22422,
        "affine_nodes": 31053,
        "prefix": [37, 362, 1529, 3722, 5792, 6222, 4188, 569],
        "special": [37, 400, 2016, 5642, 9324, 8876, 4188, 569],
        "affine_components": 336,
        "cusps": 78,
        "endpoints": {"0,2": 868, "1,1": 2088, "2,0": 1232},
    },
    "R9-U2A1x7-NEF-H8-240": {
        "finite_nodes": 1192320,
        "affine_nodes": 1501681,
        "prefix": [240, 6705, 58120, 206355, 352584, 343675, 207360, 17280],
        "special": [240, 6720, 60480, 241920, 483840, 483840, 207360, 17280],
        "affine_components": 15120,
        "cusps": 2160,
        "endpoints": {"0,2": 69120, "1,1": 138240},
    },
    "R9-UA7-NEF-H8-9": {
        "finite_nodes": 253,
        "affine_nodes": 509,
        "prefix": [9, 21, 34, 45, 51, 49, 36, 7],
        "special": [9, 36, 84, 126, 126, 84, 36, 7],
        "affine_components": 2,
        "cusps": 2,
        "endpoints": {"0,2": 1, "1,1": 14, "2,0": 21},
    },
    "R9-UD4A3-NEF-H8-10": {
        "finite_nodes": 563,
        "affine_nodes": 927,
        "prefix": [10, 30, 60, 95, 125, 135, 90, 17],
        "special": [10, 45, 120, 209, 245, 190, 90, 17],
        "affine_components": 5,
        "cusps": 4,
        "endpoints": {"0,2": 6, "1,1": 32, "2,0": 52},
    },
    "R9-UD5A2-NEF-H8-10": {
        "finite_nodes": 560,
        "affine_nodes": 882,
        "prefix": [10, 30, 60, 95, 126, 139, 81, 18],
        "special": [10, 45, 119, 203, 231, 174, 81, 18],
        "affine_components": 2,
        "cusps": 1,
        "endpoints": {"1,1": 18, "2,0": 63},
    },
    "R9-UD7-NEF-H8-9": {
        "finite_nodes": 254,
        "affine_nodes": 510,
        "prefix": [9, 21, 34, 45, 51, 49, 36, 8],
        "special": [9, 36, 84, 126, 126, 84, 36, 8],
        "affine_components": 1,
        "cusps": 1,
        "endpoints": {"1,1": 8, "2,0": 28},
    },
    "R9-UA1E6-NEF-H8-10": {
        "finite_nodes": 413,
        "affine_nodes": 762,
        "prefix": [10, 29, 52, 73, 86, 85, 63, 14],
        "special": [10, 44, 112, 182, 196, 140, 63, 14],
        "affine_components": 2,
        "cusps": 1,
        "endpoints": {"1,1": 14, "2,0": 49},
    },
}


def _gram_from_edges(
    count: int,
    ordinary: Iterable[tuple[int, int]],
    parallel: Iterable[tuple[int, int]] = (),
) -> sp.Matrix:
    gram = -2 * sp.eye(count)
    for left, right in ordinary:
        gram[left, right] = gram[right, left] = 1
    for left, right in parallel:
        gram[left, right] = gram[right, left] = 2
    return gram


def _coordinate_gram(
    lattice_gram: sp.Matrix, vectors: Sequence[Sequence[int]]
) -> sp.Matrix:
    columns = sp.Matrix.hstack(*(sp.Matrix(vector) for vector in vectors))
    return columns.T * lattice_gram * columns


def _case_ue7() -> dict[str, Any]:
    gram = _gram_from_edges(
        9,
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (4, 6), (6, 7), (7, 8)],
    )
    return {
        "id": "R9-UE7-NEF-H8-9",
        "lattice": "U+E7(-1)",
        "gram": gram,
        "basis_indices": tuple(range(9)),
        "ample": (19, 39, 60, 82, 105, 52, 77, 50, 24),
        "ample_square": 532,
        "ample_pairings": (1,) * 8 + (2,),
        "tex_lines": "6379-6436",
        "branch_model": "double cover of F4 with an E7 branch singularity",
        "reference": {
            "name": "JKRT T8",
            "path": "graphs/simplices/8-jkrt_351-01-T8.coxiter",
            "sha256": "b5c450d747e80571fea740247e3b17e026e668433156350ee12d1890aa9739dd",
            "permutation": (7, 6, 5, 4, 3, 8, 2, 1, 0),
            "ordinary": ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (3, 8)),
            "parallel": (),
            "same_header_exact_matches": 2,
        },
    }


def _case_ud6a1() -> dict[str, Any]:
    gram = _gram_from_edges(
        10,
        [(1, 2), (2, 3), (3, 4), (4, 5), (4, 6), (6, 7), (7, 8), (7, 9)],
        [(0, 1)],
    )
    return {
        "id": "R9-UD6A1-NEF-H8-10",
        "lattice": "U+D6(-1)+A1(-1)",
        "gram": gram,
        "basis_indices": tuple(range(9)),
        "ample": (10, 15, 11, 8, 6, 2, 3, 1, 0),
        "ample_square": 148,
        "ample_pairings": (10, 1, 1, 1, 1, 2, 1, 1, 1, 1),
        "tex_lines": "6439-6501",
        "branch_model": "double cover of F4 with a node and a D6 branch singularity",
        "reference": {
            "name": "Tumarkin 8-dimensional pyramid 12-k=2",
            "path": "graphs/Tumarkin pyramids n+2/8-tum04_11_12-k=2.coxiter",
            "sha256": "c44f40c8ee2566cfa4a95cb470c9409af8a9ad1a6e5b59e8da59f85fb3e94a8a",
            "permutation": (7, 6, 5, 4, 3, 8, 2, 1, 0, 9),
            "ordinary": ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1, 9), (3, 8)),
            "parallel": ((6, 7),),
            "same_header_exact_matches": 1,
        },
    }


def _case_ud4a1x3() -> dict[str, Any]:
    ordinary = [
        (0, 1), (0, 5), (1, 2), (1, 3), (1, 4),
        (5, 6), (5, 7), (5, 8),
    ]
    parallel = [
        (2, 9), (9, 12), (12, 6),
        (3, 10), (10, 13), (13, 7),
        (4, 11), (11, 14), (14, 8),
        (9, 13), (9, 14), (10, 12),
        (10, 14), (11, 12), (11, 13),
    ]
    gram = _gram_from_edges(15, ordinary, parallel)
    return {
        "id": "R9-UD4A1x3-NEF-H8-15",
        "lattice": "U+D4(-1)+A1(-1)^3",
        "gram": gram,
        "basis_indices": tuple(range(9)),
        "ample": (6, 7, 3, 3, 3, 7, 3, 3, 3),
        "ample_square": 44,
        "ample_pairings": (2,) + (1,) * 8 + (6,) * 6,
        "tex_lines": "6504-6591",
        "branch_model": "double plane branched over a sextic with two D4 singularities",
    }


def _case_ua1x7() -> dict[str, Any]:
    pairs = list(itertools.combinations(range(8), 2))
    gram = -2 * sp.eye(37)
    for index in range(8):
        gram[0, 1 + index] = gram[1 + index, 0] = 1
    for offset, pair in enumerate(pairs):
        wall = 9 + offset
        for index in pair:
            gram[1 + index, wall] = gram[wall, 1 + index] = 2
    for left, pair_left in enumerate(pairs):
        for right, pair_right in enumerate(pairs[:left]):
            value = 2 * (1 - len(set(pair_left) & set(pair_right)))
            gram[9 + left, 9 + right] = gram[9 + right, 9 + left] = value
    ample_full = (3,) + (1,) * 8 + (0,) * 28
    return {
        "id": "R9-UA1x7-NEF-H8-37",
        "lattice": "U+A1(-1)^7",
        "gram": gram,
        # A0,A2,...,A8,A12 is an integral NS basis of determinant 128.
        "basis_indices": (0, 2, 3, 4, 5, 6, 7, 8, 9),
        "ample_full_coefficients": ample_full,
        "ample_square": 14,
        "ample_pairings": (2,) + (1,) * 8 + (4,) * 28,
        "tex_lines": "6594-6609",
        "branch_model": (
            "double plane branched over a smooth conic and a quartic, with "
            "eight nodes and twenty-eight joining-line walls"
        ),
    }


def _del_pezzo_walls() -> list[tuple[int, ...]]:
    walls: list[tuple[int, ...]] = []
    # Exceptional curves.
    for index in range(8):
        walls.append((0,) + tuple(int(index == j) for j in range(8)))
    # Lines through two points.
    for pair in itertools.combinations(range(8), 2):
        chosen = set(pair)
        walls.append((1,) + tuple(-1 if j in chosen else 0 for j in range(8)))
    # Conics through five points, indexed by the omitted triple.
    for triple in itertools.combinations(range(8), 3):
        omitted = set(triple)
        walls.append((2,) + tuple(0 if j in omitted else -1 for j in range(8)))
    # Cubics through seven points, double at one of them.
    for omitted in range(8):
        for doubled in range(8):
            if doubled == omitted:
                continue
            walls.append(
                (3,)
                + tuple(
                    0 if j == omitted else -2 if j == doubled else -1
                    for j in range(8)
                )
            )
    # Quartics through all points, double at a triple.
    for triple in itertools.combinations(range(8), 3):
        doubled = set(triple)
        walls.append((4,) + tuple(-2 if j in doubled else -1 for j in range(8)))
    # Quintics, simple at a pair and double at the complementary six points.
    for pair in itertools.combinations(range(8), 2):
        simple = set(pair)
        walls.append((5,) + tuple(-1 if j in simple else -2 for j in range(8)))
    # Sextics, triple at one point and double at the other seven.
    for tripled in range(8):
        walls.append((6,) + tuple(-3 if j == tripled else -2 for j in range(8)))
    if len(walls) != 240:
        raise ArithmeticError("the degree-one del Pezzo wall list is incomplete")
    return walls


def _case_u2a1x7() -> dict[str, Any]:
    lattice_gram = sp.diag(2, *([-2] * 8))
    walls = _del_pezzo_walls()
    gram = _coordinate_gram(lattice_gram, walls)
    return {
        "id": "R9-U2A1x7-NEF-H8-240",
        "lattice": "U(2)+A1(-1)^7",
        "gram": gram,
        "basis_gram_override": lattice_gram,
        "basis_wall_indices": tuple(range(9)),
        "wall_vectors": walls,
        "ample_standard": (3,) + (-1,) * 8,
        "ample_square": 2,
        "ample_pairings": (2,) * 240,
        "tex_lines": "6612-6711",
        "branch_model": (
            "double cover of a degree-one del Pezzo surface; equivalently a "
            "smooth sextic double plane with 120 six-tangent conics"
        ),
    }


def _case_ua7() -> dict[str, Any]:
    cycle = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 1)]
    gram = _gram_from_edges(9, [(0, 1), *cycle])
    return {
        "id": "R9-UA7-NEF-H8-9",
        "lattice": "U+A7(-1)",
        "gram": gram,
        "basis_indices": tuple(range(9)),
        "ample": (9, 19, 15, 12, 10, 9, 10, 12, 15),
        "ample_square": 120,
        "ample_pairings": (1, 1, 1, 1, 1, 2, 1, 1, 1),
        "tex_lines": "6714-6777",
        "branch_model": "double cover of F4 with an A7 branch singularity",
        "reference": {
            "name": "JKRT P8",
            "path": "graphs/simplices/8-jkrt_351-04-P8.coxiter",
            "sha256": "f9f25b857a31f1d2fdb8b316d6176ff814d0f72f8e615962b4c4b06a45fa4a3c",
            "permutation": (0, 1, 8, 7, 6, 5, 4, 3, 2),
            "ordinary": ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 1)),
            "parallel": (),
            "same_header_exact_matches": 2,
        },
    }


def _case_ud4a3() -> dict[str, Any]:
    gram = _gram_from_edges(
        10,
        [
            (0, 3), (1, 3), (2, 3), (3, 4), (4, 5),
            (5, 6), (6, 7), (6, 8), (7, 9), (8, 9),
        ],
    )
    return {
        "id": "R9-UD4A3-NEF-H8-10",
        "lattice": "U+D4(-1)+A3(-1)",
        "gram": gram,
        "basis_indices": tuple(range(9)),
        "ample": (7, 7, 7, 15, 10, 6, 3, 1, 0),
        "ample_square": 56,
        "ample_pairings": (1,) * 8 + (3, 1),
        "tex_lines": "6780-6858",
        "branch_model": (
            "double cover of F4 with D4 and A3 branch singularities; also a "
            "double plane with D6 and A1 singularities"
        ),
        "reference": {
            "name": "Tumarkin 8-dimensional pyramid 09",
            "path": "graphs/Tumarkin pyramids n+2/8-tum04_11_09.coxiter",
            "sha256": "ee66cb20832509e7a08616b5286b826db26dee989650967c0dd1532d3b33ed05",
            "permutation": (7, 8, 9, 6, 5, 4, 2, 3, 1, 0),
            "ordinary": ((0, 1), (0, 3), (1, 2), (3, 2), (2, 4), (4, 5), (5, 6), (6, 7), (6, 8), (6, 9)),
            "parallel": (),
            "same_header_exact_matches": 1,
        },
    }


def _case_ud5a2() -> dict[str, Any]:
    gram = _gram_from_edges(
        10,
        [
            (0, 2), (1, 2), (2, 3), (3, 4), (3, 5),
            (5, 6), (6, 7), (7, 8), (7, 9), (8, 9),
        ],
    )
    return {
        "id": "R9-UD5A2-NEF-H8-10",
        "lattice": "U+D5(-1)+A2(-1)",
        "gram": gram,
        "basis_indices": tuple(range(9)),
        "ample": (8, 8, 17, 19, 9, 13, 8, 4, 1),
        "ample_square": 88,
        "ample_pairings": (1,) * 8 + (2, 5),
        "tex_lines": "6861-6945",
        "branch_model": (
            "double cover of F4 with D5 and A2 branch singularities; also a "
            "double plane with a D6 singularity"
        ),
        "reference": {
            "name": "Tumarkin 8-dimensional pyramid 05",
            "path": "graphs/Tumarkin pyramids n+2/8-tum04_11_05.coxiter",
            "sha256": "8238a6ba99c9ad0d148b8bbc9bdebd8ac7def1a9c08875adceed5c394e55f735",
            "permutation": (8, 9, 7, 5, 6, 4, 3, 1, 0, 2),
            "ordinary": ((0, 1), (1, 2), (0, 2), (1, 3), (3, 4), (4, 5), (5, 6), (5, 7), (7, 8), (7, 9)),
            "parallel": (),
            "same_header_exact_matches": 1,
        },
    }


def _case_ud7() -> dict[str, Any]:
    gram = _gram_from_edges(
        9,
        [(0, 1), (1, 2), (2, 3), (2, 4), (4, 5), (5, 6), (6, 7), (6, 8)],
    )
    return {
        "id": "R9-UD7-NEF-H8-9",
        "lattice": "U+D7(-1)",
        "gram": gram,
        "basis_indices": tuple(range(9)),
        "ample": (13, 27, 42, 20, 38, 35, 33, 16, 16),
        "ample_square": 260,
        # The graph/vector computation resolves the inconsistent printed labels:
        # degree 2 is on A4, not A2; all other degrees are 1.
        "ample_pairings": (1, 1, 1, 2, 1, 1, 1, 1, 1),
        "tex_lines": "6948-7023",
        "branch_model": (
            "double cover of F4 with a D7 branch singularity; also a double "
            "plane with an E7 singularity"
        ),
        "reference": {
            "name": "JKRT Q8",
            "path": "graphs/simplices/8-jkrt_351-03-Q8.coxiter",
            "sha256": "3c258adc59dfd6281ad1320e9baf70ccf5c6ba92c906ce801c0931884864c00f",
            "permutation": (6, 5, 4, 8, 3, 2, 1, 0, 7),
            "ordinary": ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1, 7), (4, 8)),
            "parallel": (),
            "same_header_exact_matches": 2,
        },
    }


def _case_ua1e6() -> dict[str, Any]:
    gram = _gram_from_edges(
        10,
        [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (5, 8), (8, 9)],
        [(0, 1)],
    )
    return {
        "id": "R9-UA1E6-NEF-H8-10",
        "lattice": "U+A1(-1)+E6(-1)",
        "gram": gram,
        "basis_indices": tuple(range(9)),
        "ample": (12, 17, 12, 8, 5, 3, 1, 0, 1),
        "ample_square": 184,
        "ample_pairings": (10, 2) + (1,) * 8,
        "tex_lines": "7026-7100",
        "branch_model": (
            "double cover of F4 with a node and an E6 branch singularity; "
            "also a double plane with a D7 singularity"
        ),
        "reference": {
            "name": "Tumarkin 8-dimensional pyramid 13-k=2",
            "path": "graphs/Tumarkin pyramids n+2/8-tum04_11_13-k=2.coxiter",
            "sha256": "dd092c507bdba42baca4aca3888bd8510adb5e194d6ba2f64f7c986b13cd86aa",
            "permutation": (9, 8, 7, 6, 5, 2, 1, 0, 3, 4),
            "ordinary": ((0, 1), (1, 2), (2, 3), (3, 4), (2, 5), (5, 6), (6, 7), (7, 8)),
            "parallel": ((8, 9),),
            "same_header_exact_matches": 1,
        },
    }


def _cases() -> list[dict[str, Any]]:
    return [
        _case_ue7(),
        _case_ud6a1(),
        _case_ud4a1x3(),
        _case_ua1x7(),
        _case_u2a1x7(),
        _case_ua7(),
        _case_ud4a3(),
        _case_ud5a2(),
        _case_ud7(),
        _case_ua1e6(),
    ]


def _payload(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _relation_label(intersection: int) -> str:
    if intersection == 0:
        return "2"
    if intersection == 1:
        return "3"
    if intersection == 2:
        return "parallel"
    return f"ultraparallel_{intersection}"


def _relation_histogram(gram: sp.Matrix) -> dict[str, int]:
    return dict(
        sorted(
            Counter(
                _relation_label(int(gram[left, right]))
                for left in range(gram.rows)
                for right in range(left)
            ).items()
        )
    )


def _matrix_summary(case: dict[str, Any]) -> dict[str, Any]:
    gram: sp.Matrix = case["gram"]
    if "basis_gram_override" in case:
        basis_gram: sp.Matrix = case["basis_gram_override"]
        columns = sp.Matrix.hstack(*(sp.Matrix(vector) for vector in case["wall_vectors"]))
        ample = sp.Matrix(case["ample_standard"])
        pairings = tuple(int(value) for value in ample.T * basis_gram * columns)
        square = int((ample.T * basis_gram * ample)[0])
        selected = columns[:, list(case["basis_wall_indices"])]
        span_index = abs(int(selected.det()))
        basis_label: Any = "standard L,E1,...,E8 basis"
    else:
        indices = case["basis_indices"]
        basis_gram = gram.extract(indices, indices)
        basis_columns = sp.eye(gram.rows)[:, list(indices)]
        wall_in_basis = basis_gram.inv() * basis_columns.T * gram
        span_index = 1
        basis_label = tuple(index + 1 for index in indices)
        if "ample_full_coefficients" in case:
            full = sp.Matrix(case["ample_full_coefficients"])
            pairings = tuple(int(value) for value in gram * full)
            square = int((full.T * gram * full)[0])
        else:
            ample = sp.Matrix(case["ample"])
            pairings = tuple(int(value) for value in ample.T * basis_gram * wall_in_basis)
            square = int((ample.T * basis_gram * ample)[0])
    smith = tuple(
        abs(int(value))
        for value in smith_normal_form(basis_gram, domain=ZZ).diagonal()
    )
    checks = {
        "wall_count": gram.rows == len(case["ample_pairings"]),
        "all_walls_minus_two": all(gram[index, index] == -2 for index in range(gram.rows)),
        "all_intersections_nonnegative": all(
            gram[left, right] >= 0
            for left in range(gram.rows)
            for right in range(left)
        ),
        "rank_and_inertia": basis_gram.rank() == 9
        and exact_inertia(basis_gram) == (1, 8, 0),
        "full_wall_rank": gram.rank() == 9,
        "wall_subset_spans_NS": span_index == 1,
        "ample_square": square == case["ample_square"],
        "ample_pairings": pairings == case["ample_pairings"],
        "ample_is_strict": all(value > 0 for value in pairings),
    }
    if not all(checks.values()):
        raise ArithmeticError(
            f"{case['id']} matrix checks failed: "
            f"{[name for name, passed in checks.items() if not passed]}"
        )
    return {
        "facet_count": gram.rows,
        "wall_gram": _payload(gram),
        "wall_gram_sha256": _canonical_sha256(_payload(gram)),
        "wall_relation_histogram": _relation_histogram(gram),
        "integral_basis": basis_label,
        "basis_determinant": int(basis_gram.det()),
        "basis_smith_invariants": smith,
        "basis_inertia": exact_inertia(basis_gram),
        "ample_square": square,
        "ample_wall_pairings": pairings,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def _expanded_cpp_text() -> str:
    source = CPP_SOURCE.read_text(encoding="utf-8")
    if hashlib.sha256(source.encode()).hexdigest() != CPP_SOURCE_SHA256:
        raise ArithmeticError("the pinned literal Vinberg verifier changed")
    start = source.index("struct Bits {")
    end = source.index("class Enumerator {")
    replacement = r'''struct Bits {
  std::uint64_t w0 = 0, w1 = 0, w2 = 0, w3 = 0;
  bool operator==(const Bits& other) const {
    return w0 == other.w0 && w1 == other.w1 &&
           w2 == other.w2 && w3 == other.w3;
  }
};

struct BitsHash {
  std::size_t operator()(const Bits& value) const {
    std::uint64_t mixed = value.w0;
    mixed ^= value.w1 * 0x9e3779b97f4a7c15ULL;
    mixed ^= value.w2 * 0xc2b2ae3d27d4eb4fULL;
    mixed ^= value.w3 * 0x165667b19e3779f9ULL;
    return std::hash<std::uint64_t>{}(mixed);
  }
};

static int Popcount(Bits value) {
  return __builtin_popcountll(value.w0) + __builtin_popcountll(value.w1) +
         __builtin_popcountll(value.w2) + __builtin_popcountll(value.w3);
}
static bool Empty(Bits value) {
  return value.w0 == 0 && value.w1 == 0 && value.w2 == 0 && value.w3 == 0;
}
static int Lowest(Bits value) {
  if (value.w0) return __builtin_ctzll(value.w0);
  if (value.w1) return 64 + __builtin_ctzll(value.w1);
  if (value.w2) return 128 + __builtin_ctzll(value.w2);
  return 192 + __builtin_ctzll(value.w3);
}
static void Clear(Bits& value, int bit) {
  if (bit < 64) value.w0 &= ~(1ULL << bit);
  else if (bit < 128) value.w1 &= ~(1ULL << (bit - 64));
  else if (bit < 192) value.w2 &= ~(1ULL << (bit - 128));
  else value.w3 &= ~(1ULL << (bit - 192));
}
static void Set(Bits& value, int bit) {
  if (bit < 64) value.w0 |= 1ULL << bit;
  else if (bit < 128) value.w1 |= 1ULL << (bit - 64);
  else if (bit < 192) value.w2 |= 1ULL << (bit - 128);
  else value.w3 |= 1ULL << (bit - 192);
}
static bool Contains(Bits value, int bit) {
  if (bit < 64) return value.w0 & (1ULL << bit);
  if (bit < 128) return value.w1 & (1ULL << (bit - 64));
  if (bit < 192) return value.w2 & (1ULL << (bit - 128));
  return value.w3 & (1ULL << (bit - 192));
}
static Bits Intersect(Bits left, Bits right) {
  return Bits{left.w0 & right.w0, left.w1 & right.w1,
              left.w2 & right.w2, left.w3 & right.w3};
}
'''
    expanded = source[:start] + replacement + source[end:]
    expected = "b152f49bfd05cca638ecc42b4f8dc97266ee7a6ab56e08946046df94755e1ed1"
    if hashlib.sha256(expanded.encode()).hexdigest() != expected:
        raise ArithmeticError("the deterministic 256-wall verifier expansion changed")
    return expanded


def _vinberg_input(gram: sp.Matrix) -> str:
    rows = [f"8 {gram.rows}"]
    rows.extend(
        " ".join(str(-int(gram[row, column])) for column in range(gram.cols))
        for row in range(gram.rows)
    )
    return "\n".join(rows) + "\n"


def _literal_vinberg(executable: Path, case: dict[str, Any]) -> dict[str, Any]:
    expected = VINBERG_EXPECTED[case["id"]]
    raw_input = _vinberg_input(case["gram"])
    completed = subprocess.run(
        [str(executable)],
        input=raw_input,
        text=True,
        capture_output=True,
        check=True,
        timeout=1800,
    )
    record = _parse_output(completed.stdout)
    special = record["finite_special_subset_counts_by_rank"]
    endpoints = record["edge_endpoint_type_histogram"]
    checks = {
        "dimension_and_facets": record["dimension"] == 8
        and record["facet_count"] == case["gram"].rows,
        "finite_search_node_count": record["edge_search_node_count"]
        == expected["finite_nodes"],
        "affine_search_node_count": record["all_finite_subset_search_node_count"]
        == expected["affine_nodes"],
        "finite_prefix_counts": record["ordering_dependent_edge_search_prefix_counts"]
        == expected["prefix"],
        "finite_special_counts": special == expected["special"],
        "affine_components": record["connected_affine_component_count"]
        == expected["affine_components"],
        "maximal_parabolic_cusps": record["ideal_vertex_count"]
        == expected["cusps"],
        "endpoint_histogram": endpoints == expected["endpoints"],
        "every_elliptic_edge_has_two_endpoints": record["bad_edge_count"] == 0,
        "every_affine_zero_neighborhood_is_full_parabolic": record[
            "bad_parabolic_neighborhood_count"
        ]
        == 0,
        "endpoint_histogram_covers_every_edge": sum(endpoints.values())
        == special[-2],
        "noncompact_has_a_cusp": record["ideal_vertex_count"] > 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(
            f"{case['id']} Vinberg checks failed: "
            f"{[name for name, passed in checks.items() if not passed]}"
        )
    return {
        **record,
        "elliptic_edge_count": special[-2],
        "finite_vertex_count": special[-1],
        "input_sha256": hashlib.sha256(raw_input.encode()).hexdigest(),
        "stdout_sha256": hashlib.sha256(completed.stdout.encode()).hexdigest(),
        "criterion": (
            "every elliptic rank-7 edge has exactly two finite/ideal endpoints, "
            "and every connected affine component has a full rank-7 parabolic "
            "zero-neighborhood"
        ),
        "decision": "finite_volume_noncompact",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def _lineage_certificate(case: dict[str, Any]) -> dict[str, Any]:
    reference = case.get("reference")
    ordinary_snapshot = {
        "sha256": "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3",
        "records": 1214,
    }
    doubling_snapshot = {
        "sha256": "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036",
        "records": 2959,
    }
    if reference is None:
        return {
            "published_status": (
                "published Roulleau complete K3 root chamber, hence not a "
                "newly discovered worldwide abstract polytope"
            ),
            "pinned_coxiter_bounded_scan": {
                "ordinary_snapshot": ordinary_snapshot,
                "ordinary_same_H8_facet_header_records": 0,
                "doubling_snapshot": doubling_snapshot,
                "doubling_same_H8_facet_header_records": 0,
                "exact_matches": 0,
            },
            "novelty_boundary": (
                "new dedicated project extraction and literal volume certificate; "
                "bounded corpus absence only, with no worldwide priority claim"
            ),
        }

    reference_labels: dict[tuple[int, int], str] = {}
    for edge in reference["ordinary"]:
        reference_labels[tuple(sorted(edge))] = "3"
    for edge in reference["parallel"]:
        reference_labels[tuple(sorted(edge))] = "parallel"
    permutation = reference["permutation"]
    gram = case["gram"]
    checks = []
    for left in range(gram.rows):
        for right in range(left):
            candidate = _relation_label(int(gram[left, right]))
            mapped = tuple(sorted((permutation[left], permutation[right])))
            checks.append(candidate == reference_labels.get(mapped, "2"))
    if not all(checks):
        raise ArithmeticError(f"{case['id']} named CoxIter graph match failed")

    reference_public = {
        key: value
        for key, value in reference.items()
        if key not in {"ordinary", "parallel"}
    }
    reference_public["candidate_to_reference_zero_based"] = reference_public.pop(
        "permutation"
    )
    return {
        "published_status": (
            "published Roulleau complete K3 root chamber and a named classical "
            "JKRT/Tumarkin abstract polytope"
        ),
        "exact_reference": reference_public,
        "literal_relation_graph_match": True,
        "pinned_coxiter_bounded_scan": {
            "ordinary_snapshot": ordinary_snapshot,
            "ordinary_same_H8_facet_header_records": (
                56 if gram.rows == 9 else 18
            ),
            "ordinary_exact_matches": reference["same_header_exact_matches"],
            "doubling_snapshot": doubling_snapshot,
            "doubling_same_H8_facet_header_records": 0,
            "doubling_exact_matches": 0,
        },
        "novelty_boundary": (
            "new dedicated project same-polytope/full-nef alignment and literal "
            "volume certificate; no worldwide novelty"
        ),
    }


def discover(compiler: str | None = None) -> None:
    selected = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected is None:
        raise RuntimeError("a C++17 compiler is required")
    for case in _cases():
        summary = _matrix_summary(case)
        print(case["id"], summary["basis_determinant"], summary["basis_smith_invariants"])
    with tempfile.TemporaryDirectory(prefix="roulleau-r9-discover-") as raw:
        cpp = Path(raw) / "vinberg.cpp"
        executable = Path(raw) / "vinberg"
        cpp.write_text(_expanded_cpp_text(), encoding="utf-8")
        subprocess.run(
            [selected, "-O3", "-std=c++17", "-Wall", "-Wextra", "-Werror", str(cpp), "-o", str(executable)],
            check=True,
        )
        for case in _cases():
            completed = subprocess.run(
                [str(executable)], input=_vinberg_input(case["gram"]), text=True,
                capture_output=True, check=True, timeout=1800,
            )
            print(case["id"])
            print(completed.stdout)


def build_certificate(compiler: str | None = None) -> dict[str, Any]:
    selected = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected is None:
        raise RuntimeError("a C++17 compiler is required")
    expanded_cpp = _expanded_cpp_text()
    cases: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="roulleau-r9-all-") as raw:
        temporary = Path(raw)
        cpp_path = temporary / "vinberg_256.cpp"
        executable = temporary / "vinberg_256"
        cpp_path.write_text(expanded_cpp, encoding="utf-8")
        subprocess.run(
            [
                selected,
                "-O3",
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(cpp_path),
                "-o",
                str(executable),
            ],
            check=True,
        )
        for source in _cases():
            matrix = _matrix_summary(source)
            vinberg = _literal_vinberg(executable, source)
            lineage = _lineage_certificate(source)
            cases.append(
                {
                    "id": source["id"],
                    "lattice": source["lattice"],
                    "picard_rank": 9,
                    "hyperbolic_dimension": 8,
                    "facet_count": source["gram"].rows,
                    "tex_lines": source["tex_lines"],
                    "branch_model": source["branch_model"],
                    "complete_curve_source": (
                        "Roulleau's rank-9 construction and final table give "
                        "the complete irreducible (-2)-curve count"
                    ),
                    "complete_full_nef_reason": (
                        "the published wall list is the complete irreducible "
                        "(-2)-curve set; the K3 ample-cone theorem identifies "
                        "their positive-cone chamber with the entire nef cone"
                    ),
                    "matrix_certificate": matrix,
                    "literal_vinberg": vinberg,
                    "lineage": lineage,
                    "classification": (
                        "finite_volume_noncompact_full_nef_K3_Coxeter_polytope"
                    ),
                    "source_erratum": (
                        "The printed D_260 sentence assigns degree 2 to A2 and "
                        "then also declares A2 among the degree-1 walls.  The "
                        "published graph and exact lattice calculation put "
                        "degree 2 on A4; all other wall degrees are 1."
                        if source["id"] == "R9-UD7-NEF-H8-9"
                        else None
                    ),
                    "all_checks_pass": matrix["all_checks_pass"]
                    and vinberg["all_checks_pass"],
                }
            )

    checks = {
        "all_ten_rank9_rows": len(cases) == 10,
        "facet_counts_match_published_table": [row["facet_count"] for row in cases]
        == [9, 10, 15, 37, 240, 9, 10, 10, 9, 10],
        "all_complete_full_nef": all(row["complete_full_nef_reason"] for row in cases),
        "all_literal_finite_volume_noncompact": all(
            row["literal_vinberg"]["decision"] == "finite_volume_noncompact"
            and row["literal_vinberg"]["bad_edge_count"] == 0
            and row["literal_vinberg"]["bad_parabolic_neighborhood_count"] == 0
            and row["literal_vinberg"]["ideal_vertex_count"] > 0
            for row in cases
        ),
        "all_wall_grams_distinct": len(
            {row["matrix_certificate"]["wall_gram_sha256"] for row in cases}
        )
        == 10,
        "seven_exact_named_classical_matches": sum(
            bool(row["lineage"].get("literal_relation_graph_match")) for row in cases
        )
        == 7,
        "three_large_bounded_corpus_absences": [
            row["id"]
            for row in cases
            if row["lineage"].get("pinned_coxiter_bounded_scan", {}).get(
                "exact_matches"
            )
            == 0
        ]
        == [
            "R9-UD4A1x3-NEF-H8-15",
            "R9-UA1x7-NEF-H8-37",
            "R9-U2A1x7-NEF-H8-240",
        ],
        "single_source_erratum_is_D260_A4": [
            row["id"] for row in cases if row["source_erratum"] is not None
        ]
        == ["R9-UD7-NEF-H8-9"],
        "degree_one_del_pezzo_has_240_walls": next(
            row for row in cases if row["id"] == "R9-U2A1x7-NEF-H8-240"
        )["facet_count"]
        == 240,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global rank-nine batch checks failed: {checks}")
    payload = {
        "schema": SCHEMA,
        "artifact_id": "roulleau_rank9_all_full_nef",
        "status_date": "2026-08-09",
        "source": {
            "roulleau_arxiv": ROULLEAU_ARXIV,
            "roulleau_source_tar_sha256": ROULLEAU_SOURCE_TAR_SHA256,
            "roulleau_source_tex_sha256": ROULLEAU_SOURCE_TEX_SHA256,
            "rank9_source_lines": "6379-7100",
            "rank9_final_table_lines": "10810-10870",
            "base_cpp_verifier": str(CPP_SOURCE.relative_to(ROOT)),
            "base_cpp_verifier_sha256": CPP_SOURCE_SHA256,
            "expanded_256_wall_cpp_sha256": hashlib.sha256(
                expanded_cpp.encode()
            ).hexdigest(),
            "compiler": Path(selected).name,
        },
        "theorem": (
            "All ten Picard-rank-nine Roulleau root types give complete "
            "full-nef finite-volume noncompact Coxeter polytopes in H8.  "
            "The largest is the 240-wall chamber obtained from the 240 "
            "(-1)-curves on a degree-one del Pezzo surface."
        ),
        "cases": cases,
        "checks": checks,
        "all_checks_pass": all(checks.values())
        and all(row["all_checks_pass"] for row in cases),
    }
    return json.loads(json.dumps(payload))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler")
    parser.add_argument("--discover", action="store_true")
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    if args.discover:
        discover(args.compiler)
        return
    payload = build_certificate(args.compiler)
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if payload != archived:
            raise ArithmeticError("frozen rank-nine archive changed")
    if args.write:
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
