#!/usr/bin/env python3
"""Certify the seven remaining rank-ten full-nef K3 Coxeter chambers.

This file starts from the literal curve configurations in Roulleau's rank-ten
atlas.  It reconstructs their wall Gram matrices in the published divisor
bases.  The final certificate also runs Vinberg's edge/cusp criterion; the
matrix reconstruction is kept here independently so that every wall and
ample-class assertion is exactly replayable.
"""

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
OUTPUT = (
    ROOT
    / "data/enumerations/roulleau_rank10_remaining_full_nef.json"
)
SCHEMA = "coxeter.roulleau_rank10_remaining_full_nef.v1"
CPP_SOURCE = (
    ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
)
CPP_SOURCE_SHA256 = (
    "45174f5461d28affc1fdb4cb55047b7a1581037d307e32800651badb3c2df045"
)
ROULLEAU_ARXIV = "https://arxiv.org/abs/2003.08985"
ROULLEAU_SOURCE_TAR_SHA256 = (
    "3e445ed01446b47ded6c007b395308de34e566b736f0ba982470bc497adeb1e9"
)
ROULLEAU_SOURCE_TEX_SHA256 = (
    "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
)


VINBERG_EXPECTED: dict[str, dict[str, Any]] = {
    "R10-UE7A1-NEF-H9-11": {
        "finite_nodes": 635,
        "affine_nodes": 1529,
        "prefix": [11, 33, 61, 89, 111, 121, 113, 80, 15],
        "special": [11, 54, 156, 294, 378, 336, 204, 80, 15],
        "affine_components": 3,
        "cusps": 2,
        "endpoints": {"0,2": 1, "1,1": 23, "2,0": 56},
    },
    "R10-UD4D4-NEF-H9-11": {
        "finite_nodes": 912,
        "affine_nodes": 1913,
        "prefix": [11, 34, 70, 115, 160, 190, 190, 125, 16],
        "special": [11, 55, 165, 330, 460, 450, 300, 125, 16],
        "affine_components": 11,
        "cusps": 10,
        "endpoints": {"0,2": 27, "1,1": 52, "2,0": 46},
    },
    "R10-UD6A1x2-NEF-H9-14": {
        "finite_nodes": 2310,
        "affine_nodes": 3536,
        "prefix": [14, 68, 200, 360, 486, 519, 394, 231, 37],
        "special": [14, 83, 280, 602, 868, 854, 566, 231, 37],
        "affine_components": 15,
        "cusps": 7,
        "endpoints": {"0,2": 15, "1,1": 99, "2,0": 117},
    },
    "R10-U2D4D4-NEF-H9-18": {
        "finite_nodes": 13162,
        "affine_nodes": 19526,
        "prefix": [18, 124, 488, 1293, 2451, 3371, 3390, 1840, 186],
        "special": [18, 145, 688, 2128, 4340, 5684, 4496, 1840, 186],
        "affine_components": 148,
        "cusps": 71,
        "endpoints": {"0,2": 630, "1,1": 746, "2,0": 464},
    },
    "R10-UD4A1x4-NEF-H9-27": {
        "finite_nodes": 15865,
        "affine_nodes": 25430,
        "prefix": [27, 199, 759, 1830, 3096, 3851, 3544, 2245, 313],
        "special": [27, 247, 1157, 3250, 5876, 6994, 5320, 2245, 313],
        "affine_components": 150,
        "cusps": 45,
        "endpoints": {"0,2": 324, "1,1": 1025, "2,0": 896},
    },
    "R10-UA1x8-NEF-H9-145": {
        "finite_nodes": 832438,
        "affine_nodes": 1186914,
        "prefix": [145, 2995, 21101, 75752, 164500, 223720, 206880, 123520, 13824],
        "special": [145, 3072, 24880, 103712, 250880, 363776, 303104, 123520, 13824],
        "affine_components": 6632,
        "cusps": 1249,
        "endpoints": {"0,2": 28800, "1,1": 65024, "2,0": 29696},
    },
    "R10-UA2E6-NEF-H9-11": {
        "finite_nodes": 913,
        "affine_nodes": 1774,
        "prefix": [11, 34, 70, 115, 161, 196, 203, 105, 17],
        "special": [11, 55, 164, 322, 434, 406, 259, 105, 17],
        "affine_components": 6,
        "cusps": 5,
        "endpoints": {"0,2": 8, "1,1": 41, "2,0": 56},
    },
}


REFERENCE_RELATIONS: dict[str, dict[str, tuple[tuple[int, int], ...]]] = {
    "graphs/Tumarkin pyramids n+2/9-tum04_11_09-k=2.coxiter": {
        "ordinary": (
            (0, 1), (1, 2), (2, 3), (3, 4), (3, 5),
            (5, 6), (6, 7), (7, 8), (8, 9),
        ),
        "parallel": ((9, 10),),
    },
    "graphs/Tumarkin pyramids n+2/9-tum04_11_06.coxiter": {
        "ordinary": (
            (0, 3), (1, 3), (2, 3), (3, 4), (4, 5),
            (5, 6), (6, 7), (7, 8), (7, 9), (7, 10),
        ),
        "parallel": (),
    },
    "graphs/Tumarkin pyramids n+2/9-tum04_11_07.coxiter": {
        "ordinary": (
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 5),
            (5, 6), (6, 7), (7, 8), (8, 6), (2, 9), (9, 10),
        ),
        "parallel": (),
    },
}


def _lattice_gram_u_root(rank: int, u_scale: int = 1) -> sp.Matrix:
    gram = sp.zeros(rank)
    gram[0, 1] = gram[1, 0] = u_scale
    for index in range(2, rank):
        gram[index, index] = -2
    return gram


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


def _case_ue7a1() -> dict[str, Any]:
    ordinary = [
        (1, 2), (2, 3), (3, 4), (4, 5), (5, 6),
        (6, 7), (6, 8), (8, 9), (9, 10),
    ]
    gram = _gram_from_edges(11, ordinary, [(0, 1)])
    return {
        "id": "R10-UE7A1-NEF-H9-11",
        "lattice": "U+E7(-1)+A1(-1)",
        "gram": gram,
        "basis_indices": tuple(range(10)),
        "ample": (15, 24, 19, 15, 12, 10, 9, 4, 5, 2),
        "ample_square": 370,
        "ample_pairings": (18,) + (1,) * 9 + (2,),
        "tex_lines": "7238-7300",
        "branch_model": (
            "double cover of F4 branched over the negative section and a "
            "curve in |3s+12f| with one node and one E7 singularity"
        ),
        "reference": {
            "name": "Tumarkin 9-dimensional pyramid with 11 facets",
            "path": "graphs/Tumarkin pyramids n+2/9-tum04_11_09-k=2.coxiter",
            "sha256": "80120b1ad7e06e5129c327cf40cd959842d897b39afbd14d6b5bacc2392d1e98",
            "candidate_to_reference_zero_based": (10, 9, 8, 7, 6, 5, 3, 4, 2, 1, 0),
        },
    }


def _case_ud4d4() -> dict[str, Any]:
    ordinary = [
        (0, 1), (1, 2), (1, 3), (1, 4), (4, 5),
        (5, 6), (6, 7), (7, 8), (7, 9), (7, 10),
    ]
    gram = _gram_from_edges(11, ordinary)
    return {
        "id": "R10-UD4D4-NEF-H9-11",
        "lattice": "U+D4(-1)^2",
        "gram": gram,
        "basis_indices": tuple(range(10)),
        "ample": (7, 15, 7, 7, 10, 6, 3, 1, 0, 0),
        "ample_square": 56,
        "ample_pairings": (1,) * 11,
        "tex_lines": "7303-7389",
        "branch_model": (
            "double cover of F4 with two D4 branch singularities; equivalently "
            "a double plane branched over a sextic with three nodes and a D6 singularity"
        ),
        "reference": {
            "name": "Tumarkin 9-dimensional pyramid with 11 facets",
            "path": "graphs/Tumarkin pyramids n+2/9-tum04_11_06.coxiter",
            "sha256": "32326979a72e9bdb5d3249da4e2b6ccf363e8c6a2502686569e531c841123ef8",
            "candidate_to_reference_zero_based": (0, 3, 1, 2, 4, 5, 6, 7, 8, 9, 10),
        },
    }


def _case_ud6a1a1() -> dict[str, Any]:
    ordinary = [
        (4, 3), (5, 3), (8, 7), (9, 7),
        (3, 2), (7, 6), (2, 1), (1, 6), (1, 0),
    ]
    parallel = [
        (11, 12), (13, 10), (4, 10), (5, 11),
        (8, 12), (9, 13), (10, 12), (11, 13),
    ]
    gram = _gram_from_edges(14, ordinary, parallel)
    return {
        "id": "R10-UD6A1x2-NEF-H9-14",
        "lattice": "U+D6(-1)+A1(-1)^2",
        "gram": gram,
        "basis_indices": tuple(range(10)),
        "ample": (7, 16, 13, 11, 5, 5, 13, 11, 5, 5),
        "ample_square": 98,
        "ample_pairings": (2,) + (1,) * 9 + (10,) * 4,
        "tex_lines": "7392-7481",
        "branch_model": (
            "double plane branched over a line and a nodal quintic; the two "
            "line/quintic crossings carry D4 singularities"
        ),
    }


def _case_u2d4d4() -> dict[str, Any]:
    ordinary = []
    parallel = []
    for offset in range(8):
        top = 1 + offset
        bottom = 9 + offset
        ordinary.extend([(0, top), (17, bottom)])
        parallel.append((top, bottom))
    gram = _gram_from_edges(18, ordinary, parallel)
    return {
        "id": "R10-U2D4D4-NEF-H9-18",
        "lattice": "U(2)+D4(-1)^2",
        "gram": gram,
        "basis_indices": tuple(list(range(8)) + [9, 17]),
        "ample": (1, 3, 0, 0, 0, 0, 0, 0, 3, 1),
        "ample_square": 8,
        "ample_pairings": (1,) * 18,
        "tex_lines": "7484-7626",
        "branch_model": (
            "double cover of F2 branched over a smooth genus-9 curve in "
            "|4s+8f| having eight even-contact fibers"
        ),
    }


def _case_ud4a1x4() -> dict[str, Any]:
    basis_gram = _gram_from_edges(
        10,
        [
            (4, 0), (4, 1), (4, 2), (4, 3), (4, 5),
            (5, 6), (6, 7), (6, 8), (6, 9),
        ],
    )
    unit = [tuple(int(i == j) for i in range(10)) for j in range(10)]
    vectors: list[tuple[int, ...]] = list(unit)
    a11 = (1, 1, 1, 1, 2, 0, -2, -1, -1, -1)
    vectors.append(a11)
    f5 = sp.Matrix((0, 0, 0, 0, 0, 1, 2, 1, 1, 1))
    for index in range(4):
        vectors.append(tuple(int(value) for value in f5 - sp.Matrix(unit[index])))
    fibers = [
        (0, 1, 1, 1, 2, 1, 0, 0, 0, 0),
        (1, 0, 1, 1, 2, 1, 0, 0, 0, 0),
        (1, 1, 0, 1, 2, 1, 0, 0, 0, 0),
        (1, 1, 1, 0, 2, 1, 0, 0, 0, 0),
    ]
    for fiber in fibers:
        for index in (7, 8, 9):
            vectors.append(
                tuple(int(value) for value in sp.Matrix(fiber) - sp.Matrix(unit[index]))
            )
    gram = _coordinate_gram(basis_gram, vectors)
    return {
        "id": "R10-UD4A1x4-NEF-H9-27",
        "lattice": "U+D4(-1)+A1(-1)^4",
        "gram": gram,
        "basis_indices": tuple(range(10)),
        "ample": (3, 3, 3, 3, 7, 3, 1, 0, 0, 0),
        "ample_square": 26,
        # The literal graph and divisor coordinates give degree 2 on A6,
        # not A5.  The prose immediately after D_26 swaps these two labels.
        "ample_pairings": (1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1) + (6,) * 16,
        "tex_lines": "7629-7761",
        "branch_model": (
            "double plane branched over two lines and a quartic; sixteen walls "
            "map to the lines joining the two four-point line/quartic intersections"
        ),
    }


def _case_ua1x8() -> dict[str, Any]:
    lattice_gram = _lattice_gram_u_root(10)
    vectors: list[tuple[int, ...]] = []
    vectors.append((-1, 1) + (0,) * 8)  # A0
    for index in range(8):
        vectors.append((1, 0) + tuple(-int(index == j) for j in range(8)))
    fiber = sp.Matrix((4, 2) + (-1,) * 8)
    for index in range(8):
        vectors.append(
            tuple(int(value) for value in fiber - sp.Matrix(vectors[1 + index]))
        )
    for subset in itertools.combinations(range(8), 3):
        chosen = set(subset)
        vectors.append((4, 4) + tuple(-2 if j in chosen else -1 for j in range(8)))
    for subset in itertools.combinations(range(8), 3):
        chosen = set(subset)
        vectors.append((2, 2) + tuple(0 if j in chosen else -1 for j in range(8)))
    for index in range(8):
        vectors.append((6, 6) + tuple(-3 if j == index else -2 for j in range(8)))
    for index in range(8):
        vectors.append((0, 0) + tuple(int(index == j) for j in range(8)))
    gram = _coordinate_gram(lattice_gram, vectors)
    return {
        "id": "R10-UA1x8-NEF-H9-145",
        "lattice": "U+A1(-1)^8",
        "gram": gram,
        "basis_vectors": tuple([0, *range(1, 9), 144]),
        "basis_indices": tuple([0, *range(1, 9), 144]),
        "ample_in_standard_basis": (7, 5) + (-2,) * 8,
        "ample_square": 6,
        "ample_pairings": (2,) + (1,) * 16 + (4,) * 128,
        "tex_lines": "7764-7969",
        "branch_model": (
            "double plane branched over a nodal sextic with eight even-contact "
            "lines and sixty-four six-tangent conics"
        ),
        "standard_lattice_gram": lattice_gram,
        "wall_vectors": vectors,
    }


def _case_ua2e6() -> dict[str, Any]:
    ordinary = [
        (0, 1), (1, 2), (4, 3), (3, 2),
        (2, 5), (5, 6), (6, 7), (7, 8),
        (8, 9), (8, 10), (9, 10),
    ]
    gram = _gram_from_edges(11, ordinary)
    return {
        "id": "R10-UA2E6-NEF-H9-11",
        "lattice": "U+A2(-1)+E6(-1)",
        "gram": gram,
        "basis_indices": tuple(range(10)),
        "ample": (10, 21, 33, 21, 10, 25, 18, 12, 7, 3),
        "ample_square": 160,
        "ample_pairings": (1,) * 10 + (10,),
        "tex_lines": "7972-8050",
        "branch_model": (
            "double cover of F4 branched over the negative section and a "
            "curve in |3s+12f| with A2 and E6 singularities"
        ),
        "reference": {
            "name": "Tumarkin 9-dimensional pyramid with 11 facets",
            "path": "graphs/Tumarkin pyramids n+2/9-tum04_11_07.coxiter",
            "sha256": "45655e8a4522a015c7e9dd339dcf58d923c9788c66444ea97b197da9597f6830",
            "candidate_to_reference_zero_based": (0, 1, 2, 9, 10, 3, 4, 5, 6, 7, 8),
        },
    }


def _cases() -> list[dict[str, Any]]:
    return [
        _case_ue7a1(),
        _case_ud4d4(),
        _case_ud6a1a1(),
        _case_u2d4d4(),
        _case_ud4a1x4(),
        _case_ua1x8(),
        _case_ua2e6(),
    ]


def _payload(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _expanded_cpp_text() -> str:
    """Expand the pinned verifier's two-word mask to three words.

    The original verifier handles at most 128 facets.  Only the bit-container
    layer changes; all finite/affine classification and Vinberg logic remains
    byte-for-byte the pinned implementation.
    """
    source = CPP_SOURCE.read_text(encoding="utf-8")
    if hashlib.sha256(source.encode()).hexdigest() != CPP_SOURCE_SHA256:
        raise ArithmeticError("the pinned literal Vinberg verifier changed")
    start = source.index("struct Bits {")
    end = source.index("class Enumerator {")
    replacement = r'''struct Bits {
  std::uint64_t lo = 0, mid = 0, hi = 0;
  bool operator==(const Bits& other) const {
    return lo == other.lo && mid == other.mid && hi == other.hi;
  }
};

struct BitsHash {
  std::size_t operator()(const Bits& value) const {
    std::uint64_t mixed = value.lo;
    mixed ^= value.mid * 0x9e3779b97f4a7c15ULL;
    mixed ^= value.hi * 0xc2b2ae3d27d4eb4fULL;
    return std::hash<std::uint64_t>{}(mixed);
  }
};

static int Popcount(Bits value) {
  return __builtin_popcountll(value.lo) + __builtin_popcountll(value.mid) +
         __builtin_popcountll(value.hi);
}
static bool Empty(Bits value) {
  return value.lo == 0 && value.mid == 0 && value.hi == 0;
}
static int Lowest(Bits value) {
  if (value.lo) return __builtin_ctzll(value.lo);
  if (value.mid) return 64 + __builtin_ctzll(value.mid);
  return 128 + __builtin_ctzll(value.hi);
}
static void Clear(Bits& value, int bit) {
  if (bit < 64) value.lo &= ~(1ULL << bit);
  else if (bit < 128) value.mid &= ~(1ULL << (bit - 64));
  else value.hi &= ~(1ULL << (bit - 128));
}
static void Set(Bits& value, int bit) {
  if (bit < 64) value.lo |= 1ULL << bit;
  else if (bit < 128) value.mid |= 1ULL << (bit - 64);
  else value.hi |= 1ULL << (bit - 128);
}
static bool Contains(Bits value, int bit) {
  if (bit < 64) return value.lo & (1ULL << bit);
  if (bit < 128) return value.mid & (1ULL << (bit - 64));
  return value.hi & (1ULL << (bit - 128));
}
static Bits Intersect(Bits left, Bits right) {
  return Bits{left.lo & right.lo, left.mid & right.mid, left.hi & right.hi};
}
'''
    expanded = source[:start] + replacement + source[end:]
    expected = "9e228e762d74e5619aa6af0990e7e25827fda7757904a7483b9fdacb619139fe"
    if hashlib.sha256(expanded.encode()).hexdigest() != expected:
        raise ArithmeticError("the deterministic 192-wall verifier expansion changed")
    return expanded


def _vinberg_input(gram: sp.Matrix) -> str:
    rows = [f"9 {gram.rows}"]
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
    )
    record = _parse_output(completed.stdout)
    special = record["finite_special_subset_counts_by_rank"]
    endpoints = record["edge_endpoint_type_histogram"]
    checks = {
        "dimension_and_facets": record["dimension"] == 9
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
            "every elliptic rank-8 edge has exactly two finite/ideal endpoints, "
            "and every connected affine component has a full rank-8 parabolic "
            "zero-neighborhood"
        ),
        "decision": "finite_volume_noncompact",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


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


def _lineage_certificate(case: dict[str, Any]) -> dict[str, Any]:
    reference = case.get("reference")
    if reference is None:
        facets = case["gram"].rows
        ordinary_headers = 0
        doubling_headers = 8 if facets == 27 else 0
        return {
            "published_status": (
                "published Roulleau complete K3 root chamber, hence not a "
                "newly discovered worldwide abstract polytope"
            ),
            "pinned_coxiter_bounded_scan": {
                "ordinary_snapshot_sha256": (
                    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
                ),
                "ordinary_records": 1214,
                "same_H9_facet_header_records": ordinary_headers,
                "doubling_snapshot_sha256": (
                    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
                ),
                "doubling_records": 2959,
                "same_H9_facet_header_records_doubling": doubling_headers,
                "same_relation_histogram_records_doubling": 0,
                "exact_matches": 0,
            },
            "novelty_boundary": (
                "new dedicated project extraction and literal volume certificate; "
                "no worldwide priority claim"
            ),
        }

    relations = REFERENCE_RELATIONS[reference["path"]]
    reference_labels: dict[tuple[int, int], str] = {}
    for edge in relations["ordinary"]:
        reference_labels[tuple(sorted(edge))] = "3"
    for edge in relations["parallel"]:
        reference_labels[tuple(sorted(edge))] = "parallel"
    permutation = reference["candidate_to_reference_zero_based"]
    checks = []
    gram = case["gram"]
    for left in range(gram.rows):
        for right in range(left):
            candidate = _relation_label(int(gram[left, right]))
            mapped = tuple(sorted((permutation[left], permutation[right])))
            target = reference_labels.get(mapped, "2")
            checks.append(candidate == target)
    if not all(checks):
        raise ArithmeticError(f"{case['id']} named Tumarkin graph match failed")
    return {
        "published_status": (
            "published Roulleau complete K3 root chamber and classical "
            "Tumarkin abstract polytope"
        ),
        "exact_reference": reference,
        "literal_relation_graph_match": True,
        "pinned_coxiter_ordinary_same_header_records": 10,
        "pinned_coxiter_ordinary_relation_histogram_matches": 1,
        "pinned_coxiter_ordinary_exact_matches": 1,
        "novelty_boundary": (
            "new dedicated project same-polytope/full-nef alignment and literal "
            "volume certificate; no worldwide novelty"
        ),
    }


def _matrix_summary(case: dict[str, Any]) -> dict[str, Any]:
    gram: sp.Matrix = case["gram"]
    basis_indices = case["basis_indices"]
    basis_gram = gram.extract(basis_indices, basis_indices)
    if case["id"] == "R10-UA1x8-NEF-H9-145":
        ample = sp.Matrix(case["ample_in_standard_basis"])
        wall_columns = sp.Matrix.hstack(
            *(sp.Matrix(vector) for vector in case["wall_vectors"])
        )
        pairings = tuple(
            int(value)
            for value in (ample.T * case["standard_lattice_gram"] * wall_columns)
        )
        square = int((ample.T * case["standard_lattice_gram"] * ample)[0])
    else:
        ample = sp.Matrix(case["ample"])
        basis_columns = sp.eye(gram.rows)[:, list(basis_indices)]
        wall_in_basis = basis_gram.inv() * basis_columns.T * gram
        pairings = tuple(int(value) for value in (ample.T * basis_gram * wall_in_basis))
        square = int((ample.T * basis_gram * ample)[0])
    smith = tuple(
        abs(int(value))
        for value in smith_normal_form(basis_gram, domain=ZZ).diagonal()
    )
    checks = {
        "wall_count": gram.rows == len(case["ample_pairings"]),
        "all_walls_minus_two": all(gram[index, index] == -2 for index in range(gram.rows)),
        "all_distinct_wall_intersections_nonnegative": all(
            gram[left, right] >= 0
            for left in range(gram.rows)
            for right in range(left)
        ),
        "basis_rank_and_inertia": basis_gram.rank() == 10
        and exact_inertia(basis_gram) == (1, 9, 0),
        "full_wall_rank": gram.rank() == 10,
        "published_ample_square": square == case["ample_square"],
        "published_ample_pairings": pairings == case["ample_pairings"],
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
        "basis_indices_one_based": tuple(index + 1 for index in basis_indices),
        "basis_determinant": int(basis_gram.det()),
        "basis_smith_invariants": smith,
        "basis_inertia": exact_inertia(basis_gram),
        "ample_square": square,
        "ample_wall_pairings": pairings,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def build_certificate(compiler: str | None = None) -> dict[str, Any]:
    selected_compiler = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected_compiler is None:
        raise RuntimeError("a C++17 compiler is required")
    expanded_cpp = _expanded_cpp_text()
    cases = []
    with tempfile.TemporaryDirectory(prefix="roulleau-r10-remaining-") as raw:
        temporary = Path(raw)
        cpp_path = temporary / "vinberg_192.cpp"
        executable = temporary / "vinberg_192"
        cpp_path.write_text(expanded_cpp, encoding="utf-8")
        subprocess.run(
            [
                selected_compiler,
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
                    "picard_rank": 10,
                    "hyperbolic_dimension": 9,
                    "facet_count": source["gram"].rows,
                    "tex_lines": source["tex_lines"],
                    "branch_model": source["branch_model"],
                    "complete_curve_source": (
                        "Roulleau's rank-10 construction and final table give "
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
                    "classification": "finite_volume_noncompact_full_nef_K3_Coxeter_polytope",
                    "source_erratum": (
                        "The printed graph and D_26 coordinates put degree 2 on "
                        "A6, not A5; all other printed degrees agree."
                        if source["id"] == "R10-UD4A1x4-NEF-H9-27"
                        else None
                    ),
                    "all_checks_pass": matrix["all_checks_pass"]
                    and vinberg["all_checks_pass"],
                }
            )

    checks = {
        "seven_remaining_rank10_rows": len(cases) == 7,
        "facet_counts_match_published_table": [row["facet_count"] for row in cases]
        == [11, 11, 14, 18, 27, 145, 11],
        "all_full_nef": all(row["complete_full_nef_reason"] for row in cases),
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
        == 7,
        "three_exact_Tumarkin_matches": sum(
            bool(row["lineage"].get("literal_relation_graph_match")) for row in cases
        )
        == 3,
        "single_source_erratum_is_A6_for_A5": [
            row["id"] for row in cases if row["source_erratum"] is not None
        ]
        == ["R10-UD4A1x4-NEF-H9-27"],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global rank-ten batch checks failed: {checks}")
    payload = {
        "schema": SCHEMA,
        "artifact_id": "roulleau_rank10_remaining_full_nef",
        "status_date": "2026-08-09",
        "source": {
            "roulleau_arxiv": ROULLEAU_ARXIV,
            "roulleau_source_tar_sha256": ROULLEAU_SOURCE_TAR_SHA256,
            "roulleau_source_tex_sha256": ROULLEAU_SOURCE_TEX_SHA256,
            "rank10_source_lines": "7238-8050",
            "rank10_final_table_lines": "10864-10904",
            "base_cpp_verifier": str(CPP_SOURCE.relative_to(ROOT)),
            "base_cpp_verifier_sha256": CPP_SOURCE_SHA256,
            "expanded_192_wall_cpp_sha256": hashlib.sha256(
                expanded_cpp.encode()
            ).hexdigest(),
            "compiler": Path(selected_compiler).name,
        },
        "theorem": (
            "The seven previously unexpanded rank-ten Roulleau root types "
            "give complete full-nef finite-volume noncompact Coxeter polytopes "
            "in H9.  Together with the separately certified U+E8 and U+D8 "
            "simplices, this closes all nine published rank-ten types positively."
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
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.compiler)
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if payload != archived:
            raise ArithmeticError("frozen remaining rank-ten archive changed")
    if args.write:
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
