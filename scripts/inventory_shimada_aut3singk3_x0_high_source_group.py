#!/usr/bin/env python3
"""Lightweight exact high-subgroup checkpoint for Shimada's ``X_0`` stock.

This deliberately stops before restricting the 3,278 walls.  It pins and
replays the complete published source, ranks every matrix in the genuine
1,440-element ample stabilizer, and proves that only one involution class can
give a proper fixed section above ``H9``.  The full weighted-wall group,
irredundancy, Coxeter-angle, and volume gates remain pending.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import sympy as sp

from scripts.certify_shimada_aut3singk3_x2_rational_involution_h10_37 import (
    _common_fixed_rank,
    _rank,
)
from scripts.ingest_shimada_aut3singk3_x2_478_stock import (
    DEFAULT_SOURCE,
    DEFAULT_TEX,
    DEFAULT_ZIP,
    FILE_SHA256,
    TEX_SHA256,
    ZIP_SHA256,
    _assignment,
    _pair,
    _preserves_gram,
    _row_times_matrix,
    _sha256,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/shimada_aut3singk3_x0_high_source_group_checkpoint.json"
SCHEMA = "coxeter.shimada.aut3singk3.x0.high_source_group_checkpoint.v1"
SURFACE_INDEX = 0
IDENTITY_INDEX = 5
HIGH_REPRESENTATIVE = 3

EXPECTED_ORBIT_SIZES = [60, 40, 180, 10, 144, 144, 240, 360, 180, 240, 240, 720, 720]
EXPECTED_WALL_NORMS = [
    Fraction(-2), Fraction(-3, 2), Fraction(-4, 3), Fraction(-1),
    Fraction(-5, 6), Fraction(-5, 6), Fraction(-2, 3), Fraction(-2, 3),
    Fraction(-1, 3), Fraction(-1, 6), Fraction(-1, 6), Fraction(-1, 6),
    Fraction(-1, 6),
]
EXPECTED_AMPLE_PAIRINGS = [2, 3, 4, 4, 5, 5, 6, 6, 6, 7, 7, 7, 7]
EXPECTED_FIXED_RANK_HISTOGRAM = {
    2: 144,
    4: 584,
    6: 450,
    8: 170,
    10: 46,
    12: 45,
    20: 1,
}


def _matrix(value: Sequence[Sequence[Fraction]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(entry.numerator, entry.denominator) for entry in row] for row in value])


def build_checkpoint(source_dir: Path, source_zip: Path, tex_source: Path) -> dict[str, Any]:
    for filename, expected in FILE_SHA256.items():
        if _sha256(source_dir / filename) != expected:
            raise ArithmeticError(f"pinned source changed: {filename}")
    if _sha256(source_zip) != ZIP_SHA256 or _sha256(tex_source) != TEX_SHA256:
        raise ArithmeticError("pinned zip or paper source changed")

    gram_raw = _assignment(source_dir / "GramS.txt", "GramS", SURFACE_INDEX)
    ample_raw = _assignment(source_dir / "ample.txt", "ample", SURFACE_INDEX)
    wall_orbits = _assignment(source_dir / "walls.txt", "walls", SURFACE_INDEX)
    group_raw = _assignment(source_dir / "AutXa.txt", "AutXa", SURFACE_INDEX)
    gram = _matrix(gram_raw)
    ample = _matrix([ample_raw])
    if not (
        gram.shape == (20, 20)
        and gram == gram.T
        and gram.det() == -36
        and exact_inertia(gram) == (1, 19, 0)
        and _pair(ample_raw, gram_raw, ample_raw) == 20
    ):
        raise ArithmeticError("the X_0 lattice or ample class changed")
    if [len(orbit) for orbit in wall_orbits] != EXPECTED_ORBIT_SIZES:
        raise ArithmeticError("the complete X_0 wall denominator changed")
    for orbit_index, orbit in enumerate(wall_orbits):
        for wall in orbit:
            if _pair(wall, gram_raw, wall) != EXPECTED_WALL_NORMS[orbit_index]:
                raise ArithmeticError(f"wall norm changed in orbit {orbit_index + 1}")
            if _pair(ample_raw, gram_raw, wall) != EXPECTED_AMPLE_PAIRINGS[orbit_index]:
                raise ArithmeticError(f"ample pairing changed in orbit {orbit_index + 1}")

    group = [_matrix(matrix) for matrix in group_raw]
    identity = sp.eye(20)
    if len(group) != 1440 or len({tuple(matrix) for matrix in group}) != 1440:
        raise ArithmeticError("the X_0 stabilizer denominator changed")
    if [index for index, matrix in enumerate(group) if matrix == identity] != [IDENTITY_INDEX]:
        raise ArithmeticError("the X_0 identity index changed")
    gram_int = [[int(value) for value in row] for row in gram_raw]
    for index, (raw, matrix) in enumerate(zip(group_raw, group, strict=True)):
        if any(value.denominator != 1 for row in raw for value in row):
            raise ArithmeticError(f"group matrix {index} is nonintegral")
        if _row_times_matrix(ample_raw, raw) != ample_raw:
            raise ArithmeticError(f"group matrix {index} does not fix a_0")
        if not _preserves_gram([[int(value) for value in row] for row in raw], gram_int):
            raise ArithmeticError(f"group matrix {index} does not preserve GramS[0]")

    fixed_ranks = [20 - _rank(matrix.T - identity) for matrix in group]
    histogram = Counter(fixed_ranks)
    if dict(sorted(histogram.items())) != EXPECTED_FIXED_RANK_HISTOGRAM:
        raise ArithmeticError("the X_0 fixed-rank histogram changed")
    arrays = np.asarray([matrix.tolist() for matrix in group], dtype=np.int64)
    maximum_entry = int(np.max(np.abs(arrays)))
    product_bound = 20 * maximum_entry**2
    if product_bound >= 2**63:
        raise ArithmeticError("the exact int64 multiplication bound failed")
    index_of = {
        tuple(int(value) for value in arrays[index].ravel()): index
        for index in range(len(group))
    }
    high = [
        index
        for index, rank in enumerate(fixed_ranks)
        if rank >= 11 and index != IDENTITY_INDEX
    ]
    identity_array = np.eye(20, dtype=np.int64)
    if not (
        len(high) == 45
        and all(fixed_ranks[index] == 12 for index in high)
        and all(np.array_equal(arrays[index] @ arrays[index], identity_array) for index in high)
    ):
        raise ArithmeticError("the X_0 high involution class changed")
    representative = arrays[HIGH_REPRESENTATIVE]
    centralizer_order = sum(
        np.array_equal(matrix @ representative, representative @ matrix)
        for matrix in arrays
    )
    if centralizer_order != 32 or 1440 // centralizer_order != 45:
        raise ArithmeticError("the X_0 high conjugacy class changed")

    high_set = set(high)
    v4_groups: set[frozenset[int]] = set()
    for position, left in enumerate(high):
        for right in high[position + 1 :]:
            product = arrays[left] @ arrays[right]
            product_index = index_of.get(tuple(int(value) for value in product.ravel()))
            if product_index is None:
                raise ArithmeticError("the X_0 list is not closed on a high pair")
            if product_index not in high_set:
                continue
            if not np.array_equal(product, arrays[right] @ arrays[left]):
                raise ArithmeticError("a high involutive product failed to commute")
            subgroup = frozenset((IDENTITY_INDEX, left, right, product_index))
            if len(subgroup) == 4:
                v4_groups.add(subgroup)
    v4_ranks = [
        _common_fixed_rank([group[index] for index in sorted(subgroup - {IDENTITY_INDEX})])
        for subgroup in v4_groups
    ]
    if len(v4_groups) != 30 or Counter(v4_ranks) != {8: 30}:
        raise ArithmeticError("the X_0 high V4 closure changed")

    return {
        "schema": SCHEMA,
        "id": "SHIMADA-AUT3SINGK3-X0-HIGH-SOURCE-GROUP-CHECKPOINT",
        "sources": {
            "paper_url": "https://arxiv.org/abs/1412.6904",
            "supplement_zip_url": "https://www.math.sci.hiroshima-u.ac.jp/shimada/preprints/Aut3SingK3s/compdata/compdata.zip",
            "supplement_zip_sha256": ZIP_SHA256,
            "paper_tex_sha256": TEX_SHA256,
            "supplement_file_sha256": FILE_SHA256,
        },
        "complete_source_stock": {
            "surface": "X_0",
            "transcendental_gram": [[6, 0], [0, 6]],
            "neron_severi_rank": 20,
            "neron_severi_determinant": -36,
            "ample_square": 20,
            "wall_count": sum(EXPECTED_ORBIT_SIZES),
            "wall_orbit_sizes": EXPECTED_ORBIT_SIZES,
            "wall_norms": [str(value) for value in EXPECTED_WALL_NORMS],
            "ample_pairings": EXPECTED_AMPLE_PAIRINGS,
            "genuine_ample_stabilizer": "Aut(X_0,a_0), order 1440",
        },
        "high_source_group_census": {
            "fixed_rank_histogram": {str(key): value for key, value in sorted(histogram.items())},
            "unique_high_class": {
                "representative_zero_based": HIGH_REPRESENTATIVE,
                "cardinality": 45,
                "centralizer_order": centralizer_order,
                "fixed_rank": 12,
                "hyperbolic_dimension": 11,
            },
            "high_closed_V4_count": len(v4_groups),
            "high_closed_V4_common_fixed_rank_histogram": {"8": 30},
            "only_C2_can_be_a_proper_fixed_rank_at_least_11_subgroup": True,
            "proof": (
                "Every element of a subgroup with common fixed rank at least 11 must "
                "be one of the 45 rank-12 involutions.  They form one class by the "
                "order-32 centralizer.  Two distinct high elements whose product is "
                "high generate one of the 30 exhaustively enumerated V4 groups, all "
                "of common fixed rank 8.  Hence the unique high C2 class is the only "
                "proper source-group candidate above H9."
            ),
        },
        "pending_gates": {
            "full_abstract_oriented_weighted_wall_group": "pending",
            "restriction_of_all_3278_walls": "pending",
            "exact_irredundancy": "pending",
            "Coxeter_pair_angle_gate": "pending",
            "Vinberg_volume_gate": "not authorized before a Coxeter survivor",
        },
        "exact_int64_audit": {
            "maximum_matrix_entry_absolute": maximum_entry,
            "single_product_entry_absolute_bound": product_bound,
            "below_signed_int64_limit": True,
        },
        "conclusion": "checkpoint only; X_0 remains the next unmined complete stock",
        "all_checks_pass": True,
    }


def validate_checkpoint(payload: dict[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or payload.get("all_checks_pass") is not True:
        raise ArithmeticError("wrong or failed X_0 checkpoint")
    if payload["complete_source_stock"]["wall_count"] != 3278:
        raise ArithmeticError("the X_0 wall count changed")
    if payload["pending_gates"]["Coxeter_pair_angle_gate"] != "pending":
        raise ArithmeticError("the X_0 checkpoint overclaims completion")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--source-zip", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--tex-source", type=Path, default=DEFAULT_TEX)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_checkpoint(args.source_dir, args.source_zip, args.tex_source)
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != encoded:
            raise ArithmeticError(f"archive mismatch: {OUTPUT}")
        print(f"archive matches: {OUTPUT}")
        return
    OUTPUT.write_text(encoded, encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "walls": 3278,
        "unique_high_class": "H11/C2",
        "restriction_and_angle_gate": "pending",
    }, indent=2))


if __name__ == "__main__":
    main()
