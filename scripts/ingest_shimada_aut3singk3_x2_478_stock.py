#!/usr/bin/env python3
"""Ingest Shimada's complete 478-wall stock for the singular K3 surface X_2.

The source is the exact computer-data supplement to arXiv:1412.6904.  For
``X_2`` the paper proves that ``D^(0)`` is an induced chamber contained in the
nef cone, that its translates tessellate the nef cone, and that the 478
primitive outward walls listed in eight ``Aut(X_2,a_2)``-orbits are complete.
The same supplement lists every one of the 720 stabilizer matrices.

This ingester deliberately keeps the complete denominator.  Later fixed-
section searches therefore do not depend on a hand-selected subset of walls
or group elements.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
from typing import Any

import sympy as sp

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/shimada_aut3singk3_x2_478_wall_stock.json"
DEFAULT_SOURCE = Path("/private/tmp/Aut3SingK3s_compdata/compdata")
DEFAULT_ZIP = Path("/private/tmp/Aut3SingK3s_compdata.zip")
DEFAULT_TEX = Path("/private/tmp/arxiv1412_source_clean/Aut3SingK3s.tex")

SCHEMA = "coxeter.shimada.aut3singk3.x2.complete_478_wall_stock.v1"
ZIP_SHA256 = "b3b943ba22a87615142a23f036351ac340c1ad4630ad720af1fe8df55e99723e"
TEX_SHA256 = "e2b37de790fc7b9d46d7d434ba28b17951ddc3fa8d795c034daa03ecad559e22"
FILE_SHA256 = {
    "GramS.txt": "b9c829e2253bdb62312d0484e622804d40b5bcf7ca6e2f1a96487d5729c84c00",
    "ample.txt": "33e29f7c6d2fc1cdb63fd10e07a5cc3192f1378c3cbceca362edbc1497b498ee",
    "walls.txt": "d9e63603667d69c71f760b3f73a887391d324c24c0c04103e14c3f0e2480eb55",
    "AutXa.txt": "d0844407e8c041177655b99b2d6af064a9077a93fa73dba33cd7a07f41111e80",
    "Isom.txt": "3b2914f2047296522818e93be68ce5eba177e09a410bd2a4e4eb07034f4426d2",
}

EXPECTED_ORBIT_SIZES = [36, 12, 40, 90, 30, 30, 120, 120]
EXPECTED_WALL_NORMS = [
    Fraction(-2), Fraction(-4, 3), Fraction(-6, 5), Fraction(-4, 5),
    Fraction(-8, 15), Fraction(-8, 15), Fraction(-2, 15), Fraction(-2, 15),
]
EXPECTED_AMPLE_PAIRINGS = [1, 2, 3, 4, 4, 4, 5, 5]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _assignment(path: Path, name: str, index: int) -> Any:
    """Parse one Magma-style rational list assignment without using eval."""
    raw = path.read_text(encoding="utf-8")
    marker = f"{name}[{index}]:="
    start = raw.find(marker)
    if start < 0:
        raise ValueError(f"missing {marker} in {path}")
    start += len(marker)
    finish = raw.find(";", start)
    if finish < 0:
        raise ValueError(f"unterminated {marker} in {path}")
    tree = ast.parse(raw[start:finish].strip(), mode="eval")

    def visit(node: ast.AST) -> Any:
        if isinstance(node, ast.List):
            return [visit(value) for value in node.elts]
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return Fraction(node.value)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -visit(node.operand)
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
            return visit(node.left) / visit(node.right)
        raise ValueError(f"unsupported syntax in {path}: {ast.dump(node)}")

    return visit(tree.body)


def _plain(value: Any) -> Any:
    if isinstance(value, Fraction):
        return value.numerator if value.denominator == 1 else str(value)
    if isinstance(value, list):
        return [_plain(item) for item in value]
    raise TypeError(type(value))


def _fraction(value: Any) -> Fraction:
    return Fraction(value)


def _pair(left: list[Fraction], gram: list[list[Fraction]], right: list[Fraction]) -> Fraction:
    return sum(
        left[i] * gram[i][j] * right[j]
        for i in range(len(left))
        for j in range(len(right))
    )


def _row_times_matrix(row: list[Fraction], matrix: list[list[Fraction]]) -> list[Fraction]:
    return [sum(row[i] * matrix[i][j] for i in range(len(row))) for j in range(len(row))]


def _preserves_gram(matrix: list[list[int]], gram: list[list[int]]) -> bool:
    rank = len(gram)
    left = [
        [sum(matrix[i][k] * gram[k][j] for k in range(rank)) for j in range(rank)]
        for i in range(rank)
    ]
    transformed = [[sum(left[i][k] * matrix[j][k] for k in range(rank)) for j in range(rank)] for i in range(rank)]
    return transformed == gram


def build_stock(source_dir: Path, source_zip: Path, tex_source: Path) -> dict[str, Any]:
    for filename, expected in FILE_SHA256.items():
        path = source_dir / filename
        if _sha256(path) != expected:
            raise ArithmeticError(f"pinned source changed: {path}")
    if _sha256(source_zip) != ZIP_SHA256 or _sha256(tex_source) != TEX_SHA256:
        raise ArithmeticError("pinned supplement zip or paper TeX source changed")

    gram = _assignment(source_dir / "GramS.txt", "GramS", 2)
    ample = _assignment(source_dir / "ample.txt", "ample", 2)
    wall_orbits = _assignment(source_dir / "walls.txt", "walls", 2)
    automorphisms = _assignment(source_dir / "AutXa.txt", "AutXa", 2)

    gram_matrix = sp.Matrix([[sp.Rational(value.numerator, value.denominator) for value in row] for row in gram])
    if not (
        gram_matrix.shape == (20, 20)
        and gram_matrix == gram_matrix.T
        and gram_matrix.det() == -15
        and exact_inertia(gram_matrix) == (1, 19, 0)
    ):
        raise ArithmeticError("unexpected X_2 Neron--Severi lattice")
    if len(ample) != 20 or _pair(ample, gram, ample) != 12:
        raise ArithmeticError("unexpected ample class")
    if [len(orbit) for orbit in wall_orbits] != EXPECTED_ORBIT_SIZES:
        raise ArithmeticError("the complete wall-orbit census changed")
    for orbit_index, orbit in enumerate(wall_orbits):
        for wall in orbit:
            if len(wall) != 20:
                raise ArithmeticError("a wall has the wrong rank")
            if _pair(wall, gram, wall) != EXPECTED_WALL_NORMS[orbit_index]:
                raise ArithmeticError(f"wall norm changed in orbit {orbit_index + 1}")
            if _pair(ample, gram, wall) != EXPECTED_AMPLE_PAIRINGS[orbit_index]:
                raise ArithmeticError(f"ample pairing changed in orbit {orbit_index + 1}")

    identity = [[Fraction(int(i == j)) for j in range(20)] for i in range(20)]
    if len(automorphisms) != 720 or len({tuple(value for row in matrix for value in row) for matrix in automorphisms}) != 720:
        raise ArithmeticError("Aut(X_2,a_2) is not the expected 720-element exact list")
    identity_indices = [index for index, matrix in enumerate(automorphisms) if matrix == identity]
    if identity_indices != [599]:
        raise ArithmeticError(f"identity index changed: {identity_indices}")
    gram_int = [[int(value) for value in row] for row in gram]
    for index, matrix in enumerate(automorphisms):
        if len(matrix) != 20 or any(len(row) != 20 for row in matrix):
            raise ArithmeticError(f"automorphism {index} has the wrong shape")
        if any(value.denominator != 1 for row in matrix for value in row):
            raise ArithmeticError(f"automorphism {index} is not integral")
        if _row_times_matrix(ample, matrix) != ample:
            raise ArithmeticError(f"automorphism {index} does not fix a_2")
        matrix_int = [[int(value) for value in row] for row in matrix]
        if not _preserves_gram(matrix_int, gram_int):
            raise ArithmeticError(f"automorphism {index} does not preserve GramS[2]")

    orbit_records = []
    for index, orbit in enumerate(wall_orbits):
        orbit_records.append({
            "orbit_one_based": index + 1,
            "size": len(orbit),
            "wall_norm": _plain(EXPECTED_WALL_NORMS[index]),
            "ample_pairing": EXPECTED_AMPLE_PAIRINGS[index],
            "primitive_outward_wall_vectors_in_published_NS_basis": _plain(orbit),
        })
    return {
        "schema": SCHEMA,
        "id": "SHIMADA-AUT3SINGK3-X2-D0-478",
        "sources": {
            "paper": "Ichiro Shimada, The automorphism groups of certain singular K3 surfaces and an Enriques surface, arXiv:1412.6904",
            "paper_url": "https://arxiv.org/abs/1412.6904",
            "supplement_page": "https://www.math.sci.hiroshima-u.ac.jp/shimada/K3andEnriques.html",
            "supplement_zip_url": "https://www.math.sci.hiroshima-u.ac.jp/shimada/preprints/Aut3SingK3s/compdata/compdata.zip",
            "supplement_zip_sha256": ZIP_SHA256,
            "paper_tex_sha256": TEX_SHA256,
            "supplement_file_sha256": FILE_SHA256,
            "selected_surface_zero_based": 2,
            "theorem_scope": {
                "transcendental_lattice": "paper display defining X_2: Gram [[2,1],[1,8]]",
                "chamber_completeness": "Theorem thm:D0: D^(0) is an induced chamber in the nef cone and its Aut(X)-translates tessellate the nef cone",
                "wall_completeness": "Theorem thm:orbits and Table table:main: complete Delta(D^(0)) in eight Aut(X_2,a_2)-orbits",
                "stabilizer": "Theorem thm:aut: Aut(X_2,a_2) is isomorphic to PGL_2(F_9), of order 720",
            },
        },
        "surface": {
            "name": "X_2",
            "kind": "complex singular K3 surface",
            "transcendental_gram": [[2, 1], [1, 8]],
            "neron_severi_rank": 20,
            "neron_severi_signature": [1, 19],
            "neron_severi_determinant": -15,
            "published_neron_severi_gram": _plain(gram),
            "published_ample_vector_a2": _plain(ample),
            "ample_square": 12,
        },
        "complete_induced_chamber": {
            "name": "D^(0)",
            "contained_in_nef_cone": True,
            "translates_tessellate_nef_cone": True,
            "wall_count": sum(EXPECTED_ORBIT_SIZES),
            "wall_orbit_count": len(EXPECTED_ORBIT_SIZES),
            "wall_orbit_sizes": EXPECTED_ORBIT_SIZES,
            "wall_orbits": orbit_records,
        },
        "ample_stabilizer": {
            "abstract_group": "PGL_2(F_9)",
            "order": 720,
            "matrix_convention": "row vectors act on the right; M G M^T = G",
            "identity_zero_based_index": 599,
            "all_integral_20_by_20_matrices_in_supplement_order": _plain(automorphisms),
        },
        "checks": {
            "pinned_sources_match": True,
            "lorentzian_NS_rank20_determinant_minus15": True,
            "ample_square_12": True,
            "complete_478_wall_orbit_table_replayed": True,
            "all_720_matrices_unique_integral_ample_fixing_isometries": True,
        },
        "all_checks_pass": True,
    }


def validate_stock(payload: dict[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or not payload.get("all_checks_pass"):
        raise ArithmeticError("wrong or failed stock schema")
    chamber = payload["complete_induced_chamber"]
    group = payload["ample_stabilizer"]
    if chamber["wall_count"] != 478 or chamber["wall_orbit_sizes"] != EXPECTED_ORBIT_SIZES:
        raise ArithmeticError("archived wall denominator changed")
    if group["order"] != 720 or len(group["all_integral_20_by_20_matrices_in_supplement_order"]) != 720:
        raise ArithmeticError("archived group denominator changed")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--source-zip", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--tex-source", type=Path, default=DEFAULT_TEX)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_stock(args.source_dir, args.source_zip, args.tex_source)
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != encoded:
            raise ArithmeticError(f"archive mismatch: {OUTPUT}")
        print(f"archive matches: {OUTPUT}")
        return
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(encoded, encoding="utf-8")
    print(json.dumps({
        "output": str(OUTPUT),
        "walls": payload["complete_induced_chamber"]["wall_count"],
        "group_order": payload["ample_stabilizer"]["order"],
        "all_checks_pass": payload["all_checks_pass"],
    }, indent=2))


if __name__ == "__main__":
    main()
