#!/usr/bin/env python3
"""Ingest Shimada's theorem-complete 1,094-wall stock for ``X_1``.

The exact source is the computational supplement to arXiv:1412.6904.  The
paper proves that ``D^(0)`` is an induced chamber in the nef cone, that its
translates tessellate the nef cone, and that the displayed wall orbits and
the 720 matrices of ``Aut(X_1,a_1)`` are complete.  This is the genuine
surface-stabilizer denominator; it is not advertised as the full abstract
weighted-wall graph group.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.ingest_shimada_aut3singk3_x2_478_stock import (
    DEFAULT_SOURCE,
    DEFAULT_TEX,
    DEFAULT_ZIP,
    FILE_SHA256,
    TEX_SHA256,
    ZIP_SHA256,
    _assignment,
    _pair,
    _plain,
    _preserves_gram,
    _row_times_matrix,
    _sha256,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/shimada_aut3singk3_x1_1094_wall_stock.json"
SCHEMA = "coxeter.shimada.aut3singk3.x1.complete_1094_wall_stock.v1"
SURFACE_INDEX = 1

EXPECTED_ORBIT_SIZES = [45, 45, 10, 30, 72, 60, 12, 40, 120, 120, 120, 180, 120, 120]
EXPECTED_WALL_NORMS = [
    Fraction(-2),
    Fraction(-1, 2),
    Fraction(-3, 2),
    Fraction(-4, 3),
    Fraction(-5, 4),
    Fraction(-1),
    Fraction(-5, 6),
    Fraction(-3, 4),
    Fraction(-7, 12),
    Fraction(-7, 12),
    Fraction(-1, 3),
    Fraction(-1, 3),
    Fraction(-1, 12),
    Fraction(-1, 12),
]
EXPECTED_AMPLE_PAIRINGS = [2, 7, 3, 4, 5, 6, 5, 6, 7, 7, 8, 8, 8, 8]


def build_stock(source_dir: Path, source_zip: Path, tex_source: Path) -> dict[str, Any]:
    for filename, expected in FILE_SHA256.items():
        path = source_dir / filename
        if _sha256(path) != expected:
            raise ArithmeticError(f"pinned source changed: {path}")
    if _sha256(source_zip) != ZIP_SHA256 or _sha256(tex_source) != TEX_SHA256:
        raise ArithmeticError("pinned supplement zip or paper TeX source changed")

    gram = _assignment(source_dir / "GramS.txt", "GramS", SURFACE_INDEX)
    ample = _assignment(source_dir / "ample.txt", "ample", SURFACE_INDEX)
    wall_orbits = _assignment(source_dir / "walls.txt", "walls", SURFACE_INDEX)
    automorphisms = _assignment(source_dir / "AutXa.txt", "AutXa", SURFACE_INDEX)

    gram_matrix = sp.Matrix(gram)
    if not (
        gram_matrix.shape == (20, 20)
        and gram_matrix == gram_matrix.T
        and gram_matrix.det() == -24
        and exact_inertia(gram_matrix) == (1, 19, 0)
    ):
        raise ArithmeticError("unexpected X_1 Neron--Severi lattice")
    if len(ample) != 20 or _pair(ample, gram, ample) != 30:
        raise ArithmeticError("unexpected X_1 ample class")
    if [len(orbit) for orbit in wall_orbits] != EXPECTED_ORBIT_SIZES:
        raise ArithmeticError("the complete X_1 wall-orbit census changed")
    for orbit_index, orbit in enumerate(wall_orbits):
        for wall in orbit:
            if len(wall) != 20:
                raise ArithmeticError("a wall has the wrong rank")
            if _pair(wall, gram, wall) != EXPECTED_WALL_NORMS[orbit_index]:
                raise ArithmeticError(f"wall norm changed in orbit {orbit_index + 1}")
            if _pair(ample, gram, wall) != EXPECTED_AMPLE_PAIRINGS[orbit_index]:
                raise ArithmeticError(f"ample pairing changed in orbit {orbit_index + 1}")

    identity = [[Fraction(int(i == j)) for j in range(20)] for i in range(20)]
    flattened = {tuple(value for row in matrix for value in row) for matrix in automorphisms}
    if len(automorphisms) != 720 or len(flattened) != 720:
        raise ArithmeticError("Aut(X_1,a_1) is not the expected 720-element list")
    identity_indices = [index for index, matrix in enumerate(automorphisms) if matrix == identity]
    if identity_indices != [7]:
        raise ArithmeticError(f"identity index changed: {identity_indices}")
    gram_int = [[int(value) for value in row] for row in gram]
    for index, matrix in enumerate(automorphisms):
        if any(value.denominator != 1 for row in matrix for value in row):
            raise ArithmeticError(f"automorphism {index} is not integral")
        if _row_times_matrix(ample, matrix) != ample:
            raise ArithmeticError(f"automorphism {index} does not fix a_1")
        if not _preserves_gram(
            [[int(value) for value in row] for row in matrix], gram_int
        ):
            raise ArithmeticError(f"automorphism {index} does not preserve GramS[1]")

    orbit_records = [
        {
            "orbit_one_based": index + 1,
            "size": len(orbit),
            "wall_norm": _plain(EXPECTED_WALL_NORMS[index]),
            "ample_pairing": EXPECTED_AMPLE_PAIRINGS[index],
            "primitive_outward_wall_vectors_in_published_NS_basis": _plain(orbit),
        }
        for index, orbit in enumerate(wall_orbits)
    ]
    return {
        "schema": SCHEMA,
        "id": "SHIMADA-AUT3SINGK3-X1-D0-1094",
        "sources": {
            "paper": "Ichiro Shimada, The automorphism groups of certain singular K3 surfaces and an Enriques surface, arXiv:1412.6904",
            "paper_url": "https://arxiv.org/abs/1412.6904",
            "supplement_page": "https://www.math.sci.hiroshima-u.ac.jp/shimada/K3andEnriques.html",
            "supplement_zip_url": "https://www.math.sci.hiroshima-u.ac.jp/shimada/preprints/Aut3SingK3s/compdata/compdata.zip",
            "supplement_zip_sha256": ZIP_SHA256,
            "paper_tex_sha256": TEX_SHA256,
            "supplement_file_sha256": FILE_SHA256,
            "selected_surface_zero_based": SURFACE_INDEX,
            "theorem_scope": {
                "transcendental_lattice": "paper display defining X_1: Gram [[2,0],[0,12]]",
                "chamber_completeness": "Theorem thm:D0: D^(0) is induced in the nef cone and its Aut(X)-translates tessellate the nef cone",
                "wall_completeness": "Theorem thm:orbits and Table table:main: complete Delta(D^(0)) in fourteen Aut(X_1,a_1)-orbits",
                "stabilizer": "Theorem thm:aut: Aut(X_1,a_1) is PGL_2(F_9), order 720",
                "denominator_caveat": "No claim that this genuine surface group is the full abstract weighted-wall graph group",
            },
        },
        "surface": {
            "name": "X_1",
            "kind": "complex singular K3 surface",
            "transcendental_gram": [[2, 0], [0, 12]],
            "neron_severi_rank": 20,
            "neron_severi_signature": [1, 19],
            "neron_severi_determinant": -24,
            "published_neron_severi_gram": _plain(gram),
            "published_ample_vector_a1": _plain(ample),
            "ample_square": 30,
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
            "identity_zero_based_index": 7,
            "all_integral_20_by_20_matrices_in_supplement_order": _plain(automorphisms),
        },
        "checks": {
            "pinned_sources_match": True,
            "lorentzian_NS_rank20_determinant_minus24": True,
            "ample_square_30": True,
            "complete_1094_wall_orbit_table_replayed": True,
            "all_720_matrices_unique_integral_ample_fixing_isometries": True,
        },
        "all_checks_pass": True,
    }


def validate_stock(payload: dict[str, Any]) -> None:
    if payload.get("schema") != SCHEMA or not payload.get("all_checks_pass"):
        raise ArithmeticError("wrong or failed stock schema")
    if payload["complete_induced_chamber"]["wall_count"] != 1094:
        raise ArithmeticError("archived wall denominator changed")
    if payload["ample_stabilizer"]["order"] != 720:
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
    print(json.dumps({"output": str(OUTPUT), "walls": 1094, "group_order": 720}, indent=2))


if __name__ == "__main__":
    main()
