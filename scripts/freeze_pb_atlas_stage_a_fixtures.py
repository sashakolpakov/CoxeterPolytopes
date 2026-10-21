#!/usr/bin/env python3
"""Freeze four arrays-only physical Coble calibration inputs for PB-Atlas A.

This is a one-time fixture extractor, not part of the recognizer.  It copies
only Picard Grams, primitive oriented wall vectors, marked boundary arrays,
and canonical classes.  Expected counts are independent replay targets.  No
stored endpoint Gram, contraction frame, incidence program, surface equation,
or output witness is copied into the numerical input.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.pb_atlas_contraction_core import INPUT_SCHEMA  # noqa: E402


OUTPUT = ROOT / "data/fixtures/pb_atlas_numerical_stage_a_inputs.json"
SCHEMA = "coxeter.pb_atlas.numerical_stage_a_fixtures.v1"

SOURCES = {
    "r11_m30_b1": ROOT / "data/enumerations/coble_one_boundary_full_ns_h10.json",
    "r12_m26_b2": ROOT / "data/enumerations/coble_type_i_two_boundary_full_ns_h11.json",
    "r15_m60_b5": ROOT / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json",
    "r18_m48_b8": ROOT / "data/enumerations/coble_type_mii_full_ns_h17.json",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _boundary_indices(gram: list[list[int]], walls: list[list[int]]) -> list[int]:
    lattice = sp.Matrix(gram)
    vectors = sp.Matrix(walls)
    wall_gram = vectors * lattice * vectors.T
    return [index for index in range(len(walls)) if wall_gram[index, index] == -4]


def build_fixtures() -> dict[str, Any]:
    one = _load(SOURCES["r11_m30_b1"])["cases"]["II"]
    two = _load(SOURCES["r12_m26_b2"])
    five = _load(SOURCES["r15_m60_b5"])
    eight = _load(SOURCES["r18_m48_b8"])

    h10_walls = one["wall_vectors"]
    h10_canonical = [-entry for entry in h10_walls[0]]
    h11_gram = two["full_NS_lattice"]["gram"]
    h11_walls = two["facets"]["NS_wall_vectors"]
    h11_boundary = _boundary_indices(h11_gram, h11_walls)
    h11_canonical = [
        -sum(h11_walls[index][column] for index in h11_boundary) // 2
        for column in range(len(h11_gram))
    ]
    h14_gram = five["full_NS_lattice"]["gram"]
    h14_walls = five["facets"]["NS_wall_vectors"]
    h17_gram = eight["full_NS"]["gram"]
    h17_walls = eight["facets"]["physical_NS_wall_vectors"]

    cases = [
        {
            "fixture_id": "r11_m30_b1",
            "calibration_alias": "H10/30",
            "input": {
                "schema": INPUT_SCHEMA,
                "picard_gram": one["lattice_gram"],
                "wall_vectors": h10_walls,
                "boundary_indices_zero_based": [0],
                "boundary_weights": [2],
                "canonical_class": h10_canonical,
            },
            "expected": {
                "wall_gram_automorphism_group_order": 24,
                "labelled_frame_counts_by_depth": [
                    1, 17, 87, 329, 948, 2094, 3472, 4417, 4455, 3291, 573
                ],
                "hirzebruch_histogram": {"F_0": 837, "F_1": 2064, "F_2": 390},
                "P2_count": 573,
            },
        },
        {
            "fixture_id": "r12_m26_b2",
            "calibration_alias": "H11/26",
            "input": {
                "schema": INPUT_SCHEMA,
                "picard_gram": h11_gram,
                "wall_vectors": h11_walls,
                "boundary_indices_zero_based": h11_boundary,
                "boundary_weights": [1] * len(h11_boundary),
                "canonical_class": h11_canonical,
            },
            "expected": {
                "wall_gram_automorphism_group_order": 8,
                "labelled_frame_counts_by_depth": [
                    1, 14, 39, 72, 125, 218, 335, 424, 486, 512, 388, 52
                ],
                "hirzebruch_histogram": {
                    "F_0": 86,
                    "F_1": 212,
                    "F_2": 54,
                    "F_3": 24,
                    "F_4": 12,
                },
                "P2_count": 52,
            },
        },
        {
            "fixture_id": "r15_m60_b5",
            "calibration_alias": "H14/60",
            "input": {
                "schema": INPUT_SCHEMA,
                "picard_gram": h14_gram,
                "wall_vectors": h14_walls,
                "boundary_indices_zero_based": _boundary_indices(h14_gram, h14_walls),
                "boundary_weights": [1] * 5,
                "canonical_class": five["full_NS_lattice"][
                    "canonical_vector_in_NS_basis"
                ],
            },
            "expected": {
                "wall_gram_automorphism_group_order": 120,
                "labelled_frame_counts_by_depth": [
                    1,
                    45,
                    355,
                    1620,
                    5200,
                    13032,
                    26880,
                    45100,
                    58640,
                    56270,
                    38136,
                    17210,
                    4615,
                    555,
                    5,
                ],
                "hirzebruch_histogram": {
                    "F_0": 105,
                    "F_1": 310,
                    "F_2": 120,
                    "F_3": 20,
                },
                "P2_count": 5,
            },
        },
        {
            "fixture_id": "r18_m48_b8",
            "calibration_alias": "H17/48",
            "input": {
                "schema": INPUT_SCHEMA,
                "picard_gram": h17_gram,
                "wall_vectors": h17_walls,
                "boundary_indices_zero_based": _boundary_indices(h17_gram, h17_walls),
                "boundary_weights": [1] * 8,
                "canonical_class": eight["full_NS"]["canonical"],
            },
            "expected": {
                "wall_gram_automorphism_group_order": 1152,
                "labelled_frame_counts_by_depth": [
                    1,
                    16,
                    216,
                    2288,
                    13628,
                    50256,
                    126664,
                    230896,
                    311238,
                    309232,
                    221416,
                    109584,
                    34940,
                    6320,
                    504,
                    16,
                    1,
                    0,
                ],
                "hirzebruch_histogram": {"F_0": 1},
                "P2_count": 0,
            },
        },
    ]
    return {
        "schema": SCHEMA,
        "scope": (
            "frozen arrays-only calibration inputs; aliases and source hashes are metadata "
            "outside the recognizer and never enter its branch parser"
        ),
        "forbidden_input_fields": [
            "surface/type name",
            "source wall labels",
            "stored endpoint Gram",
            "stored contraction frame",
            "known incidence coordinates or equations",
        ],
        "sources": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "sha256": _sha256(path),
            }
            for key, path in SOURCES.items()
        },
        "cases": cases,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    rendered = json.dumps(build_fixtures(), indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("frozen PB-Atlas Stage-A fixtures differ from exact extraction")
        print("fixture replay passed")
    elif arguments.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
