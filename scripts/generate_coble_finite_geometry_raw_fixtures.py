#!/usr/bin/env python3
"""Mechanically extract neutral raw-wall fixtures from five frozen rows.

This is deliberately separate from the recognizer.  It knows where the
calibration archives live only in order to serialize their physical Picard
Gram and primitive walls into the common input schema.  No source label or
archive digest is copied into a fixture, and the recognition code never
imports this module.
"""

from __future__ import annotations

import argparse
import copy
import itertools
import json
import math
from pathlib import Path
from typing import Any, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data/fixtures/coble_finite_geometry_raw"
RAW_SCHEMA = "coble-fg-raw-wall-input-v1"
MANIFEST_SCHEMA = "coble-fg-raw-wall-fixture-manifest-v1"


SOURCE_FIELDS = (
    (
        ROOT / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json",
        ("full_NS_lattice", "gram"),
        ("facets", "NS_wall_vectors"),
    ),
    (
        ROOT / "data/enumerations/coble_type_mii_full_ns_h17.json",
        ("full_NS", "gram"),
        ("facets", "physical_NS_wall_vectors"),
    ),
    (
        ROOT / "data/enumerations/coble_type_mi_full_ns_h11.json",
        ("surface", "NS_gram"),
        ("full_nef_chamber", "NS_wall_vectors"),
    ),
    (
        ROOT / "data/enumerations/coble_type_vi_char5_one_boundary_full_ns_h10.json",
        ("polytope", "lattice_gram"),
        ("polytope", "wall_vectors"),
    ),
    (
        ROOT / "data/enumerations/coble_type_vii_char5_one_boundary_full_ns_h10.json",
        ("polytope", "lattice_gram"),
        ("polytope", "wall_vectors"),
    ),
)


def _nested(payload: dict[str, Any], keys: Sequence[str]) -> Any:
    value: Any = payload
    for key in keys:
        value = value[key]
    return value


def _positive_fixture(path: Path, gram_keys: Sequence[str], wall_keys: Sequence[str]) -> dict[str, Any]:
    source = json.loads(path.read_text(encoding="utf-8"))
    gram = [[int(value) for value in row] for row in _nested(source, gram_keys)]
    walls = sorted(
        ([int(value) for value in row] for row in _nested(source, wall_keys)),
        key=tuple,
    )
    rank = len(gram)
    if any(len(row) != rank for row in gram) or any(len(row) != rank for row in walls):
        raise ArithmeticError("source matrix shape changed")
    if abs(int(sp.Matrix(gram).det())) != 1:
        raise ArithmeticError("physical fixture is not unimodular")
    if any(math.gcd(*(abs(value) for value in wall)) != 1 for wall in walls):
        raise ArithmeticError("source wall is not primitive")
    return {
        "schema": RAW_SCHEMA,
        "lattice_gram": gram,
        "primitive_wall_vectors": walls,
        "orientation": {"positive_half_space": "pairing(x, wall) >= 0"},
    }


def _reflect_one_root_only(payload: dict[str, Any]) -> dict[str, Any]:
    """Preserve every crude count/norm while breaking acute orientation."""

    result = copy.deepcopy(payload)
    gram = sp.Matrix(result["lattice_gram"])
    walls = result["primitive_wall_vectors"]
    roots = [
        index
        for index, wall in enumerate(walls)
        if int((sp.Matrix(wall).T * gram * sp.Matrix(wall))[0]) == -2
    ]
    left, right = next(
        (left, right)
        for left, right in itertools.combinations(roots, 2)
        if int((sp.Matrix(walls[left]).T * gram * sp.Matrix(walls[right]))[0]) == 1
    )
    replacement = [
        walls[right][column] + walls[left][column] for column in range(len(walls[right]))
    ]
    if math.gcd(*(abs(value) for value in replacement)) != 1:
        raise ArithmeticError("reflected mutation is not primitive")
    if int((sp.Matrix(replacement).T * gram * sp.Matrix(replacement))[0]) != -2:
        raise ArithmeticError("reflected mutation did not preserve root norm")
    if replacement in walls:
        raise ArithmeticError("reflected mutation collided with an existing wall")
    walls[right] = replacement
    walls.sort(key=tuple)
    return result


def _make_nonprimitive(payload: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(payload)
    result["primitive_wall_vectors"][0] = [
        2 * value for value in result["primitive_wall_vectors"][0]
    ]
    return result


def build_fixtures() -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    positives = [
        _positive_fixture(path, gram_keys, wall_keys)
        for path, gram_keys, wall_keys in SOURCE_FIELDS
    ]
    files: dict[str, dict[str, Any]] = {}
    cases = []
    for index, payload in enumerate(positives, start=1):
        filename = f"positive_{index:02d}.json"
        files[filename] = payload
        cases.append(
            {
                "case_id": f"positive_{index:02d}",
                "input": filename,
                "expected_accepted": True,
            }
        )
    for index, payload in enumerate(positives, start=1):
        filename = f"negative_norm_preserving_{index:02d}.json"
        files[filename] = _reflect_one_root_only(payload)
        cases.append(
            {
                "case_id": f"negative_norm_preserving_{index:02d}",
                "input": filename,
                "expected_accepted": False,
                "baseline": f"positive_{index:02d}.json",
                "mutation": "acute/orientation mutation: reflect one root only while preserving rank, facet count, and wall-norm histogram",
            }
        )
    files["negative_nonprimitive.json"] = _make_nonprimitive(positives[0])
    cases.append(
        {
            "case_id": "negative_nonprimitive",
            "input": "negative_nonprimitive.json",
            "expected_accepted": False,
        }
    )
    manifest = {"schema": MANIFEST_SCHEMA, "cases": cases}
    return files, manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    files, manifest = build_fixtures()
    rendered = {
        **{name: json.dumps(payload, indent=2, sort_keys=True) + "\n" for name, payload in files.items()},
        "manifest.json": json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    }
    if arguments.check:
        differing = [
            name
            for name, text in rendered.items()
            if not (arguments.output_dir / name).is_file()
            or (arguments.output_dir / name).read_text(encoding="utf-8") != text
        ]
        if differing:
            raise SystemExit(f"fixture files differ: {differing}")
    else:
        arguments.output_dir.mkdir(parents=True, exist_ok=True)
        for name, text in rendered.items():
            (arguments.output_dir / name).write_text(text, encoding="utf-8")
    print(json.dumps({"fixture_count": len(files), "manifest_cases": len(manifest["cases"])}, indent=2))


if __name__ == "__main__":
    main()
