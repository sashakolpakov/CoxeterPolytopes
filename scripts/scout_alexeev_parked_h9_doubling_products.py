#!/usr/bin/env python3
"""Scout product redoublings of the five parked rank-ten H9 chambers.

This is deliberately an exploratory exact-arithmetic script.  It constructs
finite paths, path-times-cubes, and rectangles in the Cayley graph generated
by the certified doubling walls, then checks the resulting boundary walls,
strict facet witnesses, Coxeter products, and integral reflections.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
from typing import Any, Iterable, Sequence

import sympy as sp

from scripts.certify_alexeev_h9_doubling_cubes import (
    _case_payload,
    _cone,
    _doubling_walls,
    _matrix,
    _pair,
    _primitive,
    _reflection,
    _relation_histogram,
    _relations,
    _unoriented_key,
)
from scripts.certify_alexeev_20_2_1_15_9_1_exact_cones import (
    _analyze_case,
    _parse_cdd_rays,
    _run_cdd,
)
from scripts.explore_alexeev_14_10_1_rank10_c1 import (
    _literal_vinberg,
    _time_orientation,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/alexeev_parked_rank10_h9_examples.json"
SOURCE_SHA256 = "e181b8c021fe0f37d23b53b33e4e3900559502d0f8c982a4a9e9cbcf0866e430"

CONFIGS: dict[str, dict[str, Any]] = {
    "A14_8_1_C02": {"paths": [(14, 15)], "cube": [6]},
    "A14_10_1_C01": {"paths": [], "cube": [2, 4, 10]},
    "A14_10_1_C02": {"paths": [(5, 14)], "cube": [0, 2]},
    "A16_14_1_C01": {"paths": [(0, 13), (7, 14)], "cube": []},
}


def _key(matrix: sp.Matrix) -> tuple[sp.Rational, ...]:
    return tuple(sp.Rational(value) for value in matrix)


def _load_seed(identifier: str) -> tuple[
    dict[str, Any],
    sp.Matrix,
    tuple[sp.Matrix, ...],
    sp.Matrix,
    tuple[sp.Matrix, ...],
]:
    if sha256(SOURCE.read_bytes()).hexdigest() != SOURCE_SHA256:
        raise ArithmeticError("the parked H9 source changed")
    payload = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = next(row for row in payload["cases"] if row["id"] == identifier)
    section = case["section"]
    gram = _matrix(section["saturated_fixed_lattice"]["gram"])
    active = tuple(
        row
        for row in section["restricted_mirrors"]["records"]
        if row["hyperbolic_status"] == "active_facet"
    )
    roots = tuple(
        sp.Matrix(row["primitive_integral_same_direction_coordinates"])
        for row in active
    )
    interior = sp.Matrix(
        [sp.Rational(value) for value in section["strict_fixed_interior"]["fixed_lattice_coordinates"]]
    )
    roots = tuple(-root if _pair(gram, interior, root) < 0 else root for root in roots)
    witness_by_orbit = {
        int(row["orbit_index"]): sp.Matrix(
            [sp.Rational(value) for value in row["fixed_lattice_coordinates"]]
        )
        for row in section["exact_reduction"]["strict_facet_witnesses"]
    }
    witnesses = tuple(witness_by_orbit[int(row["orbit_index"])] for row in active)
    expected = sorted(
        [index for pair in CONFIGS[identifier]["paths"] for index in pair]
        + CONFIGS[identifier]["cube"]
    )
    relations = _relations(gram, roots)
    if _doubling_walls(relations) != expected:
        raise ArithmeticError(f"{identifier}: doubling-wall stock changed")
    return section, gram, roots, interior, witnesses


def _path_transforms(
    gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    pair: tuple[int, int],
    copies: int,
    start: int = 0,
) -> tuple[sp.Matrix, ...]:
    if copies < 1 or start not in (0, 1):
        raise ValueError("copies must be positive and start must be zero or one")
    generators = tuple(_reflection(gram, roots[index]) for index in pair)
    transforms = [sp.eye(gram.rows)]
    for edge in range(copies - 1):
        transforms.append(transforms[-1] * generators[(start + edge) % 2])
    return tuple(transforms)


def _cube_transforms(
    gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    indices: Sequence[int],
) -> tuple[sp.Matrix, ...]:
    generators = tuple(_reflection(gram, roots[index]) for index in indices)
    if any(
        left * right != right * left
        for position, left in enumerate(generators)
        for right in generators[position + 1 :]
    ):
        raise ArithmeticError("the proposed cube generators do not commute")
    transforms = []
    for mask in range(1 << len(generators)):
        transform = sp.eye(gram.rows)
        for index, generator in enumerate(generators):
            if mask & (1 << index):
                transform = transform * generator
        transforms.append(transform)
    return tuple(transforms)


def _product(*factors: Sequence[sp.Matrix]) -> tuple[sp.Matrix, ...]:
    if not factors:
        raise ValueError("at least one transform factor is required")
    transforms = [sp.eye(factors[0][0].rows)]
    for factor in factors:
        transforms = [left * right for left in transforms for right in factor]
    unique = {_key(transform): transform for transform in transforms}
    if len(unique) != len(transforms):
        raise ArithmeticError("a proposed product transform set has collisions")
    return tuple(unique[key] for key in sorted(unique))


def _boundary(
    gram: sp.Matrix,
    roots: Sequence[sp.Matrix],
    interior: sp.Matrix,
    witnesses: Sequence[sp.Matrix],
    transforms: Sequence[sp.Matrix],
    generator_indices: Iterable[int],
) -> dict[str, Any]:
    transform_keys = {_key(transform) for transform in transforms}
    generators = {
        index: _reflection(gram, roots[index]) for index in generator_indices
    }
    internal: set[tuple[int, ...]] = set()
    for transform in transforms:
        for index, generator in generators.items():
            if _key(transform * generator) in transform_keys:
                internal.add(_unoriented_key(transform * roots[index]))

    raw_interior = sum(
        (transform * interior for transform in transforms),
        sp.zeros(gram.rows, 1),
    )
    boundary: dict[tuple[int, ...], tuple[sp.Matrix, sp.Matrix]] = {}
    for transform in transforms:
        for index, seed_root in enumerate(roots):
            root = _primitive(transform * seed_root)
            key = _unoriented_key(root)
            if key in internal:
                continue
            pairing = _pair(gram, raw_interior, root)
            if pairing == 0:
                raise ArithmeticError("the union interior lies on a boundary wall")
            if pairing < 0:
                root = -root
            boundary.setdefault(key, (root, transform * witnesses[index]))

    ordered = tuple(boundary[key] for key in sorted(boundary))
    boundary_roots = tuple(row[0] for row in ordered)
    boundary_witnesses = tuple(row[1] for row in ordered)
    if not (
        _pair(gram, raw_interior, raw_interior) > 0
        and all(_pair(gram, raw_interior, root) > 0 for root in boundary_roots)
        and sp.Matrix.hstack(*boundary_roots).rank() == gram.rows
    ):
        raise ArithmeticError("the union boundary is not strict and full dimensional")
    if not all(
        _pair(gram, boundary_witnesses[target], boundary_roots[target]) == 0
        and _pair(gram, boundary_witnesses[target], boundary_witnesses[target]) > 0
        and all(
            _pair(gram, boundary_witnesses[target], boundary_roots[index]) > 0
            for index in range(len(boundary_roots))
            if index != target
        )
        for target in range(len(boundary_roots))
    ):
        raise ArithmeticError("a boundary wall lacks a strict timelike witness")
    relations = _relations(gram, boundary_roots)
    for root in boundary_roots:
        _reflection(gram, root)
    return {
        "copy_count": len(transforms),
        "facet_count": len(boundary_roots),
        "relation_histogram": _relation_histogram(relations),
        "doubling_walls_one_based": [index + 1 for index in _doubling_walls(relations)],
        "roots": boundary_roots,
        "interior": raw_interior,
        "witnesses": boundary_witnesses,
    }


def _volume_audit(
    identifier: str,
    section: dict[str, Any],
    gram: sp.Matrix,
    record: dict[str, Any],
) -> dict[str, Any]:
    roots = record["roots"]
    interior = record["interior"]
    cone = _cone(gram, roots, interior)
    cdd = shutil.which("cddexec_gmp") or "cddexec_gmp"
    rays = _parse_cdd_rays(_run_cdd(cone["inequality_text"], cdd), gram.rows)
    case = _case_payload(identifier, section, gram, roots, interior)
    volume, _archive = _analyze_case(case, cone, rays)
    time = _time_orientation(case, cone, rays)
    vinberg = _literal_vinberg(case, cone, rays)
    return {
        "classification": volume["decision"]["classification"],
        "finite": volume["exact_double_description"]["finite_vertex_count"],
        "ideal": volume["exact_double_description"]["ideal_vertex_count"],
        "hyperideal": volume["exact_double_description"]["hyperideal_vertex_count"],
        "same_time_cone": time["all_extreme_rays_same_time_cone"],
        "elliptic_edges": vinberg["elliptic_edge_count"],
        "affine_components": vinberg["connected_affine_component_count"],
        "bad_edges": vinberg["bad_edge_count"],
        "bad_parabolic_neighborhoods": vinberg["bad_parabolic_neighborhood_count"],
    }


def scout(limit: int, audit_small: bool = False) -> dict[str, Any]:
    output: dict[str, Any] = {}
    for identifier, config in CONFIGS.items():
        section, gram, roots, interior, witnesses = _load_seed(identifier)
        all_generators = [
            index for pair in config["paths"] for index in pair
        ] + config["cube"]
        rows = []
        cube = _cube_transforms(gram, roots, config["cube"])
        if not config["paths"]:
            boundary = _boundary(
                gram,
                roots,
                interior,
                witnesses,
                cube,
                all_generators,
            )
            row = {
                    "shape": f"cube_{len(config['cube'])}",
                    **{
                        key: value
                        for key, value in boundary.items()
                        if key not in {"roots", "interior", "witnesses"}
                    },
                }
            if audit_small:
                row["volume_audit"] = _volume_audit(
                    f"{identifier}_{row['shape']}", section, gram, boundary
                )
            rows.append(row)
        elif len(config["paths"]) == 1:
            pair = config["paths"][0]
            for copies in range(1, limit + 1):
                for start in range(2 if copies % 2 == 0 else 1):
                    path = _path_transforms(gram, roots, pair, copies, start)
                    transforms = _product(path, cube)
                    record = _boundary(
                        gram,
                        roots,
                        interior,
                        witnesses,
                        transforms,
                        all_generators,
                    )
                    row = {
                            "shape": f"path_{copies}_phase_{start}_cube_{len(config['cube'])}",
                            **{
                                key: value
                                for key, value in record.items()
                                if key not in {"roots", "interior", "witnesses"}
                            },
                        }
                    if audit_small and copies == 2 and start == 0:
                        row["volume_audit"] = _volume_audit(
                            f"{identifier}_{row['shape']}", section, gram, record
                        )
                    rows.append(row)
        else:
            first, second = config["paths"]
            for width in range(1, limit + 1):
                for height in range(1, limit + 1):
                    path_x = _path_transforms(gram, roots, first, width, 0)
                    path_y = _path_transforms(gram, roots, second, height, 0)
                    transforms = _product(path_x, path_y)
                    record = _boundary(
                        gram,
                        roots,
                        interior,
                        witnesses,
                        transforms,
                        all_generators,
                    )
                    row = {
                            "shape": f"rectangle_{width}x{height}",
                            **{
                                key: value
                                for key, value in record.items()
                                if key not in {"roots", "interior", "witnesses"}
                            },
                        }
                    if audit_small and width == 2 and height == 2:
                        row["volume_audit"] = _volume_audit(
                            f"{identifier}_{row['shape']}", section, gram, record
                        )
                    rows.append(row)
        output[identifier] = rows
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=4)
    parser.add_argument("--audit-small", action="store_true")
    args = parser.parse_args()
    print(json.dumps(scout(args.limit, args.audit_small), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
