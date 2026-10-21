#!/usr/bin/env python3
"""Independent hostile replay of the two Alexeev H9 doubling cubes.

This audit intentionally does not import the primary doubling-cube certifier.
Starting from the pinned H9/14 factory, it independently reconstructs all 30
nonempty reflection-subgroup gluings, checks the exact convex-boundary lemma,
classifies the subset diagrams, and compares the nine representatives with the
pinned primary JSON certificate.
"""

from __future__ import annotations

from collections import Counter
from functools import reduce
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
FACTORY = ROOT / "data" / "enumerations" / "alexeev_k3_h9_14_factory.json"
PRIMARY = ROOT / "data" / "enumerations" / "alexeev_h9_doubling_cubes.json"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_h9_doubling_cubes_independent_audit.json"
)
FACTORY_SHA256 = "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775"
PRIMARY_SHA256 = "4919728a5c5ff8adebd80b250cfedfb562a040c8bd1e72e551b398b1f13e7e4c"
SCHEMA = "coxeter.alexeev_h9_doubling_cubes_independent_audit.v1"

SEEDS = {
    "T1": {
        "source_id": "A14_8_1_C01",
        "doubling_walls": (4, 7, 12, 14),
        "diagram_automorphism_order": 8,
        "type_counts": {1: 1, 2: 2, 3: 1, 4: 1},
    },
    "T2": {
        "source_id": "A18_4_0_C01",
        "doubling_walls": (1, 5, 9, 14),
        "diagram_automorphism_order": 24,
        "type_counts": {1: 1, 2: 1, 3: 1, 4: 1},
    },
}


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _rational_matrix(rows: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in rows])


def _pair(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> sp.Rational:
    return sp.Rational((left.T * gram * right)[0])


def _content(values: Iterable[int]) -> int:
    nonzero = [abs(value) for value in values if value]
    return reduce(math.gcd, nonzero)


def _primitive(vector: sp.Matrix) -> sp.Matrix:
    entries = [sp.Rational(value) for value in vector]
    denominator = int(sp.ilcm(*(entry.q for entry in entries)))
    integral = [int(entry * denominator) for entry in entries]
    divisor = _content(integral)
    return sp.Matrix([value // divisor for value in integral])


def _positive_key(vector: sp.Matrix) -> tuple[int, ...]:
    primitive = tuple(int(value) for value in _primitive(vector))
    first = next(value for value in primitive if value)
    return primitive if first > 0 else tuple(-value for value in primitive)


def _oriented(
    gram: sp.Matrix, interior: sp.Matrix, vector: sp.Matrix
) -> sp.Matrix:
    primitive = _primitive(vector)
    value = _pair(gram, interior, primitive)
    if value == 0:
        raise ArithmeticError("a purported boundary contains the averaged interior")
    return primitive if value > 0 else -primitive


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = _pair(gram, root, root)
    matrix = sp.eye(gram.rows) - sp.Rational(2, square) * root * (root.T * gram)
    if not (
        all(value.is_Integer for value in matrix)
        and matrix.T * gram * matrix == gram
        and matrix * matrix == sp.eye(gram.rows)
        and matrix.det() == -1
    ):
        raise ArithmeticError("nonintegral or invalid doubling reflection")
    return matrix


def _token_from_ratio(ratio: sp.Rational) -> str:
    known = {
        sp.Rational(0): "2",
        sp.Rational(1, 4): "3",
        sp.Rational(1, 2): "4",
        sp.Rational(3, 4): "6",
        sp.Rational(1): "parallel",
    }
    if ratio in known:
        return known[ratio]
    if ratio > 1:
        return f"dotted:{ratio}"
    raise ArithmeticError(f"non-Coxeter normalized product {ratio}")


def _relations(
    gram: sp.Matrix, roots: Sequence[sp.Matrix]
) -> tuple[tuple[str, ...], ...]:
    rows: list[tuple[str, ...]] = []
    for left, alpha in enumerate(roots):
        row = []
        for right, beta in enumerate(roots):
            if left == right:
                row.append("1")
                continue
            product = _pair(gram, alpha, beta)
            if product < 0:
                raise ArithmeticError("the reconstructed boundary is not acute")
            ratio = sp.factor(
                product**2
                / (_pair(gram, alpha, alpha) * _pair(gram, beta, beta))
            )
            row.append(_token_from_ratio(ratio))
        rows.append(tuple(row))
    return tuple(rows)


def _histogram(relations: Sequence[Sequence[str]]) -> dict[str, int]:
    return dict(
        sorted(
            Counter(
                relations[left][right]
                for left in range(len(relations))
                for right in range(left)
            ).items()
        )
    )


def _doubling_walls(relations: Sequence[Sequence[str]]) -> tuple[int, ...]:
    odd = {"3", "5", "7", "9"}
    return tuple(
        index + 1
        for index, row in enumerate(relations)
        if not any(token in odd for token in row)
    )


def _graph(relations: Sequence[Sequence[str]]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(relations)))
    for left in range(len(relations)):
        for right in range(left):
            if relations[left][right] != "2":
                graph.add_edge(left, right, relation=relations[left][right])
    return graph


EDGE_MATCH = nx.algorithms.isomorphism.categorical_edge_match("relation", None)


def _isomorphic(
    left: Sequence[Sequence[str]], right: Sequence[Sequence[str]]
) -> bool:
    return nx.is_isomorphic(_graph(left), _graph(right), edge_match=EDGE_MATCH)


def _json_hash(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return sha256(encoded).hexdigest()


def _load_seed(
    factory: dict[str, Any], seed_name: str
) -> tuple[dict[str, Any], sp.Matrix, tuple[sp.Matrix, ...], sp.Matrix, tuple[sp.Matrix, ...]]:
    spec = SEEDS[seed_name]
    construction = next(
        row for row in factory["constructions"] if row["id"] == spec["source_id"]
    )
    section = construction["section"]
    gram = _rational_matrix(section["saturated_fixed_lattice"]["gram"])
    records = tuple(
        row
        for row in section["restricted_mirrors"]["records"]
        if row["hyperbolic_status"] == "active_facet"
    )
    roots = tuple(
        sp.Matrix(row["primitive_integral_same_direction_coordinates"])
        for row in records
    )
    interior = sp.Matrix(
        [
            sp.Rational(value)
            for value in section["strict_fixed_interior"]["fixed_lattice_coordinates"]
        ]
    )
    roots = tuple(_oriented(gram, interior, root) for root in roots)
    witness_map = {
        int(row["orbit_index"]): sp.Matrix(
            [sp.Rational(value) for value in row["fixed_lattice_coordinates"]]
        )
        for row in section["exact_reduction"]["strict_facet_witnesses"]
    }
    witnesses = tuple(witness_map[int(row["orbit_index"])] for row in records)
    if not (
        gram.shape == (10, 10)
        and gram == gram.T
        and gram.det() != 0
        and len(roots) == 14
        and all(_pair(gram, interior, root) > 0 for root in roots)
        and all(
            _pair(gram, witnesses[index], roots[index]) == 0
            and _pair(gram, witnesses[index], witnesses[index]) > 0
            and all(
                _pair(gram, witnesses[index], other) > 0
                for other_index, other in enumerate(roots)
                if other_index != index
            )
            for index in range(14)
        )
    ):
        raise ArithmeticError(f"seed {seed_name} failed independent exact replay")
    return section, gram, roots, interior, witnesses


def _subgroup(generators: Sequence[sp.Matrix]) -> tuple[sp.Matrix, ...]:
    transforms = []
    rank = generators[0].rows if generators else 10
    for mask in range(1 << len(generators)):
        transform = sp.eye(rank)
        for index, generator in enumerate(generators):
            if mask & (1 << index):
                transform = generator * transform
        transforms.append(transform)
    if len({tuple(matrix) for matrix in transforms}) != 1 << len(generators):
        raise ArithmeticError("reflection subgroup order is smaller than 2^k")
    return tuple(transforms)


def _reconstruct_subset(
    gram: sp.Matrix,
    seed_roots: Sequence[sp.Matrix],
    seed_interior: sp.Matrix,
    seed_witnesses: Sequence[sp.Matrix],
    selected_zero_based: Sequence[int],
) -> dict[str, Any]:
    selected = tuple(selected_zero_based)
    generators = tuple(_reflection(gram, seed_roots[index]) for index in selected)
    if any(
        left * right != right * left
        for position, left in enumerate(generators)
        for right in generators[position + 1 :]
    ):
        raise ArithmeticError("selected orthogonal reflections do not commute")
    transforms = _subgroup(generators)
    averaged_interior = sum(
        (transform * seed_interior for transform in transforms),
        sp.zeros(gram.rows, 1),
    )
    if not (
        _pair(gram, averaged_interior, averaged_interior) > 0
        and all(generator * averaged_interior == averaged_interior for generator in generators)
    ):
        raise ArithmeticError("averaged interior is not timelike and subgroup-fixed")

    # Exact orthogonal-doubling convexity lemma: every transformed unselected
    # wall is alpha plus a nonnegative combination of the selected normals.
    coefficient_values: list[sp.Rational] = []
    for index, root in enumerate(seed_roots):
        if index in selected:
            continue
        for transform in transforms:
            image = transform * root
            coefficients = tuple(
                sp.factor(_pair(gram, image - root, seed_roots[wall]) / _pair(gram, seed_roots[wall], seed_roots[wall]))
                for wall in selected
            )
            reconstructed = root + sum(
                (coefficient * seed_roots[wall] for coefficient, wall in zip(coefficients, selected)),
                sp.zeros(gram.rows, 1),
            )
            if reconstructed != image or any(value < 0 for value in coefficients):
                raise ArithmeticError("orthogonal doubling convexity identity failed")
            coefficient_values.extend(coefficients)

    # Collect full subgroup orbits of the unselected seed walls.  Track all
    # occurrences to exclude an accidental opposite-oriented identification.
    occurrences: dict[tuple[int, ...], list[tuple[sp.Matrix, sp.Matrix]]] = {}
    for transform in transforms:
        for index, root in enumerate(seed_roots):
            if index in selected:
                continue
            image = transform * root
            occurrences.setdefault(_positive_key(image), []).append(
                (image, transform * seed_witnesses[index])
            )
    boundary_roots: list[sp.Matrix] = []
    boundary_witnesses: list[sp.Matrix] = []
    for key in sorted(occurrences):
        images = occurrences[key]
        oriented_images = [
            _oriented(gram, averaged_interior, image) for image, _ in images
        ]
        if len({tuple(image) for image in oriented_images}) != 1:
            raise ArithmeticError("opposite-oriented duplicate boundary hyperplanes")
        boundary_roots.append(oriented_images[0])
        boundary_witnesses.append(images[0][1])
    roots = tuple(boundary_roots)
    witnesses = tuple(boundary_witnesses)
    root_keys = {tuple(root) for root in roots}
    if not (
        sp.Matrix.hstack(*roots).rank() == gram.rows
        and all(_pair(gram, averaged_interior, root) > 0 for root in roots)
        and all(
            tuple(_oriented(gram, averaged_interior, generator * root)) in root_keys
            for generator in generators
            for root in roots
        )
        and all(
            _pair(gram, witnesses[index], roots[index]) == 0
            and _pair(gram, witnesses[index], witnesses[index]) > 0
            and all(
                _pair(gram, witnesses[index], other) > 0
                for other_index, other in enumerate(roots)
                if other_index != index
            )
            for index in range(len(roots))
        )
    ):
        raise ArithmeticError("boundary invariance, rank, or strict-facet replay failed")

    # Selected normals give exactly the 2^k sign sectors of the glued copies.
    sign_vectors = {
        tuple(
            1 if _pair(gram, transform * seed_interior, seed_roots[wall]) > 0 else -1
            for wall in selected
        )
        for transform in transforms
    }
    if len(sign_vectors) != len(transforms):
        raise ArithmeticError("the glued copies do not occupy all selected-wall sign sectors")

    relations = _relations(gram, roots)
    return {
        "selected_walls_one_based": [index + 1 for index in selected],
        "copy_count": len(transforms),
        "reflection_subgroup_order": len(transforms),
        "sign_sector_count": len(sign_vectors),
        "boundary_orbit_hyperplane_count": len(roots),
        "facet_count": len(roots),
        "strict_facet_witness_count": len(witnesses),
        "boundary_subgroup_invariant": True,
        "no_opposite_oriented_boundary_duplicates": True,
        "convexity_nonnegative_coefficient_minimum": str(min(coefficient_values, default=0)),
        "convexity_nonnegative_coefficient_maximum": str(max(coefficient_values, default=0)),
        "relation_histogram": _histogram(relations),
        "doubling_walls_one_based": list(_doubling_walls(relations)),
        "boundary_root_sha256": _json_hash([[str(value) for value in root] for root in roots]),
        "relation_matrix_sha256": _json_hash(relations),
        "_roots": roots,
        "_relations": relations,
    }


def _seed_automorphisms(
    relations: Sequence[Sequence[str]], doubling_zero_based: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _graph(relations), _graph(relations), edge_match=EDGE_MATCH
    )
    wall_set = set(doubling_zero_based)
    permutations = set()
    for mapping in matcher.isomorphisms_iter():
        if {mapping[wall] for wall in wall_set} != wall_set:
            raise ArithmeticError("a seed automorphism does not preserve doubling walls")
        permutations.add(tuple(mapping[wall] for wall in doubling_zero_based))
    return tuple(sorted(permutations))


def _subset_orbits(
    walls: Sequence[int], permutations: Sequence[Sequence[int]], size: int
) -> list[list[list[int]]]:
    remaining = {tuple(subset) for subset in combinations(walls, size)}
    orbits = []
    while remaining:
        representative = min(remaining)
        orbit = {
            tuple(sorted(permutation[walls.index(wall)] for wall in representative))
            for permutation in permutations
        }
        remaining -= orbit
        orbits.append([[wall + 1 for wall in subset] for subset in sorted(orbit)])
    return sorted(orbits, key=lambda orbit: orbit[0])


def _diagram_classes(rows: Sequence[dict[str, Any]]) -> list[list[list[int]]]:
    classes: list[tuple[tuple[tuple[str, ...], ...], list[list[int]]]] = []
    for row in rows:
        match = next(
            (record for record in classes if _isomorphic(row["_relations"], record[0])),
            None,
        )
        if match is None:
            classes.append((row["_relations"], [row["selected_walls_one_based"]]))
        else:
            match[1].append(row["selected_walls_one_based"])
    return sorted((sorted(subsets) for _, subsets in classes), key=lambda orbit: orbit[0])


def audit() -> dict[str, Any]:
    if _digest(FACTORY) != FACTORY_SHA256 or _digest(PRIMARY) != PRIMARY_SHA256:
        raise ArithmeticError("a pinned source or primary certificate changed")
    factory = json.loads(FACTORY.read_text(encoding="utf-8"))
    primary = json.loads(PRIMARY.read_text(encoding="utf-8"))
    seed_payloads: dict[str, Any] = {}
    all_rows: list[dict[str, Any]] = []

    for seed_name, spec in SEEDS.items():
        _, gram, roots, interior, witnesses = _load_seed(factory, seed_name)
        seed_relations = _relations(gram, roots)
        expected_walls = tuple(spec["doubling_walls"])
        if _doubling_walls(seed_relations) != expected_walls:
            raise ArithmeticError(f"{seed_name} doubling-wall stock changed")
        walls = tuple(wall - 1 for wall in expected_walls)
        if any(
            _pair(gram, roots[left], roots[right]) != 0
            for position, left in enumerate(walls)
            for right in walls[position + 1 :]
        ):
            raise ArithmeticError(f"{seed_name} doubling walls are not orthogonal")
        reflections = tuple(_reflection(gram, roots[wall]) for wall in walls)
        permutations = _seed_automorphisms(seed_relations, walls)
        if len(permutations) != spec["diagram_automorphism_order"]:
            raise ArithmeticError(f"{seed_name} seed diagram automorphism order changed")

        rows: list[dict[str, Any]] = []
        orbit_summary: dict[str, Any] = {}
        for size in range(1, 5):
            size_rows = [
                _reconstruct_subset(gram, roots, interior, witnesses, subset)
                for subset in combinations(walls, size)
            ]
            rows.extend(size_rows)
            graph_classes = _diagram_classes(size_rows)
            automorphism_orbits = _subset_orbits(walls, permutations, size)
            if graph_classes != automorphism_orbits:
                raise ArithmeticError(
                    f"{seed_name} size-{size} diagram classes differ from seed orbits"
                )
            expected_primary = sorted(
                (
                    sorted(record["selected_subsets_one_based"])
                    for record in primary["seeds"][seed_name][
                        "isomorphism_census_by_copy_count"
                    ][str(1 << size)]["types"]
                ),
                key=lambda orbit: orbit[0],
            )
            if graph_classes != expected_primary or len(graph_classes) != spec["type_counts"][size]:
                raise ArithmeticError(f"{seed_name} size-{size} primary census mismatch")
            orbit_summary[str(size)] = {
                "subset_count": len(size_rows),
                "diagram_isomorphism_type_count": len(graph_classes),
                "seed_automorphism_orbit_count": len(automorphism_orbits),
                "subset_orbits_one_based": graph_classes,
            }

        # Check each archived representative against an independently rebuilt
        # subset, including the exact oriented boundary-root list.
        by_subset = {
            tuple(row["selected_walls_one_based"]): row for row in rows
        }
        for archived in primary["seeds"][seed_name]["certified_representatives"]:
            row = by_subset[tuple(archived["selected_doubling_walls_one_based"])]
            archived_roots = tuple(
                tuple(sp.Rational(value) for value in root)
                for root in archived["roots_in_fixed_lattice_coordinates"]
            )
            rebuilt_roots = tuple(tuple(root) for root in row["_roots"])
            if not (
                archived_roots == rebuilt_roots
                and archived["facet_count"] == row["facet_count"]
                and archived["relation_histogram"] == row["relation_histogram"]
                and archived["doubling_walls_one_based"] == row["doubling_walls_one_based"]
            ):
                raise ArithmeticError(f"{archived['id']} differs from independent replay")

        public_rows = [
            {key: value for key, value in row.items() if not key.startswith("_")}
            for row in rows
        ]
        seed_payloads[seed_name] = {
            "source_id": spec["source_id"],
            "seed_facet_count": len(roots),
            "doubling_walls_one_based": list(expected_walls),
            "doubling_wall_squares": [str(_pair(gram, roots[wall], roots[wall])) for wall in walls],
            "doubling_walls_mutually_orthogonal": True,
            "integral_reflection_count": len(reflections),
            "seed_diagram_automorphism_order_on_doubling_walls": len(permutations),
            "nonempty_subset_count": len(rows),
            "subset_orbit_replay": orbit_summary,
            "subsets": public_rows,
        }
        all_rows.extend(rows)

    representatives = []
    for seed_name in SEEDS:
        for archived in primary["seeds"][seed_name]["certified_representatives"]:
            row = next(
                candidate
                for candidate in all_rows
                if candidate["selected_walls_one_based"]
                == archived["selected_doubling_walls_one_based"]
                and archived["seed_type"] == seed_name
                and candidate in [
                    internal
                    for internal in all_rows[
                        0 if seed_name == "T1" else 15 : 15 if seed_name == "T1" else 30
                    ]
                ]
            )
            representatives.append((archived["id"], row["_relations"]))
    collisions = [
        [left_id, right_id]
        for left_index, (left_id, left) in enumerate(representatives)
        for right_id, right in representatives[left_index + 1 :]
        if _isomorphic(left, right)
    ]
    if collisions:
        raise ArithmeticError(f"primary representatives collide: {collisions}")

    return {
        "schema": SCHEMA,
        "pinned_inputs": {
            "factory_path": str(FACTORY.relative_to(ROOT)),
            "factory_sha256": FACTORY_SHA256,
            "primary_certificate_path": str(PRIMARY.relative_to(ROOT)),
            "primary_certificate_sha256": PRIMARY_SHA256,
        },
        "independence": {
            "imports_primary_certifier": False,
            "calls_primary_build": False,
            "reconstructs_boundaries_from_factory_roots": True,
            "recomputes_relation_graphs": True,
        },
        "scope": {
            "orthogonal_gluing_convexity_replayed": True,
            "all_boundary_orbits_replayed": True,
            "all_nonempty_subsets_replayed": True,
            "exact_cone_or_vinberg_replayed_here": False,
            "worldwide_novelty_claimed": False,
        },
        "seeds": seed_payloads,
        "theorem": {
            "seed_count": 2,
            "total_nonempty_subset_count": 30,
            "all_subgroups_have_order_two_to_k": True,
            "all_boundaries_subgroup_invariant": True,
            "all_glued_unions_pass_exact_orthogonal_convexity_lemma": True,
            "all_facets_have_strict_timelike_witnesses": True,
            "all_subset_diagram_classes_equal_seed_automorphism_orbits": True,
            "primary_representative_count": len(representatives),
            "primary_representatives_pairwise_nonisomorphic": True,
            "primary_representative_collision_count": 0,
        },
        "decision": "independent_hostile_replay_passes",
    }


def main() -> None:
    payload = audit()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    OUTPUT.write_text(rendered, encoding="utf-8")
    print(json.dumps(payload["theorem"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
