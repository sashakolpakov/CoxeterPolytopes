#!/usr/bin/env python3
"""Exact tiny certificate for the smooth-toric full-nef Coxeter factory.

The finite classification in this file is deliberately small.  A smooth
complete toric surface with nef anticanonical class is reduced equivariantly
to P2, F0, or F2.  Reversing the contractions changes the cyclic boundary
self-intersection sequence by

    (..., a, b, ...) -> (..., a-1, -1, b-1, ...),

and nefness of -K permits the step exactly when a,b >= -1.  The resulting
closure has 16 dihedral sequence classes and stops at nine boundary rays.

For the unique nine-ray endpoint we independently work in
NS = I_(1,6), check all nine primitive wall classes, integral reflections,
Coxeter products, irredundancy, and all 84 possible rank-six active sets.
No floating point geometry or polyhedral package is used.
"""

from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from math import factorial, gcd
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "toric_full_nef_coxeter_factory.json"
)
SCHEMA = "coxeter.toric_full_nef_factory.v1"

# Minimal smooth complete toric fans compatible with D_i^2 >= -2.
# F1 is not minimal and is generated from P2 by one blowup.
SEEDS = (
    (1, 1, 1),             # P2
    (0, 0, 0, 0),          # F0
    (-2, 0, 2, 0),         # F2
)


def canonical_cycle(values: Iterable[int]) -> tuple[int, ...]:
    """Canonicalize a cyclic integer word up to the dihedral action."""

    cycle = tuple(int(value) for value in values)
    reverse = cycle[::-1]
    candidates = [
        word[index:] + word[:index]
        for word in (cycle, reverse)
        for index in range(len(cycle))
    ]
    return min(candidates)


def toric_blowups(cycle: tuple[int, ...]) -> list[tuple[int, ...]]:
    """Return every -K-nef torus-fixed blowup, modulo no symmetry."""

    result: list[tuple[int, ...]] = []
    length = len(cycle)
    for left in range(length):
        right = (left + 1) % length
        if cycle[left] < -1 or cycle[right] < -1:
            continue
        blown = list(cycle)
        if right == 0:
            blown[-1] -= 1
            blown[0] -= 1
            blown.append(-1)
        else:
            blown[left] -= 1
            blown[right] -= 1
            blown.insert(right, -1)
        result.append(canonical_cycle(blown))
    return result


def enumerate_weak_fano_cycles() -> list[tuple[int, ...]]:
    """Enumerate the complete 16-element weak-Fano fan list."""

    seen = {canonical_cycle(seed) for seed in SEEDS}
    queue = deque(sorted(seen))
    while queue:
        cycle = queue.popleft()
        for child in toric_blowups(cycle):
            if child not in seen:
                seen.add(child)
                queue.append(child)
    return sorted(seen, key=lambda cycle: (len(cycle), cycle))


def negative_wall_gram(cycle: tuple[int, ...]) -> sp.Matrix:
    """Intersection Gram of the negative invariant boundary curves."""

    indices = [index for index, square in enumerate(cycle) if square < 0]
    gram = sp.zeros(len(indices))
    for row, left in enumerate(indices):
        gram[row, row] = cycle[left]
        for column, right in enumerate(indices):
            if row == column:
                continue
            if (left - right) % len(cycle) in (1, len(cycle) - 1):
                gram[row, column] = 1
    return gram


def exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Return exact (positive, negative, zero) inertia by congruence."""

    work = sp.Matrix(matrix)
    positive = negative = zero = 0
    while work.rows:
        size = work.rows
        diagonal = next(
            (index for index in range(size) if work[index, index] != 0),
            None,
        )
        if diagonal is not None:
            order = [diagonal, *[i for i in range(size) if i != diagonal]]
            work = work.extract(order, order)
            pivot = sp.Rational(work[0, 0])
            if pivot > 0:
                positive += 1
            else:
                negative += 1
            if size == 1:
                break
            column = work[1:, 0]
            work = sp.simplify(
                work[1:, 1:] - column * column.T / pivot
            )
            continue

        off_diagonal = next(
            (
                (row, column)
                for row in range(size)
                for column in range(row + 1, size)
                if work[row, column] != 0
            ),
            None,
        )
        if off_diagonal is None:
            zero += size
            break
        row, column = off_diagonal
        order = [
            row,
            column,
            *[i for i in range(size) if i not in {row, column}],
        ]
        work = work.extract(order, order)
        block = work[:2, :2]
        positive += 1
        negative += 1
        if size == 2:
            break
        cross = work[2:, :2]
        work = sp.simplify(
            work[2:, 2:] - cross * block.inv() * cross.T
        )
    return positive, negative, zero


def coxeter_label_histogram(cycle: tuple[int, ...]) -> dict[str, int]:
    """Read labels from toric intersections of negative boundary curves."""

    indices = [index for index, square in enumerate(cycle) if square < 0]
    histogram = {"2": 0, "3": 0, "4": 0, "infinity": 0}
    for offset, left in enumerate(indices):
        for right in indices[:offset]:
            adjacent = (left - right) % len(cycle) in (
                1,
                len(cycle) - 1,
            )
            if not adjacent:
                histogram["2"] += 1
                continue
            norm_product = (-cycle[left]) * (-cycle[right])
            label = {1: "infinity", 2: "4", 4: "3"}[norm_product]
            histogram[label] += 1
    return histogram


# Endpoint coordinates are (H,E1,...,E6), with H^2=1, Ei^2=-1.
J = sp.diag(1, -1, -1, -1, -1, -1, -1)
ENDPOINT_WALL_ORDER = (
    "E1-E4",
    "E4",
    "H-E1-E2-E4",
    "E2-E5",
    "E5",
    "H-E2-E3-E5",
    "E3-E6",
    "E6",
    "H-E3-E1-E6",
)
ENDPOINT_WALLS: dict[str, tuple[int, ...]] = {
    "E1-E4": (0, 1, 0, 0, -1, 0, 0),
    "E4": (0, 0, 0, 0, 1, 0, 0),
    "H-E1-E2-E4": (1, -1, -1, 0, -1, 0, 0),
    "E2-E5": (0, 0, 1, 0, 0, -1, 0),
    "E5": (0, 0, 0, 0, 0, 1, 0),
    "H-E2-E3-E5": (1, 0, -1, -1, 0, -1, 0),
    "E3-E6": (0, 0, 0, 1, 0, 0, -1),
    "E6": (0, 0, 0, 0, 0, 0, 1),
    "H-E3-E1-E6": (1, -1, 0, -1, 0, 0, -1),
}
ENDPOINT_CANONICAL = (-3, 1, 1, 1, 1, 1, 1)
ENDPOINT_INTERIOR = (7, -2, -2, -2, -1, -1, -1)
ENDPOINT_FAN_RAYS = (
    (1, 0),
    (0, 1),
    (-1, 1),
    (-2, 1),
    (-3, 1),
    (-1, 0),
    (1, -1),
    (3, -2),
    (2, -1),
)
ENDPOINT_ANTICANONICAL_POINTS = (
    (-1, -1),  # X
    (0, -1),   # Y
    (0, 0),    # W
    (1, 2),    # Z
)


def _vector(values: Iterable[int | sp.Rational]) -> sp.Matrix:
    return sp.Matrix(list(values))


def pair(left: Iterable[int], right: Iterable[int]) -> sp.Rational:
    return sp.Rational((_vector(left).T * J * _vector(right))[0])


def square(vector: Iterable[int]) -> sp.Rational:
    return pair(vector, vector)


def endpoint_wall_matrix() -> sp.Matrix:
    return sp.Matrix.hstack(
        *(_vector(ENDPOINT_WALLS[name]) for name in ENDPOINT_WALL_ORDER)
    )


def endpoint_gram() -> sp.Matrix:
    walls = endpoint_wall_matrix()
    return walls.T * J * walls


def _det2(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def endpoint_fan_self_intersections() -> tuple[int, ...]:
    length = len(ENDPOINT_FAN_RAYS)
    return tuple(
        -_det2(
            ENDPOINT_FAN_RAYS[(index - 1) % length],
            ENDPOINT_FAN_RAYS[(index + 1) % length],
        )
        for index in range(length)
    )


def endpoint_anticanonical_lattice_points() -> tuple[tuple[int, int], ...]:
    """Enumerate the anticanonical polygon using exact proved bounds.

    Rays (1,0), (-1,0) force -1 <= a <= 1.  Ray (0,1) forces
    b >= -1, while (3,-2) and a <= 1 force b <= 2.
    """

    return tuple(
        (a, b)
        for a in range(-1, 2)
        for b in range(-1, 3)
        if all(a * x + b * y >= -1 for x, y in ENDPOINT_FAN_RAYS)
    )


def reflection_matrix(wall: Iterable[int]) -> sp.Matrix:
    root = _vector(wall)
    return sp.eye(7) - sp.Rational(2, square(root)) * root * root.T * J


def _primitive_integral(vector: sp.Matrix) -> tuple[int, ...]:
    scale = sp.ilcm(
        *[int(sp.Rational(value).q) for value in vector]
    )
    values = [int(sp.Rational(value) * scale) for value in vector]
    common = 0
    for value in values:
        common = gcd(common, abs(value))
    if common:
        values = [value // common for value in values]
    first_nonzero = next((value for value in values if value), 0)
    if first_nonzero < 0:
        values = [-value for value in values]
    return tuple(values)


def endpoint_extreme_rays() -> list[dict[str, Any]]:
    """Enumerate all endpoint rays from the 84 six-wall active sets."""

    walls = [_vector(ENDPOINT_WALLS[name]) for name in ENDPOINT_WALL_ORDER]
    rays: dict[tuple[int, ...], dict[str, Any]] = {}
    for active_subset in combinations(range(len(walls)), 6):
        equations = sp.Matrix.vstack(
            *((walls[index].T * J) for index in active_subset)
        )
        if equations.rank() != 6:
            continue
        nullspace = equations.nullspace()
        if len(nullspace) != 1:
            continue
        ray = _primitive_integral(nullspace[0])
        values = [pair(ray, wall) for wall in walls]
        if all(value <= 0 for value in values):
            ray = tuple(-value for value in ray)
            values = [-value for value in values]
        if any(value < 0 for value in values):
            continue
        ray_square = int(square(ray))
        rays[ray] = {
            "vector": list(ray),
            "square": ray_square,
            "kind": "ideal" if ray_square == 0 else "finite",
            "active_facets": [
                ENDPOINT_WALL_ORDER[index]
                for index, value in enumerate(values)
                if value == 0
            ],
        }
    return [rays[ray] for ray in sorted(rays)]


def endpoint_cusp_audits() -> list[dict[str, Any]]:
    """Literal maximal-parabolic audits for every lightlike extreme ray."""

    gram = endpoint_gram()
    names = list(ENDPOINT_WALL_ORDER)
    result: list[dict[str, Any]] = []
    for ray in endpoint_extreme_rays():
        if ray["kind"] != "ideal":
            continue
        active = [names.index(name) for name in ray["active_facets"]]
        active_gram = gram.extract(active, active)
        unseen = set(range(len(active)))
        components: list[list[int]] = []
        while unseen:
            seed = min(unseen)
            unseen.remove(seed)
            stack = [seed]
            component: list[int] = []
            while stack:
                current = stack.pop()
                component.append(current)
                neighbours = {
                    other
                    for other in unseen
                    if active_gram[current, other] != 0
                }
                unseen -= neighbours
                stack.extend(sorted(neighbours, reverse=True))
            components.append(sorted(component))
        component_rows = []
        for component in components:
            block = active_gram.extract(component, component)
            component_rows.append(
                {
                    "facets": [
                        names[active[index]] for index in component
                    ],
                    "node_count": len(component),
                    "coxeter_rank": int(block.rank()),
                    "gram_inertia": list(exact_inertia(block)),
                    "affine": exact_inertia(block)
                    == (0, len(component) - 1, 1),
                }
            )
        inactive_pairings = [
            pair(ray["vector"], ENDPOINT_WALLS[name])
            for name in names
            if name not in ray["active_facets"]
        ]
        result.append(
            {
                "ray": ray["vector"],
                "active_facet_count": len(active),
                "active_gram_rank": int(active_gram.rank()),
                "active_gram_inertia": list(exact_inertia(active_gram)),
                "components": component_rows,
                "component_rank_sum": sum(
                    row["coxeter_rank"] for row in component_rows
                ),
                "outside_pairings_strictly_positive": all(
                    value > 0 for value in inactive_pairings
                ),
                "maximal_parabolic_h6_cusp": (
                    len(active) == 7
                    and exact_inertia(active_gram) == (0, 5, 2)
                    and sorted(
                        row["coxeter_rank"] for row in component_rows
                    ) == [2, 3]
                    and all(row["affine"] for row in component_rows)
                    and all(value > 0 for value in inactive_pairings)
                ),
            }
        )
    return result


def _matrix_key(matrix: sp.Matrix, modulus: int | None = None) -> tuple[int, ...]:
    values = (int(value) for value in matrix)
    if modulus is None:
        return tuple(values)
    return tuple(value % modulus for value in values)


def _matrix_from_key(key: tuple[int, ...]) -> sp.Matrix:
    return sp.Matrix(7, 7, key)


@lru_cache(maxsize=1)
def endpoint_a2_cubed_group_keys() -> tuple[tuple[int, ...], ...]:
    """Generate the finite group on the six endpoint (-2)-roots exactly."""

    generators = [
        reflection_matrix(ENDPOINT_WALLS[name])
        for name in ENDPOINT_WALL_ORDER
        if square(ENDPOINT_WALLS[name]) == -2
    ]
    identity = _matrix_key(sp.eye(7))
    seen = {identity}
    queue = deque([identity])
    while queue:
        current = _matrix_from_key(queue.popleft())
        for generator in generators:
            candidate = _matrix_key(current * generator)
            if candidate in seen:
                continue
            seen.add(candidate)
            queue.append(candidate)
            if len(seen) > 216:
                raise ArithmeticError("the endpoint (-2)-reflection group exceeded 216")
    return tuple(sorted(seen))


def standard_cubic_line_roots() -> tuple[tuple[int, ...], ...]:
    """The 6+15+6 standard (-1)-classes on a marked cubic surface."""

    roots: set[tuple[int, ...]] = set()
    for index in range(6):
        exceptional = [0] * 7
        exceptional[index + 1] = 1
        roots.add(tuple(exceptional))
    for left, right in combinations(range(6), 2):
        line = [1] + [0] * 6
        line[left + 1] = -1
        line[right + 1] = -1
        roots.add(tuple(line))
    for omitted in range(6):
        conic = [2] + [-1] * 6
        conic[omitted + 1] = 0
        roots.add(tuple(conic))
    return tuple(sorted(roots))


def endpoint_congruence_audit() -> dict[str, Any]:
    """Audit the A2^3 gluing group and its orbit of the three (-1)-walls."""

    group_keys = endpoint_a2_cubed_group_keys()
    mod_two_image = {
        tuple(value % 2 for value in key)
        for key in group_keys
    }
    minus_two_roots = [
        ENDPOINT_WALLS[name]
        for name in ENDPOINT_WALL_ORDER
        if square(ENDPOINT_WALLS[name]) == -2
    ]
    minus_two_adjacencies = sorted(
        (left, right)
        for left, right in combinations(range(len(minus_two_roots)), 2)
        if pair(minus_two_roots[left], minus_two_roots[right]) != 0
    )
    seeds = [
        ENDPOINT_WALLS[name]
        for name in ENDPOINT_WALL_ORDER
        if square(ENDPOINT_WALLS[name]) == -1
    ]
    orbits: list[set[tuple[int, ...]]] = []
    for seed in seeds:
        vector = _vector(seed)
        orbit = {
            tuple(int(value) for value in _matrix_from_key(key) * vector)
            for key in group_keys
        }
        orbits.append(orbit)
    orbit_union = set().union(*orbits)
    cubic_lines = set(standard_cubic_line_roots())
    return {
        "minus_two_generator_count": 6,
        "finite_group_type": "W(A2)^3",
        "minus_two_diagram_is_three_disjoint_A2_edges": (
            len(minus_two_adjacencies) == 3
            and len({node for edge in minus_two_adjacencies for node in edge}) == 6
            and all(
                pair(minus_two_roots[left], minus_two_roots[right]) == 1
                for left, right in minus_two_adjacencies
            )
        ),
        "finite_group_order": len(group_keys),
        "mod_two_image_order": len(mod_two_image),
        "minus_one_seed_count": len(seeds),
        "minus_one_orbit_sizes": [len(orbit) for orbit in orbits],
        "minus_one_orbits_are_disjoint": sum(map(len, orbits)) == len(orbit_union),
        "orbit_union_count": len(orbit_union),
        "orbit_union_is_standard_6_plus_15_plus_6_cubic_lines": (
            orbit_union == cubic_lines
        ),
        "all_cubic_line_roots_have_square_minus_one": all(
            square(root) == -1 for root in cubic_lines
        ),
        "all_cubic_line_reflections_are_identity_mod_two": all(
            _matrix_key(reflection_matrix(root), 2) == _matrix_key(sp.eye(7), 2)
            for root in cubic_lines
        ),
    }


def _finite_component_order(component: list[int], gram: sp.Matrix) -> int:
    """Return the finite Coxeter order for a connected endpoint subdiagram."""

    if len(component) == 1:
        return 2
    adjacency = {
        node: [
            other
            for other in component
            if other != node and gram[node, other] != 0
        ]
        for node in component
    }
    if any(len(neighbours) > 2 for neighbours in adjacency.values()):
        raise ArithmeticError("unexpected branched endpoint subdiagram")
    endpoints = sorted(node for node in component if len(adjacency[node]) == 1)
    if len(endpoints) != 2:
        raise ArithmeticError("unexpected non-path finite endpoint subdiagram")
    ordered = [endpoints[0]]
    previous: int | None = None
    while len(ordered) < len(component):
        candidates = [
            node for node in adjacency[ordered[-1]] if node != previous
        ]
        if len(candidates) != 1:
            raise ArithmeticError("endpoint path ordering failed")
        previous, current = ordered[-1], candidates[0]
        ordered.append(current)
    labels: list[int] = []
    for left, right in zip(ordered, ordered[1:]):
        norm_product = int(gram[left, left] * gram[right, right])
        labels.append({2: 4, 4: 3}[norm_product])
    rank = len(component)
    if all(label == 3 for label in labels):
        return factorial(rank + 1)  # A_rank
    if labels.count(4) == 1 and labels[0] == 4 and all(
        label == 3 for label in labels[1:]
    ):
        return (2 ** rank) * factorial(rank)  # B_rank=C_rank
    if labels.count(4) == 1 and labels[-1] == 4 and all(
        label == 3 for label in labels[:-1]
    ):
        return (2 ** rank) * factorial(rank)  # B_rank=C_rank
    if labels == [3, 4, 3]:
        return 1152  # F4
    raise ArithmeticError(f"unrecognized finite path labels {labels}")


def _spherical_subgroup_order(subset: tuple[int, ...], gram: sp.Matrix) -> int | None:
    if not subset:
        return 1
    block = gram.extract(subset, subset)
    if exact_inertia(block) != (0, len(subset), 0):
        return None
    unseen = set(subset)
    components: list[list[int]] = []
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        stack = [seed]
        component: list[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            neighbours = {
                other for other in unseen if gram[current, other] != 0
            }
            unseen -= neighbours
            stack.extend(sorted(neighbours, reverse=True))
        components.append(sorted(component))
    order = 1
    for component in components:
        order *= _finite_component_order(component, gram)
    return order


def endpoint_solomon_euler() -> dict[str, Any]:
    """Compute the endpoint Coxeter-group Euler characteristic exactly."""

    gram = endpoint_gram()
    euler = Fraction(0)
    spherical = 0
    subset_count = 0
    for size in range(10):
        for subset in combinations(range(9), size):
            subset_count += 1
            order = _spherical_subgroup_order(subset, gram)
            if order is None:
                continue
            spherical += 1
            euler += Fraction((-1) ** size, order)
    return {
        "subset_count": subset_count,
        "spherical_subset_count": spherical,
        "euler_characteristic": str(euler),
        "gauss_bonnet_factor_absolute": "8*pi^3/15",
        "orbifold_volume": "pi^3/3240",
    }


def p6_clique_euler() -> dict[str, Any]:
    """Compute the right-angled P6 Euler characteristic from 27 line roots."""

    roots = list(standard_cubic_line_roots())
    adjacency = [
        {right for right in range(len(roots)) if right != left and pair(roots[left], roots[right]) == 0}
        for left in range(len(roots))
    ]
    counts: dict[int, int] = {0: 1}

    def extend(candidates: list[int], depth: int) -> None:
        for offset, vertex in enumerate(candidates):
            counts[depth + 1] = counts.get(depth + 1, 0) + 1
            extend(
                [
                    other
                    for other in candidates[offset + 1 :]
                    if other in adjacency[vertex]
                ],
                depth + 1,
            )

    extend(list(range(len(roots))), 0)
    euler = sum(
        (Fraction(-1, 2) ** size) * count
        for size, count in counts.items()
    )
    return {
        "clique_counts_by_size": {
            str(size): counts[size] for size in sorted(counts)
        },
        "euler_characteristic": str(euler),
        "gauss_bonnet_factor_absolute": "8*pi^3/15",
        "orbifold_volume": "pi^3/15",
    }


def endpoint_facet_witnesses() -> list[dict[str, Any]]:
    """Project one strict interior point to every wall."""

    interior = _vector(ENDPOINT_INTERIOR)
    result: list[dict[str, Any]] = []
    for name in ENDPOINT_WALL_ORDER:
        wall = _vector(ENDPOINT_WALLS[name])
        witness = interior - pair(interior, wall) / square(wall) * wall
        other_values = [
            pair(witness, ENDPOINT_WALLS[other])
            for other in ENDPOINT_WALL_ORDER
            if other != name
        ]
        result.append(
            {
                "facet": name,
                "witness": [str(value) for value in witness],
                "witness_square": str(square(witness)),
                "facet_pairing": str(pair(witness, wall)),
                "minimum_other_pairing": str(min(other_values)),
            }
        )
    return result


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _canonical_sha256(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _classification_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for cycle in enumerate_weak_fano_cycles():
        gram = negative_wall_gram(cycle)
        ray_count = len(cycle)
        picard_rank = ray_count - 2
        hyperbolic_dimension = picard_rank - 1
        rows.append(
            {
                "boundary_cycle": list(cycle),
                "fan_ray_count": ray_count,
                "picard_rank": picard_rank,
                "hyperbolic_dimension": hyperbolic_dimension,
                "anticanonical_square": 12 - ray_count,
                "negative_wall_count": sum(value < 0 for value in cycle),
                "wall_gram_rank": int(gram.rank()),
                "wall_gram_inertia": list(exact_inertia(gram)),
                "coxeter_label_histogram": coxeter_label_histogram(cycle),
            }
        )
    return rows


def build_certificate() -> dict[str, Any]:
    cycles = enumerate_weak_fano_cycles()
    endpoint_cycle = canonical_cycle((-2, -1, -2) * 3)
    assert len(cycles) == 16
    assert max(map(len, cycles)) == 9
    assert [cycle for cycle in cycles if len(cycle) == 9] == [endpoint_cycle]

    wall_matrix = endpoint_wall_matrix()
    gram = endpoint_gram()
    smith = smith_normal_form(wall_matrix, domain=ZZ)
    nonzero_smith = [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index] != 0
    ]
    rays = endpoint_extreme_rays()
    cusps = endpoint_cusp_audits()
    congruence = endpoint_congruence_audit()
    endpoint_euler = endpoint_solomon_euler()
    p6_euler = p6_clique_euler()
    witnesses = endpoint_facet_witnesses()
    reflections = [
        reflection_matrix(ENDPOINT_WALLS[name])
        for name in ENDPOINT_WALL_ORDER
    ]
    interior_pairings = [
        int(pair(ENDPOINT_INTERIOR, ENDPOINT_WALLS[name]))
        for name in ENDPOINT_WALL_ORDER
    ]

    endpoint_checks = {
        "fan_is_smooth_complete_with_expected_cycle": (
            all(
                _det2(
                    ENDPOINT_FAN_RAYS[index],
                    ENDPOINT_FAN_RAYS[(index + 1) % len(ENDPOINT_FAN_RAYS)],
                ) == 1
                for index in range(len(ENDPOINT_FAN_RAYS))
            )
            and endpoint_fan_self_intersections()
            == tuple(
                int(square(ENDPOINT_WALLS[name]))
                for name in ENDPOINT_WALL_ORDER
            )
        ),
        "anticanonical_polygon_has_xyz_equals_w_cubed_relation": (
            endpoint_anticanonical_lattice_points()
            == ENDPOINT_ANTICANONICAL_POINTS
            and tuple(
                ENDPOINT_ANTICANONICAL_POINTS[0][index]
                + ENDPOINT_ANTICANONICAL_POINTS[1][index]
                + ENDPOINT_ANTICANONICAL_POINTS[3][index]
                for index in range(2)
            )
            == tuple(3 * value for value in ENDPOINT_ANTICANONICAL_POINTS[2])
        ),
        "wall_squares_match_boundary_cycle": sorted(
            int(square(ENDPOINT_WALLS[name]))
            for name in ENDPOINT_WALL_ORDER
        ) == sorted(endpoint_cycle),
        "wall_span_is_full_and_primitive": (
            wall_matrix.rank() == 7 and nonzero_smith == [1] * 7
        ),
        "gram_has_lorentzian_quotient": exact_inertia(gram) == (1, 6, 2),
        "canonical_square_is_three": square(ENDPOINT_CANONICAL) == 3,
        "boundary_sum_is_anticanonical": tuple(
            sum(ENDPOINT_WALLS[name][index] for name in ENDPOINT_WALL_ORDER)
            for index in range(7)
        ) == tuple(-value for value in ENDPOINT_CANONICAL),
        "adjunction_holds": all(
            pair(ENDPOINT_CANONICAL, ENDPOINT_WALLS[name])
            == -2 - square(ENDPOINT_WALLS[name])
            for name in ENDPOINT_WALL_ORDER
        ),
        "strict_interior": (
            square(ENDPOINT_INTERIOR) > 0
            and all(value > 0 for value in interior_pairings)
        ),
        "all_reflections_integral_lattice_isometries": all(
            all(value.q == 1 for value in reflection)
            and reflection.det() == -1
            and reflection * reflection == sp.eye(7)
            and reflection.T * J * reflection == J
            for reflection in reflections
        ),
        "all_nine_facets_irredundant": (
            len(witnesses) == 9
            and all(row["facet_pairing"] == "0" for row in witnesses)
            and all(sp.Rational(row["minimum_other_pairing"]) > 0 for row in witnesses)
        ),
        "complete_ray_audit_has_no_hyperideal_ray": (
            len(rays) == 21
            and sum(row["kind"] == "finite" for row in rays) == 18
            and sum(row["kind"] == "ideal" for row in rays) == 3
            and all(row["square"] >= 0 for row in rays)
        ),
        "all_ideal_rays_are_literal_maximal_parabolic_cusps": (
            len(cusps) == 3
            and all(row["maximal_parabolic_h6_cusp"] for row in cusps)
        ),
        "a2_cubed_group_glues_exactly_216_endpoint_chambers": (
            congruence["minus_two_diagram_is_three_disjoint_A2_edges"]
            and congruence["finite_group_order"] == 216
            and congruence["mod_two_image_order"] == 216
            and congruence["minus_one_orbit_sizes"] == [9, 9, 9]
            and congruence["minus_one_orbits_are_disjoint"]
            and congruence[
                "orbit_union_is_standard_6_plus_15_plus_6_cubic_lines"
            ]
            and congruence["all_cubic_line_roots_have_square_minus_one"]
            and congruence[
                "all_cubic_line_reflections_are_identity_mod_two"
            ]
        ),
        "solomon_and_p6_volumes_have_ratio_216": (
            endpoint_euler["subset_count"] == 512
            and endpoint_euler["spherical_subset_count"] == 322
            and endpoint_euler["euler_characteristic"] == "-1/1728"
            and p6_euler["clique_counts_by_size"]
            == {
                "0": 1,
                "1": 27,
                "2": 216,
                "3": 720,
                "4": 1080,
                "5": 648,
                "6": 72,
            }
            and p6_euler["euler_characteristic"] == "-1/8"
        ),
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "classification": {
            "minimal_seeds": [list(canonical_cycle(seed)) for seed in SEEDS],
            "blowup_rule": "adjacent (a,b) -> (a-1,-1,b-1), allowed iff a,b>=-1",
            "dihedral_cycle_count": len(cycles),
            "fan_ray_histogram": {
                str(length): sum(len(cycle) == length for cycle in cycles)
                for length in sorted({len(cycle) for cycle in cycles})
            },
            "maximum_fan_rays": 9,
            "maximum_picard_rank": 7,
            "maximum_hyperbolic_dimension": 6,
            "unique_maximal_cycle": list(endpoint_cycle),
            "rows": _classification_rows(),
        },
        "theorem_gates": {
            "integral_full_nef_facets_iff": [
                "every negative invariant boundary curve has square -1 or -2",
                "every invariant boundary curve has square at least -2",
                "-K_X is nef",
            ],
            "coxeter_angle_rule": {
                "disjoint boundary curves": "m=2",
                "(-2)-(-2) transverse adjacency": "m=3",
                "(-1)-(-2) transverse adjacency": "m=4",
                "(-1)-(-1) transverse adjacency": "parallel, m=infinity",
            },
            "finite_volume_gate": "rho>=3 and complete rational-polyhedral nef equality",
            "compactness_gate": "no isotropic nef extremal ray",
            "toric_compactness_verdict": "noncompact for every rho>=2",
            "reflexive_polygon_bound": "B(Delta)+B(Delta_dual)=12 and B(Delta_dual)>=3, so B(Delta)<=9",
        },
        "h6_endpoint": {
            "construction": "blow up a perfect matching of three torus fixed points on the six-cycle boundary of Bl_3(P2)",
            "lattice": "I_(1,6)",
            "basis": ["H", "E1", "E2", "E3", "E4", "E5", "E6"],
            "fan_rays": [list(ray) for ray in ENDPOINT_FAN_RAYS],
            "fan_self_intersection_cycle": list(endpoint_fan_self_intersections()),
            "anticanonical_lattice_points": [
                list(point) for point in endpoint_anticanonical_lattice_points()
            ],
            "anticanonical_coordinate_order": ["X", "Y", "W", "Z"],
            "anticanonical_affine_relation": "X+Y+Z=3W on exponents",
            "anticanonical_image_equation": "X*Y*Z=W^3 in P3",
            "contracted_root_configuration": "3A2 (six (-2)-curves in three adjacent pairs)",
            "intersection_matrix": _matrix_payload(J),
            "canonical_class": list(ENDPOINT_CANONICAL),
            "strict_interior": list(ENDPOINT_INTERIOR),
            "wall_order": list(ENDPOINT_WALL_ORDER),
            "walls": {
                name: list(ENDPOINT_WALLS[name])
                for name in ENDPOINT_WALL_ORDER
            },
            "wall_gram": _matrix_payload(gram),
            "wall_gram_rank": int(gram.rank()),
            "wall_gram_inertia": list(exact_inertia(gram)),
            "wall_span_smith_nonzero": nonzero_smith,
            "coxeter_label_histogram": coxeter_label_histogram(
                tuple(int(square(ENDPOINT_WALLS[name])) for name in ENDPOINT_WALL_ORDER)
            ),
            "active_subset_trials": 84,
            "extreme_rays": rays,
            "cusp_audits": cusps,
            "congruence_and_smoothing_gluing": {
                "exact_finite_computation": congruence,
                "endpoint_solomon_audit": endpoint_euler,
                "classical_P6_right_angled_clique_audit": p6_euler,
                "external_identification_input": (
                    "Everitt--Ratcliffe--Tschantz identify the 27 standard "
                    "(-1)-root reflections with the principal level-two "
                    "reflection group Gamma_2^6 and its chamber P6"
                ),
                "group_theorem": (
                    "Gamma_2^6 is a normal subgroup of the nine-wall "
                    "endpoint reflection group with quotient W(A2)^3"
                ),
                "endpoint_group_as_congruence_preimage": (
                    "G_toric=rho_2^{-1}(rho_2(W(A2)^3)) inside the ERT "
                    "arithmetic ambient"
                ),
                "index": 216,
                "chamber_union": "P6 is the union of 216 endpoint chambers",
                "orbifold_volume_ratio_P6_over_endpoint": 216,
                "ag_interpretation": (
                    "smoothing XYZ=W^3 removes the three A2 exceptional "
                    "chains; the exact common-lattice wall crossing glues "
                    "the endpoint chambers by W(A2)^3 to the smooth-cubic "
                    "27-line nef chamber P6"
                ),
            },
            "finite_vertex_count": sum(row["kind"] == "finite" for row in rays),
            "ideal_vertex_count": sum(row["kind"] == "ideal" for row in rays),
            "hyperideal_vertex_count": sum(row["square"] < 0 for row in rays),
            "facet_witnesses": witnesses,
            "checks": endpoint_checks,
            "decision": "finite_volume_noncompact",
        },
        "scope": {
            "automatic_factory": "exactly smooth projective toric weak del Pezzo surfaces in rho>=3",
            "inverse_test": "enumerate all twelve rho>=3 toric divisor presentations and test an exact saturated lattice/wall isometry; a chosen cyclic word determines its fan, but a Gram-only shortcut is not assumed",
            "flagship_inverse_status": "scoped exact recognizer/classification, not L4 universality: the denominator is the toric-completable locus and no theorem puts an independently broad arithmetic Coxeter class inside it",
            "ag_leverage_scope": "nonzero only inside the declared smooth-toric full-nef completion species; arbitrary and compact Coxeter polytopes remain uncovered",
            "abstract_novelty": "not new: literal match to Im Hof's dimension-six twice-truncated-simplex/Napier-cycle diagram",
            "abstract_lineage": "H.-C. Im Hof, Napier cycles and hyperbolic Coxeter groups (1990), dimension-six bottom-right diagram in the compact synopsis",
            "geometric_novelty": "project-new exact toric full-nef realization; worldwide priority for the marking is not claimed",
            "outside_scope": [
                "nonintegral geometric reflection mirrors",
                "proper invariant or eigenspace sections",
                "singular toric surfaces and toric stacks",
                "divisor subsystems not equal to the complete nef cone",
            ],
        },
    }
    payload["all_checks_pass"] = all(endpoint_checks.values())
    digest_payload = dict(payload)
    payload["certificate_sha256"] = _canonical_sha256(digest_payload)
    return payload


def verify_archive() -> dict[str, Any]:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    rebuilt = build_certificate()
    if archived != rebuilt:
        raise RuntimeError("frozen toric certificate does not match exact replay")
    endpoint = rebuilt["h6_endpoint"]
    return {
        "certificate_sha256": rebuilt["certificate_sha256"],
        "cycle_count": rebuilt["classification"]["dihedral_cycle_count"],
        "maximum_hyperbolic_dimension": rebuilt["classification"]["maximum_hyperbolic_dimension"],
        "facet_count": len(endpoint["wall_order"]),
        "finite_vertex_count": endpoint["finite_vertex_count"],
        "ideal_vertex_count": endpoint["ideal_vertex_count"],
        "decision": endpoint["decision"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="write frozen JSON")
    parser.add_argument(
        "--check-archive",
        action="store_true",
        help="rebuild exactly and compare with frozen JSON",
    )
    args = parser.parse_args()

    if args.check_archive:
        print(json.dumps(verify_archive(), indent=2, sort_keys=True))
        return

    payload = build_certificate()
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "all_checks_pass": payload["all_checks_pass"],
                "certificate_sha256": payload["certificate_sha256"],
                "cycle_count": payload["classification"]["dihedral_cycle_count"],
                "maximum_hyperbolic_dimension": payload["classification"]["maximum_hyperbolic_dimension"],
                "h6_endpoint": {
                    "facets": len(payload["h6_endpoint"]["wall_order"]),
                    "finite": payload["h6_endpoint"]["finite_vertex_count"],
                    "ideal": payload["h6_endpoint"]["ideal_vertex_count"],
                    "hyperideal": payload["h6_endpoint"]["hyperideal_vertex_count"],
                },
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
