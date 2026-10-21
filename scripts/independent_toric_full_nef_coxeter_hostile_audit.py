#!/usr/bin/env python3
"""Independent tiny replay for the smooth-toric H6/9 endpoint.

This file deliberately does not import the production toric certifier.  It
rebuilds the weak-Fano boundary-cycle census from toric blowups, then checks
the endpoint in I_(1,6), including the complete affine neighborhoods of its
three cusps.  All arithmetic is exact and the largest search has 84 cases.
"""

from __future__ import annotations

from collections import Counter, deque
from itertools import combinations
import json
from math import gcd

import sympy as sp


Q = sp.diag(1, -1, -1, -1, -1, -1, -1)
WALLS = (
    (0, 1, 0, 0, -1, 0, 0),
    (0, 0, 0, 0, 1, 0, 0),
    (1, -1, -1, 0, -1, 0, 0),
    (0, 0, 1, 0, 0, -1, 0),
    (0, 0, 0, 0, 0, 1, 0),
    (1, 0, -1, -1, 0, -1, 0),
    (0, 0, 0, 1, 0, 0, -1),
    (0, 0, 0, 0, 0, 0, 1),
    (1, -1, 0, -1, 0, 0, -1),
)


def canonical_cycle(values: tuple[int, ...]) -> tuple[int, ...]:
    reverse = values[::-1]
    return min(
        word[index:] + word[:index]
        for word in (values, reverse)
        for index in range(len(values))
    )


def blowups(values: tuple[int, ...]) -> set[tuple[int, ...]]:
    children: set[tuple[int, ...]] = set()
    size = len(values)
    for left in range(size):
        right = (left + 1) % size
        if values[left] < -1 or values[right] < -1:
            continue
        child = list(values)
        if right:
            child[left] -= 1
            child[right] -= 1
            child.insert(right, -1)
        else:
            child[-1] -= 1
            child[0] -= 1
            child.append(-1)
        children.add(canonical_cycle(tuple(child)))
    return children


def cycle_census() -> set[tuple[int, ...]]:
    seeds = {
        canonical_cycle((1, 1, 1)),
        canonical_cycle((0, 0, 0, 0)),
        canonical_cycle((-2, 0, 2, 0)),
    }
    result = set(seeds)
    queue = deque(seeds)
    while queue:
        for child in blowups(queue.popleft()):
            if child not in result:
                result.add(child)
                queue.append(child)
    return result


def vector(values: tuple[int, ...]) -> sp.Matrix:
    return sp.Matrix(values)


def pair(left: tuple[int, ...], right: tuple[int, ...]) -> sp.Rational:
    return sp.Rational((vector(left).T * Q * vector(right))[0])


def primitive(values: sp.Matrix) -> tuple[int, ...]:
    denominator = sp.ilcm(*(sp.Rational(value).q for value in values))
    integers = [int(sp.Rational(value) * denominator) for value in values]
    divisor = 0
    for value in integers:
        divisor = gcd(divisor, abs(value))
    return tuple(value // divisor for value in integers)


def orient_feasible(ray: tuple[int, ...]) -> tuple[int, ...] | None:
    products = [pair(ray, wall) for wall in WALLS]
    if all(value >= 0 for value in products):
        return ray
    if all(value <= 0 for value in products):
        return tuple(-value for value in ray)
    return None


def inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Exact inertia from the signs in an LDL congruence decomposition."""

    work = sp.Matrix(matrix)
    positive = negative = zero = 0
    while work.rows:
        diagonal = next(
            (index for index in range(work.rows) if work[index, index]), None
        )
        if diagonal is None:
            edge = next(
                (
                    (row, column)
                    for row in range(work.rows)
                    for column in range(row + 1, work.rows)
                    if work[row, column]
                ),
                None,
            )
            if edge is None:
                zero += work.rows
                break
            order = [*edge, *[i for i in range(work.rows) if i not in edge]]
            work = work.extract(order, order)
            block = work[:2, :2]
            positive += 1
            negative += 1
            if work.rows == 2:
                break
            cross = work[2:, :2]
            work = work[2:, 2:] - cross * block.inv() * cross.T
            continue
        order = [diagonal, *[i for i in range(work.rows) if i != diagonal]]
        work = work.extract(order, order)
        pivot = work[0, 0]
        if pivot > 0:
            positive += 1
        else:
            negative += 1
        if work.rows == 1:
            break
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / pivot
    return positive, negative, zero


def connected_components(matrix: sp.Matrix) -> list[list[int]]:
    unseen = set(range(matrix.rows))
    components: list[list[int]] = []
    while unseen:
        stack = [unseen.pop()]
        component: list[int] = []
        while stack:
            node = stack.pop()
            component.append(node)
            neighbors = {
                other
                for other in unseen
                if matrix[node, other] != 0
            }
            unseen -= neighbors
            stack.extend(neighbors)
        components.append(sorted(component))
    return sorted(components, key=lambda item: (len(item), item))


def endpoint_replay() -> dict[str, object]:
    wall_matrix = sp.Matrix.hstack(*(vector(wall) for wall in WALLS))
    gram = wall_matrix.T * Q * wall_matrix
    rays: dict[tuple[int, ...], tuple[int, tuple[int, ...]]] = {}
    for active in combinations(range(9), 6):
        equations = sp.Matrix.vstack(
            *((vector(WALLS[index]).T * Q) for index in active)
        )
        if equations.rank() != 6:
            continue
        raw = primitive(equations.nullspace()[0])
        ray = orient_feasible(raw)
        if ray is None:
            continue
        products = [pair(ray, wall) for wall in WALLS]
        square = int(pair(ray, ray))
        rays[ray] = (
            square,
            tuple(index for index, value in enumerate(products) if value == 0),
        )

    kinds = Counter(
        "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        for square, _ in rays.values()
    )
    cusp_rows: list[dict[str, object]] = []
    for ray, (square, active) in sorted(rays.items()):
        if square:
            continue
        neighborhood = gram.extract(active, active)
        components = connected_components(neighborhood)
        component_data = []
        for component in components:
            block = neighborhood.extract(component, component)
            component_data.append(
                {
                    "nodes": len(component),
                    "rank": int(block.rank()),
                    "inertia": list(inertia(block)),
                }
            )
        cusp_rows.append(
            {
                "ray": list(ray),
                "active_wall_count": len(active),
                "neighborhood_rank": int(neighborhood.rank()),
                "neighborhood_inertia": list(inertia(neighborhood)),
                "affine_components": component_data,
            }
        )

    return {
        "wall_gram_rank": int(gram.rank()),
        "wall_gram_inertia": list(inertia(gram)),
        "active_subset_trials": 84,
        "ray_count": len(rays),
        "ray_kinds": dict(sorted(kinds.items())),
        "cusps": cusp_rows,
    }


def main() -> None:
    cycles = cycle_census()
    histogram = Counter(map(len, cycles))
    endpoint = canonical_cycle((-2, -1, -2) * 3)
    assert len(cycles) == 16
    assert histogram == {3: 1, 4: 3, 5: 2, 6: 4, 7: 2, 8: 3, 9: 1}
    assert [cycle for cycle in cycles if len(cycle) == 9] == [endpoint]

    replay = endpoint_replay()
    assert replay["wall_gram_inertia"] == [1, 6, 2]
    assert replay["ray_kinds"] == {"finite": 18, "ideal": 3}
    for cusp in replay["cusps"]:
        assert cusp["active_wall_count"] == 7
        assert cusp["neighborhood_rank"] == 5
        assert cusp["neighborhood_inertia"] == [0, 5, 2]
        assert sorted(
            (row["nodes"], row["rank"], row["inertia"])
            for row in cusp["affine_components"]
        ) == [(3, 2, [0, 2, 1]), (4, 3, [0, 3, 1])]

    print(
        json.dumps(
            {
                "classification": {
                    "cycle_count": len(cycles),
                    "fan_ray_histogram": dict(sorted(histogram.items())),
                    "unique_nine_ray_cycle": list(endpoint),
                },
                "endpoint": replay,
                "verdict": "independent exact replay passes",
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
