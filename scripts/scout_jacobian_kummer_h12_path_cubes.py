#!/usr/bin/env python3
"""Scout exact D_infinity x C2^3 products in the JK-H12 chamber.

This is an exploratory count-only lane.  It reconstructs the pinned 172-wall
carrier, verifies every proposed generator is a doubling wall, verifies the
three cube reflections commute with the path factor, and counts the exact
primitive boundary-root orbits for both path phases.
"""

from __future__ import annotations

from collections import defaultdict
import math

import scripts.certify_jacobian_kummer_h12_single_doublings as base
# All indices are zero-based.  These were found by the exact doubling-wall
# graph scout; each row is a parallel pair orthogonal to a 3-cube.
CANDIDATES = (
    ((75, 71), (35, 38, 47)),
    ((74, 70), (35, 41, 50)),
    ((73, 69), (38, 41, 55)),
    ((72, 68), (47, 50, 55)),
    ((67, 62), (30, 57, 58)),
    ((66, 63), (31, 57, 59)),
    ((65, 64), (32, 58, 59)),
    ((61, 60), (30, 31, 32)),
)


Matrix = tuple[tuple[int, ...], ...]
Vector = tuple[int, ...]


def _identity(rank: int) -> Matrix:
    return tuple(
        tuple(1 if row == column else 0 for column in range(rank))
        for row in range(rank)
    )


def _matmul(left: Matrix, right: Matrix) -> Matrix:
    rank = len(left)
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column] for middle in range(rank))
            for column in range(rank)
        )
        for row in range(rank)
    )


def _matvec(matrix: Matrix, vector: Vector) -> Vector:
    return tuple(
        sum(entry * value for entry, value in zip(row, vector, strict=True))
        for row in matrix
    )


def _reflection_matrix(gram: base.Matrix, root: base.IntegralVector) -> Matrix:
    rank = len(root)
    columns = tuple(
        base._reflection_apply(
            gram,
            root,
            tuple(1 if row == column else 0 for row in range(rank)),
        )
        for column in range(rank)
    )
    if any(value.denominator != 1 for column in columns for value in column):
        raise ArithmeticError("a proposed generator reflection is nonintegral")
    matrix = tuple(
        tuple(int(columns[column][row]) for column in range(rank))
        for row in range(rank)
    )
    if _matmul(matrix, matrix) != _identity(rank):
        raise ArithmeticError("a proposed generator reflection is not involutive")
    return matrix


def _path_transforms(
    reflections: dict[int, Matrix],
    pair: tuple[int, int],
    copies: int,
    start: int,
) -> tuple[Matrix, ...]:
    transforms = [_identity(len(next(iter(reflections.values()))))]
    for edge in range(copies - 1):
        transforms.append(
            _matmul(transforms[-1], reflections[pair[(start + edge) % 2]])
        )
    return tuple(transforms)


def _cube_transforms(
    reflections: dict[int, Matrix], cube: tuple[int, ...]
) -> tuple[Matrix, ...]:
    rank = len(next(iter(reflections.values())))
    transforms: list[Matrix] = []
    for mask in range(1 << len(cube)):
        transform = _identity(rank)
        for bit, wall in enumerate(cube):
            if mask & (1 << bit):
                transform = _matmul(transform, reflections[wall])
        transforms.append(transform)
    if len(set(transforms)) != len(transforms):
        raise ArithmeticError("a proposed cube action has a kernel")
    return tuple(transforms)


def _product(left: tuple[Matrix, ...], right: tuple[Matrix, ...]) -> tuple[Matrix, ...]:
    transforms = tuple(_matmul(a, b) for a in left for b in right)
    if len(set(transforms)) != len(transforms):
        raise ArithmeticError("a proposed product action has a collision")
    return transforms


def _primitive_key(vector: Vector) -> tuple[int, ...]:
    content = math.gcd(*(abs(value) for value in vector if value))
    primitive = tuple(value // content for value in vector)
    first = next(value for value in primitive if value)
    return primitive if first > 0 else tuple(-value for value in primitive)


def _boundary_count(
    roots: tuple[Vector, ...],
    reflections: dict[int, Matrix],
    transforms: tuple[Matrix, ...],
    generators: tuple[int, ...],
) -> int:
    transform_keys = set(transforms)
    internal: set[tuple[int, ...]] = set()
    for transform in transforms:
        for index in generators:
            if _matmul(transform, reflections[index]) in transform_keys:
                internal.add(_primitive_key(_matvec(transform, roots[index])))

    boundary: set[tuple[int, ...]] = set()
    for transform in transforms:
        for root in roots:
            key = _primitive_key(_matvec(transform, root))
            if key not in internal:
                boundary.add(key)
    return len(boundary)


def _affine_formula(values: list[int]) -> tuple[int, int] | None:
    if len(values) < 2:
        return None
    slope = values[1] - values[0]
    intercept = values[0] - slope
    if all(value == slope * index + intercept for index, value in enumerate(values, 1)):
        return slope, intercept
    return None


def main() -> None:
    _source, gram_rows, root_rows, _interior = base._load_source()
    roots = tuple(tuple(root) for root in root_rows)
    relations = base._seed_relations(gram_rows, root_rows)
    doubling = set(base._doubling_walls(relations))
    selected_stock = {index for pair, cube in CANDIDATES for index in (*pair, *cube)}
    reflections = {
        index: _reflection_matrix(gram_rows, roots[index]) for index in selected_stock
    }

    grouped: dict[tuple[tuple[int, ...], tuple[int, ...]], list[list[int]]] = defaultdict(list)
    for pair, cube in CANDIDATES:
        selected = (*pair, *cube)
        if not set(selected) <= doubling:
            raise ArithmeticError("a proposed generator is not a doubling wall")
        if not (
            base._pair(gram_rows, roots[pair[0]], roots[pair[0]]) == -4
            and base._pair(gram_rows, roots[pair[1]], roots[pair[1]]) == -4
            and base._pair(gram_rows, roots[pair[0]], roots[pair[1]]) == 4
        ):
            raise ArithmeticError("a proposed path pair is not parallel")
        cube_transforms = _cube_transforms(reflections, cube)
        path_reflections = tuple(reflections[index] for index in pair)
        if any(
            _matmul(path_reflection, cube_transform)
            != _matmul(cube_transform, path_reflection)
            for path_reflection in path_reflections
            for cube_transform in cube_transforms
        ):
            raise ArithmeticError("path and cube factors do not commute")

        phase_values: list[list[int]] = []
        for start in (0, 1):
            values: list[int] = []
            for copies in range(1, 9):
                transforms = _product(
                    _path_transforms(reflections, pair, copies, start),
                    cube_transforms,
                )
                if len(transforms) != 8 * copies:
                    raise ArithmeticError("D_infinity x C2^3 action collided")
                values.append(
                    _boundary_count(roots, reflections, transforms, selected)
                )
            phase_values.append(values)

        formulas = tuple(_affine_formula(values) for values in phase_values)
        key = (tuple(formulas), tuple(phase_values[0]))
        grouped[key].append(
            [
                [index + 1 for index in pair],
                [index + 1 for index in cube],
                phase_values[0],
                phase_values[1],
            ]
        )

    for (formulas, _values), rows in grouped.items():
        print(f"formulas={formulas}")
        for pair, cube, phase_zero, phase_one in rows:
            print(
                f"  pair={pair} cube={cube} phase0={phase_zero} phase1={phase_one}"
            )


if __name__ == "__main__":
    main()
