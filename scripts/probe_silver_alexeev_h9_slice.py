#!/usr/bin/env python3
"""Probe the silver eigenslice of Alexeev's ``(20,2,1)`` chamber.

This is deliberately a probe rather than a certificate.  It constructs the
integral self-adjoint square-root-of-two action supplied by

    saturated C5(4,3,3,3,3) + F4 + A1,

in the first twenty simple-root coordinates of Alexeev's attached 50-wall
Gram matrix.  A regular positive eigenvector is then moved into Alexeev's
chamber by exact integral root reflections.  Only after this Weyl framing do
we restrict the fifty inequalities to the ten-dimensional positive
eigenspace.

The source file is the v4 arXiv ancillary file ``all_cases_in_2el_paper.txt``.
Its digest and the exact chamber data are certified elsewhere in the
repository.  Floating point is used here only for finding a framing word and
for exploratory linear programming; every reflection matrix and the framed
endomorphism are checked over the integers.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import re
from typing import Sequence

import numpy as np
from scipy.optimize import linprog
from sympy import Matrix

from coxeter_polytopes.geometric.silver_finite import (
    silver_finite_trace_lift,
    silver_lanner_blocks,
)
from coxeter_polytopes.lattices import e8_cartan_gram


IntegerMatrix = np.ndarray


def _as_integer(matrix: Matrix | np.ndarray) -> IntegerMatrix:
    rows = matrix.tolist() if hasattr(matrix, "tolist") else matrix
    result: list[list[int]] = []
    for row in rows:
        integer_row: list[int] = []
        for entry in row:
            value = Fraction(entry)
            if value.denominator != 1:
                raise AssertionError("an expected integral matrix is fractional")
            integer_row.append(value.numerator)
        result.append(integer_row)
    return np.array(result, dtype=object)


def _alexeev_matrix(path: Path) -> Matrix:
    source = path.read_text(encoding="utf-8")
    match = re.search(
        r"^coxeter_matrix\[\(20,2,1\)\] = matrix\((\[.*\])\)$",
        source,
        flags=re.MULTILINE,
    )
    if match is None:
        raise ValueError("the (20,2,1) Coxeter matrix was not found")
    return Matrix(ast.literal_eval(match.group(1)))


def _multiplication_sqrt2(rank: int) -> IntegerMatrix:
    result = np.zeros((2 * rank, 2 * rank), dtype=object)
    for index in range(rank):
        result[2 * index, 2 * index + 1] = 2
        result[2 * index + 1, 2 * index] = 1
    return result


def _saturated_c5_action() -> tuple[IntegerMatrix, IntegerMatrix]:
    """Return the saturated C5 trace Gram and its integral sqrt(2) action."""

    free_gram = np.array(silver_lanner_blocks()[1].trace_gram, dtype=object)
    saturation = np.zeros((10, 10), dtype=object)
    for index in range(9):
        saturation[index, index] = 1
    for index in (3, 5, 9):
        saturation[index, 9] = Fraction(1, 2)
    saturation_matrix = Matrix(saturation.tolist())
    gram = _as_integer(saturation_matrix.T * Matrix(free_gram.tolist()) * saturation_matrix)
    action = _as_integer(
        saturation_matrix.inv()
        * Matrix(_multiplication_sqrt2(5).tolist())
        * saturation_matrix
    )
    return gram, action


# Columns zero and one are a U-pair.  The last eight columns generate its
# E8(-1) complement in the saturated C5 lattice.
_C5_U_E8_PREBASIS = np.array(
    [
        (-1, -1, 1, 0, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, -1, -1, 0),
        (-1, -1, 0, 1, 0, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, -1, 0, 0, 0, -1),
        (-1, -1, 0, 0, 1, 0, 0, 0, 0, 0),
        (0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
        (-1, -1, 0, 0, 0, 0, 1, 0, 0, 0),
        (0, 0, 0, 0, 0, 0, 0, 1, 0, 0),
        (-1, 0, 0, 0, 0, 0, 0, 0, 1, 0),
        (0, -1, 0, 0, 0, 0, 0, 0, 0, 1),
    ],
    dtype=object,
)

# This unimodular complement change gives the repository's standard E8 Gram.
_C5_E8_CHANGE = np.array(
    [
        (-2, 0, 0, 0, -1, 1, 0, 2),
        (-2, 0, 0, -1, 0, 2, 0, 2),
        (-1, 1, 0, 0, 0, -1, 0, -1),
        (-1, 0, 0, 0, 0, 1, -1, 1),
        (-1, 0, 1, -1, 0, 0, -1, 0),
        (0, 0, 0, 0, 0, 0, 1, 0),
        (-1, 0, 0, 0, 0, 0, -1, 1),
        (1, 0, 0, -1, 0, 0, 2, 0),
    ],
    dtype=object,
)

# The Alexeev E10 block has its affine-E8 null root in column zero and the
# adjacent overextended root in column one; the remaining columns are an E8
# complement.  The row order is [0,9,10,...,16,19] in the published basis.
_ALEXEEV_U_E8_PREBASIS = np.array(
    [
        (0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
        (3, 3, 1, 0, 0, 0, 0, 0, 0, 0),
        (2, 2, 0, 1, 0, 0, 0, 0, 0, 0),
        (4, 4, 0, 0, 1, 0, 0, 0, 0, 0),
        (6, 6, 0, 0, 0, 1, 0, 0, 0, 0),
        (5, 5, 0, 0, 0, 0, 1, 0, 0, 0),
        (4, 4, 0, 0, 0, 0, 0, 1, 0, 0),
        (3, 3, 0, 0, 0, 0, 0, 0, 1, 0),
        (2, 2, 0, 0, 0, 0, 0, 0, 0, 1),
        (1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
    ],
    dtype=object,
)


def _block_standardizer() -> IntegerMatrix:
    result = np.zeros((10, 10), dtype=object)
    result[:2, :2] = np.eye(2, dtype=object)
    result[2:, 2:] = _C5_E8_CHANGE
    return result


def _alexeev_e10_standardizer() -> IntegerMatrix:
    # The complement diagram is the standard E8 diagram after this cyclic
    # relabeling of its eight basis vectors.
    permutation = (1, 2, 3, 4, 5, 6, 7, 0)
    change = np.zeros((10, 10), dtype=object)
    change[:2, :2] = np.eye(2, dtype=object)
    for target, source in enumerate(permutation):
        change[2 + source, 2 + target] = 1
    return _ALEXEEV_U_E8_PREBASIS @ change


def silver_action_in_alexeev_basis(
    chamber: Matrix,
) -> tuple[IntegerMatrix, IntegerMatrix]:
    """Construct the integral silver action in the published 20-root basis."""

    gram = np.array(chamber[:20, :20].tolist(), dtype=object)
    standard_e8 = -np.array(e8_cartan_gram(), dtype=object)

    c5_gram, c5_action = _saturated_c5_action()
    c5_standard = _C5_U_E8_PREBASIS @ _block_standardizer()
    alexeev_standard = _alexeev_e10_standardizer()
    hyperbolic_indices = (0, *range(9, 17), 19)
    alexeev_hyperbolic = gram[np.ix_(hyperbolic_indices, hyperbolic_indices)]
    target = np.zeros((10, 10), dtype=object)
    target[:2, :2] = ((0, 1), (1, 0))
    target[2:, 2:] = standard_e8
    assert np.array_equal(c5_standard.T @ c5_gram @ c5_standard, target)
    assert np.array_equal(
        alexeev_standard.T @ alexeev_hyperbolic @ alexeev_standard,
        target,
    )
    c5_to_alexeev = _as_integer(
        Matrix(alexeev_standard.tolist()) * Matrix(c5_standard.tolist()).inv()
    )
    hyperbolic_action = _as_integer(
        Matrix(c5_to_alexeev.tolist())
        * Matrix(c5_action.tolist())
        * Matrix(c5_to_alexeev.tolist()).inv()
    )

    # The published indices 1,...,8 form an E8(-1).  Transport the explicit
    # F4 trace action through the stored unimodular F4 -> E8 isometry.
    finite_indices = tuple(range(1, 9))
    alexeev_e8 = gram[np.ix_(finite_indices, finite_indices)]
    finite_permutation = next(
        permutation
        for permutation in permutations(range(8))
        if np.array_equal(
            alexeev_e8[np.ix_(permutation, permutation)],
            standard_e8,
        )
    )
    alexeev_e8_standard = np.eye(8, dtype=object)[:, finite_permutation]
    f4_lift = silver_finite_trace_lift("F4")
    f4_to_alexeev = _as_integer(
        Matrix(alexeev_e8_standard.tolist())
        * Matrix(f4_lift.trace_isometry).inv()
    )
    finite_action = _as_integer(
        Matrix(f4_to_alexeev.tolist())
        * Matrix(_multiplication_sqrt2(4).tolist())
        * Matrix(f4_to_alexeev.tolist()).inv()
    )

    # On the last A1 trace block, changing sqrt(2)e to -sqrt(2)e changes the
    # free trace Gram's -2 off-diagonal to Alexeev's +2 off-diagonal.
    a1_indices = (17, 18)
    a1_action = np.array(((0, -2), (-1, 0)), dtype=object)

    action = np.zeros((20, 20), dtype=object)
    for indices, block in (
        (hyperbolic_indices, hyperbolic_action),
        (finite_indices, finite_action),
        (a1_indices, a1_action),
    ):
        for row, ambient_row in enumerate(indices):
            for column, ambient_column in enumerate(indices):
                action[ambient_row, ambient_column] = block[row, column]

    identity = np.eye(20, dtype=object)
    assert np.array_equal(action @ action, 2 * identity)
    assert np.array_equal(action.T @ gram, gram @ action)
    assert np.array_equal(action.T @ gram @ action, 2 * gram)
    return gram, action


def _wall_coordinates(chamber: Matrix) -> IntegerMatrix:
    basis_gram = chamber[:20, :20]
    inverse = basis_gram.inv()
    columns = []
    for index in range(chamber.rows):
        coordinates = inverse * chamber[:20, index]
        columns.append(_as_integer(coordinates).reshape(20))
    roots = np.column_stack(columns)
    gram = np.array(basis_gram.tolist(), dtype=object)
    assert np.array_equal(roots.T @ gram @ roots, np.array(chamber.tolist(), dtype=object))
    return roots


def _positive_eigenbasis(action: IntegerMatrix) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eig(np.array(action, dtype=float))
    selected = np.abs(eigenvalues - np.sqrt(2.0)) < 1e-7
    basis = np.real(eigenvectors[:, selected])
    if basis.shape != (20, 10):
        raise AssertionError("the positive silver eigenspace has the wrong dimension")
    return np.linalg.qr(basis)[0]


def _reflection(root: np.ndarray, gram: IntegerMatrix, square: int) -> IntegerMatrix:
    reflection = np.eye(20, dtype=object) - Fraction(2, square) * np.outer(
        root,
        root @ gram,
    )
    result = _as_integer(reflection)
    assert np.array_equal(result.T @ gram @ result, gram)
    return result


def weyl_frame(
    gram: IntegerMatrix,
    action: IntegerMatrix,
    roots: IntegerMatrix,
    wall_squares: Sequence[int],
    *,
    max_steps: int,
) -> tuple[IntegerMatrix, IntegerMatrix, int]:
    """Move one regular positive silver eigenvector into the root chamber."""

    plus = _positive_eigenbasis(action)
    restricted_gram = plus.T @ np.array(gram, dtype=float) @ plus
    values, vectors = np.linalg.eigh(restricted_gram)
    if values[-1] <= 0 or values[-2] >= 0:
        raise AssertionError("the positive eigenspace is not Lorentzian")
    coefficients = vectors[:, -1].copy()
    # A deterministic generic perturbation avoids landing on a mirror.
    coefficients += sum(
        ((index + 1) / 2000.0) * vectors[:, index] for index in range(9)
    )
    point = plus @ coefficients

    # The dual-to-the-first-twenty-roots vector is a strict chamber point.
    # Flip the eigenvector when necessary so both are in the same component
    # of the positive cone; without this step the reflection reduction runs
    # toward the chamber in the opposite sheet.
    chamber_point = np.linalg.solve(
        np.array(gram, dtype=float),
        np.ones(20),
    )
    if point @ np.array(gram, dtype=float) @ chamber_point < 0:
        point = -point
    if point @ np.array(gram, dtype=float) @ point <= 0:
        raise AssertionError("the framing point is not timelike")

    word: list[int] = []
    weyl = np.eye(20, dtype=object)
    real_gram = np.array(gram, dtype=float)
    real_roots = np.array(roots, dtype=float)
    normalization = np.sqrt(-np.array(wall_squares, dtype=float))
    for _ in range(max_steps):
        pairings = point @ real_gram @ real_roots
        scores = pairings / normalization
        wall = int(np.argmin(scores))
        if scores[wall] >= -1e-9:
            break
        reflection = _reflection(roots[:, wall], gram, wall_squares[wall])
        point = np.array(reflection, dtype=float) @ point
        weyl = reflection @ weyl
        word.append(wall)
    else:
        raise RuntimeError("Weyl framing did not terminate at the declared bound")

    framed_action = _as_integer(
        Matrix(weyl.tolist()) * Matrix(action.tolist()) * Matrix(weyl.tolist()).inv()
    )
    assert np.array_equal(framed_action @ framed_action, 2 * np.eye(20, dtype=object))
    assert np.array_equal(framed_action.T @ gram, gram @ framed_action)
    return framed_action, np.array(word, dtype=int), len(word)


def _projective_wall_matches(
    action: IntegerMatrix,
    roots: IntegerMatrix,
) -> list[tuple[int, int, str]]:
    lookup = {tuple(roots[:, index]): index for index in range(roots.shape[1])}
    lookup.update(
        {tuple(-entry for entry in roots[:, index]): index for index in range(roots.shape[1])}
    )
    matches = []
    for index in range(roots.shape[1]):
        image = action @ roots[:, index]
        target = lookup.get(tuple(image))
        label = "tau"
        if target is None and all(int(entry) % 2 == 0 for entry in image):
            target = lookup.get(tuple(int(entry) // 2 for entry in image))
            label = "tau/2"
        if target is not None:
            matches.append((index, target, label))
    return matches


def probe(path: Path, *, max_steps: int) -> dict[str, object]:
    chamber = _alexeev_matrix(path)
    gram, action = silver_action_in_alexeev_basis(chamber)
    roots = _wall_coordinates(chamber)
    squares = [int(chamber[index, index]) for index in range(chamber.rows)]
    framed, word, steps = weyl_frame(
        gram,
        action,
        roots,
        squares,
        max_steps=max_steps,
    )

    plus = _positive_eigenbasis(framed)
    restricted_pairings = plus.T @ np.array(gram, dtype=float) @ np.array(roots, dtype=float)
    feasibility = linprog(
        np.zeros(10),
        A_ub=-restricted_pairings.T,
        b_ub=-np.ones(chamber.rows),
        bounds=[(None, None)] * 10,
        method="highs",
    )
    restricted_signs = []
    for index in range(chamber.rows):
        root = roots[:, index]
        q_value = int(root @ gram @ root)
        b_value = int(root @ gram @ framed @ root)
        value = 2 * q_value + np.sqrt(2.0) * b_value
        restricted_signs.append(-1 if value < 0 else 1 if value > 0 else 0)

    return {
        "framing_steps": steps,
        "framing_word": word.tolist(),
        "framed_action_max_abs_entry": max(abs(int(entry)) for entry in framed.flat),
        "restricted_linear_feasible": bool(feasibility.success),
        "restricted_lp_point_square": (
            float(
                (plus @ feasibility.x)
                @ np.array(gram, dtype=float)
                @ (plus @ feasibility.x)
            )
            if feasibility.success
            else None
        ),
        "restricted_normal_sign_counts": dict(Counter(restricted_signs)),
        "projectively_matched_simple_walls": _projective_wall_matches(framed, roots),
    }


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--max-steps", type=int, default=100_000)
    args = parser.parse_args(argv)
    print(probe(args.source, max_steps=args.max_steps))


if __name__ == "__main__":
    main()
