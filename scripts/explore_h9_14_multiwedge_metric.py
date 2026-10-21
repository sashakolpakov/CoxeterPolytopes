"""Exploratory exact metric gate for the three H9/14 multiwedges.

This is deliberately a fixed-projective-realization audit.  It constructs
one integral homogeneous facet-covector matrix for each of records 43, 44,
and 47, checks its vertex complements against the frozen incidence witness,
and measures several linear right-angle constraint systems.  Failure here
does not exclude another realization of the same combinatorial type.
"""

from __future__ import annotations

from itertools import combinations
from math import gcd

from sympy import Matrix

from coxeter_polytopes.combinatorics.h9_14_record43 import (
    TARGET_MISSING_EDGES as EDGES43,
    TARGET_MISSING_FIVES as FIVES43,
    record43_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record44 import (
    TARGET_MISSING_EDGES as EDGES44,
    TARGET_MISSING_FIVES as FIVES44,
    record44_facet_complements,
)
from coxeter_polytopes.combinatorics.h9_14_record47 import (
    TARGET_MISSING_EDGES as EDGES47,
    TARGET_MISSING_FIVES as FIVES47,
    record47_facet_complements,
)


RECORD_DATA = {
    43: (EDGES43, FIVES43, record43_facet_complements, (3, 4, 5, 6, 7)),
    44: (EDGES44, FIVES44, record44_facet_complements, (2, 3, 4, 5, 6)),
    47: (EDGES47, FIVES47, record47_facet_complements, (3, 4, 5, 6, 7)),
}


def facet_covectors(record: int) -> Matrix:
    """Return rows ``ell_i`` with chamber inequalities ``ell_i(x) >= 0``.

    The first coordinate is the homogenizing coordinate.  Records 43 and 47
    use columns ``(h,x,y,z,t1,t2,s1,s2,s3,s4)``.  Record 44 uses
    ``(h,x,y,t1,t2,t3,s1,s2,s3,s4)``.  The last four variables carry the
    product four-simplex.
    """

    if record in (43, 47):
        rows = [[0] * 10 for _ in range(14)]
        rows[0] = [-1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
        rows[1] = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
        rows[2] = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
        rows[3] = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
        rows[4] = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
        rows[5] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
        rows[6] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        rows[7] = [1, 0, 0, 0, 0, 0, -1, -1, -1, -1]
        rows[8] = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
        rows[9] = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
        rows[10] = [0, 1, 0, 0, -1, -1, 0, 0, 0, 0]
        rows[11] = [10, -1, -1, -1, 0, 0, 0, 0, 0, 0]
        if record == 43:
            rows[12] = [7, 0, -1, 0, 0, 0, 0, 0, 0, 0]
            rows[13] = [16, 0, -2, -1, 0, 0, 0, 0, 0, 0]
        else:
            rows[12] = [9, 0, -1, 0, 0, 0, 0, 0, 0, 0]
            rows[13] = [9, 0, 0, -1, 0, 0, 0, 0, 0, 0]
        return Matrix(rows)

    if record == 44:
        rows = [[0] * 10 for _ in range(14)]
        rows[0] = [2, -1, 0, 0, 0, 0, 0, 0, 0, 0]
        rows[1] = [3, -1, -1, 0, 0, 0, 0, 0, 0, 0]
        rows[2] = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
        rows[3] = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
        rows[4] = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
        rows[5] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        rows[6] = [1, 0, 0, 0, 0, 0, -1, -1, -1, -1]
        rows[7] = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
        rows[8] = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
        rows[9] = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
        rows[10] = [0, 0, 1, -1, -1, -1, 0, 0, 0, 0]
        rows[11] = [2, 0, -1, 0, 0, 0, 0, 0, 0, 0]
        rows[12] = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
        rows[13] = [-1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
        return Matrix(rows)

    raise ValueError(f"unsupported H9/14 record {record}")


def _mask(indices: tuple[int, ...]) -> int:
    return sum(1 << index for index in indices)


def vertex_complements(matrix: Matrix) -> tuple[int, ...]:
    """Enumerate exact affine vertices and return their five-set complements."""

    result = []
    full = (1 << 14) - 1
    for active in combinations(range(14), 9):
        nullspace = matrix[list(active), :].nullspace()
        if len(nullspace) != 1:
            continue
        point = nullspace[0]
        if point[0] == 0:
            continue
        if point[0] < 0:
            point = -point
        point /= point[0]
        values = matrix * point
        if any(value < 0 for value in values):
            continue
        zero_set = tuple(index for index, value in enumerate(values) if value == 0)
        if zero_set != active:
            raise AssertionError("the displayed realization is not simple")
        result.append(full ^ _mask(active))
    return tuple(sorted(result))


def primitive_kernel(matrix: Matrix) -> Matrix:
    """Return a column basis of ``ker(matrix.T)`` with primitive integers."""

    columns = []
    for vector in matrix.T.nullspace():
        denominator = 1
        for entry in vector:
            denominator = denominator * entry.q // gcd(denominator, entry.q)
        values = [int(entry * denominator) for entry in vector]
        divisor = 0
        for value in values:
            divisor = gcd(divisor, abs(value))
        values = [value // divisor for value in values]
        first = next(value for value in values if value)
        if first < 0:
            values = [-value for value in values]
        columns.append(Matrix(values))
    return Matrix.hstack(*columns)


SYMMETRIC_POSITIONS = tuple(combinations(range(10), 2)) + tuple(
    (index, index) for index in range(10)
)


def orthogonality_row(left: Matrix, right: Matrix) -> list[int]:
    """Coefficients of ``left*Q*right.T`` for symmetric ``Q``."""

    result = []
    for first, second in SYMMETRIC_POSITIONS:
        if first == second:
            result.append(int(left[first] * right[first]))
        else:
            result.append(
                int(left[first] * right[second] + left[second] * right[first])
            )
    return result


def constraint_rank(matrix: Matrix, pairs: tuple[tuple[int, int], ...]) -> int:
    equations = Matrix(
        [orthogonality_row(matrix.row(i), matrix.row(j)) for i, j in pairs]
    )
    return equations.rank()


def forced_isotropic_facets(
    matrix: Matrix, pairs: tuple[tuple[int, int], ...]
) -> tuple[int, ...]:
    """Facets whose square vanishes identically under the selected zeros."""

    equations = Matrix(
        [orthogonality_row(matrix.row(i), matrix.row(j)) for i, j in pairs]
    )
    free_forms = Matrix.hstack(*equations.nullspace())
    diagonal = Matrix(
        [orthogonality_row(matrix.row(i), matrix.row(i)) for i in range(14)]
    )
    restricted = diagonal * free_forms
    return tuple(
        index
        for index in range(14)
        if all(restricted[index, column] == 0 for column in range(restricted.cols))
    )


def audit(record: int) -> dict[str, object]:
    matrix = facet_covectors(record)
    edges, fives, target_function, simplex = RECORD_DATA[record]
    missing_edges = {tuple(sorted(edge)) for edge in edges}
    meeting = tuple(
        pair for pair in combinations(range(14), 2) if pair not in missing_edges
    )
    allowed_nonright = {
        pair
        for face in fives
        for pair in combinations(sorted(face), 2)
    }
    block_supported_right = tuple(
        pair for pair in meeting if pair not in allowed_nonright
    )
    residual = tuple(index for index in range(14) if index not in simplex)
    cross = tuple((i, j) for i in simplex for j in residual)
    kernel = primitive_kernel(matrix)
    vertices = vertex_complements(matrix)
    target = tuple(target_function())
    simplex_span = matrix[list(simplex), :].rowspace()
    residual_span = matrix[list(residual), :].rowspace()
    combined = Matrix.vstack(*simplex_span, *residual_span)
    intersection_dimension = (
        len(simplex_span) + len(residual_span) - combined.rank()
    )
    return {
        "record": record,
        "facet_covectors": [list(map(int, matrix.row(i))) for i in range(14)],
        "rank": matrix.rank(),
        "kernel_basis": [list(map(int, kernel.row(i))) for i in range(14)],
        "kernel_dimension": kernel.cols,
        "vertex_count": len(vertices),
        "vertex_complements_match_archive": vertices == target,
        "meeting_pair_count": len(meeting),
        "all_meeting_right_angle_rank": constraint_rank(matrix, meeting),
        "all_meeting_right_angle_solution_dimension": 55
        - constraint_rank(matrix, meeting),
        "block_supported_right_pair_count": len(block_supported_right),
        "block_supported_right_angle_rank": constraint_rank(
            matrix, block_supported_right
        ),
        "block_supported_right_angle_solution_dimension": 55
        - constraint_rank(matrix, block_supported_right),
        "block_supported_forced_isotropic_facets": list(
            forced_isotropic_facets(matrix, block_supported_right)
        ),
        "cross_factor_pair_count": len(cross),
        "cross_factor_right_angle_rank": constraint_rank(matrix, cross),
        "cross_factor_right_angle_solution_dimension": 55
        - constraint_rank(matrix, cross),
        "simplex_covector_span_dimension": len(simplex_span),
        "residual_covector_span_dimension": len(residual_span),
        "span_intersection_dimension": intersection_dimension,
        "cross_factor_orthogonality_forces_degenerate_form": (
            intersection_dimension > 0 and combined.rank() == 10
        ),
    }


def main() -> None:
    import json

    print(json.dumps([audit(record) for record in (43, 44, 47)], indent=2))


if __name__ == "__main__":
    main()
