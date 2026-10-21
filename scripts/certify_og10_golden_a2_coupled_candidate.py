#!/usr/bin/env python3
"""Exact carrier certificate for the A2-coupled golden OG10 candidate.

The field lattice is

    [[-2*phi, 1, 0], [1, 2, 1], [0, 1, 2]] + E8.

The last positive line is coupled to the discriminant-11 binary block rather
than appended orthogonally.  The script certifies the field signatures, the
rank-22 trace lattice, an explicit primitive gluing into the OG10 Beauville
lattice, and a complete projected-wall theorem for the natural fixed H9.

The coupling defeats the *formula* behind the old ultraparallel proof: a
nonfixed norm-two pair has Coxeter label 5.  Nevertheless its fold is an
existing fixed norm-two wall.  Exact two-embedding bounds show that this is
true for every ambient norm-two wall meeting the fixed slice.  The slice is
therefore a complete chamber of the old L9 arrangement, whose archived
infinite-volume theorem rules out this carrier as a finite-volume example.
No Vinberg search is required.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence

from sympy import Matrix, Rational, eye
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from coxeter_polytopes.exact.golden import Golden, ONE, PHI, ZERO
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    integral_gram,
    is_codifferent,
    is_even_integral,
    trace_gram,
)
from coxeter_polytopes.lattices import block_diagonal, e8_cartan_gram, hyperbolic_plane
from scripts.coupled_discriminant11_candidate import multiplication_by_phi


ROOT = Path(__file__).resolve().parents[1]
OLD_FIXED_CHAMBER_CERTIFICATE = (
    ROOT / "data/candidates/og10_golden_h9_infinite_stabilizer.json"
)


def _negative(matrix: Sequence[Sequence[int]]) -> list[list[int]]:
    return [[-int(entry) for entry in row] for row in matrix]


def _sympy(matrix: Iterable[Iterable[Any]]) -> Matrix:
    return Matrix([list(row) for row in matrix])


def _integer_matrix(matrix: Matrix) -> Matrix:
    if not all(entry.is_Integer is True for entry in matrix):
        raise ArithmeticError("an expected integral matrix is nonintegral")
    return matrix.applyfunc(int)


def _rational_text(value: Any) -> str:
    numerator, denominator = value.as_numer_denom()
    numerator = int(numerator)
    denominator = int(denominator)
    return str(numerator) if denominator == 1 else f"{numerator}/{denominator}"


def _matrix_text(matrix: Matrix) -> list[list[str]]:
    return [
        [_rational_text(matrix[i, j]) for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]


def _integer_matrix_list(matrix: Matrix) -> list[list[int]]:
    integral = _integer_matrix(matrix)
    return [
        [int(integral[i, j]) for j in range(integral.cols)]
        for i in range(integral.rows)
    ]


def _golden_matrix_text(matrix: Sequence[Sequence[Golden]]) -> list[list[str]]:
    return [[str(entry) for entry in row] for row in matrix]


def _golden_vector_text(vector: Sequence[Golden]) -> list[str]:
    return [str(entry) for entry in vector]


def _golden_matmul(
    left: Sequence[Sequence[Golden]],
    right: Sequence[Sequence[Golden]],
) -> list[list[Golden]]:
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                ZERO,
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def _golden_transpose(
    matrix: Sequence[Sequence[Golden]],
) -> list[list[Golden]]:
    return [list(column) for column in zip(*matrix)]


def _golden_pair(
    left: Sequence[Golden],
    gram: Sequence[Sequence[Golden]],
    right: Sequence[Golden],
) -> Golden:
    return sum(
        (
            left[i] * gram[i][j] * right[j]
            for i in range(len(left))
            for j in range(len(right))
        ),
        ZERO,
    )


def _golden_action(
    matrix: Sequence[Sequence[Golden]],
    vector: Sequence[Golden],
) -> tuple[Golden, ...]:
    return tuple(
        sum((matrix[i][j] * vector[j] for j in range(len(vector))), ZERO)
        for i in range(len(matrix))
    )


def _inertia_tuple(matrix: Sequence[Sequence[Any]]) -> tuple[int, int, int]:
    value = inertia(matrix)
    return value.positive, value.negative, value.zero


def _smith_nonunits(matrix: Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[i, i]))
        for i in range(min(smith.rows, smith.cols))
        if abs(int(smith[i, i])) > 1
    ]


def field_gram() -> list[list[Golden]]:
    """Return the rank-11 A2-coupled field form."""

    result = [[ZERO for _ in range(11)] for _ in range(11)]
    ternary = (
        (-2 * PHI, ONE, ZERO),
        (ONE, Golden(2), ONE),
        (ZERO, ONE, Golden(2)),
    )
    for i in range(3):
        for j in range(3):
            result[i][j] = ternary[i][j]
    e8 = e8_cartan_gram()
    for i in range(8):
        for j in range(8):
            result[3 + i][3 + j] = Golden(e8[i][j])
    return result


def ternary_trace_gram() -> Matrix:
    ternary = [row[:3] for row in field_gram()[:3]]
    return _sympy(integral_gram(trace_gram(ternary, -PHI * INV_SQRT5)))


def trace_split_change() -> Matrix:
    """Split the eight golden E8 coordinates into two E8(-1) blocks."""

    change = Matrix.zeros(22)
    change[:6, :6] = eye(6)
    for index in range(8):
        old_e = 6 + 2 * index
        old_phi_e = old_e + 1
        change[old_e, 6 + index] = 1
        change[old_e, 14 + index] = -1
        change[old_phi_e, 14 + index] = 1
    return change


def expected_split_trace_gram() -> Matrix:
    negative_e8 = _negative(e8_cartan_gram())
    return _sympy(
        block_diagonal(
            _integer_matrix_list(ternary_trace_gram()),
            negative_e8,
            negative_e8,
        )
    )


def complement_gram() -> Matrix:
    return Matrix.diag(2, 30)


def gluing_vectors() -> tuple[Matrix, Matrix]:
    """Order-2 and order-10 generators of the isotropic gluing graph."""

    order_two = Matrix(
        [
            Rational(-9, 2),
            3,
            3,
            -2,
            Rational(-3, 2),
            1,
            Rational(1, 2),
            Rational(1, 2),
        ]
    )
    order_ten = Matrix(
        [
            Rational(-33, 10),
            Rational(21, 10),
            Rational(11, 5),
            Rational(-7, 5),
            Rational(-11, 10),
            Rational(7, 10),
            Rational(1, 2),
            Rational(2, 5),
        ]
    )
    return order_two, order_ten


def overlattice_basis() -> Matrix:
    """Columns form a basis of the order-20 glued rank-eight lattice."""

    return Matrix(
        [
            [1, 0, 0, 0, 0, 0, Rational(1, 2), Rational(4, 5)],
            [0, 1, 0, 0, 0, 0, Rational(1, 2), Rational(9, 10)],
            [0, 0, 1, 0, 0, 0, 0, Rational(4, 5)],
            [0, 0, 0, 1, 0, 0, 0, Rational(2, 5)],
            [0, 0, 0, 0, 1, 0, Rational(1, 2), Rational(3, 5)],
            [0, 0, 0, 0, 0, 1, Rational(1, 2), Rational(3, 10)],
            [0, 0, 0, 0, 0, 0, Rational(1, 2), 0],
            [0, 0, 0, 0, 0, 0, 0, Rational(1, 10)],
        ]
    )


def small_og10_isometry() -> Matrix:
    """Take the glued Gram to 3U + A2(-1), column convention."""

    return Matrix(
        [
            [0, 0, 0, 1, 2, 3, -5, 2],
            [0, 0, 0, 0, 1, 1, 0, 0],
            [0, 0, 1, 0, 1, 0, 0, 1],
            [0, 0, 0, 1, 1, 1, -2, 0],
            [0, 1, -1, 1, 1, 2, -2, 1],
            [0, -1, 0, -1, 0, 1, 0, 0],
            [0, 0, 1, 1, 0, -1, 0, 0],
            [1, 1, -1, 0, -1, -2, 2, 0],
        ]
    )


def standard_small_og10_gram() -> Matrix:
    a2_negative = Matrix([[-2, 1], [1, -2]])
    return Matrix.diag(
        _sympy(hyperbolic_plane()),
        _sympy(hyperbolic_plane()),
        _sympy(hyperbolic_plane()),
        a2_negative,
    )


def natural_involution() -> list[list[Golden]]:
    """The involution fixing the field hyperplane with middle coordinate 0."""

    result = [[ZERO for _ in range(11)] for _ in range(11)]
    block = (
        (ONE, ONE - PHI, ZERO),
        (ZERO, -ONE, ZERO),
        (ZERO, ONE, ONE),
    )
    for i in range(3):
        for j in range(3):
            result[i][j] = block[i][j]
    for i in range(3, 11):
        result[i][i] = ONE
    return result


def bugaenko_reversal_on_fixed_lattice() -> list[list[Golden]]:
    """Return ``A8 + 1`` in fixed coordinates ``(x,u,E8)``.

    The matrix is the independently certified nontrivial reversal of
    Bugaenko's eleven-wall diagram.  It is first displayed in the historical
    coordinate order ``(x,E8,u)`` and then permuted to this script's fixed
    coordinate order.
    """

    old_order = [
        [Golden(9, 16), ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, Golden(-2, -3), ZERO],
        [Golden(30, 50), ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE, Golden(-6, -10), ZERO],
        [Golden(60, 100), ZERO, ZERO, ZERO, ZERO, ZERO, ONE, ZERO, Golden(-12, -20), ZERO],
        [Golden(90, 150), ZERO, ZERO, ZERO, ZERO, ONE, ZERO, ZERO, Golden(-18, -30), ZERO],
        [Golden(72, 120), ZERO, ZERO, ZERO, ONE, ZERO, ZERO, ZERO, Golden(-15, -24), ZERO],
        [Golden(54, 90), ZERO, ZERO, ONE, ZERO, ZERO, ZERO, ZERO, Golden(-12, -18), ZERO],
        [Golden(36, 60), ZERO, ONE, ZERO, ZERO, ZERO, ZERO, ZERO, Golden(-8, -12), ZERO],
        [Golden(18, 30), ONE, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, Golden(-4, -6), ZERO],
        [Golden(48, 80), ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, Golden(-9, -16), ZERO],
        [ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ZERO, ONE],
    ]
    old_for_new = (0, 9, 1, 2, 3, 4, 5, 6, 7, 8)
    return [
        [old_order[old_for_new[i]][old_for_new[j]] for j in range(10)]
        for i in range(10)
    ]


def pell_fold_on_fixed_lattice() -> list[list[Golden]]:
    """Return the old complete-chamber Pell fold in ``(x,u,E8)`` order."""

    result = [[ONE if i == j else ZERO for j in range(10)] for i in range(10)]
    result[0][0] = PHI
    result[0][1] = -ONE
    result[1][0] = PHI
    result[1][1] = -PHI
    return result


def codimension_one_lifts(
    fixed_isometry: Sequence[Sequence[Golden]],
) -> tuple[list[list[Golden]], list[list[Golden]]]:
    """Return the only two field-linear isometric extensions to ``M20``.

    The fixed hyperplane is ``y=0`` and its orthogonal line is generated by
    ``nu=(1,2phi,-phi,0^8)``.  An extension must send ``nu`` to ``+nu`` or
    ``-nu``.  Since

        e_y = (-c,0,1/2,0^8) + c*nu,  c=(phi-1)/2,

    these two choices determine the full matrix.
    """

    fixed_indices = (0, 2, *range(3, 11))
    normal = (ONE, 2 * PHI, -PHI, *(ZERO for _ in range(8)))
    coefficient = (PHI - ONE) / 2
    fixed_part = (-coefficient, PHI * coefficient, *(ZERO for _ in range(8)))
    lifts: list[list[list[Golden]]] = []
    for normal_sign in (ONE, -ONE):
        lift = [[ZERO for _ in range(11)] for _ in range(11)]
        for source, full_source in enumerate(fixed_indices):
            for target, full_target in enumerate(fixed_indices):
                lift[full_target][full_source] = fixed_isometry[target][source]
        fixed_image = _golden_action(fixed_isometry, fixed_part)
        image_y = [ZERO for _ in range(11)]
        for target, full_target in enumerate(fixed_indices):
            image_y[full_target] = fixed_image[target]
        for target in range(11):
            image_y[target] += normal_sign * coefficient * normal[target]
            lift[target][1] = image_y[target]
        lifts.append(lift)
    return lifts[0], lifts[1]


def _nonintegral_entries(
    matrix: Sequence[Sequence[Golden]],
) -> list[dict[str, object]]:
    return [
        {"row": i, "column": j, "value": str(entry)}
        for i, row in enumerate(matrix)
        for j, entry in enumerate(row)
        if not entry.is_integral
    ]


def projected_wall_y_values() -> tuple[Golden, ...]:
    """Enumerate all possible middle coordinates for a wall meeting H9.

    If the fixed projection of a norm-two root has nonnegative square, then
    the identity and conjugate bounds are

        y^2 <= 4/(2+phi),  (y')^2 <= 4/(3-phi).

    They imply ``|b|<=1`` and ``|a|<=3`` for ``y=a+b*phi``.  The finite exact
    loop below is therefore exhaustive, rather than a heuristic cutoff.
    """

    identity_bound = Golden(4) / (Golden(2) + PHI)
    conjugate_bound = Golden(4) / (Golden(3) - PHI)
    return tuple(
        y
        for a in range(-3, 4)
        for b in range(-1, 2)
        for y in (Golden(a, b),)
        if (identity_bound - y * y).sign() >= 0
        and (conjugate_bound - y.conjugate() * y.conjugate()).sign() >= 0
    )


def _class_order(vector: Matrix, search_bound: int) -> int:
    for order in range(1, search_bound + 1):
        if all((order * entry).is_Integer is True for entry in vector):
            return order
    raise ArithmeticError("the discriminant-class order exceeded its bound")


def certify() -> dict[str, object]:
    h = field_gram()
    h_conjugate = [[entry.conjugate() for entry in row] for row in h]
    field_determinant = determinant(h)
    alpha = -PHI * INV_SQRT5
    trace = trace_gram(h, alpha)
    trace_integer = _sympy(integral_gram(trace))
    split = trace_split_change()
    split_trace = split.T * trace_integer * split
    expected_trace = expected_split_trace_gram()

    phi_action = _sympy(multiplication_by_phi(11))
    identity22 = eye(22)

    n_small = ternary_trace_gram()
    k_small = complement_gram()
    direct_sum = Matrix.diag(n_small, k_small)
    glue_two, glue_ten = gluing_vectors()
    glue_matrix = Matrix.hstack(eye(8), glue_two, glue_ten)
    basis = overlattice_basis()
    basis_inverse = basis.inv()
    generator_coordinates = basis_inverse * glue_matrix
    glued = _integer_matrix(basis.T * direct_sum * basis)

    n_embedding = _integer_matrix(
        basis_inverse * Matrix.vstack(eye(6), Matrix.zeros(2, 6))
    )
    k_embedding = _integer_matrix(
        basis_inverse * Matrix.vstack(Matrix.zeros(6, 2), eye(2))
    )
    n_embedding_smith = smith_normal_form(n_embedding, domain=ZZ)
    k_embedding_smith = smith_normal_form(k_embedding, domain=ZZ)

    standard_small = standard_small_og10_gram()
    small_isometry = small_og10_isometry()

    glue_two_square = (glue_two.T * direct_sum * glue_two)[0]
    glue_ten_square = (glue_ten.T * direct_sum * glue_ten)[0]
    glue_pairing = (glue_two.T * direct_sum * glue_ten)[0]

    # In A_K = Z/2 + Z/30, the gluing image is generated by (1,15)
    # and (1,12).  Its orthogonal residual is generated by (0,10),
    # represented by (0,1/3), with q=10/3 = 4/3 mod 2.
    residual_vector = Matrix([0, Rational(1, 3)])
    residual_square = (residual_vector.T * k_small * residual_vector)[0]

    j = natural_involution()
    j_square = _golden_matmul(j, j)
    j_isometry = _golden_matmul(
        _golden_transpose(j), _golden_matmul(h, j)
    )
    root = (ZERO, ONE, ZERO, *(ZERO for _ in range(8)))
    root_image = _golden_action(j, root)
    projection = tuple((left + right) / 2 for left, right in zip(root, root_image))
    normal = (ONE, 2 * PHI, -PHI, *(ZERO for _ in range(8)))
    fixed_indices = (0, 2, *range(3, 11))
    fixed_gram = [[h[i][j_index] for j_index in fixed_indices] for i in fixed_indices]

    bugaenko_reversal = bugaenko_reversal_on_fixed_lattice()
    pell_fold = pell_fold_on_fixed_lattice()
    bugaenko_lifts = codimension_one_lifts(bugaenko_reversal)
    pell_lifts = codimension_one_lifts(pell_fold)
    identity11 = [[ONE if i == k else ZERO for k in range(11)] for i in range(11)]
    identity10 = [[ONE if i == k else ZERO for k in range(10)] for i in range(10)]

    identity_y_bound = Golden(4) / (Golden(2) + PHI)
    conjugate_y_bound = Golden(4) / (Golden(3) - PHI)
    possible_y = projected_wall_y_values()
    expected_y = (-ONE, -ONE + PHI, ZERO, ONE - PHI, ONE)
    projection_square_for_y = lambda value: (
        Golden(2) - (Golden(2) + PHI) * value * value / 2
    )
    scale_for_square_one = 2 * PHI
    scale_for_square_unit = 2 * (PHI - ONE)
    folded_fixed_root = tuple(2 * PHI * entry for entry in projection)

    old_fixed_certificate = json.loads(OLD_FIXED_CHAMBER_CERTIFICATE.read_text())

    checks = {
        "field_determinant": field_determinant == -2 - 6 * PHI,
        "field_determinant_norm_minus_20": field_determinant.norm() == -20,
        "field_identity_signature_10_1": _inertia_tuple(h) == (10, 1, 0),
        "field_conjugate_signature_11_0": _inertia_tuple(h_conjugate) == (11, 0, 0),
        "trace_scalar_codifferent": is_codifferent(alpha),
        "trace_even_integral": is_even_integral(trace),
        "trace_signature_1_21": _inertia_tuple(_integer_matrix_list(trace_integer))
        == (1, 21, 0),
        "trace_determinant_minus_20": int(trace_integer.det()) == -20,
        "trace_discriminant_group_Z2_plus_Z10": _smith_nonunits(trace_integer)
        == [2, 10],
        "trace_split_change_unimodular": abs(int(split.det())) == 1,
        "trace_splits_as_ternary_plus_two_negative_E8": split_trace
        == expected_trace,
        "golden_action_polynomial": phi_action * phi_action
        == phi_action + identity22,
        "golden_action_self_adjoint": phi_action.T * trace_integer
        == trace_integer * phi_action,
        "golden_action_integral": all(
            entry.is_Integer is True for entry in phi_action
        ),
        "gluing_vectors_are_dual": all(
            entry.is_Integer is True
            for vector in (glue_two, glue_ten)
            for entry in direct_sum * vector
        ),
        "gluing_generator_orders_2_10": _class_order(glue_two, 2) == 2
        and _class_order(glue_ten, 10) == 10,
        "gluing_graph_isotropic": glue_two_square == 2
        and glue_ten_square == 2
        and glue_pairing == 2,
        "overlattice_basis_index_20": abs(basis.det()) == Rational(1, 20),
        "overlattice_basis_contains_direct_sum_and_glue": all(
            entry.is_Integer is True for entry in generator_coordinates
        ),
        "glued_lattice_even_integral": all(
            glued[i, i] % 2 == 0 for i in range(glued.rows)
        ),
        "glued_lattice_signature_3_5": _inertia_tuple(
            _integer_matrix_list(glued)
        )
        == (3, 5, 0),
        "glued_lattice_determinant_minus_3": int(glued.det()) == -3,
        "glued_discriminant_group_Z3": _smith_nonunits(glued) == [3],
        "residual_discriminant_square_10_over_3": residual_square
        == Rational(10, 3),
        "explicit_small_og10_isometry": abs(int(small_isometry.det())) == 1
        and small_isometry.T * glued * small_isometry == standard_small,
        "n_embedding_has_correct_gram": n_embedding.T * glued * n_embedding
        == n_small,
        "k_embedding_has_correct_gram": k_embedding.T * glued * k_embedding
        == k_small,
        "n_and_k_embeddings_are_orthogonal": n_embedding.T * glued * k_embedding
        == Matrix.zeros(6, 2),
        "primitive_n_embedding": all(
            abs(int(n_embedding_smith[i, i])) == 1 for i in range(6)
        ),
        "primitive_k_embedding": all(
            abs(int(k_embedding_smith[i, i])) == 1 for i in range(2)
        ),
        "trace_mod_3_nondegenerate": int(trace_integer.det()) % 3 != 0,
        "primitive_ambient_divisibility_three_picard_classes_excluded": (
            int(trace_integer.det()) % 3 != 0
            and all(abs(int(n_embedding_smith[i, i])) == 1 for i in range(6))
        ),
        "natural_involution_isometry": j_isometry == h,
        "natural_involution_order_two": j_square
        == identity11,
        "fixed_lattice_is_L9": fixed_gram
        == _field_block_l9(),
        "involution_normal_is_anti_invariant": _golden_action(j, normal)
        == tuple(-entry for entry in normal),
        "involution_normal_norm_6_plus_8phi": _golden_pair(normal, h, normal)
        == 6 + 8 * PHI,
        "root_and_image_are_norm_two": _golden_pair(root, h, root) == 2
        and _golden_pair(root_image, h, root_image) == 2,
        "root_image_coordinates": root_image
        == (ONE - PHI, -ONE, ONE, *(ZERO for _ in range(8))),
        "root_pair_is_a_5_edge": _golden_pair(root, h, root_image) == -PHI,
        "fixed_projection_is_spacelike": _golden_pair(projection, h, projection)
        == (2 - PHI) / 2,
        "m5_fold_is_existing_fixed_norm_two_root": (
            folded_fixed_root
            == (-ONE, ZERO, PHI, *(ZERO for _ in range(8)))
            and _golden_pair(folded_fixed_root, h, folded_fixed_root) == 2
            and _golden_action(j, folded_fixed_root) == folded_fixed_root
        ),
        "old_bugaenko_subspace_has_codimension_two": len(fixed_indices) == 10
        and len((0, *range(3, 11))) == 9,
        "last_root_mirror_is_tilted": h[1][2] == ONE and h[2][2] == 2,
        "projected_wall_bounds_make_enumeration_finite": (
            (Golden(4) - identity_y_bound).sign() > 0
            and (Golden(4) - conjugate_y_bound).sign() > 0
            and (Golden(2) - PHI).sign() > 0
        ),
        "projected_wall_square_formula_coefficient": (
            (PHI - ONE)
            * (PHI - ONE)
            * (Golden(6) + 8 * PHI)
            / 4
            == (Golden(2) + PHI) / 2
        ),
        "conjugate_completion_of_squares_coefficient": (
            Golden(3) / 2 - ONE / (2 * (PHI - ONE))
            == (Golden(3) - PHI) / 2
        ),
        "projected_wall_middle_coordinate_list_is_complete": possible_y
        == expected_y,
        "projected_wall_y_squares_are_exactly_two_scaling_cases": all(
            value == ZERO
            or value * value == ONE
            or value * value == Golden(2) - PHI
            for value in possible_y
        ),
        "projected_wall_cases_have_positive_projection": all(
            projection_square_for_y(value).sign() > 0 for value in possible_y
        ),
        "projected_wall_rescalings_have_norm_two": (
            scale_for_square_one
            * scale_for_square_one
            * projection_square_for_y(ONE)
            == 2
            and scale_for_square_unit
            * scale_for_square_unit
            * projection_square_for_y(PHI - ONE)
            == 2
        ),
        "projected_wall_rescalings_are_integral_on_all_coordinates": all(
            value.is_integral
            for scale in (scale_for_square_one, scale_for_square_unit)
            for value in (scale, scale / 2, scale * (ONE - PHI) / 2)
        ),
        "fixed_bugaenko_reversal_is_isometric_involution": (
            _golden_matmul(
                _golden_transpose(bugaenko_reversal),
                _golden_matmul(fixed_gram, bugaenko_reversal),
            )
            == fixed_gram
            and _golden_matmul(bugaenko_reversal, bugaenko_reversal)
            == identity10
        ),
        "both_bugaenko_reversal_lifts_are_integral_isometric_involutions": all(
            all(entry.is_integral for row in lift for entry in row)
            and _golden_matmul(
                _golden_transpose(lift), _golden_matmul(h, lift)
            )
            == h
            and _golden_matmul(lift, lift) == identity11
            for lift in bugaenko_lifts
        ),
        "fixed_pell_fold_is_isometric_involution": (
            _golden_matmul(
                _golden_transpose(pell_fold),
                _golden_matmul(fixed_gram, pell_fold),
            )
            == fixed_gram
            and _golden_matmul(pell_fold, pell_fold) == identity10
        ),
        "both_pell_fold_lifts_are_isometries_but_nonintegral": all(
            _golden_matmul(_golden_transpose(lift), _golden_matmul(h, lift))
            == h
            and _golden_matmul(lift, lift) == identity11
            and bool(_nonintegral_entries(lift))
            for lift in pell_lifts
        ),
        "archived_complete_fixed_chamber_has_infinite_volume": (
            old_fixed_certificate["verified"] is True
            and old_fixed_certificate["fixed_lattice"]["decomposition"]
            == "<-2*phi> + E8 + <2>"
            and old_fixed_certificate["group_theoretic_conclusion"][
                "P9_hyperbolic_volume"
            ]
            == "infinite"
        ),
    }
    if not all(checks.values()):
        failed = [name for name, passes in checks.items() if not passes]
        raise AssertionError(f"A2-coupled OG10 certificate failed: {failed}")

    return {
        "schema_version": "0.1.0",
        "certificate_kind": "exact_A2_coupled_golden_OG10_carrier_and_infinite_volume_screen",
        "candidate_id": "og10_golden_a2_coupled_discriminant20",
        "verified": True,
        "field": "Q(sqrt(5))",
        "ring": "Z[phi]",
        "trace_scalar": "-phi/sqrt(5)",
        "field_lattice": {
            "gram": _golden_matrix_text(h),
            "determinant": str(field_determinant),
            "determinant_norm": int(field_determinant.norm()),
            "identity_signature": [10, 1],
            "conjugate_signature": [11, 0],
            "description": "[[-2phi,1,0],[1,2,1],[0,1,2]] + E8",
        },
        "trace_lattice": {
            "ternary_trace_gram": _integer_matrix_list(n_small),
            "full_rank": 22,
            "signature": [1, 21],
            "determinant": int(trace_integer.det()),
            "smith_nonunit_invariants": _smith_nonunits(trace_integer),
            "abstract_discriminant_group": "Z/2 + Z/10",
            "split_type": "N_ternary + E8(-1) + E8(-1)",
            "multiplication_by_phi_is_integral_self_adjoint": True,
        },
        "primitive_og10_embedding": {
            "status": "proved_by_explicit_gluing_and_integral_isometry",
            "positive_complement_gram": _integer_matrix_list(k_small),
            "positive_complement_type": "<2> + <30>",
            "gluing_generator_orders": [2, 10],
            "gluing_vectors_in_N_plus_K_coordinates": [
                [_rational_text(entry) for entry in glue_two],
                [_rational_text(entry) for entry in glue_ten],
            ],
            "gluing_squares_and_pairing": [
                _rational_text(glue_two_square),
                _rational_text(glue_ten_square),
                _rational_text(glue_pairing),
            ],
            "overlattice_basis": _matrix_text(basis),
            "glued_small_gram": _integer_matrix_list(glued),
            "small_og10_isometry": _integer_matrix_list(small_isometry),
            "small_standard_type": "3U + A2(-1)",
            "residual_discriminant_group": "Z/3",
            "residual_generator_square_mod_2Z": "4/3",
            "n_embedding_smith_diagonal": [
                abs(int(n_embedding_smith[i, i])) for i in range(6)
            ],
            "k_embedding_smith_diagonal": [
                abs(int(k_embedding_smith[i, i])) for i in range(2)
            ],
            "full_ambient_type": "3U + 2E8(-1) + A2(-1)",
            "conclusion": (
                "the_trace_lattice_embeds_primitively_in_the_actual_OG10_"
                "Beauville_lattice_with_orthogonal_complement_<2>+<30>"
            ),
        },
        "mod_3": {
            "trace_determinant_coprime_to_3": True,
            "primitive_ambient_divisibility_three_picard_classes_excluded": True,
            "argument": (
                "a_primitive_ambient_divisibility_3_class_would_reduce_to_"
                "the_radical_mod_3_and_hence_to_the_radical_of_N_mod_3;_"
                "det(N)=-20_makes_that_restriction_nondegenerate"
            ),
        },
        "algebro_geometric_realization": {
            "status": "standard_period_theorem_consequence",
            "period_complement": "<2> + <30>",
            "conclusion": (
                "standard_OG10_period_surjectivity_and_the_positive_square_"
                "projectivity_criterion_realize_the_primitive_trace_lattice_"
                "as_the_full_Neron_Severi_lattice_of_a_projective_OG10_type_"
                "IHS"
            ),
            "scope": (
                "the_explicit_lattice_embedding_is_certified_here;_the_"
                "period_and_projectivity_theorems_are_external_inputs"
            ),
        },
        "local_fold_diagnostic": {
            "involution": _golden_matrix_text(j),
            "involution_fixed_lattice": "<-2phi> + E8 + <2>",
            "anti_invariant_normal": _golden_vector_text(normal),
            "anti_invariant_normal_norm": str(_golden_pair(normal, h, normal)),
            "root": _golden_vector_text(root),
            "root_image": _golden_vector_text(root_image),
            "root_pairing": str(_golden_pair(root, h, root_image)),
            "root_pair_coxeter_label": 5,
            "fixed_projection": _golden_vector_text(projection),
            "fixed_projection_norm": str(_golden_pair(projection, h, projection)),
            "folded_fixed_norm_two_root": _golden_vector_text(folded_fixed_root),
            "old_bugaenko_subspace_codimension": 2,
            "fixed_mirror_restriction": "x+u=0",
            "last_coordinate_root_functional": "y+2u",
            "certified_conclusion": (
                "the_nonfixed_root_pair_is_an_m5_pair_and_its_common_"
                "restriction_is_the_existing_fixed_norm_two_wall_with_"
                "normal_(-1,0,phi,0^8);_the_m5_pair_does_not_subdivide_the_"
                "fixed_norm_two_arrangement"
            ),
            "not_proved": (
                "the_pair_has_not_been_shown_simple_for_or_relevant_to_a_"
                "selected_complete_chamber_and_does_not_cut_the_Bugaenko_"
                "H8_relative_interior"
            ),
        },
        "projected_wall_theorem": {
            "root_coordinates": "r=(x,y,u,z)",
            "fixed_projection": (
                "p=(r+Jr)/2=(x+(1-phi)y/2,0,u+y/2,z)"
            ),
            "projection_square": "p^2=2-(2+phi)y^2/2",
            "identity_bound": str(identity_y_bound),
            "conjugate_bound": str(conjugate_y_bound),
            "conjugate_completion_of_squares": (
                "2=r'^2>=(3-phi)(y')^2/2"
            ),
            "finite_enumeration_argument": (
                "both_square_bounds_are_below_4,_so_5b^2=(y-y')^2<16_"
                "and_|b|<=1;_then_|a|<=|y|+|b|phi<2+phi<4,_so_|a|<=3"
            ),
            "identity_bound_below_4_margin": str(
                Golden(4) - identity_y_bound
            ),
            "conjugate_bound_below_4_margin": str(
                Golden(4) - conjugate_y_bound
            ),
            "possible_middle_coordinates": _golden_vector_text(possible_y),
            "rescaling_for_y_square_1": "q=2phi*p",
            "rescaling_for_y_square_2_minus_phi": "q=2(phi-1)*p",
            "rescaled_square": 2,
            "arrangement_conclusion": (
                "an_M20_norm_two_mirror_meets_the_fixed_H9_if_and_only_if_"
                "its_restriction_is_a_fixed_L9_norm_two_mirror"
            ),
            "chamber_conclusion": (
                "for_a_generic_fixed_point_v_the_unique_ambient_chamber_C_"
                "is_J_stable_and_C_intersect_H9_is_the_complete_L9_"
                "norm_two_chamber_containing_v"
            ),
            "congruence_to_archive": (
                "all_complete_L9_norm_two_chambers_are_W2(L9)_congruent_to_"
                "the_archived_infinite_volume_P9"
            ),
        },
        "fixed_symmetry_lift_diagnostics": {
            "orthogonal_normal": _golden_vector_text(normal),
            "only_two_extensions": "normal_maps_to_plus_or_minus_itself",
            "bugaenko_reversal_fixed_matrix": _golden_matrix_text(
                bugaenko_reversal
            ),
            "bugaenko_reversal_extension_status": (
                "both_extensions_are_Z[phi]_integral_isometric_involutions"
            ),
            "bugaenko_reversal_extension_y_columns": [
                _golden_vector_text([row[1] for row in lift])
                for lift in bugaenko_lifts
            ],
            "pell_fold_fixed_matrix": _golden_matrix_text(pell_fold),
            "pell_fold_extension_status": (
                "both_field_isometric_involutions_are_nonintegral"
            ),
            "pell_fold_extension_y_columns": [
                _golden_vector_text([row[1] for row in lift])
                for lift in pell_lifts
            ],
            "pell_fold_nonintegral_entries": [
                _nonintegral_entries(lift) for lift in pell_lifts
            ],
            "scope": (
                "these_lift_facts_are_secondary_diagnostics;_the_volume_"
                "decision_uses_the_complete_projected_wall_theorem_and_the_"
                "intrinsic_archived_L9_infinite_volume_theorem"
            ),
        },
        "norm_two_volume_decision": {
            "ambient_hyperbolic_dimension": 10,
            "fixed_hyperbolic_dimension": 9,
            "ambient_norm_two_chamber_volume": "infinite",
            "argument": (
                "if_the_ambient_norm_two_chamber_had_finite_volume_then_its_"
                "reflection_group_would_have_finite_index_in_the_uniform_"
                "arithmetic_group_O(M20);_the_chamber_would_be_compact,_"
                "forcing_its_complete_fixed_L9_slice_to_be_compact,_"
                "contrary_to_the_archived_infinite_volume_P9_theorem"
            ),
            "uniformity": (
                "the_conjugate_form_is_positive_definite,_so_M20_is_"
                "anisotropic_over_Q(sqrt5)_and_its_arithmetic_group_is_"
                "cocompact"
            ),
            "scope": (
                "this_decides_the_norm_two_Coxeter_chamber;_a_complete_OG10_"
                "mixed_Kahler_wall_shell_audit_for_this_new_trace_lattice_"
                "has_not_been_run"
            ),
        },
        "checks": checks,
        "proved_scope": [
            "exact_field_signatures_and_even_trace_lattice",
            "primitive_embedding_in_the_OG10_Beauville_lattice",
            "positive_rank_two_complement_<2>+<30>",
            "absence_of_primitive_ambient_divisibility_three_Picard_classes",
            "complete_classification_of_norm_two_mirrors_meeting_the_fixed_H9",
            "equality_of_the_projected_and_fixed_L9_norm_two_arrangements",
            "infinite_volume_of_every_M20_norm_two_chamber",
        ],
        "open_gates": [
            "period_and_projectivity_realization_uses_the_standard_OG10_period_theorem",
            "complete_OG10_Kahler_and_movable_wall_shell_audit_has_not_been_ported",
            "no_norm_two_Vinberg_prefix_has_been_enumerated_for_this_form",
            "finite_sidedness_of_the_infinite_volume_norm_two_chamber_is_open",
            "additional_nonreflective_OG10_Kahler_walls_could_refine_the_norm_two_chamber",
        ],
    }


def _field_block_l9() -> list[list[Golden]]:
    """The fixed form in the order (time, extra line, E8)."""

    result = [[ZERO for _ in range(10)] for _ in range(10)]
    result[0][0] = -2 * PHI
    result[1][1] = Golden(2)
    e8 = e8_cartan_gram()
    for i in range(8):
        for j in range(8):
            result[2 + i][2 + j] = Golden(e8[i][j])
    return result


def main() -> None:
    print(json.dumps(certify(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
