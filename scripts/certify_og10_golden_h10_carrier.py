#!/usr/bin/env python3
"""Exact lattice certificate for a golden H^10 carrier of OG10 type.

The certificate starts with the rank-11 ``Z[phi]``-lattice

    ([[ -2 phi, 1 ], [ 1, 2 ]] + E8) + <2>

and takes the trace pairing ``Tr(-phi/sqrt(5) h)``.  Its integral trace
lattice has signature (1, 21).  The non-unimodular part is glued to the
positive lattice <2> + <66>; an explicit integral basis identifies the
result with the Beauville lattice

    3 U + 2 E8(-1) + A2(-1)

of an irreducible holomorphic symplectic manifold of OG10 type.

This is a carrier certificate, not a Coxeter/MBM-wall certificate.  The
period-map realization and the projectivity criterion are standard external
theorems and are stated separately in the accompanying research note.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from typing import Any, Iterable

from sympy import Matrix, Rational, eye
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    change_gram_basis,
    integral_gram,
    is_codifferent,
    is_even_integral,
    trace_gram,
)
from coxeter_polytopes.lattices import block_diagonal, e8_cartan_gram, hyperbolic_plane
from scripts.coupled_discriminant11_candidate import (
    T11_REDUCED,
    field_gram as discriminant11_field_gram,
    multiplication_by_phi,
    trace_basis_change as discriminant11_trace_basis_change,
)


def _negative(matrix: list[list[int]]) -> list[list[int]]:
    return [[-entry for entry in row] for row in matrix]


def _sympy(matrix: Iterable[Iterable[Any]]) -> Matrix:
    return Matrix([list(row) for row in matrix])


def _integer_matrix(matrix: Matrix) -> Matrix:
    if not all(entry.is_Integer is True for entry in matrix):
        raise ArithmeticError("an expected integral matrix is nonintegral")
    return matrix.applyfunc(int)


def _matrix_list(matrix: Matrix) -> list[list[int]]:
    return [
        [int(matrix[i, j]) for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]


def _fraction(value: Any) -> Fraction:
    numerator, denominator = value.as_numer_denom()
    return Fraction(int(numerator), int(denominator))


def _mod_one(value: Any) -> Fraction:
    value = _fraction(value)
    return value - value.numerator // value.denominator


def _mod_two(value: Any) -> Fraction:
    value = _fraction(value)
    return value - 2 * (value // 2)


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _matrix_digest(matrix: Matrix) -> str:
    payload = json.dumps(_matrix_list(_integer_matrix(matrix)), separators=(",", ":"))
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def field_gram() -> list[list[Golden]]:
    """The discriminant-11 field form stabilized by one positive O-line."""

    old = discriminant11_field_gram()
    result = [[Golden(0) for _ in range(11)] for _ in range(11)]
    for i in range(10):
        for j in range(10):
            result[i][j] = old[i][j]
    result[10][10] = Golden(2)
    return result


def trace_basis_change() -> list[list[int]]:
    """Split the trace form as U + T11(-1) + 2E8(-1) + 2A1(-1)."""

    old = discriminant11_trace_basis_change()
    result = [[0 for _ in range(22)] for _ in range(22)]
    for i in range(20):
        for j in range(20):
            result[i][j] = old[i][j]
    # The new O-line has trace block [[-2,-2],[-2,-4]].  In the basis
    # e, phi*e-e it becomes <-2> + <-2>.
    result[20][20] = 1
    result[20][21] = -1
    result[21][21] = 1
    return result


def expected_trace_gram() -> list[list[int]]:
    negative_e8 = _negative(e8_cartan_gram())
    return block_diagonal(
        hyperbolic_plane(),
        _negative([[18, 5], [5, 2]]),
        negative_e8,
        negative_e8,
        [[-2]],
        [[-2]],
    )


def t11_reduction() -> Matrix:
    """Change [[18,5],[5,2]] to the reduced form [[2,1],[1,6]]."""

    result = eye(22)
    result[2, 2] = 0
    result[2, 3] = -1
    result[3, 2] = -1
    result[3, 3] = 2
    return result


def core_small_reordering() -> Matrix:
    """Reorder the trace lattice as (U + 2E8(-1)) + its rank-4 part."""

    source_indices = [0, 1, *range(4, 20), 2, 3, 20, 21]
    result = Matrix.zeros(22)
    for target, source in enumerate(source_indices):
        result[source, target] = 1
    return result


def core_gram() -> Matrix:
    negative_e8 = _negative(e8_cartan_gram())
    return _sympy(block_diagonal(hyperbolic_plane(), negative_e8, negative_e8))


def small_n_gram() -> Matrix:
    return Matrix(
        [
            [-2, -1, 0, 0],
            [-1, -6, 0, 0],
            [0, 0, -2, 0],
            [0, 0, 0, -2],
        ]
    )


def complement_gram() -> Matrix:
    return Matrix.diag(2, 66)


def gluing_vectors() -> tuple[Matrix, Matrix]:
    """Generators of the order-44 isotropic graph in A_N + A_T."""

    first = Matrix(
        [6, -1, Rational(11, 2), Rational(1, 2), Rational(1, 2), Rational(1, 2)]
    )
    second = Matrix(
        [Rational(-6, 11), Rational(1, 11), Rational(-1, 2), 0, 0, Rational(9, 22)]
    )
    return first, second


def overlattice_basis() -> Matrix:
    """Columns are an integral basis of the glued lattice in N+T coordinates."""

    return Matrix(
        [
            [1, 0, 0, 0, 0, Rational(3, 11)],
            [0, 1, 0, 0, 0, Rational(5, 11)],
            [0, 0, 1, 0, 0, Rational(1, 2)],
            [0, 0, 0, 1, Rational(1, 2), 0],
            [0, 0, 0, 0, Rational(1, 2), 0],
            [0, 0, 0, 0, 0, Rational(1, 22)],
        ]
    )


def small_og10_isometry() -> Matrix:
    """Take the glued rank-6 Gram matrix to 2U + A2(-1)."""

    columns = [
        [0, -1, -1, 0, 0, 2],
        [-1, -1, -1, 0, 0, 2],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, -1, 1, 0],
        [0, -1, -1, 0, 0, 1],
        [0, -1, 0, 0, 0, 1],
    ]
    return Matrix.hstack(*(Matrix(column) for column in columns))


def _class_key(vector: Matrix) -> tuple[Fraction, ...]:
    return tuple(_mod_one(entry) for entry in vector)


def _span_keys(first: Matrix, second: Matrix, first_order: int, second_order: int) -> set[tuple[Fraction, ...]]:
    return {
        _class_key(a * first + b * second)
        for a in range(first_order)
        for b in range(second_order)
    }


def _class_order(vector: Matrix, search_bound: int) -> int:
    zero = tuple(Fraction(0) for _ in vector)
    for order in range(1, search_bound + 1):
        if _class_key(order * vector) == zero:
            return order
    raise ArithmeticError("the discriminant-class order exceeded its certified bound")


def certify() -> dict[str, Any]:
    h = field_gram()
    h_conjugate = [[entry.conjugate() for entry in row] for row in h]
    alpha = -PHI * INV_SQRT5
    trace = trace_gram(h, alpha)
    trace_integer = integral_gram(trace)

    split = _sympy(trace_basis_change())
    expected = _sympy(expected_trace_gram())
    split_trace = split.T * _sympy(trace_integer) * split
    reduction = t11_reduction()
    reordering = core_small_reordering()
    full_change = split * reduction * reordering
    reduced_reordered = full_change.T * _sympy(trace_integer) * full_change

    core = core_gram()
    n_small = small_n_gram()
    n_reordered = Matrix.diag(core, n_small)

    action_original = _sympy(multiplication_by_phi(11))
    action = _integer_matrix(full_change.inv() * action_original * full_change)
    identity22 = eye(22)

    t = complement_gram()
    direct_sum = Matrix.diag(n_small, t)
    g1, g2 = gluing_vectors()
    n1, t1 = g1[:4, :], g1[4:, :]
    n2, t2 = g2[:4, :], g2[4:, :]

    n_image_keys = _span_keys(n1, n2, 2, 22)
    t_image_keys = _span_keys(t1, t2, 2, 22)
    total_image_keys = _span_keys(g1, g2, 2, 22)

    q_g1 = (g1.T * direct_sum * g1)[0]
    q_g2 = (g2.T * direct_sum * g2)[0]
    b_g1_g2 = (g1.T * direct_sum * g2)[0]

    basis = overlattice_basis()
    basis_inverse = basis.inv()
    generator_coordinates = basis_inverse * Matrix.hstack(eye(6), g1, g2)
    glued = _integer_matrix(basis.T * direct_sum * basis)

    standard_small = Matrix.diag(
        _sympy(hyperbolic_plane()),
        _sympy(hyperbolic_plane()),
        Matrix([[-2, -1], [-1, -2]]),
    )
    small_isometry = small_og10_isometry()

    n_small_embedding = _integer_matrix(
        basis_inverse * Matrix.vstack(eye(4), Matrix.zeros(2, 4))
    )
    t_embedding = _integer_matrix(
        basis_inverse * Matrix.vstack(Matrix.zeros(4, 2), eye(2))
    )
    ambient = Matrix.diag(core, glued)
    n_embedding = Matrix.zeros(24, 22)
    n_embedding[:18, :18] = eye(18)
    n_embedding[18:, 18:] = n_small_embedding
    complement_embedding = Matrix.zeros(24, 2)
    complement_embedding[18:, :] = t_embedding
    complement_embedding_smith = smith_normal_form(complement_embedding, domain=ZZ)

    full_og10_change = Matrix.diag(eye(18), small_isometry)
    og10_standard_order = Matrix.diag(core, standard_small)

    # The image of the gluing subgroup in A_T has orthogonal residual of
    # order three.  Coordinates (a,b) mean (a/2,b/66) for <2>+<66>.
    residual: list[tuple[int, int]] = []
    residual_q: list[Fraction] = []
    for a in range(2):
        for b in range(66):
            vector = Matrix([Rational(a, 2), Rational(b, 66)])
            if all(
                _fraction((vector.T * t * generator)[0]).denominator == 1
                for generator in (t1, t2)
            ):
                residual.append((a, b))
                residual_q.append(_mod_two((vector.T * t * vector)[0]))

    h_identity_signature = inertia(h)
    h_conjugate_signature = inertia(h_conjugate)
    trace_signature = inertia(trace_integer)
    n_signature = inertia(_matrix_list(n_reordered))
    ambient_signature = inertia(_matrix_list(ambient))

    smith = smith_normal_form(glued, domain=ZZ)
    smith_invariants = [
        abs(int(smith[i, i]))
        for i in range(6)
        if abs(int(smith[i, i])) > 1
    ]

    # If h(r,r)=2, the underlying divisor r and its golden multiple phi*r
    # have the two OG10 wall squares -2 and -4.  Self-adjointness gives
    # (phi*r,x)=phi*(r,x) on the phi-eigenspace, so the two ambient wall
    # hyperplanes restrict to one field-root mirror.
    root_lift_square = (alpha * 2).trace()
    phi_root_lift_square = (alpha * 2 * PHI * PHI).trace()
    root_phi_pairing = (alpha * 2 * PHI).trace()
    positive_class = Matrix([1, 1, *(0 for _ in range(20))])

    checks = {
        "field_identity_signature_10_1": (
            h_identity_signature.positive,
            h_identity_signature.negative,
            h_identity_signature.zero,
        )
        == (10, 1, 0),
        "field_conjugate_signature_11_0": (
            h_conjugate_signature.positive,
            h_conjugate_signature.negative,
            h_conjugate_signature.zero,
        )
        == (11, 0, 0),
        "trace_scalar_codifferent": is_codifferent(alpha),
        "trace_even_integral": is_even_integral(trace),
        "split_trace_lattice": split_trace == expected,
        "trace_basis_change_unimodular": abs(int(full_change.det())) == 1,
        "reduced_reordered_trace_lattice": reduced_reordered == n_reordered,
        "trace_signature_1_21": (
            trace_signature.positive,
            trace_signature.negative,
            trace_signature.zero,
        )
        == (1, 21, 0),
        "hyperbolic_phi_eigenspace_signature_1_10": (
            alpha.sign(0) == -1
            and h_identity_signature.positive == 10
            and h_identity_signature.negative == 1
        ),
        "other_eigenspace_negative_definite": (
            alpha.sign(1) == -1
            and h_conjugate_signature.positive == 11
            and h_conjugate_signature.negative == 0
        ),
        "golden_polynomial": action * action == action + identity22,
        "golden_action_self_adjoint": action.T * n_reordered == n_reordered * action,
        "golden_action_integral": all(entry.is_Integer is True for entry in action),
        "norm2_root_trace_lift_squares_and_pairing": (
            root_lift_square == -2
            and phi_root_lift_square == -4
            and root_phi_pairing == -2
        ),
        "n_discriminant_44": abs(int(n_reordered.det())) == 44,
        "n_signature_1_21": (
            n_signature.positive,
            n_signature.negative,
            n_signature.zero,
        )
        == (1, 21, 0),
        "gluing_vectors_are_dual": all(
            entry.is_Integer is True
            for vector in (g1, g2)
            for entry in direct_sum * vector
        ),
        "gluing_graph_order_44": (
            len(total_image_keys) == len(n_image_keys) == len(t_image_keys) == 44
        ),
        "gluing_generator_exact_orders_2_22": (
            _class_order(g1, 2) == 2 and _class_order(g2, 22) == 22
        ),
        "gluing_graph_isotropic": (
            _mod_two(q_g1) == 0
            and _mod_two(q_g2) == 0
            and _fraction(b_g1_g2).denominator == 1
        ),
        "overlattice_basis_index_44": abs(basis.det()) == Rational(1, 44),
        "overlattice_basis_contains_generators": all(
            entry.is_Integer is True for entry in generator_coordinates
        ),
        "glued_lattice_even": all(int(glued[i, i]) % 2 == 0 for i in range(6)),
        "glued_lattice_determinant_3": int(glued.det()) == 3,
        "glued_discriminant_group_Z3": smith_invariants == [3],
        "explicit_small_og10_isometry": (
            abs(int(small_isometry.det())) == 1
            and small_isometry.T * glued * small_isometry == standard_small
        ),
        "primitive_n_embedding": (
            n_embedding[:22, :] == eye(22)
            and n_embedding.T * ambient * n_embedding == n_reordered
        ),
        "orthogonal_complement_diag_2_66": (
            complement_embedding.T * ambient * complement_embedding == t
            and n_embedding.T * ambient * complement_embedding == Matrix.zeros(22, 2)
        ),
        "orthogonal_complement_primitive": (
            abs(int(complement_embedding_smith[0, 0])) == 1
            and abs(int(complement_embedding_smith[1, 1])) == 1
        ),
        "ambient_is_og10_beauville_lattice": (
            full_og10_change.T * ambient * full_og10_change == og10_standard_order
        ),
        "ambient_signature_3_21": (
            ambient_signature.positive,
            ambient_signature.negative,
            ambient_signature.zero,
        )
        == (3, 21, 0),
        "residual_discriminant_form_A2_minus": (
            residual == [(0, 0), (0, 22), (0, 44)]
            and residual_q == [Fraction(0), Fraction(4, 3), Fraction(4, 3)]
        ),
        "positive_class_for_projectivity": (
            (positive_class.T * n_reordered * positive_class)[0] == 2
        ),
        "period_complement_positive_rank_2": t == Matrix.diag(2, 66),
    }

    return {
        "candidate_id": "og10_golden_h10_carrier",
        "field_lattice_rank_over_Z_phi": 11,
        "field_determinant": str(determinant(h)),
        "field_signatures": {
            "identity": [
                h_identity_signature.positive,
                h_identity_signature.negative,
                h_identity_signature.zero,
            ],
            "conjugate": [
                h_conjugate_signature.positive,
                h_conjugate_signature.negative,
                h_conjugate_signature.zero,
            ],
        },
        "trace_scalar": alpha.to_json(),
        "trace_lattice": {
            "isometry_type": "U + T11_reduced(-1) + 2E8(-1) + 2A1(-1)",
            "T11_reduced": T11_REDUCED,
            "rank": 22,
            "signature": [1, 21],
            "determinant": int(n_reordered.det()),
            "golden_action_sha256": _matrix_digest(action),
            "phi_eigenspace_signature": [1, 10],
            "conjugate_eigenspace_signature": [0, 11],
        },
        "primitive_embedding": {
            "complement_gram": _matrix_list(t),
            "gluing_generators": [
                [str(entry) for entry in g1],
                [str(entry) for entry in g2],
            ],
            "gluing_graph_order": len(total_image_keys),
            "glued_rank6_gram": _matrix_list(glued),
            "glued_rank6_determinant": int(glued.det()),
            "glued_rank6_smith_invariants": smith_invariants,
            "glued_rank6_isometry_type": "2U + A2(-1)",
            "residual_T_discriminant_elements": [list(value) for value in residual],
            "residual_T_quadratics_mod_2": [
                _fraction_text(value) for value in residual_q
            ],
            "ambient_isometry_type": "3U + 2E8(-1) + A2(-1)",
        },
        "geometric_realization": {
            "external_inputs": [
                "Rapagnetta's computation of the OG10 Beauville lattice",
                "surjectivity of the marked hyperkahler period map",
                "the positive-square projectivity criterion",
            ],
            "period_complement": "<2> + <66>",
            "conclusion": (
                "there_exists_a_projective_OG10_type_IHS_with_this_"
                "Neron_Severi_lattice_and_golden_H10_carrier"
            ),
            "coxeter_wall_status": "not_yet_certified",
        },
        "root_wall_lift": {
            "field_root_norm": 2,
            "divisor_squares": [int(root_lift_square), int(phi_root_lift_square)],
            "divisor_pairing": int(root_phi_pairing),
            "wall_types": ["square_-2_divisibility_1", "square_-4_divisibility_1"],
            "restricted_hyperplanes_coincide_by_self_adjointness": True,
            "scope": (
                "every_norm_2_field_root_mirror_occurs_in_the_restricted_OG10_"
                "wall_divisor_arrangement;_no_particular_Kahler_chamber_or_"
                "completeness_statement_is_certified"
            ),
        },
        "checks": checks,
        "verified": all(checks.values()),
    }


def main() -> None:
    print(json.dumps(certify(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
