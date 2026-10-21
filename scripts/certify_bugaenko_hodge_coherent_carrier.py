"""Exact certificate for the Hodge-coherent Bugaenko K3 carrier.

The existing Bugaenko certificate acts by ``phi`` only on the rank-18
Neron--Severi lattice.  This module verifies a rational self-adjoint extension
to its rank-four K3 orthogonal complement

    T = <-2> + <2> + U.

The square-root-of-five operator ``B`` has ``B^2=5`` and ``T*B`` positive
definite.  Hence its ``+sqrt(5)`` eigenspace is the positive period plane and
``tau=(1+B)/2`` is a Hodge action satisfying ``tau^2=tau+1``.

The Hodge-theoretic and algebraicity implications are proved in the companion
report; this executable checks only the exact lattice identities on which
they rest.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from typing import Iterable, Sequence


Matrix = tuple[tuple[Fraction, ...], ...]
Quadratic = tuple[Fraction, Fraction]


def _matrix(rows: Sequence[Sequence[int | Fraction]]) -> Matrix:
    return tuple(tuple(Fraction(entry) for entry in row) for row in rows)


def _identity(size: int) -> Matrix:
    return _matrix([[int(i == j) for j in range(size)] for i in range(size)])


def _transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[j][i] for j in range(len(matrix))) for i in range(len(matrix)))


def _matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum((left[i][k] * right[k][j] for k in range(len(right))), Fraction(0))
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def _add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[i][j] + right[i][j] for j in range(len(left[0])))
        for i in range(len(left))
    )


def _scale(scalar: Fraction, matrix: Matrix) -> Matrix:
    return tuple(tuple(scalar * entry for entry in row) for row in matrix)


def _inverse(matrix: Matrix) -> Matrix:
    size = len(matrix)
    work = [list(row) + list(_identity(size)[i]) for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            raise ValueError("matrix is singular")
        work[column], work[pivot] = work[pivot], work[column]
        pivot_value = work[column][column]
        work[column] = [entry / pivot_value for entry in work[column]]
        for row in range(size):
            if row == column:
                continue
            multiplier = work[row][column]
            work[row] = [
                work[row][index] - multiplier * work[column][index]
                for index in range(2 * size)
            ]
    return tuple(tuple(row[size:]) for row in work)


def _block_diagonal(*blocks: Matrix) -> Matrix:
    size = sum(len(block) for block in blocks)
    result = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, entry in enumerate(row):
                result[offset + i][offset + j] = entry
        offset += len(block)
    return tuple(tuple(row) for row in result)


def _is_integral(matrix: Matrix) -> bool:
    return all(entry.denominator == 1 for row in matrix for entry in row)


def _quadratic(a: int | Fraction = 0, b: int | Fraction = 0) -> Quadratic:
    """Encode ``a+b*sqrt(5)`` exactly."""

    return Fraction(a), Fraction(b)


def _qadd(left: Quadratic, right: Quadratic) -> Quadratic:
    return left[0] + right[0], left[1] + right[1]


def _qscale(scalar: Quadratic, value: Quadratic) -> Quadratic:
    return (
        scalar[0] * value[0] + 5 * scalar[1] * value[1],
        scalar[0] * value[1] + scalar[1] * value[0],
    )


def _qmatvec(matrix: Matrix, vector: Sequence[Quadratic]) -> tuple[Quadratic, ...]:
    return tuple(
        (
            sum((matrix[i][j] * vector[j][0] for j in range(len(vector))), Fraction(0)),
            sum((matrix[i][j] * vector[j][1] for j in range(len(vector))), Fraction(0)),
        )
        for i in range(len(matrix))
    )


def _qvector_scale(
    scalar: Quadratic, vector: Sequence[Quadratic]
) -> tuple[Quadratic, ...]:
    return tuple(_qscale(scalar, entry) for entry in vector)


def _qpair(
    left: Sequence[Quadratic], gram: Matrix, right: Sequence[Quadratic]
) -> Quadratic:
    result = _quadratic()
    for i in range(len(left)):
        for j in range(len(right)):
            result = _qadd(
                result,
                _qscale(left[i], _qscale(_quadratic(gram[i][j]), right[j])),
            )
    return result


def _qencode(value: Quadratic) -> list[int | str]:
    return [
        entry.numerator if entry.denominator == 1 else str(entry) for entry in value
    ]


def _qvector_encode(vector: Sequence[Quadratic]) -> list[list[int | str]]:
    return [_qencode(entry) for entry in vector]


def _determinant(matrix: Matrix) -> Fraction:
    size = len(matrix)
    if size == 0:
        return Fraction(1)
    work = [list(row) for row in matrix]
    determinant = Fraction(1)
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant *= pivot_value
        for row in range(column + 1, size):
            multiplier = work[row][column] / pivot_value
            for index in range(column + 1, size):
                work[row][index] -= multiplier * work[column][index]
    return determinant


def _leading_minors(matrix: Matrix) -> tuple[Fraction, ...]:
    return tuple(
        _determinant(tuple(tuple(row[j] for j in range(size)) for row in matrix[:size]))
        for size in range(1, len(matrix) + 1)
    )


def _encode(matrix: Matrix) -> list[list[int | str]]:
    return [
        [entry.numerator if entry.denominator == 1 else str(entry) for entry in row]
        for row in matrix
    ]


def _sparse_encode(matrix: Matrix) -> list[list[int | str]]:
    return [
        [i, j, entry.numerator if entry.denominator == 1 else str(entry)]
        for i, row in enumerate(matrix)
        for j, entry in enumerate(row)
        if entry
    ]


T_GRAM = _matrix(
    [
        [-2, 0, 0, 0],
        [0, 2, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 1, 0],
    ]
)


SQRT5_ACTION = _matrix(
    [
        [-3, 2, 0, 0],
        [-2, 3, 0, 0],
        [0, 0, 0, 5],
        [0, 0, 1, 0],
    ]
)


# This commuting skew-adjoint operator makes the period construction
# explicitly CM rather than merely relying on a dimension count.
CM_ACTION = _matrix(
    [
        [0, 0, 1, -5],
        [0, 0, -1, -5],
        [-10, 10, 0, 0],
        [2, 2, 0, 0],
    ]
)


def _trace_to_standard_basis() -> Matrix:
    """Columns are the standard ``<2>+<-2>+E8(-1)^2`` basis in trace coordinates."""

    change = [[0 for _ in range(18)] for _ in range(18)]
    change[0][0] = 1
    change[0][1] = -2
    change[1][1] = 1
    for i in range(8):
        old_e = 2 + 2 * i
        old_phi_e = old_e + 1
        change[old_e][2 + i] = 1
        change[old_e][10 + i] = -1
        change[old_phi_e][10 + i] = 1
    return _matrix(change)


def _ns_standard_sqrt5_action() -> Matrix:
    """Return ``2*phi-1`` on ``N`` in its standard rank-18 basis."""

    phi = [[0 for _ in range(18)] for _ in range(18)]
    for index in range(9):
        e = 2 * index
        phi_e = e + 1
        phi[phi_e][e] = 1
        phi[e][phi_e] = 1
        phi[phi_e][phi_e] = 1
    identity = _identity(18)
    sqrt5_trace = _add(_scale(Fraction(2), _matrix(phi)), _scale(Fraction(-1), identity))
    change = _trace_to_standard_basis()
    return _matmul(_matmul(_inverse(change), sqrt5_trace), change)


def _combined_k3_basis() -> Matrix:
    """Columns span ``N+T`` inside ``U^3+E8(-1)^2``; the index is four."""

    combined = [[0 for _ in range(22)] for _ in range(22)]
    # N=<2>+<-2>+E8(-1)^2.
    combined[0][0] = combined[1][0] = 1
    combined[2][1] = 1
    combined[3][1] = -1
    for index in range(16):
        combined[6 + index][2 + index] = 1
    # T=<-2>+<2>+U.
    combined[0][18] = 1
    combined[1][18] = -1
    combined[2][19] = combined[3][19] = 1
    combined[4][20] = 1
    combined[5][21] = 1
    return _matrix(combined)


def _k3_gram() -> Matrix:
    hyperbolic_plane = _matrix([[0, 1], [1, 0]])
    e8 = [[2 if i == j else 0 for j in range(8)] for i in range(8)]
    for i, j in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)):
        e8[i][j] = e8[j][i] = -1
    e8_negative = _matrix([[-entry for entry in row] for row in e8])
    return _block_diagonal(
        hyperbolic_plane,
        hyperbolic_plane,
        hyperbolic_plane,
        e8_negative,
        e8_negative,
    )


def _ns_standard_gram() -> Matrix:
    e8 = [[2 if i == j else 0 for j in range(8)] for i in range(8)]
    for i, j in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)):
        e8[i][j] = e8[j][i] = -1
    e8_negative = _matrix([[-entry for entry in row] for row in e8])
    return _block_diagonal(_matrix([[2]]), _matrix([[-2]]), e8_negative, e8_negative)


def _discriminant_action_is_identity(action: Matrix, generator_indices: tuple[int, ...]) -> bool:
    """Check ``action(e_i/2)-e_i/2`` is integral for the stated generators."""

    for index in generator_indices:
        generator = _matrix(
            [[Fraction(int(row == index), 2)] for row in range(len(action))]
        )
        difference = _add(
            _matmul(action, generator),
            _scale(Fraction(-1), generator),
        )
        if not _is_integral(difference):
            return False
    return True


def _rank(vectors: Sequence[Sequence[Fraction]]) -> int:
    work = [list(vector) for vector in vectors]
    if not work:
        return 0
    row = 0
    for column in range(len(work[0])):
        pivot = next((index for index in range(row, len(work)) if work[index][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        pivot_value = work[row][column]
        work[row] = [entry / pivot_value for entry in work[row]]
        for index in range(len(work)):
            if index == row:
                continue
            multiplier = work[index][column]
            work[index] = [
                work[index][j] - multiplier * work[row][j]
                for j in range(len(work[index]))
            ]
        row += 1
        if row == len(work):
            break
    return row


def analyze() -> dict[str, object]:
    identity = _identity(4)
    tau = _scale(Fraction(1, 2), _add(identity, SQRT5_ACTION))
    twisted = _matmul(T_GRAM, SQRT5_ACTION)
    minors = _leading_minors(twisted)
    sqrt5 = _quadratic(0, 1)
    period_u = (
        _quadratic(2),
        _quadratic(3, 1),
        _quadratic(),
        _quadratic(),
    )
    period_v = (
        _quadratic(),
        _quadratic(),
        _quadratic(-5, -1),
        _quadratic(-1, -1),
    )
    period_u_square = _qpair(period_u, T_GRAM, period_u)
    period_v_square = _qpair(period_v, T_GRAM, period_v)
    period_cross = _qpair(period_u, T_GRAM, period_v)
    period_omega_square = (
        _qadd(period_u_square, _qscale(_quadratic(-1), period_v_square)),
        _qscale(_quadratic(2), period_cross),
    )
    period_omega_bar = _qadd(period_u_square, period_v_square)
    trace_sqrt5 = sum(
        (SQRT5_ACTION[index][index] for index in range(4)), Fraction(0)
    )
    multiplicity = 2 if trace_sqrt5 == 0 else 0
    ns_sqrt5 = _ns_standard_sqrt5_action()
    combined_basis = _combined_k3_basis()
    combined_gram = _matmul(
        _matmul(_transpose(combined_basis), _k3_gram()), combined_basis
    )
    expected_combined_gram = _block_diagonal(_ns_standard_gram(), T_GRAM)
    combined_action = _block_diagonal(ns_sqrt5, SQRT5_ACTION)
    full_sqrt5 = _matmul(
        _matmul(combined_basis, combined_action),
        _inverse(combined_basis),
    )
    full_identity = _identity(22)
    full_tau = _scale(Fraction(1, 2), _add(full_identity, full_sqrt5))
    k3_gram = _k3_gram()
    ns_discriminant_identity = _discriminant_action_is_identity(ns_sqrt5, (0, 1))
    t_discriminant_identity = _discriminant_action_is_identity(SQRT5_ACTION, (0, 1))
    cm_basis = (
        identity,
        SQRT5_ACTION,
        CM_ACTION,
        _matmul(SQRT5_ACTION, CM_ACTION),
    )
    cm_i = _scale(
        Fraction(1, 10), _matmul(SQRT5_ACTION, CM_ACTION)
    )
    cayley_t = _add(identity, SQRT5_ACTION)
    cayley_z = _matmul(
        _add(cayley_t, _scale(Fraction(-1), cm_i)),
        _inverse(_add(cayley_t, cm_i)),
    )
    cayley_z_inverse = _inverse(cayley_z)
    buskin_span = _add(
        _scale(Fraction(29, 8), _add(cayley_z, cayley_z_inverse)),
        _scale(Fraction(-15, 4), identity),
    )
    payload: dict[str, object] = {
        "schema_version": "0.2.0",
        "object": "Hodge-coherent rational extension of the Bugaenko K3 divisor action",
        "hyperbolic_dimension": 8,
        "ns_rank": 18,
        "transcendental_rank": 4,
        "transcendental_lattice": "<-2>+<2>+U",
        "transcendental_gram": _encode(T_GRAM),
        "sqrt5_action": _encode(SQRT5_ACTION),
        "tau_action": _encode(tau),
        "sqrt5_square_is_5I": _matmul(SQRT5_ACTION, SQRT5_ACTION) == _scale(Fraction(5), identity),
        "sqrt5_self_adjoint": _matmul(_transpose(SQRT5_ACTION), T_GRAM)
        == _matmul(T_GRAM, SQRT5_ACTION),
        "tau_polynomial_is_tau_plus_I": _matmul(tau, tau) == _add(tau, identity),
        "tau_self_adjoint": _matmul(_transpose(tau), T_GRAM) == _matmul(T_GRAM, tau),
        "tau_is_rational_not_integral_on_T": any(
            entry.denominator != 1 for row in tau for entry in row
        ),
        "sqrt5_twisted_gram": _encode(twisted),
        "sqrt5_twisted_leading_minors": [int(value) for value in minors],
        "sqrt5_twisted_gram_positive_definite": all(value > 0 for value in minors),
        "sqrt5_trace_zero": trace_sqrt5 == 0,
        "eigenvalue_multiplicities": {
            "+sqrt5": multiplicity,
            "-sqrt5": multiplicity,
        },
        "eigenspace_signatures": {
            "+sqrt5": [multiplicity, 0],
            "-sqrt5": [0, multiplicity],
        },
        "cm_action": _encode(CM_ACTION),
        "cm_action_commutes_with_sqrt5": _matmul(CM_ACTION, SQRT5_ACTION)
        == _matmul(SQRT5_ACTION, CM_ACTION),
        "cm_action_is_skew_adjoint": _matmul(_transpose(CM_ACTION), T_GRAM)
        == _scale(Fraction(-1), _matmul(T_GRAM, CM_ACTION)),
        "cm_action_square_is_minus_20I": _matmul(CM_ACTION, CM_ACTION)
        == _scale(Fraction(-20), identity),
        "cm_i_action_equals_BJ_over_10": _encode(cm_i),
        "cm_i_action_square_is_minus_I": _matmul(cm_i, cm_i)
        == _scale(Fraction(-1), identity),
        "cm_field_basis_rank": _rank(
            [tuple(entry for row in matrix for entry in row) for matrix in cm_basis]
        ),
        "cm_field": "Q(sqrt(5),sqrt(-20))=Q(sqrt(5),i)",
        "buskin_cayley_z_action": _encode(cayley_z),
        "buskin_cayley_z_is_T_isometry": _matmul(
            _matmul(_transpose(cayley_z), T_GRAM), cayley_z
        )
        == T_GRAM,
        "buskin_cayley_z_inverse_is_T_isometry": _matmul(
            _matmul(_transpose(cayley_z_inverse), T_GRAM), cayley_z_inverse
        )
        == T_GRAM,
        "buskin_sqrt5_span_identity": buskin_span == SQRT5_ACTION,
        "period_real_part_u": _qvector_encode(period_u),
        "period_imaginary_part_v": _qvector_encode(period_v),
        "period_Bu_is_sqrt5_u": _qmatvec(SQRT5_ACTION, period_u)
        == _qvector_scale(sqrt5, period_u),
        "period_Bv_is_sqrt5_v": _qmatvec(SQRT5_ACTION, period_v)
        == _qvector_scale(sqrt5, period_v),
        "period_Ju_is_minus_sqrt20_v": _qmatvec(CM_ACTION, period_u)
        == _qvector_scale(_quadratic(0, -2), period_v),
        "period_Jv_is_sqrt20_u": _qmatvec(CM_ACTION, period_v)
        == _qvector_scale(_quadratic(0, 2), period_u),
        "period_u_square": _qencode(period_u_square),
        "period_v_square": _qencode(period_v_square),
        "period_u_dot_v": _qencode(period_cross),
        "period_omega_square": {
            "real": _qencode(period_omega_square[0]),
            "imaginary": _qencode(period_omega_square[1]),
        },
        "period_omega_square_is_zero": period_omega_square
        == (_quadratic(), _quadratic()),
        "period_omega_dot_bar_omega": _qencode(period_omega_bar),
        "period_omega_dot_bar_omega_positive": (
            period_omega_bar[0] > 0 and period_omega_bar[1] > 0
        ),
        "integral_full_h2_order": "Z[sqrt(5)]",
        "ns_standard_sqrt5_action_sparse": _sparse_encode(ns_sqrt5),
        "sqrt5_discriminant_action_on_N_is_identity": ns_discriminant_identity,
        "sqrt5_discriminant_action_on_T_is_identity": t_discriminant_identity,
        "combined_N_plus_T_basis_index_in_k3_lattice": int(
            abs(_determinant(combined_basis))
        ),
        "combined_basis_gram_is_N_orthogonal_sum_T": combined_gram
        == expected_combined_gram,
        "full_k3_sqrt5_action_sparse": _sparse_encode(full_sqrt5),
        "full_k3_sqrt5_action_integral": _is_integral(full_sqrt5),
        "full_k3_sqrt5_square_is_5I": _matmul(full_sqrt5, full_sqrt5)
        == _scale(Fraction(5), full_identity),
        "full_k3_sqrt5_self_adjoint": _matmul(_transpose(full_sqrt5), k3_gram)
        == _matmul(k3_gram, full_sqrt5),
        "sqrt5_extends_over_unimodular_glue": (
            ns_discriminant_identity
            and t_discriminant_identity
            and _is_integral(full_sqrt5)
        ),
        "full_k3_tau_action_integral": _is_integral(full_tau),
        "full_k3_tau_polynomial_is_tau_plus_I": _matmul(full_tau, full_tau)
        == _add(full_tau, full_identity),
        "full_k3_tau_self_adjoint": _matmul(_transpose(full_tau), k3_gram)
        == _matmul(k3_gram, full_tau),
        "maximal_order_O_K_extends_integrally_on_full_h2": _is_integral(full_tau),
        "rank20_k3_h9_hodge_coherent_gate": "empty",
        "rank20_reason": "rank(T)=2 would be one-dimensional over a real quadratic field, impossible for the K3 period plane",
        "k3_hilbert_scheme_full_h2_gate": "empty_in_every_dimension",
        "k3_hilbert_scheme_reason": "b2=23 is not divisible by the quadratic field degree",
        "scope_warning": "Does not exclude NS-only eigencarriers or actions confined to a proper primitive-cohomology summand.",
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    payload["payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    return payload


def _all_exact_checks(payload: dict[str, object]) -> Iterable[bool]:
    yield payload["sqrt5_square_is_5I"] is True
    yield payload["sqrt5_self_adjoint"] is True
    yield payload["tau_polynomial_is_tau_plus_I"] is True
    yield payload["tau_self_adjoint"] is True
    yield payload["tau_is_rational_not_integral_on_T"] is True
    yield payload["sqrt5_twisted_leading_minors"] == [6, 20, 20, 100]
    yield payload["sqrt5_twisted_gram_positive_definite"] is True
    yield payload["sqrt5_trace_zero"] is True
    yield payload["eigenvalue_multiplicities"] == {"+sqrt5": 2, "-sqrt5": 2}
    yield payload["cm_action_commutes_with_sqrt5"] is True
    yield payload["cm_action_is_skew_adjoint"] is True
    yield payload["cm_action_square_is_minus_20I"] is True
    yield payload["cm_i_action_square_is_minus_I"] is True
    yield payload["cm_field_basis_rank"] == 4
    yield payload["buskin_cayley_z_is_T_isometry"] is True
    yield payload["buskin_cayley_z_inverse_is_T_isometry"] is True
    yield payload["buskin_sqrt5_span_identity"] is True
    yield payload["period_Bu_is_sqrt5_u"] is True
    yield payload["period_Bv_is_sqrt5_v"] is True
    yield payload["period_Ju_is_minus_sqrt20_v"] is True
    yield payload["period_Jv_is_sqrt20_u"] is True
    yield payload["period_omega_square_is_zero"] is True
    yield payload["period_omega_dot_bar_omega_positive"] is True
    yield payload["sqrt5_discriminant_action_on_N_is_identity"] is True
    yield payload["sqrt5_discriminant_action_on_T_is_identity"] is True
    yield payload["combined_N_plus_T_basis_index_in_k3_lattice"] == 4
    yield payload["combined_basis_gram_is_N_orthogonal_sum_T"] is True
    yield payload["full_k3_sqrt5_action_integral"] is True
    yield payload["full_k3_sqrt5_square_is_5I"] is True
    yield payload["full_k3_sqrt5_self_adjoint"] is True
    yield payload["sqrt5_extends_over_unimodular_glue"] is True
    yield payload["full_k3_tau_action_integral"] is False
    yield payload["full_k3_tau_polynomial_is_tau_plus_I"] is True
    yield payload["full_k3_tau_self_adjoint"] is True
    yield payload["maximal_order_O_K_extends_integrally_on_full_h2"] is False


if __name__ == "__main__":
    result = analyze()
    if not all(_all_exact_checks(result)):
        raise ArithmeticError("Hodge-coherent Bugaenko certificate failed")
    print(json.dumps(result, indent=2, sort_keys=True))
