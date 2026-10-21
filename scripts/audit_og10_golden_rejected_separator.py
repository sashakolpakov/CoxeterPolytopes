#!/usr/bin/env python3
"""First-principles regression for a rejected golden OG10 separator transcript.

The supplied covector had E8 entry ``8+14*phi`` while claiming projected
component ``w=(4+7*phi)*v``.  This script reconstructs ``p=F^-1 ell`` and the
full involution for both the literal and corrected covectors.  It intentionally
creates no candidate certificate: the literal transcript is malformed, and
the corrected integral carrier does not separate the archived ray.
"""

from __future__ import annotations

from fractions import Fraction
import json


Pair = tuple[Fraction, Fraction]
Vector = tuple[Pair, ...]
Matrix = tuple[tuple[Pair, ...], ...]


def pair(a: int | Fraction = 0, b: int | Fraction = 0) -> Pair:
    return Fraction(a), Fraction(b)


ZERO = pair()
ONE = pair(1)
PHI = pair(0, 1)


def add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def neg(value: Pair) -> Pair:
    return -value[0], -value[1]


def sub(left: Pair, right: Pair) -> Pair:
    return add(left, neg(right))


def mul(left: Pair, right: Pair) -> Pair:
    a, b = left
    c, d = right
    return a * c + b * d, a * d + b * c + b * d


def scale(value: Pair, scalar: int | Fraction) -> Pair:
    return value[0] * scalar, value[1] * scalar


def sum_pairs(values: object) -> Pair:
    answer = ZERO
    for value in values:
        answer = add(answer, value)
    return answer


def conjugate(value: Pair) -> Pair:
    return value[0] + value[1], -value[1]


def sign(value: Pair) -> int:
    denominator = value[0].denominator * value[1].denominator
    a = int(value[0] * denominator)
    b = int(value[1] * denominator)
    rational_part = 2 * a + b
    if rational_part == 0:
        return (b > 0) - (b < 0)
    if b == 0 or (rational_part > 0) == (b > 0):
        return (rational_part > 0) - (rational_part < 0)
    comparison = rational_part * rational_part - 5 * b * b
    if comparison > 0:
        return (rational_part > 0) - (rational_part < 0)
    return (b > 0) - (b < 0)


E8_GRAM = (
    (2, -1, 0, 0, 0, 0, 0, 0),
    (-1, 2, -1, 0, 0, 0, 0, 0),
    (0, -1, 2, -1, 0, 0, 0, -1),
    (0, 0, -1, 2, -1, 0, 0, 0),
    (0, 0, 0, -1, 2, -1, 0, 0),
    (0, 0, 0, 0, -1, 2, -1, 0),
    (0, 0, 0, 0, 0, -1, 2, 0),
    (0, 0, -1, 0, 0, 0, 0, 2),
)
ROOT_VECTOR = (-2, -4, -6, -5, -4, -3, -2, -3)

INTERIOR_SCALED = (
    pair(601, 199),
    pair(184),
    pair(364),
    pair(540),
    pair(440),
    pair(336),
    pair(228),
    pair(116),
    pair(272),
    pair(2),
)
RAY_SCALED = (
    pair(13),
    pair(0, 26),
    pair(0, 52),
    pair(0, 78),
    pair(0, 64),
    pair(0, 48),
    pair(0, 32),
    pair(0, 16),
    pair(0, 40),
    pair(0, 1),
)


def rational_inverse(matrix: tuple[tuple[int, ...], ...]) -> tuple[tuple[Fraction, ...], ...]:
    size = len(matrix)
    augmented = [
        [Fraction(matrix[i][j]) for j in range(size)]
        + [Fraction(int(i == j)) for j in range(size)]
        for i in range(size)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size) if augmented[row][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        value = augmented[column][column]
        augmented[column] = [entry / value for entry in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    augmented[row][j] - factor * augmented[column][j]
                    for j in range(2 * size)
                ]
    return tuple(
        tuple(augmented[i][size + j] for j in range(size))
        for i in range(size)
    )


def fixed_gram() -> Matrix:
    matrix = [[ZERO for _ in range(10)] for _ in range(10)]
    matrix[0][0] = pair(0, -2)
    for i in range(8):
        for j in range(8):
            matrix[i + 1][j + 1] = pair(E8_GRAM[i][j])
    matrix[9][9] = pair(2)
    return tuple(tuple(row) for row in matrix)


def solve_projection(ell: Vector) -> Vector:
    inverse = rational_inverse(E8_GRAM)
    # 1/phi=phi-1, so the first inverse entry is -1/(2phi).
    x = scale(mul(ell[0], pair(-1, 1)), Fraction(-1, 2))
    e8 = tuple(
        sum_pairs(
            scale(ell[j + 1], inverse[i][j]) for j in range(8)
        )
        for i in range(8)
    )
    return (x,) + e8 + (scale(ell[9], Fraction(1, 2)),)


def matrix_vector(matrix: Matrix, vector: Vector) -> Vector:
    return tuple(
        sum_pairs(mul(matrix[i][j], vector[j]) for j in range(len(vector)))
        for i in range(len(matrix))
    )


def dot(left: Vector, right: Vector) -> Pair:
    return sum_pairs(mul(a, b) for a, b in zip(left, right, strict=True))


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum_pairs(
                mul(left[i][k], right[k][j]) for k in range(len(right))
            )
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def transpose(matrix: Matrix) -> Matrix:
    return tuple(tuple(matrix[i][j] for i in range(len(matrix))) for j in range(len(matrix[0])))


def extended_gram(ell: Vector) -> Matrix:
    fixed = fixed_gram()
    rows = [list(row) + [ell[i]] for i, row in enumerate(fixed)]
    rows.append(list(ell) + [pair(2)])
    return tuple(tuple(row) for row in rows)


def involution(projection: Vector) -> Matrix:
    matrix = [[ZERO for _ in range(11)] for _ in range(11)]
    for i in range(10):
        matrix[i][i] = ONE
        matrix[i][10] = scale(projection[i], 2)
    matrix[10][10] = pair(-1)
    return tuple(tuple(row) for row in matrix)


def pair_json(value: Pair) -> list[int | str]:
    def scalar_json(entry: Fraction) -> int | str:
        return int(entry) if entry.denominator == 1 else str(entry)

    return [scalar_json(value[0]), scalar_json(value[1])]


def vector_json(vector: Vector) -> list[list[int | str]]:
    return [pair_json(value) for value in vector]


def case_audit(name: str, ell: Vector, claimed_w: Vector) -> dict[str, object]:
    fixed = fixed_gram()
    projection = solve_projection(ell)
    if matrix_vector(fixed, projection) != ell:
        raise AssertionError("F*p does not reconstruct ell")
    two_projection = tuple(scale(value, 2) for value in projection)
    integral = all(
        value[0].denominator == value[1].denominator == 1
        for value in two_projection
    )
    p_square = dot(projection, ell)
    q = scale(p_square, 2)
    schur_twice = sub(pair(4), q)
    identity_schur_sign = sign(schur_twice)
    conjugate_schur_sign = sign(conjugate(schur_twice))
    identity_signature = [10, 1] if identity_schur_sign > 0 else [9, 2]
    conjugate_signature = [11, 0] if conjugate_schur_sign > 0 else [10, 1]

    h = extended_gram(ell)
    j = involution(projection)
    identity_11: Matrix = tuple(
        tuple(ONE if i == k else ZERO for k in range(11)) for i in range(11)
    )
    j_squared = matmul(j, j)
    isometry = matmul(transpose(j), matmul(h, j)) == h

    on_interior = dot(ell, INTERIOR_SCALED)
    on_ray = dot(ell, RAY_SCALED)
    actual_w = projection[1:9]
    return {
        "name": name,
        "ell": vector_json(ell),
        "p": vector_json(projection),
        "two_p_is_integral": integral,
        "F_times_p_equals_ell": True,
        "actual_projected_E8_component": vector_json(actual_w),
        "claimed_w_matches_projection": actual_w == claimed_w,
        "q_equals_2_p_square": pair_json(q),
        "q_conjugate": pair_json(conjugate(q)),
        "four_minus_q": pair_json(schur_twice),
        "four_minus_q_conjugate": pair_json(conjugate(schur_twice)),
        "identity_signature_positive_negative": identity_signature,
        "conjugate_signature_positive_negative": conjugate_signature,
        "J_squared_is_identity": j_squared == identity_11,
        "J_is_H_isometry": isometry,
        "ell_on_400_times_interior": pair_json(on_interior),
        "ell_on_13_times_v61": pair_json(on_ray),
        "evaluation_signs": [sign(on_interior), sign(on_ray)],
        "cuts_v61": sign(on_interior) * sign(on_ray) == -1,
        "passes_required_signatures": identity_schur_sign > 0
        and conjugate_schur_sign > 0,
    }


def audit() -> dict[str, object]:
    gram_root = tuple(
        sum(E8_GRAM[i][j] * ROOT_VECTOR[j] for j in range(8))
        for i in range(8)
    )
    if gram_root != (0, 0, 0, 0, 0, 0, -1, 0):
        raise AssertionError("the declared E8 root convention changed")
    scalar = pair(4, 7)
    claimed_w = tuple(scale(scalar, entry) for entry in ROOT_VECTOR)
    derived_covector = tuple(
        sum_pairs(
            scale(claimed_w[j], E8_GRAM[i][j]) for j in range(8)
        )
        for i in range(8)
    )
    if derived_covector != (ZERO,) * 6 + (pair(-4, -7),) + (ZERO,):
        raise AssertionError("E8*w reconstruction changed")

    A = pair(-11, -18)
    C = pair(-1, -4)
    literal_ell = (A,) + (ZERO,) * 6 + (pair(8, 14),) + (ZERO,) + (C,)
    corrected_ell = (A,) + derived_covector + (C,)
    literal = case_audit("literal_rejected_transcript", literal_ell, claimed_w)
    corrected = case_audit("corrected_covector", corrected_ell, claimed_w)

    checks = {
        "literal_claimed_w_mismatch": literal["claimed_w_matches_projection"] is False,
        "literal_q_is_not_claimed_q": literal["q_equals_2_p_square"] != [2, -1],
        "literal_signature_gate_fails": literal["passes_required_signatures"] is False,
        "corrected_claimed_w_matches": corrected["claimed_w_matches_projection"] is True,
        "corrected_q_is_2_minus_phi": corrected["q_equals_2_p_square"] == [2, -1],
        "corrected_integral_involution": corrected["two_p_is_integral"] is True
        and corrected["J_squared_is_identity"] is True
        and corrected["J_is_H_isometry"] is True,
        "corrected_signatures_pass": corrected["passes_required_signatures"] is True,
        "corrected_does_not_cut": corrected["evaluation_signs"] == [-1, -1]
        and corrected["cuts_v61"] is False,
    }
    if not all(checks.values()):
        raise AssertionError(f"rejected-separator regression failed: {checks}")
    return {
        "status": "rejected_transcript_not_a_candidate",
        "root_cause": (
            "G*v=-e7, so w=(4+7phi)*v requires ell_E8=(-4-7phi)e7; "
            "the supplied +(8+14phi)e7 equals -2 times that covector."
        ),
        "E8_gram_times_v": list(gram_root),
        "literal": literal,
        "corrected": corrected,
        "checks": checks,
        "conclusion": (
            "The literal transcript appears to separate signs but fails the "
            "identity signature gate. The corrected integral carrier has "
            "q=2-phi and the required signatures, but both evaluations are "
            "negative, so it does not cut v61."
        ),
    }


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
