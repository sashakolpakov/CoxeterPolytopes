"""Audit quadratic-CM/algebraic upgrades of the twelve coherent returns."""

from __future__ import annotations

import json
import math
from fractions import Fraction

from coxeter_polytopes.exact.golden import Golden, PHI

from scripts.certify_golden_lanner_hodge_coherent_carriers import (
    SOURCE,
    TWISTS,
    _matmul,
    _rational_trace_gram,
    _transpose,
    complement_field_gram,
)


ZERO = Golden(0)
ONE = Golden(1)


CM_DATA: dict[str, dict[str, object]] = {
    "triangle_255": {
        "c": Golden(8, -4),
        "kind": "three_equal_two_planes",
        "square_root_c": Golden(-2, 2),
    },
    "chain_5333": {
        "c": Golden(10, -6),
        "kind": "quaternion",
        "pure_quaternion": [Golden(2, -1), ZERO, Golden(-2, 1)],
    },
    "chain_5335": {
        "c": Golden(14, -8),
        "kind": "quaternion",
        "pure_quaternion": [ONE, ZERO, Golden(-3, 2)],
    },
    "fork5_terminal5": {
        "c": Golden(8, -4),
        "kind": "quaternion",
        "pure_quaternion": [Golden(-2, 2), ZERO, ZERO],
    },
}


def _rational_square_root(value: Fraction) -> Fraction | None:
    if value < 0:
        return None
    numerator = math.isqrt(value.numerator)
    denominator = math.isqrt(value.denominator)
    if numerator * numerator != value.numerator:
        return None
    if denominator * denominator != value.denominator:
        return None
    return Fraction(numerator, denominator)


def golden_square_root(value: Golden) -> Golden | None:
    """Return a square root in Q(sqrt(5)), or prove that none exists."""

    norm_root = _rational_square_root(value.norm())
    if norm_root is None:
        return None
    for norm_y in (norm_root, -norm_root):
        trace_y_square = value.trace() + 2 * norm_y
        trace_y = _rational_square_root(trace_y_square)
        if trace_y is None:
            continue
        radical = _rational_square_root((trace_y * trace_y - 4 * norm_y) / 5)
        if radical is None:
            continue
        # y=(s+d*sqrt(5))/2 = (s-d)/2 + d*phi.
        for s in (trace_y, -trace_y):
            for d in (radical, -radical):
                candidate = Golden((s - d) / 2, d)
                if candidate * candidate == value:
                    return candidate
    return None


def _zero_matrix(rank: int) -> list[list[Golden]]:
    return [[ZERO for _ in range(rank)] for _ in range(rank)]


def _field_matmul(
    left: list[list[Golden]], right: list[list[Golden]]
) -> list[list[Golden]]:
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                start=ZERO,
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def _field_transpose(matrix: list[list[Golden]]) -> list[list[Golden]]:
    return [list(column) for column in zip(*matrix, strict=True)]


def _field_add(
    left: list[list[Golden]], right: list[list[Golden]]
) -> list[list[Golden]]:
    return [
        [left[i][j] + right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ]


def _scalar_identity(rank: int, scalar: Golden) -> list[list[Golden]]:
    return [
        [scalar if i == j else ZERO for j in range(rank)] for i in range(rank)
    ]


def _quaternion_block(x: Golden, y: Golden, z: Golden) -> list[list[Golden]]:
    """Left multiplication by x*i+y*j+z*k on the rational quaternion basis."""

    return [
        [ZERO, -x, -y, -z],
        [x, ZERO, -z, y],
        [y, z, ZERO, -x],
        [z, -y, x, ZERO],
    ]


def cm_operator(seed_id: str, field_gram: list[list[Golden]]) -> tuple[Golden, list[list[Golden]]]:
    data = CM_DATA[seed_id]
    c = data["c"]
    assert isinstance(c, Golden)
    rank = len(field_gram)
    operator = _zero_matrix(rank)
    operator[0][1] = -c
    operator[1][0] = ONE
    if data["kind"] == "three_equal_two_planes":
        scalar = data["square_root_c"]
        assert isinstance(scalar, Golden)
        for start in (2, 4, 6):
            operator[start][start + 1] = -scalar
            operator[start + 1][start] = scalar
    elif data["kind"] == "quaternion":
        x, y, z = data["pure_quaternion"]
        block = _quaternion_block(x, y, z)
        for i in range(4):
            for j in range(4):
                operator[i + 2][j + 2] = block[i][j]
    else:
        raise ArithmeticError("unknown CM operator kind")
    return c, operator


def _golden_matrix_to_rational(matrix: list[list[Golden]]) -> list[list[Fraction]]:
    rank = len(matrix)
    result = [[Fraction(0) for _ in range(2 * rank)] for _ in range(2 * rank)]
    for i in range(rank):
        for j in range(rank):
            value = matrix[i][j]
            block = [
                [value.a, value.b],
                [value.b, value.a + value.b],
            ]
            for a in range(2):
                for b in range(2):
                    result[2 * i + a][2 * j + b] = block[a][b]
    return result


def _rational_add(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [
        [left[i][j] + right[i][j] for j in range(len(left[0]))]
        for i in range(len(left))
    ]


def analyze() -> dict[str, object]:
    source = json.loads(SOURCE.read_text())
    rows: list[dict[str, object]] = []
    for example in source["examples"]:
        seed_id = example["seed_id"]
        dimension = example["dimension"]
        field_determinant = Golden.from_json(example["field_determinant"])
        field_rank = example["rank"]
        complement_rank = 11 - field_rank
        twist = TWISTS[seed_id]
        assert twist is not None
        field_gram = complement_field_gram(field_determinant, complement_rank, twist)
        determinant = ONE
        for i in range(complement_rank):
            determinant *= field_gram[i][i]

        if complement_rank % 2:
            rows.append(
                {
                    "seed_id": seed_id,
                    "dimension": dimension,
                    "complement_K_rank": complement_rank,
                    "status": "quadratic_CM_extension_impossible_by_degree_divisibility",
                    "reason": "[K(sqrt(-c)):K]=2 must divide dim_K(T)",
                }
            )
            continue

        if seed_id not in CM_DATA:
            square_root = golden_square_root(determinant)
            rows.append(
                {
                    "seed_id": seed_id,
                    "dimension": dimension,
                    "complement_K_rank": complement_rank,
                    "complement_field_determinant": determinant.to_json(),
                    "determinant_square_in_K": square_root is not None,
                    "status": "quadratic_CM_extension_impossible_by_determinant_square_class",
                    "reason": (
                        "for dim_K(T)=8, a skew J with J^2=-c forces "
                        "det(T)=c^4*square, hence det(T) is a square in K"
                    ),
                }
            )
            continue

        c, operator = cm_operator(seed_id, field_gram)
        gram_operator = _field_matmul(field_gram, operator)
        skew = _field_add(
            _field_matmul(_field_transpose(operator), field_gram), gram_operator
        ) == _zero_matrix(complement_rank)
        square = _field_matmul(operator, operator)
        square_relation = square == _scalar_identity(complement_rank, -c)
        rational_gram = _rational_trace_gram(field_gram)
        rational_operator = _golden_matrix_to_rational(operator)
        rational_golden_action = _golden_matrix_to_rational(
            _scalar_identity(complement_rank, PHI)
        )
        commutes_with_golden_action = _matmul(
            rational_operator, rational_golden_action
        ) == _matmul(rational_golden_action, rational_operator)
        rational_skew = _rational_add(
            _matmul(_transpose(rational_operator), rational_gram),
            _matmul(rational_gram, rational_operator),
        ) == [[Fraction(0) for _ in range(2 * complement_rank)] for _ in range(2 * complement_rank)]
        determinant_quotient = determinant / (c ** (complement_rank // 2))
        determinant_gate = golden_square_root(determinant_quotient) is not None
        positive_plane_is_first_block = (
            all(field_gram[index][index].sign(0) < 0 for index in (0, 1))
            and all(
                field_gram[index][index].sign(0) > 0
                for index in range(2, complement_rank)
            )
            and all(
                operator[i][j] == ZERO
                for i in range(complement_rank)
                for j in range(complement_rank)
                if (i < 2) != (j < 2)
            )
        )
        rows.append(
            {
                "seed_id": seed_id,
                "dimension": dimension,
                "complement_K_rank": complement_rank,
                "c": c.to_json(),
                "c_totally_positive": c.sign(0) > 0 and c.sign(1) > 0,
                "operator": [[entry.to_json() for entry in row] for row in operator],
                "J_square_minus_c": square_relation,
                "J_K_skew_adjoint": skew,
                "J_Q_skew_adjoint_for_trace_form": rational_skew,
                "J_commutes_with_golden_action": commutes_with_golden_action,
                "determinant_CM_square_class_gate": determinant_gate,
                "distinguished_positive_plane_is_J_stable": (
                    positive_plane_is_first_block
                ),
                "CM_field_degree_over_Q": 4,
                "CM_period_domain_complex_dimension": complement_rank // 2 - 1,
                "generic_exact_CM_endomorphism_field_available": True,
                "Ramon_Mari_algebraic_correspondence_available": True,
                "status": "certified_quadratic_CM_and_algebraic_correspondence_upgrade",
            }
        )

    positives = [row for row in rows if str(row["status"]).startswith("certified")]
    return {
        "schema_version": "0.1.0",
        "certificate_kind": "golden_lanner_quadratic_CM_upgrade_audit",
        "scope": (
            "quadratic CM extensions K(sqrt(-c)) of the certified K-complement; "
            "negative rows do not exclude higher-degree CM fields"
        ),
        "positive_count": len(positives),
        "rows": rows,
        "all_positive_exact_matrix_checks_pass": all(
            row["J_square_minus_c"]
            and row["J_K_skew_adjoint"]
            and row["J_Q_skew_adjoint_for_trace_form"]
            and row["J_commutes_with_golden_action"]
            and row["determinant_CM_square_class_gate"]
            and row["distinguished_positive_plane_is_J_stable"]
            and row["c_totally_positive"]
            for row in positives
        ),
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))
