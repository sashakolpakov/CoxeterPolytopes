#!/usr/bin/env python3
"""Exact lattice data for the OG10 Borcherds-route obstruction.

The external modular-form inputs are Wang's classification theorem and
Bruinier's Whittaker-lift obstruction principle.  This script certifies the
convention-sensitive lattice arithmetic used to apply them:

* the ordinary type-IV extension ``U + N(-1)`` has signature ``(22,2)``,
  determinant ``44``, discriminant form
  ``<2/11> + <1/2> + <1/2>``, and level ``44``;
* the minimal sign-changing field-line extension has signatures
  ``((10,2),(12,0))``; and
* its trace lattice has discriminant group
  ``Z/11 + (Z/2)^4``, order ``176``, and level ``44``.

See ``results/og10_borcherds_obstruction.md`` for theorem statements and the
two essential non-implications.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import lcm

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import inertia
from coxeter_polytopes.exact.trace_lattice import INV_SQRT5
from coxeter_polytopes.lattices import e8_cartan_gram


def _smith_nonunits(gram: Matrix) -> list[int]:
    smith = smith_normal_form(gram, domain=ZZ)
    return sorted(
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if abs(int(smith[index, index])) > 1
    )


def _level(quadratic_values: list[Fraction]) -> int:
    """Level for q:A -> Q/2Z from generator quadratic values."""

    result = 1
    for value in quadratic_values:
        candidate = 1
        while True:
            scaled = candidate * value
            if scaled.denominator == 1 and scaled.numerator % 2 == 0:
                break
            candidate += 1
        result = lcm(result, candidate)
    return result


def _signature_payload(gram: Matrix) -> dict[str, int]:
    value = inertia([[int(entry) for entry in row] for row in gram.tolist()])
    return {
        "positive": value.positive,
        "negative": value.negative,
        "zero": value.zero,
    }


def _integer_matrix_payload(gram: Matrix) -> list[list[int]]:
    return [[int(entry) for entry in row] for row in gram.tolist()]


def _fraction_payload(value: object) -> Fraction:
    return Fraction(int(value.p), int(value.q))


def certify() -> dict[str, object]:
    u = Matrix([[0, 1], [1, 0]])
    t11 = Matrix([[2, 1], [1, 6]])
    e8 = Matrix(e8_cartan_gram())
    a1 = Matrix([[2]])

    # Wang uses signature (l,2).  Since -U is isometric to U,
    # M_0=U+N(-1)=2U+T_11+2E_8+2A_1.
    ordinary = Matrix.diag(u, u, t11, e8, e8, a1, a1)
    t11_dual_generator = t11.inv() * Matrix([0, 1])
    t11_q = _fraction_payload(
        (t11_dual_generator.T * t11 * t11_dual_generator)[0]
    )
    a1_q = _fraction_payload(a1.inv()[0, 0])
    ordinary_q = [t11_q, a1_q, a1_q]

    # Add to H the even O-line with coefficient c=2(1-phi).  It is negative
    # at the identity embedding and positive at the conjugate embedding.
    field_line = 2 * (Golden(1) - PHI)
    natural_scalar = PHI * INV_SQRT5
    field_line_quadratic = natural_scalar * (Golden(1) - PHI)
    target_discriminant = natural_scalar

    # For Q_line(x)=-x^2/sqrt(5), x in O with basis (1,phi), the trace
    # bilinear form is Tr(-2xy/sqrt(5)).
    auxiliary_trace = Matrix([[0, -2], [-2, -2]])
    diagonalizing_change = Matrix([[1, 0], [-1, 1]])
    auxiliary_diagonal = diagonalizing_change.T * auxiliary_trace * diagonalizing_change

    # The H-part of the trace lattice is N(-1).  In the displayed
    # decomposition -U has again been replaced by U.
    field_trace = Matrix.diag(u, t11, e8, e8, a1, a1, auxiliary_trace)
    field_trace_q = [
        t11_q,
        a1_q,
        a1_q,
        a1_q,
        -a1_q,
    ]

    result = {
        "ordinary_type_iv": {
            "isometry_type": "2U + T11 + 2E8 + 2A1",
            "rank": ordinary.rows,
            "signature": _signature_payload(ordinary),
            "determinant": int(ordinary.det()),
            "smith_nonunits": _smith_nonunits(ordinary),
            "discriminant_group": "Z/11 + (Z/2)^2",
            "quadratic_generator_values_mod_2": [
                str(value) for value in ordinary_q
            ],
            "t11_dual_generator": [
                str(_fraction_payload(entry)) for entry in t11_dual_generator
            ],
            "level": _level(ordinary_q),
            "borcherds_input_weight": -10,
            "dual_cusp_weight": 12,
        },
        "totally_real_extension": {
            "field_line_coefficient": str(field_line),
            "field_line_signs": [field_line.sign(0), field_line.sign(1)],
            "field_signatures": [[10, 2], [12, 0]],
            "natural_quadratic_scalar": "phi/sqrt(5)",
            "field_line_quadratic": "-1/sqrt(5)",
            "field_line_identity_verified": field_line_quadratic == -INV_SQRT5,
            "norm_two_target_m": "phi/sqrt(5)",
            "norm_two_target_totally_positive": [
                target_discriminant.sign(0),
                target_discriminant.sign(1),
            ],
            "trace_auxiliary_gram": _integer_matrix_payload(auxiliary_trace),
            "trace_auxiliary_diagonal": _integer_matrix_payload(
                auxiliary_diagonal
            ),
            "trace_rank": field_trace.rows,
            "trace_signature": _signature_payload(field_trace),
            "trace_determinant": int(field_trace.det()),
            "trace_smith_nonunits": _smith_nonunits(field_trace),
            "trace_discriminant_group": "Z/11 + (Z/2)^4",
            "trace_quadratic_generator_values_mod_2": [
                str(value) for value in field_trace_q
            ],
            "trace_level": _level(field_trace_q),
            "weil_module_dimension": abs(int(field_trace.det())),
            "whittaker_input_weight": [-4, 6],
            "obstruction_cusp_weight": [6, 6],
        },
    }
    verified = (
        result["ordinary_type_iv"]["rank"] == 24
        and result["ordinary_type_iv"]["signature"]
        == {"positive": 22, "negative": 2, "zero": 0}
        and result["ordinary_type_iv"]["determinant"] == 44
        and result["ordinary_type_iv"]["smith_nonunits"] == [2, 22]
        and result["ordinary_type_iv"]["t11_dual_generator"]
        == ["-1/11", "2/11"]
        and result["ordinary_type_iv"]["level"] == 44
        and result["totally_real_extension"]["field_line_signs"] == [-1, 1]
        and result["totally_real_extension"]["field_line_identity_verified"]
        and result["totally_real_extension"]["norm_two_target_totally_positive"]
        == [1, 1]
        and auxiliary_diagonal == Matrix.diag(2, -2)
        and result["totally_real_extension"]["trace_signature"]
        == {"positive": 22, "negative": 2, "zero": 0}
        and result["totally_real_extension"]["trace_determinant"] == 176
        and result["totally_real_extension"]["trace_smith_nonunits"]
        == [2, 2, 2, 22]
        and result["totally_real_extension"]["trace_level"] == 44
    )
    return {
        "verified": verified,
        **result,
        "ordinary_conclusion": "no_reflective_Borcherds_product_by_Wang_Theorem_3_5",
        "scope": (
            "the_classification_closes_the_ordinary_type_IV_route_but_does_not_"
            "decide_the_index_of_the_real_quadratic_norm_two_reflection_group"
        ),
    }


if __name__ == "__main__":
    print(json.dumps(certify(), indent=2, sort_keys=True))
