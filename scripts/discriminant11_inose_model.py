"""Exact Inose model for the singular K3 surface of discriminant 11."""

from __future__ import annotations

import json
from fractions import Fraction

import sympy as sp


def analyze() -> dict[str, object]:
    j = Fraction(-32768)
    j_scale = Fraction(1728)
    inose_a = Fraction(64, 9)
    inose_b = Fraction(539, 27)

    # Standard Inose parameters for equal j-invariants:
    # A^3=(j/1728)^2 and B^2=(1-j/1728)^2.
    parameter_identities = (
        inose_a**3 == (j / j_scale) ** 2
        and inose_b**2 == (1 - j / j_scale) ** 2
    )

    t = sp.symbols("t")
    a4 = -3 * sp.Rational(inose_a.numerator, inose_a.denominator) * t**4
    a6 = t**5 * (
        t**2
        - 2 * sp.Rational(inose_b.numerator, inose_b.denominator) * t
        + 1
    )
    c4 = sp.factor(-48 * a4)
    discriminant = sp.factor(-16 * (4 * a4**3 + 27 * a6**2))
    expected_discriminant = -16 * t**10 * (t - 1) ** 2 * (27 * t**2 - 2102 * t + 27)

    return {
        "surface_id": "inose_discriminant11",
        "base_field": "Q",
        "cm_discriminant": -11,
        "cm_j_invariant": str(j),
        "transcendental_lattice": [[2, 1], [1, 6]],
        "inose_A": str(inose_a),
        "inose_B": str(inose_b),
        "inose_parameter_identities_verified": parameter_identities,
        "weierstrass_equation": (
            "y^2 = x^3 - (64/3)*t^4*x + "
            "t^5*(t^2 - (1078/27)*t + 1)"
        ),
        "c4": str(c4),
        "discriminant": str(discriminant),
        "discriminant_factorization_verified": sp.expand(discriminant - expected_discriminant) == 0,
        "singular_fibers": [
            {"place": "0", "kodaira_type": "II*", "euler_number": 10},
            {"place": "infinity", "kodaira_type": "II*", "euler_number": 10},
            {"place": "1", "kodaira_type": "I2", "euler_number": 2},
            {"place": "root_1(27*t^2-2102*t+27)", "kodaira_type": "I1", "euler_number": 1},
            {"place": "root_2(27*t^2-2102*t+27)", "kodaira_type": "I1", "euler_number": 1},
        ],
        "fiber_euler_sum": 24,
        "shioda_tate_mordell_weil_rank": 1,
        "geometric_status": "explicit_singular_k3_model_verified",
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, sort_keys=True))

