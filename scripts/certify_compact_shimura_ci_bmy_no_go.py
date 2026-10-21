"""Exact arithmetic for the compact type-IV complete-intersection BMY no-go.

For a compact neat type-IV quotient X of complex dimension n, write lambda
for the weight-one automorphic class.  Hirzebruch proportionality gives

    c1(TX) = -n lambda,
    c2(TX) = (n^2-n+2)/2 lambda^2.

If S is the complete intersection of n-2 divisors numerically equivalent to
k lambda, this module verifies the adjunction formulas and the factorization
of 3 c2(S)-K_S^2 used in the companion report.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb
from typing import Any


def very_ample_multiple(n: int) -> int:
    """Demailly's universal multiple after K_X is numerically n*lambda."""

    if n < 2:
        raise ValueError("the surface-section calculation requires n >= 2")
    m = 2 + comb(3 * n + 1, n)
    return 2 * n + m


def coefficients(n: int, k: int) -> dict[str, Fraction]:
    """Return coefficients after factoring out k^(n-2) * lambda^n."""

    if n < 2 or k < 1:
        raise ValueError("require n >= 2 and k >= 1")
    r = n - 2
    c2_x = Fraction(n * n - n + 2, 2)
    canonical = Fraction(n + r * k)
    c2_surface = (
        c2_x
        + Fraction(r * k * n)
        + Fraction(r * (r + 1) * k * k, 2)
    )
    defect = 3 * c2_surface - canonical * canonical
    defect_expected = (
        Fraction(n * n - 3 * n + 6, 2)
        + Fraction(n * (n - 2) * k)
        + Fraction((n - 2) * (n + 1) * k * k, 2)
    )
    # If C_i are separated mirror curves, their normal bundles are
    # lambda^{-1}; hence C_i^2=-deg_lambda(C_i).  For a reduced SNC mirror
    # boundary, log BMY adds (K.C)+2*sum(C_i^2)+sum_{i<j}(C_i.C_j), so the
    # coefficient of each individual lambda-degree is canonical-2.
    reduced_mirror_degree_coefficient = canonical - 2
    return {
        "ambient_c2": c2_x,
        "surface_canonical": canonical,
        "surface_c2": c2_surface,
        "bmy_defect": defect,
        "bmy_defect_expected": defect_expected,
        "reduced_mirror_degree_coefficient": reduced_mirror_degree_coefficient,
    }


def _encode(value: Fraction) -> int | str:
    return value.numerator if value.denominator == 1 else str(value)


def analyze() -> dict[str, Any]:
    coefficient_identity = True
    strict_positivity = True
    reduced_log_strict_positivity = True
    rows: list[dict[str, Any]] = []
    for n in range(2, 31):
        k = very_ample_multiple(n)
        data = coefficients(n, k)
        coefficient_identity &= data["bmy_defect"] == data["bmy_defect_expected"]
        strict_positivity &= data["bmy_defect"] > 0
        reduced_log_strict_positivity &= (
            data["reduced_mirror_degree_coefficient"] >= 0
        )
        if n in {2, 3, 4, 8, 9, 10}:
            rows.append(
                {
                    "dimension": n,
                    "universal_very_ample_multiple": k,
                    "canonical_coefficient": _encode(data["surface_canonical"]),
                    "c2_coefficient": _encode(data["surface_c2"]),
                    "bmy_defect_coefficient": _encode(data["bmy_defect"]),
                    "reduced_mirror_degree_coefficient": _encode(
                        data["reduced_mirror_degree_coefficient"]
                    ),
                }
            )

    bugaenko_n = 8
    bugaenko_k = very_ample_multiple(bugaenko_n)
    bugaenko = coefficients(bugaenko_n, bugaenko_k)

    # The polynomial coefficients prove positivity for every n >= 2 and k >= 1:
    # the constant has negative discriminant, and the other two are nonnegative.
    constant_discriminant = (-3) ** 2 - 4 * 1 * 6
    positivity_proof = (
        constant_discriminant < 0
        and all(n * (n - 2) >= 0 for n in range(2, 100))
        and all((n - 2) * (n + 1) >= 0 for n in range(2, 100))
    )

    return {
        "theorem": "compact_type_iv_complete_intersection_bmy_cancellation",
        "ambient_chern_classes": {
            "c1_TX": "-n*lambda",
            "c2_TX": "(n^2-n+2)/2*lambda^2",
        },
        "surface_bmy_defect": (
            "k^(n-2)*lambda^n*((n^2-3n+6)/2+n(n-2)k+"
            "(n-2)(n+1)k^2/2)"
        ),
        "coefficient_identity_verified": coefficient_identity,
        "strictly_positive_for_all_n_ge_2_k_ge_1": positivity_proof,
        "reduced_snc_mirror_log_bmy_strictly_positive_conditional_on_snc": (
            positivity_proof and reduced_log_strict_positivity
        ),
        "reduced_snc_mirror_log_bmy_defect": (
            "ordinary_defect+(n+(n-2)k-2)*sum(deg_lambda(C_i))+"
            "sum_{i<j}(C_i.C_j)"
        ),
        "sample_rows": rows,
        "bugaenko_h8": {
            "universal_very_ample_multiple": bugaenko_k,
            "canonical_coefficient": _encode(bugaenko["surface_canonical"]),
            "c2_coefficient": _encode(bugaenko["surface_c2"]),
            "bmy_defect_coefficient": _encode(bugaenko["bmy_defect"]),
            "symbolic_bmy_defect_coefficient": "27*k^2+48*k+23",
        },
        "ordinary_bmy_dimension_or_facet_bound": False,
        "unresolved_multiple_point_log_resolution_may_add_feedback": True,
        "surface_degree_inequality_stronger_than_ambient_degree_injection": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compact", action="store_true")
    args = parser.parse_args()
    print(json.dumps(analyze(), indent=None if args.compact else 2, sort_keys=True))


if __name__ == "__main__":
    main()
