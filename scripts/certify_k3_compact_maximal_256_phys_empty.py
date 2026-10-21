#!/usr/bin/env python3
"""Exact replay for the maximal (2,5,6) compact K3 PHYS obstruction.

The only non-executable input is the published finite-maximality
classification of the orientation-preserving (2,5,6) triangle group.  Every
number-field, Gram, energy-shell, and square-class calculation used after
that input is replayed exactly here.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import hashlib
from math import isqrt
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_compact_maximal_256_phys_empty.md"
ARCHIVE = ROOT / "data/enumerations/k3_compact_maximal_256_phys_empty.json"


def digest(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def is_fraction_square(value: Fraction) -> bool:
    if value < 0:
        return False
    return (
        isqrt(value.numerator) ** 2 == value.numerator
        and isqrt(value.denominator) ** 2 == value.denominator
    )


@dataclass(frozen=True)
class Q5:
    """An exact element a+b*sqrt(5)."""

    a: Fraction
    b: Fraction

    def __add__(self, other: "Q5") -> "Q5":
        return Q5(self.a + other.a, self.b + other.b)

    def __sub__(self, other: "Q5") -> "Q5":
        return Q5(self.a - other.a, self.b - other.b)

    def __mul__(self, other: "Q5") -> "Q5":
        return Q5(
            self.a * other.a + 5 * self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def norm(self) -> Fraction:
        return self.a * self.a - 5 * self.b * self.b

    def trace(self) -> Fraction:
        return 2 * self.a

    def embedding_sign(self, *, conjugate: bool = False) -> int:
        """Return the exact sign at either real embedding."""

        a = self.a
        b = -self.b if conjugate else self.b
        if b == 0:
            return (a > 0) - (a < 0)
        if a == 0:
            return (b > 0) - (b < 0)
        if (a > 0) == (b > 0):
            return (a > 0) - (a < 0)
        comparison = a * a - 5 * b * b
        if comparison == 0:
            raise AssertionError("a nonzero rational ratio cannot equal sqrt(5)")
        if comparison > 0:
            return (a > 0) - (a < 0)
        return (b > 0) - (b < 0)

    def inverse(self) -> "Q5":
        norm = self.norm()
        if norm == 0:
            raise ZeroDivisionError
        return Q5(self.a / norm, -self.b / norm)

    def __truediv__(self, other: "Q5") -> "Q5":
        return self * other.inverse()

    def coordinates(self) -> list[str]:
        return [str(self.a), str(self.b)]


def is_q5_square(value: Q5) -> bool:
    """Decide exactly whether value is a square in Q(sqrt(5))."""

    norm = value.norm()
    if not is_fraction_square(norm):
        return False
    root_norm = Fraction(isqrt(norm.numerator), isqrt(norm.denominator))
    for signed_root_norm in (root_norm, -root_norm):
        x_squared = (value.a + signed_root_norm) / 2
        if not is_fraction_square(x_squared):
            continue
        x = Fraction(
            isqrt(x_squared.numerator), isqrt(x_squared.denominator)
        )
        if x == 0:
            if value.b != 0:
                continue
            y_squared = value.a / 5
            if is_fraction_square(y_squared):
                return True
            continue
        y = value.b / (2 * x)
        if x * x + 5 * y * y == value.a:
            return True
    return False


def energy(c: int) -> Q5:
    """Return u_c=-1+c/sqrt(5)=-1+(c/5)*sqrt(5)."""

    return Q5(Fraction(-1), Fraction(c, 5))


def golden_energy_shell() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for c in range(-2, 3):
        u = energy(c)
        if u.trace() != -2:
            raise AssertionError("energy trace changed")
        if u.embedding_sign() != -1 or u.embedding_sign(conjugate=True) != -1:
            raise AssertionError("energy shell contains a nonnegative conjugate")
        numerator_norm = c * c - 5
        if numerator_norm % 3 == 0:
            raise AssertionError("energy numerator unexpectedly nonunit at (3)")
        rows.append(
            {
                "c": c,
                "element": f"-1+({c})/sqrt(5)",
                "sqrt5_coordinates": u.coordinates(),
                "field_norm": str(u.norm()),
                "numerator_norm_c2_minus_5": numerator_norm,
                "valuation_at_inert_3": 0,
            }
        )
    return rows


def square_class_census() -> dict[str, Any]:
    three = Q5(Fraction(3), Fraction(0))
    rows: list[dict[str, Any]] = []
    squares = 0
    for c in range(-2, 3):
        for e in range(-2, 3):
            quotient = energy(c) / energy(e) / three
            is_square = is_q5_square(quotient)
            squares += int(is_square)
            rows.append(
                {
                    "c": c,
                    "e": e,
                    "u_c_over_3u_e": quotient.coordinates(),
                    "is_square": is_square,
                }
            )
    if squares:
        raise AssertionError("the (2,5,6) PHYS square-class obstruction failed")
    return {
        "tested_ordered_pairs": len(rows),
        "square_count": squares,
        "rows": rows,
    }


def raw_gram_certificate() -> dict[str, Any]:
    # G=[[-6,0,3],[0,-2,phi],[3,phi,-2]], phi^2=phi+1.
    # det(G)=6*phi.  The leading minors are -6, 12, 6*phi.
    # Hence the distinguished inertia is (1,2) and the conjugate inertia
    # is (0,3).
    zero = Q5(Fraction(0), Fraction(0))
    one = Q5(Fraction(1), Fraction(0))
    phi = Q5(Fraction(1, 2), Fraction(1, 2))
    gram = [
        [Q5(Fraction(-6), Fraction(0)), zero, Q5(Fraction(3), Fraction(0))],
        [zero, Q5(Fraction(-2), Fraction(0)), phi],
        [Q5(Fraction(3), Fraction(0)), phi, Q5(Fraction(-2), Fraction(0))],
    ]

    determinant = (
        gram[0][0]
        * (gram[1][1] * gram[2][2] - gram[1][2] * gram[2][1])
        - gram[0][1]
        * (gram[1][0] * gram[2][2] - gram[1][2] * gram[2][0])
        + gram[0][2]
        * (gram[1][0] * gram[2][1] - gram[1][1] * gram[2][0])
    )
    if phi * phi != phi + one:
        raise AssertionError("golden relation changed")
    if determinant != Q5(Fraction(3), Fraction(3)):
        raise AssertionError("raw Gram determinant is not 6*phi")

    leading_pivots = [
        Q5(Fraction(-6), Fraction(0)),
        Q5(Fraction(-2), Fraction(0)),
        determinant / Q5(Fraction(12), Fraction(0)),
    ]
    distinguished_signs = [x.embedding_sign() for x in leading_pivots]
    conjugate_signs = [
        x.embedding_sign(conjugate=True) for x in leading_pivots
    ]
    if distinguished_signs != [-1, -1, 1]:
        raise AssertionError("distinguished Gram inertia changed")
    if conjugate_signs != [-1, -1, -1]:
        raise AssertionError("conjugate Gram is not negative definite")

    cosine_23 = gram[1][2] * gram[1][2] / (gram[1][1] * gram[2][2])
    cosine_31 = gram[2][0] * gram[2][0] / (gram[2][2] * gram[0][0])
    if cosine_23 != phi * phi / Q5(Fraction(4), Fraction(0)):
        raise AssertionError("order-five angle changed")
    if cosine_31 != Q5(Fraction(3, 4), Fraction(0)):
        raise AssertionError("order-six angle changed")

    two = Q5(Fraction(2), Fraction(0))
    cartan_13 = two * gram[0][2] / gram[0][0]
    cartan_31 = two * gram[2][0] / gram[2][2]
    cartan_23 = two * gram[1][2] / gram[1][1]
    cartan_32 = two * gram[2][1] / gram[2][2]
    if [cartan_13, cartan_31, cartan_23, cartan_32] != [
        Q5(Fraction(-1), Fraction(0)),
        Q5(Fraction(-3), Fraction(0)),
        Q5(Fraction(-1, 2), Fraction(-1, 2)),
        Q5(Fraction(-1, 2), Fraction(-1, 2)),
    ]:
        raise AssertionError("integral Cartan coefficients changed")

    angle_sum = Fraction(1, 2) + Fraction(1, 5) + Fraction(1, 6)
    if angle_sum != Fraction(13, 15) or not angle_sum < 1:
        raise AssertionError("the triangle is no longer compact hyperbolic")

    return {
        "basis": ["r1", "r2", "r3"],
        "gram_phi": [
            [[-6, 0], [0, 0], [3, 0]],
            [[0, 0], [-2, 0], [0, 1]],
            [[3, 0], [0, 1], [-2, 0]],
        ],
        "entry_encoding": "[a,b] means a+b*phi",
        "leading_minors": ["-6", "12", "6*phi"],
        "LDL_pivot_signs_distinguished": distinguished_signs,
        "LDL_pivot_signs_conjugate": conjugate_signs,
        "distinguished_phi_sign": 1,
        "conjugate_phi_sign": -1,
        "distinguished_inertia": [1, 2, 0],
        "conjugate_inertia": [0, 3, 0],
        "coxeter_orders": {"r1_r2": 2, "r2_r3": 5, "r3_r1": 6},
        "normalized_cosine_squares": {
            "r1_r2": "0",
            "r2_r3": "phi^2/4=cos(pi/5)^2",
            "r3_r1": "3/4=cos(pi/6)^2",
        },
        "cartan_coefficients_nonzero": {
            "C_13": "-1",
            "C_31": "-3",
            "C_23": "-phi",
            "C_32": "-phi",
        },
        "order_module_preserved": True,
        "cycle_products": ["3", "phi^2"],
        "cycle_order": "Z[phi^2]=Z[phi]=O_k",
        "compact_angle_sum": "1/2+1/5+1/6=13/15<1",
        "field": "Q(sqrt(5))",
        "field_degree": 2,
        "module_rank_over_field": 3,
        "trace_rank": 6,
    }


def build_certificate() -> dict[str, Any]:
    gram = raw_gram_certificate()
    shell = golden_energy_shell()
    census = square_class_census()
    mod_three_values = [(x * x - x - 1) % 3 for x in range(3)]
    if 0 in mod_three_values:
        raise AssertionError("3 is no longer inert in the golden field")
    if len(shell) != 5 or census["tested_ordered_pairs"] != 25:
        raise AssertionError("finite PHYS census changed")
    if gram["trace_rank"] > 20:
        raise AssertionError("this obstruction must occur below the K3 rank gate")

    return {
        "schema": "coxeter.k3-compact-maximal-256-phys-empty.v1",
        "report_sha256": sha256_file(REPORT),
        "raw_packet": gram,
        "published_maximality_input": {
            "orientation_preserving_triangle_group": "Delta(2,5,6)",
            "arithmetic": True,
            "finitely_maximal": True,
            "full_triangle_reflection_group_maximal": True,
            "source": (
                "C. Maclachlan, Triangle subgroups of hyperbolic "
                "tetrahedral groups, Pacific J. Math. 176 (1996), "
                "Table 3 and Table 4"
            ),
            "classification_boundary": (
                "Maximality is cited; all PHYS arithmetic is replayed here."
            ),
        },
        "codifferent_energy_shell": {
            "order": "O_k=Z[phi]",
            "codifferent": "O_k^vee=(1/sqrt(5))*O_k",
            "parameterization": (
                "Tr((a+b*phi)/sqrt(5))=b; Tr=-2 gives "
                "u_c=-1+c/sqrt(5)"
            ),
            "total_negativity": "|c|<sqrt(5)",
            "rows": shell,
            "cardinality": len(shell),
        },
        "intrinsic_diagonal_PHYS_criterion": {
            "square_class_sets": (
                "C_a={[u/h_a]:u in E_-2(O0)} in k^x/k^x2"
            ),
            "necessary_and_sufficient_for_diagonal_normalization": (
                "intersection_a C_a is nonempty"
            ),
            "invariant_under_facet_line_rescaling": True,
            "full_even_integral_module_also_required": True,
            "raw_maximality_controls_this_intersection": False,
            "triangle_256_pairwise_intersection_empty": True,
        },
        "inert_prime_obstruction": {
            "prime": 3,
            "minimal_polynomial_mod_3": "x^2-x-1",
            "values_at_0_1_2_mod_3": mod_three_values,
            "irreducible_mod_3": True,
            "three_is_inert": True,
            "all_energy_values_are_3_adic_units": True,
            "facet_norm_ratio_h1_over_h3": 3,
            "required_equation": "u_1/u_3=3*(a_1/a_3)^2",
            "left_valuation_parity": 0,
            "right_valuation_parity": 1,
            "contradiction": True,
        },
        "independent_exact_square_census": census,
        "compiler_verdict": {
            "rank_gate_passes": True,
            "trace_rank": 6,
            "K3_complement_rank_if_PHYS_existed": 16,
            "discriminant_length_upper_bound": 6,
            "stable_embedding_inequality": "16>=6+2",
            "EMB_would_be_automatic_after_PHYS": True,
            "PHYS": "EMPTY",
            "EMB_reached": False,
            "ROOT_reached": False,
            "K3Comp_nonempty": False,
            "raw_maximality_implies_K3Comp_nonempty": False,
            "scope": (
                "primitive full-trace, full-Picard, same-vector K3 return"
            ),
            "projected_root_or_larger_Picard_return_excluded": False,
            "alternative_AG_carrier_excluded": False,
        },
        "sharp_positive_repair": {
            "PHYS_and_EMB_must_be_assumed_or_proved": True,
            "zero_projection_roots_must_be_excluded": True,
            "wall_selective_native_saturation": (
                "For every ambient (-2) root whose nonzero restricted mirror "
                "meets P interior, its k-native reflection preserves M."
            ),
            "then_raw_maximality_removes_interior_ROOT_cuts": True,
            "facet_witnesses_still_required": True,
            "period_or_ample_gate_remaining_after_these_conditions": False,
        },
    }


def encoded_certificate() -> str:
    payload = build_certificate()
    return json.dumps(
        {"payload": payload, "payload_sha256": digest(payload)},
        indent=2,
        sort_keys=True,
    ) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-output", action="store_true")
    parser.add_argument("--check-output", action="store_true")
    args = parser.parse_args()

    encoded = encoded_certificate()
    if args.write_output:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(encoded, encoding="utf-8")
    if args.check_output:
        if not ARCHIVE.exists() or ARCHIVE.read_text(encoding="utf-8") != encoded:
            raise SystemExit("archived maximal (2,5,6) PHYS certificate differs")
    if not args.write_output and not args.check_output:
        print(encoded, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
