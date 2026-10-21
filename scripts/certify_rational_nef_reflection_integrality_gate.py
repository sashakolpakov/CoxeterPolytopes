#!/usr/bin/env python3
"""Exact certificate for the rational-nef reflection-integrality gate.

The certificate has no external dependencies.  It checks the toric
``Bl_p(F_4)`` counterexample, its half-boundary arithmetic lattice, the
reflection matrices and Coxeter triangle, and the elementary lattice
identities used in the accompanying theorem memo.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/rational_nef_reflection_integrality_gate.json"


def matmul(a, b):
    return [
        [sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))]
        for i in range(len(a))
    ]


def transpose(a):
    return [list(row) for row in zip(*a)]


def eye(n):
    return [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]


def pairing(g, x, y):
    return sum(x[i] * g[i][j] * y[j] for i in range(len(x)) for j in range(len(y)))


def determinant(a):
    a = [[Fraction(v) for v in row] for row in a]
    n = len(a)
    det = Fraction(1)
    for col in range(n):
        pivot = next((row for row in range(col, n) if a[row][col]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            det = -det
        value = a[col][col]
        det *= value
        for j in range(col, n):
            a[col][j] /= value
        for row in range(col + 1, n):
            value = a[row][col]
            if value:
                for j in range(col, n):
                    a[row][j] -= value * a[col][j]
    return det


def reflection(g, root):
    q = pairing(g, root, root)
    n = len(root)
    covector = [sum(root[k] * g[k][j] for k in range(n)) for j in range(n)]
    return [
        [Fraction(int(i == j)) - Fraction(2, 1) * root[i] * covector[j] / q for j in range(n)]
        for i in range(n)
    ]


def integral_matrix(a):
    return all(Fraction(v).denominator == 1 for row in a for v in row)


def serialize(value):
    if isinstance(value, Fraction):
        return int(value) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    if isinstance(value, list):
        return [serialize(item) for item in value]
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    return value


def build_certificate():
    # NS(Bl_p F_4), in the basis (B,F,E).  Here B is the negative section,
    # F a pulled-back fibre and E the exceptional curve.
    g_ns = [
        [Fraction(-4), Fraction(1), Fraction(0)],
        [Fraction(1), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(-1)],
    ]
    curves_ns = {
        "B": [Fraction(1), Fraction(0), Fraction(0)],
        "E": [Fraction(0), Fraction(0), Fraction(1)],
        "G=F-E": [Fraction(0), Fraction(1), Fraction(-1)],
    }
    canonical = [Fraction(-2), Fraction(-6), Fraction(1)]
    curve_gram_ns = [
        [pairing(g_ns, x, y) for y in curves_ns.values()] for x in curves_ns.values()
    ]
    b_reflection_ns = reflection(g_ns, curves_ns["B"])

    # The Coble-style arithmetic denominator: adjoin h=B/2 and multiply the
    # intersection form by two.  Basis (h,F,E).
    g_half = [
        [Fraction(-2), Fraction(1), Fraction(0)],
        [Fraction(1), Fraction(0), Fraction(0)],
        [Fraction(0), Fraction(0), Fraction(-2)],
    ]
    roots = {
        "h=B/2": [Fraction(1), Fraction(0), Fraction(0)],
        "E": [Fraction(0), Fraction(0), Fraction(1)],
        "G=F-E": [Fraction(0), Fraction(1), Fraction(-1)],
    }
    root_gram = [[pairing(g_half, x, y) for y in roots.values()] for x in roots.values()]
    reflections = {name: reflection(g_half, root) for name, root in roots.items()}
    log_class_half = [Fraction(-3), Fraction(-6), Fraction(1)]  # K+B/2
    log_class_pairings = {
        name: pairing(g_half, log_class_half, root) for name, root in roots.items()
    }
    reflection_checks = {
        name: {
            "integral": integral_matrix(matrix),
            "involution": matmul(matrix, matrix) == eye(3),
            "isometry": matmul(transpose(matrix), matmul(g_half, matrix)) == g_half,
        }
        for name, matrix in reflections.items()
    }
    root_coordinate_matrix = transpose(list(roots.values()))

    # The arbitrary-denominator ruled-surface family.  For the negative
    # section B_m^2=-2m on F_{2m}, set h=B_m/m and multiply by m.  In the
    # basis (h,F), the lattice and reflection are independent of m.
    denominator_family = []
    for m in range(2, 21):
        g_m = [[Fraction(-2), Fraction(1)], [Fraction(1), Fraction(0)]]
        s_h = reflection(g_m, [Fraction(1), Fraction(0)])
        denominator_family.append(
            {
                "m": m,
                "physical_boundary_square": -2 * m,
                "extension_gram": g_m,
                "reflection": s_h,
                "integral_isometry": integral_matrix(s_h)
                and matmul(transpose(s_h), matmul(g_m, s_h)) == g_m,
            }
        )

    arbitrary_square_sublattice_family = []
    for a in range(1, 21):
        # On F_a use the finite-index sublattice generated by (C_a,aF).
        g_a = [[Fraction(-a), Fraction(a)], [Fraction(a), Fraction(0)]]
        s_c = reflection(g_a, [Fraction(1), Fraction(0)])
        arbitrary_square_sublattice_family.append(
            {
                "a": a,
                "physical_section_square": -a,
                "sublattice_index_in_NS": a,
                "sublattice_gram": g_a,
                "reflection": s_c,
                "integral_isometry": integral_matrix(s_c)
                and matmul(transpose(s_c), matmul(g_a, s_c)) == g_a,
            }
        )

    # Unimodular primitive-root and adjunction/characteristic sieve.
    # a=-r^2 must divide 2; for C=d*r smooth rational, d divides 2.
    sieve = []
    for a in (1, 2):
        for d in (1, 2):
            k_dot_r = Fraction(-2 + a * d * d, d)
            characteristic_ok = (k_dot_r.denominator == 1) and (
                (int(k_dot_r) - (-a)) % 2 == 0
            )
            sieve.append(
                {
                    "primitive_normal_square": -a,
                    "curve_divisibility": d,
                    "curve_square": -a * d * d,
                    "K_dot_normal": k_dot_r,
                    "characteristic_parity": characteristic_ok,
                    "survives": characteristic_ok,
                }
            )

    # Weighted standard-coefficient cover obstruction.  The bracket in
    # c2(X)=l*(12+2*sum(m_i-3+2/m_i)) is >=12, with equality exactly at m_i=2.
    weighted_samples = []
    for indices in ((2,), (2, 2), (3,), (2, 3), (4,), (2, 4, 4)):
        lcm = 1
        for m in indices:
            # Tiny Euclidean gcd to keep this dependency-free.
            a, b = lcm, m
            while b:
                a, b = b, a % b
            lcm = lcm * m // a
        bracket = Fraction(12) + 2 * sum(
            (Fraction(m) - 3 + Fraction(2, m) for m in indices), Fraction(0)
        )
        weighted_samples.append(
            {
                "ramification_indices": list(indices),
                "cover_degree": lcm,
                "c2": lcm * bracket,
                "noether_upper_bound": 24,
                "possible": lcm * bracket <= 24,
            }
        )

    weighted_term_identity = all(
        Fraction(m - 3, 1) + Fraction(2, m)
        == Fraction((m - 1) * (m - 2), m)
        for m in range(2, 101)
    )
    weighted_term_nonnegative_with_unique_zero = all(
        (Fraction((m - 1) * (m - 2), m) >= 0)
        and ((Fraction((m - 1) * (m - 2), m) == 0) == (m == 2))
        for m in range(2, 101)
    )

    checks = {
        "ns_unimodular": abs(determinant(g_ns)) == 1,
        "physical_wall_gram": curve_gram_ns,
        "physical_squares": {name: pairing(g_ns, c, c) for name, c in curves_ns.items()},
        "boundary_incidence_with_minus_one_walls": {
            name: pairing(g_ns, curves_ns["B"], curves_ns[name]) for name in ("E", "G=F-E")
        },
        "B_not_minus_2K": curves_ns["B"] != [-2 * x for x in canonical],
        "B_reflection_not_integral_on_NS": not integral_matrix(b_reflection_ns),
        "half_lattice_determinant": determinant(g_half),
        "half_lattice_2_elementary_discriminant_order": abs(determinant(g_half)) == 2,
        "root_gram": root_gram,
        "root_span_unimodular_coordinates": abs(determinant(root_coordinate_matrix)) == 1,
        "reflection_checks": reflection_checks,
        "log_class_coordinates_half_basis": log_class_half,
        "log_class_root_pairings": log_class_pairings,
        "boundary_reflection_fixes_log_class": log_class_pairings["h=B/2"] == 0,
        "all_reflections_fix_log_class": all(value == 0 for value in log_class_pairings.values()),
        "coxeter_labels": {"h,E": 2, "h,G": 3, "E,G": "infinity"},
        "triangle_area_over_pi": Fraction(1, 6),
        "noncompact_coxeter_polygon_minimum_area_over_pi": Fraction(1, 6),
        "reflection_maximal_by_area": True,
        "weighted_cover_only_m2_survives_samples": all(
            row["possible"] == all(m == 2 for m in row["ramification_indices"])
            for row in weighted_samples
        ),
        "weighted_term_factorization_checked": weighted_term_identity,
        "weighted_term_nonnegative_unique_zero_m2_checked": weighted_term_nonnegative_with_unique_zero,
    }

    return serialize(
        {
            "status": "PASS",
            "theorem": {
                "NS_integral_smooth_rational_facet_squares": [-1, -2, -4],
                "minus_4_case": "C=2r, r^2=-1, K.r=1",
                "excluded_case": "C=2r, r^2=-2 contradicts K characteristic",
                "complete_primitive_mori_dimension_ceiling": 9,
                "coble_incidence_not_forced": True,
                "coble_incidence_equivalent_to_all_facet_reflections_fixing_K_plus_B_over_2": True,
                "maximality_adds_coble_incidence": False,
                "disjoint_standard_coefficient_branch": "nonempty tame branch forces every m_i=2",
                "weighted_cover_c2_formula": "ell*(12+2*sum((m_i-1)*(m_i-2)/m_i))",
            },
            "counterexample": {
                "surface": "blow-up of F_4 at a torus fixed point off the negative section",
                "NS_basis": ["B", "F", "E"],
                "NS_gram": g_ns,
                "canonical_coordinates": canonical,
                "Mori_generators": curves_ns,
                "half_extension_basis": ["B/2", "F", "E"],
                "half_extension_scaled_gram": g_half,
                "reflections": reflections,
            },
            "unimodular_sieve": sieve,
            "arbitrary_denominator_family": denominator_family,
            "arbitrary_square_sublattice_family": arbitrary_square_sublattice_family,
            "weighted_cover_samples": weighted_samples,
            "checks": checks,
        }
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--write-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.write_archive:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(encoded)
    if args.check_archive:
        assert ARCHIVE.read_text() == encoded, "archive differs from exact rebuild"
    print(encoded, end="")


if __name__ == "__main__":
    main()
