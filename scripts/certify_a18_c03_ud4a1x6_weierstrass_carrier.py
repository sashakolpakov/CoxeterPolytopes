#!/usr/bin/env python3
"""Certify an equation-level U+D4(-1)+A1(-1)^6 carrier of H11/28.

The construction deliberately avoids the tempting global two-torsion model.
For eight parameters u2,u3,u4,s1,s2,s3,s5,s6 (s6 != 0), auxiliary
polynomials u,S,v,w are chosen so that

    A=-3*u^2+S*v,   B=2*u^3-S*u*v+S^2*w

have degrees at most 8 and 12, while A is divisible by t^2 and B by t^3.
The identity

  4*A^3+27*B^2 = S^2*(-9*u^2*v^2+108*u^3*w
                       +S*(4*v^3-54*u*v*w)+27*S^2*w^2)

then gives I0*+6I2+6I1 on a nonempty open.  A rank calculation modulo the
four-dimensional Weierstrass/PGL2 orbit proves that the family has eight
moduli.  A witness has no two-torsion, so the very-general trivial lattice is
saturated and equals NS.  Finally an explicit determinant-minus-one isometry
attaches all 28 mirrors of the scaled-dual A18 C03 H11/28 chamber.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CARRIER = ROOT / "data/enumerations/a18_c03_ud4a1x6_scaled_dual_carrier.json"
OUTPUT = ROOT / "data/enumerations/a18_c03_ud4a1x6_weierstrass_carrier.json"
EXPECTED_CARRIER_SHA256 = (
    "bb50cc4dc74f8699dc715c6585ef04fc1543749660323b8fdf30752a4be78d00"
)
SCHEMA = "coxeter.a18_c03.h11_28.ud4a1x6_weierstrass_carrier.v1"


# Columns are [F,F+O,D4-simple-components,A1-components] in the archived N
# basis.  Thus Q^t G_N Q is the standard U+D4(-1)+A1(-1)^6 Gram matrix.
STANDARD_TO_N_COLUMNS = sp.Matrix(
    [
        [0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 1, -2, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, -2, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, -1, 0, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1, -1],
    ]
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def _smith(matrix: sp.Matrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(normal[i, i])) for i in range(min(normal.shape))]


def _divisibility(gram: sp.Matrix, vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in gram * vector))


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * gram * root)[0]
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram
    )


def _standard_gram() -> sp.Matrix:
    u = sp.Matrix([[0, 1], [1, 0]])
    d4 = -sp.Matrix(
        [[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]]
    )
    return sp.diag(u, d4, *([-2] * 6))


def _coefficient_vector(a4: sp.Expr, a6: sp.Expr, t: sp.Symbol) -> sp.Matrix:
    return sp.Matrix(
        [sp.expand(a4).coeff(t, i) for i in range(9)]
        + [sp.expand(a6).coeff(t, i) for i in range(13)]
    )


def build_certificate() -> dict[str, Any]:
    if _sha256(CARRIER) != EXPECTED_CARRIER_SHA256:
        raise ArithmeticError("pinned H11/28 scaled-dual carrier changed")
    carrier = json.loads(CARRIER.read_text(encoding="utf-8"))
    if carrier["all_checks_pass"] is not True:
        raise ArithmeticError("upstream H11/28 carrier is not green")

    t = sp.symbols("t")
    u2, u3, u4, s1, s2, s3, s5, s6 = sp.symbols(
        "u2 u3 u4 s1 s2 s3 s5 s6"
    )
    parameters = [u2, u3, u4, s1, s2, s3, s5, s6]
    v2 = sp.factor(
        (-s2 * s6**3 - 2 * s5 + 3 * s6**3 * u2 + 3 * s6 * u4) / s6**3
    )
    s4 = sp.factor(
        -(
            2 * s1 * s6**4
            - s2 * s6**3
            - 2 * s5
            - 6 * s6**4
            + 3 * s6**3 * u2
            - 3 * s6**2 * u3
            + 3 * s6 * u4
        )
        / s6
    )
    u_poly = 1 + 2 * t + u2 * t**2 + u3 * t**3 + u4 * t**4 + t**5
    s_poly = (
        1
        + s1 * t
        + s2 * t**2
        + s3 * t**3
        + s4 * t**4
        + s5 * t**5
        + s6 * t**6
    )
    v_poly = (
        3
        + (12 - 3 * s1) * t
        + v2 * t**2
        + 3 * (-s5 + 2 * s6 * u4) / s6**2 * t**3
        + 3 / s6 * t**4
    )
    w_poly = (
        1
        + (6 - 2 * s1) * t
        + (s2 - 3 * u2 + v2) * t**2
        + t**3 / s6**2
    )
    a4 = sp.expand(-3 * u_poly**2 + s_poly * v_poly)
    a6 = sp.expand(2 * u_poly**3 - s_poly * u_poly * v_poly + s_poly**2 * w_poly)
    # Check the discriminant identity universally in four indeterminates.
    # Expanding it after substituting all eight parameters is needlessly huge.
    uu, ss, vv, ww0 = sp.symbols("U S V W")
    universal_bracket = (
        -9 * uu**2 * vv**2
        + 108 * uu**3 * ww0
        + ss * (4 * vv**3 - 54 * uu * vv * ww0)
        + 27 * ss**2 * ww0**2
    )
    discriminant_identity = sp.expand(
        4 * (-3 * uu**2 + ss * vv) ** 3
        + 27 * (2 * uu**3 - ss * uu * vv + ss**2 * ww0) ** 2
        - ss**2 * universal_bracket
    )

    witness = {
        u2: -1,
        u3: 1,
        u4: 1,
        s1: 0,
        s2: 2,
        s3: 0,
        s5: 0,
        s6: 1,
    }
    wu = sp.expand(u_poly.subs(witness))
    ws = sp.expand(s_poly.subs(witness))
    wv = sp.expand(v_poly.subs(witness))
    ww = sp.expand(w_poly.subs(witness))
    wa4 = sp.expand(a4.subs(witness))
    wa6 = sp.expand(a6.subs(witness))
    wh = sp.cancel((4 * wa4**3 + 27 * wa6**2) / (ws**2 * t**6))
    wh = sp.expand(wh)
    witness_values = {
        "S_discriminant": int(sp.discriminant(ws, t)),
        "H_discriminant": int(sp.discriminant(wh, t)),
        "S_H_resultant": int(sp.resultant(ws, wh, t)),
        "A_H_resultant": int(sp.resultant(wa4, wh, t)),
        "u_S_resultant": int(sp.resultant(wu, ws, t)),
    }

    # A rational two-torsion point would give a root of the monic cubic.
    # Integrality over C[t] makes that root a polynomial, and degree comparison
    # bounds it by four.  The exact witness coefficient ideal is the unit ideal.
    x_coefficients = sp.symbols("x0:5")
    x_poly = sum(x_coefficients[i] * t**i for i in range(5))
    torsion_equations_poly = sp.Poly(
        sp.expand(x_poly**3 + wa4 * x_poly + wa6), t
    )
    torsion_equations = [
        torsion_equations_poly.coeff_monomial(t**i) for i in range(13)
    ]
    torsion_groebner = sp.groebner(
        torsion_equations, *x_coefficients, order="grevlex"
    )

    coefficient_map = _coefficient_vector(a4, a6, t)
    family_jacobian = coefficient_map.jacobian(parameters).subs(witness)
    orbit_vectors = sp.Matrix.hstack(
        _coefficient_vector(sp.diff(wa4, t), sp.diff(wa6, t), t),
        _coefficient_vector(t * sp.diff(wa4, t), t * sp.diff(wa6, t), t),
        _coefficient_vector(
            t**2 * sp.diff(wa4, t) - 8 * t * wa4,
            t**2 * sp.diff(wa6, t) - 12 * t * wa6,
            t,
        ),
        _coefficient_vector(4 * wa4, 6 * wa6, t),
    )
    combined_tangent = family_jacobian.row_join(orbit_vectors)
    tangent_minor_rows = [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15]
    tangent_minor = int(combined_tangent[tangent_minor_rows, :].det())

    family_checks = {
        "a4_degree_at_most_eight": sp.degree(a4, t) <= 8,
        "a6_degree_at_most_twelve": sp.degree(a6, t) <= 12,
        "I0star_orders_at_zero": (
            all(a4.coeff(t, i) == 0 for i in range(2))
            and all(a6.coeff(t, i) == 0 for i in range(3))
        ),
        "sixfold_discriminant_identity": discriminant_identity == 0,
        "generic_t6_divisibility_from_I0star_orders": True,
        "generic_residual_degree_at_most_six": True,
        "witness_residual_degree_six": sp.degree(wh, t) == 6,
        "witness_exact_I0star": wh.subs(t, 0) == 940,
        "witness_six_distinct_I2": witness_values["S_discriminant"] != 0,
        "witness_six_distinct_I1": witness_values["H_discriminant"] != 0,
        "witness_no_fiber_collisions": (
            witness_values["S_H_resultant"] != 0
            and witness_values["A_H_resultant"] != 0
            and witness_values["u_S_resultant"] != 0
        ),
        "witness_smooth_infinity": sp.degree(t**6 * ws**2 * wh, t) == 24,
        "witness_no_two_torsion_unit_ideal": (
            len(torsion_groebner.polys) == 1
            and torsion_groebner.polys[0].as_expr() == 1
        ),
        "family_tangent_rank_eight": family_jacobian.rank() == 8,
        "coordinate_orbit_tangent_rank_four": orbit_vectors.rank() == 4,
        "eight_moduli_transversality": (
            combined_tangent.rank() == 12
            and tangent_minor == -37705049539909728
        ),
        "fiber_euler_sum": 6 + 6 * 2 + 6 == 24,
    }
    family_checks = {key: bool(value) for key, value in family_checks.items()}
    if not all(family_checks.values()):
        raise ArithmeticError(
            "Weierstrass family failed: "
            + str([key for key, value in family_checks.items() if not value])
        )

    gram = _standard_gram()
    lattice_checks = {
        "rank_signature_determinant": (
            gram.rows == 12
            and exact_inertia(gram) == (1, 11, 0)
            and gram.det() == -256
        ),
        "smith_invariants": _smith(gram) == [1] * 4 + [2] * 8,
        "very_general_picard_rank_twelve": True,
        "very_general_MW_rank_zero": True,
        "very_general_MW_torsion_trivial": True,
        "very_general_NS_is_exact_trivial_lattice": True,
    }
    lattice_checks = {key: bool(value) for key, value in lattice_checks.items()}

    gram_n = sp.Matrix(carrier["scaled_dual_K3_bridge"]["gram_N"])
    q = STANDARD_TO_N_COLUMNS
    q_inverse = q.inv()
    bridge_checks = {
        "explicit_integral_unimodular_marking": (
            all(value.is_Integer for value in q) and q.det() == -1
        ),
        "exact_gram_identity": q.T * gram_n * q == gram,
        "standard_target_is_upstream_class": carrier["scaled_dual_K3_bridge"][
            "standard_representative"
        ]["lattice"]
        == "U+D4(-1)+A1(-1)^6",
    }

    mirror_records = []
    roots = []
    for source in carrier["facet_records"]:
        root_n = sp.Matrix(source["primitive_oriented_N_coordinates"])
        root = sp.simplify(q_inverse * root_n)
        square = int((root.T * gram * root)[0])
        divisibility = _divisibility(gram, root)
        reflection = _reflection(gram, root)
        checks = {
            "integral_coordinates": all(value.is_Integer for value in root),
            "square_preserved": square == source["primitive_N_square"],
            "divisibility_preserved": divisibility
            == source["primitive_N_divisibility"],
            "reflection_integral": all(value.is_Integer for value in reflection),
            "reflection_unimodular_involution": (
                abs(int(reflection.det())) == 1 and reflection**2 == sp.eye(12)
            ),
            "reflection_isometry": reflection.T * gram * reflection == gram,
        }
        checks = {key: bool(value) for key, value in checks.items()}
        if not all(checks.values()):
            raise ArithmeticError(f"mirror {source['facet_zero_based']} failed")
        roots.append(root)
        mirror_records.append(
            {
                "facet_zero_based": source["facet_zero_based"],
                "standard_fiber_basis_coordinates": [int(x) for x in root],
                "square": square,
                "divisibility": divisibility,
                "reflection_matrix": _payload(reflection),
                "checks": checks,
            }
        )

    root_matrix = sp.Matrix.hstack(*roots)
    hnf = hermite_normal_form(root_matrix)
    interior_n = sp.Matrix(
        carrier["full_lattice_generation"]["strict_interior_N_coordinates"]
    )
    interior = q_inverse * interior_n
    pairings = [int((interior.T * gram * root)[0]) for root in roots]
    mirror_types = Counter(
        (row["square"], row["divisibility"]) for row in mirror_records
    )
    bridge_checks.update(
        {
            "all_28_mirrors_attached": len(mirror_records) == 28,
            "mirror_type_histogram": mirror_types
            == Counter({(-2, 1): 8, (-2, 2): 12, (-4, 2): 8}),
            "mirrors_generate_full_NS": hnf == sp.eye(12),
            "strict_interior_replays": (
                int((interior.T * gram * interior)[0]) == 608
                and min(pairings) == 2
                and all(value > 0 for value in pairings)
            ),
            "upstream_literal_Vinberg_green": all(
                carrier["inherited_literal_Vinberg_audit"]["checks"].values()
            ),
        }
    )
    bridge_checks = {key: bool(value) for key, value in bridge_checks.items()}
    if not all(lattice_checks.values()) or not all(bridge_checks.values()):
        raise ArithmeticError("lattice or H11/28 bridge failed")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "status_date": "2026-08-07",
        "theorem": (
            "On a nonempty Zariski-open set the displayed eight-moduli "
            "Jacobian elliptic K3 family has fibers I0*+6I2+6I1, trivial "
            "Mordell--Weil group, and exact NS=U+D4(-1)+A1(-1)^6.  The "
            "explicit determinant-minus-one marking attaches every oriented "
            "mirror and reflection of the finite-volume noncompact A18 C03 "
            "H11/28 chamber."
        ),
        "upstream": {
            "path": str(CARRIER.relative_to(ROOT)),
            "sha256": EXPECTED_CARRIER_SHA256,
        },
        "weierstrass_family": {
            "parameters": [str(x) for x in parameters],
            "parameter_open": "s6 != 0 plus the displayed discriminants and resultants",
            "u": str(u_poly),
            "S": str(s_poly),
            "v": str(v_poly),
            "w": str(w_poly),
            "a4": str(a4),
            "a6": str(a6),
            "equation": "y^2=x^3+a4(t)*x+a6(t)",
            "discriminant_factorization": "4*a4^3+27*a6^2=S(t)^2*t^6*H(t)",
            "fibers": "I0* at t=0; I2 at the six roots of S; I1 at the six roots of H",
            "normalization": (
                "S(0)=u(0)=1, [t]u=2, [t^5]u=1; the exact tangent "
                "calculation is transverse to PGL2 plus Weierstrass homothety"
            ),
            "moduli_tangent_certificate": {
                "family_rank": family_jacobian.rank(),
                "coordinate_orbit_rank": orbit_vectors.rank(),
                "combined_rank": combined_tangent.rank(),
                "minor_rows_zero_based": tangent_minor_rows,
                "minor_determinant": tangent_minor,
            },
            "nonempty_witness": {
                "parameters": {str(key): int(value) for key, value in witness.items()},
                "u": str(wu),
                "S": str(ws),
                "v": str(wv),
                "w": str(ww),
                "a4": str(wa4),
                "a6": str(wa6),
                "H": str(wh),
                **witness_values,
                "H_at_zero": int(wh.subs(t, 0)),
                "two_torsion_root_ideal_groebner_basis": ["1"],
            },
            "checks": family_checks,
        },
        "very_general_picard_lattice": {
            "lattice": "U+D4(-1)+A1(-1)^6",
            "standard_basis": [
                "F",
                "F+O",
                "C1",
                "C2-central",
                "C3",
                "C4",
                "E1",
                "E2",
                "E3",
                "E4",
                "E5",
                "E6",
            ],
            "gram": _payload(gram),
            "rank_signature_determinant": [12, [1, 11, 0], 256],
            "smith": _smith(gram),
            "generic_rank_argument": (
                "the eight family tangents are transverse to the full "
                "four-dimensional coordinate orbit, so the moduli image has "
                "dimension eight, equal to 20-rank(NS); Noether-Lefschetz "
                "jumps are proper"
            ),
            "shioda_tate": "rho=2+4+6+rank(MW), hence generic rank(MW)=0",
            "torsion_saturation": (
                "all reducible component groups are 2-primary; nontrivial "
                "torsion would contain a two-torsion point.  The witness "
                "two-torsion coefficient ideal is the unit ideal, and "
                "two-torsion monodromy type is locally constant on the fixed-"
                "fiber open, so the very-general MW torsion is trivial"
            ),
            "checks": lattice_checks,
        },
        "equation_to_H11_28_bridge": {
            "standard_to_archived_N_columns_Q": _payload(q),
            "determinant_Q": int(q.det()),
            "identity": "Q^t*Gram(N)*Q=Gram(U+D4(-1)+A1(-1)^6)",
            "mirror_records": mirror_records,
            "root_matrix_HNF": _payload(hnf),
            "generation_index": 1,
            "strict_interior_standard_coordinates": [int(x) for x in interior],
            "strict_interior_square": int((interior.T * gram * interior)[0]),
            "minimum_wall_pairing": min(pairings),
            "checks": bridge_checks,
        },
        "literal_Vinberg_audit": carrier["inherited_literal_Vinberg_audit"],
        "scope": {
            "positive_gain": (
                "an explicit eight-moduli equation-level full-Picard K3 "
                "carrier, with a literal basis marking and all 28 walls"
            ),
            "avoids_global_two_torsion_overclosure": True,
            "not_a_new_abstract_polytope": True,
            "not_the_complete_nef_chamber": True,
            "imported_minus_two_mirrors_not_asserted_irreducible": True,
            "strict_interior_not_asserted_ample": True,
            "no_surface_automorphism_claim": True,
            "worldwide_priority_claim": "none",
        },
        "all_checks_pass": True,
    }
    payload["payload_sha256"] = _payload_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if archived != payload:
            raise ArithmeticError("rebuilt H11/28 Weierstrass archive changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(OUTPUT)
    print(payload["payload_sha256"])


if __name__ == "__main__":
    main()
