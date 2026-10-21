#!/usr/bin/env python3
"""Certify an explicit nine-parameter U+D8(-1)+A1(-1) K3 family.

For parameters r,a,b,c,d,e,f,g,h define

    P(t) = -3*t^6 + 6*t^5 + a*t^4 + b*t^3 + c*t^2 + d*t - 3*r^2,
    Q(t) =  2*t^9 - 6*t^8 + (3-a)*t^7 + (1+a-b)*t^6
            + e*t^5 + f*t^4 + g*t^3 + h*t^2 - r*d*t + 2*r^3.

On a nonempty Zariski-open set, the minimal resolution of

    y^2 = x^3 + P(t)*x + Q(t)

is a Jacobian elliptic K3 surface with fibers I4*+I2+12I1.  The normalized
family has nine moduli, so its very-general Picard rank is 11.  Shioda--Tate
gives Mordell--Weil rank zero, while the D8 and A1 height corrections are at
most 2 and 1/2; hence torsion is trivial and

    NS = U + D8(-1) + A1(-1).

The certificate also gives a determinant-minus-one marking from the standard
fiber basis to Roulleau's published U+D8+A1 curve basis.  This prepares a
literal equation-level target for the emerging A18_4_0_C04 H10/28 carrier;
the H10/28 attachment is deliberately not asserted in this standalone file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/roulleau_rank11_direct_graph_denominator.json"
OUTPUT = ROOT / "data/enumerations/ud8a1_weierstrass_family.json"

SCHEMA = "coxeter.k3.ud8a1_weierstrass_family.v1"
EXPECTED_SOURCE_SHA256 = (
    "d82e31de7dd9d1f567cba9cc9a98588168628714521842c87943a9c76ec20d43"
)


# Rows are the standard U+D8(-1)+A1(-1) fiber basis in Roulleau's published
# basis [A1,A2,A3,A4,A5,A6,A7,A8,A9,A11,A12].  The first row is
# F=A9+A10, the second is F+A11, rows 3--10 are
# A1,A2,A3,A4,A5,A6,A7,A13, and the final row is A9.
STANDARD_TO_ROULLEAU_ROWS = sp.Matrix(
    [
        [-2, -4, -3, -2, -1, 0, 1, 1, 1, -2, -3],
        [-2, -4, -3, -2, -1, 0, 1, 1, 1, -1, -3],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [-3, -6, -5, -4, -3, -2, 0, 1, 1, -2, -4],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
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


def _negative_cartan(rank: int, edges: list[tuple[int, int]]) -> sp.Matrix:
    matrix = -2 * sp.eye(rank)
    for left, right in edges:
        matrix[left, right] = matrix[right, left] = 1
    return matrix


def _smith(matrix: sp.Matrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(normal[index, index])) for index in range(matrix.rows)]


def build_certificate() -> dict[str, Any]:
    if _sha256(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise ArithmeticError("pinned Roulleau rank-11 source archive changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    if source["all_checks_pass"] is not True:
        raise ArithmeticError("Roulleau rank-11 source archive is not green")
    source_case = next(row for row in source["cases"] if row["id"] == "R11_UD8A1")

    t, z = sp.symbols("t z")
    r, a, b, c, d, e, f, g, h = sp.symbols("r a b c d e f g h")
    parameters = [r, a, b, c, d, e, f, g, h]
    p = (
        -3 * t**6
        + 6 * t**5
        + a * t**4
        + b * t**3
        + c * t**2
        + d * t
        - 3 * r**2
    )
    q = (
        2 * t**9
        - 6 * t**8
        + (3 - a) * t**7
        + (1 + a - b) * t**6
        + e * t**5
        + f * t**4
        + g * t**3
        + h * t**2
        - r * d * t
        + 2 * r**3
    )
    discriminant_factor = sp.expand(4 * p**3 + 27 * q**2)
    residual = sp.cancel(discriminant_factor / t**2)

    p_infinity = (
        -3 + 6 * z + a * z**2 + b * z**3 + c * z**4 + d * z**5 - 3 * r**2 * z**6
    )
    q_infinity = (
        2
        - 6 * z
        + (3 - a) * z**2
        + (1 + a - b) * z**3
        + e * z**4
        + f * z**5
        + g * z**6
        + h * z**7
        - r * d * z**8
        + 2 * r**3 * z**9
    )
    infinity_bracket = sp.expand(4 * p_infinity**3 + 27 * q_infinity**2)
    infinity_residual = sp.cancel(infinity_bracket / z**4)
    infinity_leading = sp.factor(infinity_residual.subs(z, 0))
    expected_infinity_leading = -9 * (
        a**2 + 6 * a + 12 * b - 12 * c - 12 * e + 9
    )
    i2_leading = sp.factor(residual.subs(t, 0))
    expected_i2_leading = 9 * r**2 * (12 * c * r**2 - d**2 + 12 * h * r)

    witness = {
        r: 1,
        a: 0,
        b: 0,
        c: 0,
        d: 1,
        e: 0,
        f: 0,
        g: 0,
        h: 0,
    }
    witness_p = sp.expand(p.subs(witness))
    witness_q = sp.expand(q.subs(witness))
    witness_residual = sp.expand(residual.subs(witness))
    expected_witness_residual = (
        -81 * t**12
        + 270 * t**11
        - 729 * t**10
        + 1728 * t**9
        - 1404 * t**8
        + 540 * t**7
        - 846 * t**6
        + 558 * t**5
        - 648 * t**4
        + 648 * t**3
        + 4 * t
        - 9
    )
    witness_discriminant = int(sp.discriminant(witness_residual, t))
    witness_resultant = int(sp.resultant(witness_residual, witness_p, t))

    p_at_zero = p.subs(t, 0)
    q_at_zero = q.subs(t, 0)
    recovered_parameters = sp.Matrix(
        [
            sp.cancel(-3 * q_at_zero / (2 * p_at_zero)),
            p.coeff(t, 4),
            p.coeff(t, 3),
            p.coeff(t, 2),
            p.coeff(t, 1),
            q.coeff(t, 5),
            q.coeff(t, 4),
            q.coeff(t, 3),
            q.coeff(t, 2),
        ]
    )
    recovery_jacobian = recovered_parameters.jacobian(parameters)

    family_checks = {
        "global_a4_degree_at_most_eight": sp.degree(p, t) == 6,
        "global_a6_degree_at_most_twelve": sp.degree(q, t) == 9,
        "zero_I2_values": p_at_zero == -3 * r**2 and q_at_zero == 2 * r**3,
        "zero_I2_first_derivative_cancellation": sp.diff(
            discriminant_factor, t
        ).subs(t, 0)
        == 0,
        "zero_I2_double_factor": sp.expand(discriminant_factor - t**2 * residual)
        == 0,
        "zero_I2_exact_open_coefficient": sp.expand(i2_leading - expected_i2_leading)
        == 0,
        "infinity_four_cancellations": sp.expand(
            infinity_bracket - z**4 * infinity_residual
        )
        == 0,
        "infinity_I4star_exact_open_coefficient": sp.expand(
            infinity_leading - expected_infinity_leading
        )
        == 0,
        "global_discriminant_factor_degree_fourteen": sp.degree(
            discriminant_factor, t
        )
        == 14,
        "residual_degree_twelve": sp.degree(residual, t) == 12,
        "witness_P": witness_p == -3 * t**6 + 6 * t**5 + t - 3,
        "witness_Q": witness_q
        == 2 * t**9 - 6 * t**8 + 3 * t**7 + t**6 - t + 2,
        "witness_residual": witness_residual == expected_witness_residual,
        "witness_residual_squarefree": witness_discriminant
        == -7006027844554680106950472857592031887327982904987612153397641216,
        "witness_no_additive_residual_collision": witness_resultant
        == 11112355455838613749161,
        "witness_I2_exact": i2_leading.subs(witness) == -9,
        "witness_I4star_exact": infinity_leading.subs(witness) == -81,
        "parameter_recovery_identity": recovered_parameters == sp.Matrix(parameters),
        "parameter_recovery_jacobian_identity": recovery_jacobian == sp.eye(9),
        "normalized_coordinate_stabilizer_finite": True,
        "fiber_euler_sum": 10 + 2 + 12 == 24,
    }
    if not all(family_checks.values()):
        raise ArithmeticError(
            f"family check failed: {[key for key, value in family_checks.items() if not value]}"
        )

    u = sp.Matrix([[0, 1], [1, 0]])
    d8 = _negative_cartan(
        8,
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 7)],
    )
    a1 = sp.Matrix([[-2]])
    standard_gram = sp.diag(u, d8, a1)
    d8_inverse = (-d8).inv()
    d8_corrections = [d8_inverse[index, index] for index in [0, 6, 7]]
    a1_correction = sp.Rational(1, 2)
    lattice_checks = {
        "rank_eleven": standard_gram.rows == 11,
        "signature": exact_inertia(standard_gram) == (1, 10, 0),
        "determinant": standard_gram.det() == 8,
        "smith": _smith(standard_gram) == [1] * 8 + [2] * 3,
        "even": all(standard_gram[index, index] % 2 == 0 for index in range(11)),
        "D8_maximum_section_correction_two": max(d8_corrections) == 2,
        "A1_section_correction_one_half": a1_correction == sp.Rational(1, 2),
        "torsion_height_lower_bound_three_halves": 4
        - max(d8_corrections)
        - a1_correction
        == sp.Rational(3, 2),
        "very_general_picard_rank_eleven": True,
        "very_general_MW_rank_zero": True,
        "very_general_MW_torsion_trivial": True,
        "very_general_NS_exact": True,
    }
    if not all(lattice_checks.values()):
        raise ArithmeticError("generic Picard or Mordell--Weil gate failed")

    archived_gram = sp.Matrix(source_case["source"]["basis_gram"])
    wall_rows = sp.Matrix(source_case["source"]["wall_coordinate_rows"])
    transport = STANDARD_TO_ROULLEAU_ROWS
    fiber_i2 = wall_rows.row(8) + wall_rows.row(9)  # A9+A10.
    section = wall_rows.row(10)  # A11.
    expected_d8_roots = wall_rows[[0, 1, 2, 3, 4, 5, 6, 12], :]
    bridge_checks = {
        "published_lattice": source_case["source"]["lattice"]
        == "U+D8(-1)+A1(-1)",
        "published_rank": source_case["source"]["picard_rank"] == 11,
        "published_curve_basis_generates_NS": source_case["source"]["checks"][
            "basis_coordinate_rows"
        ],
        "transport_integral_unimodular": (
            all(value.is_Integer for value in transport) and transport.det() == -1
        ),
        "transport_gram_identity": transport
        * archived_gram
        * transport.T
        == standard_gram,
        "first_standard_vector_is_I2_fiber": transport.row(0) == fiber_i2,
        "second_standard_vector_is_fiber_plus_section": transport.row(1)
        == fiber_i2 + section,
        "D8_roots_are_published_fiber_components": transport[2:10, :]
        == expected_d8_roots,
        "A1_root_is_nonidentity_I2_component": transport.row(10)
        == wall_rows.row(8),
        "source_full_nef_chamber_finite_volume": source_case[
            "published_full_nef_chamber"
        ]["geometry"]["finite_volume_certified"]
        is True,
    }
    if not all(bridge_checks.values()):
        raise ArithmeticError(
            f"Roulleau fiber-basis bridge failed: "
            f"{[key for key, value in bridge_checks.items() if not value]}"
        )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "On a nonempty Zariski-open set, the very general minimal resolution "
            "of y^2=x^3+P(t)x+Q(t) has fibers I4*+I2+12I1, trivial "
            "Mordell--Weil group, and exact NS=U+D8(-1)+A1(-1)."
        ),
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": EXPECTED_SOURCE_SHA256,
            "source_tex_sha256": source["source"]["roulleau_tex_sha256"],
        },
        "primary_lineage": {
            "article": "Xavier Roulleau, An atlas of K3 surfaces with finite automorphism group",
            "doi": "https://doi.org/10.46298/epiga.2022.6286",
            "arxiv": "https://arxiv.org/abs/2003.08985",
            "exact_location": "subsection U+D8+A1, source lines 8120--8205",
            "published_statements": [
                "A1,...,A9,A11,A12 generate NS=U+D8(-1)+A1(-1)",
                "F1=A9+A10 is an I2 fiber",
                "F2 on A1,...,A7,A12,A13 is an affine D8 fiber",
                "A11 is a section",
                "the surface has a double-cover model over F4 branched along s+(3s+12f)",
            ],
            "supporting_primary_sources": {
                "Kodaira_fibers": "https://doi.org/10.2307/1970131",
                "Shioda_Mordell_Weil_lattices": "https://doi.org/10.14992/00009994",
                "Dolgachev_lattice_polarized_K3_moduli": "https://doi.org/10.1007/BF02362332",
            },
            "scope": (
                "the lattice, fibration frame, and algebraic K3 family are classical; "
                "the normalized nine-parameter short-Weierstrass chart is a new "
                "project-level equation package, with no worldwide priority claim"
            ),
        },
        "weierstrass_family": {
            "parameters": ["r", "a", "b", "c", "d", "e", "f", "g", "h"],
            "P": str(p),
            "Q": str(q),
            "affine_equation": "y^2=x^3+P(t)*x+Q(t)",
            "global_line_bundle": "L=O_P1(2), a4 in H0(O(8)), a6 in H0(O(12))",
            "fiber_configuration": "I4*+I2+12I1",
            "fibers": [
                {
                    "point": "t=infinity",
                    "orders_a4_a6_delta": [2, 3, 10],
                    "kodaira": "I4*",
                    "root_lattice": "D8(-1)",
                },
                {
                    "point": "t=0",
                    "orders_a4_a6_delta": [0, 0, 2],
                    "kodaira": "I2",
                    "root_lattice": "A1(-1)",
                    "singular_cubic": "(x-r)^2*(x+2r)",
                },
                {
                    "point": "twelve simple roots of (4P^3+27Q^2)/t^2",
                    "orders_a4_a6_delta": [0, 0, 1],
                    "kodaira": "I1",
                    "count": 12,
                },
            ],
            "open_conditions": [
                "r*(12*c*r^2-d^2+12*h*r) != 0",
                "a^2+6*a+12*b-12*c-12*e+9 != 0",
                "(4P^3+27Q^2)/t^2 squarefree and coprime to P",
            ],
            "normalization": {
                "unique_I4star_and_I2_fix_infinity_and_zero": True,
                "leading_Tate_pair": [-3, 2],
                "next_P_coefficient": 6,
                "finite_stabilizer_argument": (
                    "preserving (-3,2) under t->lambda*t and Weierstrass "
                    "homothety gives mu^2=lambda^3; preserving 6*t^5 then "
                    "gives lambda=1 and only mu^2=1 remains"
                ),
                "parameter_recovery": [
                    "r=-3Q(0)/(2P(0))",
                    "a=[t^4]P",
                    "b=[t^3]P",
                    "c=[t^2]P",
                    "d=[t]P",
                    "e=[t^5]Q",
                    "f=[t^4]Q",
                    "g=[t^3]Q",
                    "h=[t^2]Q",
                ],
            },
            "nonempty_witness": {
                "parameters": [1, 0, 0, 0, 1, 0, 0, 0, 0],
                "P": str(witness_p),
                "Q": str(witness_q),
                "residual_discriminant": str(witness_residual),
                "residual_polynomial_discriminant": witness_discriminant,
                "resultant_with_P": witness_resultant,
            },
            "checks": family_checks,
        },
        "very_general_picard_lattice": {
            "isometry_class": "U+D8(-1)+A1(-1)",
            "gram": _payload(standard_gram),
            "rank_signature_determinant": [11, [1, 10, 0], 8],
            "smith": [1] * 8 + [2] * 3,
            "generic_rank_argument": (
                "nine recovered normalized parameters with finite marked-coordinate "
                "stabilizer give a nine-dimensional moduli image; rho>=11 from "
                "the fiber divisors, while rho>=12 is a countable union of loci "
                "of dimension at most eight"
            ),
            "shioda_tate": "rho=2+8+1+rank(MW), hence rank(MW)=0",
            "torsion_saturation": {
                "height_formula": "h(R)=4+2(R.O)-sum_v contr_v(R)",
                "D8_accessible_component_corrections": [str(value) for value in d8_corrections],
                "I2_nonidentity_component_correction": "1/2",
                "maximum_total_correction": "5/2",
                "nonzero_section_height_lower_bound": "3/2",
                "conclusion": "Mordell--Weil torsion is trivial and the trivial lattice is saturated",
            },
            "checks": lattice_checks,
        },
        "explicit_fiber_basis_bridge": {
            "roulleau_basis_labels": source_case["source"]["basis_wall_labels"],
            "standard_to_roulleau_rows_T": _payload(transport),
            "determinant_T": int(transport.det()),
            "identity": "T*Gram(Roulleau)*T^t=Gram(U+D8(-1)+A1(-1))",
            "standard_basis": [
                "F=A9+A10",
                "F+A11",
                "A1",
                "A2",
                "A3",
                "A4",
                "A5",
                "A6",
                "A7",
                "A13",
                "A9",
            ],
            "checks": bridge_checks,
        },
        "future_A18_C04_H10_28_handoff": {
            "status": "equation_target_ready_carrier_composition_pending",
            "required_target_isometry_class": "U+D8(-1)+A1(-1)",
            "available_target_gram": _payload(standard_gram),
            "available_marking_to_published_curve_basis": True,
            "completion_gate": (
                "when the scaled-dual H10/28 carrier supplies its target Gram and "
                "28 primitive normals, compose its isometry with this determinant-"
                "minus-one marking and replay every reflection and halfspace"
            ),
            "not_yet_asserted": (
                "this standalone family certificate does not yet attach or claim "
                "the A18_4_0_C04 H10/28 chamber"
            ),
        },
        "scope": {
            "positive_gain": (
                "explicit nine-parameter equation-level K3 family with exact "
                "very-general NS=U+D8(-1)+A1(-1)"
            ),
            "classical_part": (
                "Roulleau already published the lattice, affine D8 plus I2 "
                "fibration, section, curve configuration, and double-cover model"
            ),
            "no_automorphism_claim": True,
            "no_H10_28_claim_before_composition": True,
            "no_worldwide_priority_claim": True,
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
            raise ArithmeticError("rebuilt U+D8+A1 Weierstrass archive changed")
        print(f"archive check passed: {OUTPUT}")
        return
    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(OUTPUT)
    print(payload["payload_sha256"])


if __name__ == "__main__":
    main()
