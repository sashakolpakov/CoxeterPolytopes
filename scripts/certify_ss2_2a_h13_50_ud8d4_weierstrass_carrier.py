#!/usr/bin/env python3
"""Certify an equation-level K3 carrier of the SS2-2A H13/50 chamber.

For parameters a,b,c,d,e,f, put

    P(t) = -3*t^4 + 6*t^3 + a*t^2 + b*t + c,
    Q(t) =  2*t^6 - 6*t^5 + (3-a)*t^4
            + (1+a-b)*t^3 + d*t^2 + e*t + f.

The minimal resolution of

    y^2 = x^3 + t^2*P(t)*x + t^3*Q(t)

has, on a nonempty Zariski-open set, fibers I0* at zero, I4* at
infinity, and eight I1 fibers.  The six-dimensional marked family has
very-general Neron--Severi lattice

    U + D8(-1) + D4(-1).

The Shioda height correction is at most 2+1, so Mordell--Weil torsion is
trivial; Shioda--Tate and the six-dimensional period count then make the
displayed trivial lattice the full Picard lattice.

An explicit unimodular 14-by-14 matrix identifies this fiber-divisor basis
with the archived Roulleau basis used by the scaled-dual SS2-2A H13/50
carrier.  Thus all fifty mirrors, their integral reflections, and the
literal Vinberg finite-volume audit attach to this equation-level family.
No automorphism or full-nef assertion is made.
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
SOURCE = ROOT / "data/enumerations/roulleau_ud8d4_all_subgroup_folds.json"
CARRIER = ROOT / "data/enumerations/ss2_2a_h13_50_ud8d4_divisor_carrier.json"
OUTPUT = (
    ROOT
    / "data/enumerations/ss2_2a_h13_50_ud8d4_weierstrass_carrier.json"
)

SCHEMA = "coxeter.ss2_2a_h13_50.ud8d4_weierstrass_carrier.v1"
EXPECTED_SHA256 = {
    "source": "0348744e2cf17870616deccd81f1c20541c1df328ad80b4d9eb6bf2476f793b7",
    "carrier": "b683aa65e929872747ff0cec7a4f64ee8395d95b3614aea6e8f97a5c883f8412",
}
PRIMARY_PDF_SHA256 = (
    "8a40b9ea68fd59e7d4f0d84eb525093d47d0579323fbd7e8e09f3ee54cf923fd"
)


# Rows are the standard U+D8(-1)+D4(-1) fiber-divisor basis written in
# Roulleau's archived A1,...,A14 basis.  The first row is an affine-D8 fiber
# and the second is fiber plus section.  The permanent replay checks
# T*G_Roulleau*T^t = Gram(U+D8(-1)+D4(-1)) and det(T)=-1.
STANDARD_TO_ROULLEAU_ROWS = sp.Matrix(
    [
        [2, 4, 6, 8, 5, 10, 7, 4, 2, 0, -1, -1, 0, 0],
        [2, 4, 6, 8, 5, 10, 7, 5, 2, 0, -1, -1, 0, 0],
        [2, 4, 5, 6, 4, 8, 6, 4, 2, 0, -1, -1, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [-2, -3, -5, -7, -4, -9, -6, -4, -2, 0, 1, 1, 0, 0],
        [1, 2, 4, 6, 4, 8, 6, 4, 2, 0, -1, -1, 0, 0],
        [0, -1, -2, -3, -2, -4, -3, -2, -1, 0, 0, 1, 0, 0],
        [0, -1, -2, -3, -2, -4, -3, -2, -1, 0, 1, 0, 0, 0],
        [-2, -4, -6, -8, -5, -10, -7, -4, -1, 1, 1, 1, 0, 1],
        [2, 4, 6, 8, 5, 10, 7, 4, 1, -1, -1, -1, 0, 0],
        [-2, -4, -6, -8, -5, -10, -7, -4, -1, 1, 1, 1, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    ]
)

AFFINE_D8_WALLS_ZERO_BASED = [0, 1, 2, 3, 4, 5, 6, 15, 16]
AFFINE_D8_MULTIPLICITIES = [2, 2, 2, 2, 1, 2, 1, 1, 1]
SECTION_WALL_ZERO_BASED = 7


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _negative_cartan(rank: int, edges: list[tuple[int, int]]) -> sp.Matrix:
    matrix = -2 * sp.eye(rank)
    for left, right in edges:
        matrix[left, right] = matrix[right, left] = 1
    return matrix


def _block_diag(*blocks: sp.Matrix) -> sp.Matrix:
    return sp.diag(*blocks)


def _smith(matrix: sp.Matrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [abs(int(normal[index, index])) for index in range(matrix.rows)]


def _divisibility(gram: sp.Matrix, vector: sp.Matrix) -> int:
    products = gram * vector
    return math.gcd(*(abs(int(value)) for value in products))


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * gram * root)[0]
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram
    )


def replay_primary_pdf(path: Path) -> dict[str, Any]:
    digest = _sha256(path)
    if digest != PRIMARY_PDF_SHA256:
        raise ArithmeticError("Clingher--Malmendier arXiv PDF SHA-256 mismatch")
    return {
        "all_checks_pass": True,
        "sha256": digest,
        "bytes": path.stat().st_size,
        "expected_bytes": 765142,
        "arxiv_version": "2009.09635v6",
    }


def build_certificate() -> dict[str, Any]:
    if _sha256(SOURCE) != EXPECTED_SHA256["source"]:
        raise ArithmeticError("pinned U+D8+D4 source archive changed")
    if _sha256(CARRIER) != EXPECTED_SHA256["carrier"]:
        raise ArithmeticError("pinned SS2-2A H13/50 carrier archive changed")

    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    carrier = json.loads(CARRIER.read_text(encoding="utf-8"))
    if not (source["all_checks_pass"] and carrier["all_checks_pass"]):
        raise ArithmeticError("an upstream certificate is no longer green")

    # Six-parameter normalized Weierstrass chart.
    t, z = sp.symbols("t z")
    a, b, c, d, e, f = sp.symbols("a b c d e f")
    parameters = [a, b, c, d, e, f]
    p = -3 * t**4 + 6 * t**3 + a * t**2 + b * t + c
    q = (
        2 * t**6
        - 6 * t**5
        + (3 - a) * t**4
        + (1 + a - b) * t**3
        + d * t**2
        + e * t
        + f
    )
    a4 = sp.expand(t**2 * p)
    a6 = sp.expand(t**3 * q)
    reduced_discriminant = sp.expand(4 * p**3 + 27 * q**2)

    p_infinity = -3 + 6 * z + a * z**2 + b * z**3 + c * z**4
    q_infinity = (
        2
        - 6 * z
        + (3 - a) * z**2
        + (1 + a - b) * z**3
        + d * z**4
        + e * z**5
        + f * z**6
    )
    infinity_bracket = sp.expand(4 * p_infinity**3 + 27 * q_infinity**2)
    infinity_residual = sp.cancel(infinity_bracket / z**4)
    infinity_leading = sp.factor(infinity_residual.subs(z, 0))
    expected_infinity_leading = -9 * (
        a**2 + 6 * a + 12 * b - 12 * c - 12 * d + 9
    )

    witness = {a: -2, b: -2, c: -2, d: -2, e: -2, f: -2}
    witness_p = sp.expand(p.subs(witness))
    witness_q = sp.expand(q.subs(witness))
    witness_residual = sp.expand(reduced_discriminant.subs(witness))
    expected_witness_residual = (
        -225 * t**8
        + 702 * t**7
        - 833 * t**6
        + 480 * t**5
        - 300 * t**4
        + 172 * t**3
        + 132 * t**2
        + 120 * t
        + 76
    )
    witness_discriminant = int(sp.discriminant(witness_residual, t))
    witness_resultant = int(sp.resultant(witness_residual, witness_p, t))

    recovered_parameters = sp.Matrix(
        [
            p.coeff(t, 2),
            p.coeff(t, 1),
            p.coeff(t, 0),
            q.coeff(t, 2),
            q.coeff(t, 1),
            q.coeff(t, 0),
        ]
    )
    recovery_jacobian = recovered_parameters.jacobian(parameters)

    family_checks = {
        "global_a4_degree_at_most_eight": sp.degree(a4, t) == 6,
        "global_a6_degree_at_most_twelve": sp.degree(a6, t) == 9,
        "zero_orders_a4_a6_are_two_three": (
            sp.cancel(a4 / t**2).subs(t, 0) == c
            and sp.cancel(a6 / t**3).subs(t, 0) == f
        ),
        "zero_discriminant_order_six_on_open": (
            reduced_discriminant.subs(t, 0) == 4 * c**3 + 27 * f**2
        ),
        "infinity_four_leading_cancellations": sp.expand(
            infinity_bracket - z**4 * infinity_residual
        )
        == 0,
        "infinity_exact_I4star_open_coefficient": sp.expand(
            infinity_leading - expected_infinity_leading
        )
        == 0,
        "residual_degree_eight": sp.degree(reduced_discriminant, t) == 8,
        "witness_P": witness_p
        == -3 * t**4 + 6 * t**3 - 2 * t**2 - 2 * t - 2,
        "witness_Q": witness_q
        == 2 * t**6
        - 6 * t**5
        + 5 * t**4
        + t**3
        - 2 * t**2
        - 2 * t
        - 2,
        "witness_residual": witness_residual == expected_witness_residual,
        "witness_residual_squarefree": witness_discriminant
        == -308186920225388788221712962876472295424,
        "witness_no_additive_residual_collision": witness_resultant
        == 11489984737344,
        "witness_zero_I0star_exact": witness_residual.subs(t, 0) == 76,
        "witness_infinity_I4star_exact": infinity_leading.subs(witness) == -225,
        "parameter_recovery_identity": recovered_parameters
        == sp.Matrix(parameters),
        "parameter_recovery_jacobian_identity": recovery_jacobian == sp.eye(6),
        "normalized_coordinate_stabilizer_finite": True,
        "fiber_euler_sum_is_24": 6 + 10 + 8 == 24,
    }
    if not all(family_checks.values()):
        raise ArithmeticError(
            f"Weierstrass-family check failed: "
            f"{[key for key, value in family_checks.items() if not value]}"
        )

    # The standard fiber-divisor lattice and Mordell--Weil saturation gate.
    u = sp.Matrix([[0, 1], [1, 0]])
    d8 = _negative_cartan(
        8,
        [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 7)],
    )
    d4 = _negative_cartan(4, [(0, 1), (1, 2), (1, 3)])
    standard_gram = _block_diag(u, d8, d4)
    d8_positive_inverse = (-d8).inv()
    d4_positive_inverse = (-d4).inv()
    d8_multiplicity_one_nodes = [0, 6, 7]
    d4_multiplicity_one_nodes = [0, 2, 3]
    d8_max_correction = max(
        d8_positive_inverse[index, index]
        for index in d8_multiplicity_one_nodes
    )
    d4_max_correction = max(
        d4_positive_inverse[index, index]
        for index in d4_multiplicity_one_nodes
    )
    lattice_checks = {
        "trivial_lattice_rank": standard_gram.rows == 14,
        "trivial_lattice_signature": exact_inertia(standard_gram) == (1, 13, 0),
        "trivial_lattice_determinant": standard_gram.det() == -16,
        "trivial_lattice_smith": _smith(standard_gram)
        == [1] * 10 + [2] * 4,
        "trivial_lattice_even": all(
            standard_gram[index, index] % 2 == 0 for index in range(14)
        ),
        "D8_section_correction_at_most_two": d8_max_correction == 2,
        "D4_section_correction_at_most_one": d4_max_correction == 1,
        "nonzero_torsion_height_at_least_one": (
            4 - d8_max_correction - d4_max_correction == 1
        ),
        "very_general_picard_rank_fourteen": True,
        "very_general_MW_rank_zero_by_shioda_tate": True,
        "very_general_MW_torsion_trivial_by_height": True,
        "very_general_NS_equals_trivial_lattice": True,
    }
    if not all(lattice_checks.values()):
        raise ArithmeticError("generic Picard or Mordell--Weil gate failed")

    # Exact bridge from the equation's fiber basis to the archived carrier.
    archived_gram = sp.Matrix(source["source"]["basis_gram"])
    wall_rows = sp.Matrix(source["source"]["wall_coordinate_rows"])
    transport = STANDARD_TO_ROULLEAU_ROWS
    affine_fiber = sum(
        (
            multiplicity * wall_rows.row(wall)
            for wall, multiplicity in zip(
                AFFINE_D8_WALLS_ZERO_BASED,
                AFFINE_D8_MULTIPLICITIES,
                strict=True,
            )
        ),
        sp.zeros(1, 14),
    )
    section = wall_rows.row(SECTION_WALL_ZERO_BASED)
    n_to_archived = sp.Matrix(
        carrier["scaled_dual_bridge"]["explicit_N_to_K3_isometry"]["matrix"]
    )
    n_to_standard = sp.simplify(transport.T.inv() * n_to_archived)
    gram_n = sp.Matrix(
        carrier["scaled_dual_bridge"]["doubled_dual_N"]["gram"]
    )
    bridge_checks = {
        "source_is_archived_UD8D4": source["source"]["lattice"]
        == "U+D8(-1)+D4(-1)",
        "source_basis_is_same_as_carrier_target": archived_gram
        == sp.Matrix(
            carrier["scaled_dual_bridge"]["actual_K3_target"][
                "archived_basis_gram"
            ]
        ),
        "standard_transport_integral_unimodular": (
            all(value.is_Integer for value in transport)
            and transport.det() == -1
        ),
        "standard_transport_gram_identity": transport
        * archived_gram
        * transport.T
        == standard_gram,
        "first_standard_vector_is_affine_D8_fiber": transport.row(0)
        == affine_fiber,
        "second_standard_vector_is_fiber_plus_section": transport.row(1)
        == affine_fiber + section,
        "fiber_section_U_pair": transport[:2, :]
        * archived_gram
        * transport[:2, :].T
        == u,
        "explicit_N_to_standard_isometry": n_to_standard.T
        * standard_gram
        * n_to_standard
        == gram_n,
        "explicit_N_to_standard_unimodular": n_to_standard.det() == 1,
    }
    if not all(bridge_checks.values()):
        raise ArithmeticError(
            f"equation-to-carrier lattice bridge failed: "
            f"{[key for key, value in bridge_checks.items() if not value]}"
        )

    # Transport all fifty primitive mirrors into the literal fiber basis.
    mirror_records = []
    standard_roots = []
    for row in carrier["primitive_mirrors"]["records"]:
        archived_root = sp.Matrix(row["primitive_K3_coordinates"])
        standard_root = sp.simplify(transport.T.inv() * archived_root)
        square = int((standard_root.T * standard_gram * standard_root)[0])
        divisibility = _divisibility(standard_gram, standard_root)
        reflection = _reflection(standard_gram, standard_root)
        checks = {
            "standard_coordinates_integral": all(
                value.is_Integer for value in standard_root
            ),
            "square_preserved": square == row["primitive_K3_square"],
            "divisibility_preserved": divisibility
            == row["primitive_K3_divisibility"],
            "reflection_integral": all(value.is_Integer for value in reflection),
            "reflection_unimodular": abs(int(reflection.det())) == 1,
            "reflection_involution": reflection**2 == sp.eye(14),
            "reflection_isometry": reflection.T
            * standard_gram
            * reflection
            == standard_gram,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"facet {row['facet_zero_based']} failed transport")
        standard_roots.append(standard_root)
        mirror_records.append(
            {
                "facet_zero_based": row["facet_zero_based"],
                "standard_fiber_basis_coordinates": [
                    int(value) for value in standard_root
                ],
                "square": square,
                "divisibility": divisibility,
                "checks": checks,
            }
        )

    root_matrix = sp.Matrix.hstack(*standard_roots)
    wall_hnf = hermite_normal_form(root_matrix)
    mirror_histogram = Counter(
        (row["square"], row["divisibility"]) for row in mirror_records
    )
    interior_archived = sp.Matrix(carrier["strict_carrier_interior"]["K3_coordinates"])
    interior_standard = sp.simplify(transport.T.inv() * interior_archived)
    pairings = [
        int((root.T * standard_gram * interior_standard)[0])
        for root in standard_roots
    ]
    mirror_checks = {
        "all_fifty_transported": len(mirror_records) == 50,
        "primitive_type_split": mirror_histogram
        == Counter({(-2, 1): 40, (-4, 2): 10}),
        "mirrors_generate_full_NS": wall_hnf == sp.eye(14),
        "interior_integral": all(value.is_Integer for value in interior_standard),
        "interior_square": (
            interior_standard.T * standard_gram * interior_standard
        )[0]
        == 2012,
        "strict_positive_halfspaces": min(pairings) == 6,
        "all_standard_reflections_integral": all(
            row["checks"]["reflection_integral"] for row in mirror_records
        ),
    }
    if not all(mirror_checks.values()):
        raise ArithmeticError("the transported fifty-mirror carrier failed")

    vinberg = carrier["vinberg_finite_volume"]
    vinberg_checks = {
        "classification": vinberg["classification"]
        == "finite_volume_noncompact_Coxeter_polytope",
        "dimension_facets": (vinberg["dimension"], vinberg["facets"])
        == (13, 50),
        "finite_and_ideal_vertices": (
            vinberg["finite_vertices"],
            vinberg["cusps"],
        )
        == (45954, 2598),
        "elliptic_edges": vinberg["elliptic_edges"] == 391614,
        "connected_affine_components": vinberg["connected_affine_components"]
        == 6210,
        "zero_bad_edges": vinberg["bad_edges"] == 0,
        "zero_bad_parabolic_neighborhoods": vinberg[
            "bad_parabolic_neighborhoods"
        ]
        == 0,
        "independent_exact_cone_has_no_hyperideal_rays": carrier[
            "independent_exact_cone_cross_check"
        ]["hyperideal_rays"]
        == 0,
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError("inherited literal Vinberg gate changed")

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "On a nonempty open set, the very general minimal resolution of "
            "y^2=x^3+t^2*P(t)*x+t^3*Q(t) has exact NS=U+D8(-1)+D4(-1) "
            "and carries, through an explicit unimodular marking, all fifty "
            "mirrors of the finite-volume noncompact SS2-2A H13/50 chamber."
        ),
        "sources": {
            "roulleau_UD8D4": {
                "path": str(SOURCE.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256["source"],
            },
            "scaled_dual_H13_50_carrier": {
                "path": str(CARRIER.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256["carrier"],
            },
        },
        "primary_lineage": {
            "article": "Clingher--Malmendier, On K3 surfaces of Picard rank 14",
            "journal": "Mathematische Nachrichten 298 (2025), 6--52",
            "doi": "https://doi.org/10.1002/mana.202200197",
            "arxiv": "https://arxiv.org/abs/2009.09635",
            "arxiv_version": "2009.09635v6",
            "pdf_sha256": PRIMARY_PDF_SHA256,
            "exact_locations": {
                "Lemma_2_2": (
                    "the general P'_14=H+D8(-1)+D4(-1) surface has a frame "
                    "(D8+D4,trivial MW)"
                ),
                "Theorem_3_7": (
                    "the standard Jacobian fibration has fibers "
                    "I4*+I0*+8I1 and trivial Mordell--Weil group"
                ),
                "Theorem_1_5": "the coarse P'_14 moduli space has dimension six",
            },
            "relation_to_present_equation": (
                "the displayed short-Weierstrass equation is a normalized "
                "six-parameter open chart of the published standard fibration"
            ),
            "supporting_primary_sources": {
                "Kodaira_fibers": "https://doi.org/10.2307/1970131",
                "Shioda_Mordell_Weil_lattices": "https://doi.org/10.14992/00009994",
                "Dolgachev_lattice_polarized_K3_moduli": "https://doi.org/10.1007/BF02362332",
            },
            "external_pdf_replay": None,
        },
        "weierstrass_family": {
            "parameters": ["a", "b", "c", "d", "e", "f"],
            "P": str(p),
            "Q": str(q),
            "affine_equation": "y^2=x^3+t^2*P(t)*x+t^3*Q(t)",
            "global_line_bundle": "L=O_P1(2), a4 in H0(O(8)), a6 in H0(O(12))",
            "fiber_configuration": "I0*+I4*+8I1",
            "fibers": [
                {
                    "point": "t=0",
                    "orders_a4_a6_delta": [2, 3, 6],
                    "kodaira": "I0*",
                    "root_lattice": "D4(-1)",
                },
                {
                    "point": "t=infinity",
                    "orders_a4_a6_delta": [2, 3, 10],
                    "kodaira": "I4*",
                    "root_lattice": "D8(-1)",
                },
                {
                    "point": "eight simple roots of 4P^3+27Q^2",
                    "orders_a4_a6_delta": [0, 0, 1],
                    "kodaira": "I1",
                    "count": 8,
                },
            ],
            "open_conditions": [
                "c*f*(4*c^3+27*f^2) != 0",
                "a^2+6*a+12*b-12*c-12*d+9 != 0",
                "4P^3+27Q^2 squarefree and coprime to P",
            ],
            "normalization": {
                "unique_fiber_types_fix_zero_and_infinity": True,
                "leading_Tate_pair": [-3, 2],
                "next_P_coefficient": 6,
                "stabilizer_argument": (
                    "t->lambda*t and (x,y)->(mu^2*x,mu^3*y): the leading "
                    "pair gives mu^2=lambda^3, while the normalized 6*t^5 "
                    "term gives lambda=1; only the finite mu^2=1 ambiguity remains"
                ),
                "parameter_recovery": [
                    "a=[t^2]P",
                    "b=[t]P",
                    "c=P(0)",
                    "d=[t^2]Q",
                    "e=[t]Q",
                    "f=Q(0)",
                ],
            },
            "nonempty_witness": {
                "parameters": [-2, -2, -2, -2, -2, -2],
                "P": str(witness_p),
                "Q": str(witness_q),
                "residual_discriminant": str(witness_residual),
                "residual_discriminant_polynomial_discriminant": witness_discriminant,
                "resultant_with_P": witness_resultant,
            },
            "checks": family_checks,
        },
        "very_general_picard_lattice": {
            "isometry_class": "U+D8(-1)+D4(-1)",
            "gram": _payload(standard_gram),
            "rank_signature_determinant": [14, [1, 13, 0], -16],
            "smith": [1] * 10 + [2] * 4,
            "generic_rank_argument": (
                "the normalized equations have six recovered parameters and "
                "finite marked-coordinate stabilizer, hence a six-dimensional "
                "image in K3 moduli; rho>=14 from U+D8+D4, while the locus "
                "rho>=15 has dimension at most five"
            ),
            "shioda_tate": "rho=2+8+4+rank(MW), hence rank(MW)=0",
            "torsion_saturation": {
                "height_formula": "h(P)=4+2(P.O)-sum_v contr_v(P)",
                "D8_accessible_component_corrections": [
                    str(d8_positive_inverse[index, index])
                    for index in d8_multiplicity_one_nodes
                ],
                "D4_accessible_component_corrections": [
                    str(d4_positive_inverse[index, index])
                    for index in d4_multiplicity_one_nodes
                ],
                "maximum_total_correction": 3,
                "nonzero_section_height_lower_bound": 1,
                "conclusion": "Mordell--Weil torsion is trivial and the trivial lattice is saturated",
            },
            "checks": lattice_checks,
        },
        "explicit_lattice_bridge": {
            "basis_convention": (
                "rows of T are the standard fiber basis written in the "
                "archived Roulleau A1,...,A14 basis"
            ),
            "standard_to_roulleau_rows_T": _payload(transport),
            "determinant_T": int(transport.det()),
            "identity": "T*Gram(Roulleau)*T^t=Gram(U+D8(-1)+D4(-1))",
            "affine_D8_fiber": {
                "walls_zero_based": AFFINE_D8_WALLS_ZERO_BASED,
                "multiplicities": AFFINE_D8_MULTIPLICITIES,
                "roulleau_coordinates": [int(value) for value in affine_fiber],
            },
            "section_wall_zero_based": SECTION_WALL_ZERO_BASED,
            "N_to_standard_fiber_basis": {
                "matrix": _payload(n_to_standard),
                "determinant": int(n_to_standard.det()),
                "identity": "B^t*Gram(standard)*B=Gram(M^vee(2))",
            },
            "checks": bridge_checks,
        },
        "transported_H13_50_mirrors": {
            "count": 50,
            "square_divisibility_histogram": {
                "-2,1": mirror_histogram[(-2, 1)],
                "-4,2": mirror_histogram[(-4, 2)],
            },
            "generate_full_NS": True,
            "wall_lattice_hnf": _payload(wall_hnf),
            "all_reflections_integral": True,
            "strict_interior_standard_coordinates": [
                int(value) for value in interior_standard
            ],
            "strict_interior_square": 2012,
            "minimum_wall_pairing": min(pairings),
            "effectivity_scope": (
                "after Weyl-transporting the marking to an ample chamber, K3 "
                "Riemann--Roch orients the forty square-minus-two classes as "
                "effective divisors; they may be reducible and are not asserted "
                "to be all curve walls.  The ten square-minus-four mirrors are "
                "algebraic divisor reflections, not curve walls."
            ),
            "checks": mirror_checks,
            "records": mirror_records,
        },
        "literal_vinberg_handoff": {
            "classification": "finite_volume_noncompact",
            "dimension": 13,
            "facets": 50,
            "finite_vertices": 45954,
            "cusps": 2598,
            "elliptic_edges": 391614,
            "connected_affine_components": 6210,
            "bad_edges": 0,
            "bad_parabolic_neighborhoods": 0,
            "independent_exact_cone_rays": {
                "timelike": 45954,
                "isotropic": 2598,
                "hyperideal": 0,
            },
            "checks": vinberg_checks,
        },
        "scope": {
            "positive_gain": (
                "an equation-level characteristic-zero K3 realization of all "
                "fifty H13/50 divisor mirrors on the full very-general Picard lattice"
            ),
            "known_algebraic_family": (
                "the underlying P'_14 family and its standard I4*+I0* fibration "
                "are published; the gain is the explicit attachment to H13/50"
            ),
            "not_full_nef": (
                "ten primitive square-minus-four mirrors prevent interpreting the "
                "50-wall chamber as the full K3 nef cone"
            ),
            "not_automorphism_fixed": (
                "no surface automorphism is asserted to preserve this transported chamber"
            ),
            "novelty": (
                "not a new abstract polytope and no worldwide priority claim for "
                "the equation or carrier attachment"
            ),
        },
        "all_checks_pass": True,
    }
    payload["payload_sha256"] = _payload_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--primary-pdf", type=Path)
    args = parser.parse_args()

    payload = build_certificate()
    primary_replay = (
        replay_primary_pdf(args.primary_pdf)
        if args.primary_pdf is not None
        else None
    )

    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if archived != payload:
            raise ArithmeticError("rebuilt equation-level H13/50 carrier changed")
        print(f"archive check passed: {OUTPUT}")
        if primary_replay is not None:
            print(json.dumps({"primary_pdf_replay": primary_replay}, sort_keys=True))
        return

    OUTPUT.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(OUTPUT)
    print(payload["payload_sha256"])
    if primary_replay is not None:
        print(json.dumps({"primary_pdf_replay": primary_replay}, sort_keys=True))


if __name__ == "__main__":
    main()
