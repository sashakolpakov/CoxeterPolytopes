#!/usr/bin/env python3
"""Exact equation-to-polytope certificate for Vinberg's X_3 chamber.

Two very different algebraic presentations of the same singular K3 surface
are certified:

1. the extremal Inose fibration

       y^2 = x^3 + t^5 (t-1)^2,

   with fibers II*, II*, IV and trivial Mordell--Weil group; and

2. the generalized Kummer/order-three quotient model, in a constant-rescaled
   Weierstrass form

       y^2 = x^3 + (t^3-1)^4,

   with three IV* fibers and Mordell--Weil torsion Z/3.

Both have Neron--Severi lattice

    U + 2 E8(-1) + A2(-1),

so a marking transports the already complete 36-wall Vinberg system to their
divisor lattices.  The script verifies the Weierstrass factorizations, the
index-three E6^3 glue, and the exact handoff to the archived complete-wall and
finite-volume certificate in data/enumerations/vinberg_x3_h19.json.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.certify_vinberg_x3_h19 import build_wall_gram
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
POLYTOPE_ARCHIVE = ROOT / "data" / "enumerations" / "vinberg_x3_h19.json"
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "vinberg_x3_h19_explicit_carriers.json"
)


def negative_cartan(rank: int, edges: list[tuple[int, int]]) -> sp.Matrix:
    gram = -2 * sp.eye(rank)
    for left, right in edges:
        gram[left, right] = gram[right, left] = 1
    return gram


U = sp.Matrix([[0, 1], [1, 0]])
A2 = negative_cartan(2, [(0, 1)])
E6 = negative_cartan(
    6,
    [(0, 1), (1, 2), (2, 3), (3, 4), (2, 5)],
)
E8 = negative_cartan(
    8,
    [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)],
)


def matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def rational_matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [str(sp.factor(matrix[row, column])) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def normalized_smith(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return sorted(
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index] != 0
    )


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload_digest(payload: dict[str, Any]) -> str:
    copy = dict(payload)
    copy.pop("payload_sha256", None)
    encoded = json.dumps(copy, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def order_mod_two(value: sp.Rational) -> sp.Rational:
    return sp.Rational(value - 2 * sp.floor(value / 2))


def first_discriminant_norm_mod_two(gram: sp.Matrix) -> sp.Rational:
    inverse = gram.inv()
    for column in range(inverse.cols):
        vector = inverse[:, column]
        if any(sp.denom(value) != 1 for value in vector):
            return order_mod_two(sp.factor((vector.T * gram * vector)[0]))
    raise AssertionError("unimodular matrix has no discriminant generator")


def e6_diagonal_glue() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return the E6^3 trivial lattice, glue basis, and glued Gram matrix.

    The first dual fundamental weight has square -4/3.  Its diagonal sum in
    the three E6 factors has square -4 and gives the isotropic order-three
    class met by a nonzero 3-torsion section.  The basis matrix has determinant
    1/3, so the resulting even integral lattice is exactly the index-three
    primitive closure of U+E6(-1)^3.
    """

    trivial = sp.diag(U, E6, E6, E6)
    fundamental_weight = E6.inv()[:, 0]
    glue = sp.Matrix(
        [0, 0]
        + list(fundamental_weight)
        + list(fundamental_weight)
        + list(fundamental_weight)
    )
    integral_glue = 3 * glue

    # The first two E6 coordinates of 3*glue are -4 and -5.  Together with
    # e_0+e_1 they form a unimodular 2-by-2 minor (determinant one).  This
    # gives a literal basis of Z^20 + Z*glue after dividing the first column
    # by three.
    standard = sp.eye(20)
    unimodular_columns = [integral_glue, standard[:, 2] + standard[:, 3]]
    unimodular_columns.extend(
        standard[:, index] for index in range(20) if index not in (2, 3)
    )
    basis = sp.Matrix.hstack(*unimodular_columns)
    assert basis.det() == 1
    basis[:, 0] = glue
    glued = sp.simplify(basis.T * trivial * basis)
    return trivial, basis, glued


def build_certificate() -> dict[str, Any]:
    archived = json.loads(POLYTOPE_ARCHIVE.read_text(encoding="utf-8"))
    wall_gram = build_wall_gram()

    t, u, x = sp.symbols("t u x")

    # The II*+II*+IV Inose model.
    inose_a6 = t**5 * (t - 1) ** 2
    inose_delta = sp.factor(-432 * inose_a6**2)
    inose_a6_infinity = sp.factor(u**12 * inose_a6.subs(t, 1 / u))
    inose_delta_infinity = sp.factor(
        u**24 * inose_delta.subs(t, 1 / u)
    )

    # The 3 IV* quotient model, after a harmless constant complex rescaling
    # of y^2=x^3-(t^3-1)^4.
    quotient_a6 = (t**3 - 1) ** 4
    quotient_delta = sp.factor(-432 * quotient_a6**2)
    quotient_a6_infinity = sp.factor(u**12 * quotient_a6.subs(t, 1 / u))
    quotient_delta_infinity = sp.factor(
        u**24 * quotient_delta.subs(t, 1 / u)
    )
    torsion_y = (t**3 - 1) ** 2

    target_gram = sp.diag(U, E8, E8, A2)
    trivial_e6, e6_glue_basis, glued_e6 = e6_diagonal_glue()

    target_smith = normalized_smith(target_gram)
    glued_smith = normalized_smith(glued_e6)
    target_q = first_discriminant_norm_mod_two(target_gram)
    glued_q = first_discriminant_norm_mod_two(glued_e6)

    source_checks = {
        "archive_all_checks_pass": archived["all_checks_pass"] is True,
        "wall_gram_reconstructs_archive": wall_gram
        == sp.Matrix(archived["wall_gram"]),
        "complete_wall_count_36": archived["facet_count"] == 36,
        "hyperbolic_dimension_19": archived["hyperbolic_dimension"] == 19,
        "wall_gram_rank_20": wall_gram.rank() == 20,
        "wall_gram_inertia_1_19_16": exact_inertia(wall_gram)
        == (1, 19, 16),
        "wall_gram_smith_1^19_3": normalized_smith(wall_gram)
        == [1] * 19 + [3],
        "vinberg_every_edge_has_two_endpoints": archived[
            "vinberg_endpoint_audit"
        ]["bad_edge_count"]
        == 0
        and archived["vinberg_endpoint_audit"]["extension_histogram"]
        == {"2": 87046},
        "vinberg_no_endpoint_generated_nonedges": archived[
            "vinberg_endpoint_audit"
        ]["finite_generated_nonedge_count"]
        == 0
        and archived["vinberg_endpoint_audit"][
            "ideal_generated_nonedge_count"
        ]
        == 0,
        "finite_volume_noncompact": archived["volume"]["finite_volume"]
        is True
        and archived["volume"]["compact"] is False,
    }

    checks = {
        "inose_discriminant_factorization": inose_delta
        == -432 * t**10 * (t - 1) ** 4,
        "inose_infinity_a6": sp.expand(
            inose_a6_infinity - u**5 * (1 - u) ** 2
        )
        == 0,
        "inose_infinity_delta": sp.expand(
            inose_delta_infinity + 432 * u**10 * (1 - u) ** 4
        )
        == 0,
        "inose_euler_sum_24": 10 + 10 + 4 == 24,
        "quotient_discriminant_factorization": sp.expand(
            quotient_delta + 432 * (t**3 - 1) ** 8
        )
        == 0,
        "quotient_three_distinct_fiber_places": sp.resultant(
            t**3 - 1, 3 * t**2, t
        )
        != 0,
        "quotient_infinity_smooth": sp.expand(
            quotient_a6_infinity - (1 - u**3) ** 4
        )
        == 0
        and sp.expand(quotient_delta_infinity + 432 * (1 - u**3) ** 8)
        == 0,
        "quotient_euler_sum_24": 3 * 8 == 24,
        "displayed_section_lies_on_quotient_model": sp.expand(
            torsion_y**2 - quotient_a6
        )
        == 0,
        "displayed_section_is_nonzero_3_torsion_by_triple_tangent": sp.expand(
            torsion_y**2 - (x**3 + quotient_a6)
        )
        == -x**3,
        "target_rank_signature_determinant": target_gram.rank() == 20
        and exact_inertia(target_gram) == (1, 19, 0)
        and target_gram.det() == -3,
        "target_smith_1^19_3": target_smith == [1] * 19 + [3],
        "e6_trivial_rank_signature_determinant": trivial_e6.rank() == 20
        and exact_inertia(trivial_e6) == (1, 19, 0)
        and trivial_e6.det() == -27,
        "diagonal_e6_glue_has_order_3": e6_glue_basis.det()
        == sp.Rational(1, 3),
        "diagonal_e6_glue_is_even_integral": all(
            sp.denom(value) == 1 for value in glued_e6
        )
        and all(glued_e6[index, index] % 2 == 0 for index in range(20)),
        "glued_e6_rank_signature_determinant": glued_e6.rank() == 20
        and exact_inertia(glued_e6) == (1, 19, 0)
        and glued_e6.det() == -3,
        "glued_e6_smith_1^19_3": glued_smith == [1] * 19 + [3],
        "target_and_glued_discriminant_forms_match": target_q
        == glued_q
        == sp.Rational(4, 3),
        **source_checks,
    }

    six_fibrations = [
        {
            "id": "E1",
            "equation": "y^2=x^3+t^5*(t-1)^2",
            "reducible_kodaira_fibers": ["II*", "II*", "IV"],
            "frame_root_lattice": "E8(-1)^2+A2(-1)",
            "mordell_weil": "0",
            "directly_certified_here": True,
        },
        {
            "id": "E2",
            "equation": "y^2=x^3-t^3*x^2+t^3*x",
            "reducible_kodaira_fibers": ["I6*", "III*"],
            "frame_root_lattice": "D10(-1)+E7(-1)",
            "mordell_weil": "Z + Z/2",
            "directly_certified_here": False,
        },
        {
            "id": "E3",
            "equation": "y^2=x^3+(t^3-2)*x^2+x",
            "reducible_kodaira_fibers": ["I12*", "I3"],
            "frame_root_lattice": "D16(-1)+A2(-1)",
            "mordell_weil": "Z/2",
            "directly_certified_here": False,
        },
        {
            "id": "E4",
            "equation": "y^2=x^3+t^4*x^2-8*t^2*x+16",
            "reducible_kodaira_fibers": ["I18"],
            "frame_root_lattice": "A17(-1)",
            "mordell_weil": "Z + Z/3",
            "directly_certified_here": False,
        },
        {
            "id": "E5",
            "equation": "y^2=x^3+(t^4-8*t)*x^2+16*t^2*x",
            "reducible_kodaira_fibers": ["I12", "I3*"],
            "frame_root_lattice": "A11(-1)+D7(-1)",
            "mordell_weil": "Z/4",
            "directly_certified_here": False,
        },
        {
            "id": "E6",
            "equation": "y^2-((t+1)^2/t^4)*y=x^3",
            "reducible_kodaira_fibers": ["IV*", "IV*", "IV*"],
            "frame_root_lattice": "E6(-1)^3",
            "mordell_weil": "Z/3",
            "directly_certified_here": False,
            "directly_certified_constant_rescaled_quotient_model": (
                "y^2=x^3+(t^3-1)^4"
            ),
        },
    ]

    payload: dict[str, Any] = {
        "artifact_status": (
            "exact_explicit_algebro_geometric_carriers_for_known_H19_chamber"
        ),
        "novelty_status": (
            "new repository-level equation/abelian-quotient/divisor-lattice "
            "realization of the known Vinberg X3 chamber; neither the abstract "
            "polytope nor the published equations are claimed new worldwide"
        ),
        "positive_theorem": {
            "surface": "Vinberg's singular K3 surface X3",
            "neron_severi": "U+2E8(-1)+A2(-1)",
            "transcendental": "A2",
            "hyperbolic_dimension": 19,
            "coxeter_facets": 36,
            "finite_volume": True,
            "compact": False,
            "explicit_carrier_count_certified_here": 2,
            "published_elliptic_fibrations_up_to_automorphism": 6,
        },
        "inose_carrier": {
            "affine_equation": "y^2=x^3+t^5*(t-1)^2",
            "homogeneous_equation": (
                "Y^2=X^3+T^5*(T-S)^2*S^5, deg(X)=4, deg(Y)=6"
            ),
            "c4": "0",
            "discriminant": str(inose_delta),
            "fibers": [
                {
                    "place": "t=0",
                    "valuations_a6_delta": [5, 10],
                    "kodaira": "II*",
                    "root_lattice": "E8(-1)",
                    "euler": 10,
                },
                {
                    "place": "t=1",
                    "valuations_a6_delta": [2, 4],
                    "kodaira": "IV",
                    "root_lattice": "A2(-1)",
                    "euler": 4,
                },
                {
                    "place": "t=infinity",
                    "valuations_a6_delta": [5, 10],
                    "kodaira": "II*",
                    "root_lattice": "E8(-1)",
                    "euler": 10,
                },
            ],
            "trivial_lattice": "U+2E8(-1)+A2(-1)",
            "rank": 20,
            "determinant": -3,
            "mordell_weil": "0",
            "ns_exactness": (
                "rank 20 reaches the complex K3 Picard bound; an index m "
                "same-rank integral overclosure would require m^2 dividing 3, "
                "so m=1"
            ),
        },
        "abelian_quotient_carrier": {
            "abelian_surface": "E_omega x E_omega, E_omega: v^2=u^3-1",
            "group": (
                "<alpha x alpha^2>, alpha(v,u)=(v,omega*u), omega^3=1"
            ),
            "quotient_singularities": "9 rational double points of type A2",
            "surface": "minimal resolution of the quotient",
            "weierstrass_equation_after_constant_complex_rescaling": (
                "y^2=x^3+(t^3-1)^4"
            ),
            "source_equation_before_constant_rescaling": (
                "y^2=x^3-(t^3-s^3)^4"
            ),
            "fibers": {
                "count": 3,
                "kodaira": "IV*",
                "root_lattice_each": "E6(-1)",
                "euler_each": 8,
            },
            "three_torsion_sections": [
                "(x,y)=(0,(t^3-1)^2)",
                "(x,y)=(0,-(t^3-1)^2)",
            ],
            "trivial_lattice": "U+E6(-1)^3",
            "trivial_lattice_determinant": -27,
            "primitive_closure_index": 3,
            "primitive_closure": "U+(E6(-1)^3)'",
            "primitive_closure_determinant": -3,
            "primitive_closure_gram": matrix_payload(glued_e6),
            "diagonal_glue_basis_in_trivial_coordinates": rational_matrix_payload(
                e6_glue_basis
            ),
            "discriminant_generator_q_mod_2": str(glued_q),
            "identified_ns": "U+2E8(-1)+A2(-1)",
            "transcendental": "A2",
        },
        "six_fibration_atlas": {
            "status": (
                "published by Brandhorst--Zach, Theorem 2.5; E1 and a "
                "constant-rescaled E6 quotient model are independently "
                "certified in this artifact"
            ),
            "models": six_fibrations,
        },
        "lattice_identification": {
            "target_gram": matrix_payload(target_gram),
            "rank_signature_determinant": [20, [1, 19, 0], -3],
            "smith": target_smith,
            "discriminant_form_nonzero_q_mod_2": str(target_q),
            "identification_reason": (
                "the Inose fiber components give the target lattice literally; "
                "the diagonal E6^3 primitive closure has the same signature and "
                "discriminant form and is the unique even indefinite lattice in "
                "this genus (equivalently both models are X3)"
            ),
        },
        "hyperbolic_handoff": {
            "polytope_archive": str(POLYTOPE_ARCHIVE.relative_to(ROOT)),
            "polytope_archive_sha256": sha256_file(POLYTOPE_ARCHIVE),
            "wall_gram_reconstructed_not_copied": True,
            "dimension": archived["hyperbolic_dimension"],
            "facets": archived["facet_count"],
            "root_norm_counts": archived["lattice"]["root_norm_counts"],
            "finite_vertices": archived["vertices_and_cusps"][
                "finite_vertex_count"
            ],
            "cusps": archived["vertices_and_cusps"]["ideal_vertex_count"],
            "elliptic_rank18_edges": archived["vinberg_endpoint_audit"][
                "elliptic_rank18_edge_count"
            ],
            "endpoint_type_histogram": archived["vinberg_endpoint_audit"][
                "endpoint_type_histogram"
            ],
            "bad_edges": archived["vinberg_endpoint_audit"]["bad_edge_count"],
            "finite_volume": archived["volume"]["finite_volume"],
            "compact": archived["volume"]["compact"],
            "completeness": archived["completeness"],
        },
        "algebro_geometric_interpretation": {
            "divisor_carrier": (
                "a marking NS(X3)=U+2E8(-1)+A2(-1) transports all 36 "
                "primitive wall normals to divisor classes on either explicit "
                "surface model"
            ),
            "long_walls": (
                "24 square-minus-two divisibility-one walls are the rational-"
                "curve walls of the Vinberg/OSCAR chamber"
            ),
            "short_walls": (
                "12 square-minus-six divisibility-three integral divisor "
                "classes supply the additional Vinberg/Borcherds cuts; they "
                "are not asserted to be irreducible curves"
            ),
            "nef_cone_relation": (
                "Brandhorst--Zach compute a 36-wall chamber contained in the "
                "nef cone of this explicit X3; Vinberg--Looijenga supply the "
                "complete reflection and finite-volume theorem"
            ),
        },
        "sources": [
            {
                "authors": "Simon Brandhorst and Matthias Zach",
                "title": "Elliptic Fibrations on Vinberg's Most Algebraic K3 Surface",
                "arxiv": "2311.11766",
                "locator": "Sections 2.1-2.6, especially Theorem 2.5",
            },
            {
                "author": "È. B. Vinberg",
                "title": "The Two Most Algebraic K3 Surfaces",
                "doi": "10.1007/BF01456933",
                "year": 1983,
            },
            {
                "author": "Eduard Looijenga",
                "title": "The period map for cubic fourfolds",
                "doi": "10.1007/s00222-009-0178-6",
                "locator": "Section 1, Properties 1.4 and Figure 1",
            },
            {
                "author": "Alice Garbagnati Prieto Montañez",
                "title": "Automorphisms on Algebraic Varieties: K3 Surfaces, Fano Varieties and Cubic Fourfolds",
                "locator": (
                    "Chapter 2, Propositions 2.1.2, 2.2.3, 2.2.4 and equations "
                    "(2.1)-(2.3)"
                ),
                "url": "https://amsdottorato.unibo.it/10149/",
            },
        ],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["payload_sha256"] = payload_digest(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if args.check_archive:
        if not OUTPUT.exists():
            raise SystemExit(f"missing archive: {OUTPUT}")
        frozen = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if frozen != payload:
            raise SystemExit("frozen explicit-X3 carrier certificate differs")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "all_checks_pass": payload["all_checks_pass"],
                "carriers": payload["positive_theorem"][
                    "explicit_carrier_count_certified_here"
                ],
                "dimension": payload["positive_theorem"][
                    "hyperbolic_dimension"
                ],
                "facets": payload["positive_theorem"]["coxeter_facets"],
                "finite_volume": payload["positive_theorem"]["finite_volume"],
                "payload_sha256": payload["payload_sha256"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
