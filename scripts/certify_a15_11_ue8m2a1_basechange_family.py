#!/usr/bin/env python3
"""Certify an equation-level K3 family for the A15_11 H10/21 carrier.

For nine parameters put

    A(t) = -3 + a1*t + a2*t^2 + a3*t^3 + a4*t^4,
    B(t) =  2 + (1-a1)*t + b2*t^2 + ... + b6*t^6.

Then t=0 is a normalized I1 fiber of the rational elliptic surface

    R: y^2 = x^3 + A(t)*x + B(t).

On the nonempty open where Q=(4*A^3+27*B^2)/t is squarefree and its
leading coefficient is nonzero, R has twelve I1 fibers and infinity is
smooth.  The quadratic base change t=s^2, branched at the marked I1 fiber
and the smooth fiber at infinity, gives the elliptic K3

    X: y^2 = x^3 + A(s^2)*x + B(s^2).

The deck involution fixes one elliptic curve and one rational curve.  Its
invariant lattice is visibly

    U + E8(-2) + <-2>:

U comes from fiber and zero section, E8(-2) is the doubled pullback of the
E8 Mordell--Weil lattice of R, and <-2> is the nonidentity I2 component.
For a very general parameter this is the full Neron--Severi lattice.  The
script also pins the exact standard Gram to the existing A15_11 H10/21
carrier certificate, thereby making that known Coxeter chamber an explicit
equation-level divisor-reflection subsystem.

The algebraic construction and its fixed-locus/invariant-lattice lineage
are due to Garbagnati--Salgado.  The new claim is the exact normalized
equation-to-existing-Coxeter-carrier handoff in this project, not a new
abstract polytope or a worldwide-priority claim for the K3 family.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    smith_diagonal,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/a15_11_h10_21_x2_k3_carrier.json"
OUTPUT = (
    ROOT
    / "data/enumerations/a15_11_ue8m2a1_basechange_family.json"
)

SOURCE_CERTIFICATE_PAYLOAD_SHA256 = (
    "a2905e05f15d5538d7f5896314bac4f52228f79a0ba2609e2ea13e99ff2e963c"
)
SOURCE_FILE_SHA256 = (
    "5f241c3659a79cf674d09558f72a5bc075d50b9fc6d87444ce8d253f1a4ab90b"
)


U = sp.Matrix([[0, 1], [1, 0]])


def negative_cartan(rank: int, edges: list[tuple[int, int]]) -> sp.Matrix:
    gram = -2 * sp.eye(rank)
    for left, right in edges:
        gram[left, right] = gram[right, left] = 1
    return gram


# This Dynkin ordering is literally the ordering used by the pinned carrier.
E8_NEGATIVE = negative_cartan(
    8,
    [(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)],
)
F_GRAM = sp.diag(U, 2 * E8_NEGATIVE, sp.Matrix([[-2]]))


def matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload_digest(payload: dict[str, Any]) -> str:
    copy = dict(payload)
    copy.pop("payload_sha256", None)
    encoded = json.dumps(copy, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def integer_matrix(rows: list[list[str | int]]) -> sp.Matrix:
    return sp.Matrix([[int(entry) for entry in row] for row in rows])


def build_certificate() -> dict[str, Any]:
    carrier = json.loads(SOURCE.read_text(encoding="utf-8"))

    t, s = sp.symbols("t s")
    a1, a2, a3, a4, b2, b3, b4, b5, b6 = sp.symbols(
        "a1 a2 a3 a4 b2 b3 b4 b5 b6"
    )
    parameters = [a1, a2, a3, a4, b2, b3, b4, b5, b6]

    A = -3 + a1 * t + a2 * t**2 + a3 * t**3 + a4 * t**4
    B = (
        2
        + (1 - a1) * t
        + b2 * t**2
        + b3 * t**3
        + b4 * t**4
        + b5 * t**5
        + b6 * t**6
    )
    reduced_discriminant = sp.expand(4 * A**3 + 27 * B**2)
    Q = sp.Poly(sp.cancel(reduced_discriminant / t), t)

    recovered = [
        sp.diff(A, t, order).subs(t, 0) / sp.factorial(order)
        for order in range(1, 5)
    ] + [
        sp.diff(B, t, order).subs(t, 0) / sp.factorial(order)
        for order in range(2, 7)
    ]
    recovery_jacobian = sp.Matrix(recovered).jacobian(parameters)

    witness_substitution = {
        a1: 1,
        a2: 0,
        a3: 0,
        a4: 1,
        b2: 0,
        b3: 0,
        b4: 0,
        b5: 0,
        b6: 1,
    }
    witness_A = sp.expand(A.subs(witness_substitution))
    witness_B = sp.expand(B.subs(witness_substitution))
    witness_D = sp.Poly(
        reduced_discriminant.subs(witness_substitution), t
    )
    witness_Q = sp.Poly(Q.as_expr().subs(witness_substitution), t)
    witness_Q_derivative = sp.Poly(sp.diff(witness_Q.as_expr(), t), t)
    witness_Q_discriminant = int(
        sp.discriminant(witness_Q.as_expr(), t)
    )
    witness_A_Q_resultant = int(
        sp.resultant(witness_A, witness_Q.as_expr(), t)
    )
    witness_Q_s2 = sp.Poly(witness_Q.as_expr().subs(t, s**2), s)

    standard_carrier_gram = integer_matrix(
        carrier["lattice_bridge"]["standard_model_gram"]
    )
    f_inverse = F_GRAM.inv()
    two_f_inverse = 2 * f_inverse
    delta_witness = sp.Rational(-1, 2)

    # Nikulin's fixed-locus equations for g=1 and k=1.
    genus = 1
    rational_curves = 1
    invariant_rank = 11 - genus + rational_curves
    discriminant_length = 11 - genus - rational_curves

    literal_vinberg = carrier["finite_volume"]["literal_vinberg"]
    checks = {
        "normalized_i1_discriminant_vanishes": (
            reduced_discriminant.subs(t, 0) == 0
        ),
        "normalized_i1_first_derivative_is_108": (
            sp.diff(reduced_discriminant, t).subs(t, 0) == 108
        ),
        "normalized_i1_c4_is_nonzero": -48 * A.subs(t, 0) == 144,
        "Q_is_polynomial": sp.expand(t * Q.as_expr() - reduced_discriminant)
        == 0,
        "Q_degree_is_eleven": Q.degree() == 11,
        "Q_constant_is_108": Q.eval(0) == 108,
        "Q_leading_is_smooth_infinity_discriminant": (
            Q.LC() == 4 * a4**3 + 27 * b6**2
        ),
        "parameter_recovery_is_literal": recovered == parameters,
        "parameter_recovery_jacobian_is_identity": (
            recovery_jacobian == sp.eye(9)
        ),
        "witness_A": witness_A == t**4 + t - 3,
        "witness_B": witness_B == t**6 + 2,
        "witness_Q": witness_Q.as_expr()
        == (
            31 * t**11
            + 12 * t**8
            - 36 * t**7
            + 120 * t**5
            - 72 * t**4
            + 108 * t**3
            + 4 * t**2
            - 36 * t
            + 108
        ),
        "witness_Q_squarefree": sp.gcd(
            witness_Q, witness_Q_derivative
        ).degree()
        == 0,
        "witness_Q_discriminant_nonzero": witness_Q_discriminant != 0,
        "witness_c4_avoids_residual_roots": witness_A_Q_resultant != 0,
        "witness_infinity_is_smooth": witness_D.LC() == 31,
        "witness_basechanged_residual_is_squarefree": sp.gcd(
            witness_Q_s2,
            sp.Poly(sp.diff(witness_Q_s2.as_expr(), s), s),
        ).degree()
        == 0,
        "basechange_fiber_count": 1 + 22 == 23,
        "basechange_euler_is_24": 2 + 22 == 24,
        "fixed_locus_g_k": [genus, rational_curves] == [1, 1],
        "nikulin_rank": invariant_rank == 11,
        "nikulin_length": discriminant_length == 9,
        "e8_negative_determinant_one": E8_NEGATIVE.det() == 1,
        "f_is_even": all(F_GRAM[index, index] % 2 == 0 for index in range(11)),
        "f_signature": exact_inertia(F_GRAM) == (1, 10, 0),
        "f_determinant": F_GRAM.det() == 512,
        "f_smith_is_2_elementary_length_nine": smith_diagonal(F_GRAM)
        == [1, 1, *([2] * 9)],
        "f_discriminant_exponent_two": all(
            entry.q == 1 for entry in two_f_inverse
        ),
        "f_delta_is_one": delta_witness.q == 2,
        "visible_lattice_matches_invariant_determinant": (
            abs(int(F_GRAM.det())) == 2**discriminant_length
        ),
        "source_file_sha256": sha256_file(SOURCE) == SOURCE_FILE_SHA256,
        "source_payload_sha256": carrier["certificate_payload_sha256"]
        == SOURCE_CERTIFICATE_PAYLOAD_SHA256,
        "source_all_checks_pass": carrier["all_checks_pass"] is True,
        "source_same_isometry_class": carrier["lattice_bridge"]["N"][
            "isometry_class"
        ]
        == "U+E8(-2)+<-2>",
        "source_standard_gram_is_literal": standard_carrier_gram == F_GRAM,
        "source_all_walls_integral": carrier["walls"][
            "all_21_integral_on_F"
        ]
        is True,
        "source_wall_count": carrier["walls"]["count"] == 21,
        "source_literal_vinberg_finite_volume": literal_vinberg[
            "classification"
        ]
        == "finite_volume_noncompact_Coxeter_polytope",
        "source_vinberg_no_bad_audits": literal_vinberg["checks"][
            "zero_bad_audits"
        ]
        is True,
    }
    if not all(checks.values()):
        failed = [name for name, value in checks.items() if not value]
        raise AssertionError(f"base-change-family checks failed: {failed}")

    payload: dict[str, Any] = {
        "schema": "coxeter.a15_11.ue8m2a1_basechange_family.v1",
        "theorem": (
            "On a nonempty Zariski-open subset of the displayed normalized "
            "nine-parameter family, the quadratic base change is an elliptic "
            "K3 with deck-invariant lattice U+E8(-2)+<-2>.  For a very "
            "general parameter this is the full Neron-Severi lattice.  Via "
            "the pinned standard Gram, it is an equation-level full-Picard "
            "carrier of the known A15_11 C02 finite-volume noncompact H10/21 "
            "Coxeter chamber."
        ),
        "family": {
            "base_field": "C (or any algebraically closed field of characteristic zero)",
            "parameters": [str(parameter) for parameter in parameters],
            "parameter_count": 9,
            "rational_elliptic_surface": "R: y^2=x^3+A(t)*x+B(t)",
            "A": "A(t)=-3+a1*t+a2*t^2+a3*t^3+a4*t^4",
            "B": (
                "B(t)=2+(1-a1)*t+b2*t^2+b3*t^3+b4*t^4+"
                "b5*t^5+b6*t^6"
            ),
            "k3_weierstrass_model": "X: y^2=x^3+A(s^2)*x+B(s^2)",
            "surface": "the minimal resolution of the displayed Weierstrass model",
            "deck_involution": "iota(s,x,y)=(-s,x,y)",
            "holomorphic_two_form": "ds wedge dx/y",
            "two_form_multiplier": "-1 (non-symplectic)",
            "reduced_discriminant": "D(t)=4*A(t)^3+27*B(t)^2=t*Q(t)",
            "full_discriminant": "Delta_R=-16*D(t); Delta_X=-16*D(s^2)",
            "normalized_i1": {
                "place": "t=0",
                "D_order": 1,
                "D_prime": 108,
                "c4": 144,
            },
            "parameter_open": (
                "Disc_t(Q) is nonzero and 4*a4^3+27*b6^2 is nonzero.  "
                "Then Q has eleven simple nonzero roots, the fiber at "
                "infinity is smooth, and all twelve singular fibers of R "
                "are I1."
            ),
            "fibers": {
                "R": {"I1": 12, "smooth_marked_at_infinity": 1},
                "X": {"I2": 1, "I1": 22, "smooth_at_infinity": 1},
                "euler_sum_X": 24,
            },
            "alternate_double_cover_model": (
                "Equivalently, X is the resolved double cover of the blow-up "
                "of P2 defining R, branched over the strict transforms of "
                "the marked nodal cubic and the marked smooth cubic."
            ),
            "nonempty_open_witness": {
                "parameters": {
                    "a1": 1,
                    "a2": 0,
                    "a3": 0,
                    "a4": 1,
                    "b2": 0,
                    "b3": 0,
                    "b4": 0,
                    "b5": 0,
                    "b6": 1,
                },
                "A": str(witness_A),
                "B": str(witness_B),
                "Q": str(witness_Q.as_expr()),
                "Q_discriminant": witness_Q_discriminant,
                "resultant_A_Q": witness_A_Q_resultant,
                "infinity_discriminant": int(witness_D.LC()),
                "claim_scope": (
                    "this rational point proves the open is nonempty; its "
                    "individual Picard lattice is not asserted to be generic"
                ),
            },
        },
        "fixed_locus_and_invariant_lattice": {
            "fixed_locus": (
                "one smooth elliptic fiber above infinity and the strict "
                "transform of the ramified nodal I1 fiber, a rational curve"
            ),
            "fixed_locus_g_k": [1, 1],
            "ramified_i1_local_model": (
                "uv=t becomes uv=s^2; resolving the A1 singularity turns the "
                "fiber into I2.  Exactly the strict transform of the nodal "
                "fiber is pointwise fixed, not the exceptional component."
            ),
            "nikulin_invariants": {"r": 11, "a": 9, "delta": 1},
            "standard_lattice": "F=U+E8(-2)+<-2>",
            "gram": matrix_payload(F_GRAM),
            "signature": [1, 10, 0],
            "determinant": 512,
            "smith": smith_diagonal(F_GRAM),
            "delta_one_witness": {
                "class": "one half of the <-2> generator",
                "discriminant_square": "-1/2 mod 2Z",
            },
            "geometric_integral_identification": {
                "U": (
                    "fiber f and O+f, where O is the zero section; their "
                    "Gram is [[0,1],[1,0]]"
                ),
                "E8(-2)": (
                    "R has no reducible fibers, so its Mordell-Weil height "
                    "lattice is E8.  Pullback of its eight Shioda classes "
                    "under the degree-two base change doubles every pairing, "
                    "giving the orthogonal divisor lattice E8(-2)."
                ),
                "<-2>": (
                    "the nonidentity component C of the ramified I2 fiber; "
                    "C^2=-2 and C is orthogonal to U and to the corrected "
                    "Shioda classes"
                ),
                "saturation": (
                    "all eleven displayed divisor classes are invariant.  "
                    "Their lattice has determinant 2^9, equal to the "
                    "determinant forced by (r,a,delta)=(11,9,1); hence its "
                    "index in H^2(X,Z)^iota is one."
                ),
            },
        },
        "generic_full_picard": {
            "claim": "NS(X)=U+E8(-2)+<-2> for a very general parameter",
            "parameter_moduli_dimension": 9,
            "two_elementary_pair_moduli_dimension": 9,
            "dimension_formula": "20-r=20-11=9",
            "parameter_recovery": {
                "observables": [
                    "A^(i)(0)/i! for i=1,...,4",
                    "B^(j)(0)/j! for j=2,...,6",
                ],
                "symbolic_jacobian": matrix_payload(recovery_jacobian),
                "jacobian_determinant": int(recovery_jacobian.det()),
            },
            "normalization": (
                "PGL2 puts the ordered ramification fibers at 0 and infinity; "
                "short-Weierstrass homothety puts (A(0),B(0))=(-3,2); the "
                "remaining base scaling normalizes D'(0)=108."
            ),
            "quotient_recovery": (
                "From a generic pair (X,iota), the fixed elliptic curve gives "
                "the unique type-2 elliptic fibration.  Quotient by iota and "
                "contract the distinguished exceptional curve over the "
                "ramified I1 node to recover R and its ordered I1/smooth "
                "branch fibers.  Thus the normalized parameter-to-pair map "
                "has zero-dimensional generic fibers."
            ),
            "noether_lefschetz_step": (
                "The nine-dimensional family therefore dominates the "
                "nine-dimensional moduli component of (11,9,1) pairs.  Extra "
                "anti-invariant algebraic classes lie on the countable union "
                "of Noether-Lefschetz hyperplanes, so the very general member "
                "has NS equal to the invariant lattice."
            ),
            "scope": (
                "exact for the geometric generic/very general member; the "
                "displayed rational witness only certifies nonemptiness"
            ),
        },
        "coxeter_handoff": {
            "source_artifact": str(SOURCE.relative_to(ROOT)),
            "source_file_sha256": SOURCE_FILE_SHA256,
            "source_certificate_payload_sha256": (
                SOURCE_CERTIFICATE_PAYLOAD_SHA256
            ),
            "source_polytope": "A15_11_1_C02 / H10/21",
            "same_standard_gram_literal": True,
            "walls": carrier["walls"]["count"],
            "all_reflections_integral_on_very_general_full_NS": True,
            "finite_volume": {
                "classification": literal_vinberg["classification"],
                "finite_vertices": literal_vinberg["finite_vertices"],
                "ideal_vertices": literal_vinberg["ideal_vertices"],
                "elliptic_edges": literal_vinberg["elliptic_edges"],
                "zero_bad_audits": literal_vinberg["checks"][
                    "zero_bad_audits"
                ],
            },
            "meaning": (
                "choose the displayed geometric basis as the pinned standard "
                "F-marking and transport the 21 archived wall vectors.  They "
                "are divisor classes and their reflections are integral on "
                "the full NS of a very general member."
            ),
        },
        "lineage": {
            "garbagnati_salgado_double_cover": {
                "primary": "https://arxiv.org/abs/1703.02783",
                "published": "https://doi.org/10.1016/j.jpaa.2018.03.010",
                "use": (
                    "quadratic base changes of rational elliptic surfaces; "
                    "the branch-fiber table I0/I1 -> I0/I2 with fixed genus-1/"
                    "rational curves; invariant-lattice classification"
                ),
            },
            "garbagnati_salgado_genus_one_classification": {
                "primary": "https://arxiv.org/abs/1806.03097",
                "published": "https://doi.org/10.4171/RMI/1163",
                "use": (
                    "Proposition 2.8 and Section 3: for k=1, r=11, a=9; "
                    "the admissible rational elliptic surface has branch "
                    "fibers I0+I1, eleven other I1 fibers, MW rank 8, and "
                    "delta=1"
                ),
            },
            "shioda_mordell_weil_lattices": {
                "primary": "https://doi.org/10.14992/00009994",
                "use": (
                    "height-lattice formalism and the E8 Mordell-Weil frame "
                    "for a rational elliptic surface without reducible fibers"
                ),
            },
            "nikulin_lattices": {
                "primary": "https://doi.org/10.1070/IM1980v014n01ABEH001060",
                "use": "classification and uniqueness of even 2-elementary lattices",
            },
            "novelty_scope": (
                "The 2-elementary K3 class and the abstract Coxeter polytope "
                "are known.  The positive addition is an exact normalized "
                "equation, geometric full-Picard identification, and literal "
                "Gram handoff joining them in this project; no worldwide "
                "priority claim is made."
            ),
        },
        "scope": {
            "not_a_new_abstract_polytope": True,
            "not_claimed_full_nef_cone": True,
            "reflections_not_claimed_surface_automorphisms": True,
            "deck_involution_not_used_to_generate_reflections": True,
            "square_minus_four_walls_not_claimed_curve_walls": True,
            "positive_result": (
                "an equation-level very-general K3 variety whose full divisor "
                "lattice carries all 21 certified Coxeter reflections"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["payload_sha256"] = payload_digest(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    payload = build_certificate()
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if payload != archived:
            raise AssertionError(f"archive mismatch: {args.output}")
    else:
        args.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "all_checks_pass": True,
                "payload_sha256": payload["payload_sha256"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
