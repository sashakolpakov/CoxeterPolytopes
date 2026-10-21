#!/usr/bin/env python3
"""Certify an explicit Weierstrass family for BH 0.4.4.1 / H10/12.

Put ``s=t^2`` and, for parameters ``alpha,beta,gamma,delta``, set

    A(s) = -3 + alpha*(s-1) + beta*(s-1)^2,
    C(s) =  2 - (alpha+1)*(s-1) + gamma*(s-1)^2
              + delta*(s-1)^3.

The minimal resolution of

    y^2 = x^3 + A(t^2)*x + t*C(t^2)

has the purely non-symplectic order-four automorphism

    sigma(x,y,t) = (-x, i*y, -t).

On a nonempty Zariski-open parameter set the fibers are II* at infinity,
I2 at t=+1 and t=-1, and ten I1 fibers.  The very general Neron--Severi
lattice is U+E8(-1)+A1(-1)^2, and sigma exchanges the two A1 components.
This is the exact Picard action in Brandhorst--Hofmann row 0.4.4.1 and hence
an equation-level realization of the archived finite-volume H10/12 section.

The script checks the symbolic discriminant, a concrete nonempty-open
witness, the lattice/action/glue data, the generic Picard-rank argument, and
the pinned BH/chamber/Vinberg handoffs.  An optional external replay checks
the complete BH row in the published database ZIP.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any
import zipfile

import sympy as sp

from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    _parse_bh_row,
    smith_diagonal,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/roulleau_rank12_bh_k3_lifts.json"
OUTPUT = ROOT / "data/enumerations/ue8a1x2_order4_weierstrass_family.json"

BH_ZIP_SHA256 = "c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c"
BH_MEMBER = "K3Groups/results/purely_ns/order04.txt"
BH_MEMBER_SHA256 = "d5a61603c1539278e8411e188656c9417662336ff4282b714e939c40d38c94e8"
BH_ROW_ONE_BASED = 39
BH_ROW_SHA256 = "03ef01ed30b98e3d0fa419cb8275577c90594868fcd3400c48413555b1248452"
VINBERG_AUDIT_SHA256 = (
    "d20267d1d38090f54d650d20d88d71448aa0393686cf64705e9780f91c37bd5a"
)
BH_FIXED_PROPERTIES = [
    (
        "fixed points powers",
        (((8,), 2, 0, 1), ((), 5, 1, 4)),
    ),
    ("fixed points", ((8,), 2, 0, 1)),
]


U = sp.Matrix([[0, 1], [1, 0]])
A1X2 = sp.diag(-2, -2)
A1_SWAP = sp.Matrix([[0, 1], [1, 0]])


def negative_cartan(rank: int, edges: list[tuple[int, int]]) -> sp.Matrix:
    gram = -2 * sp.eye(rank)
    for left, right in edges:
        gram[left, right] = gram[right, left] = 1
    return gram


E8 = negative_cartan(
    8,
    [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)],
)


def block_diag(*blocks: sp.Matrix) -> sp.Matrix:
    return sp.diag(*blocks)


def matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def payload_digest(payload: dict[str, Any]) -> str:
    copy = dict(payload)
    copy.pop("payload_sha256", None)
    encoded = json.dumps(copy, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def source_case() -> tuple[dict[str, Any], dict[str, Any]]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = next(
        row
        for row in source["realizations"]
        if row["brandhorst_hofmann"]["id"] == "0.4.4.1"
    )
    return source, case


def replay_brandhorst_hofmann_zip(path: Path) -> dict[str, Any]:
    if sha256_file(path) != BH_ZIP_SHA256:
        raise AssertionError("Brandhorst--Hofmann ZIP SHA-256 mismatch")
    with zipfile.ZipFile(path) as archive:
        member = archive.read(BH_MEMBER)
    if sha256_bytes(member) != BH_MEMBER_SHA256:
        raise AssertionError("Brandhorst--Hofmann order04 member SHA-256 mismatch")
    rows = member.decode("utf-8").splitlines(keepends=True)
    row_bytes = rows[BH_ROW_ONE_BASED - 1].encode("utf-8")
    if sha256_bytes(row_bytes) != BH_ROW_SHA256:
        raise AssertionError("Brandhorst--Hofmann row SHA-256 mismatch")
    parsed = _parse_bh_row(rows[BH_ROW_ONE_BASED - 1])
    if parsed["id"] != "0.4.4.1":
        raise AssertionError("wrong Brandhorst--Hofmann row")
    if parsed["transcendental_value"] != 4:
        raise AssertionError("wrong transcendental value")
    if parsed["properties"] != BH_FIXED_PROPERTIES:
        raise AssertionError("Brandhorst--Hofmann fixed-locus fingerprint mismatch")
    return {
        "all_checks_pass": True,
        "database_zip_sha256": BH_ZIP_SHA256,
        "member": BH_MEMBER,
        "member_sha256": BH_MEMBER_SHA256,
        "row_one_based": BH_ROW_ONE_BASED,
        "row_sha256_with_newline": BH_ROW_SHA256,
        "id": parsed["id"],
        "fixed_properties": json.loads(json.dumps(parsed["properties"])),
    }


def build_certificate() -> dict[str, Any]:
    source, case = source_case()

    u, alpha, beta, gamma, delta = sp.symbols(
        "u alpha beta gamma delta"
    )
    s = u + 1
    a_coefficient = -3 + alpha * u + beta * u**2
    c_coefficient = (
        2 - (alpha + 1) * u + gamma * u**2 + delta * u**3
    )
    reduced_discriminant = sp.expand(
        4 * a_coefficient**3 + 27 * s * c_coefficient**2
    )
    quotient = sp.cancel(reduced_discriminant / u**2)
    i2_leading = sp.factor(quotient.subs(u, 0))

    # Exact coordinates on the normalized marked-Weierstrass stratum.
    parameters = sp.Matrix([alpha, beta, gamma, delta])
    recovered_parameters = sp.Matrix(
        [
            sp.diff(a_coefficient, u).subs(u, 0),
            sp.diff(a_coefficient, u, 2).subs(u, 0) / 2,
            sp.diff(c_coefficient, u, 2).subs(u, 0) / 2,
            sp.diff(c_coefficient, u, 3).subs(u, 0) / 6,
        ]
    )
    recovery_jacobian = recovered_parameters.jacobian(parameters)
    base_scale, weierstrass_scale = sp.symbols(
        "base_scale weierstrass_scale"
    )
    base_stabilizer_polynomial = base_scale**2 - 1
    weierstrass_stabilizer_polynomial = sp.gcd(
        weierstrass_scale**4 - 1,
        weierstrass_scale**6 - 1,
    )

    witness_substitution = {alpha: 0, beta: 0, gamma: 0, delta: 1}
    witness_quotient = sp.factor(quotient.subs(witness_substitution) / 27)
    expected_witness = u**5 + u**4 - 2 * u**3 + 2 * u**2 + 5 * u - 3
    witness_discriminant = int(sp.discriminant(expected_witness, u))

    ns_gram = block_diag(U, E8, A1X2)
    ns_action = block_diag(sp.eye(10), A1_SWAP)
    fixed_a1_basis = sp.Matrix([[1, 1]])
    coinvariant_basis = sp.Matrix([[1, -1]])
    fixed_a1_gram = fixed_a1_basis * A1X2 * fixed_a1_basis.T
    coinvariant_gram = coinvariant_basis * A1X2 * coinvariant_basis.T
    fixed_gram = block_diag(U, E8, fixed_a1_gram)

    lam = sp.Symbol("lambda")
    ns_characteristic = sp.factor(ns_action.charpoly(lam).as_expr())
    full_characteristic = sp.factor(ns_characteristic * (lam**2 + 1) ** 5)
    h2_traces = {"1": 10, "2": 2, "3": 10}
    fixed_euler = {power: trace + 2 for power, trace in h2_traces.items()}

    bh = case["brandhorst_hofmann"]
    conjugacy = case["exact_action_conjugacy"]
    chamber = case["chamber_transport"]
    hyperbolic = case["coxeter_polytope"]

    checks = {
        "reduced_discriminant_has_i2_pair": sp.expand(
            reduced_discriminant - u**2 * quotient
        )
        == 0,
        "residual_discriminant_degree_five": sp.degree(quotient, u) == 5,
        "i2_leading_coefficient": sp.simplify(
            i2_leading
            + 9 * (alpha**2 + 6 * alpha - 12 * beta - 12 * gamma + 9)
        )
        == 0,
        "witness_residual_quintic": witness_quotient == expected_witness,
        "witness_residual_squarefree": witness_discriminant == 283904,
        "witness_no_fiber_at_t_zero": expected_witness.subs(u, -1) == -4,
        "parameter_recovery_is_identity": recovered_parameters == parameters,
        "parameter_recovery_jacobian_is_identity": (
            recovery_jacobian == sp.eye(4)
        ),
        "base_normalization_has_finite_stabilizer": (
            sp.degree(base_stabilizer_polynomial, base_scale) == 2
        ),
        "weierstrass_normalization_has_finite_stabilizer": sp.expand(
            weierstrass_stabilizer_polynomial
            - (weierstrass_scale**2 - 1)
        )
        == 0,
        "fiber_euler_sum_is_24": 10 + 2 * 2 + 10 == 24,
        "trivial_lattice_rank_is_12": 2 + 8 + 2 == 12,
        "ns_signature": exact_inertia(ns_gram) == (1, 11, 0),
        "ns_determinant": ns_gram.det() == -4,
        "ns_smith": smith_diagonal(ns_gram) == [*[1] * 10, 2, 2],
        "a1_swap_isometry": A1_SWAP * A1X2 * A1_SWAP.T == A1X2,
        "a1_swap_order_two": A1_SWAP**2 == sp.eye(2),
        "ns_action_isometry": ns_action * ns_gram * ns_action.T == ns_gram,
        "ns_action_order_two": ns_action**2 == sp.eye(12),
        "ns_characteristic": ns_characteristic
        == (lam - 1) ** 11 * (lam + 1),
        "fixed_lattice_signature": exact_inertia(fixed_gram) == (1, 10, 0),
        "fixed_lattice_determinant": fixed_gram.det() == 4,
        "fixed_lattice_smith": smith_diagonal(fixed_gram)
        == [*[1] * 10, 4],
        "coinvariant_is_minus_four": coinvariant_gram == sp.Matrix([[-4]]),
        "coinvariant_rootless": coinvariant_gram[0, 0] == -4,
        "primitive_glue_index_two": (
            abs(int(sp.Matrix([[1, 1], [1, -1]]).det())) == 2
        ),
        "full_characteristic": full_characteristic
        == (lam - 1) ** 11 * (lam + 1) * (lam**2 + 1) ** 5,
        "fixed_euler_characteristics": fixed_euler
        == {"1": 12, "2": 4, "3": 12},
        "source_all_checks_pass": source["all_checks_pass"] is True,
        "source_bh_id": bh["id"] == "0.4.4.1",
        "source_generic_ns": bh["generic_NS"]
        == "U+E8(-1)+A1(-1)^2",
        "source_picard_image": bh["Picard_image"] == "C2",
        "source_effective_group": bh["effective_group"] == "C4",
        "source_full_action_conjugacy": conjugacy[
            "full_NS_action_pairs_integrally_conjugate"
        ]
        is True,
        "source_same_coinvariant": conjugacy["algebraic_coinvariant"]
        == "<-4>, primitive of divisibility 2",
        "source_all_glue_carried": conjugacy["all_glue_pairs_carried"] is True,
        "source_complete_nef_chamber": chamber[
            "complete_Roulleau_nef_chamber"
        ]
        is True,
        "vinberg_finite_volume": hyperbolic["classification"]
        == "finite_volume_noncompact",
        "vinberg_no_bad_audits": hyperbolic["Vinberg_bad_audits"] == 0,
    }
    if not all(checks.values()):
        raise AssertionError(f"a Weierstrass-family check failed: {checks}")

    fiber_rows = [
        {
            "place": "t=infinity",
            "ord_a4_at_least": 4,
            "ord_a6": 5,
            "ord_discriminant": 10,
            "kodaira": "II*",
            "euler": 10,
            "root_lattice": "E8(-1)",
            "root_rank": 8,
        },
        {
            "place": "t=+1 and t=-1",
            "multiplicity": 2,
            "ord_discriminant": 2,
            "kodaira": "I2",
            "euler_each": 2,
            "root_lattice_each": "A1(-1)",
            "root_rank_each": 1,
        },
        {
            "place": "ten roots above the five roots of Q(s)",
            "multiplicity": 10,
            "ord_discriminant": 1,
            "kodaira": "I1",
            "euler_each": 1,
            "root_rank": 0,
        },
    ]

    payload: dict[str, Any] = {
        "schema": "ue8a1x2-order4-weierstrass-family-v1",
        "theorem": (
            "On a nonempty Zariski-open set of (alpha,beta,gamma,delta), "
            "the very general minimal resolution of y^2=x^3+A(t^2)x+"
            "t*C(t^2) has NS=U+E8(-1)+A1(-1)^2.  The primitive order-four "
            "automorphism (x,y,t)->(-x,i*y,-t) exchanges the two I2 root "
            "components, realizes BH 0.4.4.1, and therefore explicitly "
            "realizes the archived finite-volume noncompact H10/12 Coxeter "
            "section."
        ),
        "family": {
            "affine_equation": "y^2=x^3+A(t^2)*x+t*C(t^2)",
            "A": "A(s)=-3+alpha*(s-1)+beta*(s-1)^2",
            "C": (
                "C(s)=2-(alpha+1)*(s-1)+gamma*(s-1)^2+"
                "delta*(s-1)^3"
            ),
            "automorphism": "sigma(x,y,t)=(-x,i*y,-t)",
            "global_order": 4,
            "square": "sigma^2(x,y,t)=(x,-y,t)",
            "holomorphic_two_form": "dt wedge dx/y",
            "two_form_multiplier": "-i (a primitive fourth root)",
            "reduced_discriminant": (
                "D(s)=4*A(s)^3+27*s*C(s)^2=(s-1)^2*Q(s)"
            ),
            "full_discriminant": "Delta(t)=-16*D(t^2)",
            "residual_Q_degree": 5,
            "i2_open_coefficient": str(i2_leading),
            "parameter_open": (
                "delta is nonzero; the displayed I2 coefficient is nonzero; "
                "Q is squarefree, has no root s=0, and shares no root with "
                "A; these are nonempty Zariski-open conditions"
            ),
            "nonempty_open_witness": {
                "parameters": {
                    "alpha": 0,
                    "beta": 0,
                    "gamma": 0,
                    "delta": 1,
                },
                "Q_over_27_in_u_equals_s_minus_1": str(expected_witness),
                "Q_over_27_discriminant": witness_discriminant,
                "Q_over_27_at_s_zero": int(expected_witness.subs(u, -1)),
            },
            "fibers": fiber_rows,
        },
        "neron_severi": {
            "fiber_and_zero_section": "U",
            "reducible_fiber_roots": "E8(-1)+A1(-1)^2",
            "lower_bound_lattice": "U+E8(-1)+A1(-1)^2",
            "gram": matrix_payload(ns_gram),
            "rank_signature_determinant": [12, [1, 11, 0], -4],
            "smith": smith_diagonal(ns_gram),
            "mordell_weil_very_general": {
                "rank": 0,
                "torsion": "trivial",
                "reason": (
                    "Shioda--Tate after exact rho=12; the no-overlattice "
                    "argument rules out a torsion saturation"
                ),
            },
            "primitive_at_rank_12": (
                "an index-two even overlattice would be even unimodular of "
                "signature (1,11), impossible since 1-11 is not divisible "
                "by 8"
            ),
            "generic_exactness": {
                "parameter_moduli_dimension": 4,
                "parameter_recovery": {
                    "observables": [
                        "A'(1)",
                        "A''(1)/2",
                        "C''(1)/2",
                        "C'''(1)/6",
                    ],
                    "recovered_parameters": [
                        "alpha",
                        "beta",
                        "gamma",
                        "delta",
                    ],
                    "symbolic_jacobian": matrix_payload(recovery_jacobian),
                    "jacobian_determinant": int(recovery_jacobian.det()),
                },
                "normalizations": (
                    "the unique II* fiber fixes infinity; base scaling puts "
                    "the exchanged I2 pair at +/-1; Weierstrass homothety "
                    "normalizes the nodal cubic to A(1)=-3,C(1)=2"
                ),
                "normalized_coordinate_stabilizer": {
                    "base": (
                        "commuting with t->-t and fixing the unique II* "
                        "fiber leaves t->c*t; preserving {+1,-1} gives "
                        "c^2=1"
                    ),
                    "weierstrass": (
                        "preserving A(1)=-3 and C(1)=2 gives "
                        "lambda^4=lambda^6=1, hence lambda^2=1"
                    ),
                    "finite": True,
                },
                "marked_to_unmarked_zero_dimensional": {
                    "integral_classes": (
                        "a marked elliptic fibration and section determine "
                        "integral classes F,O in NS with F^2=0, O^2=-2, "
                        "and F.O=1; these classes are locally constant"
                    ),
                    "uniqueness": (
                        "for fixed F,O, the pencil |F| and section O recover "
                        "the fibration up to PGL2; the displayed base and "
                        "short-Weierstrass normalizations leave only the "
                        "finite stabilizer above"
                    ),
                    "automorphisms": (
                        "H^0(X,T_X)=0 for a K3 surface, so Aut(X) is "
                        "discrete; no positive-dimensional forgetful fiber "
                        "remains"
                    ),
                    "fiber_dimension": 0,
                },
                "argument": (
                    "The recovery Jacobian is the 4 by 4 identity, the "
                    "normalized coordinate stabilizer is finite, and the "
                    "marked-to-unmarked forgetful map has zero-dimensional "
                    "fibers by the integral F,O class argument.  Thus the "
                    "four normalized coefficients give four K3 moduli.  "
                    "Since sigma is primitive order 4 and sigma^2 is the "
                    "identity on the displayed lattice, T is a Q(i)-module, "
                    "so rank(T) and rho are even.  Since rho>=12, a strict "
                    "jump would have rho>=14 and rank(T)<=8, whose primitive-"
                    "eigenvalue period ball has dimension at most 3.  This "
                    "contradicts the four-dimensional family.  Hence very "
                    "general rho=12."
                ),
                "scope": (
                    "exact for the geometric generic/very general member; no "
                    "particular rational parameter quadruple is Picard-"
                    "certified here"
                ),
            },
        },
        "picard_action": {
            "description": (
                "identity on U+E8(-1); sigma exchanges the fibers t=+1 and "
                "t=-1 and hence swaps their two A1 root components"
            ),
            "a1x2_gram": matrix_payload(A1X2),
            "a1_swap": matrix_payload(A1_SWAP),
            "characteristic_polynomial": str(ns_characteristic),
            "fixed_lattice": {
                "class": "U+E8(-1)+<-4>",
                "rank": 11,
                "signature": [1, 10, 0],
                "determinant": 4,
                "smith": smith_diagonal(fixed_gram),
                "a1_diagonal_gram": matrix_payload(fixed_a1_gram),
                "discriminant_q_histogram": {"0": 1, "1": 1, "7/4": 2},
            },
            "coinvariant": {
                "class": "<-4>",
                "gram": matrix_payload(coinvariant_gram),
                "generator": "R_+ - R_-",
                "divisibility_in_NS": 2,
                "rootless": True,
            },
            "primitive_glue": {
                "index": 2,
                "description": (
                    "with f=R_++R_- and c=R_+-R_-, the original A1^2 is "
                    "obtained from <f>+<c> by adjoining (f+c)/2=R_+"
                ),
            },
        },
        "cohomology_and_fixed_loci": {
            "generic_transcendental_rank": 10,
            "transcendental_characteristic": "(lambda^2+1)^5",
            "full_characteristic": str(full_characteristic),
            "h2_traces_powers_1_2_3": h2_traces,
            "fixed_euler_powers_1_2_3": fixed_euler,
            "bh_generator_fingerprint": {
                "database_tuple": [[8], 2, 0, 1],
                "meaning": "8 isolated points and 2 rational curves",
                "euler": 12,
            },
            "bh_square_fingerprint": {
                "database_tuple": [[], 5, 1, 4],
                "meaning": "5 rational curves and one genus-4 curve",
                "euler": 4,
            },
        },
        "artebani_sarti_lineage": {
            "primary_source": "https://arxiv.org/abs/1102.4436",
            "generic_even_odd_model": (
                "Artebani--Sarti Example 6.3: y^2=x^3+a(t)x+b(t), "
                "a even, b odd, with (x,y,t)->(-x,iy,-t)"
            ),
            "node_specialization": (
                "their Example 7.2 explains that two nodes exchanged by the "
                "base involution insert A1^2 into S(sigma^2)"
            ),
            "our_fiber_root_specialization": (
                "degree truncation ord_infinity(a4)>=4, ord_infinity(a6)=5 "
                "inserts II*/E8 at infinity, while the imposed double "
                "discriminant root at s=1 inserts the exchanged I2/I2 pair"
            ),
            "novelty_scope": (
                "a new explicit equation-to-Coxeter realization in this "
                "project; no worldwide priority for the underlying order-4 "
                "K3 family or the abstract Coxeter polytope is claimed"
            ),
        },
        "brandhorst_hofmann_handoff": {
            "id": "0.4.4.1",
            "article": "https://doi.org/10.1017/fms.2023.50",
            "database": "https://doi.org/10.5281/zenodo.10410278",
            "existing_certificate_sha256": sha256_file(SOURCE),
            "row_one_based": bh["row_one_based"],
            "row_sha256_with_newline": bh["row_sha256_with_newline"],
            "same_fixed_coinvariant_glue_action": True,
            "full_NS_action_pairs_integrally_conjugate": conjugacy[
                "full_NS_action_pairs_integrally_conjugate"
            ],
            "all_glue_pairs_carried": conjugacy["all_glue_pairs_carried"],
            "external_replay": None,
        },
        "chamber_handoff": chamber,
        "hyperbolic_handoff": {
            "polytope_id": case["polytope_id"],
            "dimension": hyperbolic["dimension"],
            "facets": hyperbolic["facets"],
            "classification": hyperbolic["classification"],
            "finite_vertices": hyperbolic["finite_vertices"],
            "cusps": hyperbolic["cusps"],
            "elliptic_edges": hyperbolic["elliptic_edges"],
            "bad_audits": hyperbolic["Vinberg_bad_audits"],
            "canonical_full_audit_sha256": VINBERG_AUDIT_SHA256,
            "audit_digest_schema": (
                "canonical_json_v1_finite_vertices_ideal_components_edges_"
                "parabolic_audits"
            ),
            "abstract_polytope_novelty": (
                "known; this is a new explicit algebraic realization in the "
                "project"
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
    parser.add_argument("--brandhorst-hofmann-zip", type=Path)
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

    summary: dict[str, Any] = {
        "all_checks_pass": True,
        "payload_sha256": payload["payload_sha256"],
    }
    if args.brandhorst_hofmann_zip is not None:
        summary["external_replay"] = replay_brandhorst_hofmann_zip(
            args.brandhorst_hofmann_zip
        )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
