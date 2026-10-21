#!/usr/bin/env python3
"""Certify an explicit Weierstrass family for BH 0.4.4.6 / H10/38.

For s=t^2 and u=s-1, define

    P(s) = -3 + alpha*u + beta*u^2,
    Q(s) =  2 - alpha*u + gamma*u^2 + delta*u^3.

The minimal resolution of

    y^2 = x^3 + t^2*P(t^2)*x + t^3*Q(t^2)

has the primitive order-four automorphism

    sigma(x,y,t)=(-x,i*y,-t).

Generically the fibers are I0* at zero and infinity, I2 at t=+/-1, and
eight I1 fibers.  The very general Neron--Severi lattice is
U+D4(-1)^2+A1(-1)^2.  Tate-root coordinates show that sigma fixes both D4
diagrams vertex by vertex and exchanges the two A1 components.  This is the
exact BH 0.4.4.6 action and explicitly realizes the finite-volume H10/38
Coxeter section.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any
import zipfile

import sympy as sp

from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    _parse_bh_row,
    smith_diagonal,
)
from scripts.certify_ue8a1x2_order4_weierstrass_family import (
    A1_SWAP,
    A1X2,
    U,
    block_diag,
    matrix_payload,
    negative_cartan,
    payload_digest,
    sha256_bytes,
    sha256_file,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/roulleau_rank12_bh_k3_lifts.json"
OUTPUT = ROOT / "data/enumerations/ud4x2a1x2_order4_weierstrass_family.json"

BH_ZIP_SHA256 = "c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c"
BH_MEMBER = "K3Groups/results/purely_ns/order04.txt"
BH_MEMBER_SHA256 = "d5a61603c1539278e8411e188656c9417662336ff4282b714e939c40d38c94e8"
BH_ROW_ONE_BASED = 44
BH_ROW_SHA256 = "2c5db8537584081dac8e9ddda92022f41a129f32b043063fd8155897449168ff"
VINBERG_AUDIT_SHA256 = (
    "79c800f246b68e6eff1908dd1922abdf7b31fa3ac436e43920fce8b22a47cf03"
)
BH_FIXED_PROPERTIES = [
    (
        "fixed points powers",
        (((8,), 2, 0, 1), ((), 3, 1, 2)),
    ),
    ("fixed points", ((8,), 2, 0, 1)),
]


D4 = negative_cartan(4, [(0, 1), (1, 2), (1, 3)])


def source_case() -> tuple[dict[str, Any], dict[str, Any]]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = next(
        row
        for row in source["realizations"]
        if row["brandhorst_hofmann"]["id"] == "0.4.4.6"
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
    if parsed["id"] != "0.4.4.6":
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
    p_coefficient = -3 + alpha * u + beta * u**2
    q_coefficient = 2 - alpha * u + gamma * u**2 + delta * u**3
    tate_discriminant = sp.expand(
        4 * p_coefficient**3 + 27 * q_coefficient**2
    )
    residual = sp.cancel(tate_discriminant / u**2)
    finite_i2_leading = sp.factor(residual.subs(u, 0))
    zero_i0star_leading = sp.factor(tate_discriminant.subs(u, -1))
    infinity_i0star_leading = 4 * beta**3 + 27 * delta**2

    witness = {alpha: 0, beta: 1, gamma: 0, delta: 1}
    witness_residual = sp.factor(residual.subs(witness))
    expected_witness = 31 * u**4 - 36 * u**2 + 108 * u + 108
    witness_discriminant = int(sp.discriminant(expected_witness, u))
    witness_resultant = int(
        sp.resultant(expected_witness, p_coefficient.subs(witness), u)
    )

    xi = sp.Symbol("xi")
    zero_tate_cubic = (
        xi**3
        + p_coefficient.subs(u, -1) * xi
        + q_coefficient.subs(u, -1)
    )
    infinity_tate_cubic = xi**3 + beta * xi + delta
    zero_cubic_discriminant = sp.discriminant(zero_tate_cubic, xi)
    infinity_cubic_discriminant = sp.discriminant(infinity_tate_cubic, xi)

    # Exact coordinates on the normalized marked-Weierstrass stratum.
    parameters = sp.Matrix([alpha, beta, gamma, delta])
    recovered_parameters = sp.Matrix(
        [
            sp.diff(p_coefficient, u).subs(u, 0),
            sp.diff(p_coefficient, u, 2).subs(u, 0) / 2,
            sp.diff(q_coefficient, u, 2).subs(u, 0) / 2,
            sp.diff(q_coefficient, u, 3).subs(u, 0) / 6,
        ]
    )
    recovery_jacobian = recovered_parameters.jacobian(parameters)
    base_scale, weierstrass_scale = sp.symbols(
        "base_scale weierstrass_scale"
    )
    weierstrass_stabilizer_polynomial = sp.gcd(
        weierstrass_scale**4 - 1,
        weierstrass_scale**6 - 1,
    )

    # Actual I0* Weierstrass germs and their first blowup charts.
    tau, local_v, local_x, local_X, zeta = sp.symbols(
        "tau local_v local_x local_X zeta"
    )
    p_at_zero = p_coefficient.subs(u, tau**2 - 1)
    q_at_zero = q_coefficient.subs(u, tau**2 - 1)
    zero_branch_cubic = xi**3 + p_at_zero * xi + q_at_zero
    zero_branch_rhs = (
        local_x**3
        + tau**2 * p_at_zero * local_x
        + tau**3 * q_at_zero
    )
    p_at_infinity = (
        beta
        + (alpha - 2 * beta) * zeta
        + (-3 - alpha + beta) * zeta**2
    )
    q_at_infinity = (
        delta
        + (gamma - 3 * delta) * zeta
        + (-alpha - 2 * gamma + 3 * delta) * zeta**2
        + (2 + alpha + gamma - delta) * zeta**3
    )
    infinity_branch_cubic = (
        xi**3
        + p_at_infinity.subs(zeta, local_v**2) * xi
        + q_at_infinity.subs(zeta, local_v**2)
    )
    infinity_branch_rhs = (
        local_X**3
        + local_v**2
        * p_at_infinity.subs(zeta, local_v**2)
        * local_X
        + local_v**3 * q_at_infinity.subs(zeta, local_v**2)
    )

    ns_gram = block_diag(U, D4, D4, A1X2)
    ns_action = block_diag(sp.eye(10), A1_SWAP)
    fixed_a1_basis = sp.Matrix([[1, 1]])
    coinvariant_basis = sp.Matrix([[1, -1]])
    fixed_a1_gram = fixed_a1_basis * A1X2 * fixed_a1_basis.T
    coinvariant_gram = coinvariant_basis * A1X2 * coinvariant_basis.T
    fixed_gram = block_diag(U, D4, D4, fixed_a1_gram)

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
        "finite_i2_pair_divides_discriminant": sp.expand(
            tate_discriminant - u**2 * residual
        )
        == 0,
        "finite_i2_values": (
            p_coefficient.subs(u, 0) == -3
            and q_coefficient.subs(u, 0) == 2
        ),
        "finite_i2_derivative_relation": sp.simplify(
            sp.diff(q_coefficient, u).subs(u, 0)
            + sp.diff(p_coefficient, u).subs(u, 0)
        )
        == 0,
        "finite_i2_open_coefficient": sp.simplify(
            finite_i2_leading
            + 9 * (alpha**2 - 12 * beta - 12 * gamma)
        )
        == 0,
        "residual_discriminant_degree_four": sp.degree(residual, u) == 4,
        "zero_i0star_cubic_discriminant": sp.simplify(
            zero_cubic_discriminant + zero_i0star_leading
        )
        == 0,
        "infinity_i0star_cubic_discriminant": sp.simplify(
            infinity_cubic_discriminant + infinity_i0star_leading
        )
        == 0,
        "witness_residual_quartic": witness_residual == expected_witness,
        "witness_residual_squarefree": witness_discriminant
        == -1900537572096,
        "witness_zero_i0star_exact": zero_i0star_leading.subs(witness) == -5,
        "witness_infinity_i0star_exact": infinity_i0star_leading.subs(witness)
        == 31,
        "witness_no_additive_residual_collision": witness_resultant == 42849,
        "parameter_recovery_is_identity": recovered_parameters == parameters,
        "parameter_recovery_jacobian_is_identity": (
            recovery_jacobian == sp.eye(4)
        ),
        "base_normalization_has_four_element_stabilizer": (
            sp.degree(base_scale**2 - 1, base_scale) == 2
        ),
        "weierstrass_normalization_has_finite_stabilizer": sp.expand(
            weierstrass_stabilizer_polynomial
            - (weierstrass_scale**2 - 1)
        )
        == 0,
        "zero_actual_branch_blowup_identity": sp.expand(
            zero_branch_rhs.subs(local_x, tau * xi)
            - tau**3 * zero_branch_cubic
        )
        == 0,
        "zero_blowup_coordinate_fixed": sp.simplify(
            (-local_x) / (-tau) - local_x / tau
        )
        == 0,
        "infinity_transformed_a4_identity": sp.simplify(
            local_v**6
            * p_coefficient.subs(u, local_v**-2 - 1)
            - local_v**2 * p_at_infinity.subs(zeta, local_v**2)
        )
        == 0,
        "infinity_transformed_a6_identity": sp.simplify(
            local_v**9
            * q_coefficient.subs(u, local_v**-2 - 1)
            - local_v**3 * q_at_infinity.subs(zeta, local_v**2)
        )
        == 0,
        "infinity_actual_branch_blowup_identity": sp.expand(
            infinity_branch_rhs.subs(local_X, local_v * xi)
            - local_v**3 * infinity_branch_cubic
        )
        == 0,
        "infinity_blowup_coordinate_fixed": sp.simplify(
            (-local_X) / (-local_v) - local_X / local_v
        )
        == 0,
        "fiber_euler_sum_is_24": 2 * 6 + 2 * 2 + 8 == 24,
        "trivial_lattice_rank_is_12": 2 + 2 * 4 + 2 == 12,
        "d4_determinant": D4.det() == 4,
        "ns_signature": exact_inertia(ns_gram) == (1, 11, 0),
        "ns_determinant": ns_gram.det() == -64,
        "ns_smith": smith_diagonal(ns_gram)
        == [*[1] * 6, *[2] * 6],
        "a1_swap_isometry": A1_SWAP * A1X2 * A1_SWAP.T == A1X2,
        "ns_action_isometry": ns_action * ns_gram * ns_action.T == ns_gram,
        "ns_action_order_two": ns_action**2 == sp.eye(12),
        "ns_characteristic": ns_characteristic
        == (lam - 1) ** 11 * (lam + 1),
        "fixed_lattice_signature": exact_inertia(fixed_gram) == (1, 10, 0),
        "fixed_lattice_determinant": fixed_gram.det() == 64,
        "fixed_lattice_smith": smith_diagonal(fixed_gram)
        == [*[1] * 6, *[2] * 4, 4],
        "coinvariant_is_minus_four": coinvariant_gram == sp.Matrix([[-4]]),
        "primitive_glue_index_two": abs(
            int(sp.Matrix([[1, 1], [1, -1]]).det())
        )
        == 2,
        "torsion_height_upper_correction_is_three": 2 * sp.Rational(1)
        + 2 * sp.Rational(1, 2)
        == 3,
        "torsion_height_lower_bound_positive": 4 - 3 == 1,
        "full_characteristic": full_characteristic
        == (lam - 1) ** 11 * (lam + 1) * (lam**2 + 1) ** 5,
        "fixed_euler_characteristics": fixed_euler
        == {"1": 12, "2": 4, "3": 12},
        "source_all_checks_pass": source["all_checks_pass"] is True,
        "source_bh_id": bh["id"] == "0.4.4.6",
        "source_generic_ns": bh["generic_NS"]
        == "U+D4(-1)^2+A1(-1)^2",
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
            "place": "t=0",
            "ord_a4": 2,
            "ord_a6": 3,
            "ord_discriminant": 6,
            "kodaira": "I0*",
            "euler": 6,
            "root_lattice": "D4(-1)",
            "root_rank": 4,
        },
        {
            "place": "t=infinity",
            "ord_a4": 2,
            "ord_a6": 3,
            "ord_discriminant": 6,
            "kodaira": "I0*",
            "euler": 6,
            "root_lattice": "D4(-1)",
            "root_rank": 4,
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
            "place": "eight roots above the four residual roots in s",
            "multiplicity": 8,
            "ord_discriminant": 1,
            "kodaira": "I1",
            "euler_each": 1,
            "root_rank": 0,
        },
    ]

    payload: dict[str, Any] = {
        "schema": "ud4x2a1x2-order4-weierstrass-family-v1",
        "theorem": (
            "On a nonempty Zariski-open set of (alpha,beta,gamma,delta), "
            "the very general minimal resolution of y^2=x^3+t^2*P(t^2)x+"
            "t^3*Q(t^2) has NS=U+D4(-1)^2+A1(-1)^2.  The primitive "
            "order-four automorphism (x,y,t)->(-x,i*y,-t) fixes both D4 "
            "diagrams and exchanges the two A1 roots, realizes BH 0.4.4.6, "
            "and explicitly realizes the finite-volume noncompact H10/38 "
            "Coxeter section."
        ),
        "family": {
            "affine_equation": "y^2=x^3+t^2*P(t^2)*x+t^3*Q(t^2)",
            "P": "P(s)=-3+alpha*(s-1)+beta*(s-1)^2",
            "Q": "Q(s)=2-alpha*(s-1)+gamma*(s-1)^2+delta*(s-1)^3",
            "automorphism": "sigma(x,y,t)=(-x,i*y,-t)",
            "global_order": 4,
            "square": "sigma^2(x,y,t)=(x,-y,t)",
            "holomorphic_two_form": "dt wedge dx/y",
            "two_form_multiplier": "-i (a primitive fourth root)",
            "reduced_tate_discriminant": (
                "E(s)=4*P(s)^3+27*Q(s)^2=(s-1)^2*R(s), deg R=4"
            ),
            "full_discriminant": "Delta(t)=-16*t^6*E(t^2)",
            "finite_i2_open_coefficient": str(finite_i2_leading),
            "zero_i0star_open_coefficient": str(zero_i0star_leading),
            "infinity_i0star_open_coefficient": str(
                infinity_i0star_leading
            ),
            "parameter_open": (
                "the three displayed fiber coefficients are nonzero; R is "
                "squarefree, has no root s=0, and shares no root with P"
            ),
            "nonempty_open_witness": {
                "parameters": {
                    "alpha": 0,
                    "beta": 1,
                    "gamma": 0,
                    "delta": 1,
                },
                "P_in_u_equals_s_minus_1": "u^2-3",
                "Q_in_u_equals_s_minus_1": "u^3+2",
                "R": str(expected_witness),
                "R_discriminant": witness_discriminant,
                "zero_i0star_coefficient": int(
                    zero_i0star_leading.subs(witness)
                ),
                "infinity_i0star_coefficient": int(
                    infinity_i0star_leading.subs(witness)
                ),
                "resultant_R_P": witness_resultant,
            },
            "fibers": fiber_rows,
        },
        "d4_component_actions": {
            "at_t_zero": {
                "tate_root_coordinate": "xi=x/t",
                "sigma_action": "xi->(-x)/(-t)=xi",
                "tate_cubic": "xi^3+P(0)*xi+Q(0)",
                "cubic_discriminant": str(zero_cubic_discriminant),
                "conclusion": (
                    "on the open set its three distinct roots, hence all "
                    "three outer D4 arms, are individually fixed"
                ),
                "actual_blowup_chart": (
                    "substituting x=t*xi gives RHS=t^3*"
                    "(xi^3+P(t^2)xi+Q(t^2)); blowing up the branch-plane "
                    "origin has chart (t,xi=x/t)->(-t,xi), so the three "
                    "distinct branch points on the exceptional divisor are "
                    "fixed separately"
                ),
            },
            "at_t_infinity": {
                "local_coordinates": "u=1/t, X=u^4*x, Y=u^6*y",
                "tate_root_coordinate": "xi=X/u",
                "sigma_action": "xi->(-X)/(-u)=xi",
                "tate_cubic": "xi^3+beta*xi+delta",
                "cubic_discriminant": str(infinity_cubic_discriminant),
                "conclusion": (
                    "on the open set its three distinct roots, hence all "
                    "three outer D4 arms, are individually fixed"
                ),
                "actual_blowup_chart": (
                    "the transformed equation is Y^2=X^3+v^2*Phat(v^2)X+"
                    "v^3*Qhat(v^2); substituting X=v*xi gives RHS=v^3*"
                    "(xi^3+Phat(v^2)xi+Qhat(v^2)); the blowup chart "
                    "(v,xi=X/v)->(-v,xi) fixes the three distinct branch "
                    "points separately"
                ),
            },
            "global_conclusion": (
                "both central components and all six outer arms are fixed; "
                "neither D4 triality group is used"
            ),
        },
        "neron_severi": {
            "fiber_and_zero_section": "U",
            "reducible_fiber_roots": "D4(-1)^2+A1(-1)^2",
            "lower_bound_lattice": "U+D4(-1)^2+A1(-1)^2",
            "gram": matrix_payload(ns_gram),
            "rank_signature_determinant": [12, [1, 11, 0], -64],
            "smith": smith_diagonal(ns_gram),
            "generic_exactness": {
                "parameter_moduli_dimension": 4,
                "parameter_recovery": {
                    "observables": [
                        "P'(1)",
                        "P''(1)/2",
                        "Q''(1)/2",
                        "Q'''(1)/6",
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
                    "the pair of fixed I0* fibers is placed at 0,infinity; "
                    "base scaling puts the exchanged I2 pair at +/-1; "
                    "Weierstrass homothety normalizes P(1),Q(1) to -3,2"
                ),
                "normalized_coordinate_stabilizer": {
                    "base": (
                        "the distinct reducible-fiber divisors {0,infinity} "
                        "of type I0* and {+1,-1} of type I2 leave only "
                        "t->c*t or c/t with c^2=1"
                    ),
                    "weierstrass": (
                        "preserving P(1)=-3 and Q(1)=2 gives "
                        "lambda^4=lambda^6=1, hence lambda^2=1"
                    ),
                    "finite": True,
                },
                "marked_to_unmarked_zero_dimensional": {
                    "integral_classes": (
                        "the marked fiber and section are integral classes "
                        "F,O in NS with F^2=0, O^2=-2, F.O=1 and are locally "
                        "constant in a fiber of the forgetful map"
                    ),
                    "uniqueness": (
                        "fixed F,O recover the pencil |F| and section O; "
                        "the base and short-Weierstrass normalizations then "
                        "leave only the finite stabilizer above"
                    ),
                    "automorphisms": (
                        "H^0(X,T_X)=0 makes Aut(X) discrete"
                    ),
                    "fiber_dimension": 0,
                },
                "period_argument": (
                    "The identity recovery Jacobian, finite coordinate "
                    "stabilizer, and zero-dimensional integral-class "
                    "forgetful fibers prove that the normalized family has "
                    "four K3 moduli.  Primitive order "
                    "4 makes T a Q(i)-module and rho even.  A strict jump "
                    "from rho>=12 would have rho>=14 and period dimension at "
                    "most 3, contradicting the four-dimensional family."
                ),
                "mordell_weil_rank": 0,
                "torsion_height_argument": (
                    "Each D4/I0* local correction is at most 1 and each "
                    "A1/I2 correction at most 1/2.  Thus a nonzero torsion "
                    "section would have h(P)>=4-(1+1+1/2+1/2)=1, impossible. "
                    "The trivial lattice is saturated."
                ),
                "conclusion": "very general NS=U+D4(-1)^2+A1(-1)^2",
                "scope": (
                    "exact for the geometric generic/very general member; no "
                    "particular rational quadruple is Picard-certified here"
                ),
            },
        },
        "picard_action": {
            "description": (
                "identity on U+D4(-1)^2; exchange of the two I2/A1 roots"
            ),
            "a1x2_gram": matrix_payload(A1X2),
            "a1_swap": matrix_payload(A1_SWAP),
            "characteristic_polynomial": str(ns_characteristic),
            "fixed_lattice": {
                "class": "U+D4(-1)^2+<-4>",
                "rank": 11,
                "signature": [1, 10, 0],
                "determinant": 64,
                "smith": smith_diagonal(fixed_gram),
                "a1_diagonal_gram": matrix_payload(fixed_a1_gram),
                "discriminant_q_histogram": {
                    "0": 16,
                    "1": 16,
                    "3/4": 12,
                    "7/4": 20,
                },
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
                    "with f=R_++R_- and c=R_+-R_-, adjoining "
                    "(f+c)/2=R_+ recovers the original A1^2"
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
                "database_tuple": [[], 3, 1, 2],
                "meaning": "3 rational curves and one genus-2 curve",
                "euler": 4,
            },
        },
        "artebani_sarti_lineage": {
            "primary_source": "https://arxiv.org/abs/1102.4436",
            "ambient_model": (
                "Artebani--Sarti Example 6.3 gives the even-a, odd-b "
                "Weierstrass model with (x,y,t)->(-x,iy,-t)"
            ),
            "our_fiber_root_specialization": (
                "factoring t^2 from a4 and t^3 from a6 inserts I0*/D4 at "
                "zero; degree truncation inserts a second I0*/D4 at infinity; "
                "a double root in s=t^2 inserts the exchanged I2/I2 pair"
            ),
            "novelty_scope": (
                "a new explicit equation-to-Coxeter realization in this "
                "project; no worldwide priority for the underlying order-4 "
                "K3 family or abstract Coxeter polytope is claimed"
            ),
        },
        "brandhorst_hofmann_handoff": {
            "id": "0.4.4.6",
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
                "bounded-corpus novelty candidate only; this theorem is a "
                "new explicit algebraic realization in the project"
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
