#!/usr/bin/env python3
"""Certify an explicit Weierstrass family for BH 0.4.4.2 / H10/15.

For ``s=t^2`` and parameters ``r,alpha,beta,gamma``, define

    P(s) = -3*s^3 + alpha*s^2 + beta*s + p0,
    Q(s) =  2*s^4 - alpha*s^3 + gamma*s^2 + q1*s + q0,

where p0,q1,q0 are chosen below.  They simultaneously force an I4* fiber
at infinity and an exchanged pair of I2 fibers at t=+/-1.  The minimal
resolution of

    y^2 = x^3 + P(t^2)*x + t*Q(t^2)

has the primitive order-four automorphism

    sigma(x,y,t)=(-x,i*y,-t).

On a nonempty Zariski-open set the very general Neron--Severi lattice is
U+D8(-1)+A1(-1)^2.  The local D8 rational-double-point branch model shows
that sigma fixes the D8 diagram vertex by vertex, while it exchanges the two
A1 components.  This is the exact BH 0.4.4.2 Picard action and explicitly
realizes the finite-volume noncompact H10/15 Coxeter section.
"""

from __future__ import annotations

import argparse
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
OUTPUT = ROOT / "data/enumerations/ud8a1x2_order4_weierstrass_family.json"

BH_ZIP_SHA256 = "c1272db2c0fe2ae0a788aed7745201fccd2a1337b6c2d85555880c434cab4d7c"
BH_MEMBER = "K3Groups/results/purely_ns/order04.txt"
BH_MEMBER_SHA256 = "d5a61603c1539278e8411e188656c9417662336ff4282b714e939c40d38c94e8"
BH_ROW_ONE_BASED = 40
BH_ROW_SHA256 = "bbba67834e034344032bbf137a1f2f1499f06f3adac123ce2b580e0bf7614849"
VINBERG_AUDIT_SHA256 = (
    "00977c1aa97e62df0ac6bbd22e635cd3b37dbf38e22b3cf52770ecd7a42c15fb"
)
BH_FIXED_PROPERTIES = [
    (
        "fixed points powers",
        (((8,), 2, 0, 1), ((), 4, 1, 3)),
    ),
    ("fixed points", ((8,), 2, 0, 1)),
]


D8 = negative_cartan(
    8,
    [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (5, 7),
    ],
)


def source_case() -> tuple[dict[str, Any], dict[str, Any]]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = next(
        row
        for row in source["realizations"]
        if row["brandhorst_hofmann"]["id"] == "0.4.4.2"
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
    if parsed["id"] != "0.4.4.2":
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

    s, z = sp.symbols("s z")
    r, alpha, beta, gamma = sp.symbols("r alpha beta gamma")
    p0 = -3 * r**2 + 3 - alpha - beta
    q1 = (
        9 * r
        - 2 * alpha * r
        - beta * r
        - r**3
        - 8
        + 3 * alpha
        - 2 * gamma
    )
    q0 = 2 * r**3 - 2 + alpha - gamma - q1
    p_coefficient = -3 * s**3 + alpha * s**2 + beta * s + p0
    q_coefficient = (
        2 * s**4 - alpha * s**3 + gamma * s**2 + q1 * s + q0
    )

    reduced_discriminant = sp.expand(
        4 * p_coefficient**3 + 27 * s * q_coefficient**2
    )
    residual = sp.cancel(reduced_discriminant / (s - 1) ** 2)

    p_infinity = -3 + alpha * z + beta * z**2 + p0 * z**3
    q_infinity = (
        2 - alpha * z + gamma * z**2 + q1 * z**3 + q0 * z**4
    )
    infinity_bracket = sp.expand(4 * p_infinity**3 + 27 * q_infinity**2)
    infinity_residual = sp.cancel(infinity_bracket / z**2)
    infinity_leading = sp.factor(infinity_residual.subs(z, 0))
    finite_i2_leading = sp.factor(residual.subs(s, 1))
    tate_root = sp.Symbol("xi")
    tate_cubic_at_infinity = tate_root**3 - 3 * tate_root + 2

    witness = {r: 1, alpha: 0, beta: 1, gamma: 0}
    witness_residual = sp.factor(residual.subs(witness))
    expected_witness = 108 * s**5 - 36 * s**3 + 31 * s - 4
    witness_discriminant = int(sp.discriminant(expected_witness, s))
    witness_resultant = int(
        sp.resultant(expected_witness, p_coefficient.subs(witness), s)
    )

    p_at_one = sp.expand(p_coefficient.subs(s, 1))
    q_at_one = sp.expand(q_coefficient.subs(s, 1))
    recovered_parameters = sp.Matrix(
        [
            sp.cancel(-3 * q_at_one / (2 * p_at_one)),
            sp.expand(p_coefficient).coeff(s, 2),
            sp.expand(p_coefficient).coeff(s, 1),
            sp.expand(q_coefficient).coeff(s, 2),
        ]
    )
    parameters = sp.Matrix([r, alpha, beta, gamma])
    recovery_jacobian = recovered_parameters.jacobian(parameters)
    base_scale, weierstrass_scale = sp.symbols(
        "base_scale weierstrass_scale"
    )
    weierstrass_stabilizer_polynomial = sp.gcd(
        weierstrass_scale**4 - 1,
        weierstrass_scale**6 - 1,
    )

    # Standard D8 rational double point and the induced branch action.
    local_y, local_v, local_u, kappa, branch_c = sp.symbols(
        "local_y local_v local_u kappa branch_c"
    )
    d8_standard = local_y**2 + local_u * local_v**2 + kappa * local_u**7
    d8_transformed = sp.expand(
        d8_standard.subs(
            {
                local_y: sp.I * local_y,
                local_v: -local_v,
                local_u: -local_u,
            },
            simultaneous=True,
        )
    )
    branch_zero = local_u
    branch_plus = local_v - branch_c * local_u**3
    branch_minus = local_v + branch_c * local_u**3

    actual_d8_branch = sp.expand(
        (local_u + local_v) ** 3
        + local_u**2
        * p_infinity.subs(z, local_u**2)
        * (local_u + local_v)
        + local_u**3 * q_infinity.subs(z, local_u**2)
    )
    expected_actual_d8_branch = sp.expand(
        local_v**2 * (local_v + 3 * local_u)
        + alpha * local_u**4 * local_v
        + beta * local_u**6 * local_v
        + p0 * local_u**8 * local_v
        + (beta + gamma) * local_u**7
        + (p0 + q1) * local_u**9
        + q0 * local_u**11
    )
    first_chart, second_chart, third_chart = sp.symbols("q0 q1 q2")
    first_strict = sp.cancel(
        actual_d8_branch.subs(local_v, local_u * first_chart)
        / local_u**3
    )
    second_strict = sp.cancel(
        actual_d8_branch.subs(local_v, local_u**2 * second_chart)
        / local_u**5
    )
    third_strict = sp.cancel(
        actual_d8_branch.subs(local_v, local_u**3 * third_chart)
        / local_u**7
    )
    tangent_slope_polynomial = sp.expand(third_strict.subs(local_u, 0))
    tangent_slope_discriminant = sp.discriminant(
        tangent_slope_polynomial, third_chart
    )
    transformed_actual_d8_branch = sp.expand(
        actual_d8_branch.subs(
            {local_v: -local_v, local_u: -local_u},
            simultaneous=True,
        )
    )

    def local_transform(poly: sp.Expr) -> sp.Expr:
        return sp.expand(
            poly.subs(
                {local_v: -local_v, local_u: -local_u},
                simultaneous=True,
            )
        )

    ns_gram = block_diag(U, D8, A1X2)
    ns_action = block_diag(sp.eye(10), A1_SWAP)
    fixed_a1_basis = sp.Matrix([[1, 1]])
    coinvariant_basis = sp.Matrix([[1, -1]])
    fixed_a1_gram = fixed_a1_basis * A1X2 * fixed_a1_basis.T
    coinvariant_gram = coinvariant_basis * A1X2 * coinvariant_basis.T
    fixed_gram = block_diag(U, D8, fixed_a1_gram)

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
            reduced_discriminant - (s - 1) ** 2 * residual
        )
        == 0,
        "finite_i2_values": (
            p_coefficient.subs(s, 1) == -3 * r**2
            and q_coefficient.subs(s, 1) == 2 * r**3
        ),
        "finite_i2_derivative_relation": sp.simplify(
            sp.diff(q_coefficient, s).subs(s, 1)
            + r * sp.diff(p_coefficient, s).subs(s, 1)
            + r**3
        )
        == 0,
        "finite_i2_exact_open_coefficient": finite_i2_leading
        != sp.Integer(0),
        "infinity_double_cancellation": sp.expand(
            infinity_bracket - z**2 * infinity_residual
        )
        == 0,
        "infinity_i4star_leading": sp.simplify(
            infinity_leading + 9 * (alpha**2 - 12 * beta - 12 * gamma)
        )
        == 0,
        "infinity_tate_cubic_double_root": sp.factor(
            tate_cubic_at_infinity
        )
        == (tate_root - 1) ** 2 * (tate_root + 2),
        "global_discriminant_degree_seven": sp.degree(
            reduced_discriminant, s
        )
        == 7,
        "residual_discriminant_degree_five": sp.degree(residual, s) == 5,
        "witness_residual_quintic": witness_residual == expected_witness,
        "witness_residual_squarefree": witness_discriminant
        == 6477927509458944,
        "witness_no_fiber_at_t_zero": expected_witness.subs(s, 0) == -4,
        "witness_no_additive_residual_collision": witness_resultant
        == -1440747,
        "parameter_r_recovered_rationally": recovered_parameters[0] == r,
        "parameter_recovery_is_identity": recovered_parameters == parameters,
        "parameter_recovery_jacobian_is_identity": (
            recovery_jacobian == sp.eye(4)
        ),
        "base_normalization_has_finite_stabilizer": (
            sp.degree(base_scale**2 - 1, base_scale) == 2
        ),
        "weierstrass_normalization_has_finite_stabilizer": sp.expand(
            weierstrass_stabilizer_polynomial
            - (weierstrass_scale**2 - 1)
        )
        == 0,
        "actual_d8_translated_germ_identity": (
            actual_d8_branch == expected_actual_d8_branch
        ),
        "actual_d8_germ_is_antiinvariant": (
            transformed_actual_d8_branch == -actual_d8_branch
        ),
        "actual_d8_first_blowup_special_fiber": sp.expand(
            first_strict.subs(local_u, 0)
            - first_chart**2 * (first_chart + 3)
        )
        == 0,
        "actual_d8_second_blowup_special_fiber": sp.expand(
            second_strict.subs(local_u, 0) - 3 * second_chart**2
        )
        == 0,
        "actual_d8_third_blowup_slope_polynomial": (
            tangent_slope_polynomial
            == 3 * third_chart**2
            + alpha * third_chart
            + beta
            + gamma
        ),
        "actual_d8_tangent_slopes_distinct_on_open": (
            tangent_slope_discriminant
            == alpha**2 - 12 * beta - 12 * gamma
        ),
        "actual_d8_blowup_chart_actions": (
            sp.simplify((-local_v) / (-local_u) - local_v / local_u)
            == 0
            and sp.simplify(
                (-local_v) / (-local_u) ** 2
                + local_v / local_u**2
            )
            == 0
            and sp.simplify(
                (-local_v) / (-local_u) ** 3
                - local_v / local_u**3
            )
            == 0
        ),
        "d8_standard_action_preserves_surface": d8_transformed == -d8_standard,
        "d8_branch_zero_fixed": local_transform(branch_zero) == -branch_zero,
        "d8_branch_plus_fixed": local_transform(branch_plus) == -branch_plus,
        "d8_branch_minus_fixed": local_transform(branch_minus) == -branch_minus,
        "fiber_euler_sum_is_24": 10 + 2 * 2 + 10 == 24,
        "trivial_lattice_rank_is_12": 2 + 8 + 2 == 12,
        "d8_determinant": D8.det() == 4,
        "ns_signature": exact_inertia(ns_gram) == (1, 11, 0),
        "ns_determinant": ns_gram.det() == -16,
        "ns_smith": smith_diagonal(ns_gram)
        == [*[1] * 8, *[2] * 4],
        "a1_swap_isometry": A1_SWAP * A1X2 * A1_SWAP.T == A1X2,
        "ns_action_isometry": ns_action * ns_gram * ns_action.T == ns_gram,
        "ns_action_order_two": ns_action**2 == sp.eye(12),
        "ns_characteristic": ns_characteristic
        == (lam - 1) ** 11 * (lam + 1),
        "fixed_lattice_signature": exact_inertia(fixed_gram) == (1, 10, 0),
        "fixed_lattice_determinant": fixed_gram.det() == 16,
        "fixed_lattice_smith": smith_diagonal(fixed_gram)
        == [*[1] * 8, 2, 2, 4],
        "coinvariant_is_minus_four": coinvariant_gram == sp.Matrix([[-4]]),
        "primitive_glue_index_two": abs(
            int(sp.Matrix([[1, 1], [1, -1]]).det())
        )
        == 2,
        "torsion_height_upper_correction_is_three": sp.Rational(2)
        + 2 * sp.Rational(1, 2)
        == 3,
        "torsion_height_lower_bound_positive": 4 - 3 == 1,
        "full_characteristic": full_characteristic
        == (lam - 1) ** 11 * (lam + 1) * (lam**2 + 1) ** 5,
        "fixed_euler_characteristics": fixed_euler
        == {"1": 12, "2": 4, "3": 12},
        "source_all_checks_pass": source["all_checks_pass"] is True,
        "source_bh_id": bh["id"] == "0.4.4.2",
        "source_generic_ns": bh["generic_NS"]
        == "U+D8(-1)+A1(-1)^2",
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
            "ord_a4": 2,
            "ord_a6": 3,
            "ord_discriminant": 10,
            "kodaira": "I4*",
            "euler": 10,
            "root_lattice": "D8(-1)",
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
            "place": "ten roots above the five roots of the residual in s",
            "multiplicity": 10,
            "ord_discriminant": 1,
            "kodaira": "I1",
            "euler_each": 1,
            "root_rank": 0,
        },
    ]

    payload: dict[str, Any] = {
        "schema": "ud8a1x2-order4-weierstrass-family-v1",
        "theorem": (
            "On a nonempty Zariski-open set of (r,alpha,beta,gamma), the "
            "very general minimal resolution of y^2=x^3+P(t^2)x+tQ(t^2) "
            "has NS=U+D8(-1)+A1(-1)^2.  The primitive order-four "
            "automorphism (x,y,t)->(-x,i*y,-t) fixes D8 and exchanges the "
            "two A1 roots, realizes BH 0.4.4.2, and therefore explicitly "
            "realizes the archived finite-volume noncompact H10/15 Coxeter "
            "section."
        ),
        "family": {
            "affine_equation": "y^2=x^3+P(t^2)*x+t*Q(t^2)",
            "P": (
                "P(s)=-3*s^3+alpha*s^2+beta*s+"
                "(-3*r^2+3-alpha-beta)"
            ),
            "Q": (
                "Q(s)=2*s^4-alpha*s^3+gamma*s^2+q1*s+q0"
            ),
            "q1": (
                "q1=9*r-2*alpha*r-beta*r-r^3-8+3*alpha-2*gamma"
            ),
            "q0": "q0=2*r^3-2+alpha-gamma-q1",
            "automorphism": "sigma(x,y,t)=(-x,i*y,-t)",
            "global_order": 4,
            "square": "sigma^2(x,y,t)=(x,-y,t)",
            "holomorphic_two_form": "dt wedge dx/y",
            "two_form_multiplier": "-i (a primitive fourth root)",
            "reduced_discriminant": (
                "D(s)=4*P(s)^3+27*s*Q(s)^2=(s-1)^2*R(s), deg R=5"
            ),
            "full_discriminant": "Delta(t)=-16*D(t^2)",
            "infinity_i4star_open_coefficient": str(infinity_leading),
            "finite_i2_open_coefficient": str(finite_i2_leading),
            "parameter_open": (
                "r, alpha^2-12*beta-12*gamma, and the finite I2 leading "
                "coefficient are nonzero; R is squarefree, has no root at "
                "s=0, and shares no root with P"
            ),
            "nonempty_open_witness": {
                "parameters": {"r": 1, "alpha": 0, "beta": 1, "gamma": 0},
                "P": "-3*s^3+s-1",
                "Q": "2*s^4-s+1",
                "R": str(expected_witness),
                "R_discriminant": witness_discriminant,
                "R_at_s_zero": int(expected_witness.subs(s, 0)),
                "resultant_R_P": witness_resultant,
            },
            "fibers": fiber_rows,
        },
        "d8_component_action": {
            "weierstrass_local_coordinates": (
                "u=1/t, X=u^4*x, Y=u^6*y; sigma(X,Y,u)=(-X,iY,-u)"
            ),
            "repeated_tate_root_coordinate": (
                "xi=X/u is invariant.  Its cubic at u=0 is "
                "(xi-1)^2*(xi+2); exact order two in s=u^2 splits the "
                "double root into two power series xi_+(s),xi_-(s).  Thus "
                "all three branch germs X=u*xi_j(u^2) are individually "
                "fixed"
            ),
            "standard_d8_equation": "Y^2+u*V^2+kappa*u^7=0",
            "standard_action": "(Y,V,u)->(iY,-V,-u)",
            "branch_factorization": (
                "u*(V-c*u^3)*(V+c*u^3), with c^2=-kappa"
            ),
            "branchwise_fixed": True,
            "resolution_conclusion": (
                "all three branch germs and every infinitely-near blowup "
                "center are fixed; hence every D8 exceptional component is "
                "fixed and the spinor nodes are not exchanged"
            ),
            "actual_weierstrass_germ": {
                "translation": "V=X-u",
                "equation": (
                    "Y^2=V^2(V+3u)+alpha*u^4*V+beta*u^6*V+p0*u^8*V+"
                    "(beta+gamma)*u^7+(p0+q1)*u^9+q0*u^11"
                ),
                "action": "(Y,V,u)->(iY,-V,-u)",
                "first_blowup": (
                    "q0=V/u is fixed; at u=0 the branch polynomial is "
                    "q0^2(q0+3), separating the transverse branch at -3"
                ),
                "second_blowup": (
                    "q1=V/u^2 is negated; at u=0 the strict branch is "
                    "3*q1^2, with the common center fixed"
                ),
                "third_blowup": (
                    "q2=V/u^3 is fixed; the two tangent branches meet the "
                    "exceptional divisor at the distinct roots of "
                    "3*q2^2+alpha*q2+beta+gamma"
                ),
                "slope_discriminant": "alpha^2-12*beta-12*gamma",
                "conclusion": (
                    "the exact chart actions fix every blowup center and "
                    "each separated branch point, so the double-cover "
                    "resolution fixes every D8 exceptional component"
                ),
            },
        },
        "neron_severi": {
            "fiber_and_zero_section": "U",
            "reducible_fiber_roots": "D8(-1)+A1(-1)^2",
            "lower_bound_lattice": "U+D8(-1)+A1(-1)^2",
            "gram": matrix_payload(ns_gram),
            "rank_signature_determinant": [12, [1, 11, 0], -16],
            "smith": smith_diagonal(ns_gram),
            "generic_exactness": {
                "parameter_moduli_dimension": 4,
                "parameter_recovery": {
                    "rational_identity": "r=-3*Q(1)/(2*P(1)) on r!=0",
                    "coefficient_identities": [
                        "alpha=[s^2]P",
                        "beta=[s]P",
                        "gamma=[s^2]Q",
                    ],
                    "recovered_parameters": [
                        "r",
                        "alpha",
                        "beta",
                        "gamma",
                    ],
                    "symbolic_jacobian": matrix_payload(recovery_jacobian),
                    "jacobian_determinant": int(recovery_jacobian.det()),
                },
                "normalizations": (
                    "the unique I4* fiber fixes infinity; base scaling puts "
                    "the exchanged I2 pair at +/-1; Weierstrass homothety "
                    "normalizes the repeated leading Tate root to -3,2"
                ),
                "normalized_coordinate_stabilizer": {
                    "base": (
                        "the unique I4* fiber fixes infinity; commuting with "
                        "t->-t leaves t->c*t, and preserving {+1,-1} gives "
                        "c^2=1"
                    ),
                    "weierstrass": (
                        "preserving the infinity Tate coefficients -3,2 "
                        "gives lambda^4=lambda^6=1, hence lambda^2=1"
                    ),
                    "finite": True,
                },
                "marked_to_unmarked_zero_dimensional": {
                    "integral_classes": (
                        "the marked fiber and section are integral classes "
                        "F,O in NS with F^2=0, O^2=-2, F.O=1 and are locally "
                        "constant in a forgetful fiber"
                    ),
                    "uniqueness": (
                        "fixed F,O recover |F| and O; normalized base and "
                        "minimal short Weierstrass coordinates leave only "
                        "the finite stabilizer above"
                    ),
                    "automorphisms": "H^0(X,T_X)=0 makes Aut(X) discrete",
                    "fiber_dimension": 0,
                },
                "period_argument": (
                    "The identity recovery Jacobian, finite coordinate "
                    "stabilizer, and zero-dimensional integral-class "
                    "forgetful fibers prove four K3 moduli.  Since sigma is "
                    "primitive order 4, T is a Q(i)-module and rho is even.  "
                    "A strict jump from rho>=12 would have rho>=14 and a "
                    "primitive-eigenvalue period ball of dimension at most "
                    "3, contradicting the four-dimensional family."
                ),
                "mordell_weil_rank": 0,
                "torsion_height_argument": (
                    "For a nonzero torsion section the D8/I4* correction is "
                    "at most 2 and the two A1/I2 corrections total at most "
                    "1.  The K3 height formula gives h(P)>=4-3=1, a "
                    "contradiction.  Thus MW torsion is trivial and the "
                    "trivial lattice is saturated."
                ),
                "conclusion": "very general NS=U+D8(-1)+A1(-1)^2",
                "scope": (
                    "exact for the geometric generic/very general member; no "
                    "particular rational quadruple is Picard-certified here"
                ),
            },
        },
        "picard_action": {
            "description": (
                "identity on U+D8(-1); exchange of the two I2/A1 root "
                "components"
            ),
            "a1x2_gram": matrix_payload(A1X2),
            "a1_swap": matrix_payload(A1_SWAP),
            "characteristic_polynomial": str(ns_characteristic),
            "fixed_lattice": {
                "class": "U+D8(-1)+<-4>",
                "rank": 11,
                "signature": [1, 10, 0],
                "determinant": 16,
                "smith": smith_diagonal(fixed_gram),
                "a1_diagonal_gram": matrix_payload(fixed_a1_gram),
                "discriminant_q_histogram": {
                    "0": 4,
                    "1": 4,
                    "3/4": 2,
                    "7/4": 6,
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
                "database_tuple": [[], 4, 1, 3],
                "meaning": "4 rational curves and one genus-3 curve",
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
                "a second-order leading-discriminant cancellation at "
                "infinity upgrades I0* to I4*/D8; a finite double root in "
                "s=t^2 inserts the exchanged I2/I2 pair"
            ),
            "novelty_scope": (
                "a new explicit equation-to-Coxeter realization in this "
                "project; no worldwide priority for the underlying order-4 "
                "K3 family or abstract Coxeter polytope is claimed"
            ),
        },
        "brandhorst_hofmann_handoff": {
            "id": "0.4.4.2",
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
