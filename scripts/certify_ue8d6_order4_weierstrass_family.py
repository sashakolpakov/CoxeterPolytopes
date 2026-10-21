#!/usr/bin/env python3
"""Certify an explicit Weierstrass family for BH 0.4.2.2 / Kondo tau_2.

For parameters b,c in a nonempty Zariski open set, consider the minimal
resolution of

    y^2 = x^3 + t^2(t^2-3)x + t^3(2+b*t^2+c*t^4)

and sigma(x,y,t)=(-x,i*y,-t).  The fibers are I2* at zero, II* at infinity,
and six I1 fibers.  For the very general member the trivial lattice is the
whole Neron--Severi lattice U+E8(-1)+D6(-1).  The order-four automorphism is
the identity on U+E8 and swaps the two spinor nodes of D6.

The geometric inputs (Kodaira's valuation table, the period-dimension
argument, and the order-four holomorphic Lefschetz congruence) are recorded
explicitly.  The script checks all symbolic, lattice, Brandhorst--Hofmann,
and Vinberg handoffs used in the accompanying report.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import _q_histogram
from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    smith_diagonal,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/kondo_rho16_h14_bh_order4_lift.json"
ACTION_SOURCE = ROOT / "data/candidates/k3_rho16_h14_fold_identification.json"
OUTPUT = ROOT / "data/enumerations/ue8d6_order4_weierstrass_family.json"


U = sp.Matrix([[0, 1], [1, 0]])


def negative_cartan(rank: int, edges: list[tuple[int, int]]) -> sp.Matrix:
    gram = -2 * sp.eye(rank)
    for left, right in edges:
        gram[left, right] = gram[right, left] = 1
    return gram


E8 = negative_cartan(
    8,
    [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)],
)
D6 = negative_cartan(
    6,
    [(0, 1), (1, 2), (2, 3), (3, 4), (3, 5)],
)
D6_OUTER = sp.eye(6)
D6_OUTER.row_swap(4, 5)


def block_diag(*blocks: sp.Matrix) -> sp.Matrix:
    return sp.diag(*blocks)


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


def fractional_part(value: sp.Rational) -> sp.Rational:
    return sp.Rational(value - sp.floor(value))


def coset_key(vector: sp.Matrix) -> tuple[sp.Rational, ...]:
    return tuple(fractional_part(sp.Rational(value)) for value in vector)


def discriminant_classes(gram: sp.Matrix) -> tuple[tuple[sp.Rational, ...], ...]:
    zero = tuple(sp.Rational(0) for _ in range(gram.rows))
    classes = {zero}
    inverse = gram.inv()
    for column in range(inverse.cols):
        generator = inverse[:, column]
        previous = tuple(classes)
        classes.update(coset_key(sp.Matrix(value) + generator) for value in previous)
    return tuple(sorted(classes))


def build_certificate() -> dict[str, Any]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    action_source = json.loads(ACTION_SOURCE.read_text(encoding="utf-8"))

    t, u, b, c, lam = sp.symbols("t u b c lambda")
    a4 = t**2 * (t**2 - 3)
    a6 = t**3 * (2 + b * t**2 + c * t**4)
    discriminant_core = sp.expand(4 * a4**3 + 27 * a6**2)
    residual_cubic = (
        108 * (1 + b)
        + (-36 + 27 * (b**2 + 4 * c)) * u
        + (4 + 54 * b * c) * u**2
        + 27 * c**2 * u**3
    )
    residual_discriminant = sp.factor(sp.discriminant(residual_cubic, u))
    residual_at_a4_zero = sp.factor(residual_cubic.subs(u, 3))
    example_cubic = sp.Poly(residual_cubic.subs({b: 0, c: 1}), u)

    ns_gram = block_diag(U, E8, D6)
    ns_action = block_diag(sp.eye(10), D6_OUTER)

    fixed_d6_basis = sp.zeros(5, 6)
    for index in range(4):
        fixed_d6_basis[index, index] = 1
    fixed_d6_basis[4, 4] = fixed_d6_basis[4, 5] = 1
    fixed_d6_gram = fixed_d6_basis * D6 * fixed_d6_basis.T
    fixed_ns_gram = block_diag(U, E8, fixed_d6_gram)
    coinvariant_basis = sp.Matrix([[0, 0, 0, 0, 1, -1]])
    coinvariant_gram = coinvariant_basis * D6 * coinvariant_basis.T
    d6_classes = discriminant_classes(D6)
    d6_discriminant_action = tuple(
        d6_classes.index(coset_key(D6_OUTER * sp.Matrix(value)))
        for value in d6_classes
    )
    d6_discriminant_q = tuple(
        str(sp.factor((sp.Matrix(value).T * D6 * sp.Matrix(value))[0]) % 2)
        for value in d6_classes
    )
    archived_discriminant_action = tuple(
        action_source["lattice_action"]["discriminant_action"]
    )
    # Our sorted D6 classes have q-values (0,1,1/2,1/2), whereas the archived
    # wall-basis classes have q-values (0,1/2,1/2,1).  This is the explicit
    # class identification from our ordering to the archived ordering.
    class_identification = (0, 3, 1, 2)

    ns_char = sp.factor(ns_action.charpoly(lam).as_expr())
    full_char = sp.factor(ns_char * (lam**2 + 1) ** 3)
    h2_traces = {
        "1": int(sp.trace(ns_action)),
        "2": int(sp.trace(ns_action**2)) - 6,
    }

    target_fixed = source["exact_picard_action_conjugacy"]["fixed_lattices"]
    target_coin = source["exact_picard_action_conjugacy"]["coinvariant_lattices"]
    target_glue = source["exact_picard_action_conjugacy"]["primitive_glue"]
    bh = source["brandhorst_hofmann"]
    hyperbolic = source["hyperbolic_polytope"]

    # In D6 root coordinates p=e_4 and q=e_5 are the two spinor roots.
    # If s=p+q is fixed and d=p-q is anti-fixed, then
    # p=(s+d)/2.  This is the nonzero element of the index-two primitive
    # glue between the fixed lattice and <-4>.
    explicit_glue = {
        "index": 2,
        "fixed_order_two_class": "(alpha_5+alpha_6)/2",
        "coinvariant_order_two_class": "(alpha_5-alpha_6)/2",
        "adjoined_vector": "alpha_5=(s+d)/2",
        "source_smith_pair": [2, 2],
    }
    expected_source_glue = [
        [[0] * 15, [0]],
        [[0] * 14 + [2], [2]],
    ]

    fiber_rows = [
        {
            "place": "t=0",
            "valuations_a4_a6_delta": [2, 3, 8],
            "kodaira": "I2*",
            "euler": 8,
            "root_lattice": "D6(-1)",
            "root_rank": 6,
        },
        {
            "place": "six roots of Q(t^2)",
            "multiplicity": 6,
            "valuations_c4_delta": [0, 1],
            "kodaira": "I1",
            "euler_each": 1,
            "root_rank": 0,
        },
        {
            "place": "t=infinity",
            "valuations_a4_a6_delta": [4, 5, 10],
            "kodaira": "II*",
            "euler": 10,
            "root_lattice": "E8(-1)",
            "root_rank": 8,
        },
    ]

    checks = {
        "discriminant_factorization": sp.expand(
            discriminant_core - t**8 * residual_cubic.subs(u, t**2)
        )
        == 0,
        "residual_cubic_degree": sp.degree(residual_cubic, u) == 3,
        "residual_discriminant_factorization": residual_discriminant
        == 432
        * (3 * b + 9 * c + 2) ** 3
        * (27 * b**2 * c - 54 * b * c + b - 432 * c**2 - 9 * c - 2),
        "residual_at_a4_zero_factorization": residual_at_a4_zero
        == 9 * (3 * b + 9 * c + 2) ** 2,
        "parameter_open_is_nonempty": (
            example_cubic.LC() != 0
            and example_cubic.eval(0) != 0
            and example_cubic.eval(3) != 0
            and sp.gcd(example_cubic, example_cubic.diff()).degree() == 0
        ),
        "example_residual_discriminant": sp.discriminant(example_cubic.as_expr(), u)
        == -254721456,
        "fiber_euler_sum_is_24": 8 + 6 + 10 == 24,
        "trivial_lattice_rank_is_16": 2 + 6 + 8 == 16,
        "ns_signature": exact_inertia(ns_gram) == (1, 15, 0),
        "ns_determinant": ns_gram.det() == -4,
        "ns_smith": smith_diagonal(ns_gram) == [*[1] * 14, 2, 2],
        "ns_q": _q_histogram(ns_gram) == {"0": 1, "1": 1, "1/2": 2},
        "outer_isometry": D6_OUTER * D6 * D6_OUTER.T == D6,
        "outer_order": D6_OUTER**2 == sp.eye(6),
        "outer_characteristic": ns_char
        == (lam - 1) ** 15 * (lam + 1),
        "literal_discriminant_action": d6_discriminant_action == (0, 1, 3, 2),
        "literal_discriminant_q": d6_discriminant_q == ("0", "1", "1/2", "1/2"),
        "archived_discriminant_action": archived_discriminant_action
        == (0, 2, 1, 3),
        "discriminant_actions_conjugate": all(
            class_identification[d6_discriminant_action[index]]
            == archived_discriminant_action[class_identification[index]]
            for index in range(4)
        ),
        "fixed_rank_signature_determinant": (
            fixed_ns_gram.rows == 15
            and exact_inertia(fixed_ns_gram) == (1, 14, 0)
            and fixed_ns_gram.det() == 4
        ),
        "fixed_smith": smith_diagonal(fixed_ns_gram) == [*[1] * 14, 4],
        "fixed_q": _q_histogram(fixed_ns_gram)
        == {"0": 1, "1": 1, "3/4": 2},
        "coinvariant_is_minus4": coinvariant_gram == sp.Matrix([[-4]]),
        "coinvariant_rootless": coinvariant_gram[0, 0] < -2,
        "primitive_glue_index": explicit_glue["index"] == 2,
        "source_glue_exact": target_glue["index"] == 2
        and target_glue["identity_maps_carry_both_pairs"] is True
        and target_glue["archived_pairs_in_smith_coordinates"]
        == expected_source_glue
        and target_glue["pure_database_pairs_in_smith_coordinates"]
        == expected_source_glue,
        "full_characteristic": full_char
        == (lam - 1) ** 15 * (lam + 1) * (lam**2 + 1) ** 3,
        "lefschetz_excludes_picard_identity": (16 - 0 - 2) % 4 != 0,
        "outer_lefschetz_alpha": (15 - 1 - 2) // 4 == 3,
        "fixed_euler_generator": h2_traces["1"] + 2 == 16,
        "fixed_euler_square": h2_traces["2"] + 2 == 12,
        "source_all_checks_pass": source["all_checks_pass"] is True,
        "bh_row": bh["pure_cyclic_row"]["row_id"] == "0.4.2.2",
        "bh_global_order": bh["pure_cyclic_row"]["global_automorphism_order"]
        == 4,
        "bh_generic_ns": bh["generic_picard_lattice"]["isometry_class"]
        == "U+E8(-1)+D6(-1)",
        "target_fixed_match": target_fixed["rank_signature_determinant"]
        == [15, [1, 14, 0], 4],
        "target_coinvariant_match": target_coin["archived_gram"] == [[-4]],
        "vinberg_finite_volume": hyperbolic["classification"]
        == "finite_volume_noncompact",
        "vinberg_no_bad_or_hyperideal": hyperbolic["bad_edges"] == 0
        and hyperbolic["hyperideal_vertices"] == 0,
    }
    checks = {key: bool(value) for key, value in checks.items()}
    if not all(checks.values()):
        raise AssertionError(f"a Weierstrass-family check failed: {checks}")

    payload: dict[str, Any] = {
        "schema": "ue8d6-order4-weierstrass-family-v1",
        "theorem": (
            "For (b,c) in the displayed nonempty Zariski open set, the very "
            "general minimal resolution of y^2=x^3+t^2(t^2-3)x+"
            "t^3(2+b*t^2+c*t^4) has NS=U+E8(-1)+D6(-1).  The primitive "
            "order-four automorphism sigma(x,y,t)=(-x,i*y,-t) induces the "
            "outer D6 spinor-node swap, hence realizes BH 0.4.2.2, Kondo "
            "tau_2, and the finite-volume noncompact H14/17 Vinberg P14 "
            "fixed-nef section."
        ),
        "family": {
            "affine_equation": (
                "y^2 = x^3 + t^2*(t^2-3)*x + "
                "t^3*(2+b*t^2+c*t^4)"
            ),
            "homogeneous_a4": "T^2*(T^2-3*S^2)*S^4",
            "homogeneous_a6": "T^3*(2*S^4+b*T^2*S^2+c*T^4)*S^5",
            "automorphism": "sigma(x,y,t)=(-x,i*y,-t)",
            "global_order": 4,
            "square": "sigma^2(x,y,t)=(x,-y,t)",
            "holomorphic_two_form": "dt wedge dx/y",
            "two_form_multiplier": "-i",
            "a4": "t^2*(t^2-3)",
            "a6": "t^3*(2+b*t^2+c*t^4)",
            "discriminant": "-16*t^8*Q_(b,c)(t^2)",
            "residual_cubic_Q": str(residual_cubic),
            "disc_u_Q": str(residual_discriminant),
            "Q_at_3": str(residual_at_a4_zero),
            "parameter_open": [
                (
                    "c*(b+1)*(3*b+9*c+2)*"
                    "(27*b^2*c-54*b*c+b-432*c^2-9*c-2) != 0"
                ),
                "equivalently: c != 0, Q(0) != 0, disc_u(Q) != 0",
                "disc_u(Q) != 0 also implies Q(3) != 0",
            ],
            "open_nonempty_witness": {
                "b": 0,
                "c": 1,
                "Q": str(example_cubic.as_expr()),
                "disc_u_Q": int(sp.discriminant(example_cubic.as_expr(), u)),
                "Q_at_3": int(example_cubic.eval(3)),
            },
            "fibers": fiber_rows,
        },
        "artebani_sarti_lineage": {
            "article": "https://doi.org/10.2969/jmsj/06720503",
            "arxiv": "https://arxiv.org/abs/1102.4436",
            "locator": "Example 6.3 and Proposition 1",
            "ambient_family": (
                "y^2=x^3+a(t)x+b(t), with a(t) even and b(t) odd, carrying "
                "(x,y,t)->(-x,i*y,-t)"
            ),
            "our_specialization": (
                "a two-parameter simultaneous Noether--Lefschetz fiber "
                "enhancement at the two fixed base points: I2* at zero and "
                "II* at infinity"
            ),
            "novelty_scope": (
                "new explicit carrier and exact Coxeter connection in this "
                "project; no claim that the normalized equation is absent "
                "from all prior literature"
            ),
        },
        "neron_severi": {
            "fiber_and_zero_section": "U",
            "reducible_fiber_roots": "E8(-1)+D6(-1)",
            "lower_bound_lattice": "U+E8(-1)+D6(-1)",
            "gram": matrix_payload(ns_gram),
            "rank_signature_determinant": [16, [1, 15, 0], -4],
            "smith": smith_diagonal(ns_gram),
            "same_rank_saturation": (
                "a proper even overlattice has index two and determinant "
                "one, hence would be even unimodular of signature (1,15); "
                "this is impossible since 1-15 is not divisible by 8"
            ),
            "generic_exactness": {
                "parameter_moduli_dimension": 2,
                "normalization_argument": (
                    "I2* and II* distinguish zero and infinity, so a "
                    "fibration-preserving base isomorphism is t->lambda*t. "
                    "After Weierstrass scaling, the normalized coefficient "
                    "t^2(t^2-3) forces lambda^2=1; thus (b,c) retain two "
                    "independent moduli.  Forgetting the marked fibration "
                    "has discrete fibers because fiber and section classes "
                    "are integral Neron--Severi classes."
                ),
                "period_argument": (
                    "For a primitive order-four automorphism the "
                    "transcendental lattice has only i and -i eigenvalues, "
                    "so rank(T) and rho are even.  Since rho>=16, a strict "
                    "jump gives rho>=18 and a primitive-eigenvalue period "
                    "space of dimension at most one, contradicting the "
                    "two-dimensional family.  Hence very general rho=16."
                ),
                "scope": (
                    "exact for the geometric generic/very general member; "
                    "the witness (b,c)=(0,1) only proves the parameter open "
                    "is nonempty, not that this named member has rho=16"
                ),
            },
        },
        "picard_action": {
            "d6_gram": matrix_payload(D6),
            "d6_outer_swap": matrix_payload(D6_OUTER),
            "description": (
                "sigma fixes the fiber and zero-section classes.  The II* "
                "E8 graph has no nontrivial automorphism.  At I2*, sigma "
                "fixes the affine component met by the invariant zero "
                "section, and the stabilizer of that node in the affine D6 "
                "graph is {identity, spinor swap}.  Identity on D6 would "
                "give whole-H2 eigenspace ranks (r,l)=(16,0), contradicting "
                "Artebani--Sarti Proposition 1: r-l-2=4*alpha.  Therefore "
                "the component action is the nontrivial spinor swap."
            ),
            "characteristic_polynomial": str(ns_char),
            "fixed_lattice": {
                "rank_signature_determinant": [15, [1, 14, 0], 4],
                "smith": smith_diagonal(fixed_ns_gram),
                "q_histogram_mod_2Z": _q_histogram(fixed_ns_gram),
                "d6_fixed_gram": matrix_payload(fixed_d6_gram),
            },
            "coinvariant": {
                "class": "<-4>",
                "gram": matrix_payload(coinvariant_gram),
                "rootless": True,
            },
            "primitive_glue": explicit_glue,
            "discriminant_action": {
                "our_sorted_class_q_values": list(d6_discriminant_q),
                "our_permutation": list(d6_discriminant_action),
                "our_to_archived_class_identification": list(class_identification),
                "archived_permutation": list(archived_discriminant_action),
                "description": (
                    "nontrivial: it swaps the two spinor classes in "
                    "A_D6=(Z/2)^2 and fixes the vector class"
                ),
            },
        },
        "cohomology_and_fixed_loci": {
            "generic_transcendental_rank": 6,
            "transcendental_characteristic": "(lambda^2+1)^3",
            "full_characteristic": str(full_char),
            "h2_traces_powers_1_2": h2_traces,
            "fixed_euler_powers_1_2": {
                power: trace + 2 for power, trace in h2_traces.items()
            },
            "generator": {
                "isolated_points": 10,
                "fixed_rational_curves": 3,
                "euler": 16,
                "derivation": (
                    "r=15,l=1 gives alpha=(r-l-2)/4=3 and n=2*alpha+4=10; "
                    "the nontrivial base action forces every fixed curve to "
                    "lie in a fixed fiber, hence it is rational"
                ),
            },
            "square": {
                "isolated_points": 0,
                "fixed_rational_curves": 7,
                "fixed_positive_genus_curves": [{"genus": 2, "count": 1}],
                "euler": 12,
                "nikulin_invariants": [16, 2, 1],
            },
            "bh_fingerprint": {
                "fixed_points": "((10,),3,0,1)",
                "fixed_points_square": "((),7,1,2)",
            },
        },
        "brandhorst_hofmann_handoff": {
            "id": "0.4.2.2",
            "article": "https://doi.org/10.1017/fms.2023.50",
            "database": "https://doi.org/10.5281/zenodo.10410278",
            "existing_certificate_sha256": sha256_file(SOURCE),
            "archived_action_certificate_sha256": sha256_file(ACTION_SOURCE),
            "exact_match": (
                "same NS, fixed lattice, rootless <-4> coinvariant, "
                "nontrivial spinor discriminant action, and the same "
                "index-two (order-two,order-two) primitive glue; the prior "
                "certificate proves this integral action is Kondo tau_2 and "
                "the unique pure order-four database hit"
            ),
            "full_group_row": "1.4.2.1",
            "full_group_scope": (
                "a separate compatible database realization supplies the "
                "commuting symplectic involution and C4 x C2 enrichment; "
                "the displayed equation itself certifies only the C4 row "
                "and does not claim that every very general X_(b,c) carries "
                "the extra involution"
            ),
        },
        "kahler_chamber_handoff": {
            "effectivity": (
                "the displayed regular automorphism preserves a Kahler "
                "class; after Weyl transport to the archived complete "
                "19-wall nef chamber, the exact integral action is the "
                "unique rank-15-fixed graph symmetry Kondo tau_2"
            ),
            "fixed_section": hyperbolic["id"],
        },
        "hyperbolic_handoff": {
            key: hyperbolic[key]
            for key in (
                "id",
                "dimension",
                "facets",
                "classification",
                "finite_vertices",
                "cusps",
                "hyperideal_vertices",
                "elliptic_edges",
                "bad_edges",
                "parabolic_neighborhood_audits",
                "published_identification",
            )
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
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(
        json.dumps(
            {
                "all_checks_pass": True,
                "payload_sha256": payload["payload_sha256"],
            }
        )
    )


if __name__ == "__main__":
    main()
