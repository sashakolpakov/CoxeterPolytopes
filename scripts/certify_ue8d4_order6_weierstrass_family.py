#!/usr/bin/env python3
"""Certify an explicit Weierstrass family for BH 0.6.3.1 / D4 triality.

For pairwise distinct nonzero ``1,a,b,c`` consider the minimal Weierstrass K3

    y^2 = x^3 + t^3 (t-1)(t-a)(t-b)(t-c)

and ``zeta(x,y,t) = (omega*x,-y,t)``, where ``omega^3=1`` and ``omega != 1``.
The fibers are I0* at zero, four II fibers, and II* at infinity.  The I0*
components give D4 and zeta cyclically permutes its three nonidentity arms.

This script checks the exact lattice/action, moduli-dimension, fixed-locus,
Brandhorst--Hofmann, and Vinberg handoffs used in the accompanying report.
The algebro-geometric inputs (Kodaira's valuation table and the standard
period-dimension argument) are recorded explicitly rather than disguised as
computer calculations.
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
SOURCE = (
    ROOT
    / "data/enumerations/roulleau_ue8d4_generic_ns_k3_automorphisms.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/ue8d4_order6_weierstrass_family.json"
)


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
D4 = negative_cartan(4, [(0, 1), (1, 2), (1, 3)])

# Row-action convention, identical to the archived R14_UE8D4_C3_S3 matrix.
D4_TRIALITY = sp.Matrix(
    [
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
    ]
)


def block_diag(*blocks: sp.Matrix) -> sp.Matrix:
    return sp.diag(*blocks)


def matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload_digest(payload: dict[str, Any]) -> str:
    copy = dict(payload)
    copy.pop("payload_sha256", None)
    encoded = json.dumps(copy, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def source_case() -> tuple[dict[str, Any], dict[str, Any]]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    case = next(
        row
        for row in source["cases"]
        if row["brandhorst_hofmann_row"]["id"] == "0.6.3.1"
    )
    return source, case


def build_certificate() -> dict[str, Any]:
    source, case = source_case()
    ns_gram = block_diag(U, E8, D4)
    ns_action = block_diag(sp.eye(10), D4_TRIALITY)
    fixed_d4_basis = sp.Matrix([[0, 1, 0, 0], [1, 0, 1, 1]])
    fixed_d4_gram = fixed_d4_basis * D4 * fixed_d4_basis.T
    coinvariant_basis = sp.Matrix([[1, 0, -1, 0], [0, 0, 1, -1]])
    coinvariant_gram = coinvariant_basis * D4 * coinvariant_basis.T

    fiber_rows = [
        {
            "place": "t=0",
            "ord_a6": 3,
            "ord_discriminant": 6,
            "kodaira": "I0*",
            "euler": 6,
            "root_lattice": "D4(-1)",
            "root_rank": 4,
        },
        {
            "place": "four roots of q4",
            "multiplicity": 4,
            "ord_a6": 1,
            "ord_discriminant": 2,
            "kodaira": "II",
            "euler_each": 2,
            "root_rank": 0,
        },
        {
            "place": "t=infinity",
            "ord_a6": 5,
            "ord_discriminant": 10,
            "kodaira": "II*",
            "euler": 10,
            "root_lattice": "E8(-1)",
            "root_rank": 8,
        },
    ]

    bh_match = case["exact_picard_action_match"]
    hyperbolic = case["hyperbolic_polytope"]
    expected_action = sp.Matrix(
        bh_match["equivariant_orthogonal_decomposition"]["d4_action"]
    )

    phi3 = sp.Symbol("lambda") ** 2 + sp.Symbol("lambda") + 1
    phi6 = sp.Symbol("lambda") ** 2 - sp.Symbol("lambda") + 1
    ns_char = sp.factor(ns_action.charpoly().as_expr())
    full_char = sp.factor(ns_char * phi6**4)
    traces = {
        str(power): int(sp.trace(ns_action**power))
        + 4
        * (
            1
            if power % 6 in (1, 5)
            else -1
            if power % 6 in (2, 4)
            else -2
            if power % 6 == 3
            else 2
        )
        for power in (1, 2, 3)
    }

    checks = {
        "fiber_euler_sum_is_24": 6 + 4 * 2 + 10 == 24,
        "trivial_lattice_rank_is_14": 2 + 4 + 8 == 14,
        "ns_signature": exact_inertia(ns_gram) == (1, 13, 0),
        "ns_determinant": ns_gram.det() == -4,
        "ns_smith": smith_diagonal(ns_gram) == [*[1] * 12, 2, 2],
        "triality_isometry": D4_TRIALITY * D4 * D4_TRIALITY.T == D4,
        "triality_order": D4_TRIALITY**3 == sp.eye(4),
        "triality_matches_bh": D4_TRIALITY == expected_action,
        "ns_characteristic": ns_char
        == (sp.Symbol("lambda") - 1) ** 12 * phi3,
        "fixed_d4_determinant": fixed_d4_gram.det() == 3,
        "fixed_ns_rank_determinant": (
            fixed_d4_gram.rows + 10 == 12 and -fixed_d4_gram.det() == -3
        ),
        "coinvariant_is_A2_minus2": coinvariant_gram
        == sp.Matrix([[-4, 2], [2, -4]]),
        "full_characteristic": full_char
        == (sp.Symbol("lambda") - 1) ** 12 * phi3 * phi6**4,
        "fixed_euler_characteristics": traces == {"1": 15, "2": 7, "3": 6},
        "source_all_checks_pass": source["all_checks_pass"] is True,
        "bh_generic_ns": case["brandhorst_hofmann_row"]["generic_picard_lattice"][
            "isometry_class"
        ]
        == "U+E8(-1)+D4(-1)",
        "bh_picard_action_order": bh_match["picard_action_order"] == 3,
        "vinberg_finite_volume": hyperbolic["classification"]
        == "finite_volume_noncompact",
    }
    if not all(checks.values()):
        raise AssertionError(f"a Weierstrass-family check failed: {checks}")

    payload: dict[str, Any] = {
        "schema": "ue8d4-order6-weierstrass-family-v1",
        "theorem": (
            "The very general minimal resolution of "
            "y^2=x^3+t^3(t-1)(t-a)(t-b)(t-c), with pairwise distinct "
            "nonzero 1,a,b,c, has NS=U+E8(-1)+D4(-1).  The order-six "
            "automorphism (x,y,t)->(omega*x,-y,t) induces D4 triality and "
            "therefore explicitly realizes BH 0.6.3.1 and the archived "
            "finite-volume noncompact H11/13 Coxeter section."
        ),
        "family": {
            "affine_equation": "y^2 = x^3 + t^3*(t-1)*(t-a)*(t-b)*(t-c)",
            "homogeneous_a6": (
                "U^3*(U-V)*(U-aV)*(U-bV)*(U-cV)*V^5"
            ),
            "parameter_open": (
                "a,b,c are nonzero, pairwise distinct, and different from 1"
            ),
            "automorphism": "zeta(x,y,t)=(omega*x,-y,t), omega^3=1, omega!=1",
            "global_order": 6,
            "holomorphic_two_form": "dt wedge dx/y",
            "two_form_multiplier": "-omega (a primitive sixth root)",
            "base_action": "identity",
            "discriminant": "-432*[t^3*(t-1)*(t-a)*(t-b)*(t-c)]^2",
            "fibers": fiber_rows,
        },
        "dillies_lineage": {
            "primary_source": "https://arxiv.org/abs/0912.5228",
            "classification_row": 11,
            "pre_collision_match_scope": (
                "Dillies table 1 row 11 gives the fixed-locus data for zeta "
                "and zeta^2; it does not determine the zeta^3 fixed locus"
            ),
            "row_11_generic_root_multiplicities": [5, 2, 1, 1, 1, 1, 1],
            "our_root_multiplicities": [5, 3, 1, 1, 1, 1],
            "noether_lefschetz_collision": "II + IV -> I0*",
            "dimension_before_collision": 4,
            "dimension_after_collision": 3,
            "zeta3_collision_change": {
                "fixed_three_section": "normalization of x^3=-p12(t)",
                "before_collision_full_branch_points": 7,
                "before_collision_genus": 5,
                "after_collision_full_branch_points": 5,
                "after_collision_genus": 3,
                "riemann_hurwitz": "2g-2=-6+2r for r full branch points",
            },
            "fixed_locus_generator": {
                "bh_tuple": [[5, 8], 2, 0, 1],
                "meaning": "5 type-(2,5) points, 8 type-(3,4) points, 2 rational curves",
                "status": (
                    "detailed BH 0.6.3.1 fingerprint, compatible with the "
                    "Dillies row-11 local data; this script directly checks "
                    "its Euler characteristic"
                ),
            },
            "fixed_locus_square": {
                "bh_tuple": [[5], 3, 1, 2],
                "meaning": "5 points, 3 rational curves, one genus-2 curve",
                "status": (
                    "detailed BH 0.6.3.1 fingerprint, compatible with the "
                    "Dillies row-11 order-three fixed-locus data; this script "
                    "directly checks its Euler characteristic"
                ),
            },
            "fixed_locus_cube": {
                "bh_tuple": [[], 6, 1, 3],
                "meaning": "6 rational curves and one genus-3 curve",
                "status": (
                    "detailed BH 0.6.3.1 fingerprint; the genus-three change "
                    "is independently obtained from the cyclic three-section "
                    "by Riemann--Hurwitz, while this script directly checks "
                    "the total Euler characteristic"
                ),
            },
        },
        "neron_severi": {
            "fiber_and_zero_section": "U",
            "reducible_fiber_roots": "E8(-1)+D4(-1)",
            "lower_bound_lattice": "U+E8(-1)+D4(-1)",
            "gram": matrix_payload(ns_gram),
            "rank_signature_determinant": [14, [1, 13, 0], -4],
            "smith": smith_diagonal(ns_gram),
            "primitive_at_rank_14": (
                "an index-two even overlattice would be even unimodular of "
                "signature (1,13), impossible since 1-13 is not divisible by 8"
            ),
            "generic_exactness": {
                "parameter_moduli_dimension": 3,
                "argument": (
                    "I0* and II* distinguish 0 and infinity; the four II "
                    "positions modulo t-scaling have dimension 3.  Forgetting "
                    "the marked fibration has discrete fibers.  Since zeta is "
                    "primitive order 6, T is a Q(zeta_6)-module, so rank(T) "
                    "and rho=22-rank(T) are even.  Since rho>=14, any strict "
                    "jump would therefore have rho>=16.  Then rank(T)<=6 and "
                    "the primitive-eigenvalue period domain would have "
                    "dimension at most 2, contradicting the 3-dimensional "
                    "family.  Hence very general rho=14."
                ),
                "scope": (
                    "exact for the geometric generic/very general member; no "
                    "specific rational triple (a,b,c) is certified here"
                ),
            },
        },
        "picard_action": {
            "description": (
                "identity on U+E8(-1); the I0* Tate cubic has three roots "
                "cycled by x->omega*x, giving the outer D4 3-cycle"
            ),
            "d4_gram": matrix_payload(D4),
            "d4_triality": matrix_payload(D4_TRIALITY),
            "characteristic_polynomial": str(ns_char),
            "fixed_lattice": {
                "class": "U+E8(-1)+A2(-1)",
                "rank": 12,
                "determinant": -3,
                "d4_fixed_gram": matrix_payload(fixed_d4_gram),
            },
            "coinvariant": {
                "class": "A2(-2)",
                "gram": matrix_payload(coinvariant_gram),
                "rootless": True,
            },
        },
        "cohomology_and_fixed_loci": {
            "generic_transcendental_rank": 8,
            "transcendental_characteristic": "(lambda^2-lambda+1)^4",
            "full_characteristic": str(full_char),
            "h2_traces_powers_1_2_3": traces,
            "fixed_euler_powers_1_2_3": {
                power: trace + 2 for power, trace in traces.items()
            },
        },
        "brandhorst_hofmann_handoff": {
            "id": "0.6.3.1",
            "article": "https://doi.org/10.1017/fms.2023.50",
            "database": "https://doi.org/10.5281/zenodo.10410278",
            "existing_certificate_sha256": sha256_file(SOURCE),
            "archived_type_id": bh_match["archived_type_id"],
            "same_d4_action_entrywise": True,
            "same_discriminant_action": bh_match["discriminant_form_and_action"][
                "same_q_action_class"
            ],
        },
        "hyperbolic_handoff": {
            "dimension": hyperbolic["dimension"],
            "facets": hyperbolic["facets"],
            "classification": hyperbolic["classification"],
            "finite_vertices": hyperbolic["finite_vertices"],
            "cusps": hyperbolic["cusps"],
            "elliptic_edges": hyperbolic["elliptic_edges"],
            "bad_edges": hyperbolic["bad_edges"],
            "abstract_polytope_novelty": "known; this is a new explicit realization in the project",
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
    print(json.dumps({"all_checks_pass": True, "payload_sha256": payload["payload_sha256"]}))


if __name__ == "__main__":
    main()
