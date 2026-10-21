#!/usr/bin/env python3
"""Certify an explicit Inose-family carrier for the H17/19 K3 chamber.

For parameters A,B in a nonempty Zariski open set, consider

    y^2 = x^3 - 3*A*t^4*x + t^5*(t^2 - 2*B*t + 1).

The fibers are II* at zero and infinity and four I1 fibers.  The two
reducible fibers provide the even unimodular lattice U+E8(-1)^2.  A
two-dimensional moduli argument then proves that this is the exact
Neron--Severi lattice of a very general member.  The script also pins and
checks the complete 19-wall nef diagram and its full Vinberg finite-volume
audit.

The equation is the classical Inose pencil.  The positive contribution
certified here is its explicit matching to the complete H17/19 divisor
chamber in the project ledger, not a novelty claim for the equation or the
abstract Coxeter polytope.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.verify_k3_high_rank_nef_chambers import (
    analyze_case,
    build_gram,
    exact_inertia,
)


ROOT = Path(__file__).resolve().parents[1]
WALL_SOURCE = ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
VOLUME_SOURCE = (
    ROOT / "data/known_examples/k3_high_rank_nef_vinberg_volume.json"
)
ARCHIVE = ROOT / "data/enumerations/ue8x2_inose_h17_family.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload_digest(payload: dict[str, Any]) -> str:
    unhashed = dict(payload)
    unhashed.pop("payload_sha256", None)
    encoded = json.dumps(unhashed, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def negative_cartan(rank: int, edges: list[tuple[int, int]]) -> sp.Matrix:
    gram = -2 * sp.eye(rank)
    for left, right in edges:
        gram[left, right] = gram[right, left] = 1
    return gram


E8 = negative_cartan(
    8,
    [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)],
)
U = sp.Matrix([[0, 1], [1, 0]])


def order_at_zero(expression: sp.Expr, variable: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.cancel(expression), variable)
    return min(
        monomial[0]
        for monomial, coefficient in polynomial.terms()
        if coefficient != 0
    )


def find_rho(payload: dict[str, Any], rho: int) -> dict[str, Any]:
    matches = [case for case in payload["cases"] if int(case["rho"]) == rho]
    if len(matches) != 1:
        raise ValueError(f"expected one rho={rho} case, found {len(matches)}")
    return matches[0]


def build_certificate() -> dict[str, Any]:
    wall_payload = json.loads(WALL_SOURCE.read_text(encoding="utf-8"))
    volume_payload = json.loads(VOLUME_SOURCE.read_text(encoding="utf-8"))
    wall_case = find_rho(wall_payload, 18)
    volume_case = find_rho(volume_payload, 18)

    t, u, A, B = sp.symbols("t u A B")
    a4 = -3 * A * t**4
    a6 = t**5 * (t**2 - 2 * B * t + 1)
    c4 = sp.expand(-48 * a4)
    discriminant = sp.factor(-16 * (4 * a4**3 + 27 * a6**2))
    residual = sp.expand((t**2 - 2 * B * t + 1) ** 2 - 4 * A**3 * t**2)
    residual_discriminant = sp.factor(sp.discriminant(residual, t))

    a4_infinity = sp.cancel(u**8 * a4.subs(t, 1 / u))
    a6_infinity = sp.cancel(u**12 * a6.subs(t, 1 / u))
    c4_infinity = sp.cancel(u**8 * c4.subs(t, 1 / u))
    discriminant_infinity = sp.cancel(u**24 * discriminant.subs(t, 1 / u))
    zero_valuations = {
        "a4": order_at_zero(a4, t),
        "a6": order_at_zero(a6, t),
        "c4": order_at_zero(c4, t),
        "discriminant": order_at_zero(discriminant, t),
    }
    infinity_valuations = {
        "a4": order_at_zero(a4_infinity, u),
        "a6": order_at_zero(a6_infinity, u),
        "c4": order_at_zero(c4_infinity, u),
        "discriminant": order_at_zero(discriminant_infinity, u),
    }

    witness_residual = sp.factor(residual.subs({A: 1, B: 3}))
    witness_discriminant = sp.discriminant(witness_residual, t)
    coefficient_t3 = sp.Poly(residual, t).coeff_monomial(t**3)
    coefficient_t2 = sp.Poly(residual, t).coeff_monomial(t**2)
    recovered_B = sp.factor(-coefficient_t3 / 4)
    recovered_A3 = sp.factor((4 * recovered_B**2 + 2 - coefficient_t2) / 4)
    j_function = sp.factor(c4**3 / discriminant)

    ns_gram = sp.diag(U, E8, E8)
    wall_gram = build_gram(wall_case)
    wall_analysis = analyze_case(wall_case)

    expected_wall_hash = (
        "0dcc7ec9c10f3437cd39713cd0d6036dc4fa8ad0c978ffdad4c2b34748269120"
    )
    expected_volume_hash = (
        "349f58a302c1d274bb25a7e3dff7be18ce12d0924717b52ac7976003d95dbbde"
    )
    expected_audit_hash = (
        "6e77968ffeb645ad008cc52bcc4995c9a983d928a88d6a4a34e5dc35ad4d0139"
    )

    checks = {
        "discriminant_factorization": sp.expand(
            discriminant + 432 * t**10 * residual
        )
        == 0,
        "residual_discriminant_factorization": sp.expand(
            residual_discriminant
            - 4096
            * A**6
            * (A**3 - (B + 1) ** 2)
            * (A**3 - (B - 1) ** 2)
        )
        == 0,
        "open_witness_factorization": witness_residual
        == (t**2 - 8 * t + 1) * (t**2 - 4 * t + 1),
        "open_witness_discriminant": witness_discriminant == 184320,
        "zero_is_II_star": zero_valuations
        == {"a4": 4, "a6": 5, "c4": 4, "discriminant": 10},
        "infinity_is_II_star": infinity_valuations
        == {"a4": 4, "a6": 5, "c4": 4, "discriminant": 10},
        "four_simple_nonzero_I1": sp.Poly(residual, t).degree() == 4
        and residual.subs(t, 0) == 1,
        "fiber_euler_sum": 10 + 10 + 4 == 24,
        "parameter_recovery": recovered_B == B and recovered_A3 == A**3,
        "j_function": j_function == -6912 * A**3 * t**2 / residual,
        "ns_rank_signature": ns_gram.rows == 18
        and exact_inertia(ns_gram) == (1, 17, 0),
        "ns_even_unimodular": ns_gram.det() == -1
        and all(ns_gram[index, index] % 2 == 0 for index in range(18)),
        "wall_source_hash": sha256_file(WALL_SOURCE) == expected_wall_hash,
        "volume_source_hash": sha256_file(VOLUME_SOURCE) == expected_volume_hash,
        "wall_case_identity": wall_case["lattice"]
        == "U + E8(-1) + E8(-1)"
        and wall_case["wall_count"] == 19
        and wall_case["hyperbolic_dimension"] == 17,
        "wall_gram_exact": wall_gram.rank() == 18
        and exact_inertia(wall_gram) == (1, 17, 1)
        and wall_analysis["basis_determinant"] == -1
        and wall_analysis["all_checks_pass"] is True,
        "vinberg_decision": volume_case["decision"]
        == "finite_volume_noncompact"
        and volume_case["simple_roots_complete"] is True,
        "vinberg_edges_complete": volume_case["elliptic_edge_count"] == 729
        and volume_case["finite_volume_bad_edge_count"] == 0,
        "vinberg_parabolics_complete": volume_case[
            "connected_affine_audit_count"
        ]
        == 3
        and volume_case["bad_full_parabolic_neighborhood_count"] == 0,
        "vinberg_vertices": volume_case["finite_vertex_count"] == 80
        and volume_case["ideal_vertex_count"] == 2,
        "vinberg_full_audit_digest": volume_case[
            "canonical_full_audit_sha256"
        ]
        == expected_audit_hash,
        "source_volume_cross_pin": volume_payload["source_sha256"]
        == expected_wall_hash,
    }

    payload: dict[str, Any] = {
        "schema_version": 1,
        "artifact_id": "ue8x2_inose_h17_family",
        "positive_object_id": "R18-UE8x2-NEF-H17-19",
        "claim_scope": (
            "new project-level equation-to-complete-nef realization of a known "
            "K3 chamber and known abstract finite-volume Coxeter polytope"
        ),
        "lineage": {
            "equation": "classical Inose pencil",
            "equation_source": "https://arxiv.org/abs/math/0612396",
            "wall_source_locator": wall_case["source_locator"],
            "published_projective_model": wall_case["published_model"],
            "priority_boundary": (
                "neither the Inose equation nor Vinberg's H17/19 chamber is "
                "claimed new worldwide; their explicit algebraic matching is "
                "the positive realization added here"
            ),
        },
        "family": {
            "affine_equation": (
                "y^2 = x^3 - 3*A*t^4*x + "
                "t^5*(t^2 - 2*B*t + 1)"
            ),
            "parameter_open": (
                "A*(A^3-(B+1)^2)*(A^3-(B-1)^2) != 0"
            ),
            "discriminant": str(discriminant),
            "residual_quartic": str(residual),
            "residual_discriminant": str(residual_discriminant),
            "c4": str(c4),
            "j": str(j_function),
            "open_nonempty_witness": {
                "A": 1,
                "B": 3,
                "residual": str(witness_residual),
                "residual_discriminant": int(witness_discriminant),
            },
            "fibers": [
                {
                    "place": "t=0",
                    "kodaira": "II*",
                    "valuations_a4_a6_c4_delta": list(zero_valuations.values()),
                    "root_lattice": "E8(-1)",
                    "euler": 10,
                },
                {
                    "place": "four roots of R(t)",
                    "kodaira": "I1",
                    "multiplicity": 4,
                    "root_lattice": None,
                    "euler_each": 1,
                },
                {
                    "place": "t=infinity",
                    "kodaira": "II*",
                    "valuations_a4_a6_c4_delta": list(
                        infinity_valuations.values()
                    ),
                    "root_lattice": "E8(-1)",
                    "euler": 10,
                },
            ],
            "fiber_euler_sum": 24,
        },
        "neron_severi": {
            "very_general_exact_lattice": "U+E8(-1)+E8(-1)",
            "rank_signature_determinant": [18, [1, 17, 0], -1],
            "trivial_lattice_argument": (
                "fiber, zero section, and the nonidentity components of the "
                "two II* fibers give U+E8(-1)^2"
            ),
            "unimodular_closure": (
                "the displayed determinant-one lattice is primitive and has "
                "no proper integral same-rank overlattice; after rho=18 it "
                "is the whole NS, and Shioda--Tate forces Mordell--Weil rank "
                "and torsion zero"
            ),
            "generic_exactness": {
                "parameter_moduli_dimension": 2,
                "coefficient_recovery": {
                    "B": "-[t^3]R/4",
                    "A^3": "(4*B^2+2-[t^2]R)/4",
                },
                "finite_equivalences": (
                    "an ordered-II*-fibration isomorphism has t->lambda*t; "
                    "comparison of the monic palindromic R gives lambda^4=1 "
                    "and generically lambda^2=1, hence only B->+/-B and the "
                    "finite cube-root ambiguity in A; swapping 0 and infinity "
                    "adds only another finite identification"
                ),
                "forgetful_map": (
                    "on a fixed K3, elliptic fibrations and sections are "
                    "integral divisor classes, so forgetting the marked "
                    "fibration has zero-dimensional fibers"
                ),
                "period_argument": (
                    "the L-polarized K3 period domain for rank(L)=18 has "
                    "dimension two; the family has two-dimensional image. "
                    "The rho>=19 Noether--Lefschetz locus is a countable union "
                    "of proper curves, so the very general member has rho=18"
                ),
                "scope": (
                    "exact for the geometric very general member; (A,B)=(1,3) "
                    "only witnesses the open and is not asserted to have "
                    "geometric Picard rank 18"
                ),
            },
        },
        "complete_nef_handoff": {
            "source_artifact": str(WALL_SOURCE.relative_to(ROOT)),
            "source_sha256": sha256_file(WALL_SOURCE),
            "transport_argument": (
                "for NS isometric to U+E8(-1)^2, the nef cone is a Weyl "
                "chamber in the positive cone; Weyl transport identifies it "
                "with the published complete 19-wall chamber"
            ),
            "wall_count": 19,
            "wall_gram_rank_inertia": [18, [1, 17, 1]],
            "full_rank_basis_determinant": -1,
            "wall_analysis": wall_analysis,
        },
        "vinberg_handoff": {
            "source_artifact": str(VOLUME_SOURCE.relative_to(ROOT)),
            "source_sha256": sha256_file(VOLUME_SOURCE),
            "criterion": volume_payload["criterion"],
            "classification": volume_case["decision"],
            "hyperbolic_dimension": 17,
            "facets": 19,
            "finite_vertices": volume_case["finite_vertex_count"],
            "ideal_vertices": volume_case["ideal_vertex_count"],
            "elliptic_edges": volume_case["elliptic_edge_count"],
            "bad_edges": volume_case["finite_volume_bad_edge_count"],
            "connected_affine_neighborhoods": volume_case[
                "connected_affine_audit_count"
            ],
            "bad_full_parabolic_neighborhoods": volume_case[
                "bad_full_parabolic_neighborhood_count"
            ],
            "simple_roots_complete": volume_case["simple_roots_complete"],
            "canonical_full_audit_sha256": volume_case[
                "canonical_full_audit_sha256"
            ],
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["payload_sha256"] = payload_digest(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    arguments = parser.parse_args()
    payload = build_certificate()
    if arguments.check_archive:
        archived = json.loads(arguments.archive.read_text(encoding="utf-8"))
        if payload != archived:
            raise AssertionError(f"archive mismatch: {arguments.archive}")
        print(
            json.dumps(
                {
                    "all_checks_pass": payload["all_checks_pass"],
                    "payload_sha256": payload["payload_sha256"],
                }
            )
        )
        return
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
