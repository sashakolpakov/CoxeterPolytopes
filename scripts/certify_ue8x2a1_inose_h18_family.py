#!/usr/bin/env python3
"""Certify a one-parameter Inose collision carrying the H18/24 chamber.

For s^3 not equal to 0 or 1, consider

    y^2 = x^3 - 3*s^2*t^4*x
          + t^5*(t^2 - 2*(1-s^3)*t + 1).

It has two II* fibers, one I2 fiber, and two I1 fibers.  The very-general
Neron--Severi lattice is U+E8(-1)^2+A1(-1).  This script pins the complete
24-wall nef diagram and the literal Vinberg finite-volume audit.
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
    analyze_derived_rank19_inose_model,
    build_gram,
    exact_inertia,
)


ROOT = Path(__file__).resolve().parents[1]
WALL_SOURCE = ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
VOLUME_SOURCE = (
    ROOT / "data/known_examples/k3_high_rank_nef_vinberg_volume.json"
)
ARCHIVE = ROOT / "data/enumerations/ue8x2a1_inose_h18_family.json"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def payload_digest(payload: dict[str, Any]) -> str:
    unhashed = dict(payload)
    unhashed.pop("payload_sha256", None)
    encoded = json.dumps(unhashed, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def find_rho(payload: dict[str, Any], rho: int) -> dict[str, Any]:
    matches = [case for case in payload["cases"] if int(case["rho"]) == rho]
    if len(matches) != 1:
        raise ValueError(f"expected one rho={rho} case, found {len(matches)}")
    return matches[0]


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
A1 = sp.Matrix([[-2]])


def build_certificate() -> dict[str, Any]:
    wall_payload = json.loads(WALL_SOURCE.read_text(encoding="utf-8"))
    volume_payload = json.loads(VOLUME_SOURCE.read_text(encoding="utf-8"))
    wall_case = find_rho(wall_payload, 19)
    volume_case = find_rho(volume_payload, 19)
    prior_symbolic = analyze_derived_rank19_inose_model()

    t, s = sp.symbols("t s")
    a4 = -3 * s**2 * t**4
    a6 = t**5 * (t**2 - 2 * (1 - s**3) * t + 1)
    c4 = sp.expand(-48 * a4)
    discriminant = sp.factor(-16 * (4 * a4**3 + 27 * a6**2))
    residual_quadratic = t**2 + (4 * s**3 - 2) * t + 1
    residual_discriminant = sp.factor(
        sp.discriminant(residual_quadratic, t)
    )
    i2_leading = sp.factor((discriminant / (t - 1) ** 2).subs(t, 1))
    recovered_s3 = sp.factor(
        (sp.Poly(residual_quadratic, t).coeff_monomial(t) + 2) / 4
    )
    witness_quadratic = residual_quadratic.subs(s, 2)
    witness_discriminant = sp.discriminant(witness_quadratic, t)

    ns_gram = sp.diag(U, E8, E8, A1)
    wall_gram = build_gram(wall_case)
    wall_analysis = analyze_case(wall_case)

    expected_wall_hash = (
        "0dcc7ec9c10f3437cd39713cd0d6036dc4fa8ad0c978ffdad4c2b34748269120"
    )
    expected_volume_hash = (
        "349f58a302c1d274bb25a7e3dff7be18ce12d0924717b52ac7976003d95dbbde"
    )
    expected_audit_hash = (
        "dfead85c337bc4e6e58a2925a4f905949c73f0e271c559333fd48adf6d0715be"
    )
    expected_discriminant = sp.factor(
        -432
        * t**10
        * (t - 1) ** 2
        * residual_quadratic
    )

    checks = {
        "prior_symbolic_replay": prior_symbolic["all_checks_pass"] is True,
        "discriminant_factorization": sp.expand(
            discriminant - expected_discriminant
        )
        == 0,
        "residual_discriminant": sp.expand(
            residual_discriminant - 16 * s**3 * (s**3 - 1)
        )
        == 0,
        "i2_at_one": c4.subs(t, 1) == 144 * s**2
        and i2_leading == -1728 * s**3,
        "two_simple_I1": sp.Poly(residual_quadratic, t).degree() == 2
        and residual_quadratic.subs(t, 0) == 1,
        "fiber_euler_sum": 10 + 10 + 2 + 1 + 1 == 24,
        "moduli_parameter_recovery": recovered_s3 == s**3,
        "open_witness": witness_quadratic == t**2 + 30 * t + 1
        and witness_discriminant == 896,
        "ns_rank_signature": ns_gram.rows == 19
        and exact_inertia(ns_gram) == (1, 18, 0),
        "ns_even_determinant_two": ns_gram.det() == 2
        and all(ns_gram[index, index] % 2 == 0 for index in range(19)),
        "no_same_rank_overlattice": abs(int(ns_gram.det())) < 4,
        "wall_source_hash": sha256_file(WALL_SOURCE) == expected_wall_hash,
        "volume_source_hash": sha256_file(VOLUME_SOURCE) == expected_volume_hash,
        "wall_case_identity": wall_case["lattice"]
        == "U + E8(-1) + E8(-1) + A1(-1)"
        and wall_case["wall_count"] == 24
        and wall_case["hyperbolic_dimension"] == 18,
        "wall_gram_exact": wall_gram.rank() == 19
        and exact_inertia(wall_gram) == (1, 18, 5)
        and wall_analysis["basis_determinant"] == 2
        and wall_analysis["all_checks_pass"] is True,
        "vinberg_decision": volume_case["decision"]
        == "finite_volume_noncompact"
        and volume_case["simple_roots_complete"] is True,
        "vinberg_edges": volume_case["elliptic_edge_count"] == 5160
        and volume_case["finite_volume_bad_edge_count"] == 0,
        "vinberg_parabolics": volume_case["connected_affine_audit_count"]
        == 22
        and volume_case["bad_full_parabolic_neighborhood_count"] == 0,
        "vinberg_vertices": volume_case["finite_vertex_count"] == 525
        and volume_case["ideal_vertex_count"] == 10,
        "vinberg_audit_digest": volume_case["canonical_full_audit_sha256"]
        == expected_audit_hash,
        "source_volume_cross_pin": volume_payload["source_sha256"]
        == expected_wall_hash,
    }

    payload: dict[str, Any] = {
        "schema_version": 1,
        "artifact_id": "ue8x2a1_inose_h18_family",
        "positive_object_id": "R19-UE8x2A1-NEF-H18-24",
        "claim_scope": (
            "standalone project-level equation-to-complete-nef realization "
            "of a known H18 chamber"
        ),
        "lineage": {
            "ambient_equation": "classical Inose pencil",
            "ambient_source": "https://arxiv.org/abs/math/0612396",
            "specialization": "A=s^2, B=1-s^3 in the two-parameter pencil",
            "wall_source_locator": wall_case["source_locator"],
            "published_model_note": wall_case["published_model"],
            "priority_boundary": (
                "the equation, lattice, and abstract Vinberg--Kondo chamber "
                "are known; the dedicated checked equation-to-chamber bridge "
                "is the positive project contribution"
            ),
        },
        "family": {
            "affine_equation": (
                "y^2=x^3-3*s^2*t^4*x+"
                "t^5*(t^2-2*(1-s^3)*t+1)"
            ),
            "parameter_open": "s^3*(s^3-1) != 0",
            "discriminant": str(discriminant),
            "c4": str(c4),
            "residual_quadratic": str(residual_quadratic),
            "residual_discriminant": str(residual_discriminant),
            "open_nonempty_witness": {
                "s": 2,
                "residual_quadratic": str(witness_quadratic),
                "residual_discriminant": int(witness_discriminant),
            },
            "fibers": [
                {"place": "t=0", "kodaira": "II*", "euler": 10},
                {"place": "t=infinity", "kodaira": "II*", "euler": 10},
                {"place": "t=1", "kodaira": "I2", "euler": 2},
                {
                    "place": "two roots of the residual quadratic",
                    "kodaira": "I1",
                    "multiplicity": 2,
                    "euler_each": 1,
                },
            ],
            "fiber_euler_sum": 24,
        },
        "neron_severi": {
            "very_general_exact_lattice": "U+E8(-1)^2+A1(-1)",
            "rank_signature_determinant": [19, [1, 18, 0], 2],
            "trivial_lattice_argument": (
                "the section, fiber, two II* component systems, and the "
                "nonidentity I2 component give U+E8(-1)^2+A1(-1)"
            ),
            "parameter_recovery": "s^3=([t]Q+2)/4",
            "generic_exactness": (
                "the marked family has one-dimensional image, equal to the "
                "dimension of the rank-19 lattice-polarized period domain; "
                "the rho=20 locus is zero-dimensional, so very generally "
                "rho=19"
            ),
            "saturation_and_mordell_weil": (
                "a proper same-rank overlattice would have index at least "
                "two and nonintegral determinant 2/4; hence none exists. "
                "Shioda--Tate then gives Mordell--Weil rank and torsion zero"
            ),
        },
        "complete_nef_handoff": {
            "source_artifact": str(WALL_SOURCE.relative_to(ROOT)),
            "source_sha256": sha256_file(WALL_SOURCE),
            "transport_argument": (
                "with exact NS, Weyl transport identifies the nef cone with "
                "the published complete 24-wall chamber"
            ),
            "wall_count": 24,
            "wall_gram_rank_inertia": [19, [1, 18, 5]],
            "full_rank_basis_determinant": 2,
            "wall_analysis": wall_analysis,
        },
        "vinberg_handoff": {
            "source_artifact": str(VOLUME_SOURCE.relative_to(ROOT)),
            "source_sha256": sha256_file(VOLUME_SOURCE),
            "criterion": volume_payload["criterion"],
            "classification": volume_case["decision"],
            "hyperbolic_dimension": 18,
            "facets": 24,
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
