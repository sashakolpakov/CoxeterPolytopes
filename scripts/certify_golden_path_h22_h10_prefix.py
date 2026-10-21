#!/usr/bin/env python3
"""Freeze the exact H<=10 Vinberg prefix of the rank-23 golden path host."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden, ONE, PHI
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.search.golden_path_host import (
    COEFFICIENT,
    DETERMINANT,
    DIMENSION,
    HEIGHT_BOUND,
    MIRROR_NORM_PROOF,
    NATIVE_RANK,
    P101,
    P5,
    WEYL_SQUARE,
    canonical_mirror_norms,
    cartan_matrix,
    golden_path_h10_audit,
    lattice_gram,
    pair,
    strict_interior_point,
    weyl_vector,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/golden_path_h22_h10_prefix.json"


def _golden(value: Golden | int) -> dict[str, str]:
    return Golden.coerce(value).to_json()


def _vector(vector) -> list[dict[str, str]]:
    return [_golden(value) for value in vector]


def _matrix(matrix) -> list[list[dict[str, str]]]:
    return [_vector(row) for row in matrix]


def _candidate(candidate) -> dict:
    return {
        "height": _golden(candidate.height),
        "norm": _golden(candidate.norm),
        "control": _golden(candidate.control),
        "vector": _vector(candidate.vector),
        "simple_A22_pairings": _vector(candidate.simple_pairings),
        "reflection_coefficients_on_basis": _vector(candidate.reflection_coefficients),
        "primitive_unit_control": abs(candidate.control.norm()) == 1,
    }


def _geometry(geometry) -> dict:
    return {
        "facet_count": geometry.facet_count,
        "finite_vertex_count": len(geometry.finite_vertices),
        "cusp_count": len(geometry.cusps),
        "hyperideal_ray_count": len(geometry.hyperideal_rays),
        "opposite_sheet_ray_count": len(geometry.opposite_sheet_rays),
        "finite_vertices": [list(face) for face in geometry.finite_vertices],
        "cusps": [list(face) for face in geometry.cusps],
        "hyperideal_rays": [list(face) for face in geometry.hyperideal_rays],
        "finite_volume": not geometry.hyperideal_rays and not geometry.opposite_sheet_rays,
        "compact": (
            not geometry.cusps
            and not geometry.hyperideal_rays
            and not geometry.opposite_sheet_rays
        ),
    }


def _inertia(value) -> dict[str, int]:
    return {
        "positive": value.positive,
        "negative": value.negative,
        "zero": value.zero,
    }


def build_payload() -> dict:
    audit = golden_path_h10_audit()
    gram = lattice_gram()
    rho = weyl_vector()
    witness = audit.hyperideal_witness
    controls = [
        {
            "height": _golden(row.height),
            "norm": _golden(row.norm),
            "control": _golden(row.control),
            "spatial_target": _golden(row.spatial_target),
            "reflection_target": _golden(row.reflection_target),
            "alphabet_size": row.alphabet_size,
            "recursion_nodes": row.recursion_nodes,
            "energy_solution_count": row.energy_solution_count,
            "integral_candidate_count": row.integral_candidate_count,
        }
        for row in audit.controls
    ]
    decisions = [
        {
            "candidate": _candidate(row.candidate),
            "accepted": row.accepted,
            "first_positive_wall": row.first_positive_wall,
            "first_positive_pairing": (
                None
                if row.first_positive_pairing is None
                else _golden(row.first_positive_pairing)
            ),
        }
        for row in audit.decisions
    ]
    return {
        "schema": "coxeter.golden_path_h22_h10_prefix.v1",
        "scope": {
            "field": "Q(sqrt(5))",
            "ring": "Z[phi]",
            "native_rank": NATIVE_RANK,
            "hyperbolic_dimension": DIMENSION,
            "height_bound_inclusive": _golden(HEIGHT_BOUND),
            "claim": (
                "Exact host arithmetic, complete primitive mirror-norm candidate "
                "alphabet, and complete chronological Vinberg prefix through H<=10."
            ),
            "not_claimed": [
                "Vinberg termination above height 10",
                "finite index of the reflection subgroup",
                "a compact or finite-volume Coxeter chamber",
                "golden-field nonexistence in H22--H29",
            ],
        },
        "host": {
            "coefficient_c": _golden(COEFFICIENT),
            "cartan_definition": (
                "A00=...=A22,22=2; A01=-1; A10=-c; "
                "Ai,i+1=Ai+1,i=-1 for 1<=i<22"
            ),
            "cartan_matrix": _matrix(cartan_matrix()),
            "symmetrizer_diagonal": [_golden(COEFFICIENT)]
            + [_golden(ONE) for _ in range(22)],
            "symmetric_gram": _matrix(gram),
            "gram_equals_symmetrizer_times_cartan": True,
            "determinant": _golden(determinant(gram)),
            "expected_determinant": _golden(DETERMINANT),
            "determinant_factorization": {
                "identity": "-84-172phi=-4(2+phi)(4+13phi)",
                "p5_generator": _golden(P5),
                "p101_canonical_totally_positive_generator": _golden(P101),
                "p101_association": (
                    "phi(4+13phi)=13+17phi=phi^2(9+4phi)"
                ),
                "p5_absolute_norm": int(P5.norm()),
                "p101_absolute_norm": int(P101.norm()),
            },
            "identity_inertia": _inertia(inertia(gram, embedding=0)),
            "conjugate_inertia": _inertia(inertia(gram, embedding=1)),
            "integral_reflections_in_basis": True,
            "minimal_vinberg_field": {
                "field": "Q(sqrt(5))",
                "proof": (
                    "The two-cycle Cartan product A01*A10=c=2+2phi is "
                    "nonrational, while all Cartan data lie in Q(sqrt(5))."
                ),
            },
            "anisotropic": {
                "over_Qsqrt5": True,
                "proof": (
                    "A nonzero K-isotropic vector would remain isotropic after "
                    "conjugation, contradicting positive definiteness at the "
                    "second real embedding."
                ),
                "arithmetic_orthogonal_host_cocompact": True,
                "reflection_subgroup_finite_index": "OPEN",
            },
        },
        "weyl_height": {
            "rho": _vector(rho),
            "rho_square": _golden(pair(rho, rho)),
            "expected_rho_square": _golden(WEYL_SQUARE),
            "normalization": (
                "B*rho=rho^2*e0, so B^-1*e0=rho/rho^2; "
                "H(r)=-rho^2*x^2/q for r=(x,y), q=B(r,r)."
            ),
            "strict_prefix_interior": _vector(strict_interior_point()),
        },
        "mirror_norms": {
            "candidate_count": len(canonical_mirror_norms()),
            "canonical_candidates": [_golden(value) for value in canonical_mirror_norms()],
            "local_exponents": {
                "inert_2": [0, 1, 2, 3],
                "ramified_p5": [0, 1],
                "split_p101": [0, 1],
            },
            "proof": MIRROR_NORM_PROOF,
            "scope_word": "candidate alphabet; empty local classes are retained",
        },
        "enumeration": {
            "control_count": len(audit.controls),
            "sum_recursion_nodes": sum(row.recursion_nodes for row in audit.controls),
            "sum_energy_solutions": sum(row.energy_solution_count for row in audit.controls),
            "integral_candidate_count": len(audit.candidates),
            "control_stream_sha256": audit.control_stream_sha256,
            "candidate_stream_sha256": audit.candidate_stream_sha256,
            "method": (
                "Exact two-place control boxes; z_i^2<=2T alphabet; "
                "A_d inverse Schur pruning for every prefix; exact A22 inverse "
                "reconstruction; basis reflection-integrality; chronological acute gate."
            ),
            "controls": controls,
            "decisions": decisions,
        },
        "accepted_prefix": {
            "height_zero_wall_count": 22,
            "positive_wall_count": len(audit.accepted_positive_roots),
            "total_wall_count": len(audit.walls),
            "accepted_positive_roots": [
                _candidate(row) for row in audit.accepted_positive_roots
            ],
            "first_positive_wall_is_seed_hyperplane": True,
            "seed_hyperplane_explanation": (
                "The q=4 root (phi-1)e0 is a unit multiple of e0, hence "
                "the original terminal dotted path hyperplane, not a new mirror."
            ),
            "new_wall_count": 1,
            "new_wall_description": (
                "q=2, x=1, y1=1+phi, y2=0, y3=...=y22=-phi; "
                "pairings with the seed,e2,e3,e22 are -2,-1,-phi,-phi."
            ),
            "rejected_simultaneous_root": (
                "q=2 root e0+e1 is rejected because its pairing with "
                "the earlier (phi-1)e0 seed wall is +2phi."
            ),
            "walls": [_vector(wall) for wall in audit.walls],
        },
        "geometry": {
            "23_wall_seed_cone": _geometry(audit.seed_geometry),
            "24_wall_H10_prefix_cone": _geometry(audit.geometry),
            "explicit_hyperideal_extreme_ray": {
                "active_facets": list(witness.active_facets),
                "ray": _vector(witness.ray),
                "norm": _golden(witness.norm),
                "wall_pairings": _vector(witness.wall_pairings),
            },
            "bounded_verdict": (
                "The complete H<=10 prefix cone is infinite-volume: it has "
                "41 exact hyperideal extreme rays.  Later walls may cut them."
            ),
        },
        "engine_lineage": {
            "existing_pinned_alvin_inspected": True,
            "alvin_used_for_this_host": False,
            "reason": (
                "The pinned AlVin command interface accepts diagonal forms. "
                "This host is a non-diagonal integral O_K lattice, and no "
                "O_K-unimodular diagonalization preserving its lattice was "
                "constructed.  Feeding only a K-diagonalization would change "
                "the enumerated integral lattice.  The bounded prefix therefore "
                "uses the repository's native exact two-place/Schur/DD engines."
            ),
        },
        "verdict": {
            "new_compact_example": False,
            "new_exact_Vinberg_wall": True,
            "complete_through_height_10": True,
            "reflectivity": "OPEN",
            "golden_H22_H29_nonexistence": "OPEN",
        },
    }


def encoded_payload() -> bytes:
    return (json.dumps(build_payload(), indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    encoded = encoded_payload()
    if arguments.verify:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != encoded:
            raise SystemExit("FAIL: frozen golden path H10 JSON differs from recomputation")
        status = "PASS"
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(encoded)
        status = "WROTE"
    payload = build_payload()
    print(
        json.dumps(
            {
                "status": status,
                "output": str(OUTPUT),
                "sha256": sha256(encoded).hexdigest(),
                "controls": payload["enumeration"]["control_count"],
                "candidates": payload["enumeration"]["integral_candidate_count"],
                "accepted_walls": payload["accepted_prefix"]["total_wall_count"],
                "hyperideal_rays": payload["geometry"]["24_wall_H10_prefix_cone"][
                    "hyperideal_ray_count"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
