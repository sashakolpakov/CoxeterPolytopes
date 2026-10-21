#!/usr/bin/env python3
"""Certify a very-general full-Picard K3 carrier for K3-D4-H12/17.

The archived order-three D4-triality fold already supplies a complete
17-wall finite-volume H12 chamber.  This composition certificate does not
try to realize triality as an automorphism of a new surface.  Instead it

* freezes the saturated rank-13 fixed lattice and all 17 wall coordinates;
* gives a determinant-one decomposition as
  U+E8(-1)+A2(-1)+A1(-1);
* checks every primitive reflection and index-one wall span;
* gives an explicit primitive embedding in the K3 lattice; and
* replays the literal Vinberg endpoint and affine-cusp gates.

Consequently a seven-dimensional very-general lattice-polarized K3 family
has this lattice as its full Neron--Severi group and carries the chamber as
a mixed divisor-reflection subsystem.  This is not a full-nef or
surface-automorphism assertion.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "data/candidates/k3_rho15_h12_cyclic_fold_audit.json"
PRODUCTION = ROOT / "data/enumerations/k3_rank15_order3_h12_fold.json"
LINEAGE = (
    ROOT
    / "data/enumerations/k3_rank15_order3_h12_scharlau_kirschmer_audit.json"
)
OUTPUT = (
    ROOT / "data/enumerations/k3_d4_h12_full_picard_k3_carrier.json"
)

SCHEMA = "coxeter.k3_d4.h12_17.full_picard_k3_carrier.v1"
EXPECTED_SHA256 = {
    "candidate": "a93178362a063fc96c721e3dfd721655cf2093cde3ce7cdb4ef60ccc4e22a28d",
    "production": "54fea0698d48241cd50d6654ccc7ad687e227e9fa9d8908e011e7455782e097b",
    "lineage": "8263b22ff7691abb0fa0a75c18c4bad73ec67a60ab33e4283ff6c27f37a8e613",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(index)
        ]
        for index in range(matrix.rows)
    ]


def _vector_payload(vector: sp.Matrix) -> list[int | str]:
    return [
        int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
        for value in vector
    ]


def _smith(matrix: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(diagonal[index, index]))
        for index in range(min(diagonal.shape))
        if diagonal[index, index]
    ]


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector))


def _divisibility(gram: sp.Matrix, vector: sp.Matrix) -> int:
    return _content(gram * vector)


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = sp.Rational((root.T * gram * root)[0])
    return sp.simplify(
        sp.eye(gram.rows)
        - sp.Rational(2, square) * root * (root.T * gram)
    )


def _primitive_k3_embedding(target: sp.Matrix) -> dict[str, Any]:
    """Embed U+A2(-1)+A1(-1)+E8(-1) primitively in Lambda_K3."""

    U = sp.Matrix([[0, 1], [1, 0]])
    E8m = target[5:13, 5:13]
    ambient = sp.diag(U, U, U, E8m, E8m)
    embedding = sp.zeros(22, 13)

    # Target U -> first ambient U.
    embedding[0, 0] = 1
    embedding[1, 1] = 1
    # Target A2 uses adjacent simple roots 0,1 of the second E8.
    embedding[14, 2] = 1
    embedding[15, 3] = 1
    # Target A1 uses simple root 7, orthogonal to those two roots.
    embedding[21, 4] = 1
    # Target E8 -> first ambient E8.
    for index in range(8):
        embedding[6 + index, 5 + index] = 1

    unit_rows = [0, 1, 14, 15, 21, 6, 7, 8, 9, 10, 11, 12, 13]
    orthogonal_vectors = (embedding.T * ambient).nullspace()
    orthogonal_basis = sp.Matrix.hstack(*orthogonal_vectors)
    orthogonal_gram = sp.simplify(
        orthogonal_basis.T * ambient * orthogonal_basis
    )
    checks = {
        "ambient_even_unimodular_signature_3_19": (
            ambient.det() == -1
            and exact_inertia(ambient) == (3, 19, 0)
            and all(ambient[index, index] % 2 == 0 for index in range(22))
        ),
        "E8_block_unimodular_negative_definite": (
            E8m.det() == 1 and exact_inertia(E8m) == (0, 8, 0)
        ),
        "A2_and_A1_are_orthogonal_basis_subdiagram": (
            E8m.extract([0, 1, 7], [0, 1, 7])
            == target.extract([2, 3, 4], [2, 3, 4])
        ),
        "embedding_isometric": embedding.T * ambient * embedding == target,
        "embedding_primitive_Smith_all_one": _smith(embedding) == [1] * 13,
        "visible_unit_minor": embedding.extract(unit_rows, range(13)).det() == 1,
        "orthogonal_complement_rank_signature": (
            len(orthogonal_vectors) == 9
            and orthogonal_basis.rank() == 9
            and exact_inertia(orthogonal_gram) == (2, 7, 0)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(
            [key for key, value in checks.items() if not value]
        )
    return {
        "ambient_lattice": "Lambda_K3=U^3+E8(-1)^2",
        "embedding_description": (
            "U enters the first ambient U; the target E8(-1) enters the "
            "first ambient E8; A2(-1)+A1(-1) is the primitive basis-node "
            "subdiagram {0,1}+{7} in the second ambient E8"
        ),
        "embedding_matrix_columns": _matrix_payload(embedding),
        "visible_unit_minor_rows_zero_based": unit_rows,
        "orthogonal_complement_rational_basis_columns": _matrix_payload(
            orthogonal_basis
        ),
        "orthogonal_complement_gram": _matrix_payload(orthogonal_gram),
        "orthogonal_complement_signature": [2, 7, 0],
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "candidate": CANDIDATE,
        "production": PRODUCTION,
        "lineage": LINEAGE,
    }
    for key, path in paths.items():
        actual = _sha256(path)
        if actual != EXPECTED_SHA256[key]:
            raise ArithmeticError(
                f"pinned {key} archive changed: {actual}"
            )

    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    production = json.loads(PRODUCTION.read_text(encoding="utf-8"))
    lineage = json.loads(LINEAGE.read_text(encoding="utf-8"))
    if not all(candidate["checks"].values()):
        raise ArithmeticError("the independent fold certificate is not green")
    if production["all_checks_pass"] is not True:
        raise ArithmeticError("the production fold certificate is not green")
    if not all(lineage["checks"].values()):
        raise ArithmeticError("the lineage/full-root certificate is not green")

    gram_fixed = sp.Matrix(candidate["matrices"]["fixed_lattice_gram_13"])
    roots_fixed = sp.Matrix(
        candidate["matrices"]["folded_normal_coordinates_13x17"]
    )
    orbit_gram = sp.Matrix(candidate["matrices"]["orbit_gram_17"])

    # The lineage matrix C satisfies C^T*(-G)*C=U+A2(+)+A1(+)+E8(+).
    # Negating the first U basis vector changes -U back to U, yielding the
    # negative-definite root blocks in the K3 convention.
    decomposition = sp.Matrix(
        lineage["exact_lattice_identity"]["candidate_decomposition_matrix"]
    )
    sign_change = sp.diag(-1, *([1] * 12))
    standard_to_fixed = decomposition * sign_change
    gram_standard = sp.simplify(
        standard_to_fixed.T * gram_fixed * standard_to_fixed
    )

    U = sp.Matrix([[0, 1], [1, 0]])
    A2m = sp.Matrix([[-2, 1], [1, -2]])
    A1m = sp.Matrix([[-2]])
    E8m = -sp.Matrix(
        [
            [2 if i == j else -1 if (min(i, j), max(i, j)) in {
                (0, 1), (1, 2), (2, 3), (3, 4),
                (4, 5), (4, 6), (6, 7),
            } else 0 for j in range(8)]
            for i in range(8)
        ]
    )
    expected_standard = sp.diag(U, A2m, A1m, E8m)
    fixed_to_standard = standard_to_fixed.inv()
    roots_standard = sp.simplify(fixed_to_standard * roots_fixed)

    strict_fixed = sp.Matrix(
        candidate["order_three_H12_fold"]["strict_interior"][
            "basis_coordinates"
        ]
    )
    strict_standard = sp.simplify(fixed_to_standard * strict_fixed)
    expected_pairings = candidate["order_three_H12_fold"]["strict_interior"][
        "wall_pairings"
    ]

    records: list[dict[str, Any]] = []
    types: Counter[tuple[int, int]] = Counter()
    pairings: list[int] = []
    for index in range(roots_standard.cols):
        root = roots_standard[:, index]
        source_root = roots_fixed[:, index]
        reflection = _reflection(gram_standard, root)
        square = int((root.T * gram_standard * root)[0])
        divisibility = _divisibility(gram_standard, root)
        pairing = int((strict_standard.T * gram_standard * root)[0])
        checks = {
            "transport_integral": all(value.is_Integer for value in root),
            "primitive": _content(root) == 1,
            "square_preserved": (
                square == int((source_root.T * gram_fixed * source_root)[0])
            ),
            "orientation_strict": pairing > 0,
            "reflection_integral": all(value.is_Integer for value in reflection),
            "reflection_unimodular": reflection.det() == -1,
            "reflection_involutive": reflection**2 == sp.eye(13),
            "reflection_preserves_gram": (
                reflection.T * gram_standard * reflection == gram_standard
            ),
            "reflection_negates_root": reflection * root == -root,
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"wall {index + 1} failed: "
                f"{[key for key, value in checks.items() if not value]}"
            )
        types[(square, divisibility)] += 1
        pairings.append(pairing)
        records.append(
            {
                "facet_one_based": index + 1,
                "ambient_wall_orbit": candidate["order_three_H12_fold"][
                    "folded_reflections"
                ][index]["ambient_wall_orbit"],
                "primitive_fixed_basis_coordinates": _vector_payload(source_root),
                "primitive_standard_Picard_coordinates": _vector_payload(root),
                "square": square,
                "divisibility_in_full_Picard": divisibility,
                "strict_interior_pairing": pairing,
                "checks": checks,
            }
        )

    transported_gram = sp.simplify(
        roots_standard.T * gram_standard * roots_standard
    )
    span_smith = _smith(roots_standard.T)
    lattice_checks = {
        "saturated_fixed_rank13_H12": (
            gram_fixed.rows == 13
            and exact_inertia(gram_fixed) == (1, 12, 0)
            and candidate["checks"]["fixed_lattice_saturated_rank13_H12"]
        ),
        "determinant6_smith_1x12_6": (
            gram_fixed.det() == gram_standard.det() == 6
            and _smith(gram_standard) == [1] * 12 + [6]
        ),
        "even": all(gram_standard[index, index] % 2 == 0 for index in range(13)),
        "determinant_one_standard_decomposition": (
            abs(int(standard_to_fixed.det())) == 1
            and gram_standard == expected_standard
        ),
        "all_17_wall_transports": all(
            all(record["checks"].values()) for record in records
        ),
        "facet_Gram_preserved": transported_gram == orbit_gram,
        "walls_span_full_Picard_index_one": (
            roots_standard.rank() == 13 and span_smith == [1] * 13
        ),
        "mirror_type_histogram": dict(types)
        == {(-6, 3): 2, (-2, 1): 13, (-2, 2): 2},
        "strict_square1110": (
            (strict_standard.T * gram_standard * strict_standard)[0] == 1110
        ),
        "strict_pairings_literal": pairings == expected_pairings,
    }
    if not all(lattice_checks.values()):
        raise ArithmeticError(
            [key for key, value in lattice_checks.items() if not value]
        )

    embedding = _primitive_k3_embedding(gram_standard)

    literal = candidate["literal_vinberg_coxiter_prop_1_10"]
    production_literal = production["literal_coxiter_prop_1_10"]
    volume_checks = {
        "independent_finite_volume_noncompact": literal[
            "finite_volume_noncompact"
        ],
        "93_finite_3_ideal_0_hyperideal": (
            candidate["order_three_H12_fold"]["finite_vertex_count"] == 93
            and candidate["order_three_H12_fold"]["ideal_vertex_count"] == 3
            and candidate["order_three_H12_fold"]["hyperideal_vertex_count"] == 0
        ),
        "literal_623_edges_each_two_endpoints": (
            literal["elliptic_edge_count"] == 623
            and literal["extension_histogram"] == {"2": 623}
            and literal["bad_edges"] == []
        ),
        "all_affine_neighborhoods_full_cusps": literal[
            "full_zero_neighbourhoods_match_ideal_vertices"
        ],
        "production_literal_agrees": (
            production_literal["spherical_rank_11_edge_count"] == 623
            and production_literal["spherical_rank_12_vertex_count"] == 93
            and production_literal["parabolic_rank_11_vertex_count"] == 3
            and production_literal["edge_extension_histogram"] == {"2": 623}
            and production_literal["bad_edges"] == []
        ),
        "optimized_production_agrees": (
            production["optimized_vinberg"]["elliptic_edge_count"] == 623
            and production["optimized_vinberg"]["finite_vertex_count"] == 93
            and production["optimized_vinberg"]["ideal_vertex_count"] == 3
            and production["optimized_vinberg"]["finite_volume_bad_edges"] == []
        ),
    }
    if not all(volume_checks.values()):
        raise ArithmeticError(
            [key for key, value in volume_checks.items() if not value]
        )

    full_group_checks = {
        "lineage_all_checks_green": all(lineage["checks"].values()),
        "complete_simple_system": (
            lineage["full_reflective_root_audit"][
                "finite_volume_stopping_step"
            ]["conclusion"]
            == (
                "The 17 candidate roots are the complete simple-root system "
                "of the full reflection group W(S)."
            )
        ),
        "known_Scharlau_Kirschmer_entry": (
            lineage["source_freeze"]["selected_entry_one_based"] == 6
        ),
        "not_new_abstract_chamber": (
            lineage["novelty_conclusion"][
                "new_abstract_hyperbolic_coxeter_polytope"
            ]
            is False
        ),
    }

    sources = {
        key: {
            "path": str(path.relative_to(ROOT)),
            "sha256": EXPECTED_SHA256[key],
        }
        for key, path in paths.items()
    }
    checks = {
        "source_hashes_pinned": all(
            _sha256(path) == EXPECTED_SHA256[key]
            for key, path in paths.items()
        ),
        "lattice_and_walls_green": all(lattice_checks.values()),
        "primitive_K3_embedding_green": all(embedding["checks"].values()),
        "literal_volume_green": all(volume_checks.values()),
        "full_reflection_lineage_green": all(full_group_checks.values()),
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "A seven-dimensional very-general family of projective K3 "
            "surfaces has full NS=U+E8(-1)+A2(-1)+A1(-1) and carries the "
            "known finite-volume noncompact K3-D4-H12/17 chamber as 17 "
            "primitive integral divisor mirrors generating full Picard."
        ),
        "sources": sources,
        "saturated_fixed_lattice_handoff": {
            "source_ambient_Picard": "U+E8(-1)+D4(-1)+A1(-1)",
            "source_action": "order-three D4 triality on the divisor lattice",
            "source_action_surface_scope": candidate[
                "surface_realization_boundary"
            ],
            "fixed_lattice_gram": _matrix_payload(gram_fixed),
            "saturation_check": candidate["checks"][
                "fixed_lattice_saturated_rank13_H12"
            ],
            "standard_to_fixed_determinant_one_matrix": _matrix_payload(
                standard_to_fixed
            ),
        },
        "full_Picard_lattice": {
            "name": "S=U+E8(-1)+A2(-1)+A1(-1)",
            "ordered_standard_decomposition": "U+A2(-1)+A1(-1)+E8(-1)",
            "rank": 13,
            "signature": [1, 12, 0],
            "determinant": 6,
            "smith_invariants": _smith(gram_standard),
            "standard_gram": _matrix_payload(gram_standard),
            "primitive_K3_embedding": embedding,
            "period_consequence": {
                "K3_moduli_dimension": 7,
                "very_general_full_Picard": "NS(X)=S",
                "reason": (
                    "the displayed primitive embedding defines a nonempty "
                    "S-polarized K3 period domain; away from the countable "
                    "Noether--Lefschetz hyperplanes no extra algebraic class occurs"
                ),
            },
            "checks": lattice_checks,
        },
        "transported_chamber": {
            "id": "K3-D4-H12",
            "hyperbolic_dimension": 12,
            "facet_count": 17,
            "wall_span_rank": roots_standard.rank(),
            "wall_span_smith": span_smith,
            "wall_span_index_in_full_Picard": math.prod(span_smith),
            "primitive_mirror_type_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(types.items())
            },
            "strict_interior_coordinates": _vector_payload(strict_standard),
            "strict_interior_square": 1110,
            "strict_wall_pairings": pairings,
            "minimum_wall_pairing": min(pairings),
            "wall_records": records,
            "transported_facet_gram_sha256": _payload_sha256(
                _matrix_payload(transported_gram)
            ),
        },
        "finite_volume": {
            "criterion": (
                "literal CoxIter Proposition 1.10 / Vinberg edge-endpoint "
                "criterion plus all affine zero-neighbourhood cusp tests"
            ),
            "classification": "finite_volume_noncompact_Coxeter_polytope",
            "finite_vertices": 93,
            "cusps": 3,
            "hyperideal_extreme_rays": 0,
            "elliptic_rank11_edges": 623,
            "endpoint_type_histogram": literal["endpoint_type_histogram"],
            "checks": volume_checks,
        },
        "reflection_group_lineage": {
            "full_chamber": True,
            "index": 1,
            "known_source": "Scharlau--Kirschmer List12/Gram12 entry 6",
            "identity_scope": lineage["identity_hierarchy"][
                "highest_defensible_statement"
            ],
            "checks": full_group_checks,
        },
        "scope": {
            "positive_gain": (
                "new project-level very-general complex-K3 full-Picard "
                "divisor-reflection carrier of a known chamber"
            ),
            "not_new_abstract_polytope": True,
            "not_full_nef": (
                "the two primitive square-minus-six mirrors are mixed divisor "
                "cuts, so the 17-wall reflection chamber is not asserted to "
                "be the generic K3 nef cone"
            ),
            "triality_not_surface_automorphism": (
                "the order-three graph/lattice action is used only to derive "
                "the fixed lattice and walls; it is not claimed as an "
                "automorphism of the very-general carrier"
            ),
            "reflections_not_surface_automorphisms": (
                "integral Picard reflections are not asserted to extend to "
                "Hodge isometries or K3 automorphisms"
            ),
            "effectivity": (
                "after Weyl re-marking, Riemann--Roch orients the fifteen "
                "square-minus-two mirrors effectively; irreducibility is not claimed"
            ),
            "equation_level": False,
            "worldwide_priority_claim": "none",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["payload_sha256"] = _payload_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    payload = build_certificate()
    if not payload["all_checks_pass"]:
        raise ArithmeticError("K3-D4 H12 carrier certificate failed")
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(
            encoding="utf-8"
        ) != rendered:
            raise ArithmeticError("archived K3-D4 H12 carrier changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={payload['payload_sha256']}")


if __name__ == "__main__":
    main()
