#!/usr/bin/env python3
"""Certify an equation-level full-Picard K3 carrier for X2 H10/37.

This is a composition certificate.  It pins

* the exact 37-wall rational-involution section of Shimada's X_2 chamber;
* the determinant-one identification of its fixed divisor lattice with
  F=U+E8(-2)+<-2>; and
* the explicit nine-parameter quadratic-base-change K3 family whose very
  general member has full Neron--Severi lattice F.

All 37 primitive wall normals are transported to the geometric standard
F-basis, checked to span F, and checked reflection-by-reflection.  The
finite-volume statement is inherited only after replaying every stored
literal Vinberg endpoint and affine-neighbourhood check.  No full-nef or
surface-automorphism assertion is made for the very-general carrier.
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
X2 = (
    ROOT
    / "data/enumerations/shimada_aut3singk3_x2_rational_involution_h10_37.json"
)
X2_STOCK = (
    ROOT / "data/enumerations/shimada_aut3singk3_x2_478_wall_stock.json"
)
BRIDGE = ROOT / "data/enumerations/a15_11_h10_21_x2_k3_carrier.json"
FAMILY = ROOT / "data/enumerations/a15_11_ue8m2a1_basechange_family.json"
OUTPUT = (
    ROOT
    / "data/enumerations/x2_h10_37_ue8m2a1_full_picard_carrier.json"
)

SCHEMA = "coxeter.x2.h10_37.ue8m2a1_full_picard_carrier.v1"
EXPECTED_SHA256 = {
    "x2": "12914344403b7d3d6b69a2866c4150f0a0420a64914e4bc8584fa372799599a7",
    "x2_stock": "10d1c767d21d7194c6ccb603c5dff29520201954b4a07fd3c5f214d452ea2f74",
    "bridge": "5f241c3659a79cf674d09558f72a5bc075d50b9fc6d87444ce8d253f1a4ab90b",
    "family": "1a5902092e3150bd380af19916ec5a6733392404a6cf2ecdba523f82b488532c",
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
    """Embed F=U+E8(-2)+<-2> primitively in Lambda_K3."""

    U = sp.Matrix([[0, 1], [1, 0]])
    E8m = sp.simplify(target[2:10, 2:10] / 2)
    ambient = sp.diag(U, U, U, E8m, E8m)

    embedding = sp.zeros(22, 11)
    embedding[0, 0] = 1
    embedding[1, 1] = 1
    for index in range(8):
        embedding[6 + index, 2 + index] = 1
        embedding[14 + index, 2 + index] = 1
    # e-f in the second U is primitive and has square -2.
    embedding[2, 10] = 1
    embedding[3, 10] = -1

    complement = sp.zeros(22, 11)
    complement[4, 0] = 1
    complement[5, 1] = 1
    for index in range(8):
        complement[6 + index, 2 + index] = 1
        complement[14 + index, 2 + index] = -1
    # e+f is orthogonal to e-f and has square +2.
    complement[2, 10] = 1
    complement[3, 10] = 1

    complement_gram = sp.diag(U, target[2:10, 2:10], sp.Matrix([[2]]))
    combined = embedding.row_join(complement)
    checks = {
        "ambient_even_unimodular_signature_3_19": (
            ambient.det() == -1
            and exact_inertia(ambient) == (3, 19, 0)
            and all(ambient[index, index] % 2 == 0 for index in range(22))
        ),
        "embedding_isometric": embedding.T * ambient * embedding == target,
        "embedding_primitive_Smith_all_one": _smith(embedding) == [1] * 11,
        "visible_unit_minor": (
            embedding.extract(
                [0, 1, 6, 7, 8, 9, 10, 11, 12, 13, 2], range(11)
            ).det()
            == 1
        ),
        "complement_isometric": (
            complement.T * ambient * complement == complement_gram
        ),
        "orthogonal": embedding.T * ambient * complement == sp.zeros(11),
        "gluing_index_equals_discriminant_order": (
            abs(int(combined.det())) == abs(int(target.det())) == 512
        ),
        "opposite_discriminant_form_visible": (
            target[:10, :10] == complement_gram[:10, :10]
            and target[10, 10] == -2
            and complement_gram[10, 10] == 2
            and _smith(target) == _smith(complement_gram)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(
            [key for key, value in checks.items() if not value]
        )
    return {
        "ambient_lattice": "Lambda_K3=U^3+E8(-1)^2",
        "embedding_description": (
            "U is the first ambient U; E8(-2) is diagonal in the two "
            "E8(-1) factors; <-2> is generated by e-f in the second U"
        ),
        "embedding_matrix_columns": _matrix_payload(embedding),
        "orthogonal_complement": "U+E8(-2)+<2>",
        "orthogonal_complement_gram": _matrix_payload(complement_gram),
        "orthogonal_complement_matrix_columns": _matrix_payload(complement),
        "combined_embedding_complement_index": abs(int(combined.det())),
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "x2": X2,
        "x2_stock": X2_STOCK,
        "bridge": BRIDGE,
        "family": FAMILY,
    }
    for key, path in paths.items():
        actual = _sha256(path)
        if actual != EXPECTED_SHA256[key]:
            raise ArithmeticError(
                f"pinned {key} archive changed: {actual}"
            )

    x2 = json.loads(X2.read_text(encoding="utf-8"))
    stock = json.loads(X2_STOCK.read_text(encoding="utf-8"))
    bridge = json.loads(BRIDGE.read_text(encoding="utf-8"))
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    if not all(
        document["all_checks_pass"] for document in (x2, stock, family)
    ):
        raise ArithmeticError("an upstream theorem archive is not green")
    if bridge["all_checks_pass"] is not True:
        raise ArithmeticError("the pinned lattice bridge is not green")

    fixed = x2["saturated_fixed_divisor_lattice"]
    gram_fixed = sp.Matrix(fixed["gram"])
    gram_standard = sp.Matrix(
        family["fixed_locus_and_invariant_lattice"]["gram"]
    )
    standard_to_fixed = sp.Matrix(
        bridge["lattice_bridge"]["standard_to_actual_X2_fixed_Q"]
    )
    fixed_to_standard = standard_to_fixed.inv()

    fixed_rows = sp.Matrix(fixed["basis_rows_in_published_NS_basis"])
    gram_x2 = sp.Matrix(stock["surface"]["published_neron_severi_gram"])
    ample_x2 = sp.Matrix(stock["surface"]["published_ample_vector_a2"])
    ample_fixed = fixed_rows.T.gauss_jordan_solve(ample_x2)[0]
    ample_standard = sp.simplify(fixed_to_standard * ample_fixed)

    wall_records: list[dict[str, Any]] = []
    standard_roots: list[sp.Matrix] = []
    types: Counter[tuple[int, int]] = Counter()
    pairings: list[int] = []
    for source in fixed["mirror_records"]:
        root_fixed = sp.Matrix(source["fixed_basis_coordinates"])
        root_standard = sp.simplify(fixed_to_standard * root_fixed)
        reflection = _reflection(gram_standard, root_standard)
        square = int((root_standard.T * gram_standard * root_standard)[0])
        divisibility = _divisibility(gram_standard, root_standard)
        pairing = int(
            (ample_standard.T * gram_standard * root_standard)[0]
        )
        checks = {
            "transport_integral": all(value.is_Integer for value in root_standard),
            "primitive": _content(root_standard) == 1,
            "square_preserved": square == source["square"],
            "divisibility_preserved": divisibility == source["divisibility"],
            "orientation_preserved_by_square12_interior": pairing > 0,
            "reflection_integral": all(value.is_Integer for value in reflection),
            "reflection_unimodular": reflection.det() == -1,
            "reflection_involutive": reflection**2 == sp.eye(11),
            "reflection_preserves_gram": (
                reflection.T * gram_standard * reflection == gram_standard
            ),
            "reflection_negates_root": reflection * root_standard == -root_standard,
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"wall {source['facet_one_based']} failed: "
                f"{[key for key, value in checks.items() if not value]}"
            )
        types[(square, divisibility)] += 1
        pairings.append(pairing)
        standard_roots.append(root_standard)
        wall_records.append(
            {
                "facet_one_based": source["facet_one_based"],
                "primitive_fixed_basis_coordinates": source[
                    "fixed_basis_coordinates"
                ],
                "primitive_geometric_standard_F_coordinates": _vector_payload(
                    root_standard
                ),
                "square": square,
                "divisibility_in_full_Picard": divisibility,
                "square12_interior_pairing": pairing,
                "checks": checks,
            }
        )

    roots_by_rows = sp.Matrix.hstack(*standard_roots).T
    transported_gram = sp.Matrix(
        [
            [
                (left.T * gram_standard * right)[0]
                for right in standard_roots
            ]
            for left in standard_roots
        ]
    )
    expected_gram = sp.Matrix(x2["rational_involution_fold"]["facet_gram"])
    span_smith = _smith(roots_by_rows)

    volume = x2["literal_vinberg_volume_certificate"]
    volume_checks = {
        "classification_exact": (
            volume["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
        ),
        "every_elliptic_edge_has_two_endpoints": volume["checks"][
            "every_edge_has_exactly_two_finite_or_ideal_endpoints"
        ],
        "every_affine_zero_neighborhood_is_full_parabolic": volume["checks"][
            "all_full_affine_zero_neighborhoods_are_rank_9_parabolic"
        ],
        "zero_bad_edges": volume["bad_edge_count"] == 0,
        "zero_bad_parabolic_neighborhoods": (
            volume["bad_parabolic_neighborhood_count"] == 0
        ),
        "noncompact_cusps_exist": volume["ideal_vertex_count"] == 261,
        "complete_rank9_rank10_enumeration": volume["checks"][
            "all_elliptic_rank_9_edges_and_rank_10_vertices_enumerated"
        ],
    }
    if not all(volume_checks.values()):
        raise ArithmeticError("the inherited literal Vinberg gate changed")

    family_checks = {
        "explicit_family_is_green": family["all_checks_pass"],
        "nine_parameters": family["family"]["parameter_count"] == 9,
        "nonempty_open_witness": bool(
            family["family"]["nonempty_open_witness"]
        ),
        "very_general_full_NS_literal": (
            family["generic_full_picard"]["claim"]
            == "NS(X)=U+E8(-2)+<-2> for a very general parameter"
        ),
        "same_standard_gram": (
            gram_standard == sp.Matrix(fixed["standard_model_gram"])
        ),
        "same_nikulin_invariants": (
            family["fixed_locus_and_invariant_lattice"]["nikulin_invariants"]
            == fixed["nikulin_invariants"]
            == {"r": 11, "a": 9, "delta": 1}
        ),
        "full_moduli_dimension": (
            family["generic_full_picard"]["parameter_moduli_dimension"] == 9
            and family["generic_full_picard"][
                "two_elementary_pair_moduli_dimension"
            ]
            == 9
        ),
    }
    if not all(family_checks.values()):
        raise ArithmeticError("the explicit-family handoff changed")

    transport_checks = {
        "fixed_lattice_signature_1_10": (
            exact_inertia(gram_fixed) == (1, 10, 0)
            and exact_inertia(gram_standard) == (1, 10, 0)
        ),
        "fixed_lattice_even_determinant512": (
            gram_fixed.det() == gram_standard.det() == 512
            and all(gram_standard[index, index] % 2 == 0 for index in range(11))
        ),
        "unimodular_standard_to_fixed_isometry": (
            abs(int(standard_to_fixed.det())) == 1
            and standard_to_fixed.T * gram_fixed * standard_to_fixed
            == gram_standard
        ),
        "published_ample_lies_in_fixed_lattice": (
            fixed_rows.T * ample_fixed == ample_x2
        ),
        "square12_positive_interior": (
            (ample_standard.T * gram_standard * ample_standard)[0] == 12
            and min(pairings) > 0
        ),
        "all_37_transport_checks": all(
            all(record["checks"].values()) for record in wall_records
        ),
        "transported_Gram_literal": transported_gram == expected_gram,
        "all_37_roots_span_full_Picard": (
            roots_by_rows.rank() == 11 and span_smith == [1] * 11
        ),
        "mirror_type_histogram": dict(types)
        == {(-4, 2): 10, (-2, 1): 6, (-2, 2): 21},
    }
    if not all(transport_checks.values()):
        raise ArithmeticError(
            [key for key, value in transport_checks.items() if not value]
        )

    embedding = _primitive_k3_embedding(gram_standard)
    sources = {
        key: {
            "path": str(path.relative_to(ROOT)),
            "sha256": EXPECTED_SHA256[key],
        }
        for key, path in paths.items()
    }
    checks = {
        "all_source_hashes_pinned": all(
            _sha256(path) == EXPECTED_SHA256[key]
            for key, path in paths.items()
        ),
        "family_handoff_green": all(family_checks.values()),
        "wall_transport_green": all(transport_checks.values()),
        "primitive_K3_embedding_green": all(embedding["checks"].values()),
        "literal_Vinberg_green": all(volume_checks.values()),
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "Every very-general marked member of the explicit normalized "
            "nine-parameter quadratic-base-change K3 family, whose full "
            "Neron--Severi lattice is F=U+E8(-2)+<-2>, carries the known "
            "X2-RAT-C2-H10/37 finite-volume noncompact chamber as 37 "
            "primitive integral divisor mirrors spanning full NS."
        ),
        "sources": sources,
        "explicit_family": {
            "equation": family["family"]["k3_weierstrass_model"],
            "coefficient_A": family["family"]["A"],
            "coefficient_B": family["family"]["B"],
            "parameter_count": family["family"]["parameter_count"],
            "parameter_open": family["family"]["parameter_open"],
            "nonempty_open_witness": family["family"][
                "nonempty_open_witness"
            ],
            "fibers": family["family"]["fibers"]["X"],
            "very_general_full_Picard": family["generic_full_picard"][
                "claim"
            ],
            "geometric_basis": family[
                "fixed_locus_and_invariant_lattice"
            ]["geometric_integral_identification"],
            "checks": family_checks,
        },
        "full_Picard_lattice": {
            "name": "F=U+E8(-2)+<-2>",
            "rank": 11,
            "signature": [1, 10, 0],
            "determinant": 512,
            "smith_invariants": _smith(gram_standard),
            "nikulin_invariants": fixed["nikulin_invariants"],
            "standard_geometric_gram": _matrix_payload(gram_standard),
            "standard_to_X2_fixed_isometry": _matrix_payload(
                standard_to_fixed
            ),
            "primitive_K3_embedding": embedding,
            "K3_moduli_dimension": 9,
        },
        "transported_H10_37_chamber": {
            "id": "X2-RAT-C2-H10-37",
            "hyperbolic_dimension": 10,
            "facet_count": 37,
            "walls_span_full_Picard": True,
            "wall_span_rank": roots_by_rows.rank(),
            "wall_span_smith": span_smith,
            "wall_span_index": math.prod(span_smith),
            "primitive_mirror_type_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(types.items())
            },
            "transported_facet_gram_sha256": _payload_sha256(
                _matrix_payload(transported_gram)
            ),
            "square12_strict_interior_coordinates": _vector_payload(
                ample_standard
            ),
            "square12_strict_interior_pairings": pairings,
            "minimum_wall_pairing": min(pairings),
            "wall_records": wall_records,
            "checks": transport_checks,
        },
        "finite_volume": {
            "criterion": volume["criterion"],
            "classification": volume["classification"],
            "finite_vertices": volume["finite_vertex_count"],
            "cusps": volume["ideal_vertex_count"],
            "elliptic_rank9_edges": volume["elliptic_edge_count"],
            "bad_edges": volume["bad_edge_count"],
            "bad_parabolic_neighborhoods": volume[
                "bad_parabolic_neighborhood_count"
            ],
            "checks": volume_checks,
        },
        "shared_carrier_gain": {
            "same_explicit_family_already_carries": "A15_11-C2-H10-21",
            "newly_attached_chamber": "X2-RAT-C2-H10-37",
            "meaning": (
                "one explicit nine-dimensional very-general full-Picard "
                "complex-K3 family now carries two distinct certified H10 "
                "mixed divisor-reflection chambers, with 21 and 37 facets"
            ),
        },
        "scope": {
            "positive_gain": (
                "new project-level equation and full-Picard deformation "
                "realization of the already certified H10/37 chamber"
            ),
            "known_abstract_polytope": True,
            "not_new_abstract_polytope": True,
            "not_full_nef": (
                "the ten square-minus-four mirrors are divisor cuts; no "
                "claim identifies the 37-wall chamber with the full K3 nef cone"
            ),
            "not_surface_automorphisms": (
                "integral reflections of full NS are not asserted to extend "
                "to Hodge isometries or K3 surface automorphisms"
            ),
            "effectivity": (
                "after Weyl re-marking, Riemann--Roch orients the 27 "
                "square-minus-two mirrors effectively; irreducibility is not claimed"
            ),
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
        raise ArithmeticError("X2 H10/37 full-Picard carrier failed")
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(
            encoding="utf-8"
        ) != rendered:
            raise ArithmeticError("archived X2 H10/37 carrier changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={payload['payload_sha256']}")


if __name__ == "__main__":
    main()
