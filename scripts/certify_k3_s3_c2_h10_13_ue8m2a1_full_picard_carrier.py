#!/usr/bin/env python3
"""Certify the common equation-family carrier for K3-S3-C2-H10/13.

This composition certificate pins an explicit determinant-one isometry from
the saturated fixed lattice of the published rank-19 K3 transposition fold to

    F = U + E8(-2) + <-2>.

It transports all thirteen primitive facet normals to the geometric basis of
the normalized nine-parameter quadratic-base-change K3 family, verifies that
they span full Picard, and reruns the literal and optimized finite-covolume
criteria from the facet Gram.  The very-general carrier is a divisor-
reflection realization; the special S3 action and quotient-nef interpretation
are not asserted to deform, and the chamber is not asserted to be the full nef
cone of a very-general member.
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

from coxeter_polytopes.verify.vinberg_volume import (
    verify_vinberg_volume_criterion,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import (
    _enumerate_vertices,
    _literal_coxiter_prop_1_10_audit,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT / "data/enumerations/k3_rank18_19_all_subgroup_completion.json"
)
FAMILY = ROOT / "data/enumerations/a15_11_ue8m2a1_basechange_family.json"
OUTPUT = (
    ROOT
    / "data/enumerations/"
    "k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.json"
)

SCHEMA = "coxeter.k3_s3_c2.h10_13.ue8m2a1_full_picard_carrier.v1"
EXPECTED_SHA256 = {
    "source": "148ce106b968709c94b31ac8bed13b17269861ea6051eeba68b5bc8b89e7e5c4",
    "family": "1a5902092e3150bd380af19916ec5a6733392404a6cf2ecdba523f82b488532c",
}


# Columns are the geometric standard F basis expressed in the archived fixed
# basis.  The first two columns arise from the displayed isotropic U pair;
# the remaining nine come from an exact positive-definite complement isometry.
STANDARD_TO_FIXED = sp.Matrix(
    [
        [0, 0, 1, 0, 2, 0, -2, 1, 0, 0, 1],
        [0, -1, 2, 0, 4, 0, -4, 2, 0, 0, 0],
        [0, 0, -1, 0, -2, 0, 2, -1, 0, 0, 0],
        [-1, -1, 1, 0, 3, 0, -3, 2, 0, -1, -1],
        [-1, -1, 1, 0, 2, 0, -2, 1, 1, -2, -1],
        [-2, -2, 0, 0, 1, 0, -1, 2, 0, -4, -2],
        [-2, -2, -1, 0, 0, 0, 0, 2, 0, -4, -2],
        [-3, -3, -3, 1, -1, 0, 0, 3, 0, -6, -3],
        [-2, -2, -2, 0, -2, 1, 0, 2, 0, -4, -2],
        [-1, -1, -1, 0, -2, 0, 1, 1, 0, -2, -1],
        [1, 1, 1, 0, 0, 0, 0, -1, 0, 2, 1],
    ]
)

ISOTROPIC_U_VECTOR = sp.Matrix([0, 0, 0, -1, -1, -2, -2, -3, -2, -1, 1])
ISOTROPIC_U_PARTNER = sp.Matrix(
    [0, -1, 0, -1, -1, -2, -2, -3, -2, -1, 1]
)

# A deliberately small strict chamber-interior class in the archived fixed
# basis.  Its thirteen pairings and positive square are recomputed below.
FIXED_INTERIOR = sp.Matrix([9, 20, 10, 19, 14, 15, 12, 16, 13, 11, -5])


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
    """Embed F primitively into Lambda_K3 and display its complement."""

    U = sp.Matrix([[0, 1], [1, 0]])
    E8m = sp.simplify(target[2:10, 2:10] / 2)
    ambient = sp.diag(U, U, U, E8m, E8m)

    embedding = sp.zeros(22, 11)
    embedding[0, 0] = 1
    embedding[1, 1] = 1
    for index in range(8):
        embedding[6 + index, 2 + index] = 1
        embedding[14 + index, 2 + index] = 1
    embedding[2, 10] = 1
    embedding[3, 10] = -1

    complement = sp.zeros(22, 11)
    complement[4, 0] = 1
    complement[5, 1] = 1
    for index in range(8):
        complement[6 + index, 2 + index] = 1
        complement[14 + index, 2 + index] = -1
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
    }
    if not all(checks.values()):
        raise ArithmeticError([key for key, value in checks.items() if not value])
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


def _fresh_volume_replay(
    facet_gram: sp.Matrix, source_geometry: dict[str, Any]
) -> dict[str, Any]:
    """Rerun both exact finite-volume lanes from the 13-wall Gram."""

    dimension = 10
    facets = facet_gram.rows
    vinberg_gram = [
        [-int(facet_gram[left, right]) for right in range(facets)]
        for left in range(facets)
    ]
    literal = _literal_coxiter_prop_1_10_audit(vinberg_gram, dimension)
    optimized = verify_vinberg_volume_criterion(
        vinberg_gram, dimension, simple_roots_complete=True
    )
    vertices, _, _ = _enumerate_vertices(
        facet_gram, list(range(facets)), []
    )

    finite = [vertex for vertex in vertices if vertex["type"] == "finite"]
    ideal = [vertex for vertex in vertices if vertex["type"] == "ideal"]
    hyperideal = [
        vertex for vertex in vertices if vertex["type"] == "hyperideal"
    ]
    literal_finite = set(literal["spherical_vertices"])
    literal_ideal = set(literal["parabolic_vertices"])
    optimized_finite = set(optimized.finite_vertices)
    optimized_ideal = {endpoint.facets for endpoint in optimized.ideal_vertices}
    ray_finite = {
        tuple(index - 1 for index in vertex["active_facets"])
        for vertex in finite
    }
    ray_ideal = {
        tuple(index - 1 for index in vertex["active_facets"])
        for vertex in ideal
    }

    endpoint_histogram = {
        f"{finite_count},{ideal_count}": multiplicity
        for (finite_count, ideal_count), multiplicity in sorted(
            literal["endpoint_type_histogram"].items()
        )
    }
    literal_parabolics_one_based = sorted(
        [index + 1 for index in vertex]
        for vertex in literal["parabolic_vertices"]
    )
    archived_literal = source_geometry["literal_coxiter_prop_1_10"]
    archived_optimized = source_geometry["optimized_vinberg"]
    archived_rays = source_geometry["exact_cone"]

    ideal_link_checks = []
    for vertex in ideal:
        active = [index - 1 for index in vertex["active_facets"]]
        active_gram = facet_gram.extract(active, active)
        inertia = exact_inertia(active_gram)
        ideal_link_checks.append(
            {
                "active_facets": vertex["active_facets"],
                "inertia": list(inertia),
                "rank": active_gram.rank(),
                "negative_semidefinite_rank9": (
                    inertia[0] == 0
                    and inertia[1] == 9
                    and active_gram.rank() == 9
                ),
            }
        )

    checks = {
        "fresh_extreme_ray_census_34_4_0": (
            len(finite), len(ideal), len(hyperideal)
        )
        == (34, 4, 0),
        "fresh_literal_208_edges": len(literal["spherical_edges"]) == 208,
        "fresh_literal_34_finite_vertices": len(literal_finite) == 34,
        "fresh_literal_4_parabolic_vertices": len(literal_ideal) == 4,
        "fresh_literal_every_edge_has_two_endpoints": (
            literal["extension_histogram"] == {2: 208}
            and not literal["bad_edges"]
        ),
        "fresh_literal_endpoint_histogram": endpoint_histogram
        == {"0,2": 5, "1,1": 66, "2,0": 137},
        "fresh_literal_finite_covolume": literal["finite_covolume"] is True,
        "fresh_optimized_finite_volume_noncompact": (
            optimized.decision == "finite_volume_noncompact"
            and optimized.finite_volume is True
            and optimized.compact is False
        ),
        "fresh_optimized_counts": (
            optimized.finite_vertex_count,
            optimized.ideal_vertex_count,
            optimized.elliptic_edge_count,
        )
        == (34, 4, 208),
        "fresh_optimized_zero_bad_edges": not optimized.finite_volume_bad_edges,
        "fresh_all_affine_neighborhood_audits_good": (
            len(optimized.parabolic_neighborhood_audits) == 7
            and all(
                audit.is_rank_dimension_minus_one_parabolic
                for audit in optimized.parabolic_neighborhood_audits
            )
        ),
        "literal_optimized_and_extreme_endpoint_sets_agree": (
            literal_finite == optimized_finite == ray_finite
            and literal_ideal == optimized_ideal == ray_ideal
        ),
        "all_four_ideal_links_are_rank9_parabolic": (
            len(ideal_link_checks) == 4
            and all(
                record["negative_semidefinite_rank9"]
                for record in ideal_link_checks
            )
        ),
        "fresh_literal_matches_pinned_archive": (
            archived_literal["spherical_rank_n_minus_one_count"] == 208
            and archived_literal["spherical_rank_n_count"] == 34
            and archived_literal["parabolic_rank_n_minus_one_count"] == 4
            and archived_literal["extension_histogram"] == {"2": 208}
            and archived_literal["endpoint_type_histogram"]
            == endpoint_histogram
            and sorted(archived_literal["parabolic_vertices"])
            == literal_parabolics_one_based
            and archived_literal["bad_edges"] == []
            and archived_literal["finite_covolume"] is True
        ),
        "fresh_optimized_matches_pinned_archive": (
            archived_optimized["decision"] == optimized.decision
            and archived_optimized["finite_vertex_count"] == 34
            and archived_optimized["ideal_vertex_count"] == 4
            and archived_optimized["elliptic_edge_count"] == 208
            and archived_optimized["finite_volume_bad_edges"] == []
            and archived_optimized["parabolic_neighborhood_audit_count"] == 7
        ),
        "fresh_rays_match_pinned_archive_counts": (
            archived_rays["finite_vertex_count"] == len(finite)
            and archived_rays["ideal_vertex_count"] == len(ideal)
            and archived_rays["hyperideal_vertex_count"] == len(hyperideal)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError([key for key, value in checks.items() if not value])

    return {
        "criterion": (
            "fresh literal CoxIter Proposition 1.10 all-subset extension "
            "enumeration, independent optimized Vinberg full-neighborhood "
            "audit, and exact projective extreme-ray enumeration"
        ),
        "classification": optimized.decision,
        "finite_vertices": len(finite),
        "cusps": len(ideal),
        "hyperideal_vertices": len(hyperideal),
        "elliptic_rank9_edges": len(literal["spherical_edges"]),
        "endpoint_type_histogram": endpoint_histogram,
        "literal_bad_edges": len(literal["bad_edges"]),
        "optimized_bad_edges": len(optimized.finite_volume_bad_edges),
        "parabolic_neighborhood_audits": len(
            optimized.parabolic_neighborhood_audits
        ),
        "parabolic_vertices_one_based": literal_parabolics_one_based,
        "ideal_link_checks": ideal_link_checks,
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    paths = {"source": SOURCE, "family": FAMILY}
    for key, path in paths.items():
        actual = _sha256(path)
        if actual != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed: {actual}")

    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    family = json.loads(FAMILY.read_text(encoding="utf-8"))
    if not source["all_checks_pass"] or not family["all_checks_pass"]:
        raise ArithmeticError("an upstream theorem archive is not green")

    geometry = source["survivor_geometry"]
    if geometry["id"] != "K3-S3-C2-H10-13":
        raise ArithmeticError("the pinned survivor changed")
    fixed = geometry["fixed_lattice"]
    gram_fixed = sp.Matrix(fixed["gram"])
    gram_standard = sp.Matrix(
        family["fixed_locus_and_invariant_lattice"]["gram"]
    )
    fixed_to_standard = STANDARD_TO_FIXED.inv()

    bridge_checks = {
        "fixed_signature_1_10": exact_inertia(gram_fixed) == (1, 10, 0),
        "standard_signature_1_10": exact_inertia(gram_standard) == (1, 10, 0),
        "both_even_determinant512": (
            gram_fixed.det() == gram_standard.det() == 512
            and all(gram_fixed[index, index] % 2 == 0 for index in range(11))
            and all(gram_standard[index, index] % 2 == 0 for index in range(11))
        ),
        "standard_to_fixed_integral_determinant_one": (
            all(value.is_Integer for value in STANDARD_TO_FIXED)
            and STANDARD_TO_FIXED.det() == 1
        ),
        "explicit_Gram_isometry": (
            STANDARD_TO_FIXED.T * gram_fixed * STANDARD_TO_FIXED
            == gram_standard
        ),
        "inverse_integral_determinant_one": (
            all(value.is_Integer for value in fixed_to_standard)
            and fixed_to_standard.det() == 1
        ),
        "displayed_U_pair_is_isotropic_unimodular": (
            (ISOTROPIC_U_VECTOR.T * gram_fixed * ISOTROPIC_U_VECTOR)[0]
            == 0
            and (
                ISOTROPIC_U_PARTNER.T
                * gram_fixed
                * ISOTROPIC_U_PARTNER
            )[0]
            == 0
            and (
                ISOTROPIC_U_VECTOR.T
                * gram_fixed
                * ISOTROPIC_U_PARTNER
            )[0]
            == 1
            and STANDARD_TO_FIXED[:, 0] == ISOTROPIC_U_VECTOR
            and STANDARD_TO_FIXED[:, 1] == ISOTROPIC_U_PARTNER
        ),
        "Smith_invariants_1_1_2x9": (
            _smith(gram_fixed) == _smith(gram_standard) == [1, 1] + [2] * 9
        ),
    }
    if not all(bridge_checks.values()):
        raise ArithmeticError(
            [key for key, value in bridge_checks.items() if not value]
        )

    interior_standard = fixed_to_standard * FIXED_INTERIOR
    roots_standard: list[sp.Matrix] = []
    wall_records: list[dict[str, Any]] = []
    type_histogram: Counter[tuple[int, int]] = Counter()
    pairings: list[int] = []
    for record in fixed["mirror_records"]:
        if record["role"] != "facet":
            continue
        root_fixed = sp.Matrix(record["fixed_basis_coordinates"])
        root_standard = sp.simplify(fixed_to_standard * root_fixed)
        square = int((root_standard.T * gram_standard * root_standard)[0])
        divisibility = _divisibility(gram_standard, root_standard)
        reflection = _reflection(gram_standard, root_standard)
        pairing = int(
            (interior_standard.T * gram_standard * root_standard)[0]
        )
        checks = {
            "transport_integral": all(value.is_Integer for value in root_standard),
            "primitive": _content(root_standard) == 1,
            "square_preserved": square == record["square"],
            "strict_interior_orientation": pairing > 0,
            "reflection_integral": all(value.is_Integer for value in reflection),
            "reflection_unimodular": reflection.det() == -1,
            "reflection_involutive": reflection**2 == sp.eye(11),
            "reflection_preserves_Gram": (
                reflection.T * gram_standard * reflection == gram_standard
            ),
            "reflection_negates_root": reflection * root_standard == -root_standard,
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"wall {record['orbit_index']} failed: "
                f"{[key for key, value in checks.items() if not value]}"
            )
        roots_standard.append(root_standard)
        pairings.append(pairing)
        type_histogram[(square, divisibility)] += 1
        wall_records.append(
            {
                "facet_one_based": record["orbit_index"],
                "source_walls": record["source_walls"],
                "primitive_fixed_basis_coordinates": record[
                    "fixed_basis_coordinates"
                ],
                "primitive_geometric_standard_F_coordinates": _vector_payload(
                    root_standard
                ),
                "square": square,
                "divisibility_in_full_Picard": divisibility,
                "square358_interior_pairing": pairing,
                "checks": checks,
            }
        )

    roots_by_rows = sp.Matrix.hstack(*roots_standard).T
    span_smith = _smith(roots_by_rows)
    transported_facet_gram = sp.Matrix(
        [
            [
                int((left.T * gram_standard * right)[0])
                for right in roots_standard
            ]
            for left in roots_standard
        ]
    )
    expected_facet_gram = sp.Matrix(geometry["facet_gram"])
    interior_square = int(
        (interior_standard.T * gram_standard * interior_standard)[0]
    )
    transport_checks = {
        "all_13_roots_transport": len(roots_standard) == 13,
        "all_wall_checks_green": all(
            all(record["checks"].values()) for record in wall_records
        ),
        "transported_facet_Gram_literal": (
            transported_facet_gram == expected_facet_gram
        ),
        "all_13_roots_span_full_Picard": (
            roots_by_rows.rank() == 11 and span_smith == [1] * 11
        ),
        "mirror_type_histogram": dict(type_histogram)
        == {(-4, 2): 9, (-2, 1): 2, (-2, 2): 2},
        "square358_strict_interior": (
            interior_square == 358
            and pairings == [2, 2, 2, 2, 2, 6, 2, 2, 2, 2, 2, 2, 38]
        ),
    }
    if not all(transport_checks.values()):
        raise ArithmeticError(
            [key for key, value in transport_checks.items() if not value]
        )

    volume = _fresh_volume_replay(expected_facet_gram, geometry)
    family_checks = {
        "explicit_family_green": family["all_checks_pass"] is True,
        "nine_parameters": family["family"]["parameter_count"] == 9,
        "nonempty_open_witness": bool(
            family["family"]["nonempty_open_witness"]
        ),
        "very_general_full_NS_literal": (
            family["generic_full_picard"]["claim"]
            == "NS(X)=U+E8(-2)+<-2> for a very general parameter"
        ),
        "full_moduli_dimension_nine": (
            family["generic_full_picard"]["parameter_moduli_dimension"] == 9
            and family["generic_full_picard"][
                "two_elementary_pair_moduli_dimension"
            ]
            == 9
        ),
        "same_standard_geometric_Gram": (
            gram_standard
            == sp.Matrix(family["fixed_locus_and_invariant_lattice"]["gram"])
        ),
        "same_Nikulin_invariants": (
            family["fixed_locus_and_invariant_lattice"]["nikulin_invariants"]
            == {"r": 11, "a": 9, "delta": 1}
        ),
    }
    if not all(family_checks.values()):
        raise ArithmeticError(
            [key for key, value in family_checks.items() if not value]
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
        "determinant_one_bridge_green": all(bridge_checks.values()),
        "wall_transport_green": all(transport_checks.values()),
        "fresh_finite_volume_replay_green": all(volume["checks"].values()),
        "explicit_family_genericization_green": all(family_checks.values()),
        "primitive_K3_embedding_green": all(embedding["checks"].values()),
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "Every very-general marked member of the explicit normalized "
            "nine-parameter quadratic-base-change K3 family, whose full "
            "Neron--Severi lattice is F=U+E8(-2)+<-2>, carries the known "
            "K3-S3-C2-H10/13 finite-volume noncompact chamber as thirteen "
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
            "very_general_full_Picard": family["generic_full_picard"]["claim"],
            "Noether_Lefschetz_argument": family["generic_full_picard"][
                "noether_lefschetz_step"
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
            "nikulin_invariants": {"r": 11, "a": 9, "delta": 1},
            "standard_geometric_gram": _matrix_payload(gram_standard),
            "standard_to_archived_fixed_GL11Z": _matrix_payload(
                STANDARD_TO_FIXED
            ),
            "archived_fixed_to_standard_GL11Z": _matrix_payload(
                fixed_to_standard
            ),
            "isometry_determinant": int(STANDARD_TO_FIXED.det()),
            "isometry_construction": {
                "isotropic_U_vector_in_fixed_basis": _vector_payload(
                    ISOTROPIC_U_VECTOR
                ),
                "isotropic_U_partner_in_fixed_basis": _vector_payload(
                    ISOTROPIC_U_PARTNER
                ),
                "definite_complement": "E8(-2)+<-2>",
                "method": (
                    "split the displayed integral U, then apply an exact "
                    "rank-nine positive-definite qfisom certificate"
                ),
            },
            "primitive_K3_embedding": embedding,
            "K3_moduli_dimension": 9,
            "checks": bridge_checks,
        },
        "transported_H10_13_chamber": {
            "id": "K3-S3-C2-H10-13",
            "hyperbolic_dimension": 10,
            "facet_count": 13,
            "known_full_reflection_chamber": True,
            "walls_span_full_Picard": True,
            "wall_span_rank": roots_by_rows.rank(),
            "wall_span_smith": span_smith,
            "wall_span_index": math.prod(span_smith),
            "primitive_mirror_type_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(type_histogram.items())
            },
            "transported_facet_gram_sha256": _payload_sha256(
                _matrix_payload(transported_facet_gram)
            ),
            "square358_strict_interior_fixed_coordinates": _vector_payload(
                FIXED_INTERIOR
            ),
            "square358_strict_interior_standard_F_coordinates": _vector_payload(
                interior_standard
            ),
            "square358_strict_interior_pairings": pairings,
            "minimum_wall_pairing": min(pairings),
            "wall_records": wall_records,
            "checks": transport_checks,
        },
        "finite_volume": volume,
        "shared_carrier_gain": {
            "one_explicit_family_carries": [
                "K3-S3-C2-H10-13",
                "A15_11-C2-H10-21",
                "X2-RAT-C2-H10-37",
            ],
            "facet_counts": [13, 21, 37],
            "meaning": (
                "one explicit nine-dimensional very-general full-Picard "
                "complex-K3 family now carries three distinct certified H10 "
                "mixed divisor-reflection chambers"
            ),
        },
        "scope": {
            "positive_gain": (
                "new project-level equation and very-general full-Picard "
                "deformation realization of the already known H10/13 chamber"
            ),
            "known_abstract_polytope": True,
            "not_new_abstract_polytope": True,
            "special_surface_geometry_kept_separate": (
                "Kondo's S3 action and the quotient-nef interpretation live "
                "on the special rank-19 source K3 and are not asserted on a "
                "very-general base-change member"
            ),
            "not_full_nef": (
                "the nine square-minus-four mirrors are integral divisor "
                "reflections; no claim identifies this chamber with the full "
                "nef cone of the very-general K3"
            ),
            "not_surface_automorphisms": (
                "the thirteen integral Picard reflections are not asserted "
                "to extend to Hodge isometries or K3 surface automorphisms"
            ),
            "effectivity": (
                "after Weyl re-marking, Riemann--Roch orients the four "
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
        raise ArithmeticError("K3-S3-C2 H10/13 full-Picard carrier failed")
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(
            encoding="utf-8"
        ) != rendered:
            raise ArithmeticError("archived K3-S3-C2 H10/13 carrier changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={payload['payload_sha256']}")


if __name__ == "__main__":
    main()
