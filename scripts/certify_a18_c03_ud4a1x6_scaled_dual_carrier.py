#!/usr/bin/env python3
"""Certify a full-Picard K3 scaled-dual carrier of A18 C03 H11/28.

The A18 C03 fixed lattice has a unique reflection-stable odd index-four
closure M of determinant 16.  Its doubled dual N=M^vee(2) is even,
2-elementary, and has invariants (r,a,delta)=(12,8,1).  An explicit
discriminant graph glue with T=U+U(2)+A1(-1)^6 produces the even unimodular
K3 lattice and proves that N embeds primitively.  Period surjectivity then
gives projective complex K3 surfaces with full Picard lattice N.  All 28
oriented halfspaces and reflections transport to N, their primitive mirrors
generate N, and the literal upstream Vinberg audit is replayed.
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
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _q_histogram,
)
from scripts.certify_rank12_five_chamber_simultaneous_carrier import (
    _binary_discriminant_isometry,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = (
    ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json"
)
STABLE = (
    ROOT
    / "data/enumerations/alexeev_remaining_high_fixed_section_stable_overlattices.json"
)
VINBERG = ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_vinberg.json"
OUTPUT = ROOT / "data/enumerations/a18_c03_ud4a1x6_scaled_dual_carrier.json"

SCHEMA = "coxeter.a18_c03.ud4a1x6_scaled_dual_k3_carrier.v1"
EXPECTED_SHA256 = {
    "construction": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
    "stable": "882204e5457e96a3930399201dd061bd6bc35839497693a384d0b931f10924a0",
    "vinberg": "d5e56b83f66182a0c52758ee2be8d09217c66f0e4ffe4bfc82bc73b867f2ab0c",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(row)
        ]
        for row in range(matrix.rows)
    ]


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _smith(matrix: sp.Matrix) -> list[int]:
    normal = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(normal[index, index]))
        for index in range(min(normal.rows, normal.cols))
    ]


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector))


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * gram * root)[0]
    return sp.simplify(
        sp.eye(gram.rows) - sp.Rational(2, square) * root * root.T * gram
    )


def _histogram_payload(counter: Counter[tuple[int, ...]]) -> dict[str, int]:
    return {
        ",".join(str(value) for value in key): count
        for key, count in sorted(counter.items())
    }


def _standard_picard_gram() -> sp.Matrix:
    """The standard representative U + D4(-1) + A1(-1)^6."""

    hyperbolic_plane = sp.Matrix([[0, 1], [1, 0]])
    d4_negative = -sp.Matrix(
        [[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]]
    )
    return sp.diag(hyperbolic_plane, d4_negative, sp.diag(*([-2] * 6)))


def _complement_gram() -> sp.Matrix:
    """A signature-(2,8) anti-discriminant complement."""

    hyperbolic_plane = sp.Matrix([[0, 1], [1, 0]])
    doubled_plane = sp.Matrix([[0, 2], [2, 0]])
    return sp.diag(hyperbolic_plane, doubled_plane, sp.diag(*([-2] * 6)))


def _explicit_k3_glue(gram_n: sp.Matrix) -> dict[str, Any]:
    """Glue N and T along a displayed anti-isometry of discriminant forms."""

    complement = _complement_gram()
    anti_isometry = _binary_discriminant_isometry(gram_n, -complement)
    if anti_isometry is None:
        raise ArithmeticError("N and the proposed complement are not anti-isometric")

    moduli_n, left_inverse_n, _ = _discriminant_data(gram_n)
    moduli_t, left_inverse_minus_t, _ = _discriminant_data(-complement)
    indices_n = [index for index, modulus in enumerate(moduli_n) if modulus > 1]
    indices_t = [index for index, modulus in enumerate(moduli_t) if modulus > 1]
    if not (
        len(indices_n) == len(indices_t) == 8
        and all(moduli_n[index] == 2 for index in indices_n)
        and all(moduli_t[index] == 2 for index in indices_t)
    ):
        raise ArithmeticError("the expected elementary discriminant bases changed")

    glue_vectors: list[sp.Matrix] = []
    generators = [sp.eye(22).col(index) for index in range(22)]
    for source_index, image in zip(
        indices_n, anti_isometry["basis_images_as_columns"]
    ):
        source_snf_vector = sp.zeros(12, 1)
        source_snf_vector[source_index] = 1
        source_dual = gram_n.inv() * left_inverse_n * source_snf_vector

        target_snf_vector = sp.zeros(10, 1)
        for bit, target_index in zip(image, indices_t):
            target_snf_vector[target_index] = bit
        target_dual = (
            (-complement).inv() * left_inverse_minus_t * target_snf_vector
        )
        glue_vector = source_dual.col_join(target_dual)
        glue_vectors.append(glue_vector)
        generators.append(glue_vector)

    denominator = math.lcm(
        *(int(value.q) for vector in generators for value in vector)
    )
    integer_generators = sp.Matrix.hstack(
        *(denominator * vector for vector in generators)
    ).applyfunc(int)
    glue_basis = hermite_normal_form(integer_generators) / denominator
    split_gram = sp.diag(gram_n, complement)
    k3_gram = sp.simplify(glue_basis.T * split_gram * glue_basis)

    n_embedding = sp.simplify(glue_basis.inv() * sp.eye(22)[:, :12])
    t_embedding = sp.simplify(glue_basis.inv() * sp.eye(22)[:, 12:])
    checks = {
        "anti_isometry_all_256_classes": (
            anti_isometry["all_elements_checked"] == 256
            and anti_isometry["quadratic_form_preserved_mod_2Z"]
            and anti_isometry["bilinear_form_preserved_mod_Z"]
            and anti_isometry["bijective"]
        ),
        "eight_order_two_glue_generators": (
            denominator == 2 and len(glue_vectors) == 8
        ),
        "glue_index_256": abs(glue_basis.det()) == sp.Rational(1, 256),
        "even_unimodular_signature_3_19": (
            all(value.is_Integer for value in k3_gram)
            and all(int(k3_gram[index, index]) % 2 == 0 for index in range(22))
            and k3_gram.det() == -1
            and exact_inertia(k3_gram) == (3, 19, 0)
        ),
        "N_embedding_integral_and_isometric": (
            all(value.is_Integer for value in n_embedding)
            and n_embedding.T * k3_gram * n_embedding == gram_n
        ),
        "T_embedding_integral_and_isometric": (
            all(value.is_Integer for value in t_embedding)
            and t_embedding.T * k3_gram * t_embedding == complement
        ),
        "N_and_T_orthogonal": (
            n_embedding.T * k3_gram * t_embedding == sp.zeros(12, 10)
        ),
        "both_embeddings_primitive": (
            _smith(n_embedding) == [1] * 12
            and _smith(t_embedding) == [1] * 10
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(
            f"explicit K3 glue failed: "
            f"{[key for key, value in checks.items() if not value]}"
        )
    return {
        "complement": {
            "lattice": "T=U+U(2)+A1(-1)^6",
            "gram": _payload(complement),
            "signature": [2, 8, 0],
            "determinant": int(complement.det()),
            "smith_invariants": _smith(complement),
            "discriminant_q_histogram": _q_histogram(complement),
        },
        "discriminant_anti_isometry": anti_isometry,
        "glue_vectors_in_N_plus_T_coordinates": [
            [
                int(value)
                if sp.Rational(value).q == 1
                else str(sp.Rational(value))
                for value in vector
            ]
            for vector in glue_vectors
        ],
        "glued_lattice_basis_columns_in_N_plus_T": _payload(glue_basis),
        "glued_even_unimodular_gram": _payload(k3_gram),
        "N_embedding_columns_in_glued_K3_basis": _payload(n_embedding),
        "T_embedding_columns_in_glued_K3_basis": _payload(t_embedding),
        "classification": (
            "the unique even unimodular lattice of signature (3,19), "
            "hence the K3 lattice U^3+E8(-1)^2"
        ),
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "construction": CONSTRUCTION,
        "stable": STABLE,
        "vinberg": VINBERG,
    }
    for key, path in paths.items():
        if _sha256(path) != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed")

    construction = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    stable = json.loads(STABLE.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    case = next(row for row in construction["cases"] if row["id"] == "A18_4_0_C03")
    stable_case = next(row for row in stable["cases"] if row["id"] == "A18_4_0_C03")
    vinberg_case = next(row for row in vinberg["cases"] if row["id"] == "A18_4_0_C03")

    gram_l = sp.Matrix(case["saturated_fixed_lattice"]["gram"])
    closure = stable_case["unique_maximal_reflection_stable_overlattice"]
    closure_basis = sp.Matrix(closure["basis_columns_in_original_fixed_lattice"])
    gram_m = sp.simplify(closure_basis.T * gram_l * closure_basis)
    gram_n = sp.simplify(2 * gram_m.inv())
    standard_picard = _standard_picard_gram()

    standard_discriminant_isometry = _binary_discriminant_isometry(
        gram_n, standard_picard
    )
    if standard_discriminant_isometry is None:
        raise ArithmeticError("N missed its proposed standard 2-elementary genus")
    k3_glue = _explicit_k3_glue(gram_n)

    lattice_checks = {
        "all_three_upstream_certificates_green": (
            construction["all_checks_pass"]
            and stable["all_checks_pass"]
            and vinberg["all_checks_pass"]
            and all(case["checks"].values())
        ),
        "source_fixed_lattice_rank12_det256_even": (
            exact_inertia(gram_l) == (1, 11, 0)
            and gram_l.det() == -256
            and _smith(gram_l) == [1] * 6 + [2] * 4 + [4] * 2
            and all(int(gram_l[index, index]) % 2 == 0 for index in range(12))
        ),
        "unique_stable_odd_index_four_closure": (
            closure["unique"] is True
            and closure["index_over_original_lattice"] == 4
            and abs(closure_basis.det()) == sp.Rational(1, 4)
            and gram_m == sp.Matrix(closure["gram"])
            and gram_m.det() == -16
            and exact_inertia(gram_m) == (1, 11, 0)
            and _smith(gram_m) == [1] * 8 + [2] * 4
            and any(int(gram_m[index, index]) % 2 for index in range(12))
        ),
        "doubled_dual_even_2_elementary_12_8_1": (
            all(value.is_Integer for value in gram_n)
            and gram_n.det() == -256
            and exact_inertia(gram_n) == (1, 11, 0)
            and _smith(gram_n) == [1] * 4 + [2] * 8
            and all(int(gram_n[index, index]) % 2 == 0 for index in range(12))
            and _q_histogram(gram_n)
            == {"0": 64, "1": 64, "1/2": 56, "3/2": 72}
        ),
        "exact_standard_discriminant_form": (
            standard_picard.det() == -256
            and exact_inertia(standard_picard) == (1, 11, 0)
            and _smith(standard_picard) == [1] * 4 + [2] * 8
            and _q_histogram(standard_picard) == _q_histogram(gram_n)
            and standard_discriminant_isometry["all_elements_checked"] == 256
            and standard_discriminant_isometry["bijective"]
        ),
        "explicit_primitive_K3_glue": all(k3_glue["checks"].values()),
        "nikulin_sufficient_codimension_gate": (22 - 12 == 8 + 2),
    }
    if not all(lattice_checks.values()):
        raise ArithmeticError(
            f"A18 C03 K3 lattice gate failed: "
            f"{[key for key, value in lattice_checks.items() if not value]}"
        )

    active_orbits = set(case["exact_reduction"]["active_facet_orbit_indices"])
    active_records = [
        row
        for row in case["restricted_mirrors"]["records"]
        if row["orbit_index"] in active_orbits
    ]
    if len(active_records) != 28:
        raise ArithmeticError("active C03 facet list changed")

    inverse_closure_basis = closure_basis.inv()
    records: list[dict[str, Any]] = []
    roots_m: list[sp.Matrix] = []
    roots_n: list[sp.Matrix] = []
    closure_types: Counter[tuple[int, int]] = Counter()
    target_types: Counter[tuple[int, int]] = Counter()
    source_to_target_types: Counter[tuple[int, ...]] = Counter()

    for index, source in enumerate(active_records):
        root_l = sp.Matrix(source["primitive_integral_same_direction_coordinates"])
        orientation = int(sp.Rational(source["presentation_scale_from_primitive"]))
        if orientation not in (-1, 1):
            raise ArithmeticError("unexpected non-unit C03 presentation scale")
        raw_root_m = inverse_closure_basis * root_l
        if any(not value.is_Integer for value in raw_root_m):
            raise ArithmeticError(f"facet {index} is absent from the stable closure")
        closure_content = _content(raw_root_m)
        root_m = orientation * (raw_root_m / closure_content).applyfunc(int)
        div_m = _content(gram_m * root_m)
        root_n = (gram_m * root_m / div_m).applyfunc(int)

        square_l = int((root_l.T * gram_l * root_l)[0])
        div_l = _content(gram_l * root_l)
        square_m = int((root_m.T * gram_m * root_m)[0])
        square_n = int((root_n.T * gram_n * root_n)[0])
        div_n = _content(gram_n * root_n)
        reflection_m = _reflection(gram_m, root_m)
        reflection_n = _reflection(gram_n, root_n)
        factor_from_m = sp.Rational(2, div_m)
        factor_from_l = sp.Rational(2, div_m * closure_content)

        checks = {
            "source_square_archived": square_l == int(source["square"]),
            "source_pairing_column_and_orientation": (
                gram_l * (orientation * root_l)
                == sp.Matrix(source["fixed_pairing_column"])
            ),
            "primitive_in_M_and_N": _content(root_m) == _content(root_n) == 1,
            "scaled_dual_square_law": square_n * div_m**2 == 2 * square_m,
            "exact_positive_halfspace_transport": (
                factor_from_m > 0
                and factor_from_l > 0
                and root_n.T * gram_n * gram_m
                == factor_from_m * root_m.T * gram_m
            ),
            "closure_reflection_integral": all(
                value.is_Integer for value in reflection_m
            ),
            "contragredient_reflection": reflection_n == reflection_m.T,
            "target_reflection_integral": all(
                value.is_Integer for value in reflection_n
            ),
            "target_reflection_unimodular_involution_isometry": (
                reflection_n.det() == -1
                and reflection_n**2 == sp.eye(12)
                and reflection_n.T * gram_n * reflection_n == gram_n
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"facet {index} failed: "
                f"{[key for key, value in checks.items() if not value]}"
            )

        closure_types[(square_m, div_m)] += 1
        target_types[(square_n, div_n)] += 1
        source_to_target_types[
            (square_l, div_l, closure_content, square_m, div_m, square_n, div_n)
        ] += 1
        roots_m.append(root_m)
        roots_n.append(root_n)
        records.append(
            {
                "facet_zero_based": index,
                "source_orbit_index_one_based": source["orbit_index"],
                "source_L_coordinates_unoriented": [int(value) for value in root_l],
                "source_presentation_orientation": orientation,
                "source_L_square": square_l,
                "source_L_divisibility": div_l,
                "coordinates_in_M_before_primitive_reduction_unoriented": [
                    int(value) for value in raw_root_m
                ],
                "M_coordinate_content": closure_content,
                "primitive_oriented_M_coordinates": [int(value) for value in root_m],
                "primitive_M_square": square_m,
                "primitive_M_divisibility": div_m,
                "primitive_oriented_N_coordinates": [int(value) for value in root_n],
                "primitive_N_square": square_n,
                "primitive_N_divisibility": div_n,
                "positive_halfspace_factor_from_M": str(factor_from_m),
                "positive_halfspace_factor_from_original_L": str(factor_from_l),
                "N_reflection_matrix": _payload(reflection_n),
                "checks": checks,
            }
        )

    pairwise_transport = all(
        (roots_n[left].T * gram_n * roots_n[right])[0]
        * _content(gram_m * roots_m[left])
        * _content(gram_m * roots_m[right])
        == 2 * (roots_m[left].T * gram_m * roots_m[right])[0]
        for left in range(28)
        for right in range(28)
    )
    roots_n_matrix = sp.Matrix.hstack(*roots_n).applyfunc(int)
    hnf_n = hermite_normal_form(roots_n_matrix)
    carrier = sum(roots_n, sp.zeros(12, 1))
    carrier_pairings = [
        int((carrier.T * gram_n * root)[0]) for root in roots_n
    ]
    generation_checks = {
        "all_28_active_facets_transport": len(records) == 28,
        "exact_pairwise_positive_Gram_rescaling": pairwise_transport,
        "primitive_N_roots_generate_full_Picard_lattice": hnf_n == sp.eye(12),
        "sum_is_strict_timelike_interior": (
            int((carrier.T * gram_n * carrier)[0]) == 608
            and min(carrier_pairings) == 2
            and all(value > 0 for value in carrier_pairings)
        ),
        "closure_type_histogram": closure_types
        == Counter({(-1, 1): 12, (-2, 1): 8, (-4, 2): 8}),
        "target_type_histogram": target_types
        == Counter({(-2, 2): 12, (-4, 2): 8, (-2, 1): 8}),
    }
    if not all(generation_checks.values()):
        raise ArithmeticError(
            f"generation/interior audit failed: "
            f"{[key for key, value in generation_checks.items() if not value]}"
        )

    vinberg_checks = {
        "literal_case_green": vinberg_case["all_checks_pass"] is True,
        "finite_volume_noncompact": (
            vinberg_case["finite_volume"] is True
            and vinberg_case["compact"] is False
            and vinberg_case["classification"]
            == "finite_volume_noncompact_Coxeter_polytope"
        ),
        "every_elliptic_edge_has_two_endpoints": (
            vinberg_case["bad_edge_count"] == 0
            and vinberg_case["elliptic_edge_count"] == 8851
        ),
        "all_parabolic_neighborhoods_good": (
            vinberg_case["bad_parabolic_neighborhood_count"] == 0
            and vinberg_case["connected_affine_component_count"] == 186
        ),
        "finite_and_ideal_vertices": (
            vinberg_case["finite_vertex_count"] == 1252
            and vinberg_case["ideal_vertex_count"] == 64
        ),
    }
    if not all(vinberg_checks.values()):
        raise ArithmeticError(f"inherited Vinberg audit failed: {vinberg_checks}")

    checks = {
        "lattice_and_K3_embedding": all(lattice_checks.values()),
        "all_facet_records": all(all(row["checks"].values()) for row in records),
        "full_generation_and_interior": all(generation_checks.values()),
        "literal_Vinberg_finite_volume": all(vinberg_checks.values()),
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-07",
        "id": "A18-C03-H11-28-UD4A1X6-SCALED-DUAL-K3-CARRIER",
        "upstream": {
            key: {
                "path": str(path.relative_to(ROOT)),
                "sha256": EXPECTED_SHA256[key],
            }
            for key, path in paths.items()
        },
        "theorem": (
            "The finite-volume noncompact A18 C03 H11/28 chamber has a "
            "second full-Picard complex-K3 divisor-lattice carrier. After "
            "its unique reflection-stable odd index-four closure M, the "
            "doubled dual N=M^vee(2) is the even 2-elementary lattice of "
            "type (12,8,1), embeds primitively in the K3 lattice by the "
            "displayed discriminant graph glue, and occurs as NS of a "
            "projective K3; all 28 positive halfspaces and reflections "
            "transport and their primitive mirrors generate NS."
        ),
        "stable_closure": {
            "operation": "unique reflection-stable odd index-four closure before duality",
            "basis_columns_in_original_fixed_Q_space": _payload(closure_basis),
            "gram": _payload(gram_m),
            "determinant": int(gram_m.det()),
            "absolute_determinant": abs(int(gram_m.det())),
            "smith_invariants": _smith(gram_m),
            "parity": "odd",
            "primitive_mirror_types_square_divisibility": _histogram_payload(
                closure_types
            ),
        },
        "scaled_dual_K3_bridge": {
            "operation": "N=M^vee(2), with Gram 2*Gram(M)^(-1)",
            "gram_N": _payload(gram_n),
            "N_determinant": int(gram_n.det()),
            "N_absolute_determinant": abs(int(gram_n.det())),
            "N_smith_invariants": _smith(gram_n),
            "N_parity": "even",
            "N_signature": [1, 11, 0],
            "Nikulin_invariants": {"r": 12, "a": 8, "delta": 1},
            "discriminant_group": "(Z/2)^8",
            "discriminant_quadratic_form": (
                "q_D4(-1) + <-1/2>^6; equivalently the form of "
                "U+D4(-1)+A1(-1)^6"
            ),
            "discriminant_q_mod_2_histogram": _q_histogram(gram_n),
            "standard_representative": {
                "lattice": "U+D4(-1)+A1(-1)^6",
                "gram": _payload(standard_picard),
                "exact_discriminant_form_isometry": standard_discriminant_isometry,
                "uniqueness_criterion": (
                    "even indefinite 2-elementary classification; "
                    "rank 12 >= length(A_N)+2 = 10"
                ),
            },
            "primitive_K3_embedding": k3_glue,
            "period_consequence": {
                "statement": (
                    "By K3 period surjectivity, a generic period in the "
                    "signature-(2,8) orthogonal complement avoids every "
                    "extra integral Noether-Lefschetz hyperplane; the "
                    "resulting projective complex K3 has NS exactly N."
                ),
                "N_polarized_moduli_dimension": 8,
                "full_Picard_rank": 12,
                "not_merely_a_sublattice": True,
            },
            "primitive_target_mirror_types_square_divisibility": (
                _histogram_payload(target_types)
            ),
            "source_to_closure_to_target_type_histogram": (
                _histogram_payload(source_to_target_types)
            ),
            "positive_halfspace_identity": (
                "For oriented primitive q in M, p=Gram(M)q/div_M(q) and "
                "Phi(x)=Gram(M)x satisfy (Phi(x),p)_N="
                "(2/div_M(q))*(x,q)_M."
            ),
            "checks": lattice_checks,
        },
        "facet_records": records,
        "full_lattice_generation": {
            "N_root_matrix_HNF": _payload(hnf_n),
            "generation_index": 1,
            "strict_interior_N_coordinates": [int(value) for value in carrier],
            "strict_interior_square": int((carrier.T * gram_n * carrier)[0]),
            "minimum_wall_pairing": min(carrier_pairings),
            "checks": generation_checks,
        },
        "inherited_literal_Vinberg_audit": {
            "criterion": vinberg_case["criterion"],
            "dimension": vinberg_case["dimension"],
            "facet_count": vinberg_case["facet_count"],
            "finite_vertex_count": vinberg_case["finite_vertex_count"],
            "ideal_vertex_count": vinberg_case["ideal_vertex_count"],
            "elliptic_edge_count": vinberg_case["elliptic_edge_count"],
            "connected_affine_component_count": vinberg_case[
                "connected_affine_component_count"
            ],
            "bad_edge_count": vinberg_case["bad_edge_count"],
            "bad_parabolic_neighborhood_count": vinberg_case[
                "bad_parabolic_neighborhood_count"
            ],
            "classification": vinberg_case["classification"],
            "checks": vinberg_checks,
        },
        "lineage_and_scope": {
            "existing_geometric_carrier": (
                "a special Picard-rank-20 K3 with symplectic Nikulin "
                "involution and fixed rank-12 lattice L of determinant 256"
            ),
            "new_positive_gain": (
                "a second full-Picard rank-12 K3 carrier on the scaled-dual "
                "2-elementary lattice N of determinant 256, with an explicit "
                "primitive K3-lattice glue and all 28 walls transported"
            ),
            "arithmetically_different_from_original_carrier": (
                "L and N have the same determinant but different Smith "
                "invariants: L has 1^6 2^4 4^2, whereas N has 1^4 2^8"
            ),
            "known_abstract_polytope": True,
            "not_a_new_abstract_polytope": True,
            "not_the_complete_K3_nef_chamber": (
                "the imported arrangement has eight primitive square-minus-four "
                "mirrors and is only certified as a divisor-reflection chamber"
            ),
            "not_an_automorphism_fixed_chamber_claim": True,
            "minus_two_divisors_not_asserted_irreducible": True,
            "no_named_or_equation_family_claim": True,
            "worldwide_priority_claim": "none",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    payload = build_certificate()
    if not payload["all_checks_pass"]:
        raise ArithmeticError("A18 C03 scaled-dual K3 carrier certificate failed")
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists():
            raise ArithmeticError(f"missing archive: {args.output}")
        if args.output.read_text(encoding="utf-8") != serialized:
            raise ArithmeticError("archived A18 C03 scaled-dual payload changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={_payload_sha256(payload)}")


if __name__ == "__main__":
    main()
