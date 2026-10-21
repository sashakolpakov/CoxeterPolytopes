#!/usr/bin/env python3
"""Certify generic full-Picard K3 carriers for the two high A20 folds.

The known A20-C2a-H12/19 and A20-C2b-H11/20 chambers first occur as
automorphism-fixed sections on one singular Picard-rank-20 K3.  Here that
special geometry is used only to reconstruct primitive embeddings of the two
saturated fixed kernels in the K3 lattice.  A very-general period then has
full Picard equal to the corresponding chamber lattice.

For both chambers this script transports every oriented primitive wall,
checks index-one wall span and every integral reflection, exhausts the
reflection-stable discriminant-module overclosures, and pins the independent
literal Vinberg and GMP exact-cone volume certificates.  No automorphism,
full-nef, or equation-family claim is made for the generic carriers.
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

from scripts.certify_a18_c03_c04_special_k3_symplectic_lifts import (
    saturated_row_kernel,
)
from scripts.certify_a20_c2a_c2b_k3_automorphism_lifts import (
    discriminant_action,
    glue_k3_lattice,
    matrix_digest,
    two_elementary_basis,
)
from scripts.certify_alexeev_remaining_stable_overlattices import (
    _overlattice_basis,
    _stable_subgroups,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
INVENTORY = (
    ROOT
    / "data/enumerations/"
    "alexeev_20_2_1_15_9_1_wall_symmetry_inventory.json"
)
CONSTRUCTION = (
    ROOT
    / "data/enumerations/"
    "alexeev_20_2_1_15_9_1_fixed_section_census.json"
)
EXACT_CONES = (
    ROOT
    / "data/enumerations/"
    "alexeev_20_2_1_15_9_1_exact_cone_census.json"
)
VINBERG = (
    ROOT
    / "data/enumerations/"
    "alexeev_20_2_1_15_9_1_literal_vinberg_census.json"
)
SPECIAL_LIFT = (
    ROOT / "data/enumerations/a20_c2a_c2b_k3_automorphism_lifts.json"
)
OUTPUT = (
    ROOT / "data/enumerations/a20_pair_full_picard_k3_carriers.json"
)

SCHEMA = "coxeter.a20_pair.full_picard_k3_carriers.v1"
EXPECTED_SHA256 = {
    "inventory": "20a257e60b768571a7e7c8aa64e60350bdf7940f5d221cfe031fc2e6ad879e09",
    "construction": "0294897befa7357f485d45df17d5e3e203656a3605c9c3231ab362d3a30754a6",
    "exact_cones": "0279262bdef519c18793d9f565bb91a730e0dc6bad5fb2ed2edc2a9ec371b9d1",
    "vinberg": "1ca27c1287a05824fcbede3d13a0351ea9f3bc7889463bb024229faaf847a283",
    "special_lift": "8ad3566d70dd3b253d27510a3fb0d8daf0271afa20154656f6dc2f7fa28a3314",
}

CONFIGURATIONS = {
    "A20_2_1_C01": {
        "ledger_id": "A20-C2a-H12-19",
        "dimension": 12,
        "facets": 19,
        "rank": 13,
        "determinant": 256,
        "smith": [1] * 6 + [2] * 6 + [4],
        "transcendental_action": sp.Matrix([[0, -1], [1, 0]]),
        "moduli_dimension": 7,
        "mirror_types": {(-2, 1): 7, (-4, 2): 12},
        "interior_square": 32068,
        "interior_pairings": [
            4,
            8,
            8,
            8,
            8,
            8,
            8,
            180,
            4,
            180,
            4,
            4,
            4,
            4,
            112,
            396,
            176,
            780,
            780,
        ],
        "finite_vertices": 527,
        "cusps": 25,
        "edges": 3846,
        "endpoint_histogram": {"ff": 2589, "fi": 1146, "ii": 111},
        "affine_components": 59,
    },
    "A20_2_1_C02": {
        "ledger_id": "A20-C2b-H11-20",
        "dimension": 11,
        "facets": 20,
        "rank": 12,
        "determinant": -1024,
        "smith": [1] * 2 + [2] * 10,
        "transcendental_action": -sp.eye(2),
        "moduli_dimension": 8,
        "mirror_types": {(-2, 1): 3, (-2, 2): 6, (-4, 2): 11},
        "interior_square": 5386,
        "interior_pairings": [
            4,
            2,
            4,
            4,
            2,
            58,
            4,
            58,
            4,
            4,
            4,
            2,
            2,
            88,
            230,
            88,
            218,
            260,
            454,
            540,
        ],
        "finite_vertices": 311,
        "cusps": 29,
        "edges": 2229,
        "endpoint_histogram": {"ff": 1312, "fi": 797, "ii": 120},
        "affine_components": 64,
    },
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


def _primitive_integer(vector: sp.Matrix) -> sp.Matrix:
    denominator = int(sp.ilcm(*(sp.Rational(value).q for value in vector)))
    integral = sp.Matrix([int(sp.Rational(value) * denominator) for value in vector])
    content = _content(integral)
    return integral / content


def _stable_overclosure_audit(
    gram: sp.Matrix, reflections: list[sp.Matrix]
) -> dict[str, Any]:
    raw = _stable_subgroups(gram, reflections)
    groups = raw.pop("groups_internal")
    left_inverse = raw.pop("left_smith_inverse")
    maxima = raw["maximal_stable_subgroup_positions_zero_based"]
    if not (
        raw["stable_integral_subgroup_orders"] == [1, 2]
        and raw["stable_even_subgroup_orders"] == [1]
        and len(maxima) == 1
    ):
        raise ArithmeticError("A20 stable-overlattice denominator changed")

    position = maxima[0]
    maximum = groups[position]
    basis = _overlattice_basis(
        gram,
        left_inverse,
        raw["stable_subgroups"][position]["cosets_in_smith_coordinates"],
    )
    maximum_gram = sp.simplify(basis.T * gram * basis)
    actions = [sp.simplify(basis.inv() * reflection * basis) for reflection in reflections]
    checks = {
        "exhaustive_integral_orders_1_2": raw[
            "stable_integral_subgroup_orders"
        ]
        == [1, 2],
        "only_even_order_is_trivial": raw["stable_even_subgroup_orders"] == [1],
        "unique_maximum_has_index_two": len(maximum) == 2,
        "index_square_determinant_law": (
            maximum_gram.det() == gram.det() / 4
            and abs(sp.Rational(1, basis.det())) == 2
        ),
        "maximum_is_odd": any(
            maximum_gram[index, index] % 2 for index in range(gram.rows)
        ),
        "all_reflections_preserve_maximum": all(
            all(value.is_Integer for value in action)
            and abs(int(action.det())) == 1
            and action.T * maximum_gram * action == maximum_gram
            for action in actions
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(
            [key for key, value in checks.items() if not value]
        )
    return {
        "exhaustive_discriminant_module_audit": raw,
        "unique_maximal_integral_overclosure": {
            "index": 2,
            "basis_columns_in_original_lattice_Q": _matrix_payload(basis),
            "gram": _matrix_payload(maximum_gram),
            "determinant": int(maximum_gram.det()),
            "smith_invariants": _smith(maximum_gram),
            "parity": "odd",
        },
        "K3_admissible_conclusion": (
            "the original fixed lattice is the unique even reflection-stable "
            "overclosure; the sole proper stable integral maximum is odd and "
            "cannot be a K3 Neron--Severi lattice"
        ),
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    paths = {
        "inventory": INVENTORY,
        "construction": CONSTRUCTION,
        "exact_cones": EXACT_CONES,
        "vinberg": VINBERG,
        "special_lift": SPECIAL_LIFT,
    }
    for key, path in paths.items():
        actual = _sha256(path)
        if actual != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} archive changed: {actual}")

    inventory_data = json.loads(INVENTORY.read_text(encoding="utf-8"))
    construction_data = json.loads(CONSTRUCTION.read_text(encoding="utf-8"))
    cone_data = json.loads(EXACT_CONES.read_text(encoding="utf-8"))
    vinberg_data = json.loads(VINBERG.read_text(encoding="utf-8"))
    lift_data = json.loads(SPECIAL_LIFT.read_text(encoding="utf-8"))
    if not lift_data["all_checks_pass"]:
        raise ArithmeticError("the special A20 K3 lift archive is not green")

    inventory = next(
        row for row in inventory_data["examples"] if row["key"] == [20, 2, 1]
    )
    n_gram = sp.Matrix(inventory["ambient"]["gram"])
    if not (
        n_gram.shape == (20, 20)
        and n_gram.det() == -4
        and exact_inertia(n_gram) == (1, 19, 0)
    ):
        raise ArithmeticError("the ambient A20 Picard lattice changed")

    discriminant_basis = two_elementary_basis(n_gram, 2)
    t_gram = sp.diag(2, 2)
    identification = sp.Matrix([[0, 1], [1, 1]])
    sources = {
        key: {
            "path": str(path.relative_to(ROOT)),
            "sha256": EXPECTED_SHA256[key],
        }
        for key, path in paths.items()
    }

    cases = []
    common_k3_gram: sp.Matrix | None = None
    common_k3_basis: sp.Matrix | None = None
    for case_id, config in CONFIGURATIONS.items():
        source = next(
            row for row in construction_data["cases"] if row["id"] == case_id
        )
        cone = next(row for row in cone_data["cases"] if row["id"] == case_id)
        vinberg = next(
            row for row in vinberg_data["cases"] if row["id"] == case_id
        )
        special = next(
            row for row in lift_data["cases"] if row["id"] == case_id
        )

        gram = sp.Matrix(source["saturated_fixed_lattice"]["gram"])
        fixed_basis = sp.Matrix(
            source["saturated_fixed_lattice"]["basis_primal_rows"]
        )
        action = sp.Matrix(
            source["case"]["representative_nonidentity_primal_actions"][0]
        )
        t_action = config["transcendental_action"]
        k3_basis, k3_gram, k3_action = glue_k3_lattice(
            n_gram,
            action,
            discriminant_basis,
            t_gram,
            t_action,
            identification,
        )
        if common_k3_gram is None:
            common_k3_gram = k3_gram
            common_k3_basis = k3_basis
        elif not (
            k3_gram == common_k3_gram and k3_basis == common_k3_basis
        ):
            raise ArithmeticError("the two A20 kernels left the common K3 lattice")

        extended_fixed_basis = fixed_basis.row_join(sp.zeros(config["rank"], 2))
        embedding_rows = sp.simplify(extended_fixed_basis * k3_basis.inv())
        complement_rows = saturated_row_kernel(k3_gram * embedding_rows.T)
        complement_gram = sp.simplify(
            complement_rows * k3_gram * complement_rows.T
        )
        combined = embedding_rows.col_join(complement_rows)
        embedding_checks = {
            "marked_K3_lattice_even_unimodular_signature_3_19": (
                k3_gram.det() == -1
                and exact_inertia(k3_gram) == (3, 19, 0)
                and all(k3_gram[index, index] % 2 == 0 for index in range(22))
            ),
            "embedding_integral": all(value.is_Integer for value in embedding_rows),
            "embedding_isometric": embedding_rows * k3_gram * embedding_rows.T
            == gram,
            "embedding_primitive_Smith_all_one": _smith(embedding_rows)
            == [1] * config["rank"],
            "orthogonal_complement_saturated": _smith(complement_rows)
            == [1] * (22 - config["rank"]),
            "orthogonal_complement_signature": exact_inertia(complement_gram)
            == (2, 20 - config["rank"], 0),
            "orthogonal_complement_determinant_order": abs(
                int(complement_gram.det())
            )
            == abs(int(gram.det())),
            "orthogonal": complement_rows * k3_gram * embedding_rows.T
            == sp.zeros(22 - config["rank"], config["rank"]),
            "combined_gluing_index_equals_discriminant_order": abs(
                int(combined.det())
            )
            == abs(int(gram.det())),
            "special_lift_reconstructed": (
                matrix_digest(k3_gram)
                == lift_data["common_marked_k3_lattice"]["glued_gram_sha256"]
                and matrix_digest(k3_basis)
                == lift_data["common_marked_k3_lattice"]["glued_basis_sha256"]
                and matrix_digest(k3_action) == special["k3_action_sha256"]
            ),
            "fixed_kernel_saturated_in_ambient_NS": (
                source["saturated_fixed_lattice"]["constraint_clearing_denominator"]
                == 1
            ),
        }
        if not all(embedding_checks.values()):
            raise ArithmeticError(
                f"{case_id} embedding failed: "
                f"{[key for key, value in embedding_checks.items() if not value]}"
            )

        active = set(source["exact_reduction"]["active_facet_orbit_indices"])
        interior_rational = sp.Matrix(
            [
                sp.Rational(value)
                for value in source["strict_fixed_interior"][
                    "fixed_lattice_coordinates"
                ]
            ]
        )
        interior = _primitive_integer(interior_rational)
        roots = []
        wall_records = []
        types: Counter[tuple[int, int]] = Counter()
        pairings = []
        reflections = []
        for record in source["restricted_mirrors"]["records"]:
            if int(record["orbit_index"]) not in active:
                continue
            sign = 1 if sp.Rational(record["presentation_scale_from_primitive"]) > 0 else -1
            root = sign * sp.Matrix(
                record["primitive_integral_same_direction_coordinates"]
            )
            square = int((root.T * gram * root)[0])
            divisibility = _divisibility(gram, root)
            pairing = int((interior.T * gram * root)[0])
            reflection = _reflection(gram, root)
            checks = {
                "primitive": _content(root) == 1,
                "square_preserved": square == int(record["square"]),
                "strict_positive_orientation": pairing > 0,
                "reflection_integral": all(value.is_Integer for value in reflection),
                "reflection_unimodular": reflection.det() == -1,
                "reflection_involutive": reflection**2 == sp.eye(config["rank"]),
                "reflection_preserves_Gram": reflection.T * gram * reflection
                == gram,
                "reflection_negates_root": reflection * root == -root,
            }
            if not all(checks.values()):
                raise ArithmeticError(
                    f"{case_id} wall {record['orbit_index']} failed"
                )
            roots.append(root)
            reflections.append(reflection)
            pairings.append(pairing)
            types[(square, divisibility)] += 1
            wall_records.append(
                {
                    "facet_one_based": len(wall_records) + 1,
                    "source_orbit_index": record["orbit_index"],
                    "source_walls": record["source_wall_indices"],
                    "orientation_sign_from_presentation_scale": sign,
                    "primitive_fixed_lattice_coordinates": _vector_payload(root),
                    "square": square,
                    "divisibility_in_full_Picard": divisibility,
                    "strict_interior_pairing": pairing,
                    "checks": checks,
                }
            )

        root_rows = sp.Matrix.hstack(*roots).T
        span_smith = _smith(root_rows)
        facet_gram = sp.Matrix(
            [
                [int((left.T * gram * right)[0]) for right in roots]
                for left in roots
            ]
        )
        type_payload = {
            f"{square},{divisibility}": count
            for (square, divisibility), count in sorted(types.items())
        }
        expected_types = {
            f"{square},{divisibility}": count
            for (square, divisibility), count in sorted(
                config["mirror_types"].items()
            )
        }
        wall_checks = {
            "complete_active_wall_count": len(roots) == config["facets"],
            "facet_Gram_literal": facet_gram
            == sp.Matrix(source["restricted_mirrors"]["facet_gram"]),
            "wall_span_index_one": (
                root_rows.rank() == config["rank"]
                and span_smith == [1] * config["rank"]
            ),
            "mirror_type_histogram": type_payload == expected_types,
            "primitive_strict_interior": _content(interior) == 1,
            "strict_interior_square": int((interior.T * gram * interior)[0])
            == config["interior_square"],
            "strict_interior_pairings": pairings == config["interior_pairings"],
            "all_reflections_green": all(
                all(record["checks"].values()) for record in wall_records
            ),
        }
        if not all(wall_checks.values()):
            raise ArithmeticError(
                f"{case_id} wall system failed: "
                f"{[key for key, value in wall_checks.items() if not value]}"
            )

        overclosures = _stable_overclosure_audit(gram, reflections)
        literal = vinberg["literal_endpoint_audit"]
        exact = cone["exact_double_description"]
        volume_checks = {
            "literal_classification_finite_volume_noncompact": (
                vinberg["decision"]["classification"]
                == "finite_volume_noncompact"
                and vinberg["decision"]["finite_volume"] is True
                and vinberg["decision"]["compact"] is False
            ),
            "literal_counts": (
                literal["finite_vertex_count"],
                literal["ideal_vertex_count"],
                literal["edge_count"],
            )
            == (
                config["finite_vertices"],
                config["cusps"],
                config["edges"],
            ),
            "literal_endpoint_histogram": literal[
                "edge_endpoint_type_histogram"
            ]
            == config["endpoint_histogram"],
            "every_literal_edge_has_two_endpoints": literal[
                "edge_endpoint_count_histogram"
            ]
            == {"2": config["edges"]},
            "zero_literal_failures": (
                literal["bad_edge_count"] == 0
                and literal["bad_parabolic_neighborhood_count"] == 0
                and all(vinberg["checks"].values())
            ),
            "connected_affine_component_count": vinberg[
                "diagram_enumeration"
            ]["connected_affine_component_count"]
            == config["affine_components"],
            "independent_GMP_counts_and_no_hyperideal": (
                exact["finite_vertex_count"] == config["finite_vertices"]
                and exact["ideal_vertex_count"] == config["cusps"]
                and exact["hyperideal_vertex_count"] == 0
                and exact["extreme_ray_count"]
                == config["finite_vertices"] + config["cusps"]
                and all(cone["checks"].values())
            ),
            "same_complete_facet_count_in_both_lanes": (
                vinberg["facet_count"] == cone["facet_count"] == config["facets"]
            ),
        }
        if not all(volume_checks.values()):
            raise ArithmeticError(
                f"{case_id} volume handoff failed: "
                f"{[key for key, value in volume_checks.items() if not value]}"
            )

        case_checks = {
            "source_fixed_lattice_invariants": (
                gram.rows == config["rank"]
                and int(gram.det()) == config["determinant"]
                and _smith(gram) == config["smith"]
                and exact_inertia(gram) == (1, config["dimension"], 0)
            ),
            "primitive_embedding_green": all(embedding_checks.values()),
            "wall_system_green": all(wall_checks.values()),
            "stable_overclosure_denominator_green": all(
                overclosures["checks"].values()
            ),
            "volume_green": all(volume_checks.values()),
        }
        cases.append(
            {
                "id": case_id,
                "ledger_id": config["ledger_id"],
                "theorem": (
                    f"a very-general marked K3 in a {config['moduli_dimension']}-"
                    f"dimensional {config['ledger_id']} lattice-polarized period "
                    "family has full Picard equal to the chamber lattice and "
                    "carries every archived primitive facet reflection"
                ),
                "carrier_lattice": {
                    "name": (
                        "U+E8(-2)+<-2>^2"
                        if case_id == "A20_2_1_C02"
                        else "saturated A20-C2a rank-13 fixed lattice"
                    ),
                    "gram": _matrix_payload(gram),
                    "rank": gram.rows,
                    "signature": list(exact_inertia(gram)),
                    "determinant": int(gram.det()),
                    "smith_invariants": _smith(gram),
                    "even": all(
                        gram[index, index] % 2 == 0 for index in range(gram.rows)
                    ),
                    "primitive_embedding_in_Lambda_K3": {
                        "ambient_marked_K3_gram": _matrix_payload(k3_gram),
                        "embedding_matrix_rows": _matrix_payload(embedding_rows),
                        "orthogonal_complement_matrix_rows": _matrix_payload(
                            complement_rows
                        ),
                        "orthogonal_complement_gram": _matrix_payload(
                            complement_gram
                        ),
                        "orthogonal_complement_signature": list(
                            exact_inertia(complement_gram)
                        ),
                        "combined_gluing_index": abs(int(combined.det())),
                        "checks": embedding_checks,
                    },
                    "very_general_period_family": {
                        "dimension": config["moduli_dimension"],
                        "full_Picard": True,
                        "claim": "NS(X) equals the displayed carrier lattice",
                        "reason": (
                            "period surjectivity and avoidance of the countable "
                            "union of extra algebraic hyperplanes in the "
                            "orthogonal period domain"
                        ),
                    },
                },
                "oriented_chamber": {
                    "hyperbolic_dimension": config["dimension"],
                    "facet_count": config["facets"],
                    "wall_records": wall_records,
                    "wall_span_rank": root_rows.rank(),
                    "wall_span_smith": span_smith,
                    "wall_span_index": math.prod(span_smith),
                    "primitive_mirror_type_histogram": type_payload,
                    "facet_gram_sha256": _payload_sha256(
                        _matrix_payload(facet_gram)
                    ),
                    "strict_primitive_interior_coordinates": _vector_payload(
                        interior
                    ),
                    "strict_interior_square": config["interior_square"],
                    "strict_interior_pairings": pairings,
                    "minimum_wall_pairing": min(pairings),
                    "checks": wall_checks,
                },
                "reflection_stable_overclosures": overclosures,
                "finite_volume": {
                    "classification": "finite_volume_noncompact",
                    "proof_lanes": [
                        "literal Vinberg all-edge/full-affine-neighborhood audit",
                        "independent GMP rational double-description",
                    ],
                    "finite_vertices": config["finite_vertices"],
                    "cusps": config["cusps"],
                    "hyperideal_vertices": 0,
                    "elliptic_rank_dimension_minus_one_edges": config["edges"],
                    "endpoint_type_histogram": config["endpoint_histogram"],
                    "connected_affine_components": config["affine_components"],
                    "bad_edges": 0,
                    "bad_parabolic_neighborhoods": 0,
                    "checks": volume_checks,
                },
                "scope": {
                    "positive_gain": (
                        "new project-level very-general full-Picard period "
                        "carrier of an already certified chamber"
                    ),
                    "not_new_abstract_polytope": True,
                    "not_equation_level": True,
                    "not_full_nef": (
                        "square-minus-four mirrors are mixed divisor cuts; no "
                        "claim identifies the chamber with the generic K3 nef cone"
                    ),
                    "not_generic_surface_automorphism": (
                        "the source order-four automorphism or involution is a "
                        "special-period embedding witness and is not asserted on "
                        "the returned very-general K3"
                    ),
                    "not_reflection_automorphisms": (
                        "integral Picard reflections are not asserted to extend "
                        "to Hodge isometries or surface automorphisms"
                    ),
                },
                "checks": case_checks,
                "all_checks_pass": all(case_checks.values()),
            }
        )

    if common_k3_gram is None or common_k3_basis is None:
        raise ArithmeticError("no A20 cases were constructed")
    checks = {
        "all_source_hashes_pinned": all(
            _sha256(path) == EXPECTED_SHA256[key]
            for key, path in paths.items()
        ),
        "both_cases_green": len(cases) == 2
        and all(case["all_checks_pass"] for case in cases),
        "common_marked_K3_reconstructed": (
            matrix_digest(common_k3_gram)
            == lift_data["common_marked_k3_lattice"]["glued_gram_sha256"]
            and matrix_digest(common_k3_basis)
            == lift_data["common_marked_k3_lattice"]["glued_basis_sha256"]
        ),
        "both_unique_even_stable_carriers_are_original_kernels": all(
            case["reflection_stable_overclosures"][
                "exhaustive_discriminant_module_audit"
            ]["stable_even_subgroup_orders"]
            == [1]
            for case in cases
        ),
    }
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "The A20-C2a-H12/19 and A20-C2b-H11/20 chambers each admit a "
            "very-general full-Picard complex-K3 divisor-reflection carrier. "
            "Their saturated fixed kernels are primitively embedded in one "
            "marked K3 lattice; every wall spans and reflects the full carrier."
        ),
        "sources": sources,
        "common_special_period_witness": {
            "ambient_Picard": "the even 2-elementary (20,2,1) lattice",
            "same_singular_K3": True,
            "source_actions": [
                "non-symplectic order-four automorphism for C01",
                "non-symplectic involution for C02",
            ],
            "role_in_this_theorem": (
                "primitive-embedding witness only; the generic carrier theorem "
                "does not deform either automorphism"
            ),
            "common_marked_Lambda_K3_gram_sha256": matrix_digest(
                common_k3_gram
            ),
        },
        "cases": cases,
        "paired_gain": {
            "new_full_Picard_period_carriers": 2,
            "geometric_slot_count_change": 0,
            "reason": (
                "both chambers already had genuine special-period K3 "
                "automorphism carriers; the gain is deformation to full Picard"
            ),
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
        raise ArithmeticError("A20 paired full-Picard carriers failed")
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(
            encoding="utf-8"
        ) != rendered:
            raise ArithmeticError("archived A20 paired carrier changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={payload['payload_sha256']}")


if __name__ == "__main__":
    main()
