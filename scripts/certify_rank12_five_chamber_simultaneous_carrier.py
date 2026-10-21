#!/usr/bin/env python3
"""Certify one order-four K3 carrier for all five rank-12 H10 chambers.

The Brandhorst--Hofmann row 0.4.4.1 has generic Neron--Severi lattice

    N = U + E8(-1) + A1(-1)^2

and its order-four generator induces an involution on N.  In the pinned BH
basis, five invariant full-rank divisor lattices are obtained by imposing the
following coordinate parities:

    H10/12: none
    H10/15: coordinate 12 even
    H10/16: coordinate 10 even
    H10/38: coordinates 7 and 12 even
    H10/45: coordinates 10 and 12 even.

The chamber action is identified independently in each sublattice by the
complete rank-11 involution census of its root-chamber group.  Dedicated
H10/16 and H10/45 action/glue packages are retained as corroboration, not as
the simultaneity proof.  This script performs the common exact lattice and
Vinberg replay and records the simultaneous-carrier theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_alexeev_remaining_stable_overlattices import (
    _discriminant_data,
    _q_histogram,
)
from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    saturated_row_kernel,
    smith_diagonal,
)
from scripts.certify_roulleau_rank12_bh_k3_lifts import NS_GRAMS, PURE_ACTION
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
GRAPH_SOURCE = ROOT / "data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json"
D4_SOURCE = ROOT / "data/enumerations/roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json"
BH_SOURCE = ROOT / "data/enumerations/roulleau_rank12_bh_k3_lifts.json"
WEIERSTRASS_SOURCE = ROOT / "data/enumerations/ue8a1x2_order4_weierstrass_family.json"
H10_16_CARRIER_SOURCE = (
    ROOT
    / "data/enumerations/r12_h10_16_order4_divisor_overlattice_carrier.json"
)
H10_45_CARRIER_SOURCE = (
    ROOT / "data/enumerations/r12_h10_45_bh_divisor_sublattice_carrier.json"
)
OUTPUT = ROOT / "data/enumerations/rank12_five_chamber_simultaneous_carrier.json"


CASES: tuple[dict[str, Any], ...] = (
    {
        "short": "H10/12",
        "id": "R12-E8A1x2-C2-H10-12",
        "even_coordinates_one_based": [],
        "index": 1,
        "determinant": -4,
        "smith": [*[1] * 10, 2, 2],
        "q": {"0": 1, "1": 1, "3/2": 2},
        "fixed_q": {"0": 1, "1": 1, "7/4": 2},
        "facets": 12,
        "finite_vertices": 17,
        "cusps": 2,
        "elliptic_edges": 99,
        "parabolic_audits": 3,
        "source_kind": "graph",
        "stock": "R12_UE8A1x2",
        "class_index": 2,
        "chamber_census_key": "R12_UE8A1x2",
        "chamber_class_index": 2,
        "action_type": "standard_A1_swap_on_U_E8_A1x2",
    },
    {
        "short": "H10/15",
        "id": "R12-D8A1x2-C2a-H10-15",
        "even_coordinates_one_based": [12],
        "index": 2,
        "determinant": -16,
        "smith": [*[1] * 8, *[2] * 4],
        "q": {"0": 4, "1": 4, "1/2": 2, "3/2": 6},
        "fixed_q": {"0": 4, "1": 4, "3/4": 2, "7/4": 6},
        "facets": 15,
        "finite_vertices": 42,
        "cusps": 10,
        "elliptic_edges": 294,
        "parabolic_audits": 17,
        "source_kind": "graph",
        "stock": "R12_UD8A1x2",
        "class_index": 2,
        "chamber_census_key": "R12_UD8A1x2",
        "chamber_class_index": 2,
        "action_type": "standard_A1_swap_on_U_D8_A1x2",
    },
    {
        "short": "H10/16",
        "id": "R12-D8A1x2-C2b-H10-16",
        "even_coordinates_one_based": [10],
        "index": 2,
        "determinant": -16,
        "smith": [*[1] * 8, *[2] * 4],
        "q": {"0": 4, "1": 4, "1/2": 2, "3/2": 6},
        "fixed_q": {
            "0": 2,
            "1": 2,
            "1/2": 2,
            "1/4": 2,
            "3/2": 2,
            "5/4": 2,
            "7/4": 4,
        },
        "facets": 16,
        "finite_vertices": 70,
        "cusps": 6,
        "elliptic_edges": 426,
        "parabolic_audits": 14,
        "source_kind": "graph",
        "stock": "R12_UD8A1x2",
        "class_index": 3,
        "chamber_census_key": "R12_UD8A1x2",
        "chamber_class_index": 3,
        "action_type": "alternative_D8_spinor_swap_pair",
    },
    {
        "short": "H10/38",
        "id": "R12-D4x2A1x2-p67-H10-38",
        "even_coordinates_one_based": [7, 12],
        "index": 4,
        "determinant": -64,
        "smith": [*[1] * 6, *[2] * 6],
        "q": {"0": 16, "1": 16, "1/2": 12, "3/2": 20},
        "fixed_q": {"0": 16, "1": 16, "3/4": 12, "7/4": 20},
        "facets": 38,
        "finite_vertices": 526,
        "cusps": 141,
        "elliptic_edges": 4518,
        "parabolic_audits": 307,
        "source_kind": "d4",
        "source_id": "R12-D4x2A1x2-p67-H10-38",
        "chamber_census_key": "R12_UD4x2A1x2",
        "chamber_class_index": 4,
        "action_type": "standard_A1_swap_on_U_D4x2_A1x2",
    },
    {
        "short": "H10/45",
        "id": "R12-D4x2A1x2-S6transposition-H10-45",
        "even_coordinates_one_based": [10, 12],
        "index": 4,
        "determinant": -64,
        "smith": [*[1] * 6, *[2] * 6],
        "q": {"0": 16, "1": 16, "1/2": 12, "3/2": 20},
        "fixed_q": {
            "0": 8,
            "1": 8,
            "1/2": 8,
            "1/4": 8,
            "3/2": 8,
            "3/4": 4,
            "5/4": 8,
            "7/4": 12,
        },
        "facets": 45,
        "finite_vertices": 1048,
        "cusps": 93,
        "elliptic_edges": 7332,
        "parabolic_audits": 302,
        "source_kind": "d4",
        "source_id": "R12-D4x2A1x2-S6transposition-H10-45",
        "chamber_census_key": "R12_UD4x2A1x2",
        "chamber_class_index": 12,
        "action_type": "alternative_D4_outer_swap_pair",
    },
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(value) for value in matrix.row(i)] for i in range(matrix.rows)]


def sublattice_basis(even_coordinates_one_based: list[int]) -> sp.Matrix:
    basis = sp.eye(12)
    for coordinate in even_coordinates_one_based:
        basis[coordinate - 1, coordinate - 1] = 2
    return basis


def _xor(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a ^ b for a, b in zip(left, right))


def _binary_discriminant_isometry(
    source_gram: sp.Matrix, target_gram: sp.Matrix
) -> dict[str, Any] | None:
    """Find an exact isometry between small elementary 2-discriminant forms.

    A q-histogram is only a fingerprint, not an isometry certificate.  All
    lattices in this theorem have elementary 2-discriminant groups of length
    at most six, so a basis backtrack can preserve q and every bilinear
    pairing exactly and then replay all group elements.
    """

    source_moduli, _, source_discriminant_gram = _discriminant_data(source_gram)
    target_moduli, _, target_discriminant_gram = _discriminant_data(target_gram)
    source_indices = [
        index for index, modulus in enumerate(source_moduli) if modulus > 1
    ]
    target_indices = [
        index for index, modulus in enumerate(target_moduli) if modulus > 1
    ]
    if not (
        len(source_indices) == len(target_indices)
        and all(source_moduli[index] == 2 for index in source_indices)
        and all(target_moduli[index] == 2 for index in target_indices)
    ):
        return None
    dimension = len(source_indices)
    reduced_elements = list(itertools.product((0, 1), repeat=dimension))

    def full_vector(
        reduced: tuple[int, ...], indices: list[int], rank: int
    ) -> sp.Matrix:
        vector = sp.zeros(rank, 1)
        for bit, index in zip(reduced, indices):
            vector[index] = bit
        return vector

    source_basis = [
        tuple(int(position == index) for position in range(dimension))
        for index in range(dimension)
    ]

    def q_value(
        reduced: tuple[int, ...],
        indices: list[int],
        discriminant_gram: sp.Matrix,
        rank: int,
    ) -> sp.Rational:
        vector = full_vector(reduced, indices, rank)
        return sp.Rational((vector.T * discriminant_gram * vector)[0]) % 2

    def b_value(
        left: tuple[int, ...],
        right: tuple[int, ...],
        indices: list[int],
        discriminant_gram: sp.Matrix,
        rank: int,
    ) -> sp.Rational:
        left_vector = full_vector(left, indices, rank)
        right_vector = full_vector(right, indices, rank)
        return sp.Rational(
            (left_vector.T * discriminant_gram * right_vector)[0]
        ) % 1

    source_q = {
        vector: q_value(
            vector, source_indices, source_discriminant_gram, source_gram.rows
        )
        for vector in reduced_elements
    }
    target_q = {
        vector: q_value(
            vector, target_indices, target_discriminant_gram, target_gram.rows
        )
        for vector in reduced_elements
    }
    source_b = {
        (left, right): b_value(
            left,
            right,
            source_indices,
            source_discriminant_gram,
            source_gram.rows,
        )
        for left in source_basis
        for right in source_basis
    }

    def span(vectors: list[tuple[int, ...]]) -> set[tuple[int, ...]]:
        result = {tuple(0 for _ in range(dimension))}
        for vector in vectors:
            result |= {_xor(element, vector) for element in tuple(result)}
        return result

    images: list[tuple[int, ...]] = []

    def backtrack(index: int) -> bool:
        if index == dimension:
            return True
        forbidden = span(images)
        source_vector = source_basis[index]
        for candidate in reduced_elements:
            if candidate in forbidden or target_q[candidate] != source_q[source_vector]:
                continue
            if any(
                b_value(
                    candidate,
                    previous_image,
                    target_indices,
                    target_discriminant_gram,
                    target_gram.rows,
                )
                != source_b[(source_vector, source_basis[previous_index])]
                for previous_index, previous_image in enumerate(images)
            ):
                continue
            images.append(candidate)
            if backtrack(index + 1):
                return True
            images.pop()
        return False

    if not backtrack(0):
        return None

    image_table: dict[tuple[int, ...], tuple[int, ...]] = {}
    for element in reduced_elements:
        image = tuple(0 for _ in range(dimension))
        for bit, basis_image in zip(element, images):
            if bit:
                image = _xor(image, basis_image)
        image_table[element] = image
    if not (
        len(set(image_table.values())) == 2**dimension
        and all(source_q[element] == target_q[image] for element, image in image_table.items())
    ):
        raise AssertionError("a reconstructed discriminant isometry failed replay")
    return {
        "group": f"(Z/2)^{dimension}",
        "basis_images_as_columns": [list(image) for image in images],
        "all_elements_checked": 2**dimension,
        "quadratic_form_preserved_mod_2Z": True,
        "bilinear_form_preserved_mod_Z": True,
        "bijective": True,
    }


def source_lattice_gram(
    specification: dict[str, Any], graph: dict[str, Any], d4: dict[str, Any]
) -> sp.Matrix:
    if specification["source_kind"] == "graph":
        stock = next(
            row for row in graph["stocks"] if row["stock_key"] == specification["stock"]
        )
        return sp.Matrix(stock["basis_gram"])
    return sp.Matrix(d4["coordinate_reconstruction"]["basis_gram"])


def _lift_wall_involution(
    wall_rows: sp.Matrix, gram: sp.Matrix, permutation: list[int]
) -> sp.Matrix:
    rank = gram.rows
    action = wall_rows[[permutation[index] for index in range(rank)], :]
    expected = wall_rows[
        [permutation[index] for index in range(wall_rows.rows)], :
    ]
    if not (
        all(value.is_Integer for value in action)
        and abs(int(action.det())) == 1
        and action**2 == sp.eye(rank)
        and wall_rows * action == expected
        and action * gram * action.T == gram
    ):
        raise AssertionError("a chamber involution failed its exact lattice lift")
    return action


def rank_eleven_chamber_involution_census(
    graph: dict[str, Any], d4: dict[str, Any]
) -> dict[str, Any]:
    """Recompute every rank-11 involution class in the three chamber groups."""

    families: dict[str, Any] = {}
    for stock_key in ("R12_UE8A1x2", "R12_UD8A1x2"):
        stock = next(row for row in graph["stocks"] if row["stock_key"] == stock_key)
        gram = sp.Matrix(stock["basis_gram"])
        wall_rows = sp.Matrix(stock["wall_coordinate_rows"])
        entries = []
        for class_record in stock["classes"]:
            if not (
                class_record["order"] == 2 and class_record["fixed_rank"] == 11
            ):
                continue
            generator = [
                value - 1
                for value in class_record["representative_generators_one_based"][0]
            ]
            action = _lift_wall_involution(wall_rows, gram, generator)
            fixed_basis = saturated_row_kernel(action - sp.eye(12))
            fixed_gram = fixed_basis * gram * fixed_basis.T
            entries.append(
                {
                    "class_index": class_record["class_index"],
                    "conjugacy_class_size": class_record["conjugacy_class_size"],
                    "fixed_rank": fixed_basis.rows,
                    "fixed_smith": smith_diagonal(fixed_gram),
                    "fixed_discriminant_q_histogram_mod_two": _q_histogram(
                        fixed_gram
                    ),
                }
            )
        families[stock_key] = {
            "complete_chamber_group_order": stock["weighted_diagram_group"]["order"],
            "complete_subgroup_class_count": stock["subgroup_class_count"],
            "rank_eleven_involution_class_count": len(entries),
            "classes": entries,
        }

    gram = sp.Matrix(d4["coordinate_reconstruction"]["basis_gram"])
    wall_rows = sp.Matrix(d4["coordinate_reconstruction"]["wall_coordinate_rows"])
    entries = []
    for class_record in d4["subgroup_census"]["classes"]:
        if not (class_record["order"] == 2 and class_record["fixed_rank"] == 11):
            continue
        generator = class_record["representative_generators_zero_based"][0]
        action = _lift_wall_involution(wall_rows, gram, generator)
        fixed_basis = saturated_row_kernel(action - sp.eye(12))
        fixed_gram = fixed_basis * gram * fixed_basis.T
        entries.append(
            {
                "class_index": class_record["class_index"],
                "conjugacy_class_size": class_record["conjugacy_class_size"],
                "fixed_rank": fixed_basis.rows,
                "fixed_smith": smith_diagonal(fixed_gram),
                "fixed_discriminant_q_histogram_mod_two": _q_histogram(fixed_gram),
            }
        )
    families["R12_UD4x2A1x2"] = {
        "complete_chamber_group_order": d4["weighted_wall_group"]["order"],
        "complete_subgroup_class_count": d4["subgroup_census"]["class_count"],
        "rank_eleven_involution_class_count": len(entries),
        "classes": entries,
    }

    expected_indices = {
        "R12_UE8A1x2": [2],
        "R12_UD8A1x2": [2, 3],
        "R12_UD4x2A1x2": [4, 12],
    }
    checks = {
        "upstream_complete_censuses_green": graph["all_checks_pass"]
        and d4["all_checks_pass"],
        "all_rank_eleven_involution_classes_recovered": all(
            [entry["class_index"] for entry in families[key]["classes"]]
            == expected
            for key, expected in expected_indices.items()
        ),
        "fixed_fingerprints_pairwise_distinct_within_each_chamber_group": all(
            len(
                {
                    (
                        tuple(entry["fixed_smith"]),
                        tuple(
                            sorted(
                                entry[
                                    "fixed_discriminant_q_histogram_mod_two"
                                ].items()
                            )
                        ),
                    )
                    for entry in family["classes"]
                }
            )
            == family["rank_eleven_involution_class_count"]
            for family in families.values()
        ),
    }
    if not all(checks.values()):
        raise AssertionError(f"the rank-11 chamber-action census changed: {checks}")
    return {
        "scope": (
            "all conjugacy classes of order-two subgroups with fixed rank 11 "
            "in each complete weighted root-chamber group"
        ),
        "families": families,
        "checks": checks,
        "all_checks_pass": True,
    }


def source_geometry(
    specification: dict[str, Any], graph: dict[str, Any], d4: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    if specification["source_kind"] == "graph":
        stock = next(
            row for row in graph["stocks"] if row["stock_key"] == specification["stock"]
        )
        case = next(
            row
            for row in stock["survivors"]
            if row["class_index"] == specification["class_index"]
        )
        geometry = case["geometry"]
        vinberg = geometry["vinberg_volume_criterion"]
        summary = {
            "dimension": geometry["hyperbolic_dimension"],
            "facets": geometry["facet_count"],
            "finite_vertices": geometry["finite_vertex_count"],
            "cusps": geometry["ideal_vertex_count"],
            "elliptic_edges": vinberg["elliptic_edge_count"],
            "parabolic_audits": vinberg["parabolic_neighborhood_audit_count"],
            "bad_audits": len(vinberg["finite_volume_bad_edges"]),
            "classification": vinberg["decision"],
            "complete_audit_green": (
                case["all_checks_pass"]
                and vinberg["simple_roots_complete"]
                and bool(vinberg["canonical_full_audit_sha256"])
            ),
        }
        return case["fixed_lattice"], summary
    case = next(row for row in d4["cases"] if row["id"] == specification["source_id"])
    vinberg = case["vinberg_finite_volume_audit"]
    summary = {
        "dimension": vinberg["dimension"],
        "facets": vinberg["facet_count"],
        "finite_vertices": vinberg["finite_vertex_count"],
        "cusps": vinberg["ideal_vertex_count"],
        "elliptic_edges": vinberg["elliptic_edge_count"],
        "parabolic_audits": vinberg["connected_affine_component_count"],
        "bad_audits": vinberg["bad_edge_count"]
        + vinberg["bad_parabolic_neighborhood_count"],
        "classification": "finite_volume_noncompact",
        "complete_audit_green": vinberg["all_checks_pass"]
        and vinberg["finite_volume"]
        and not vinberg["compact"],
    }
    return case["fixed_lattice"], summary


def build_certificate() -> dict[str, Any]:
    graph = json.loads(GRAPH_SOURCE.read_text(encoding="utf-8"))
    d4 = json.loads(D4_SOURCE.read_text(encoding="utf-8"))
    bh = json.loads(BH_SOURCE.read_text(encoding="utf-8"))
    weierstrass = json.loads(WEIERSTRASS_SOURCE.read_text(encoding="utf-8"))
    h10_16_carrier = json.loads(
        H10_16_CARRIER_SOURCE.read_text(encoding="utf-8")
    )
    h10_45_carrier = json.loads(
        H10_45_CARRIER_SOURCE.read_text(encoding="utf-8")
    )
    if not all(
        source.get("all_checks_pass")
        for source in (
            graph,
            d4,
            bh,
            weierstrass,
            h10_16_carrier,
            h10_45_carrier,
        )
    ):
        raise AssertionError("an upstream rank-12 carrier source is no longer green")

    chamber_census = rank_eleven_chamber_involution_census(graph, d4)
    ambient_gram = sp.Matrix(NS_GRAMS["0.4.4.1"])
    ambient_action = sp.Matrix(PURE_ACTION)
    ambient_fixed_basis = saturated_row_kernel(ambient_action - sp.eye(12))
    ambient_coinvariant_basis = saturated_row_kernel(
        ambient_gram * ambient_fixed_basis.T
    )
    ambient_coinvariant_gram = sp.simplify(
        ambient_coinvariant_basis
        * ambient_gram
        * ambient_coinvariant_basis.T
    )
    if not (
        exact_inertia(ambient_gram) == (1, 11, 0)
        and ambient_gram.det() == -4
        and ambient_action * ambient_gram * ambient_action.T == ambient_gram
        and ambient_action**2 == sp.eye(12)
        and ambient_coinvariant_gram == sp.Matrix([[-4]])
        and weierstrass["family"]["global_order"] == 4
        and weierstrass["brandhorst_hofmann_handoff"]["id"] == "0.4.4.1"
        and weierstrass["brandhorst_hofmann_handoff"][
            "full_NS_action_pairs_integrally_conjugate"
        ]
    ):
        raise AssertionError("the pinned BH row-39 Picard action changed")

    records = []
    for specification in CASES:
        basis = sublattice_basis(specification["even_coordinates_one_based"])
        gram = sp.simplify(basis * ambient_gram * basis.T)
        action = sp.simplify(basis * ambient_action * basis.inv())
        fixed_basis = saturated_row_kernel(action - sp.eye(12))
        fixed_gram = sp.simplify(fixed_basis * gram * fixed_basis.T)
        coinvariant_basis = saturated_row_kernel(gram * fixed_basis.T)
        coinvariant_gram = sp.simplify(
            coinvariant_basis * gram * coinvariant_basis.T
        )
        source_fixed, geometry = source_geometry(specification, graph, d4)
        source_gram = source_lattice_gram(specification, graph, d4)
        discriminant_isometry = _binary_discriminant_isometry(source_gram, gram)
        discriminant_length = len(
            [value for value in smith_diagonal(gram) if value > 1]
        )
        family_census = chamber_census["families"][
            specification["chamber_census_key"]
        ]
        fixed_smith = smith_diagonal(fixed_gram)
        fixed_q = _q_histogram(fixed_gram)
        matching_chamber_classes = [
            entry
            for entry in family_census["classes"]
            if entry["fixed_smith"] == fixed_smith
            and entry["fixed_discriminant_q_histogram_mod_two"] == fixed_q
        ]
        expected_geometry = {
            "dimension": 10,
            "facets": specification["facets"],
            "finite_vertices": specification["finite_vertices"],
            "cusps": specification["cusps"],
            "elliptic_edges": specification["elliptic_edges"],
            "parabolic_audits": specification["parabolic_audits"],
            "bad_audits": 0,
            "classification": "finite_volume_noncompact",
            "complete_audit_green": True,
        }
        source_q = source_fixed["discriminant_q_histogram_mod_two"]
        source_fixed_smith = source_fixed["smith_invariants"]
        if specification["short"] == "H10/16":
            exact_pair_corroboration = h10_16_carrier[
                "direct_row39_divisor_sublattice"
            ]["exact_action_pair_conjugacy_to_archived_C2b"][
                "full_integral_lattice_action_pairs_conjugate"
            ]
            exact_pair_source: str | None = h10_16_carrier["schema"]
        elif specification["short"] == "H10/45":
            row39_h10_45 = next(
                row
                for row in h10_45_carrier["carriers"]
                if row["brandhorst_hofmann"]["id"] == "0.4.4.1"
            )
            exact_pair_corroboration = row39_h10_45["exact_H10_45_action_pair"][
                "full_integral_lattice_action_pairs_conjugate"
            ]
            exact_pair_source = h10_45_carrier["schema"]
        else:
            exact_pair_corroboration = None
            exact_pair_source = None
        checks = {
            "basis_index": abs(int(basis.det())) == specification["index"],
            "full_rank_even_lorentzian": (
                gram.shape == (12, 12)
                and exact_inertia(gram) == (1, 11, 0)
                and all(gram[i, i] % 2 == 0 for i in range(12))
            ),
            "ambient_determinant": gram.det() == specification["determinant"],
            "ambient_smith": smith_diagonal(gram) == specification["smith"],
            "ambient_discriminant_q": _q_histogram(gram) == specification["q"],
            "source_discriminant_form_exactly_isometric": discriminant_isometry
            is not None,
            "source_lattice_identified_by_indefinite_uniqueness": (
                discriminant_isometry is not None
                and gram.rows >= discriminant_length + 2
                and exact_inertia(source_gram) == exact_inertia(gram)
            ),
            "action_integral": all(sp.Rational(value).q == 1 for value in action),
            "one_ambient_action_stabilizes_sublattice": basis * ambient_action
            == action * basis,
            "action_is_involution": action**2 == sp.eye(12),
            "action_isometry": action * gram * action.T == gram,
            "action_character": str(sp.factor(action.charpoly().as_expr()))
            == "(lambda - 1)**11*(lambda + 1)",
            "fixed_rank_signature": (
                fixed_gram.rows == 11 and exact_inertia(fixed_gram) == (1, 10, 0)
            ),
            "fixed_discriminant_q": (
                fixed_q == source_q == specification["fixed_q"]
            ),
            "fixed_smith": fixed_smith == source_fixed_smith,
            "algebraic_coinvariant_is_rootless_minus_four": coinvariant_gram
            == sp.Matrix([[-4]]),
            "primitive_fixed_coinvariant_glue_index_two": abs(
                int(fixed_gram.det() * coinvariant_gram.det() / gram.det())
            )
            == 4,
            "unique_complete_chamber_group_fingerprint": len(
                matching_chamber_classes
            )
            == 1,
            "fingerprint_forces_archived_chamber_class": (
                len(matching_chamber_classes) == 1
                and matching_chamber_classes[0]["class_index"]
                == specification["chamber_class_index"]
            ),
            "literal_vinberg_gate": geometry == expected_geometry,
        }
        if not all(checks.values()):
            raise AssertionError(
                f"a five-chamber carrier gate changed for {specification['short']}: {checks}"
            )
        records.append(
            {
                "polytope": specification["short"],
                "id": specification["id"],
                "action_type": specification["action_type"],
                "target_specific_exact_pair_corroboration": {
                    "source": exact_pair_source,
                    "available": exact_pair_corroboration is not None,
                    "passes": exact_pair_corroboration,
                    "role": (
                        "corroboration only; the common-ample chamber census, "
                        "not a target-specific Weyl conjugator, proves simultaneity"
                    ),
                },
                "sublattice_index_in_UE8A1x2": specification["index"],
                "even_coordinates_one_based_in_pinned_row39_basis": specification[
                    "even_coordinates_one_based"
                ],
                "basis": matrix_payload(basis),
                "gram": matrix_payload(gram),
                "action": matrix_payload(action),
                "fixed_gram": matrix_payload(fixed_gram),
                "coinvariant_gram": matrix_payload(coinvariant_gram),
                "source_lattice_discriminant_form_isometry": discriminant_isometry,
                "complete_chamber_group_fingerprint_match": {
                    "census_key": specification["chamber_census_key"],
                    "matching_classes": matching_chamber_classes,
                    "forced_class_index": specification["chamber_class_index"],
                },
                "geometry": geometry,
                "checks": checks,
                "all_checks_pass": True,
            }
        )

    simultaneity_hypotheses = {
        "one_actual_order_four_K3_automorphism": (
            weierstrass["all_checks_pass"]
            and weierstrass["family"]["global_order"] == 4
            and weierstrass["brandhorst_hofmann_handoff"]["id"] == "0.4.4.1"
            and weierstrass["brandhorst_hofmann_handoff"][
                "full_NS_action_pairs_integrally_conjugate"
            ]
        ),
        "one_common_Picard_involution": all(
            row["checks"]["one_ambient_action_stabilizes_sublattice"]
            for row in records
        ),
        "all_five_divisor_lattices_full_rank_finite_index": all(
            row["checks"]["basis_index"]
            and row["checks"]["full_rank_even_lorentzian"]
            for row in records
        ),
        "all_five_source_lattices_exactly_identified": all(
            row["checks"]["source_discriminant_form_exactly_isometric"]
            and row["checks"]["source_lattice_identified_by_indefinite_uniqueness"]
            for row in records
        ),
        "common_real_positive_cone": all(
            len(row["gram"]) == ambient_gram.rows for row in records
        ),
        "common_algebraic_coinvariant_rootless_minus_four": (
            ambient_coinvariant_gram == sp.Matrix([[-4]])
            and all(
                row["checks"]["algebraic_coinvariant_is_rootless_minus_four"]
                for row in records
            )
        ),
        "rank_eleven_chamber_involution_censuses_complete": chamber_census[
            "all_checks_pass"
        ],
        "each_common_ample_chamber_action_uniquely_forced": all(
            row["checks"]["unique_complete_chamber_group_fingerprint"]
            and row["checks"]["fingerprint_forces_archived_chamber_class"]
            for row in records
        ),
        "all_five_literal_Vinberg_gates": all(
            row["checks"]["literal_vinberg_gate"] for row in records
        ),
    }
    simultaneity_conclusion = all(simultaneity_hypotheses.values())

    certificate = {
        "schema": "coxeter.rank12.five_chamber_simultaneous_k3_carrier.v2",
        "artifact_id": "R12-FIVE-CHAMBER-SIMULTANEOUS-K3-CARRIER",
        "theorem": (
            "The purely non-symplectic order-four UE8A1x2 K3 action carries "
            "five invariant full-rank divisor-subsystem chambers.  Their fixed "
            "sections are exactly the finite-volume noncompact H10/12, H10/15, "
            "H10/16, H10/38, and H10/45 Coxeter polytopes."
        ),
        "carrier": {
            "brandhorst_hofmann_id": "0.4.4.1",
            "generic_NS": "U+E8(-1)+A1(-1)^2",
            "automorphism_order": 4,
            "Picard_action_order": 2,
            "explicit_family_schema": weierstrass["schema"],
            "equation": "y^2=x^3+A(t^2)x+t*C(t^2)",
        },
        "finite_index_divisor_carrier_criterion": {
            "statement": (
                "Let g be one actual finite-order automorphism of a projective K3 "
                "surface and let M_i be finitely many g-stable, full-rank divisor "
                "sublattices of NS.  Average one ample class over <g> to obtain a "
                "single invariant ample h.  The complete M_i-root chamber containing "
                "h is g-stable.  If the complete chamber-group census uniquely "
                "identifies the induced action from its fixed-lattice fingerprint, "
                "then all corresponding fixed Coxeter sections occur simultaneously. "
                "Riemann--Roch makes the h-positive simple roots effective; extra NS "
                "roots may only refine these divisor-subsystem chambers."
            ),
            "common_coinvariant": "<-4>",
            "verified_hypotheses": simultaneity_hypotheses,
            "proof_steps": [
                "Choose any ample class H on the one explicit row-39 K3 and set h=sum_{j=0}^3 (sigma^j)^*H; h is ample and sigma-invariant.",
                "For each full-rank stable M_i, h lies in a unique complete M_i-root chamber, and sigma preserves that chamber because it fixes h.",
                "Conjugate that chamber to the archived complete chamber.  The induced chamber involution has the fixed Smith/q fingerprint computed from sigma|M_i.",
                "The exhaustive chamber-group census has exactly one rank-11 involution class with that fingerprint, forcing H10/12, H10/15, H10/16, H10/38, or H10/45 respectively.",
                "Every h-positive square-minus-two simple root is effective by K3 Riemann--Roch, simultaneously for all five systems.",
            ],
            "common_invariant_ample_class_exists": simultaneity_hypotheses[
                "one_actual_order_four_K3_automorphism"
            ]
            and simultaneity_hypotheses["one_common_Picard_involution"],
            "all_five_divisor_chambers_can_be_oriented_around_that_class": simultaneity_conclusion,
            "simultaneous_carrier_conclusion": simultaneity_conclusion,
            "full_nef_scope": (
                "Only H10/12 is the fixed section of the full NS nef chamber in "
                "this simultaneous carrier.  The other four are effective "
                "finite-index divisor-subsystem chambers."
            ),
        },
        "distinct_hyperbolic_polytope_count": 5,
        "newly_geometric_slots": ["H10/16", "H10/45"],
        "complete_rank_eleven_chamber_involution_census": chamber_census,
        "records": records,
        "source_sha256": {
            "graph": sha256(GRAPH_SOURCE),
            "d4": sha256(D4_SOURCE),
            "bh_lifts": sha256(BH_SOURCE),
            "weierstrass": sha256(WEIERSTRASS_SOURCE),
            "h10_16_carrier": sha256(H10_16_CARRIER_SOURCE),
            "h10_45_carrier": sha256(H10_45_CARRIER_SOURCE),
        },
        "novelty_scope": (
            "new project-level simultaneous algebro-geometric carrier theorem; "
            "worldwide novelty of the abstract H10/16 and H10/45 diagrams remains open"
        ),
        "all_checks_pass": simultaneity_conclusion
        and all(row["all_checks_pass"] for row in records),
    }
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stdout", action="store_true")
    arguments = parser.parse_args()
    certificate = build_certificate()
    payload = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.stdout:
        print(payload, end="")
    else:
        OUTPUT.write_text(payload, encoding="utf-8")
        print(OUTPUT)


if __name__ == "__main__":
    main()
