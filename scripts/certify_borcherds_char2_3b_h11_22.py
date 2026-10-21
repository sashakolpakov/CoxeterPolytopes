#!/usr/bin/env python3
"""Certify the geometric order-3 fixed H11/22 chamber in Borcherds' H21 stock.

The ambient 210-wall Gram is Alexeev's transcription of the Borcherds /
Dolgachev--Kondo characteristic-two chamber.  A fixed representative of the
ATLAS 3b wall-diagram class is replayed exactly.  The verifier derives its
integral Neron--Severi action, saturated fixed lattice, complete restricted
facet set, Coxeter labels, and Vinberg edge/cusp certificate.
"""

from __future__ import annotations

import argparse
import base64
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_alexeev_high_dim_reflection_examples import (
    SOURCE_SHA256,
    _extract_matrix,
)
from scripts.certify_shimada_disc16_d5_h10_30_vinberg_faces import (
    DiagramEnumerator,
)
from scripts.certify_shimada_disc16_d5_h14_fixed_lattice import (
    _smith_decomposition,
)
from scripts.search_k3_high_rank_involution_folds import (
    _coxeter_audit,
    _cycles,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _orbit_sum_matrix


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = Path(
    "/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "borcherds_char2_3b_h11_22_certificate.json"
)
SCHEMA = "coxeter.borcherds_char2_3b_h11_22.v1"

# One representative of the unique order-3, fixed-rank-12 class.  Values are
# the zero-based images of the 210 published walls, packed as unsigned bytes.
PERMUTATION_BASE64 = (
    "HQEyAwQFBiARHwoZCBYAFRAMDTgjGhIJCxgPGy4OAhcxHCI3QRMtJihwXyxGJyEqMAce"
    "koRHdhQlVGdSYSkzP1FzSaQ7lytsgFh9bY+ba3qGYkRNXbI0rEKcXjyWOYUvNltAg3uJ"
    "oYKxi1e6NcNDmD18cSSfVWBjb2mqjnJ+SmWQtTp3VlrB0J1/o3mKtGSZSKs+iGjAoKh4T"
    "MdTuJO3zVy+gYxuxLOVRdFPxWrGqcyllHXJtqKNwlmaTr2vz2bLp1CeS7CRyrmHpq3OvH"
    "S/u8iu"
)
EXPECTED = {
    "atlas_class": "3b",
    "class_size": 672,
    "order": 3,
    "fixed_rank": 12,
    "dimension": 11,
    "wall_orbits": 78,
    "facets": 22,
    "isotropic_orbits": 50,
    "timelike_orbits": 6,
    "elliptic_face_vector": [
        22, 231, 1480, 6320, 18704, 38986,
        57182, 57785, 38140, 14580, 2304,
    ],
    "ideal_vertices": 72,
}


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _digest(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _mask_digest(masks: set[int] | tuple[int, ...]) -> str:
    return _digest([hex(mask) for mask in sorted(masks)])


def _permutation() -> tuple[int, ...]:
    values = tuple(base64.b64decode(PERMUTATION_BASE64))
    if len(values) != 210 or set(values) != set(range(210)):
        raise ArithmeticError("archived 3b permutation is malformed")
    return values


def _order(permutation: tuple[int, ...]) -> int:
    power = tuple(range(len(permutation)))
    for exponent in range(1, 100):
        power = tuple(permutation[power[index]] for index in range(len(power)))
        if power == tuple(range(len(power))):
            return exponent
    raise ArithmeticError("permutation order exceeded the safety bound")


def _ambient_action(
    wall_gram: sp.Matrix, permutation: tuple[int, ...], order: int = 3
) -> tuple[sp.Matrix, sp.Matrix]:
    basis_gram = wall_gram[:22, :22]
    if basis_gram.det() != -4:
        raise ArithmeticError("the first 22 Borcherds walls no longer form NS")
    inverse = basis_gram.inv()
    columns = [
        inverse * wall_gram[:22, permutation[index]] for index in range(22)
    ]
    action = sp.Matrix.hstack(*columns)
    if not (
        all(value.is_Integer for value in action)
        and abs(int(action.det())) == 1
        and action.T * basis_gram * action == basis_gram
        and action**order == sp.eye(22)
    ):
        raise ArithmeticError("the 3b wall permutation did not lift integrally")
    for wall in range(210):
        coordinates = inverse * wall_gram[:22, wall]
        image_coordinates = action * coordinates
        expected = inverse * wall_gram[:22, permutation[wall]]
        if image_coordinates != expected:
            raise ArithmeticError("the NS action fails on a published wall")
    return basis_gram, action


def _fixed_lattice(
    basis_gram: sp.Matrix,
    action: sp.Matrix,
    *,
    expected_rank: int = 12,
    expected_inertia: tuple[int, int, int] = (1, 11, 0),
    expected_determinant: int = -243,
    expected_smith: list[int] | None = None,
) -> tuple[sp.Matrix, sp.Matrix, list[int]]:
    diagonal, _, right = _smith_decomposition(action - sp.eye(22))
    zero_columns = tuple(
        column
        for column in range(22)
        if all(diagonal[row, column] == 0 for row in range(22))
    )
    if len(zero_columns) != expected_rank:
        raise ArithmeticError("the saturated fixed rank changed")
    fixed_basis = right[:, list(zero_columns)]
    if not (
        action * fixed_basis == fixed_basis
        and all(value.is_Integer for value in fixed_basis)
    ):
        raise ArithmeticError("Smith fixed basis failed replay")
    fixed_gram = fixed_basis.T * basis_gram * fixed_basis
    smith = smith_normal_form(fixed_gram, domain=ZZ)
    smith_absolute = [
        abs(int(smith[index, index])) for index in range(fixed_gram.rows)
    ]
    if not (
        exact_inertia(fixed_gram) == expected_inertia
        and int(fixed_gram.det()) == expected_determinant
        and smith_absolute == (expected_smith or [1] * 7 + [3] * 5)
    ):
        raise ArithmeticError("the fixed lattice invariants changed")
    return fixed_basis, fixed_gram, smith_absolute


def _facet_lattice_records(
    wall_gram: sp.Matrix,
    basis_gram: sp.Matrix,
    fixed_basis: sp.Matrix,
    fixed_gram: sp.Matrix,
    cycles: list[list[int]],
    facet_indices: list[int],
    expected_histogram: Counter[tuple[int, int, int]] | None = None,
) -> list[dict[str, Any]]:
    inverse = basis_gram.inv()
    root_coordinates = [inverse * wall_gram[:22, index] for index in range(210)]
    records = []
    for orbit_index in facet_indices:
        ambient = sum(
            (root_coordinates[label - 1] for label in cycles[orbit_index]),
            sp.zeros(22, 1),
        )
        solution_set = sp.linsolve((fixed_basis, ambient))
        if solution_set is sp.EmptySet:
            raise ArithmeticError("a facet orbit sum left the fixed lattice")
        coordinates = sp.Matrix(next(iter(solution_set)))
        if any(value.free_symbols or not value.is_Integer for value in coordinates):
            raise ArithmeticError("a facet direction is not integral")
        common = math.gcd(*(abs(int(value)) for value in coordinates))
        primitive = coordinates / common
        square = int((primitive.T * fixed_gram * primitive)[0])
        divisibility = math.gcd(
            *(abs(int(value)) for value in fixed_gram * primitive)
        )
        reflection = (
            sp.eye(fixed_gram.rows)
            - sp.Rational(2, square) * primitive * primitive.T * fixed_gram
        )
        integral = all(value.is_Integer for value in reflection)
        if not (
            integral
            and reflection.T * fixed_gram * reflection == fixed_gram
            and reflection**2 == sp.eye(fixed_gram.rows)
        ):
            raise ArithmeticError("a facet reflection is not integral")
        records.append(
            {
                "orbit_index_one_based": orbit_index + 1,
                "source_wall_labels": cycles[orbit_index],
                "orbit_size": len(cycles[orbit_index]),
                "orbit_sum_scale_from_primitive": common,
                "primitive_fixed_coordinates": [int(value) for value in primitive],
                "primitive_square": square,
                "primitive_divisibility": divisibility,
                "reflection_integral": integral,
            }
        )
    histogram = Counter(
        (
            row["orbit_size"],
            row["primitive_square"],
            row["primitive_divisibility"],
        )
        for row in records
    )
    if histogram != (expected_histogram or Counter({(1, -2, 1): 12, (3, -6, 3): 10})):
        raise ArithmeticError("the mixed reflective-root facet types changed")
    return records


def _vinberg(
    positive_gram: sp.Matrix,
    *,
    dimension: int = 11,
    expected_face_vector: list[int] | None = None,
    expected_ideal_vertices: int = 72,
    expected_endpoint_types: Counter[tuple[int, int]] | None = None,
) -> dict[str, Any]:
    facet_count = positive_gram.rows
    enumerator = DiagramEnumerator(positive_gram, dimension=dimension)
    finite_by_order = {
        order: enumerator.finite_subdiagrams(order)
        for order in range(1, dimension + 1)
    }
    face_vector = [len(finite_by_order[order]) for order in range(1, dimension + 1)]
    if face_vector != (expected_face_vector or EXPECTED["elliptic_face_vector"]):
        raise ArithmeticError("the exact elliptic face vector changed")
    edges = finite_by_order[dimension - 1]
    finite_vertices = finite_by_order[dimension]

    affine = enumerator.connected_affine_subdiagrams()
    full_neighborhood: dict[int, int] = {}
    for component_mask in affine:
        component = tuple(
            index for index in range(facet_count) if component_mask & (1 << index)
        )
        zero = tuple(
            facet
            for facet in range(facet_count)
            if facet not in component
            and all(positive_gram[facet, node] == 0 for node in component)
        )
        full_neighborhood[component_mask] = component_mask | sum(
            1 << index for index in zero
        )

    ideal_components: dict[int, tuple[int, ...]] = {}
    bad_affine = []
    for component_mask, neighborhood in full_neighborhood.items():
        components = enumerator._components(neighborhood, all_nonorthogonal=True)
        valid = (
            all(component in affine for component in components)
            and sum(component.bit_count() - 1 for component in components)
            == dimension - 1
            and all(full_neighborhood[component] == neighborhood for component in components)
        )
        if valid:
            ideal_components[neighborhood] = components
        else:
            bad_affine.append(component_mask)
    if bad_affine:
        raise ArithmeticError("an affine component has a bad full zero-neighborhood")
    if len(ideal_components) != expected_ideal_vertices:
        raise ArithmeticError("the ideal-vertex count changed")

    finite_by_edge: dict[int, list[int]] = defaultdict(list)
    for vertex_mask in finite_vertices:
        remaining = vertex_mask
        while remaining:
            bit = remaining & -remaining
            remaining ^= bit
            finite_by_edge[vertex_mask ^ bit].append(vertex_mask)

    endpoint_types: Counter[tuple[int, int]] = Counter()
    bad_edges = []
    for edge_mask in edges:
        finite_count = len(finite_by_edge.get(edge_mask, ()))
        ideal_count = 0
        for neighborhood, components in ideal_components.items():
            if edge_mask & ~neighborhood:
                continue
            if all((component & ~edge_mask).bit_count() == 1 for component in components):
                ideal_count += 1
        endpoint_types[(finite_count, ideal_count)] += 1
        if finite_count + ideal_count != 2:
            bad_edges.append(edge_mask)
    expected_types = expected_endpoint_types or Counter(
        {(2, 0): 10965, (1, 1): 3414, (0, 2): 201}
    )
    if bad_edges or endpoint_types != expected_types:
        raise ArithmeticError("Vinberg edge endpoints changed")

    return {
        "criterion": "Vinberg edge-endpoint plus full affine zero-neighborhood",
        "elliptic_face_vector_rank_1_through_11": face_vector,
        "elliptic_edge_count": len(edges),
        "finite_vertex_count": len(finite_vertices),
        "ideal_vertex_count": len(ideal_components),
        "connected_affine_component_count": len(affine),
        "edge_endpoint_type_histogram": {
            f"{finite},{ideal}": count
            for (finite, ideal), count in sorted(endpoint_types.items())
        },
        "bad_edge_count": len(bad_edges),
        "bad_parabolic_neighborhood_count": len(bad_affine),
        "elliptic_edge_masks_sha256": _mask_digest(set(edges)),
        "finite_vertex_masks_sha256": _mask_digest(set(finite_vertices)),
        "ideal_vertex_masks_sha256": _mask_digest(set(ideal_components)),
        "finite_volume": True,
        "compact": False,
        "decision": "finite_volume_noncompact",
    }


def build_certificate(source_path: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    source_raw = source_path.read_bytes()
    source_sha256 = hashlib.sha256(source_raw).hexdigest()
    if source_sha256 != SOURCE_SHA256:
        raise ArithmeticError("the pinned Alexeev ancillary source hash changed")
    wall_gram = _extract_matrix(source_raw.decode("utf-8"), (22, 2, 0))
    if not (
        wall_gram.shape == (210, 210)
        and wall_gram.rank() == 22
        and exact_inertia(wall_gram) == (1, 21, 188)
    ):
        raise ArithmeticError("the ambient H21/210 Gram changed")

    permutation = _permutation()
    if _order(permutation) != 3:
        raise ArithmeticError("the selected wall action is not order three")
    if any(
        wall_gram[left, right]
        != wall_gram[permutation[left], permutation[right]]
        for left in range(210)
        for right in range(210)
    ):
        raise ArithmeticError("the selected action does not preserve all 210 walls")
    cycles = _cycles(permutation)
    orbit_sums = _orbit_sum_matrix(210, cycles)
    orbit_gram = orbit_sums.T * wall_gram * orbit_sums
    if not (
        len(cycles) == EXPECTED["wall_orbits"]
        and sorted(map(len, cycles)) == [1] * 12 + [3] * 66
        and orbit_gram.rank() == EXPECTED["fixed_rank"]
        and exact_inertia(orbit_gram) == (1, 11, 66)
    ):
        raise ArithmeticError("the fixed wall-orbit Gram changed")

    facet_indices = [
        index for index in range(orbit_gram.rows) if orbit_gram[index, index] < 0
    ]
    isotropic = [
        index for index in range(orbit_gram.rows) if orbit_gram[index, index] == 0
    ]
    timelike = [
        index for index in range(orbit_gram.rows) if orbit_gram[index, index] > 0
    ]
    if (
        len(facet_indices), len(isotropic), len(timelike)
    ) != (
        EXPECTED["facets"], EXPECTED["isotropic_orbits"], EXPECTED["timelike_orbits"]
    ):
        raise ArithmeticError("the causal wall-orbit census changed")

    all_ones = sp.ones(orbit_gram.rows, 1)
    ample_pairings = orbit_gram * all_ones
    ample_square = int((all_ones.T * orbit_gram * all_ones)[0])
    if ample_square != 137214 or min(ample_pairings) != 99:
        raise ArithmeticError("the strict fixed interior witness changed")
    facet_witnesses = []
    for facet in facet_indices:
        norm = sp.Rational(orbit_gram[facet, facet])
        witness = all_ones.copy()
        witness[facet] -= sp.Rational(ample_pairings[facet], norm)
        pairings = orbit_gram * witness
        square = sp.Rational((witness.T * orbit_gram * witness)[0])
        if not (
            square > 0
            and pairings[facet] == 0
            and all(value > 0 for index, value in enumerate(pairings) if index != facet)
        ):
            raise ArithmeticError("a restricted wall is not a genuine facet")
        facet_witnesses.append(
            {
                "orbit_index_one_based": facet + 1,
                "witness_square": str(square),
                "minimum_other_wall_pairing": str(
                    min(value for index, value in enumerate(pairings) if index != facet)
                ),
            }
        )

    facet_gram = orbit_gram.extract(facet_indices, facet_indices)
    coxeter = _coxeter_audit(facet_gram, range(facet_gram.rows), [cycles[i] for i in facet_indices])
    expected_edges = {
        "2": 186, "3": 35, "4": 0, "6": 10,
        "infinity": 0, "ultraparallel": 0,
    }
    if not (
        coxeter["all_pairs_are_exact_coxeter_relations"]
        and coxeter["edge_counts"] == expected_edges
    ):
        raise ArithmeticError("the exact Coxeter pair gate changed")

    basis_gram, action = _ambient_action(wall_gram, permutation)
    fixed_basis, fixed_gram, fixed_smith = _fixed_lattice(basis_gram, action)
    facet_lattice = _facet_lattice_records(
        wall_gram, basis_gram, fixed_basis, fixed_gram, cycles, facet_indices
    )
    facet_payload = _matrix_payload(facet_gram)
    vinberg = _vinberg(-facet_gram)

    checks = {
        "pinned_source_hash": True,
        "ambient_H21_210_gram": True,
        "order_three_gram_automorphism": True,
        "integral_unimodular_NS_action": True,
        "saturated_fixed_lattice": True,
        "strict_fixed_interior": True,
        "all_22_facets_irredundant": True,
        "all_facet_pairs_exact_coxeter": True,
        "all_facet_reflections_integral": True,
        "vinberg_edge_endpoint_criterion": vinberg["finite_volume"],
        "vinberg_full_parabolic_neighborhoods": (
            vinberg["bad_parabolic_neighborhood_count"] == 0
        ),
    }
    return {
        "schema": SCHEMA,
        "id": "SS2-3B-H11-22",
        "source": {
            "path": str(source_path),
            "sha256": source_sha256,
            "matrix_key": "(22,2,0)",
            "lineage": "Borcherds H21/210 chamber; Dolgachev--Kondo K3 realization; Alexeev matrix transcription",
        },
        "ambient": {
            "lattice": "U + D20(-1)",
            "rank": 22,
            "wall_count": 210,
            "wall_gram_rank": 22,
            "wall_gram_inertia": [1, 21, 188],
            "basis_gram": _matrix_payload(basis_gram),
        },
        "action": {
            "atlas_class": EXPECTED["atlas_class"],
            "full_diagram_group_class_size": EXPECTED["class_size"],
            "order": 3,
            "permutation_images_one_based": [value + 1 for value in permutation],
            "cycle_length_histogram": {"1": 12, "3": 66},
            "NS_action_columns": _matrix_payload(action),
            "trace": int(action.trace()),
            "geometric_subgroup_status": "contained: every odd-order element lies in the index-two Dolgachev--Kondo geometric subgroup",
        },
        "fixed_lattice": {
            "definition": "NS(X) intersect Fix_Q(g)",
            "basis_columns_in_ambient_basis": _matrix_payload(fixed_basis),
            "gram": _matrix_payload(fixed_gram),
            "signature_positive_negative": [1, 11],
            "determinant": int(fixed_gram.det()),
            "smith_diagonal_absolute": fixed_smith,
            "identified_class": "II_(11,1)(3^-5); Scharlau--Kirschmer List11 entry 27",
        },
        "restricted_chamber": {
            "wall_orbit_count": len(cycles),
            "negative_facet_count": len(facet_indices),
            "isotropic_inactive_orbit_count": len(isotropic),
            "timelike_inactive_orbit_count": len(timelike),
            "strict_interior_square": ample_square,
            "strict_interior_minimum_wall_pairing": int(min(ample_pairings)),
            "facet_orbit_indices_one_based": [index + 1 for index in facet_indices],
            "facet_gram": facet_payload,
            "facet_gram_sha256": _digest(facet_payload),
            "facet_witnesses": facet_witnesses,
            "primitive_facet_lattice_records": facet_lattice,
            "primitive_facet_type_histogram": {
                "orbit1_square-2_div1": 12,
                "orbit3_square-6_div3": 10,
            },
            "coxeter_edge_counts": coxeter["edge_counts"],
        },
        "vinberg": vinberg,
        "decision": {
            "hyperbolic_dimension": 11,
            "facet_count": 22,
            "finite_volume": True,
            "compact": False,
            "classification": "finite_volume_noncompact_Coxeter_polytope",
            "abstract_novelty": "known: Scharlau--Kirschmer List11 entry 27",
            "algebro_geometric_novelty": "priority-open explicit K3 order-3 fixed-divisor realization",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", nargs="?", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.source)
    if args.check_archive:
        archived = json.loads(args.output.read_text(encoding="utf-8"))
        if payload != archived:
            raise ArithmeticError("rebuilt certificate differs from archive")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "id": payload["id"],
        "decision": payload["decision"]["classification"],
        "elliptic_edges": payload["vinberg"]["elliptic_edge_count"],
        "finite_vertices": payload["vinberg"]["finite_vertex_count"],
        "ideal_vertices": payload["vinberg"]["ideal_vertex_count"],
        "facet_gram_sha256": payload["restricted_chamber"]["facet_gram_sha256"],
        "all_checks_pass": payload["all_checks_pass"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
