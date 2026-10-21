#!/usr/bin/env python3
"""Certify the geometric class-2d fixed H10/47 Borcherds chamber."""

from __future__ import annotations

import argparse
import base64
from collections import Counter
import hashlib
import json
from pathlib import Path

import sympy as sp

from scripts.certify_alexeev_high_dim_reflection_examples import (
    SOURCE_SHA256,
    _extract_matrix,
)
from scripts.certify_borcherds_char2_3b_h11_22 import (
    DEFAULT_SOURCE,
    _ambient_action,
    _digest,
    _facet_lattice_records,
    _fixed_lattice,
    _matrix_payload,
    _order,
    _vinberg,
)
from scripts.search_k3_high_rank_involution_folds import _coxeter_audit, _cycles
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _orbit_sum_matrix


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT / "data/enumerations/borcherds_char2_2d_h10_47_certificate.json"
)
SCHEMA = "coxeter.borcherds_char2_2d_h10_47.v1"
PERMUTATION_BASE64 = (
    "AQAfJxUeCR0aBg0oIgoYET8PMmpABCEcDiYIMRcHBQItFgyauG8ZAwvHYsmCIDBsLhsS"
    "rLxynEiRUVKr0KKxEBS5QldQvnhgN4h9Vp1NsqhEOTprWZtLQ3lUz2jIgc29R2cqwXuT"
    "ZmFbrRNTL8yGJXS1Nbtwy5+URli2ZIBKqZl8XSyYjZBupEnKjsKjhIqwhTizZXfDxbSDf"
    "yNVNkzGdqavPYyHpaCqT36nOzNp0aGPPk6Sl3F6xCRBv3M0X0W6wGOLlbeWnilcK4l1b"
    "V7OWjyu"
)
FACE_VECTOR = [47, 540, 3230, 12215, 31212, 54582, 64050, 47760, 19910, 3136]
ENDPOINT_TYPES = Counter({(2, 0): 12050, (1, 1): 7260, (0, 2): 600})
FACET_TYPES = Counter(
    {
        (2, -4, 4): 20,
        (2, -4, 2): 16,
        (1, -4, 4): 6,
        (2, -2, 2): 5,
    }
)


def _permutation() -> tuple[int, ...]:
    result = tuple(base64.b64decode(PERMUTATION_BASE64))
    if len(result) != 210 or set(result) != set(range(210)):
        raise ArithmeticError("archived 2d permutation is malformed")
    return result


def build_certificate(source_path: Path = DEFAULT_SOURCE) -> dict:
    raw = source_path.read_bytes()
    source_hash = hashlib.sha256(raw).hexdigest()
    if source_hash != SOURCE_SHA256:
        raise ArithmeticError("the pinned Alexeev source hash changed")
    wall_gram = _extract_matrix(raw.decode("utf-8"), (22, 2, 0))
    if not (
        wall_gram.shape == (210, 210)
        and wall_gram.rank() == 22
        and exact_inertia(wall_gram) == (1, 21, 188)
    ):
        raise ArithmeticError("the ambient H21/210 Gram changed")

    permutation = _permutation()
    if _order(permutation) != 2:
        raise ArithmeticError("the selected class-2d action is not an involution")
    if any(
        wall_gram[left, right]
        != wall_gram[permutation[left], permutation[right]]
        for left in range(210)
        for right in range(210)
    ):
        raise ArithmeticError("class 2d does not preserve the wall Gram")
    cycles = _cycles(permutation)
    sums = _orbit_sum_matrix(210, cycles)
    orbit_gram = sums.T * wall_gram * sums
    if not (
        len(cycles) == 108
        and sorted(map(len, cycles)) == [1] * 6 + [2] * 102
        and orbit_gram.rank() == 11
        and exact_inertia(orbit_gram) == (1, 10, 97)
    ):
        raise ArithmeticError("the class-2d fixed wall Gram changed")

    facet_indices = [i for i in range(orbit_gram.rows) if orbit_gram[i, i] < 0]
    isotropic = [i for i in range(orbit_gram.rows) if orbit_gram[i, i] == 0]
    timelike = [i for i in range(orbit_gram.rows) if orbit_gram[i, i] > 0]
    if (len(facet_indices), len(isotropic), len(timelike)) != (47, 15, 46):
        raise ArithmeticError("the class-2d causal orbit census changed")

    interior = sp.ones(orbit_gram.rows, 1)
    interior_pairings = orbit_gram * interior
    interior_square = int((interior.T * orbit_gram * interior)[0])
    if interior_square != 137214 or min(interior_pairings) != 198:
        raise ArithmeticError("the class-2d strict interior changed")
    facet_witnesses = []
    for facet in facet_indices:
        witness = interior.copy()
        witness[facet] -= sp.Rational(
            interior_pairings[facet], orbit_gram[facet, facet]
        )
        pairings = orbit_gram * witness
        square = sp.Rational((witness.T * orbit_gram * witness)[0])
        if not (
            square > 0
            and pairings[facet] == 0
            and all(value > 0 for index, value in enumerate(pairings) if index != facet)
        ):
            raise ArithmeticError("a class-2d facet lacks a strict witness")
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
    coxeter = _coxeter_audit(
        facet_gram, range(facet_gram.rows), [cycles[i] for i in facet_indices]
    )
    expected_edges = {
        "2": 490,
        "3": 30,
        "4": 20,
        "6": 0,
        "infinity": 191,
        "ultraparallel": 350,
    }
    if not (
        coxeter["all_pairs_are_exact_coxeter_relations"]
        and coxeter["edge_counts"] == expected_edges
    ):
        raise ArithmeticError("the class-2d Coxeter gate changed")

    basis_gram, action = _ambient_action(wall_gram, permutation, order=2)
    if int(action.trace()) != 0:
        raise ArithmeticError("the class-2d NS trace changed")
    fixed_basis, fixed_gram, fixed_smith = _fixed_lattice(
        basis_gram,
        action,
        expected_rank=11,
        expected_inertia=(1, 10, 0),
        expected_determinant=4096,
        expected_smith=[2] * 10 + [4],
    )
    facet_lattice = _facet_lattice_records(
        wall_gram,
        basis_gram,
        fixed_basis,
        fixed_gram,
        cycles,
        facet_indices,
        expected_histogram=FACET_TYPES,
    )
    facet_payload = _matrix_payload(facet_gram)
    vinberg = _vinberg(
        -facet_gram,
        dimension=10,
        expected_face_vector=FACE_VECTOR,
        expected_ideal_vertices=216,
        expected_endpoint_types=ENDPOINT_TYPES,
    )
    # The shared enumerator was first written for the H11/22 certificate and
    # retains that dimension in its legacy field name.  Keep this archive
    # dimensionally honest without perturbing the already frozen H11 schema.
    vinberg["elliptic_face_vector_rank_1_through_10"] = vinberg.pop(
        "elliptic_face_vector_rank_1_through_11"
    )
    checks = {
        "pinned_source_hash": True,
        "ambient_H21_210_gram": True,
        "order_two_gram_automorphism": True,
        "integral_unimodular_NS_action": True,
        "saturated_fixed_lattice": True,
        "strict_fixed_interior": True,
        "all_47_facets_irredundant": True,
        "all_facet_pairs_exact_coxeter": True,
        "all_facet_reflections_integral": True,
        "vinberg_edge_endpoint_criterion": vinberg["finite_volume"],
        "vinberg_full_parabolic_neighborhoods": (
            vinberg["bad_parabolic_neighborhood_count"] == 0
        ),
    }
    return {
        "schema": SCHEMA,
        "id": "SS2-2D-H10-47",
        "source": {
            "path": str(source_path),
            "sha256": source_hash,
            "matrix_key": "(22,2,0)",
        },
        "action": {
            "atlas_class": "2d",
            "full_diagram_group_class_size": 1008,
            "order": 2,
            "permutation_images_one_based": [value + 1 for value in permutation],
            "cycle_length_histogram": {"1": 6, "2": 102},
            "NS_action_columns": _matrix_payload(action),
            "trace": int(action.trace()),
            "geometric_subgroup_status": "contained in the Dolgachev--Kondo index-two geometric subgroup",
        },
        "fixed_lattice": {
            "definition": "NS(X) intersect Fix_Q(g)",
            "basis_columns_in_ambient_basis": _matrix_payload(fixed_basis),
            "gram": _matrix_payload(fixed_gram),
            "signature_positive_negative": [1, 10],
            "determinant": int(fixed_gram.det()),
            "smith_diagonal_absolute": fixed_smith,
            "rescaled_class": "2 times f_2^10; Scharlau--Kirschmer List10 entry 2",
        },
        "restricted_chamber": {
            "wall_orbit_count": len(cycles),
            "negative_facet_count": len(facet_indices),
            "isotropic_inactive_orbit_count": len(isotropic),
            "timelike_inactive_orbit_count": len(timelike),
            "strict_interior_square": interior_square,
            "strict_interior_minimum_wall_pairing": int(min(interior_pairings)),
            "facet_orbit_indices_one_based": [i + 1 for i in facet_indices],
            "facet_gram": facet_payload,
            "facet_gram_sha256": _digest(facet_payload),
            "facet_witnesses": facet_witnesses,
            "primitive_facet_lattice_records": facet_lattice,
            "primitive_facet_type_histogram": {
                "orbit1_square-4_div4": 6,
                "orbit2_square-2_div2": 5,
                "orbit2_square-4_div2": 16,
                "orbit2_square-4_div4": 20,
            },
            "coxeter_edge_counts": coxeter["edge_counts"],
        },
        "vinberg": vinberg,
        "lineage": {
            "maximal_chamber_facet_count": 13,
            "maximal_chamber_ideal_vertex_count": 3,
            "solomon_euler_characteristic": "-527/120960",
            "maximal_group_euler_characteristic": "-31/11147673600",
            "reflection_subgroup_index": 1566720,
            "covolume": "pi^5*527/3572100",
        },
        "decision": {
            "hyperbolic_dimension": 10,
            "facet_count": 47,
            "finite_volume": True,
            "compact": False,
            "classification": "finite_volume_noncompact_Coxeter_polytope",
            "abstract_novelty": "finite-index subgroup chamber in the classical f_2^10 commensurability class",
            "algebro_geometric_novelty": "priority-open explicit K3 involution fixed-divisor realization",
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
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt H10/47 certificate differs from archive")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "id": payload["id"],
        "elliptic_edges": payload["vinberg"]["elliptic_edge_count"],
        "finite_vertices": payload["vinberg"]["finite_vertex_count"],
        "ideal_vertices": payload["vinberg"]["ideal_vertex_count"],
        "index": payload["lineage"]["reflection_subgroup_index"],
        "facet_gram_sha256": payload["restricted_chamber"]["facet_gram_sha256"],
        "all_checks_pass": payload["all_checks_pass"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
