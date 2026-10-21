#!/usr/bin/env python3
"""Independent exact audit of the rho=15 cyclic K3 divisor-wall folds.

This certifier imports none of the production K3 fold or Vinberg-volume
implementations.  It exhausts the weighted wall-graph automorphism group,
separates the order-three and order-six orbit partitions, and reconstructs
the viable H12 cone by exact rational arithmetic.  It then enumerates every
extreme ray, every Vinberg--CoxIter endpoint, every connected affine full
zero-neighbourhood, and the complete face lattice.

The generic exact helpers are shared with the independent H14 audit; that
helper file and this file are both content-hashed in the certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import sympy as sp

from scripts.independent_certify_k3_rho16_h14_fold import (
    _LiteralDiagramAudit,
    _build_ambient_gram,
    _canonical_sha256,
    _coxeter_label,
    _cycles,
    _discriminant_action_trivial,
    _enumerate_extreme_rays,
    _exact_inertia,
    _facet_audit,
    _is_integral_matrix,
    _matrix_payload,
    _matrix_sha256,
    _orbit_incidence,
    _reflection_matrix,
    _sha256_file,
    _smith_diagonal,
    _weighted_gram_automorphisms,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
PRODUCTION_PATH = ROOT / "data/enumerations/k3_rank15_order3_h12_fold.json"
CERTIFICATE_PATH = ROOT / "data/candidates/k3_rho15_h12_cyclic_fold_audit.json"
HELPER_PATH = ROOT / "scripts/independent_certify_k3_rho16_h14_fold.py"

ORDER_THREE_CYCLES = tuple(
    [(index,) for index in range(1, 14)]
    + [(14, 15, 16), (17, 18, 19), (20,), (21,)]
)

EXPECTED_F_VECTOR = [
    96,
    623,
    2054,
    4433,
    6807,
    7680,
    6444,
    4017,
    1835,
    596,
    130,
    17,
    1,
]


def _permutation_from_cycles(
    size: int, cycles: Sequence[Sequence[int]]
) -> tuple[int, ...]:
    permutation = list(range(size))
    seen: set[int] = set()
    for raw_cycle in cycles:
        cycle = tuple(int(label) - 1 for label in raw_cycle)
        if not cycle or seen.intersection(cycle):
            raise AssertionError("invalid or overlapping cycle list")
        seen.update(cycle)
        for position, source in enumerate(cycle):
            permutation[source] = cycle[(position + 1) % len(cycle)]
    if seen != set(range(size)):
        raise AssertionError("cycle list does not cover every wall")
    return tuple(permutation)


def _compose(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    return tuple(int(left[int(right[index])]) for index in range(len(left)))


def _element_order(permutation: Sequence[int]) -> int:
    identity = tuple(range(len(permutation)))
    power = identity
    for order in range(1, 25):
        power = _compose(permutation, power)
        if power == identity:
            return order
    raise AssertionError("graph automorphism order exceeds group size bound")


def _exact_integer_rank(rows: Sequence[Sequence[int]], columns: int) -> int:
    """Fraction Gaussian rank used only by the independent face census."""

    matrix = [[Fraction(int(entry)) for entry in row] for row in rows]
    rank = 0
    for column in range(columns):
        pivot = next(
            (
                row
                for row in range(rank, len(matrix))
                if matrix[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        for row in range(rank + 1, len(matrix)):
            if not matrix[row][column]:
                continue
            coefficient = matrix[row][column] / pivot_value
            for other in range(column, columns):
                matrix[row][other] -= coefficient * matrix[rank][other]
        rank += 1
    return rank


def _face_lattice_audit(
    vertices: Sequence[dict[str, Any]], pairing_rows: sp.Matrix, dimension: int
) -> dict[str, Any]:
    """Enumerate every face from exact vertex--facet incidence closure."""

    facet_count = pairing_rows.rows
    vertex_count = len(vertices)
    active_masks: list[int] = []
    for vertex in vertices:
        mask = 0
        for label in vertex["active_facets"]:
            mask |= 1 << (int(label) - 1)
        active_masks.append(mask)

    facet_vertex_masks: list[int] = []
    for facet in range(facet_count):
        mask = 0
        for vertex, active in enumerate(active_masks):
            if active >> facet & 1:
                mask |= 1 << vertex
        facet_vertex_masks.append(mask)

    all_vertices = (1 << vertex_count) - 1
    intersection_vertex_masks = [0] * (1 << facet_count)
    intersection_vertex_masks[0] = all_vertices
    for facet_mask in range(1, 1 << facet_count):
        bit = facet_mask & -facet_mask
        facet = bit.bit_length() - 1
        intersection_vertex_masks[facet_mask] = (
            intersection_vertex_masks[facet_mask ^ bit]
            & facet_vertex_masks[facet]
        )
    face_vertex_masks = sorted(
        {mask for mask in intersection_vertex_masks if mask}
    )

    integer_rows = [
        [int(pairing_rows[row, column]) for column in range(pairing_rows.cols)]
        for row in range(pairing_rows.rows)
    ]
    rank_cache: dict[int, int] = {}
    records: list[tuple[int, int, int]] = []
    f_vector = [0] * (dimension + 1)
    for vertex_mask in face_vertex_masks:
        closure = (1 << facet_count) - 1
        remaining = vertex_mask
        while remaining:
            bit = remaining & -remaining
            closure &= active_masks[bit.bit_length() - 1]
            remaining ^= bit
        if closure not in rank_cache:
            rank_cache[closure] = _exact_integer_rank(
                [
                    integer_rows[facet]
                    for facet in range(facet_count)
                    if closure >> facet & 1
                ],
                pairing_rows.cols,
            )
        face_dimension = dimension - rank_cache[closure]
        if not 0 <= face_dimension <= dimension:
            raise AssertionError("face dimension outside the polytope")
        f_vector[face_dimension] += 1
        records.append((vertex_mask, closure, face_dimension))

    return {
        "method": (
            "all 2^17 facet intersections, deduplicated by their nonempty "
            "sets of exact extreme rays; dimensions are 12 minus the exact "
            "rational rank of the common active pairing rows"
        ),
        "facet_subset_count": 1 << facet_count,
        "nonempty_face_count": len(face_vertex_masks),
        "distinct_active_closure_count": len(rank_cache),
        "f_vector_vertices_through_polytope": f_vector,
        "alternating_euler_sum": sum(
            (-1) ** index * value for index, value in enumerate(f_vector)
        ),
        "canonical_full_face_lattice_sha256": _canonical_sha256(records),
    }


def _production_agreement(
    production: dict[str, Any],
    orbit_gram: sp.Matrix,
    vertices: Sequence[dict[str, Any]],
    endpoint_audit: dict[str, Any],
    relation_counts: Counter[str],
) -> dict[str, bool]:
    independent_vertex_types = {
        (vertex["type"], tuple(vertex["active_facets"])) for vertex in vertices
    }
    production_vertex_types = {
        (vertex["type"], tuple(vertex["active_facets"]))
        for vertex in production["vertices"]
    }
    production_counts = {
        key: value
        for key, value in production["coxeter"]["edge_counts"].items()
        if value
    }
    return {
        "orbit_gram_equal": production["orbit_gram"] == _matrix_payload(orbit_gram),
        "extreme_ray_type_and_active_sets_equal": (
            production_vertex_types == independent_vertex_types
        ),
        "finite_ideal_hyperideal_counts_equal": (
            production["finite_vertex_count"],
            production["ideal_vertex_count"],
            production["hyperideal_vertex_count"],
        )
        == (
            endpoint_audit["finite_vertex_count"],
            endpoint_audit["ideal_vertex_count"],
            0,
        ),
        "coxeter_relation_counts_equal": production_counts
        == dict(sorted(relation_counts.items())),
        "literal_edge_and_endpoint_counts_equal": (
            production["literal_coxiter_prop_1_10"][
                "spherical_rank_11_edge_count"
            ],
            production["literal_coxiter_prop_1_10"][
                "spherical_rank_12_vertex_count"
            ],
            production["literal_coxiter_prop_1_10"][
                "parabolic_rank_11_vertex_count"
            ],
        )
        == (
            endpoint_audit["elliptic_edge_count"],
            endpoint_audit["finite_vertex_count"],
            endpoint_audit["ideal_vertex_count"],
        ),
    }


def build_certificate() -> dict[str, Any]:
    source = json.loads(SOURCE_PATH.read_text())
    production = json.loads(PRODUCTION_PATH.read_text())
    case = next(case for case in source["cases"] if int(case["rho"]) == 15)
    ambient = _build_ambient_gram(case)
    basis = ambient[:15, :15]
    wall_coordinates = basis.inv() * ambient[:15, :]
    if not _is_integral_matrix(wall_coordinates):
        raise AssertionError("published walls are not integral in the basis")

    automorphisms = _weighted_gram_automorphisms(ambient)
    order_histogram = Counter(_element_order(element) for element in automorphisms)
    order_three = tuple(
        element for element in automorphisms if _element_order(element) == 3
    )
    order_six = tuple(
        element for element in automorphisms if _element_order(element) == 6
    )

    proposed = _permutation_from_cycles(21, ORDER_THREE_CYCLES)
    if proposed not in order_three:
        raise AssertionError("advertised order-three action is absent")
    proposed_inverse = _compose(proposed, proposed)
    order_three_partitions = {
        frozenset(frozenset(cycle) for cycle in _cycles(element))
        for element in order_three
    }

    order_six_records: list[dict[str, Any]] = []
    for element in order_six:
        cycles = _cycles(element)
        incidence = _orbit_incidence(21, cycles)
        gram = incidence.T * ambient * incidence
        order_six_records.append(
            {
                "cycles": [list(cycle) for cycle in cycles],
                "orbit_count": len(cycles),
                "orbit_gram_rank": gram.rank(),
                "orbit_gram_inertia": list(_exact_inertia(gram)),
                "orbit_normal_squares": [int(entry) for entry in gram.diagonal()],
                "valid_spacelike_Coxeter_wall_partition": all(
                    entry < 0 for entry in gram.diagonal()
                ),
            }
        )

    incidence = _orbit_incidence(21, ORDER_THREE_CYCLES)
    orbit_gram = incidence.T * ambient * incidence
    proposed_action = wall_coordinates[:, list(proposed[:15])]
    orbit_coordinates = wall_coordinates * incidence

    # The first 13 basis walls are fixed pointwise and form a coordinate
    # direct summand, hence a saturated basis of the rank-13 fixed lattice.
    fixed_gram = basis[:13, :13]
    fixed_normals = orbit_coordinates[:13, :]
    pairing_rows = fixed_normals.T * fixed_gram
    if pairing_rows != orbit_gram[:, :13]:
        raise AssertionError("fixed pairings disagree with the orbit Gram")

    interior = fixed_gram.inv() * sp.Matrix([2] * 13)
    interior_pairings = pairing_rows * interior
    interior_square = int((interior.T * fixed_gram * interior)[0])

    vertices = _enumerate_extreme_rays(fixed_gram, pairing_rows)
    finite_vertices = tuple(v for v in vertices if v["type"] == "finite")
    ideal_vertices = tuple(v for v in vertices if v["type"] == "ideal")
    hyperideal_vertices = tuple(v for v in vertices if v["type"] == "hyperideal")
    facets = _facet_audit(vertices, fixed_gram, pairing_rows)

    reflections: list[dict[str, Any]] = []
    for index in range(fixed_normals.cols):
        normal = fixed_normals[:, index]
        matrix = _reflection_matrix(normal, fixed_gram)
        reflections.append(
            {
                "facet": index + 1,
                "ambient_wall_orbit": list(ORDER_THREE_CYCLES[index]),
                "normal_coordinates": [int(entry) for entry in normal],
                "normal_square": int((normal.T * fixed_gram * normal)[0]),
                "matrix": _matrix_payload(matrix),
                "matrix_sha256": _matrix_sha256(matrix),
                "integral": _is_integral_matrix(matrix),
                "involution": matrix * matrix == sp.eye(13),
                "isometry": matrix.T * fixed_gram * matrix == fixed_gram,
            }
        )

    positive_wall_gram = -orbit_gram
    relation_counts: Counter[str] = Counter()
    relation_edges: list[dict[str, Any]] = []
    for left in range(17):
        for right in range(left + 1, 17):
            label = _coxeter_label(
                int(positive_wall_gram[left, left]),
                int(positive_wall_gram[right, right]),
                int(positive_wall_gram[left, right]),
            )
            relation_counts[label] += 1
            if label != "2":
                relation_edges.append(
                    {
                        "facets": [left + 1, right + 1],
                        "label": label,
                        "normalized_square_ratio": str(
                            sp.Rational(positive_wall_gram[left, right]) ** 2
                            / (
                                positive_wall_gram[left, left]
                                * positive_wall_gram[right, right]
                            )
                        ),
                    }
                )

    diagram = _LiteralDiagramAudit(positive_wall_gram)
    endpoint_audit = diagram.literal_endpoint_audit(12)
    endpoint_finite_sets = {
        tuple(facets)
        for facets in endpoint_audit["elliptic_rank_n_finite_vertices"]
    }
    endpoint_ideal_sets = {
        tuple(vertex["facets"])
        for vertex in endpoint_audit[
            "parabolic_rank_n_minus_one_ideal_vertices"
        ]
    }
    extreme_finite_sets = {
        tuple(vertex["active_facets"]) for vertex in finite_vertices
    }
    extreme_ideal_sets = {
        tuple(vertex["active_facets"]) for vertex in ideal_vertices
    }

    face_lattice = _face_lattice_audit(vertices, pairing_rows, 12)
    production_agreement = _production_agreement(
        production, orbit_gram, vertices, endpoint_audit, relation_counts
    )

    identity = sp.eye(15)
    action_discriminant_trivial = _discriminant_action_trivial(
        proposed_action, basis
    )
    checks = {
        "ambient_complete_wall_lattice_exact": (
            ambient.rank() == 15
            and _exact_inertia(ambient) == (1, 14, 6)
            and int(basis.det()) == 8
            and _smith_diagonal(basis) == [1] * 12 + [2, 2, 2]
        ),
        "full_weighted_graph_automorphism_census": (
            len(automorphisms) == 12
            and order_histogram == Counter({1: 1, 2: 7, 3: 2, 6: 2})
        ),
        "inverse_order_three_pair_has_one_orbit_partition": (
            len(order_three) == 2
            and proposed_inverse in order_three
            and len(order_three_partitions) == 1
        ),
        "order_six_partition_is_not_H12_or_all_spacelike": (
            len(order_six_records) == 2
            and all(record["orbit_count"] == 9 for record in order_six_records)
            and all(record["orbit_gram_rank"] == 7 for record in order_six_records)
            and all(
                record["orbit_gram_inertia"] == [1, 6, 2]
                for record in order_six_records
            )
            and all(
                record["orbit_normal_squares"]
                == [-4] * 6 + [-2, 0, 8]
                for record in order_six_records
            )
            and not any(
                record["valid_spacelike_Coxeter_wall_partition"]
                for record in order_six_records
            )
        ),
        "order_three_NS_action_integral_order_three_isometry": (
            _is_integral_matrix(proposed_action)
            and proposed_action**3 == identity
            and proposed_action != identity
            and proposed_action.T * basis * proposed_action == basis
            and proposed_action * orbit_coordinates == orbit_coordinates
        ),
        "fixed_lattice_saturated_rank13_H12": (
            (proposed_action - identity).rank() == 2
            and orbit_coordinates.rank() == 13
            and _exact_inertia(fixed_gram) == (1, 12, 0)
            and int(fixed_gram.det()) == 6
            and _smith_diagonal(fixed_gram) == [1] * 12 + [6]
        ),
        "orbit_gram_rank_inertia_and_norms": (
            orbit_gram.rank() == 13
            and _exact_inertia(orbit_gram) == (1, 12, 4)
            and list(orbit_gram.diagonal())
            == [-2] * 13 + [-6, -6, -2, -2]
        ),
        "all_17_complete_fixed_nef_inequalities": (
            len(ORDER_THREE_CYCLES) == 17
            and {
                label for orbit in ORDER_THREE_CYCLES for label in orbit
            }
            == set(range(1, 22))
            and pairing_rows.rank() == 13
        ),
        "strict_integral_timelike_interior": (
            _is_integral_matrix(interior)
            and interior_square == 1110
            and [int(entry) for entry in interior_pairings]
            == [2] * 13 + [54, 54, 58, 58]
        ),
        "all_17_facets_irredundant_with_relative_interiors": all(
            facet["irredundant_with_strict_relative_interior"]
            for facet in facets
        ),
        "all_17_folded_reflections_integral": all(
            reflection["integral"]
            and reflection["involution"]
            and reflection["isometry"]
            for reflection in reflections
        ),
        "exact_coxeter_relation_census": relation_counts
        == Counter(
            {
                "2": 116,
                "3": 12,
                "6": 2,
                "infinity": 3,
                "ultraparallel": 3,
            }
        ),
        "extreme_rays_93_finite_3_ideal_no_hyperideal": (
            len(finite_vertices) == 93
            and len(ideal_vertices) == 3
            and not hyperideal_vertices
        ),
        "literal_623_edges_each_have_two_endpoints": (
            endpoint_audit["elliptic_edge_count"] == 623
            and endpoint_audit["extension_histogram"] == {"2": 623}
            and endpoint_audit["endpoint_type_histogram"]
            == {"0,2": 3, "1,1": 124, "2,0": 496}
            and not endpoint_audit["bad_edges"]
        ),
        "literal_endpoints_equal_extreme_rays": (
            endpoint_finite_sets == extreme_finite_sets
            and endpoint_ideal_sets == extreme_ideal_sets
        ),
        "all_affine_zero_neighbourhoods_are_full_cusps": (
            endpoint_audit["connected_affine_component_count"] == 8
            and endpoint_audit["distinct_full_zero_neighbourhood_count"] == 3
            and endpoint_audit["full_zero_neighbourhoods_match_ideal_vertices"]
            and all(
                audit["is_full_rank_dimension_minus_one_parabolic"]
                for audit in endpoint_audit[
                    "connected_affine_zero_neighbourhood_audits"
                ]
            )
        ),
        "finite_volume_noncompact": endpoint_audit[
            "finite_volume_noncompact"
        ],
        "complete_face_lattice_f_vector": (
            face_lattice["facet_subset_count"] == 131072
            and face_lattice["nonempty_face_count"] == sum(EXPECTED_F_VECTOR)
            and face_lattice["f_vector_vertices_through_polytope"]
            == EXPECTED_F_VECTOR
            and face_lattice["alternating_euler_sum"] == 1
        ),
        "independent_and_production_certificates_agree": all(
            production_agreement.values()
        ),
        "order_three_discriminant_action_is_nontrivial": (
            not action_discriminant_trivial and proposed_action**3 == identity
        ),
    }
    if not all(checks.values()):
        raise AssertionError(f"rho15 H12 independent audit failed: {checks}")

    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "artifact_id": "independent_k3_rho15_h12_cyclic_fold_audit",
        "certificate_kind": (
            "exact_divisor_wall_cyclic_fold_finite_volume_and_face_lattice"
        ),
        "source": {
            "path": str(SOURCE_PATH.relative_to(ROOT)),
            "sha256": _sha256_file(SOURCE_PATH),
            "ambient_lattice": case["lattice"],
            "published_model": case["published_model"],
            "source_locator": case["source_locator"],
            "complete_wall_count": int(case["wall_count"]),
        },
        "weighted_wall_graph_cyclic_census": {
            "automorphism_count": len(automorphisms),
            "group_identification": "S3 x C2",
            "element_order_histogram": {
                str(order): count for order, count in sorted(order_histogram.items())
            },
            "order_three_elements": [
                [list(cycle) for cycle in _cycles(element)]
                for element in order_three
            ],
            "order_three_elements_are_inverse_with_same_partition": True,
            "order_six_elements": order_six_records,
            "order_six_conclusion": (
                "The order-six cyclic partitions have only nine orbits and "
                "fixed H6 signature; their orbit sums include one isotropic "
                "and one timelike vector, so they do not define the H12 fold."
            ),
        },
        "order_three_H12_fold": {
            "wall_orbits": [list(cycle) for cycle in ORDER_THREE_CYCLES],
            "fixed_lattice_rank": 13,
            "fixed_lattice_inertia": list(_exact_inertia(fixed_gram)),
            "fixed_lattice_determinant": int(fixed_gram.det()),
            "fixed_lattice_smith_diagonal": _smith_diagonal(fixed_gram),
            "likely_abstract_lattice_identification": (
                "U + E8(-1) + A2(-1) + A1(-1) (rank/signature/determinant "
                "and displayed basis Gram agree; no uniqueness claim needed)"
            ),
            "hyperbolic_dimension": 12,
            "facet_count": 17,
            "orbit_gram_rank": orbit_gram.rank(),
            "orbit_gram_inertia": list(_exact_inertia(orbit_gram)),
            "orbit_normal_squares": [int(entry) for entry in orbit_gram.diagonal()],
            "complete_inequality_argument": (
                "The source supplies all 21 nef inequalities. On a fixed "
                "class the three inequalities in either 3-cycle are equal "
                "and collapse to the orbit-sum inequality; the 15 singleton "
                "inequalities remain. The 17 orbits cover every source wall."
            ),
            "strict_interior": {
                "basis_coordinates": [int(entry) for entry in interior],
                "square": interior_square,
                "wall_pairings": [int(entry) for entry in interior_pairings],
            },
            "facet_audits": list(facets),
            "folded_reflections": reflections,
            "coxeter_relation_counts": dict(sorted(relation_counts.items())),
            "nonright_coxeter_edges": relation_edges,
            "vertices": list(vertices),
            "finite_vertex_count": len(finite_vertices),
            "ideal_vertex_count": len(ideal_vertices),
            "hyperideal_vertex_count": len(hyperideal_vertices),
        },
        "literal_vinberg_coxiter_prop_1_10": endpoint_audit,
        "complete_face_lattice": face_lattice,
        "production_certificate_crosscheck": {
            "path": str(PRODUCTION_PATH.relative_to(ROOT)),
            "sha256": _sha256_file(PRODUCTION_PATH),
            "agreement": production_agreement,
        },
        "surface_realization_boundary": {
            "action_trivial_on_NS_discriminant": action_discriminant_trivial,
            "status": (
                "not the published general-surface automorphism action; "
                "special-period Torelli realization is not needed or claimed"
            ),
            "role": (
                "Surface-automorphism realization is not a gate: the complete "
                "divisor-wall configuration plus its integral graph/lattice "
                "symmetry already defines the exact fixed-cone polytope."
            ),
        },
        "novelty_audit": {
            "exact_new_example_claim": False,
            "status": "strong_novelty_candidate_but_literature_exhaustion_incomplete",
            "coxiter_repository": {
                "url": "https://github.com/rgugliel/CoxIter/tree/master/graphs",
                "audited_commit": "b800be48d8240fbabd2d1c9542e61054584ad7fe",
                "dimension_12_file_count": 11,
                "dimension_12_facet_counts": {
                    "13": 6,
                    "14": 3,
                    "16": 1,
                    "18": 1,
                    "17": 0,
                },
                "conclusion": (
                    "No distributed CoxIter H12 diagram has 17 facets, so "
                    "none is isomorphic to this diagram."
                ),
            },
            "warning": (
                "CoxIter's distributed examples are not an exhaustive "
                "classification of finite-volume H12 Coxeter polytopes. "
                "Absence there is evidence, not a proof of global novelty."
            ),
        },
        "matrices": {
            "ambient_curve_gram_21": _matrix_payload(ambient),
            "published_basis_gram_15": _matrix_payload(basis),
            "order_three_NS_action_15": _matrix_payload(proposed_action),
            "orbit_incidence_21x17": _matrix_payload(incidence),
            "orbit_gram_17": _matrix_payload(orbit_gram),
            "fixed_lattice_gram_13": _matrix_payload(fixed_gram),
            "folded_normal_coordinates_13x17": _matrix_payload(fixed_normals),
            "positive_normal_coxeter_gram_17": _matrix_payload(
                positive_wall_gram
            ),
        },
        "matrix_sha256": {
            "ambient_curve_gram_21": _matrix_sha256(ambient),
            "published_basis_gram_15": _matrix_sha256(basis),
            "order_three_NS_action_15": _matrix_sha256(proposed_action),
            "orbit_gram_17": _matrix_sha256(orbit_gram),
            "fixed_lattice_gram_13": _matrix_sha256(fixed_gram),
            "positive_normal_coxeter_gram_17": _matrix_sha256(
                positive_wall_gram
            ),
            "all_extreme_rays": _canonical_sha256(vertices),
            "all_edge_extensions": _canonical_sha256(
                endpoint_audit["edge_extensions"]
            ),
        },
        "checks": checks,
        "source_sha256": {
            "certifier": _sha256_file(Path(__file__).resolve()),
            "independent_generic_helper": _sha256_file(HELPER_PATH),
            "k3_high_rank_nef_chambers": _sha256_file(SOURCE_PATH),
            "production_crosscheck": _sha256_file(PRODUCTION_PATH),
        },
        "commands": [
            "PYTHONPATH=src:. python3 scripts/independent_certify_k3_rho15_h12_cyclic_fold.py --write-certificate",
            "PYTHONPATH=src:. python3 scripts/independent_certify_k3_rho15_h12_cyclic_fold.py --verify-certificate",
            "PYTHONPATH=src:. pytest -q tests/test_k3_rho15_h12_cyclic_fold_audit.py",
        ],
        "conclusion": {
            "abstract_hyperbolic_coxeter_polytope": True,
            "dimension": 12,
            "finite_volume": True,
            "compact": False,
            "finite_vertices": 93,
            "ideal_vertices": 3,
            "facets": 17,
            "f_vector": EXPECTED_F_VECTOR,
            "positive_example": True,
            "global_novelty_proved": False,
            "summary": (
                "The order-three divisor-wall fold is an exact finite-volume "
                "noncompact H12 Coxeter polytope. The order-six partitions "
                "are invalid for H12. No match occurs in the distributed "
                "CoxIter corpus, but global literature novelty remains open."
            ),
        },
    }
    payload["certificate_payload_sha256"] = _canonical_sha256(payload)
    return payload


def verify_certificate(payload: dict[str, Any]) -> dict[str, bool]:
    body = dict(payload)
    digest = body.pop("certificate_payload_sha256", None)
    fresh = build_certificate()
    checks = {
        "payload_digest_matches": digest == _canonical_sha256(body),
        "source_hashes_match": payload.get("source_sha256")
        == fresh.get("source_sha256"),
        "full_independent_recomputation_matches": payload == fresh,
        "all_exact_checks_are_true": all(payload["checks"].values()),
    }
    if not all(checks.values()):
        raise AssertionError(f"rho15 H12 certificate verification failed: {checks}")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-certificate", action="store_true")
    parser.add_argument("--verify-certificate", action="store_true")
    arguments = parser.parse_args()
    if arguments.write_certificate:
        payload = build_certificate()
        CERTIFICATE_PATH.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )
        print(f"wrote {CERTIFICATE_PATH}")
        print(payload["certificate_payload_sha256"])
        return
    payload = json.loads(CERTIFICATE_PATH.read_text())
    print(json.dumps(verify_certificate(payload), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
