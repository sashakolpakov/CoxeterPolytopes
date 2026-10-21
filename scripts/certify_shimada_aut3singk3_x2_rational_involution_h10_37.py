#!/usr/bin/env python3
"""Certify the X_2 rational-involution fixed chamber in H^10.

The denominator is Shimada's theorem-complete 478-wall induced chamber
``D^(0)`` and all 720 matrices of ``Aut(X_2,a_2)``.  The certificate does
four separate things:

* exhausts every proper subgroup which can have fixed rank at least 11;
* restricts all 478 inequalities to both involution classes;
* applies Vinberg's literal edge/endpoint and full-parabolic-neighbourhood
  tests to the unique Coxeter survivor; and
* identifies its fixed lattice and its finite-index position in the known
  full reflection group of the 2-elementary lattice (11,9,1).

Matrices in Shimada's supplement act on row vectors on the right.  Fixed
lattice coordinate matrices below use the usual column-vector convention
only when constructing reflection matrices.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Any, Iterable, Sequence

import numpy as np
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from scripts.certify_alexeev_remaining_fixed_sections_lineage import _solomon
from scripts.certify_alexeev_remaining_stable_overlattices import (
    _q_histogram,
    _stable_subgroups,
)
from scripts.certify_borcherds_char2_high_fixed_sections_vinberg import (
    _parse_output,
)
from scripts.certify_higher_dimensional_rational_lineage import _e8_cartan
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
STOCK = ROOT / "data/enumerations/shimada_aut3singk3_x2_478_wall_stock.json"
K3_FULL = ROOT / "data/enumerations/k3_rank18_19_all_subgroup_completion.json"
CPP_SOURCE = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
OUTPUT = (
    ROOT
    / "data/enumerations/shimada_aut3singk3_x2_rational_involution_h10_37.json"
)

SCHEMA = "coxeter.shimada.aut3singk3.x2.rational_involution_h10_37.v1"
STOCK_SHA256 = "10d1c767d21d7194c6ccb603c5dff29520201954b4a07fd3c5f214d452ea2f74"
K3_FULL_SHA256 = "148ce106b968709c94b31ac8bed13b17269861ea6051eeba68b5bc8b89e7e5c4"
CPP_SHA256 = "45174f5461d28affc1fdb4cb55047b7a1581037d307e32800651badb3c2df045"

IDENTITY_INDEX = 599
RATIONAL_REPRESENTATIVE = 596
SYMPLECTIC_REPRESENTATIVE = 470

EXPECTED_FIXED_RANK_HISTOGRAM = {
    3: 324,
    4: 144,
    6: 90,
    8: 80,
    11: 36,
    12: 45,
    20: 1,
}

EXPECTED_VINBERG = {
    "dimension": 10,
    "facet_count": 37,
    "edge_search_node_count": 254182,
    "all_finite_subset_search_node_count": 306198,
    "ordering_dependent_edge_search_prefix_counts": [
        37, 470, 3115, 12422, 32580, 56242, 67559, 54170, 24350, 3236
    ],
    "finite_special_subset_counts_by_rank": [
        37, 490, 3415, 14490, 39842, 72352, 85585, 62400, 24350, 3236
    ],
    "connected_affine_component_count": 813,
    "ideal_vertex_count": 261,
    "bad_parabolic_neighborhood_count": 0,
    "bad_edge_count": 0,
    "edge_endpoint_type_histogram": {"0,2": 2360, "1,1": 11620, "2,0": 10370},
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _matrix(payload: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _fraction_text(value: int | Fraction | sp.Rational) -> str:
    value = sp.Rational(value)
    return str(int(value)) if value.q == 1 else str(value)


def _payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [_fraction_text(matrix[row, col]) for col in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _integer_payload(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.is_Integer is not True for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return [[int(value) for value in row] for row in matrix.tolist()]


def _rank(matrix: sp.Matrix) -> int:
    return DomainMatrix.from_list_sympy(
        matrix.rows, matrix.cols, matrix.tolist()
    ).rank()


def _smith(matrix: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(matrix, domain=ZZ).diagonal()
    return [abs(int(value)) for value in diagonal]


def _primitive_integer_row(vector: sp.Matrix) -> tuple[int, ...]:
    values = [sp.Rational(value) for value in vector]
    denominator = math.lcm(*(int(value.q) for value in values))
    integers = [int(value * denominator) for value in values]
    content = math.gcd(*(abs(value) for value in integers))
    if content == 0:
        raise ArithmeticError("zero restricted wall")
    # The sign is the oriented half-space sign and must be preserved.
    return tuple(value // content for value in integers)


def _matrix_key(matrix: sp.Matrix) -> tuple[int, ...]:
    if any(value.is_Integer is not True for value in matrix):
        raise ArithmeticError("nonintegral group matrix")
    return tuple(int(value) for value in matrix)


def _smith_decomposition(matrix: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    diagonal, left, right = smith_normal_decomp(
        DomainMatrix.from_Matrix(matrix).convert_to(ZZ)
    )
    return diagonal.to_Matrix(), left.to_Matrix(), right.to_Matrix()


def _saturated_fixed_basis(
    generators: Sequence[sp.Matrix], ambient_rank: int
) -> tuple[sp.Matrix, dict[str, Any]]:
    identity = sp.eye(ambient_rank)
    equations = generators[0].T - identity
    for generator in generators[1:]:
        equations = equations.col_join(generator.T - identity)
    diagonal, _, right = _smith_decomposition(equations)
    equation_rank = _rank(equations)
    basis = right[:, equation_rank:].T
    reduced = basis.lll(delta=sp.Rational(3, 4))
    _, pivots = basis.rref()
    pivot_columns = list(pivots)
    change = reduced[:, pivot_columns] * basis[:, pivot_columns].inv()
    if not (
        change * basis == reduced
        and abs(int(change.det())) == 1
        and all(reduced * generator == reduced for generator in generators)
    ):
        raise ArithmeticError("failed to recover the saturated fixed kernel")
    return reduced, {
        "equation_matrix_rank": equation_rank,
        "nonzero_smith_diagonal": [
            int(diagonal[i, i])
            for i in range(min(diagonal.shape))
            if diagonal[i, i]
        ],
        "kernel_rank": reduced.rows,
        "lll_basis_change_determinant": int(change.det()),
        "saturated_smith_kernel": True,
    }


def _load_denominator() -> tuple[
    dict[str, Any], sp.Matrix, sp.Matrix, list[sp.Matrix], list[sp.Matrix]
]:
    if _sha256(STOCK) != STOCK_SHA256:
        raise ArithmeticError("the frozen complete X_2 wall stock changed")
    stock = json.loads(STOCK.read_text(encoding="utf-8"))
    if not (
        stock["all_checks_pass"] is True
        and stock["complete_induced_chamber"]["wall_count"] == 478
        and stock["ample_stabilizer"]["order"] == 720
    ):
        raise ArithmeticError("the X_2 source denominator is incomplete")
    gram = _matrix(stock["surface"]["published_neron_severi_gram"])
    ample = sp.Matrix(1, 20, stock["surface"]["published_ample_vector_a2"])
    walls = [
        sp.Matrix(1, 20, wall)
        for orbit in stock["complete_induced_chamber"]["wall_orbits"]
        for wall in orbit["primitive_outward_wall_vectors_in_published_NS_basis"]
    ]
    group = [
        _matrix(matrix)
        for matrix in stock["ample_stabilizer"][
            "all_integral_20_by_20_matrices_in_supplement_order"
        ]
    ]
    return stock, gram, ample, walls, group


def _common_fixed_rank(generators: Sequence[sp.Matrix]) -> int:
    equations = generators[0].T - sp.eye(20)
    for generator in generators[1:]:
        equations = equations.col_join(generator.T - sp.eye(20))
    return 20 - _rank(equations)


def _subgroup_denominator(group: Sequence[sp.Matrix]) -> dict[str, Any]:
    identity = sp.eye(20)
    fixed_ranks = [20 - _rank(matrix.T - identity) for matrix in group]
    histogram = Counter(fixed_ranks)
    if dict(sorted(histogram.items())) != EXPECTED_FIXED_RANK_HISTOGRAM:
        raise ArithmeticError("the 720-element fixed-rank census changed")

    arrays = np.asarray(
        [[[int(value) for value in row] for row in matrix.tolist()] for matrix in group],
        dtype=np.int64,
    )
    maximum_entry = int(np.max(np.abs(arrays)))
    product_bound = 20 * maximum_entry**2
    if product_bound >= 2**63:
        raise ArithmeticError("int64 group multiplication is not exact")
    index_of = {
        tuple(int(value) for value in arrays[index].ravel()): index
        for index in range(len(group))
    }
    if len(index_of) != 720:
        raise ArithmeticError("duplicate stabilizer matrices")

    high = [
        index
        for index, rank in enumerate(fixed_ranks)
        if rank >= 11 and index != IDENTITY_INDEX
    ]
    if len(high) != 81:
        raise ArithmeticError("high involution denominator changed")
    identity_array = np.eye(20, dtype=np.int64)
    if not all(np.array_equal(arrays[index] @ arrays[index], identity_array) for index in high):
        raise ArithmeticError("a high nonidentity element is not an involution")

    high_set = set(high)
    v4_groups: set[frozenset[int]] = set()
    for position, left in enumerate(high):
        for right in high[position + 1 :]:
            product = arrays[left] @ arrays[right]
            product_index = index_of.get(tuple(int(value) for value in product.ravel()))
            if product_index is None:
                raise ArithmeticError("the published 720 matrices are not closed")
            if product_index not in high_set:
                continue
            if not np.array_equal(product, arrays[right] @ arrays[left]):
                raise ArithmeticError("an involutive high product did not commute")
            subgroup = frozenset((IDENTITY_INDEX, left, right, product_index))
            if len(subgroup) == 4:
                v4_groups.add(subgroup)

    v4_rows = []
    for subgroup in sorted(v4_groups, key=lambda row: tuple(sorted(row))):
        nonidentity = sorted(subgroup - {IDENTITY_INDEX})
        rank = _common_fixed_rank([group[index] for index in nonidentity])
        v4_rows.append(
            {
                "elements_zero_based": sorted(subgroup),
                "nonidentity_element_fixed_ranks": sorted(
                    fixed_ranks[index] for index in nonidentity
                ),
                "common_fixed_rank": rank,
            }
        )
    v4_rank_histogram = Counter(row["common_fixed_rank"] for row in v4_rows)
    v4_pattern_histogram = Counter(
        tuple(row["nonidentity_element_fixed_ranks"]) for row in v4_rows
    )
    if not (
        len(v4_rows) == 120
        and v4_rank_histogram == {7: 90, 8: 30}
        and v4_pattern_histogram == {(11, 11, 12): 90, (12, 12, 12): 30}
    ):
        raise ArithmeticError("the high-element V4 closure changed")

    checks = {
        "all_720_elements_ranked_exactly": sum(histogram.values()) == 720,
        "source_involution_type_counts_replayed": (
            fixed_ranks.count(11) == 36 and fixed_ranks.count(12) == 45
        ),
        "all_81_high_nonidentity_elements_are_involutions": len(high) == 81,
        "all_high_closed_V4_subgroups_enumerated": len(v4_rows) == 120,
        "every_high_closed_V4_drops_below_rank_11": all(
            row["common_fixed_rank"] < 11 for row in v4_rows
        ),
        "only_C2_can_be_a_proper_fixed_rank_at_least_11_subgroup": True,
    }
    if not all(checks.values()):
        raise ArithmeticError("proper high-subgroup closure failed")
    return {
        "fixed_rank_histogram": {str(key): value for key, value in sorted(histogram.items())},
        "high_nonidentity_element_count": len(high),
        "high_involution_classes": [
            {
                "surface_type": "rational",
                "source_cardinality": 36,
                "fixed_rank": 11,
                "hyperbolic_dimension": 10,
                "representative_zero_based": RATIONAL_REPRESENTATIVE,
            },
            {
                "surface_type": "symplectic",
                "source_cardinality": 45,
                "fixed_rank": 12,
                "hyperbolic_dimension": 11,
                "representative_zero_based": SYMPLECTIC_REPRESENTATIVE,
            },
        ],
        "v4_subgroup_count": len(v4_rows),
        "v4_common_fixed_rank_histogram": {
            str(key): value for key, value in sorted(v4_rank_histogram.items())
        },
        "v4_element_rank_pattern_histogram": {
            ",".join(map(str, key)): value
            for key, value in sorted(v4_pattern_histogram.items())
        },
        "v4_subgroups": v4_rows,
        "exact_int64_audit": {
            "maximum_matrix_entry_absolute": maximum_entry,
            "single_product_entry_absolute_bound": product_bound,
            "below_signed_int64_limit": True,
        },
        "proof_of_exhaustion": (
            "If Fix(H) has rank at least 11, every element of H has individual "
            "fixed rank at least 11.  The 81 possibilities are all involutions.  "
            "For two distinct such elements their product is again high; hence it "
            "is an involution, the two commute, and they generate one of the 120 "
            "enumerated V4 groups.  Every such V4 has common fixed rank 7 or 8.  "
            "Therefore a proper high subgroup is C2, in exactly the two classes above."
        ),
        "checks": checks,
    }


def _relation_audit(facets: Sequence[sp.Matrix], gram: sp.Matrix) -> dict[str, Any]:
    counts: Counter[str] = Counter()
    forbidden = []
    for right in range(len(facets)):
        for left in range(right):
            left_norm = sp.Rational((facets[left] * gram * facets[left].T)[0])
            right_norm = sp.Rational((facets[right] * gram * facets[right].T)[0])
            pairing = sp.Rational((facets[left] * gram * facets[right].T)[0])
            ratio = sp.factor(pairing**2 / (left_norm * right_norm))
            if pairing == 0:
                label = "2"
            elif pairing < 0:
                label = "forbidden"
            elif ratio == sp.Rational(1, 4):
                label = "3"
            elif ratio == sp.Rational(1, 2):
                label = "4"
            elif ratio == sp.Rational(3, 4):
                label = "6"
            elif ratio == 1:
                label = "infinity"
            elif ratio > 1:
                label = "ultraparallel"
            else:
                label = "forbidden"
            counts[label] += 1
            if label == "forbidden" and len(forbidden) < 12:
                forbidden.append(
                    {
                        "facet_pair_one_based": [left + 1, right + 1],
                        "inner_product": _fraction_text(pairing),
                        "squares": [_fraction_text(left_norm), _fraction_text(right_norm)],
                        "normalized_square_ratio": _fraction_text(ratio),
                    }
                )
    for label in ("2", "3", "4", "6", "infinity", "ultraparallel", "forbidden"):
        counts[label] += 0
    return {
        "pair_count": len(facets) * (len(facets) - 1) // 2,
        "coxeter_edge_counts": dict(counts),
        "all_pairs_coxeter": counts["forbidden"] == 0,
        "first_forbidden_pairs": forbidden,
    }


def _fold(
    representative: int,
    surface_type: str,
    gram: sp.Matrix,
    ample: sp.Matrix,
    walls: Sequence[sp.Matrix],
    group: Sequence[sp.Matrix],
) -> dict[str, Any]:
    action = group[representative]
    if action**2 != sp.eye(20):
        raise ArithmeticError("fold representative is not an involution")
    wall_index = {tuple(sp.Rational(value) for value in wall): i for i, wall in enumerate(walls)}
    permutation = []
    for wall in walls:
        image = tuple(sp.Rational(value) for value in wall * action)
        if image not in wall_index:
            raise ArithmeticError("the involution does not permute all 478 walls")
        permutation.append(wall_index[image])

    seen: set[int] = set()
    wall_orbits: list[list[int]] = []
    for start in range(len(walls)):
        if start in seen:
            continue
        orbit = []
        current = start
        while current not in orbit:
            orbit.append(current)
            seen.add(current)
            current = permutation[current]
        if current != start or len(orbit) not in {1, 2}:
            raise ArithmeticError("unexpected wall cycle")
        wall_orbits.append(orbit)

    unique: dict[tuple[int, ...], dict[str, Any]] = {}
    for orbit in wall_orbits:
        restriction = sum((walls[index] for index in orbit), sp.zeros(1, 20))
        primitive = _primitive_integer_row(restriction)
        record = unique.setdefault(
            primitive,
            {
                "primitive_fixed_vector_in_published_NS_basis": list(primitive),
                "source_wall_orbits_one_based": [],
            },
        )
        record["source_wall_orbits_one_based"].append([index + 1 for index in orbit])

    rays = [sp.Matrix(1, 20, primitive) for primitive in unique]
    records = list(unique.values())
    for record, ray in zip(records, rays):
        norm = sp.Rational((ray * gram * ray.T)[0])
        ample_pairing = sp.Rational((ample * gram * ray.T)[0])
        if ample_pairing <= 0 or ray * action != ray:
            raise ArithmeticError("a restricted half-space lost its orientation or fixedness")
        record["source_orbit_multiplicity"] = len(record["source_wall_orbits_one_based"])
        record["square"] = _fraction_text(norm)
        record["ample_pairing"] = _fraction_text(ample_pairing)
        record["lorentz_role"] = (
            "hyperbolic_facet_candidate"
            if norm < 0
            else "ideal_support_nonfacet"
            if norm == 0
            else "timelike_support_nonfacet"
        )

    facet_positions = [
        index
        for index, ray in enumerate(rays)
        if (ray * gram * ray.T)[0] < 0
    ]
    facets = [rays[index] for index in facet_positions]
    witnesses = []
    for position in facet_positions:
        target = rays[position]
        target_norm = sp.Rational((target * gram * target.T)[0])
        witness = ample - sp.Rational(
            (ample * gram * target.T)[0], target_norm
        ) * target
        slacks = [sp.Rational((witness * gram * ray.T)[0]) for ray in rays]
        others = [value for index, value in enumerate(slacks) if index != position]
        if not (
            (witness * gram * witness.T)[0] > 0
            and slacks[position] == 0
            and min(others) > 0
        ):
            raise ArithmeticError("a restricted negative wall is not a strict facet")
        witnesses.append(
            {
                "unique_restriction_index_one_based": position + 1,
                "witness_in_published_NS_basis": [
                    _fraction_text(value) for value in witness
                ],
                "witness_square": _fraction_text((witness * gram * witness.T)[0]),
                "minimum_other_restriction_pairing": _fraction_text(min(others)),
            }
        )
        records[position]["lorentz_role"] = "strict_hyperbolic_facet"

    facet_gram = sp.Matrix(
        [[(left * gram * right.T)[0] for right in facets] for left in facets]
    )
    positive_wall_gram = -facet_gram
    relations = _relation_audit(facets, gram)
    sign_histogram = Counter(
        "negative" if (ray * gram * ray.T)[0] < 0
        else "zero" if (ray * gram * ray.T)[0] == 0
        else "positive"
        for ray in rays
    )
    for sign in ("negative", "zero", "positive"):
        sign_histogram[sign] += 0
    multiplicity_histogram = Counter(
        len(record["source_wall_orbits_one_based"]) for record in records
    )
    result = {
        "surface_involution_type": surface_type,
        "representative_zero_based": representative,
        "fixed_rank": 20 - _rank(action.T - sp.eye(20)),
        "hyperbolic_dimension": 19 - _rank(action.T - sp.eye(20)),
        "complete_source_wall_count": len(walls),
        "wall_orbit_count": len(wall_orbits),
        "wall_orbit_size_histogram": {
            str(key): value
            for key, value in sorted(Counter(map(len, wall_orbits)).items())
        },
        "positive_ray_deduplicated_restriction_count": len(rays),
        "duplicate_multiplicity_histogram": {
            str(key): value for key, value in sorted(multiplicity_histogram.items())
        },
        "restriction_square_sign_histogram": dict(sign_histogram),
        "strict_facet_count": len(facets),
        "strict_facet_witnesses": witnesses,
        "unique_restrictions": records,
        "facet_gram": _payload(facet_gram),
        "positive_wall_gram": _payload(positive_wall_gram),
        "positive_wall_gram_sha256": _canonical_sha256(_payload(positive_wall_gram)),
        "facet_gram_rank": _rank(facet_gram),
        "positive_wall_gram_inertia": list(exact_inertia(positive_wall_gram)),
        "facet_norm_histogram": {
            str(key): value
            for key, value in sorted(
                Counter(int(facet_gram[i, i]) for i in range(len(facets))).items()
            )
        },
        "pair_gate": relations,
        "_facet_rays": facets,
        "_all_rays": rays,
        "_facet_gram_matrix": facet_gram,
        "_positive_wall_gram_matrix": positive_wall_gram,
    }
    return result


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    norm = sp.Rational((root.T * gram * root)[0])
    return sp.eye(gram.rows) - sp.Rational(2, norm) * root * (root.T * gram)


def _fixed_lattice_audit(
    action: sp.Matrix,
    ambient_gram: sp.Matrix,
    facets: Sequence[sp.Matrix],
) -> dict[str, Any]:
    basis, kernel = _saturated_fixed_basis((action,), 20)
    fixed_gram = basis * ambient_gram * basis.T
    _, pivot_columns = basis.rref()
    pivot_columns = list(pivot_columns)
    pivot_inverse = basis[:, pivot_columns].inv()
    coordinates = sp.Matrix.vstack(
        *(facet[:, pivot_columns] * pivot_inverse for facet in facets)
    )
    if not all(value.is_Integer for value in coordinates):
        raise ArithmeticError("a facet root left the saturated fixed lattice")
    coordinate_smith = _smith(coordinates)

    reflections = []
    mirror_records = []
    for index in range(coordinates.rows):
        root = coordinates.row(index).T
        norm = sp.Rational((root.T * fixed_gram * root)[0])
        divisor = math.gcd(
            *(abs(int(value)) for value in (root.T * fixed_gram))
        )
        content = math.gcd(*(abs(int(value)) for value in root))
        reflection = _reflection(fixed_gram, root)
        integral = all(value.is_Integer for value in reflection)
        if not (
            integral
            and abs(int(reflection.det())) == 1
            and reflection.T * fixed_gram * reflection == fixed_gram
        ):
            raise ArithmeticError("a facet reflection is not integral on Fix(iota)")
        reflections.append(reflection)
        mirror_records.append(
            {
                "facet_one_based": index + 1,
                "fixed_basis_coordinates": [int(value) for value in root],
                "coordinate_content": content,
                "square": int(norm),
                "divisibility": divisor,
                "reflection_integral": integral,
            }
        )

    stable = _stable_subgroups(fixed_gram, reflections)
    stable.pop("groups_internal")
    stable.pop("left_smith_inverse")
    if not (
        stable["stable_integral_subgroup_orders"] == [1]
        and stable["stable_even_subgroup_orders"] == [1]
    ):
        raise ArithmeticError("the reflection-stable overclosure census changed")

    smith = _smith(fixed_gram)
    inverse = fixed_gram.inv()
    delta = int(any(sp.Rational(inverse[i, i]).q == 2 for i in range(11)))
    standard = sp.diag(
        sp.Matrix([[0, 1], [1, 0]]),
        -2 * _e8_cartan(),
        sp.Matrix([[-2]]),
    )
    scaled_dual = sp.simplify(2 * inverse)
    target_checks = {
        "rank_11_signature_1_10": fixed_gram.rows == 11
        and exact_inertia(fixed_gram) == (1, 10, 0),
        "even": all(fixed_gram[i, i] % 2 == 0 for i in range(11)),
        "determinant_512": fixed_gram.det() == 512,
        "smith_1^2_2^9": smith == [1, 1] + [2] * 9,
        "two_elementary_delta_1": all(2 * value in ZZ for value in inverse)
        and delta == 1,
        "standard_model_has_same_Nikulin_invariants": (
            standard.det() == 512
            and _smith(standard) == [1, 1] + [2] * 9
            and exact_inertia(standard) == (1, 10, 0)
        ),
        "all_37_roots_primitive": all(
            record["coordinate_content"] == 1 for record in mirror_records
        ),
        "all_37_reflections_integral": all(
            record["reflection_integral"] for record in mirror_records
        ),
        "facet_roots_generate_entire_fixed_lattice": coordinate_smith == [1] * 11,
        "no_proper_integral_reflection_stable_overclosure": stable[
            "stable_integral_subgroup_orders"
        ]
        == [1],
        "scaled_dual_integral_odd_determinant_4": (
            all(value.is_Integer for value in scaled_dual)
            and any(scaled_dual[i, i] % 2 for i in range(11))
            and abs(int(scaled_dual.det())) == 4
            and _smith(scaled_dual) == [1] * 9 + [2, 2]
        ),
    }
    if not all(target_checks.values()):
        raise ArithmeticError(f"fixed-lattice classification changed: {target_checks}")

    return {
        "basis_rows_in_published_NS_basis": _integer_payload(basis),
        "smith_kernel_audit": kernel,
        "gram": _integer_payload(fixed_gram),
        "rank": 11,
        "inertia": [1, 10, 0],
        "determinant": 512,
        "even": True,
        "smith_diagonal_absolute": smith,
        "nikulin_invariants": {"r": 11, "a": 9, "delta": 1},
        "classification": (
            "By uniqueness of indefinite even 2-elementary lattices with "
            "these invariants: U + E8(-2) + <-2>."
        ),
        "standard_model_gram": _integer_payload(standard),
        "discriminant_q_mod_2_histogram": _q_histogram(fixed_gram),
        "mirror_type_histogram": {
            key: value
            for key, value in sorted(
                Counter(
                    f"{row['square']},div{row['divisibility']}"
                    for row in mirror_records
                ).items()
            )
        },
        "mirror_records": mirror_records,
        "facet_mirror_coordinate_smith_diagonal": coordinate_smith,
        "facet_reflection_actions_sha256": _canonical_sha256(
            [_payload(reflection) for reflection in reflections]
        ),
        "exhaustive_reflection_stable_overlattice_audit": stable,
        "scaled_dual_lineage": {
            "operation": "F^vee(2), Gram 2 F^{-1}; not an overlattice",
            "gram": _integer_payload(scaled_dual),
            "determinant": 4,
            "parity": "odd",
            "smith_diagonal_absolute": [1] * 9 + [2, 2],
            "scharlau_kirschmer_class": "List10 entry 6, genus I_{10,1}(2^2), after sign convention",
        },
        "checks": target_checks,
        "_gram_matrix": fixed_gram,
    }


def _run_vinberg(
    positive_wall_gram: sp.Matrix,
    cpp_source: Path = CPP_SOURCE,
    compiler: str | None = None,
) -> dict[str, Any]:
    if _sha256(cpp_source) != CPP_SHA256:
        raise ArithmeticError("the literal Vinberg verifier source changed")
    selected = compiler or shutil.which("clang++") or shutil.which("g++")
    if selected is None:
        raise RuntimeError("a C++17 compiler is required")
    if any(value.is_Integer is not True for value in positive_wall_gram):
        raise ArithmeticError("Vinberg input is not integral")
    rows = [f"10 {positive_wall_gram.rows}"]
    rows.extend(
        " ".join(str(int(value)) for value in positive_wall_gram.row(index))
        for index in range(positive_wall_gram.rows)
    )
    stream = "\n".join(rows) + "\n"
    with tempfile.TemporaryDirectory(prefix="shimada-x2-h10-37-vinberg-") as raw:
        executable = Path(raw) / "vinberg"
        subprocess.run(
            [
                selected,
                "-O3",
                "-std=c++17",
                "-Wall",
                "-Wextra",
                "-Werror",
                str(cpp_source),
                "-o",
                str(executable),
            ],
            check=True,
        )
        completed = subprocess.run(
            [str(executable)],
            input=stream,
            text=True,
            capture_output=True,
            check=True,
            timeout=3600,
        )
    record = _parse_output(completed.stdout)
    special = record["finite_special_subset_counts_by_rank"]
    endpoints = record["edge_endpoint_type_histogram"]
    checks = {
        "literal_output_matches_frozen_exact_census": record == EXPECTED_VINBERG,
        "all_finite_special_subsets_enumerated": (
            record["all_finite_subset_search_node_count"] == 1 + sum(special)
        ),
        "all_elliptic_rank_9_edges_and_rank_10_vertices_enumerated": (
            special[-2] == 24350 and special[-1] == 3236
        ),
        "all_full_affine_zero_neighborhoods_are_rank_9_parabolic": (
            record["bad_parabolic_neighborhood_count"] == 0
        ),
        "every_edge_has_exactly_two_finite_or_ideal_endpoints": (
            record["bad_edge_count"] == 0
            and sum(endpoints.values()) == special[-2]
            and sum(
                (int(key.split(",")[0]) + int(key.split(",")[1])) * value
                for key, value in endpoints.items()
            )
            == 2 * special[-2]
        ),
        "finite_vertex_edge_incidence": (
            2 * endpoints["2,0"] + endpoints["1,1"] == 10 * special[-1]
        ),
        "cusps_exist": record["ideal_vertex_count"] == 261,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"literal Vinberg audit changed: {checks}")
    return {
        **record,
        "elliptic_edge_count": special[-2],
        "finite_vertex_count": special[-1],
        "literal_stdout": completed.stdout,
        "criterion": (
            "Vinberg's edge-endpoint criterion plus exhaustive full zero-"
            "neighbourhood tests for every irreducible affine component"
        ),
        "finite_volume": True,
        "compact": False,
        "classification": "finite_volume_noncompact_Coxeter_polytope",
        "checks": checks,
    }


def _lineage(
    candidate_positive_gram: sp.Matrix,
    fixed_lattice: dict[str, Any],
) -> dict[str, Any]:
    if _sha256(K3_FULL) != K3_FULL_SHA256:
        raise ArithmeticError("the frozen K3-S3 full-chamber comparator changed")
    k3 = json.loads(K3_FULL.read_text(encoding="utf-8"))
    full = k3["survivor_geometry"]
    if not (
        full["id"] == "K3-S3-C2-H10-13"
        and full["fixed_lattice"]["determinant"] == 512
        and full["fixed_lattice"]["smith_diagonal_absolute"] == [1, 1] + [2] * 9
        and full["fixed_lattice"]["even"] is True
        and full["fixed_lattice"]["facet_mirrors_generate_saturated_fixed_lattice"]
    ):
        raise ArithmeticError("the (11,9,1) full-chamber comparator changed")

    full_positive_gram = -_matrix(full["facet_gram"])
    candidate_solomon = _solomon(_payload(candidate_positive_gram))
    full_solomon = _solomon(_payload(full_positive_gram))
    candidate_chi = Fraction(candidate_solomon["rational_euler_characteristic"])
    full_chi = Fraction(full_solomon["rational_euler_characteristic"])
    index = candidate_chi / full_chi
    if index.denominator != 1 or index != 241920:
        raise ArithmeticError("the exact full-reflection-group index changed")

    relation_transcription_hash = _canonical_sha256(
        full["normalized_coxeter_relation_matrix"]
    )
    checks = {
        "same_unique_even_2_elementary_lattice_11_9_1": (
            fixed_lattice["nikulin_invariants"] == {"r": 11, "a": 9, "delta": 1}
        ),
        "published_full_Wr_diagram_has_13_facets": len(full["facet_gram"]) == 13,
        "full_chamber_has_34_finite_vertices_4_cusps": (
            full["exact_cone"]["finite_vertex_count"] == 34
            and full["exact_cone"]["ideal_vertex_count"] == 4
        ),
        "full_chamber_Solomon_characteristic": str(full_chi)
        == "-527/33443020800",
        "candidate_Solomon_characteristic": str(candidate_chi) == "-527/138240",
        "exact_index_241920": index == 241920,
        "candidate_reflections_lie_in_full_root_reflection_group": all(
            row["reflection_integral"] for row in fixed_lattice["mirror_records"]
        ),
        "no_new_integral_or_arithmetic_lattice_class": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"full-group lineage changed: {checks}")
    return {
        "fixed_lattice_identification": (
            "Both lattices are the unique indefinite even 2-elementary lattice "
            "(r,a,delta)=(11,9,1), hence are integrally isometric to "
            "U+E8(-2)+<-2>."
        ),
        "published_full_reflection_chamber": {
            "id": "K3-S3-C2-H10-13",
            "facet_count": 13,
            "finite_vertex_count": 34,
            "cusp_count": 4,
            "source": (
                "Valery Alexeev and Philip Engel, Compactifications of moduli "
                "spaces of K3 surfaces with a nonsymplectic involution, "
                "arXiv:2208.10383, finite-covolume theorem and Coxeter figure "
                "for (11,9,1)"
            ),
            "source_url": "https://arxiv.org/abs/2208.10383",
            "source_tex_sha256": "76b4c279d6ae490c8f76ef05f6b980d8536b7b00463ec965270ff445ecbd48c0",
            "figure_relation_transcription_sha256": relation_transcription_hash,
            "scharlau_kirschmer": {
                "list": "List10",
                "entry_one_based": 6,
                "genus_symbol": "I_{10,1}(2^2)",
                "facet_count": 13,
                "face_vector": [13, 76, 265, 615, 1001, 1168, 975, 565, 208, 34],
                "cusp_count": 4,
                "list_sha256": "9e627c4e1303015fd96acdcf6ccbeb264ebadba49abfcd761705f54d211b2e21",
                "gram_sha256": "b8186e1fa2771c517a6040998f3b234568c3edde3ec8b7bb8412865d7710af18",
            },
        },
        "solomon": {
            "candidate": candidate_solomon,
            "full_Wr_chamber": full_solomon,
            "index_candidate_in_full_Wr": int(index),
            "index_formula": "chi(W_candidate)/chi(W_r)",
        },
        "reflection_stable_overlattice_status": (
            "closed: F itself is the only integral overclosure invariant under "
            "all 37 displayed reflections"
        ),
        "bounded_literal_diagram_search": {
            "pre_catch_ledger_sha256": "3482a0faa91604b5bbdba633076f3ce2873763113737250db34755a20ec245b1",
            "pre_catch_same_header_count": 0,
            "coxiter_repository": "https://github.com/rgugl/CoxIter",
            "coxiter_commit": "b800be48d8240fbabd2d1c9542e61054584ad7fe",
            "ordinary_file_count": 1214,
            "ordinary_H10_37_header_count": 0,
            "ordinary_archive_sha256": "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3",
            "doubling_file_count": 2959,
            "doubling_H10_37_header_count": 0,
            "doubling_archive_sha256": "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036",
            "scope_warning": (
                "Header absence is exact only in these pinned project and CoxIter "
                "snapshots.  It is not a worldwide priority theorem."
            ),
        },
        "novelty_status": {
            "new_arithmetic_lattice_class": False,
            "new_full_reflection_group": False,
            "new_abstract_polytope_worldwide": "open",
            "new_algebro_geometric_realization_candidate": True,
            "precise_claim": (
                "A newly extracted 37-facet finite-index chamber from a genuine "
                "rational automorphism of Shimada's singular K3 surface X_2; its "
                "fixed lattice, containing full chamber, and index 241920 are known."
            ),
        },
        "checks": checks,
    }


def _public(record: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in record.items() if not key.startswith("_")}


def build_certificate(
    cpp_source: Path = CPP_SOURCE,
    compiler: str | None = None,
) -> dict[str, Any]:
    stock, gram, ample, walls, group = _load_denominator()
    subgroups = _subgroup_denominator(group)
    rational = _fold(
        RATIONAL_REPRESENTATIVE, "rational", gram, ample, walls, group
    )
    symplectic = _fold(
        SYMPLECTIC_REPRESENTATIVE, "symplectic", gram, ample, walls, group
    )

    expected_rational = {
        "fixed_rank": 11,
        "wall_orbit_count": 242,
        "positive_ray_deduplicated_restriction_count": 207,
        "restriction_square_sign_histogram": {
            "negative": 37,
            "zero": 45,
            "positive": 125,
        },
        "strict_facet_count": 37,
        "facet_norm_histogram": {"-4": 10, "-2": 27},
    }
    expected_symplectic = {
        "fixed_rank": 12,
        "wall_orbit_count": 266,
        "positive_ray_deduplicated_restriction_count": 266,
        "restriction_square_sign_histogram": {
            "negative": 90,
            "zero": 0,
            "positive": 176,
        },
        "strict_facet_count": 90,
        "facet_norm_histogram": {
            "-120": 12,
            "-30": 24,
            "-20": 26,
            "-12": 8,
            "-4": 16,
            "-2": 4,
        },
    }
    for key, value in expected_rational.items():
        if rational[key] != value:
            raise ArithmeticError(f"rational fold changed at {key}")
    for key, value in expected_symplectic.items():
        if symplectic[key] != value:
            raise ArithmeticError(f"symplectic fold changed at {key}")
    if rational["pair_gate"]["coxeter_edge_counts"] != {
        "2": 405,
        "3": 15,
        "4": 70,
        "6": 0,
        "infinity": 121,
        "ultraparallel": 55,
        "forbidden": 0,
    }:
        raise ArithmeticError("the H10/37 Coxeter relation census changed")
    if symplectic["pair_gate"]["coxeter_edge_counts"] != {
        "2": 1152,
        "3": 48,
        "4": 8,
        "6": 32,
        "infinity": 155,
        "ultraparallel": 2258,
        "forbidden": 352,
    }:
        raise ArithmeticError("the symplectic non-Coxeter obstruction changed")

    fixed = _fixed_lattice_audit(
        group[RATIONAL_REPRESENTATIVE], gram, rational["_facet_rays"]
    )
    vinberg = _run_vinberg(
        rational["_positive_wall_gram_matrix"], cpp_source, compiler
    )
    lineage = _lineage(rational["_positive_wall_gram_matrix"], fixed)

    checks = {
        "complete_478_wall_720_group_source_replayed": True,
        "all_proper_fixed_rank_at_least_11_subgroups_exhausted": all(
            subgroups["checks"].values()
        ),
        "rational_class_is_strict_integral_Coxeter_H10_37": (
            rational["pair_gate"]["all_pairs_coxeter"] is True
            and rational["strict_facet_count"] == 37
            and all(fixed["checks"].values())
        ),
        "symplectic_class_decisively_non_Coxeter": (
            symplectic["pair_gate"]["coxeter_edge_counts"]["forbidden"] == 352
        ),
        "literal_Vinberg_finite_volume_noncompact": (
            vinberg["finite_volume"] is True
            and vinberg["compact"] is False
            and all(vinberg["checks"].values())
        ),
        "fixed_lattice_and_full_group_lineage_closed": all(
            lineage["checks"].values()
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"composite X_2 certificate failed: {checks}")

    fixed_public = _public(fixed)
    return {
        "schema": SCHEMA,
        "id": "SHIMADA-X2-RATIONAL-C2-H10-37",
        "status_date": "2026-08-07",
        "source": {
            "stock_path": str(STOCK.relative_to(ROOT)),
            "stock_sha256": STOCK_SHA256,
            "paper": stock["sources"]["paper"],
            "paper_url": stock["sources"]["paper_url"],
            "supplement_zip_url": stock["sources"]["supplement_zip_url"],
            "supplement_zip_sha256": stock["sources"]["supplement_zip_sha256"],
            "theorem_scope": stock["sources"]["theorem_scope"],
            "involution_type_source": (
                "Shimada Proposition prop:InvolsD0: X_2 has exactly 45 "
                "symplectic and 36 rational involutions in Aut(X_2,a_2)."
            ),
        },
        "method": {
            "algebro_geometric_carrier": (
                "the saturated invariant Neron--Severi lattice and the restriction "
                "of Shimada's theorem-complete induced nef chamber"
            ),
            "subgroup_quantifier": (
                "all subgroups capable of fixed rank at least 11 in the complete "
                "720-element genuine surface-automorphism group"
            ),
            "volume_criterion": (
                "literal Vinberg finite-special-subdiagram edge/endpoint audit and "
                "full affine zero-neighbourhood cusp audit"
            ),
        },
        "complete_high_subgroup_denominator": subgroups,
        "rational_involution_fold": _public(rational),
        "symplectic_involution_negative_screen": _public(symplectic),
        "saturated_fixed_divisor_lattice": fixed_public,
        "literal_vinberg_volume_certificate": vinberg,
        "arithmetic_and_literal_lineage": lineage,
        "decision": {
            "classification": "finite_volume_noncompact_hyperbolic_Coxeter_polytope",
            "hyperbolic_dimension": 10,
            "facet_count": 37,
            "finite_vertex_count": 3236,
            "cusp_count": 261,
            "elliptic_edge_count": 24350,
            "fixed_lattice": "U+E8(-2)+<-2>",
            "index_in_full_reflection_group": 241920,
            "genuine_surface_automorphism": True,
            "worldwide_literal_priority": "open",
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpp-source", type=Path, default=CPP_SOURCE)
    parser.add_argument("--compiler")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.cpp_source, args.compiler)
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != encoded:
            raise ArithmeticError(f"archive mismatch: {args.output}")
        print(f"archive matches: {args.output}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded, encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "classification": payload["decision"]["classification"],
                "dimension": payload["decision"]["hyperbolic_dimension"],
                "facets": payload["decision"]["facet_count"],
                "finite_vertices": payload["decision"]["finite_vertex_count"],
                "cusps": payload["decision"]["cusp_count"],
                "full_group_index": payload["decision"]["index_in_full_reflection_group"],
                "all_checks_pass": payload["all_checks_pass"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
