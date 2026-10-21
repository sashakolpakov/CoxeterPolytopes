#!/usr/bin/env python3
"""Certify two exact all-length H14 parallel-doubling path families.

The frozen H14/100 and H14/101 chambers have the same explicit pair of
parallel doubling walls.  This audit uses only their integral wall vectors
and the ambient Gram matrix.  In particular, it never invokes cdd, CoxIter,
or an extreme-ray enumeration.

For each parent, the two endpoint reflections generate an infinite dihedral
group.  The parabolic translation T satisfies (T-I)^3=0, so every possible
collision between transported side walls reduces to an integer quadratic
equation.  Solving those equations for all z in Z gives an all-N, rather than
bounded, facet formula.  A small direct reconstruction is retained only as a
redundant regression check.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable, Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
SINGLE = ROOT / "data/enumerations/coble_type_vi_h14_single_doublings.json"
SOURCE = (
    ROOT
    / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json"
)
SECOND = (
    ROOT
    / "data/enumerations/coble_type_vi_h14_second_generation_doublings.json"
)
OUTPUT = (
    ROOT
    / "data/enumerations/coble_type_vi_h14_parallel_path_families.json"
)
REPORT = ROOT / "results/coble_type_vi_h14_parallel_path_families.md"

SCHEMA = "coxeter.coble_type_vi_h14_parallel_path_families.v1"
PINNED = {
    "single_doublings": "d32c9dbd96d1bb97c49736bd0d3d4fb6b967eb3754a771fc61d557d9b3a7d88d",
    "type_vi_source": "8247ad85f7970d3bf1f08fe15517b70376409c17cc2784b8f14b129aaecc2ede",
    "second_generation": "11da08bf3e663807fd76a0a022288a0f15140d2fc342c890066b595c706244a7",
}
PARENTS = ("H14/100", "H14/101")
EXPECTED = {
    "H14/100": {
        "facet_count": 100,
        "neither": 77,
        "finite_collision_count": 91,
        "formula": "81*N+19-floor(N/2)",
        "second_generation_id": "COBLE-VI3-D4A02-H14-180",
        "second_generation_facets": 180,
    },
    "H14/101": {
        "facet_count": 101,
        "neither": 78,
        "finite_collision_count": 92,
        "formula": "82*N+19-floor(N/2)",
        "second_generation_id": "COBLE-VI3-D4B04-H14-182",
        "second_generation_facets": 182,
    },
}

IntegralVector = tuple[int, ...]
IntegralMatrix = tuple[tuple[int, ...], ...]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_hash(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _integral_matrix(matrix: sp.Matrix) -> IntegralMatrix:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("an integral matrix acquired a denominator")
    return tuple(tuple(int(value) for value in matrix.row(row)) for row in range(matrix.rows))


def _matvec(matrix: IntegralMatrix, vector: IntegralVector) -> IntegralVector:
    return tuple(
        sum(value * coordinate for value, coordinate in zip(row, vector, strict=True))
        for row in matrix
    )


def _scale(coefficient: int, vector: IntegralVector) -> IntegralVector:
    return tuple(coefficient * value for value in vector)


def _add(*vectors: IntegralVector) -> IntegralVector:
    return tuple(sum(values) for values in zip(*vectors, strict=True))


def _content(vector: IntegralVector) -> int:
    return math.gcd(*(abs(value) for value in vector))


def _primitive(vector: Sequence[int | Fraction]) -> IntegralVector:
    denominator = math.lcm(*(value.denominator if isinstance(value, Fraction) else 1 for value in vector))
    integers = tuple(
        int(value * denominator) if isinstance(value, Fraction) else value * denominator
        for value in vector
    )
    content = _content(integers)
    return tuple(value // content for value in integers)


def _unoriented_key(vector: IntegralVector) -> IntegralVector:
    content = _content(vector)
    values = tuple(value // content for value in vector)
    first = next(value for value in values if value)
    return tuple(-value for value in values) if first < 0 else values


def _pair(
    gram: IntegralMatrix, left: IntegralVector, right: IntegralVector
) -> int:
    return sum(
        left[row] * gram[row][column] * right[column]
        for row in range(len(left))
        for column in range(len(right))
    )


def _reflection(gram: sp.Matrix, root: IntegralVector) -> sp.Matrix:
    vector = sp.Matrix(root)
    square = (vector.T * gram * vector)[0]
    return sp.eye(gram.rows) - sp.Rational(2, square) * vector * (vector.T * gram)


def _smith_diagonal(vectors: Sequence[IntegralVector]) -> list[int]:
    matrix = sp.Matrix.hstack(*(sp.Matrix(vector) for vector in vectors))
    return [
        abs(int(value))
        for value in smith_normal_form(matrix, domain=ZZ).diagonal()
    ]


def _integer_roots_of_scalar_polynomial(
    quadratic: int, linear: int, constant: int
) -> tuple[int, ...]:
    if quadratic == 0:
        if linear == 0 or (-constant) % linear:
            return ()
        return ((-constant) // linear,)
    discriminant = linear * linear - 4 * quadratic * constant
    if discriminant < 0:
        return ()
    square_root = math.isqrt(discriminant)
    if square_root * square_root != discriminant:
        return ()
    denominator = 2 * quadratic
    roots: set[int] = set()
    for numerator in (-linear + square_root, -linear - square_root):
        if numerator % denominator == 0:
            roots.add(numerator // denominator)
    return tuple(sorted(roots))


def _translation_power_apply(
    nilpotent: IntegralMatrix,
    nilpotent_square: IntegralMatrix,
    vector: IntegralVector,
    exponent: int,
) -> IntegralVector:
    return _add(
        vector,
        _scale(exponent, _matvec(nilpotent, vector)),
        _scale(
            exponent * (exponent - 1) // 2,
            _matvec(nilpotent_square, vector),
        ),
    )


def _solve_translation_equation(
    nilpotent: IntegralMatrix,
    nilpotent_square: IntegralMatrix,
    vector: IntegralVector,
    target: IntegralVector,
) -> tuple[int, ...] | None:
    """Solve T^z(vector)=target for every z in Z; None means every z."""

    first = _matvec(nilpotent, vector)
    second = _matvec(nilpotent_square, vector)
    coefficients = tuple(
        (
            second[index],
            2 * first[index] - second[index],
            2 * (vector[index] - target[index]),
        )
        for index in range(len(vector))
    )
    pivot = next(
        (coefficient for coefficient in coefficients if coefficient != (0, 0, 0)),
        None,
    )
    if pivot is None:
        return None
    candidates = _integer_roots_of_scalar_polynomial(*pivot)
    return tuple(
        exponent
        for exponent in candidates
        if _translation_power_apply(
            nilpotent, nilpotent_square, vector, exponent
        )
        == target
    )


def _chain_image(
    chamber_index: int,
    vector: IntegralVector,
    right_reflection: IntegralMatrix,
    nilpotent: IntegralMatrix,
    nilpotent_square: IntegralMatrix,
) -> IntegralVector:
    exponent = chamber_index // 2
    source = vector if chamber_index % 2 == 0 else _matvec(right_reflection, vector)
    return _translation_power_apply(
        nilpotent, nilpotent_square, source, exponent
    )


def _facet_formula(parent: str, chamber_count: int) -> int:
    if chamber_count < 1:
        raise ValueError("a finite path must contain at least one chamber")
    slope = 81 if parent == "H14/100" else 82
    return slope * chamber_count + 19 - chamber_count // 2


def _direct_boundary(
    chamber_count: int,
    roots: Sequence[IntegralVector],
    left: int,
    right: int,
    right_reflection: IntegralMatrix,
    nilpotent: IntegralMatrix,
    nilpotent_square: IntegralMatrix,
) -> tuple[IntegralVector, ...]:
    nonend = tuple(
        index for index in range(len(roots)) if index not in {left, right}
    )
    boundary: dict[IntegralVector, IntegralVector] = {}
    for chamber_index in range(chamber_count):
        for root_index in nonend:
            image = _chain_image(
                chamber_index,
                roots[root_index],
                right_reflection,
                nilpotent,
                nilpotent_square,
            )
            boundary.setdefault(_unoriented_key(image), image)
    last = chamber_count - 1
    right_seed = right if last % 2 == 0 else left
    outer_right = _chain_image(
        last,
        roots[right_seed],
        right_reflection,
        nilpotent,
        nilpotent_square,
    )
    for endpoint in (roots[left], outer_right):
        key = _unoriented_key(endpoint)
        if key in boundary:
            raise ArithmeticError("an outer endpoint collided with a side wall")
        boundary[key] = endpoint
    return tuple(boundary[key] for key in sorted(boundary))


def _seed_support_audit(
    gram: sp.Matrix, roots: Sequence[IntegralVector], interior_payload: Sequence[Any]
) -> dict[str, Any]:
    """Replay the projected strict witness for every parent facet cheaply."""

    interior = sp.Matrix([sp.Rational(str(value)) for value in interior_payload])
    root_matrix = sp.Matrix.hstack(*(sp.Matrix(root) for root in roots))
    wall_gram = root_matrix.T * gram * root_matrix
    interior_pairs = list(interior.T * gram * root_matrix)
    interior_square = (interior.T * gram * interior)[0]
    minimum_slack: sp.Rational | None = None
    minimum_square: sp.Rational | None = None
    for target in range(len(roots)):
        coefficient = interior_pairs[target] / wall_gram[target, target]
        witness_square = interior_square - (
            interior_pairs[target] ** 2 / wall_gram[target, target]
        )
        slacks = [
            interior_pairs[other] - coefficient * wall_gram[target, other]
            for other in range(len(roots))
        ]
        if not (
            slacks[target] == 0
            and witness_square > 0
            and all(value > 0 for index, value in enumerate(slacks) if index != target)
        ):
            raise ArithmeticError("a projected parent-facet witness failed")
        local = min(value for index, value in enumerate(slacks) if index != target)
        minimum_slack = local if minimum_slack is None else min(minimum_slack, local)
        minimum_square = (
            witness_square
            if minimum_square is None
            else min(minimum_square, witness_square)
        )
    return {
        "facet_count": len(roots),
        "projected_strict_witness_count": len(roots),
        "minimum_off_target_slack": str(minimum_slack),
        "minimum_witness_square": str(minimum_square),
        "all_parent_facets_have_strict_timelike_witnesses": True,
        "all_N_transport_reason": (
            "transport a strict seed-facet witness by the chamber isometry; "
            "the witness lies away from both path glues unless it is the chosen "
            "outer endpoint, so every distinct boundary mirror supports a facet"
        ),
    }


def _collision_audit(
    roots: Sequence[IntegralVector],
    left: int,
    right: int,
    right_reflection: IntegralMatrix,
    nilpotent: IntegralMatrix,
    nilpotent_square: IntegralMatrix,
    categories: dict[str, list[int]],
    expected_finite_count: int,
) -> dict[str, Any]:
    nonend = tuple(
        index for index in range(len(roots)) if index not in {left, right}
    )
    branches: dict[str, Callable[[IntegralVector], IntegralVector]] = {
        "T^z": lambda vector: vector,
        "T^z*S_b": lambda vector: _matvec(right_reflection, vector),
    }
    finite_histogram: Counter[tuple[str, int, int]] = Counter()
    identity_histogram: Counter[str] = Counter()
    finite_count = 0
    identity_count = 0
    bad_finite: list[list[Any]] = []
    bad_identities: list[list[Any]] = []
    common = set(categories["orthogonal_to_both_one_based"])

    for branch_name, branch in branches.items():
        for source_index in nonend:
            vector = branch(roots[source_index])
            for target_index in nonend:
                for sign in (1, -1):
                    solutions = _solve_translation_equation(
                        nilpotent,
                        nilpotent_square,
                        vector,
                        _scale(sign, roots[target_index]),
                    )
                    if solutions is None:
                        identity_count += 1
                        identity_histogram[branch_name] += 1
                        if not (
                            source_index == target_index
                            and sign == 1
                            and source_index + 1 in common
                        ):
                            bad_identities.append(
                                [branch_name, source_index + 1, target_index + 1, sign]
                            )
                    else:
                        for exponent in solutions:
                            finite_count += 1
                            finite_histogram[(branch_name, exponent, sign)] += 1
                            if not (source_index == target_index and sign == 1):
                                bad_finite.append(
                                    [
                                        branch_name,
                                        source_index + 1,
                                        target_index + 1,
                                        sign,
                                        exponent,
                                    ]
                                )

    expected_histogram = Counter(
        {
            ("T^z", 0, 1): len(nonend) - 14,
            ("T^z*S_b", 0, 1): 4,
            ("T^z*S_b", -1, 1): 3,
        }
    )
    if not (
        finite_count == expected_finite_count
        and finite_histogram == expected_histogram
        and not bad_finite
        and identity_count == 28
        and identity_histogram == Counter({"T^z": 14, "T^z*S_b": 14})
        and not bad_identities
    ):
        raise ArithmeticError("the exhaustive side-wall collision census changed")

    endpoint_side_count = 0
    for branch in branches.values():
        for source_index in (left, right):
            vector = branch(roots[source_index])
            for target_index in nonend:
                for sign in (1, -1):
                    solutions = _solve_translation_equation(
                        nilpotent,
                        nilpotent_square,
                        vector,
                        _scale(sign, roots[target_index]),
                    )
                    if solutions is None:
                        endpoint_side_count += 1
                    else:
                        endpoint_side_count += len(solutions)
    if endpoint_side_count:
        raise ArithmeticError("an endpoint orbit met a side-wall orbit")

    endpoint_equations: dict[str, list[dict[str, int]]] = {}
    for name, vector in {
        "T^z(b)=plus_or_minus_a": roots[right],
        "T^z*S_b(a)=plus_or_minus_a": _matvec(
            right_reflection, roots[left]
        ),
    }.items():
        records: list[dict[str, int]] = []
        for sign in (1, -1):
            solutions = _solve_translation_equation(
                nilpotent,
                nilpotent_square,
                vector,
                _scale(sign, roots[left]),
            )
            if solutions is None:
                raise ArithmeticError("an unexpected all-z endpoint identity appeared")
            records.extend({"z": exponent, "sign": sign} for exponent in solutions)
        endpoint_equations[name] = records
    if not (
        endpoint_equations["T^z(b)=plus_or_minus_a"] == []
        and endpoint_equations["T^z*S_b(a)=plus_or_minus_a"]
        == [{"z": -1, "sign": -1}]
    ):
        raise ArithmeticError("the outer-endpoint collision equations changed")

    return {
        "domain": (
            f"all z in Z, two parity branches, both signs, and all "
            f"{len(nonend)}^2 ordered side-wall pairs"
        ),
        "method": (
            "T^z=I+zK+z(z-1)K^2/2 with K^3=0; solve one nonzero "
            "integer linear/quadratic coordinate and verify every coordinate"
        ),
        "bounded_exponent_scan_used": False,
        "finite_solution_count": finite_count,
        "finite_solution_histogram": {
            f"{branch};z={exponent};sign={sign}": count
            for (branch, exponent, sign), count in sorted(finite_histogram.items())
        },
        "all_z_identity_count": identity_count,
        "all_z_identity_histogram": dict(sorted(identity_histogram.items())),
        "identities_are_exactly_the_14_common_orthogonals": True,
        "cross_seed_wall_collision_count": len(bad_finite),
        "negative_root_collision_count": 0,
        "endpoint_to_side_collision_count": endpoint_side_count,
        "outer_endpoint_equations": endpoint_equations,
        "all_exact_checks_pass": True,
    }


def _second_generation_anchor(
    parent: str,
    gram: sp.Matrix,
    roots: Sequence[IntegralVector],
    interior_payload: Sequence[Any],
    selected: int,
    right_reflection: sp.Matrix,
    second: dict[str, Any],
) -> dict[str, Any]:
    """Rebuild the N=2 wall list and match its frozen second-generation row."""

    interior = sp.Matrix([sp.Rational(str(value)) for value in interior_payload])
    doubled_interior = interior + right_reflection * interior
    boundary: dict[IntegralVector, IntegralVector] = {}
    for copy in range(2):
        for index, root_tuple in enumerate(roots):
            if index == selected:
                continue
            root = sp.Matrix(root_tuple)
            image = root if copy == 0 else right_reflection * root
            if (doubled_interior.T * gram * image)[0] < 0:
                image = -image
            vector = tuple(int(value) for value in image)
            boundary.setdefault(_unoriented_key(vector), vector)

    expected = EXPECTED[parent]
    row = next(
        case
        for case in second["cases"]
        if case["id"] == expected["second_generation_id"]
    )
    digest = _json_hash([list(vector) for vector in boundary.values()])
    if not (
        selected + 1 == 20
        and len(boundary) == expected["second_generation_facets"]
        and row["parent"] == parent
        and row["selected_parent_wall_one_based"] == 20
        and row["facet_count"] == len(boundary)
        and row["arithmetic_wall_vectors_sha256"] == digest
    ):
        raise ArithmeticError("the N=2 path no longer matches its frozen child")
    return {
        "N": 2,
        "frozen_second_generation_id": row["id"],
        "selected_parent_wall_one_based": 20,
        "facet_count": len(boundary),
        "ordered_arithmetic_wall_vectors_sha256": digest,
        "exact_coordinate_list_match": True,
    }


def _certify_parent(
    parent: str,
    payload: dict[str, Any],
    gram: sp.Matrix,
    gram_integral: IntegralMatrix,
    ns_gram: sp.Matrix,
    second: dict[str, Any],
    direct_max_n: int,
) -> dict[str, Any]:
    construction = payload["construction"]
    roots = tuple(
        tuple(int(value) for value in vector)
        for vector in construction["wall_vectors_in_arithmetic_basis"]
    )
    physical_roots = tuple(
        tuple(int(value) for value in vector)
        for vector in construction["physical_NS_divisor_carrier"][
            "wall_vectors_in_NS_basis"
        ]
    )
    explicit = construction["redoubling_frontier"]["explicit_redoublable_pair"]
    left, right = explicit["wall_indices_zero_based"]
    expected = EXPECTED[parent]
    if not (
        len(roots) == expected["facet_count"]
        and (left, right) == (11, 19)
        and tuple(explicit["left_root"]) == roots[left]
        and tuple(explicit["right_root"]) == roots[right]
        and left in construction["redoubling_frontier"]["doubling_walls_zero_based"]
        and right in construction["redoubling_frontier"]["doubling_walls_zero_based"]
    ):
        raise ArithmeticError(f"the frozen {parent} endpoint pair changed")

    left_root = roots[left]
    right_root = roots[right]
    left_reflection_sympy = _reflection(gram, left_root)
    right_reflection_sympy = _reflection(gram, right_root)
    translation_sympy = right_reflection_sympy * left_reflection_sympy
    nilpotent_sympy = translation_sympy - sp.eye(gram.rows)
    nilpotent_square_sympy = nilpotent_sympy**2
    left_reflection = _integral_matrix(left_reflection_sympy)
    right_reflection = _integral_matrix(right_reflection_sympy)
    translation = _integral_matrix(translation_sympy)
    nilpotent = _integral_matrix(nilpotent_sympy)
    nilpotent_square = _integral_matrix(nilpotent_square_sympy)
    variable = translation_sympy.charpoly().gen
    half_sum = tuple(
        Fraction(a + b, 2) for a, b in zip(left_root, right_root, strict=True)
    )
    null_direction = _primitive(half_sum)
    if not (
        _pair(gram_integral, left_root, left_root) == -2
        and _pair(gram_integral, right_root, right_root) == -2
        and _pair(gram_integral, left_root, right_root) == 2
        and _pair(gram_integral, null_direction, null_direction) == 0
        and _content(null_direction) == 1
        and left_reflection_sympy.T * gram * left_reflection_sympy == gram
        and right_reflection_sympy.T * gram * right_reflection_sympy == gram
        and left_reflection_sympy**2 == right_reflection_sympy**2 == sp.eye(gram.rows)
        and sp.expand(translation_sympy.charpoly().as_expr())
        == sp.expand((variable - 1) ** gram.rows)
        and nilpotent_sympy.rank() == 2
        and nilpotent_square_sympy.rank() == 1
        and nilpotent_sympy**3 == sp.zeros(gram.rows)
        and nilpotent_square_sympy != sp.zeros(gram.rows)
    ):
        raise ArithmeticError(f"the {parent} dihedral translation gate failed")

    nonend = tuple(index for index in range(len(roots)) if index not in {left, right})
    orthogonal_left = {
        index for index in nonend if _pair(gram_integral, left_root, roots[index]) == 0
    }
    orthogonal_right = {
        index for index in nonend if _pair(gram_integral, right_root, roots[index]) == 0
    }
    category_sets = {
        "orthogonal_to_both_one_based": orthogonal_left & orthogonal_right,
        "orthogonal_to_a_only_one_based": orthogonal_left - orthogonal_right,
        "orthogonal_to_b_only_one_based": orthogonal_right - orthogonal_left,
        "orthogonal_to_neither_one_based": set(nonend)
        - (orthogonal_left | orthogonal_right),
    }
    categories = {
        key: sorted(index + 1 for index in indices)
        for key, indices in category_sets.items()
    }
    category_counts = {key: len(value) for key, value in categories.items()}
    if category_counts != {
        "orthogonal_to_both_one_based": 14,
        "orthogonal_to_a_only_one_based": 3,
        "orthogonal_to_b_only_one_based": 4,
        "orthogonal_to_neither_one_based": expected["neither"],
    }:
        raise ArithmeticError(f"the {parent} endpoint categories changed")

    collision = _collision_audit(
        roots,
        left,
        right,
        right_reflection,
        nilpotent,
        nilpotent_square,
        categories,
        expected["finite_collision_count"],
    )
    support = _seed_support_audit(
        gram, roots, construction["strict_interior"]
    )

    arithmetic_smith = _smith_diagonal(
        tuple(root for index, root in enumerate(roots) if index != right)
    )
    physical_smith = _smith_diagonal(
        tuple(root for index, root in enumerate(physical_roots) if index != right)
    )
    physical_left_reflection = _reflection(ns_gram, physical_roots[left])
    physical_right_reflection = _reflection(ns_gram, physical_roots[right])
    if not (
        arithmetic_smith == physical_smith == [1] * 15
        and _integral_matrix(physical_left_reflection)
        and _integral_matrix(physical_right_reflection)
        and physical_left_reflection.T * ns_gram * physical_left_reflection == ns_gram
        and physical_right_reflection.T * ns_gram * physical_right_reflection == ns_gram
    ):
        raise ArithmeticError(f"the {parent} full divisor-carrier gate failed")

    direct: list[dict[str, Any]] = []
    for chamber_count in range(1, direct_max_n + 1):
        boundary = _direct_boundary(
            chamber_count,
            roots,
            left,
            right,
            right_reflection,
            nilpotent,
            nilpotent_square,
        )
        expected_count = _facet_formula(parent, chamber_count)
        if not (
            len(boundary) == expected_count
            and all(_content(vector) == 1 for vector in boundary)
        ):
            raise ArithmeticError(
                f"the direct {parent}, N={chamber_count} boundary audit failed"
            )
        direct.append(
            {
                "N": chamber_count,
                "facet_count": len(boundary),
                "canonical_unoriented_wall_set_sha256": _json_hash(
                    [list(_unoriented_key(vector)) for vector in boundary]
                ),
                "all_walls_primitive_integral": True,
            }
        )

    anchor = _second_generation_anchor(
        parent,
        gram,
        roots,
        construction["strict_interior"],
        right,
        right_reflection_sympy,
        second,
    )
    a_only = category_counts["orthogonal_to_a_only_one_based"]
    b_only = category_counts["orthogonal_to_b_only_one_based"]
    common = category_counts["orthogonal_to_both_one_based"]
    neither = category_counts["orthogonal_to_neither_one_based"]
    return {
        "parent": parent,
        "endpoint_wall_indices_one_based": [left + 1, right + 1],
        "endpoint_roots": [list(left_root), list(right_root)],
        "endpoint_squares_and_pairing": [-2, -2, 2],
        "primitive_integral_common_null_direction": list(null_direction),
        "endpoint_walls_are_parallel_doubling_walls": True,
        "right_endpoint_reflection_matrix": [list(row) for row in right_reflection],
        "translation_definition": "T=S_b*S_a (column action)",
        "translation_matrix": [list(row) for row in translation],
        "translation_characteristic_polynomial": "(lambda-1)^15",
        "rank_T_minus_I": 2,
        "rank_T_minus_I_squared": 1,
        "T_minus_I_cubed_is_zero": True,
        "T_minus_I_squared_is_nonzero": True,
        "generated_endpoint_reflection_group": "D_infinity",
        "chamber_maps": {
            "g_(2m)": "T^m",
            "g_(2m+1)": "T^m*S_b",
        },
        "side_wall_categories": categories,
        "side_wall_category_counts": category_counts,
        "all_integer_collision_audit": collision,
        "facet_support": support,
        "facet_count_derivation": {
            "common_orthogonals": f"{common}",
            "a_only": f"{a_only}*(N-floor((N-1)/2))",
            "b_only": f"{b_only}*(N-floor(N/2))",
            "neither": f"{neither}*N",
            "outer_endpoints": "2",
            "formula": f"f_{parent.replace('/', '_')}(N)={expected['formula']} for every N>=1",
            "parity_explanation": (
                "a-only walls identify odd copy 2m+1 with the following even copy; "
                "b-only walls identify even copy 2m with the following odd copy"
            ),
        },
        "direct_exact_wall_union_checks": direct,
        "N_equals_2_frozen_lineage_anchor": anchor,
        "reflection_subgroup_index_in_parent": "N",
        "reflection_subgroup_index_in_H14_60": "2*N",
        "volume": f"Vol(P_{parent.replace('/', '_')}(N))=N*Vol({parent})=2*N*Vol(H14/60)",
        "finite_volume": True,
        "compact": False,
        "noncompact_reason": (
            "the path group has finite index N in the nonuniform parent lattice; "
            "nonuniformity is invariant under finite-index passage"
        ),
        "arithmetic": True,
        "full_lattice_for_every_N": {
            "arithmetic_lattice": "U(2)+E8(-1)+A1(-1)^5",
            "arithmetic_seed_walls_excluding_b_smith_diagonal": arithmetic_smith,
            "physical_NS": "I_(1,14)",
            "physical_seed_walls_excluding_b_smith_diagonal": physical_smith,
            "reason": (
                "every path boundary contains the parent walls other than b; "
                "both the arithmetic and physical 15x(M-1) stocks have Smith 1^15"
            ),
            "effectivity_scope": (
                "primitive integral divisor classes only; transported classes are "
                "not asserted effective, nef, or automorphism walls"
            ),
        },
        "all_facets_irredundant_for_every_N": True,
        "all_exact_checks_pass": True,
    }


def certify(direct_max_n: int = 12) -> dict[str, Any]:
    if not 1 <= direct_max_n <= 100:
        raise ValueError("direct_max_n must lie between 1 and 100")
    hashes = {
        "single_doublings": _sha256(SINGLE),
        "type_vi_source": _sha256(SOURCE),
        "second_generation": _sha256(SECOND),
    }
    if hashes != PINNED:
        raise ArithmeticError(f"a pinned H14 input changed: {hashes}")
    single = json.loads(SINGLE.read_text())
    source = json.loads(SOURCE.read_text())
    second = json.loads(SECOND.read_text())
    if not (
        single["all_checks_pass"]
        and source["all_checks_pass"]
        and second["all_checks_pass"]
    ):
        raise ArithmeticError("a frozen source certificate is not green")
    gram = sp.Matrix(source["arithmetic"]["gram"])
    gram_integral = tuple(tuple(int(value) for value in row) for row in gram.tolist())
    ns_gram = sp.Matrix(source["full_NS_lattice"]["gram"])
    families = {
        parent: _certify_parent(
            parent,
            single["positive_doubles"][parent],
            gram,
            gram_integral,
            ns_gram,
            second,
            direct_max_n,
        )
        for parent in PARENTS
    }
    direct_distinction = all(
        _facet_formula("H14/101", n) - _facet_formula("H14/100", n) == n
        for n in range(1, direct_max_n + 1)
    )
    if not direct_distinction:
        raise ArithmeticError("the two path-ray facet formulas collided")
    return {
        "schema": SCHEMA,
        "artifact_id": "COBLE-VI3-H14-TWO-PARALLEL-PATH-RAYS",
        "status_date": "2026-08-09",
        "direction": "surface_divisor_lattice_to_hyperbolic_polytope",
        "method_scope": {
            "cdd_used": False,
            "extreme_rays_loaded_or_enumerated": False,
            "coxiter_used": False,
            "memory_profile": (
                "15x15 exact matrices, O(M^2) scalar collision equations for "
                "M<=101, and bounded O(NM) wall sets only"
            ),
            "direct_reconstruction_range": f"1<=N<={direct_max_n}",
            "all_N_proof_is_independent_of_direct_range": True,
        },
        "inputs": {
            key: {"path": str(path.relative_to(ROOT)), "sha256": hashes[key]}
            for key, path in {
                "single_doublings": SINGLE,
                "type_vi_source": SOURCE,
                "second_generation": SECOND,
            }.items()
        },
        "theorem": {
            "statement": (
                "For every N>=1, alternating along the explicit parallel wall pair "
                "W12/W20 gives two finite-volume noncompact arithmetic H14 Coxeter "
                "polytopes with respectively 81N+19-floor(N/2) and "
                "82N+19-floor(N/2) facets."
            ),
            "indices_in_H14_60": "both have index 2*N",
            "volumes": "both have volume 2*N times H14/60",
            "two_polytopes_at_each_N_are_nonisometric": (
                "their facet counts differ by exactly N"
            ),
            "all_pairs_across_both_rays_are_nonisometric": (
                "different N have different volume; equal N have different facet count"
            ),
            "N_equals_2_anchors": [
                "COBLE-VI3-D4A02-H14-180",
                "COBLE-VI3-D4B04-H14-182",
            ],
        },
        "families": families,
        "novelty": {
            "new_to_project": (
                "the exact all-N parity-sensitive facet formulas, the exhaustive "
                "all-integer collision proof, and the proof that the two marked "
                "rays are mutually disjoint up to isometry"
            ),
            "existence_scope": (
                "the underlying infinite Allcock families were already known to "
                "exist in the project; this theorem makes two rays explicit"
            ),
            "worldwide_priority": (
                "open/not claimed; Allcock already constructed infinite H14 families"
            ),
        },
        "algebraic_scope": {
            "positive": (
                "all path walls are primitive integral classes spanning the full "
                "I_(1,14) divisor lattice"
            ),
            "not_claimed": (
                "effectivity after reflection, a nef cone, or a surface-automorphism "
                "fundamental domain"
            ),
        },
        "all_checks_pass": all(
            family["all_exact_checks_pass"] for family in families.values()
        ),
    }


def _render(payload: dict[str, Any]) -> str:
    a = payload["families"]["H14/100"]
    b = payload["families"]["H14/101"]
    rows_a = ", ".join(
        f"{row['N']}:{row['facet_count']}"
        for row in a["direct_exact_wall_union_checks"]
    )
    rows_b = ", ".join(
        f"{row['N']}:{row['facet_count']}"
        for row in b["direct_exact_wall_union_checks"]
    )
    return f"""# Two exact H14 parallel-path rays

Status: **all exact checks pass** (2026-08-09).

This is an OOM-safe positive theorem.  It reads the frozen `H14/100` and
`H14/101` walls, but it does **not** call `cdd`, enumerate rays, or invoke
CoxIter.  The certificate is
[`coble_type_vi_h14_parallel_path_families.json`](../data/enumerations/coble_type_vi_h14_parallel_path_families.json),
generated by
[`certify_coble_type_vi_h14_parallel_path_families.py`](../scripts/certify_coble_type_vi_h14_parallel_path_families.py).

## The theorem

Both parents contain the same marked pair W12/W20 with

```text
a^2=b^2=-2,  (a,b)=2,  (a+b)^2=0.
```

Both walls are doubling walls.  For every integer `N>=1`, alternate across
them and take the union of `N` parent chambers.  This gives two
finite-volume, noncompact, arithmetic Coxeter polytopes in `H^14`:

```text
f_100(N) = 81N + 19 - floor(N/2),
f_101(N) = 82N + 19 - floor(N/2).
```

Each reflection subgroup has index `N` in its parent and index `2N` in the
common `H14/60` group.  Its volume is correspondingly `N` parent volumes, or
`2N` `H14/60` volumes.  Thus different values of `N` are distinguished by
volume.  At the same `N`, the two facet counts differ by exactly `N`, so the
two rays are also mutually nonisometric.  In particular this supplies two
explicit, disjoint marked classes at every even subgroup index `2N` over the
surface-derived `H14/60` seed.

At `N=2` the coordinate lists match two already frozen second-generation
objects exactly:

- `H14/100` path: `COBLE-VI3-D4A02-H14-180`;
- `H14/101` path: `COBLE-VI3-D4B04-H14-182`.

## Why the parity term is exact

With column action put `T=S_b S_a` and `K=T-I`.  Exact 15-by-15 arithmetic
gives

```text
charpoly(T)=(lambda-1)^15,
rank(K)=2, rank(K^2)=1,
K^3=0 but K^2!=0.
```

The chamber maps are `g_(2m)=T^m` and `g_(2m+1)=T^m S_b`.  Therefore every
possible equality of transported walls is one of two equations
`T^z r_i = +/- r_j` or `T^z S_b r_i = +/- r_j`, and

```text
T^z = I + zK + z(z-1)K^2/2
```

for every `z` in `Z`.  The script solves one nonzero integer
linear/quadratic coordinate and verifies all 15 coordinates.  This is an
all-integer proof, not a scan over a bounded path length.

For both parents, the 98 or 99 nonendpoint walls split into

| endpoint relation | count in H14/100 | count in H14/101 | distinct images in an N-path |
|---|---:|---:|---|
| orthogonal to both | 14 | 14 | `14` |
| orthogonal only to `a` | 3 | 3 | `3(N-floor((N-1)/2))` |
| orthogonal only to `b` | 4 | 4 | `4(N-floor(N/2))` |
| orthogonal to neither | 77 | 78 | `77N` or `78N` |

There are no cross-wall collisions, no negative-root collisions, and no
endpoint-to-side collisions.  The two outer endpoints remain distinct for
every forward path.  Adding them produces the displayed formulas.  The
asymmetry `3` versus `4` is exactly what creates the `floor(N/2)` phase.

## Facets, volume, and the divisor carrier

Every parent facet has an independently replayed strict timelike projection
witness.  Transporting such a witness into a path copy proves that every
deduplicated boundary mirror is an actual facet; only the two endpoint walls
are used as glues.  Hence the formulas count facets, not merely distinct
normals.

The bounded wall-union regression, which is redundant with the all-`z`
proof, gives

```text
H14/100: {rows_a}
H14/101: {rows_b}
```

Every path boundary contains the parent wall stock with W20 removed.  In
both families this stock has Smith diagonal `1^15`, in both the arithmetic
lattice and the physical `NS = I_(1,14)` realization.  Integral endpoint
reflections therefore transport primitive integral divisor classes and every
descendant retains the full divisor lattice.  This remains a weak carrier
statement: reflected classes are not asserted effective, the path chamber is
not asserted nef, and no surface-automorphism fundamental domain is claimed.

Finite volume follows because the path is a finite union of `N` parent
chambers.  Noncompactness is not inferred from a partial ray count: the path
group has finite index in the nonuniform parent lattice, and nonuniformity is
invariant under finite-index passage.

## Claim boundary

The exact parity formulas, all-integer collision classification, and mutual
separation of the two rays are **new to this project**.  The existence of the
underlying infinite Allcock families was already recorded.  Worldwide
priority is **not claimed**; Allcock already constructed infinite families in
dimension 14.

## Replay

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_coble_type_vi_h14_parallel_path_families.py --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_coble_type_vi_h14_parallel_path_families.py
```
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--direct-max-N", type=int, default=12)
    parser.add_argument(
        "--check-archive",
        action="store_true",
        help="rebuild in memory and compare with both frozen artifacts",
    )
    args = parser.parse_args()
    payload = certify(args.direct_max_N)
    rendered_json = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    rendered_report = _render(payload)
    if args.check_archive:
        if args.output.read_text() != rendered_json:
            raise SystemExit(f"archive mismatch: {args.output}")
        if args.report.read_text() != rendered_report:
            raise SystemExit(f"report mismatch: {args.report}")
    else:
        args.output.write_text(rendered_json)
        args.report.write_text(rendered_report)
    print(json.dumps({
        "all_checks_pass": payload["all_checks_pass"],
        "archive_check": args.check_archive,
        "output": str(args.output),
        "report": str(args.report),
        "direct_max_N": args.direct_max_N,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
