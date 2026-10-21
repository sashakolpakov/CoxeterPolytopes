#!/usr/bin/env python3
"""Gram-only hostile audit of the A18_4_0_C01 T2 lattice precompiler.

The sole mathematical input is the pinned 14 by 14 facet Gram row in
``alexeev_k3_h9_14_factory.json``.  This script does not read or import the
C01/C07 strong-promotion package or a production lattice-precompiler module.

It independently enumerates primitive minimal even crystallographic normal
scalings, the two one-/two-sheet-compatible lanes, every reflection-stable
even overlattice in those lanes, and the finite IC0 obstruction gates.  For
the unique survivor it gives literal U + E8(-2), rank-18 glue, involution,
24-root, and fold-back matrices.  No datum-specific IC2 identification with
the existing C01 chamber/action is inferred from equal abstract invariants.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict, deque
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/alexeev_k3_h9_14_factory.json"
OUTPUT = (
    ROOT
    / "data/enumerations/a18_c01_t2_lattice_precompiler_hostile_audit.json"
)

SCHEMA = "coxeter.a18_c01.t2_lattice_precompiler.hostile_audit.v1"
STATUS_DATE = "2026-08-08"
SOURCE_ROW_ID = "A18_4_0_C01"
PINNED_FACET_GRAM_SHA256 = (
    "6979b16de7583c20277077a142df863e9f9a7209c443e9252505234b15d66c56"
)

Matrix = sp.Matrix
BitBasis = tuple[int, ...]


def _sha256_json(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _fraction(value: sp.Rational | Fraction | int) -> str:
    value = sp.Rational(value)
    return str(int(value)) if value.q == 1 else f"{value.p}/{value.q}"


def _payload_matrix(matrix: Matrix) -> list[list[str]]:
    return [[_fraction(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _integer_matrix(matrix: Matrix) -> Matrix:
    if not all(sp.denom(value) == 1 for value in matrix):
        raise ArithmeticError("matrix is not integral")
    return matrix.applyfunc(int)


def _smith_absolute(matrix: Matrix) -> list[int]:
    smith = smith_normal_form(_integer_matrix(matrix), domain=ZZ)
    return [
        abs(int(smith[i, i]))
        for i in range(min(smith.rows, smith.cols))
        if smith[i, i]
    ]


def _content(vector: Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector if value))


def _inertia(matrix: Matrix) -> tuple[int, int, int]:
    """Exact symmetric Gaussian elimination with 1x1/2x2 pivots."""

    work = Matrix(matrix)
    positive = negative = zero = 0
    while work.rows:
        diagonal = next((i for i in range(work.rows) if work[i, i]), None)
        if diagonal is not None:
            if diagonal:
                work.row_swap(0, diagonal)
                work.col_swap(0, diagonal)
            pivot = work[0, 0]
            positive += 1 if pivot > 0 else 0
            negative += 1 if pivot < 0 else 0
            if work.rows == 1:
                break
            column = work[1:, 0]
            work = (work[1:, 1:] - column * column.T / pivot).applyfunc(sp.cancel)
            continue
        off = next(
            (
                (i, j)
                for i in range(work.rows)
                for j in range(i + 1, work.cols)
                if work[i, j]
            ),
            None,
        )
        if off is None:
            zero += work.rows
            break
        i, j = off
        work.row_swap(0, i)
        work.col_swap(0, i)
        if j == 0:
            j = i
        work.row_swap(1, j)
        work.col_swap(1, j)
        block = work[:2, :2]
        # With every diagonal zero, det(block) is strictly negative.
        if block.det() >= 0:
            raise ArithmeticError("unexpected exact inertia pivot")
        positive += 1
        negative += 1
        if work.rows == 2:
            break
        cross = work[2:, :2]
        work = (work[2:, 2:] - cross * block.inv() * cross.T).applyfunc(sp.cancel)
    return positive, negative, zero


def _load_pinned_gram() -> Matrix:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    rows = [row for row in source["constructions"] if row.get("id") == SOURCE_ROW_ID]
    if len(rows) != 1:
        raise ArithmeticError("the pinned C01 construction row is not unique")
    # Deliberately read no lattice/action data from the row: the Gram is the
    # only mathematical input to this hostile precompiler.
    raw = rows[0]["section"]["restricted_mirrors"]["facet_gram"]
    if _sha256_json(raw) != PINNED_FACET_GRAM_SHA256:
        raise ArithmeticError("the pinned C01 facet Gram changed")
    gram = Matrix([[int(value) for value in row] for row in raw])
    if gram.shape != (14, 14) or gram != gram.T or gram.rank() != 10:
        raise ArithmeticError("the pinned facet Gram shape/rank changed")
    if _inertia(gram) != (1, 9, 4):
        raise ArithmeticError("the pinned facet Gram inertia changed")
    return gram


def _gf2_basis(vectors: Iterable[int]) -> BitBasis:
    pivots: dict[int, int] = {}
    for value in vectors:
        value = int(value)
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
                continue
            pivots[pivot] = value
            for other in tuple(pivots):
                if other != pivot and ((pivots[other] >> pivot) & 1):
                    pivots[other] ^= value
            break
    return tuple(pivots[pivot] for pivot in sorted(pivots, reverse=True))


def _gf2_span(basis: Sequence[int]) -> set[int]:
    span = {0}
    for value in basis:
        span |= {old ^ value for old in tuple(span)}
    return span


def _gf2_rank(vectors: Iterable[int]) -> int:
    return len(_gf2_basis(vectors))


def _column_bits(matrix: Matrix, column: int) -> int:
    return sum((int(matrix[row, column]) & 1) << row for row in range(matrix.rows))


def _minimum_even_factor(coefficients: Sequence[sp.Rational]) -> sp.Rational:
    """Smallest q>0 with q*c integral for every supplied rational c."""

    denominator = math.lcm(*(int(value.q) for value in coefficients if value))
    divisor = math.gcd(
        *(
            abs(int(value.p)) * (denominator // int(value.q))
            for value in coefficients
            if value
        )
    )
    return sp.Rational(denominator, divisor)


def _normal_scaling_rays(gram: Matrix) -> list[dict[str, Any]]:
    """Enumerate all crystallographic rays, fixing the first short scale."""

    short = [i for i in range(14) if gram[i, i] == -2]
    long = [i for i in range(14) if gram[i, i] == -4]
    if short != [0, 4, 8, 13] or len(long) != 10:
        raise ArithmeticError("the short/long facet partition changed")

    adjacency = {index: set() for index in long}
    for left in long:
        for right in long:
            if left < right and gram[left, right]:
                adjacency[left].add(right)
                adjacency[right].add(left)
    components: list[list[int]] = []
    unseen = set(long)
    while unseen:
        stack = [min(unseen)]
        unseen.remove(stack[0])
        component: list[int] = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor in sorted(adjacency[vertex]):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    components.sort()
    if list(map(len, components)) != [5, 5]:
        raise ArithmeticError("the two m3 long components changed")
    if not all(
        sum(gram[node, neighbor] != 0 for neighbor in component) == 1
        for node in short
        for component in components
    ):
        raise ArithmeticError("a short node no longer meets each long component once")

    exponent_rows: set[tuple[int, ...]] = set()
    # Along an m3 edge the exponents agree.  Along an m4 short--long edge,
    # e_long-e_short is 0 or -1.  Gauge-fix e(short[0])=0.
    for first_exponent, second_exponent in itertools.product((0, -1), repeat=2):
        allowed: dict[int, set[int]] = {}
        for node in short:
            allowed[node] = {first_exponent, first_exponent + 1} & {
                second_exponent,
                second_exponent + 1,
            }
        if 0 not in allowed[short[0]]:
            continue
        for choices in itertools.product(*(sorted(allowed[node]) for node in short[1:])):
            exponents = [0] * 14
            for node in components[0]:
                exponents[node] = first_exponent
            for node in components[1]:
                exponents[node] = second_exponent
            for node, exponent in zip(short[1:], choices, strict=True):
                exponents[node] = exponent
            exponent_rows.add(tuple(exponents))

    rows: list[dict[str, Any]] = []
    for exponents in sorted(exponent_rows):
        relative = [sp.Rational(2) ** exponent for exponent in exponents]
        constraints: list[sp.Rational] = []
        for i in range(14):
            for j in range(i, 14):
                if not gram[i, j]:
                    continue
                coefficient = sp.Rational(gram[i, j]) * relative[i] * relative[j]
                constraints.append(abs(coefficient / 2 if i == j else coefficient))
        common_square = _minimum_even_factor(constraints)
        diagonal = sp.diag(*relative)
        scaled = (common_square * diagonal * gram * diagonal).applyfunc(sp.cancel)
        scaled = _integer_matrix(scaled)
        cartan = Matrix(
            14,
            14,
            lambda i, j: 0
            if not scaled[i, j]
            else sp.Rational(2 * scaled[i, j], scaled[i, i]),
        )
        if not all(sp.denom(value) == 1 for value in cartan):
            raise ArithmeticError("a purported scaling is not crystallographic")
        if not all(scaled[i, i] % 2 == 0 for i in range(14)):
            raise ArithmeticError("a purported scaling is not even")
        norms = [-int(scaled[i, i]) for i in range(14)]
        rows.append(
            {
                "exponents_base_2_with_W1_gauge_zero": list(exponents),
                "minimal_common_squared_factor": _fraction(common_square),
                "squared_scale_factors": [
                    _fraction(common_square * value * value) for value in relative
                ],
                "normal_norms_absolute": norms,
                "one_two_sheet_compatible": set(norms) <= {2, 4},
                "scaled_gram": scaled,
            }
        )
    if len(rows) != 18:
        raise ArithmeticError("the complete primitive scaling count changed")
    return rows


def _generated_lattice(scaled_gram: Matrix) -> dict[str, Matrix | list[int]]:
    """Build the saturated coordinate lattice generated by all 14 normals."""

    basis_indices = list(scaled_gram.rref()[1])
    basis_gram = scaled_gram.extract(basis_indices, basis_indices)
    if len(basis_indices) != 10 or not basis_gram.det():
        raise ArithmeticError("could not choose a nonsingular normal basis")
    normal_coordinates = basis_gram.inv() * scaled_gram.extract(basis_indices, range(14))
    denominator = math.lcm(*(int(sp.denom(value)) for value in normal_coordinates))
    column_lattice = _integer_matrix(normal_coordinates * denominator)
    coordinate_basis = sp.Rational(1, denominator) * hermite_normal_form(column_lattice)
    lattice_gram = _integer_matrix(coordinate_basis.T * basis_gram * coordinate_basis)
    normals = _integer_matrix(coordinate_basis.inv() * normal_coordinates)
    if normals.T * lattice_gram * normals != scaled_gram:
        raise ArithmeticError("the generated lattice does not recover its Gram")
    return {
        "basis_indices": basis_indices,
        "gram": lattice_gram,
        "normals": normals,
    }


def _quotient_mod_two(gram: Matrix) -> dict[str, Any]:
    """Represent an elementary 2-discriminant group as F2^a."""

    rank = gram.rows
    image_basis = _gf2_basis(_column_bits(gram, column) for column in range(rank))

    def reduce_value(value: int) -> int:
        for basis_value in image_basis:
            pivot = basis_value.bit_length() - 1
            if (value >> pivot) & 1:
                value ^= basis_value
        return value

    quotient_basis: list[int] = []
    known = {0}
    coordinate = {0: 0}
    for bit_index in range(rank):
        value = reduce_value(1 << bit_index)
        if value in known:
            continue
        old = tuple(known)
        bit = 1 << len(quotient_basis)
        quotient_basis.append(value)
        for previous in old:
            new = previous ^ value
            known.add(new)
            coordinate[new] = coordinate[previous] ^ bit
    classes = sorted(known)
    if len(classes) != 1 << len(quotient_basis):
        raise ArithmeticError("failed to coordinatize a 2-discriminant group")
    by_coordinate: list[int] = [0] * len(classes)
    for value, index in coordinate.items():
        by_coordinate[index] = value
    inverse = gram.inv()
    q4: list[int] = []
    for value in by_coordinate:
        covector = Matrix([(value >> i) & 1 for i in range(rank)])
        doubled_q = 2 * (covector.T * inverse * covector)[0]
        if sp.denom(doubled_q) != 1:
            raise ArithmeticError("the discriminant exponent exceeds two")
        q4.append(int(doubled_q) % 4)
    return {
        "reduce_value": reduce_value,
        "coordinate": coordinate,
        "basis_classes": tuple(quotient_basis),
        "class_by_coordinate": tuple(by_coordinate),
        "q4": tuple(q4),
    }


def _reflection_matrices(gram: Matrix, normals: Matrix) -> tuple[Matrix, ...]:
    identity = sp.eye(gram.rows)
    reflections: list[Matrix] = []
    for column in range(normals.cols):
        root = normals[:, column]
        square = (root.T * gram * root)[0]
        reflection = identity - sp.Rational(2, square) * root * (root.T * gram)
        reflections.append(_integer_matrix(reflection))
    return tuple(reflections)


def _stable_even_subgroups(gram: Matrix, normals: Matrix) -> tuple[list[BitBasis], dict[str, Any]]:
    quotient = _quotient_mod_two(gram)
    basis_classes = quotient["basis_classes"]
    reduce_value = quotient["reduce_value"]
    coordinate = quotient["coordinate"]
    rank = gram.rows
    action_generators: list[tuple[int, ...]] = []
    for reflection in _reflection_matrices(gram, normals):
        dual_action = reflection.inv().T
        images: list[int] = []
        for value in basis_classes:
            covector = Matrix([(value >> i) & 1 for i in range(rank)])
            image = dual_action * covector
            raw = sum((int(image[i]) & 1) << i for i in range(rank))
            images.append(coordinate[reduce_value(raw)])
        action_generators.append(tuple(images))

    def act(generator: Sequence[int], value: int) -> int:
        image = 0
        for index, basis_image in enumerate(generator):
            if (value >> index) & 1:
                image ^= basis_image
        return image

    def invariant_closure(seed: Sequence[int]) -> BitBasis:
        basis = _gf2_basis(seed)
        while True:
            expanded = _gf2_basis(
                (
                    *basis,
                    *(
                        act(generator, value)
                        for generator in action_generators
                        for value in basis
                    ),
                )
            )
            if expanded == basis:
                return basis
            basis = expanded

    submodules: set[BitBasis] = {()}
    queue: deque[BitBasis] = deque([()])
    group_size = 1 << len(basis_classes)
    while queue:
        basis = queue.popleft()
        span = _gf2_span(basis)
        for value in range(1, group_size):
            if value in span:
                continue
            extension = invariant_closure((*basis, value))
            if extension not in submodules:
                submodules.add(extension)
                queue.append(extension)
    isotropic = [
        basis
        for basis in submodules
        if all(quotient["q4"][value] == 0 for value in _gf2_span(basis))
    ]
    isotropic.sort(key=lambda basis: (len(basis), basis))
    quotient["action_generators"] = tuple(action_generators)
    quotient["invariant_submodule_count"] = len(submodules)
    return isotropic, quotient


def _overlattice(
    gram: Matrix,
    normals: Matrix,
    subgroup: BitBasis,
    quotient: dict[str, Any],
) -> dict[str, Any]:
    inverse = gram.inv()
    extras: list[Matrix] = []
    for coordinate in subgroup:
        value = quotient["class_by_coordinate"][coordinate]
        covector = Matrix([(value >> i) & 1 for i in range(gram.rows)])
        extras.append(inverse * covector)
    generators = sp.eye(gram.rows)
    if extras:
        generators = generators.row_join(Matrix.hstack(*extras))
    denominator = math.lcm(*(int(sp.denom(value)) for value in generators))
    basis = sp.Rational(1, denominator) * hermite_normal_form(
        _integer_matrix(generators * denominator)
    )
    over_gram = _integer_matrix(basis.T * gram * basis)
    over_normals = _integer_matrix(basis.inv() * normals)
    if over_normals.T * over_gram * over_normals != normals.T * gram * normals:
        raise ArithmeticError("an overlattice changed the marked normal Gram")
    if not all(over_gram[i, i] % 2 == 0 for i in range(over_gram.rows)):
        raise ArithmeticError("an enumerated overlattice is not even")
    reflections = _reflection_matrices(over_gram, over_normals)
    if len(reflections) != 14:
        raise ArithmeticError("an overlattice lost a marked reflection")
    if not all(_content(over_normals[:, i]) == 1 for i in range(14)):
        raise ArithmeticError("a marked normal became imprimitive")

    smith = _smith_absolute(over_gram)
    if not all(value in (1, 2) for value in smith):
        raise ArithmeticError("an even stable overlattice ceased to be 2-elementary")
    discr = _quotient_mod_two(over_gram)
    length = sum(value == 2 for value in smith)
    if length != len(discr["basis_classes"]):
        raise ArithmeticError("Smith and mod-two discriminant lengths disagree")

    norms = [-int((over_normals[:, i].T * over_gram * over_normals[:, i])[0]) for i in range(14)]
    long_indices = [i for i, norm in enumerate(norms) if norm == 4]
    long_classes = [
        sum((int(over_normals[row, index]) & 1) << row for row in range(10))
        for index in long_indices
    ]
    repeated: list[dict[str, int]] = []
    for left_position, left in enumerate(long_indices):
        for right_position in range(left_position):
            right = long_indices[right_position]
            if long_classes[left_position] == long_classes[right_position]:
                repeated.append(
                    {
                        "left_wall_one_based": right + 1,
                        "right_wall_one_based": left + 1,
                        "folded_pairing": int(
                            (over_normals[:, right].T * over_gram * over_normals[:, left])[0]
                        ),
                    }
                )
    bilinear_rows: list[int] = []
    radical_facets: list[int] = []
    for left_position, left in enumerate(long_indices):
        row_bits = 0
        for right_position, right in enumerate(long_indices):
            pairing = int(
                (over_normals[:, left].T * over_gram * over_normals[:, right])[0]
            )
            if pairing % 2:
                raise ArithmeticError("a long class is not dual")
            row_bits |= ((pairing // 2) & 1) << right_position
        bilinear_rows.append(row_bits)
        if not row_bits:
            radical_facets.append(left + 1)
    long_span_rank = _gf2_rank(long_classes)
    long_radical_dimension = long_span_rank - _gf2_rank(bilinear_rows)
    if long_radical_dimension < 0:
        raise ArithmeticError("invalid long-class radical dimension")

    return {
        "subgroup_basis": subgroup,
        "glue_generators_in_root_lattice_basis": extras,
        "basis_in_root_lattice_tensor_Q": basis,
        "gram": over_gram,
        "normals": over_normals,
        "index": 1 << len(subgroup),
        "determinant_absolute": abs(int(over_gram.det())),
        "smith_diagonal_absolute": smith,
        "discriminant_length": length,
        "discriminant_delta": int(any(value % 2 for value in discr["q4"])),
        "discriminant_q4_histogram": dict(sorted(Counter(discr["q4"]).items())),
        "normal_norms_absolute": norms,
        "long_wall_indices_one_based": [index + 1 for index in long_indices],
        "long_class_span_rank": long_span_rank,
        "long_class_radical_dimension": long_radical_dimension,
        "radical_long_wall_classes_one_based": radical_facets,
        "repeated_long_class_pairs": repeated,
        "dual_rank_eight_length_gate": length <= 8,
        "no_long_facet_class_is_radical": not radical_facets,
        "no_repeated_long_classes": not repeated,
    }


def _public_overlattice(row: dict[str, Any], number: int) -> dict[str, Any]:
    return {
        "number": number,
        "glue_dimension": len(row["subgroup_basis"]),
        "index": row["index"],
        "glue_generators_in_root_lattice_basis": [
            [_fraction(value) for value in vector]
            for vector in row["glue_generators_in_root_lattice_basis"]
        ],
        "determinant_absolute": row["determinant_absolute"],
        "smith_diagonal_absolute": row["smith_diagonal_absolute"],
        "discriminant_length": row["discriminant_length"],
        "discriminant_delta": row["discriminant_delta"],
        "discriminant_q4_histogram": {
            str(key): value for key, value in row["discriminant_q4_histogram"].items()
        },
        "normal_norms_absolute": row["normal_norms_absolute"],
        "long_wall_indices_one_based": row["long_wall_indices_one_based"],
        "long_class_span_rank": row["long_class_span_rank"],
        "long_class_radical_dimension": row["long_class_radical_dimension"],
        "radical_long_wall_classes_one_based": row[
            "radical_long_wall_classes_one_based"
        ],
        "repeated_long_class_pairs": row["repeated_long_class_pairs"],
        "gates": {
            "dual_rank_eight_discriminant_length_at_most_eight": row[
                "dual_rank_eight_length_gate"
            ],
            "no_long_facet_class_in_long_class_radical": row[
                "no_long_facet_class_is_radical"
            ],
            "no_repeated_long_classes": row["no_repeated_long_classes"],
        },
    }


def _candidate_certificate(candidate: dict[str, Any], pinned_gram: Matrix) -> dict[str, Any]:
    """Give exact U+E8(-2), rank-18 glue, action, and fold-back witnesses."""

    gram = candidate["gram"]
    normals = candidate["normals"]
    # Independently found unimodular splitting witness.  Its verification,
    # not a classification label, proves the claimed direct sum.
    hyperbolic_left = Matrix([1, -2, -1, -3, -2, 1, 2, 0, 2, 3])
    hyperbolic_right = Matrix([2, -2, -1, -3, -2, 1, 2, 0, 2, 3])
    complement_columns: list[Matrix] = []
    for parameter in range(8):
        vector = sp.zeros(10, 1)
        if parameter < 5:
            vector[[1, 2, 3, 5, 6][parameter]] = 1
        elif parameter == 5:
            vector[7] = 1
            vector[4] = -2
        elif parameter == 6:
            vector[8] = 1
            vector[4] = 1
        else:
            vector[0] = 1
            vector[9] = 2
            vector[4] = -3
        complement_columns.append(vector)
    splitting = Matrix.hstack(
        hyperbolic_left,
        hyperbolic_right,
        *complement_columns,
    )
    if abs(int(splitting.det())) != 1:
        raise ArithmeticError("the U+E8(-2) splitting is not unimodular")
    block = _integer_matrix(splitting.T * gram * splitting)
    hyperbolic = Matrix([[0, 1], [1, 0]])
    if block[:2, :2] != hyperbolic or block[:2, 2:] != sp.zeros(2, 8):
        raise ArithmeticError("the displayed U summand did not split")
    if any(block[i, j] % 2 for i in range(2, 10) for j in range(2, 10)):
        raise ArithmeticError("the complementary summand is not twice integral")
    e8_negative = _integer_matrix(block[2:, 2:] / 2)
    if not (
        e8_negative.det() == 1
        and _inertia(e8_negative) == (0, 8, 0)
        and all(e8_negative[i, i] % 2 == 0 for i in range(8))
    ):
        raise ArithmeticError("the complement is not even unimodular E8(-1)")
    block_normals = _integer_matrix(splitting.inv() * normals)

    # Select exact E8(-1) roots in the ten long discriminant classes.  Eight
    # are already roots; two exact congruent representatives are supplied and
    # checked.  No production root/action data enters this choice.
    exceptional_anti = {
        7: Matrix([0, 0, 0, 0, 0, 1, 2, 0]),
        12: Matrix([-3, -2, -4, 1, 2, -1, 2, 2]),
    }
    anti_roots: dict[int, Matrix] = {}
    direct_sum_gram = sp.diag(block, 2 * e8_negative)
    direct_sum_roots: list[Matrix] = []
    root_orbits: list[list[int]] = []
    glue_vectors: list[Matrix] = []
    for wall in range(14):
        fixed_normal = block_normals[:, wall]
        norm = -int((fixed_normal.T * block * fixed_normal)[0])
        if norm == 2:
            source_root = Matrix.vstack(fixed_normal, sp.zeros(8, 1))
            direct_sum_roots.append(source_root)
            root_orbits.append([len(direct_sum_roots) - 1])
            continue
        if norm != 4:
            raise ArithmeticError("the candidate has an invalid sheet norm")
        e8_component = fixed_normal[2:, 0]
        anti_root = exceptional_anti.get(wall, e8_component)
        if (anti_root.T * e8_negative * anti_root)[0] != -2:
            raise ArithmeticError("an anti representative is not an E8 root")
        if any(int(anti_root[i] - e8_component[i]) % 2 for i in range(8)):
            raise ArithmeticError("an anti root has the wrong discriminant class")
        if any(int(fixed_normal[i]) % 2 for i in range(2)):
            raise ArithmeticError("a long normal has a nondivisible U coordinate")
        anti_roots[wall] = anti_root
        fixed_part = Matrix.vstack(fixed_normal, sp.zeros(8, 1))
        anti_part = Matrix.vstack(sp.zeros(10, 1), anti_root)
        plus = (fixed_part + anti_part) / 2
        minus = (fixed_part - anti_part) / 2
        # These are deliberately half-integral in S direct-sum E8(-2); they
        # become integral generators of the rank-18 glued lattice below.
        direct_sum_roots.extend((plus, minus))
        root_orbits.append([len(direct_sum_roots) - 2, len(direct_sum_roots) - 1])
        glue_vectors.append(plus)

    long_fixed_classes = [
        sum((int(block_normals[row, wall]) & 1) << (row - 2) for row in range(2, 10))
        for wall in sorted(anti_roots)
    ]
    long_anti_classes = [
        sum((int(anti_roots[wall][row]) & 1) << row for row in range(8))
        for wall in sorted(anti_roots)
    ]
    if long_fixed_classes != long_anti_classes:
        raise ArithmeticError("the fixed/anti long-class markings disagree")
    if _gf2_rank(long_fixed_classes) != 6:
        raise ArithmeticError("the candidate long-class glue rank changed")

    generators = sp.eye(18).row_join(Matrix.hstack(*glue_vectors))
    denominator = math.lcm(*(int(sp.denom(value)) for value in generators))
    ambient_basis = sp.Rational(1, denominator) * hermite_normal_form(
        _integer_matrix(generators * denominator)
    )
    ambient_index = sp.Rational(1, abs(ambient_basis.det()))
    if ambient_index != 64:
        raise ArithmeticError("the rank-18 long-class glue index changed")
    ambient_gram = _integer_matrix(ambient_basis.T * direct_sum_gram * ambient_basis)
    ambient_roots = _integer_matrix(
        ambient_basis.inv() * Matrix.hstack(*direct_sum_roots)
    )
    if not all(ambient_gram[i, i] % 2 == 0 for i in range(18)):
        raise ArithmeticError("the rank-18 glued lattice is not even")
    if _inertia(ambient_gram) != (1, 17, 0):
        raise ArithmeticError("the rank-18 glued lattice has wrong signature")
    ambient_smith = _smith_absolute(ambient_gram)
    ambient_discr = _quotient_mod_two(ambient_gram)
    if not (
        ambient_gram.det() == -16
        and ambient_smith == [1] * 14 + [2] * 4
        and not any(value % 2 for value in ambient_discr["q4"])
    ):
        raise ArithmeticError("the rank-18 (18,4,0) invariants changed")

    direct_involution = sp.diag(sp.eye(10), -sp.eye(8))
    involution = _integer_matrix(ambient_basis.inv() * direct_involution * ambient_basis)
    if not (
        involution * involution == sp.eye(18)
        and involution.T * ambient_gram * involution == ambient_gram
        and involution.trace() == 2
    ):
        raise ArithmeticError("the rank-18 involution failed")
    if not all(
        involution * ambient_roots[:, orbit[0]] == ambient_roots[:, orbit[-1]]
        for orbit in root_orbits
    ):
        raise ArithmeticError("the source roots do not have the promised orbits")

    root_gram = ambient_roots.T * ambient_gram * ambient_roots
    if not all(root_gram[i, i] == -2 for i in range(root_gram.rows)):
        raise ArithmeticError("a lifted source root is not square minus two")
    off_diagonal = [
        int(root_gram[i, j])
        for i in range(root_gram.rows)
        for j in range(i)
    ]
    if min(off_diagonal) < 0 or set(off_diagonal) - {0, 1}:
        raise ArithmeticError("the lifted source roots are not acute")
    root_hnf = hermite_normal_form(ambient_roots)
    if root_hnf.shape != (18, 18) or abs(int(root_hnf.det())) != 1:
        raise ArithmeticError("the 24 lifted roots do not span the ambient lattice")

    folded_roots = [
        ambient_roots[:, orbit[0]]
        if len(orbit) == 1
        else ambient_roots[:, orbit[0]] + ambient_roots[:, orbit[1]]
        for orbit in root_orbits
    ]
    folded = Matrix.hstack(*folded_roots)
    if folded.T * ambient_gram * folded != pinned_gram:
        raise ArithmeticError("rank-18 fold-back does not recover the pinned Gram")

    fixed_basis_direct = Matrix.vstack(sp.eye(10), sp.zeros(8, 10))
    anti_basis_direct = Matrix.vstack(sp.zeros(10, 8), sp.eye(8))
    fixed_basis = _integer_matrix(ambient_basis.inv() * fixed_basis_direct)
    anti_basis = _integer_matrix(ambient_basis.inv() * anti_basis_direct)
    if not (
        fixed_basis.T * ambient_gram * fixed_basis == block
        and anti_basis.T * ambient_gram * anti_basis == 2 * e8_negative
        and _smith_absolute(fixed_basis) == [1] * 10
        and _smith_absolute(anti_basis) == [1] * 8
    ):
        raise ArithmeticError("the fixed/anti eigensublattices are not primitive")
    if not (
        involution * fixed_basis == fixed_basis
        and involution * anti_basis == -anti_basis
    ):
        raise ArithmeticError("the eigensublattice signs changed")

    return {
        "fixed_lattice": {
            "isometry_class": "U+E8(-2)",
            "splitting_basis_columns_in_candidate_basis": _payload_matrix(splitting),
            "block_gram": _payload_matrix(block),
            "U_gram": _payload_matrix(hyperbolic),
            "E8_minus_one_gram": _payload_matrix(e8_negative),
            "splitting_determinant_absolute": abs(int(splitting.det())),
            "E8_minus_one_even_unimodular_negative_definite": True,
        },
        "long_class_glue": {
            "long_wall_indices_one_based": [wall + 1 for wall in sorted(anti_roots)],
            "fixed_long_classes_as_bits": long_fixed_classes,
            "anti_long_classes_as_bits": long_anti_classes,
            "class_span_rank": _gf2_rank(long_fixed_classes),
            "anti_E8_root_representatives": {
                str(wall + 1): [int(value) for value in anti_roots[wall]]
                for wall in sorted(anti_roots)
            },
            "glue_index": int(ambient_index),
        },
        "rank18": {
            "abstract_2_elementary_invariants": [18, 4, 0],
            "gram": _payload_matrix(ambient_gram),
            "determinant": int(ambient_gram.det()),
            "smith_diagonal_absolute": ambient_smith,
            "signature": [1, 17],
            "ambient_basis_columns_in_fixed_plus_anti_tensor_Q": _payload_matrix(
                ambient_basis
            ),
            "involution": _payload_matrix(involution),
            "involution_order": 2,
            "involution_trace": int(involution.trace()),
            "fixed_rank": 10,
            "anti_rank": 8,
            "fixed_basis_columns": _payload_matrix(fixed_basis),
            "anti_basis_columns": _payload_matrix(anti_basis),
            "fixed_and_anti_sublattices_primitive": True,
            "anti_lattice": "E8(-2)",
            "anti_lattice_rootless": True,
        },
        "fold_back": {
            "source_root_count": ambient_roots.cols,
            "source_roots_span_rank18_with_index_one": True,
            "source_roots": _payload_matrix(ambient_roots),
            "source_root_gram_pairing_histogram": {
                str(key): value
                for key, value in sorted(Counter(off_diagonal).items())
            },
            "source_root_orbits_one_based": [
                [index + 1 for index in orbit] for orbit in root_orbits
            ],
            "one_sheet_orbit_count": sum(len(orbit) == 1 for orbit in root_orbits),
            "two_sheet_orbit_count": sum(len(orbit) == 2 for orbit in root_orbits),
            "folded_roots": _payload_matrix(folded),
            "folded_gram_sha256": _sha256_json(_payload_matrix(folded.T * ambient_gram * folded)),
            "folded_gram_equals_pinned_facet_gram": True,
            "all_source_roots_square_minus_two": True,
            "all_source_root_pairings_acute_zero_or_one": True,
            "involution_has_promised_root_orbits": True,
        },
    }


def audit() -> dict[str, Any]:
    pinned_gram = _load_pinned_gram()
    scaling_rows = _normal_scaling_rays(pinned_gram)
    compatible = [row for row in scaling_rows if row["one_two_sheet_compatible"]]
    if len(compatible) != 2:
        raise ArithmeticError("the one-/two-sheet lane count changed")
    compatible.sort(key=lambda row: row["normal_norms_absolute"])

    lane_payloads: list[dict[str, Any]] = []
    internal_overlattices: list[tuple[str, dict[str, Any]]] = []
    for row in compatible:
        label = (
            "displayed"
            if row["scaled_gram"] == pinned_gram
            else "fully_reversed"
        )
        generated = _generated_lattice(row["scaled_gram"])
        lattice_gram = generated["gram"]
        normals = generated["normals"]
        subgroups, quotient = _stable_even_subgroups(lattice_gram, normals)
        overlattices = [
            _overlattice(lattice_gram, normals, subgroup, quotient)
            for subgroup in subgroups
        ]
        expected_count = 16 if label == "displayed" else 5
        if len(overlattices) != expected_count:
            raise ArithmeticError(f"the {label} stable-even count changed")
        for over in overlattices:
            internal_overlattices.append((label, over))
        lane_payloads.append(
            {
                "label": label,
                "normal_norms_absolute": row["normal_norms_absolute"],
                "scaled_gram": _payload_matrix(row["scaled_gram"]),
                "root_lattice": {
                    "gram": _payload_matrix(lattice_gram),
                    "normal_coordinates": _payload_matrix(normals),
                    "determinant_absolute": abs(int(lattice_gram.det())),
                    "smith_diagonal_absolute": _smith_absolute(lattice_gram),
                    "signature": [1, 9],
                    "invariant_discriminant_submodule_count": quotient[
                        "invariant_submodule_count"
                    ],
                },
                "reflection_stable_even_overlattice_count": len(overlattices),
                "overlattices": [
                    _public_overlattice(over, number)
                    for number, over in enumerate(overlattices, 1)
                ],
            }
        )

    survivors = [
        (label, over)
        for label, over in internal_overlattices
        if over["dual_rank_eight_length_gate"]
        and over["no_long_facet_class_is_radical"]
        and over["no_repeated_long_classes"]
    ]
    if len(survivors) != 1 or survivors[0][0] != "displayed":
        raise ArithmeticError("the finite IC0 obstruction survivor is not unique")
    candidate = survivors[0][1]
    if not (
        candidate["determinant_absolute"] == 256
        and candidate["smith_diagonal_absolute"] == [1, 1] + [2] * 8
        and candidate["discriminant_length"] == 8
        and candidate["discriminant_delta"] == 0
        and candidate["long_class_span_rank"] == 6
    ):
        raise ArithmeticError("the unique candidate invariants changed")
    certificate = _candidate_certificate(candidate, pinned_gram)

    displayed = next(row for row in lane_payloads if row["label"] == "displayed")
    reversed_lane = next(row for row in lane_payloads if row["label"] == "fully_reversed")
    displayed_base = displayed["overlattices"][0]
    reversed_base = reversed_lane["overlattices"][0]
    if not (
        displayed_base["discriminant_length"] == 10
        and not displayed_base["gates"][
            "dual_rank_eight_discriminant_length_at_most_eight"
        ]
        and reversed_base["radical_long_wall_classes_one_based"] == [1, 5, 9, 14]
        and all(
            pair["folded_pairing"] == 0
            for row in lane_payloads
            for over in row["overlattices"]
            for pair in over["repeated_long_class_pairs"]
        )
    ):
        raise ArithmeticError("an advertised obstruction witness changed")

    checks = {
        "only_pinned_C01_facet_Gram_used_as_mathematical_input": True,
        "facet_Gram_rank10_and_inertia_1_9_4": True,
        "all_18_primitive_minimal_even_crystallographic_scaling_rays_enumerated": len(
            scaling_rows
        )
        == 18,
        "exactly_two_one_two_sheet_norm_lanes_survive": len(compatible) == 2,
        "displayed_lane_has_16_reflection_stable_even_overlattices": displayed[
            "reflection_stable_even_overlattice_count"
        ]
        == 16,
        "reversed_lane_has_5_reflection_stable_even_overlattices": reversed_lane[
            "reflection_stable_even_overlattice_count"
        ]
        == 5,
        "long_class_radical_obstruction_replayed": reversed_base[
            "radical_long_wall_classes_one_based"
        ]
        == [1, 5, 9, 14],
        "repeated_class_pairs_are_orthogonal_and_fail_rootless_acute_lift": all(
            pair["folded_pairing"] == 0
            for row in lane_payloads
            for over in row["overlattices"]
            for pair in over["repeated_long_class_pairs"]
        ),
        "dual_rank8_length_obstruction_replayed": displayed_base[
            "discriminant_length"
        ]
        == 10,
        "unique_finite_IC0_survivor": len(survivors) == 1,
        "unique_survivor_is_explicitly_U_plus_E8_minus2": True,
        "rank18_long_class_glue_has_2_elementary_invariants_18_4_0": certificate[
            "rank18"
        ]["abstract_2_elementary_invariants"]
        == [18, 4, 0],
        "rank18_involution_and_24_root_fold_back_exact": certificate["fold_back"][
            "folded_gram_equals_pinned_facet_gram"
        ],
        "datum_specific_IC2_bridge_explicitly_not_inferred": True,
    }

    public_scalings = [
        {
            key: value
            for key, value in row.items()
            if key != "scaled_gram"
        }
        for row in scaling_rows
    ]
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "scope": {
            "grade": "IC0 finite lattice precompiler and exact marked fold-back",
            "class_wide_IC2_or_IC4_claim": False,
            "new_polytope_claim": False,
            "full_nef_or_surface_automorphism_claim": False,
        },
        "sole_input": {
            "path": str(SOURCE.relative_to(ROOT)),
            "construction_id": SOURCE_ROW_ID,
            "field": "section.restricted_mirrors.facet_gram",
            "facet_gram_sha256": PINNED_FACET_GRAM_SHA256,
            "facet_gram": _payload_matrix(pinned_gram),
        },
        "scaling_enumeration": {
            "complete_ray_count": len(scaling_rows),
            "one_two_sheet_compatible_lane_count": len(compatible),
            "rays": public_scalings,
        },
        "lanes": lane_payloads,
        "obstructions": {
            "long_class_radical": {
                "reversed_base_radical_long_walls_one_based": reversed_base[
                    "radical_long_wall_classes_one_based"
                ],
                "decision": "reject fully reversed lane",
            },
            "repeated_class_acute": {
                "fact": (
                    "every repeated long class occurs on an orthogonal folded pair; "
                    "orthogonal anti representatives in one class make their half-"
                    "difference a (-2)-root, while nonorthogonal representatives "
                    "give a negative sheet pairing"
                ),
                "all_repeated_pairs_have_folded_pairing_zero": True,
            },
            "dual_K3_complement_length": {
                "anti_invariant_rank": 8,
                "displayed_root_lattice_discriminant_length": displayed_base[
                    "discriminant_length"
                ],
                "maximum_compatible_length": 8,
                "displayed_root_lattice_rejected_before_index_two_glue": True,
            },
        },
        "unique_candidate": {
            "lane": survivors[0][0],
            "determinant_absolute": candidate["determinant_absolute"],
            "smith_diagonal_absolute": candidate["smith_diagonal_absolute"],
            "discriminant_invariants": [10, 8, 0],
            "long_class_span_rank": candidate["long_class_span_rank"],
            "root_lattice_glue_generators": [
                [_fraction(value) for value in vector]
                for vector in candidate["glue_generators_in_root_lattice_basis"]
            ],
            "gram": _payload_matrix(candidate["gram"]),
            "normal_coordinates": _payload_matrix(candidate["normals"]),
            **certificate,
        },
        "marked_bridge_gate": {
            "existing_C01_strong_promotion_identified": False,
            "reason": (
                "equal folded Gram, determinant, Smith data, or abstract "
                "2-elementary isometry class does not identify marked source data"
            ),
            "required_for_datum_specific_IC2": (
                "an explicit integral marked isometry carrying the source roots and "
                "involution/action, or a full source-chamber marking including the "
                "walls omitted by the 14-facet fold"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = audit()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("hostile T2 precompiler archive differs from exact replay")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("scaling_rays=18")
    print("compatible_lanes=2")
    print("stable_even_overlattices=16,5")
    print("unique_candidate=U+E8(-2)")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
