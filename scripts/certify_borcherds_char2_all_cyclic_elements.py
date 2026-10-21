#!/usr/bin/env python3
"""Exhaust cyclic element classes of the Borcherds characteristic-two stock.

The 42 square-minus-two walls in the published 210-wall Gram matrix are the
point--line Levi graph of PG(2,4).  We reconstruct PGammaL(3,4).2 on those
walls, lift every generator uniquely through the 168 distinct long-wall
pairing signatures, and ask GAP for every conjugacy class.  Every class is
then lifted to the rank-22 Neron--Severi lattice and subjected to the fixed
rank gate.  High classes receive the exact restricted-facet, Coxeter, and
integral fixed-lattice audits used by the individual 2d and 3b certificates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any, Iterable

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_alexeev_high_dim_reflection_examples import (
    SOURCE_SHA256,
    _extract_matrix,
)
from scripts.certify_borcherds_char2_3b_h11_22 import (
    DEFAULT_SOURCE,
    _digest,
    _matrix_payload,
    _permutation as _class_3b_permutation,
)
from scripts.certify_borcherds_char2_2d_h10_47 import (
    _permutation as _class_2d_permutation,
)
from scripts.certify_borcherds_char2_2a_h13_50_2c_h14_98 import (
    CASES as HIGH_INVOLUTION_CASES,
    _permutation as _high_involution_permutation,
)
from scripts.certify_shimada_disc16_d5_h14_fixed_lattice import (
    _smith_decomposition,
)
from scripts.search_k3_high_rank_involution_folds import _coxeter_audit, _cycles
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia
from scripts.verify_k3_quotient_nef_slices import _orbit_sum_matrix


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "borcherds_char2_all_cyclic_element_completion.json"
)
SCHEMA = "coxeter.borcherds_char2_all_cyclic_element_completion.v1"

# A pinned incidence-preserving identification from our canonical PG(2,4)
# ordering (21 points, then 21 lines) to the labels of the 42 short walls in
# Alexeev's transcription.  It is verified against the Gram on every replay.
CANONICAL_TO_SOURCE_SHORT = (
    0, 29, 10, 14, 4, 2, 28, 12, 6, 38, 17, 27, 50, 33, 40, 45, 8, 16,
    39, 30, 46, 26, 49, 24, 3, 31, 48, 18, 7, 9, 23, 25, 5, 13, 21, 15,
    34, 63, 22, 1, 32, 11,
)


def _gf4_add(left: int, right: int) -> int:
    return left ^ right


def _gf4_mul(left: int, right: int) -> int:
    """Multiply polynomial-basis elements modulo x^2+x+1."""
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left & 4:
            left ^= 7
    return result


GF4_INVERSE = {1: 1, 2: 3, 3: 2}


def _projective_normalize(vector: Iterable[int]) -> tuple[int, int, int]:
    values = tuple(vector)
    for value in values:
        if value:
            scale = GF4_INVERSE[value]
            return tuple(_gf4_mul(scale, entry) for entry in values)  # type: ignore[return-value]
    raise ValueError("the zero vector has no projective class")


def _projective_points() -> tuple[tuple[int, int, int], ...]:
    result = {
        _projective_normalize((first, second, third))
        for first in range(4)
        for second in range(4)
        for third in range(4)
        if first or second or third
    }
    if len(result) != 21:
        raise ArithmeticError("PG(2,4) point count changed")
    return tuple(sorted(result))


POINTS = _projective_points()


def _gf4_sum(values: Iterable[int]) -> int:
    result = 0
    for value in values:
        result = _gf4_add(result, value)
    return result


def _matrix_vector(
    matrix: tuple[tuple[int, int, int], ...], vector: tuple[int, int, int]
) -> tuple[int, int, int]:
    return tuple(
        _gf4_sum(_gf4_mul(matrix[row][column], vector[column]) for column in range(3))
        for row in range(3)
    )  # type: ignore[return-value]


def _matrix_inverse(
    matrix: tuple[tuple[int, int, int], ...]
) -> tuple[tuple[int, int, int], ...]:
    rows = [
        list(matrix[row]) + [int(row == column) for column in range(3)]
        for row in range(3)
    ]
    for column in range(3):
        pivot = next(row for row in range(column, 3) if rows[row][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = GF4_INVERSE[rows[column][column]]
        rows[column] = [_gf4_mul(scale, entry) for entry in rows[column]]
        for row in range(3):
            if row == column or not rows[row][column]:
                continue
            scale = rows[row][column]
            rows[row] = [
                _gf4_add(entry, _gf4_mul(scale, pivot_entry))
                for entry, pivot_entry in zip(rows[row], rows[column])
            ]
    return tuple(
        tuple(rows[row][3 + column] for column in range(3))
        for row in range(3)
    )


def _transpose(
    matrix: tuple[tuple[int, int, int], ...]
) -> tuple[tuple[int, int, int], ...]:
    return tuple(tuple(matrix[column][row] for column in range(3)) for row in range(3))


def _linear_levi_permutation(
    matrix: tuple[tuple[int, int, int], ...]
) -> tuple[int, ...]:
    inverse_transpose = _transpose(_matrix_inverse(matrix))
    point_images = [
        POINTS.index(_projective_normalize(_matrix_vector(matrix, point)))
        for point in POINTS
    ]
    line_images = [
        21 + POINTS.index(
            _projective_normalize(_matrix_vector(inverse_transpose, line))
        )
        for line in POINTS
    ]
    return tuple(point_images + line_images)


def _frobenius_levi_permutation() -> tuple[int, ...]:
    result = []
    for offset in (0, 21):
        result.extend(
            offset
            + POINTS.index(
                _projective_normalize(tuple(_gf4_mul(entry, entry) for entry in vector))
            )
            for vector in POINTS
        )
    return tuple(result)


def _duality_levi_permutation() -> tuple[int, ...]:
    return tuple(index + 21 if index < 21 else index - 21 for index in range(42))


def _canonical_generators() -> dict[str, tuple[int, ...]]:
    transvection = ((1, 1, 0), (0, 1, 0), (0, 0, 1))
    coordinate_cycle = ((0, 1, 0), (0, 0, 1), (1, 0, 0))
    projective_diagonal = ((2, 0, 0), (0, 1, 0), (0, 0, 1))
    return {
        "linear_transvection": _linear_levi_permutation(transvection),
        "linear_coordinate_cycle": _linear_levi_permutation(coordinate_cycle),
        "projective_diagonal": _linear_levi_permutation(projective_diagonal),
        "field_frobenius": _frobenius_levi_permutation(),
        "point_line_duality": _duality_levi_permutation(),
    }


def _short_indices(wall_gram: sp.Matrix) -> tuple[int, ...]:
    result = tuple(index for index in range(210) if wall_gram[index, index] == -2)
    if len(result) != 42:
        raise ArithmeticError("the short-wall count changed")
    return result


def _verify_incidence_identification(wall_gram: sp.Matrix) -> None:
    if set(CANONICAL_TO_SOURCE_SHORT) != set(_short_indices(wall_gram)):
        raise ArithmeticError("the pinned PG(2,4) short-wall labels changed")
    for point_index, point in enumerate(POINTS):
        for line_index, line in enumerate(POINTS):
            incident = _gf4_sum(
                _gf4_mul(point[column], line[column]) for column in range(3)
            ) == 0
            pairing = wall_gram[
                CANONICAL_TO_SOURCE_SHORT[point_index],
                CANONICAL_TO_SOURCE_SHORT[21 + line_index],
            ]
            if pairing != int(incident):
                raise ArithmeticError("the pinned PG(2,4) incidence map changed")


def _lift_levi_permutation(
    wall_gram: sp.Matrix, canonical_permutation: tuple[int, ...]
) -> tuple[int, ...]:
    short = _short_indices(wall_gram)
    short_position = {label: index for index, label in enumerate(short)}
    source_map = {
        CANONICAL_TO_SOURCE_SHORT[index]:
        CANONICAL_TO_SOURCE_SHORT[canonical_permutation[index]]
        for index in range(42)
    }
    short_permutation = tuple(short_position[source_map[label]] for label in short)
    inverse_short = [0] * 42
    for source, target in enumerate(short_permutation):
        inverse_short[target] = source

    long = tuple(index for index in range(210) if wall_gram[index, index] == -4)
    signature_to_long = {
        tuple(int(wall_gram[index, label]) for label in short): index
        for index in long
    }
    if len(signature_to_long) != 168:
        raise ArithmeticError("long walls no longer have unique short-wall signatures")

    result = list(range(210))
    for position, label in enumerate(short):
        result[label] = short[short_permutation[position]]
    for label in long:
        target_signature = tuple(
            int(wall_gram[label, short[inverse_short[position]]])
            for position in range(42)
        )
        result[label] = signature_to_long[target_signature]
    permutation = tuple(result)
    if set(permutation) != set(range(210)):
        raise ArithmeticError("a Levi action did not lift bijectively")
    if any(
        wall_gram[left, right]
        != wall_gram[permutation[left], permutation[right]]
        for left in range(210)
        for right in range(210)
    ):
        raise ArithmeticError("a lifted Levi action fails the full weighted Gram")
    return permutation


def _gap_permutation(permutation: tuple[int, ...]) -> str:
    return "PermList([" + ",".join(str(value + 1) for value in permutation) + "])"


def _gap_classes(
    generators: dict[str, tuple[int, ...]], gap_executable: str = "gap"
) -> tuple[dict[str, int], list[dict[str, Any]]]:
    names = tuple(generators)
    definitions = ",\n".join(_gap_permutation(generators[name]) for name in names)
    program = f"""
SizeScreen([100000,100000]);;
Reset(GlobalRandomSource,1);;
Reset(GlobalMersenneTwister,1);;
gens := [{definitions}];;
G := Group(gens);;
H := Group([gens[1],gens[2],gens[3],gens[5]]);;
Print("GROUP|",Size(G),"|",Size(H),"\\n");
classes := ConjugacyClasses(G);;
Print("COUNT|",Length(classes),"\\n");
for index in [1..Length(classes)] do
  representative := Representative(classes[index]);
  Print("CLASS|",index,"|",Order(representative),"|",Size(classes[index]),"|",
        representative in H,"|",
        JoinStringsWithSeparator(List(ListPerm(representative,210),String),","),"\\n");
od;
"""
    completed = subprocess.run(
        [gap_executable, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    group_data: dict[str, int] = {}
    classes: list[dict[str, Any]] = []
    for line in completed.stdout.splitlines():
        fields = line.split("|")
        if fields[0] == "GROUP":
            group_data = {
                "full_order": int(fields[1]),
                "geometric_subgroup_order": int(fields[2]),
            }
        elif fields[0] == "COUNT":
            group_data["class_count"] = int(fields[1])
        elif fields[0] == "CLASS":
            if len(fields) != 6:
                raise ArithmeticError(
                    "malformed GAP class line: " + repr(line)
                    + "\ncomplete GAP stdout:\n" + completed.stdout
                    + "\ncomplete GAP stderr:\n" + completed.stderr
                )
            permutation = tuple(int(value) - 1 for value in fields[5].split(","))
            classes.append(
                {
                    "_gap_index": int(fields[1]),
                    "order": int(fields[2]),
                    "class_size": int(fields[3]),
                    "geometric": fields[4] == "true",
                    "permutation": permutation,
                }
            )
    if group_data != {
        "full_order": 241920,
        "geometric_subgroup_order": 120960,
        "class_count": len(classes),
    }:
        raise ArithmeticError(f"unexpected GAP group data: {group_data}")
    if sum(row["class_size"] for row in classes) != 241920:
        raise ArithmeticError("GAP conjugacy classes do not partition the full group")
    return group_data, classes


def _permutation_order(permutation: tuple[int, ...]) -> int:
    return math.lcm(*(len(cycle) for cycle in _cycles(permutation)))


def _cycle_histogram(permutation: tuple[int, ...]) -> dict[str, int]:
    histogram: dict[str, int] = {}
    for cycle in _cycles(permutation):
        key = str(len(cycle))
        histogram[key] = histogram.get(key, 0) + 1
    return dict(sorted(histogram.items(), key=lambda item: int(item[0])))


def _ambient_action(
    wall_gram: sp.Matrix, permutation: tuple[int, ...]
) -> tuple[sp.Matrix, sp.Matrix]:
    basis_gram = wall_gram[:22, :22]
    if basis_gram.det() != -4:
        raise ArithmeticError("the first 22 walls no longer form an ambient basis")
    inverse = basis_gram.inv()
    action = sp.Matrix.hstack(
        *(inverse * wall_gram[:22, permutation[index]] for index in range(22))
    )
    order = _permutation_order(permutation)
    if not (
        all(value.is_Integer for value in action)
        and abs(int(action.det())) == 1
        and action.T * basis_gram * action == basis_gram
        and action**order == sp.eye(22)
    ):
        raise ArithmeticError("a wall class did not lift to an integral NS action")
    for wall in range(210):
        coordinates = inverse * wall_gram[:22, wall]
        expected = inverse * wall_gram[:22, permutation[wall]]
        if action * coordinates != expected:
            raise ArithmeticError("an NS class action failed on a published wall")
    return basis_gram, action


def _fixed_lattice(
    basis_gram: sp.Matrix, action: sp.Matrix
) -> tuple[sp.Matrix, sp.Matrix, list[int]]:
    diagonal, _, right = _smith_decomposition(action - sp.eye(22))
    zero_columns = tuple(
        column
        for column in range(22)
        if all(diagonal[row, column] == 0 for row in range(22))
    )
    fixed_basis = right[:, list(zero_columns)]
    if not (
        action * fixed_basis == fixed_basis
        and all(value.is_Integer for value in fixed_basis)
    ):
        raise ArithmeticError("the Smith fixed basis failed")
    fixed_gram = fixed_basis.T * basis_gram * fixed_basis
    smith = smith_normal_form(fixed_gram, domain=ZZ)
    smith_absolute = [
        abs(int(smith[index, index])) for index in range(fixed_gram.rows)
    ]
    if exact_inertia(fixed_gram) != (1, fixed_gram.rows - 1, 0):
        raise ArithmeticError("a finite cyclic fixed lattice is not hyperbolic")
    return fixed_basis, fixed_gram, smith_absolute


def _class_name(row: dict[str, Any]) -> str:
    prime_dictionary = {
        (2, 315, 14): "2a",
        (2, 280, 10): "2b",
        (2, 360, 15): "2c",
        (2, 1008, 11): "2d",
        (3, 2240, 10): "3a",
        (3, 672, 12): "3b",
        (3, 1920, 8): "3c",
    }
    prime_name = prime_dictionary.get(
        (row["order"], row["class_size"], row["fixed_rank"])
    )
    if prime_name is not None:
        return prime_name
    index = row.get("stable_index", row.get("_gap_index", 0))
    return f"o{row['order']}-c{int(index):02d}"


def _pinned_high_representative(class_name: str) -> tuple[int, ...] | None:
    """Use the already frozen representative for exact cross-certificate hashes."""
    if class_name == "2d":
        return _class_2d_permutation()
    if class_name == "3b":
        return _class_3b_permutation()
    if class_name in HIGH_INVOLUTION_CASES:
        return _high_involution_permutation(HIGH_INVOLUTION_CASES[class_name])
    return None


def _facet_lattice_records(
    wall_gram: sp.Matrix,
    basis_gram: sp.Matrix,
    fixed_basis: sp.Matrix,
    fixed_gram: sp.Matrix,
    cycles: list[list[int]],
    facet_indices: list[int],
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
    return records


def _screen_high_class(
    wall_gram: sp.Matrix,
    basis_gram: sp.Matrix,
    action: sp.Matrix,
    permutation: tuple[int, ...],
) -> dict[str, Any]:
    cycles = _cycles(permutation)
    sums = _orbit_sum_matrix(210, cycles)
    orbit_gram = sums.T * wall_gram * sums
    fixed_rank = 22 - (action - sp.eye(22)).rank()
    if not (
        orbit_gram.rank() == fixed_rank
        and exact_inertia(orbit_gram) == (1, fixed_rank - 1, len(cycles) - fixed_rank)
    ):
        raise ArithmeticError("the wall-orbit Gram does not span the fixed space")
    facets = [index for index in range(len(cycles)) if orbit_gram[index, index] < 0]
    isotropic = [index for index in range(len(cycles)) if orbit_gram[index, index] == 0]
    timelike = [index for index in range(len(cycles)) if orbit_gram[index, index] > 0]

    interior = sp.ones(len(cycles), 1)
    pairings = orbit_gram * interior
    interior_square = int((interior.T * orbit_gram * interior)[0])
    if not (interior_square > 0 and min(pairings) > 0):
        raise ArithmeticError("the invariant all-wall sum is not a strict interior point")
    witnesses = []
    for facet in facets:
        witness = interior.copy()
        witness[facet] -= sp.Rational(pairings[facet], orbit_gram[facet, facet])
        witness_pairings = orbit_gram * witness
        witness_square = sp.Rational((witness.T * orbit_gram * witness)[0])
        genuine = (
            witness_square > 0
            and witness_pairings[facet] == 0
            and all(
                value > 0
                for index, value in enumerate(witness_pairings)
                if index != facet
            )
        )
        witnesses.append(
            {
                "orbit_index_one_based": facet + 1,
                "genuine": genuine,
                "witness_square": str(witness_square),
                "minimum_other_pairing": str(
                    min(
                        value
                        for index, value in enumerate(witness_pairings)
                        if index != facet
                    )
                ),
            }
        )

    facet_gram = orbit_gram.extract(facets, facets)
    coxeter = _coxeter_audit(facet_gram, range(len(facets)), [cycles[i] for i in facets])
    fixed_basis, fixed_gram, fixed_smith = _fixed_lattice(basis_gram, action)
    lattice_records = _facet_lattice_records(
        wall_gram, basis_gram, fixed_basis, fixed_gram, cycles, facets
    )
    facet_payload = _matrix_payload(facet_gram)
    all_genuine = all(record["genuine"] for record in witnesses)
    all_reflective = all(record["reflection_integral"] for record in lattice_records)
    return {
        "wall_orbit_count": len(cycles),
        "orbit_length_histogram": {
            str(length): sum(len(cycle) == length for cycle in cycles)
            for length in sorted({len(cycle) for cycle in cycles})
        },
        "negative_facet_count": len(facets),
        "isotropic_inactive_orbit_count": len(isotropic),
        "timelike_inactive_orbit_count": len(timelike),
        "strict_interior_square": interior_square,
        "strict_interior_minimum_pairing": int(min(pairings)),
        "all_negative_restrictions_are_genuine_facets": all_genuine,
        "facet_witnesses": witnesses,
        "facet_orbit_indices_one_based": [index + 1 for index in facets],
        "facet_gram": facet_payload,
        "facet_gram_sha256": _digest(facet_payload),
        "coxeter": coxeter,
        "fixed_lattice": {
            "basis_columns_in_ambient_basis": _matrix_payload(fixed_basis),
            "gram": _matrix_payload(fixed_gram),
            "rank": fixed_gram.rows,
            "signature_positive_negative": [1, fixed_gram.rows - 1],
            "determinant": int(fixed_gram.det()),
            "smith_diagonal_absolute": fixed_smith,
        },
        "primitive_facet_lattice_records": lattice_records,
        "all_facet_reflections_integral": all_reflective,
        "coxeter_reflective_candidate": (
            all_genuine
            and coxeter["all_pairs_are_exact_coxeter_relations"]
            and all_reflective
        ),
    }


def build_certificate(
    source_path: Path = DEFAULT_SOURCE, gap_executable: str = "gap"
) -> dict[str, Any]:
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
        raise ArithmeticError("the ambient Borcherds Gram changed")
    _verify_incidence_identification(wall_gram)

    levi_generators = _canonical_generators()
    generators = {
        name: _lift_levi_permutation(wall_gram, permutation)
        for name, permutation in levi_generators.items()
    }
    group, classes = _gap_classes(generators, gap_executable)
    basis_gram = wall_gram[:22, :22]
    class_rows = []
    for source_row in classes:
        gap_permutation = source_row.pop("permutation")
        source_row.pop("_gap_index")
        if _permutation_order(gap_permutation) != source_row["order"]:
            raise ArithmeticError("GAP and direct permutation orders disagree")
        _, gap_action = _ambient_action(wall_gram, gap_permutation)
        fixed_rank = 22 - (gap_action - sp.eye(22)).rank()
        preliminary = {**source_row, "fixed_rank": fixed_rank}
        class_name = _class_name(preliminary)
        pinned_representative = _pinned_high_representative(class_name)
        permutation = pinned_representative or gap_permutation
        if _permutation_order(permutation) != source_row["order"]:
            raise ArithmeticError("a pinned class representative has the wrong order")
        _, action = _ambient_action(wall_gram, permutation)
        selected_fixed_rank = 22 - (action - sp.eye(22)).rank()
        if selected_fixed_rank != fixed_rank:
            raise ArithmeticError("a pinned representative left its GAP class fingerprint")
        fixed_rank = 22 - (action - sp.eye(22)).rank()
        row = {
            **source_row,
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "trace_on_NS": int(action.trace()),
            "fixed_short_wall_count": sum(
                gap_permutation[index] == index
                for index in _short_indices(wall_gram)
            ),
            "cycle_length_histogram_210": _cycle_histogram(gap_permutation),
        }
        row["class_name"] = class_name
        if pinned_representative is not None:
            row["pinned_representative_permutation_images_one_based"] = [
                value + 1 for value in permutation
            ]
            row["pinned_representative_permutation_sha256"] = _digest(
                [value + 1 for value in permutation]
            )
        if fixed_rank >= 11:
            row["high_rank_screen"] = _screen_high_class(
                wall_gram, basis_gram, action, permutation
            )
        class_rows.append(row)

    def stable_key(row: dict[str, Any]) -> tuple[Any, ...]:
        return (
            row["order"],
            row["class_size"],
            row["geometric"],
            row["fixed_rank"],
            row["trace_on_NS"],
            row["fixed_short_wall_count"],
            tuple(
                (int(length), count)
                for length, count in row["cycle_length_histogram_210"].items()
            ),
        )

    keys = [stable_key(row) for row in class_rows]
    if len(set(keys)) != len(class_rows):
        raise ArithmeticError("the invariant class fingerprints are not unique")
    class_rows.sort(key=stable_key)
    for stable_index, row in enumerate(class_rows, start=1):
        row["stable_index"] = stable_index
        row["class_name"] = _class_name(row)
        fingerprint = {
            "order": row["order"],
            "class_size": row["class_size"],
            "geometric": row["geometric"],
            "fixed_rank": row["fixed_rank"],
            "trace_on_NS": row["trace_on_NS"],
            "fixed_short_wall_count": row["fixed_short_wall_count"],
            "cycle_length_histogram_210": row["cycle_length_histogram_210"],
        }
        row["invariant_class_fingerprint_sha256"] = _digest(fingerprint)

    high = [row for row in class_rows if row["fixed_rank"] >= 11]
    nonidentity_high = [row for row in high if row["order"] > 1]
    current_casts = {
        "2a": {
            "id": "SS2-2A-H13-50",
            "facet_gram_sha256": HIGH_INVOLUTION_CASES["2a"]["facet_gram_sha256"],
        },
        "2c": {
            "id": "SS2-2C-H14-98",
            "facet_gram_sha256": HIGH_INVOLUTION_CASES["2c"]["facet_gram_sha256"],
        },
        "2d": {
            "id": "SS2-2D-H10-47",
            "facet_gram_sha256": "454e34ef2358603578a4360145ac6bdd7abfe6e0b539344323ffa5e354038881",
        },
        "3b": {
            "id": "SS2-3B-H11-22",
            "facet_gram_sha256": "2944c4ce6ebab9a5fb1b52c9a0f7c2a105656863036c3e2badbe52b15c9763cb",
        },
    }
    for row in nonidentity_high:
        cast = current_casts[row["class_name"]]
        screen = row["high_rank_screen"]
        screen["deduplication"] = {
            "current_cast_id": cast["id"],
            "current_cast_facet_gram_sha256": cast["facet_gram_sha256"],
            "exact_facet_gram_match": (
                screen["facet_gram_sha256"] == cast["facet_gram_sha256"]
            ),
            "conclusion": "same pinned chamber, not an additional cyclic cast",
        }
        if not screen["deduplication"]["exact_facet_gram_match"]:
            raise ArithmeticError("a high cyclic class failed current-cast deduplication")
    return {
        "schema": SCHEMA,
        "source": {
            "path": str(source_path),
            "sha256": source_hash,
            "matrix_key": "(22,2,0)",
        },
        "ambient": {
            "lattice": "U + D20(-1)",
            "rank": 22,
            "signature_positive_negative": [1, 21],
            "short_wall_count": 42,
            "long_wall_count": 168,
            "wall_count": 210,
        },
        "reconstruction": {
            "group": "PGammaL(3,4).2 = PSL(3,4).D12",
            **group,
            "geometric_subgroup": "PGL(3,4).2",
            "geometric_subgroup_index": 2,
            "generator_names": list(generators),
            "generator_permutations_sha256": {
                name: _digest([value + 1 for value in permutation])
                for name, permutation in generators.items()
            },
            "all_generators_preserve_full_weighted_gram": True,
            "long_wall_short_signature_count": 168,
        },
        "rank_gate": {
            "threshold_fixed_rank": 11,
            "threshold_hyperbolic_dimension": 10,
            "high_class_count_including_identity": len(high),
            "nonidentity_high_class_count": len(nonidentity_high),
            "high_class_names_including_identity": [
                row["class_name"] for row in high
            ],
            "nonidentity_high_class_names": [
                row["class_name"] for row in nonidentity_high
            ],
            "conclusion": (
                "all nonidentity cyclic sections in H10 or above are exactly "
                "the four current casts 2a, 2c, 2d, and 3b"
            ),
        },
        "classes": class_rows,
        "checks": {
            "pinned_source_hash": True,
            "ambient_H21_210_gram": True,
            "PG2F4_incidence_identification": True,
            "unique_long_wall_lifts": True,
            "full_group_order_241920": group["full_order"] == 241920,
            "geometric_index_two_subgroup_order_120960": (
                group["geometric_subgroup_order"] == 120960
            ),
            "all_element_classes_exhausted": (
                sum(row["class_size"] for row in class_rows) == 241920
            ),
            "every_class_lifts_integrally_to_NS": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", nargs="?", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--gap", default="gap")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.source, args.gap)
    if args.check_archive:
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt all-cyclic certificate differs from archive")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "class_count": len(payload["classes"]),
                "high_classes": payload["rank_gate"][
                    "nonidentity_high_class_names"
                ],
                "candidate_classes": [
                    row["class_name"]
                    for row in payload["classes"]
                    if row["order"] > 1
                    and row.get("high_rank_screen", {}).get(
                        "coxeter_reflective_candidate", False
                    )
                ],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
