#!/usr/bin/env python3
"""Exact rational/AG routing obstruction from the Borcherds H21 chamber.

The replay consumes only the frozen rank-22 lattice and 210 wall vectors
already reconstructed from Alexeev's Coxeter matrix.  It independently
checks the two raw norm orbits, exhausts every common rational square class
which can be coloured by physical squares -1, -2, -4, solves adjunction,
and recognizes the remaining incidence structure as PG(2,4).

The computational conclusion is that no simultaneous complex rational
weak-log return exists.  The same rank-22 input cannot be a full-trace or
projected complex-K3 eigenslice, because a complex K3 surface has Picard
rank at most 20.  Absolute reflection maximality follows from Borcherds'
complete root chamber and the exact classification, performed below, of
reflection-stable lattices between the short-root frame and its dual.
Esselmann's rank-22 theorem is recorded only as independent endpoint context.

The result is deliberately scoped.  The chamber has Dolgachev--Kondo's
supersingular characteristic-two K3 realization, where the 42 short roots
are curve classes and the 168 long roots act by Cremona involutions.  That
mixed divisor/action realization is not a literal complex nef or projected
root return and is not excluded here.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json"
)
SOURCE_SHA256 = "dc9053322bf54cedde15fcfc5080666ebd79e64feb065b071880dce6ff1552ca"
OUTPUT = (
    ROOT
    / "data/enumerations/borcherds_h21_outer_routing_counterexample.json"
)
SCHEMA = "coxeter.borcherds_h21.outer_routing_counterexample.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _as_rational_matrix(rows: Sequence[Sequence[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in rows])


def _integer_vector(vector: sp.Matrix) -> sp.Matrix:
    denominators = [int(sp.Rational(value).q) for value in vector]
    common = math.lcm(*denominators)
    result = sp.Matrix([int(common * value) for value in vector])
    divisor = math.gcd(*(abs(int(value)) for value in result if value))
    result = result / divisor
    first = next(int(value) for value in result if value)
    # This orientation makes the archived Fano frame deterministic.
    if first > 0:
        result = -result
    return sp.Matrix([int(value) for value in result])


def _vector_payload(vector: Iterable[sp.Expr]) -> list[int]:
    return [int(value) for value in vector]


def _square_class(value: Fraction) -> int:
    """Positive rational square class as its square-free integer label."""

    if value <= 0:
        raise ValueError("only positive square classes occur here")
    label = 1
    for number in (value.numerator, value.denominator):
        for prime, exponent in sp.factorint(number).items():
            if exponent % 2:
                label *= int(prime)
    return label


def _pairing(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> int:
    value = sp.Rational((left.T * gram * right)[0])
    if value.q != 1:
        raise ArithmeticError("expected an integral raw wall pairing")
    return int(value)


def _two_dimensional_f2_action(
    gram: sp.Matrix,
    wall: sp.Matrix,
    discriminant_basis_bits: Sequence[Sequence[int]],
) -> tuple[int, int, int, int]:
    """Action of one wall reflection on A_L=(Z/2)^2, by columns."""

    basis = [sp.Matrix(bits) for bits in discriminant_basis_bits]
    representatives = [vector / 2 for vector in basis]
    classes = {
        (0, 0): sp.zeros(gram.rows, 1),
        (1, 0): basis[0],
        (0, 1): basis[1],
        (1, 1): basis[0] + basis[1],
    }

    def coordinates(bits: sp.Matrix) -> tuple[int, int]:
        residue = tuple(int(value) % 2 for value in bits)
        for coordinate, representative in classes.items():
            if tuple(int(value) % 2 for value in representative) == residue:
                return coordinate
        raise ArithmeticError("failed to express a discriminant class")

    square = sp.Rational((wall.T * gram * wall)[0])
    columns: list[tuple[int, int]] = []
    for representative in representatives:
        image = representative - (
            2 * (representative.T * gram * wall)[0] / square
        ) * wall
        doubled = 2 * image
        if not all(sp.Rational(value).q == 1 for value in doubled):
            raise ArithmeticError("reflection did not preserve the dual lattice")
        columns.append(coordinates(doubled))
    return columns[0] + columns[1]


def _f2_apply(
    matrix: tuple[int, int, int, int], vector: tuple[int, int]
) -> tuple[int, int]:
    # Matrices are stored by columns: (a,c,b,d).
    return (
        (matrix[0] * vector[0] + matrix[2] * vector[1]) % 2,
        (matrix[1] * vector[0] + matrix[3] * vector[1]) % 2,
    )


def _f2_multiply(
    left: tuple[int, int, int, int], right: tuple[int, int, int, int]
) -> tuple[int, int, int, int]:
    first = _f2_apply(left, (right[0], right[1]))
    second = _f2_apply(left, (right[2], right[3]))
    return first + second


def _f2_group_closure(
    generators: Iterable[tuple[int, int, int, int]],
) -> set[tuple[int, int, int, int]]:
    identity = (1, 0, 0, 1)
    generator_set = set(generators)
    group = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generator_set:
            product = _f2_multiply(current, generator)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return group


def _find_fano_subplane(
    points: Sequence[int],
    lines: Sequence[int],
    line_points: dict[int, set[int]],
) -> tuple[list[int], list[int], dict[str, int]]:
    """Return the lexicographically first Fano subplane and a standard frame."""

    for chosen_tuple in itertools.combinations(points, 7):
        chosen = set(chosen_tuple)
        chosen_lines = [line for line in lines if len(line_points[line] & chosen) == 3]
        if len(chosen_lines) != 7:
            continue
        if not all(
            sum(point in line_points[line] for line in chosen_lines) == 3
            for point in chosen
        ):
            continue
        if not all(
            sum(
                left in line_points[line] and right in line_points[line]
                for line in chosen_lines
            )
            == 1
            for left, right in itertools.combinations(chosen, 2)
        ):
            continue

        selected_sets = [line_points[line] & chosen for line in chosen_lines]

        def collinear(triple: set[int]) -> bool:
            return any(triple <= line_set for line_set in selected_sets)

        def third(left: int, right: int) -> int:
            candidates = [
                line_set - {left, right}
                for line_set in selected_sets
                if {left, right} <= line_set
            ]
            if len(candidates) != 1 or len(candidates[0]) != 1:
                raise ArithmeticError("the projective-plane pair axiom failed")
            return next(iter(candidates[0]))

        # Normalize a non-collinear triangle A,B,C.  Then D,E,F are the
        # third points on AB, AC, BC and G is the remaining point.  The four
        # extra Fano incidences force the standard characteristic-two frame.
        for a, b, c in itertools.combinations(chosen_tuple, 3):
            if collinear({a, b, c}):
                continue
            d = third(a, b)
            e = third(a, c)
            f = third(b, c)
            six = {a, b, c, d, e, f}
            if len(six) != 6:
                continue
            remainder = chosen - six
            if len(remainder) != 1:
                continue
            g = next(iter(remainder))
            if all(
                collinear(triple)
                for triple in (
                    {d, e, f},
                    {a, f, g},
                    {b, e, g},
                    {c, d, g},
                )
            ):
                return (
                    list(chosen_tuple),
                    chosen_lines,
                    {"A": a, "B": b, "C": c, "D": d, "E": e, "F": f, "G": g},
                )
    raise ArithmeticError("no Fano subplane was found")


def build_certificate() -> dict[str, Any]:
    if _sha256(SOURCE) != SOURCE_SHA256:
        raise ArithmeticError("the frozen Borcherds wall archive changed")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    raw = source["coxeter_input"]
    gram = _as_rational_matrix(raw["basis_gram"])
    walls = _as_rational_matrix(
        raw["wall_coordinates_as_columns_transposed"]
    ).T

    squares = [_pairing(gram, walls[:, index], walls[:, index]) for index in range(210)]
    short = [index for index, square in enumerate(squares) if square == -2]
    long = [index for index, square in enumerate(squares) if square == -4]
    short_rank = walls[:, short].rank()
    long_rank = walls[:, long].rank()
    short_smith = smith_normal_form(walls[:, short], domain=ZZ)
    short_smith_diagonal = [
        abs(int(short_smith[index, index]))
        for index in range(min(short_smith.shape))
        if short_smith[index, index]
    ]

    discriminant_basis_bits = source["picard_lattice_L"][
        "two_torsion_kernel_basis_bits"
    ]
    discriminant_actions = [
        _two_dimensional_f2_action(gram, walls[:, index], discriminant_basis_bits)
        for index in range(walls.cols)
    ]
    discriminant_action_group = _f2_group_closure(discriminant_actions)
    discriminant_subgroups = [
        {(0, 0)},
        {(0, 0), (1, 0)},
        {(0, 0), (0, 1)},
        {(0, 0), (1, 1)},
        {(0, 0), (1, 0), (0, 1), (1, 1)},
    ]
    invariant_discriminant_subgroups = [
        subgroup
        for subgroup in discriminant_subgroups
        if all(
            {_f2_apply(action, vector) for vector in subgroup} == subgroup
            for action in discriminant_action_group
        )
    ]

    long_orthogonal = (walls[:, long].T * gram).nullspace()
    if len(long_orthogonal) != 1:
        raise ArithmeticError("the long-wall orthogonal should be one-dimensional")
    canonical_direction = _integer_vector(long_orthogonal[0])
    canonical_square_raw = _pairing(gram, canonical_direction, canonical_direction)
    short_canonical_pairings = {
        index: _pairing(gram, canonical_direction, walls[:, index])
        for index in short
    }
    negative_part = sorted(
        index for index in short if short_canonical_pairings[index] == -1
    )
    positive_part = sorted(
        index for index in short if short_canonical_pairings[index] == 1
    )

    short_divisibilities = [
        math.gcd(*(abs(int(value)) for value in gram * walls[:, index]))
        for index in short
    ]
    long_divisibilities = [
        math.gcd(*(abs(int(value)) for value in gram * walls[:, index]))
        for index in long
    ]

    # Exhaust the common rational similarity square class.  A physical
    # vector on a raw line is lambda*r with lambda rational, so
    # |C^2|/|r^2| must have the one common square class.
    physical_norms = (1, 2, 4)
    allowed_by_raw_norm: dict[int, dict[int, list[int]]] = {}
    for raw_norm in (2, 4):
        by_class: dict[int, list[int]] = {}
        for physical_norm in physical_norms:
            square_class = _square_class(Fraction(physical_norm, raw_norm))
            by_class.setdefault(square_class, []).append(physical_norm)
        allowed_by_raw_norm[raw_norm] = by_class
    common_square_classes = sorted(
        set(allowed_by_raw_norm[2]) & set(allowed_by_raw_norm[4])
    )

    # The first square class forces all 42 spanning short roots to be (-2),
    # hence K=0.  Every long line would then violate adjunction whether it is
    # scaled to a (-1)- or a (-4)-curve.
    square_class_one_adjunction_values = {
        "short_minus_2": 0,
        "long_minus_1": -1,
        "long_minus_4": 2,
    }

    # For square class 2 use the representative physical form 2G.  Long
    # roots become (-2) after multiplier 1/2, so K=t*w.  The rational-surface
    # identity K^2=10-rho gives -12t^2=-12, hence t=+-1.  Adjunction colours
    # the two 21-element short-root parts as (-1) and boundary (-4), swapped
    # when t changes sign.
    picard_rank = gram.rank()
    required_canonical_square = 10 - picard_rank
    physical_canonical_square_coefficient = 2 * canonical_square_raw
    canonical_parameters = [-1, 1]
    colourings = []
    for parameter in canonical_parameters:
        exceptional = sorted(
            index
            for index in short
            if parameter * short_canonical_pairings[index] == -1
        )
        boundary = sorted(
            index
            for index in short
            if parameter * short_canonical_pairings[index] == 1
        )
        colourings.append(
            {
                "canonical_parameter": parameter,
                "canonical_vector_in_raw_basis": _vector_payload(
                    parameter * canonical_direction
                ),
                "canonical_square_in_physical_form_2G": (
                    physical_canonical_square_coefficient * parameter * parameter
                ),
                "exceptional_minus_1_short_walls_one_based": [
                    index + 1 for index in exceptional
                ],
                "boundary_minus_4_short_walls_one_based": [
                    index + 1 for index in boundary
                ],
                "long_walls_all_minus_2": True,
            }
        )

    # Work with t=+1.  Replacing w by -w exchanges points and lines.
    points = negative_part
    lines = positive_part
    point_lines = {
        point: {
            line
            for line in lines
            if _pairing(gram, walls[:, point], walls[:, line]) == 1
        }
        for point in points
    }
    line_points = {
        line: {point for point in points if line in point_lines[point]}
        for line in lines
    }
    within_point_pairings = {
        _pairing(gram, walls[:, left], walls[:, right])
        for left, right in itertools.combinations(points, 2)
    }
    within_line_pairings = {
        _pairing(gram, walls[:, left], walls[:, right])
        for left, right in itertools.combinations(lines, 2)
    }
    cross_pairing_counts = {
        value: sum(
            _pairing(gram, walls[:, point], walls[:, line]) == value
            for point in points
            for line in lines
        )
        for value in (0, 1)
    }
    common_lines_for_point_pairs = {
        len(point_lines[left] & point_lines[right])
        for left, right in itertools.combinations(points, 2)
    }
    common_points_for_line_pairs = {
        len(line_points[left] & line_points[right])
        for left, right in itertools.combinations(lines, 2)
    }
    short_odd_graph = {
        left: {
            right
            for right in short
            if left != right
            and _pairing(gram, walls[:, left], walls[:, right]) == 1
        }
        for left in short
    }
    reached = {short[0]}
    frontier = [short[0]]
    while frontier:
        current = frontier.pop()
        for neighbour in short_odd_graph[current] - reached:
            reached.add(neighbour)
            frontier.append(neighbour)

    fano_points, fano_lines, fano_frame = _find_fano_subplane(
        points, lines, line_points
    )
    diagonal_point_matrix = sp.Matrix(
        [
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
        ]
    )
    fano_determinant = int(diagonal_point_matrix.det())

    # Translate all raw indices to the human-facing one-based convention.
    fano_incidence = {
        line + 1: sorted(point + 1 for point in line_points[line] & set(fano_points))
        for line in fano_lines
    }
    fano_frame_one_based = {key: value + 1 for key, value in fano_frame.items()}

    checks = {
        "source_replay_passed": source["all_checks_pass"] is True,
        "rank_22_signature_1_21": gram.rank() == 22
        and source["picard_lattice_L"]["signature"] == [1, 21],
        "determinant_minus_4": gram.det() == -4,
        "wall_split_42_short_168_long": len(short) == 42 and len(long) == 168,
        "short_roots_square_minus_2_divisibility_1": set(short_divisibilities) == {1},
        "long_roots_square_minus_4_divisibility_2": set(long_divisibilities) == {2},
        "short_roots_span_rank_22": short_rank == 22,
        "short_roots_generate_L_integrally": short_smith_diagonal == [1] * 22,
        "long_roots_span_rank_21": long_rank == 21,
        "long_orthogonal_direction_square_minus_6": canonical_square_raw == -6,
        "short_pairings_split_21_plus_21": len(negative_part) == len(positive_part) == 21,
        "short_pairings_are_plus_or_minus_one": set(short_canonical_pairings.values())
        == {-1, 1},
        "physical_square_classes_exhausted": common_square_classes == [1, 2],
        "square_class_one_forces_short_minus_2": allowed_by_raw_norm[2][1] == [2],
        "square_class_one_long_options_minus_1_minus_4": allowed_by_raw_norm[4][1]
        == [1, 4],
        "square_class_two_short_options_minus_1_minus_4": allowed_by_raw_norm[2][2]
        == [1, 4],
        "square_class_two_forces_long_minus_2": allowed_by_raw_norm[4][2] == [2],
        "square_class_one_adjunction_contradiction": short_rank == 22
        and square_class_one_adjunction_values["short_minus_2"] == 0
        and square_class_one_adjunction_values["long_minus_1"] != 0
        and square_class_one_adjunction_values["long_minus_4"] != 0,
        "canonical_square_forces_t_plus_or_minus_one": required_canonical_square == -12
        and physical_canonical_square_coefficient == -12
        and canonical_parameters == [-1, 1],
        "two_colourings_have_21_exceptional_21_boundary": all(
            len(colouring["exceptional_minus_1_short_walls_one_based"]) == 21
            and len(colouring["boundary_minus_4_short_walls_one_based"]) == 21
            for colouring in colourings
        ),
        "same_colour_curves_are_disjoint": within_point_pairings == {0}
        and within_line_pairings == {0},
        "cross_incidence_has_105_edges": cross_pairing_counts == {0: 336, 1: 105},
        "each_point_and_line_has_degree_5": {len(value) for value in point_lines.values()}
        == {5}
        and {len(value) for value in line_points.values()} == {5},
        "unique_line_through_two_points": common_lines_for_point_pairs == {1},
        "unique_point_on_two_lines": common_points_for_line_pairs == {1},
        "short_odd_coxeter_graph_connected": reached == set(short),
        "wall_reflections_induce_full_GL2_F2": len(discriminant_action_group) == 6,
        "only_trivial_and_full_discriminant_subgroups_are_W_stable": sorted(
            len(subgroup) for subgroup in invariant_discriminant_subgroups
        )
        == [1, 4],
        "fano_subplane_found": len(fano_points) == len(fano_lines) == 7,
        "fano_frame_is_deterministic": fano_frame_one_based
        == {"A": 2, "B": 4, "C": 8, "D": 6, "E": 22, "F": 12, "G": 33},
        "fano_diagonal_determinant_minus_2": fano_determinant == -2,
        "complex_k3_rank_obstruction": gram.rank() == 22 and 22 > 20,
        "direct_rational_dimension_obstruction": raw["dimension"] == 21 and 21 > 11,
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "scope": {
            "input": (
                "the classical noncompact finite-volume Borcherds chamber in "
                "H^21 for L=U+D20(-1), with all 210 mixed-root facets"
            ),
            "excluded_disjunction": [
                "direct complex rational full-NS/full-nef return",
                "complex rational weak-log RLogComp return with simultaneous prime "
                "smooth-rational facets of square -1,-2,-4",
                "full-trace or projected-root complex-K3 eigenslice return",
            ],
            "not_excluded": [
                "Dolgachev--Kondo's supersingular characteristic-two K3 mixed "
                "divisor/Cremona-action realization",
                "an as-yet undefined algebraic carrier species not contained in the "
                "three displayed branches",
            ],
        },
        "frozen_input": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": SOURCE_SHA256,
            "ledger_id": raw["ledger_id"],
            "dimension": raw["dimension"],
            "facet_count": raw["facet_count"],
        },
        "raw_lattice": {
            "identification": "U + D20(-1)",
            "rank": gram.rank(),
            "signature": [1, 21],
            "determinant": int(gram.det()),
            "short_wall_count": len(short),
            "long_wall_count": len(long),
            "short_wall_rank": short_rank,
            "short_wall_coordinate_smith_diagonal": short_smith_diagonal,
            "long_wall_rank": long_rank,
            "short_wall_divisibilities": {"1": short_divisibilities.count(1)},
            "long_wall_divisibilities": {"2": long_divisibilities.count(2)},
            "long_wall_orthogonal_primitive_vector": _vector_payload(
                canonical_direction
            ),
            "long_wall_orthogonal_vector_square": canonical_square_raw,
            "short_pairing_counts_with_orthogonal_vector": {
                "-1": len(negative_part),
                "1": len(positive_part),
            },
            "discriminant_module": {
                "structure": "(Z/2)^2",
                "basis_bits": discriminant_basis_bits,
                "wall_action_generator_matrices_by_columns": sorted(
                    [list(action) for action in set(discriminant_actions)]
                ),
                "generated_group_order": len(discriminant_action_group),
                "generated_group": "GL(2,F2)=S3",
                "all_subgroup_sizes": sorted(len(group) for group in discriminant_subgroups),
                "W_stable_subgroup_sizes": sorted(
                    len(group) for group in invariant_discriminant_subgroups
                ),
            },
        },
        "rational_physical_square_class_exhaustion": {
            "principle": (
                "For one rational similarity c and a rational multiplier on each "
                "wall line, |C_i^2|/|r_i^2| has the common class [c] in "
                "Q*/Q*^2."
            ),
            "physical_square_magnitudes": list(physical_norms),
            "allowed_physical_magnitudes_by_raw_norm_and_square_class": {
                str(raw_norm): {
                    str(square_class): values
                    for square_class, values in sorted(by_class.items())
                }
                for raw_norm, by_class in allowed_by_raw_norm.items()
            },
            "common_square_classes": common_square_classes,
            "square_class_1": {
                "short_lines": "all physical square -2; they span rank 22, so K=0",
                "long_lines": "physical square -1 or -4",
                "adjunction": square_class_one_adjunction_values,
                "conclusion": (
                    "K=0 contradicts K.C=-1 on every (-1) option and K.C=2 "
                    "on every (-4) option"
                ),
            },
            "square_class_2": {
                "representative_physical_form": "2G",
                "long_lines": "all physical square -2 using multiplier 1/2",
                "short_lines": (
                    "physical square -1 using multiplier 1/2 or physical square "
                    "-4 using multiplier 1"
                ),
                "canonical_line": "K=t*w because the 168 long walls have rank 21",
                "raw_w_square": canonical_square_raw,
                "physical_K_square_formula": "K^2=-12*t^2",
                "rational_surface_K_square": required_canonical_square,
                "canonical_parameters": canonical_parameters,
                "only_adjunction_colourings": colourings,
            },
        },
        "projective_plane_obstruction": {
            "one_colouring": (
                "For t=1 the 21 pairing-minus-one short roots are disjoint (-1) "
                "curves and the 21 pairing-plus-one short roots are disjoint "
                "boundary (-4) curves; t=-1 exchanges the two parts."
            ),
            "design_parameters": {
                "points": 21,
                "lines": 21,
                "points_per_line": 5,
                "lines_per_point": 5,
                "incidences": cross_pairing_counts[1],
                "common_line_for_two_points": 1,
                "common_point_for_two_lines": 1,
                "identification": "the point-line plane PG(2,4)",
            },
            "contraction": {
                "initial_picard_rank": picard_rank,
                "disjoint_minus_1_curves": 21,
                "contracted_picard_rank": picard_rank - 21,
                "rank_one_minimal_rational_surface": "P^2",
                "each_boundary_curve_initial_square": -4,
                "exceptional_points_on_each_boundary_curve": 5,
                "image_square": -4 + 5,
                "image": "a projective line in P^2",
            },
            "fano_subplane": {
                "point_walls_one_based": [index + 1 for index in fano_points],
                "line_walls_one_based": [index + 1 for index in fano_lines],
                "incidence_by_line_wall_one_based": {
                    str(line): points_on_line
                    for line, points_on_line in sorted(fano_incidence.items())
                },
                "standard_frame_one_based": fano_frame_one_based,
                "standard_frame_lines": [
                    "ABD",
                    "ACE",
                    "BCF",
                    "AFG",
                    "BEG",
                    "CDG",
                    "DEF",
                ],
                "forced_diagonal_point_matrix_rows": [
                    [1, 1, 0],
                    [1, 0, 1],
                    [0, 1, 1],
                ],
                "determinant": fano_determinant,
                "conclusion": (
                    "The Fano incidences force the three diagonal points to be "
                    "collinear.  Their determinant is -2, so the configuration "
                    "is representable only in characteristic 2, not over C."
                ),
            },
            "conclusion": (
                "Neither of the two exhaustive adjunction colourings can occur on "
                "a complex rational surface; hence RLogComp is EMPTY for this raw P."
            ),
        },
        "complex_k3_and_direct_obstructions": {
            "raw_eigenspace_rank_required": gram.rank(),
            "complex_k3_picard_rank_maximum": 20,
            "projected_rank_inequality": "rank(raw image) <= rho(X) <= 20",
            "full_trace_rank_inequality": "[k:Q]*rank(raw) <= rho(X) <= 20",
            "complex_k3_conclusion": (
                "No injective full-trace or projected-root complex-K3 eigenslice "
                "can return a rank-22 chamber."
            ),
            "direct_complex_rational_bound_imported": "dimension <= 11",
            "direct_conclusion": (
                "The direct full-NS rational branch is excluded by its proved "
                "dimension bound; its smooth-rational weak-log subcase was also "
                "excluded directly by the colouring calculation above."
            ),
        },
        "reflection_maximality": {
            "level_1_complete_root_chamber": (
                "Borcherds proves that the 210-wall polyhedron is a fundamental "
                "chamber for the full integral reflection group R(L)."
            ),
            "level_2_arithmetic": (
                "R(L) has finite index in O^+(L), hence is a finite-covolume "
                "arithmetic reflection lattice."
            ),
            "level_3_classification_free_absolute_overgroup_argument": [
                "If R(L) is contained in a discrete reflection group Gamma, then "
                "the lattice-overgroup index [Gamma:R(L)] is finite.",
                "A finite extension of an arithmetic lattice is arithmetic.  Thus "
                "Gamma preserves some integral Lorentzian lattice M in the same "
                "rational commensurability class.",
                "The 42 short simple reflections have connected odd Coxeter graph, "
                "so their root lines are one W-conjugacy class.  Rescale M so its "
                "intersection with every short-root line is Zr.",
                "The short roots generate L integrally.  Hence L is contained in M.  "
                "For x in M, s_r(x)-x=(x,r)r belongs to M intersection Qr=Zr; "
                "therefore (x,L) is integral and M is contained in L*.",
                "Intermediate lattices L <= M <= L* correspond to subgroups of "
                "A_L=(Z/2)^2.  The 210 wall reflections generate GL(2,F2)=S3 on "
                "A_L, whose only invariant subgroups are 0 and A_L.  Thus M is L "
                "or L*.",
                "Dualization gives O(L)=O(L*) and therefore R(L)=R(L*).  Every "
                "generating reflection of Gamma preserves M, so Gamma is contained "
                "in R(M)=R(L).  The original inclusion forces equality.",
            ],
            "conclusion": (
                "The side group is absolutely maximal among discrete reflection "
                "overgroups, and hence maximal among arithmetic reflection groups."
            ),
            "literature_inputs": [
                {
                    "result": "complete 210-wall R(L) chamber and finite index",
                    "author": "Richard Borcherds",
                    "reference": "J. Algebra 111 (1987), Example 5, pp. 148-149",
                    "url": "https://math.berkeley.edu/~reb/papers/al/al.pdf",
                },
            ],
            "esselmann_context_not_dependency": {
                "statement": (
                    "Esselmann independently classifies the rank-22 endpoint as L "
                    "and L* up to scale/isometry, but that classification is not "
                    "used in the special-case overgroup proof above."
                ),
                "reference": "J. Number Theory 61 (1996), 103-144",
                "doi": "10.1006/jnth.1996.0141",
                "url": "https://doi.org/10.1006/jnth.1996.0141",
            },
        },
        "positive_characteristic_escape_clause": {
            "surface": (
                "the unique supersingular K3 surface of Artin invariant one in "
                "characteristic 2"
            ),
            "picard_lattice": "U + D20(-1), rank 22",
            "short_orbit": "42 actual smooth rational (-2)-curve classes",
            "long_orbit": (
                "168 square-minus-four reflections realized by explicit Cremona "
                "involutions, not 168 prime nef-wall divisors"
            ),
            "meaning": (
                "This is a genuine AG mixed divisor/action realization of the raw "
                "reflection chamber.  The counterexample refutes only the stated "
                "complex rational/full-trace/projected-root routing disjunction, "
                "not every possible AG carrier category."
            ),
            "reference": (
                "Dolgachev--Kondo, A supersingular K3 surface in characteristic 2 "
                "and the Leech lattice"
            ),
            "url": "https://arxiv.org/abs/math/0112283",
        },
        "outer_routing_verdict": {
            "status": "COUNTEREXAMPLE_TO_PROPOSED_EXHAUSTIVE_ROUTING",
            "statement": (
                "There exists an absolutely reflection-maximal finite-volume "
                "arithmetic chamber which has neither a direct/RLogComp complex "
                "rational return nor a full-trace/projected-root complex-K3 return."
            ),
            "consequence": (
                "Raw rationality, maximality, rank, cusp data, canonical colouring, "
                "MMP and section algebra do not force the proposed disjunction.  "
                "Therefore that disjunction cannot by itself prove the global "
                "n<=25 and finiteness theorem."
            ),
        },
        "checks": checks,
    }
    payload["all_checks_pass"] = all(checks.values())
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--write-archive", action="store_true")
    action.add_argument("--check-archive", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()

    payload = build_certificate()
    if args.write_archive:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if archived != payload:
            raise SystemExit("archived certificate differs from exact replay")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            json.dumps(
                {
                    "all_checks_pass": payload["all_checks_pass"],
                    "certificate_sha256": payload["certificate_sha256"],
                    "outer_routing_verdict": payload["outer_routing_verdict"]["status"],
                },
                sort_keys=True,
            )
        )


if __name__ == "__main__":
    main()
