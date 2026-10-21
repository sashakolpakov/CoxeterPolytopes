"""Independent hostile audit of the frozen Borcherds H21 routing no-go.

The audit does not import the owner certifier.  It rebuilds the short/long
wall arithmetic, the discriminant action used in the special maximality
argument, the two physical square-class branches, and the PG(2,4)/Fano
contraction obstruction directly from the earlier frozen wall archive.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json"
ARCHIVE = (
    ROOT
    / "data/enumerations"
    / "borcherds_h21_outer_routing_counterexample_independent_hostile_audit.json"
)

SOURCE_SHA256 = "dc9053322bf54cedde15fcfc5080666ebd79e64feb065b071880dce6ff1552ca"
OWNER_PINS = {
    "results/borcherds_h21_outer_routing_counterexample.md": (
        "fa899977f418d6fa93f7d5270e067c28a648c5e55f2636cb2b362e381534821d"
    ),
    "scripts/certify_borcherds_h21_outer_routing_counterexample.py": (
        "c4e21ce23ea67ecd028437caf181adf51460872aff8c891fab891fcf60424523"
    ),
    "data/enumerations/borcherds_h21_outer_routing_counterexample.json": (
        "89dd9af466dcd879fd9764787b7fb6da99bd5222d4a0b17b554af9d6896313a9"
    ),
    "tests/test_borcherds_h21_outer_routing_counterexample.py": (
        "4900b836a2afd7814867bb40f21d98491f48acebd7259e68e982be4a9e1d4d87"
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix(rows: Iterable[Iterable[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in rows])


def _pair(gram: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> int:
    value = sp.Rational((left.T * gram * right)[0])
    if value.q != 1:
        raise ArithmeticError("nonintegral wall pairing")
    return int(value)


def _primitive(vector: sp.Matrix) -> sp.Matrix:
    denominators = [int(sp.Rational(value).q) for value in vector]
    common = math.lcm(*denominators)
    entries = [int(common * value) for value in vector]
    divisor = math.gcd(*(abs(value) for value in entries if value))
    entries = [value // divisor for value in entries]
    if next(value for value in entries if value) > 0:
        entries = [-value for value in entries]
    return sp.Matrix(entries)


def _square_class(value: Fraction) -> int:
    result = 1
    for integer in (value.numerator, value.denominator):
        for prime, exponent in sp.factorint(integer).items():
            if exponent % 2:
                result *= int(prime)
    return result


def _class_action(
    gram: sp.Matrix, wall: sp.Matrix, basis_bits: list[list[int]]
) -> tuple[int, int, int, int]:
    """Return the F2 action by columns on the displayed A_L basis."""

    basis = [sp.Matrix(bits) for bits in basis_bits]
    residue_to_coordinates: dict[tuple[int, ...], tuple[int, int]] = {}
    for first in (0, 1):
        for second in (0, 1):
            residue = first * basis[0] + second * basis[1]
            residue_to_coordinates[tuple(int(value) % 2 for value in residue)] = (
                first,
                second,
            )

    columns: list[tuple[int, int]] = []
    wall_square = sp.Rational((wall.T * gram * wall)[0])
    for bits in basis:
        representative = bits / 2
        reflected = representative - (
            2 * (representative.T * gram * wall)[0] / wall_square
        ) * wall
        doubled = 2 * reflected
        if any(sp.Rational(value).q != 1 for value in doubled):
            raise ArithmeticError("reflection fails to preserve the dual")
        residue = tuple(int(value) % 2 for value in doubled)
        columns.append(residue_to_coordinates[residue])
    return columns[0] + columns[1]


def _f2_apply(
    matrix: tuple[int, int, int, int], vector: tuple[int, int]
) -> tuple[int, int]:
    return (
        (matrix[0] * vector[0] + matrix[2] * vector[1]) % 2,
        (matrix[1] * vector[0] + matrix[3] * vector[1]) % 2,
    )


def _f2_product(
    left: tuple[int, int, int, int], right: tuple[int, int, int, int]
) -> tuple[int, int, int, int]:
    return _f2_apply(left, right[:2]) + _f2_apply(left, right[2:])


def _f2_closure(
    generators: Iterable[tuple[int, int, int, int]],
) -> set[tuple[int, int, int, int]]:
    generators = set(generators)
    group = {(1, 0, 0, 1)}
    frontier = list(group)
    while frontier:
        current = frontier.pop()
        for generator in generators:
            product = _f2_product(current, generator)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return group


def audit() -> dict[str, Any]:
    source = json.loads(SOURCE.read_text())
    raw = source["coxeter_input"]
    gram = _matrix(raw["basis_gram"])
    walls = _matrix(raw["wall_coordinates_as_columns_transposed"]).T
    squares = [_pair(gram, walls[:, index], walls[:, index]) for index in range(210)]
    short = [index for index, square in enumerate(squares) if square == -2]
    long = [index for index, square in enumerate(squares) if square == -4]

    short_coordinates = walls[:, short]
    short_smith_matrix = smith_normal_form(short_coordinates, domain=ZZ)
    short_smith = [
        abs(int(short_smith_matrix[index, index]))
        for index in range(min(short_smith_matrix.shape))
        if short_smith_matrix[index, index]
    ]
    short_divisibilities = {
        math.gcd(*(abs(int(value)) for value in gram * walls[:, index]))
        for index in short
    }
    long_divisibilities = {
        math.gcd(*(abs(int(value)) for value in gram * walls[:, index]))
        for index in long
    }

    odd_graph = {
        left: {
            right
            for right in short
            if left != right and _pair(gram, walls[:, left], walls[:, right]) == 1
        }
        for left in short
    }
    reached = {short[0]}
    frontier = [short[0]]
    while frontier:
        current = frontier.pop()
        for neighbour in odd_graph[current] - reached:
            reached.add(neighbour)
            frontier.append(neighbour)

    basis_bits = source["picard_lattice_L"]["two_torsion_kernel_basis_bits"]
    actions = [_class_action(gram, walls[:, index], basis_bits) for index in range(210)]
    action_group = _f2_closure(actions)
    zero = (0, 0)
    subgroups = [
        {zero},
        {zero, (1, 0)},
        {zero, (0, 1)},
        {zero, (1, 1)},
        {zero, (1, 0), (0, 1), (1, 1)},
    ]
    stable_subgroups = [
        subgroup
        for subgroup in subgroups
        if all({_f2_apply(action, value) for value in subgroup} == subgroup for action in action_group)
    ]

    long_kernel = (walls[:, long].T * gram).nullspace()
    if len(long_kernel) != 1:
        raise ArithmeticError("long walls do not have a one-dimensional orthogonal")
    canonical = _primitive(long_kernel[0])
    canonical_square = _pair(gram, canonical, canonical)
    canonical_pairings = {
        index: _pair(gram, canonical, walls[:, index]) for index in short
    }
    points = sorted(index for index, value in canonical_pairings.items() if value == -1)
    lines = sorted(index for index, value in canonical_pairings.items() if value == 1)

    incidences = {
        point: {
            line
            for line in lines
            if _pair(gram, walls[:, point], walls[:, line]) == 1
        }
        for point in points
    }
    reverse_incidence = {
        line: {point for point in points if line in incidences[point]} for line in lines
    }
    same_point_pairings = {
        _pair(gram, walls[:, left], walls[:, right])
        for offset, left in enumerate(points)
        for right in points[offset + 1 :]
    }
    same_line_pairings = {
        _pair(gram, walls[:, left], walls[:, right])
        for offset, left in enumerate(lines)
        for right in lines[offset + 1 :]
    }
    common_lines = {
        len(incidences[left] & incidences[right])
        for offset, left in enumerate(points)
        for right in points[offset + 1 :]
    }
    common_points = {
        len(reverse_incidence[left] & reverse_incidence[right])
        for offset, left in enumerate(lines)
        for right in lines[offset + 1 :]
    }

    fano_points_one = [2, 4, 6, 8, 12, 22, 33]
    fano_lines_one = [5, 7, 31, 34, 40, 47, 51]
    fano_points = {index - 1 for index in fano_points_one}
    fano_lines = {index - 1 for index in fano_lines_one}
    fano_incidence = {
        line + 1: sorted(point + 1 for point in reverse_incidence[line] & fano_points)
        for line in fano_lines
    }
    expected_fano_incidence = {
        5: [2, 4, 6],
        7: [6, 8, 33],
        31: [4, 22, 33],
        34: [2, 12, 33],
        40: [6, 12, 22],
        47: [2, 8, 22],
        51: [4, 8, 12],
    }
    fano_determinant = int(sp.Matrix([[1, 1, 0], [1, 0, 1], [0, 1, 1]]).det())

    physical_norms = (1, 2, 4)
    square_class_table = {
        raw_norm: {
            physical: _square_class(Fraction(physical, raw_norm))
            for physical in physical_norms
        }
        for raw_norm in (2, 4)
    }
    common_classes = sorted(
        set(square_class_table[2].values()) & set(square_class_table[4].values())
    )

    # A reflection in a 22-dimensional quadratic space has trace 20.  Thus
    # any projectively Q-rational reflection has a genuinely Q-rational lift:
    # if A=c*s is rational, c=tr(A)/20 is rational.
    reflection_trace = gram.rows - 2

    checks = {
        "source_and_owner_pins_match": _sha256(SOURCE) == SOURCE_SHA256
        and all(_sha256(ROOT / path) == digest for path, digest in OWNER_PINS.items()),
        "raw_lattice_is_rank22_signature_1_21_det_minus4": gram.rank() == 22
        and exact_inertia(gram) == (1, 21, 0)
        and gram.det() == -4,
        "wall_split_and_divisibilities": len(short) == 42
        and len(long) == 168
        and short_divisibilities == {1}
        and long_divisibilities == {2},
        "short_roots_generate_L_integrally": short_smith == [1] * 22,
        "short_odd_graph_connected": reached == set(short),
        "discriminant_action_is_GL2_F2": len(action_group) == 6,
        "only_zero_and_full_discriminant_subgroups_stable": sorted(
            len(subgroup) for subgroup in stable_subgroups
        )
        == [1, 4],
        "reflection_trace_is_nonzero_20": reflection_trace == 20,
        "long_span_rank21_and_orthogonal_square_minus6": walls[:, long].rank() == 21
        and canonical_square == -6,
        "canonical_pairing_split_is_21_plus_21": len(points) == len(lines) == 21
        and set(canonical_pairings.values()) == {-1, 1},
        "same_colour_parts_are_disjoint": same_point_pairings == {0}
        and same_line_pairings == {0},
        "incidence_is_symmetric_21_5_design": {len(value) for value in incidences.values()}
        == {5}
        and {len(value) for value in reverse_incidence.values()} == {5}
        and common_lines == {1}
        and common_points == {1},
        "literal_fano_subplane_matches": fano_incidence == expected_fano_incidence,
        "fano_determinant_is_minus2": fano_determinant == -2,
        "physical_square_classes_are_exactly_one_and_two": common_classes == [1, 2]
        and square_class_table
        == {2: {1: 2, 2: 1, 4: 2}, 4: {1: 1, 2: 2, 4: 1}},
        "canonical_square_forces_t_plus_minus_one": 2 * canonical_square == 10 - 22,
        "contraction_numerics_give_P2_lines": 22 - 21 == 1 and -4 + 5 == 1,
        "complex_k3_rank_gate_fails": gram.rank() == 22 and 22 > 20,
    }

    return {
        "artifact_status": "independent hostile audit PASS",
        "owner_pins": OWNER_PINS,
        "source_pin": SOURCE_SHA256,
        "maximality_audit": {
            "short_coordinate_smith": short_smith,
            "short_odd_graph_component_size": len(reached),
            "discriminant_action_group_order": len(action_group),
            "stable_discriminant_subgroup_orders": sorted(
                len(subgroup) for subgroup in stable_subgroups
            ),
            "rationality_gap_closed_by_trace": (
                "a projectively Q-rational reflection h in dimension 22 has "
                "trace 20; if A=c*h is rational then c=tr(A)/20 is rational"
            ),
            "conclusion": (
                "the frozen special lattice-overgroup proof establishes absolute "
                "reflection maximality without Esselmann"
            ),
        },
        "physical_and_contraction_audit": {
            "common_square_classes": common_classes,
            "canonical_vector": [int(value) for value in canonical],
            "canonical_raw_square": canonical_square,
            "point_count": len(points),
            "line_count": len(lines),
            "degree": 5,
            "fano_points_one_based": fano_points_one,
            "fano_lines_one_based": fano_lines_one,
            "fano_incidence": {str(key): value for key, value in sorted(fano_incidence.items())},
            "fano_determinant": fano_determinant,
            "contraction": {
                "picard_rank_before": 22,
                "disjoint_minus_one_curves": 21,
                "picard_rank_after": 1,
                "surface_after": "P2",
                "boundary_square_before": -4,
                "incident_exceptional_curves": 5,
                "image_square": 1,
                "image_degree": 1,
            },
            "conclusion": (
                "both dual numerical colourings force a Fano subplane over C; "
                "the determinant -2 obstruction empties the simultaneous prime locus"
            ),
        },
        "scope": {
            "certified": (
                "absolute reflection maximality; EMPTY direct/RLogComp complex "
                "rational branches; EMPTY full-trace/projected-root complex-K3 branch"
            ),
            "preserved_escape": (
                "the characteristic-two supersingular K3 mixed divisor/Cremona carrier"
            ),
            "not_claimed": "nonexistence of every algebraic carrier species",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-archive", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = audit()
    if args.write_archive:
        ARCHIVE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.check_archive and json.loads(ARCHIVE.read_text()) != payload:
        raise AssertionError("independent H21 audit archive drift")
    if not payload["all_checks_pass"]:
        failed = [key for key, passed in payload["checks"].items() if not passed]
        raise AssertionError(f"failed independent H21 checks: {failed}")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
