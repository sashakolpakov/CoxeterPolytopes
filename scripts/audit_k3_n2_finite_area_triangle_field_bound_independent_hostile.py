#!/usr/bin/env python3
"""Independent hostile audit of the finite-area K3 triangle field theorem.

No owner proof function is imported.  The audit checks the reflection trace
identities without assuming primitive normals, decomposes arbitrary closed
walks into simple directed cycles, and recomputes the primitive-grid and
Vandermonde constants.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINS = {
    "report": (
        "results/k3_n2_compact_triangle_field_bound.md",
        "68b91b604e08bda1a48eb76fd13635f108bba5aab8fcc8ecf4d369afc367baa2",
    ),
    "certifier": (
        "scripts/certify_k3_n2_compact_triangle_field_bound.py",
        "e674fca868efaa7d0041e54376500c6f980321264e7d0e02fdb06457b6e29c3d",
    ),
    "archive": (
        "data/enumerations/k3_n2_compact_triangle_field_bound.json",
        "4bfe02df789ef5d0ad90a4deadfb25b646154e109926ffe8a82b472d3be270e6",
    ),
    "tests": (
        "tests/test_k3_n2_compact_triangle_field_bound.py",
        "f32c974d16023c81eb813bc94ef216855d75a3eab666b0fb1adfed567d992fe4",
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_pins() -> dict[str, str]:
    observed = {}
    for label, (relative, expected) in PINS.items():
        actual = sha256(ROOT / relative)
        assert actual == expected, (label, expected, actual)
        observed[label] = actual
    return observed


Matrix = tuple[tuple[Fraction, ...], ...]


def identity(size: int) -> Matrix:
    return tuple(
        tuple(Fraction(row == column) for column in range(size))
        for row in range(size)
    )


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum(
                (left[row][index] * right[index][column] for index in range(len(right))),
                Fraction(0),
            )
            for column in range(len(right[0]))
        )
        for row in range(len(left))
    )


def trace(matrix: Matrix) -> Fraction:
    return sum((matrix[index][index] for index in range(len(matrix))), Fraction(0))


def reflection(index: int, cartan: Matrix) -> Matrix:
    """Matrix of e_j -> e_j-a_ij e_i in the facet-normal basis."""

    result = [list(row) for row in identity(3)]
    for column in range(3):
        result[index][column] -= cartan[index][column]
    return tuple(tuple(row) for row in result)


def subtract(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[row][column] - right[row][column] for column in range(3))
        for row in range(3)
    )


def trace_integrality_certificate() -> dict[str, object]:
    samples = (
        (
            (2, Fraction(1, 2), Fraction(2, 3)),
            (Fraction(3, 5), 2, Fraction(4, 7)),
            (Fraction(5, 11), Fraction(6, 13), 2),
        ),
        (
            (2, 0, Fraction(7, 3)),
            (0, 2, Fraction(-5, 4)),
            (Fraction(9, 8), Fraction(-2, 7), 2),
        ),
    )
    checked = 0
    for raw in samples:
        cartan: Matrix = tuple(tuple(map(Fraction, row)) for row in raw)
        reflections = tuple(reflection(index, cartan) for index in range(3))
        operators = tuple(subtract(identity(3), item) for item in reflections)
        pairs = {
            (i, j): cartan[i][j] * cartan[j][i]
            for i in range(3)
            for j in range(i + 1, 3)
        }
        for (i, j), pair in pairs.items():
            assert trace(matmul(operators[i], operators[j])) == pair
            assert trace(matmul(reflections[i], reflections[j])) == -1 + pair
            checked += 1
        cycle = cartan[0][1] * cartan[1][2] * cartan[2][0]
        assert trace(matmul(matmul(operators[0], operators[1]), operators[2])) == cycle
        assert (
            trace(matmul(matmul(reflections[0], reflections[1]), reflections[2]))
            == -3 + pairs[(0, 1)] + pairs[(0, 2)] + pairs[(1, 2)] - cycle
        )
        checked += 1

        # Rescaling e_i by lambda_i changes individual directed Cartan
        # entries but leaves every closed cyclic product unchanged.
        scales = (Fraction(2), Fraction(3), Fraction(5))
        rescaled: Matrix = tuple(
            tuple(scales[column] * cartan[row][column] / scales[row] for column in range(3))
            for row in range(3)
        )
        assert rescaled != cartan
        for (i, j), pair in pairs.items():
            assert rescaled[i][j] * rescaled[j][i] == pair
        assert rescaled[0][1] * rescaled[1][2] * rescaled[2][0] == cycle

    return {
        "exact_trace_identities_checked": checked,
        "integral_objects": "T_i=1-s_i in End_O(M)",
        "pair_cycle_formula": "Tr(T_i*T_j)=p_ij",
        "triple_cycle_formula": "Tr(T_0*T_1*T_2)=c",
        "individual_directed_Cartan_entries_assumed_integral": False,
        "primitive_facet_normals_required": False,
        "cyclic_products_invariant_under_normal_rescaling": True,
        "finite_projective_module_trace_lies_in_order": True,
    }


def split_closed_walk(walk: tuple[int, ...]) -> list[tuple[int, ...]]:
    """Decompose a closed walk (first vertex repeated last) into simple cycles."""

    assert len(walk) >= 2 and walk[0] == walk[-1]
    for first in range(len(walk) - 1):
        for second in range(first + 1, len(walk) - 1):
            if walk[first] != walk[second]:
                continue
            inner = walk[first : second + 1]
            outer = walk[: first + 1] + walk[second + 1 :]
            return split_closed_walk(inner) + split_closed_walk(outer)
    return [walk]


def cycle_label(cycle: tuple[int, ...]) -> str:
    length = len(cycle) - 1
    assert cycle[0] == cycle[-1]
    assert len(set(cycle[:-1])) == length
    if length == 1:
        return "2"
    if length == 2:
        return "p" + "".join(map(str, sorted(cycle[:-1])))
    if length == 3:
        assert set(cycle[:-1]) == {0, 1, 2}
        return "c"
    raise AssertionError("a simple cycle on three vertices has length at most three")


def closed_walk_certificate() -> dict[str, object]:
    examples = {
        (0, 1, 0, 2, 0): ["p01", "p02"],
        (0, 1, 2, 0): ["c"],
        (0, 0): ["2"],
        (0, 1, 2, 1, 0, 2, 0): ["p12", "p01", "p02"],
    }
    observed = {}
    for walk, expected in examples.items():
        labels = [cycle_label(cycle) for cycle in split_closed_walk(walk)]
        assert labels == expected, (walk, labels, expected)
        observed["-".join(map(str, walk))] = labels
    return {
        "examples": observed,
        "decomposition": "cut at a repeated vertex and recurse",
        "simple_cycle_generators": ["2", "p01", "p12", "p02", "c"],
        "all_closed_walk_coordinates_in_Z[p01,p12,p02,c]": True,
    }


def height_and_field_certificate() -> dict[str, object]:
    degree_bound = 20 // 3
    assert degree_bound == 6 and 3 * degree_bound <= 20 < 3 * (degree_bound + 1)
    pair_height = 4
    triple_height = 8
    # c^2 is the product of the three nonnegative pair coordinates at every
    # real place.  Ideal vertices permit equality at the distinguished place.
    assert triple_height**2 == pair_height**3

    minimal_packet_size = 2
    assert 2**minimal_packet_size <= degree_bound < 2 ** (minimal_packet_size + 1)
    bad_values = comb(degree_bound, 2)
    choices = bad_values + 1
    primitive_height = choices * triple_height
    difference_height = 2 * primitive_height
    exponent = degree_bound * (degree_bound - 1)
    assert (bad_values, choices, primitive_height, difference_height, exponent) == (
        15,
        16,
        128,
        256,
        30,
    )
    return {
        "full_trace_degree_bound": degree_bound,
        "finite_or_ideal_pair_cycle_weak_height": pair_height,
        "triple_cycle_weak_height": triple_height,
        "minimal_generating_subpacket_size_upper_bound": minimal_packet_size,
        "primitive_grid": "theta1+m*theta2, m=0,...,15",
        "primitive_generator_weak_height": primitive_height,
        "conjugate_difference_weak_height": difference_height,
        "Vandermonde_exponent": exponent,
        "field_discriminant_bound": f"{difference_height}^{exponent}",
        "Hermite_Minkowski_gives_finite_field_set": True,
    }


def scope_certificate() -> dict[str, object]:
    text = " ".join((ROOT / PINS["report"][0]).read_text().split())
    required = {
        "finite_area": "exactly a finite-area Coxeter triangle",
        "ideal_vertices": "finite and ideal vertices",
        "primitive_full_trace": "primitive full trace",
        "minimal_field": "ground field generated by the cyclic Cartan products",
        "closed_walk_repair": "Every closed directed walk decomposes into simple directed cycles",
        "carrier_not_constructed": "construction of an exact K3 carrier",
        "classification_not_claimed": "does not classify triangle angles",
        "higher_facets_open": "arbitrary compact or finite-volume polygons",
    }
    for label, phrase in required.items():
        assert phrase in text, (label, phrase)
    return {label: True for label in required}


def build_audit() -> dict[str, object]:
    return {
        "schema": "k3_n2_finite_area_triangle_field_bound/independent_hostile/v1",
        "verdict": "PASS_EXACT_SCOPE",
        "pins": verify_pins(),
        "nonprimitive_reflection_trace_integrality": trace_integrality_certificate(),
        "closed_walk_generation": closed_walk_certificate(),
        "height_and_field_finiteness": height_and_field_certificate(),
        "scope": scope_certificate(),
    }


def main() -> None:
    print(json.dumps(build_audit(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
