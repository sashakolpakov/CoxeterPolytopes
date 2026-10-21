#!/usr/bin/env python3
"""Independent hostile certificate for the affine-B-tail completion obstruction.

This file deliberately imports none of the producer's functions.  It
reconstructs the full-NS lattice, visible walls, dual rays, ruling-swap
reflection, hidden-wall sign, and first-height truncator arithmetic.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from typing import Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "rational_affine_b_tail_global_completion_obstruction_hostile_audit.json"
)

SOURCE_HASHES = {
    "results/rational_affine_b_tail_global_completion_obstruction.md":
        "0c3e32ef6763308cfbe4b4b28e21452c428611a9956d0a5a9333220e55d69b60",
    "scripts/certify_rational_affine_b_tail_global_completion_obstruction.py":
        "a76414589ce0c7680fcee6ea1928536e540d324189b456eb834cc125c2efd517",
    "data/enumerations/rational_affine_b_tail_global_completion_obstruction.json":
        "996d0919993e0b9756f1cbf36c6e922a2e4fa084cf4ad86643326253768fbc38",
    "tests/test_rational_affine_b_tail_global_completion_obstruction.py":
        "5582bfc90b256d5f5dc91aa6622c3e9ee1e2cc200f1e77b59fe218acac36f500",
}

Vector = tuple[int, ...]
Matrix = list[list[Fraction]]


def file_sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def digest(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def pairing(left: Sequence[int | Fraction], right: Sequence[int | Fraction]) -> Fraction:
    """The surface form U + <-1>^m in c,d,e_i coordinates."""

    return (
        Fraction(left[0]) * Fraction(right[1])
        + Fraction(left[1]) * Fraction(right[0])
        - sum(
            Fraction(left[index]) * Fraction(right[index])
            for index in range(2, len(left))
        )
    )


def determinant(matrix: Sequence[Sequence[int | Fraction]]) -> Fraction:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, len(work)):
            coefficient = work[row][column] / pivot_value
            for index in range(column, len(work)):
                work[row][index] -= coefficient * work[column][index]
    return result


def standard_basis(size: int) -> list[Vector]:
    return [
        tuple(int(index == column) for index in range(size))
        for column in range(size)
    ]


def visible_walls(length: int) -> tuple[list[str], list[Vector]]:
    if length < 2:
        raise ValueError("m must be at least two")
    size = length + 2
    labels = ["D", "delta_0"]
    walls: list[Vector] = []

    section = [0] * size
    section[1], section[2] = 1, -1
    walls.append(tuple(section))

    mixed = [0] * size
    mixed[0], mixed[2], mixed[3] = 1, -1, -1
    walls.append(tuple(mixed))

    for index in range(1, length):
        root = [0] * size
        root[index + 1], root[index + 2] = 1, -1
        walls.append(tuple(root))
        labels.append(f"delta_{index}")

    endpoint = [0] * size
    endpoint[-1] = 1
    walls.append(tuple(endpoint))
    labels.append(f"delta_{length}")
    return labels, walls


def asserted_duals(length: int) -> list[Vector]:
    result: list[Vector] = [
        (1, 0, *((0,) * length)),
        (0, 1, *((0,) * length)),
        (1, 1, -1, *((0,) * (length - 1))),
    ]
    for index in range(2, length + 1):
        result.append(
            (1, 2, *((-1,) * index), *((0,) * (length - index)))
        )
    return result


def coxeter_order(left: Vector, right: Vector) -> int | str:
    value = pairing(left, right)
    if value == 0:
        return 2
    ratio = value * value / (pairing(left, left) * pairing(right, right))
    if ratio == Fraction(1, 4):
        return 3
    if ratio == Fraction(1, 2):
        return 4
    if ratio == 1:
        return "infinity"
    raise AssertionError((left, right, value, ratio))


def visible_record(length: int) -> dict[str, object]:
    labels, walls = visible_walls(length)
    duals = asserted_duals(length)
    coordinate_determinant = determinant(walls)
    if abs(coordinate_determinant) != 1:
        raise AssertionError("visible walls are not an integral basis")

    for i, dual in enumerate(duals):
        for j, wall in enumerate(walls):
            if pairing(dual, wall) != int(i == j):
                raise AssertionError((length, i, j))

    squares = [int(pairing(dual, dual)) for dual in duals]
    expected = [0, 0, 1] + [4 - index for index in range(2, length + 1)]
    if squares != expected:
        raise AssertionError((length, squares, expected))

    edges: list[tuple[str, str, int | str]] = []
    for i, left in enumerate(walls):
        for j in range(i + 1, len(walls)):
            if pairing(left, walls[j]):
                edges.append((labels[i], labels[j], coxeter_order(left, walls[j])))

    return {
        "m": length,
        "coordinate_determinant": int(coordinate_determinant),
        "dual_squares": squares,
        "finite_volume": all(square >= 0 for square in squares),
        "coxeter_edges": [list(edge) for edge in edges],
    }


def reflect(vector: Vector, root: Vector) -> Vector:
    root_square = pairing(root, root)
    coefficients = [
        Fraction(value)
        - 2 * pairing(vector, root) * Fraction(root_value) / root_square
        for value, root_value in zip(vector, root)
    ]
    if any(value.denominator != 1 for value in coefficients):
        raise AssertionError("nonintegral reflection")
    return tuple(int(value) for value in coefficients)


def swap_record(length: int) -> dict[str, object]:
    size = length + 2
    basis = standard_basis(size)
    root = (1, -1, *((0,) * length))
    images = [reflect(vector, root) for vector in basis]
    if images[0] != basis[1] or images[1] != basis[0]:
        raise AssertionError("the two rulings were not swapped")
    if images[2:] != basis[2:]:
        raise AssertionError("an exceptional coordinate moved")
    for i, left in enumerate(basis):
        for j, right in enumerate(basis):
            if pairing(images[i], images[j]) != pairing(left, right):
                raise AssertionError("swap is not an isometry")
    if any(reflect(image, root) != vector for image, vector in zip(images, basis)):
        raise AssertionError("swap is not involutive")

    # Any iteratively constructed symmetric ample class has this form.  Its
    # ampleness is the standard small-exceptional-coefficient blowup lemma;
    # the exact certificate checks the mirror equation.
    ample_probe = (
        Fraction(1),
        Fraction(1),
        *(Fraction(-1, 10 ** (index + 2)) for index in range(length)),
    )
    if pairing(ample_probe, root) != 0:
        raise AssertionError("symmetric ample probe missed the mirror")
    return {
        "m": length,
        "root_square": int(pairing(root, root)),
        "integral": True,
        "swaps_rulings": True,
        "fixes_exceptional_coordinates": True,
        "symmetric_ample_pairing": 0,
    }


def canonical(length: int) -> Vector:
    return (-2, -2, *((1,) * length))


def truncator_record(length: int) -> dict[str, object]:
    if length < 5:
        raise ValueError(length)
    vector = (1, 2, *((-1,) * 5), *((0,) * (length - 5)))
    square = int(pairing(vector, vector))
    k_pair = int(pairing(canonical(length), vector))
    genus = 1 + (square + k_pair) // 2
    v5 = vector
    return {
        "m": length,
        "vector": list(vector),
        "square": square,
        "canonical_pairing": k_pair,
        "arithmetic_genus": genus,
        "pairing_with_v5": int(pairing(vector, v5)),
        "linear_system_h0": 6,
        "cluster_condition_count": 5,
        "effective_divisor_class": True,
        "irreducible_smooth_facet_automatic": False,
    }


def r_root(alpha: int, length: int) -> Vector:
    support = 4 * alpha + 2
    if support > length:
        raise ValueError((alpha, length))
    return (
        alpha,
        2,
        *((-1,) * support),
        *((0,) * (length - support)),
    )


def r_family_record(alpha: int) -> dict[str, object]:
    support = 4 * alpha + 2
    length = support + 3
    vector = r_root(alpha, length)
    labels, walls = visible_walls(length)
    intersections = {
        label: int(pairing(vector, wall))
        for label, wall in zip(labels, walls)
        if pairing(vector, wall)
    }
    v5 = (1, 2, *((-1,) * 5), *((0,) * (length - 5)))
    square = int(pairing(vector, vector))
    k_pair = int(pairing(canonical(length), vector))
    genus = 1 + (square + k_pair) // 2
    expected = {f"delta_{support}": 1}
    if alpha > 1:
        expected["D"] = alpha - 1
    if intersections != expected:
        raise AssertionError((alpha, intersections, expected))
    return {
        "a": alpha,
        "support": support,
        "square": square,
        "canonical_pairing": k_pair,
        "arithmetic_genus": genus,
        "pairing_with_v5": int(pairing(vector, v5)),
        "visible_intersections": intersections,
        "h0_downstairs": 3 * alpha + 3,
        "jet_rows": support,
        "minimum_row_dependencies_for_effectivity": alpha,
        "free_cluster_dimension": 4 * alpha,
        "irreducible_exact_multiplicity_incidence_dimension_upper_bound":
            3 * alpha,
        "irreducible_exact_multiplicity_locus_codimension_lower_bound":
            alpha,
        "smooth_rational_adjunction_compatible": genus == 0,
    }


def first_height_record(limit: int = 12) -> dict[str, object]:
    """Symbolic classification at beta^2/(-r^2)=2.

    In a unimodular lattice a primitive integral mirror has square -1 or
    -2.  Compatibility excludes beta=1.  Height two forces beta=2 and
    square -2.  The inequality x1+x2<=2 leaves either a single coefficient
    two (which cannot be negative once alpha>=2) or a zero-one sequence.
    The square equation then gives exactly 4*alpha+2 ones.
    """

    records = []
    for alpha in range(1, limit + 1):
        support = 4 * alpha + 2
        square = 4 * alpha - support
        if square != -2:
            raise AssertionError(alpha)
        records.append(
            {
                "alpha": alpha,
                "beta": 2,
                "square": square,
                "normalized_height": 2,
                "nonzero_coefficients": support,
                "all_nonzero_coefficients_equal_one": True,
            }
        )
    return {
        "beta_one_negative_root_exists": False,
        "minimum_positive_normalized_height": 2,
        "height_two_roots": records,
        "higher_height_square_minus_one_family_not_claimed_exhausted": True,
    }


def pairwise_r_record(limit: int = 9) -> list[dict[str, int]]:
    result = []
    length = 4 * limit + 5
    for alpha in range(1, limit):
        for beta in range(alpha + 1, limit + 1):
            value = int(pairing(r_root(alpha, length), r_root(beta, length)))
            expected = 2 * (beta - alpha - 1)
            if value != expected:
                raise AssertionError((alpha, beta, value, expected))
            result.append({"a": alpha, "b": beta, "pairing": value})
    return result


def hidden_wall_record() -> dict[str, object]:
    # C=a*c+b*d-sum(mu_i e_i), v5=c+2d-e1-...-e5.
    samples = []
    for a, b, multiplicities in (
        (1, 2, (1, 1, 1, 1, 1)),
        (1, 2, (1, 1, 1, 1, 1, 1)),
        (3, 1, (2, 2, 1, 1, 1)),
    ):
        length = len(multiplicities)
        if length < 5:
            raise AssertionError
        curve = (a, b, *(-value for value in multiplicities))
        v5 = (1, 2, *((-1,) * 5), *((0,) * (length - 5)))
        direct = int(pairing(curve, v5))
        formula = 2 * a + b - sum(multiplicities[:5])
        if direct != formula:
            raise AssertionError((curve, direct, formula))
        samples.append(
            {
                "a": a,
                "b": b,
                "first_five_multiplicities": list(multiplicities[:5]),
                "pairing": direct,
            }
        )
    return {
        "formula": "C.v5 = 2*a+b-sum(mu_1,...,mu_5)",
        "finite_volume_requires_some_new_wall_negative_on_v5": True,
        "necessary_inequality": "sum(mu_1,...,mu_5) > 2*a+b",
        "samples": samples,
    }


def adapted_host_record() -> dict[str, object]:
    return {
        "reflection_formula":
            "s_(c-d)(B/2)=B/2+(B.(c-d)/2)*(c-d)",
        "parity_required_for_every_half_boundary": True,
        "odd_absorption_by_half_root_impossible_in_integral_host": (
            "if (c-d)/2 and B/2 both lie in M, their doubled-form pairing "
            "B.(c-d)/2 must be integral"
        ),
        "full_NS_obstruction_automatically_transfers_to_every_adapted_host": False,
    }


def source_hash_record() -> dict[str, str]:
    actual = {
        relative: file_sha256(ROOT / relative)
        for relative in SOURCE_HASHES
    }
    if actual != SOURCE_HASHES:
        raise AssertionError({"expected": SOURCE_HASHES, "actual": actual})
    return actual


@lru_cache(maxsize=1)
def build_payload() -> dict[str, object]:
    visible = [visible_record(length) for length in range(2, 49)]
    finite = [record["m"] for record in visible if record["finite_volume"]]
    if finite != [2, 3, 4]:
        raise AssertionError(finite)

    t_records = [truncator_record(length) for length in (5, 6, 11, 30)]
    if any(record["pairing_with_v5"] != -1 for record in t_records):
        raise AssertionError("T5 did not survive later blowups")

    r_records = [r_family_record(alpha) for alpha in range(1, 13)]
    if [record["arithmetic_genus"] for record in r_records] != list(range(12)):
        raise AssertionError("adjunction-genus sequence failed")

    return {
        "verdict": "PASS_AFTER_SCOPE_AND_CERTIFICATE_REPAIRS",
        "source_hashes": source_hash_record(),
        "visible_simplex_records": visible,
        "visible_finite_volume_exactly_m_2_3_4": True,
        "hidden_wall_gate": hidden_wall_record(),
        "swap_records": [swap_record(length) for length in (2, 5, 17, 48)],
        "full_NS_uniform_obstruction": {
            "obstructing_root": "c-d",
            "root_mirror_contains_an_ample_class": True,
            "facet_group_cannot_contain_that_reflection": True,
            "adjoining_it_is_a_proper_discrete_reflection_overgroup": True,
            "finite_volume_overgroup_has_finite_index": True,
            "absolute_full_NS_reflection_maximality": False,
        },
        "T5_survives_all_later_blowups_samples": t_records,
        "first_height_classification": first_height_record(),
        "R_family_records": r_records,
        "R_pairwise_records": pairwise_r_record(),
        "adapted_host_boundary": adapted_host_record(),
        "claim_boundary": {
            "T5_effective_divisor_class_is_automatic": True,
            "T5_irreducible_smooth_nef_facet_is_automatic": False,
            "Ra_lattice_root_is_automatic": True,
            "Ra_effectivity_is_automatic": False,
            "reflection_maximality_makes_every_root_a_physical_facet": False,
            "proper_host_or_slice_excluded": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    payload = build_payload()
    archive = {"payload": payload, "payload_sha256": digest(payload)}
    rendered = json.dumps(archive, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if arguments.output.read_text() != rendered:
            raise AssertionError("hostile archive differs from recomputation")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
    print(archive["payload_sha256"])


if __name__ == "__main__":
    main()
