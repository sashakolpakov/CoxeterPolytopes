#!/usr/bin/env python3
"""Polytope-first OG10 carrier audit for the Borcherds D4 H21/210 chamber.

The input is only Alexeev's frozen copy of the 210-by-210 Coxeter Gram
matrix for the classical ``(22,2,0)`` chamber.  The replay reconstructs its
rank-22 lattice, enumerates all discriminant gluings into the OG10
Beauville lattice, and then applies the OG10 wall-divisor list.

The outcome has two deliberately separate parts.

* A primitive full-Picard OG10 period carrier exists and its embedding orbit
  is unique.
* The displayed 210-wall Coxeter chamber cannot be a complete Kahler/nef
  chamber: an additional square-minus-four, divisibility-one wall passes
  through a strictly interior point.

This is a POLYTOPE -> VARIETY computation.  It is not a nonexistence theorem
for the classical hyperbolic polytope and imports no surface-first wall list.
"""

from __future__ import annotations

import argparse
from collections import Counter
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

from scripts.certify_alexeev_high_dim_reflection_examples import (
    SOURCE_SHA256,
    _extract_matrix,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = Path(
    "/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt"
)
OUTPUT = (
    ROOT
    / "data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json"
)
SCHEMA = "coxeter.borcherds_d4_h21.og10_polytope_first.v1"
KEY = (22, 2, 0)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _rational(value: sp.Expr) -> str:
    value = sp.Rational(value)
    return str(int(value)) if value.q == 1 else f"{value.p}/{value.q}"


def _matrix_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [_rational(matrix[row, col]) for col in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _vector_payload(vector: Iterable[sp.Expr]) -> list[str]:
    return [_rational(value) for value in vector]


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index]
    ]


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Return exact rational inertia without floating-point eigenvalues."""

    work = sp.Matrix(matrix)
    positive = negative = 0
    while work.rows:
        pivot = next((i for i in range(work.rows) if work[i, i]), None)
        if pivot is None:
            pair = next(
                (
                    (i, j)
                    for i in range(work.rows)
                    for j in range(i)
                    if work[i, j]
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            i, j = pair
            transform = sp.eye(work.rows)
            transform[:, i] = sp.eye(work.rows)[:, i] + sp.eye(work.rows)[:, j]
            work = transform.T * work * transform
            pivot = i
        order = [pivot, *range(pivot), *range(pivot + 1, work.rows)]
        work = work.extract(order, order)
        value = work[0, 0]
        if value > 0:
            positive += 1
        else:
            negative += 1
        if work.rows == 1:
            break
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / value
    return positive, negative, 0


def _kernel_mod_prime(matrix: sp.Matrix, prime: int) -> list[list[int]]:
    """RREF kernel basis over F_prime, represented by least residues."""

    rows, cols = matrix.shape
    work = [
        [int(matrix[row, col]) % prime for col in range(cols)]
        for row in range(rows)
    ]
    pivot_columns: list[int] = []
    pivot_row = 0
    for column in range(cols):
        candidate = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if candidate is None:
            continue
        work[pivot_row], work[candidate] = work[candidate], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [
            (entry * inverse) % prime for entry in work[pivot_row]
        ]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                (left - multiple * right) % prime
                for left, right in zip(work[row], work[pivot_row])
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == rows:
            break
    free_columns = [
        column for column in range(cols) if column not in pivot_columns
    ]
    basis: list[list[int]] = []
    for free in free_columns:
        vector = [0] * cols
        vector[free] = 1
        for row, pivot in enumerate(pivot_columns):
            vector[pivot] = (-work[row][free]) % prime
        if any(
            sum(int(matrix[row, col]) * vector[col] for col in range(cols))
            % prime
            for row in range(rows)
        ):
            raise ArithmeticError("finite-field kernel reconstruction failed")
        basis.append(vector)
    return basis


def _xor(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    return tuple(a ^ b for a, b in zip(left, right))


def _two_torsion(
    gram: sp.Matrix,
) -> tuple[list[list[int]], list[dict[str, Any]]]:
    basis = _kernel_mod_prime(gram, 2)
    if len(basis) != 2:
        raise ArithmeticError("expected a two-dimensional 2-torsion kernel")
    zero = (0,) * gram.rows
    first, second = (tuple(vector) for vector in basis)
    bit_vectors = [first, second, _xor(first, second)]
    classes: list[dict[str, Any]] = []
    for bits in bit_vectors:
        representative = sp.Matrix(bits) / 2
        dual_pairings = gram * representative
        square = sp.Rational((representative.T * gram * representative)[0])
        classes.append(
            {
                "bits": list(bits),
                "representative": _vector_payload(representative),
                "dual_pairings_integral": all(value.q == 1 for value in dual_pairings),
                "square": _rational(square),
                "q_mod_2": _rational(square % 2),
                "order_two": _xor(bits, bits) == zero,
            }
        )
    return basis, classes


def _gluing_graph_census(
    l_nonzero: Sequence[Sequence[int]],
    k_nonzero: Sequence[Sequence[int]],
) -> tuple[int, int]:
    """Count ordered generator presentations and distinct graph subgroups."""

    presentations = 0
    subgroups: set[tuple[tuple[int, ...], ...]] = set()
    zero_l = (0,) * len(l_nonzero[0])
    zero_k = (0,) * len(k_nonzero[0])
    for left_basis in itertools.permutations(l_nonzero, 2):
        if tuple(left_basis[0]) == tuple(left_basis[1]):
            continue
        for right_basis in itertools.permutations(k_nonzero, 2):
            if tuple(right_basis[0]) == tuple(right_basis[1]):
                continue
            presentations += 1
            graph = (
                zero_l + zero_k,
                tuple(left_basis[0]) + tuple(right_basis[0]),
                tuple(left_basis[1]) + tuple(right_basis[1]),
                _xor(left_basis[0], left_basis[1])
                + _xor(right_basis[0], right_basis[1]),
            )
            subgroups.add(tuple(sorted(graph)))
    return presentations, len(subgroups)


def _divisibility(gram: sp.Matrix, vector: sp.Matrix) -> int:
    values = [abs(int(value)) for value in gram * vector]
    return math.gcd(*values)


def _reduced_even_binary_det_twelve() -> list[tuple[int, int, int]]:
    forms: list[tuple[int, int, int]] = []
    for a in range(1, 13):
        for c in range(a, 13):
            for b in range(-a, a + 1):
                if 4 * a * c - b * b != 12:
                    continue
                if (abs(b) == a or a == c) and b < 0:
                    continue
                forms.append((a, b, c))
    return forms


def _binary_form_audit(form: tuple[int, int, int]) -> dict[str, Any]:
    a, b, c = form
    gram = sp.Matrix([[2 * a, b], [b, 2 * c]])
    classes = []
    for bits in ((1, 0), (0, 1), (1, 1)):
        representative = sp.Matrix(bits) / 2
        square = sp.Rational((representative.T * gram * representative)[0])
        classes.append(_rational(square % 2))
    return {
        "reduced_form_abc": list(form),
        "gram": _matrix_payload(gram),
        "determinant": int(gram.det()),
        "two_torsion_q_mod_2": classes,
        "compatible_with_L_two_primary_form": classes == ["1", "1", "1"],
    }


def _a2_two_automorphisms(gram: sp.Matrix) -> tuple[int, int]:
    nonzero = ((1, 0), (0, 1), (1, 1))
    isometries: list[sp.Matrix] = []
    permutations: set[tuple[int, ...]] = set()
    for entries in itertools.product(range(-2, 3), repeat=4):
        transform = sp.Matrix(2, 2, entries)
        if abs(int(transform.det())) != 1 or transform.T * gram * transform != gram:
            continue
        isometries.append(transform)
        permutations.add(
            tuple(
                nonzero.index(
                    tuple(int(value) % 2 for value in transform * sp.Matrix(bits))
                )
                for bits in nonzero
            )
        )
    return len(isometries), len(permutations)


def _counter_payload(counter: Counter[tuple[int, int]]) -> dict[str, int]:
    return {
        f"square_{square}_div_{divisibility}": count
        for (square, divisibility), count in sorted(counter.items())
    }


def build_certificate(source_path: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    raw = source_path.read_bytes()
    source_sha256 = hashlib.sha256(raw).hexdigest()
    source = raw.decode("utf-8")
    coxeter_gram = _extract_matrix(source, KEY)
    if coxeter_gram.shape != (210, 210):
        raise ArithmeticError("the frozen Coxeter Gram is not 210 by 210")

    # The first 22 walls form the canonical nonsingular block in the source.
    lattice_gram = coxeter_gram[:22, :22]
    wall_coordinates = lattice_gram.inv() * coxeter_gram[:22, :]
    coordinate_integrality = all(value.q == 1 for value in wall_coordinates)
    gram_reconstruction = (
        wall_coordinates.T * lattice_gram * wall_coordinates == coxeter_gram
    )

    l_basis, l_classes = _two_torsion(lattice_gram)
    complement_gram = sp.Matrix([[4, -2], [-2, 4]])
    k_basis, k_classes = _two_torsion(complement_gram)
    ordered_gluings, distinct_gluings = _gluing_graph_census(
        [item["bits"] for item in l_classes],
        [item["bits"] for item in k_classes],
    )

    # Choose the graph sending the finite-field kernel basis of L to the
    # standard kernel basis of A2(2).  Its two glue vectors, together with
    # the 22 standard L vectors, are a basis of the index-four overlattice.
    direct_gram = sp.diag(lattice_gram, complement_gram)
    glue_columns = []
    for left_bits, right_bits in zip(l_basis, k_basis):
        glue_columns.append(
            sp.Matrix.vstack(sp.Matrix(left_bits) / 2, sp.Matrix(right_bits) / 2)
        )
    direct_l_basis = sp.Matrix.vstack(sp.eye(22), sp.zeros(2, 22))
    overlattice_basis = sp.Matrix.hstack(direct_l_basis, *glue_columns)
    ambient_gram = sp.simplify(overlattice_basis.T * direct_gram * overlattice_basis)
    l_embedding = sp.simplify(overlattice_basis.inv() * direct_l_basis)

    ambient_discriminant_kernel = _kernel_mod_prime(ambient_gram, 3)
    if len(ambient_discriminant_kernel) != 1:
        raise ArithmeticError("ambient discriminant is not cyclic of order three")
    ambient_discriminant_representative = (
        sp.Matrix(ambient_discriminant_kernel[0]) / 3
    )
    ambient_discriminant_square = sp.Rational(
        (
            ambient_discriminant_representative.T
            * ambient_gram
            * ambient_discriminant_representative
        )[0]
    )

    ambient_walls = l_embedding * wall_coordinates
    l_wall_types: Counter[tuple[int, int]] = Counter()
    ambient_wall_types: Counter[tuple[int, int]] = Counter()
    for column in range(ambient_walls.cols):
        l_wall = wall_coordinates[:, column]
        ambient_wall = ambient_walls[:, column]
        square = int((l_wall.T * lattice_gram * l_wall)[0])
        l_wall_types[(square, _divisibility(lattice_gram, l_wall))] += 1
        ambient_wall_types[
            (square, _divisibility(ambient_gram, ambient_wall))
        ] += 1

    reduced_forms = _reduced_even_binary_det_twelve()
    complement_audits = [_binary_form_audit(form) for form in reduced_forms]
    a2_isometry_order, a2_discriminant_action_order = _a2_two_automorphisms(
        complement_gram
    )

    # A P-first hidden-wall witness.  The positive vector x pairs strictly
    # positively with every one of the 210 archived facets, while the
    # additional primitive vector v is orthogonal to x.
    hidden_root = sp.Matrix([-1, 1, *([0] * 20)])
    interior_vector = lattice_gram.inv() * sp.ones(22, 1)
    facet_pairings = [
        sp.Rational((interior_vector.T * lattice_gram * wall_coordinates[:, j])[0])
        for j in range(210)
    ]
    hidden_ambient = l_embedding * hidden_root
    hidden_is_listed = any(
        wall_coordinates[:, column] in (hidden_root, -hidden_root)
        for column in range(wall_coordinates.cols)
    )

    no_divisibility_three_argument = {
        "determinant_mod_3": int(lattice_gram.det()) % 3,
        "gram_invertible_mod_3": int(lattice_gram.det()) % 3 != 0,
        "argument": (
            "If a primitive c in L had ambient divisibility 3, then Gc=0 "
            "mod 3 because L is a primitive sublattice.  Since det(G)=-4 "
            "is a unit mod 3, c=0 mod 3, contradicting primitivity."
        ),
        "conclusion": (
            "No primitive class in this Picard lattice has ambient "
            "divisibility 3; the OG10 (-6,3) and (-24,3) wall types cannot "
            "occur inside L."
        ),
    }

    checks = {
        "source_digest_pinned": source_sha256 == SOURCE_SHA256,
        "coxeter_gram_symmetric": coxeter_gram == coxeter_gram.T,
        "coxeter_gram_rank_22": coxeter_gram.rank() == 22,
        "coxeter_inertia_1_21_188": _exact_inertia(coxeter_gram) == (1, 21, 188),
        "leading_lattice_det_minus_4": lattice_gram.det() == -4,
        "leading_lattice_inertia_1_21": _exact_inertia(lattice_gram) == (1, 21, 0),
        "all_wall_coordinates_integral": coordinate_integrality,
        "all_210_pairings_reconstructed": gram_reconstruction,
        "L_discriminant_is_2_by_2": _smith_diagonal(lattice_gram) == [1] * 20 + [2, 2],
        "L_nonzero_discriminant_values_all_minus_1": [
            item["square"] for item in l_classes
        ]
        == ["-1", "-1", "-1"],
        "K_is_A2_of_scale_2": complement_gram.det() == 12,
        "K_nonzero_two_torsion_values_all_plus_1": [
            item["square"] for item in k_classes
        ]
        == ["1", "1", "1"],
        "six_discriminant_graphs": ordered_gluings == 36 and distinct_gluings == 6,
        "overlattice_basis_index_four": abs(overlattice_basis.det()) == sp.Rational(1, 4),
        "ambient_even_integral": all(value.q == 1 for value in ambient_gram)
        and all(ambient_gram[i, i] % 2 == 0 for i in range(24)),
        "ambient_signature_3_21": _exact_inertia(ambient_gram) == (3, 21, 0),
        "ambient_det_minus_3": ambient_gram.det() == -3,
        "ambient_smith_1_23_3": _smith_diagonal(ambient_gram) == [1] * 23 + [3],
        "ambient_discriminant_q_4_over_3": ambient_discriminant_square % 2
        == sp.Rational(4, 3),
        "L_embedding_integral_and_primitive": all(value.q == 1 for value in l_embedding)
        and _smith_diagonal(l_embedding) == [1] * 22,
        "wall_types_in_L_are_42_short_168_long": l_wall_types
        == Counter({(-2, 1): 42, (-4, 2): 168}),
        "all_210_are_OG10_wall_types": ambient_wall_types
        == Counter({(-2, 1): 42, (-4, 1): 168}),
        "rank_two_complements_exhausted": reduced_forms == [(1, 0, 3), (2, 2, 2)],
        "only_A2_scale_2_is_glue_compatible": [
            item["compatible_with_L_two_primary_form"]
            for item in complement_audits
        ]
        == [False, True],
        "A2_scale_2_acts_as_full_S3_on_two_torsion": a2_isometry_order == 12
        and a2_discriminant_action_order == 6,
        "no_divisibility_three_classes_in_L": int(lattice_gram.det()) % 3 != 0,
        "interior_vector_positive": (interior_vector.T * lattice_gram * interior_vector)[0]
        == 581,
        "all_210_facet_pairings_strict": min(facet_pairings) == 1
        and max(facet_pairings) == 421
        and all(value > 0 for value in facet_pairings),
        "hidden_root_is_new_OG10_wall": not hidden_is_listed
        and (hidden_root.T * lattice_gram * hidden_root)[0] == -4
        and _divisibility(ambient_gram, hidden_ambient) == 1,
        "hidden_wall_crosses_strict_interior": (
            interior_vector.T * lattice_gram * hidden_root
        )[0]
        == 0,
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "directional_accounting": {
            "direction": "POLYTOPE -> OG10 VARIETY -> POLYTOPE-CONE CRITERION",
            "first_input": (
                "the frozen 210-by-210 Coxeter Gram matrix; no variety or "
                "surface wall list is used to construct the candidate"
            ),
            "positive_result": (
                "the full-Picard OG10 period-carrier moduli problem is nonempty"
            ),
            "negative_result": (
                "the moduli problem requiring these 210 facets to be the "
                "complete Kahler/nef chamber is empty"
            ),
            "scope_warning": (
                "This does not obstruct the Borcherds hyperbolic polytope, "
                "nor does failure of the OG10 carrier exhaust other carrier species."
            ),
        },
        "source": {
            "matrix_key": list(KEY),
            "source_file_name": source_path.name,
            "source_sha256": source_sha256,
            "expected_sha256": SOURCE_SHA256,
            "url": "https://arxiv.org/abs/2209.09110",
        },
        "coxeter_input": {
            "ledger_id": "BOR-D4-H21-210",
            "dimension": 21,
            "facet_count": 210,
            "gram_rank": int(coxeter_gram.rank()),
            "gram_inertia_positive_negative_zero": list(_exact_inertia(coxeter_gram)),
            "basis_facet_indices_one_based": list(range(1, 23)),
            "basis_gram": _matrix_payload(lattice_gram),
            "basis_gram_sha256": _canonical_sha256(_matrix_payload(lattice_gram)),
            "wall_coordinates_as_columns_transposed": _matrix_payload(
                wall_coordinates.T
            ),
            "wall_coordinates_sha256": _canonical_sha256(
                _matrix_payload(wall_coordinates.T)
            ),
            "wall_norm_counts": {"-2": 42, "-4": 168},
            "full_gram_reconstructed": gram_reconstruction,
        },
        "picard_lattice_L": {
            "identification": "U + D20(-1)",
            "rank": 22,
            "signature": [1, 21],
            "determinant": int(lattice_gram.det()),
            "smith_diagonal": _smith_diagonal(lattice_gram),
            "two_torsion_kernel_basis_bits": l_basis,
            "nonzero_discriminant_classes": l_classes,
            "wall_type_counts": _counter_payload(l_wall_types),
        },
        "primitive_OG10_embedding": {
            "positive_complement": "A2(2)",
            "complement_gram": _matrix_payload(complement_gram),
            "complement_two_torsion_kernel_basis_bits": k_basis,
            "complement_nonzero_two_torsion_classes": k_classes,
            "ordered_generator_presentations": ordered_gluings,
            "distinct_isotropic_graph_subgroups": distinct_gluings,
            "chosen_glue_vectors_in_L_plus_K": [
                _vector_payload(vector) for vector in glue_columns
            ],
            "overlattice_basis_determinant": _rational(overlattice_basis.det()),
            "ambient_gram": _matrix_payload(ambient_gram),
            "ambient_gram_sha256": _canonical_sha256(_matrix_payload(ambient_gram)),
            "ambient_signature": [3, 21],
            "ambient_determinant": int(ambient_gram.det()),
            "ambient_smith_diagonal": _smith_diagonal(ambient_gram),
            "ambient_discriminant_generator": _vector_payload(
                ambient_discriminant_representative
            ),
            "ambient_discriminant_square": _rational(ambient_discriminant_square),
            "ambient_discriminant_q_mod_2": _rational(
                ambient_discriminant_square % 2
            ),
            "ambient_identification": (
                "3U + 2E8(-1) + A2(-1), the OG10 Beauville lattice"
            ),
            "L_embedding_columns": _matrix_payload(l_embedding),
            "L_embedding_smith_diagonal": _smith_diagonal(l_embedding),
            "ambient_wall_type_counts": _counter_payload(ambient_wall_types),
            "period_statement": (
                "Period-map surjectivity gives marked OG10 manifolds with "
                "exact Picard lattice L; the orthogonal positive rank-two "
                "complement makes this a zero-dimensional period stratum."
            ),
        },
        "embedding_orbit_exhaustion": {
            "determinant_formula": (
                "The 2-primary discriminant of L must be killed completely. "
                "Thus the glue has order 4 and det(K)=3*4^2/4=12."
            ),
            "reduced_even_positive_binary_forms_of_determinant_12": complement_audits,
            "compatible_complement_count": sum(
                item["compatible_with_L_two_primary_form"]
                for item in complement_audits
            ),
            "A2_2_integral_isometry_group_order": a2_isometry_order,
            "induced_two_torsion_permutation_group_order": a2_discriminant_action_order,
            "conclusion": (
                "A2(2) is the unique complement and its isometry group is "
                "surjective on the six anti-isometries; hence there is one "
                "primitive embedding orbit."
            ),
        },
        "OG10_wall_divisor_gate": {
            "published_wall_types": [
                {"square": -2, "divisibility": 1},
                {"square": -4, "divisibility": 1},
                {"square": -6, "divisibility": 3},
                {"square": -24, "divisibility": 3},
            ],
            "all_210_input_facets_are_wall_divisors": True,
            "no_divisibility_three_argument": no_divisibility_three_argument,
            "criterion": (
                "A full Kahler/nef promotion is possible iff the proposed cone "
                "is a component of the positive cone after removing every "
                "OG10 wall-divisor hyperplane.  One crossing hidden wall is a "
                "complete obstruction."
            ),
        },
        "hidden_wall_certificate": {
            "hidden_root_in_L_basis": _vector_payload(hidden_root),
            "hidden_root_square": int(
                (hidden_root.T * lattice_gram * hidden_root)[0]
            ),
            "hidden_root_ambient_divisibility": _divisibility(
                ambient_gram, hidden_ambient
            ),
            "hidden_root_absent_from_210_facets": not hidden_is_listed,
            "strict_interior_vector_in_L_basis": _vector_payload(interior_vector),
            "interior_square": int(
                (interior_vector.T * lattice_gram * interior_vector)[0]
            ),
            "input_facet_pairing_minimum": _rational(min(facet_pairings)),
            "input_facet_pairing_maximum": _rational(max(facet_pairings)),
            "hidden_root_pairing": _rational(
                (interior_vector.T * lattice_gram * hidden_root)[0]
            ),
            "conclusion": (
                "The square-minus-four divisibility-one OG10 wall v^perp "
                "contains a point strictly interior to all 210 input halfspaces. "
                "Therefore the 210-wall chamber is not a complete Kahler/nef cone."
            ),
        },
        "literature_inputs": [
            {
                "url": "https://arxiv.org/abs/math/0606409",
                "use": "Rapagnetta's OG10 Beauville lattice and b2=24",
            },
            {
                "url": "https://arxiv.org/abs/2010.12511",
                "use": (
                    "Mongardi--Onorati's OG10 wall-divisor classification and "
                    "Kahler-cone component theorem"
                ),
            },
            {
                "url": "https://arxiv.org/abs/2209.09110",
                "use": "the pinned (22,2,0) Coxeter Gram matrix",
            },
        ],
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", nargs="?", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    payload = build_certificate(args.source)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise AssertionError("archived OG10 carrier audit has drifted")
        print(f"archive check passed: {OUTPUT}")
        print(f"certificate_sha256={payload['certificate_sha256']}")
        return
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)
        print(payload["certificate_sha256"])
        return
    print(rendered, end="")


if __name__ == "__main__":
    main()
