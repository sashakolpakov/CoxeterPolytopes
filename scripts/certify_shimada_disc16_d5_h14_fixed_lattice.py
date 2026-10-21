"""Certify the saturated fixed lattice and mirrors of disc16 D5 H14.

The involution and all facet orbit sums come from the independent primary-data
screen.  The row-coordinate fixed lattice is the integral kernel of
``M^T-I``.  A Smith decomposition supplies a saturated basis; Euclidean LLL
is used only for a smaller unimodular basis of the same lattice.

Each facet normal is converted from Shimada's dual-basis coordinates into the
fixed-lattice basis, replaced by the primitive integral vector defining the
same mirror, and its reflection is replayed as an integral isometry.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from scripts.independent_certify_shimada_disc16_d5_h14 import (
    SUMMARY_OUTPUT as VOLUME_CERTIFICATE,
)
from scripts.independent_screen_shimada_disc16_d5_h14 import (
    OUTPUT as SCREEN_CERTIFICATE,
    ROOT,
)
from scripts.search_jacobian_kummer_chamber_folds import (
    _matrix_payload,
    _rational_text,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h14_fixed_lattice_certificate.json"
)


def _matrix(payload: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(entry) for entry in row] for row in payload])


def _smith_decomposition(matrix: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    diagonal, left, right = smith_normal_decomp(
        DomainMatrix.from_Matrix(matrix).convert_to(ZZ)
    )
    return diagonal.to_Matrix(), left.to_Matrix(), right.to_Matrix()


def _gf2_nullspace(matrix: sp.Matrix) -> tuple[tuple[int, ...], ...]:
    rows, columns = matrix.shape
    work = [
        [int(matrix[row, column]) % 2 for column in range(columns)]
        for row in range(rows)
    ]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(columns):
        swap = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if swap is None:
            continue
        work[pivot_row], work[swap] = work[swap], work[pivot_row]
        for row in range(rows):
            if row != pivot_row and work[row][column]:
                work[row] = [
                    left ^ right
                    for left, right in zip(work[row], work[pivot_row], strict=True)
                ]
        pivots.append(column)
        pivot_row += 1
    free = [column for column in range(columns) if column not in pivots]
    basis: list[tuple[int, ...]] = []
    for free_column in free:
        vector = [0] * columns
        vector[free_column] = 1
        for row, pivot in reversed(list(enumerate(pivots))):
            vector[pivot] = sum(
                work[row][column] * vector[column] for column in free
            ) % 2
        basis.append(tuple(vector))
    return tuple(basis)


def _mod_rational(value: sp.Rational, modulus: int) -> sp.Rational:
    value = sp.Rational(value)
    return sp.Rational(value - modulus * sp.floor(value / modulus))


def _primitive(vector: sp.Matrix) -> tuple[sp.Matrix, sp.Rational]:
    denominators = [int(sp.Rational(value).q) for value in vector]
    denominator = int(sp.ilcm(*denominators))
    integers = [int(sp.Rational(value) * denominator) for value in vector]
    common = math.gcd(*(abs(value) for value in integers))
    if common == 0:
        raise ArithmeticError("zero mirror normal")
    primitive = [value // common for value in integers]
    first = next(value for value in primitive if value)
    if first < 0:
        primitive = [-value for value in primitive]
    primitive_matrix = sp.Matrix(primitive)
    scale = next(
        sp.Rational(vector[index], primitive_matrix[index])
        for index in range(vector.rows)
        if primitive_matrix[index]
    )
    if vector != scale * primitive_matrix:
        raise ArithmeticError("failed to primitive-normalize a mirror")
    return primitive_matrix, scale


def build_certificate(
    screen_path: Path = SCREEN_CERTIFICATE,
    volume_path: Path = VOLUME_CERTIFICATE,
) -> dict[str, Any]:
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    volume = json.loads(volume_path.read_text(encoding="utf-8"))
    if not (
        screen["facet_count"] == 60
        and screen["representative"]["fixed_rank"] == 15
        and volume["decision"]["finite_volume"] is True
        and volume["exact_double_description"]["hyperideal_vertex_count"] == 0
    ):
        raise ArithmeticError("the fixed-lattice inputs are not the certified H14 chamber")

    ambient_gram = _matrix(screen["ambient_lattice"]["gram"])
    ambient_dual_gram = ambient_gram.inv()
    involution = _matrix(screen["representative"]["raw_primal_row_matrix"])
    constraint = involution.T - sp.eye(20)
    smith_constraint, _, right = _smith_decomposition(constraint)
    constraint_diagonal = tuple(smith_constraint[index, index] for index in range(20))
    zero_columns = tuple(
        column
        for column in range(20)
        if all(smith_constraint[row, column] == 0 for row in range(20))
    )
    if not (
        constraint_diagonal[:5] == (1, 1, 1, 1, 1)
        and all(value == 0 for value in constraint_diagonal[5:])
        and zero_columns == tuple(range(5, 20))
    ):
        raise ArithmeticError("the fixed-kernel Smith form changed")
    smith_kernel_basis = right[:, list(zero_columns)].T
    reduced_basis = smith_kernel_basis.lll(delta=sp.Rational(3, 4))
    _, pivot_columns = smith_kernel_basis.rref()
    pivot_columns = tuple(pivot_columns)
    basis_change = (
        reduced_basis[:, list(pivot_columns)]
        * smith_kernel_basis[:, list(pivot_columns)].inv()
    )
    if not (
        basis_change * smith_kernel_basis == reduced_basis
        and abs(int(basis_change.det())) == 1
        and reduced_basis * involution == reduced_basis
        and reduced_basis.rank() == 15
    ):
        raise ArithmeticError("LLL did not preserve the saturated fixed lattice")

    fixed_gram = reduced_basis * ambient_gram * reduced_basis.T
    fixed_snf = smith_normal_form(fixed_gram, domain=ZZ)
    fixed_snf_diagonal = tuple(abs(int(fixed_snf[index, index])) for index in range(15))
    if not (
        fixed_gram.det() == 128
        and exact_inertia(fixed_gram) == (1, 14, 0)
        and all(int(fixed_gram[index, index]) % 2 == 0 for index in range(15))
        and fixed_snf_diagonal == (1,) * 8 + (2,) * 7
    ):
        raise ArithmeticError("the saturated fixed-lattice invariants changed")

    discriminant_bits = _gf2_nullspace(fixed_gram)
    if len(discriminant_bits) != 7:
        raise ArithmeticError("the mod-2 dual quotient has wrong dimension")
    discriminant_columns = sp.Matrix.hstack(
        *(sp.Matrix(vector) for vector in discriminant_bits)
    )
    if any(int(value) % 2 for value in fixed_gram * discriminant_columns):
        raise ArithmeticError("a proposed half-vector is not in the dual lattice")
    discriminant_form = (
        discriminant_columns.T * fixed_gram * discriminant_columns / 4
    )
    discriminant_form_mod = [
        [
            _rational_text(
                _mod_rational(
                    discriminant_form[left, right], 2 if left == right else 1
                )
            )
            for right in range(7)
        ]
        for left in range(7)
    ]
    delta = int(
        any(
            sp.Rational(discriminant_form[index, index]).q == 2
            for index in range(7)
        )
    )
    if delta != 1:
        raise ArithmeticError("the fixed lattice has unexpected 2-elementary parity")

    walls = tuple(
        sp.Matrix(1, 20, wall)
        for wall in screen["chamber"]["published_primitive_dual_walls"]
    )
    wall_orbits = tuple(
        tuple(index - 1 for index in orbit) for orbit in screen["wall_orbits"]
    )
    facet_orbits = tuple(index - 1 for index in screen["facet_orbit_indices"])
    mirror_records: list[dict[str, Any]] = []
    mirror_types: Counter[tuple[int, int]] = Counter()
    presentation_types: Counter[tuple[str, int, int, str]] = Counter()
    reflection_matrices: list[list[list[str]]] = []
    for facet_position, orbit_index in enumerate(facet_orbits):
        dual_normal = sum(
            (walls[index] for index in wall_orbits[orbit_index]),
            sp.zeros(1, 20),
        )
        pairing_column = reduced_basis * dual_normal.T
        fixed_coordinates = fixed_gram.inv() * pairing_column
        primitive, scale = _primitive(fixed_coordinates)
        actual_primal_coordinates = dual_normal * ambient_dual_gram
        if fixed_coordinates.T * reduced_basis != actual_primal_coordinates:
            raise ArithmeticError("a facet normal was not recovered in the fixed span")
        square = int((primitive.T * fixed_gram * primitive)[0])
        divisibility = math.gcd(
            *(abs(int(value)) for value in fixed_gram * primitive)
        )
        reflection = (
            sp.eye(15)
            - sp.Rational(2, square)
            * (fixed_gram * primitive)
            * primitive.T
        )
        if not (
            all(value.is_Integer for value in reflection)
            and reflection**2 == sp.eye(15)
            and reflection * fixed_gram * reflection.T == fixed_gram
            and reflection.det() == -1
        ):
            raise ArithmeticError("a facet reflection is not an integral isometry")
        presentation_square = sp.Rational(
            (dual_normal * ambient_dual_gram * dual_normal.T)[0]
        )
        if presentation_square != scale**2 * square:
            raise ArithmeticError("primitive and presented mirror squares disagree")
        mirror_types[(square, divisibility)] += 1
        presentation_types[
            (_rational_text(presentation_square), square, divisibility, _rational_text(scale))
        ] += 1
        reflection_matrices.append(_matrix_payload(reflection))
        mirror_records.append(
            {
                "facet_index": facet_position + 1,
                "wall_orbit_index": orbit_index + 1,
                "source_wall_indices": [
                    index + 1 for index in wall_orbits[orbit_index]
                ],
                "presented_dual_normal": [int(value) for value in dual_normal],
                "presented_square": _rational_text(presentation_square),
                "primitive_fixed_coordinates": [int(value) for value in primitive],
                "presented_equals_scale_times_primitive": _rational_text(scale),
                "primitive_square": square,
                "primitive_divisibility": divisibility,
                "reflection_integral": True,
            }
        )

    expected_mirror_types = Counter({(-2, 2): 35, (-2, 1): 15, (-4, 2): 10})
    if mirror_types != expected_mirror_types:
        raise ArithmeticError("the primitive mirror-type histogram changed")

    return {
        "schema": "shimada_disc16_d5_h14_fixed_lattice_certificate_v1",
        "source": {
            "screen_certificate": str(screen_path.relative_to(ROOT)),
            "volume_certificate": str(volume_path.relative_to(ROOT)),
            "primary_disc16_sha256": screen["source"]["sha256"],
            "representative": "unique fixed-rank-15 element of the D5 stabilizer",
        },
        "saturated_fixed_lattice": {
            "definition": "Fix(M) = ker_Z(M^T-I) for Shimada row action M",
            "constraint_smith_diagonal": [int(value) for value in constraint_diagonal],
            "rank": 15,
            "saturated": True,
            "saturation_reason": (
                "zero-column basis in a unimodular Smith right transformation; "
                "equivalently, an integral kernel is primitive"
            ),
            "smith_kernel_basis_rows_in_NS": _matrix_payload(smith_kernel_basis),
            "reduced_basis_rows_in_NS": _matrix_payload(reduced_basis),
            "reduced_to_smith_basis_change": _matrix_payload(basis_change),
            "basis_change_determinant": int(basis_change.det()),
            "gram": _matrix_payload(fixed_gram),
            "inertia": list(exact_inertia(fixed_gram)),
            "determinant": int(fixed_gram.det()),
            "even": True,
            "smith_diagonal": list(fixed_snf_diagonal),
            "discriminant_group": "(Z/2Z)^7",
            "two_elementary_invariants": {"r": 15, "a": 7, "delta": delta},
            "standard_isometry_type": "U(2) + E8(-1) + A1(-1)^5",
            "identification_reason": (
                "standard uniqueness classification of indefinite even "
                "2-elementary lattices with signature (1,14) and invariants (15,7,1)"
            ),
        },
        "discriminant_form": {
            "generator_model": (
                "each bit vector y represents the dual class y/2 in the "
                "reduced fixed-lattice basis"
            ),
            "generator_bit_vectors": [list(vector) for vector in discriminant_bits],
            "quadratic_bilinear_matrix_mod_Z_with_diagonal_mod_2Z": discriminant_form_mod,
            "all_generators_have_order_two": True,
            "delta": delta,
        },
        "mirrors": {
            "facet_count": len(mirror_records),
            "primitive_type_histogram": {
                f"square_{square}_div_{divisibility}": count
                for (square, divisibility), count in sorted(mirror_types.items())
            },
            "presentation_type_histogram": {
                (
                    f"presented_square_{presented}__primitive_square_{square}"
                    f"_div_{divisibility}__scale_{scale}"
                ): count
                for (presented, square, divisibility, scale), count in sorted(
                    presentation_types.items()
                )
            },
            "all_reflections_integral": True,
            "all_reflections_preserve_fixed_gram": True,
            "records": mirror_records,
            "reflection_matrices_in_reduced_basis": reflection_matrices,
        },
        "reflectivity": {
            "integral_reflection_group": True,
            "finite_covolume": True,
            "fixed_lattice_reflective": True,
            "reason": (
                "the 60 integral mirrors bound the independently certified "
                "finite-volume H14 Coxeter chamber"
            ),
            "does_not_establish_literature_novelty": True,
        },
        "checks": {
            "fixed_kernel_saturated": True,
            "reduced_basis_is_unimodularly_equivalent_to_smith_basis": True,
            "fixed_gram_exact": True,
            "discriminant_form_exact": True,
            "all_60_mirrors_primitive_normalized": True,
            "all_60_reflections_integral_involutive_isometries": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen", type=Path, default=SCREEN_CERTIFICATE)
    parser.add_argument("--volume", type=Path, default=VOLUME_CERTIFICATE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate(arguments.screen, arguments.volume)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "rank": result["saturated_fixed_lattice"]["rank"],
                "determinant": result["saturated_fixed_lattice"]["determinant"],
                "smith_diagonal": result["saturated_fixed_lattice"]["smith_diagonal"],
                "two_elementary_invariants": result["saturated_fixed_lattice"]["two_elementary_invariants"],
                "standard_isometry_type": result["saturated_fixed_lattice"]["standard_isometry_type"],
                "primitive_mirror_types": result["mirrors"]["primitive_type_histogram"],
                "all_reflections_integral": result["mirrors"]["all_reflections_integral"],
                "reflective": result["reflectivity"]["fixed_lattice_reflective"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
