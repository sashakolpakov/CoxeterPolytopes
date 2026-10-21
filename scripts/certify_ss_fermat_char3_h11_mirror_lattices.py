#!/usr/bin/env python3
"""Certify H11 mirror types and the minimal stable lattice for V4-H11."""

from __future__ import annotations

from collections import Counter
import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

from scripts.certify_ss_fermat_char3_survivor_fixed_lattices import (
    _matrix_payload,
    _primitive_integer_column,
)
from scripts.screen_ss_fermat_quartic_char3_c2_survivors import wall_cast
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/enumerations"
FIXED = DATA / "ss_fermat_quartic_char3_survivor_fixed_lattices.json"
C2_CLASSES = DATA / "ss_fermat_quartic_char3_c2_element_class_screen.json"
C2_SCREEN = DATA / "ss_fermat_quartic_char3_c2_survivor_screen.json"
V4_SCREEN = DATA / "ss_fermat_quartic_char3_v4_fold_screen.json"
C2_VOLUME = DATA / "ss_fermat_quartic_char3_c2_h11_vinberg_certificate.json"
V4_VOLUME = DATA / "ss_fermat_quartic_char3_v4_h11_exact_cone_certificate.json"
OUTPUT = DATA / "ss_fermat_quartic_char3_h11_mirror_lattices.json"


def _fraction_text(value: sp.Rational) -> str:
    value = sp.Rational(value)
    return str(value.p) if value.q == 1 else f"{value.p}/{value.q}"


def _reflection(root: sp.Matrix, gram: sp.Matrix) -> sp.Matrix:
    square = sp.Rational((root.T * gram * root)[0])
    return sp.eye(gram.rows) - sp.Rational(2, 1) / square * root * (root.T * gram)


def _reflection_record(
    index: int,
    root: sp.Matrix,
    gram: sp.Matrix,
    *,
    coordinate_key: str = "fixed_basis_coordinates",
) -> tuple[dict[str, Any], sp.Matrix]:
    square = int((root.T * gram * root)[0])
    divisibility = math.gcd(*(abs(int(value)) for value in gram * root))
    reflection = _reflection(root, gram)
    denominator = int(sp.ilcm(*[sp.Rational(value).q for value in reflection]))
    rational_valid = (
        reflection.T * gram * reflection == gram
        and reflection * reflection == sp.eye(gram.rows)
        and reflection * root == -root
    )
    integral = denominator == 1
    if not rational_valid:
        raise ArithmeticError("a mirror reflection failed exact rational replay")
    return (
        {
            "restricted_ray_index": index,
            coordinate_key: [int(value) for value in root],
            "square": square,
            "divisibility": divisibility,
            "reflection_matrix_denominator": denominator,
            "reflection_integral": integral,
            "reflection_rational_involutive_isometry": rational_valid,
        },
        reflection,
    )


def _roots_in_fixed_basis(
    fixed_case: dict[str, Any],
    indexed_ambient_dual_normals: list[tuple[int, list[int] | tuple[int, ...]]],
) -> tuple[list[dict[str, Any]], list[sp.Matrix], list[sp.Matrix]]:
    gram = sp.Matrix(fixed_case["gram"])
    basis = sp.Matrix(fixed_case["basis_rows_in_published_NS_basis"])
    records: list[dict[str, Any]] = []
    roots: list[sp.Matrix] = []
    reflections: list[sp.Matrix] = []
    for index, ambient_normal in indexed_ambient_dual_normals:
        root = _primitive_integer_column(
            gram.inv() * basis * sp.Matrix(ambient_normal)
        )
        record, reflection = _reflection_record(index, root, gram)
        records.append(record)
        roots.append(root)
        reflections.append(reflection)
    if len({tuple(root) for root in roots}) != len(roots):
        raise ArithmeticError("two facets produced the same primitive mirror")
    return records, roots, reflections


def _histogram(records: list[dict[str, Any]]) -> dict[str, int]:
    counts = Counter(
        (
            int(record["square"]),
            int(record["divisibility"]),
            bool(record["reflection_integral"]),
        )
        for record in records
    )
    return {
        f"{square},{divisibility},{'integral' if integral else 'nonintegral'}": count
        for (square, divisibility, integral), count in sorted(counts.items())
    }


def _minimal_stable_overlattice(
    gram: sp.Matrix,
    roots: list[sp.Matrix],
    reflections: list[sp.Matrix],
) -> dict[str, Any]:
    identity = sp.eye(gram.rows)
    initially_nonintegral = [
        index
        for index, reflection in enumerate(reflections)
        if any(value.is_Integer is not True for value in reflection)
    ]
    generators = identity.row_join(
        sp.Matrix.hstack(*(reflections[index] for index in initially_nonintegral))
    )
    denominator = int(sp.ilcm(*[sp.Rational(value).q for value in generators]))
    integral_generators = sp.Matrix(
        generators.rows,
        generators.cols,
        [int(value * denominator) for value in generators],
    )
    # Columns of B generate L + sum_s s(L).  Every stable lattice containing
    # L must contain this one-step closure.
    basis = hermite_normal_form(integral_generators) / denominator
    basis_inverse = basis.inv()
    stable_matrices = [basis_inverse * reflection * basis for reflection in reflections]
    if any(
        value.is_Integer is not True
        for reflection in stable_matrices
        for value in reflection
    ):
        raise ArithmeticError("the one-step reflection closure is not stable")
    if any(value.is_Integer is not True for value in basis_inverse):
        raise ArithmeticError("the closure does not contain the fixed divisor lattice")

    stable_gram = basis.T * gram * basis
    if any(value.is_Integer is not True for value in stable_gram):
        raise ArithmeticError("minimal stable overmodule is not integral")
    smith = smith_normal_form(stable_gram, domain=sp.ZZ)
    smith_diagonal = [abs(int(smith[index, index])) for index in range(gram.rows)]
    if exact_inertia(stable_gram) != (1, 11, 0):
        raise ArithmeticError("stable lattice lost the H11 signature")

    stable_records: list[dict[str, Any]] = []
    for source_record, root in zip(
        [None] * len(roots), roots, strict=True
    ):
        del source_record
        stable_root = _primitive_integer_column(basis_inverse * root)
        record, _ = _reflection_record(
            0,
            stable_root,
            stable_gram,
            coordinate_key="stable_lattice_coordinates",
        )
        stable_records.append(record)
    stable_types = Counter(
        (record["square"], record["divisibility"])
        for record in stable_records
    )
    return {
        "construction": "L_stable = L_fixed + sum over facet reflections s(L_fixed)",
        "minimality_reason": (
            "every reflection-stable lattice containing L_fixed contains the displayed "
            "one-step closure; the closure itself is stable under every facet reflection"
        ),
        "initial_nonintegral_reflection_count": len(initially_nonintegral),
        "initial_nonintegral_facet_positions": [
            index + 1 for index in initially_nonintegral
        ],
        "basis_columns_in_fixed_lattice_rational_coordinates": [
            [_fraction_text(basis[row, column]) for column in range(basis.cols)]
            for row in range(basis.rows)
        ],
        "basis_denominator": int(
            sp.ilcm(*[sp.Rational(value).q for value in basis])
        ),
        "overlattice_index": int(sp.Rational(1, 1) / abs(basis.det())),
        "gram": _matrix_payload(stable_gram),
        "signature": [1, 11],
        "even": all(int(stable_gram[index, index]) % 2 == 0 for index in range(12)),
        "determinant_absolute": abs(int(stable_gram.det())),
        "smith_diagonal_absolute": smith_diagonal,
        "all_328_reflections_integral_on_stable_lattice": True,
        "primitive_square_divisibility_histogram": {
            f"{square},{divisibility}": count
            for (square, divisibility), count in sorted(stable_types.items())
        },
        "mirrors_in_stable_lattice": stable_records,
    }


def build_certificate() -> dict[str, Any]:
    fixed = json.loads(FIXED.read_text(encoding="utf-8"))
    by_id = {case["id"]: case for case in fixed["cases"]}

    cyclic = json.loads(C2_CLASSES.read_text(encoding="utf-8"))
    c2_record = next(
        record
        for record in cyclic["high_fixed_rank_casts"]
        if int(record["hyperbolic_dimension"]) == 11
    )
    c2_cast = wall_cast(np.asarray(c2_record["representative_matrix"], dtype=np.int64))
    c2_reduction = json.loads(C2_SCREEN.read_text(encoding="utf-8"))
    c2_screen = next(
        cast for cast in c2_reduction["casts"] if cast["hyperbolic_dimension"] == 11
    )
    c2_normals = [
        (index, c2_cast["rays"][index - 1])
        for index in c2_screen["facet_restricted_ray_indices"]
    ]
    c2_records, _, _ = _roots_in_fixed_basis(by_id["char3-C2-H11"], c2_normals)

    v4 = json.loads(V4_SCREEN.read_text(encoding="utf-8"))
    v4_cast = next(
        cast
        for cast in v4["casts"]
        if int(cast["subgroup"]["hyperbolic_dimension"]) == 11
    )
    v4_normals = [
        (
            int(facet["restricted_ray_index"]),
            facet["primitive_normal_covector"],
        )
        for facet in v4_cast["strict_facet_gate"]["certificates"]
    ]
    v4_records, v4_roots, v4_reflections = _roots_in_fixed_basis(
        by_id["char3-V4-H11"], v4_normals
    )
    stable = _minimal_stable_overlattice(
        sp.Matrix(by_id["char3-V4-H11"]["gram"]),
        v4_roots,
        v4_reflections,
    )

    c2_volume = json.loads(C2_VOLUME.read_text(encoding="utf-8"))
    v4_volume = json.loads(V4_VOLUME.read_text(encoding="utf-8"))
    checks = {
        "c2_has_206_distinct_primitive_mirrors": len(c2_records) == 206,
        "all_c2_reflections_integral_on_fixed_lattice": all(
            record["reflection_integral"] for record in c2_records
        ),
        "c2_type_histogram_expected": _histogram(c2_records)
        == {
            "-4,2,integral": 30,
            "-2,1,integral": 12,
            "-2,2,integral": 164,
        },
        "v4_has_328_distinct_primitive_mirrors": len(v4_records) == 328,
        "v4_fixed_lattice_has_exactly_36_nonintegral_reflections": sum(
            not record["reflection_integral"] for record in v4_records
        )
        == 36,
        "v4_nonintegral_type_is_exactly_minus24_divisibility6": _histogram(v4_records).get(
            "-24,6,nonintegral"
        )
        == 36,
        "minimal_stable_overlattice_has_index4": stable["overlattice_index"] == 4,
        "minimal_stable_overlattice_is_even_integral_det48": (
            stable["even"] and stable["determinant_absolute"] == 48
        ),
        "all_v4_reflections_integral_on_minimal_stable_overlattice": stable[
            "all_328_reflections_integral_on_stable_lattice"
        ],
        "both_h11_chambers_are_finite_covolume": (
            c2_volume["finite_volume"]
            and v4_volume["decision"]["finite_volume"]
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"H11 mirror-lattice checks failed: {checks}")
    return {
        "schema": "ss_fermat_quartic_char3_h11_mirror_lattices_v1",
        "c2_h11": {
            "fixed_lattice_id": "char3-C2-H11",
            "facet_count": len(c2_records),
            "primitive_square_divisibility_integrality_histogram": _histogram(
                c2_records
            ),
            "all_reflections_integral_on_fixed_lattice": True,
            "mirrors": c2_records,
        },
        "v4_h11": {
            "fixed_lattice_id": "char3-V4-H11",
            "facet_count": len(v4_records),
            "primitive_square_divisibility_integrality_histogram": _histogram(
                v4_records
            ),
            "integral_reflection_count_on_fixed_lattice": sum(
                record["reflection_integral"] for record in v4_records
            ),
            "nonintegral_reflection_count_on_fixed_lattice": sum(
                not record["reflection_integral"] for record in v4_records
            ),
            "mirrors_in_fixed_lattice": v4_records,
            "minimal_reflection_stable_overlattice": stable,
            "arithmeticity": {
                "field": "Q",
                "reason": (
                    "all facet reflections preserve the integral Lorentzian stable "
                    "overlattice and the Coxeter group has finite covolume"
                ),
                "arithmetic": True,
            },
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    result = build_certificate()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archive differs: {args.output}")
    else:
        args.output.write_text(rendered, encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
