#!/usr/bin/env python3
"""Certify the integral lattice behind the 172-facet Jacobian--Kummer H12 fold.

The input is the exact trace-nine involution and its 172 irredundant wall
orbits archived in ``jacobian_kummer_h12_interim_facet_gram.json``.  This
script independently returns to Shimada's published ``GramS16`` and wall
vectors and certifies:

* a saturated rank-13 basis of the integral fixed lattice;
* signature (1,12), determinant 256, and Smith invariants 1^8,2^2,4^3;
* the discriminant form ``u(2) + <3/4> + <3/4> + <7/4>`` by five explicit
  dual-lattice representatives, without depending on a noncanonical Smith
  basis;
* the identification ``U(2) + D5(-1)^2 + <-4>`` via the standard indefinite
  even-lattice uniqueness criterion; and
* integrality of all 172 facet reflections and their primitive
  (square, divisibility) histogram.

This is a lattice/reflection certificate.  The logically separate finite-
volume proof is archived by ``verify_jacobian_kummer_h12_vinberg.py`` and
``independent_certify_jacobian_kummer_h12.py``.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from itertools import product
import json
from math import gcd
from pathlib import Path
import re
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "data"
    / "external"
    / "dolgachev_shimada_15_nodal_2019"
    / "15nodalcompdata.txt"
)
CHAMBER = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_interim_facet_gram.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_fixed_lattice_certificate.json"
)


def _extract_assignment(text: str, name: str) -> str:
    """Extract one top-level bracketed GAP assignment."""

    marker = f"{name}:="
    start = text.index(marker) + len(marker)
    while text[start].isspace():
        start += 1
    if text[start] != "[":
        raise ValueError(f"assignment {name} is not a bracketed list")
    depth = 0
    quote: str | None = None
    for end in range(start, len(text)):
        character = text[end]
        if quote is not None:
            if character == quote and text[end - 1] != "\\":
                quote = None
            continue
        if character in {'"', "'"}:
            quote = character
        elif character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
            if depth == 0:
                return text[start : end + 1]
    raise ValueError(f"unterminated assignment {name}")


def _parse_published_data() -> tuple[sp.Matrix, tuple[sp.Matrix, ...]]:
    text = SOURCE.read_text(encoding="utf-8").replace("\\\n", "")
    gram = sp.Matrix(ast.literal_eval(_extract_assignment(text, "GramS16")))
    walls_text = _extract_assignment(text, "D16WallRecs")
    wall_vectors = tuple(
        sp.Matrix(
            1,
            17,
            [sp.Rational(entry.strip()) for entry in vector.split(",")],
        )
        for vector in re.findall(r"\bvector\s*:=\s*\[([^\]]+)\]", walls_text)
    )
    if gram.shape != (17, 17) or len(wall_vectors) != 316:
        raise ArithmeticError("published S16/D16 data have unexpected dimensions")
    return gram, wall_vectors


def _integer_payload(matrix: sp.Matrix) -> list[list[int]]:
    if not all(sp.Rational(entry).q == 1 for entry in matrix):
        raise ArithmeticError("matrix is not integral")
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _rational_text(value: sp.Rational) -> str:
    value = sp.Rational(value)
    return str(value.p) if value.q == 1 else f"{value.p}/{value.q}"


def _rational_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [_rational_text(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _primitive_integer_column(vector: sp.Matrix) -> sp.Matrix:
    denominators = [sp.Rational(entry).q for entry in vector]
    denominator = int(sp.ilcm(*denominators))
    entries = [int(sp.Rational(entry) * denominator) for entry in vector]
    content = 0
    for entry in entries:
        content = gcd(content, abs(entry))
    if content == 0:
        raise ArithmeticError("zero vector has no primitive normalization")
    return sp.Matrix([entry // content for entry in entries])


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Exact rational congruence elimination."""

    work = sp.Matrix(matrix)
    positive = negative = 0
    while work.rows:
        pivot = next((i for i in range(work.rows) if work[i, i]), None)
        if pivot is None:
            pair = next(
                (
                    (i, j)
                    for i in range(work.rows)
                    for j in range(i + 1, work.rows)
                    if work[i, j]
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            i, j = pair
            order = [i, j] + [k for k in range(work.rows) if k not in pair]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0
        if pivot:
            order = [pivot] + [k for k in range(work.rows) if k != pivot]
            work = work.extract(order, order)
        diagonal = sp.Rational(work[0, 0])
        positive += int(bool(diagonal > 0))
        negative += int(bool(diagonal < 0))
        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / diagonal
    return positive, negative, 0


def _fixed_basis(involution: sp.Matrix) -> sp.Matrix:
    rows: list[list[int]] = []
    for vector in (involution.T - sp.eye(17)).nullspace():
        primitive = _primitive_integer_column(vector)
        rows.append([int(entry) for entry in primitive])
    basis = sp.Matrix(rows)
    if basis.shape != (13, 17) or basis * involution != basis:
        raise ArithmeticError("failed to reconstruct the rank-13 fixed space")
    return basis


def _target_gram() -> sp.Matrix:
    u2 = sp.Matrix([[0, 2], [2, 0]])
    d5 = sp.diag(*([-2] * 5))
    for left, right in ((0, 1), (1, 2), (2, 3), (2, 4)):
        d5[left, right] = d5[right, left] = 1
    return sp.diag(u2, d5, d5, sp.Matrix([[-4]]))


def _fractional_part(value: sp.Rational) -> sp.Rational:
    value = sp.Rational(value)
    return value - sp.floor(value)


def _reduced_form(gram: sp.Matrix, representatives: sp.Matrix) -> sp.Matrix:
    form = representatives.T * gram * representatives
    reduced = sp.zeros(form.rows)
    for row in range(form.rows):
        for column in range(form.cols):
            modulus = 2 if row == column else 1
            value = sp.Rational(form[row, column])
            reduced[row, column] = value - modulus * sp.floor(value / modulus)
    return reduced


def _class_order(vector: sp.Matrix) -> int:
    for order in range(1, 17):
        if all(sp.Rational(order * entry).q == 1 for entry in vector):
            return order
    raise ArithmeticError("dual class order unexpectedly exceeds 16")


def _distinct_class_count(
    representatives: sp.Matrix, orders: Iterable[int]
) -> int:
    classes: set[tuple[sp.Rational, ...]] = set()
    for coefficients in product(*(range(order) for order in orders)):
        vector = representatives * sp.Matrix(coefficients)
        classes.add(tuple(_fractional_part(entry) for entry in vector))
    return len(classes)


def _fixed_discriminant_representatives() -> sp.Matrix:
    """Columns u,v,a,b,c in L*/L, reduced modulo the fixed basis."""

    return sp.Matrix(
        [
            [0, sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(3, 4), 0],
            [0, 0, sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(3, 4)],
            [0, 0, sp.Rational(1, 2), sp.Rational(1, 4), sp.Rational(1, 4)],
            [sp.Rational(1, 2), sp.Rational(1, 2), 0, sp.Rational(1, 4), sp.Rational(1, 4)],
            [0, 0, 0, sp.Rational(1, 2), 0],
            [0, 0, 0, sp.Rational(1, 2), 0],
            [sp.Rational(1, 2), sp.Rational(1, 2), 0, sp.Rational(3, 4), sp.Rational(3, 4)],
            [0, sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(3, 4)],
            [0, 0, sp.Rational(3, 4), sp.Rational(3, 4), 0],
            [sp.Rational(1, 2), 0, 0, sp.Rational(1, 2), sp.Rational(3, 4)],
            [sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2), 0, 0],
            [0, 0, sp.Rational(1, 4), sp.Rational(1, 2), 0],
            [sp.Rational(1, 2), 0, 0, sp.Rational(1, 2), sp.Rational(1, 2)],
        ]
    )


def _target_discriminant_representatives(target: sp.Matrix) -> sp.Matrix:
    representatives = sp.zeros(13, 5)
    representatives[0, 0] = sp.Rational(1, 2)
    representatives[1, 1] = sp.Rational(1, 2)
    d5_inverse = target[2:7, 2:7].inv()
    representatives[2:7, 2] = d5_inverse[:, 3]
    representatives[7:12, 3] = d5_inverse[:, 3]
    representatives[12, 4] = sp.Rational(1, 4)
    return representatives


def build_certificate() -> dict[str, Any]:
    chamber = json.loads(CHAMBER.read_text(encoding="utf-8"))
    gram, published_walls = _parse_published_data()
    involution = sp.Matrix(chamber["group_class"]["representative_matrix"])
    basis = _fixed_basis(involution)
    fixed_gram = basis * gram * basis.T

    # A unit maximal minor proves that the row lattice is primitive.  Since it
    # already spans the entire rational fixed space, it is the saturated
    # integral fixed lattice.
    unit_minor_columns = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 13, 14)
    unit_minor = int(basis[:, unit_minor_columns].det())
    if abs(unit_minor) != 1:
        raise ArithmeticError("advertised fixed-basis unit minor is not a unit")

    smith = smith_normal_form(fixed_gram, domain=sp.ZZ)
    smith_diagonal = [abs(int(smith[i, i])) for i in range(13)]
    target = _target_gram()
    target_smith = smith_normal_form(target, domain=sp.ZZ)
    target_smith_diagonal = [abs(int(target_smith[i, i])) for i in range(13)]

    fixed_dual = _fixed_discriminant_representatives()
    target_dual = _target_discriminant_representatives(target)
    expected_orders = [2, 2, 4, 4, 4]
    expected_form = sp.Matrix(
        [
            [0, sp.Rational(1, 2), 0, 0, 0],
            [sp.Rational(1, 2), 0, 0, 0, 0],
            [0, 0, sp.Rational(3, 4), 0, 0],
            [0, 0, 0, sp.Rational(3, 4), 0],
            [0, 0, 0, 0, sp.Rational(7, 4)],
        ]
    )

    fixed_orders = [_class_order(fixed_dual[:, i]) for i in range(5)]
    target_orders = [_class_order(target_dual[:, i]) for i in range(5)]
    fixed_class_count = _distinct_class_count(fixed_dual, fixed_orders)
    target_class_count = _distinct_class_count(target_dual, target_orders)
    fixed_form = _reduced_form(fixed_gram, fixed_dual)
    target_form = _reduced_form(target, target_dual)

    mirror_histogram: Counter[tuple[int, int]] = Counter()
    source_breakdown: Counter[tuple[str, int, int]] = Counter()
    reflection_checks: list[bool] = []
    for wall_orbit in chamber["facet_source_wall_orbits"]:
        orbit_sum = sum(
            (published_walls[int(index) - 1] for index in wall_orbit),
            sp.zeros(1, 17),
        )
        coordinates = fixed_gram.inv() * basis * gram * orbit_sum.T
        primitive = _primitive_integer_column(coordinates)
        square = int((primitive.T * fixed_gram * primitive)[0])
        divisibility = 0
        for entry in fixed_gram * primitive:
            divisibility = gcd(divisibility, abs(int(entry)))
        mirror_histogram[(square, divisibility)] += 1
        source_square = _rational_text((orbit_sum * gram * orbit_sum.T)[0])
        source_breakdown[(source_square, square, divisibility)] += 1

        reflection = (
            sp.eye(13)
            - sp.Rational(2, square)
            * (fixed_gram * primitive)
            * primitive.T
        )
        reflection_checks.append(
            all(sp.Rational(entry).q == 1 for entry in reflection)
            and reflection**2 == sp.eye(13)
            and reflection * fixed_gram * reflection.T == fixed_gram
            and primitive.T * reflection == -primitive.T
        )

    checks = {
        "representative_is_trace_nine_involution": (
            involution.trace() == 9 and involution**2 == sp.eye(17)
        ),
        "fixed_basis_is_pointwise_fixed": basis * involution == basis,
        "fixed_basis_is_saturated_by_unit_minor": abs(unit_minor) == 1,
        "fixed_lattice_has_signature_1_12": _exact_inertia(fixed_gram) == (1, 12, 0),
        "fixed_lattice_is_even": all(int(fixed_gram[i, i]) % 2 == 0 for i in range(13)),
        "fixed_lattice_has_determinant_256": int(fixed_gram.det()) == 256,
        "fixed_smith_diagonal_is_1_8_2_2_4_3": smith_diagonal
        == [1] * 8 + [2] * 2 + [4] * 3,
        "five_fixed_dual_representatives_lie_in_dual": all(
            sp.Rational(entry).q == 1 for entry in fixed_gram * fixed_dual
        ),
        "fixed_dual_orders_are_2_2_4_4_4": fixed_orders == expected_orders,
        "fixed_dual_representatives_generate_all_256_classes": fixed_class_count
        == 256,
        "fixed_discriminant_form_is_canonical_form": fixed_form == expected_form,
        "target_has_same_signature_determinant_and_smith_diagonal": (
            _exact_inertia(target) == (1, 12, 0)
            and int(target.det()) == 256
            and target_smith_diagonal == smith_diagonal
        ),
        "target_dual_representatives_generate_all_256_classes": (
            target_orders == expected_orders and target_class_count == 256
        ),
        "target_discriminant_form_is_same_canonical_form": target_form
        == expected_form,
        "all_172_primitive_mirror_reflections_are_integral_isometries": (
            len(reflection_checks) == 172 and all(reflection_checks)
        ),
        "primitive_mirror_histogram_matches": mirror_histogram
        == Counter({(-2, 1): 16, (-4, 2): 36, (-4, 4): 96, (-8, 4): 24}),
    }

    return {
        "schema": "jacobian_kummer_h12_fixed_lattice_certificate_v1",
        "source": {
            "surface": "general genus-2 Jacobian Kummer K3 surface Y16",
            "paper": "Dolgachev--Shimada, 15-nodal quartic surfaces. Part II",
            "paper_arxiv": "https://arxiv.org/abs/1908.05390",
            "paper_doi": "https://doi.org/10.1007/s12215-019-00464-7",
            "computational_data": str(SOURCE.relative_to(ROOT)),
            "chamber_certificate": str(CHAMBER.relative_to(ROOT)),
        },
        "fixed_lattice": {
            "basis_rows_in_published_S16_basis": _integer_payload(basis),
            "unit_maximal_minor_columns_zero_based": list(unit_minor_columns),
            "unit_maximal_minor_determinant": unit_minor,
            "gram": _integer_payload(fixed_gram),
            "rank": 13,
            "signature": [1, 12, 0],
            "determinant": int(fixed_gram.det()),
            "smith_diagonal": smith_diagonal,
            "discriminant_group": "(Z/2)^2 + (Z/4)^3",
            "discriminant_length": 5,
            "discriminant_generator_names": ["u", "v", "a", "b", "c"],
            "discriminant_generator_orders": fixed_orders,
            "dual_generator_columns_mod_Z": _rational_payload(fixed_dual),
            "reduced_discriminant_form_gram": _rational_payload(fixed_form),
            "canonical_discriminant_form": "u(2) + <3/4> + <3/4> + <7/4>",
        },
        "target_identification": {
            "lattice": "U(2) + D5(-1)^2 + <-4>",
            "gram": _integer_payload(target),
            "reduced_discriminant_form_gram": _rational_payload(target_form),
            "uniqueness_criterion": {
                "reference": "Nikulin, Integral symmetric bilinear forms, Math. USSR Izv. 14 (1980)",
                "doi": "https://doi.org/10.1070/IM1980v014n01ABEH001060",
                "indefinite_even": True,
                "rank": 13,
                "discriminant_length": 5,
                "rank_at_least_length_plus_two": True,
                "conclusion": "the two lattices are isometric",
            },
        },
        "facet_reflections": {
            "facet_count": len(reflection_checks),
            "primitive_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(mirror_histogram.items())
            },
            "source_orbit_sum_square_to_primitive_type": {
                f"{source_square} -> ({square},{divisibility})": count
                for (source_square, square, divisibility), count in sorted(
                    source_breakdown.items()
                )
            },
            "all_reflections_integral_involutive_isometries": all(reflection_checks),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archived fixed-lattice certificate is stale")
        print("archived Jacobian--Kummer H12 fixed-lattice certificate verified")
    else:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)


if __name__ == "__main__":
    main()
