#!/usr/bin/env python3
"""Certify the integral lattice of the Jacobian--Kummer H10 class_0007 fold.

The fixed space in the discovery certificate is rational.  This verifier
returns to Dolgachev--Shimada's integral Neron--Severi lattice, computes the
saturated integral kernel of the involution with a Smith decomposition, and
certifies the intrinsic lattice and all facet reflections.

The resulting lattice is

    U + E8(-2) + <-4>.

Its discriminant form is exhibited as the orthogonal sum of the plus-type
eight-dimensional even quadratic space over F_2 and <-1/4>.  Nikulin's
indefinite uniqueness criterion then identifies the two integral lattices.
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
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp


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
    / "jacobian_kummer_h10_class0007_screen.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h10_class0007_fixed_lattice_certificate.json"
)


def _extract_assignment(text: str, name: str) -> str:
    marker = f"{name}:="
    start = text.index(marker) + len(marker)
    while text[start].isspace():
        start += 1
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
    walls = tuple(
        sp.Matrix(
            1,
            17,
            [sp.Rational(entry.strip()) for entry in vector.split(",")],
        )
        for vector in re.findall(r"\bvector\s*:=\s*\[([^\]]+)\]", walls_text)
    )
    if gram.shape != (17, 17) or len(walls) != 316:
        raise ArithmeticError("unexpected published S16/D16 dimensions")
    return gram, walls


def _smith_decomposition(matrix: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    domain_matrix = DomainMatrix.from_list_sympy(
        *matrix.shape, matrix.tolist()
    ).convert_to(ZZ)
    diagonal, left, right = smith_normal_decomp(domain_matrix)
    return diagonal.to_Matrix(), left.to_Matrix(), right.to_Matrix()


def _saturated_fixed_basis(involution: sp.Matrix) -> sp.Matrix:
    equation = involution.T - sp.eye(17)
    diagonal, _, right = _smith_decomposition(equation)
    rank = equation.rank()
    basis = right[:, rank:].T
    if any(diagonal[i, i] == 0 for i in range(rank)):
        raise ArithmeticError("Smith rank disagrees with rational rank")
    if basis.shape != (11, 17) or basis * involution != basis:
        raise ArithmeticError("failed to obtain the saturated fixed kernel")
    return basis


def _e8_cartan() -> sp.Matrix:
    return sp.Matrix(
        [
            [2, 0, -1, 0, 0, 0, 0, 0],
            [0, 2, 0, -1, 0, 0, 0, 0],
            [-1, 0, 2, -1, 0, 0, 0, 0],
            [0, -1, -1, 2, -1, 0, 0, 0],
            [0, 0, 0, -1, 2, -1, 0, 0],
            [0, 0, 0, 0, -1, 2, -1, 0],
            [0, 0, 0, 0, 0, -1, 2, -1],
            [0, 0, 0, 0, 0, 0, -1, 2],
        ]
    )


def _target_gram() -> sp.Matrix:
    return sp.diag(
        sp.Matrix([[0, 1], [1, 0]]),
        -2 * _e8_cartan(),
        sp.Matrix([[-4]]),
    )


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
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
            order = list(pair) + [i for i in range(work.rows) if i not in pair]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0
        if pivot:
            order = [pivot] + [i for i in range(work.rows) if i != pivot]
            work = work.extract(order, order)
        value = sp.Rational(work[0, 0])
        positive += int(bool(value > 0))
        negative += int(bool(value < 0))
        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / value
    return positive, negative, 0


def _primitive_integer_column(vector: sp.Matrix) -> sp.Matrix:
    denominator = int(sp.ilcm(*[sp.Rational(entry).q for entry in vector]))
    entries = [int(sp.Rational(entry) * denominator) for entry in vector]
    content = 0
    for entry in entries:
        content = gcd(content, abs(entry))
    if content == 0:
        raise ArithmeticError("zero vector")
    return sp.Matrix([entry // content for entry in entries])


def _integer_payload(matrix: sp.Matrix) -> list[list[int]]:
    if not all(sp.Rational(entry).q == 1 for entry in matrix):
        raise ArithmeticError("matrix is not integral")
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _rational_payload(matrix: sp.Matrix) -> list[list[str]]:
    return [[str(sp.Rational(matrix[i, j])) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _discriminant_basis(gram: sp.Matrix) -> tuple[list[int], sp.Matrix, sp.Matrix]:
    diagonal, _, right = _smith_decomposition(gram)
    invariants = [abs(int(diagonal[i, i])) for i in range(gram.rows)]
    indices = [i for i, value in enumerate(invariants) if value > 1]
    orders = [invariants[i] for i in indices]
    representatives = sp.Matrix.hstack(
        *[right[:, i] / invariants[i] for i in indices]
    )
    if not all(sp.Rational(entry).q == 1 for entry in gram * representatives):
        raise ArithmeticError("advertised discriminant representatives are not dual")
    form = representatives.T * gram * representatives
    reduced = sp.zeros(form.rows)
    for i in range(form.rows):
        for j in range(form.cols):
            modulus = 2 if i == j else 1
            reduced[i, j] = form[i, j] - modulus * sp.floor(form[i, j] / modulus)
    return orders, representatives, reduced


def _finite_form_split(form: sp.Matrix, orders: list[int]) -> dict[str, Any]:
    elements = list(product(*(range(order) for order in orders)))

    def raw(left: tuple[int, ...], right: tuple[int, ...]) -> sp.Rational:
        return sp.Rational(
            sum(
                left[i] * form[i, j] * right[j]
                for i in range(len(orders))
                for j in range(len(orders))
            )
        )

    def quadratic(element: tuple[int, ...]) -> sp.Rational:
        value = raw(element, element)
        return value - 2 * sp.floor(value / 2)

    def pairing(left: tuple[int, ...], right: tuple[int, ...]) -> sp.Rational:
        value = raw(left, right)
        return value - sp.floor(value)

    def element_order(element: tuple[int, ...]) -> int:
        for candidate in (1, 2, 4):
            if all(candidate * entry % order == 0 for entry, order in zip(element, orders)):
                return candidate
        raise ArithmeticError("unexpected discriminant element order")

    # Select the lexicographically first <-1/4> generator whose orthogonal
    # complement is the even plus-type quadratic space over F_2.  Smith
    # generators are noncanonical, so the fixed and target lattices need not
    # use the same coefficient tuple.
    generator: tuple[int, ...] | None = None
    complement: list[tuple[int, ...]] = []
    value_histogram: Counter[str] = Counter()
    for candidate in elements:
        if element_order(candidate) != 4 or quadratic(candidate) != sp.Rational(7, 4):
            continue
        candidate_complement = [
            element for element in elements if pairing(element, candidate) == 0
        ]
        candidate_histogram = Counter(
            str(quadratic(element)) for element in candidate_complement
        )
        if (
            len(candidate_complement) == 256
            and all(element_order(element) <= 2 for element in candidate_complement)
            and candidate_histogram == Counter({"0": 136, "1": 120})
        ):
            generator = candidate
            complement = candidate_complement
            value_histogram = candidate_histogram
            break
    if generator is None:
        generator = tuple(0 for _ in orders)
    checks = {
        "generator_has_order_4": element_order(generator) == 4,
        "generator_has_form_minus_one_quarter": quadratic(generator) == sp.Rational(7, 4),
        "orthogonal_complement_has_order_256": len(complement) == 256,
        "orthogonal_complement_is_elementary_2": all(
            element_order(element) <= 2 for element in complement
        ),
        "orthogonal_complement_is_even_plus_type": value_histogram
        == Counter({"0": 136, "1": 120}),
    }
    return {
        "order_four_generator_coefficients": list(generator),
        "order_four_generator_form": "7/4",
        "orthogonal_complement_order": len(complement),
        "orthogonal_complement_quadratic_value_histogram": dict(value_histogram),
        "canonical_form": "q_plus(F2^8) + <-1/4>",
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def build_certificate() -> dict[str, Any]:
    chamber = json.loads(CHAMBER.read_text(encoding="utf-8"))
    gram, published_walls = _parse_published_data()
    involution = sp.Matrix(
        [[sp.Rational(entry) for entry in row] for row in chamber["group_class"]["representative_matrix"]]
    )
    basis = _saturated_fixed_basis(involution)
    fixed_gram = basis * gram * basis.T
    target = _target_gram()

    fixed_smith = [
        abs(int(value))
        for value in smith_normal_form(fixed_gram, domain=sp.ZZ).diagonal()
    ]
    target_smith = [
        abs(int(value))
        for value in smith_normal_form(target, domain=sp.ZZ).diagonal()
    ]
    fixed_orders, fixed_dual, fixed_form = _discriminant_basis(fixed_gram)
    target_orders, _, target_form = _discriminant_basis(target)
    fixed_split = _finite_form_split(fixed_form, fixed_orders)
    target_split = _finite_form_split(target_form, target_orders)

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
        source_breakdown[(str((orbit_sum * gram * orbit_sum.T)[0]), square, divisibility)] += 1
        reflection = (
            sp.eye(11)
            - sp.Rational(2, square) * (fixed_gram * primitive) * primitive.T
        )
        reflection_checks.append(
            all(sp.Rational(entry).q == 1 for entry in reflection)
            and reflection**2 == sp.eye(11)
            and reflection * fixed_gram * reflection.T == fixed_gram
            and primitive.T * reflection == -primitive.T
        )

    expected_smith = [1, 1] + [2] * 8 + [4]
    checks = {
        "representative_is_trace_five_involution": (
            involution.trace() == 5 and involution**2 == sp.eye(17)
        ),
        "fixed_basis_is_saturated_smith_kernel": basis * involution == basis,
        "fixed_signature_is_1_10": _exact_inertia(fixed_gram) == (1, 10, 0),
        "fixed_lattice_is_even": all(fixed_gram[i, i] % 2 == 0 for i in range(11)),
        "fixed_determinant_is_1024": int(fixed_gram.det()) == 1024,
        "fixed_smith_diagonal_matches": fixed_smith == expected_smith,
        "target_signature_determinant_smith_match": (
            _exact_inertia(target) == (1, 10, 0)
            and int(target.det()) == 1024
            and target_smith == expected_smith
        ),
        "fixed_discriminant_form_has_canonical_split": fixed_split["all_checks_pass"],
        "target_discriminant_form_has_same_canonical_split": target_split["all_checks_pass"],
        "all_84_mirror_reflections_are_integral_isometries": (
            len(reflection_checks) == 84 and all(reflection_checks)
        ),
        "primitive_mirror_histogram_matches": mirror_histogram
        == Counter({(-2, 1): 8, (-4, 2): 33, (-4, 4): 43}),
    }

    return {
        "schema": "jacobian_kummer_h10_class0007_fixed_lattice_certificate_v1",
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
            "gram": _integer_payload(fixed_gram),
            "rank": 11,
            "signature": [1, 10, 0],
            "determinant": int(fixed_gram.det()),
            "smith_diagonal": fixed_smith,
            "discriminant_group": "(Z/2)^8 + Z/4",
            "discriminant_generator_orders": fixed_orders,
            "dual_generator_columns_mod_Z": _rational_payload(fixed_dual),
            "reduced_discriminant_form_gram": _rational_payload(fixed_form),
            "discriminant_form_split": fixed_split,
        },
        "target_identification": {
            "lattice": "U + E8(-2) + <-4>",
            "gram": _integer_payload(target),
            "smith_diagonal": target_smith,
            "reduced_discriminant_form_gram": _rational_payload(target_form),
            "discriminant_form_split": target_split,
            "uniqueness_criterion": {
                "reference": "Nikulin, Integral symmetric bilinear forms, Theorem 1.14.2",
                "doi": "https://doi.org/10.1070/IM1980v014n01ABEH001060",
                "indefinite_even": True,
                "rank": 11,
                "discriminant_length": 9,
                "rank_at_least_length_plus_two": True,
                "conclusion": "the fixed and target lattices are isometric",
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
                for (source_square, square, divisibility), count in sorted(source_breakdown.items())
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
            raise SystemExit("archived H10 fixed-lattice certificate is stale")
        print("archived Jacobian--Kummer H10 fixed-lattice certificate verified")
    else:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)


if __name__ == "__main__":
    main()
