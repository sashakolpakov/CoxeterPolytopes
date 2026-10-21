#!/usr/bin/env python3
"""Independent hostile replay of the Type-II downstairs anti-Gram compiler.

The production certificate uses a negative breadth-first switching gauge,
project congruence code, and PARI/GP short-vector enumeration.  This replay
instead parses Shimada's ``12B`` record directly, uses a positive depth-first
gauge, computes inertia from exact characteristic polynomials and Sturm root
counts, constructs the full half-sum image from all eight binary kernel
cosets in one HNF, and enumerates short vectors by exact rational LDL recursion.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
import hashlib
import itertools
import json
from math import isqrt, lcm
from pathlib import Path
import tarfile
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form


ROOT = Path(__file__).resolve().parents[1]
SHIMADA_DATA = (
    ROOT
    / "data"
    / "external"
    / "shimada_l10l26note"
    / "L10L26compdata2.txt"
)
PRODUCTION_SCRIPT = (
    ROOT / "scripts" / "certify_enriques_type_ii_downstairs_anti_gram_census.py"
)
PRODUCTION = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_ii_downstairs_anti_gram_census.json"
)
PRODUCTION_REPORT = (
    ROOT / "results" / "enriques_type_ii_downstairs_anti_gram_census.md"
)
PRODUCTION_TESTS = (
    ROOT / "tests" / "test_enriques_type_ii_downstairs_anti_gram_census.py"
)
BRANDHORST_SHIMADA = (
    ROOT
    / "data"
    / "external"
    / "brandhorst_shimada_aut_enriques"
    / "AutEnrVol.pdf"
)
BRANDHORST_SHIMADA_README = BRANDHORST_SHIMADA.with_name("README.md")
MARTIN_PDF = (
    ROOT / "data" / "external" / "martin_enriques_2019" / "Martin2019.pdf"
)
MARTIN_README = MARTIN_PDF.with_name("README.md")
MARTIN_SOURCE = MARTIN_PDF.with_name("1703.08419-source.tar")
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "enriques_type_ii_downstairs_anti_gram_census_hostile_audit.json"
)

EXPECTED_SHA256 = {
    "shimada": "218c8188e467bde29369cdeda3771d767c3f70105d8fb56def8d698b90e75a6c",
    "production_script": "0fa0fcb8b4dff6ca2f9ed0b3c667680695a758332631851b3cc3348ab1618436",
    "production_archive": "04d4ac3b122d3ef0cf0b9feeb90e03f06e0692aa6a982df9d52a89df090e72e0",
    "production_report": "936da9d7dd0266e5f98e587f60db76b15e4805aad80444d0fb2cd1ce0f60711a",
    "production_tests": "567e80e28f171fcfe31fb5359f31037f265261f25093d951d40f41c202c6b9c2",
    "brandhorst_shimada_pdf": "e11f06d376996bbb2aa3221c6e67ac82a9fa99e9e56096ebc35279cf206e610b",
    "brandhorst_shimada_readme": "90e12733841c6048547a1849c605716c61935506cfcee0d8611f0eb40be7ffea",
    "martin_pdf": "7831fbf14ed90f4a30eec68a989dce660724e6384bd1da3f108d373929b905a9",
    "martin_readme": "51cae8a572170fbe79bb33ac1a3a7a44ad202e4b0d3b6861ffb6d2fa59c65704",
    "martin_source": "77f30c710719d5f30106604710a644410f956861903a71b9a21792acd2588776",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _vector_hash(vectors: Iterable[Iterable[int]]) -> str:
    canonical = sorted(tuple(int(value) for value in vector) for vector in vectors)
    return hashlib.sha256(
        json.dumps(canonical, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _matrix_key(matrix: Iterable[Iterable[int]]) -> tuple[int, ...]:
    rows = [list(map(int, row)) for row in matrix]
    return tuple(
        rows[left][right]
        for left in range(len(rows))
        for right in range(left, len(rows))
    )


def _key_set_hash(keys: Iterable[tuple[int, ...]]) -> str:
    return hashlib.sha256(
        json.dumps(sorted(keys), separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _extract_balanced_list(record: str, field: str) -> list[Any]:
    marker = f"{field} :="
    field_position = record.find(marker)
    if field_position < 0:
        raise ArithmeticError(f"could not find Shimada field {field}")
    start = record.find("[", field_position + len(marker))
    if start < 0:
        raise ArithmeticError(f"could not find opening bracket for Shimada field {field}")
    depth = 0
    for position in range(start, len(record)):
        if record[position] == "[":
            depth += 1
        elif record[position] == "]":
            depth -= 1
            if depth == 0:
                return ast.literal_eval(record[start : position + 1])
    raise ArithmeticError(f"unterminated Shimada field {field}")


def _type_ii_shimada_record() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Parse only literal integer arrays in record 12B, without project parsers."""

    text = SHIMADA_DATA.read_text(encoding="utf-8")
    name_marker = 'name := "12B"'
    name_position = text.find(name_marker)
    if name_position < 0 or text.find(name_marker, name_position + 1) >= 0:
        raise ArithmeticError("Shimada input has no unique 12B record")
    record_start = text.rfind("rec(", 0, name_position)
    record_end = text.find("), rec(", name_position)
    if record_start < 0 or record_end < 0:
        raise ArithmeticError("could not bound Shimada 12B record")
    record = text[record_start:record_end]
    gram_l10 = sp.Matrix(_extract_balanced_list(record, "GramL10"))
    configuration = sp.Matrix(_extract_balanced_list(record, "configuration"))
    walls = sp.Matrix(_extract_balanced_list(record, "walls"))
    return gram_l10, configuration, walls


def _switch(matrix: list[list[int]], signs: tuple[int, ...]) -> list[list[int]]:
    return [
        [
            signs[left] * int(matrix[left][right]) * signs[right]
            for right in range(len(matrix))
        ]
        for left in range(len(matrix))
    ]


def _dfs_positive_canonical(
    matrix: list[list[int]],
) -> tuple[tuple[int, ...], list[list[int]], list[int]]:
    """Independent gauge: force sorted DFS spanning-tree edges positive."""

    size = len(matrix)
    adjacency = [
        [other for other in range(size) if other != vertex and matrix[vertex][other]]
        for vertex in range(size)
    ]
    signs: list[int | None] = [None] * size

    def visit(vertex: int) -> None:
        for other in sorted(adjacency[vertex]):
            if signs[other] is None:
                edge_sign = 1 if matrix[vertex][other] > 0 else -1
                signs[other] = int(signs[vertex]) * edge_sign
                visit(other)

    for root in range(size):
        if signs[root] is None:
            signs[root] = 1
            visit(root)
    integral_signs = [int(value) for value in signs]
    canonical = _switch(matrix, tuple(integral_signs))
    return _matrix_key(canonical), canonical, integral_signs


def _direct_switch_orbit(matrix: list[list[int]]) -> set[tuple[int, ...]]:
    size = len(matrix)
    return {
        _matrix_key(_switch(matrix, signs + (1,)))
        for signs in itertools.product((-1, 1), repeat=size - 1)
    }


def _sturm_inertia(matrix: list[list[int]]) -> tuple[int, int, int]:
    """Exact inertia from the characteristic polynomial and Sturm counts."""

    gram = sp.Matrix(matrix)
    variable = sp.Symbol("lambda")
    polynomial = gram.charpoly(variable).as_poly()
    zero = 0
    zero_factor = sp.Poly(variable, variable, domain=sp.ZZ)
    while polynomial.eval(0) == 0:
        polynomial = polynomial.exquo(zero_factor)
        zero += 1
    positive = 0
    negative = 0
    for factor, multiplicity in sp.factor_list(polynomial)[1]:
        factor = sp.Poly(factor, variable)
        negative += int(factor.count_roots(-sp.oo, 0)) * int(multiplicity)
        positive += int(factor.count_roots(0, sp.oo)) * int(multiplicity)
    if positive + negative + zero != gram.rows:
        raise ArithmeticError("Sturm inertia did not account for all eigenvalues")
    return positive, negative, zero


def _histogram(counter: Counter[tuple[int, int, int]]) -> dict[str, int]:
    return {
        f"p{positive}_n{negative}_z{zero}": count
        for (positive, negative, zero), count in sorted(counter.items())
    }


def _gf2_rank(rows: Iterable[Iterable[int]]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        value = sum((int(entry) & 1) << index for index, entry in enumerate(row))
        while value:
            pivot = value.bit_length() - 1
            if pivot in pivots:
                value ^= pivots[pivot]
            else:
                pivots[pivot] = value
                break
    return len(pivots)


def _binary_left_kernel(walls: list[list[int]]) -> list[tuple[int, ...]]:
    kernel = []
    for coefficients in itertools.product((0, 1), repeat=len(walls)):
        if all(
            sum(coefficients[row] * walls[row][column] for row in range(len(walls)))
            % 2
            == 0
            for column in range(len(walls[0]))
        ):
            kernel.append(coefficients)
    return kernel


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=sp.ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index]
    ]


def _ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def _ldl_upper(matrix: sp.Matrix) -> tuple[list[sp.Rational], sp.Matrix]:
    """Exact ``matrix = U.T*diag(d)*U`` with unit upper-triangular U."""

    if matrix.rows != matrix.cols or matrix != matrix.T:
        raise ArithmeticError("LDL input is not symmetric square")
    upper = sp.eye(matrix.rows)
    diagonal: list[sp.Rational] = []
    for column in range(matrix.rows):
        pivot = sp.Rational(
            matrix[column, column]
            - sum(
                diagonal[index] * upper[index, column] ** 2
                for index in range(column)
            )
        )
        if pivot <= 0:
            raise ArithmeticError("LDL input is not positive definite")
        diagonal.append(pivot)
        for target in range(column + 1, matrix.rows):
            upper[column, target] = sp.Rational(
                matrix[column, target]
                - sum(
                    diagonal[index]
                    * upper[index, column]
                    * upper[index, target]
                    for index in range(column)
                ),
                pivot,
            )
    return diagonal, upper


def _exact_short_vectors(
    matrix: sp.Matrix, bound: int
) -> tuple[list[tuple[int, ...]], dict[str, Any]]:
    """Complete exact ellipsoid enumeration, independent of PARI qfminim."""

    diagonal, upper = _ldl_upper(matrix)
    vector = [0] * matrix.rows
    vectors: list[tuple[int, ...]] = []
    leaf_count = 0

    def recurse(index: int, remaining: sp.Rational) -> None:
        nonlocal leaf_count
        if index < 0:
            leaf_count += 1
            if any(vector):
                column = sp.Matrix(vector)
                exact = int((column.T * matrix * column)[0])
                if not 0 < exact <= bound:
                    raise ArithmeticError("LDL recursion admitted an out-of-bound vector")
                vectors.append(tuple(vector))
            return
        center = sp.Rational(
            sum(
                upper[index, later] * vector[later]
                for later in range(index + 1, matrix.rows)
            )
        )
        allowance = sp.Rational(remaining / diagonal[index])
        if allowance < 0:
            return
        center_numerator = int(center.p)
        center_denominator = int(center.q)
        radius_numerator = int(allowance.p) * center_denominator**2
        radius_denominator = int(allowance.q)
        radius = isqrt(radius_numerator // radius_denominator)
        low = _ceil_div(-radius - center_numerator, center_denominator)
        high = (radius - center_numerator) // center_denominator
        for value in range(low, high + 1):
            term = diagonal[index] * (sp.Rational(value) + center) ** 2
            if term <= remaining:
                vector[index] = value
                recurse(index - 1, remaining - term)
        vector[index] = 0

    recurse(matrix.rows - 1, sp.Rational(bound))
    vectors = sorted(set(vectors))
    if leaf_count != len(vectors) + 1:
        raise ArithmeticError("LDL recursion produced duplicate leaves")
    return vectors, {
        "engine": "independent exact rational LDL ellipsoid recursion",
        "bound": bound,
        "nonzero_vector_count": len(vectors),
        "including_zero_leaf_count": leaf_count,
        "vector_set_sha256": _vector_hash(vectors),
    }


def _canonical_mod_sign(vector: Iterable[int]) -> tuple[int, ...]:
    values = tuple(int(value) for value in vector)
    first = next(value for value in values if value)
    return values if first > 0 else tuple(-value for value in values)


def _dynkin_arms(cartan: sp.Matrix) -> list[int]:
    if any(cartan[index, index] != 2 for index in range(cartan.rows)):
        raise ArithmeticError("simple-root Gram has a non-two diagonal")
    adjacency = [[] for _ in range(cartan.rows)]
    for left in range(cartan.rows):
        for right in range(left + 1, cartan.rows):
            value = int(cartan[left, right])
            if value == -1:
                adjacency[left].append(right)
                adjacency[right].append(left)
            elif value != 0:
                raise ArithmeticError("simple-root Gram is not simply laced")
    branch = [index for index, neighbours in enumerate(adjacency) if len(neighbours) == 3]
    if len(branch) != 1 or sum(map(len, adjacency)) // 2 != cartan.rows - 1:
        raise ArithmeticError("simple-root graph is not a three-arm tree")
    origin = branch[0]
    arms = []
    for neighbour in adjacency[origin]:
        previous, current, length = origin, neighbour, 1
        while len(adjacency[current]) == 2:
            following = next(vertex for vertex in adjacency[current] if vertex != previous)
            previous, current = current, following
            length += 1
        if len(adjacency[current]) != 1:
            raise ArithmeticError("simple-root arm did not terminate at a leaf")
        arms.append(length)
    return sorted(arms)


def _classify_anti_lattice(gram: sp.Matrix) -> dict[str, Any]:
    positive = -gram
    roots_two, enumeration_two = _exact_short_vectors(positive, 2)
    if roots_two:
        if len(roots_two) != 18:
            raise ArithmeticError("non-rootless residue did not have eighteen roots")
        scaled = positive / 2
        if any(value.q != 1 for value in scaled):
            raise ArithmeticError("cubic scaled form is not integral")
        representatives = sorted({_canonical_mod_sign(vector) for vector in roots_two})
        root_basis = sp.Matrix(representatives)
        divided_gram = root_basis * scaled * root_basis.T
        if (
            len(representatives) != 9
            or abs(int(root_basis.det())) != 1
            or divided_gram != sp.eye(9)
        ):
            raise ArithmeticError("eighteen roots do not identify the cubic lattice")
        return {
            "isometry_class": "<-2>^9",
            "rootless_through_norm_two": False,
            "signed_minus_two_vector_count": len(roots_two),
            "minus_two_vectors_modulo_sign_count": len(representatives),
            "minus_two_vector_representatives": [list(vector) for vector in representatives],
            "unimodular_root_basis_determinant_absolute": abs(int(root_basis.det())),
            "divided_root_basis_gram": _rows(divided_gram),
            "minimum_enumeration": enumeration_two,
        }

    roots_four, enumeration_four = _exact_short_vectors(positive, 4)
    if len(roots_four) != 144 or any(
        int((sp.Matrix(vector).T * positive * sp.Matrix(vector))[0]) != 4
        for vector in roots_four
    ):
        raise ArithmeticError("rootless residue did not have 144 norm-four vectors")
    scaled = positive / 2
    if any(value.q != 1 for value in scaled):
        raise ArithmeticError("D9 scaled form is not integral")
    max_coordinate = max(abs(value) for vector in roots_four for value in vector)
    weights = [1]
    for _ in range(1, gram.rows):
        weights.append(max_coordinate * sum(weights) + 1)
    if any(sum(value * weight for value, weight in zip(root, weights)) == 0 for root in roots_four):
        raise ArithmeticError("hostile positive-system functional hit a root")
    positive_roots = {
        root
        for root in roots_four
        if sum(value * weight for value, weight in zip(root, weights)) > 0
    }
    simple_roots = sorted(
        root
        for root in positive_roots
        if not any(
            tuple(root[index] - summand[index] for index in range(gram.rows))
            in positive_roots
            for summand in positive_roots
        )
    )
    simple_basis = sp.Matrix(simple_roots)
    cartan = simple_basis * scaled * simple_basis.T
    arms = _dynkin_arms(cartan)
    if (
        len(simple_roots) != 9
        or abs(int(simple_basis.det())) != 1
        or arms != [1, 1, 6]
    ):
        raise ArithmeticError("norm-four roots do not identify D9")
    return {
        "isometry_class": "D9(-2)",
        "rootless_through_norm_two": True,
        "signed_minus_two_vector_count": 0,
        "signed_norm_four_vector_count": len(roots_four),
        "positive_system_functional_weights": weights,
        "simple_root_count": len(simple_roots),
        "simple_root_basis_rows": [list(vector) for vector in simple_roots],
        "simple_root_basis_determinant_absolute": abs(int(simple_basis.det())),
        "divided_simple_root_gram": _rows(cartan),
        "dynkin_arm_lengths": arms,
        "minimum_enumeration": enumeration_two,
        "norm_four_enumeration": enumeration_four,
    }


def _full_half_sum_lattice(
    gram_rows: list[list[int]], kernel: list[tuple[int, ...]]
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Compute the image of every congruence coset in one exact HNF."""

    gram = sp.Matrix(gram_rows)
    if gram.rank() != 9:
        raise ArithmeticError("full-half-sum construction expected rank nine")
    basis_indices = next(
        indices
        for indices in itertools.combinations(tuple(reversed(range(12))), 9)
        if gram.extract(indices, indices).det()
    )
    delta_basis_gram = gram.extract(basis_indices, basis_indices)
    coordinates = gram[:, list(basis_indices)] * delta_basis_gram.inv()
    if any(value.q != 1 for value in coordinates):
        raise ArithmeticError("reverse-lex delta basis did not give integral coordinates")
    if coordinates * delta_basis_gram * coordinates.T != gram:
        raise ArithmeticError("rank-nine factorization did not reconstruct anti-Gram")

    # 2*e_i lies in the congruence lattice and maps to delta_i.  Every one
    # of the eight kernel classes maps to k*C/2.  Using all classes avoids
    # selecting the production kernel basis and avoids its first HNF.
    generators = [list(coordinates.row(index)) for index in range(12)]
    coset_images = []
    fraction_classes = set()
    for coefficients in kernel:
        image = sp.Matrix([coefficients]) * coordinates / 2
        coset_images.append(image)
        twice = 2 * image
        if any(value.q != 1 for value in twice):
            raise ArithmeticError("binary half-sum image has non-half-integral coordinate")
        fraction_classes.add(tuple(int(value) % 2 for value in twice))
        if any(coefficients):
            generators.append(list(image.row(0)))
    denominator = lcm(*(int(value.q) for row in generators for value in row))
    integer_generators = sp.Matrix(
        [[int(value * denominator) for value in row] for row in generators]
    )
    image_hnf = hermite_normal_form(integer_generators.T).T
    if image_hnf.shape != (9, 9):
        raise ArithmeticError("one-stage image HNF did not have rank nine")
    image_basis = image_hnf / denominator
    inverse_basis = image_basis.inv()
    if any(value.q != 1 for value in sp.Matrix(generators) * inverse_basis):
        raise ArithmeticError("image HNF failed to contain a direct generator")
    expected_covolume = sp.Rational(1, len(fraction_classes))
    if abs(image_basis.det()) != expected_covolume:
        raise ArithmeticError("HNF covolume disagrees with direct fractional cosets")

    anti_gram = image_basis * delta_basis_gram * image_basis.T
    if anti_gram != anti_gram.T or any(value.q != 1 for value in anti_gram):
        raise ArithmeticError("full half-sum image does not have integral symmetric Gram")
    if any(int(anti_gram[index, index]) % 2 for index in range(9)):
        raise ArithmeticError("full half-sum image is not even")
    classification = _classify_anti_lattice(anti_gram)
    result = {
        "formula": "A={sum z_i*delta_i/2 : z in Z^12 and zW=0 mod 2}",
        "construction": "all eight binary kernel cosets mapped directly, then one exact HNF",
        "different_from_production_two_stage_hnf": True,
        "reverse_lex_delta_basis_indices_one_based": [index + 1 for index in basis_indices],
        "direct_generator_count": len(generators),
        "common_denominator_before_hnf": denominator,
        "fractional_coset_count_mod_Z9": len(fraction_classes),
        "basis_covolume_in_reverse_lex_delta_span": str(expected_covolume),
        "gram": _rows(anti_gram),
        "rank": anti_gram.rank(),
        "absolute_determinant": abs(int(anti_gram.det())),
        "smith_diagonal": _smith_diagonal(anti_gram),
        "even": True,
        "classification": classification,
    }
    internal = {
        "delta_coordinates": coordinates,
        "image_basis": image_basis,
        "anti_gram": anti_gram,
        "coset_images": coset_images,
    }
    return result, internal


def _binary_witnesses(
    gram: list[list[int]], kernel: list[tuple[int, ...]]
) -> list[list[int]]:
    witnesses = []
    for vector in kernel:
        numerator = sum(
            vector[left] * gram[left][right] * vector[right]
            for left in range(12)
            for right in range(12)
        )
        if numerator == -8:
            witnesses.append(list(vector))
    return witnesses


def _ternary_witnesses(
    gram: list[list[int]],
    walls: list[list[int]],
    kernel: list[tuple[int, ...]],
) -> list[dict[str, Any]]:
    """Exhaust {-1,0,1}^12 using its eight possible support parities."""

    witnesses = []
    for support_vector in kernel:
        support = [index for index, value in enumerate(support_vector) if value]
        if not support:
            continue
        # Fix the first nonzero coefficient to +1 to quotient global sign.
        for remaining_signs in itertools.product((-1, 1), repeat=len(support) - 1):
            coefficients = [0] * 12
            coefficients[support[0]] = 1
            for index, sign in zip(support[1:], remaining_signs):
                coefficients[index] = sign
            numerator = sum(
                coefficients[left] * gram[left][right] * coefficients[right]
                for left in support
                for right in support
            )
            if numerator != -8:
                continue
            downstairs_sum = [
                sum(coefficients[row] * walls[row][column] for row in range(12))
                for column in range(10)
            ]
            if any(value % 2 for value in downstairs_sum):
                raise ArithmeticError("kernel-supported ternary witness failed parity")
            witnesses.append(
                {
                    "coefficients": coefficients,
                    "support_one_based": [index + 1 for index in support],
                    "downstairs_sum": downstairs_sum,
                    "four_times_square": numerator,
                    "square": -2,
                }
            )
    witnesses.sort(key=lambda item: (len(item["support_one_based"]), item["coefficients"]))
    return witnesses


def _coefficient_hash(witnesses: list[dict[str, Any]]) -> str:
    return hashlib.sha256(
        json.dumps(
            [item["coefficients"] for item in witnesses], separators=(",", ":")
        ).encode("utf-8")
    ).hexdigest()


def _transport_coefficients(
    coefficients: Iterable[int], signs: Iterable[int]
) -> tuple[int, ...]:
    transported = tuple(int(value) * int(sign) for value, sign in zip(coefficients, signs))
    return _canonical_mod_sign(transported)


def _provenance_checks() -> dict[str, Any]:
    brandhorst_readme = BRANDHORST_SHIMADA_README.read_text(encoding="utf-8")
    brandhorst_clauses = [
        "Definition 1.1 and Propositions 1.2--1.3",
        "Table 1.1 row 184",
        "Definition 1.4(ii)",
        "Section 7.4",
        "`(D9,D9)`-generic, Number 184",
        "`12B` chamber",
    ]
    if any(clause not in brandhorst_readme for clause in brandhorst_clauses):
        raise ArithmeticError("Brandhorst--Shimada pinned provenance note changed")
    with tarfile.open(MARTIN_SOURCE, mode="r:*") as archive:
        member = archive.getmember("EnriquesFiniteAut.tex")
        extracted = archive.extractfile(member)
        if extracted is None:
            raise ArithmeticError("Martin source archive lost its main TeX file")
        martin_tex = " ".join(extracted.read().decode("utf-8").split())
    martin_clauses = [
        r"y^2 + \beta (s^2+s) xy + \beta^2 (s^2+s)^2 y = x^3 + \beta (s^2+s) x^2",
        r"J(\sigma): s \mapsto -s-1",
        r"translation by the section $N^- = (0,0)$",
        r"fixed point free if and only if $\beta \neq -64$",
        r"Enriques surfaces of type $\II$ are parametrized by $\bbA^1 - \{0,-64\}$",
    ]
    if any(" ".join(clause.split()) not in martin_tex for clause in martin_clauses):
        raise ArithmeticError("Martin Type-II equation/existence provenance changed")
    production_report = " ".join(
        PRODUCTION_REPORT.read_text(encoding="utf-8").split()
    )
    if "No new abstract Coxeter polytope is claimed." not in production_report:
        raise ArithmeticError("production scope statement changed")
    return {
        "brandhorst_shimada": {
            "pdf_sha256_pinned": True,
            "readme_sha256_pinned": True,
            "hostile_local_text_extraction_checked": True,
            "checked_locations": [
                "Definition 1.1 and Propositions 1.2--1.3",
                "Table 1.1 row 184",
                "Definition 1.4(ii)",
                "Section 7.4",
            ],
            "checked_statement": (
                "M_D9 has anti complement D9(-2); a very general Type-II "
                "surface is (D9,D9)-generic No.184 with nef chamber 12B"
            ),
        },
        "martin": {
            "pdf_and_source_sha256_pinned": True,
            "source_tex_checked_in_memory": True,
            "checked_statement": (
                "the Type-II K3-cover equation, involution, and fixed-point-free "
                "range beta != 0,-64 give independent existence lineage"
            ),
            "no_lift_gram_imported": True,
        },
    }


def analyze() -> dict[str, Any]:
    input_paths = {
        "shimada": SHIMADA_DATA,
        "production_script": PRODUCTION_SCRIPT,
        "production_archive": PRODUCTION,
        "production_report": PRODUCTION_REPORT,
        "production_tests": PRODUCTION_TESTS,
        "brandhorst_shimada_pdf": BRANDHORST_SHIMADA,
        "brandhorst_shimada_readme": BRANDHORST_SHIMADA_README,
        "martin_pdf": MARTIN_PDF,
        "martin_readme": MARTIN_README,
        "martin_source": MARTIN_SOURCE,
    }
    actual_hashes = {name: _sha256(path) for name, path in input_paths.items()}
    if actual_hashes != EXPECTED_SHA256:
        raise AssertionError(f"a pinned Type-II artifact changed: {actual_hashes}")
    production = json.loads(PRODUCTION.read_text(encoding="utf-8"))
    if not production["all_checks_pass"]:
        raise ArithmeticError("production archive does not claim its own checks pass")
    provenance = _provenance_checks()

    gram_l10, configuration_matrix, walls_matrix = _type_ii_shimada_record()
    if (
        gram_l10.shape != (10, 10)
        or configuration_matrix.shape != (12, 12)
        or walls_matrix.shape != (12, 10)
        or walls_matrix * gram_l10 * walls_matrix.T != configuration_matrix
    ):
        raise ArithmeticError("independent Shimada 12B reconstruction failed")
    configuration = _rows(configuration_matrix)
    walls = _rows(walls_matrix)
    edges = [
        (left, right)
        for left in range(12)
        for right in range(left + 1, 12)
        if configuration[left][right] == 1
    ]
    if len(edges) != 15 or any(
        configuration[left][right] not in (0, 1)
        for left in range(12)
        for right in range(left + 1, 12)
    ):
        raise ArithmeticError("independent Type-II edge alphabet changed")
    support = {0}
    while True:
        enlarged = support | {
            right if left in support else left
            for left, right in edges
            if left in support or right in support
        }
        if enlarged == support:
            break
        support = enlarged
    if support != set(range(12)):
        raise ArithmeticError("Type-II support graph is not connected")

    multiplicities: Counter[tuple[int, ...]] = Counter()
    representatives: dict[tuple[int, ...], list[list[int]]] = {}
    for values in itertools.product((-2, 2), repeat=len(edges)):
        gram = [[-4 if left == right else 0 for right in range(12)] for left in range(12)]
        for (left, right), value in zip(edges, values):
            gram[left][right] = gram[right][left] = value
        canonical_key, canonical, _ = _dfs_positive_canonical(gram)
        multiplicities[canonical_key] += 1
        representatives[canonical_key] = canonical
    if sum(multiplicities.values()) != 32768 or len(representatives) != 16:
        raise ArithmeticError("independent DFS switching quotient changed")

    production_records = production["exact_inertia_census"]["orbit_records"]
    production_in_hostile_gauge: dict[tuple[int, ...], dict[str, Any]] = {}
    production_to_hostile_signs: dict[tuple[int, ...], list[int]] = {}
    for record in production_records:
        key, _, signs = _dfs_positive_canonical(record["canonical_gram"])
        if key in production_in_hostile_gauge:
            raise ArithmeticError("production has duplicate hostile switching classes")
        production_in_hostile_gauge[key] = record
        production_to_hostile_signs[key] = signs
    if set(production_in_hostile_gauge) != set(representatives):
        raise ArithmeticError("production and hostile switching quotients differ")

    weighted: Counter[tuple[int, int, int]] = Counter()
    unweighted: Counter[tuple[int, int, int]] = Counter()
    inertia_by_key: dict[tuple[int, ...], tuple[int, int, int]] = {}
    for key, representative in representatives.items():
        direct_orbit = _direct_switch_orbit(representative)
        if len(direct_orbit) != 2048 or multiplicities[key] != 2048:
            raise ArithmeticError("complete Type-II switch orbit is not size 2048")
        if any(
            _dfs_positive_canonical(_switch(representative, signs + (1,)))[0] != key
            for signs in itertools.product((-1, 1), repeat=11)
        ):
            raise ArithmeticError("positive DFS gauge is not switching invariant")
        signature = _sturm_inertia(representative)
        inertia_by_key[key] = signature
        weighted[signature] += multiplicities[key]
        unweighted[signature] += 1
        production_record = production_in_hostile_gauge[key]
        production_inertia = production_record["inertia"]
        if signature != (
            production_inertia["positive"],
            production_inertia["negative"],
            production_inertia["zero"],
        ) or production_record["raw_assignment_multiplicity"] != multiplicities[key]:
            raise ArithmeticError("hostile Sturm inertia differs from production")
    expected_weighted = {
        "p0_n9_z3": 4096,
        "p1_n10_z1": 12288,
        "p1_n9_z2": 2048,
        "p2_n10_z0": 14336,
    }
    expected_unweighted = {
        "p0_n9_z3": 2,
        "p1_n10_z1": 6,
        "p1_n9_z2": 1,
        "p2_n10_z0": 7,
    }
    if _histogram(weighted) != expected_weighted or _histogram(unweighted) != expected_unweighted:
        raise ArithmeticError("independent Type-II inertia histograms changed")

    kernel = _binary_left_kernel(walls)
    if len(kernel) != 8 or _gf2_rank(walls) != 9:
        raise ArithmeticError("independent wall congruence kernel changed")
    kernel_hash = _vector_hash(kernel)
    nonpositive_keys = [key for key, value in inertia_by_key.items() if value == (0, 9, 3)]
    if len(nonpositive_keys) != 2:
        raise ArithmeticError("hostile signature gate did not leave two rank-nine classes")
    hostile_lattices = {
        key: _full_half_sum_lattice(representatives[key], kernel)
        for key in nonpositive_keys
    }
    hostile_survivor_keys = [
        key
        for key, (lattice, _) in hostile_lattices.items()
        if lattice["classification"]["rootless_through_norm_two"]
    ]
    if len(hostile_survivor_keys) != 1:
        raise ArithmeticError("independent full lattices did not select one survivor")
    survivor_key = hostile_survivor_keys[0]
    killed_key = next(key for key in nonpositive_keys if key != survivor_key)
    survivor_lattice, survivor_internal = hostile_lattices[survivor_key]
    killed_lattice, killed_internal = hostile_lattices[killed_key]
    production_survivor = production_in_hostile_gauge[survivor_key][
        "generated_anti_lattice"
    ]
    production_killed = production_in_hostile_gauge[killed_key][
        "generated_anti_lattice"
    ]
    if not (
        survivor_lattice["classification"]["isometry_class"] == "D9(-2)"
        and survivor_lattice["absolute_determinant"] == 2048
        and survivor_lattice["smith_diagonal"] == [2] * 8 + [8]
        and killed_lattice["classification"]["isometry_class"] == "<-2>^9"
        and killed_lattice["absolute_determinant"] == 512
        and killed_lattice["smith_diagonal"] == [2] * 9
    ):
        raise ArithmeticError("independent full anti-lattice classification changed")
    for hostile, archived in (
        (survivor_lattice, production_survivor),
        (killed_lattice, production_killed),
    ):
        if (
            hostile["absolute_determinant"] != archived["absolute_determinant"]
            or hostile["smith_diagonal"] != archived["smith_diagonal"]
            or hostile["classification"]["isometry_class"]
            != archived["scaled_lattice_classification"]["isometry_class"]
        ):
            raise ArithmeticError("hostile anti-lattice invariant differs from production")

    hostile_survivor_binary = _binary_witnesses(
        representatives[survivor_key], kernel
    )
    hostile_killed_binary = _binary_witnesses(representatives[killed_key], kernel)
    production_survivor_binary = _binary_witnesses(
        production_in_hostile_gauge[survivor_key]["canonical_gram"], kernel
    )
    production_killed_binary = _binary_witnesses(
        production_in_hostile_gauge[killed_key]["canonical_gram"], kernel
    )
    if not (
        production_survivor_binary == []
        and production_killed_binary == []
        and hostile_survivor_binary == []
        and len(hostile_killed_binary) == 3
    ):
        raise ArithmeticError("binary half-sum gauge-dependence census changed")
    survivor_ternary = _ternary_witnesses(representatives[survivor_key], walls, kernel)
    killed_ternary = _ternary_witnesses(representatives[killed_key], walls, kernel)
    if survivor_ternary or len(killed_ternary) != 36 or {
        len(item["support_one_based"]) for item in killed_ternary
    } != {9}:
        raise ArithmeticError("hostile ternary witness census changed")

    # Recompute the production-gauge witness hash, then transport the complete
    # 36-element set into the hostile gauge rather than trusting three samples.
    production_killed_record = production_in_hostile_gauge[killed_key]
    production_gauge_witnesses = _ternary_witnesses(
        production_killed_record["canonical_gram"], walls, kernel
    )
    archived_witness_search = production_killed_record["ternary_half_sum_root_search"]
    if (
        _coefficient_hash(production_gauge_witnesses)
        != archived_witness_search["all_witness_coefficients_sha256"]
        or len(production_gauge_witnesses) != 36
    ):
        raise ArithmeticError("production-gauge ternary witness archive did not replay")
    killed_signs = production_to_hostile_signs[killed_key]
    transported_production_witnesses = {
        _transport_coefficients(item["coefficients"], killed_signs)
        for item in production_gauge_witnesses
    }
    hostile_witness_set = {
        tuple(item["coefficients"]) for item in killed_ternary
    }
    if transported_production_witnesses != hostile_witness_set:
        raise ArithmeticError("complete ternary witness set did not transport between gauges")

    # Map all 36 coefficient witnesses into the exact HNF lattice.  They must
    # cover exactly the nine cubic root pairs, four presentations per pair.
    root_fibres: Counter[tuple[int, ...]] = Counter()
    coordinate_map = killed_internal["delta_coordinates"] / 2
    image_basis_inverse = killed_internal["image_basis"].inv()
    for item in killed_ternary:
        lattice_coordinates = (
            sp.Matrix([item["coefficients"]])
            * coordinate_map
            * image_basis_inverse
        )
        if any(value.q != 1 for value in lattice_coordinates):
            raise ArithmeticError("ternary witness is not integral in hostile HNF basis")
        root_fibres[_canonical_mod_sign(lattice_coordinates)] += 1
    enumerated_root_pairs = {
        tuple(vector)
        for vector in killed_lattice["classification"][
            "minus_two_vector_representatives"
        ]
    }
    if set(root_fibres) != enumerated_root_pairs or set(root_fibres.values()) != {4}:
        raise ArithmeticError("36 coefficient witnesses do not cover 9 roots fourfold")

    # Exhaust every nontrivial sheet switch and every killed witness directly.
    sheet_switch_count = 0
    witness_transport_checks = 0
    killed_gram = representatives[killed_key]
    for initial in itertools.product((-1, 1), repeat=11):
        signs = initial + (1,)
        switched_gram = _switch(killed_gram, signs)
        sheet_switch_count += 1
        for item in killed_ternary:
            transported = [
                value * sign for value, sign in zip(item["coefficients"], signs)
            ]
            if any(
                sum(transported[row] * walls[row][column] for row in range(12)) % 2
                for column in range(10)
            ):
                raise ArithmeticError("sheet transport changed congruence parity")
            square_times_four = sum(
                transported[left]
                * switched_gram[left][right]
                * transported[right]
                for left in range(12)
                for right in range(12)
            )
            if square_times_four != -8:
                raise ArithmeticError("sheet transport changed a root norm")
            witness_transport_checks += 1
    if sheet_switch_count != 2048 or witness_transport_checks != 2048 * 36:
        raise ArithmeticError("sheet-invariance exhaustion denominator changed")

    checks = {
        "all_32768_anti_grams_generated_from_independently_parsed_12B": True,
        "positive_DFS_gauge_gives_16_complete_orbits_of_size_2048": True,
        "all_16_exact_Sturm_inertias_match_production": True,
        "exactly_two_negative_semidefinite_rank_9_orbits": True,
        "all_eight_mod_2_kernel_cosets_used_in_one_stage_image_HNF": True,
        "independent_LDL_enumeration_identifies_D9_minus2_and_cubic_lattices": True,
        "binary_shortcut_is_recorded_as_sheet_gauge_dependent": True,
        "all_36_ternary_coefficients_map_four_to_each_of_nine_cubic_roots": True,
        "all_2048_sheet_switches_transport_all_36_witnesses": True,
        "Brandhorst_Shimada_and_Martin_provenance_pins_checked": True,
        "no_new_surface_or_polytope_claim_preserved": True,
    }
    return {
        "schema": "enriques_type_ii_downstairs_anti_gram_census_hostile_audit_v1",
        "verdict": "PASS/FREEZE",
        "pinned_inputs": actual_hashes,
        "independent_downstairs_parse": {
            "record": "12B",
            "parser": "record-local balanced-bracket literal parser; no project GAP parser",
            "wall_count": 12,
            "m1_edge_count": len(edges),
            "support_connected": True,
            "wall_rank_over_Q": walls_matrix.rank(),
            "wall_rank_mod_2": _gf2_rank(walls),
            "binary_left_kernel_size": len(kernel),
            "binary_left_kernel_dimension": 3,
            "all_kernel_vectors_sha256": kernel_hash,
        },
        "independent_switching_census": {
            "canonical_gauge": "positive sorted depth-first spanning-tree edges",
            "different_from_production_negative_BFS_gauge": True,
            "raw_assignment_count": sum(multiplicities.values()),
            "orbit_count": len(representatives),
            "orbit_size_histogram": {"2048": 16},
            "cycle_rank": len(edges) - 12 + 1,
            "canonical_key_set_sha256": _key_set_hash(representatives),
            "all_direct_orbits_replayed": True,
        },
        "independent_exact_inertia": {
            "engine": "factorized exact characteristic polynomial plus Sturm root counts",
            "raw_assignment_histogram": _histogram(weighted),
            "switch_orbit_histogram": _histogram(unweighted),
            "negative_semidefinite_rank_9_orbit_count": len(nonpositive_keys),
        },
        "full_half_sum_lattices": {
            "method": (
                "factor through a reverse-lex rank-nine delta basis; map 2Z^12 "
                "and every one of the eight parity-kernel cosets; take one exact HNF"
            ),
            "survivor_hostile_key_sha256": _key_set_hash([survivor_key]),
            "survivor": survivor_lattice,
            "eliminated_hostile_key_sha256": _key_set_hash([killed_key]),
            "eliminated": killed_lattice,
        },
        "ternary_root_audit": {
            "enumeration_method": (
                "exhaust all signings of all eight allowed support-parity vectors; "
                "this is exactly {-1,0,1}^12 modulo global sign and zW=0 mod 2"
            ),
            "survivor_witness_count_modulo_sign": len(survivor_ternary),
            "eliminated_witness_count_modulo_sign": len(killed_ternary),
            "binary_shortcut_gauge_dependence": {
                "production_negative_BFS_survivor_count": len(
                    production_survivor_binary
                ),
                "production_negative_BFS_eliminated_count": len(
                    production_killed_binary
                ),
                "hostile_positive_DFS_survivor_count": len(hostile_survivor_binary),
                "hostile_positive_DFS_eliminated_count": len(hostile_killed_binary),
                "hostile_positive_DFS_eliminated_binary_witnesses": (
                    hostile_killed_binary
                ),
                "conclusion": (
                    "0/1 coefficient searches depend on sheet labels; only the full "
                    "congruence lattice and signed witness transport are invariant"
                ),
            },
            "eliminated_witness_coefficients_sha256_hostile_gauge": _coefficient_hash(
                killed_ternary
            ),
            "production_gauge_witness_hash_replayed": _coefficient_hash(
                production_gauge_witnesses
            ),
            "minimum_support_size": min(
                len(item["support_one_based"]) for item in killed_ternary
            ),
            "first_three_hostile_gauge_witnesses": killed_ternary[:3],
            "distinct_lattice_root_pairs": len(root_fibres),
            "coefficient_presentations_per_root_pair": {
                str(multiplicity): count
                for multiplicity, count in sorted(
                    Counter(root_fibres.values()).items()
                )
            },
            "complete_set_transports_from_production_gauge": True,
        },
        "sheet_invariance": {
            "sheet_switches_checked": sheet_switch_count,
            "witnesses_per_switch": len(killed_ternary),
            "direct_parity_and_norm_checks": witness_transport_checks,
            "coefficient_transport": "z -> zE for D -> EDE",
            "survivor_rootlessness_invariant_by_the_same_integral_lattice_isometry": True,
        },
        "geometric_provenance": provenance,
        "frozen_conclusion": (
            "The Type-II 12B two-sheet compiler is sound: its 32,768 anti-Grams "
            "have 16 switching classes; signature leaves two rank-nine classes; "
            "the full half-sum lattice is D9(-2) for exactly one and <-2>^9 for "
            "the other, whose 36 ternary presentations are genuine sheet-invariant roots."
        ),
        "scope": {
            "positive_result": (
                "a second exact surface-to-wall compiler and a forced realization "
                "of the known very general Type-II 12B carrier"
            ),
            "new_surface_claimed": False,
            "new_abstract_hyperbolic_Coxeter_polytope_claimed": False,
            "period_or_ample_cone_reconstructed_from_downstairs_alone": False,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(analyze(), indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived hostile audit differs: {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(hashlib.sha256(rendered.encode("utf-8")).hexdigest())


if __name__ == "__main__":
    main()
