#!/usr/bin/env python3
"""Certify a very-general full-Picard K3 carrier for GHESS-C2-H10/27.

The generic-Hessian source package already proves that the unique high-rank
transposition fold is the known finite-volume H10/27 chamber.  This
composition certificate adds the missing full-Picard period carrier:

* reconstruct the saturated determinant-256 fixed lattice and all 27 walls;
* identify its complete discriminant form;
* glue the Hessian Picard lattice explicitly into an even unimodular
  signature-(3,19) lattice and check that the fixed lattice is primitive;
* prove that the 27 primitive mirrors generate full Picard and that every
  reflection is integral; and
* rerun both literal Vinberg and exact GMP double-description volume gates.

The result is a period-level mixed divisor-reflection carrier.  It is not an
equation family, a full-nef statement, a surface-automorphism statement, or
a new abstract polytope.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Iterable

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_generic_hessian_transposition_h10_27 import (
    _fold_data,
    _picard_model,
    _vinberg,
)
from scripts.ingest_generic_hessian_quartic_84_stock import (
    FORMAL_LABELS,
    _basis_coordinates,
    _vector,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/generic_hessian_transposition_h10_27.json"
STOCK = ROOT / "data/enumerations/generic_hessian_quartic_84_wall_stock.json"
CPP = ROOT / "scripts/certify_borcherds_char2_high_fixed_sections_vinberg.cpp"
OUTPUT = ROOT / "data/enumerations/ghess_h10_27_full_picard_k3_carrier.json"

SCHEMA = "coxeter.ghess.h10_27.full_picard_complex_k3_carrier.v1"
EXPECTED_SHA256 = {
    "source": "c637c54385b33f1f116fff8b9d3cb833c08f48e22b78b8cd809b6d25a04532fb",
    "stock": "29ef2b4e33f804dfccb4f36bce2eb6a1cbb13f8ba2501c3612ff372a8397c3c6",
    "vinberg_cpp": "45174f5461d28affc1fdb4cb55047b7a1581037d307e32800651badb3c2df045",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _text(value: int | Fraction | sp.Rational) -> str:
    value = Fraction(value)
    return str(value.numerator) if value.denominator == 1 else str(value)


def _matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [
            int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
            for value in matrix.row(index)
        ]
        for index in range(matrix.rows)
    ]


def _vector_payload(vector: Iterable[Any]) -> list[int | str]:
    return [
        int(value) if sp.Rational(value).q == 1 else str(sp.Rational(value))
        for value in vector
    ]


def _smith(matrix: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(diagonal[index, index]))
        for index in range(min(diagonal.shape))
        if diagonal[index, index]
    ]


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(value)) for value in vector if value))


def _divisibility(gram: sp.Matrix, vector: sp.Matrix) -> int:
    return _content(gram * vector)


def _reflection(gram: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = sp.Rational((root.T * gram * root)[0])
    return sp.simplify(
        sp.eye(gram.rows)
        - sp.Rational(2, square) * root * (root.T * gram)
    )


def _mod_one(value: int | Fraction | sp.Rational) -> Fraction:
    value = Fraction(value)
    return value - value.numerator // value.denominator


def _mod_two(value: int | Fraction | sp.Rational) -> Fraction:
    value = Fraction(value)
    return value - 2 * (value // 2)


def _add_classes(left: tuple[Fraction, ...], right: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(_mod_one(a + b) for a, b in zip(left, right))


def _multiply_class(coefficient: int, value: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    result = (Fraction(0),) * len(value)
    for _ in range(coefficient):
        result = _add_classes(result, value)
    return result


def _class_order(value: tuple[Fraction, ...]) -> int:
    result = (Fraction(0),) * len(value)
    for order in range(1, 100):
        result = _add_classes(result, value)
        if not any(result):
            return order
    raise ArithmeticError("discriminant-class order search did not stop")


def _generated_classes(generators: list[tuple[Fraction, ...]]) -> set[tuple[Fraction, ...]]:
    result = {(Fraction(0),) * len(generators[0])}
    for generator in generators:
        result = {
            _add_classes(value, _multiply_class(index, generator))
            for value in result
            for index in range(_class_order(generator))
        }
    return result


def _finite_q(gram: sp.Matrix, value: tuple[Fraction, ...]) -> Fraction:
    vector = sp.Matrix(value)
    return _mod_two((vector.T * gram * vector)[0])


def _finite_b(
    gram: sp.Matrix,
    left: tuple[Fraction, ...],
    right: tuple[Fraction, ...],
) -> Fraction:
    return _mod_one((sp.Matrix(left).T * gram * sp.Matrix(right))[0])


def _u_split(gram: sp.Matrix) -> dict[str, Any]:
    """Give a determinant-one decomposition S=U+K_H."""

    isotropic_left = sp.Matrix([-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0])
    isotropic_right = sp.Matrix([-1, -1, 0, -1, -1, -1, -1, -1, -1, -1, 0])

    # Free coordinates x_2,...,x_10 obey x_2+x_3+x_4=0 mod 2.
    parity_basis = [
        (1, 1, 0) + (0,) * 6,
        (1, 0, 1) + (0,) * 6,
        (0, 1, 1) + (0,) * 6,
    ]
    for index in range(3, 9):
        row = [0] * 9
        row[index] = 1
        parity_basis.append(tuple(row))

    orthogonal_columns = []
    for free in parity_basis:
        vector = [0] * 11
        for index, value in zip(range(2, 11), free):
            vector[index] = value
        vector[0] = 2 * vector[2] - vector[10]
        vector[1] = (
            (vector[2] + vector[3] + vector[4]) // 2
            + vector[5]
            + vector[6]
            - vector[7]
            - vector[8]
            - vector[10]
        )
        orthogonal_columns.append(sp.Matrix(vector))

    change = isotropic_left.row_join(isotropic_right).row_join(
        sp.Matrix.hstack(*orthogonal_columns)
    )
    split = sp.simplify(change.T * gram * change)
    U = sp.Matrix([[0, 1], [1, 0]])
    negative = split[2:, 2:]
    checks = {
        "determinant_one_change": abs(int(change.det())) == 1,
        "literal_U_block": split[:2, :2] == U,
        "orthogonal_direct_sum": (
            split[:2, 2:] == sp.zeros(2, 9)
            and split[2:, :2] == sp.zeros(9, 2)
        ),
        "negative_definite_rank9": exact_inertia(negative) == (0, 9, 0),
        "negative_determinant_minus256": negative.det() == -256,
    }
    if not all(checks.values()):
        raise ArithmeticError([key for key, value in checks.items() if not value])
    return {
        "identification": "S_H=U+K_H, with K_H the displayed even negative-definite rank-nine lattice",
        "standard_to_source_fixed_basis": _matrix_payload(change),
        "K_H_gram": _matrix_payload(negative),
        "checks": checks,
    }


def _discriminant_form(gram: sp.Matrix) -> dict[str, Any]:
    """Identify q_S as <3/4> + u(2)^2 + v(2)."""

    half = Fraction(1, 2)
    quarter = Fraction(1, 4)
    cyclic = (0, half, 0, 0, 0, 0, 0, quarter, quarter, quarter, 0)
    elementary = [
        (0, 0, 0, 0, 0, 0, 0, 0, half, half, 0),
        (0, 0, 0, 0, 0, 0, 0, half, 0, half, 0),
        (0, 0, 0, 0, 0, 0, half, 0, 0, half, 0),
        (0, 0, 0, 0, 0, half, 0, 0, 0, half, 0),
        (0, 0, 0, half, half, 0, 0, 0, 0, 0, 0),
        (0, 0, half, 0, half, 0, 0, 0, 0, 0, 0),
    ]
    elementary = [tuple(Fraction(value) for value in row) for row in elementary]

    def total(*indices: int) -> tuple[Fraction, ...]:
        result = (Fraction(0),) * 11
        for index in indices:
            result = _add_classes(result, elementary[index])
        return result

    normal = [
        tuple(Fraction(value) for value in cyclic),
        total(3, 5),
        total(2, 4, 5),
        total(1, 3, 5),
        total(0, 2, 4, 5),
        total(5),
        total(4),
    ]
    orders = [_class_order(value) for value in normal]
    q_values = [_finite_q(gram, value) for value in normal]
    b_matrix = sp.Matrix(
        [[_finite_b(gram, left, right) for right in normal] for left in normal]
    )
    expected_b = sp.diag(
        sp.Matrix([[Fraction(3, 4)]]),
        sp.Matrix([[0, half], [half, 0]]),
        sp.Matrix([[0, half], [half, 0]]),
        sp.Matrix([[0, half], [half, 0]]),
    )
    expected_q = [Fraction(3, 4), 0, 0, 0, 0, 1, 1]
    dual_checks = [
        all(value.is_Integer for value in gram * sp.Matrix(generator))
        for generator in normal
    ]
    checks = {
        "all_generators_lie_in_dual": all(dual_checks),
        "orders_4_2x6": orders == [4, 2, 2, 2, 2, 2, 2],
        "generators_exhaust_all_256_classes": len(_generated_classes(normal)) == 256,
        "quadratic_values_are_normal_form": q_values == expected_q,
        "bilinear_form_is_orthogonal_normal_form": b_matrix == expected_b,
        "smith_group_agrees": _smith(gram) == [1] * 4 + [2] * 6 + [4],
    }
    if not all(checks.values()):
        raise ArithmeticError([key for key, value in checks.items() if not value])
    return {
        "group": "A_S=(Z/4)+(Z/2)^6",
        "quadratic_form": "q_S=<3/4> on Z/4 + u(2)^2 + v(2)",
        "convention": (
            "q takes values in Q/2Z; u(2) has q(e)=q(f)=0 and "
            "b(e,f)=1/2; v(2) has q(e)=q(f)=1 and b(e,f)=1/2"
        ),
        "normal_generators_in_source_dual_coordinates": [
            _vector_payload(value) for value in normal
        ],
        "generator_orders": orders,
        "generator_q_values_mod_2Z": [_text(value) for value in q_values],
        "generator_bilinear_matrix_mod_Z": _matrix_payload(b_matrix),
        "checks": checks,
    }


def _explicit_k3_embedding(
    picard_gram: sp.Matrix,
    fixed_basis: sp.Matrix,
    fixed_gram: sp.Matrix,
) -> dict[str, Any]:
    """Glue Pic(Hess) to U+U(2)+A2(-2), then include the fixed kernel."""

    U = sp.Matrix([[0, 1], [1, 0]])
    transcendental = sp.diag(
        U,
        2 * U,
        -2 * sp.Matrix([[2, -1], [-1, 2]]),
    )
    direct_sum = sp.diag(picard_gram, transcendental)
    half, third, sixth = Fraction(1, 2), Fraction(1, 3), Fraction(1, 6)
    picard_generators = [
        (0, 0, 0, 0, third, 2 * third, half, sixth, 5 * sixth, half, 0, 2 * third, third, 0, third, 2 * third),
        (0, 0, 0, half, half, half, 0, 0, 0, half, 0, 0, 0, 0, 0, 0),
        (0, half, half, 0, 0, half, 0, 0, half, 0, 0, 0, 0, 0, 0, 0),
        (half, 0, half, 0, half, 0, 0, half, 0, 0, 0, 0, 0, 0, 0, 0),
    ]
    transcendental_generators = [
        (0, 0, 0, half, third, 2 * third),
        (0, 0, half, 0, 0, 0),
        (0, 0, half, half, 0, half),
        (0, 0, half, half, half, 0),
    ]
    p_generators = [tuple(Fraction(value) for value in row) for row in picard_generators]
    t_generators = [tuple(Fraction(value) for value in row) for row in transcendental_generators]
    orders = [_class_order(value) for value in p_generators]
    p_q = [_finite_q(picard_gram, value) for value in p_generators]
    t_q = [_finite_q(transcendental, value) for value in t_generators]
    t_orders = [_class_order(value) for value in t_generators]
    p_b = sp.Matrix(
        [[_finite_b(picard_gram, left, right) for right in p_generators] for left in p_generators]
    )
    t_b = sp.Matrix(
        [[_finite_b(transcendental, left, right) for right in t_generators] for left in t_generators]
    )

    glue_vectors = [
        sp.Matrix(left + right)
        for left, right in zip(p_generators, t_generators)
    ]
    retained_standard_indices = [
        index for index in range(22) if index not in (0, 1, 3, 7)
    ]
    glue_basis = sp.eye(22)[:, retained_standard_indices].row_join(
        sp.Matrix.hstack(*glue_vectors)
    )
    ambient = sp.simplify(glue_basis.T * direct_sum * glue_basis)
    direct_to_glued = sp.simplify(glue_basis.inv())
    picard_embedding = direct_to_glued[:, :16]
    fixed_embedding = sp.simplify(picard_embedding * fixed_basis)
    rational_complement = sp.Matrix.hstack(
        *(fixed_embedding.T * ambient).nullspace()
    )
    complement_gram = sp.simplify(rational_complement.T * ambient * rational_complement)

    checks = {
        "Hessian_Picard_signature_1_15_determinant_minus48": (
            exact_inertia(picard_gram) == (1, 15, 0)
            and picard_gram.det() == -48
        ),
        "transcendental_is_U_plus_U2_plus_A2minus2": (
            exact_inertia(transcendental) == (2, 4, 0)
            and transcendental.det() == 48
        ),
        "discriminant_generator_orders_6_2_2_2": orders == [6, 2, 2, 2],
        "all_generators_lie_in_their_dual_lattices": (
            all(
                all(value.is_Integer for value in picard_gram * sp.Matrix(generator))
                for generator in p_generators
            )
            and all(
                all(value.is_Integer for value in transcendental * sp.Matrix(generator))
                for generator in t_generators
            )
        ),
        "transcendental_generator_orders_6_2_2_2": t_orders == [6, 2, 2, 2],
        "Picard_generators_exhaust_48_classes": (
            len(_generated_classes(p_generators)) == 48
        ),
        "transcendental_generators_exhaust_48_classes": (
            len(_generated_classes(t_generators)) == 48
        ),
        "quadratic_forms_are_opposites": all(
            _mod_two(left + right) == 0 for left, right in zip(p_q, t_q)
        ),
        "bilinear_forms_are_opposites": all(
            _mod_one(p_b[i, j] + t_b[i, j]) == 0
            for i in range(4)
            for j in range(4)
        ),
        "glue_basis_index_48": abs(glue_basis.det()) == sp.Rational(1, 48),
        "ambient_integral_even": (
            all(value.is_Integer for value in ambient)
            and all(int(ambient[index, index]) % 2 == 0 for index in range(22))
        ),
        "ambient_unimodular_signature_3_19": (
            abs(int(ambient.det())) == 1
            and exact_inertia(ambient) == (3, 19, 0)
        ),
        "Picard_embedding_integral_primitive_isometric": (
            all(value.is_Integer for value in picard_embedding)
            and _smith(picard_embedding) == [1] * 16
            and picard_embedding.T * ambient * picard_embedding == picard_gram
        ),
        "fixed_embedding_integral_primitive_isometric": (
            all(value.is_Integer for value in fixed_embedding)
            and _smith(fixed_embedding) == [1] * 11
            and fixed_embedding.T * ambient * fixed_embedding == fixed_gram
        ),
        "fixed_orthogonal_complement_signature_2_9": (
            rational_complement.rank() == 11
            and exact_inertia(complement_gram) == (2, 9, 0)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError([key for key, value in checks.items() if not value])
    return {
        "criterion": "explicit full-discriminant graph gluing, not a stable-range shortcut",
        "Hessian_Picard_gram": _matrix_payload(picard_gram),
        "orthogonal_transcendental_lattice": "T_H=U+U(2)+A2(-2)",
        "orthogonal_transcendental_gram": _matrix_payload(transcendental),
        "Picard_discriminant_generators": [_vector_payload(value) for value in p_generators],
        "transcendental_discriminant_generators": [_vector_payload(value) for value in t_generators],
        "generator_orders": orders,
        "Picard_generator_q_values": [_text(value) for value in p_q],
        "transcendental_generator_q_values": [_text(value) for value in t_q],
        "Picard_generator_bilinear_matrix": _matrix_payload(p_b),
        "transcendental_generator_bilinear_matrix": _matrix_payload(t_b),
        "graph_glue_basis_columns_in_Picard_plus_T_Q_space": _matrix_payload(glue_basis),
        "glued_even_unimodular_gram": _matrix_payload(ambient),
        "ambient_identification": (
            "the unique even unimodular lattice of signature (3,19), hence Lambda_K3"
        ),
        "Hessian_Picard_embedding_columns_in_glued_basis": _matrix_payload(picard_embedding),
        "fixed_lattice_embedding_columns_in_glued_basis": _matrix_payload(fixed_embedding),
        "fixed_orthogonal_complement_rational_basis_columns": _matrix_payload(rational_complement),
        "fixed_orthogonal_complement_rational_gram": _matrix_payload(complement_gram),
        "fixed_orthogonal_complement_signature": [2, 9, 0],
        "checks": checks,
    }


def _primitive_integer_ray(tokens: Iterable[str]) -> tuple[int, ...]:
    values = tuple(Fraction(token) for token in tokens)
    denominator = math.lcm(*(value.denominator for value in values))
    integers = tuple(
        value.numerator * (denominator // value.denominator) for value in values
    )
    divisor = math.gcd(*(abs(value) for value in integers))
    if divisor == 0:
        raise ArithmeticError("cddlib emitted the zero ray")
    return tuple(value // divisor for value in integers)


def _exact_cone(gram: sp.Matrix, roots: sp.Matrix) -> dict[str, Any]:
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for the exact cone replay")
    inequality = roots.T * gram
    lines = ["H-representation", "begin", f"{inequality.rows} {inequality.cols + 1} rational"]
    lines.extend(
        "0 " + " ".join(_text(value) for value in inequality.row(index))
        for index in range(inequality.rows)
    )
    lines.append("end")
    completed = subprocess.run(
        [executable, "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        timeout=600,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(completed.stdout[-4000:] + completed.stderr[-4000:])
    output_lines = completed.stdout.splitlines()
    begin = output_lines.index("begin")
    count_text, width_text, kind = output_lines[begin + 1].split()
    count, width = int(count_text), int(width_text)
    if kind != "rational" or width != 12:
        raise ArithmeticError("cddlib emitted the wrong generator format")

    causal: Counter[str] = Counter()
    incidence: Counter[tuple[str, int]] = Counter()
    facet_incidence: Counter[int] = Counter()
    seen = set()
    for line in output_lines[begin + 2 : begin + 2 + count]:
        tokens = line.split()
        if len(tokens) != width or Fraction(tokens[0]) != 0:
            raise ArithmeticError("cddlib emitted a non-ray row")
        ray = _primitive_integer_ray(tokens[1:])
        if ray in seen:
            raise ArithmeticError("cddlib emitted duplicate projective rays")
        seen.add(ray)
        vector = sp.Matrix(ray)
        slacks = inequality * vector
        if any(value < 0 for value in slacks):
            raise ArithmeticError("an exact ray violates a chamber inequality")
        active = [index for index, value in enumerate(slacks) if value == 0]
        square = int((vector.T * gram * vector)[0])
        ray_type = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        causal[ray_type] += 1
        incidence[(ray_type, len(active))] += 1
        facet_incidence.update(active)

    expected_incidence = Counter(
        {("finite", 10): 895, ("ideal", 11): 6, ("ideal", 12): 58, ("ideal", 13): 3}
    )
    checks = {
        "all_962_extreme_rays_distinct": len(seen) == 962,
        "895_finite_67_ideal_0_hyperideal": causal == {"finite": 895, "ideal": 67},
        "active_facet_histogram_exact": incidence == expected_incidence,
        "every_wall_is_an_actual_cone_facet": set(facet_incidence) == set(range(27)),
    }
    if not all(checks.values()):
        raise ArithmeticError([key for key, value in checks.items() if not value])
    return {
        "engine": "cddlib GMP rational double description",
        "extreme_ray_count": len(seen),
        "finite_vertex_count": causal["finite"],
        "ideal_vertex_count": causal["ideal"],
        "hyperideal_extreme_ray_count": causal["hyperideal"],
        "active_facet_count_histogram": {
            f"{ray_type}_{number}": value
            for (ray_type, number), value in sorted(incidence.items())
        },
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    paths = {"source": SOURCE, "stock": STOCK, "vinberg_cpp": CPP}
    for key, path in paths.items():
        actual = _sha256(path)
        if actual != EXPECTED_SHA256[key]:
            raise ArithmeticError(f"pinned {key} input changed: {actual}")
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    if source["all_checks_pass"] is not True:
        raise ArithmeticError("the generic-Hessian source certificate is not green")

    picard = _picard_model()
    fold = _fold_data(picard)
    fixed_gram = picard["fixed_gram"]
    if _matrix_payload(fixed_gram) != source["fixed_divisor_lattice"]["gram"]:
        raise ArithmeticError("the reconstructed fixed Gram changed")

    roots = []
    records = []
    types: Counter[tuple[int, int]] = Counter()
    source_records = source["restricted_chamber"]["primitive_fixed_mirror_records"]
    for facet, orbit_index in enumerate(fold["negative"]):
        quotient_coordinates = _basis_coordinates(
            fold["roots"][orbit_index],
            picard["pivots"],
            picard["quotient_gram_inverse"],
        )
        picard_coordinates = picard["picard_basis"].inv() * quotient_coordinates
        fixed_coordinates = picard["fixed_basis"].gauss_jordan_solve(picard_coordinates)[0]
        content = _content(fixed_coordinates)
        root = fixed_coordinates / content
        reflection = _reflection(fixed_gram, root)
        square = int((root.T * fixed_gram * root)[0])
        divisibility = _divisibility(fixed_gram, root)
        checks = {
            "integral": all(value.is_Integer for value in root),
            "primitive": _content(root) == 1,
            "source_content_agrees": content == source_records[facet]["orbit_sum_content_in_fixed_lattice"],
            "source_square_agrees": square == source_records[facet]["primitive_square"],
            "source_divisibility_agrees": divisibility == source_records[facet]["primitive_divisibility"],
            "reflection_integral": all(value.is_Integer for value in reflection),
            "reflection_unimodular": reflection.det() == -1,
            "reflection_involutive": reflection**2 == sp.eye(11),
            "reflection_preserves_full_Picard": reflection.T * fixed_gram * reflection == fixed_gram,
            "reflection_negates_root": reflection * root == -root,
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"wall {facet + 1}: {[key for key, value in checks.items() if not value]}"
            )
        roots.append(root)
        types[(square, divisibility)] += 1
        records.append(
            {
                "facet_one_based": facet + 1,
                "source_wall_orbit_one_based": [
                    value + 1 for value in fold["cycles"][orbit_index]
                ],
                "primitive_fixed_Picard_coordinates": _vector_payload(root),
                "orbit_sum_content": content,
                "square": square,
                "divisibility_in_full_Picard": divisibility,
                "checks": checks,
            }
        )
    root_matrix = sp.Matrix.hstack(*roots)
    primitive_facet_gram = sp.simplify(root_matrix.T * fixed_gram * root_matrix)

    ample_formal = _vector((kind, label, 1) for kind, label in FORMAL_LABELS)
    ample_quotient = _basis_coordinates(
        ample_formal, picard["pivots"], picard["quotient_gram_inverse"]
    )
    ample_picard = picard["picard_basis"].inv() * ample_quotient
    ample_fixed = picard["fixed_basis"].gauss_jordan_solve(ample_picard)[0]
    ample_pairings = [
        int(value) for value in (ample_fixed.T * fixed_gram * root_matrix)
    ]
    for record, pairing in zip(records, ample_pairings):
        record["strict_interior_pairing"] = pairing

    split = _u_split(fixed_gram)
    discriminant = _discriminant_form(fixed_gram)
    embedding = _explicit_k3_embedding(
        picard["picard_gram"], picard["fixed_basis"], fixed_gram
    )
    vinberg = _vinberg(primitive_facet_gram)
    exact_cone = _exact_cone(fixed_gram, root_matrix)

    lattice_checks = {
        "saturated_fixed_rank11_signature_1_10": (
            picard["fixed_basis"].cols == 11
            and source["fixed_divisor_lattice"]["saturation_audit"]["saturated_smith_kernel"]
            and exact_inertia(fixed_gram) == (1, 10, 0)
        ),
        "even_determinant256_smith_1x4_2x6_4": (
            fixed_gram.det() == 256
            and _smith(fixed_gram) == [1] * 4 + [2] * 6 + [4]
            and all(int(fixed_gram[index, index]) % 2 == 0 for index in range(11))
        ),
        "all_27_primitive_integral_reflections": all(
            all(record["checks"].values()) for record in records
        ),
        "mirror_type_histogram": dict(types)
        == {(-2, 1): 8, (-4, 2): 10, (-4, 4): 9},
        "walls_span_full_Picard_index_one": (
            root_matrix.rank() == 11 and _smith(root_matrix.T) == [1] * 11
        ),
        "strict_square20_interior": (
            (ample_fixed.T * fixed_gram * ample_fixed)[0] == 20
            and min(ample_pairings) == 1
        ),
        "facet_relation_matrix_lineage_unchanged": (
            fold["facet_gram"]
            == sp.Matrix(source["restricted_chamber"]["facet_gram"])
            and
            source["lineage"]["normalized_relation_matrices_match_entrywise"]
            and source["lineage"]["known_id"] == "D16-D4-C7-H10-27"
        ),
    }
    volume_checks = {
        "literal_Vinberg_895_finite": vinberg["finite_vertex_count"] == 895,
        "literal_Vinberg_67_cusps": vinberg["ideal_vertex_count"] == 67,
        "literal_Vinberg_6221_edges": vinberg["elliptic_edge_count"] == 6221,
        "literal_Vinberg_zero_edge_or_cusp_failures": (
            vinberg["bad_edge_count"] == 0
            and vinberg["bad_parabolic_neighborhood_count"] == 0
        ),
        "literal_Vinberg_finite_volume_noncompact": (
            vinberg["finite_volume"] and not vinberg["compact"]
        ),
        "exact_cone_895_67_0": (
            exact_cone["finite_vertex_count"] == 895
            and exact_cone["ideal_vertex_count"] == 67
            and exact_cone["hyperideal_extreme_ray_count"] == 0
        ),
    }

    sources = {
        key: {"path": str(path.relative_to(ROOT)), "sha256": EXPECTED_SHA256[key]}
        for key, path in paths.items()
    }
    checks = {
        "source_hashes_pinned": all(
            _sha256(path) == EXPECTED_SHA256[key] for key, path in paths.items()
        ),
        "fixed_lattice_and_walls_green": all(lattice_checks.values()),
        "discriminant_form_green": all(discriminant["checks"].values()),
        "explicit_primitive_K3_embedding_green": all(embedding["checks"].values()),
        "literal_Vinberg_and_exact_cone_green": all(volume_checks.values()),
    }
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "theorem": (
            "A nine-dimensional very-general family of projective complex K3 "
            "surfaces has full NS=S_H, the saturated rank-11 determinant-256 "
            "generic-Hessian transposition-fixed lattice, and carries all 27 "
            "primitive mirrors of the known finite-volume noncompact H10/27 "
            "chamber; those mirrors generate full Picard."
        ),
        "sources": sources,
        "full_Picard_lattice": {
            "name": "S_H (generic-Hessian transposition-fixed divisor lattice)",
            "rank": 11,
            "signature": [1, 10, 0],
            "determinant": 256,
            "smith_invariants": _smith(fixed_gram),
            "gram_in_source_fixed_basis": _matrix_payload(fixed_gram),
            "standard_identification": split,
            "discriminant_form": discriminant,
            "primitive_K3_embedding": embedding,
            "period_consequence": {
                "K3_moduli_dimension": 9,
                "very_general_full_Picard": "NS(X)=S_H",
                "reason": (
                    "the displayed primitive embedding defines a nonempty "
                    "S_H-polarized K3 period domain; outside the countable "
                    "Noether--Lefschetz hyperplanes there is no extra algebraic class"
                ),
            },
            "checks": lattice_checks,
        },
        "transported_chamber": {
            "id": "GHESS-C2-H10-27",
            "known_abstract_id": "D16-D4-C7-H10-27",
            "hyperbolic_dimension": 10,
            "facet_count": 27,
            "primitive_mirror_type_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(types.items())
            },
            "wall_span_rank": root_matrix.rank(),
            "wall_span_smith": _smith(root_matrix.T),
            "wall_span_index_in_full_Picard": math.prod(_smith(root_matrix.T)),
            "strict_interior_coordinates": _vector_payload(ample_fixed),
            "strict_interior_square": int((ample_fixed.T * fixed_gram * ample_fixed)[0]),
            "strict_wall_pairings": ample_pairings,
            "minimum_wall_pairing": min(ample_pairings),
            "wall_records": records,
            "primitive_facet_gram": _matrix_payload(primitive_facet_gram),
            "primitive_facet_gram_sha256": _payload_sha256(
                _matrix_payload(primitive_facet_gram)
            ),
        },
        "finite_volume": {
            "classification": "finite_volume_noncompact_Coxeter_polytope",
            "literal_Vinberg": vinberg,
            "exact_GMP_cone": exact_cone,
            "checks": volume_checks,
        },
        "lineage": {
            "normalized_diagram": "D16-D4-C7-H10-27",
            "known_full_chamber_facets": source["lineage"]["known_full_reflection_chamber_facets"],
            "known_reflection_subgroup_index": source["lineage"]["known_reflection_subgroup_index"],
            "new_abstract_polytope": False,
            "new_arithmetic_class": False,
        },
        "scope": {
            "positive_gain": (
                "new project-level very-general complex-K3 full-Picard "
                "divisor-reflection carrier of the known H10/27 chamber"
            ),
            "not_equation_level": True,
            "not_full_nef": (
                "the nineteen primitive square-minus-four mirrors are mixed "
                "divisor cuts, not K3 negative-curve walls"
            ),
            "not_surface_automorphism": (
                "the source transposition and the Picard reflections are not "
                "claimed to extend to automorphisms of the very-general K3"
            ),
            "not_new_abstract_polytope": True,
            "worldwide_priority_claim": "none",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["payload_sha256"] = _payload_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    if not payload["all_checks_pass"]:
        raise ArithmeticError("GHESS H10/27 full-Picard carrier failed")
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise ArithmeticError("archived GHESS H10/27 carrier changed")
        print(f"verified {args.output}")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"wrote {args.output}")
    print(f"payload_sha256={payload['payload_sha256']}")


if __name__ == "__main__":
    main()
