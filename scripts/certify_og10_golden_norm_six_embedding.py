#!/usr/bin/env python3
"""Exact primitive OG10 embedding of the simple golden norm-six carrier.

The field lattice is H6(k) with k=(0^9,1).  Its trace lattice N has rank 22,
signature (1,21), determinant -36, and discriminant group (Z/6)^2.  This
certificate glues N to T=<2>+<6> along an explicit order-12 isotropic graph.
The resulting even lattice is identified by a literal 24-by-24 isometry with
the OG10 Beauville lattice 3U+2E8(-1)+A2(-1).  The selected trace lift of the
norm-six field root is primitive of square -6 and ambient divisibility 3.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import reduce
import hashlib
from itertools import product
import json
from math import gcd
from pathlib import Path
from typing import Any, Iterable

from sympy import Matrix, Rational, eye
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.exact.matrix import determinant, inertia
from coxeter_polytopes.exact.trace_lattice import (
    INV_SQRT5,
    integral_gram,
    is_even_integral,
    trace_gram,
)
from coxeter_polytopes.lattices import e8_cartan_gram, hyperbolic_plane
from coxeter_polytopes.search.og10_golden_h9_vinberg import h9_field_gram
from scripts.certify_og10_golden_rank_two_reduction import _append_block


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_PATH = (
    ROOT / "data/candidates/og10_golden_norm_six_embedding.json"
)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_json_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _fraction(value: Any) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    numerator, denominator = value.as_numer_denom()
    return Fraction(int(numerator), int(denominator))


def _mod_one(value: Any) -> Fraction:
    value = _fraction(value)
    return value - value.numerator // value.denominator


def _mod_two(value: Any) -> Fraction:
    value = _fraction(value)
    return value - 2 * (value // 2)


def _scalar_json(value: Any) -> int | str:
    value = _fraction(value)
    if value.denominator == 1:
        return value.numerator
    return f"{value.numerator}/{value.denominator}"


def _matrix_json(matrix: Matrix) -> list[list[int | str]]:
    return [
        [_scalar_json(matrix[i, j]) for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]


def _matrix_digest(matrix: Matrix) -> str:
    encoded = json.dumps(_matrix_json(matrix), separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _integer_matrix(matrix: Matrix) -> Matrix:
    if not all(entry.is_Integer is True for entry in matrix):
        raise AssertionError("expected an integral matrix")
    return matrix.applyfunc(int)


def _class_key(vector: Matrix) -> tuple[Fraction, ...]:
    return tuple(_mod_one(entry) for entry in vector)


def _class_order(vector: Matrix, bound: int = 12) -> int:
    zero = tuple(Fraction(0) for _ in vector)
    for order in range(1, bound + 1):
        if _class_key(order * vector) == zero:
            return order
    raise AssertionError("discriminant class order exceeded the bound")


def _span_keys(first: Matrix, second: Matrix) -> set[tuple[Fraction, ...]]:
    return {
        _class_key(a * first + b * second)
        for a in range(2)
        for b in range(6)
    }


def _smith_diagonal(matrix: Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[i, i]))
        for i in range(min(smith.rows, smith.cols))
        if smith[i, i] != 0
    ]


def _gcd_entries(vector: Matrix) -> int:
    return reduce(gcd, (abs(int(entry)) for entry in vector), 0)


def field_and_trace_grams() -> tuple[list[list[Golden]], Matrix]:
    field = _append_block([list(row) for row in h9_field_gram()], [[Golden(6)]])
    field[9][10] = field[10][9] = Golden(3)
    trace = Matrix(integral_gram(trace_gram(field, -PHI * INV_SQRT5)))
    return field, trace


def small_trace_gram() -> Matrix:
    return Matrix(
        [
            [2, 4, 0, 0, 0, 0],
            [4, 6, 0, 0, 0, 0],
            [0, 0, -2, -2, -3, -3],
            [0, 0, -2, -4, -3, -6],
            [0, 0, -3, -3, -6, -6],
            [0, 0, -3, -6, -6, -12],
        ]
    )


def trace_split_change() -> Matrix:
    """Put N in the order 2E8(-1) plus (x,u,e) trace coordinates."""

    columns = []
    for index in range(8):
        column = Matrix.zeros(22, 1)
        column[2 + 2 * index] = 1
        columns.append(column)
    for index in range(8):
        column = Matrix.zeros(22, 1)
        column[2 + 2 * index] = -1
        column[2 + 2 * index + 1] = 1
        columns.append(column)
    for position in (0, 1, 18, 19, 20, 21):
        column = Matrix.zeros(22, 1)
        column[position] = 1
        columns.append(column)
    return Matrix.hstack(*columns)


def complement_gram() -> Matrix:
    return Matrix.diag(2, 6)


def chosen_gluing_vectors() -> tuple[Matrix, Matrix]:
    """Order-2 and order-6 generators in A_N plus A_T."""

    first = Matrix([0, Rational(1, 2), 0, 0, 0, 0, Rational(1, 2), 0])
    second = Matrix(
        [
            Rational(1, 2),
            0,
            0,
            0,
            Rational(-1, 3),
            Rational(1, 3),
            0,
            Rational(1, 6),
        ]
    )
    return first, second


def small_overlattice_basis() -> Matrix:
    """Columns form the glued rank-eight lattice in N_small plus T."""

    first, second = chosen_gluing_vectors()
    return Matrix.hstack(
        Matrix.vstack(eye(6), Matrix.zeros(2, 6)),
        first,
        second,
    )


def small_og10_isometry() -> Matrix:
    """Take the glued small Gram to 3U plus A2(-1)."""

    return Matrix(
        [
            [0, 1, -1, 1, 2, 1, 1, 2],
            [0, 0, 0, -1, 0, 0, 0, 0],
            [0, 0, 0, 0, -1, -1, 0, -3],
            [0, 0, 0, 0, 0, 1, 1, 2],
            [0, 0, 0, 0, 0, 0, 0, 1],
            [0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 1, 1, 0, 0, 0, 0],
            [1, -1, 0, 0, -4, -2, -2, -4],
        ]
    )


def standard_og10_gram() -> Matrix:
    u = Matrix(hyperbolic_plane())
    negative_e8 = -Matrix(e8_cartan_gram())
    a2_negative = Matrix([[-2, -1], [-1, -2]])
    return Matrix.diag(u, u, u, negative_e8, negative_e8, a2_negative)


def full_embedding_data(trace: Matrix) -> dict[str, Matrix]:
    negative_e8 = -Matrix(e8_cartan_gram())
    small_n = small_trace_gram()
    t = complement_gram()
    split = trace_split_change()
    split_target = Matrix.diag(negative_e8, negative_e8, small_n)
    if split.T * trace * split != split_target or abs(int(split.det())) != 1:
        raise AssertionError("trace splitting failed")

    direct_sum = Matrix.diag(trace, t)
    reordered_direct_change = Matrix.diag(split, eye(2))
    small_basis = small_overlattice_basis()
    full_gluing_basis = reordered_direct_change * Matrix.diag(eye(16), small_basis)
    ambient = _integer_matrix(full_gluing_basis.T * direct_sum * full_gluing_basis)

    current_isometry = Matrix.diag(eye(16), small_og10_isometry())
    # Reorder 2E8(-1)+3U+A2(-1) to 3U+2E8(-1)+A2(-1).
    source_order = [*range(16, 22), *range(0, 16), 22, 23]
    permutation = Matrix.zeros(24)
    for target, source in enumerate(source_order):
        permutation[source, target] = 1
    ambient_isometry = _integer_matrix(current_isometry * permutation)
    standard = standard_og10_gram()
    if ambient_isometry.T * ambient * ambient_isometry != standard:
        raise AssertionError("the explicit 24x24 OG10 isometry failed")

    ambient_inverse = full_gluing_basis.inv()
    n_inclusion = Matrix.vstack(eye(22), Matrix.zeros(2, 22))
    t_inclusion = Matrix.vstack(Matrix.zeros(22, 2), eye(2))
    n_embedding = _integer_matrix(ambient_inverse * n_inclusion)
    t_embedding = _integer_matrix(ambient_inverse * t_inclusion)
    combined_standard_basis = full_gluing_basis * ambient_isometry
    return {
        "trace_split_change": split,
        "direct_sum_gram": direct_sum,
        "full_gluing_basis": full_gluing_basis,
        "ambient_gram": ambient,
        "ambient_isometry": ambient_isometry,
        "standard_og10_gram": standard,
        "combined_standard_basis": combined_standard_basis,
        "n_embedding": n_embedding,
        "t_embedding": t_embedding,
    }


def discriminant_graph_audit() -> dict[str, object]:
    n = small_trace_gram()
    t = complement_gram()
    direct = Matrix.diag(n, t)
    n_classes = [
        Matrix(
            [
                Rational(a, 2),
                Rational(b, 2),
                0,
                0,
                Rational(c, 3),
                Rational(d, 3),
            ]
        )
        for a, b, c, d in product(range(2), range(2), range(3), range(3))
    ]
    if len({_class_key(vector) for vector in n_classes}) != 36:
        raise AssertionError("the explicit A_N representatives are incomplete")
    if any(
        not all(entry.is_Integer is True for entry in n * vector)
        for vector in n_classes
    ):
        raise AssertionError("an A_N representative is not dual")

    q_n = lambda vector: _mod_two((vector.T * n * vector)[0])
    b_n = lambda left, right: _mod_one((left.T * n * right)[0])
    d_class = Matrix([0, 0, 0, 0, Rational(1, 3), 0])

    candidates: list[tuple[Matrix, Matrix]] = []
    selected_candidates: list[tuple[Matrix, Matrix]] = []
    for image_order_two in n_classes:
        if _class_order(image_order_two) != 2:
            continue
        if q_n(image_order_two) != Fraction(3, 2):
            continue
        for image_order_six in n_classes:
            if _class_order(image_order_six) != 6:
                continue
            if q_n(image_order_six) != Fraction(11, 6):
                continue
            if b_n(image_order_two, image_order_six) != 0:
                continue
            if len(_span_keys(image_order_two, image_order_six)) != 12:
                continue
            candidates.append((image_order_two, image_order_six))
            if (
                b_n(d_class, image_order_two) == 0
                and b_n(d_class, image_order_six) == 0
            ):
                selected_candidates.append((image_order_two, image_order_six))

    first, second = chosen_gluing_vectors()
    first_n, first_t = first[:6, :], first[6:, :]
    second_n, second_t = second[:6, :], second[6:, :]
    chosen_key = (_class_key(first_n), _class_key(second_n))
    selected_keys = {
        (_class_key(left), _class_key(right))
        for left, right in selected_candidates
    }
    if chosen_key not in selected_keys:
        raise AssertionError("chosen graph does not preserve the selected D class")

    q_first = _fraction((first.T * direct * first)[0])
    q_second = _fraction((second.T * direct * second)[0])
    b_first_second = _fraction((first.T * direct * second)[0])
    n_span = _span_keys(first_n, second_n)
    t_span = _span_keys(first_t, second_t)
    total_span = _span_keys(first, second)
    selected_sign_orbit = {
        (_class_key(first_n), _class_key(second_n)),
        (_class_key(first_n), _class_key(-second_n)),
    }
    return {
        "n_discriminant_smith": _smith_diagonal(n),
        "n_discriminant_group": "Z/6 + Z/6",
        "n_discriminant_class_count": len(n_classes),
        "t_discriminant_smith": _smith_diagonal(t),
        "t_discriminant_group": "Z/2 + Z/6",
        "t_discriminant_class_count": 12,
        "all_injective_anti_isometry_generator_pairs": len(candidates),
        "pairs_preserving_selected_D_over_3": len(selected_candidates),
        "selected_pairs_form_one_complement_sign_orbit": (
            selected_keys == selected_sign_orbit
        ),
        "all_generator_images": [
            {
                "order_two_image": _matrix_json(left.T)[0],
                "order_six_image": _matrix_json(right.T)[0],
                "selected_D_survives": (
                    b_n(d_class, left) == 0 and b_n(d_class, right) == 0
                ),
            }
            for left, right in candidates
        ],
        "chosen_order_two_generator": _matrix_json(first.T)[0],
        "chosen_order_six_generator": _matrix_json(second.T)[0],
        "chosen_generator_orders": [_class_order(first), _class_order(second)],
        "chosen_generator_squares": [
            _scalar_json(q_first),
            _scalar_json(q_second),
        ],
        "chosen_generator_pairing": _scalar_json(b_first_second),
        "chosen_graph_image_orders": {
            "N": len(n_span),
            "T": len(t_span),
            "N_plus_T": len(total_span),
        },
        "chosen_graph_isotropic": (
            _mod_two(q_first) == 0
            and _mod_two(q_second) == 0
            and b_first_second.denominator == 1
        ),
        "selected_D_over_3_square_mod_2Z": _scalar_json(q_n(d_class)),
        "selected_D_over_3_orthogonal_to_graph": (
            b_n(d_class, first_n) == 0 and b_n(d_class, second_n) == 0
        ),
    }


def _source_hashes() -> dict[str, str]:
    paths = {
        "certifier": Path(__file__).resolve(),
        "rank_two_norm_six_compiler": ROOT
        / "scripts/certify_og10_golden_rank_two_reduction.py",
        "golden_arithmetic": ROOT / "src/coxeter_polytopes/exact/golden.py",
        "trace_lattice": ROOT / "src/coxeter_polytopes/exact/trace_lattice.py",
        "h9_fixed_field_gram": ROOT
        / "src/coxeter_polytopes/search/og10_golden_h9_vinberg.py",
        "standard_lattices": ROOT / "src/coxeter_polytopes/lattices.py",
    }
    return {name: _sha256_file(path) for name, path in paths.items()}


def build_certificate() -> dict[str, object]:
    field, trace = field_and_trace_grams()
    data = full_embedding_data(trace)
    graph = discriminant_graph_audit()
    ambient = data["ambient_gram"]
    n_embedding = data["n_embedding"]
    t_embedding = data["t_embedding"]

    selected_d = Matrix([0] * 20 + [1, 0])
    golden_partner = Matrix([0] * 20 + [-1, 1])
    selected_ambient = n_embedding * selected_d
    partner_ambient = n_embedding * golden_partner
    selected_pairing_row = ambient * selected_ambient
    partner_pairing_row = ambient * partner_ambient
    positive_class = Matrix([1] + [0] * 21)
    positive_ambient = n_embedding * positive_class
    selected_pair_gram = Matrix.hstack(selected_d, golden_partner).T * trace * (
        Matrix.hstack(selected_d, golden_partner)
    )

    n_smith = _smith_diagonal(n_embedding)
    t_smith = _smith_diagonal(t_embedding)
    ambient_smith = _smith_diagonal(ambient)
    trace_inertia = inertia([[int(entry) for entry in row] for row in trace.tolist()])
    ambient_inertia = inertia(
        [[int(entry) for entry in row] for row in ambient.tolist()]
    )
    field_identity_inertia = inertia(field)
    field_conjugate_inertia = inertia(
        [[entry.conjugate() for entry in row] for row in field]
    )

    small_direct = Matrix.diag(small_trace_gram(), complement_gram())
    small_basis = small_overlattice_basis()
    first, second = chosen_gluing_vectors()
    generator_coordinates = small_basis.inv() * Matrix.hstack(eye(8), first, second)
    checks = {
        "field_signatures_10_1_and_11_0": (
            (field_identity_inertia.positive, field_identity_inertia.negative)
            == (10, 1)
            and (
                field_conjugate_inertia.positive,
                field_conjugate_inertia.negative,
            )
            == (11, 0)
        ),
        "field_norm_six_root_defines_an_integral_reflection": (
            field[10][10] == Golden(6)
            and all((entry / 3).is_integral for entry in field[10])
        ),
        "trace_even_rank22_signature_1_21_determinant_minus36": (
            is_even_integral(trace.tolist())
            and trace.rows == 22
            and (trace_inertia.positive, trace_inertia.negative) == (1, 21)
            and int(trace.det()) == -36
        ),
        "discriminant_graph_enumeration_complete": (
            graph["n_discriminant_class_count"] == 36
            and graph["t_discriminant_class_count"] == 12
            and graph["all_injective_anti_isometry_generator_pairs"] == 4
            and graph["pairs_preserving_selected_D_over_3"] == 2
            and graph["selected_pairs_form_one_complement_sign_orbit"] is True
        ),
        "chosen_graph_is_order12_isotropic_and_injective": (
            graph["chosen_graph_isotropic"] is True
            and graph["chosen_graph_image_orders"]
            == {"N": 12, "T": 12, "N_plus_T": 12}
            and graph["selected_D_over_3_orthogonal_to_graph"] is True
        ),
        "small_overlattice_basis_index12_contains_all_generators": (
            abs(small_basis.det()) == Rational(1, 12)
            and all(entry.is_Integer is True for entry in generator_coordinates)
        ),
        "ambient_even_signature_3_21_determinant_minus3": (
            all(int(ambient[i, i]) % 2 == 0 for i in range(24))
            and (ambient_inertia.positive, ambient_inertia.negative) == (3, 21)
            and int(ambient.det()) == -3
            and ambient_smith == [1] * 23 + [3]
        ),
        "explicit_24x24_isometry_to_standard_og10": (
            abs(int(data["ambient_isometry"].det())) == 1
            and data["ambient_isometry"].T
            * ambient
            * data["ambient_isometry"]
            == data["standard_og10_gram"]
        ),
        "N_embedding_is_primitive_and_isometric": (
            n_smith == [1] * 22
            and n_embedding.T * ambient * n_embedding == trace
        ),
        "complement_embedding_is_primitive_orthogonal_and_diag_2_6": (
            t_smith == [1, 1]
            and t_embedding.T * ambient * t_embedding == complement_gram()
            and n_embedding.T * ambient * t_embedding == Matrix.zeros(22, 2)
        ),
        "selected_D_is_primitive_square_minus6_divisibility3": (
            _gcd_entries(selected_ambient) == 1
            and (selected_ambient.T * ambient * selected_ambient)[0] == -6
            and _gcd_entries(selected_pairing_row) == 3
        ),
        "golden_partner_is_primitive_square_minus6_divisibility1": (
            _gcd_entries(partner_ambient) == 1
            and (partner_ambient.T * ambient * partner_ambient)[0] == -6
            and _gcd_entries(partner_pairing_row) == 1
        ),
        "paired_trace_lifts_are_orthogonal_square_minus6": (
            selected_pair_gram == Matrix.diag(-6, -6)
        ),
        "positive_square2_class_proves_projectivity": (
            (positive_ambient.T * ambient * positive_ambient)[0] == 2
        ),
        "positive_rank2_period_complement": complement_gram() == Matrix.diag(2, 6),
    }
    if not all(checks.values()):
        raise AssertionError(f"norm-six embedding checks failed: {checks}")

    matrix_payload = {
        name: _matrix_json(matrix)
        for name, matrix in {
            "trace_gram_22": trace,
            "trace_split_change_22": data["trace_split_change"],
            "small_trace_gram_6": small_trace_gram(),
            "complement_gram_2": complement_gram(),
            "small_overlattice_basis_8": small_basis,
            "small_glued_gram_8": small_basis.T * small_direct * small_basis,
            "small_og10_isometry_8": small_og10_isometry(),
            "full_gluing_basis_in_N_plus_T_24": data["full_gluing_basis"],
            "ambient_gram_24": ambient,
            "ambient_to_standard_og10_isometry_24": data["ambient_isometry"],
            "standard_og10_gram_24": data["standard_og10_gram"],
            "combined_standard_basis_in_N_plus_T_24": data[
                "combined_standard_basis"
            ],
            "primitive_N_embedding_24x22": n_embedding,
            "primitive_T_embedding_24x2": t_embedding,
        }.items()
    }
    matrix_digests = {
        name: _matrix_digest(matrix)
        for name, matrix in {
            "trace_gram_22": trace,
            "full_gluing_basis_in_N_plus_T_24": data["full_gluing_basis"],
            "ambient_gram_24": ambient,
            "ambient_to_standard_og10_isometry_24": data["ambient_isometry"],
            "combined_standard_basis_in_N_plus_T_24": data[
                "combined_standard_basis"
            ],
            "primitive_N_embedding_24x22": n_embedding,
            "primitive_T_embedding_24x2": t_embedding,
        }.items()
    }

    payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "candidate_id": "og10_golden_norm_six_simple_carrier_embedding",
        "certificate_kind": "exact_primitive_OG10_lattice_embedding",
        "field_carrier": {
            "definition": "H6(k)=[[F,3k],[3k^t,6]], F=<-2phi>+E8+<2>",
            "normalized_covector": "k=(0^9,1)",
            "field_rank": 11,
            "field_determinant": str(determinant(field)),
            "field_signatures": {
                "identity": [10, 1],
                "conjugate": [11, 0],
            },
        },
        "trace_lattice": {
            "rank": 22,
            "signature": [1, 21],
            "determinant": -36,
            "smith_diagonal": _smith_diagonal(trace),
            "split_type": "2E8(-1) plus explicit rank-six lattice",
        },
        "discriminant_form_gluing": graph,
        "gluing": {
            "complement": "<2>+<6>",
            "graph_order": 12,
            "overlattice_index": 12,
            "full_gluing_basis_determinant": _scalar_json(
                data["full_gluing_basis"].det()
            ),
            "ambient_determinant": int(ambient.det()),
            "ambient_smith_diagonal": ambient_smith,
        },
        "primitive_embeddings": {
            "N_embedding_smith_diagonal": n_smith,
            "T_embedding_smith_diagonal": t_smith,
            "orthogonal_complement_gram": _matrix_json(
                t_embedding.T * ambient * t_embedding
            ),
        },
        "selected_trace_lifts": {
            "field_reflection": {
                "formula": "s_e(x)=x-H6(x,e)*e/3",
                "integrality_reason": "H6(e,L) is contained in 3*Z[phi]",
                "restricted_mirror": (
                    "On the identity eigenspace, B(e,x)=(-phi/sqrt(5))*"
                    "H6(e,x), so D^perp restricts to the field mirror e^perp."
                ),
            },
            "D_equals_e": {
                "ambient_coordinates": _matrix_json(selected_ambient.T)[0],
                "square": int((selected_ambient.T * ambient * selected_ambient)[0]),
                "ambient_pairing_row": _matrix_json(selected_pairing_row.T)[0],
                "ambient_divisibility": _gcd_entries(selected_pairing_row),
                "primitive": _gcd_entries(selected_ambient) == 1,
                "discriminant_value": "-2/3 = 4/3 mod 2Z",
            },
            "E_equals_phi_minus_1_times_e": {
                "ambient_coordinates": _matrix_json(partner_ambient.T)[0],
                "square": int((partner_ambient.T * ambient * partner_ambient)[0]),
                "ambient_pairing_row": _matrix_json(partner_pairing_row.T)[0],
                "ambient_divisibility": _gcd_entries(partner_pairing_row),
                "primitive": _gcd_entries(partner_ambient) == 1,
            },
        },
        "period_projectivity_and_MBM": {
            "period_complement": "T=<2>+<6> is positive definite of rank two",
            "explicit_period_line": "sigma=t1+i*t2/sqrt(3), so sigma^2=0 and (sigma,conjugate(sigma))=4",
            "full_NS_statement": (
                "Period surjectivity for OG10 applied to this primitive "
                "embedding realizes N as the full Picard lattice."
            ),
            "positive_class": {
                "trace_coordinates": _matrix_json(positive_class.T)[0],
                "ambient_coordinates": _matrix_json(positive_ambient.T)[0],
                "square": 2,
            },
            "projectivity_input": (
                "A marked OG10 manifold with a positive integral (1,1) "
                "class is projective."
            ),
            "MBM_input": (
                "Mongardi-Onorati, Proposition 3.1 (arXiv:2010.12511): "
                "a multiple of a primitive OG10 class of square -6 and "
                "divisibility 3 is stably prime-exceptional. Their "
                "Proposition 5.4 says the primitive class is a wall divisor; "
                "this is the cohomological/MBM wall statement used here."
            ),
            "selected_D_is_an_OG10_MBM_wall_class": True,
        },
        "matrices": matrix_payload,
        "matrix_sha256": matrix_digests,
        "checks": checks,
        "source_sha256": _source_hashes(),
        "commands": [
            "PYTHONPATH=src:. python3 scripts/certify_og10_golden_norm_six_embedding.py --write-certificate",
            "PYTHONPATH=src:. python3 scripts/certify_og10_golden_norm_six_embedding.py --verify-certificate",
            "PYTHONPATH=src:. pytest -q tests/test_og10_golden_norm_six_embedding.py",
        ],
        "conclusion": (
            "The simple norm-six golden trace lattice embeds primitively as "
            "the full Picard lattice of a projective OG10-type manifold, and "
            "the selected trace lift e is a primitive ambient (-6,3) MBM class."
        ),
        "scope": (
            "This is a positive algebro-geometric carrier and wall theorem. "
            "It does not certify a finite Coxeter chamber, compactness, finite "
            "volume, or a new hyperbolic Coxeter polytope."
        ),
    }
    payload["certificate_payload_sha256"] = _canonical_json_sha256(payload)
    return payload


def verify_certificate(payload: dict[str, object]) -> dict[str, bool]:
    body = dict(payload)
    digest = body.pop("certificate_payload_sha256", None)
    recomputed = build_certificate()
    checks = {
        "payload_digest_matches": digest == _canonical_json_sha256(body),
        "source_hashes_match": payload.get("source_sha256") == _source_hashes(),
        "full_exact_recomputation_matches": payload == recomputed,
        "all_frozen_checks_are_true": all(payload["checks"].values()),
    }
    if not all(checks.values()):
        raise AssertionError(f"norm-six embedding verification failed: {checks}")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-certificate", action="store_true")
    parser.add_argument("--verify-certificate", action="store_true")
    args = parser.parse_args()
    if args.write_certificate:
        payload = build_certificate()
        CERTIFICATE_PATH.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {CERTIFICATE_PATH}")
        print(payload["certificate_payload_sha256"])
        return
    payload = json.loads(CERTIFICATE_PATH.read_text())
    print(json.dumps(verify_certificate(payload), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
