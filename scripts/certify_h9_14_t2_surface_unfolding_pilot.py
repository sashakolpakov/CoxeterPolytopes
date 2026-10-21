#!/usr/bin/env python3
"""Surface-free finite unfolding compiler for the T2/C01 H9/14 Gram."""

from __future__ import annotations

import argparse
import ast
from collections import Counter, deque
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
FACTORY = ROOT / "data/enumerations/alexeev_k3_h9_14_factory.json"
PROMOTION = ROOT / "data/enumerations/a18_c01_c07_k3_nef_promotions.json"
ANCILLARY = Path("/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt")
OUTPUT = ROOT / "data/enumerations/h9_14_t2_surface_unfolding_pilot.json"

FACTORY_SHA256 = "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775"
PROMOTION_SHA256 = "2ab7acc379e48ac55944e92b8d614c0992d7c482af1300c670f2bf8e899ba1b1"
ANCILLARY_SHA256 = "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
FACET_GRAM_SHA256 = "6979b16de7583c20277077a142df863e9f9a7209c443e9252505234b15d66c56"
ACTION_SHA256 = "288839c3574fb374424ee7500e25dbb56030cf1515cf8f5428a3dce98e2739be"
SCHEMA = "coxeter.h9_14_t2_surface_unfolding_pilot.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json_sha256(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _matrix_rows(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _smith(matrix: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(diagonal[i, i]))
        for i in range(min(diagonal.shape))
        if diagonal[i, i]
    ]


def _bit_vector(value: int, size: int) -> sp.Matrix:
    return sp.Matrix([(value >> index) & 1 for index in range(size)])


def _vector_bit(vector: Sequence[Any]) -> int:
    return sum((int(vector[index]) & 1) << index for index in range(len(vector)))


def _reduce_basis(vectors: Iterable[int], size: int) -> tuple[int, ...]:
    basis = [0] * size
    for original in vectors:
        value = int(original)
        while value:
            pivot = value.bit_length() - 1
            if basis[pivot]:
                value ^= basis[pivot]
                continue
            basis[pivot] = value
            for lower in range(pivot):
                if basis[lower] and ((basis[pivot] >> lower) & 1):
                    basis[pivot] ^= basis[lower]
            for upper in range(pivot + 1, size):
                if basis[upper] and ((basis[upper] >> pivot) & 1):
                    basis[upper] ^= basis[pivot]
            break
    return tuple(value for value in basis if value)


def _span(basis: Sequence[int]) -> tuple[int, ...]:
    values = [0]
    for vector in basis:
        values += [value ^ vector for value in values]
    return tuple(values)


def _gf2_rank_rows(rows: Sequence[Sequence[int]]) -> int:
    return len(
        _reduce_basis(
            (_vector_bit(row) for row in rows), len(rows[0]) if rows else 0
        )
    )


def _relation(value: int, basis: Sequence[int]) -> int:
    for coefficients in range(1 << len(basis)):
        image = 0
        for index, vector in enumerate(basis):
            if (coefficients >> index) & 1:
                image ^= vector
        if image == value:
            return coefficients
    raise ArithmeticError("GF(2) vector left its advertised span")


def _inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    work = sp.Matrix(matrix)
    positive = negative = zero = 0
    while work.rows:
        diagonal = next((index for index in range(work.rows) if work[index, index]), None)
        if diagonal is not None:
            work.row_swap(0, diagonal)
            work.col_swap(0, diagonal)
            pivot = work[0, 0]
            positive += 1 if pivot > 0 else 0
            negative += 1 if pivot < 0 else 0
            if work.rows == 1:
                break
            column = work[1:, 0]
            work = (work[1:, 1:] - column * column.T / pivot).applyfunc(sp.cancel)
            continue
        off_diagonal = next(
            (
                (left, right)
                for left in range(work.rows)
                for right in range(left + 1, work.cols)
                if work[left, right]
            ),
            None,
        )
        if off_diagonal is None:
            zero += work.rows
            break
        left, right = off_diagonal
        work.row_swap(0, left)
        work.col_swap(0, left)
        if right == 0:
            right = left
        work.row_swap(1, right)
        work.col_swap(1, right)
        block = work[:2, :2]
        if block.det() >= 0:
            raise ArithmeticError("unexpected exact inertia pivot")
        positive += 1
        negative += 1
        if work.rows == 2:
            break
        cross = work[2:, :2]
        work = (work[2:, 2:] - cross * block.inv() * cross.T).applyfunc(sp.cancel)
    return positive, negative, zero


def _load_factory_input() -> tuple[dict[str, Any], sp.Matrix]:
    if _sha256(FACTORY) != FACTORY_SHA256:
        raise ArithmeticError("the pinned H9/14 factory archive changed")
    factory = json.loads(FACTORY.read_text(encoding="utf-8"))
    row = next(item for item in factory["constructions"] if item["id"] == "A18_4_0_C01")
    raw_gram = row["section"]["restricted_mirrors"]["facet_gram"]
    if _json_sha256(raw_gram) != FACET_GRAM_SHA256:
        raise ArithmeticError("the pinned T2/C01 facet Gram changed")
    gram = sp.Matrix([[int(value) for value in values] for values in raw_gram])
    return row, gram


def _load_comparison_inputs() -> dict[str, Any]:
    if _sha256(PROMOTION) != PROMOTION_SHA256:
        raise ArithmeticError("the comparison-only C01 promotion archive changed")
    if _sha256(ANCILLARY) != ANCILLARY_SHA256:
        raise ArithmeticError("the comparison-only Alexeev ancillary changed")
    return json.loads(PROMOTION.read_text(encoding="utf-8"))


def _coxeter_data(gram: sp.Matrix) -> tuple[dict[int, list[tuple[int, int]]], nx.Graph]:
    edges = {3: [], 4: []}
    graph = nx.Graph()
    graph.add_nodes_from(range(gram.rows))
    histogram: Counter[int] = Counter()
    for left in range(gram.rows):
        for right in range(left + 1, gram.rows):
            if gram[left, right] == 0:
                order = 2
            else:
                normalized = Fraction(
                    4 * int(gram[left, right]) ** 2,
                    int(gram[left, left] * gram[right, right]),
                )
                order = {Fraction(1): 3, Fraction(2): 4}[normalized]
                edges[order].append((left, right))
            histogram[order] += 1
            if order != 2:
                graph.add_edge(left, right, order=order)
    if histogram != Counter({2: 75, 3: 8, 4: 8}):
        raise ArithmeticError("the normalized T2 relation inventory changed")
    return edges, graph


def _normal_scalings(
    edges: dict[int, list[tuple[int, int]]], graph: nx.Graph
) -> tuple[list[dict[str, Any]], tuple[tuple[int, ...], ...]]:
    parent = list(range(14))

    def find(value: int) -> int:
        if parent[value] != value:
            parent[value] = find(parent[value])
        return parent[value]

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for edge in edges[3]:
        union(*edge)
    components: list[int] = []
    for vertex in range(14):
        root = find(vertex)
        if root not in components:
            components.append(root)
    component_index = {root: index for index, root in enumerate(components)}
    ratio_edges = [
        (component_index[find(left)], component_index[find(right)])
        for left, right in edges[4]
    ]
    solutions: set[tuple[int, ...]] = set()
    for signs in itertools.product((-1, 1), repeat=len(ratio_edges)):
        adjacency: list[list[tuple[int, int]]] = [[] for _ in components]
        for (left, right), sign in zip(ratio_edges, signs, strict=True):
            adjacency[left].append((right, sign))
            adjacency[right].append((left, -sign))
        heights: list[int | None] = [None] * len(components)
        heights[0] = 0
        stack = [0]
        consistent = True
        while stack:
            left = stack.pop()
            for right, change in adjacency[left]:
                value = int(heights[left]) + change
                if heights[right] is None:
                    heights[right] = value
                    stack.append(right)
                elif heights[right] != value:
                    consistent = False
        if consistent and all(value is not None for value in heights):
            minimum = min(int(value) for value in heights)
            normalized = [int(value) - minimum for value in heights]
            solutions.add(
                tuple(
                    2 * 2 ** normalized[component_index[find(vertex)]]
                    for vertex in range(14)
                )
            )
    scalings = tuple(sorted(solutions))
    automorphisms = tuple(
        tuple(mapping[index] for index in range(14))
        for mapping in nx.algorithms.isomorphism.GraphMatcher(
            graph,
            graph,
            edge_match=lambda left, right: left["order"] == right["order"],
        ).isomorphisms_iter()
    )
    if len(automorphisms) != 48:
        raise ArithmeticError("the exact T2 diagram automorphism group changed")
    unseen = set(scalings)
    orbit_number: dict[tuple[int, ...], int] = {}
    orbit_sizes: list[int] = []
    while unseen:
        representative = min(unseen)
        orbit = {
            tuple(representative[permutation[index]] for index in range(14))
            for permutation in automorphisms
        } & set(scalings)
        number = len(orbit_sizes) + 1
        for scaling in orbit:
            orbit_number[scaling] = number
        orbit_sizes.append(len(orbit))
        unseen -= orbit
    if len(scalings) != 18 or sorted(orbit_sizes) != [1, 1, 2, 4, 4, 6]:
        raise ArithmeticError("the primitive normal-scaling census changed")
    rows = []
    for number, squares in enumerate(scalings, start=1):
        histogram = Counter(squares)
        rows.append(
            {
                "scaling_number": number,
                "normal_square_absolutes": list(squares),
                "normal_square_histogram": {
                    str(key): value for key, value in sorted(histogram.items())
                },
                "diagram_automorphism_orbit_number": orbit_number[squares],
                "two_sheet_minus_two_compatible": set(squares) <= {2, 4},
                "norm_eight_two_sheet_obstruction": 8 in squares,
            }
        )
    return rows, scalings


def _gram_for_scaling(
    squares: Sequence[int], edges: dict[int, list[tuple[int, int]]]
) -> sp.Matrix:
    gram = sp.diag(*[-value for value in squares])
    for left, right in edges[3]:
        if squares[left] != squares[right]:
            raise ArithmeticError("an m3 edge acquired unequal root lengths")
        gram[left, right] = gram[right, left] = squares[left] // 2
    for left, right in edges[4]:
        low, high = sorted((squares[left], squares[right]))
        if high != 2 * low:
            raise ArithmeticError("an m4 edge lost its crystallographic ratio")
        gram[left, right] = gram[right, left] = high // 2
    return gram


def _root_span(gram: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, list[int]]:
    indices: list[int] | None = None
    for candidate in itertools.combinations(range(14), 10):
        block = gram.extract(candidate, candidate)
        if block.det():
            indices = list(candidate)
            break
    if indices is None:
        raise ArithmeticError("the normalized Gram has no rank-ten principal block")
    span_gram = gram.extract(indices, indices)
    coordinates = span_gram.inv() * gram.extract(indices, range(14))
    if not all(value.q == 1 for value in coordinates):
        raise ArithmeticError("the fourteen roots do not have an integral root-span basis")
    coordinates = coordinates.applyfunc(int)
    if coordinates[:, indices] != sp.eye(10):
        raise ArithmeticError("the selected root-span coordinates changed")
    return span_gram, coordinates, indices


def _discriminant_vectors(gram: sp.Matrix) -> tuple[int, ...]:
    size = gram.rows
    row_masks = [
        sum((int(gram[row, column]) & 1) << column for column in range(size))
        for row in range(size)
    ]
    return tuple(
        value
        for value in range(1 << size)
        if all((value & mask).bit_count() % 2 == 0 for mask in row_masks)
    )


def _discriminant_q(gram: sp.Matrix, value: int) -> Fraction:
    vector = _bit_vector(value, gram.rows)
    numerator = int((vector.T * gram * vector)[0])
    return Fraction(numerator % 8, 4)


def _two_elementary_data(gram: sp.Matrix) -> dict[str, Any]:
    smith = _smith(gram)
    if any(value not in (1, 2) for value in smith):
        raise ArithmeticError("a fixed-lattice candidate is not 2-elementary")
    vectors = _discriminant_vectors(gram)
    length = int(math.log2(len(vectors)))
    q_histogram = Counter(str(_discriminant_q(gram, value)) for value in vectors)
    integral_values = []
    for value in vectors:
        vector = _bit_vector(value, gram.rows)
        numerator = int((vector.T * gram * vector)[0])
        integral_values.append(numerator % 4 == 0)
    return {
        "smith_diagonal_absolute": smith,
        "discriminant_length": length,
        "delta": 0 if all(integral_values) else 1,
        "q_mod_2_histogram": {
            key: value for key, value in sorted(q_histogram.items())
        },
    }


def _reflection_stable_isotropic_subspaces(
    gram: sp.Matrix, root_coordinates: sp.Matrix, root_squares: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    size = gram.rows
    discriminant = _discriminant_vectors(gram)
    isotropic: dict[int, bool] = {}
    for value in discriminant:
        vector = _bit_vector(value, size)
        isotropic[value] = int((vector.T * gram * vector)[0]) % 8 == 0
    root_bits = [_vector_bit(root_coordinates[:, index]) for index in range(14)]
    actions: list[tuple[int, int]] = []
    for index, square in enumerate(root_squares):
        if square == 2:
            continue
        if square != 4:
            raise ArithmeticError("fixed-lattice enumeration received a non-two-sheet root")
        functional = 0
        root = root_coordinates[:, index]
        for coordinate in range(size):
            basis = sp.eye(size)[:, coordinate]
            coefficient = int((basis.T * gram * root)[0]) // 2
            if coefficient % 2:
                functional |= 1 << coordinate
        actions.append((functional, root_bits[index]))

    def act(value: int, action: tuple[int, int]) -> int:
        functional, root = action
        return value ^ (root if (value & functional).bit_count() % 2 else 0)

    def closure(basis: tuple[int, ...]) -> tuple[int, ...]:
        while True:
            enlarged = _reduce_basis(
                basis + tuple(act(value, action) for value in basis for action in actions),
                size,
            )
            if enlarged == basis:
                return basis
            basis = enlarged

    zero: tuple[int, ...] = ()
    seen = {zero}
    queue = deque([zero])
    while queue:
        basis = queue.popleft()
        current = set(_span(basis))
        for value in discriminant:
            if value in current:
                continue
            candidate = closure(_reduce_basis(basis + (value,), size))
            if candidate in seen:
                continue
            if all(isotropic[element] for element in _span(candidate)):
                seen.add(candidate)
                queue.append(candidate)
    return tuple(sorted(seen, key=lambda basis: (len(basis), basis)))


def _overlattice(
    span_gram: sp.Matrix,
    root_coordinates: sp.Matrix,
    subgroup: Sequence[int],
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    generators = 2 * sp.eye(10)
    for value in subgroup:
        generators = generators.row_join(_bit_vector(value, 10))
    basis_in_span = hermite_normal_form(generators.applyfunc(int)) / 2
    gram = sp.simplify(basis_in_span.T * span_gram * basis_in_span)
    roots = sp.simplify(basis_in_span.inv() * root_coordinates)
    if not all(value.q == 1 for value in gram) or not all(value.q == 1 for value in roots):
        raise ArithmeticError("an even glue subgroup did not produce an integral overlattice")
    return basis_in_span, gram.applyfunc(int), roots.applyfunc(int)


def _reflection_integral(gram: sp.Matrix, root: sp.Matrix, square: int) -> bool:
    reflection = sp.eye(gram.rows) + sp.Rational(2, square) * root * (root.T * gram)
    return all(value.q == 1 for value in reflection)


def _long_space(
    folded_gram: sp.Matrix, roots: sp.Matrix
) -> tuple[list[int], list[int], list[int], list[list[int]], int]:
    long_indices = [index for index in range(14) if folded_gram[index, index] == -4]
    classes = [_vector_bit(roots[:, index]) for index in long_indices]
    basis_classes: list[int] = []
    basis_indices: list[int] = []
    for index, value in zip(long_indices, classes, strict=True):
        if len(_reduce_basis(basis_classes + [value], 10)) > len(basis_classes):
            basis_classes.append(value)
            basis_indices.append(index)
    polar = [
        [int(folded_gram[left, right] // 2) % 2 for right in basis_indices]
        for left in basis_indices
    ]
    radical_dimension = len(basis_indices) - _gf2_rank_rows(polar)
    return long_indices, classes, basis_indices, polar, radical_dimension


def _e8_data() -> tuple[sp.Matrix, dict[int, tuple[tuple[int, ...], ...]]]:
    gram = 2 * sp.eye(8)
    for left, right in ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7)):
        gram[left, right] = gram[right, left] = -1
    roots: set[tuple[int, ...]] = {(1, 0, 0, 0, 0, 0, 0, 0)}
    queue = deque(roots)
    while queue:
        root = sp.Matrix(queue.popleft())
        for index in range(8):
            image = list(root)
            image[index] -= int((root.T * gram[:, index])[0])
            candidate = tuple(int(value) for value in image)
            if candidate not in roots:
                roots.add(candidate)
                queue.append(candidate)
    by_class: dict[int, list[tuple[int, ...]]] = {}
    for root in roots:
        by_class.setdefault(_vector_bit(root), []).append(root)
    frozen = {key: tuple(sorted(values)) for key, values in by_class.items()}
    if not (
        gram.det() == 1
        and len(roots) == 240
        and len(frozen) == 120
        and all(len(values) == 2 and values[0] == tuple(-x for x in values[1]) for values in frozen.values())
    ):
        raise ArithmeticError("the exact E8 root/class census changed")
    return gram, frozen


def _e8_embedding(
    polar: Sequence[Sequence[int]],
    e8_gram: sp.Matrix,
    root_classes: dict[int, tuple[tuple[int, ...], ...]],
) -> tuple[int, ...] | None:
    dimension = len(polar)
    candidates = tuple(sorted(root_classes))

    def pairing(left: int, right: int) -> int:
        return int(
            (_bit_vector(left, 8).T * e8_gram * _bit_vector(right, 8))[0]
        ) % 2

    chosen = [candidates[0]]  # the exact E8 Weyl orbit above is all 120 root classes

    def search(index: int) -> tuple[int, ...] | None:
        if index == dimension:
            return tuple(chosen)
        for candidate in candidates:
            if len(_reduce_basis(chosen + [candidate], 8)) != len(chosen) + 1:
                continue
            if all(pairing(candidate, chosen[old]) == polar[index][old] for old in range(index)):
                chosen.append(candidate)
                answer = search(index + 1)
                if answer is not None:
                    return answer
                chosen.pop()
        return None

    return search(1)


def _content(vector: sp.Matrix) -> int:
    entries = [abs(int(value)) for value in vector if value]
    return math.gcd(*entries) if entries else 0


def _fixed_candidates(
    scaling_number: int,
    folded_gram: sp.Matrix,
    e8_gram: sp.Matrix,
    e8_root_classes: dict[int, tuple[tuple[int, ...], ...]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Enumerate all marked, reflection-stable, even fixed overlattices."""

    span_gram, span_roots, basis_indices = _root_span(folded_gram)
    squares = [-int(folded_gram[index, index]) for index in range(14)]
    subgroups = _reflection_stable_isotropic_subspaces(
        span_gram, span_roots, squares
    )
    public: list[dict[str, Any]] = []
    internal: list[dict[str, Any]] = []
    for number, subgroup in enumerate(subgroups, start=1):
        basis, gram, roots = _overlattice(span_gram, span_roots, subgroup)
        if not all(gram[index, index] % 2 == 0 for index in range(10)):
            raise ArithmeticError("an enumerated fixed overlattice is not even")
        if roots.T * gram * roots != folded_gram:
            raise ArithmeticError("a fixed overlattice changed the folded Gram")
        if not all(_content(roots[:, index]) == 1 for index in range(14)):
            raise ArithmeticError("a marked facet root became imprimitive")
        if not all(
            _reflection_integral(gram, roots[:, index], squares[index])
            for index in range(14)
        ):
            raise ArithmeticError("a marked reflection ceased to be integral")

        discr = _two_elementary_data(gram)
        long_indices, classes, long_basis_indices, polar, radical = _long_space(
            folded_gram, roots
        )
        repeated = []
        for left_position, left in enumerate(long_indices):
            for right_position in range(left_position):
                right = long_indices[right_position]
                if classes[left_position] == classes[right_position]:
                    repeated.append(
                        {
                            "left_facet_one_based": right + 1,
                            "right_facet_one_based": left + 1,
                            "folded_pairing": int(folded_gram[right, left]),
                            "class_hex": hex(classes[left_position]),
                        }
                    )
        distinct = not repeated
        dimension = len(long_basis_indices)
        symplectic_subspace_bound = radical <= 8 - dimension
        embedding = (
            _e8_embedding(polar, e8_gram, e8_root_classes)
            if distinct and symplectic_subspace_bound
            else None
        )
        minimal_rank18_discriminant_length = (
            discr["discriminant_length"] + 8 - 2 * dimension
        )
        weak_k3_length_gate = minimal_rank18_discriminant_length <= 4
        survives = distinct and embedding is not None and weak_k3_length_gate
        record = {
            "scaling_number": scaling_number,
            "candidate_number": number,
            "glue_dimension": len(subgroup),
            "glue_basis_little_endian_hex": [hex(value) for value in subgroup],
            "index_over_root_span": 1 << len(subgroup),
            "determinant": int(gram.det()),
            "smith_diagonal_absolute": discr["smith_diagonal_absolute"],
            "signature": list(_inertia(gram)[:2]),
            "discriminant_length": discr["discriminant_length"],
            "discriminant_delta": discr["delta"],
            "q_in_Q_mod_2Z_histogram": discr["q_mod_2_histogram"],
            "long_facet_indices_one_based": [value + 1 for value in long_indices],
            "long_class_span_dimension": dimension,
            "long_class_basis_facets_one_based": [
                value + 1 for value in long_basis_indices
            ],
            "long_class_polar_matrix_mod_two": polar,
            "long_class_radical_dimension": radical,
            "zero_long_class_facets_one_based": [
                index + 1
                for index, value in zip(long_indices, classes, strict=True)
                if value == 0
            ],
            "repeated_long_class_pairs": repeated,
            "minimal_rank18_discriminant_length": minimal_rank18_discriminant_length,
            "gates": {
                "long_classes_distinct": distinct,
                "symplectic_subspace_bound": symplectic_subspace_bound,
                "embeddable_as_nonsingular_E8_classes": embedding is not None,
                "rank4_K3_complement_length": weak_k3_length_gate,
            },
            "E8_class_obstruction": (
                "repeated_class"
                if not distinct
                else (
                    "radical_dimension_exceeds_8_minus_dimension"
                    if not symplectic_subspace_bound
                    else ("none" if embedding is not None else "embedding_search_failed")
                )
            ),
            "survives_minimal_two_sheet_weak_K3_gate": survives,
        }
        public.append(record)
        internal.append(
            {
                **record,
                "subgroup": tuple(subgroup),
                "basis_in_root_span": basis,
                "gram": gram,
                "roots": roots,
                "root_span_gram": span_gram,
                "root_span_roots": span_roots,
                "root_span_basis_indices": basis_indices,
                "long_classes": classes,
                "e8_embedding_basis_classes": embedding,
            }
        )
    return public, internal


def _integral(matrix: sp.Matrix, message: str) -> sp.Matrix:
    if not all(value.q == 1 for value in matrix):
        raise ArithmeticError(message)
    return matrix.applyfunc(int)


def _k3_glue(
    n_gram: sp.Matrix, involution: sp.Matrix
) -> dict[str, Any]:
    """Glue the compiled rank-18 lattice to U(2)^2 and extend the involution."""

    source_classes = _reduce_basis(_discriminant_vectors(n_gram), 18)
    if len(source_classes) != 4:
        raise ArithmeticError("the compiled rank-18 discriminant length is not four")
    t_gram = sp.Matrix(
        [[0, 2, 0, 0], [2, 0, 0, 0], [0, 0, 0, 2], [0, 0, 2, 0]]
    )

    def q_value(gram: sp.Matrix, value: int) -> Fraction:
        vector = _bit_vector(value, gram.rows)
        return Fraction(int((vector.T * gram * vector)[0]) % 8, 4)

    source_q = {
        coefficients: q_value(
            n_gram,
            next(
                value
                for value in _span(source_classes)
                if _relation(value, source_classes) == coefficients
            ),
        )
        for coefficients in range(16)
    }
    target_basis: tuple[int, ...] | None = None
    for candidate in itertools.permutations(range(1, 16), 4):
        if len(_reduce_basis(candidate, 4)) != 4:
            continue
        if all(
            q_value(
                t_gram,
                next(
                    value
                    for value in _span(candidate)
                    if _relation(value, candidate) == coefficients
                ),
            )
            == source_q[coefficients]
            for coefficients in range(16)
        ):
            target_basis = candidate
            break
    if target_basis is None:
        raise ArithmeticError("no discriminant anti-isometry to U(2)^2 was found")

    direct_gram = sp.diag(n_gram, t_gram)
    glue_columns = []
    for source, target in zip(source_classes, target_basis, strict=True):
        glue_columns.append(
            sp.Matrix.vstack(_bit_vector(source, 18), _bit_vector(target, 4))
        )
    generators = (2 * sp.eye(22)).row_join(sp.Matrix.hstack(*glue_columns))
    basis = hermite_normal_form(generators.applyfunc(int)) / 2
    k3_gram = _integral(
        sp.simplify(basis.T * direct_gram * basis),
        "the K3 graph glue is not integral",
    )
    direct_action = sp.diag(involution, sp.eye(4))
    action = _integral(
        sp.simplify(basis.inv() * direct_action * basis),
        "the involution did not extend over the K3 graph glue",
    )
    if not (
        all(k3_gram[index, index] % 2 == 0 for index in range(22))
        and k3_gram.det() == -1
        and _inertia(k3_gram) == (3, 19, 0)
        and action.T * k3_gram * action == k3_gram
        and action**2 == sp.eye(22)
    ):
        raise ArithmeticError("the glued K3 lattice/action checks failed")
    n_embedding = _integral(
        basis.inv() * sp.Matrix.vstack(sp.eye(18), sp.zeros(4, 18)),
        "the rank-18 lattice is not integrally embedded in the K3 glue",
    )
    t_embedding = _integral(
        basis.inv() * sp.Matrix.vstack(sp.zeros(18, 4), sp.eye(4)),
        "the rank-four complement is not integrally embedded in the K3 glue",
    )
    if _smith(n_embedding) != [1] * 18 or _smith(t_embedding) != [1] * 4:
        raise ArithmeticError("a K3 graph-glue summand is not primitive")
    symbol = sp.Symbol("lambda")
    expected_charpoly = (symbol - 1) ** 14 * (symbol + 1) ** 8
    if sp.factor(action.charpoly(symbol).as_expr()) != expected_charpoly:
        raise ArithmeticError("the extended involution characteristic polynomial changed")
    return {
        "complement_lattice": "U(2)+U(2)",
        "complement_gram": _matrix_rows(t_gram),
        "source_discriminant_basis_little_endian_hex": [
            hex(value) for value in source_classes
        ],
        "target_discriminant_basis_little_endian_hex": [
            hex(value) for value in target_basis
        ],
        "graph_glue_index": 16,
        "k3_lattice_rank_signature_determinant": [22, 3, 19, -1],
        "k3_gram": _matrix_rows(k3_gram),
        "k3_basis_columns_in_N_plus_T_tensor_Q": [
            [str(basis[row, column]) for column in range(22)]
            for row in range(22)
        ],
        "extended_involution": _matrix_rows(action),
        "extended_involution_characteristic_polynomial": "(lambda-1)^14*(lambda+1)^8",
        "N_and_T_embeddings_primitive": True,
    }


def _construct_survivor(
    candidate: dict[str, Any],
    folded_gram: sp.Matrix,
    e8_gram: sp.Matrix,
    e8_root_classes: dict[int, tuple[tuple[int, ...], ...]],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Construct the minimal rank-18 two-sheet unfolding and K3 embedding."""

    fixed_gram: sp.Matrix = candidate["gram"]
    fixed_roots: sp.Matrix = candidate["roots"]
    long_indices = [index for index in range(14) if folded_gram[index, index] == -4]
    long_basis_indices = [
        value - 1 for value in candidate["long_class_basis_facets_one_based"]
    ]
    fixed_basis_classes = [
        _vector_bit(fixed_roots[:, index]) for index in long_basis_indices
    ]
    target_basis_classes: tuple[int, ...] = candidate[
        "e8_embedding_basis_classes"
    ]
    if len(target_basis_classes) != len(fixed_basis_classes):
        raise ArithmeticError("the long-class E8 embedding changed dimension")

    anti_roots: dict[int, sp.Matrix] = {}
    target_classes: dict[int, int] = {}
    for index in long_indices:
        fixed_class = _vector_bit(fixed_roots[:, index])
        relation = _relation(fixed_class, fixed_basis_classes)
        target = 0
        for position, basis_class in enumerate(target_basis_classes):
            if (relation >> position) & 1:
                target ^= basis_class
        representatives = e8_root_classes.get(target)
        if representatives is None:
            raise ArithmeticError("a long class did not map to an E8 root class")
        anti_roots[index] = sp.Matrix(max(representatives))
        target_classes[index] = target

    anti_gram = -2 * e8_gram
    direct_gram = sp.diag(fixed_gram, anti_gram)
    direct_roots: list[sp.Matrix] = []
    root_orbits: list[list[int]] = []
    glue_columns: list[sp.Matrix] = []
    for facet in range(14):
        fixed = sp.Matrix.vstack(fixed_roots[:, facet], sp.zeros(8, 1))
        if folded_gram[facet, facet] == -2:
            direct_roots.append(fixed)
            root_orbits.append([len(direct_roots) - 1])
            continue
        anti = sp.Matrix.vstack(sp.zeros(10, 1), anti_roots[facet])
        plus = (fixed + anti) / 2
        minus = (fixed - anti) / 2
        direct_roots.extend((plus, minus))
        root_orbits.append([len(direct_roots) - 2, len(direct_roots) - 1])
        glue_columns.append(plus)

    generators = (2 * sp.eye(18)).row_join(
        sp.Matrix.hstack(*[(2 * vector).applyfunc(int) for vector in glue_columns])
    )
    basis = hermite_normal_form(generators) / 2
    n_gram = _integral(
        sp.simplify(basis.T * direct_gram * basis),
        "the rank-18 graph glue is not integral",
    )
    roots = _integral(
        sp.simplify(basis.inv() * sp.Matrix.hstack(*direct_roots)),
        "a lifted source root is not integral in the graph glue",
    )
    involution = _integral(
        sp.simplify(
            basis.inv() * sp.diag(sp.eye(10), -sp.eye(8)) * basis
        ),
        "the sheet involution is not integral",
    )
    fixed_embedding = _integral(
        basis.inv() * sp.Matrix.vstack(sp.eye(10), sp.zeros(8, 10)),
        "the fixed lattice does not embed integrally",
    )
    anti_embedding = _integral(
        basis.inv() * sp.Matrix.vstack(sp.zeros(10, 8), sp.eye(8)),
        "the coinvariant lattice does not embed integrally",
    )
    if not (
        n_gram.det() == -16
        and _smith(n_gram) == [1] * 14 + [2] * 4
        and _inertia(n_gram) == (1, 17, 0)
        and all(n_gram[index, index] % 2 == 0 for index in range(18))
        and involution**2 == sp.eye(18)
        and involution.T * n_gram * involution == n_gram
        and _smith(fixed_embedding) == [1] * 10
        and _smith(anti_embedding) == [1] * 8
    ):
        raise ArithmeticError("the rank-18 lattice/involution invariant package failed")
    if not (
        fixed_embedding.T * n_gram * fixed_embedding == fixed_gram
        and anti_embedding.T * n_gram * anti_embedding == anti_gram
        and involution * fixed_embedding == fixed_embedding
        and involution * anti_embedding == -anti_embedding
    ):
        raise ArithmeticError("the primitive fixed/coinvariant identifications failed")

    root_gram = roots.T * n_gram * roots
    off_diagonal = [
        int(root_gram[left, right])
        for left in range(root_gram.rows)
        for right in range(left)
    ]
    if not (
        all(root_gram[index, index] == -2 for index in range(root_gram.rows))
        and set(off_diagonal) <= {0, 1}
        and min(off_diagonal) >= 0
    ):
        raise ArithmeticError("the 24 source-root lift is not acute and simply laced")
    root_hnf = hermite_normal_form(roots)
    if root_hnf.shape != (18, 18) or abs(int(root_hnf.det())) != 1:
        raise ArithmeticError("the 24 lifted roots do not span the rank-18 lattice")
    folded_roots = sp.Matrix.hstack(
        *[
            roots[:, orbit[0]]
            if len(orbit) == 1
            else roots[:, orbit[0]] + roots[:, orbit[1]]
            for orbit in root_orbits
        ]
    )
    if folded_roots.T * n_gram * folded_roots != folded_gram:
        raise ArithmeticError("the lifted roots do not fold back to the input Gram")
    if not all(
        involution * roots[:, orbit[0]] == roots[:, orbit[-1]]
        for orbit in root_orbits
    ):
        raise ArithmeticError("the involution does not realize the sheet orbits")

    n_discr = _two_elementary_data(n_gram)
    if not (
        n_discr["discriminant_length"] == 4
        and n_discr["delta"] == 0
        and n_discr["q_mod_2_histogram"] == {"0": 10, "1": 6}
    ):
        raise ArithmeticError("the compiled rank-18 discriminant form changed")
    k3 = _k3_glue(n_gram, involution)
    source_root_groups = [
        [index + 1 for index in orbit] for orbit in root_orbits
    ]
    public = {
        "fixed_lattice": {
            "rank_signature": [10, 1, 9],
            "determinant": int(fixed_gram.det()),
            "smith_diagonal_absolute": _smith(fixed_gram),
            "two_elementary_invariants": [10, 8, 0],
            "isometry_class": "U+E8(-2)",
            "classification_reason": (
                "unique indefinite even 2-elementary lattice of signature (1,9) "
                "and invariants (r,a,delta)=(10,8,0)"
            ),
            "root_span_glue_basis_little_endian_hex": [
                hex(value) for value in candidate["subgroup"]
            ],
            "gram": _matrix_rows(fixed_gram),
            "facet_root_coordinates": _matrix_rows(fixed_roots),
        },
        "E8_two_sheet_data": {
            "long_class_span_dimension": len(fixed_basis_classes),
            "fixed_basis_facets_one_based": [value + 1 for value in long_basis_indices],
            "target_basis_classes_little_endian_hex": [
                f"{value:02x}" for value in target_basis_classes
            ],
            "all_target_classes_little_endian_hex": [
                f"{target_classes[index]:02x}" for index in long_indices
            ],
            "anti_E8_root_representatives_by_facet": {
                str(index + 1): [int(value) for value in anti_roots[index]]
                for index in long_indices
            },
            "graph_glue_index": 64,
            "absolute_anti_pairings_equal_folded_long_pairings": all(
                abs(int((anti_roots[left].T * anti_gram * anti_roots[right])[0]))
                == int(folded_gram[left, right])
                for left in long_indices
                for right in long_indices
                if left < right
            ),
        },
        "rank18_unfolding": {
            "rank_signature_determinant": [18, 1, 17, -16],
            "smith_diagonal_absolute": _smith(n_gram),
            "two_elementary_invariants": [18, 4, 0],
            "gram": _matrix_rows(n_gram),
            "basis_columns_in_fixed_plus_E8_minus2_tensor_Q": [
                [str(basis[row, column]) for column in range(18)]
                for row in range(18)
            ],
            "involution": _matrix_rows(involution),
            "involution_characteristic_polynomial": "(lambda-1)^10*(lambda+1)^8",
            "fixed_and_coinvariant_embeddings_primitive": True,
            "fixed_embedding_columns": _matrix_rows(fixed_embedding),
            "coinvariant_embedding_columns": _matrix_rows(anti_embedding),
            "coinvariant_isometry_class": "E8(-2)",
            "coinvariant_has_no_minus_two_vectors": True,
        },
        "source_roots_and_fold": {
            "source_root_count": roots.cols,
            "source_root_coordinates": _matrix_rows(roots),
            "source_root_orbits_one_based": source_root_groups,
            "one_sheet_orbit_count": sum(len(orbit) == 1 for orbit in root_orbits),
            "two_sheet_orbit_count": sum(len(orbit) == 2 for orbit in root_orbits),
            "source_roots_span_N_with_index_one": True,
            "off_diagonal_pairing_histogram": {
                str(key): value for key, value in sorted(Counter(off_diagonal).items())
            },
            "all_roots_square_minus_two_and_pairings_zero_or_one": True,
            "folded_root_coordinates": _matrix_rows(folded_roots),
            "folded_gram_equals_selected_scaling": True,
        },
        "weak_K3_embedding_and_involution_extension": k3,
    }
    internal = {
        "n_gram": n_gram,
        "roots": roots,
        "root_orbits": root_orbits,
        "involution": involution,
    }
    return public, internal


def _ancillary_matrix(key: tuple[int, int, int]) -> sp.Matrix:
    source = ANCILLARY.read_text(encoding="utf-8")
    marker = f"coxeter_matrix[({key[0]},{key[1]},{key[2]})] = matrix(["
    start = source.index(marker) + len(marker) - 1
    depth = 0
    end = start
    while end < len(source):
        if source[end] == "[":
            depth += 1
        elif source[end] == "]":
            depth -= 1
            if depth == 0:
                end += 1
                break
        end += 1
    rows = ast.literal_eval(source[start:end])
    matrix = sp.Matrix(rows)
    if matrix.shape != (48, 48) or matrix != matrix.T:
        raise ArithmeticError("the comparison-only (18,4,0) wall Gram changed")
    return matrix


def _marked_archive_bridge(
    row: dict[str, Any],
    promotion: dict[str, Any],
    compiled: dict[str, Any],
    input_gram: sp.Matrix,
) -> dict[str, Any]:
    """Identify the compiled marking with the archived C01 source datum."""

    wall_gram = _ancillary_matrix((18, 4, 0))
    archived_gram = wall_gram[:18, :18]
    wall_rows = _integral(
        wall_gram[:, :18] * archived_gram.inv(),
        "the 48 ancillary walls are not integral in the leading-wall basis",
    )
    if wall_rows * archived_gram * wall_rows.T != wall_gram:
        raise ArithmeticError("the ancillary 48-wall coordinates do not reconstruct")

    raw_action = row["section"]["case"]["representative_nonidentity_primal_actions"][0]
    if _json_sha256(raw_action) != ACTION_SHA256:
        raise ArithmeticError("the archived C01 action hash changed")
    archived_action = sp.Matrix(
        [[sp.Rational(value) for value in values] for values in raw_action]
    )
    promotion_case = next(
        case for case in promotion["cases"] if case["id"] == "A18_4_0_C01"
    )
    if not (
        promotion["all_checks_pass"] is True
        and promotion["source_file_sha256"]["factory"] == FACTORY_SHA256
        and promotion_case["ambient_NS"]["action_sha256"] == ACTION_SHA256
    ):
        raise ArithmeticError("the comparison-only strong-promotion chain changed")

    active_records = [
        record
        for record in row["section"]["restricted_mirrors"]["records"]
        if record["hyperbolic_status"] == "active_facet"
    ]
    archived_groups = [record["source_wall_indices"] for record in active_records]
    if len(archived_groups) != 14 or [len(group) for group in archived_groups].count(2) != 10:
        raise ArithmeticError("the active C01 source-root orbit inventory changed")
    folded_archived = sp.Matrix.hstack(
        *[
            sum((wall_rows[index - 1, :].T for index in group), sp.zeros(18, 1))
            for group in archived_groups
        ]
    )
    if folded_archived.T * archived_gram * folded_archived != input_gram:
        raise ArithmeticError("the archived source roots do not fold to the input Gram")

    n_gram: sp.Matrix = compiled["n_gram"]
    roots: sp.Matrix = compiled["roots"]
    root_orbits: list[list[int]] = compiled["root_orbits"]
    involution: sp.Matrix = compiled["involution"]
    pivot_columns = list(roots.rref()[1])
    if len(pivot_columns) != 18:
        raise ArithmeticError("the compiled roots lost rank")
    compiler_basis = roots[:, pivot_columns]
    integral_isometries: list[dict[str, Any]] = []
    for mask in range(1 << 10):
        source_columns: list[sp.Matrix] = []
        pair_number = 0
        for group in archived_groups:
            ordered = list(group)
            if len(group) == 2:
                if (mask >> pair_number) & 1:
                    ordered.reverse()
                pair_number += 1
            source_columns.extend(wall_rows[index - 1, :].T for index in ordered)
        source_roots = sp.Matrix.hstack(*source_columns)
        if source_roots[:, pivot_columns].det() == 0:
            continue
        mapping = sp.simplify(source_roots[:, pivot_columns] * compiler_basis.inv())
        if not all(value.q == 1 for value in mapping):
            continue
        mapping = mapping.applyfunc(int)
        if not (
            abs(int(mapping.det())) == 1
            and mapping * roots == source_roots
            and mapping.T * archived_gram * mapping == n_gram
            and mapping * involution == archived_action.T * mapping
        ):
            continue
        integral_isometries.append(
            {
                "sheet_swap_mask_little_endian_hex": hex(mask),
                "determinant": int(mapping.det()),
                "mapping": mapping,
            }
        )
    if len(integral_isometries) != 4:
        raise ArithmeticError("the exact marked compiler/archive bridge count changed")
    witness = integral_isometries[0]
    if witness["sheet_swap_mask_little_endian_hex"] != "0x0":
        raise ArithmeticError("the canonical marked bridge is no longer the unswapped one")
    mapping_rows = _matrix_rows(witness["mapping"])
    mapping_hash = _json_sha256(mapping_rows)
    expected_mapping_hash = (
        "15f683beded0f7366fe8be371298e5c17bf42238e36c2bec7833398acf6d6e85"
    )
    if mapping_hash != expected_mapping_hash:
        raise ArithmeticError("the canonical marked bridge matrix changed")
    return {
        "role": "comparison-only; not an input to the surface-free uniqueness theorem",
        "inputs": {
            "promotion_path": str(PROMOTION.relative_to(ROOT)),
            "promotion_sha256": PROMOTION_SHA256,
            "ancillary_path": str(ANCILLARY),
            "ancillary_sha256": ANCILLARY_SHA256,
            "archived_action_canonical_json_sha256": ACTION_SHA256,
        },
        "active_source_wall_orbits_one_based": archived_groups,
        "sheet_assignment_count_checked": 1 << 10,
        "integral_marked_isometry_count": len(integral_isometries),
        "integral_isometry_swap_masks_little_endian_hex": [
            record["sheet_swap_mask_little_endian_hex"]
            for record in integral_isometries
        ],
        "canonical_witness": {
            "sheet_swap_mask_little_endian_hex": witness[
                "sheet_swap_mask_little_endian_hex"
            ],
            "determinant": witness["determinant"],
            "compiler_to_archived_coordinate_matrix": mapping_rows,
            "matrix_canonical_json_sha256": mapping_hash,
            "maps_all_24_source_roots_exactly": True,
            "intertwines_involutions_exactly": True,
            "is_integral_unimodular_isometry": True,
        },
        "strong_promotion_chain": {
            "promotion_all_checks_pass": True,
            "promotion_pins_same_factory_sha256": True,
            "promotion_C01_action_sha256_equals_bridge_action": True,
            "datum_specific_existing_C01_promotion_reached_without_coordinate_gap": True,
        },
    }


def build_certificate() -> dict[str, Any]:
    """Replay the finite compiler first, then attach the comparison bridge."""

    row, input_gram = _load_factory_input()
    if not (
        input_gram.shape == (14, 14)
        and input_gram == input_gram.T
        and input_gram.rank() == 10
        and _inertia(input_gram) == (1, 9, 4)
    ):
        raise ArithmeticError("the surface-free input Gram shape/signature changed")
    edges, diagram = _coxeter_data(input_gram)
    scaling_rows, scalings = _normal_scalings(edges, diagram)
    compatible_internal: list[tuple[dict[str, Any], sp.Matrix]] = []
    for public, squares in zip(scaling_rows, scalings, strict=True):
        scaled_gram = _gram_for_scaling(squares, edges)
        if scaled_gram.rank() != 10 or _inertia(scaled_gram) != (1, 9, 4):
            raise ArithmeticError("a crystallographic scaling changed hyperbolic rank")
        public["scaled_gram_integer_canonical_json_sha256"] = _json_sha256(
            _matrix_rows(scaled_gram)
        )
        public["equals_pinned_input_gram_numerically"] = scaled_gram == input_gram
        if public["two_sheet_minus_two_compatible"]:
            compatible_internal.append((public, scaled_gram))
    if not (
        len(compatible_internal) == 2
        and sum(row["norm_eight_two_sheet_obstruction"] for row in scaling_rows) == 16
    ):
        raise ArithmeticError("the two-sheet scaling boundary changed")

    e8_gram, e8_root_classes = _e8_data()
    lanes: list[dict[str, Any]] = []
    all_candidates: list[dict[str, Any]] = []
    for scaling, folded_gram in compatible_internal:
        lane = "archived_orientation" if folded_gram == input_gram else "dual_orientation"
        public_candidates, internal_candidates = _fixed_candidates(
            scaling["scaling_number"], folded_gram, e8_gram, e8_root_classes
        )
        expected_count = 16 if lane == "archived_orientation" else 5
        expected_histogram = (
            {0: 1, 1: 7, 2: 7, 3: 1}
            if lane == "archived_orientation"
            else {0: 1, 1: 3, 2: 1}
        )
        if len(public_candidates) != expected_count or Counter(
            item["glue_dimension"] for item in public_candidates
        ) != expected_histogram:
            raise ArithmeticError(f"the complete {lane} fixed-overlattice census changed")
        for item in internal_candidates:
            item["lane"] = lane
            all_candidates.append(item)
        lanes.append(
            {
                "label": lane,
                "scaling_number": scaling["scaling_number"],
                "normal_square_absolutes": scaling["normal_square_absolutes"],
                "folded_gram": _matrix_rows(folded_gram),
                "root_span": {
                    "basis_facets_one_based": [
                        value + 1
                        for value in internal_candidates[0]["root_span_basis_indices"]
                    ],
                    "determinant": int(
                        internal_candidates[0]["root_span_gram"].det()
                    ),
                    "smith_diagonal_absolute": _smith(
                        internal_candidates[0]["root_span_gram"]
                    ),
                },
                "reflection_stable_even_overlattice_count": len(public_candidates),
                "glue_dimension_histogram": {
                    str(key): value for key, value in expected_histogram.items()
                },
                "candidates": public_candidates,
            }
        )
    survivors = [
        candidate
        for candidate in all_candidates
        if candidate["survives_minimal_two_sheet_weak_K3_gate"]
    ]
    if len(survivors) != 1:
        raise ArithmeticError("the minimal two-sheet weak-K3 survivor is not unique")
    survivor = survivors[0]
    if not (
        survivor["lane"] == "archived_orientation"
        and survivor["candidate_number"] == 8
        and survivor["subgroup"] == (0x28A,)
        and survivor["determinant"] == -256
        and survivor["smith_diagonal_absolute"] == [1, 1] + [2] * 8
        and survivor["discriminant_length"] == 8
        and survivor["discriminant_delta"] == 0
        and survivor["long_class_span_dimension"] == 6
    ):
        raise ArithmeticError("the unique fixed-lattice survivor invariants changed")
    archived_gram = next(
        gram
        for scaling, gram in compatible_internal
        if scaling["scaling_number"] == survivor["scaling_number"]
    )
    survivor_public, survivor_internal = _construct_survivor(
        survivor, archived_gram, e8_gram, e8_root_classes
    )

    # Only after the surface-free theorem is complete do we read the action,
    # ancillary wall matrix, and previously certified strong-promotion archive.
    promotion = _load_comparison_inputs()
    bridge = _marked_archive_bridge(
        row, promotion, survivor_internal, input_gram
    )
    literal = row["literal_vinberg"]
    if not (
        row["exact_cone"]["decision"]["classification"]
        == "finite_volume_noncompact"
        and literal["literal_vinberg_passes"] is True
        and literal["finite_endpoint_count"] == 90
        and literal["ideal_endpoint_count"] == 20
        and literal["bad_edge_count"] == 0
        and literal["bad_parabolic_neighborhood_count"] == 0
    ):
        raise ArithmeticError("the known top-level Vinberg comparison changed")

    archived_lane = next(lane for lane in lanes if lane["label"] == "archived_orientation")
    dual_lane = next(lane for lane in lanes if lane["label"] == "dual_orientation")
    archived_root_span = archived_lane["candidates"][0]
    dual_root_span = dual_lane["candidates"][0]
    if not (
        archived_root_span["long_class_span_dimension"] == 7
        and archived_root_span["long_class_radical_dimension"] == 3
        and archived_root_span["E8_class_obstruction"]
        == "radical_dimension_exceeds_8_minus_dimension"
        and dual_root_span["gates"]["embeddable_as_nonsingular_E8_classes"]
        and dual_root_span["minimal_rank18_discriminant_length"] > 4
    ):
        raise ArithmeticError("the advertised two primary obstruction witnesses changed")

    checks = {
        "surface_free_input_is_only_the_pinned_normalized_facet_Gram": True,
        "normalized_relation_census_is_75_m2_8_m3_8_m4": True,
        "diagram_automorphism_group_has_order_48": True,
        "all_18_primitive_even_crystallographic_scalings_enumerated": True,
        "exactly_two_one_or_two_sheet_scaling_lanes": True,
        "all_16_and_5_reflection_stable_even_overlattices_enumerated": True,
        "E8_class_radical_repetition_and_K3_length_gates_applied_separately": True,
        "unique_fixed_lattice_survivor_is_U_plus_E8_minus2": True,
        "minimal_rank18_two_sheet_unfolding_is_even_and_exact": True,
        "fixed_and_coinvariant_sublattices_are_primitive": True,
        "24_compiled_roots_are_minus_two_acute_and_fold_exactly": True,
        "weak_primitive_K3_embedding_and_involution_extension_are_explicit": True,
        "exact_marked_bridge_to_existing_C01_datum_is_closed": True,
        "known_top_level_literal_Vinberg_archive_is_finite_volume_noncompact": True,
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-08",
        "surface_free_logical_input": {
            "path": str(FACTORY.relative_to(ROOT)),
            "construction_id": "A18_4_0_C01",
            "field": "section.restricted_mirrors.facet_gram",
            "factory_file_sha256": FACTORY_SHA256,
            "facet_gram_string_canonical_json_sha256": FACET_GRAM_SHA256,
            "facet_gram": row["section"]["restricted_mirrors"]["facet_gram"],
            "rank_and_inertia": [10, 1, 9, 4],
            "normalized_coxeter_relation_counts": {"m2": 75, "m3": 8, "m4": 8},
        },
        "finite_surface_free_compiler": {
            "scaling_ray_count": len(scaling_rows),
            "diagram_automorphism_group_order": 48,
            "diagram_automorphism_orbit_sizes": sorted(
                Counter(row["diagram_automorphism_orbit_number"] for row in scaling_rows).values()
            ),
            "two_sheet_compatible_scaling_count": len(compatible_internal),
            "norm_eight_obstructed_scaling_count": sum(
                row["norm_eight_two_sheet_obstruction"] for row in scaling_rows
            ),
            "scalings": scaling_rows,
            "lanes": lanes,
            "total_fixed_overlattice_candidate_count": len(all_candidates),
            "unique_survivor_count": len(survivors),
            "theorem": (
                "Among all primitive minimal even crystallographic normal scalings "
                "of the pinned Gram, all reflection-stable even fixed overlattices, "
                "and minimal E8(-2) one-/two-sheet unfoldings satisfying the weak "
                "rank-four K3-complement length gate, exactly one marked candidate "
                "survives; its fixed lattice is U+E8(-2)."
            ),
        },
        "unique_survivor": {
            "lane": survivor["lane"],
            "scaling_number": survivor["scaling_number"],
            "candidate_number": survivor["candidate_number"],
            **survivor_public,
        },
        "comparison_only_marked_bridge": bridge,
        "known_geometry_from_archives": {
            "role": "validates this one known folded datum; not used in the compiler uniqueness proof",
            "exact_cone_decision": row["exact_cone"]["decision"],
            "literal_vinberg": literal,
            "classification": "finite_volume_noncompact",
        },
        "claim_boundary": {
            "proved_grade": (
                "IC0 finite lattice-precompiler uniqueness, exact rank-18 fold-back, "
                "weak K3 lattice embedding, involution extension, and marked bridge "
                "to the existing datum-specific promotion"
            ),
            "new_hyperbolic_polytope_claim": False,
            "class_wide_surface_correspondence_claim": False,
            "class_wide_IC2_claim": False,
            "class_wide_IC4_claim": False,
        },
        "open_gates": {
            "IC2_surface_geometry": [
                "effectivity and common ample orientation",
                "irreducibility, or an explicit reducible-wall interpretation",
                "extremality and complete-root/hidden-wall separation",
            ],
            "period_rootlessness": (
                "E8(-2) is rootless and generic periods in U(2)^2 exist, but this "
                "pilot makes no new class-wide Torelli or moduli-stratum theorem"
            ),
            "IC4_polytope_geometry": [
                "reverse extraction for every retained surface stratum",
                "exact redundancy removal",
                "complete Vinberg elliptic-edge endpoint and parabolic-neighborhood audit",
            ],
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("H9/14 T2 surface-unfolding archive differs from exact replay")
        print("archive_check=PASS")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("surface_free_scalings=18")
    print("stable_even_fixed_overlattices=16,5")
    print("unique_candidate=U+E8(-2)")
    print("marked_archive_isometries=4")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
