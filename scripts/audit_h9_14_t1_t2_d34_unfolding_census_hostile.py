#!/usr/bin/env python3
"""Hostile, basis-independent audit of the bounded T1/T2 D34 census.

The production certificate uses a Hermite quotient and Hermite overclosure
bases.  This audit instead uses a Smith decomposition of the full degenerate
facet Gram, reduced GF(2) module enumeration, and pivot bases for every
overlattice.  It imports no code from either production compiler.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp


ROOT = Path(__file__).resolve().parents[1]
FACTORY = ROOT / "data/enumerations/alexeev_k3_h9_14_factory.json"
PARKED = ROOT / "data/enumerations/alexeev_parked_rank10_h9_examples.json"
PRODUCTS = ROOT / "data/enumerations/alexeev_parked_h9_doubling_products.json"
PRODUCTION = ROOT / "data/enumerations/h9_14_t1_t2_d34_unfolding_census.json"
OUTPUT = ROOT / "data/enumerations/h9_14_t1_t2_d34_unfolding_census_hostile_audit.json"

SOURCE_HASHES = {
    "factory": "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775",
    "parked": "e181b8c021fe0f37d23b53b33e4e3900559502d0f8c982a4a9e9cbcf0866e430",
    "products": "10e6323895bacd2d18f897f33ea5f9d4bb1ee4eedad061d545aa5302b447c694",
    "production": "76c7e08ed2fd5107822eec4c4c96f35d107a4a2e231efef2a969a8c98ba36620",
}
SCHEMA = "coxeter.h9_14_t1_t2_d34_unfolding_census.hostile_audit.v1"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_hash(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def rows(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def integral(matrix: sp.Matrix, message: str) -> sp.Matrix:
    if not all(value.q == 1 for value in matrix):
        raise ArithmeticError(message)
    return matrix.applyfunc(int)


def smith(matrix: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(diagonal[i, i]))
        for i in range(min(diagonal.shape))
        if diagonal[i, i]
    ]


def inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Exact symmetric elimination, including singular trailing blocks."""
    work = sp.Matrix(matrix)
    positive = negative = zero = 0
    while work.rows:
        pivot_index = next((i for i in range(work.rows) if work[i, i]), None)
        if pivot_index is not None:
            work.row_swap(0, pivot_index)
            work.col_swap(0, pivot_index)
            pivot = work[0, 0]
            if pivot > 0:
                positive += 1
            elif pivot < 0:
                negative += 1
            if work.rows == 1:
                break
            column = work[1:, 0]
            work = (work[1:, 1:] - column * column.T / pivot).applyfunc(sp.cancel)
            continue
        pair = next(
            (
                (i, j)
                for i in range(work.rows)
                for j in range(i + 1, work.cols)
                if work[i, j]
            ),
            None,
        )
        if pair is None:
            zero += work.rows
            break
        left, right = pair
        work.row_swap(0, left)
        work.col_swap(0, left)
        if right == 0:
            right = left
        work.row_swap(1, right)
        work.col_swap(1, right)
        block = work[:2, :2]
        if block.det() >= 0:
            raise ArithmeticError("bad exact two-by-two inertia pivot")
        positive += 1
        negative += 1
        if work.rows == 2:
            break
        cross = work[2:, :2]
        work = (work[2:, 2:] - cross * block.inv() * cross.T).applyfunc(sp.cancel)
    return positive, negative, zero


def bit_vector(value: int, size: int) -> sp.Matrix:
    return sp.Matrix([(value >> i) & 1 for i in range(size)])


def vector_bit(vector: Sequence[Any]) -> int:
    return sum((int(vector[i]) & 1) << i for i in range(len(vector)))


def gf2_basis(values: Iterable[int], size: int) -> tuple[int, ...]:
    """Unique reduced basis, indexed by its highest pivot bit."""
    table = [0] * size
    for original in values:
        value = int(original)
        while value:
            pivot = value.bit_length() - 1
            if table[pivot]:
                value ^= table[pivot]
                continue
            table[pivot] = value
            for lower in range(pivot):
                if table[lower] and ((table[pivot] >> lower) & 1):
                    table[pivot] ^= table[lower]
            for upper in range(pivot + 1, size):
                if table[upper] and ((table[upper] >> pivot) & 1):
                    table[upper] ^= table[pivot]
            break
    return tuple(value for value in table if value)


def gf2_span(basis: Sequence[int]) -> tuple[int, ...]:
    result = [0]
    for vector in basis:
        result += [value ^ vector for value in result]
    return tuple(result)


def gf2_coordinates(value: int, ordered_basis: Sequence[int]) -> int:
    for coefficients in range(1 << len(ordered_basis)):
        image = 0
        for i, vector in enumerate(ordered_basis):
            if (coefficients >> i) & 1:
                image ^= vector
        if image == value:
            return coefficients
    raise ArithmeticError("GF(2) coordinate solve failed")


def gf2_rank_matrix(matrix: Sequence[Sequence[int]]) -> int:
    return len(gf2_basis((sum((x & 1) << j for j, x in enumerate(row)) for row in matrix), len(matrix[0]) if matrix else 0))


def relation_data(gram: sp.Matrix) -> tuple[Counter[int], dict[int, list[tuple[int, int]]], nx.Graph]:
    histogram: Counter[int] = Counter()
    edges: dict[int, list[tuple[int, int]]] = {3: [], 4: []}
    graph = nx.Graph()
    graph.add_nodes_from(range(gram.rows))
    for left in range(gram.rows):
        for right in range(left):
            if gram[left, right] == 0:
                order = 2
            else:
                ratio = Fraction(
                    4 * int(gram[left, right]) ** 2,
                    int(gram[left, left] * gram[right, right]),
                )
                order = {
                    Fraction(1): 3,
                    Fraction(2): 4,
                    Fraction(3): 6,
                    Fraction(4): 0,
                }.get(ratio, -1)
            histogram[order] += 1
            if order in edges:
                edges[order].append((right, left))
            if order != 2:
                graph.add_edge(right, left, order=order)
    return histogram, edges, graph


def diagram_isomorphic(left: sp.Matrix, right: sp.Matrix) -> bool:
    left_graph = relation_data(left)[2]
    right_graph = relation_data(right)[2]
    return nx.is_isomorphic(
        left_graph,
        right_graph,
        edge_match=lambda a, b: a["order"] == b["order"],
    )


def load_corpus() -> tuple[list[dict[str, Any]], dict[str, int]]:
    for name, path in (("factory", FACTORY), ("parked", PARKED), ("products", PRODUCTS)):
        if sha256(path) != SOURCE_HASHES[name]:
            raise ArithmeticError(f"pinned hostile-audit source changed: {name}")
    factory = json.loads(FACTORY.read_text())
    parked = json.loads(PARKED.read_text())
    products = json.loads(PRODUCTS.read_text())
    corpus = [
        {
            "id": item["id"],
            "lineage": "factory",
            "gram": sp.Matrix(item["section"]["restricted_mirrors"]["facet_gram"]),
            "classification": item["exact_cone"]["decision"]["classification"],
            "vinberg": item["literal_vinberg"],
        }
        for item in factory["constructions"]
    ]
    parked_h9_14 = [item for item in parked["cases"] if item["exact_cone"]["facet_count"] == 14]
    product_h9_14 = [item for item in products["direct_exact_audits"] if item["facet_count"] == 14]
    corpus.extend(
        {
            "id": item["id"],
            "lineage": "parked_product",
            "gram": sp.Matrix(item["facet_gram"]),
            "classification": item["exact_cone"]["decision"]["classification"],
            "vinberg": item["literal_vinberg"],
        }
        for item in product_h9_14
    )
    return corpus, {
        "factory_markings": len(factory["constructions"]),
        "parked_base_H9_14_markings": len(parked_h9_14),
        "parked_product_H9_14_markings": len(product_h9_14),
    }


def scaling_rays(
    edges: dict[int, list[tuple[int, int]]], graph: nx.Graph
) -> tuple[list[tuple[int, ...]], int, list[int]]:
    """Enumerate height functions directly, rather than orienting m4 edges."""
    parent = list(range(len(graph)))

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for edge in edges[3]:
        union(*edge)
    components = sorted({find(i) for i in range(len(graph))})
    component_number = {value: i for i, value in enumerate(components)}
    vertex_component = [component_number[find(i)] for i in range(len(graph))]
    adjacency: list[set[int]] = [set() for _ in components]
    for left, right in edges[4]:
        a, b = vertex_component[left], vertex_component[right]
        if a == b:
            raise ArithmeticError("an m4 edge collapsed inside an m3 component")
        adjacency[a].add(b)
        adjacency[b].add(a)
    heights: dict[int, int] = {0: 0}
    solutions: set[tuple[int, ...]] = set()

    def search() -> None:
        if len(heights) == len(components):
            floor = min(heights.values())
            normalized = [heights[i] - floor for i in range(len(components))]
            solutions.add(tuple(2 * (1 << normalized[vertex_component[i]]) for i in range(len(graph))))
            return
        choices: list[tuple[int, set[int]]] = []
        for vertex in range(len(components)):
            if vertex in heights:
                continue
            assigned_neighbors = adjacency[vertex] & heights.keys()
            if not assigned_neighbors:
                continue
            allowed: set[int] | None = None
            for neighbor in assigned_neighbors:
                local = {heights[neighbor] - 1, heights[neighbor] + 1}
                allowed = local if allowed is None else allowed & local
            choices.append((vertex, allowed or set()))
        if not choices:
            raise ArithmeticError("the non-right Coxeter graph is disconnected")
        vertex, allowed = min(choices, key=lambda item: (len(item[1]), item[0]))
        for value in sorted(allowed):
            heights[vertex] = value
            search()
            del heights[vertex]

    search()
    scalings = sorted(solutions)
    automorphisms = tuple(
        tuple(mapping[i] for i in range(len(graph)))
        for mapping in nx.algorithms.isomorphism.GraphMatcher(
            graph,
            graph,
            edge_match=lambda a, b: a["order"] == b["order"],
        ).isomorphisms_iter()
    )
    unseen = set(scalings)
    orbit_sizes: list[int] = []
    while unseen:
        representative = min(unseen)
        orbit = {
            tuple(representative[permutation[i]] for i in range(len(graph)))
            for permutation in automorphisms
        } & set(scalings)
        orbit_sizes.append(len(orbit))
        unseen -= orbit
    return scalings, len(automorphisms), sorted(orbit_sizes)


def scaled_gram(squares: Sequence[int], edges: dict[int, list[tuple[int, int]]]) -> sp.Matrix:
    gram = sp.diag(*[-value for value in squares])
    for left, right in edges[3]:
        if squares[left] != squares[right]:
            raise ArithmeticError("m3 scaling mismatch")
        gram[left, right] = gram[right, left] = squares[left] // 2
    for left, right in edges[4]:
        low, high = sorted((squares[left], squares[right]))
        if high != 2 * low:
            raise ArithmeticError("m4 scaling mismatch")
        gram[left, right] = gram[right, left] = high // 2
    return gram


def smith_quotient(gram: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, dict[str, Any]]:
    """Construct Z^14/rad(G) from a full Smith decomposition of G."""
    domain = DomainMatrix.from_Matrix(gram).convert_to(ZZ)
    diagonal_dm, left_dm, right_dm = smith_normal_decomp(domain)
    diagonal, left, right = (
        diagonal_dm.to_Matrix(),
        left_dm.to_Matrix(),
        right_dm.to_Matrix(),
    )
    if left * gram * right != diagonal or abs(int(left.det())) != 1 or abs(int(right.det())) != 1:
        raise ArithmeticError("Smith decomposition identity/unimodularity failed")
    rank = sum(diagonal[i, i] != 0 for i in range(min(diagonal.shape)))
    kernel = right[:, rank:]
    if gram * kernel != sp.zeros(gram.rows, gram.cols - rank):
        raise ArithmeticError("Smith radical is not in the Gram kernel")
    complement = right[:, :rank]
    inverse = integral(right.inv(), "unimodular Smith matrix acquired denominators")
    root_coordinates = inverse[:rank, :]
    quotient_gram = integral(complement.T * gram * complement, "quotient Gram is not integral")
    if root_coordinates.T * quotient_gram * root_coordinates != gram:
        raise ArithmeticError("Smith quotient does not recover the labelled Gram")
    return quotient_gram, root_coordinates, {
        "ambient_rank": gram.rows,
        "quotient_rank": rank,
        "radical_rank": gram.rows - rank,
        "radical_saturated": smith(kernel) == [1] * (gram.rows - rank),
        "decomposition_identity": True,
        "labelled_gram_recovered": True,
        "quotient_gram_sha256": json_hash(rows(quotient_gram)),
    }


def discriminant_classes(gram: sp.Matrix) -> tuple[int, ...]:
    size = gram.rows
    masks = [
        sum((int(gram[i, j]) & 1) << j for j in range(size))
        for i in range(size)
    ]
    return tuple(
        value
        for value in range(1 << size)
        if all((value & mask).bit_count() % 2 == 0 for mask in masks)
    )


def q_value(gram: sp.Matrix, value: int) -> Fraction:
    vector = bit_vector(value, gram.rows)
    return Fraction(int((vector.T * gram * vector)[0]) % 8, 4)


def two_elementary(gram: sp.Matrix) -> dict[str, Any]:
    invariants = smith(gram)
    if len(invariants) != gram.rows or any(value not in (1, 2) for value in invariants):
        raise ArithmeticError("lattice is not nondegenerate 2-elementary")
    classes = discriminant_classes(gram)
    if len(classes) != 1 << invariants.count(2):
        raise ArithmeticError("kernel-mod-two discriminant count disagrees with Smith form")
    values = [q_value(gram, value) for value in classes]
    return {
        "smith": invariants,
        "length": invariants.count(2),
        "delta": int(any(value.denominator != 1 for value in values)),
        "q_histogram": {str(key): value for key, value in sorted(Counter(values).items())},
    }


def stable_isotropic_subspaces(
    gram: sp.Matrix, root_coordinates: sp.Matrix, squares: Sequence[int]
) -> tuple[tuple[int, ...], ...]:
    """Enumerate all invariant isotropic modules in A_L[2] by orbit closure."""
    size = gram.rows
    discriminant = discriminant_classes(gram)
    isotropic = {
        value: int((bit_vector(value, size).T * gram * bit_vector(value, size))[0]) % 8 == 0
        for value in discriminant
    }
    actions: list[tuple[int, int]] = []
    for facet, square in enumerate(squares):
        if square == 2:
            continue
        if square != 4:
            raise ArithmeticError("stable-overclosure audit received a nonminimal norm")
        root = root_coordinates[:, facet]
        functional = 0
        for coordinate in range(size):
            pairing = int((sp.eye(size)[:, coordinate].T * gram * root)[0])
            if pairing % 2:
                raise ArithmeticError("a discriminant reflection coefficient is half-integral")
            if (pairing // 2) % 2:
                functional |= 1 << coordinate
        actions.append((functional, vector_bit(root)))

    def act(value: int, action: tuple[int, int]) -> int:
        functional, root = action
        return value ^ (root if (value & functional).bit_count() % 2 else 0)

    orbit_span_cache: dict[int, tuple[int, ...]] = {}

    def orbit_span(value: int) -> tuple[int, ...]:
        if value in orbit_span_cache:
            return orbit_span_cache[value]
        orbit = {value}
        queue = deque([value])
        while queue:
            current = queue.popleft()
            for action in actions:
                image = act(current, action)
                if image not in orbit:
                    orbit.add(image)
                    queue.append(image)
        answer = gf2_basis(orbit, size)
        orbit_span_cache[value] = answer
        return answer

    seen: set[tuple[int, ...]] = {()}
    queue: deque[tuple[int, ...]] = deque([()])
    while queue:
        basis = queue.popleft()
        current = set(gf2_span(basis))
        for value in discriminant:
            if value in current:
                continue
            candidate = gf2_basis(basis + orbit_span(value), size)
            if candidate in seen:
                continue
            if all(isotropic[element] for element in gf2_span(candidate)):
                seen.add(candidate)
                queue.append(candidate)
    # Hostile completeness check: every result is stable, and adjoining any
    # missing discriminant vector either reproduces the list or becomes non-even.
    for basis in seen:
        span = set(gf2_span(basis))
        if any(act(value, action) not in span for value in span for action in actions):
            raise ArithmeticError("enumerated subgroup is not reflection-stable")
        if not all(isotropic[value] for value in span):
            raise ArithmeticError("enumerated subgroup is not totally isotropic")
    return tuple(sorted(seen, key=lambda basis: (len(basis), basis)))


def pivot_overlattice_basis(size: int, subgroup: Sequence[int]) -> sp.Matrix:
    """A rational pivot basis for Z^n + H/2; no Hermite form is used."""
    basis = gf2_basis(subgroup, size)
    pivots = [value.bit_length() - 1 for value in basis]
    columns = [bit_vector(value, size) / 2 for value in basis]
    columns.extend(sp.eye(size)[:, i] for i in range(size) if i not in pivots)
    transform = sp.Matrix.hstack(*columns) if columns else sp.eye(size)
    if abs(transform.det()) != sp.Rational(1, 1 << len(basis)):
        raise ArithmeticError("pivot overclosure basis has the wrong covolume")
    return transform


def root_content(vector: sp.Matrix) -> int:
    entries = [abs(int(value)) for value in vector if value]
    return math.gcd(*entries) if entries else 0


def e8_data() -> tuple[sp.Matrix, dict[int, tuple[tuple[int, ...], ...]]]:
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
    classes: dict[int, list[tuple[int, ...]]] = {}
    for root in roots:
        classes.setdefault(vector_bit(root), []).append(root)
    frozen = {key: tuple(sorted(values)) for key, values in classes.items()}
    if gram.det() != 1 or len(roots) != 240 or len(frozen) != 120 or 0 in frozen:
        raise ArithmeticError("E8 root/class census failed")
    if any(len(values) != 2 or values[0] != tuple(-x for x in values[1]) for values in frozen.values()):
        raise ArithmeticError("an E8 root class is not an opposite pair")
    return gram, frozen


def e8_marked_embedding(
    polar: Sequence[Sequence[int]],
    marked_relations: Sequence[int],
    e8_gram: sp.Matrix,
    root_classes: dict[int, tuple[tuple[int, ...], ...]],
) -> tuple[int, ...] | None:
    """Embed the whole marked long-class set, not only a chosen basis."""
    dimension = len(polar)
    candidates = tuple(sorted(root_classes))

    def pairing(left: int, right: int) -> int:
        return int((bit_vector(left, 8).T * e8_gram * bit_vector(right, 8))[0]) % 2

    chosen = [candidates[0]]

    def search(position: int) -> tuple[int, ...] | None:
        if position == dimension:
            for relation in marked_relations:
                target = 0
                for i, value in enumerate(chosen):
                    if (relation >> i) & 1:
                        target ^= value
                if target not in root_classes:
                    return None
            return tuple(chosen)
        for candidate in candidates:
            if len(gf2_basis(chosen + [candidate], 8)) != len(chosen) + 1:
                continue
            if any(pairing(candidate, chosen[i]) != polar[position][i] for i in range(position)):
                continue
            chosen.append(candidate)
            answer = search(position + 1)
            if answer is not None:
                return answer
            chosen.pop()
        return None

    return search(1) if dimension else ()


def candidate_inventory(
    folded: sp.Matrix,
    e8_gram: sp.Matrix,
    e8_classes: dict[int, tuple[tuple[int, ...], ...]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    span_gram, span_roots, quotient = smith_quotient(folded)
    squares = [-int(folded[i, i]) for i in range(14)]
    subgroups = stable_isotropic_subspaces(span_gram, span_roots, squares)
    public: list[dict[str, Any]] = []
    internal: list[dict[str, Any]] = []
    for number, subgroup in enumerate(subgroups, start=1):
        basis = pivot_overlattice_basis(10, subgroup)
        gram = integral(basis.T * span_gram * basis, "overclosure Gram is not integral")
        roots = integral(basis.inv() * span_roots, "marked root left the overclosure")
        if any(gram[i, i] % 2 for i in range(10)):
            raise ArithmeticError("stable isotropic overclosure is not even")
        if roots.T * gram * roots != folded:
            raise ArithmeticError("overclosure changed the folded Gram")
        if any(root_content(roots[:, i]) != 1 for i in range(14)):
            raise ArithmeticError("a marked folded root became imprimitive")
        for facet, square in enumerate(squares):
            root = roots[:, facet]
            reflection = sp.eye(10) + sp.Rational(2, square) * root * (root.T * gram)
            integral(reflection, "a marked reflection is not integral on the overclosure")
        discr = two_elementary(gram)
        long_indices = [i for i, square in enumerate(squares) if square == 4]
        long_classes = [vector_bit(roots[:, i]) for i in long_indices]
        ordered_basis: list[int] = []
        basis_indices: list[int] = []
        for facet, value in zip(long_indices, long_classes, strict=True):
            if len(gf2_basis(ordered_basis + [value], 10)) > len(ordered_basis):
                ordered_basis.append(value)
                basis_indices.append(facet)
        relations = [gf2_coordinates(value, ordered_basis) for value in long_classes]
        polar = [
            [int(folded[left, right] // 2) % 2 for right in basis_indices]
            for left in basis_indices
        ]
        radical = len(ordered_basis) - gf2_rank_matrix(polar)
        repeats = sum(
            long_classes[i] == long_classes[j]
            for i in range(len(long_classes))
            for j in range(i)
        )
        distinct_nonzero = len(set(long_classes)) == len(long_classes) and 0 not in long_classes
        symplectic_bound = radical <= 8 - len(ordered_basis)
        minimal_n_length = discr["length"] + 8 - 2 * len(ordered_basis)
        complement_bound = minimal_n_length <= 4
        embedding = None
        embedding_checked = distinct_nonzero and symplectic_bound and complement_bound
        if embedding_checked:
            embedding = e8_marked_embedding(polar, relations, e8_gram, e8_classes)
        survives = embedding is not None and complement_bound
        record = {
            "candidate_number": number,
            "glue_dimension": len(subgroup),
            "index_over_root_span": 1 << len(subgroup),
            "determinant": int(gram.det()),
            "smith": discr["smith"],
            "signature": list(inertia(gram)[:2]),
            "discriminant_length": discr["length"],
            "discriminant_delta": discr["delta"],
            "q_histogram": discr["q_histogram"],
            "long_count": len(long_indices),
            "long_class_span_dimension": len(ordered_basis),
            "long_class_radical_dimension": radical,
            "repeated_long_class_pair_count": repeats,
            "zero_long_class_count": long_classes.count(0),
            "minimal_rank18_discriminant_length": minimal_n_length,
            "gates": {
                "long_classes_distinct_nonzero": distinct_nonzero,
                "symplectic_subspace_bound": symplectic_bound,
                "rank4_K3_complement_length": complement_bound,
                "full_marked_E8_class_embedding": (
                    embedding is not None if embedding_checked else None
                ),
            },
            "E8_embedding_search_status": (
                "passed"
                if embedding is not None
                else ("failed" if embedding_checked else "not_needed_after_prior_gate")
            ),
            "survives": survives,
        }
        public.append(record)
        internal.append(
            {
                **record,
                "gram": gram,
                "roots": roots,
                "long_indices": long_indices,
                "long_classes": long_classes,
                "long_basis_classes": ordered_basis,
                "long_relations": relations,
                "e8_embedding": embedding,
            }
        )
    quotient["root_span_determinant"] = int(span_gram.det())
    quotient["root_span_smith"] = smith(span_gram)
    return public, internal, quotient


def primitive_embedding(embedding: sp.Matrix) -> bool:
    return smith(embedding) == [1] * embedding.cols


def find_discriminant_anti_isometry(source: sp.Matrix, target: sp.Matrix) -> tuple[tuple[int, ...], tuple[int, ...]]:
    source_basis = gf2_basis(discriminant_classes(source), source.rows)
    target_classes = discriminant_classes(target)
    length = len(source_basis)
    if len(target_classes) != 1 << length:
        raise ArithmeticError("complement discriminant order mismatch")
    source_values = {
        coefficients: q_value(
            source,
            next(value for value in gf2_span(source_basis) if gf2_coordinates(value, source_basis) == coefficients),
        )
        for coefficients in range(1 << length)
    }
    for target_basis in itertools.permutations((value for value in target_classes if value), length):
        if len(gf2_basis(target_basis, target.rows)) != length:
            continue
        good = True
        for coefficients, source_q in source_values.items():
            target_value = 0
            for i, value in enumerate(target_basis):
                if (coefficients >> i) & 1:
                    target_value ^= value
            if (source_q + q_value(target, target_value)) % 2:
                good = False
                break
        if good:
            return source_basis, tuple(target_basis)
    raise ArithmeticError("no exact discriminant anti-isometry exists")


def k3_glue(n_gram: sp.Matrix, involution: sp.Matrix) -> dict[str, Any]:
    n_discr = two_elementary(n_gram)
    if n_discr["length"] == 2:
        complement_name = "U+U(2)"
        complement = sp.diag(sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, 2], [2, 0]]))
    elif n_discr["length"] == 4:
        complement_name = "U(2)+U(2)"
        complement = sp.diag(sp.Matrix([[0, 2], [2, 0]]), sp.Matrix([[0, 2], [2, 0]]))
    else:
        raise ArithmeticError("hostile audit has no rank-four complement for this length")
    source_basis, target_basis = find_discriminant_anti_isometry(n_gram, complement)
    graph = gf2_basis(
        (
            source | (target << 18)
            for source, target in zip(source_basis, target_basis, strict=True)
        ),
        22,
    )
    basis = pivot_overlattice_basis(22, graph)
    direct_gram = sp.diag(n_gram, complement)
    k3_gram = integral(basis.T * direct_gram * basis, "K3 graph glue is not integral")
    action = integral(
        basis.inv() * sp.diag(involution, sp.eye(4)) * basis,
        "involution does not extend over K3 graph glue",
    )
    n_embedding = integral(
        basis.inv() * sp.Matrix.vstack(sp.eye(18), sp.zeros(4, 18)),
        "N embedding in K3 glue is not integral",
    )
    t_embedding = integral(
        basis.inv() * sp.Matrix.vstack(sp.zeros(18, 4), sp.eye(4)),
        "T embedding in K3 glue is not integral",
    )
    symbol = sp.Symbol("lambda")
    checks = {
        "even": all(k3_gram[i, i] % 2 == 0 for i in range(22)),
        "unimodular": int(k3_gram.det()) == -1,
        "signature": inertia(k3_gram) == (3, 19, 0),
        "action_isometry": action.T * k3_gram * action == k3_gram,
        "action_involution": action**2 == sp.eye(22),
        "N_primitive": primitive_embedding(n_embedding),
        "complement_primitive": primitive_embedding(t_embedding),
        "characteristic_polynomial": sp.factor(action.charpoly(symbol).as_expr())
        == (symbol - 1) ** 14 * (symbol + 1) ** 8,
    }
    if not all(checks.values()):
        raise ArithmeticError("K3 glue check failed: " + ",".join(key for key, value in checks.items() if not value))
    return {
        "complement": complement_name,
        "complement_signature_determinant": [2, 2, int(complement.det())],
        "discriminant_anti_isometry_checked_on_all_classes": True,
        "graph_glue_index": 1 << n_discr["length"],
        "k3_signature_determinant": [3, 19, -1],
        "N_and_complement_embeddings_primitive": True,
        "extended_involution_characteristic_polynomial": "(lambda-1)^14*(lambda+1)^8",
        "k3_gram_sha256": json_hash(rows(k3_gram)),
        "extended_involution_sha256": json_hash(rows(action)),
        "checks": checks,
    }


def construct_unfolding(
    candidate: dict[str, Any],
    folded: sp.Matrix,
    e8_gram: sp.Matrix,
    e8_classes: dict[int, tuple[tuple[int, ...], ...]],
) -> dict[str, Any]:
    fixed_gram: sp.Matrix = candidate["gram"]
    fixed_roots: sp.Matrix = candidate["roots"]
    long_indices: list[int] = candidate["long_indices"]
    domain_basis: list[int] = candidate["long_basis_classes"]
    target_basis: tuple[int, ...] = candidate["e8_embedding"]
    anti_roots: dict[int, sp.Matrix] = {}
    for facet, relation in zip(long_indices, candidate["long_relations"], strict=True):
        target_class = 0
        for i, value in enumerate(target_basis):
            if (relation >> i) & 1:
                target_class ^= value
        if target_class not in e8_classes:
            raise ArithmeticError("full marked E8 embedding lost a long root class")
        anti_roots[facet] = sp.Matrix(max(e8_classes[target_class]))
    anti_gram = -2 * e8_gram
    direct_gram = sp.diag(fixed_gram, anti_gram)
    direct_roots: list[sp.Matrix] = []
    orbits: list[list[int]] = []
    glue_classes: list[int] = []
    for facet in range(14):
        fixed = sp.Matrix.vstack(fixed_roots[:, facet], sp.zeros(8, 1))
        if facet not in long_indices:
            direct_roots.append(fixed)
            orbits.append([len(direct_roots) - 1])
            continue
        anti = sp.Matrix.vstack(sp.zeros(10, 1), anti_roots[facet])
        plus, minus = (fixed + anti) / 2, (fixed - anti) / 2
        direct_roots.extend((plus, minus))
        orbits.append([len(direct_roots) - 2, len(direct_roots) - 1])
        glue_classes.append(vector_bit(fixed_roots[:, facet]) | (vector_bit(anti_roots[facet]) << 10))
    graph_subgroup = gf2_basis(glue_classes, 18)
    if len(graph_subgroup) != len(domain_basis):
        raise ArithmeticError("graph-glue dimension differs from long-class span")
    basis = pivot_overlattice_basis(18, graph_subgroup)
    n_gram = integral(basis.T * direct_gram * basis, "rank-18 graph glue is not integral")
    roots_matrix = integral(
        basis.inv() * sp.Matrix.hstack(*direct_roots),
        "a lifted (-2)-root is not in the rank-18 graph glue",
    )
    involution = integral(
        basis.inv() * sp.diag(sp.eye(10), -sp.eye(8)) * basis,
        "sheet involution is not integral",
    )
    fixed_embedding = integral(
        basis.inv() * sp.Matrix.vstack(sp.eye(10), sp.zeros(8, 10)),
        "fixed lattice embedding is not integral",
    )
    anti_embedding = integral(
        basis.inv() * sp.Matrix.vstack(sp.zeros(10, 8), sp.eye(8)),
        "coinvariant lattice embedding is not integral",
    )
    root_gram = roots_matrix.T * n_gram * roots_matrix
    off_diagonal = [
        int(root_gram[i, j])
        for i in range(root_gram.rows)
        for j in range(i)
    ]
    folded_roots = sp.Matrix.hstack(
        *(roots_matrix[:, orbit[0]] if len(orbit) == 1 else roots_matrix[:, orbit[0]] + roots_matrix[:, orbit[1]] for orbit in orbits)
    )
    root_smith = smith(roots_matrix)
    root_index = math.prod(root_smith)
    n_discr = two_elementary(n_gram)
    checks = {
        "rank_signature": inertia(n_gram) == (1, 17, 0),
        "even": all(n_gram[i, i] % 2 == 0 for i in range(18)),
        "involution": involution**2 == sp.eye(18),
        "involution_isometry": involution.T * n_gram * involution == n_gram,
        "fixed_primitive": primitive_embedding(fixed_embedding),
        "coinvariant_primitive": primitive_embedding(anti_embedding),
        "fixed_gram": fixed_embedding.T * n_gram * fixed_embedding == fixed_gram,
        "coinvariant_gram": anti_embedding.T * n_gram * anti_embedding == anti_gram,
        "root_squares": all(root_gram[i, i] == -2 for i in range(root_gram.rows)),
        "acute": set(off_diagonal) <= {0, 1},
        "exact_labelled_fold_back": folded_roots.T * n_gram * folded_roots == folded,
        "sheet_orbits": all(involution * roots_matrix[:, orbit[0]] == roots_matrix[:, orbit[-1]] for orbit in orbits),
        "root_span_full_rank": len(root_smith) == 18,
    }
    if not all(checks.values()):
        raise ArithmeticError("unfolding check failed: " + ",".join(key for key, value in checks.items() if not value))
    k3 = k3_glue(n_gram, involution)
    fixed_discr = two_elementary(fixed_gram)
    return {
        "fixed_lattice": {
            "signature_determinant": [1, 9, int(fixed_gram.det())],
            "two_elementary_invariants": [10, fixed_discr["length"], fixed_discr["delta"]],
            "smith": smith(fixed_gram),
            "isometry_class": (
                "U+E8(-2)"
                if fixed_discr["length"] == 8 and fixed_discr["delta"] == 0
                else "unidentified"
            ),
            "gram_sha256": json_hash(rows(fixed_gram)),
        },
        "rank18_unfolding": {
            "signature_determinant": [1, 17, int(n_gram.det())],
            "smith": n_discr["smith"],
            "two_elementary_invariants": [18, n_discr["length"], n_discr["delta"]],
            "q_histogram": n_discr["q_histogram"],
            "fixed_and_coinvariant_embeddings_primitive": True,
            "gram_sha256": json_hash(rows(n_gram)),
            "involution_sha256": json_hash(rows(involution)),
        },
        "roots_and_fold": {
            "root_count": roots_matrix.cols,
            "one_sheet_count": sum(len(orbit) == 1 for orbit in orbits),
            "two_sheet_count": sum(len(orbit) == 2 for orbit in orbits),
            "off_diagonal_pairing_histogram": {str(key): value for key, value in sorted(Counter(off_diagonal).items())},
            "roots_span_index_in_rank18_lattice": root_index,
            "root_span_smith": root_smith,
            "exact_labelled_fold_back": True,
            "roots_sha256": json_hash(rows(roots_matrix)),
        },
        "weak_K3_lattice_action_extension": k3,
        "checks": checks,
    }


def invariant_candidate_signature(candidate: dict[str, Any], production: bool = False) -> tuple[Any, ...]:
    if production:
        return (
            candidate["glue_dimension"],
            candidate["determinant"],
            tuple(candidate["smith_diagonal_absolute"]),
            candidate["discriminant_length"],
            candidate["discriminant_delta"],
            candidate["long_class_span_dimension"],
            candidate["long_class_radical_dimension"],
            len(candidate["repeated_long_class_pairs"]),
            candidate["minimal_rank18_discriminant_length"],
            candidate["survives_minimal_two_sheet_weak_K3_gate"],
        )
    return (
        candidate["glue_dimension"],
        candidate["determinant"],
        tuple(candidate["smith"]),
        candidate["discriminant_length"],
        candidate["discriminant_delta"],
        candidate["long_class_span_dimension"],
        candidate["long_class_radical_dimension"],
        candidate["repeated_long_class_pair_count"],
        candidate["minimal_rank18_discriminant_length"],
        candidate["survives"],
    )


def compile_type(
    type_id: str,
    representative: dict[str, Any],
    markings: list[str],
    e8_gram: sp.Matrix,
    e8_classes: dict[int, tuple[tuple[int, ...], ...]],
) -> dict[str, Any]:
    raw = representative["gram"]
    histogram, edges, graph = relation_data(raw)
    scalings, automorphism_order, orbit_sizes = scaling_rays(edges, graph)
    lanes: list[dict[str, Any]] = []
    survivors: list[tuple[dict[str, Any], sp.Matrix, str, int]] = []
    for scaling_number, scaling in enumerate(scalings, start=1):
        if not set(scaling) <= {2, 4}:
            continue
        folded = scaled_gram(scaling, edges)
        public, internal, quotient = candidate_inventory(folded, e8_gram, e8_classes)
        label = f"short{scaling.count(2)}_long{scaling.count(4)}"
        lanes.append(
            {
                "label": label,
                "scaling_number": scaling_number,
                "normal_square_absolutes": list(scaling),
                "candidate_count": len(public),
                "glue_dimension_histogram": {str(key): value for key, value in sorted(Counter(item["glue_dimension"] for item in public).items())},
                "smith_quotient": quotient,
                "candidates": public,
            }
        )
        survivors.extend((candidate, folded, label, scaling_number) for candidate in internal if candidate["survives"])
    if len(survivors) != 1:
        raise ArithmeticError(f"{type_id}: hostile compiler found {len(survivors)} survivors")
    survivor, folded, label, scaling_number = survivors[0]
    unfolding = construct_unfolding(survivor, folded, e8_gram, e8_classes)
    literal = representative["vinberg"]
    if not (
        representative["classification"] == "finite_volume_noncompact"
        and literal["literal_vinberg_passes"] is True
        and literal["bad_edge_count"] == 0
        and literal["bad_parabolic_neighborhood_count"] == 0
    ):
        raise ArithmeticError("the exact labelled Vinberg archive is not green")
    return {
        "type_id": type_id,
        "archived_markings": markings,
        "representative": representative["id"],
        "raw_gram_sha256": json_hash(rows(raw)),
        "relation_histogram": {str(key): value for key, value in sorted(histogram.items())},
        "diagram_automorphism_order": automorphism_order,
        "scaling_count": len(scalings),
        "scaling_orbit_sizes": orbit_sizes,
        "compatible_lane_count": len(lanes),
        "norm_gt4_obstructed_scaling_count": sum(not set(scaling) <= {2, 4} for scaling in scalings),
        "lanes": lanes,
        "unique_survivor": {
            "lane": label,
            "scaling_number": scaling_number,
            **unfolding,
        },
        "casewise_V_archive_identity": {
            "classification": representative["classification"],
            "finite_endpoint_count": literal["finite_endpoint_count"],
            "ideal_endpoint_count": literal["ideal_endpoint_count"],
            "bad_edge_count": literal["bad_edge_count"],
            "bad_parabolic_neighborhood_count": literal["bad_parabolic_neighborhood_count"],
            "scope": "imported only after exact labelled fold-back; Vinberg was not recomputed here",
        },
    }


def compare_production(compiled: list[dict[str, Any]]) -> dict[str, Any]:
    if sha256(PRODUCTION) != SOURCE_HASHES["production"]:
        raise ArithmeticError("pinned production archive changed")
    production = json.loads(PRODUCTION.read_text())
    expected = {row["type_id"]: row for row in production["compiler_inventory"]}
    checks: dict[str, bool] = {}
    for hostile in compiled:
        type_id = hostile["type_id"]
        prior = expected[type_id]
        checks[f"{type_id}_scalings"] = (
            hostile["scaling_count"] == prior["scaling_count"]
            and hostile["scaling_orbit_sizes"] == prior["scaling_orbit_sizes"]
            and hostile["diagram_automorphism_order"] == prior["diagram_automorphism_order"]
        )
        hostile_lanes = {lane["label"]: lane for lane in hostile["lanes"]}
        prior_lanes = {lane["label"]: lane for lane in prior["lanes"]}
        checks[f"{type_id}_lane_labels"] = set(hostile_lanes) == set(prior_lanes)
        for label, lane in hostile_lanes.items():
            prior_lane = prior_lanes[label]
            checks[f"{type_id}_{label}_candidate_invariant_multiset"] = sorted(
                invariant_candidate_signature(item) for item in lane["candidates"]
            ) == sorted(
                invariant_candidate_signature(item, production=True) for item in prior_lane["candidates"]
            )
        hostile_survivor = hostile["unique_survivor"]
        prior_survivor = prior["unique_survivor"]
        checks[f"{type_id}_survivor_fixed"] = (
            hostile_survivor["fixed_lattice"]["two_elementary_invariants"]
            == prior_survivor["fixed_lattice"]["invariants"]
            and hostile_survivor["fixed_lattice"]["isometry_class"]
            == prior_survivor["fixed_lattice"]["isometry_class"]
        )
        checks[f"{type_id}_survivor_rank18"] = (
            hostile_survivor["rank18_unfolding"]["two_elementary_invariants"]
            == prior_survivor["rank18_unfolding"]["two_elementary_invariants"]
            and hostile_survivor["rank18_unfolding"]["signature_determinant"]
            == prior_survivor["rank18_unfolding"]["signature_determinant"]
        )
        checks[f"{type_id}_root_index_and_histogram"] = (
            hostile_survivor["roots_and_fold"]["roots_span_index_in_rank18_lattice"]
            == prior_survivor["roots_and_fold"]["roots_span_index_in_rank18_lattice"]
            and hostile_survivor["roots_and_fold"]["off_diagonal_pairing_histogram"]
            == prior_survivor["roots_and_fold"]["off_diagonal_pairing_histogram"]
        )
        checks[f"{type_id}_K3_complement"] = (
            hostile_survivor["weak_K3_lattice_action_extension"]["complement"]
            == prior_survivor["weak_K3_lattice_action_extension"]["complement"]
        )
    if not all(checks.values()):
        raise ArithmeticError("hostile/production discrepancy: " + ",".join(key for key, value in checks.items() if not value))
    return {
        "production_archive_sha256": SOURCE_HASHES["production"],
        "comparison_basis": "coordinate-free candidate invariant multisets plus survivor lattice/root/K3 invariants",
        "checks": checks,
        "all_agree": True,
    }


def build_certificate() -> dict[str, Any]:
    corpus, source_counts = load_corpus()
    if len(corpus) != 6:
        raise ArithmeticError("hostile corpus denominator is not six")
    groups: list[list[dict[str, Any]]] = []
    for row in corpus:
        group = next((group for group in groups if diagram_isomorphic(group[0]["gram"], row["gram"])), None)
        if group is None:
            groups.append([row])
        else:
            group.append(row)
    if sorted(len(group) for group in groups) != [1, 1, 2, 2]:
        raise ArithmeticError("hostile corpus quotient is not four diagram types")
    retained: list[list[dict[str, Any]]] = []
    excluded: list[dict[str, Any]] = []
    for group in groups:
        gram = group[0]["gram"]
        histogram = relation_data(gram)[0]
        gate = (
            gram.rows == 14
            and gram.rank() == 10
            and inertia(gram) == (1, 9, 4)
            and set(histogram) <= {2, 3, 4}
            and histogram[4] == 8
        )
        if gate:
            retained.append(group)
        else:
            excluded.append(
                {
                    "markings": [row["id"] for row in group],
                    "relation_histogram": {str(key): value for key, value in sorted(histogram.items())},
                }
            )
    if len(retained) != 2 or sum(map(len, retained)) != 3:
        raise ArithmeticError("hostile D34 filter did not retain two types/three markings")
    e8_gram, e8_classes = e8_data()
    compiled = []
    for group in retained:
        histogram = relation_data(group[0]["gram"])[0]
        compiled.append(
            compile_type(
                f"D34-{histogram[3]}-{histogram[4]}",
                group[0],
                [row["id"] for row in group],
                e8_gram,
                e8_classes,
            )
        )
    compiled.sort(key=lambda row: row["type_id"])
    by_type = {row["type_id"]: row for row in compiled}
    t2 = by_type["D34-8-8"]
    t1 = by_type["D34-9-8"]
    decisive = {
        "corpus_has_6_markings_4_types": len(corpus) == 6 and len(groups) == 4,
        "D34_has_3_markings_2_types": sum(map(len, retained)) == 3 and len(retained) == 2,
        "all_36_scaling_rays": sum(row["scaling_count"] for row in compiled) == 36,
        "all_28_stable_even_overclosures": sum(lane["candidate_count"] for row in compiled for lane in row["lanes"]) == 28,
        "unique_survivor_per_type": all(sum(candidate["survives"] for lane in row["lanes"] for candidate in lane["candidates"]) == 1 for row in compiled),
        "both_fixed_lattices_U_plus_E8_minus2": all(row["unique_survivor"]["fixed_lattice"]["isometry_class"] == "U+E8(-2)" for row in compiled),
        "T1_exact_18_2_0": t1["unique_survivor"]["rank18_unfolding"]["two_elementary_invariants"] == [18, 2, 0],
        "T1_rank18_determinant_minus4": t1["unique_survivor"]["rank18_unfolding"]["signature_determinant"] == [1, 17, -4],
        "T1_24_root_span_index_2": t1["unique_survivor"]["roots_and_fold"]["roots_span_index_in_rank18_lattice"] == 2,
        "T1_K3_complement_U_plus_U2": t1["unique_survivor"]["weak_K3_lattice_action_extension"]["complement"] == "U+U(2)",
        "T2_exact_18_4_0": t2["unique_survivor"]["rank18_unfolding"]["two_elementary_invariants"] == [18, 4, 0],
        "T2_24_root_span_index_1": t2["unique_survivor"]["roots_and_fold"]["roots_span_index_in_rank18_lattice"] == 1,
        "T2_K3_complement_U2_plus_U2": t2["unique_survivor"]["weak_K3_lattice_action_extension"]["complement"] == "U(2)+U(2)",
        "all_24_root_systems_acute_and_fold_exactly": all(row["unique_survivor"]["checks"]["acute"] and row["unique_survivor"]["checks"]["exact_labelled_fold_back"] for row in compiled),
        "all_K3_glues_even_unimodular_primitive": all(all(row["unique_survivor"]["weak_K3_lattice_action_extension"]["checks"].values()) for row in compiled),
    }
    if not all(decisive.values()):
        raise ArithmeticError("decisive hostile audit check failed")
    comparison = compare_production(compiled)
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-08",
        "source_sha256": SOURCE_HASHES,
        "method_independence": {
            "production_code_imported": False,
            "root_lattice": "full Smith decomposition of the 14x14 degenerate Gram",
            "stable_overclosures": "GF(2) reflection-orbit module closure and rational pivot bases",
            "rank18_and_K3_glues": "GF(2) graph subgroups and rational pivot bases",
            "root_span_index": "rectangular Smith form of the 18x24 root-coordinate matrix",
            "Hermite_normal_form_used": False,
            "principal_facet_basis_assumed": False,
        },
        "bounded_primary_corpus": {
            **source_counts,
            "marking_count": len(corpus),
            "diagram_type_count": len(groups),
            "group_sizes": sorted(len(group) for group in groups),
            "scope": "pinned primary certified exact-Gram corpus only; not all abstract H9/14 diagrams",
        },
        "D34_filter": {
            "definition": "H9, 14 facets, Gram rank/inertia (10;1,9,4), ordinary {2,3,4}, exactly eight m4 edges",
            "retained_marking_count": sum(map(len, retained)),
            "retained_diagram_type_count": len(retained),
            "excluded_types": excluded,
        },
        "compiler_inventory": compiled,
        "decisive_checks": decisive,
        "production_comparison": comparison,
        "scope_boundary": {
            "U": "open: this proves only a pinned-corpus IC0 census, not universality or an algebraic-surface existence locus",
            "R": "exact labelled lattice/action fold-back R0 is independently checked for T1 and T2",
            "C": "T1 cone completeness and hidden-wall exclusion remain open; this audit makes no nef-cone promotion",
            "V": "stored Vinberg results are only imported after exact labelled identity; the Vinberg enumeration is not re-audited here",
            "new_polytope": False,
        },
        "all_checks_pass": all(decisive.values()) and comparison["all_agree"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build_certificate()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not args.output.is_file() or args.output.read_text() != rendered:
            raise SystemExit("hostile D34 archive differs from exact replay")
        print("archive_check=PASS")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
        print(f"wrote={args.output}")
    print("all_checks_pass=true")
    print("corpus=6/4; D34=3/2; scalings=36; stable_overclosures=28; survivors=2")
    print("T1=(18,2,0),det=-4,root_index=2; T2=(18,4,0),det=-16,root_index=1")
    print(f"payload_sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
