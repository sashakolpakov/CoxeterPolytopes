#!/usr/bin/env python3
"""Independent exact audit of the proposed rho=16, H14 K3 wall fold.

The proposed wall permutation fixes 1,...,15 and swaps (16,18),(17,19).
This script deliberately imports neither the quotient-nef fold verifier nor
the Vinberg-volume implementation.  It rebuilds the published 19-wall Gram
matrix, the induced Neron--Severi isometry, the 17 orbit inequalities, all
extreme rays, and CoxIter Proposition 1.10 by exact rational arithmetic.

The metric fold is a finite-volume noncompact H14 Coxeter polytope.  The
same computation also detects the algebraic flaw: the proposed isometry acts
nontrivially on A_NS=(Z/2)^2, so it cannot be induced by an involution of a
K3 surface whose full Neron--Severi lattice is the advertised rank-16
lattice.  Accordingly this certificate makes no K3 quotient or novelty claim.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache, reduce
import hashlib
from itertools import combinations
import json
from math import gcd
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = ROOT / "data/known_examples/k3_high_rank_nef_chambers.json"
CERTIFICATE_PATH = ROOT / "data/candidates/k3_rho16_h14_fold_audit.json"

PROPOSED_CYCLES = tuple(
    [(index,) for index in range(1, 16)] + [(16, 18), (17, 19)]
)
PUBLISHED_NONTRIVIAL_DECK_CYCLES = (
    (1, 15),
    (2, 14),
    (3, 12),
    (4, 13),
    (5, 11),
    (6, 10),
    (7, 9),
    (8,),
    (16, 17),
    (18, 19),
)


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int | str]]:
    answer: list[list[int | str]] = []
    for row in range(matrix.rows):
        payload_row: list[int | str] = []
        for column in range(matrix.cols):
            value = sp.Rational(matrix[row, column])
            payload_row.append(
                int(value) if value.q == 1 else f"{int(value.p)}/{int(value.q)}"
            )
        answer.append(payload_row)
    return answer


def _matrix_sha256(matrix: sp.Matrix) -> str:
    return _canonical_sha256(_matrix_payload(matrix))


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Exact rational congruence elimination, independent of eigenvalues."""

    work = sp.Matrix(matrix)
    positive = 0
    negative = 0
    while work.rows:
        pivot = next((i for i in range(work.rows) if work[i, i] != 0), None)
        if pivot is None:
            pair = next(
                (
                    (i, j)
                    for i in range(work.rows)
                    for j in range(i + 1, work.rows)
                    if work[i, j] != 0
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            first, second = pair
            order = [first, second] + [
                i for i in range(work.rows) if i not in pair
            ]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0
        if pivot:
            order = [pivot] + [i for i in range(work.rows) if i != pivot]
            work = work.extract(order, order)
        diagonal = sp.Rational(work[0, 0])
        if diagonal > 0:
            positive += 1
        else:
            negative += 1
        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / diagonal
    return positive, negative, 0


def _build_ambient_gram(case: dict[str, Any]) -> sp.Matrix:
    size = int(case["wall_count"])
    gram = sp.diag(*([-2] * size))
    seen: set[tuple[int, int]] = set()
    labels: set[int] = set()
    for value_text, raw_pairs in case["intersections"].items():
        value = int(value_text)
        if value <= 0:
            raise AssertionError("curve intersections must be positive")
        for raw_left, raw_right in raw_pairs:
            left, right = sorted((int(raw_left), int(raw_right)))
            if not (1 <= left < right <= size) or (left, right) in seen:
                raise AssertionError("malformed or duplicate curve intersection")
            seen.add((left, right))
            labels.update((left, right))
            gram[left - 1, right - 1] = value
            gram[right - 1, left - 1] = value
    if labels != set(range(1, size + 1)):
        raise AssertionError("the source does not cover exactly walls 1,...,19")
    return gram


def _permutation_from_cycles(
    size: int, cycles: Sequence[Sequence[int]]
) -> tuple[int, ...]:
    permutation = list(range(size))
    seen: set[int] = set()
    for raw_cycle in cycles:
        cycle = tuple(int(label) - 1 for label in raw_cycle)
        if not cycle or seen.intersection(cycle):
            raise AssertionError("invalid or overlapping cycle list")
        if any(index < 0 or index >= size for index in cycle):
            raise AssertionError("cycle label outside the wall set")
        seen.update(cycle)
        for position, source in enumerate(cycle):
            permutation[source] = cycle[(position + 1) % len(cycle)]
    if seen != set(range(size)):
        raise AssertionError("cycles do not cover every wall")
    return tuple(permutation)


def _cycles(permutation: Sequence[int]) -> tuple[tuple[int, ...], ...]:
    seen: set[int] = set()
    answer: list[tuple[int, ...]] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        orbit: list[int] = []
        current = start
        while current not in seen:
            seen.add(current)
            orbit.append(current + 1)
            current = int(permutation[current])
        answer.append(tuple(orbit))
    return tuple(answer)


def _weighted_gram_automorphisms(gram: sp.Matrix) -> tuple[tuple[int, ...], ...]:
    """Exhaust the complete weighted graph by relation-preserving backtracking."""

    size = gram.rows
    relations = tuple(
        tuple(int(gram[left, right]) for right in range(size))
        for left in range(size)
    )
    signatures = tuple(tuple(sorted(row)) for row in relations)
    candidates = tuple(
        tuple(target for target in range(size) if signatures[target] == signatures[source])
        for source in range(size)
    )
    order = tuple(sorted(range(size), key=lambda source: (len(candidates[source]), source)))
    image = [-1] * size
    used: set[int] = set()
    answers: list[tuple[int, ...]] = []

    def visit(depth: int) -> None:
        if depth == size:
            permutation = tuple(image)
            if any(
                relations[left][right]
                != relations[permutation[left]][permutation[right]]
                for left in range(size)
                for right in range(size)
            ):
                raise AssertionError("false graph automorphism")
            answers.append(permutation)
            return
        source = order[depth]
        for target in candidates[source]:
            if target in used:
                continue
            if any(
                image[other] >= 0
                and relations[source][other]
                != relations[target][image[other]]
                for other in range(size)
            ):
                continue
            image[source] = target
            used.add(target)
            visit(depth + 1)
            used.remove(target)
            image[source] = -1

    visit(0)
    return tuple(sorted(answers))


def _orbit_incidence(size: int, cycles: Sequence[Sequence[int]]) -> sp.Matrix:
    incidence = sp.zeros(size, len(cycles))
    for column, orbit in enumerate(cycles):
        for label in orbit:
            incidence[int(label) - 1, column] = 1
    return incidence


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index] != 0
    ]


def _positive_primitive(vector: sp.Matrix) -> sp.Matrix:
    values = [sp.Rational(entry) for entry in vector]
    denominator = int(sp.ilcm(*[entry.q for entry in values]))
    integers = [int(entry * denominator) for entry in values]
    common = reduce(gcd, (abs(entry) for entry in integers), 0)
    if common == 0:
        raise AssertionError("zero projective vector")
    return sp.Matrix([entry // common for entry in integers])


def _enumerate_extreme_rays(
    fixed_gram: sp.Matrix, pairing_rows: sp.Matrix
) -> tuple[dict[str, Any], ...]:
    rank = fixed_gram.rows
    facets = pairing_rows.rows
    vertices: dict[tuple[int, ...], dict[str, Any]] = {}
    for active_candidate in combinations(range(facets), rank - 1):
        equations = pairing_rows.extract(active_candidate, range(rank))
        if equations.rank() != rank - 1:
            continue
        nullspace = equations.nullspace()
        if len(nullspace) != 1:
            raise AssertionError("codimension-one ray has wrong nullity")
        oriented: sp.Matrix | None = None
        for vector in (
            _positive_primitive(nullspace[0]),
            _positive_primitive(-nullspace[0]),
        ):
            pairings = pairing_rows * vector
            if all(bool(entry >= 0) for entry in pairings):
                oriented = vector
                break
        if oriented is None:
            continue
        pairings = pairing_rows * oriented
        active = tuple(
            index for index in range(facets) if pairings[index] == 0
        )
        if pairing_rows.extract(active, range(rank)).rank() != rank - 1:
            raise AssertionError("recorded active set does not span the ray equations")
        square = int((oriented.T * fixed_gram * oriented)[0])
        record = {
            "basis_coordinates": [int(entry) for entry in oriented],
            "active_facets": [index + 1 for index in active],
            "square": square,
            "type": (
                "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
            ),
        }
        key = tuple(int(entry) for entry in oriented)
        if key in vertices and vertices[key] != record:
            raise AssertionError("inconsistent duplicate extreme ray")
        vertices[key] = record
    return tuple(
        vertices[key]
        for key in sorted(vertices)
    )


def _facet_audit(
    vertices: Sequence[dict[str, Any]],
    fixed_gram: sp.Matrix,
    pairing_rows: sp.Matrix,
) -> tuple[dict[str, Any], ...]:
    answer: list[dict[str, Any]] = []
    for facet in range(pairing_rows.rows):
        equality_vectors = [
            sp.Matrix(vertex["basis_coordinates"])
            for vertex in vertices
            if facet + 1 in vertex["active_facets"]
        ]
        equality_rank = (
            sp.Matrix.hstack(*equality_vectors).rank() if equality_vectors else 0
        )
        witness = sum(equality_vectors, sp.zeros(fixed_gram.rows, 1))
        pairings = pairing_rows * witness
        other_pairings = [
            int(pairings[index])
            for index in range(pairing_rows.rows)
            if index != facet
        ]
        square = int((witness.T * fixed_gram * witness)[0])
        answer.append(
            {
                "facet": facet + 1,
                "equality_extreme_ray_count": len(equality_vectors),
                "equality_span_rank": equality_rank,
                "relative_interior_witness": [int(entry) for entry in witness],
                "relative_interior_witness_square": square,
                "target_pairing": int(pairings[facet]),
                "minimum_other_pairing": min(other_pairings),
                "irredundant_with_strict_relative_interior": (
                    equality_rank == fixed_gram.rows - 1
                    and pairings[facet] == 0
                    and square > 0
                    and min(other_pairings) > 0
                ),
            }
        )
    return tuple(answer)


def _coxeter_label(diagonal_left: int, diagonal_right: int, product: int) -> str:
    if diagonal_left <= 0 or diagonal_right <= 0 or product > 0:
        raise AssertionError("not an acute positive-normal Gram pair")
    ratio = Fraction(product * product, diagonal_left * diagonal_right)
    labels = {
        Fraction(0): "2",
        Fraction(1, 4): "3",
        Fraction(1, 2): "4",
        Fraction(3, 4): "6",
        Fraction(1): "infinity",
    }
    if ratio in labels:
        return labels[ratio]
    if ratio > 1:
        return "ultraparallel"
    raise AssertionError(f"non-Coxeter normalized product {ratio}")


class _LiteralDiagramAudit:
    """Fast exact rational principal-diagram classifier keyed by bit masks."""

    def __init__(self, wall_gram: sp.Matrix):
        self.gram = tuple(
            tuple(int(wall_gram[row, column]) for column in range(wall_gram.cols))
            for row in range(wall_gram.rows)
        )
        self.size = wall_gram.rows

    @lru_cache(maxsize=None)
    def pivots(self, mask: int) -> tuple[Fraction, ...]:
        indices = [index for index in range(self.size) if mask >> index & 1]
        matrix = [
            [Fraction(self.gram[left][right]) for right in indices]
            for left in indices
        ]
        pivots: list[Fraction] = []
        for pivot_index in range(len(matrix)):
            pivot = matrix[pivot_index][pivot_index]
            pivots.append(pivot)
            if pivot == 0:
                return tuple(pivots)
            for row in range(pivot_index + 1, len(matrix)):
                for column in range(row, len(matrix)):
                    value = (
                        matrix[row][column]
                        - matrix[row][pivot_index]
                        * matrix[pivot_index][column]
                        / pivot
                    )
                    matrix[row][column] = value
                    matrix[column][row] = value
        return tuple(pivots)

    @lru_cache(maxsize=None)
    def positive_definite(self, mask: int) -> bool:
        pivots = self.pivots(mask)
        return len(pivots) == mask.bit_count() and all(pivot > 0 for pivot in pivots)

    @lru_cache(maxsize=None)
    def components(self, mask: int) -> tuple[int, ...]:
        unseen = {index for index in range(self.size) if mask >> index & 1}
        components: list[int] = []
        while unseen:
            start = min(unseen)
            unseen.remove(start)
            component = 1 << start
            frontier = [start]
            while frontier:
                vertex = frontier.pop()
                neighbours = {
                    other
                    for other in unseen
                    if self.gram[vertex][other] != 0
                }
                unseen.difference_update(neighbours)
                for other in neighbours:
                    component |= 1 << other
                frontier.extend(neighbours)
            components.append(component)
        return tuple(components)

    @lru_cache(maxsize=None)
    def connected_affine(self, mask: int) -> bool:
        if mask.bit_count() < 2 or len(self.components(mask)) != 1:
            return False
        pivots = self.pivots(mask)
        if (
            len(pivots) != mask.bit_count()
            or pivots[-1] != 0
            or any(pivot <= 0 for pivot in pivots[:-1])
        ):
            return False
        return all(
            self.positive_definite(mask ^ (1 << vertex))
            for vertex in range(self.size)
            if mask >> vertex & 1
        )

    @lru_cache(maxsize=None)
    def parabolic(self, mask: int) -> bool:
        components = self.components(mask)
        return bool(components) and all(
            self.connected_affine(component) for component in components
        )

    def parabolic_rank(self, mask: int) -> int:
        if not self.parabolic(mask):
            raise AssertionError("rank requested for a non-parabolic diagram")
        return sum(component.bit_count() - 1 for component in self.components(mask))

    def mask(self, vertices: Iterable[int]) -> int:
        return sum(1 << int(vertex) for vertex in vertices)

    def labels(self, mask: int) -> list[int]:
        return [index + 1 for index in range(self.size) if mask >> index & 1]

    def literal_endpoint_audit(self, dimension: int) -> dict[str, Any]:
        edge_masks = tuple(
            mask
            for vertices in combinations(range(self.size), dimension - 1)
            if self.positive_definite(mask := self.mask(vertices))
        )
        finite_masks = tuple(
            mask
            for vertices in combinations(range(self.size), dimension)
            if self.positive_definite(mask := self.mask(vertices))
        )
        ideal_masks: list[int] = []
        for order in range(dimension, self.size + 1):
            for vertices in combinations(range(self.size), order):
                mask = self.mask(vertices)
                if self.parabolic(mask) and self.parabolic_rank(mask) == dimension - 1:
                    ideal_masks.append(mask)

        edge_extensions: list[dict[str, Any]] = []
        endpoint_type_histogram: Counter[tuple[int, int]] = Counter()
        extension_histogram: Counter[int] = Counter()
        for edge in edge_masks:
            finite = tuple(vertex for vertex in finite_masks if edge & vertex == edge)
            ideal = tuple(
                vertex
                for vertex in ideal_masks
                if edge & vertex == edge
                and all(
                    (component & ~edge).bit_count() == 1
                    for component in self.components(vertex)
                )
            )
            endpoint_type_histogram[(len(finite), len(ideal))] += 1
            extension_histogram[len(finite) + len(ideal)] += 1
            edge_extensions.append(
                {
                    "edge": self.labels(edge),
                    "finite_endpoints": [self.labels(vertex) for vertex in finite],
                    "ideal_endpoints": [self.labels(vertex) for vertex in ideal],
                    "endpoint_count": len(finite) + len(ideal),
                }
            )

        connected_affine_masks = tuple(
            mask
            for mask in range(1, 1 << self.size)
            if self.connected_affine(mask)
        )
        neighbourhood_audits: list[dict[str, Any]] = []
        full_neighbourhoods: set[int] = set()
        for affine in connected_affine_masks:
            zero_mask = self.mask(
                vertex
                for vertex in range(self.size)
                if not (affine >> vertex & 1)
                and all(
                    self.gram[vertex][node] == 0
                    for node in range(self.size)
                    if affine >> node & 1
                )
            )
            neighbourhood = affine | zero_mask
            full = (
                self.parabolic(neighbourhood)
                and self.parabolic_rank(neighbourhood) == dimension - 1
            )
            if full:
                full_neighbourhoods.add(neighbourhood)
            neighbourhood_audits.append(
                {
                    "affine_component": self.labels(affine),
                    "full_zero_neighbourhood": self.labels(neighbourhood),
                    "neighbourhood_components": [
                        self.labels(component)
                        for component in self.components(neighbourhood)
                    ],
                    "is_full_rank_dimension_minus_one_parabolic": full,
                }
            )

        return {
            "dimension": dimension,
            "elliptic_rank_n_minus_one_edges": [
                self.labels(mask) for mask in edge_masks
            ],
            "elliptic_rank_n_finite_vertices": [
                self.labels(mask) for mask in finite_masks
            ],
            "parabolic_rank_n_minus_one_ideal_vertices": [
                {
                    "facets": self.labels(mask),
                    "affine_components": [
                        self.labels(component) for component in self.components(mask)
                    ],
                }
                for mask in ideal_masks
            ],
            "edge_extensions": edge_extensions,
            "elliptic_edge_count": len(edge_masks),
            "finite_vertex_count": len(finite_masks),
            "ideal_vertex_count": len(ideal_masks),
            "extension_histogram": {
                str(key): value for key, value in sorted(extension_histogram.items())
            },
            "endpoint_type_histogram": {
                f"{finite},{ideal}": value
                for (finite, ideal), value in sorted(endpoint_type_histogram.items())
            },
            "bad_edges": [
                audit["edge"]
                for audit in edge_extensions
                if audit["endpoint_count"] != 2
            ],
            "connected_affine_zero_neighbourhood_audits": neighbourhood_audits,
            "connected_affine_component_count": len(connected_affine_masks),
            "distinct_full_zero_neighbourhood_count": len(full_neighbourhoods),
            "full_zero_neighbourhoods_match_ideal_vertices": (
                full_neighbourhoods == set(ideal_masks)
            ),
            "finite_volume_noncompact": (
                bool(finite_masks or ideal_masks)
                and bool(ideal_masks)
                and all(audit["endpoint_count"] == 2 for audit in edge_extensions)
                and all(
                    audit["is_full_rank_dimension_minus_one_parabolic"]
                    for audit in neighbourhood_audits
                )
                and full_neighbourhoods == set(ideal_masks)
            ),
        }


def _reflection_matrix(normal: sp.Matrix, gram: sp.Matrix) -> sp.Matrix:
    square = sp.Rational((normal.T * gram * normal)[0])
    return sp.eye(gram.rows) - sp.Rational(2, 1) / square * normal * normal.T * gram


def _is_integral_matrix(matrix: sp.Matrix) -> bool:
    return all(sp.Rational(entry).q == 1 for entry in matrix)


def _discriminant_action_trivial(action: sp.Matrix, gram: sp.Matrix) -> bool:
    return _is_integral_matrix((action - sp.eye(action.rows)) * gram.inv())


def build_certificate() -> dict[str, Any]:
    source = json.loads(SOURCE_PATH.read_text())
    case = next(case for case in source["cases"] if int(case["rho"]) == 16)
    ambient = _build_ambient_gram(case)
    basis = ambient[:16, :16]
    wall_coordinates = basis.inv() * ambient[:16, :]
    if not _is_integral_matrix(wall_coordinates):
        raise AssertionError("the 19 walls are not integral in the published basis")

    proposed_permutation = _permutation_from_cycles(19, PROPOSED_CYCLES)
    permutation_matrix = sp.zeros(19)
    for source_index, target_index in enumerate(proposed_permutation):
        permutation_matrix[target_index, source_index] = 1
    proposed_action = wall_coordinates[:, list(proposed_permutation[:16])]
    incidence = _orbit_incidence(19, PROPOSED_CYCLES)
    orbit_coordinates_16 = wall_coordinates * incidence
    orbit_gram = incidence.T * ambient * incidence

    # The proposed action fixes the first fifteen published basis vectors.
    # They are therefore a saturated basis of the full fixed lattice.
    fixed_gram = basis[:15, :15]
    fixed_normals = orbit_coordinates_16[:15, :]
    pairing_rows = fixed_normals.T * fixed_gram
    if pairing_rows != orbit_gram[:, :15]:
        raise AssertionError("fixed wall pairing rows disagree with the orbit Gram")

    reflection_records: list[dict[str, Any]] = []
    for index in range(fixed_normals.cols):
        normal = fixed_normals[:, index]
        reflection = _reflection_matrix(normal, fixed_gram)
        reflection_records.append(
            {
                "facet": index + 1,
                "ambient_wall_orbit": list(PROPOSED_CYCLES[index]),
                "normal_coordinates": [int(entry) for entry in normal],
                "normal_square": int((normal.T * fixed_gram * normal)[0]),
                "matrix": _matrix_payload(reflection),
                "matrix_sha256": _matrix_sha256(reflection),
                "integral": _is_integral_matrix(reflection),
                "involution": reflection * reflection == sp.eye(15),
                "isometry": reflection.T * fixed_gram * reflection == fixed_gram,
            }
        )

    interior = fixed_gram.inv() * sp.Matrix([2] * 15)
    interior_pairings = pairing_rows * interior
    interior_square = int((interior.T * fixed_gram * interior)[0])
    vertices = _enumerate_extreme_rays(fixed_gram, pairing_rows)
    facets = _facet_audit(vertices, fixed_gram, pairing_rows)
    finite_vertices = tuple(vertex for vertex in vertices if vertex["type"] == "finite")
    ideal_vertices = tuple(vertex for vertex in vertices if vertex["type"] == "ideal")
    hyperideal_vertices = tuple(
        vertex for vertex in vertices if vertex["type"] == "hyperideal"
    )

    positive_wall_gram = -orbit_gram
    relation_counts: Counter[str] = Counter()
    relation_edges: list[dict[str, Any]] = []
    for left in range(17):
        for right in range(left + 1, 17):
            label = _coxeter_label(
                int(positive_wall_gram[left, left]),
                int(positive_wall_gram[right, right]),
                int(positive_wall_gram[left, right]),
            )
            relation_counts[label] += 1
            if label != "2":
                relation_edges.append(
                    {"facets": [left + 1, right + 1], "label": label}
                )

    diagram = _LiteralDiagramAudit(positive_wall_gram)
    endpoint_audit = diagram.literal_endpoint_audit(14)
    finite_endpoint_sets = {
        tuple(facets)
        for facets in endpoint_audit["elliptic_rank_n_finite_vertices"]
    }
    ideal_endpoint_sets = {
        tuple(endpoint["facets"])
        for endpoint in endpoint_audit["parabolic_rank_n_minus_one_ideal_vertices"]
    }
    extreme_finite_sets = {
        tuple(vertex["active_facets"]) for vertex in finite_vertices
    }
    extreme_ideal_sets = {
        tuple(vertex["active_facets"]) for vertex in ideal_vertices
    }

    graph_automorphisms = _weighted_gram_automorphisms(ambient)
    deck_permutation = _permutation_from_cycles(
        19, PUBLISHED_NONTRIVIAL_DECK_CYCLES
    )
    deck_action = wall_coordinates[:, list(deck_permutation[:16])]
    dual_witness = basis.inv()[:, 3]
    dual_witness_difference = (
        proposed_action - sp.eye(16)
    ) * dual_witness
    discriminant_smith = _smith_diagonal(basis)
    proposed_discriminant_trivial = _discriminant_action_trivial(
        proposed_action, basis
    )
    deck_discriminant_trivial = _discriminant_action_trivial(deck_action, basis)

    checks = {
        "source_is_full_rho16_nef_configuration": (
            source["scope"].startswith("Full nef chambers")
            and int(case["rho"]) == 16
            and int(case["wall_count"]) == 19
        ),
        "ambient_gram_rank_signature_and_basis": (
            ambient.rank() == 16
            and _exact_inertia(ambient) == (1, 15, 3)
            and int(basis.det()) == -4
            and _smith_diagonal(basis) == [1] * 14 + [2, 2]
        ),
        "proposed_cycles_are_an_involutive_gram_symmetry": (
            permutation_matrix.T * ambient * permutation_matrix == ambient
            and all(
                proposed_permutation[proposed_permutation[index]] == index
                for index in range(19)
            )
        ),
        "induced_NS_action_is_integral_involutive_isometry": (
            _is_integral_matrix(proposed_action)
            and proposed_action * proposed_action == sp.eye(16)
            and proposed_action.T * basis * proposed_action == basis
        ),
        "fixed_lattice_is_saturated_rank15_lorentzian": (
            proposed_action[:15, :15] == sp.eye(15)
            and (proposed_action - sp.eye(16)).rank() == 1
            and orbit_coordinates_16.rank() == 15
            and proposed_action * orbit_coordinates_16 == orbit_coordinates_16
            and _exact_inertia(fixed_gram) == (1, 14, 0)
            and int(fixed_gram.det()) == 4
        ),
        "orbit_gram_rank_signature_and_norms": (
            orbit_gram.rank() == 15
            and _exact_inertia(orbit_gram) == (1, 14, 2)
            and list(orbit_gram.diagonal()) == [-2] * 15 + [-4, -4]
        ),
        "complete_fixed_nef_inequalities_are_all_17_orbits": (
            len(PROPOSED_CYCLES) == 17
            and set(label for orbit in PROPOSED_CYCLES for label in orbit)
            == set(range(1, 20))
            and pairing_rows.rows == 17
            and pairing_rows.rank() == 15
        ),
        "strict_chamber_interior": (
            _is_integral_matrix(interior)
            and interior_square == 1540
            and [int(entry) for entry in interior_pairings]
            == [2] * 15 + [44, 44]
        ),
        "all_17_walls_irredundant_with_strict_relative_interiors": all(
            record["irredundant_with_strict_relative_interior"]
            for record in facets
        ),
        "all_folded_reflections_integral_involutive_isometries": all(
            record["integral"] and record["involution"] and record["isometry"]
            for record in reflection_records
        ),
        "exact_acute_coxeter_relations": (
            relation_counts
            == Counter({"2": 119, "3": 14, "4": 2, "infinity": 1})
        ),
        "extreme_rays_are_89_finite_5_ideal_no_hyperideal": (
            len(finite_vertices) == 89
            and len(ideal_vertices) == 5
            and not hyperideal_vertices
        ),
        "literal_coxiter_704_edges_have_exactly_two_endpoints": (
            endpoint_audit["elliptic_edge_count"] == 704
            and endpoint_audit["extension_histogram"] == {"2": 704}
            and not endpoint_audit["bad_edges"]
        ),
        "literal_endpoint_sets_match_extreme_rays": (
            finite_endpoint_sets == extreme_finite_sets
            and ideal_endpoint_sets == extreme_ideal_sets
        ),
        "all_affine_zero_neighbourhoods_are_full_parabolic_cusps": (
            endpoint_audit["connected_affine_component_count"] == 8
            and endpoint_audit["distinct_full_zero_neighbourhood_count"] == 5
            and endpoint_audit["full_zero_neighbourhoods_match_ideal_vertices"]
            and all(
                audit["is_full_rank_dimension_minus_one_parabolic"]
                for audit in endpoint_audit[
                    "connected_affine_zero_neighbourhood_audits"
                ]
            )
        ),
        "abstract_fold_is_finite_volume_noncompact": (
            endpoint_audit["finite_volume_noncompact"] is True
        ),
        "weighted_graph_automorphism_census_is_V4": (
            len(graph_automorphisms) == 4
            and proposed_permutation in graph_automorphisms
            and deck_permutation in graph_automorphisms
        ),
        "proposed_action_is_nontrivial_on_discriminant": (
            discriminant_smith == [1] * 14 + [2, 2]
            and proposed_discriminant_trivial is False
            and not _is_integral_matrix(dual_witness_difference)
        ),
        "published_deck_action_is_trivial_on_discriminant": (
            deck_discriminant_trivial is True
        ),
        "no_K3_involution_with_full_NS_realizes_proposed_action": (
            proposed_discriminant_trivial is False
            and discriminant_smith == [1] * 14 + [2, 2]
        ),
    }
    if not all(checks.values()):
        raise AssertionError(f"rho16 H14 fold audit failed: {checks}")

    graph_cycle_census = sorted(
        (list(map(list, _cycles(permutation))) for permutation in graph_automorphisms),
        key=lambda cycles: (len(cycles), cycles),
    )
    payload: dict[str, Any] = {
        "schema_version": "1.0.0",
        "artifact_id": "independent_k3_rho16_h14_fold_audit",
        "certificate_kind": (
            "exact_abstract_fold_plus_algebraic_realization_obstruction"
        ),
        "source": {
            "path": str(SOURCE_PATH.relative_to(ROOT)),
            "sha256": _sha256_file(SOURCE_PATH),
            "artifact_id": source["artifact_id"],
            "artifact_status": source["artifact_status"],
            "scope": source["scope"],
            "rho16_locator": case["source_locator"],
            "published_model": case["published_model"],
            "published_general_automorphism_group": case[
                "general_automorphism_group"
            ],
        },
        "proposed_fold": {
            "cycles": [list(cycle) for cycle in PROPOSED_CYCLES],
            "ambient_picard_rank": 16,
            "fixed_rank": 15,
            "hyperbolic_dimension": 14,
            "ambient_wall_count": 19,
            "folded_facet_count": 17,
            "status": "abstract_lattice_and_complete_halfspace_fold_only",
        },
        "ambient_lattice": {
            "name": case["lattice"],
            "wall_gram_rank": ambient.rank(),
            "wall_gram_inertia": list(_exact_inertia(ambient)),
            "basis_determinant": int(basis.det()),
            "basis_smith_diagonal": discriminant_smith,
            "discriminant_group": "Z/2 + Z/2",
        },
        "fixed_nef_fold": {
            "orbit_gram_rank": orbit_gram.rank(),
            "orbit_gram_inertia": list(_exact_inertia(orbit_gram)),
            "orbit_normal_squares": [int(entry) for entry in orbit_gram.diagonal()],
            "completeness_argument": (
                "The source supplies all 19 nef inequalities. On a fixed class, "
                "the two inequalities in each 2-cycle have equal pairing and "
                "collapse to their orbit-sum inequality; singleton inequalities "
                "remain unchanged. The 17 cycles cover all 19 walls."
            ),
            "strict_interior": {
                "basis_coordinates": [int(entry) for entry in interior],
                "square": interior_square,
                "wall_pairings": [int(entry) for entry in interior_pairings],
            },
            "facet_audits": list(facets),
            "folded_reflections": reflection_records,
            "coxeter_relation_counts": dict(sorted(relation_counts.items())),
            "nonright_coxeter_edges": relation_edges,
            "vertices": list(vertices),
            "finite_vertex_count": len(finite_vertices),
            "ideal_vertex_count": len(ideal_vertices),
            "hyperideal_vertex_count": len(hyperideal_vertices),
        },
        "literal_vinberg_coxiter_prop_1_10": endpoint_audit,
        "algebraic_realization_audit": {
            "weighted_graph_automorphism_count": len(graph_automorphisms),
            "weighted_graph_group": "V4",
            "weighted_graph_cycle_census": graph_cycle_census,
            "published_nontrivial_D2_deck_cycles": [
                list(cycle) for cycle in PUBLISHED_NONTRIVIAL_DECK_CYCLES
            ],
            "proposed_action_trivial_on_A_NS": proposed_discriminant_trivial,
            "published_deck_action_trivial_on_A_NS": deck_discriminant_trivial,
            "dual_discriminant_witness": {
                "definition": "fourth dual basis vector B^{-1}e_4",
                "coordinates": _matrix_payload(dual_witness.T)[0],
                "J_witness_minus_witness": _matrix_payload(
                    dual_witness_difference.T
                )[0],
                "difference_is_integral": _is_integral_matrix(
                    dual_witness_difference
                ),
            },
            "involution_gluing_obstruction": (
                "If an involution f of a K3 surface with full NS=N induces J, "
                "then f acts on the period by epsilon=+/-1. The opposite "
                "eigenspace in T=N^perp is rational and orthogonal to the "
                "period, hence would contain an extra integral (1,1) class. "
                "Thus f|T=epsilon*I. Since A_T is 2-elementary, +/-I acts "
                "trivially on A_T; unimodular K3 gluing forces J to act "
                "trivially on A_N. The explicit dual witness disproves this."
            ),
            "conclusion": (
                "No involution of a K3 surface whose full Neron-Severi lattice "
                "is U+E8(-1)+D6(-1) can induce the proposed small wall action."
            ),
            "special_period_warning": (
                "A realization after a Picard-rank jump would no longer inherit "
                "the certified completeness of the 19-wall nef configuration."
            ),
        },
        "matrices": {
            "ambient_curve_gram_19": _matrix_payload(ambient),
            "published_basis_gram_16": _matrix_payload(basis),
            "wall_coordinates_in_published_basis_16x19": _matrix_payload(
                wall_coordinates
            ),
            "proposed_NS_action_16": _matrix_payload(proposed_action),
            "orbit_incidence_19x17": _matrix_payload(incidence),
            "orbit_gram_17": _matrix_payload(orbit_gram),
            "fixed_lattice_gram_15": _matrix_payload(fixed_gram),
            "folded_normal_coordinates_15x17": _matrix_payload(fixed_normals),
            "positive_normal_coxeter_gram_17": _matrix_payload(
                positive_wall_gram
            ),
        },
        "matrix_sha256": {
            "ambient_curve_gram_19": _matrix_sha256(ambient),
            "published_basis_gram_16": _matrix_sha256(basis),
            "proposed_NS_action_16": _matrix_sha256(proposed_action),
            "orbit_gram_17": _matrix_sha256(orbit_gram),
            "fixed_lattice_gram_15": _matrix_sha256(fixed_gram),
            "positive_normal_coxeter_gram_17": _matrix_sha256(
                positive_wall_gram
            ),
            "all_extreme_rays": _canonical_sha256(vertices),
            "all_edge_extensions": _canonical_sha256(
                endpoint_audit["edge_extensions"]
            ),
        },
        "checks": checks,
        "source_sha256": {
            "certifier": _sha256_file(Path(__file__).resolve()),
            "k3_high_rank_nef_chambers": _sha256_file(SOURCE_PATH),
        },
        "commands": [
            "PYTHONPATH=src:. python3 scripts/independent_certify_k3_rho16_h14_fold.py --write-certificate",
            "PYTHONPATH=src:. python3 scripts/independent_certify_k3_rho16_h14_fold.py --verify-certificate",
            "PYTHONPATH=src:. pytest -q tests/test_k3_rho16_h14_fold_audit.py",
        ],
        "conclusion": {
            "abstract_hyperbolic_coxeter_polytope": True,
            "finite_volume": True,
            "compact": False,
            "algebraic_K3_quotient_nef_realization": False,
            "new_hyperbolic_example_claim": False,
            "summary": (
                "The proposed diagram fold is an exact finite-volume "
                "noncompact H14 Coxeter polytope, but the specified small "
                "involution cannot be realized on a K3 with this full NS."
            ),
        },
    }
    payload["certificate_payload_sha256"] = _canonical_sha256(payload)
    return payload


def verify_certificate(payload: dict[str, Any]) -> dict[str, bool]:
    body = dict(payload)
    digest = body.pop("certificate_payload_sha256", None)
    fresh = build_certificate()
    checks = {
        "payload_digest_matches": digest == _canonical_sha256(body),
        "source_hashes_match": payload.get("source_sha256")
        == fresh.get("source_sha256"),
        "full_independent_recomputation_matches": payload == fresh,
        "all_exact_checks_are_true": all(payload["checks"].values()),
    }
    if not all(checks.values()):
        raise AssertionError(f"rho16 H14 certificate verification failed: {checks}")
    return checks


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-certificate", action="store_true")
    parser.add_argument("--verify-certificate", action="store_true")
    arguments = parser.parse_args()
    if arguments.write_certificate:
        payload = build_certificate()
        CERTIFICATE_PATH.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )
        print(f"wrote {CERTIFICATE_PATH}")
        print(payload["certificate_payload_sha256"])
        return
    payload = json.loads(CERTIFICATE_PATH.read_text())
    print(json.dumps(verify_certificate(payload), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
