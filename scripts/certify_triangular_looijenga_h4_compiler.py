#!/usr/bin/env python3
"""Surface-free compiler for triangular full-nef Looijenga carriers in H4.

The raw input is an integral Coxeter wall system in I_(1,4), not a surface.
Adjunction recovers the canonical class from the wall norms.  The compiler
then exhausts every bounded anticanonical triangle, plane/infinitely-near
blowup certificate, period-root branch, numerical exceptional class, and
hidden boundary/exceptional inequality.  The surviving complete nef cones
are classified by their weighted Gram matrices and checked by exact rays.

The finite arithmetic is proved and replayed here.  The passage from a
nonempty period/chamber stratum to an actual marked Looijenga pair uses the
period-surjectivity, root, and ample-cone theorems of Friedman and of
Gross--Hacking--Keel, which are listed as external theorem inputs in the
archive and companion report.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, permutations, product
import json
from math import gcd, isqrt
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "triangular_looijenga_h4_compiler.json"
)
SCHEMA = "coxeter.triangular_looijenga_h4_compiler.v1"

# Coordinates are (H,E1,E2,E3,E4), with H^2=1 and Ei^2=-1.
J = sp.diag(1, -1, -1, -1, -1)
K = (-3, 1, 1, 1, 1)
MINUS_K = tuple(-entry for entry in K)

BOUNDARY_SEQUENCE_ORDER = (
    (-2, -2, 3),
    (-2, -1, 2),
    (-2, 0, 1),
    (-1, -1, 1),
    (-1, 0, 0),
)

# In a blowdown basis (h;e1,...,e4), a row (a;m1,...,m4) denotes
# a*h-sum(m_i*e_i).  A negative m_i records an exceptional boundary
# component.  The last template has e1-e4 and hence one infinitely-near
# blowup.  Rows and exceptional labels are canonicalized before comparison.
PLANE_BLOWUP_TEMPLATES: dict[tuple[int, int, int], dict[str, Any]] = {
    (-1, 0, 0): {
        "rows": ((1, 0, 0, 0, 1), (1, 0, 0, 1, 0), (1, 1, 1, 0, 0)),
        "model": "P2 triangle; smooth-boundary blowups distributed (2,1,1)",
    },
    (-1, -1, 1): {
        "rows": ((1, 0, 0, 0, 0), (1, 0, 0, 1, 1), (1, 1, 1, 0, 0)),
        "model": "P2 triangle; smooth-boundary blowups distributed (2,2,0)",
    },
    (-2, 0, 1): {
        "rows": ((1, 0, 0, 0, 0), (1, 0, 0, 0, 1), (1, 1, 1, 1, 0)),
        "model": "P2 triangle; smooth-boundary blowups distributed (3,1,0)",
    },
    (-2, -1, 2): {
        "rows": ((0, -1, 0, 0, 0), (1, 1, 0, 1, 1), (2, 1, 1, 0, 0)),
        "model": (
            "line+conic cubic; blow up one line-conic node so its exceptional "
            "curve is the third boundary component, then make three smooth-boundary blowups"
        ),
    },
    (-2, -2, 3): {
        "rows": ((0, -1, 0, 0, 1), (1, 1, 1, 1, 0), (2, 1, 0, 0, 0)),
        "model": (
            "line+conic cubic; blow up one node, two further smooth points on the "
            "line, and one point of the exceptional boundary away from its nodes "
            "(the last point is infinitely near)"
        ),
    },
}


Vector = tuple[int, int, int, int, int]
Triangle = tuple[Vector, Vector, Vector]
WallSet = tuple[Vector, ...]


def _v(values: Iterable[int | sp.Rational]) -> sp.Matrix:
    return sp.Matrix(list(values))


def pair(left: Iterable[int], right: Iterable[int]) -> sp.Rational:
    left_values = tuple(left)
    right_values = tuple(right)
    return sp.Rational(
        left_values[0] * right_values[0]
        - sum(
            left_values[index] * right_values[index]
            for index in range(1, 5)
        )
    )


def square(vector: Iterable[int]) -> sp.Rational:
    return pair(vector, vector)


def _rank(vectors: Sequence[Vector]) -> int:
    if not vectors:
        return 0
    return sp.Matrix.hstack(*(_v(vector) for vector in vectors)).rank()


def _add(*vectors: Vector) -> Vector:
    return tuple(sum(entries) for entries in zip(*vectors, strict=True))  # type: ignore[return-value]


def _sub(left: Vector, right: Vector) -> Vector:
    return tuple(a - b for a, b in zip(left, right, strict=True))  # type: ignore[return-value]


def _canonical_sha256(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _rational_payload(value: sp.Rational) -> int | str:
    value = sp.Rational(value)
    if value.q == 1:
        return int(value)
    return f"{int(value.p)}/{int(value.q)}"


def _vector_payload(vector: Iterable[sp.Rational]) -> list[int | str]:
    return [_rational_payload(value) for value in vector]


def _primitive(vector: sp.Matrix) -> Vector:
    denominators = [int(sp.Rational(value).q) for value in vector]
    scale = int(sp.ilcm(*denominators))
    values = [int(sp.Rational(value) * scale) for value in vector]
    common = 0
    for value in values:
        common = gcd(common, abs(value))
    values = [value // common for value in values]
    return tuple(values)  # type: ignore[return-value]


def _is_characteristic(vector: Vector) -> bool:
    basis = [tuple(int(i == j) for i in range(5)) for j in range(5)]
    return all((pair(vector, basis[index]) - J[index, index]) % 2 == 0 for index in range(5))


def _lattice_index(vectors: Sequence[Vector]) -> int:
    """Index of the lattice generated by a full-rank vector family in Z^5."""

    if _rank(vectors) != 5:
        return 0
    columns = sp.Matrix.hstack(*(_v(vector) for vector in vectors))
    common = 0
    for indices in combinations(range(len(vectors)), 5):
        common = gcd(common, abs(int(columns[:, indices].det())))
    return common


def _gcd_maximal_minors(columns: sp.Matrix) -> int:
    rank = columns.cols
    common = 0
    for rows in combinations(range(columns.rows), rank):
        common = gcd(common, abs(int(columns.extract(rows, range(rank)).det())))
    return common


def _gcd_pair_minors(left: Vector, right: Vector) -> int:
    common = 0
    for i, j in combinations(range(5), 2):
        common = gcd(common, abs(left[i] * right[j] - left[j] * right[i]))
    return common


def _recover_canonical(walls: Sequence[Vector]) -> tuple[sp.Rational, ...] | None:
    if _rank(walls) != 5:
        return None
    equations = sp.Matrix.vstack(*((_v(wall).T * J) for wall in walls))
    targets = sp.Matrix([-2 - square(wall) for wall in walls])
    solution_set = sp.linsolve((equations, targets))
    if solution_set is sp.EmptySet or len(solution_set.args) != 1:
        return None
    solution = tuple(sp.Rational(value) for value in solution_set.args[0])
    if any(value.free_symbols for value in solution):
        return None
    return solution


def numerical_exceptional_classes() -> list[Vector]:
    """All e with e^2=K.e=-1; Cauchy forces the H coefficient to be 0 or 1."""

    result: list[Vector] = []
    for degree in (0, 1):
        for coefficients in product(range(-2, 3), repeat=4):
            vector: Vector = (degree, *coefficients)
            if square(vector) == -1 and pair(K, vector) == -1:
                result.append(vector)
    return sorted(result)


def canonical_roots() -> list[Vector]:
    """All roots of K^perp; its positive opposite has the A4 norm form."""

    result: list[Vector] = []
    # If u=(a,b,c,d,-3a-b-c-d) and -u^2=2, the inverse-form bound gives
    # |a|,|b|,|c|,|d|<2.  Thus this is an exact, not heuristic, box.
    for a, b, c, d in product(range(-1, 2), repeat=4):
        vector: Vector = (a, b, c, d, -3 * a - b - c - d)
        if square(vector) == -2:
            result.append(vector)
    return sorted(result)


def boundary_component_candidates() -> list[Vector]:
    """All bounded triangle components, using the exact K^perp norm bound."""

    result: set[Vector] = set()
    for component_square in range(-2, 4):
        canonical_pairing = -component_square - 2
        # u=5d-(K.d)K is in K^perp and
        # -u^2=5(s^2-s+4), at most 50 for -2<=s<=3.
        target_u_square = 5 * (-component_square * component_square + component_square - 4)
        for a in range(-6, 7):
            for b, c, d in product(range(-7, 8), repeat=3):
                e = -3 * a - b - c - d
                if abs(e) > 7:
                    continue
                u: Vector = (a, b, c, d, e)
                if square(u) != target_u_square:
                    continue
                numerators = [
                    u[index] + canonical_pairing * K[index]
                    for index in range(5)
                ]
                if any(value % 5 for value in numerators):
                    continue
                component: Vector = tuple(value // 5 for value in numerators)  # type: ignore[assignment]
                if (
                    square(component) == component_square
                    and pair(K, component) == canonical_pairing
                ):
                    result.add(component)
    return sorted(result)


def anticanonical_triangles(components: Sequence[Vector]) -> list[Triangle]:
    lookup = set(components)
    result: set[Triangle] = set()
    for left_index, left in enumerate(components):
        for middle in components[left_index:]:
            if pair(left, middle) != 1:
                continue
            right: Vector = tuple(
                MINUS_K[index] - left[index] - middle[index]
                for index in range(5)
            )  # type: ignore[assignment]
            if right not in lookup:
                continue
            if pair(left, right) != 1 or pair(middle, right) != 1:
                continue
            result.add(tuple(sorted((left, middle, right))))
    return sorted(result)


def exceptional_blowdown_bases(exceptionals: Sequence[Vector]) -> list[tuple[Vector, tuple[Vector, ...]]]:
    result: list[tuple[Vector, tuple[Vector, ...]]] = []
    for exceptional_tuple in combinations(exceptionals, 4):
        if any(pair(left, right) != 0 for left, right in combinations(exceptional_tuple, 2)):
            continue
        numerator = _sub(_add(*exceptional_tuple), K)
        if any(value % 3 for value in numerator):
            continue
        line: Vector = tuple(value // 3 for value in numerator)  # type: ignore[assignment]
        if square(line) != 1 or any(pair(line, exceptional) != 0 for exceptional in exceptional_tuple):
            continue
        result.append((line, exceptional_tuple))
    return sorted(result)


def _canonical_plane_rows(rows: Sequence[tuple[int, ...]]) -> tuple[tuple[int, ...], ...]:
    return min(
        tuple(
            sorted(
                (row[0], *(row[1 + index] for index in exceptional_permutation))
                for row in rows
            )
        )
        for exceptional_permutation in permutations(range(4))
    )


def blowdown_template_witnesses(
    triangle: Triangle,
    bases: Sequence[tuple[Vector, tuple[Vector, ...]]],
) -> list[dict[str, Any]]:
    sequence = tuple(sorted(int(square(component)) for component in triangle))
    template = _canonical_plane_rows(PLANE_BLOWUP_TEMPLATES[sequence]["rows"])
    witnesses: list[dict[str, Any]] = []
    for line, exceptionals in bases:
        rows = tuple(
            (
                int(pair(line, component)),
                *(int(pair(exceptional, component)) for exceptional in exceptionals),
            )
            for component in triangle
        )
        if _canonical_plane_rows(rows) != template:
            continue
        witnesses.append(
            {
                "line": list(line),
                "exceptionals": [list(vector) for vector in exceptionals],
                "rows": [list(row) for row in rows],
            }
        )
    return witnesses


def _in_integer_span(vector: Vector, basis: Sequence[Vector]) -> bool:
    if not basis:
        return False
    columns = sp.Matrix.hstack(*(_v(item) for item in basis))
    if columns.rank() != sp.Matrix.hstack(columns, _v(vector)).rank():
        return False
    solution = next(iter(sp.linsolve((columns, _v(vector)))))
    return all(sp.Rational(value).q == 1 for value in solution)


def simple_root_systems(roots: Sequence[Vector]) -> list[tuple[Vector, ...]]:
    """Every oriented simple system of every closed subsystem of a rank-2 root set."""

    result: set[tuple[Vector, ...]] = {()}
    result.update((root,) for root in roots)
    for left, right in combinations(roots, 2):
        if _rank((left, right)) != 2 or pair(left, right) not in (0, 1):
            continue
        columns = sp.Matrix.hstack(_v(left), _v(right))
        good = True
        for root in roots:
            if not _in_integer_span(root, (left, right)):
                continue
            coordinates = next(iter(sp.linsolve((columns, _v(root)))))
            if not (
                all(value >= 0 for value in coordinates)
                or all(value <= 0 for value in coordinates)
            ):
                good = False
                break
        if good:
            result.add(tuple(sorted((left, right))))
    return sorted(result)


def _saturated_kernel_basis(triangle: Triangle) -> tuple[Vector, Vector]:
    """Find a primitive Z-basis of D^perp; rank two and a bound-three witness suffices."""

    for bound in range(1, 4):
        vectors = [
            vector
            for vector in product(range(-bound, bound + 1), repeat=5)
            if any(vector) and all(pair(vector, component) == 0 for component in triangle)
        ]
        vectors.sort(key=lambda vector: (-int(square(vector)), sum(abs(value) for value in vector), vector))
        for left, right in combinations(vectors, 2):
            if _gcd_pair_minors(left, right) == 1:
                return left, right  # type: ignore[return-value]
    raise AssertionError("no saturated D^perp basis found in the certified bound")


def _coordinates(vector: Vector, basis: Sequence[Vector]) -> tuple[int, ...]:
    columns = sp.Matrix.hstack(*(_v(item) for item in basis))
    solution = next(iter(sp.linsolve((columns, _v(vector)))))
    if any(sp.Rational(value).q != 1 for value in solution):
        raise AssertionError("nonintegral coordinates in a saturated basis")
    return tuple(int(value) for value in solution)


def period_profile(
    triangle: Triangle,
    roots: Sequence[Vector],
    simple_roots: Sequence[Vector],
    basis: tuple[Vector, Vector] | None = None,
) -> dict[str, Any]:
    if basis is None:
        basis = _saturated_kernel_basis(triangle)
    coordinates = [_coordinates(root, basis) for root in simple_roots]
    rank = len(simple_roots)
    if rank == 0:
        root_pairs = len(roots) // 2
        stratum = (
            "(C*)^2 minus 1 root hypertorus"
            if root_pairs == 1
            else f"(C*)^2 minus {root_pairs} root hypertori"
        )
        quotient = "Z^2"
        primitive = True
        forbidden_exponents: list[int] = []
    elif rank == 1:
        a, b = coordinates[0]
        primitive = gcd(abs(a), abs(b)) == 1
        intended = {simple_roots[0], tuple(-value for value in simple_roots[0])}
        forbidden_exponents = sorted(
            {
                abs(a * y - b * x)
                for root in roots
                if root not in intended
                for x, y in [_coordinates(root, basis)]
            }
            - {0}
        )
        stratum = "C*" if not forbidden_exponents else "C* minus {1}"
        quotient = "Z"
    else:
        coordinate_matrix = sp.Matrix.hstack(*(_v(value) for value in coordinates))
        index = abs(int(coordinate_matrix.det()))
        primitive = index == 1
        forbidden_exponents = []
        stratum = "point" if index == 1 else f"{index} torsion characters"
        quotient = "0" if index == 1 else f"finite group of order {index}"
    return {
        "lambda_basis": [list(vector) for vector in basis],
        "lambda_gram": [
            [int(pair(left, right)) for right in basis]
            for left in basis
        ],
        "simple_coordinates": [list(value) for value in coordinates],
        "root_count": len(roots),
        "root_subsystem_rank": rank,
        "quotient": quotient,
        "root_span_primitive": primitive,
        "forbidden_free_exponents": forbidden_exponents,
        "exact_stratum": stratum,
        "nonempty": primitive,
    }


def complete_negative_walls(
    triangle: Triangle,
    simple_roots: Sequence[Vector],
    exceptionals: Sequence[Vector],
) -> WallSet:
    negative_boundary = [component for component in triangle if square(component) < 0]
    irreducible_exceptionals = [
        exceptional
        for exceptional in exceptionals
        if exceptional in negative_boundary
        or (
            all(pair(exceptional, component) >= 0 for component in negative_boundary)
            and all(pair(exceptional, root) >= 0 for root in simple_roots)
        )
    ]
    return tuple(sorted(set(negative_boundary + list(simple_roots) + irreducible_exceptionals)))


def _monoid_decompositions(walls: WallSet, max_weight: int = 6) -> dict[Vector, tuple[int, ...]]:
    zero: Vector = (0, 0, 0, 0, 0)
    result: dict[Vector, tuple[int, ...]] = {zero: (0,) * len(walls)}
    # Increasing total weight makes the recorded decompositions minimal.
    for weight in range(1, max_weight + 1):
        for indices in combinations_with_replacement(range(len(walls)), weight):
            vector = zero
            coefficients = [0] * len(walls)
            for index in indices:
                vector = _add(vector, walls[index])
                coefficients[index] += 1
            result.setdefault(vector, tuple(coefficients))
    return result


def _coxeter_relation(left: Vector, right: Vector) -> str:
    inner = int(pair(left, right))
    left_norm = -int(square(left))
    right_norm = -int(square(right))
    product_norm = left_norm * right_norm
    if inner == 0:
        return "m2"
    if inner * inner == product_norm:
        return "parallel"
    if inner * inner > product_norm:
        return "ultraparallel"
    if left_norm == right_norm == 2 and inner == 1:
        return "m3"
    if {left_norm, right_norm} == {1, 2} and inner == 1:
        return "m4"
    return "not_coxeter"


def _gram(walls: Sequence[Vector]) -> sp.Matrix:
    return sp.Matrix([[pair(left, right) for right in walls] for left in walls])


def _canonical_gram(walls: WallSet) -> tuple[tuple[int, ...], tuple[int, ...]]:
    gram = _gram(walls)
    norm_groups: dict[int, list[int]] = defaultdict(list)
    for index, wall in enumerate(walls):
        norm_groups[int(square(wall))].append(index)
    best_flat: tuple[int, ...] | None = None
    best_order: tuple[int, ...] | None = None
    for blocks in product(
        *(list(permutations(norm_groups[norm])) for norm in sorted(norm_groups))
    ):
        order = tuple(index for block in blocks for index in block)
        flat = tuple(int(gram[left, right]) for left in order for right in order)
        if best_flat is None or flat < best_flat:
            best_flat = flat
            best_order = order
    assert best_flat is not None and best_order is not None
    return best_flat, best_order


def _automorphism_order(walls: WallSet) -> int:
    gram = _gram(walls)
    norm_groups: dict[int, list[int]] = defaultdict(list)
    for index, wall in enumerate(walls):
        norm_groups[int(square(wall))].append(index)
    count = 0
    for blocks in product(
        *(list(permutations(norm_groups[norm])) for norm in sorted(norm_groups))
    ):
        order = tuple(index for block in blocks for index in block)
        if all(
            gram[left, right] == gram[order[left], order[right]]
            for left in range(len(walls))
            for right in range(len(walls))
        ):
            count += 1
    return count


def _strict_interior(walls: WallSet) -> Vector:
    minus_two = [wall for wall in walls if square(wall) == -2]
    for bound in range(1, 8):
        for perturbation in product(range(-bound, bound + 1), repeat=5):
            if not all(pair(perturbation, wall) > 0 for wall in minus_two):
                continue
            for scale in range(1, 101):
                candidate: Vector = tuple(
                    scale * MINUS_K[index] + perturbation[index]
                    for index in range(5)
                )  # type: ignore[assignment]
                if square(candidate) > 0 and all(pair(candidate, wall) > 0 for wall in walls):
                    return candidate
    raise AssertionError("strict interior search exceeded certified bound")


def _facet_witnesses(walls: WallSet, interior: Vector) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for index, wall in enumerate(walls):
        witness = _v(interior) - pair(interior, wall) / square(wall) * _v(wall)
        values = [sp.Rational((_v(other).T * J * witness)[0]) for other in walls]
        result.append(
            {
                "facet": index,
                "witness": _vector_payload(witness),
                "square": _rational_payload(sp.Rational((witness.T * J * witness)[0])),
                "minimum_other_pairing": _rational_payload(
                    min(value for other_index, value in enumerate(values) if other_index != index)
                ),
            }
        )
    return result


def _extreme_rays(walls: WallSet) -> list[dict[str, Any]]:
    result: dict[Vector, dict[str, Any]] = {}
    for indices in combinations(range(len(walls)), 4):
        equations = sp.Matrix.vstack(*((_v(walls[index]).T * J) for index in indices))
        nullspace = equations.nullspace()
        if len(nullspace) != 1:
            continue
        ray = _primitive(nullspace[0])
        values = [pair(ray, wall) for wall in walls]
        if not all(value >= 0 for value in values):
            ray = tuple(-value for value in ray)  # type: ignore[assignment]
            values = [-value for value in values]
        if not all(value >= 0 for value in values):
            continue
        ray_square = int(square(ray))
        result[ray] = {
            "vector": list(ray),
            "square": ray_square,
            "kind": "finite" if ray_square > 0 else "ideal" if ray_square == 0 else "hyperideal",
            "active_facets": [index for index, value in enumerate(values) if value == 0],
        }
    return [result[ray] for ray in sorted(result)]


def _reflection_integral(wall: Vector) -> bool:
    root = _v(wall)
    reflection = sp.eye(5) - sp.Rational(2, square(wall)) * root * root.T * J
    return (
        all(sp.Rational(value).q == 1 for value in reflection)
        and reflection.T * J * reflection == J
        and reflection * reflection == sp.eye(5)
        and reflection.det() == -1
    )


def _period_label(root_count: int, rank: int) -> str:
    if rank == 0:
        pairs = root_count // 2
        return (
            "(C*)^2 minus 1 root hypertorus"
            if pairs == 1
            else f"(C*)^2 minus {pairs} root hypertori"
        )
    if rank == 1:
        return "C*" if root_count == 2 else "C* minus {1}"
    return "point"


def recognize_raw(walls: Sequence[Sequence[int]], type_catalog: Sequence[dict[str, Any]] | None = None) -> dict[str, Any]:
    """Recognize a raw marked wall system without accepting a supplied surface."""

    normalized: WallSet = tuple(sorted(tuple(int(value) for value in wall) for wall in walls))  # type: ignore[assignment]
    if any(len(wall) != 5 for wall in normalized):
        return {"carrier": False, "reason": "walls_are_not_vectors_in_I_1_4"}
    if any(gcd(*[abs(value) for value in wall]) != 1 for wall in normalized):
        return {"carrier": False, "reason": "nonprimitive_wall"}
    if len(normalized) > 8:
        return {"carrier": False, "reason": "facet_bound_exceeded"}
    if any(square(wall) not in (-1, -2) for wall in normalized):
        return {"carrier": False, "reason": "wall_norm_outside_minus1_minus2"}
    if _rank(normalized) != 5:
        return {"carrier": False, "reason": "walls_do_not_span_H4_lattice"}
    if any(pair(left, right) < 0 for left, right in combinations(normalized, 2)):
        return {"carrier": False, "reason": "wall_orientation_is_not_acute"}
    if any(_coxeter_relation(left, right) == "not_coxeter" for left, right in combinations(normalized, 2)):
        return {"carrier": False, "reason": "non_coxeter_intersection"}
    canonical = _recover_canonical(normalized)
    if canonical is None or any(value.q != 1 for value in canonical):
        return {"carrier": False, "reason": "adjunction_canonical_is_not_integral"}
    canonical_integer: Vector = tuple(int(value) for value in canonical)  # type: ignore[assignment]
    if square(canonical_integer) != 5 or not _is_characteristic(canonical_integer):
        return {"carrier": False, "reason": "adjunction_canonical_fails_characteristic_square_five_gate"}
    if _lattice_index(normalized) != 1:
        return {"carrier": False, "reason": "surviving_full_nef_walls_do_not_generate_I_1_4"}
    flat, _ = _canonical_gram(normalized)
    if type_catalog is None:
        type_catalog = build_certificate()["surviving_gram_types"]
    for row in type_catalog:
        archived_flat = tuple(value for matrix_row in row["canonical_gram"] for value in matrix_row)
        if flat == archived_flat:
            return {
                "carrier": True,
                "type_id": row["type_id"],
                "canonical": list(canonical_integer),
                "surface_strata": row["surface_strata"],
                "decision": row["cone_certificate"]["decision"],
            }
    return {
        "carrier": False,
        "reason": "unlisted_complete_nef_gram",
        "interpretation": "triangular Looijenga carrier obstruction only; not a polytope nonexistence theorem",
    }


def build_certificate() -> dict[str, Any]:
    exceptionals = numerical_exceptional_classes()
    roots = canonical_roots()
    components = boundary_component_candidates()
    triangles = anticanonical_triangles(components)
    blowdown_bases = exceptional_blowdown_bases(exceptionals)

    triangle_sequence_counts: Counter[tuple[int, int, int]] = Counter()
    blowdown_witness_profile: Counter[tuple[tuple[int, int, int], int]] = Counter()
    triangle_data: dict[Triangle, dict[str, Any]] = {}
    all_templates_valid = True
    all_kernel_bases_saturated = True
    for triangle in triangles:
        sequence = tuple(sorted(int(square(component)) for component in triangle))
        triangle_sequence_counts[sequence] += 1
        witnesses = blowdown_template_witnesses(triangle, blowdown_bases)
        blowdown_witness_profile[(sequence, len(witnesses))] += 1
        all_templates_valid &= bool(witnesses)
        boundary_roots = tuple(
            root for root in roots if all(pair(root, component) == 0 for component in triangle)
        )
        kernel_basis = _saturated_kernel_basis(triangle)
        all_kernel_bases_saturated &= _gcd_maximal_minors(
            sp.Matrix.hstack(*(_v(vector) for vector in kernel_basis))
        ) == 1
        triangle_data[triangle] = {
            "sequence": sequence,
            "roots": boundary_roots,
            "kernel_basis": kernel_basis,
            "template_witness_count": len(witnesses),
        }

    branch_profile: Counter[tuple[tuple[int, int, int], int, int]] = Counter()
    excluded_over_facet_bound = 0
    survivor_instances: list[dict[str, Any]] = []
    wall_decorations: dict[WallSet, list[dict[str, Any]]] = defaultdict(list)
    monoid_cache: dict[WallSet, dict[Vector, tuple[int, ...]]] = {}
    maximum_hidden_wall_weight = 0
    hidden_weight_profile: Counter[tuple[str, int]] = Counter()
    all_hidden_inequalities_generated = True
    all_period_strata_nonempty = True
    all_branch_walls_coxeter = True
    all_survivor_wall_lattices_unimodular = True
    all_canonical_vectors_recovered = True

    for triangle in triangles:
        data = triangle_data[triangle]
        boundary_roots: tuple[Vector, ...] = data["roots"]
        for simple_roots in simple_root_systems(boundary_roots):
            walls = complete_negative_walls(triangle, simple_roots, exceptionals)
            branch_profile[(data["sequence"], len(simple_roots), len(walls))] += 1
            if len(walls) > 8:
                excluded_over_facet_bound += 1
                continue
            if _rank(walls) != 5:
                raise AssertionError("a purported full-nef wall family lost full rank")
            period = period_profile(
                triangle,
                boundary_roots,
                simple_roots,
                data["kernel_basis"],
            )
            all_period_strata_nonempty &= bool(period["nonempty"])
            all_branch_walls_coxeter &= all(
                pair(left, right) >= 0 and _coxeter_relation(left, right) != "not_coxeter"
                for left, right in combinations(walls, 2)
            )
            all_survivor_wall_lattices_unimodular &= _lattice_index(walls) == 1
            recovered = _recover_canonical(walls)
            all_canonical_vectors_recovered &= recovered == tuple(sp.Rational(value) for value in K)

            if walls not in monoid_cache:
                monoid_cache[walls] = _monoid_decompositions(walls)
            monoid = monoid_cache[walls]
            for kind, targets in (("exceptional", exceptionals), ("boundary", triangle)):
                for target in targets:
                    if target not in monoid:
                        all_hidden_inequalities_generated = False
                        continue
                    weight = sum(monoid[target])
                    maximum_hidden_wall_weight = max(maximum_hidden_wall_weight, weight)
                    hidden_weight_profile[(kind, weight)] += 1

            decoration = {
                "triangle": [list(component) for component in triangle],
                "boundary_sequence": list(data["sequence"]),
                "simple_roots": [list(root) for root in simple_roots],
                "period_stratum": period["exact_stratum"],
            }
            survivor_instances.append(
                {
                    "walls": walls,
                    "triangle": triangle,
                    "sequence": data["sequence"],
                    "simple_roots": simple_roots,
                    "root_count": len(boundary_roots),
                    "period": period,
                }
            )
            wall_decorations[walls].append(decoration)

    gram_classes: dict[tuple[int, ...], list[WallSet]] = defaultdict(list)
    for walls in wall_decorations:
        flat, _ = _canonical_gram(walls)
        gram_classes[flat].append(walls)

    type_ids = {
        (6, 3): "T6a",
        (6, 2): "T6b",
        (7, 2): "T7",
        (8, 1): "T8",
    }
    type_rows: list[dict[str, Any]] = []
    all_type_facets_irredundant = True
    all_type_rays_nonspacelike = True
    all_reflections_integral = True
    for flat, marked_sets in sorted(
        gram_classes.items(),
        key=lambda item: (len(item[1][0]), -sum(square(wall) == -2 for wall in item[1][0]), item[0]),
    ):
        representative = min(marked_sets)
        _, canonical_order = _canonical_gram(representative)
        canonical_walls: WallSet = tuple(representative[index] for index in canonical_order)
        canonical_gram = _gram(canonical_walls)
        facet_count = len(canonical_walls)
        minus_two_count = sum(square(wall) == -2 for wall in canonical_walls)
        type_id = type_ids[(facet_count, minus_two_count)]
        decorations = [
            decoration
            for wall_set in marked_sets
            for decoration in wall_decorations[wall_set]
        ]
        stratum_profile = Counter(
            (
                tuple(decoration["boundary_sequence"]),
                len(decoration["simple_roots"]),
                decoration["period_stratum"],
            )
            for decoration in decorations
        )
        interior = _strict_interior(canonical_walls)
        facets = _facet_witnesses(canonical_walls, interior)
        rays = _extreme_rays(canonical_walls)
        ray_profile = Counter(ray["kind"] for ray in rays)
        relations = Counter(
            _coxeter_relation(left, right)
            for left, right in combinations(canonical_walls, 2)
        )
        facets_irredundant = all(
            sp.Rational(row["square"]) > 0 and sp.Rational(row["minimum_other_pairing"]) > 0
            for row in facets
        )
        rays_nonspacelike = ray_profile.get("hyperideal", 0) == 0
        reflections_integral = all(_reflection_integral(wall) for wall in canonical_walls)
        all_type_facets_irredundant &= facets_irredundant
        all_type_rays_nonspacelike &= rays_nonspacelike
        all_reflections_integral &= reflections_integral
        type_rows.append(
            {
                "type_id": type_id,
                "facet_count": facet_count,
                "wall_norm_profile": {
                    "minus_1": facet_count - minus_two_count,
                    "minus_2": minus_two_count,
                },
                "marked_wall_set_count": len(marked_sets),
                "decoration_branch_count": len(decorations),
                "representative_walls": [list(wall) for wall in canonical_walls],
                "canonical_gram": _matrix_payload(canonical_gram),
                "weighted_gram_automorphism_order": _automorphism_order(canonical_walls),
                "coxeter_relation_profile": dict(sorted(relations.items())),
                "surface_strata": [
                    {
                        "boundary_sequence": list(sequence),
                        "simple_root_rank": rank,
                        "exact_period_stratum": stratum,
                        "decoration_branch_count": count,
                    }
                    for (sequence, rank, stratum), count in sorted(stratum_profile.items())
                ],
                "cone_certificate": {
                    "strict_interior": list(interior),
                    "facet_witnesses": facets,
                    "extreme_rays": rays,
                    "vertex_count": len(rays),
                    "finite_vertex_count": ray_profile.get("finite", 0),
                    "ideal_vertex_count": ray_profile.get("ideal", 0),
                    "hyperideal_ray_count": ray_profile.get("hyperideal", 0),
                    "decision": "finite_volume_noncompact" if rays_nonspacelike and ray_profile.get("ideal", 0) else "unresolved",
                },
                "checks": {
                    "walls_generate_I_1_4": _lattice_index(canonical_walls) == 1,
                    "canonical_recovers_K": _recover_canonical(canonical_walls) == tuple(sp.Rational(value) for value in K),
                    "all_facets_irredundant": facets_irredundant,
                    "all_extreme_rays_nonspacelike": rays_nonspacelike,
                    "all_reflections_integral": reflections_integral,
                },
            }
        )

    root_type_profile = {
        str(list(sequence)): {
            "triangle_count": triangle_sequence_counts[sequence],
            "D_perp_root_count": len(next(data["roots"] for triangle, data in triangle_data.items() if data["sequence"] == sequence)),
        }
        for sequence in BOUNDARY_SEQUENCE_ORDER
    }
    branch_profile_rows = [
        {
            "boundary_sequence": list(sequence),
            "simple_root_rank": rank,
            "complete_wall_count": wall_count,
            "branch_count": count,
            "inside_facet_bound": wall_count <= 8,
        }
        for (sequence, rank, wall_count), count in sorted(branch_profile.items())
    ]

    checks = {
        "ambient_is_odd_unimodular_I_1_4": J.det() == 1 and K == (-3, 1, 1, 1, 1),
        "canonical_is_characteristic_of_square_5": _is_characteristic(K) and square(K) == 5,
        "canonical_perp_is_negative_definite_A4": len(roots) == 20,
        "numerical_exceptional_count_is_10": len(exceptionals) == 10,
        "bounded_component_count_is_70": len(components) == 70,
        "bounded_triangle_count_is_165": len(triangles) == 165,
        "triangle_sequence_profile_is_60_60_20_15_10": [
            triangle_sequence_counts[sequence] for sequence in BOUNDARY_SEQUENCE_ORDER
        ] == [60, 60, 20, 15, 10],
        "all_triangles_have_plane_or_infinitely_near_blowup_template": all_templates_valid,
        "exceptional_blowdown_basis_count_is_5": len(blowdown_bases) == 5,
        "all_D_perp_bases_are_saturated": all_kernel_bases_saturated,
        "branch_count_is_885": sum(branch_profile.values()) == 885,
        "exactly_25_branches_exceed_eight_facets": excluded_over_facet_bound == 25,
        "surviving_decoration_branch_count_is_860": len(survivor_instances) == 860,
        "surviving_marked_wall_set_count_is_260": len(wall_decorations) == 260,
        "all_period_character_strata_nonempty": all_period_strata_nonempty,
        "all_hidden_exceptional_and_boundary_inequalities_are_wall_sums": all_hidden_inequalities_generated,
        "maximum_hidden_wall_sum_weight_is_6": maximum_hidden_wall_weight == 6,
        "all_surviving_wall_systems_are_crystallographic_coxeter": all_branch_walls_coxeter,
        "all_surviving_wall_systems_generate_I_1_4": all_survivor_wall_lattices_unimodular,
        "canonical_recovered_from_every_surviving_raw_wall_system": all_canonical_vectors_recovered,
        "weighted_gram_type_count_is_4": len(type_rows) == 4,
        "marked_wall_set_profile_is_120_60_60_20": [row["marked_wall_set_count"] for row in type_rows] == [120, 60, 60, 20],
        "all_type_facets_irredundant": all_type_facets_irredundant,
        "all_type_extreme_rays_nonspacelike": all_type_rays_nonspacelike,
        "all_type_reflections_integral": all_reflections_integral,
        "all_four_types_are_finite_volume_noncompact": all(
            row["cone_certificate"]["decision"] == "finite_volume_noncompact" for row in type_rows
        ),
        "every_square_minus_two_class_in_D_perp_is_a_deformation_root": all(
            any(
                root in instance["simple_roots"]
                for instance in survivor_instances
                if instance["triangle"] == triangle
            )
            for triangle, data in triangle_data.items()
            for root in data["roots"]
        ),
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "kind": "surface_free_triangular_looijenga_full_nef_carrier_compiler",
        "status": "finite_exhaustive_iff_in_C_triangle_4_8_with_cited_surface_realization_inputs",
        "class_definition": {
            "name": "C_triangle(4,8)",
            "ambient_lattice": "I_(1,4)",
            "hyperbolic_dimension": 4,
            "maximum_facets": 8,
            "primitive_wall_squares": [-2, -1],
            "requirements": [
                "full rational wall span",
                "acute crystallographic Coxeter intersection alphabet",
                "no supplied surface or supplied boundary decoration",
            ],
        },
        "termination_proof": {
            "canonical_recovery": "full wall span makes K.w=-2-w^2 a unique rational system",
            "negative_definite_shell": (
                "For u in K^perp, write u=(a,b,c,d,-3a-b-c-d).  Then "
                "-u^2=z^T Q z with Q=[[8,3,3,3],[3,2,1,1],[3,1,2,1],[3,1,1,2]] positive definite."
            ),
            "inverse_form": [
                ["4/5", "-3/5", "-3/5", "-3/5"],
                ["-3/5", "6/5", "1/5", "1/5"],
                ["-3/5", "1/5", "6/5", "1/5"],
                ["-3/5", "1/5", "1/5", "6/5"],
            ],
            "component_bound": (
                "For d^2=s in [-2,3] and K.d=-s-2, u=5d-(K.d)K has "
                "u^2=5(-s^2+s-4), so -u^2<=50.  Cauchy in Q gives "
                "|u0|<=6 and |u1|,...,|u4|<=7; the compiler enumerates this exact box."
            ),
            "root_bound": "The same form at -u^2=2 gives the exact 20-root A4 shell.",
            "exceptional_bound": (
                "For e=aH+sum(c_i Ei), e^2=K.e=-1 gives Cauchy inequality "
                "5a^2-6a-3<=0, hence a is 0 or 1 and exactly ten classes occur."
            ),
        },
        "surface_free_census": {
            "canonical": list(K),
            "canonical_square": int(square(K)),
            "root_count_in_K_perp": len(roots),
            "numerical_exceptional_count": len(exceptionals),
            "numerical_exceptionals": [list(vector) for vector in exceptionals],
            "bounded_component_count": len(components),
            "component_square_profile": [
                {"square": component_square, "count": count}
                for component_square, count in sorted(
                    Counter(int(square(component)) for component in components).items()
                )
            ],
            "triangle_count": len(triangles),
            "triangle_sequence_profile": [
                {"boundary_sequence": list(sequence), "count": triangle_sequence_counts[sequence]}
                for sequence in BOUNDARY_SEQUENCE_ORDER
            ],
            "blowdown_basis_count": len(blowdown_bases),
            "blowdown_templates": [
                {
                    "boundary_sequence": list(sequence),
                    "plane_rows": [list(row) for row in PLANE_BLOWUP_TEMPLATES[sequence]["rows"]],
                    "model": PLANE_BLOWUP_TEMPLATES[sequence]["model"],
                    "triangle_count": triangle_sequence_counts[sequence],
                    "witness_count_profile": [
                        {"witnesses_per_triangle": witnesses, "triangle_count": count}
                        for (candidate_sequence, witnesses), count in sorted(blowdown_witness_profile.items())
                        if candidate_sequence == sequence
                    ],
                }
                for sequence in BOUNDARY_SEQUENCE_ORDER
            ],
            "root_type_profile": root_type_profile,
        },
        "branch_census": {
            "all_branch_count": sum(branch_profile.values()),
            "excluded_over_eight_facets": excluded_over_facet_bound,
            "surviving_decoration_branches": len(survivor_instances),
            "surviving_marked_wall_sets": len(wall_decorations),
            "rows": branch_profile_rows,
            "hidden_wall_sum_maximum_weight": maximum_hidden_wall_weight,
            "hidden_wall_sum_weight_profile": [
                {"kind": kind, "weight": weight, "target_count": count}
                for (kind, weight), count in sorted(hidden_weight_profile.items())
            ],
        },
        "surviving_gram_types": type_rows,
        "inverse_recognition_theorem": {
            "statement": (
                "A raw wall system P in C_triangle(4,8) has a marked rational full-nef "
                "Looijenga-pair carrier with an anticanonical triangle iff the adjunction "
                "canonical gate and the finite decoration/period/C_D++ compiler are nonempty. "
                "Equivalently in this census, its integral weighted Gram system is one of "
                "T6a,T6b,T7,T8."
            ),
            "forward_direction": (
                "A carrier supplies K, its triangle, its exact period-root subsystem, and all "
                "negative curves; the finite exhaustion therefore places it in one compiler row."
            ),
            "reverse_direction": (
                "A blowup template constructs the deformation type; the nonempty exact character "
                "stratum and lifted Weyl chamber construct the requested nodal branch; exhaustive "
                "nonnegative wall sums prove equality with C_D++, so Friedman/GHK give an actual "
                "marked pair whose full nef cone is P."
            ),
            "scope_warning": (
                "Failure is a triangular-Looijenga carrier obstruction, never by itself an "
                "abstract hyperbolic-polytope nonexistence theorem.  This is a total detector "
                "for the chosen carrier species, not U_total asserting that every raw P has "
                "a carrier, and it transfers no dimension bound or polytope no-go theorem."
            ),
        },
        "external_theorem_inputs": [
            {
                "source": "Gross--Hacking--Keel, Moduli of surfaces with an anti-canonical cycle",
                "url": "https://arxiv.org/abs/1211.6367",
                "uses": [
                    "toric models after corner blowups (Definition 1.2 and Lemma 1.3)",
                    "root criterion and nodal-versus-period-trivial roots (Theorem 3.1 and Proposition 3.4)",
                    "marked moduli as the lifted period torus (Theorem 6.1)",
                ],
            },
            {
                "source": "Friedman, On the ample cone of a rational surface with an anticanonical cycle",
                "url": "https://arxiv.org/abs/1207.7012",
                "uses": [
                    "period surjectivity (Theorem 1.3)",
                    "complete negative-curve and numerical-exceptional wall formula (Lemmas 2.3--2.4)",
                    "Looijenga-root criterion (Theorem 2.14)",
                ],
            },
            {
                "source": "Friedman, On the geometry of anticanonical pairs",
                "url": "https://arxiv.org/abs/1502.02560",
                "uses": [
                    "period surjectivity (Theorem 3.17)",
                    "root/effectivity criterion (Theorem 6.6)",
                ],
            },
        ],
        "lineage": {
            "T8": "recovers the previously frozen eight-wall H4 inverse Looijenga carrier",
            "T6a_T6b_T7": (
                "project-new positive full-nef Looijenga carrier types and explicit arithmetic "
                "H4 Coxeter realizations; worldwide abstract-polytope novelty is not claimed "
                "without a dedicated lineage search"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()) and all(all(row["checks"].values()) for row in type_rows),
    }
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def verify_archive() -> dict[str, Any]:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    fresh = build_certificate()
    if archived != fresh:
        raise AssertionError("archived triangular Looijenga compiler has drifted")
    if not fresh["all_checks_pass"]:
        failed = [name for name, passed in fresh["checks"].items() if not passed]
        raise AssertionError(f"failed checks: {failed}")
    return {
        "certificate_sha256": fresh["certificate_sha256"],
        "triangle_count": fresh["surface_free_census"]["triangle_count"],
        "decoration_branch_count": fresh["branch_census"]["surviving_decoration_branches"],
        "marked_wall_set_count": fresh["branch_census"]["surviving_marked_wall_sets"],
        "gram_types": [row["type_id"] for row in fresh["surviving_gram_types"]],
        "decisions": [row["cone_certificate"]["decision"] for row in fresh["surviving_gram_types"]],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    if args.check_archive:
        print(json.dumps(verify_archive(), indent=2, sort_keys=True))
        return
    result = build_certificate()
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
