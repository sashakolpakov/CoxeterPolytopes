#!/usr/bin/env python3
"""Exact integral-triality census for Bredon's 24 ADEG polytopes.

The input is the Coxeter-diagram table in N. Bredon, *On ADEG-polyhedra in
hyperbolic spaces*, arXiv:2507.05153v1, Theorem 4.2 and Table 2.  Every
unlabelled edge has label three and every labelled edge has label six.

An integral ``D4``-triality unfolding must assign square 2 or square 6 to
each wall.  Label-three edges preserve the length and label-six edges switch
it.  Hence, after contracting all label-three components, the label-six
quotient must be loop-free and bipartite.  For each passing diagram this
script constructs both long/short orientations, unfolds every square-six
wall into three orthogonal square-two walls, and computes the exact rank,
inertia, saturated wall lattice and discriminant data.

This is a carrier *precompiler*.  A K3 return additionally requires a
primitive K3 embedding and a complete-root/no-hidden-wall proof.  The script
does not promote a passing unfolding to a nef-cone realization.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/bredon_adeg24_triality_census.json"
SCHEMA = "coxeter.bredon_adeg24_triality_census.v1"
SOURCE = {
    "author": "Naomi Bredon",
    "title": "On ADEG-polyhedra in hyperbolic spaces",
    "arxiv": "2507.05153v1",
    "location": "Theorem 4.2 and Table 2",
}


Edge = tuple[int, int]


@dataclass(frozen=True)
class Diagram:
    id: str
    dimension: int
    vertices: int
    m3: tuple[Edge, ...]
    m6: tuple[Edge, ...]
    source_order: int


def path_edges(vertices: int) -> tuple[Edge, ...]:
    return tuple((index, index + 1) for index in range(vertices - 1))


def labelled_path(
    identifier: str, dimension: int, labels: tuple[int, ...], source_order: int
) -> Diagram:
    edges = path_edges(len(labels) + 1)
    return Diagram(
        identifier,
        dimension,
        len(labels) + 1,
        tuple(edge for edge, label in zip(edges, labels, strict=True) if label == 3),
        tuple(edge for edge, label in zip(edges, labels, strict=True) if label == 6),
        source_order,
    )


DIAGRAMS: tuple[Diagram, ...] = (
    # H^2: four ADEG triangles.
    Diagram("ADEG-N2-01-path-66", 2, 3, (), ((0, 1), (1, 2)), 1),
    Diagram("ADEG-N2-02-triangle-633", 2, 3, ((1, 2), (2, 0)), ((0, 1),), 2),
    Diagram("ADEG-N2-03-triangle-663", 2, 3, ((2, 0),), ((0, 1), (1, 2)), 3),
    Diagram("ADEG-N2-04-triangle-666", 2, 3, (), ((0, 1), (1, 2), (2, 0)), 4),
    # H^3: seven ADEG tetrahedra.
    labelled_path("ADEG-N3-01-path-633", 3, (6, 3, 3), 5),
    labelled_path("ADEG-N3-02-path-636", 3, (6, 3, 6), 6),
    labelled_path("ADEG-N3-03-path-363", 3, (3, 6, 3), 7),
    Diagram("ADEG-N3-04-star-6", 3, 4, ((1, 2), (1, 3)), ((0, 1),), 8),
    Diagram(
        "ADEG-N3-05-triangle-leaf-6",
        3,
        4,
        ((1, 2), (1, 3), (2, 3)),
        ((0, 1),),
        9,
    ),
    Diagram(
        "ADEG-N3-06-square-one-6",
        3,
        4,
        ((1, 2), (2, 3), (3, 0)),
        ((0, 1),),
        10,
    ),
    Diagram(
        "ADEG-N3-07-square-opposite-6",
        3,
        4,
        ((1, 2), (3, 0)),
        ((0, 1), (2, 3)),
        11,
    ),
    # H^5: two pyramids followed by the two Im Hof doubly-truncated simplices.
    labelled_path("ADEG-N5-01-path-633336", 5, (6, 3, 3, 3, 3, 6), 12),
    Diagram(
        "ADEG-N5-02-path-triangle",
        5,
        7,
        ((1, 2), (2, 3), (3, 4), (4, 5), (4, 6), (5, 6)),
        ((0, 1),),
        13,
    ),
    Diagram(
        "ADEG-N5-03-octagon-alternating",
        5,
        8,
        ((1, 2), (3, 4), (5, 6), (7, 0)),
        ((0, 1), (2, 3), (4, 5), (6, 7)),
        14,
    ),
    Diagram(
        "ADEG-N5-04-octagon-clustered",
        5,
        8,
        ((1, 2), (4, 5), (6, 7), (7, 0)),
        ((0, 1), (2, 3), (3, 4), (5, 6)),
        15,
    ),
    # H^6: the unique ADEG pyramid.
    Diagram(
        "ADEG-N6-01-diamond-tail",
        6,
        8,
        ((0, 1), (1, 2), (2, 3), (3, 0), (2, 4), (4, 5), (5, 6)),
        ((6, 7),),
        16,
    ),
    # H^7: five ADEG pyramids.
    Diagram(
        "ADEG-N7-01-pentagon-tail",
        7,
        9,
        ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (1, 5), (5, 6), (6, 7)),
        ((7, 8),),
        17,
    ),
    Diagram(
        "ADEG-N7-02-path-star",
        7,
        9,
        ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 7), (5, 8)),
        ((0, 1),),
        18,
    ),
    Diagram(
        "ADEG-N7-03-three-armed-path",
        7,
        10,
        ((1, 2), (2, 3), (3, 4), (4, 5), (3, 7), (7, 8)),
        ((0, 1), (5, 6), (8, 9)),
        19,
    ),
    Diagram(
        "ADEG-N7-04-triangle-branch-path",
        7,
        10,
        ((1, 2), (2, 3), (3, 4), (4, 5), (3, 7), (7, 8), (7, 9), (8, 9)),
        ((0, 1), (5, 6)),
        20,
    ),
    Diagram(
        "ADEG-N7-05-two-triangles-branch",
        7,
        10,
        (
            (0, 1), (1, 2), (2, 0), (1, 3), (3, 4),
            (4, 5), (4, 6), (5, 6), (3, 7), (7, 8),
        ),
        ((8, 9),),
        21,
    ),
    # H^9: one pyramid and Bredon's exceptional P_star.
    Diagram(
        "ADEG-N9-01-branched-path",
        9,
        11,
        ((1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (6, 9), (9, 10)),
        ((0, 1),),
        22,
    ),
    Diagram(
        "ADEG-N9-02-P-star",
        9,
        14,
        tuple(
            edge
            for a, b, c in ((0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 10, 11))
            for edge in ((12, a), (b, c), (c, 13))
        ),
        ((0, 1), (3, 4), (6, 7), (9, 10)),
        23,
    ),
    # H^11: the unique ADEG pyramid.
    Diagram(
        "ADEG-N11-01-branched-path",
        11,
        13,
        (
            (0, 1), (1, 2), (2, 3), (3, 4), (5, 1), (4, 6),
            (6, 8), (8, 7), (0, 9), (7, 10), (10, 11),
        ),
        ((11, 12),),
        24,
    ),
)


def connected_components(vertices: int, edges: Iterable[Edge]) -> tuple[tuple[int, ...], ...]:
    adjacency = [set() for _ in range(vertices)]
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    unseen = set(range(vertices))
    components: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        stack = [start]
        component: set[int] = set()
        while stack:
            current = stack.pop()
            if current in component:
                continue
            component.add(current)
            stack.extend(adjacency[current] - component)
        unseen -= component
        components.append(tuple(sorted(component)))
    return tuple(components)


def triality_orientations(diagram: Diagram) -> tuple[dict[str, Any], list[set[int]]]:
    components = connected_components(diagram.vertices, diagram.m3)
    component_of = {
        vertex: index
        for index, component in enumerate(components)
        for vertex in component
    }
    quotient = [set() for _ in components]
    loops: list[Edge] = []
    for left, right in diagram.m6:
        a, b = component_of[left], component_of[right]
        if a == b:
            loops.append((left, right))
        else:
            quotient[a].add(b)
            quotient[b].add(a)
    colors: dict[int, int] = {}
    odd_cycle = False
    for start in range(len(components)):
        if start in colors:
            continue
        colors[start] = 0
        stack = [start]
        while stack:
            current = stack.pop()
            for neighbour in quotient[current]:
                expected = 1 - colors[current]
                if neighbour in colors and colors[neighbour] != expected:
                    odd_cycle = True
                elif neighbour not in colors:
                    colors[neighbour] = expected
                    stack.append(neighbour)
    passes = not loops and not odd_cycle
    orientations: list[set[int]] = []
    if passes:
        # Every ADEG Coxeter graph is connected, so the bipartition has two
        # global long/short orientations.  Assert this rather than silently
        # undercounting a disconnected source transcription.
        all_edges = diagram.m3 + diagram.m6
        if len(connected_components(diagram.vertices, all_edges)) != 1:
            raise AssertionError(f"{diagram.id}: source diagram is disconnected")
        for long_color in (0, 1):
            orientations.append(
                {
                    vertex
                    for vertex in range(diagram.vertices)
                    if colors[component_of[vertex]] == long_color
                }
            )
    return (
        {
            "passes": passes,
            "m3_components": [list(component) for component in components],
            "m6_loops_after_contraction": [list(edge) for edge in loops],
            "m6_quotient_bipartite": not odd_cycle,
            "criterion": "contract m3 components; m6 quotient is loop-free and bipartite",
        },
        orientations,
    )


def integral_gram(diagram: Diagram, long_vertices: set[int]) -> sp.Matrix:
    gram = sp.zeros(diagram.vertices)
    for vertex in range(diagram.vertices):
        gram[vertex, vertex] = 6 if vertex in long_vertices else 2
    for left, right in diagram.m3:
        if (left in long_vertices) != (right in long_vertices):
            raise AssertionError("m3 edge changed root length")
        value = -3 if left in long_vertices else -1
        gram[left, right] = gram[right, left] = value
    for left, right in diagram.m6:
        if (left in long_vertices) == (right in long_vertices):
            raise AssertionError("m6 edge did not switch root length")
        gram[left, right] = gram[right, left] = -3
    return gram


def normalized_symmetric_gram(diagram: Diagram) -> sp.Matrix:
    """Return twice Bredon's normalized Gram over ``Q(sqrt(3))``."""

    gram = 2 * sp.eye(diagram.vertices)
    for left, right in diagram.m3:
        gram[left, right] = gram[right, left] = -1
    for left, right in diagram.m6:
        gram[left, right] = gram[right, left] = -sp.sqrt(3)
    return gram


def exact_quadratic_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Congruence inertia with exact signs in ``Q(sqrt(3))``."""

    work = sp.Matrix(matrix)
    positive = negative = 0
    while work.rows:
        pivot = next(
            (index for index in range(work.rows) if sp.simplify(work[index, index]) != 0),
            None,
        )
        if pivot is None:
            pair = next(
                (
                    (left, right)
                    for left in range(work.rows)
                    for right in range(left + 1, work.rows)
                    if sp.simplify(work[left, right]) != 0
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            left, right = pair
            order = [left, right] + [
                index for index in range(work.rows) if index not in pair
            ]
            work = work.extract(order, order)
            change = sp.eye(work.rows)
            change[1, 0] = 1
            work = change.T * work * change
            pivot = 0
        if pivot:
            order = [pivot] + [index for index in range(work.rows) if index != pivot]
            work = work.extract(order, order)
        diagonal = sp.simplify(work[0, 0])
        sign = sp.sign(diagonal)
        if sign == 1:
            positive += 1
        elif sign == -1:
            negative += 1
        else:
            raise AssertionError(f"undecided exact Q(sqrt(3)) sign: {diagonal}")
        if work.rows == 1:
            return positive, negative, 0
        column = work[1:, 0]
        work = (work[1:, 1:] - column * column.T / diagonal).applyfunc(sp.simplify)
    return positive, negative, 0


def unfold(diagram: Diagram, long_vertices: set[int], folded: sp.Matrix) -> sp.Matrix:
    vertices = [
        (vertex, phase)
        for vertex in range(diagram.vertices)
        for phase in (range(3) if vertex in long_vertices else range(1))
    ]
    position = {vertex: index for index, vertex in enumerate(vertices)}
    gram = 2 * sp.eye(len(vertices))
    for left, right in diagram.m3:
        if left in long_vertices:
            for phase in range(3):
                i, j = position[left, phase], position[right, phase]
                gram[i, j] = gram[j, i] = -1
        else:
            i, j = position[left, 0], position[right, 0]
            gram[i, j] = gram[j, i] = -1
    for left, right in diagram.m6:
        long = left if left in long_vertices else right
        short = right if left in long_vertices else left
        for phase in range(3):
            i, j = position[long, phase], position[short, 0]
            gram[i, j] = gram[j, i] = -1
    orbit_sum = sp.zeros(diagram.vertices, len(vertices))
    for vertex in range(diagram.vertices):
        phases = range(3) if vertex in long_vertices else range(1)
        for phase in phases:
            orbit_sum[vertex, position[vertex, phase]] = 1
    if orbit_sum * gram * orbit_sum.T != folded:
        raise AssertionError(f"{diagram.id}: orbit sums do not recover folded Gram")
    return gram


def independent_principal_indices(gram: sp.Matrix) -> tuple[int, ...]:
    rank = gram.rank()
    # The RREF pivot set usually works and avoids a combinatorial search.
    pivots = tuple(gram.rref()[1])
    if len(pivots) == rank and gram.extract(pivots, pivots).det() != 0:
        return pivots
    for indices in itertools.combinations(range(gram.rows), rank):
        if gram.extract(indices, indices).det() != 0:
            return indices
    raise AssertionError("symmetric Gram has no full-rank principal minor")


def saturated_wall_lattice(gram: sp.Matrix) -> tuple[sp.Matrix, dict[str, Any]]:
    indices = independent_principal_indices(gram)
    basis_gram = gram.extract(indices, indices)
    coordinates = gram[:, indices] * basis_gram.inv()
    denominator = sp.ilcm(*[sp.denom(value) for value in coordinates])
    integer_coordinates = (denominator * coordinates).applyfunc(int)
    hnf_columns = hermite_normal_form(integer_coordinates.T)
    if hnf_columns.shape != (len(indices), len(indices)):
        raise AssertionError("wall coordinate lattice did not have full rank")
    basis_change = hnf_columns.T / denominator
    saturated_gram = sp.simplify(basis_change * basis_gram * basis_change.T)
    if any(value.q != 1 for value in saturated_gram):
        raise AssertionError("the root-generated quotient lattice is not integral")
    saturated_gram = saturated_gram.applyfunc(int)
    if any(saturated_gram[i, i] % 2 for i in range(saturated_gram.rows)):
        raise AssertionError("the unfolded root lattice is not even")
    snf = smith_normal_form(saturated_gram, domain=ZZ)
    smith = [abs(int(snf[i, i])) for i in range(snf.rows)]
    discriminant_length = sum(value > 1 for value in smith)
    return saturated_gram, {
        "independent_root_indices": list(indices),
        "root_coordinate_denominator": int(denominator),
        "basis_change_from_independent_roots": [
            [str(value) for value in row] for row in basis_change.tolist()
        ],
        "determinant": int(saturated_gram.det()),
        "smith_diagonal": smith,
        "discriminant_length": discriminant_length,
    }


def analyze() -> dict[str, Any]:
    if len(DIAGRAMS) != 24 or [row.source_order for row in DIAGRAMS] != list(range(1, 25)):
        raise AssertionError("the primary 24-row table is not represented exactly once")
    rows: list[dict[str, Any]] = []
    for diagram in DIAGRAMS:
        edges = diagram.m3 + diagram.m6
        if any(left == right for left, right in edges):
            raise AssertionError(f"{diagram.id}: loop in source graph")
        if len({tuple(sorted(edge)) for edge in edges}) != len(edges):
            raise AssertionError(f"{diagram.id}: duplicate source edge")
        normalized = normalized_symmetric_gram(diagram)
        normalized_inertia = exact_quadratic_inertia(normalized)
        expected_normalized = (
            diagram.dimension,
            1,
            diagram.vertices - diagram.dimension - 1,
        )
        if normalized.rank() != diagram.dimension + 1 or normalized_inertia != expected_normalized:
            raise AssertionError(
                f"{diagram.id}: primary-table transcription has rank/inertia "
                f"{normalized.rank(), normalized_inertia}, expected "
                f"{diagram.dimension + 1, expected_normalized}"
            )
        criterion, orientations = triality_orientations(diagram)
        orientation_rows = []
        for orientation_index, long_vertices in enumerate(orientations):
            folded = integral_gram(diagram, long_vertices)
            folded_inertia = exact_inertia(folded)
            expected_folded = (
                diagram.dimension,
                1,
                diagram.vertices - diagram.dimension - 1,
            )
            if folded_inertia != expected_folded:
                raise AssertionError(
                    f"{diagram.id}: transcription has folded inertia {folded_inertia}, "
                    f"expected {expected_folded}"
                )
            unfolded = unfold(diagram, long_vertices, folded)
            unfolded_inertia = exact_inertia(unfolded)
            lattice_gram, lattice = saturated_wall_lattice(unfolded)
            lattice_rank = lattice_gram.rows
            unfolded_lorentzian = unfolded_inertia == (
                lattice_rank - 1,
                1,
                unfolded.rows - lattice_rank,
            )
            gap = 22 - lattice_rank
            stable_k3 = (
                unfolded_lorentzian
                and
                lattice_rank <= 19
                and gap >= int(lattice["discriminant_length"]) + 2
            )
            orientation_rows.append(
                {
                    "orientation": orientation_index,
                    "long_vertices": sorted(long_vertices),
                    "short_vertices": sorted(set(range(diagram.vertices)) - long_vertices),
                    "folded_integral_gram": [
                        [int(value) for value in row] for row in folded.tolist()
                    ],
                    "folded_rank_inertia": [folded.rank(), list(folded_inertia)],
                    "unfolded_root_count": unfolded.rows,
                    "unfolded_rank_inertia": [unfolded.rank(), list(unfolded_inertia)],
                    "unfolded_lorentzian": unfolded_lorentzian,
                    "unfolded_gram": [
                        [int(value) for value in row] for row in unfolded.tolist()
                    ],
                    "root_generated_lattice": lattice,
                    "nikulin_stable_range_primitive_k3_embedding": stable_k3,
                    "strong_same_P_k3_return_certified": False,
                    "strong_return_missing_gate": (
                        "primitive embedding outside stable range and complete-root/no-hidden-wall equality"
                        if not stable_k3
                        else "complete-root/no-hidden-wall equality"
                    ),
                }
            )
        rows.append(
            {
                "id": diagram.id,
                "source_table_order": diagram.source_order,
                "dimension": diagram.dimension,
                "facet_count": diagram.vertices,
                "m3_edges": [list(edge) for edge in diagram.m3],
                "m6_edges": [list(edge) for edge in diagram.m6],
                "folded_gram_rank_expected": diagram.dimension + 1,
                "normalized_symmetric_rank_inertia": [
                    normalized.rank(),
                    list(normalized_inertia),
                ],
                "triality_criterion": criterion,
                "integral_triality_orientation_count": len(orientation_rows),
                "orientations": orientation_rows,
            }
        )
    passing = [row for row in rows if row["triality_criterion"]["passes"]]
    high = [row for row in passing if row["dimension"] >= 5]
    stable = [
        (row["id"], orientation["orientation"])
        for row in passing
        for orientation in row["orientations"]
        if orientation["nikulin_stable_range_primitive_k3_embedding"]
    ]
    payload = {
        "schema": SCHEMA,
        "source": SOURCE,
        "scope": (
            "all 24 published ADEG diagrams; exact integral square-2/square-6 "
            "D4-triality precompiler only, not a complete K3-nef return theorem"
        ),
        "counts": {
            "published_diagrams": len(rows),
            "integral_triality_pass": len(passing),
            "integral_triality_fail": len(rows) - len(passing),
            "passing_dimension_at_least_5": len(high),
            "orientation_packages": sum(len(row["orientations"]) for row in rows),
            "stable_range_k3_embedding_packages": len(stable),
            "strong_same_P_k3_returns": 0,
        },
        "passing_ids": [row["id"] for row in passing],
        "failing_ids": [row["id"] for row in rows if row not in passing],
        "stable_range_packages": [list(item) for item in stable],
        "rows": rows,
        "status": (
            "exact finite carrier precompiler; positive K3-nef claims await "
            "complete-root and hidden-wall certification"
        ),
    }
    payload["canonical_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    arguments = parser.parse_args()
    payload = analyze()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(encoded, encoding="utf-8")
    if arguments.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != encoded:
            raise SystemExit("archive does not match exact replay")
    if not arguments.quiet:
        print(encoded, end="")


if __name__ == "__main__":
    main()
