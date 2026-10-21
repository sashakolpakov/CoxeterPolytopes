#!/usr/bin/env python3
"""Certify the full-NS Type-MII Coble chamber in hyperbolic dimension 17.

Kondo's characteristic-three Type-MII surface is the blow-up of
``P1 x P1`` at its sixteen ``F3``-rational grid points.  This script builds
the Picard lattice and every visible curve from that model, enumerates all
primitive cusps of the complete Coble--Mukai chamber, and checks that none
has an integral one-boundary exceptional lift.  It then identifies the
complete 48-wall arrangement, entry by entry and by an integral marked
isometry, with Alexeev's published ``(18,4,0)`` Coxeter chamber.

This is a surface-to-polytope certificate and a marked input for the future
nine-class inverse compiler.  It does not assert an unmarked
polytope-to-Coble existence theorem or novelty of the abstract chamber.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
from typing import Any, Sequence

import sympy as sp
from sympy.matrices.normalforms import (
    hermite_normal_form,
    smith_normal_decomp,
    smith_normal_form,
)
from sympy.polys.domains import ZZ

from scripts.certify_coble_one_boundary_full_ns_h10 import _cdd_rays
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/coble_type_mii_full_ns_h17.json"
ALEXEEV = (
    ROOT
    / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json"
)
SCHEMA = "coxeter.coble_type_mii_full_ns_h17.v1"

PINNED = {
    "kondo_source_tar": "6989648564a80167f2daa35da64b90a063a0a090b9f23ff843b7fabb9d96a98d",
    "kondo_tex": "8e23df4e77a2666a7a4f9a634c1ec2f36066353bf0bdcd0a72b11a0a70bfb137",
    "alexeev_source": "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457",
    "alexeev_archive": "c339b56a6d0b68c0939dab85ba6390c7429a30e56eeafec822ac0d5373ec0e91",
}

# Our wall order is B1,...,B4,B'1,...,B'4,E11,...,E44 and then the
# graph curves C_sigma in lexicographic S4 order.  Entries are zero-based
# positions in Alexeev's ancillary (18,4,0) matrix.
MII_TO_ALEXEEV = (
    0, 29, 5, 13, 11, 3, 15, 7,
    23, 24, 17, 18, 10, 2, 16, 8, 21, 4, 31, 6, 12, 22, 14, 30,
    1, 41, 28, 32, 42, 20, 9, 38, 19, 39, 36, 25,
    34, 40, 26, 46, 37, 47, 44, 27, 43, 35, 33, 45,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("an integral matrix became nonintegral")
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _rational_rows(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(value) if value.q == 1 else str(value) for value in matrix.row(row)]
        for row in range(matrix.rows)
    ]


def _half_integral_basis(
    generators: Sequence[sp.Matrix], rank: int = 18
) -> sp.Matrix:
    raw = sp.Matrix.hstack(*(2 * vector for vector in generators))
    doubled = sp.Matrix(
        [[int(raw[row, column]) for column in range(raw.cols)] for row in range(raw.rows)]
    )
    if doubled.rows != rank or doubled.rank() != rank:
        raise ArithmeticError("the half-integral generators do not have full rank")
    if any(value.q != 1 for value in doubled):
        raise ArithmeticError("a proposed half-integral generator has denominator >2")
    return hermite_normal_form(doubled) / 2


def _surface_model() -> dict[str, Any]:
    """Build Bl_(P1(F3)xP1(F3))(P1xP1) and its 48 negative curves."""

    gram = sp.zeros(18)
    gram[0, 1] = gram[1, 0] = 1
    for index in range(2, 18):
        gram[index, index] = -1
    identity = sp.eye(18)
    f1, f2 = identity[:, 0], identity[:, 1]
    exceptional = {
        (row, column): identity[:, 2 + 4 * row + column]
        for row in range(4)
        for column in range(4)
    }
    boundaries = [
        f1
        - sum(
            (exceptional[row, column] for column in range(4)),
            sp.zeros(18, 1),
        )
        for row in range(4)
    ] + [
        f2
        - sum(
            (exceptional[row, column] for row in range(4)),
            sp.zeros(18, 1),
        )
        for column in range(4)
    ]
    permutations = list(itertools.permutations(range(4)))
    graph_curves = [
        f1
        + f2
        - sum(
            (exceptional[row, permutation[row]] for row in range(4)),
            sp.zeros(18, 1),
        )
        for permutation in permutations
    ]
    exceptionals = [
        exceptional[row, column] for row in range(4) for column in range(4)
    ]
    physical_walls = [*boundaries, *exceptionals, *graph_curves]
    physical_matrix = sp.Matrix.hstack(*physical_walls)
    physical_gram = physical_matrix.T * gram * physical_matrix
    canonical = -2 * f1 - 2 * f2 + sum(exceptionals, sp.zeros(18, 1))
    strict_interior = 5 * f1 + 5 * f2 - sum(exceptionals, sp.zeros(18, 1))

    labels = [*(f"B{index + 1}" for index in range(4))]
    labels.extend(f"B'{index + 1}" for index in range(4))
    labels.extend(
        f"E{row + 1}{column + 1}"
        for row in range(4)
        for column in range(4)
    )
    labels.extend(
        "C_" + "".join(str(value + 1) for value in permutation)
        for permutation in permutations
    )

    relation_histogram: Counter[str] = Counter()
    for left in range(48):
        for right in range(left):
            product = physical_gram[left, right]
            ratio = sp.cancel(
                product**2
                / (physical_gram[left, left] * physical_gram[right, right])
            )
            relation = {
                sp.Rational(0): "m2",
                sp.Rational(1, 4): "m3",
                sp.Rational(1, 2): "m4",
                sp.Rational(1): "parallel",
            }.get(ratio)
            if relation is None:
                raise ArithmeticError(f"a non-Coxeter wall ratio appeared: {ratio}")
            relation_histogram[relation] += 1

    # Orthogonal projection of a strict chamber interior to any wall is a
    # timelike relative-interior point of that facet because the wall Gram
    # is acute.  This gives an exact irredundancy witness for every wall.
    facet_witnesses = []
    for index, wall in enumerate(physical_walls):
        square = (wall.T * gram * wall)[0]
        pairing = (strict_interior.T * gram * wall)[0]
        witness = strict_interior - pairing / square * wall
        witness_pairings = [
            (witness.T * gram * other)[0] for other in physical_walls
        ]
        if not (
            witness_pairings[index] == 0
            and all(
                value > 0
                for other, value in enumerate(witness_pairings)
                if other != index
            )
            and (witness.T * gram * witness)[0] > 0
        ):
            raise ArithmeticError("a wall failed its exact facet witness")
        facet_witnesses.append(witness)

    # Pic~=<Pic,B_i/2>.  Only seven half-boundary classes are independent,
    # since sum B_i/2=-K.  Doubling the form gives an integral reflective
    # lattice L with Nikulin invariants (18,4,0).
    pic_tilde_basis = _half_integral_basis(
        [*identity.columnspace(), *(boundary / 2 for boundary in boundaries)]
    )
    arithmetic_gram = 2 * pic_tilde_basis.T * gram * pic_tilde_basis
    arithmetic_normals = [
        *(boundary / 2 for boundary in boundaries),
        *exceptionals,
        *graph_curves,
    ]
    arithmetic_walls = pic_tilde_basis.inv() * sp.Matrix.hstack(
        *arithmetic_normals
    )
    arithmetic_wall_gram = (
        arithmetic_walls.T * arithmetic_gram * arithmetic_walls
    )
    smith = [
        abs(int(value))
        for value in smith_normal_form(arithmetic_gram, domain=ZZ).diagonal()
    ]
    dual = arithmetic_gram.inv()
    delta = int(any(sp.denom(dual[index, index]) != 1 for index in range(18)))
    wall_span = [
        abs(int(value))
        for value in smith_normal_form(arithmetic_walls, domain=ZZ).diagonal()
    ]
    root_types: Counter[tuple[int, int]] = Counter()
    integral_reflections = []
    for index in range(48):
        wall = arithmetic_walls[:, index]
        square = int((wall.T * arithmetic_gram * wall)[0])
        divisibility = math.gcd(
            *(abs(int(value)) for value in arithmetic_gram * wall)
        )
        root_types[(square, divisibility)] += 1
        reflection = (
            sp.eye(18)
            - sp.Rational(2, square)
            * wall
            * (wall.T * arithmetic_gram)
        )
        integral_reflections.append(
            all(value.q == 1 for value in reflection)
            and reflection.T * arithmetic_gram * reflection == arithmetic_gram
            and reflection**2 == sp.eye(18)
        )

    ns_reflections = []
    for wall in physical_walls:
        square = (wall.T * gram * wall)[0]
        reflection = sp.eye(18) - sp.Rational(2, square) * wall * (wall.T * gram)
        ns_reflections.append(all(value.q == 1 for value in reflection))

    checks = {
        "NS_is_odd_unimodular_I_1_17": gram.det() == -1
        and exact_inertia(gram) == (1, 17, 0)
        and any(int(gram[index, index]) % 2 for index in range(18)),
        "anti_bicanonical_boundary": sum(boundaries, sp.zeros(18, 1))
        == -2 * canonical,
        "eight_disjoint_minus_four_boundaries": sp.Matrix.hstack(
            *boundaries
        ).T
        * gram
        * sp.Matrix.hstack(*boundaries)
        == -4 * sp.eye(8),
        "canonical_square_minus_eight": (canonical.T * gram * canonical)[0] == -8,
        "wall_norm_split": Counter(physical_gram.diagonal())
        == Counter({-4: 8, -1: 16, -2: 24}),
        "physical_walls_span_NS": physical_matrix.rank() == 18
        and [
            abs(int(value))
            for value in smith_normal_form(physical_matrix, domain=ZZ).diagonal()
        ]
        == [1] * 18,
        "acute": all(
            physical_gram[left, right] >= 0
            for left in range(48)
            for right in range(left)
        ),
        "coxeter_relation_histogram": relation_histogram
        == Counter({"m2": 796, "m3": 128, "m4": 96, "parallel": 108}),
        "strict_interior_square_34": (strict_interior.T * gram * strict_interior)[0]
        == 34,
        "strict_interior_pairings": Counter(
            (strict_interior.T * gram * wall)[0] for wall in physical_walls
        )
        == Counter({1: 24, 6: 24}),
        "all_48_facets_irredundant": len(facet_witnesses) == 48,
        "extended_index_128": abs(pic_tilde_basis.det()) == sp.Rational(1, 128),
        "arithmetic_lattice_18_4_0": arithmetic_gram.det() == -16
        and exact_inertia(arithmetic_gram) == (1, 17, 0)
        and all(int(arithmetic_gram[index, index]) % 2 == 0 for index in range(18))
        and smith == [1] * 14 + [2] * 4
        and delta == 0,
        "arithmetic_wall_types": root_types
        == Counter({(-2, 1): 24, (-4, 2): 24}),
        "arithmetic_walls_span_L": wall_span == [1] * 18,
        "all_48_arithmetic_reflections_integral": all(integral_reflections),
        "only_boundaries_fail_on_NS": [
            index for index, value in enumerate(ns_reflections) if not value
        ]
        == list(range(8)),
        "physical_and_arithmetic_Grams_are_same_normalized_arrangement": all(
            physical_gram[left, right] ** 2
            * arithmetic_wall_gram[left, left]
            * arithmetic_wall_gram[right, right]
            == arithmetic_wall_gram[left, right] ** 2
            * physical_gram[left, left]
            * physical_gram[right, right]
            for left in range(48)
            for right in range(48)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Type-MII surface model failed: {checks}")

    return {
        "gram": gram,
        "f1": f1,
        "f2": f2,
        "exceptional": exceptional,
        "boundaries": boundaries,
        "permutations": permutations,
        "graph_curves": graph_curves,
        "physical_walls": physical_walls,
        "physical_gram": physical_gram,
        "canonical": canonical,
        "strict_interior": strict_interior,
        "facet_witnesses": facet_witnesses,
        "labels": labels,
        "relation_histogram": relation_histogram,
        "pic_tilde_basis": pic_tilde_basis,
        "arithmetic_gram": arithmetic_gram,
        "arithmetic_walls": arithmetic_walls,
        "arithmetic_wall_gram": arithmetic_wall_gram,
        "root_types": root_types,
        "checks": checks,
    }


def _cusp_component_type(
    component: Sequence[int], root_gram: sp.Matrix
) -> str:
    if len(component) == 2 and root_gram[component[0], component[1]] == 2:
        return "2*A1~"
    degrees = []
    for vertex in component:
        degrees.append(
            sum(
                root_gram[vertex, other] == 1
                for other in component
                if other != vertex
            )
        )
    if len(component) >= 3 and degrees == [2] * len(component):
        return f"A{len(component) - 1}~"
    raise ArithmeticError("an active cusp component is not affine A type")


def _components(nodes: Sequence[int], root_gram: sp.Matrix) -> list[list[int]]:
    unseen = set(nodes)
    result = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        queue = [start]
        while queue:
            vertex = queue.pop()
            adjacent = {
                other
                for other in unseen
                if root_gram[vertex, other] != 0
            }
            unseen.difference_update(adjacent)
            component.update(adjacent)
            queue.extend(adjacent)
        result.append(sorted(component))
    return result


def _cm_chamber(surface: dict[str, Any], executable: str) -> dict[str, Any]:
    """Build CM(S), enumerate all 202 primitive cusps, and test every lift."""

    gram: sp.Matrix = surface["gram"]
    boundaries: list[sp.Matrix] = surface["boundaries"]
    exceptional: dict[tuple[int, int], sp.Matrix] = surface["exceptional"]
    graph_curves: list[sp.Matrix] = surface["graph_curves"]
    pic_basis: sp.Matrix = surface["pic_tilde_basis"]

    constraint = 2 * sp.Matrix.hstack(*boundaries).T * gram * pic_basis
    if any(value.q != 1 for value in constraint):
        raise ArithmeticError("the CM orthogonality constraint became nonintegral")
    diagonal, left_transform, right_transform = smith_normal_decomp(
        constraint, domain=ZZ
    )
    rank = constraint.rank()
    cm_basis = pic_basis * right_transform[:, rank:]
    cm_gram = cm_basis.T * gram * cm_basis

    cm_roots = []
    root_labels = []
    for row in range(4):
        for column in range(4):
            cm_roots.append(
                boundaries[row] / 2
                + boundaries[4 + column] / 2
                + 2 * exceptional[row, column]
            )
            root_labels.append(f"e{row + 1}{column + 1}")
    cm_roots.extend(graph_curves)
    root_labels.extend(surface["labels"][24:])
    root_matrix = sp.Matrix.hstack(*cm_roots)
    root_coordinates = cm_basis.gauss_jordan_solve(root_matrix)[0]
    root_gram = root_matrix.T * gram * root_matrix
    root_span = [
        abs(int(value))
        for value in smith_normal_form(root_coordinates, domain=ZZ).diagonal()
    ]

    inequalities = root_coordinates.T * cm_gram
    rays = sorted(_cdd_rays(inequalities, executable))
    counts: Counter[str] = Counter()
    cusps = []
    for ray in rays:
        vector = sp.Matrix(ray)
        square = (vector.T * cm_gram * vector)[0]
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        counts[kind] += 1
        if square == 0:
            cusps.append(tuple(int(value) for value in ray))

    cusp_records = []
    for cusp in cusps:
        fiber = cm_basis * sp.Matrix(cusp)
        pairings = inequalities * sp.Matrix(cusp)
        active = [index for index, value in enumerate(pairings) if value == 0]
        membership = []
        for boundary in boundaries:
            candidate = fiber - boundary / 2
            membership.append(all(value.q == 1 for value in candidate))
        cusp_records.append(
            {
                "CM_coordinates": list(cusp),
                "active_root_indices_zero_based": active,
                "active_root_labels": [root_labels[index] for index in active],
                "single_boundary_integral_lift_membership": membership,
            }
        )

    # The natural S4 x S4 and factor-swap action has four cusp orbits.  This
    # independently recovers exactly Kondo's four maximal parabolic types.
    identity4 = tuple(range(4))

    def ambient_action(
        row_permutation: tuple[int, ...] = identity4,
        column_permutation: tuple[int, ...] = identity4,
        transpose: bool = False,
    ) -> sp.Matrix:
        action = sp.zeros(18)
        if not transpose:
            action[:, 0] = surface["f1"]
            action[:, 1] = surface["f2"]
            for row in range(4):
                for column in range(4):
                    action[:, 2 + 4 * row + column] = exceptional[
                        row_permutation[row], column_permutation[column]
                    ]
        else:
            action[:, 0] = surface["f2"]
            action[:, 1] = surface["f1"]
            for row in range(4):
                for column in range(4):
                    action[:, 2 + 4 * row + column] = exceptional[
                        column_permutation[column], row_permutation[row]
                    ]
        return action

    ambient_generators = []
    for index in range(3):
        adjacent = list(identity4)
        adjacent[index], adjacent[index + 1] = adjacent[index + 1], adjacent[index]
        ambient_generators.append(ambient_action(tuple(adjacent), identity4))
        ambient_generators.append(ambient_action(identity4, tuple(adjacent)))
    ambient_generators.append(ambient_action(transpose=True))
    cm_generators = [
        cm_basis.gauss_jordan_solve(action * cm_basis)[0]
        for action in ambient_generators
    ]
    if any(
        any(value.q != 1 for value in generator) or abs(generator.det()) != 1
        for generator in cm_generators
    ):
        raise ArithmeticError("a grid symmetry failed to preserve CM integrally")

    lookup = {cusp: index for index, cusp in enumerate(cusps)}
    unseen = set(range(len(cusps)))
    cusp_orbits = []
    while unseen:
        seed = min(unseen)
        orbit = {seed}
        queue: deque[int] = deque([seed])
        while queue:
            current = queue.popleft()
            vector = sp.Matrix(cusps[current])
            for generator in cm_generators:
                image = tuple(int(value) for value in generator * vector)
                target = lookup.get(image)
                if target is None:
                    raise ArithmeticError("a grid symmetry left the cusp set")
                if target not in orbit:
                    orbit.add(target)
                    queue.append(target)
        unseen.difference_update(orbit)
        representative = min(orbit)
        active = cusp_records[representative]["active_root_indices_zero_based"]
        components = _components(active, root_gram)
        types = sorted(
            (_cusp_component_type(component, root_gram) for component in components),
            key=lambda value: (len(value), value),
        )
        cusp_orbits.append(
            {
                "size": len(orbit),
                "representative_cusp_index_zero_based": representative,
                "affine_components": types,
                "active_root_count": len(active),
            }
        )
    cusp_orbits.sort(key=lambda record: (record["size"], record["affine_components"]))
    orbit_signature = {
        (record["size"], tuple(record["affine_components"]))
        for record in cusp_orbits
    }
    expected_orbits = {
        (72, ("A7~", "2*A1~")),
        (96, ("A2~", "A5~", "2*A1~")),
        (18, ("A3~", "A3~", "2*A1~", "2*A1~")),
        (16, ("A2~", "A2~", "A2~", "A2~")),
    }

    checks = {
        "constraint_smith_diagonal": [
            abs(int(diagonal[index, index])) for index in range(8)
        ]
        == [2, 2, 2, 2, 2, 2, 2, 4],
        "CM_rank_ten": cm_basis.shape == (18, 10),
        "CM_is_E10": cm_gram.det() == -1
        and exact_inertia(cm_gram) == (1, 9, 0)
        and all(int(cm_gram[index, index]) % 2 == 0 for index in range(10)),
        "forty_roots_norm_minus_two": root_gram.rank() == 10
        and Counter(root_gram.diagonal()) == Counter({-2: 40}),
        "forty_roots_span_CM": all(value.q == 1 for value in root_coordinates)
        and root_span == [1] * 10,
        "CM_exact_cone_counts": counts
        == Counter({"finite": 4_752, "ideal": 202, "hyperideal": 0}),
        "all_202_primitive_cusps_recorded": len(cusp_records) == 202,
        "no_single_boundary_cusp_lifts": all(
            not any(record["single_boundary_integral_lift_membership"])
            for record in cusp_records
        ),
        "four_Kondo_cusp_orbits": orbit_signature == expected_orbits,
        "grid_symmetry_order": 2 * math.factorial(4) ** 2 == 1_152,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Type-MII CM chamber failed: {checks}")

    return {
        "cm_basis": cm_basis,
        "cm_gram": cm_gram,
        "root_coordinates": root_coordinates,
        "root_gram": root_gram,
        "root_labels": root_labels,
        "ray_count": len(rays),
        "finite_vertex_count": counts["finite"],
        "ideal_vertex_count": counts["ideal"],
        "hyperideal_vertex_count": counts["hyperideal"],
        "primitive_ray_list_sha256": hashlib.sha256(
            json.dumps(rays, separators=(",", ":")).encode()
        ).hexdigest(),
        "cusps": cusp_records,
        "cusp_orbits": cusp_orbits,
        "checks": checks,
    }


def _lineage(surface: dict[str, Any]) -> dict[str, Any]:
    """Give an entrywise and integral marked match to Alexeev (18,4,0)."""

    if _sha256(ALEXEEV) != PINNED["alexeev_archive"]:
        raise ArithmeticError("the pinned Alexeev reconstruction archive changed")
    archive = json.loads(ALEXEEV.read_text(encoding="utf-8"))
    record = next(
        row
        for row in archive["inventory"]["examples"]
        if row["key"] == [18, 4, 0]
    )
    alexeev_gram = sp.Matrix(
        [[sp.Rational(value) for value in row] for row in record["ambient"]["gram"]]
    )
    alexeev_walls = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in record["walls"]["coordinates_in_leading_root_basis"]
        ]
    )
    alexeev_wall_gram = alexeev_walls * alexeev_gram * alexeev_walls.T
    our_wall_gram: sp.Matrix = surface["arithmetic_wall_gram"]
    if not all(
        our_wall_gram[left, right]
        == alexeev_wall_gram[MII_TO_ALEXEEV[left], MII_TO_ALEXEEV[right]]
        for left in range(48)
        for right in range(48)
    ):
        raise ArithmeticError("the MII and Alexeev weighted Grams differ")

    inverse = [0] * 48
    for source, target in enumerate(MII_TO_ALEXEEV):
        inverse[target] = source
    our_walls: sp.Matrix = surface["arithmetic_walls"]
    marking = sp.Matrix.hstack(*(our_walls[:, inverse[index]] for index in range(18)))
    all_wall_marking = all(
        marking * alexeev_walls[index, :].T == our_walls[:, inverse[index]]
        for index in range(48)
    )
    checks = {
        "Alexeev_archive_green": archive["all_checks_pass"],
        "Alexeev_48_roots": record["walls"]["count"] == 48,
        "Alexeev_lattice_18_4_0": record["ambient"]["determinant"] == -16
        and record["ambient"]["inertia"] == [1, 17, 0],
        "all_2304_weighted_Gram_entries_equal": True,
        "integral_marking_unimodular": marking.det() == 1,
        "marked_lattice_Gram_equal": marking.T
        * surface["arithmetic_gram"]
        * marking
        == alexeev_gram,
        "all_48_marked_roots_equal": all_wall_marking,
        "weighted_diagram_group_order_1152": record["weighted_diagram_group"]["order"]
        == 1_152,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Type-MII lineage failed: {checks}")
    return {
        "exact_identification": "Alexeev/Keum--Kondo (18,4,0), H17/48",
        "MII_to_Alexeev_zero_based_wall_permutation": list(MII_TO_ALEXEEV),
        "Alexeev_basis_to_MII_arithmetic_basis": _rows(marking),
        "published_Vinberg_data": {
            "complete_simple_root_count": 48,
            "Vinberg_algorithm_termination_steps": 2,
            "maximal_parabolic_subdiagram_count": 5_244,
            "maximal_parabolic_orbit_count": 17,
            "diagram_group": "C2 semidirect (S4 x S4)",
            "diagram_group_order": 1_152,
            "orthogonal_group_decomposition": "O+(L)=Aut(K4,4) semidirect W_r",
            "decision": "finite_volume_noncompact",
        },
        "lineage": [
            "Keum--Kondo (2001) described the same (18,4,0) Coxeter diagram using Borcherds/Leech methods.",
            "Alexeev (arXiv:2209.09110v4) gave the K4,4^(2)+perfect-matching description, proved reflectivity, and enumerated all maximal parabolics.",
            "Kondo (arXiv:2102.00591v2) constructed the characteristic-three Type-MII Coble surface and its eight boundaries, sixteen exceptional curves, twenty-four (1,1)-curves, and forty complete CM roots.",
            "The present certificate supplies the explicit full-NS marked isometry and exhaustive no-extra-cusp-lift audit; it does not claim a new abstract chamber.",
        ],
        "priority": {
            "new_abstract_polytope": False,
            "new_arithmetic_lattice": False,
            "worldwide_new_surface_realization": "not claimed",
            "project_gain": (
                "first exact full-nef Type-MII/H17 bridge in the project and "
                "the positive top stratum of the nine-class Coble program"
            ),
        },
        "checks": checks,
    }


def build_certificate(cdd_executable: str | None = None) -> dict[str, Any]:
    executable = cdd_executable or shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for the exact CM cusp enumeration")
    surface = _surface_model()
    cm = _cm_chamber(surface, executable)
    lineage = _lineage(surface)

    # Negative-curve exhaustion is now literal.  Kondo proves that the 40 CM
    # roots are all effective irreducible roots.  Sixteen are the boundary-pair
    # (-1)-roots and twenty-four are actual (-2)-curves.  Any other (-1)-curve
    # has B.E=2 by adjunction; if it meets one boundary twice then
    # f=E+B_i/2 is a primitive CM cusp.  All 202 cusps and all eight integral
    # lifts were tested above and none survives.
    negative_curve_exhaustion = {
        "Kondo_complete_effective_irreducible_CM_roots": 40,
        "minus_one_roots_meeting_two_boundaries": 16,
        "actual_minus_two_curves": 24,
        "primitive_CM_cusps_tested": 202,
        "boundary_tests_per_cusp": 8,
        "total_single_boundary_lift_tests": 1_616,
        "integral_single_boundary_lifts": 0,
        "boundary_curves": 8,
        "complete_negative_curve_count": 48,
        "argument": [
            "Kondo: the only negative curves on a Coble surface are boundaries, (-1)-curves, and (-2)-curves.",
            "Adjunction gives B.E=2 for every (-1)-curve.",
            "The (1,1) boundary-incidence case is exactly one of the sixteen complete (-1)-roots.",
            "The (2) case gives the primitive isotropic CM class f=E+B_i/2 and hence an extreme ideal ray of the CM chamber.",
            "Exact cdd enumeration gives all 202 primitive CM cusps; no f-B_i/2 belongs to NS in any of 1,616 tests.",
            "Therefore the eight boundaries, sixteen exceptional curves, and twenty-four graph curves are all negative curves and cut out the full nef chamber.",
        ],
        "full_nef_conclusion": True,
    }

    volume = {
        "decision": "finite_volume_noncompact",
        "method_1_structural": (
            "complete negative-curve exhaustion identifies the chamber with the "
            "full nef cone; Hodge index places every nef ray in the closed "
            "future cone"
        ),
        "method_2_exact_Vinberg_lineage": (
            "the integral marked isometry identifies all 48 simple roots with "
            "Alexeev's complete (18,4,0) diagram; O+(L)=Aut(K4,4) semidirect "
            "W_r proves finite covolume, and 5,244 maximal parabolic diagrams "
            "in 17 symmetry orbits prove noncompactness"
        ),
        "finite_index_of_reflection_group_in_O_plus": 1_152,
        "maximal_parabolic_subdiagram_count": 5_244,
        "maximal_parabolic_orbit_count": 17,
        "compact": False,
        "finite_volume": True,
        "raw_full_vertex_enumeration_required": False,
    }

    checks = {
        "surface_model": all(surface["checks"].values()),
        "CM_cusp_exhaustion": all(cm["checks"].values()),
        "complete_negative_curves": negative_curve_exhaustion["full_nef_conclusion"],
        "full_H17_Coxeter_chamber": len(surface["physical_walls"]) == 48,
        "all_facets_irredundant": surface["checks"]["all_48_facets_irredundant"],
        "arithmetic_reflection_lattice": surface["checks"][
            "all_48_arithmetic_reflections_integral"
        ],
        "finite_volume_noncompact": volume["finite_volume"] and not volume["compact"],
        "exact_known_chamber_lineage": all(lineage["checks"].values()),
        "not_new_abstract_polytope": not lineage["priority"]["new_abstract_polytope"],
        "no_inverse_theorem_overclaim": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global Type-MII H17 checks failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "direction": "SURFACE_TO_POLYTOPE_with_marked_inverse_input_only",
        "claim": {
            "theorem": (
                "The projectivized full Neron--Severi nef cone of Kondo's "
                "unique characteristic-three Type-MII Coble surface is the "
                "48-facet finite-volume noncompact Alexeev/Keum--Kondo "
                "(18,4,0) arithmetic Coxeter polytope in H17."
            ),
            "not_claimed": [
                "a new abstract Coxeter polytope",
                "a new arithmetic reflection lattice",
                "worldwide priority for the full-NS interpretation",
                "an unmarked polytope-to-Coble existence theorem",
                "completion of the nine-class raw-P inverse compiler",
            ],
        },
        "AG_native_stratum": {
            "classification": "Kondo Type MII",
            "characteristic": 3,
            "isomorphism_classes": 1,
            "boundary_component_count": 8,
            "K_square": -8,
            "Picard_rank": 18,
            "hyperbolic_dimension": 17,
            "covering_K3": "supersingular Kum(E x E), Artin invariant 1",
            "automorphism_group": "(S4 x S4) semidirect C2",
            "automorphism_group_order": 1_152,
            "R_invariant": "(D8 + A2^2, (Z/2)^2)",
            "Coble_class_endpoint": (
                "unique b=8 top stratum; it realizes the sharp AG-native "
                "dimension endpoint H17 inside the nine finite-Aut Coble classes"
            ),
        },
        "sources": {
            "kondo": {
                "title": "Coble surfaces with finite automorphism group",
                "arxiv": "2102.00591v2",
                "source_tar_sha256": PINNED["kondo_source_tar"],
                "main_tex_sha256": PINNED["kondo_tex"],
                "used": [
                    "Theorem 1.1/Table 2",
                    "Propositions NefCone, Vinberg2, Halphen",
                    "Example TypeMII, Lemma parabolicMII, Proposition MII",
                    "Section 6 Type MII and uniqueness",
                ],
            },
            "alexeev": {
                "title": "Reflective hyperbolic 2-elementary lattices, K3 surfaces and hyperkahler varieties",
                "arxiv": "2209.09110v4",
                "ancillary_sha256": PINNED["alexeev_source"],
                "archive": str(ALEXEEV.relative_to(ROOT)),
                "archive_sha256": PINNED["alexeev_archive"],
                "used": [
                    "Theorems 1 and 2",
                    "section (18,4,0)",
                    "complete 48-root ancillary Gram",
                    "5,244 maximal parabolics in 17 diagram-group orbits",
                ],
            },
            "keum_kondo": {
                "title": "The automorphism groups of Kummer surfaces associated with the product of two elliptic curves",
                "journal": "Trans. Amer. Math. Soc. 353 (2001), 1469-1487",
                "doi": "10.1090/S0002-9947-00-02631-3",
                "scope": "prior Borcherds--Leech description of the same (18,4,0) diagram",
            },
        },
        "full_NS": {
            "basis": ["F1", "F2", *(f"E{row + 1}{column + 1}" for row in range(4) for column in range(4))],
            "gram": _rows(surface["gram"]),
            "canonical": _rows(surface["canonical"].T)[0],
            "boundary_vectors": _rows(sp.Matrix.hstack(*surface["boundaries"]).T),
            "visible_curves_span_NS_with_index": 1,
            "identification": "odd unimodular I_(1,17)",
        },
        "curve_exhaustion": negative_curve_exhaustion,
        "Coble_Mukai_chamber": {
            "basis_in_NS_coordinates": _rational_rows(cm["cm_basis"]),
            "gram": _rows(cm["cm_gram"]),
            "identification": "E10 = II_(1,9)",
            "root_labels": cm["root_labels"],
            "root_vectors_in_CM_basis": _rows(cm["root_coordinates"].T),
            "root_gram": _rows(cm["root_gram"]),
            "exact_double_description": {
                "extreme_ray_count": cm["ray_count"],
                "finite_vertex_count": cm["finite_vertex_count"],
                "ideal_vertex_count": cm["ideal_vertex_count"],
                "hyperideal_vertex_count": cm["hyperideal_vertex_count"],
                "primitive_ray_list_sha256": cm["primitive_ray_list_sha256"],
            },
            "primitive_cusps": cm["cusps"],
            "cusp_orbits_under_C2_semidirect_S4xS4": cm["cusp_orbits"],
        },
        "facets": {
            "labels": surface["labels"],
            "physical_NS_wall_vectors": _rows(
                sp.Matrix.hstack(*surface["physical_walls"]).T
            ),
            "physical_wall_gram": _rows(surface["physical_gram"]),
            "norm_histogram": {"-4": 8, "-2": 24, "-1": 16},
            "relation_histogram": dict(sorted(surface["relation_histogram"].items())),
            "strict_interior": _rows(surface["strict_interior"].T)[0],
            "strict_interior_square": 34,
            "facet_witnesses": _rational_rows(
                sp.Matrix.hstack(*surface["facet_witnesses"]).T
            ),
            "all_facets_irredundant": True,
        },
        "arithmetic": {
            "lattice": "L=2*Pic~(S), Nikulin (18,4,0) = U + D8(-1)^2",
            "basis_in_NS_coordinates": _rational_rows(surface["pic_tilde_basis"]),
            "gram": _rows(surface["arithmetic_gram"]),
            "wall_vectors": _rows(surface["arithmetic_walls"].T),
            "wall_gram": _rows(surface["arithmetic_wall_gram"]),
            "primitive_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(surface["root_types"].items())
            },
            "walls_span_L": True,
            "all_reflections_integral": True,
            "boundary_reflections_fail_on_NS_but_are_integral_on_L": True,
        },
        "finite_volume": volume,
        "lineage": lineage,
        "inverse_program_status": {
            "marked_boundary_fold_available": True,
            "raw_P_exhaustive_compiler_available": False,
            "inverse_theorem_credit": False,
            "use": (
                "positive top-stratum template and exact return target for the "
                "future nine-class Coble P->surface compiler"
            ),
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--cdd", default=None)
    args = parser.parse_args()
    certificate = build_certificate(args.cdd)
    if args.check_archive:
        if not args.output.exists() or json.loads(
            args.output.read_text(encoding="utf-8")
        ) != certificate:
            raise SystemExit("archived Type-MII H17 certificate differs")
        print("archived Type-MII H17 certificate matches")
        return
    args.output.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(args.output)


if __name__ == "__main__":
    main()
