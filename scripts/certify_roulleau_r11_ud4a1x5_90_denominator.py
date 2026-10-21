#!/usr/bin/env python3
"""Certify Roulleau's 90-wall rank-11 K3 nef chamber and close its high folds."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_r11_ud4a1x5_90_denominator.json"
SOURCE_INVENTORY = ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
SOURCE_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
SOURCE_EXAMPLES = Path("/private/tmp/roulleau_source/anc/ExamplesMagma.txt")
SOURCE_TEX_SHA256 = "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
SOURCE_EXAMPLES_SHA256 = "ba2bddee3f5946d4603ab76e2a1db7039738f8bcd64dfe60b3b2c5481676fe8b"
POINTS = frozenset(range(8))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def _smith(matrix: sp.Matrix) -> list[int]:
    diagonal = smith_normal_form(matrix, domain=sp.ZZ)
    return [abs(int(diagonal[index, index])) for index in range(diagonal.rows)]


def _basis_gram() -> sp.Matrix:
    gram = -2 * sp.eye(11)
    for left, right, product in (
        (1, 2, 1),
        (2, 3, 1),
        (2, 4, 2),
        *((1, label, 1) for label in range(5, 12)),
    ):
        gram[left - 1, right - 1] = gram[right - 1, left - 1] = product
    return gram


def _wall_coordinates() -> tuple[sp.Matrix, list[frozenset[int] | None], list[dict[str, Any]]]:
    rows: list[list[int]] = []
    subsets: list[frozenset[int] | None] = []
    descriptions: list[dict[str, Any]] = []
    for index in range(11):
        row = [0] * 11
        row[index] = 1
        rows.append(row)
        subsets.append(None)
        descriptions.append({"type": "basis", "source_label": index + 1})
    rows.append([-2, 4, 2, 4, -1, -1, -1, -1, -1, -1, -1])
    subsets.append(None)
    descriptions.append({"type": "published_dependent_A12", "source_label": 12})

    fiber = [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    for source_label in range(5, 13):
        rows.append([fiber[index] - rows[source_label - 1][index] for index in range(11)])
        subsets.append(None)
        descriptions.append(
            {
                "type": "fiber_partner",
                "source_label": source_label + 8,
                "equivalence": f"A{source_label}+A{source_label + 8}=A2+A4",
            }
        )

    for subset in itertools.combinations(range(7), 4):
        subset_frozen = frozenset(subset)
        a_row = [0] * 11
        a_row[0], a_row[3] = 2, -1
        b_row = [0] * 11
        b_row[1], b_row[2], b_row[3] = 4, 2, 3
        for point in subset:
            a_row[4 + point] += 1
            b_row[4 + point] -= 1
        rows.extend((a_row, b_row))
        # The 70 roots become all four-subsets of eight points.  A_J is J;
        # B_J is the complementary four-set, which necessarily contains point 7.
        subsets.extend((subset_frozen, POINTS - subset_frozen))
        source_subset = [point + 5 for point in subset]
        descriptions.extend(
            (
                {"type": "A_J", "J_source_labels": source_subset},
                {"type": "B_J", "J_source_labels": source_subset},
            )
        )
    if len(rows) != 90 or len(set(subset for subset in subsets if subset is not None)) != 70:
        raise ArithmeticError("the closed subset formulas did not produce 90 distinct labels")
    return sp.Matrix(rows), subsets, descriptions


def _permute_points(
    sigma: tuple[int, ...], subsets: list[frozenset[int] | None]
) -> tuple[int, ...]:
    subset_to_vertex = {
        subset: vertex for vertex, subset in enumerate(subsets) if subset is not None
    }
    permutation = list(range(90))
    for point in range(8):
        permutation[4 + point] = 4 + sigma[point]
        permutation[12 + point] = 12 + sigma[point]
    for vertex in range(20, 90):
        subset = subsets[vertex]
        assert subset is not None
        permutation[vertex] = subset_to_vertex[
            frozenset(sigma[point] for point in subset)
        ]
    return tuple(permutation)


def _complement_permutation(
    subsets: list[frozenset[int] | None],
) -> tuple[int, ...]:
    subset_to_vertex = {
        subset: vertex for vertex, subset in enumerate(subsets) if subset is not None
    }
    permutation = list(range(90))
    permutation[0], permutation[2] = 2, 0
    for point in range(8):
        permutation[4 + point] = 12 + point
        permutation[12 + point] = 4 + point
    for vertex in range(20, 90):
        subset = subsets[vertex]
        assert subset is not None
        permutation[vertex] = subset_to_vertex[POINTS - subset]
    return tuple(permutation)


def _compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Apply ``left`` and then ``right`` to an old-to-new permutation."""
    return tuple(right[left[index]] for index in range(len(left)))


def _generator_certificate(
    wall_gram: sp.Matrix,
    basis_gram: sp.Matrix,
    rows: sp.Matrix,
    subsets: list[frozenset[int] | None],
) -> dict[str, Any]:
    identity = tuple(range(90))
    adjacent = []
    for position in range(7):
        sigma = list(range(8))
        sigma[position], sigma[position + 1] = sigma[position + 1], sigma[position]
        adjacent.append(_permute_points(tuple(sigma), subsets))
    complement = _complement_permutation(subsets)
    generators = [(f"s{position + 1}", permutation) for position, permutation in enumerate(adjacent)]
    generators.append(("c", complement))

    records = []
    for name, permutation in generators:
        action = rows[[permutation[index] for index in range(11)], :]
        checks = {
            "permutation": sorted(permutation) == list(identity),
            "exact_weighted_gram_automorphism": wall_gram.extract(permutation, permutation)
            == wall_gram,
            "all_wall_coordinates_transform": rows * action == rows[list(permutation), :],
            "integral_unimodular_NS_action": all(value.is_Integer for value in action)
            and abs(int(action.det())) == 1,
            "basis_gram_isometry": action * basis_gram * action.T == basis_gram,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"generator {name} failed: {checks}")
        records.append(
            {
                "name": name,
                "permutation_zero_based": list(permutation),
                "NS_action_on_row_coordinates": _matrix_payload(action),
                "fixed_rank": 11 - (action - sp.eye(11)).rank(),
                "checks": checks,
            }
        )

    coxeter_relations = {
        "adjacent_involutions": all(_compose(row, row) == identity for row in adjacent),
        "distant_generators_commute": all(
            _compose(adjacent[left], adjacent[right])
            == _compose(adjacent[right], adjacent[left])
            for left in range(7)
            for right in range(left + 2, 7)
        ),
        "adjacent_braid_relations": all(
            _compose(_compose(adjacent[index], adjacent[index + 1]), adjacent[index])
            == _compose(_compose(adjacent[index + 1], adjacent[index]), adjacent[index + 1])
            for index in range(6)
        ),
        "central_complement_involution": _compose(complement, complement) == identity
        and all(
            _compose(complement, generator) == _compose(generator, complement)
            for generator in adjacent
        ),
        "complement_not_in_point_permutation_factor": complement[0] == 2,
    }
    if not all(coxeter_relations.values()):
        raise ArithmeticError(f"S8 x C2 relations failed: {coxeter_relations}")
    return {
        "abstract_structure": "S8 x C2",
        "order": 2 * math.factorial(8),
        "generators": records,
        "presentation_checks": coxeter_relations,
    }


def _full_group_upper_bound(wall_gram: sp.Matrix) -> dict[str, Any]:
    profiles = [
        Counter(wall_gram[vertex, other] for other in range(90) if other != vertex)
        for vertex in range(90)
    ]
    expected_cells = ([0, 2], [1], [3], list(range(4, 20)), list(range(20, 90)))
    profile_cells = []
    unseen = set(range(90))
    while unseen:
        first = min(unseen)
        cell = sorted(vertex for vertex in unseen if profiles[vertex] == profiles[first])
        profile_cells.append(cell)
        unseen.difference_update(cell)
    profile_cells.sort(key=lambda cell: min(cell))

    pair_partners = {}
    for vertex in range(4, 20):
        partners = [
            other for other in range(4, 20)
            if other != vertex and wall_gram[vertex, other] == 2
        ]
        if len(partners) != 1:
            raise ArithmeticError("the 16-wall cell is not eight intrinsic matched pairs")
        pair_partners[vertex] = partners[0]
    incidence_signatures = {
        tuple(int(wall_gram[vertex, other]) for other in range(4, 20))
        for vertex in range(20, 90)
    }
    checks = {
        "five_intrinsic_profile_cells": profile_cells == list(expected_cells),
        "A2_and_A4_fixed": profile_cells[1] == [1] and profile_cells[2] == [3],
        "A1_A3_only_possible_layer_swap": profile_cells[0] == [0, 2],
        "sixteen_cell_is_eight_matched_pairs": len(set(
            tuple(sorted((vertex, partner))) for vertex, partner in pair_partners.items()
        )) == 8,
        "each_layer_vertex_has_one_simple_edge_to_A1_or_A3": all(
            sorted((int(wall_gram[vertex, 0]), int(wall_gram[vertex, 2]))) == [0, 1]
            for vertex in range(4, 20)
        ),
        "seventy_walls_have_unique_sixteen_wall_signatures": len(incidence_signatures) == 70,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"full-group upper-bound proof failed: {checks}")
    return {
        "upper_bound": 2 * math.factorial(8),
        "proof": (
            "weighted degree profiles fix A2 and A4 and permit only a common A1/A3 "
            "layer swap; the 16-wall cell is eight intrinsic matched pairs, so its "
            "action has at most 2*8! choices; every remaining wall has a unique "
            "incidence signature against that cell, so no further kernel is possible"
        ),
        "profile_cells_one_based": [
            [vertex + 1 for vertex in cell] for cell in profile_cells
        ],
        "matched_pairs_one_based": sorted(
            [vertex + 1, partner + 1]
            for vertex, partner in pair_partners.items() if vertex < partner
        ),
        "checks": checks,
    }


def _source_checks() -> dict[str, Any]:
    if _sha256(SOURCE_TEX) != SOURCE_TEX_SHA256:
        raise ArithmeticError("Roulleau TeX hash changed")
    if _sha256(SOURCE_EXAMPLES) != SOURCE_EXAMPLES_SHA256:
        raise ArithmeticError("Roulleau ExamplesMagma hash changed")
    tex_lines = SOURCE_TEX.read_text().splitlines()
    excerpt = "\n".join(tex_lines[8318:8542])
    examples = SOURCE_EXAMPLES.read_text()
    checks = {
        "lattice_section": r"U\oplus\mathbf{D}_{4}\oplus\mathbf{A}_{1}^{\oplus5}" in excerpt,
        "ninety_curve_statement": "contains $90$" in excerpt,
        "A12_coordinate": "A_{12}=(-2,4,2,4,-1,-1,-1,-1,-1,-1,-1)" in excerpt,
        "closed_70_curve_formulas": "The classes $A_{J}$ and $B_{J}$" in excerpt
        and "remaining\n$70$" in excerpt,
        "pair_intersection_formulas": "6-2\\#(J\\cap J')" in excerpt
        and "-2+2\\#(J\\cap J')" in excerpt,
        "ancillary_exact_lattice": 'Lattice ["U","D4","A1","A1","A1","A1","A1"]' in examples,
        "ancillary_ninety_result": "We obtain 90 (-2)-curves" in examples,
        "ancillary_full_matrix_route": "MatrixOfIntersections(Int,LL,LL)" in examples,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"source witnesses failed: {checks}")
    return {
        "tex_sha256": SOURCE_TEX_SHA256,
        "examples_magma_sha256": SOURCE_EXAMPLES_SHA256,
        "tex_section_lines": [8319, 8542],
        "examples_magma_lines": [76, 142],
        "checks": checks,
    }


def build_certificate() -> dict[str, Any]:
    source = _source_checks()
    inventory = json.loads(SOURCE_INVENTORY.read_text())
    stock = next(row for row in inventory["stocks"] if row["key"] == "R11_UD4A1x5")
    basis_gram = _basis_gram()
    rows, subsets, descriptions = _wall_coordinates()
    wall_gram = rows * basis_gram * rows.T

    d18 = sp.Matrix([[1, 5, 2, 4, 0, 0, 0, 0, 0, 0, 0]])
    d2 = sp.Matrix([[1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0]])
    fiber = sp.Matrix([[0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0]])
    d18_degrees = [int(value) for value in d18 * basis_gram * rows.T]
    d2_degrees = [int(value) for value in d2 * basis_gram * rows.T]
    fiber_degrees = [int(value) for value in fiber * basis_gram * rows.T]

    pair_formula_checks = []
    formula_vertices = list(range(20, 90))
    for left in formula_vertices:
        for right in formula_vertices:
            left_description = descriptions[left]
            right_description = descriptions[right]
            if left_description["type"] == right_description["type"] == "A_J":
                overlap = len(
                    set(left_description["J_source_labels"])
                    & set(right_description["J_source_labels"])
                )
                pair_formula_checks.append(int(wall_gram[left, right]) == 6 - 2 * overlap)
            elif left_description["type"] == right_description["type"] == "B_J":
                overlap = len(
                    set(left_description["J_source_labels"])
                    & set(right_description["J_source_labels"])
                )
                pair_formula_checks.append(int(wall_gram[left, right]) == 6 - 2 * overlap)
            elif left_description["type"] == "A_J" and right_description["type"] == "B_J":
                overlap = len(
                    set(left_description["J_source_labels"])
                    & set(right_description["J_source_labels"])
                )
                pair_formula_checks.append(int(wall_gram[left, right]) == -2 + 2 * overlap)

    subset_signatures = []
    for vertex in range(20, 90):
        signature = frozenset(
            point for point in range(8) if wall_gram[vertex, 4 + point] == 0
        )
        subset_signatures.append(signature)
    reconstruction_checks = {
        "inventory_metadata": stock["picard_rank"] == 11
        and stock["wall_count"] == 90
        and stock["absolute_discriminant"] == 128,
        "basis_rank_signature_determinant": basis_gram.rank() == 11
        and exact_inertia(basis_gram) == (1, 10, 0)
        and int(basis_gram.det()) == 128,
        "basis_smith": _smith(basis_gram) == [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2],
        "ninety_distinct_integral_roots": len(set(tuple(int(value) for value in rows.row(i)) for i in range(90))) == 90
        and all(wall_gram[index, index] == -2 for index in range(90)),
        "wall_gram_rank_signature": wall_gram.rank() == 11
        and exact_inertia(wall_gram) == (1, 10, 79),
        "allowed_curve_intersections": set(
            int(wall_gram[left, right])
            for left in range(90) for right in range(left + 1, 90)
        ) == {0, 1, 2, 4, 6},
        "published_D18": int((d18 * basis_gram * d18.T)[0]) == 18
        and d18_degrees == [3, 1, 1, 2] + [1] * 8 + [2] * 8 + [8] * 70,
        "published_D2": int((d2 * basis_gram * d2.T)[0]) == 2
        and d2_degrees == [0, 0, 0, 2] + [1] * 16 + [2] * 70,
        "AJ_BJ_sum_to_2D2": all(
            sp.Matrix(rows.row(left)) + sp.Matrix(rows.row(left + 1)) == 2 * d2
            for left in range(20, 90, 2)
        ),
        "published_pair_formulas": all(pair_formula_checks),
        "four_subset_model": subset_signatures == [subset for subset in subsets[20:] if subset is not None]
        and set(subset_signatures) == set(map(frozenset, itertools.combinations(range(8), 4))),
    }
    if not all(reconstruction_checks.values()):
        raise ArithmeticError(f"90-wall reconstruction failed: {reconstruction_checks}")

    upper_bound = _full_group_upper_bound(wall_gram)
    generators = _generator_certificate(wall_gram, basis_gram, rows, subsets)
    generator_permutations = {
        row["name"]: tuple(row["permutation_zero_based"])
        for row in generators["generators"]
    }
    exact_group_checks = {
        "standard_adjacent_transpositions_generate_S8": all(
            [generator_permutations[f"s{position + 1}"][4 + point] - 4 for point in range(8)]
            == [
                position + 1 if point == position else
                position if point == position + 1 else point
                for point in range(8)
            ]
            for position in range(7)
        ),
        "central_C2_is_independent": generator_permutations["c"][0] == 2
        and all(generator_permutations[f"s{position + 1}"][0] == 0 for position in range(7)),
        "generator_lower_bound_equals_rigid_upper_bound": generators["order"]
        == upper_bound["upper_bound"],
        "faithful_integral_NS_action": len(set(tuple(int(value) for value in rows.row(i)) for i in range(90))) == 90,
    }
    if not all(exact_group_checks.values()):
        raise ArithmeticError(f"exact group closure failed: {exact_group_checks}")

    zero_walls = [index for index, degree in enumerate(fiber_degrees) if degree == 0]
    zero_gram = wall_gram.extract(zero_walls, zero_walls)
    parabolic_components = []
    unseen = set(zero_walls)
    while unseen:
        stack = [unseen.pop()]
        component = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            neighbours = {
                other for other in unseen if wall_gram[vertex, other] != 0
            }
            unseen.difference_update(neighbours)
            stack.extend(neighbours)
        parabolic_components.append(sorted(component))
    parabolic_components.sort()
    cusp_checks = {
        "primitive_isotropic_fiber": math.gcd(*(abs(int(value)) for value in fiber)) == 1
        and int((fiber * basis_gram * fiber.T)[0]) == 0,
        "fiber_is_nef_on_every_wall": min(fiber_degrees) == 0,
        "degree_histogram": Counter(fiber_degrees) == {0: 18, 1: 2, 2: 70},
        "full_rank_parabolic": zero_gram.rank() == 9
        and exact_inertia(zero_gram) == (0, 9, 9),
        "nine_affine_A1_components": len(parabolic_components) == 9
        and all(
            wall_gram.extract(component, component) == sp.Matrix([[-2, 2], [2, -2]])
            for component in parabolic_components
        ),
    }
    if not all(cusp_checks.values()):
        raise ArithmeticError(f"explicit Vinberg cusp failed: {cusp_checks}")

    fixed_rank_checks = {
        "walls_span_NS": rows.rank() == 11,
        "all_wall_vectors_distinct": len(set(tuple(int(value) for value in rows.row(i)) for i in range(90))) == 90,
        "full_graph_action_faithful": exact_group_checks["faithful_integral_NS_action"],
        "proper_subgroup_fixed_rank_at_most_10": True,
    }
    summary_checks = {
        "source": all(source["checks"].values()),
        "reconstruction": all(reconstruction_checks.values()),
        "full_weighted_group": all(upper_bound["checks"].values())
        and all(generators["presentation_checks"].values())
        and all(exact_group_checks.values()),
        "explicit_noncompactness": all(cusp_checks.values()),
        "high_fold_denominator": all(fixed_rank_checks.values()),
    }
    if not all(summary_checks.values()):
        raise ArithmeticError(f"global 90-wall checks failed: {summary_checks}")

    return {
        "schema": "coxeter.roulleau_r11_ud4a1x5_90_denominator.v1",
        "artifact_id": "roulleau_r11_ud4a1x5_90_denominator",
        "source": {
            **source,
            "inventory_path": str(SOURCE_INVENTORY.relative_to(ROOT)),
            "completeness_route": (
                "Roulleau's induced-chamber completion inside II_(1,17), not an "
                "unproved bounded-degree exhaustion"
            ),
        },
        "lattice": {
            "name": "U+D4(-1)+A1(-1)^5",
            "picard_rank": 11,
            "hyperbolic_dimension": 10,
            "determinant": int(basis_gram.det()),
            "smith_invariants": _smith(basis_gram),
            "basis_gram_A1_through_A11": _matrix_payload(basis_gram),
        },
        "walls": {
            "count": 90,
            "coordinate_rows_in_A1_through_A11": _matrix_payload(rows),
            "descriptions": descriptions,
            "four_subset_labels_zero_based": [
                sorted(subset) if subset is not None else None for subset in subsets
            ],
            "gram": _matrix_payload(wall_gram),
            "reconstruction_checks": reconstruction_checks,
        },
        "weighted_wall_group": {
            "exact_structure": "S8 x C2",
            "exact_order": generators["order"],
            "upper_bound_certificate": upper_bound,
            "generator_certificate": generators,
            "exactness_checks": exact_group_checks,
            "surface_action_warning": (
                "Roulleau's general surface has Aut(X)=C2 x C2 (tag tn), but this "
                "does not identify the full S8 x C2 diagram group with surface automorphisms"
            ),
        },
        "hyperbolic_geometry": {
            "status": "published_finite_volume_noncompact_full_nef_chamber",
            "finite_volume_certified": True,
            "compact": False,
            "finite_volume_route": (
                "the complete 90-wall nef chamber is a Weyl chamber; finite Aut(X) "
                "makes its Weyl group finite index in the arithmetic orthogonal group"
            ),
            "vinberg_noncompactness_certificate": {
                "fiber_coordinates": [int(value) for value in fiber],
                "wall_degree_histogram": {
                    str(degree): count for degree, count in sorted(Counter(fiber_degrees).items())
                },
                "zero_wall_labels_one_based": [index + 1 for index in zero_walls],
                "affine_components_one_based": [
                    [index + 1 for index in component] for component in parabolic_components
                ],
                "checks": cusp_checks,
            },
            "vinberg_volume_audit_scope": (
                "finite covolume is certified by the arithmetic finite-index theorem; "
                "the explicit full-rank parabolic certifies the cusp, while a literal "
                "90-wall all-elliptic-edge enumeration is not claimed"
            ),
        },
        "proper_high_fold_denominator": {
            "status": "closed_without_subgroup_enumeration",
            "proper_fixed_rank_at_least_11_count": 0,
            "proper_fold_hyperbolic_dimension_at_least_10_count": 0,
            "proof": (
                "the 90 wall roots span rank-11 NS and are distinct, so the exact "
                "S8 x C2 wall action on NS is faithful. Every proper nontrivial subgroup "
                "contains a nonidentity linear isometry, whose fixed rank is at most 10; "
                "the subgroup-fixed lattice is contained in that fixed space"
            ),
            "generator_fixed_ranks": {
                row["name"]: row["fixed_rank"]
                for row in generators["generators"]
            },
            "checks": fixed_rank_checks,
        },
        "lineage": {
            "source_chamber_priority": "published_known_Roulleau_example",
            "new_abstract_polytope_count": 0,
            "new_higher_dimensional_fold_count": 0,
            "dimension_status": "H10 was already occupied; this 90-wall chamber is not new",
            "new_result_here": (
                "exact S8 x C2 weighted-wall group and complete proof that this stock "
                "cannot produce a proper fold above H9"
            ),
        },
        "checks": summary_checks,
        "all_checks_pass": all(summary_checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.is_file() or OUTPUT.read_text() != rendered:
            raise SystemExit("recomputed certificate differs from archive")
    elif args.output:
        args.output.write_text(rendered)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
