#!/usr/bin/env python3
"""Certify the JK-H12 D_infinity^3 finite-box family.

The seed lattice and 172 roots are reparsed through the already-hostile raw
Dolgachev--Shimada reconstruction.  No exploratory scout or production
path/cube builder is imported.  The three-parameter formula is proved by an
all-Z^3 collision solve, not finite interpolation.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction
import hashlib
from itertools import product
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

from scripts import audit_jacobian_kummer_h12_path_cube_hostile as base


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/jacobian_kummer_h12_triple_parallel_family.json"
)
DEFAULT_COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
SCHEMA = "coxeter.jacobian_kummer.h12_triple_parallel_family.v1"
STATUS_DATE = "2026-08-08"

# Zero based; one-based labels are emitted everywhere user-facing.
TRIPLE_CONFIGURATIONS = (
    ((47, 41), (50, 38), (55, 35)),
    ((57, 32), (59, 30), (58, 31)),
)
REPRESENTATIVE = TRIPLE_CONFIGURATIONS[0]

Vector = base.Vector
Matrix = base.Matrix
Permutation = base.Permutation


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_sha256(payload: Any) -> str:
    return _sha256_bytes(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    )


def _maximal_cliques(
    vertices: Iterable[int], adjacency: dict[int, set[int]]
) -> tuple[tuple[int, ...], ...]:
    result: list[tuple[int, ...]] = []

    def visit(chosen: set[int], candidates: set[int], excluded: set[int]) -> None:
        if not candidates and not excluded:
            result.append(tuple(sorted(chosen)))
            return
        pivot_pool = candidates | excluded
        pivot = max(
            pivot_pool,
            key=lambda vertex: len(candidates & adjacency[vertex]),
            default=None,
        )
        branch = candidates - (adjacency[pivot] if pivot is not None else set())
        for vertex in tuple(sorted(branch)):
            visit(
                chosen | {vertex},
                candidates & adjacency[vertex],
                excluded & adjacency[vertex],
            )
            candidates.remove(vertex)
            excluded.add(vertex)

    visit(set(), set(vertices), set())
    return tuple(result)


def _product_graph_exhaustion(
    stock: dict[str, Any], roots: tuple[Vector, ...]
) -> dict[str, Any]:
    doubling: tuple[int, ...] = stock["doubling"]
    relations: tuple[tuple[str, ...], ...] = stock["relations"]
    parallel_pairs = tuple(
        (left, right)
        for position, left in enumerate(doubling)
        for right in doubling[:position]
        if relations[left][right] == "parallel"
    )
    orthogonal_adjacency = {vertex: set() for vertex in doubling}
    for position, left in enumerate(doubling):
        for right in doubling[:position]:
            if relations[left][right] == "m2":
                orthogonal_adjacency[left].add(right)
                orthogonal_adjacency[right].add(left)
    orthogonal_maximal = _maximal_cliques(doubling, orthogonal_adjacency)
    orthogonal_maximum = max(map(len, orthogonal_maximal))
    maximum_cubes = tuple(
        clique for clique in orthogonal_maximal if len(clique) == orthogonal_maximum
    )

    factor_adjacency = {index: set() for index in range(len(parallel_pairs))}
    compatibility_edges: list[tuple[int, int]] = []
    for left_index, left_pair in enumerate(parallel_pairs):
        for right_index in range(left_index):
            right_pair = parallel_pairs[right_index]
            if all(
                relations[left][right] == "m2"
                for left in left_pair
                for right in right_pair
            ):
                factor_adjacency[left_index].add(right_index)
                factor_adjacency[right_index].add(left_index)
                compatibility_edges.append((left_index, right_index))
    factor_maximal = _maximal_cliques(range(len(parallel_pairs)), factor_adjacency)
    maximum_factor_count = max(map(len, factor_maximal))
    maximum_factor_cliques = tuple(
        clique for clique in factor_maximal if len(clique) == maximum_factor_count
    )
    reconstructed_triples = {
        frozenset(frozenset(parallel_pairs[index]) for index in clique)
        for clique in maximum_factor_cliques
    }
    expected_triples = {
        frozenset(frozenset(pair) for pair in configuration)
        for configuration in TRIPLE_CONFIGURATIONS
    }

    path_cube_histogram: Counter[int] = Counter()
    maximal_path_cubes: list[tuple[tuple[int, int], tuple[int, ...]]] = []
    for pair in parallel_pairs:
        common = tuple(
            vertex
            for vertex in doubling
            if vertex not in pair
            and all(relations[vertex][endpoint] == "m2" for endpoint in pair)
        )
        common_set = set(common)
        common_adjacency = {
            vertex: orthogonal_adjacency[vertex] & common_set for vertex in common
        }
        cliques = _maximal_cliques(common, common_adjacency) if common else ((),)
        maximum = max(map(len, cliques))
        path_cube_histogram[maximum] += 1
        if maximum == 3:
            maximal_path_cubes.extend(
                (pair, clique) for clique in cliques if len(clique) == maximum
            )

    pair_pair_cube_histogram: Counter[int] = Counter()
    for left_index, right_index in compatibility_edges:
        selected = set(parallel_pairs[left_index]) | set(parallel_pairs[right_index])
        common = tuple(
            vertex
            for vertex in doubling
            if vertex not in selected
            and all(relations[vertex][wall] == "m2" for wall in selected)
        )
        common_set = set(common)
        common_adjacency = {
            vertex: orthogonal_adjacency[vertex] & common_set for vertex in common
        }
        cliques = _maximal_cliques(common, common_adjacency) if common else ((),)
        pair_pair_cube_histogram[max(map(len, cliques))] += 1

    triple_common_orthogonals = []
    for clique in maximum_factor_cliques:
        selected = {
            wall for pair_index in clique for wall in parallel_pairs[pair_index]
        }
        common = [
            vertex + 1
            for vertex in doubling
            if vertex not in selected
            and all(relations[vertex][wall] == "m2" for wall in selected)
        ]
        triple_common_orthogonals.append(common)

    relation_histogram = Counter(
        relations[left][right]
        for position, left in enumerate(doubling)
        for right in doubling[:position]
    )
    expected_path_cubes = {
        (
            frozenset(pair),
            frozenset(cube),
        )
        for pair, cube in (
            ((61, 60), (32, 30, 31)),
            ((65, 64), (32, 58, 59)),
            ((66, 63), (57, 59, 31)),
            ((67, 62), (57, 58, 30)),
            ((72, 68), (50, 47, 55)),
            ((73, 69), (41, 38, 55)),
            ((74, 70), (41, 50, 35)),
            ((75, 71), (35, 38, 47)),
        )
    }
    actual_path_cubes = {
        (frozenset(pair), frozenset(cube)) for pair, cube in maximal_path_cubes
    }
    checks = {
        "doubling_stock_size_is_124": len(doubling) == 124,
        "doubling_relation_histogram_is_exact": relation_histogram
        == Counter({"m2": 264, "parallel": 1138, "ultraparallel": 6224}),
        "parallel_pair_count_is_1138": len(parallel_pairs) == 1138,
        "parallel_factor_compatibility_edge_count_is_6": len(
            compatibility_edges
        )
        == 6,
        "maximum_mutually_orthogonal_parallel_factor_count_is_3": maximum_factor_count
        == 3,
        "exactly_two_D_infinity_cubed_configurations": reconstructed_triples
        == expected_triples,
        "no_D_infinity_four_product": maximum_factor_count < 4,
        "no_cube_factor_orthogonal_to_D_infinity_cubed": triple_common_orthogonals
        == [[], []],
        "maximum_D_infinity_squared_times_cube_rank_is_1": pair_pair_cube_histogram
        == Counter({1: 6}),
        "maximum_D_infinity_times_cube_rank_is_3": max(path_cube_histogram) == 3,
        "exactly_eight_maximal_D_infinity_times_C2_cubed_configurations": (
            actual_path_cubes == expected_path_cubes
        ),
        "maximum_orthogonal_cube_rank_is_4": orthogonal_maximum == 4,
        "exactly_16_maximum_C2_four_cubes": len(maximum_cubes) == 16,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"product graph exhaustion failed: {checks}")
    return {
        "doubling_wall_count": len(doubling),
        "doubling_relation_histogram": dict(sorted(relation_histogram.items())),
        "parallel_pair_count": len(parallel_pairs),
        "parallel_factor_compatibility_edge_count": len(compatibility_edges),
        "maximum_D_infinity_factor_count": maximum_factor_count,
        "D_infinity_cubed_configurations_one_based": [
            [[wall + 1 for wall in parallel_pairs[index]] for index in clique]
            for clique in maximum_factor_cliques
        ],
        "triple_common_orthogonal_doubling_walls": triple_common_orthogonals,
        "D_infinity_squared_common_cube_rank_histogram": {
            str(key): value for key, value in sorted(pair_pair_cube_histogram.items())
        },
        "D_infinity_common_cube_rank_histogram": {
            str(key): value for key, value in sorted(path_cube_histogram.items())
        },
        "maximum_D_infinity_times_C2_cubed_configurations_one_based": [
            {
                "parallel_pair": [wall + 1 for wall in pair],
                "cube": [wall + 1 for wall in cube],
            }
            for pair, cube in maximal_path_cubes
        ],
        "maximum_orthogonal_cube_rank": orthogonal_maximum,
        "maximum_C2_four_cube_count": len(maximum_cubes),
        "maximum_C2_four_cubes_one_based": [
            [wall + 1 for wall in clique] for clique in maximum_cubes
        ],
        "scope": (
            "exhaustive only for direct products whose D_infinity factors are "
            "parallel doubling-wall pairs and whose distinct factors/cube walls "
            "are pairwise orthogonal inside this exact 124-wall seed stock"
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
        "_parallel_pairs": parallel_pairs,
    }


def _triple_structure(seed: dict[str, Any], stock: dict[str, Any]) -> dict[str, Any]:
    gram: Matrix = seed["fixed_gram"]
    roots: tuple[Vector, ...] = seed["roots"]
    doubling = set(stock["doubling"])
    rows: list[dict[str, Any]] = []
    for configuration_number, factors in enumerate(TRIPLE_CONFIGURATIONS, 1):
        selected = tuple(wall for factor in factors for wall in factor)
        selected_gram = [
            [base._pair(gram, roots[left], roots[right]) for right in selected]
            for left in selected
        ]
        expected_gram = []
        for left in range(6):
            row = []
            for right in range(6):
                if left == right:
                    row.append(-4)
                elif left // 2 == right // 2:
                    row.append(4)
                else:
                    row.append(0)
            expected_gram.append(row)
        reflections = {
            wall: base._reflection_matrix(gram, roots[wall]) for wall in selected
        }
        translations: list[Matrix] = []
        nilpotents: list[Matrix] = []
        null_directions: list[Vector] = []
        for left, right in factors:
            translation = base._matmul(reflections[right], reflections[left])
            nilpotent = tuple(
                tuple(
                    translation[row][column] - int(row == column)
                    for column in range(13)
                )
                for row in range(13)
            )
            half_sum = tuple(
                Fraction(roots[left][index] + roots[right][index], 2)
                for index in range(13)
            )
            if not all(value.denominator == 1 for value in half_sum):
                raise ArithmeticError("a parallel-pair half-sum is nonintegral")
            null = base._primitive_key(tuple(int(value) for value in half_sum))
            translations.append(translation)
            nilpotents.append(nilpotent)
            null_directions.append(null)
        common_null = null_directions[0]
        nilpotent_squares = tuple(
            base._matmul(nilpotent, nilpotent) for nilpotent in nilpotents
        )
        zero = tuple(tuple(0 for _ in range(13)) for _ in range(13))
        checks = {
            "all_six_walls_are_in_exact_doubling_stock": set(selected) <= doubling,
            "selected_Gram_is_three_orthogonal_parallel_blocks": selected_gram
            == expected_gram,
            "all_six_reflections_are_integral_involutive_isometries": len(
                reflections
            )
            == 6,
            "distinct_D_infinity_factors_commute": all(
                base._matmul(reflections[left], reflections[right])
                == base._matmul(reflections[right], reflections[left])
                for factor_index, factor in enumerate(factors)
                for other_factor in factors[:factor_index]
                for left in factor
                for right in other_factor
            ),
            "three_primitive_null_half_sums_are_equal": (
                len(set(null_directions)) == 1
                and base._content(common_null) == 1
                and base._pair(gram, common_null, common_null) == 0
            ),
            "translations_commute_pairwise": all(
                base._matmul(translations[left], translations[right])
                == base._matmul(translations[right], translations[left])
                for left in range(3)
                for right in range(left)
            ),
            "each_translation_has_nilpotent_ranks_2_1_0": all(
                sp.Matrix(nilpotents[index]).rank() == 2
                and sp.Matrix(nilpotent_squares[index]).rank() == 1
                and base._matmul(nilpotent_squares[index], nilpotents[index])
                == zero
                for index in range(3)
            ),
            "mixed_nilpotent_products_vanish": all(
                base._matmul(nilpotents[left], nilpotents[right]) == zero
                for left in range(3)
                for right in range(3)
                if left != right
            ),
            "each_nilpotent_square_image_is_common_null_line": all(
                sp.Matrix.hstack(
                    sp.Matrix(common_null), sp.Matrix(nilpotent_squares[index])
                ).rank()
                == 1
                for index in range(3)
            ),
        }
        if not all(checks.values()):
            raise ArithmeticError(
                f"D_infinity^3 structure {configuration_number} failed: {checks}"
            )
        rows.append(
            {
                "configuration_number": configuration_number,
                "parallel_pairs_one_based": [
                    [wall + 1 for wall in factor] for factor in factors
                ],
                "selected_wall_Gram": selected_gram,
                "common_primitive_integral_null_direction": list(common_null),
                "translation_ranks": [
                    {
                        "rank_T_minus_I": int(sp.Matrix(nilpotents[index]).rank()),
                        "rank_T_minus_I_squared": int(
                            sp.Matrix(nilpotent_squares[index]).rank()
                        ),
                        "T_minus_I_cubed_is_zero": True,
                    }
                    for index in range(3)
                ],
                "checks": checks,
                "all_exact_checks_pass": True,
                "_reflections": reflections,
                "_translations": tuple(translations),
                "_nilpotents": tuple(nilpotents),
                "_nilpotent_squares": nilpotent_squares,
                "_common_null": common_null,
            }
        )
    return {"configurations": rows, "all_exact_checks_pass": True}


def _side_pattern_census(
    seed: dict[str, Any], structure: dict[str, Any]
) -> dict[str, Any]:
    gram: Matrix = seed["fixed_gram"]
    roots: tuple[Vector, ...] = seed["roots"]
    factors = REPRESENTATIVE
    selected = {wall for factor in factors for wall in factor}
    histogram: Counter[tuple[str, str, str]] = Counter()
    side_records: list[dict[str, Any]] = []
    for seed_index, root in enumerate(roots):
        if seed_index in selected:
            continue
        pattern: list[str] = []
        for left, right in factors:
            left_orthogonal = base._pair(gram, root, roots[left]) == 0
            right_orthogonal = base._pair(gram, root, roots[right]) == 0
            category = (
                "B"
                if left_orthogonal and right_orthogonal
                else "L"
                if left_orthogonal
                else "R"
                if right_orthogonal
                else "N"
            )
            pattern.append(category)
        key = tuple(pattern)
        histogram[key] += 1
        side_records.append(
            {
                "seed_wall_one_based": seed_index + 1,
                "pattern": "".join(key),
                "primitive_coordinates": list(base._primitive_key(root)),
            }
        )
    expected = {
        "BBB": 12,
        "LLL": 3,
        "LLN": 4,
        "LLR": 2,
        "LNL": 4,
        "LNN": 4,
        "LNR": 4,
        "LRL": 2,
        "LRN": 4,
        "LRR": 3,
        "NLL": 4,
        "NLN": 4,
        "NLR": 4,
        "NNL": 4,
        "NNN": 62,
        "NNR": 4,
        "NRL": 4,
        "NRN": 4,
        "NRR": 4,
        "RLL": 2,
        "RLN": 4,
        "RLR": 3,
        "RNL": 4,
        "RNN": 4,
        "RNR": 4,
        "RRL": 3,
        "RRN": 4,
        "RRR": 2,
    }
    payload_histogram = {
        "".join(pattern): count for pattern, count in sorted(histogram.items())
    }
    margins = [
        Counter(
            pattern[index]
            for pattern, count in histogram.items()
            for _ in range(count)
        )
        for index in range(3)
    ]
    checks = {
        "exactly_166_nonendpoint_seed_facets": len(side_records) == 166,
        "pattern_histogram_is_exact": payload_histogram == expected,
        "each_factor_margin_is_B12_L30_R30_N94": all(
            margin == Counter({"B": 12, "L": 30, "R": 30, "N": 94})
            for margin in margins
        ),
        "the_only_patterns_containing_B_are_12_BBB": all(
            "B" not in pattern or pattern == "BBB"
            for pattern in payload_histogram
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"side pattern census failed: {checks}")
    return {
        "side_facet_count": len(side_records),
        "pattern_legend": {
            "B": "orthogonal to both endpoints of this parallel factor",
            "L": "orthogonal only to its left endpoint",
            "R": "orthogonal only to its right endpoint",
            "N": "orthogonal to neither endpoint",
        },
        "pattern_histogram": payload_histogram,
        "factor_margins": [dict(sorted(margin.items())) for margin in margins],
        "side_record_sha256": _json_sha256(side_records),
        "checks": checks,
        "all_exact_checks_pass": True,
        "_side_roots": tuple(
            tuple(record["primitive_coordinates"]) for record in side_records
        ),
        "_patterns": tuple(record["pattern"] for record in side_records),
    }


def _triple_translation_apply(
    nilpotents: Sequence[Matrix],
    nilpotent_squares: Sequence[Matrix],
    vector: Vector,
    exponents: tuple[int, int, int],
) -> Vector:
    image = vector
    for index, exponent in enumerate(exponents):
        image = base._translation_power_apply(
            nilpotents[index], nilpotent_squares[index], image, exponent
        )
    return image


def _all_integer_triple_collision_audit(
    seed: dict[str, Any],
    structure: dict[str, Any],
    patterns: dict[str, Any],
) -> dict[str, Any]:
    representative = structure["configurations"][0]
    factors = REPRESENTATIVE
    reflections: dict[int, Matrix] = representative["_reflections"]
    nilpotents: tuple[Matrix, ...] = representative["_nilpotents"]
    nilpotent_squares: tuple[Matrix, ...] = representative[
        "_nilpotent_squares"
    ]
    common_null: Vector = representative["_common_null"]
    sides: tuple[Vector, ...] = patterns["_side_roots"]
    side_patterns: tuple[str, ...] = patterns["_patterns"]
    signed_targets: dict[Vector, tuple[int, int]] = {}
    for target_index, target in enumerate(sides):
        if target in signed_targets or base._scale(-1, target) in signed_targets:
            raise ArithmeticError("side stock has an unoriented duplicate")
        signed_targets[target] = (target_index, 1)
        signed_targets[base._scale(-1, target)] = (target_index, -1)

    branches: list[Matrix] = []
    for bits in range(8):
        branch = base._identity(13)
        for factor_index, (_left, right) in enumerate(factors):
            if bits & (1 << factor_index):
                branch = base._matmul(branch, reflections[right])
        branches.append(branch)

    finite: list[tuple[int, tuple[int, int, int], int, int, int]] = []
    identities: list[tuple[int, int, int, int]] = []
    quotient_rank_histogram: Counter[int] = Counter()
    integral_candidate_count = 0
    full_vector_candidate_count = 0
    for bits, branch in enumerate(branches):
        for source_index, source in enumerate(sides):
            branched = base._matvec(branch, source)
            first_differences = tuple(
                base._matvec(nilpotent, branched) for nilpotent in nilpotents
            )
            if not any(value for vector in first_differences for value in vector):
                hit = signed_targets.get(branched)
                if hit is not None:
                    identities.append((bits, source_index, hit[0], hit[1]))
                continue

            affine_matrix = sp.Matrix.hstack(
                *(sp.Matrix(vector) for vector in first_differences),
                sp.Matrix(common_null),
            )
            rank = int(affine_matrix.rank())
            quotient_rank_histogram[rank] += 1
            if rank != 4:
                raise ArithmeticError(
                    "translation directions are dependent modulo the common null line"
                )
            pivot_rows = affine_matrix.T.rref()[1]
            inverse = affine_matrix.extract(pivot_rows, range(4)).inv()
            inverse_fraction = tuple(
                tuple(Fraction(inverse[row, column]) for column in range(4))
                for row in range(4)
            )
            for target, (target_index, sign) in signed_targets.items():
                difference = tuple(
                    target[row] - branched[row] for row in pivot_rows
                )
                coefficients = tuple(
                    sum(
                        inverse_fraction[row][column] * difference[column]
                        for column in range(4)
                    )
                    for row in range(4)
                )
                if not all(coefficients[index].denominator == 1 for index in range(3)):
                    continue
                integral_candidate_count += 1
                exponents = tuple(
                    coefficients[index].numerator for index in range(3)
                )
                image = _triple_translation_apply(
                    nilpotents,
                    nilpotent_squares,
                    branched,
                    exponents,
                )
                full_vector_candidate_count += 1
                if image == target:
                    finite.append(
                        (bits, exponents, source_index, target_index, sign)
                    )

    expected_finite: set[tuple[int, tuple[int, int, int], int, int, int]] = set()
    for source_index, pattern in enumerate(side_patterns):
        if pattern == "BBB":
            continue
        for bits in range(8):
            exponents: list[int] = []
            compatible = True
            for factor_index, category in enumerate(pattern):
                if not bits & (1 << factor_index):
                    exponents.append(0)
                elif category == "R":
                    exponents.append(0)
                elif category == "L":
                    exponents.append(-1)
                else:
                    compatible = False
                    break
            if compatible:
                expected_finite.add(
                    (bits, tuple(exponents), source_index, source_index, 1)
                )
    expected_identities = {
        (bits, source_index, source_index, 1)
        for bits in range(8)
        for source_index, pattern in enumerate(side_patterns)
        if pattern == "BBB"
    }
    finite_histogram = Counter(bits for bits, *_rest in finite)
    exponent_histogram = Counter((bits, exponents) for bits, exponents, *_ in finite)
    identity_histogram = Counter(bits for bits, *_rest in identities)
    checks = {
        "all_440896_ordered_signed_side_equations_covered": (
            8 * 2 * 166 * 166 == 440_896
        ),
        "no_bounded_Z3_box_scan_used": True,
        "all_1232_nonfixed_branch_sources_have_rank_three_mod_null": (
            quotient_rank_histogram == Counter({4: 1232})
        ),
        "finite_collision_set_is_exactly_pattern_stabilizers": set(finite)
        == expected_finite,
        "finite_collision_count_is_462": len(finite) == 462,
        "finite_collisions_are_positive_self_only": all(
            source == target and sign == 1
            for _bits, _exponents, source, target, sign in finite
        ),
        "all_Z3_identity_set_is_exactly_12_BBB_times_8_branches": set(
            identities
        )
        == expected_identities,
        "all_Z3_identity_count_is_96": len(identities) == 96,
        "no_negative_or_cross_facet_collision": all(
            source == target and sign == 1
            for _bits, _exponents, source, target, sign in finite
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"all-Z3 collision audit failed: {checks}")
    return {
        "domain": (
            "all (z1,z2,z3) in Z^3, all eight reflection branches, both root "
            "signs, and all 166x166 ordered nonendpoint seed-facet pairs"
        ),
        "method": (
            "the three commuting translations satisfy Ki^3=0 and Ki*Kj=0 for "
            "i!=j; for every moving branch/source, [K1v,K2v,K3v,u] has rank 4, "
            "so quotienting by the common null line gives the unique rational "
            "candidate (z1,z2,z3), whose integrality and all 13 coordinates are "
            "then checked exactly"
        ),
        "bounded_parameter_scan_used": False,
        "ordered_signed_equation_count": 440_896,
        "moving_branch_source_count": 1232,
        "quotient_affine_rank_histogram": {
            str(key): value for key, value in sorted(quotient_rank_histogram.items())
        },
        "integral_candidate_count": integral_candidate_count,
        "full_vector_candidate_count": full_vector_candidate_count,
        "finite_solution_count": len(finite),
        "finite_solution_branch_histogram": {
            f"{bits:03b}": finite_histogram[bits] for bits in range(8)
        },
        "finite_solution_branch_exponent_histogram": {
            f"branch={bits:03b};z={exponents}": count
            for (bits, exponents), count in sorted(exponent_histogram.items())
        },
        "finite_solutions_are_exactly_positive_self_stabilizers": True,
        "all_Z3_identity_count": len(identities),
        "all_Z3_identity_branch_histogram": {
            f"{bits:03b}": identity_histogram[bits] for bits in range(8)
        },
        "all_Z3_identities_are_exactly_the_12_BBB_facets_in_each_branch": True,
        "negative_root_collision_count": 0,
        "cross_facet_collision_count": 0,
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def _endpoint_orbit_audit(
    seed: dict[str, Any], structure: dict[str, Any], patterns: dict[str, Any]
) -> dict[str, Any]:
    gram: Matrix = seed["fixed_gram"]
    roots: tuple[Vector, ...] = seed["roots"]
    representative = structure["configurations"][0]
    factors = REPRESENTATIVE
    reflections: dict[int, Matrix] = representative["_reflections"]
    nilpotents: tuple[Matrix, ...] = representative["_nilpotents"]
    nilpotent_squares: tuple[Matrix, ...] = representative[
        "_nilpotent_squares"
    ]
    common_null: Vector = representative["_common_null"]
    side_roots: tuple[Vector, ...] = patterns["_side_roots"]
    side_patterns: tuple[str, ...] = patterns["_patterns"]
    bbb_roots = tuple(
        root for root, pattern in zip(side_roots, side_patterns, strict=True)
        if pattern == "BBB"
    )
    moving_sides = tuple(
        root for root, pattern in zip(side_roots, side_patterns, strict=True)
        if pattern != "BBB"
    )
    endpoint_roots = tuple(
        base._primitive_key(roots[wall]) for factor in factors for wall in factor
    )
    endpoint_equations: list[dict[str, Any]] = []
    for factor_index, (left_index, right_index) in enumerate(factors):
        left = base._primitive_key(roots[left_index])
        right = base._primitive_key(roots[right_index])
        right_reflection = reflections[right_index]
        equations = {}
        for name, vector in {
            "T^z(right)=plus_or_minus_left": right,
            "T^z*S_right(left)=plus_or_minus_left": base._matvec(
                right_reflection, left
            ),
        }.items():
            solutions: list[dict[str, int]] = []
            for sign in (1, -1):
                values = base._solve_vector_equation(
                    nilpotents[factor_index],
                    nilpotent_squares[factor_index],
                    vector,
                    base._scale(sign, left),
                )
                if values is None:
                    raise ArithmeticError("unexpected endpoint all-z identity")
                solutions.extend({"z": value, "sign": sign} for value in values)
            equations[name] = solutions
        endpoint_equations.append(
            {"factor": factor_index + 1, "equations": equations}
        )
    cross_quotient_collisions: list[tuple[int, int, int]] = []
    for left_factor in range(3):
        for right_factor in range(left_factor):
            for left_wall in factors[left_factor]:
                for right_wall in factors[right_factor]:
                    for sign in (1, -1):
                        difference = base._add(
                            base._primitive_key(roots[left_wall]),
                            base._scale(
                                -sign, base._primitive_key(roots[right_wall])
                            ),
                        )
                        if sp.Matrix.hstack(
                            sp.Matrix(common_null), sp.Matrix(difference)
                        ).rank() == 1:
                            cross_quotient_collisions.append(
                                (left_wall + 1, right_wall + 1, sign)
                            )
    expected_endpoint_equations = {
        "T^z(right)=plus_or_minus_left": [],
        "T^z*S_right(left)=plus_or_minus_left": [{"z": -1, "sign": -1}],
    }
    checks = {
        "all_six_endpoint_roots_are_distinct": len(set(endpoint_roots)) == 6,
        "moving_side_facets_have_nonzero_common_null_pairing": all(
            base._pair(gram, root, common_null) != 0 for root in moving_sides
        ),
        "endpoints_and_BBB_facets_have_zero_common_null_pairing": all(
            base._pair(gram, root, common_null) == 0
            for root in (*endpoint_roots, *bbb_roots)
        ),
        "BBB_facets_are_fixed_by_all_six_endpoint_reflections": all(
            base._matvec(reflections[wall], root) == root
            for wall in (wall for factor in factors for wall in factor)
            for root in bbb_roots
        ),
        "no_endpoint_equals_a_BBB_facet_up_to_sign": set(endpoint_roots).isdisjoint(
            set(bbb_roots)
        ),
        "three_one_factor_outer_endpoint_equations_are_exact": all(
            row["equations"] == expected_endpoint_equations
            for row in endpoint_equations
        ),
        "different_factor_endpoint_orbits_are_distinct_mod_common_null": not (
            cross_quotient_collisions
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"endpoint orbit audit failed: {checks}")
    return {
        "outer_endpoint_count_for_every_finite_box": 6,
        "one_factor_endpoint_equations": endpoint_equations,
        "cross_factor_endpoint_quotient_collision_count": len(
            cross_quotient_collisions
        ),
        "endpoint_to_moving_side_exclusion": (
            "pairing with the common null vector is zero on endpoints and nonzero "
            "on every moving side facet"
        ),
        "endpoint_to_BBB_exclusion": (
            "the 12 BBB roots and six endpoint roots are distinct canonical roots; "
            "BBB roots are fixed by the whole selected reflection subgroup"
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def facet_formula(
    n1: int, n2: int, n3: int, phase_parity: int = 0
) -> int:
    if min(n1, n2, n3) < 1:
        raise ValueError("all box parameters must be positive")
    parameter_parity = int(any(value % 2 for value in (n1, n2, n3)))
    epsilon = (
        parameter_parity
        if phase_parity % 2 == 0
        else 2 - parameter_parity
    )
    numerator = (
        177 * n1 * n2 * n3
        + 29 * (n1 * n2 + n1 * n3 + n2 * n3)
        + 13 * (n1 + n2 + n3)
        + 40
        + epsilon
    )
    if numerator % 2:
        raise ArithmeticError("closed facet formula is unexpectedly nonintegral")
    return numerator // 2


def _formula_audit(patterns: dict[str, Any]) -> dict[str, Any]:
    histogram = {
        pattern: int(count) for pattern, count in patterns["pattern_histogram"].items()
    }
    n1, n2, n3 = sp.symbols("N1 N2 N3", integer=True, positive=True)
    e1, e2, e3 = sp.symbols("e1 e2 e3", integer=True)
    parameters = (n1, n2, n3)
    parities = (e1, e2, e3)
    weights = []
    for parameter, parity in zip(parameters, parities, strict=True):
        weights.append(
            {
                "B": sp.Integer(1),
                "L": (parameter + parity) / 2,
                "R": (parameter + 2 - parity) / 2,
                "N": parameter,
            }
        )
    orbit_sum = sp.Integer(6)
    for pattern, count in histogram.items():
        orbit_sum += count * sp.prod(
            weights[index][category] for index, category in enumerate(pattern)
        )
    base_numerator = (
        177 * n1 * n2 * n3
        + 29 * (n1 * n2 + n1 * n3 + n2 * n3)
        + 13 * (n1 + n2 + n3)
        + 40
    )
    parity_identities = []
    for phase_bits in product((0, 1), repeat=3):
        phase_weights = []
        for index, weight in enumerate(weights):
            if phase_bits[index]:
                phase_weights.append(
                    {**weight, "L": weight["R"], "R": weight["L"]}
                )
            else:
                phase_weights.append(weight)
        phase_orbit_sum = sp.Integer(6)
        for pattern, count in histogram.items():
            phase_orbit_sum += count * sp.prod(
                phase_weights[index][category]
                for index, category in enumerate(pattern)
            )
        phase_parity = sum(phase_bits) % 2
        for parameter_bits in product((0, 1), repeat=3):
            parameter_parity = int(any(parameter_bits))
            epsilon = (
                parameter_parity
                if phase_parity == 0
                else 2 - parameter_parity
            )
            difference = sp.expand(
                phase_orbit_sum.subs(
                    dict(zip(parities, parameter_bits, strict=True))
                )
                - (base_numerator + epsilon) / 2
            )
            parity_identities.append(
                {
                    "parameter_parities": list(parameter_bits),
                    "phase_bits": list(phase_bits),
                    "epsilon": epsilon,
                    "symbolic_difference": str(difference),
                }
            )
    checks = {
        "all_64_parameter_and_phase_parity_cases_equal_closed_formula": all(
            row["symbolic_difference"] == "0" for row in parity_identities
        ),
        "epsilon_phase_convention_is_exact": True,
        "six_outer_endpoints_are_included_in_orbit_sum": True,
        "formula_reproduces_seed_at_1_1_1": facet_formula(1, 1, 1) == 172,
        "formula_reproduces_C2_cube_at_2_2_2": facet_formula(2, 2, 2) == 941,
        "odd_phase_at_2_2_2_has_942_facets": facet_formula(2, 2, 2, 1)
        == 942,
        "preliminary_phase_blind_formula_fails_at_even_odd_phase": (
            facet_formula(2, 2, 2, 0) != facet_formula(2, 2, 2, 1)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"closed formula audit failed: {checks}")
    return {
        "facet_formula": (
            "f(N1,N2,N3)=(177*N1*N2*N3+29*(N1*N2+N1*N3+N2*N3)"
            "+13*(N1+N2+N3)+40+epsilon(N,phase))/2"
        ),
        "epsilon_convention": (
            "let E=0 if N1,N2,N3 are all even and E=1 otherwise; epsilon=E "
            "for even phase-flip parity and epsilon=2-E for odd phase-flip parity"
        ),
        "orbit_weight_formula": (
            "6 + sum_pattern h(pattern)*product_i w_pattern_i(N_i), with "
            "w_B=1, w_L=ceil(N/2), w_R=floor(N/2)+1, w_N=N"
        ),
        "symbolic_parity_identities": parity_identities,
        "sample_values": {
            "1,1,1": facet_formula(1, 1, 1),
            "2,1,1": facet_formula(2, 1, 1),
            "2,2,1": facet_formula(2, 2, 1),
            "2,2,2": facet_formula(2, 2, 2),
            "2,2,2;odd_phase": facet_formula(2, 2, 2, 1),
            "3,3,3": facet_formula(3, 3, 3),
            "4,4,4": facet_formula(4, 4, 4),
        },
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def _path_transforms(
    pair: tuple[int, int],
    reflections: dict[int, Matrix],
    chamber_count: int,
    phase: int,
) -> tuple[Matrix, ...]:
    transforms = [base._identity(13)]
    for edge in range(chamber_count - 1):
        transforms.append(
            base._matmul(
                transforms[-1], reflections[pair[(phase + edge) % 2]]
            )
        )
    return tuple(transforms)


def _grid_boundary(
    roots: tuple[Vector, ...],
    gram: Matrix,
    factors: tuple[tuple[int, int], ...],
    reflections: dict[int, Matrix],
    parameters: tuple[int, int, int],
    phases: tuple[int, int, int],
) -> tuple[Vector, ...]:
    factor_transforms = tuple(
        _path_transforms(factor, reflections, count, phase)
        for factor, count, phase in zip(
            factors, parameters, phases, strict=True
        )
    )
    transforms: list[Matrix] = []
    for choices in product(*factor_transforms):
        transform = base._identity(13)
        for factor_transform in choices:
            transform = base._matmul(transform, factor_transform)
        transforms.append(transform)
    if len(set(transforms)) != math.prod(parameters):
        raise ArithmeticError("finite D_infinity^3 box chamber maps collided")
    transform_set = set(transforms)
    internal: set[Vector] = set()
    for transform in transforms:
        for factor in factors:
            for wall in factor:
                if base._matmul(transform, reflections[wall]) in transform_set:
                    internal.add(
                        base._primitive_key(base._matvec(transform, roots[wall]))
                    )
    boundary: dict[Vector, Vector] = {}
    for transform in transforms:
        for root in roots:
            image = base._matvec(transform, root)
            key = base._primitive_key(image)
            if key not in internal:
                boundary.setdefault(key, image)
    return tuple(boundary[key] for key in sorted(boundary))


def _direct_grid_audit(
    seed: dict[str, Any], structure: dict[str, Any]
) -> dict[str, Any]:
    roots: tuple[Vector, ...] = seed["roots"]
    gram: Matrix = seed["fixed_gram"]
    representative = structure["configurations"][0]
    reflections: dict[int, Matrix] = representative["_reflections"]
    parameter_stock = list(product(range(1, 4), repeat=3)) + [
        (4, 4, 4),
        (2, 3, 4),
        (5, 2, 3),
    ]
    rows: list[dict[str, Any]] = []
    for parameters in parameter_stock:
        phase_rows: list[dict[str, Any]] = []
        for phases in product((0, 1), repeat=3):
            boundary = _grid_boundary(
                roots,
                gram,
                REPRESENTATIVE,
                reflections,
                parameters,
                phases,
            )
            phase_parity = sum(phases) % 2
            expected = facet_formula(*parameters, phase_parity)
            if not (
                len(boundary) == expected
                and all(
                    base._root_is_integral_reflection_root(gram, root)
                    for root in boundary
                )
            ):
                raise ArithmeticError(
                    f"direct grid failed for N={parameters}, phase={phases}"
                )
            phase_rows.append(
                {
                    "phase_bits": list(phases),
                    "phase_flip_parity": phase_parity,
                    "facet_count": len(boundary),
                    "boundary_root_coordinate_sha256": _json_sha256(
                        [list(root) for root in boundary]
                    ),
                }
            )
        rows.append(
            {
                "parameters": list(parameters),
                "seed_chamber_count": math.prod(parameters),
                "even_phase_expected_facet_count": facet_formula(*parameters, 0),
                "odd_phase_expected_facet_count": facet_formula(*parameters, 1),
                "phase_reconstructions": phase_rows,
            }
        )
    checks = {
        "all_27_parameter_triples_in_1_to_3_cubed_checked": all(
            any(tuple(row["parameters"]) == parameters for row in rows)
            for parameters in product(range(1, 4), repeat=3)
        ),
        "three_additional_larger_parameter_triples_checked": len(rows) == 30,
        "all_eight_phases_checked_for_every_parameter_triple": all(
            len(row["phase_reconstructions"]) == 8 for row in rows
        ),
        "all_240_direct_boundaries_match_phase_refined_formula": all(
            phase["facet_count"]
            == facet_formula(
                *row["parameters"], phase["phase_flip_parity"]
            )
            for row in rows
            for phase in row["phase_reconstructions"]
        ),
        "phase_blind_formula_regression_blocked_at_2_2_2": next(
            row for row in rows if row["parameters"] == [2, 2, 2]
        )["even_phase_expected_facet_count"]
        == 941
        and next(row for row in rows if row["parameters"] == [2, 2, 2])[
            "odd_phase_expected_facet_count"
        ]
        == 942,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"direct grid audit failed: {checks}")
    return {"rows": rows, "checks": checks, "all_exact_checks_pass": True}


def _smith_diagonal(vectors: Sequence[Vector]) -> list[int]:
    matrix = sp.Matrix.hstack(*(sp.Matrix(vector) for vector in vectors))
    smith = smith_normal_form(matrix, domain=sp.ZZ)
    return [abs(int(smith[index, index])) for index in range(matrix.rows)]


def _full_lattice_core_audit(
    seed: dict[str, Any], patterns: dict[str, Any]
) -> dict[str, Any]:
    roots: tuple[Vector, ...] = seed["roots"]
    side_roots: tuple[Vector, ...] = patterns["_side_roots"]
    rows = []
    for endpoint_choices in product((0, 1), repeat=3):
        core = side_roots + tuple(
            roots[REPRESENTATIVE[index][choice]]
            for index, choice in enumerate(endpoint_choices)
        )
        diagonal = _smith_diagonal(core)
        rows.append(
            {
                "retained_endpoint_choices": list(endpoint_choices),
                "core_root_count": len(core),
                "Smith_diagonal": diagonal,
            }
        )
    checks = {
        "all_eight_phase_cores_have_169_roots": all(
            row["core_root_count"] == 169 for row in rows
        ),
        "all_eight_phase_cores_span_full_rank_13_lattice": all(
            row["Smith_diagonal"] == [1] * 13 for row in rows
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"full lattice core audit failed: {checks}")
    return {
        "phase_cores": rows,
        "persistence_reason": (
            "every finite box boundary contains the 166 nonendpoint seed roots "
            "and one outer endpoint root from each factor; all eight possible "
            "169-root cores have Smith diagonal 1^13"
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def _configuration_key(
    factors: Iterable[Iterable[int]],
) -> frozenset[frozenset[int]]:
    return frozenset(frozenset(factor) for factor in factors)


def _symmetry_orbit_audit(seed: dict[str, Any]) -> dict[str, Any]:
    source = seed["source"]
    walls = seed["walls"]
    fold = seed["fold"]
    involution_permutation: Permutation = seed["involution_permutation"]
    generators: tuple[sp.Matrix, ...] = source["stabilizer_generators"]
    generator_permutations = tuple(
        base._wall_permutation(matrix, walls) for matrix in generators
    )
    identity = tuple(range(316))
    words: dict[Permutation, tuple[int, ...]] = {identity: ()}
    queue: deque[Permutation] = deque((identity,))
    while queue:
        element = queue.popleft()
        for generator_index, generator in enumerate(generator_permutations):
            image = base._compose_permutations(element, generator)
            if image not in words:
                words[image] = words[element] + (generator_index,)
                queue.append(image)
    centralizer = tuple(
        (element, word)
        for element, word in words.items()
        if base._compose_permutations(element, involution_permutation)
        == base._compose_permutations(involution_permutation, element)
    )
    all_orbits = tuple(
        frozenset(int(index) - 1 for index in orbit) for orbit in fold["wall_orbits"]
    )
    orbit_index = {orbit: index for index, orbit in enumerate(all_orbits)}
    facet_orbit_indices = tuple(
        int(index) - 1 for index in fold["facet_orbit_indices"]
    )
    facet_position = {
        orbit_index_value: position
        for position, orbit_index_value in enumerate(facet_orbit_indices)
    }
    target_lookup = {
        _configuration_key(configuration): index
        for index, configuration in enumerate(TRIPLE_CONFIGURATIONS)
    }
    hits: Counter[int] = Counter()
    shortest: dict[int, tuple[int, ...]] = {}
    stabilizer_actions: Counter[tuple[tuple[int, int], ...]] = Counter()
    shortest_action_words: dict[
        tuple[tuple[int, int], ...], tuple[int, ...]
    ] = {}
    for element, word in centralizer:
        orbit_action = tuple(
            orbit_index[
                frozenset(element[wall] for wall in orbit)
            ]
            for orbit in all_orbits
        )
        facet_action = tuple(
            facet_position[orbit_action[index]] for index in facet_orbit_indices
        )
        image_factors = tuple(
            tuple(facet_action[wall] for wall in factor)
            for factor in REPRESENTATIVE
        )
        target = target_lookup.get(_configuration_key(image_factors))
        if target is not None:
            hits[target] += 1
            if target not in shortest or len(word) < len(shortest[target]):
                shortest[target] = word
        if _configuration_key(image_factors) != _configuration_key(REPRESENTATIVE):
            continue
        action: list[tuple[int, int]] = []
        for image_factor in image_factors:
            match = next(
                (
                    (factor_index, int(image_factor == tuple(reversed(factor))))
                    for factor_index, factor in enumerate(REPRESENTATIVE)
                    if set(image_factor) == set(factor)
                ),
                None,
            )
            if match is None:
                raise ArithmeticError("configuration stabilizer factor action failed")
            action.append(match)
        signature = tuple(action)
        stabilizer_actions[signature] += 1
        if signature not in shortest_action_words or len(word) < len(
            shortest_action_words[signature]
        ):
            shortest_action_words[signature] = word

    witness_rows = []
    for target in range(2):
        word = shortest[target]
        matrix = sp.eye(17)
        for generator_index in word:
            matrix = matrix * generators[generator_index]
        if not (
            matrix * seed["ambient_gram"] * matrix.T == seed["ambient_gram"]
            and seed["ample"] * matrix == seed["ample"]
            and matrix * seed["involution"] == seed["involution"] * matrix
        ):
            raise ArithmeticError("explicit triple-product conjugacy word failed")
        witness_rows.append(
            {
                "target_configuration_number": target + 1,
                "published_generator_word_one_based": [
                    index + 1 for index in word
                ],
                "ambient_matrix_sha256": base._matrix_sha256(matrix),
            }
        )
    action_rows = [
        {
            "factor_images_zero_based_and_endpoint_flip_bits": [
                list(pair) for pair in signature
            ],
            "centralizer_element_count": stabilizer_actions[signature],
            "shortest_published_generator_word_one_based": [
                index + 1 for index in shortest_action_words[signature]
            ],
        }
        for signature in sorted(stabilizer_actions)
    ]
    factor_permutations = {
        tuple(target for target, _flip in signature)
        for signature in stabilizer_actions
    }
    endpoint_flip_parities = {
        sum(flip for _target, flip in signature) % 2
        for signature in stabilizer_actions
    }
    checks = {
        "published_chamber_group_order_is_23040": len(words) == 23_040,
        "trace_nine_involution_centralizer_order_is_768": len(centralizer) == 768,
        "two_triple_products_form_one_orbit_with_384_maps_each": hits
        == Counter({0: 384, 1: 384}),
        "configuration_stabilizer_order_is_384": sum(stabilizer_actions.values())
        == 384,
        "stabilizer_realizes_24_labelled_factor_endpoint_actions": len(
            stabilizer_actions
        )
        == 24,
        "every_factor_permutation_is_realized": len(factor_permutations) == 6,
        "centralizer_stabilizer_realizes_exactly_even_endpoint_flip_parity": endpoint_flip_parities
        == {0},
        "each_labelled_action_has_16_lifts": set(stabilizer_actions.values())
        == {16},
        "explicit_exact_word_to_both_configurations": len(witness_rows) == 2,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"triple symmetry audit failed: {checks}")
    return {
        "published_chamber_group_order": len(words),
        "involution_centralizer_order": len(centralizer),
        "configuration_orbit_size": len(hits),
        "maps_from_representative_to_each_configuration": {
            str(index + 1): hits[index] for index in range(2)
        },
        "configuration_stabilizer_order": sum(stabilizer_actions.values()),
        "explicit_configuration_words": witness_rows,
        "labelled_factor_endpoint_actions": action_rows,
        "phase_and_parameter_consequence": (
            "the two listed triples are conjugate and all six permutations of "
            "the three box parameters are lattice-isometric; the certified chamber "
            "centralizer realizes exactly the even endpoint-flip subgroup, producing "
            "two phase-parity orbits rather than silently identifying all 48 labelled "
            "factor/endpoint actions"
        ),
        "checks": checks,
        "all_exact_checks_pass": True,
    }


def _public(record: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in record.items() if not key.startswith("_")}


def build(coxiter_archive: Path = DEFAULT_COXITER_ARCHIVE) -> dict[str, Any]:
    seed = base._reconstruct_seed()
    stock = base._doubling_stock(seed["fixed_gram"], seed["roots"])
    graph = _product_graph_exhaustion(stock, seed["roots"])
    structure = _triple_structure(seed, stock)
    patterns = _side_pattern_census(seed, structure)
    collision = _all_integer_triple_collision_audit(seed, structure, patterns)
    endpoints = _endpoint_orbit_audit(seed, structure, patterns)
    formula = _formula_audit(patterns)
    direct = _direct_grid_audit(seed, structure)
    lattice_core = _full_lattice_core_audit(seed, patterns)
    symmetry = _symmetry_orbit_audit(seed)
    volume_carrier = base._volume_and_carrier_audit(seed)
    coxiter_raw = base._bounded_coxiter_audit(coxiter_archive)
    coxiter = {
        **coxiter_raw,
        "family_minimum_facet_count": 172,
        "minimum_nonseed_descendant_facet_count": 296,
        "literal_seed_or_descendant_headers_absent": True,
        "scope_warning": (
            "These two pinned CoxIter lanes are not a publication-complete "
            "denominator; literal header absence is not worldwide novelty."
        ),
    }

    checks = {
        "raw_hostile_seed_reconstruction_green": all(seed["checks"].values()),
        "exact_124_wall_doubling_stock_green": all(stock["checks"].values()),
        "scoped_product_graph_exhaustion_green": all(graph["checks"].values()),
        "both_D_infinity_cubed_structures_green": all(
            all(row["checks"].values()) for row in structure["configurations"]
        ),
        "side_pattern_census_green": all(patterns["checks"].values()),
        "all_Z3_collision_proof_green": all(collision["checks"].values()),
        "six_endpoint_orbit_proof_green": all(endpoints["checks"].values()),
        "closed_formula_symbolic_proof_green": all(formula["checks"].values()),
        "direct_240_boundary_reconstructions_green": all(direct["checks"].values()),
        "full_Picard_lattice_core_persistence_green": all(
            lattice_core["checks"].values()
        ),
        "two_configuration_symmetry_orbit_green": all(symmetry["checks"].values()),
        "finite_volume_noncompact_carrier_gate_green": all(
            volume_carrier["checks"].values()
        ),
        "bounded_CoxIter_header_replay_green": all(coxiter["checks"].values()),
        "worldwide_priority_explicitly_open": coxiter["worldwide_priority"]
        == "open",
    }
    if not all(checks.values()):
        raise ArithmeticError(f"triple-parallel family failed: {checks}")

    structure_public = {
        "configurations": [_public(row) for row in structure["configurations"]],
        "all_exact_checks_pass": True,
    }
    patterns_public = _public(patterns)
    graph_public = _public(graph)
    theorem = {
        "hyperbolic_dimension": 12,
        "parameter_domain": "all integers N1,N2,N3>=1",
        "phase_parameter": (
            "sigma in {0,1}, the parity of the three endpoint-orientation flips"
        ),
        "selected_subgroup": "D_infinity^3",
        "seed_chamber_count_and_reflection_subgroup_index": "N1*N2*N3",
        "volume": "Vol(P_N1_N2_N3)=N1*N2*N3*Vol(P_seed)",
        "finite_volume": True,
        "compact": False,
        "facet_formula": formula["facet_formula"],
        "epsilon_convention": formula["epsilon_convention"],
        "phase_classification": (
            "the certified centralizer has two phase-parity orbits; when all Ni are "
            "even their facet counts differ by one, proving the two resulting "
            "polytopes nonisometric; when some Ni is odd the counts coincide and no "
            "wider isometry identification is claimed here"
        ),
        "pairwise_nonisometric_when_volume_or_facet_counts_differ": True,
        "full_fixed_Picard_lattice_for_every_parameter_triple": True,
        "primitive_integral_reflection_roots_for_every_parameter_triple": True,
        "Coxeter_reason": (
            "the index-Ni reflection subgroup of each commuting D_infinity "
            "factor has an Ni-chamber interval; their direct product box is a "
            "Coxeter fundamental chamber, and the all-Z3 audit gives its exact "
            "exterior-wall identifications"
        ),
        "finite_volume_reason": (
            "the box is the union of exactly N1*N2*N3 interiors-disjoint copies "
            "of the independently finite-volume seed chamber"
        ),
        "noncompact_reason": (
            "the box contains a closed noncompact seed chamber and retains the "
            "common ideal direction"
        ),
        "carrier": volume_carrier,
        "claim_scope": (
            "abstract finite-volume noncompact Coxeter polytopes on a weak full-"
            "Picard K3 divisor-reflection carrier; no full-nef or surface-"
            "automorphism claim"
        ),
    }
    return {
        "schema": SCHEMA,
        "status_date": STATUS_DATE,
        "all_checks_pass": all(checks.values()),
        "independence": {
            "seed_source": (
                "prior hostile raw Dolgachev--Shimada/fold reconstruction, not a "
                "production path/cube builder"
            ),
            "exploratory_scout_imported_or_invoked": False,
            "production_path_cube_builder_imported_or_invoked": False,
            "finite_parameter_interpolation_used": False,
        },
        "seed": {
            "id": "JK-H12-172",
            "dimension": 12,
            "facet_count": 172,
            "fixed_rank": 13,
            "fixed_signature": [1, 12, 0],
            "fixed_determinant": 256,
            "root_coordinate_sha256": _json_sha256(
                [list(root) for root in seed["roots"]]
            ),
            "raw_reconstruction_checks": seed["checks"],
        },
        "scoped_maximum_product_theorem": graph_public,
        "triple_parallel_structures": structure_public,
        "side_pattern_census": patterns_public,
        "all_integer_Z3_collision_audit": collision,
        "endpoint_orbit_audit": endpoints,
        "formula_audit": formula,
        "direct_grid_reconstructions": direct,
        "full_lattice_core_audit": lattice_core,
        "symmetry_orbit_audit": symmetry,
        "theorem": theorem,
        "bounded_coxiter_corpus_replay": coxiter,
        "claim_boundary": {
            "new_project_family": True,
            "worldwide_priority": "open",
            "weak_full_Picard_divisor_carrier": True,
            "full_nef_claimed": False,
            "surface_automorphism_claimed": False,
            "scoped_product_maximum_not_arbitrary_redoubling_classification": True,
        },
        "checks": checks,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coxiter-archive", type=Path, default=DEFAULT_COXITER_ARCHIVE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    payload = build(arguments.coxiter_archive)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.exists():
            raise SystemExit(f"archive missing: {arguments.output}")
        if arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive differs from exact rebuild")
        print("archive_check=PASS")
    else:
        arguments.output.write_text(rendered, encoding="utf-8")
        print(f"wrote={arguments.output}")
    print("all_checks_pass=true")
    print("family=D_infinity^3")
    print(f"payload_sha256={_sha256_bytes(rendered.encode())}")


if __name__ == "__main__":
    main()
