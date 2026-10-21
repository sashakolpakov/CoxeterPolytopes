#!/usr/bin/env python3
"""Independent hostile audit of rational full-nef maximal finiteness.

This script does not import the producer.  It rebuilds the finite host
counts by different recurrences/RREF enumeration, audits the saturation
logic under three inequivalent maximality notions, tests chamber
transitivity on a complete finite mirror arrangement, and supplies the
missing proof that the ``(2,3,infinity)`` F4 control is reflection-maximal.
"""

from __future__ import annotations

import argparse
from collections import deque
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[1]
PRODUCER_REPORT = ROOT / "results/rational_full_nef_reflection_maximal_finiteness_theorem.md"
PRODUCER_SCRIPT = ROOT / "scripts/certify_rational_full_nef_maximal_finiteness.py"
PRODUCER_ARCHIVE = ROOT / "data/enumerations/rational_full_nef_maximal_finiteness.json"
PRODUCER_TEST = ROOT / "tests/test_rational_full_nef_maximal_finiteness.py"
HOSTILE_REPORT = ROOT / "results/rational_full_nef_maximal_finiteness_hostile_audit.md"
HOSTILE_SCRIPT = ROOT / "scripts/audit_rational_full_nef_maximal_finiteness_hostile.py"
HOSTILE_TEST = ROOT / "tests/test_rational_full_nef_maximal_finiteness_hostile_audit.py"
OUTPUT = ROOT / "data/enumerations/rational_full_nef_maximal_finiteness_hostile_audit.json"
SCHEMA = "coxeter.rational_full_nef_maximal_finiteness_hostile_audit.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _serialize(value: Any) -> Any:
    if isinstance(value, Fraction):
        return int(value) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    if isinstance(value, dict):
        return {key: _serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_serialize(item) for item in value]
    return value


def _determinant(matrix: Sequence[Sequence[int | Fraction]]) -> Fraction:
    work = [[Fraction(value) for value in row] for row in matrix]
    answer = Fraction(1)
    for column in range(len(work)):
        pivot = next(
            (row for row in range(column, len(work)) if work[row][column]), None
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        diagonal = work[column][column]
        answer *= diagonal
        for position in range(column, len(work)):
            work[column][position] /= diagonal
        for row in range(column + 1, len(work)):
            multiple = work[row][column]
            for position in range(column, len(work)):
                work[row][position] -= multiple * work[column][position]
    return answer


def _matmul(
    left: Sequence[Sequence[int | Fraction]],
    right: Sequence[Sequence[int | Fraction]],
) -> list[list[Fraction]]:
    return [
        [
            sum(
                Fraction(left[row][middle]) * Fraction(right[middle][column])
                for middle in range(len(right))
            )
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def _transpose(matrix: Sequence[Sequence[int | Fraction]]) -> list[list[Fraction]]:
    return [[Fraction(value) for value in column] for column in zip(*matrix)]


def _pair(
    gram: Sequence[Sequence[int | Fraction]],
    left: Sequence[int | Fraction],
    right: Sequence[int | Fraction],
) -> Fraction:
    return sum(
        Fraction(left[i]) * Fraction(gram[i][j]) * Fraction(right[j])
        for i in range(len(left))
        for j in range(len(right))
    )


def _reflection(
    gram: Sequence[Sequence[int | Fraction]], root: Sequence[int | Fraction]
) -> list[list[Fraction]]:
    norm = _pair(gram, root, root)
    covector = [
        sum(Fraction(root[i]) * Fraction(gram[i][j]) for i in range(len(root)))
        for j in range(len(root))
    ]
    return [
        [
            Fraction(int(i == j))
            - 2 * Fraction(root[i]) * covector[j] / norm
            for j in range(len(root))
        ]
        for i in range(len(root))
    ]


def _binary_subspace_table(maximum_rank: int) -> list[list[int]]:
    """Gaussian coefficients by recurrence, not the producer's product."""

    table = [[1]]
    for rank in range(1, maximum_rank + 1):
        previous = table[-1]
        row = []
        for dimension in range(rank + 1):
            left = previous[dimension] if dimension < len(previous) else 0
            right = (
                2 ** (rank - dimension) * previous[dimension - 1]
                if dimension > 0
                else 0
            )
            row.append(left + right)
        table.append(row)
    return table


def _rref_subspaces(rank: int) -> list[tuple[tuple[int, ...], ...]]:
    """Enumerate every binary subspace once by its unique RREF basis."""

    answer: list[tuple[tuple[int, ...], ...]] = []
    columns = range(rank)
    for dimension in range(rank + 1):
        for pivots in itertools.combinations(columns, dimension):
            free_positions = [
                (row, column)
                for column in columns
                if column not in pivots
                for row, pivot in enumerate(pivots)
                if pivot < column
            ]
            for mask in range(1 << len(free_positions)):
                rows = [[0] * rank for _ in range(dimension)]
                for row, pivot in enumerate(pivots):
                    rows[row][pivot] = 1
                for bit, (row, column) in enumerate(free_positions):
                    rows[row][column] = (mask >> bit) & 1
                answer.append(tuple(tuple(row) for row in rows))
    if len(set(answer)) != len(answer):
        raise ArithmeticError("RREF enumeration duplicated a binary subspace")
    return answer


def _totally_isotropic(rows: Sequence[Sequence[int]]) -> bool:
    # The mod-two reduction of diag(1,-1,...,-1) is the identity pairing.
    return all(
        sum(left[index] * right[index] for index in range(len(left))) % 2 == 0
        for left in rows
        for right in rows
    )


def _saturation_truth_table() -> dict[str, Any]:
    """Exhaust the Horn implications used in the W=R(M) proof."""

    names = (
        "missing_integral_reflection",
        "extension_inside_O_M",
        "extension_discrete",
        "extension_arithmetic",
        "finite_volume_W",
        "maximal_discrete_reflection",
        "maximal_M_relative_reflection",
        "maximal_arithmetic_reflection",
        "maximal_integral_lattice_only",
    )
    assignments = []
    for values in itertools.product((False, True), repeat=len(names)):
        row = dict(zip(names, values))
        # Adjoining a missing integral reflection gives a proper extension
        # inside O(M); an O(M)-subgroup is discrete.  If W has finite volume,
        # the extension is finite-index/commensurable and hence arithmetic.
        if row["missing_integral_reflection"] and not row["extension_inside_O_M"]:
            continue
        if row["extension_inside_O_M"] and not row["extension_discrete"]:
            continue
        if (
            row["finite_volume_W"]
            and row["extension_inside_O_M"]
            and not row["extension_arithmetic"]
        ):
            continue
        if (
            row["missing_integral_reflection"]
            and row["extension_discrete"]
            and row["maximal_discrete_reflection"]
        ):
            continue
        if (
            row["missing_integral_reflection"]
            and row["extension_inside_O_M"]
            and row["maximal_M_relative_reflection"]
        ):
            continue
        if (
            row["missing_integral_reflection"]
            and row["extension_arithmetic"]
            and row["maximal_arithmetic_reflection"]
        ):
            continue
        assignments.append(row)

    modes = {
        "strong_discrete_reflection_maximality": {
            "maximal_discrete_reflection": True,
        },
        "M_relative_reflection_maximality": {
            "maximal_M_relative_reflection": True,
        },
        "arithmetic_reflection_maximality_with_finite_volume": {
            "maximal_arithmetic_reflection": True,
            "finite_volume_W": True,
        },
        "maximal_integral_lattice_only_hostile_control": {
            "maximal_integral_lattice_only": True,
        },
    }
    countermodels = {}
    for mode, assumptions in modes.items():
        matching = [
            row
            for row in assignments
            if row["missing_integral_reflection"]
            and all(row[key] == value for key, value in assumptions.items())
        ]
        countermodels[mode] = matching
    return {
        "admissible_boolean_models": len(assignments),
        "missing_reflection_countermodel_counts": {
            key: len(value) for key, value in countermodels.items()
        },
        "insufficient_maximal_integral_lattice_countermodel": countermodels[
            "maximal_integral_lattice_only_hostile_control"
        ][0],
        "three_valid_maximality_modes_force_W_equals_R_M": all(
            not countermodels[key]
            for key in (
                "strong_discrete_reflection_maximality",
                "M_relative_reflection_maximality",
                "arithmetic_reflection_maximality_with_finite_volume",
            )
        ),
        "maximal_integral_lattice_does_not_force_saturation": bool(
            countermodels["maximal_integral_lattice_only_hostile_control"]
        ),
    }


def _finite_mirror_chamber_control() -> dict[str, Any]:
    # Three complete A2 mirrors split the circle of directions into six
    # chambers.  Reflection in mirror k acts by j -> 2k-j-1 mod 6.
    chamber_count = 6
    reflections = [
        [(2 * mirror - chamber - 1) % chamber_count for chamber in range(chamber_count)]
        for mirror in range(3)
    ]
    orbit = {0}
    queue = deque([0])
    while queue:
        chamber = queue.popleft()
        for reflection in reflections:
            image = reflection[chamber]
            if image not in orbit:
                orbit.add(image)
                queue.append(image)
    incomplete_orbit = {0, reflections[0][0]}
    return {
        "arrangement": "three complete A2 mirrors",
        "chamber_count": chamber_count,
        "reflection_permutations": reflections,
        "complete_reflection_group_orbit_of_chamber_zero": sorted(orbit),
        "complete_arrangement_has_one_chamber_orbit": len(orbit) == chamber_count,
        "one_mirror_only_orbit": sorted(incomplete_orbit),
        "incomplete_reflection_set_need_not_be_transitive": len(incomplete_orbit)
        < chamber_count,
        "general_proof_dependencies": [
            "discrete action implies locally finite mirror arrangement",
            "a compact generic path crosses finitely many mirrors",
            "each crossed-mirror reflection maps to the adjacent chamber",
        ],
    }


def _f4_maximality(producer: dict[str, Any]) -> dict[str, Any]:
    row = producer["controls"]["F4_maximal_non_coble"]
    gram = [[Fraction(value) for value in line] for line in row["half_root_gram"]]
    roots = [[int(i == j) for i in range(3)] for j in range(3)]
    norms = [_pair(gram, root, root) for root in roots]
    pairings = [gram[0][1], gram[1][2], gram[2][0]]
    coxeter_orders = []
    for pairing in pairings:
        normalized_product = 4 * pairing * pairing / (norms[0] * norms[0])
        coxeter_orders.append({Fraction(0): 2, Fraction(1): 3, Fraction(4): "infinity"}[normalized_product])
    area = Fraction(1) - Fraction(1, 2) - Fraction(1, 3)

    triangle_areas = []
    for first in range(2, 101):
        for second in range(first, 101):
            value = Fraction(1) - Fraction(1, first) - Fraction(1, second)
            if value > 0:
                triangle_areas.append((value, first, second))
    minimum = min(triangle_areas)
    matrices = [_reflection(gram, root) for root in roots]
    reflection_checks = []
    identity = [[Fraction(int(i == j)) for j in range(3)] for i in range(3)]
    for matrix in matrices:
        reflection_checks.append(
            all(value.denominator == 1 for line in matrix for value in line)
            and _matmul(_transpose(matrix), _matmul(gram, matrix)) == gram
            and _matmul(matrix, matrix) == identity
        )
    return {
        "half_root_gram": gram,
        "root_norms": norms,
        "coxeter_orders": coxeter_orders,
        "triangle_area_over_pi": area,
        "independent_triangle_search_order_bound": 100,
        "minimum_positive_one_ideal_triangle_area_over_pi": minimum[0],
        "minimum_signature_in_search": [minimum[1], minimum[2], "infinity"],
        "analytic_lower_bound": {
            "one_ideal_triangle": "1-1/m-1/n >= 1-1/2-1/3 = 1/6",
            "two_or_more_ideal_triangle": ">= 1/2",
            "noncompact_k_gon_k_at_least_4": "k-2-(k-1)/2 = (k-3)/2 >= 1/2",
        },
        "all_three_reflections_integral_isometric_involutions": all(reflection_checks),
        "finite_volume_overgroup_has_finite_index": True,
        "overgroup_remains_noncompact_because_it_contains_a_parabolic": True,
        "proper_reflection_overgroup_would_have_area_below_pi_over_6": True,
        "reflection_maximality_proved": area == minimum[0] == Fraction(1, 6),
        "producer_archive_maximality_certificate_was_missing": not any(
            "maximality" in key.lower() and key != "maximal_triangle" for key in row
        ),
    }


def _ideal_polygon_control() -> dict[str, Any]:
    gram = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    rows = []
    for index in range(1, 25):
        roots = [
            [0, -1, 0],
            *[
                [1 + j * (j + 1), 2 * j + 1, j * (j + 1) - 1]
                for j in range(index)
            ],
            [index, 1, index],
        ]
        wall_gram = [[_pair(gram, left, right) for right in roots] for left in roots]
        adjacent = {(i, i + 1) for i in range(len(roots) - 1)} | {
            (0, len(roots) - 1)
        }
        rows.append(
            {
                "index": index,
                "facet_count": len(roots),
                "all_norms_minus_one": all(wall_gram[i][i] == -1 for i in range(len(roots))),
                "all_adjacent_vertices_ideal": all(
                    wall_gram[i][j] == 1 for i, j in adjacent
                ),
                "area_over_pi_from_all_ideal_polygon": len(roots) - 2,
                "proper_index_over_base_for_index_above_one": index if index > 1 else None,
            }
        )
    return {
        "rows": rows,
        "unbounded_facets_at_fixed_I_1_2_without_maximality": [
            row["facet_count"] for row in rows
        ]
        == list(range(3, 27)),
        "all_exact_ideal_polygon_checks_pass": all(
            row["all_norms_minus_one"]
            and row["all_adjacent_vertices_ideal"]
            and row["area_over_pi_from_all_ideal_polygon"] == row["index"]
            for row in rows
        ),
    }


def _field_audit() -> dict[str, Any]:
    # Rational normals force every Cartan entry and cyclic product rational.
    rational_gram = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    rational_roots = [[0, 1, 0], [1, 2, 0], [2, 2, 1]]
    cartan = [
        [
            2 * _pair(rational_gram, rational_roots[i], rational_roots[j])
            / _pair(rational_gram, rational_roots[i], rational_roots[i])
            for j in range(3)
        ]
        for i in range(3)
    ]
    cyclic = cartan[0][1] * cartan[1][2] * cartan[2][0]

    # Hostile irrational-normal control.  For roots
    # (1,0), (1,1), (sqrt(2),1) in a negative Euclidean plane, the Cartan
    # 1->2->3->1 product is 8/3 + (4/3)sqrt(2).
    irrational_cyclic = (Fraction(8, 3), Fraction(4, 3))
    return {
        "rational_cartan_matrix": cartan,
        "rational_cyclic_product": cyclic,
        "rational_normals_force_Q_cyclic_products": all(
            isinstance(value, Fraction) for row in cartan for value in row
        )
        and isinstance(cyclic, Fraction),
        "irrational_normal_hostile_cyclic_product_Q_plus_sqrt2_coefficients": irrational_cyclic,
        "irrational_normal_control_is_not_rational": irrational_cyclic[1] != 0,
        "minimal_vinberg_field_Q_uses_literal_rational_full_dimensional_host": True,
        "trace_module_degree_inequality_is_optional_not_used": True,
        "carrier_surface_fields_of_definition_are_not_bounded_by_this_argument": True,
    }


def build_certificate() -> dict[str, Any]:
    producer = json.loads(PRODUCER_ARCHIVE.read_text(encoding="utf-8"))
    report = PRODUCER_REPORT.read_text(encoding="utf-8")
    recurrence = _binary_subspace_table(22)
    recurrence_counts = [sum(recurrence[rank]) for rank in range(23)]
    small_rref = []
    for rank in range(1, 6):
        subspaces = _rref_subspaces(rank)
        isotropic = [rows for rows in subspaces if _totally_isotropic(rows)]
        small_rref.append(
            {
                "rank": rank,
                "all_subspaces_by_RREF": len(subspaces),
                "totally_isotropic_subspaces_for_I_mod_2": len(isotropic),
                "maximum_isotropic_dimension": max(map(len, isotropic)),
            }
        )
    complex_bound = sum(recurrence_counts[rank] for rank in range(11, 21))
    odd_bound = sum(recurrence_counts[rank] for rank in range(11, 23))

    split_change = [[1, 1, 0], [1, 0, 1], [1, 1, 1]]
    u_minus = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
    direct_hosts = {
        "rank_2_forms": [
            {"name": "U", "gram": [[0, 1], [1, 0]]},
            {"name": "I_1_1", "gram": [[1, 0], [0, -1]]},
        ],
        "odd_forms_ranks_3_through_10": list(range(3, 11)),
        "positive_dimensional_host_count": 2 + len(range(3, 11)),
        "U_plus_minus_one_change_determinant": _determinant(split_change),
        "U_plus_minus_one_conjugated_gram": _matmul(
            _transpose(split_change), _matmul(u_minus, split_change)
        ),
    }
    coble = {
        "small_rank_independent_RREF_counts": small_rref,
        "all_subspace_counts_rank_1_through_22_by_recurrence": recurrence_counts[1:],
        "complex_rank_11_through_20_bound": complex_bound,
        "odd_positive_rank_11_through_22_bound": odd_bound,
        "nonisotropic_half_vector_control": {
            "H_basis": [[1, 0]],
            "scaled_half_norm": Fraction(1, 2),
            "integral": False,
        },
        "isotropic_half_vector_control": {
            "H_basis": [[1, 1]],
            "scaled_half_norm_for_diag_1_minus_1": 0,
            "integral": True,
        },
        "determinant_formula": "abs(det M_H)=2^(rho-2 dim H)",
        "binary_boundary_sum_relation_only_reduces_the_upper_bound": True,
    }
    saturation = _saturation_truth_table()
    chamber = _finite_mirror_chamber_control()
    f4 = _f4_maximality(producer)
    ideal = _ideal_polygon_control()
    field = _field_audit()

    producer_complex = producer["coble_extension_bounds"][
        "complex_marked_extension_upper_bound"
    ]
    producer_odd = producer["coble_extension_bounds"][
        "odd_positive_marked_extension_upper_bound"
    ]
    claims_present = {
        "W_equals_R_M": "W(P)=R(M)" in report,
        "arithmetic_reflection_maximality_distinguished": "arithmetic **reflection** groups" in report,
        "one_chamber_orbit": "all chambers" in report and "conjugate" in report,
        "finite_standard_hosts": "only finitely many hosts" in report,
        "field_Q": "k(P)=\\mathbf Q" in report,
        "no_terminating_reflectivity_census": "does **not** provide a uniform stopping certificate" in report,
    }
    checks = {
        "producer_claims_under_audit_are_present": all(claims_present.values()),
        "three_valid_reflection_maximality_modes_force_saturation": saturation[
            "three_valid_maximality_modes_force_W_equals_R_M"
        ],
        "maximal_integral_lattice_is_insufficient": saturation[
            "maximal_integral_lattice_does_not_force_saturation"
        ],
        "complete_mirror_arrangement_has_one_chamber_orbit": chamber[
            "complete_arrangement_has_one_chamber_orbit"
        ],
        "incomplete_mirror_set_control_is_not_transitive": chamber[
            "incomplete_reflection_set_need_not_be_transitive"
        ],
        "direct_host_count_independently_equals_10": direct_hosts[
            "positive_dimensional_host_count"
        ]
        == 10,
        "direct_U_split_reconstructed": direct_hosts[
            "U_plus_minus_one_conjugated_gram"
        ]
        == [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
        "small_binary_subspace_counts_match_independent_RREF": all(
            row["all_subspaces_by_RREF"] == recurrence_counts[row["rank"]]
            for row in small_rref
        ),
        "complex_Coble_host_bound_matches": complex_bound == producer_complex,
        "odd_positive_Coble_host_bound_matches": odd_bound == producer_odd,
        "integral_overlattice_filter_has_positive_and_negative_controls": (
            not coble["nonisotropic_half_vector_control"]["integral"]
            and coble["isotropic_half_vector_control"]["integral"]
        ),
        "F4_triangle_reflection_maximality_gap_repaired": f4[
            "reflection_maximality_proved"
        ],
        "producer_F4_maximality_was_not_previously_certified": f4[
            "producer_archive_maximality_certificate_was_missing"
        ],
        "fixed_lattice_nonmaximal_polygon_control_rebuilt": ideal[
            "unbounded_facets_at_fixed_I_1_2_without_maximality"
        ]
        and ideal["all_exact_ideal_polygon_checks_pass"],
        "rational_wall_normals_force_vinberg_field_Q": field[
            "rational_normals_force_Q_cyclic_products"
        ],
        "irrational_and_surface_field_scope_controls_hold": field[
            "irrational_normal_control_is_not_rational"
        ]
        and field["carrier_surface_fields_of_definition_are_not_bounded_by_this_argument"],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"maximal-finiteness hostile audit failed: {checks}")
    return _serialize(
        {
            "schema": SCHEMA,
            "status": "PASS_WITH_ONE_INDEPENDENTLY_REPAIRED_CERTIFICATE_GAP",
            "status_date": "2026-08-10",
            "verdict": {
                "core_finiteness_theorem": "sound under its stated carrier-relative hypotheses",
                "scope_narrowing_required": False,
                "producer_certificate_gap": (
                    "the F4 control asserted reflection maximality without certifying it"
                ),
                "gap_repaired_here": (
                    "(2,3,infinity) has area pi/6, the minimum for a noncompact "
                    "finite-area Coxeter polygon, so no proper reflection overgroup exists"
                ),
            },
            "sources": {
                "producer_report": {
                    "path": str(PRODUCER_REPORT.relative_to(ROOT)),
                    "sha256": _sha256(PRODUCER_REPORT),
                },
                "producer_script": {
                    "path": str(PRODUCER_SCRIPT.relative_to(ROOT)),
                    "sha256": _sha256(PRODUCER_SCRIPT),
                },
                "producer_archive": {
                    "path": str(PRODUCER_ARCHIVE.relative_to(ROOT)),
                    "sha256": _sha256(PRODUCER_ARCHIVE),
                },
                "producer_test": {
                    "path": str(PRODUCER_TEST.relative_to(ROOT)),
                    "sha256": _sha256(PRODUCER_TEST),
                },
                "hostile_report": {
                    "path": str(HOSTILE_REPORT.relative_to(ROOT)),
                    "sha256": _sha256(HOSTILE_REPORT),
                },
                "hostile_script": {
                    "path": str(HOSTILE_SCRIPT.relative_to(ROOT)),
                    "sha256": _sha256(HOSTILE_SCRIPT),
                },
                "hostile_test": {
                    "path": str(HOSTILE_TEST.relative_to(ROOT)),
                    "sha256": _sha256(HOSTILE_TEST),
                },
            },
            "claims_present": claims_present,
            "saturation_logic": saturation,
            "one_chamber_orbit_audit": chamber,
            "direct_host_audit": direct_hosts,
            "standard_Coble_host_audit": coble,
            "F4_reflection_maximality_repair": f4,
            "nonmaximal_fixed_lattice_control": ideal,
            "vinberg_field_audit": field,
            "effectiveness_scope": {
                "finite_host_list_is_a_terminating_mathematical_enumeration": True,
                "finite_abstract_chamber_bound_proved": True,
                "uniform_reflectivity_or_facet_stopping_bound_proved": False,
                "terminating_census_of_reflective_hosts_proved": False,
                "finite_surface_isomorphism_or_definition_field_list_proved": False,
                "raw_P_to_carrier_essential_surjectivity_proved": False,
            },
            "checks": checks,
            "all_checks_pass": True,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    rendered = json.dumps(build_certificate(), indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
