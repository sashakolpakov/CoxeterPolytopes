#!/usr/bin/env python3
"""Certify the first generated PB-Atlas Stage-B2b geometric component.

The calibration packet is selected only by its proximity-forest invariant.
The program enumerates the finite free-chart cover, generates every degree-one
root determinant from the divisor table, solves the surviving incidence ideal,
and lets the cluster-collision test choose the good component.  It then
constructs every remaining facet section and the unique bi-anticanonical
boundary through the public Stage-B2a API.

No stored point coordinates, old incidence equation, surface equation, Coble
type, or accepted component is read by the construction.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Mapping, Sequence

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_pb_atlas_cluster_chart_stage_b2a import (
    _all_p2_packets,
    _matrix_from_record,
    _unique_packet_by_forest,
)
from scripts.pb_atlas_cluster_chart_stage_b2a import (
    compile_p2_cluster_chart,
    compile_p2_divisor_evaluation_matrix,
)
from scripts.pb_atlas_geometry_stage_b2b import (
    compile_closed_section_geometry,
    closed_linear_system_rank,
)


OUTPUT = ROOT / "data/enumerations/pb_atlas_h10_geometry_stage_b2b.json"
SCHEMA = "coxeter.pb_atlas.h10_geometry_stage_b2b.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _root_incidence_rows(packet: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    rows = [
        row
        for row in packet["divisor_table"]
        if not row["is_contracted_wall"]
        and row["all_multiplicities_nonnegative"]
        and row["self_intersection"] == -2
        and row["plane_degree"] == 1
    ]
    if len(rows) != 9:
        raise ArithmeticError("the structural plane packet does not have nine line roots")
    return rows


def _matrix_determinant(
    packet: Mapping[str, Any], row: Mapping[str, Any], options: Mapping[str, Any]
) -> sp.Expr:
    record = compile_p2_divisor_evaluation_matrix(
        packet,
        degree=1,
        multiplicities=row["multiplicities_chronological"],
        chart_options=options,
    )
    matrix = _matrix_from_record(record)
    if matrix.shape != (3, 3):
        raise ArithmeticError("a line-root incidence matrix is not 3 by 3")
    return sp.factor(matrix.det())


def _incidence_chart_census(packet: Mapping[str, Any]) -> dict[str, Any]:
    centers = packet["reverse_blowup_program"]["center_records"]
    free = [
        index for index, center in enumerate(centers) if center["center_kind"] == "free"
    ]
    rows = _root_incidence_rows(packet)
    records: list[dict[str, Any]] = []
    viable: list[dict[str, Any]] = []
    for modes_tuple in itertools.product(("A", "B"), repeat=len(free)):
        modes = dict(zip(free, modes_tuple))
        options = {"free_chart_mode_by_center": modes}
        determinants = [
            (int(row["wall_index_zero_based"]), _matrix_determinant(packet, row, options))
            for row in rows
        ]
        unit_obstructions = [
            wall
            for wall, determinant in determinants
            if not determinant.free_symbols and determinant != 0
        ]
        record = {
            "free_chart_modes": {str(key): value for key, value in modes.items()},
            "root_determinants": {
                str(wall): str(determinant) for wall, determinant in determinants
            },
            "constant_nonzero_obstruction_walls": unit_obstructions,
            "incidence_possible_on_chart": not unit_obstructions,
        }
        records.append(record)
        if not unit_obstructions:
            viable.append({"record": record, "options": options, "equations": [value for _wall, value in determinants if value != 0]})
    if len(records) != 8 or len(viable) != 1:
        raise ArithmeticError("the finite free-chart incidence census changed")

    survivor = viable[0]
    compiled = compile_p2_cluster_chart(packet, survivor["options"])
    variables = [sp.Symbol(name) for name in compiled["chart"]["coordinate_variables"]]
    groebner = sp.groebner(survivor["equations"], *variables, order="lex")
    solutions = sp.solve(survivor["equations"], variables, dict=True)
    if len(solutions) != 2 or any(set(solution) != set(variables) for solution in solutions):
        raise ArithmeticError("the root incidence ideal did not produce two closed components")
    return {
        "all_chart_records": records,
        "viable_chart_count": len(viable),
        "surviving_chart_options": survivor["options"],
        "coordinate_variables": [str(value) for value in variables],
        "groebner_basis": [str(sp.factor(value.as_expr())) for value in groebner.polys],
        "solutions": solutions,
    }


def _normalized_projective(values: Sequence[sp.Expr]) -> tuple[sp.Rational, ...]:
    rational = [sp.Rational(value) for value in values]
    for value in rational:
        if value:
            return tuple(sp.cancel(entry / value) for entry in rational)
    raise ArithmeticError("a specialized proper point is the zero vector")


def _proper_collisions(
    packet: Mapping[str, Any],
    options: Mapping[str, Any],
    solution: Mapping[sp.Symbol, sp.Expr],
) -> list[list[int]]:
    chart = compile_p2_cluster_chart(packet, options)["chart"]
    points: list[tuple[int, tuple[sp.Rational, ...]]] = []
    for record in chart["center_charts"]:
        if record["center_kind"] != "proper":
            continue
        coordinates = [
            sp.sympify(value).subs(solution, simultaneous=True)
            for value in record["homogeneous_center_coordinates"]
        ]
        points.append(
            (
                int(record["center_position_zero_based"]),
                _normalized_projective(coordinates),
            )
        )
    return [
        [left_index, right_index]
        for (left_index, left), (right_index, right) in itertools.combinations(points, 2)
        if left == right
    ]


def _selected_chart_open_audit(
    packet: Mapping[str, Any],
    options: Mapping[str, Any],
    solution: Mapping[sp.Symbol, sp.Expr],
) -> dict[str, Any]:
    chart = compile_p2_cluster_chart(packet, options)["chart"]
    proper: dict[int, tuple[sp.Rational, ...]] = {}
    affine_patch_values: dict[str, str] = {}
    for record in chart["center_charts"]:
        if record["center_kind"] != "proper":
            continue
        position = int(record["center_position_zero_based"])
        coordinates = tuple(
            sp.Rational(sp.sympify(value).subs(solution, simultaneous=True))
            for value in record["homogeneous_center_coordinates"]
        )
        proper[position] = coordinates
        patch = int(record["proper_affine_patch_coordinate_zero_based"])
        affine_patch_values[str(position)] = str(coordinates[patch])

    gauge = [
        int(value)
        for value in chart["pgl3_gauge"][
            "selected_proper_center_positions_zero_based"
        ]
    ]
    gauge_determinants = {
        "_".join(str(value) for value in triple): str(
            sp.Matrix.hstack(*(sp.Matrix(proper[value]) for value in triple)).det()
        )
        for triple in itertools.combinations(gauge, 3)
    }
    passes = (
        all(sp.Rational(value) != 0 for value in affine_patch_values.values())
        and all(sp.Rational(value) != 0 for value in gauge_determinants.values())
    )
    return {
        "selected_affine_patch_coordinate_values": affine_patch_values,
        "selected_PGL3_gauge_triple_determinants": gauge_determinants,
        "every_declared_chart_open_factor_nonzero": passes,
    }


def _singular_absolute_irreducible(polynomial: str, degree: int) -> dict[str, Any]:
    if degree == 1:
        return {
            "engine": "degree-one theorem",
            "absolute_factor_count_with_unit": 2,
            "absolute_component_count": 1,
            "sole_nontrivial_factor_multiplicity": 1,
            "absolutely_irreducible": True,
        }
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for the absolute-factorization gate")
    singular_polynomial = (
        polynomial.replace("X", "x").replace("Y", "y").replace("Z", "z")
    )
    program = f'''
ring r=0,(x,y,z),dp;
poly f={singular_polynomial};
LIB "absfact.lib";
def ABSR=absFactorize(f);
setring ABSR;
print("ABS_FACTOR_COUNT_WITH_UNIT"); print(size(absolute_factors[1]));
print("ABS_COMPONENT_COUNT"); print(absolute_factors[4]);
print("ABS_SOLE_NONTRIVIAL_FACTOR_MULTIPLICITY"); print(absolute_factors[2][2]);
exit;
'''
    completed = subprocess.run(
        [executable, "-q", "--no-rc"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    parsed: dict[str, int] = {}
    for index, line in enumerate(lines[:-1]):
        if line == "ABS_FACTOR_COUNT_WITH_UNIT":
            parsed["absolute_factor_count_with_unit"] = int(lines[index + 1])
        elif line == "ABS_COMPONENT_COUNT":
            parsed["absolute_component_count"] = int(lines[index + 1])
        elif line == "ABS_SOLE_NONTRIVIAL_FACTOR_MULTIPLICITY":
            parsed["sole_nontrivial_factor_multiplicity"] = int(lines[index + 1])
    irreducible = parsed == {
        "absolute_factor_count_with_unit": 2,
        "absolute_component_count": 1,
        "sole_nontrivial_factor_multiplicity": 1,
    }
    if not irreducible:
        raise ArithmeticError(
            f"a generated facet polynomial is not absolutely irreducible: {parsed}"
        )
    return {"engine": "Singular absfact over Q", **parsed, "absolutely_irreducible": True}


def _singular_engine_record() -> dict[str, str]:
    executable = shutil.which("Singular")
    if executable is None:
        raise RuntimeError("Singular is required for the absolute-factorization gate")
    completed = subprocess.run(
        [executable, "--version"],
        text=True,
        capture_output=True,
        check=True,
        timeout=15,
    )
    first_line = completed.stdout.splitlines()[0].strip()
    return {
        "executable": executable,
        "version_first_line": first_line,
        "library": "absfact.lib",
        "calculation_domain": "absolute factorization of homogeneous polynomials over Q",
    }


def _closed_component_geometry(
    packet: Mapping[str, Any],
    options: Mapping[str, Any],
    solution: Mapping[sp.Symbol, sp.Expr],
) -> dict[str, Any]:
    coordinates = {str(key): value for key, value in solution.items()}
    minus_K = closed_linear_system_rank(
        packet,
        system="minus_K",
        chart_options=options,
        closed_coordinates=coordinates,
    )
    minus_two_K = closed_linear_system_rank(
        packet,
        system="minus_two_K",
        chart_options=options,
        closed_coordinates=coordinates,
    )
    boundary = compile_closed_section_geometry(
        packet,
        degree=6,
        multiplicities=[2] * len(packet["reverse_blowup_program"]["center_records"]),
        chart_options=options,
        closed_coordinates=coordinates,
    )
    boundary["absolute_factorization"] = _singular_absolute_irreducible(
        boundary["homogeneous_polynomial"], 6
    )

    boundary_root_rows = [
        row
        for row in packet["divisor_table"]
        if not row["is_contracted_wall"]
        and row["anti_bicanonical_boundary_weight"]
    ]
    if len(boundary_root_rows) != 1:
        raise ArithmeticError("the packet does not have exactly one physical boundary wall")
    boundary_root = boundary_root_rows[0]
    if not (
        boundary_root["anti_bicanonical_boundary_weight"] == 2
        and boundary_root["plane_degree"] == 3
        and boundary_root["multiplicities_chronological"]
        == [1] * len(packet["reverse_blowup_program"]["center_records"])
        and boundary_root["self_intersection"] == -1
        and boundary_root["canonical_intersection"] == 1
    ):
        raise ArithmeticError("the physical boundary is not twice its declared wall root")
    boundary.update(
        {
            "wall_index_zero_based": int(boundary_root["wall_index_zero_based"]),
            "physical_curve_is_twice_wall_root": True,
            "wall_root_plane_degree": int(boundary_root["plane_degree"]),
            "wall_root_self_intersection": int(boundary_root["self_intersection"]),
            "physical_curve_self_intersection": -4,
            "physical_curve_canonical_intersection": 2,
        }
    )

    divisor_rows = [
        row
        for row in packet["divisor_table"]
        if not row["is_contracted_wall"]
        and row["all_multiplicities_nonnegative"]
        and not row["anti_bicanonical_boundary_weight"]
    ]
    curves = []
    for row in divisor_rows:
        geometry = compile_closed_section_geometry(
            packet,
            degree=int(row["plane_degree"]),
            multiplicities=row["multiplicities_chronological"],
            chart_options=options,
            closed_coordinates=coordinates,
        )
        geometry["absolute_factorization"] = _singular_absolute_irreducible(
            geometry["homogeneous_polynomial"], int(row["plane_degree"])
        )
        if row["riemann_roch"]["arithmetic_genus"] != 0:
            raise ArithmeticError("a generated nonboundary facet does not have genus zero")
        curves.append(
            {
                "wall_index_zero_based": int(row["wall_index_zero_based"]),
                "self_intersection": int(row["self_intersection"]),
                "canonical_intersection": int(row["canonical_intersection"]),
                "arithmetic_genus": 0,
                "smooth_geometrically_integral_rational_by_exact_multiplicity_and_genus_zero": True,
                **geometry,
            }
        )

    contracted = [
        row for row in packet["divisor_table"] if row["is_contracted_wall"]
    ]
    center_records = packet["reverse_blowup_program"]["center_records"]
    contracted_wall_indices = {
        int(row["wall_index_zero_based"]) for row in contracted
    }
    exceptional_source_indices = {
        int(center["source_contracted_wall_index_zero_based"])
        for center in center_records
    }
    if (
        len(contracted) != len(center_records)
        or len(contracted_wall_indices) != len(contracted)
        or contracted_wall_indices != exceptional_source_indices
    ):
        raise ArithmeticError("the contracted walls are not exactly the blowup exceptional curves")
    contracted_by_index = {
        int(row["wall_index_zero_based"]): row for row in contracted
    }
    exceptional_coefficient_checks = []
    for center in center_records:
        wall_index = int(center["source_contracted_wall_index_zero_based"])
        row = contracted_by_index[wall_index]
        strict = [
            int(value)
            for value in center[
                "strict_exceptional_coordinates_in_chronological_total_basis"
            ]
        ]
        exceptional_coefficients = [
            int(value) for value in row["exceptional_coefficients_chronological"]
        ]
        multiplicities = [
            int(value) for value in row["multiplicities_chronological"]
        ]
        passes = (
            exceptional_coefficients == strict
            and multiplicities == [-value for value in strict]
        )
        exceptional_coefficient_checks.append(
            {
                "center_position_zero_based": int(
                    center["center_position_zero_based"]
                ),
                "wall_index_zero_based": wall_index,
                "strict_exceptional_coordinates": strict,
                "row_exceptional_coefficients": exceptional_coefficients,
                "row_multiplicities": multiplicities,
                "coefficientwise_match": passes,
            }
        )
    if not all(row["coefficientwise_match"] for row in exceptional_coefficient_checks):
        raise ArithmeticError("a contracted wall is not its source strict exceptional curve")
    noncontracted_indices = {
        int(curve["wall_index_zero_based"]) for curve in curves
    }
    boundary_index = int(boundary_root["wall_index_zero_based"])
    all_table_indices = {
        int(row["wall_index_zero_based"]) for row in packet["divisor_table"]
    }
    partitioned_indices = contracted_wall_indices | noncontracted_indices | {
        boundary_index
    }
    if (
        len(all_table_indices) != len(packet["divisor_table"])
        or contracted_wall_indices & noncontracted_indices
        or boundary_index in contracted_wall_indices | noncontracted_indices
        or partitioned_indices != all_table_indices
    ):
        raise ArithmeticError("the constructed facets do not partition the supplied wall list")
    physical_wall_count = len(contracted) + len(curves) + 1
    return {
        "closed_coordinates": {str(key): str(value) for key, value in solution.items()},
        "minus_K": minus_K,
        "minus_two_K": minus_two_K,
        "physical_boundary": boundary,
        "noncontracted_nonboundary_facets": curves,
        "contracted_exceptional_facets": {
            "count": len(contracted),
            "wall_indices_zero_based": sorted(contracted_wall_indices),
            "source_wall_indices_from_reverse_blowups_zero_based": sorted(
                exceptional_source_indices
            ),
            "coefficientwise_source_exceptional_checks": exceptional_coefficient_checks,
            "every_source_exceptional_matches_coefficientwise": True,
            "smooth_irreducible_by_reverse_blowup_construction": True,
        },
        "physical_wall_count": physical_wall_count,
        "constructed_wall_indices_partition_supplied_wall_list": True,
        "every_proposed_facet_effective_irreducible": (
            physical_wall_count == len(packet["divisor_table"])
            and all(curve["all_checks_pass"] for curve in curves)
            and all(
                curve["absolute_factorization"]["absolutely_irreducible"]
                for curve in curves
            )
            and boundary["all_checks_pass"]
            and boundary["absolute_factorization"]["absolutely_irreducible"]
        ),
        "strong_Coble": (
            minus_K["kernel_dimension"] == 0
            and minus_two_K["kernel_dimension"] == 1
            and boundary["every_declared_multiplicity_is_exact"]
            and boundary["absolute_factorization"]["absolutely_irreducible"]
        ),
    }


def build_certificate() -> dict[str, Any]:
    source = _unique_packet_by_forest(
        _all_p2_packets(), [2, 2, 2, 1, 1, 1, 1]
    )
    packet = source["packet"]
    incidence = _incidence_chart_census(packet)
    options = incidence.pop("surviving_chart_options")
    component_records = []
    good = []
    for solution in incidence.pop("solutions"):
        collisions = _proper_collisions(packet, options, solution)
        chart_open = _selected_chart_open_audit(packet, options, solution)
        record = {
            "coordinates": {str(key): str(value) for key, value in solution.items()},
            "proper_center_collision_pairs_zero_based": collisions,
            "collision_free": not collisions,
            "selected_chart_open_audit": chart_open,
        }
        component_records.append(record)
        if not collisions and chart_open["every_declared_chart_open_factor_nonzero"]:
            good.append(solution)
    if len(good) != 1:
        raise ArithmeticError("collision avoidance did not select one incidence component")
    geometry = _closed_component_geometry(packet, options, good[0])
    if geometry["physical_wall_count"] != 30:
        raise ArithmeticError("the generated physical facet count changed")

    centres = packet["reverse_blowup_program"]["center_records"]
    satellite_centres = [
        int(center["center_position_zero_based"])
        for center in centres
        if center["center_kind"] == "satellite"
    ]
    if satellite_centres:
        raise ArithmeticError(
            "the H10 calibration unexpectedly contains satellite centres; "
            "their collision gate is not part of this certificate"
        )
    sibling_free_collision_debt = any(
        sum(
            child["center_kind"] == "free"
            and child["parent_center_positions_zero_based"] == [parent]
            for child in centres
        )
        > 1
        for parent in range(len(centres))
    )
    if sibling_free_collision_debt:
        raise ArithmeticError("the calibration has an unchecked same-exceptional collision")

    checks = {
        "packet_selected_by_proximity_forest_not_name": True,
        "all_eight_free_chart_modes_exhausted": len(incidence["all_chart_records"]) == 8,
        "one_viable_chart_on_the_deterministic_PGL3_open": incidence["viable_chart_count"] == 1,
        "root_incidence_has_exactly_two_closed_components": len(component_records) == 2,
        "collision_open_selects_exactly_one_component": len(good) == 1,
        "selected_incidence_component_lies_in_declared_chart_open": all(
            record["selected_chart_open_audit"][
                "every_declared_chart_open_factor_nonzero"
            ]
            for record in component_records
            if record["collision_free"]
        ),
        "anticanonical_system_empty": geometry["minus_K"]["kernel_dimension"] == 0,
        "antibicanonical_system_unique": geometry["minus_two_K"]["kernel_dimension"] == 1,
        "strong_Coble_constructed": geometry["strong_Coble"],
        "all_thirty_facets_constructed_effective_irreducible": geometry[
            "every_proposed_facet_effective_irreducible"
        ],
        "contracted_wall_set_equals_reverse_blowup_exceptional_set": (
            geometry["contracted_exceptional_facets"]["wall_indices_zero_based"]
            == geometry["contracted_exceptional_facets"][
                "source_wall_indices_from_reverse_blowups_zero_based"
            ]
        ),
        "no_satellite_or_same_parent_free_collision_debt": (
            not satellite_centres and not sibling_free_collision_debt
        ),
        "no_old_coordinate_or_surface_equation_read": True,
        "no_same_P_or_maximality_claim_in_this_certificate": True,
    }
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-10",
        "status": "PROVED_GENERATED_CLOSED_COMPONENT_STRONG_COBLE_AND_ALL_FACET_EFFECTIVITY",
        "input_policy": {
            "reads": [
                "one Stage-B1 numerical P2 packet",
                "its chronological cluster and divisor table",
                "Stage-B2a exact chart and jet compiler",
            ],
            "does_not_read": [
                "surface or Coble type name",
                "stored point coordinates",
                "stored incidence equation or accepted factor",
                "stored curve or boundary equation",
                "Kondo or Nikulin-Vinberg classification",
            ],
            "calibration_selector": "unique P2 packet with proximity component sizes [2,2,2,1,1,1,1]",
        },
        "generated_incidence": {
            **incidence,
            "component_records": component_records,
            "collision_avoidance_is_generated_not_hard_coded": True,
        },
        "geometric_component": geometry,
        "exact_algebra_engines": {
            "incidence_and_linear_algebra": "SymPy exact QQ arithmetic",
            "absolute_factorization": _singular_engine_record(),
        },
        "dimension_bridge": {
            "general_statement": (
                "Any raw packet for which this characteristic-zero finite gate constructs a strong Coble surface "
                "of Picard rank n+1 has a smooth canonical K3 double cover containing the adapted rank-(n+1) lattice; "
                "therefore n<=19 over C."
            ),
            "current_compiler_base_field": "Q, with geometric conclusion over C",
            "odd_characteristic_analogue": (
                "An independently constructed strong Coble carrier in characteristic not two has n<=21; "
                "the present QQ/Singular gate does not construct such a carrier."
            ),
            "does_not_require_same_nef_cone_for_rank_bound": True,
            "reflection_maximal_inputs_get_same_P_from_the_independent_exact_return_theorem": True,
            "current_calibration_picard_rank": 11,
            "current_calibration_hyperbolic_dimension": 10,
            "current_calibration_not_claimed_to_be_reflection_maximal_here": True,
        },
        "scope_boundary": {
            "proved": "one generated closed component, all physical facets, strong-Coble existence, and the class-relative K3 rank implication",
            "not_proved": [
                "a second packet passing the same complete B2b grammar",
                "all Stage-B2a charts or ruled endpoints",
                "reflection maximality of this calibration",
                "global raw arithmetic-chamber essential surjectivity",
            ],
            "rating": "strict generated B2b positive and reusable dimension bridge; not global L4",
        },
        "sources": {
            "scripts/pb_atlas_cluster_chart_stage_b2a.py": _sha256(
                ROOT / "scripts/pb_atlas_cluster_chart_stage_b2a.py"
            ),
            "scripts/pb_atlas_geometry_stage_b2b.py": _sha256(
                ROOT / "scripts/pb_atlas_geometry_stage_b2b.py"
            ),
            "data/enumerations/pb_atlas_numerical_cluster_stage_b1.json": _sha256(
                ROOT / "data/enumerations/pb_atlas_numerical_cluster_stage_b1.json"
            ),
            "results/canonical_orthogonal_k3_cover_dimension_theorem.md": _sha256(
                ROOT / "results/canonical_orthogonal_k3_cover_dimension_theorem.md"
            ),
            "results/reflection_maximal_coble_exact_nef_return.md": _sha256(
                ROOT / "results/reflection_maximal_coble_exact_nef_return.md"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-output", type=Path)
    arguments = parser.parse_args()
    certificate = build_certificate()
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
    if arguments.check_output is not None:
        if arguments.check_output.read_text(encoding="utf-8") != rendered:
            raise SystemExit("frozen Stage-B2b archive differs from exact replay")
    if not arguments.write and arguments.check_output is None:
        print(rendered, end="")


if __name__ == "__main__":
    main()
