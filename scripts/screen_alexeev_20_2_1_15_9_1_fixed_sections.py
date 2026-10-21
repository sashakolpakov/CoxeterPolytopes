"""Exact strict-facet/Coxeter screen for three Alexeev subgroup sections.

Input is the complete weighted-diagram inventory for ``(20,2,1)`` and
``(15,9,1)``.  Every proper subgroup class with fixed rank at least 11 is
screened.  The ambient integral lattice, its saturated subgroup-fixed
lattice, the rational restricted mirror configuration, and the geometric
realization of the diagram subgroup are kept as distinct objects.

A case is rejected only by a forbidden exact product between independently
strict-witnessed hyperbolic facets.  A zero-obstruction case is printed as a
``SURVIVOR`` and must pass a separate literal Vinberg/dual-cone volume gate.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
from scripts.certify_shimada_disc16_d5_full_subgroup_fishing import (
    _dependence_on_retained_facets,
    _pair_gate,
)
from scripts.certify_shimada_disc16_d5_h14_fixed_lattice import (
    _smith_decomposition,
)
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    CONFIGS,
    OUTPUT as INVENTORY,
    SCHEMA as INVENTORY_SCHEMA,
    SOURCE,
    SOURCE_SHA256,
    _ambient_reconstruction,
    _lift_actions,
)
from scripts.screen_shimada_disc16_atlas_archived_high_sections import (
    _audit_witness,
    _primitive,
    _strict_timelike_witness,
)
from scripts.screen_shimada_disc16_fresh_weighted_high_sections import (
    _exact_face_ray_facet_witness,
    _reflection_integrality,
    _simplicial_exact_facet_witness,
)
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _exact_nonnegative_dependence,
    _matrix_payload,
    _orbits,
    _rational_text,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_fixed_section_census.json"
)
SCHEMA = "coxeter.alexeev_20_2_1_15_9_1_fixed_section_census.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _matrix(payload: list[list[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _saturated_fixed_lattice(
    subgroup_actions: tuple[sp.Matrix, ...],
    ambient_gram: sp.Matrix,
    expected_rank: int,
) -> dict[str, Any]:
    rank = ambient_gram.rows
    nonidentity = tuple(action for action in subgroup_actions if action != sp.eye(rank))
    if not nonidentity:
        raise ArithmeticError("a proper subgroup has no nonidentity action")
    constraints = sp.Matrix.vstack(*(action.T - sp.eye(rank) for action in nonidentity))
    _, independent_rows = constraints.T.rref()
    reduced = constraints[list(independent_rows), :]
    if reduced.rank() != rank - expected_rank:
        raise ArithmeticError("fixed-space rank differs from the inventory")
    denominator = int(sp.ilcm(*(sp.Rational(value).q for value in reduced)))
    integral = denominator * reduced
    if any(not value.is_Integer for value in integral):
        raise ArithmeticError("fixed constraints did not clear over Z")
    diagonal, _, right = _smith_decomposition(integral)
    zero_columns = tuple(
        column
        for column in range(rank)
        if all(diagonal[row, column] == 0 for row in range(diagonal.rows))
    )
    if len(zero_columns) != expected_rank:
        raise ArithmeticError("Smith kernel has the wrong fixed rank")
    smith_basis = right[:, list(zero_columns)].T
    basis = smith_basis.lll(delta=sp.Rational(3, 4))
    _, pivot_columns = smith_basis.rref()
    change = basis[:, list(pivot_columns)] * smith_basis[:, list(pivot_columns)].inv()
    if not (
        basis.rank() == expected_rank
        and basis * constraints.T == sp.zeros(expected_rank, constraints.rows)
        and all(basis * action == basis for action in subgroup_actions)
        and change * smith_basis == basis
        and abs(int(change.det())) == 1
    ):
        raise ArithmeticError("saturated fixed lattice failed exact replay")
    gram = basis * ambient_gram * basis.T
    if exact_inertia(gram) != (1, expected_rank - 1, 0):
        raise ArithmeticError("fixed lattice is not Lorentzian")
    smith = smith_normal_form(gram, domain=ZZ)
    return {
        "basis": basis,
        "gram": gram,
        "constraint_rank": reduced.rows,
        "constraint_clearing_denominator": denominator,
        "smith_constraint_diagonal": [
            int(diagonal[index, index])
            for index in range(min(diagonal.rows, diagonal.cols))
        ],
        "determinant": int(gram.det()),
        "smith_diagonal_absolute": [
            abs(int(smith[index, index])) for index in range(expected_rank)
        ],
    }


def _fixed_coordinates(vector: sp.Matrix, fixed_basis: sp.Matrix) -> sp.Matrix:
    solution = sp.linsolve((fixed_basis.T, vector.T))
    if solution is sp.EmptySet:
        raise ArithmeticError("a fixed ambient vector is absent from the fixed basis span")
    coordinates = sp.Matrix(tuple(next(iter(solution))))
    if coordinates.T * fixed_basis != vector:
        raise ArithmeticError("fixed coordinates failed replay")
    return coordinates


def _screen_case(
    *,
    key: tuple[int, int, int],
    example: dict[str, Any],
    class_record: dict[str, Any],
    wall_gram: sp.Matrix,
    source_interior: sp.Matrix | None = None,
    minimum_fixed_rank: int = 11,
) -> dict[str, Any]:
    rank = int(example["ambient"]["rank"])
    ambient_gram, wall_rows = _ambient_reconstruction(wall_gram, rank)
    group = tuple(
        tuple(map(int, element))
        for element in example["weighted_diagram_group"]["elements"]
    )
    actions = _lift_actions(group, wall_rows, ambient_gram)
    subgroup: Subgroup = frozenset(
        tuple(map(int, element)) for element in class_record["representative_elements"]
    )
    if not (
        len(subgroup) == int(class_record["order"])
        and subgroup.issubset(actions)
        and int(class_record["fixed_rank"]) >= minimum_fixed_rank
    ):
        raise ArithmeticError("selected fixed-rank subgroup changed")
    subgroup_actions = tuple(actions[element] for element in sorted(subgroup))
    fixed = _saturated_fixed_lattice(
        subgroup_actions, ambient_gram, int(class_record["fixed_rank"])
    )
    fixed_basis: sp.Matrix = fixed["basis"]
    fixed_gram: sp.Matrix = fixed["gram"]
    fixed_rank = fixed_basis.rows

    # The dual vector pairing one with every leading basis root is a strict
    # timelike interior point for the two originally pinned Alexeev chambers.
    # Other complete Alexeev chambers may use a leading root basis whose dual
    # sum is spacelike, even though (for example) the sum of all chamber roots
    # is a strict timelike interior point.  Permit that independently audited
    # row vector to be supplied without changing the original default.
    source_interior_was_supplied = source_interior is not None
    if source_interior is None:
        source_interior = sp.ones(1, rank) * ambient_gram.inv()
    if source_interior.shape != (1, rank):
        raise ArithmeticError("source interior has the wrong ambient shape")
    source_pairings = tuple(
        sp.Rational(value) for value in source_interior * ambient_gram * wall_rows.T
    )
    if not (
        (source_interior * ambient_gram * source_interior.T)[0] > 0
        and all(value > 0 for value in source_pairings)
    ):
        raise ArithmeticError("canonical Alexeev chamber interior is not strict timelike")
    averaged_interior = sum(
        (source_interior * action for action in subgroup_actions),
        sp.zeros(1, rank),
    ) / len(subgroup_actions)
    if not (
        all(averaged_interior * action == averaged_interior for action in subgroup_actions)
        and (averaged_interior * ambient_gram * averaged_interior.T)[0] > 0
        and all(
            value > 0
            for value in averaged_interior * ambient_gram * wall_rows.T
        )
    ):
        raise ArithmeticError("averaged subgroup-fixed chamber interior failed")
    ample = _fixed_coordinates(averaged_interior, fixed_basis)

    wall_orbits = _orbits(subgroup, wall_rows.rows)
    orbit_normals = tuple(
        sum((wall_rows.row(index) for index in orbit), sp.zeros(1, rank))
        for orbit in wall_orbits
    )
    if any(
        normal * action != normal
        for normal in orbit_normals
        for action in subgroup_actions
    ):
        raise ArithmeticError("a folded wall sum is not subgroup fixed")
    stacked = sp.Matrix.vstack(*orbit_normals)
    pairing_matrix = sp.Matrix.hstack(
        *(fixed_basis * ambient_gram * normal.T for normal in orbit_normals)
    )
    if pairing_matrix.rank() != fixed_rank:
        raise ArithmeticError("folded wall inequalities do not span the fixed dual")
    inverse_fixed_gram = fixed_gram.inv()
    orbit_gram = pairing_matrix.T * inverse_fixed_gram * pairing_matrix
    ambient_orbit_gram = stacked * ambient_gram * stacked.T
    if not (
        orbit_gram == ambient_orbit_gram
        and exact_inertia(orbit_gram)
        == (1, fixed_rank - 1, len(wall_orbits) - fixed_rank)
        and all(value > 0 for value in ample.T * pairing_matrix)
    ):
        raise ArithmeticError("ambient and fixed folded-mirror geometry disagree")

    coordinates = pairing_matrix.T
    conic_facets = []
    conic_redundant = []
    for target in range(len(wall_orbits)):
        destination = (
            conic_redundant
            if _exact_nonnegative_dependence(coordinates, target) is not None
            else conic_facets
        )
        destination.append(target)
    redundancy_certificates = []
    for target in conic_redundant:
        dependence = _dependence_on_retained_facets(
            coordinates, target, tuple(conic_facets)
        )
        redundancy_certificates.append(
            {
                "orbit_index": target + 1,
                "terms": [
                    {
                        "facet_orbit_index": support + 1,
                        "coefficient": _rational_text(coefficient),
                    }
                    for support, coefficient in dependence
                ],
            }
        )

    active_facets: list[int] = []
    causal_inactive: list[int] = []
    unresolved: list[int] = []
    witnesses = []

    def accept_witness(target: int, result: tuple[sp.Matrix, dict[str, Any]]) -> None:
        witness, record = result
        square, pairings = _audit_witness(
            witness, target, fixed_gram, pairing_matrix
        )
        active_facets.append(target)
        witnesses.append(
            record
            | {
                "orbit_index": target + 1,
                "fixed_lattice_coordinates": [
                    _rational_text(sp.Rational(value)) for value in witness
                ],
                "witness_square": _rational_text(square),
                "minimum_strict_pairing": _rational_text(
                    min(
                        value
                        for index, value in enumerate(pairings)
                        if index != target
                    )
                ),
            }
        )

    for target in conic_facets:
        if orbit_gram[target, target] >= 0:
            causal_inactive.append(target)
            continue
        result = _strict_timelike_witness(
            target, fixed_gram, pairing_matrix, orbit_gram, ample
        )
        if result is None:
            result = _simplicial_exact_facet_witness(
                target, fixed_gram, pairing_matrix
            )
        if result is None:
            unresolved.append(target)
        else:
            accept_witness(target, result)

    preliminary = _pair_gate(orbit_gram.extract(active_facets, active_facets))
    if unresolved and preliminary["non_coxeter_pair_count"] == 0:
        still_unresolved = []
        for target in unresolved:
            result = _exact_face_ray_facet_witness(
                target, fixed_gram, pairing_matrix
            )
            if result is None:
                still_unresolved.append(target)
            else:
                accept_witness(target, result)
        unresolved = still_unresolved

    active_facets.sort()
    witnesses.sort(key=lambda row: row["orbit_index"])
    facet_gram = orbit_gram.extract(active_facets, active_facets)
    pair_gate = _pair_gate(facet_gram)
    if pair_gate["non_coxeter_pair_count"]:
        classification = "not_coxeter"
    elif unresolved:
        classification = "unresolved_facet_causality"
    else:
        classification = "coxeter_candidate_requires_volume_audit"

    mirror_records = []
    reflection_flags = {}
    for orbit_index, (orbit, normal) in enumerate(
        zip(wall_orbits, orbit_normals, strict=True)
    ):
        fixed_normal = inverse_fixed_gram * pairing_matrix[:, orbit_index]
        primitive, scale = _primitive(fixed_normal)
        norm = sp.Rational(orbit_gram[orbit_index, orbit_index])
        integral_reflection = None
        if norm < 0:
            integral_reflection, _ = _reflection_integrality(primitive, fixed_gram)
            reflection_flags[orbit_index] = integral_reflection
        mirror_records.append(
            {
                "orbit_index": orbit_index + 1,
                "source_wall_indices": [index + 1 for index in orbit],
                "ambient_orbit_sum_primal_coordinates": [int(value) for value in normal],
                "fixed_pairing_column": [
                    _rational_text(sp.Rational(value))
                    for value in pairing_matrix[:, orbit_index]
                ],
                "rational_fixed_normal_coordinates": [
                    _rational_text(sp.Rational(value)) for value in fixed_normal
                ],
                "primitive_integral_same_direction_coordinates": [
                    int(value) for value in primitive
                ],
                "presentation_scale_from_primitive": _rational_text(scale),
                "square": _rational_text(norm),
                "reflection_integral_on_saturated_fixed_lattice": integral_reflection,
                "conic_status": "facet" if orbit_index in conic_facets else "redundant",
                "hyperbolic_status": (
                    "active_facet"
                    if orbit_index in active_facets
                    else "causally_inactive_nonnegative_normal_square"
                    if orbit_index in causal_inactive
                    else "unresolved"
                    if orbit_index in unresolved
                    else "conically_redundant"
                ),
            }
        )

    case_id = f"A{key[0]}_{key[1]}_{key[2]}_C{int(class_record['class_index']):02d}"
    return {
        "id": case_id,
        "ambient_key": list(key),
        "case": {
            "class_index": int(class_record["class_index"]),
            "subgroup_order": len(subgroup),
            "subgroup_conjugacy_class_size": int(class_record["conjugacy_class_size"]),
            "cyclic": bool(class_record["cyclic"]),
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "source_wall_count": wall_rows.rows,
            "wall_orbit_count": len(wall_orbits),
            "wall_orbit_size_histogram": {
                str(size): count
                for size, count in sorted(Counter(map(len, wall_orbits)).items())
            },
            "representative_wall_permutations": [
                list(element) for element in sorted(subgroup)
            ],
            "representative_nonidentity_primal_actions": [
                _matrix_payload(actions[element])
                for element in sorted(subgroup)
                if actions[element] != sp.eye(rank)
            ],
        },
        "geometry_scope": {
            "diagram_subgroup_used": True,
            "geometric_automorphism_membership_certified": False,
            "ambient_lattice_and_saturated_fixed_lattice_distinguished": True,
            "fixed_divisor_lattice_and_restricted_mirror_lattice_distinguished": True,
            "finite_volume_gate_performed": False,
        },
        "saturated_fixed_lattice": {
            "basis_primal_rows": _matrix_payload(fixed_basis),
            "gram": _matrix_payload(fixed_gram),
            "determinant": fixed["determinant"],
            "smith_diagonal_absolute": fixed["smith_diagonal_absolute"],
            "constraint_rank": fixed["constraint_rank"],
            "constraint_clearing_denominator": fixed[
                "constraint_clearing_denominator"
            ],
            "smith_constraint_diagonal": fixed["smith_constraint_diagonal"],
            "inertia": list(exact_inertia(fixed_gram)),
        },
        "strict_fixed_interior": {
            "source_method": (
                "caller-supplied exact strict timelike interior"
                if source_interior_was_supplied
                else "pairing one with every leading ambient root"
            ),
            "fixed_lattice_coordinates": [
                _rational_text(sp.Rational(value)) for value in ample
            ],
            "source_square": _rational_text(
                sp.Rational((source_interior * ambient_gram * source_interior.T)[0])
            ),
            "averaged_fixed_square": _rational_text(
                sp.Rational((ample.T * fixed_gram * ample)[0])
            ),
            "minimum_folded_pairing": _rational_text(
                min(sp.Rational(value) for value in ample.T * pairing_matrix)
            ),
        },
        "exact_reduction": {
            "conic_facet_count": len(conic_facets),
            "conic_redundant_count": len(conic_redundant),
            "conic_facet_orbit_indices": [index + 1 for index in conic_facets],
            "redundancy_certificates": redundancy_certificates,
            "active_hyperbolic_facet_count": len(active_facets),
            "active_facet_orbit_indices": [index + 1 for index in active_facets],
            "causally_inactive_orbit_indices": [index + 1 for index in causal_inactive],
            "unresolved_orbit_indices": [index + 1 for index in unresolved],
            "strict_facet_witnesses": witnesses,
        },
        "restricted_mirrors": {
            "records": mirror_records,
            "all_active_reflections_integral_on_saturated_fixed_lattice": all(
                reflection_flags[index] for index in active_facets
            ),
            "facet_gram": _matrix_payload(facet_gram),
            "positive_wall_gram": _matrix_payload(-facet_gram),
        },
        "coxeter_gate": pair_gate,
        "decision": {
            "classification": classification,
            "finite_volume_claimed": False,
        },
        "checks": {
            "ambient_and_fixed_mirror_grams_agree": True,
            "saturated_integral_fixed_lattice_replayed": True,
            "every_removed_orbit_is_an_exact_nonnegative_sum_of_final_conic_facets": True,
            "every_active_facet_has_an_exact_strict_timelike_witness": True,
            "every_active_facet_pair_checked_exactly": True,
        },
    }


def build(
    source: Path = SOURCE,
    inventory_path: Path = INVENTORY,
) -> dict[str, Any]:
    if _sha256(source) != SOURCE_SHA256:
        raise ArithmeticError("Alexeev ancillary source hash changed")
    inventory_hash = _sha256(inventory_path)
    inventory = json.loads(inventory_path.read_text())
    if not (
        inventory["schema"] == INVENTORY_SCHEMA
        and inventory["source"]["sha256"] == SOURCE_SHA256
    ):
        raise ArithmeticError("weighted-diagram inventory has incompatible scope")
    source_text = source.read_text(encoding="utf-8")
    cases = []
    for example in inventory["examples"]:
        key = tuple(map(int, example["key"]))
        if key not in CONFIGS:
            raise ArithmeticError("unexpected ambient key in inventory")
        wall_gram = _extract_matrix(source_text, key)
        for class_record in example["subgroup_classes"]:
            if not class_record["proper_fixed_section_above_h9"]:
                continue
            case = _screen_case(
                key=key,
                example=example,
                class_record=class_record,
                wall_gram=wall_gram,
            )
            cases.append(case)
            print(
                json.dumps(
                    {
                        "id": case["id"],
                        "dimension": case["case"]["hyperbolic_dimension"],
                        "facets": case["exact_reduction"][
                            "active_hyperbolic_facet_count"
                        ],
                        "bad_pairs": case["coxeter_gate"][
                            "non_coxeter_pair_count"
                        ],
                        "decision": case["decision"]["classification"],
                    }
                ),
                flush=True,
            )
            if case["decision"]["classification"] == "coxeter_candidate_requires_volume_audit":
                print(f"SURVIVOR {case['id']}", flush=True)

    classification_histogram = Counter(
        case["decision"]["classification"] for case in cases
    )
    dimension_histogram = Counter(case["case"]["hyperbolic_dimension"] for case in cases)
    survivors = [
        case["id"]
        for case in cases
        if case["decision"]["classification"]
        == "coxeter_candidate_requires_volume_audit"
    ]
    return {
        "schema": SCHEMA,
        "source": {
            "primary_path": str(source),
            "primary_sha256": SOURCE_SHA256,
            "inventory_path": str(inventory_path.relative_to(ROOT)),
            "inventory_sha256": inventory_hash,
        },
        "scope": {
            "ambient_keys": [list(key) for key in CONFIGS],
            "excluded_parallel_assignment": [22, 2, 0],
            "proper_subgroup_fixed_rank_cutoff": 11,
            "diagram_symmetry_not_conflated_with_geometric_automorphism": True,
            "finite_volume_test_deferred_until_exact_zero_angle_obstruction": True,
        },
        "summary": {
            "case_count": len(cases),
            "dimension_histogram": {
                str(key): value for key, value in sorted(dimension_histogram.items())
            },
            "classification_histogram": dict(sorted(classification_histogram.items())),
            "total_certified_non_coxeter_pair_count": sum(
                case["coxeter_gate"]["non_coxeter_pair_count"] for case in cases
            ),
            "survivor_count_requiring_volume_audit": len(survivors),
        },
        "survivors_requiring_volume_audit": survivors,
        "cases": cases,
    }


def main() -> None:
    payload = build()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
