"""Exact reduction and one-star traversal of the physical D_55 eigensection.

The marked chamber is reached from D_68 through physical D_86, crossing
D_86 wall (20, 0).  The verifier restricts all 41 primitive walls to the
fixed golden eigenspace, proves an irredundant 11-facet description over
Q(phi), and crosses every internal facet once.

Representative indices do not identify physical chambers.  Return and
collision checks below therefore compare complete transported wall sets.
Nothing here asserts that the global restricted chamber graph has closed.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Sequence

from sympy import Matrix

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.geometric.golden_halfspaces import (
    audit_golden_halfspace_certificate,
    positive_proportional_classes,
)
from coxeter_polytopes.geometric.k3 import K3_TRACE_SCALAR, field_inner_product
from scripts.analyze_shimada_disc11_d86 import D86_FACET_WITNESSES
from scripts.analyze_shimada_disc11_eigen_nef import (
    _crossing_point,
    _point_record,
    _restricted_coefficients,
)
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    D68_INTERNAL_FACET_WITNESSES,
    DISC11_SOURCE_SHA256,
    analyze_conjugated_roots,
    expected_shimada_gram,
    transported_identity_eigenbasis,
)


def _g(a: int, b: int = 0, denominator: int = 1) -> Golden:
    return Golden(Fraction(a, denominator), Fraction(b, denominator))


D55_STATE_SHA256 = "f9164f790643d027483fae91f1fccf21e7b04a372a0c7c38675ffa23edfbdf0a"
D55_PHYSICAL_WALL_ROWS_SHA256 = (
    "3aafb6ae4a58df9712b4edfd0bddbecc858cc859bb4a08ace14e674562819103"
)
D55_PHYSICAL_WALL_SET_SHA256 = (
    "90b630a04c3e6fa81ed2a9130abe387afe761b96965cc6345321e3727d9de245"
)
D55_RESTRICTED_ROWS_SHA256 = (
    "fa14c05fe1f69a18f4974cfb27a69dc1d14967408970e6eda311ff4c47ec8344"
)


D55_POSITIVE_PROPORTIONAL_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0),),
    ((0, 1), (6, 0)),
    ((1, 0), (3, 0)),
    ((1, 1), (8, 1)),
    ((2, 0),),
    ((3, 1),),
    ((4, 0), (5, 0)),
    ((4, 1), (7, 1)),
    ((5, 1),),
    ((7, 0), (8, 0)),
    ((9, 0), (10, 0)),
    ((11, 0),),
    ((11, 1), (21, 0)),
    ((12, 0),),
    ((12, 1),),
    ((13, 0),),
    ((13, 1),),
    ((14, 0),),
    ((15, 0),),
    ((15, 1),),
    ((16, 0),),
    ((16, 1),),
    ((17, 0),),
    ((18, 0),),
    ((18, 1),),
    ((19, 0),),
    ((20, 0),),
    ((21, 1),),
    ((22, 0),),
    ((22, 1),),
    ((23, 0),),
    ((23, 1),),
    ((24, 0),),
)


D55_ROOT_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 1), (6, 0)),
    ((1, 0), (3, 0)),
    ((1, 1), (8, 1)),
    ((4, 0), (5, 0)),
    ((4, 1), (7, 1)),
    ((7, 0), (8, 0)),
    ((9, 0), (10, 0)),
)
D55_INTERNAL_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((12, 1),),
    ((13, 1),),
    ((14, 0),),
    ((16, 1),),
)
D55_FACET_CLASSES = D55_ROOT_FACET_CLASSES + D55_INTERNAL_FACET_CLASSES
D55_FACET_REPRESENTATIVES = tuple(item[0] for item in D55_FACET_CLASSES)


# Exact timelike points on the relative interiors of the eleven faces.
D55_FACET_WITNESSES: tuple[
    tuple[tuple[tuple[int, int], ...], tuple[Golden, ...]], ...
] = (
    (((0, 1), (6, 0)), (
        _g(-69, 35), _g(53, -35), _g(-7, 3, 3), _g(1, 0, 3),
        _g(0), _g(-11, 3, 3), _g(-19, 3, 3), _g(-8, 1),
        _g(-35, 6, 3), _g(-8),
    )),
    (((1, 0), (3, 0)), (
        _g(-49, 23), _g(37, -25), _g(-7, 3, 3), _g(1, 0, 3),
        _g(0), _g(-5, 0, 3), _g(-13, 0, 3), _g(-6),
        _g(-29, 3, 3), _g(-6, -1),
    )),
    (((1, 1), (8, 1)), (
        _g(52, -38), _g(-6, 3), _g(-8, 4, 3), _g(2, -1, 3),
        _g(0), _g(-13, 5, 3), _g(-29, 10, 3), _g(-12, 4),
        _g(-55, 20, 3), _g(-15, 5),
    )),
    (((4, 0), (5, 0)), (
        _g(-149, 87), _g(119, -75), _g(-3, 1), _g(-1), _g(-2),
        _g(-6, 2), _g(-11, 4), _g(-13, 5), _g(-19, 8), _g(-16, 6),
    )),
    (((4, 1), (7, 1)), (
        _g(-89, 49), _g(70, -45), _g(-4, 1), _g(-3), _g(-5),
        _g(-7, 1), _g(-10, 2), _g(-10, 2), _g(-14, 4), _g(-12, 2),
    )),
    (((7, 0), (8, 0)), (
        _g(33, -26), _g(-1), _g(0, -1), _g(3, -3), _g(3, -4),
        _g(-1, -2), _g(-6), _g(-8, 1), _g(-14, 4), _g(-10, 1),
    )),
    (((9, 0), (10, 0)), (
        _g(-121, 67), _g(95, -61), _g(-7, 3, 3), _g(1, 0, 3),
        _g(0), _g(-11, 3, 3), _g(-28, 9, 3), _g(-11, 3),
        _g(-50, 15, 3), _g(-13, 3),
    )),
    (((12, 1),), (
        _g(-57, 29), _g(44, -29), _g(3, -3), _g(6, -5),
        _g(11, -9), _g(7, -7), _g(2, -5), _g(0, -4),
        _g(-6, -1), _g(-2, -4),
    )),
    (((13, 1),), (
        _g(-61, -33, 5), _g(19, -33, 5), _g(43, -71, 15),
        _g(44, -43, 15), _g(0), _g(-16, -28, 15),
        _g(-47, -56, 15), _g(-11, -33, 5), _g(-79, -112, 15),
        _g(-21, -33, 5),
    )),
    (((14, 0),), (
        _g(-124, 70), _g(97, -62), _g(-2), _g(1, -2), _g(1, -3),
        _g(-3, -1), _g(-8, 1), _g(-10, 2), _g(-16, 5), _g(-12, 2),
    )),
    (((16, 1),), (
        _g(11, -15), _g(-12, 5), _g(-2, 0, 3), _g(-4, 3, 3),
        _g(0), _g(-1, -3, 3), _g(-5, -6, 3), _g(0, -4),
        _g(-7, -12, 3), _g(-2, -4),
    )),
)


F01, F10, F11, F40, F41, F70, F90, F121, F131, F140, F161 = (
    D55_FACET_REPRESENTATIVES
)


# Fixed independent supports; coefficients are recovered and audited exactly.
D55_REDUNDANCY_SUPPORTS: tuple[
    tuple[tuple[int, int], tuple[tuple[int, int], ...]], ...
] = (
    ((0, 0), (F01, F10, F11, F40, F41, F90, F121)),
    ((2, 0), (F01, F10, F11, F40, F41, F70, F121, F131, F161)),
    ((3, 1), (F01, F10, F11, F40, F41, F90, F121, F161)),
    ((5, 1), (F121, F161)),
    ((11, 0), (F70, F131, F140, F161)),
    ((11, 1), (F70, F131, F161)),
    ((12, 0), (F70, F131, F140, F161)),
    ((13, 0), (F01, F10, F11, F40, F41, F70, F90, F121, F131, F161)),
    ((15, 0), (F01, F10, F11, F40, F41, F70, F90, F121, F131, F161)),
    ((15, 1), (F11, F70, F121, F131, F140, F161)),
    ((16, 0), (F01, F10, F11, F40, F41, F70, F90, F121, F161)),
    ((17, 0), (F01, F10, F11, F40, F41, F90, F121)),
    ((18, 0), (F10, F11, F40, F41, F70, F121, F131, F140, F161)),
    ((18, 1), (F01, F10, F11, F40, F41, F70, F90, F121, F131, F161)),
    ((19, 0), (F01, F10, F11, F40, F41, F90, F121, F161)),
    ((20, 0), (F01, F10, F11, F40, F41, F70, F90, F121, F131, F161)),
    ((21, 1), (F10, F11, F40, F41, F70, F121, F131, F140, F161)),
    ((22, 0), (F01, F10, F11, F40, F41, F90, F121, F161)),
    ((22, 1), (F01, F10, F11, F40, F41, F70, F90, F121, F161)),
    ((23, 0), (F11, F70, F121, F131, F140, F161)),
    ((23, 1), (F01, F10, F11, F40, F41, F90, F121, F161)),
    ((24, 0), (F10, F11, F40, F41, F70, F121, F131, F140, F161)),
)


D55_INTERNAL_CROSSINGS = (
    ((12, 1), 20, 3, "5aecf9e9ef2e330b1d822471efa130f1cc1af31d9694cfed4cfdd15fac135611"),
    ((13, 1), 14, 2, "24cbab8eba94d9ecb0c569f85e22b3569f4739dff28803cdd8209ecd259058e9"),
    ((14, 0), 16, 3, "a26f1a20c5b38cb4fc00ff55da413e40a3ff719d8565d5f430eb25ac9b28f953"),
    ((16, 1), 86, 2, "3e8dccdc9aedb61f0e32c9b8ffb5f61e893986a545a0d6108814f362494eb3b9"),
)


def _canonical_sha(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _pair(row: Sequence[Golden], point: Sequence[Golden]) -> Golden:
    return sum((a * b for a, b in zip(row, point)), Golden(0))


def _physical_wall_set(atlas: ShimadaDisc11Atlas, state) -> frozenset[tuple[int, ...]]:
    return frozenset(wall for orbit in atlas.state_wall_orbits(state) for wall in orbit)


def _physical_wall_set_sha(atlas: ShimadaDisc11Atlas, state) -> str:
    return _canonical_sha(sorted(_physical_wall_set(atlas, state)))


def _solve_fixed_support(
    target: tuple[Golden, ...],
    generators: tuple[tuple[Golden, ...], ...],
) -> tuple[Golden, ...]:
    """Solve one archived independent conic support over Q(phi)."""

    variable_count = len(generators)
    work = [
        [generators[column][coordinate] for column in range(variable_count)]
        + [target[coordinate]]
        for coordinate in range(len(target))
    ]
    pivot_row = 0
    pivot_rows: dict[int, int] = {}
    for column in range(variable_count):
        pivot = next(
            (row for row in range(pivot_row, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            raise ValueError("an archived D_55 redundancy support is dependent")
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(len(work)):
            if row == pivot_row or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][index] - factor * work[pivot_row][index]
                for index in range(variable_count + 1)
            ]
        pivot_rows[column] = pivot_row
        pivot_row += 1
    if any(not any(row[:-1]) and row[-1] for row in work):
        raise ValueError("an archived D_55 redundancy support is inconsistent")
    solution = tuple(work[pivot_rows[column]][-1] for column in range(variable_count))
    reconstructed = tuple(
        sum(
            (solution[index] * generators[index][coordinate] for index in range(variable_count)),
            Golden(0),
        )
        for coordinate in range(len(target))
    )
    if reconstructed != target:
        raise ArithmeticError("exact D_55 redundancy reconstruction failed")
    return solution


def analyze_d55(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")

    atlas = ShimadaDisc11Atlas(path)
    basis = transported_identity_eigenbasis()
    d68 = atlas.representative_state(68)
    d86 = atlas.neighboring_state(d68, 11, 1)
    if d86 is None:
        raise ValueError("the physical D_86 state is unavailable")
    d55 = atlas.neighboring_state(d86, 20, 0)
    if d55 is None or d55.index != 55:
        raise ValueError("D_86 wall (20, 0) did not reach D_55")

    physical_orbits = atlas.state_wall_orbits(d55)
    physical_rows = [
        [[orbit_index, member_index], list(wall)]
        for orbit_index, orbit in enumerate(physical_orbits)
        for member_index, wall in enumerate(orbit)
    ]
    coefficients = {
        (orbit_index, member_index): _restricted_coefficients(wall, basis)
        for orbit_index, orbit in enumerate(physical_orbits)
        for member_index, wall in enumerate(orbit)
    }
    restricted_rows = [
        [[*label], [entry.to_json() for entry in row]]
        for label, row in sorted(coefficients.items())
    ]
    class_records = positive_proportional_classes(coefficients)
    classes = tuple(record.members for record in class_records)

    dependencies = {}
    dependence_records = []
    dependences_verified = True
    for target, support in D55_REDUNDANCY_SUPPORTS:
        solution = _solve_fixed_support(
            coefficients[target], tuple(coefficients[label] for label in support)
        )
        terms = tuple(zip(solution, support))
        verified = all(value.sign() > 0 for value in solution)
        dependences_verified &= verified
        dependencies[target] = terms
        dependence_records.append(
            {
                "target": list(target),
                "support": [list(label) for label in support],
                "coefficients": [value.to_json() for value in solution],
                "all_coefficients_positive": verified,
                "identity_verified": True,
            }
        )

    witness_records = []
    witnesses_verified = True
    for target_class, witness in D55_FACET_WITNESSES:
        target = set(target_class)
        pairings = {label: _pair(row, witness) for label, row in coefficients.items()}
        q_square = K3_TRACE_SCALAR * field_inner_product(field_gram(), witness, witness)
        target_zero = all(pairings[label] == 0 for label in target)
        others_positive = all(
            value.sign() > 0 for label, value in pairings.items() if label not in target
        )
        verified = target_zero and others_positive and q_square.sign() > 0
        witnesses_verified &= verified
        witness_records.append(
            {
                "wall_class": [list(label) for label in target_class],
                "q_square": q_square.to_json(),
                "target_zero": target_zero,
                "all_others_positive": others_positive,
                "minimum_other_pairing": min(
                    value for label, value in pairings.items() if label not in target
                ).to_json(),
            }
        )

    d86_face = next(
        witness for wall_class, witness in D86_FACET_WITNESSES
        if wall_class == ((20, 0),)
    )
    d86_interior = _crossing_point(D68_INTERNAL_FACET_WITNESSES[1][1], 2)
    d55_interior = tuple(
        d86_face[index] - Golden(1) / 16 * (d86_interior[index] - d86_face[index])
        for index in range(10)
    )
    generic_audit = audit_golden_halfspace_certificate(
        coefficients,
        facet_classes=D55_FACET_CLASSES,
        facet_witnesses={
            wall_class[0]: witness for wall_class, witness in D55_FACET_WITNESSES
        },
        redundancies=dependencies,
        interior=d55_interior,
        quadratic_gram=tuple(
            tuple(K3_TRACE_SCALAR * entry for entry in row) for row in field_gram()
        ),
    )

    chamber = atlas.chamber(55)
    gram = Matrix(expected_shimada_gram())
    inverse_gram = gram.inv()
    root_labels = []
    roots = []
    for orbit_index, adjacency in enumerate(chamber.adjacencies):
        if adjacency.kind != "_minustwo":
            continue
        for member_index, wall in enumerate(physical_orbits[orbit_index]):
            primal = Matrix(1, 20, wall) * inverse_gram
            if any(not entry.is_Integer for entry in primal):
                raise ValueError("a D_55 root wall recovered a nonintegral root")
            root = tuple(int(entry) for entry in primal)
            if (Matrix(1, 20, root) * gram * Matrix(20, 1, root))[0] != -2:
                raise ValueError("a D_55 root wall did not have square -2")
            root_labels.append((orbit_index, member_index))
            roots.append(root)
    root_filter = analyze_conjugated_roots(tuple(roots))
    surviving_root_labels = {
        root_labels[index] for index in root_filter["surviving_zero_based_indices"]
    }
    root_facet_labels = {label for wall_class in D55_ROOT_FACET_CLASSES for label in wall_class}

    witness_by_label = {
        wall_class[0]: witness
        for wall_class, witness in D55_FACET_WITNESSES
        if wall_class in D55_INTERNAL_FACET_CLASSES
    }
    known_states = {"D68": d68, "D86": d86, "D55": d55}
    crossing_records = []
    crossing_states = []
    crossings_verified = True
    for label, expected_target, denominator, expected_wall_hash in D55_INTERNAL_CROSSINGS:
        neighbor = atlas.neighboring_state(d55, *label)
        if neighbor is None:
            raise ValueError("a D_55 internal facet became a K3 boundary")
        face = witness_by_label[label]
        crossing_point = tuple(
            face[index] - Golden(1) / denominator * (d55_interior[index] - face[index])
            for index in range(10)
        )
        record = _point_record(atlas, neighbor, crossing_point, basis)
        wall_hash = _physical_wall_set_sha(atlas, neighbor)
        physical_matches = [
            name
            for name, state in known_states.items()
            if _physical_wall_set(atlas, neighbor) == _physical_wall_set(atlas, state)
        ]
        verified = (
            neighbor.index == expected_target
            and wall_hash == expected_wall_hash
            and record["q_square_positive"]
            and record["negative_wall_pairings"] == 0
            and record["zero_wall_pairings"] == 0
        )
        crossings_verified &= verified
        crossing_states.append(neighbor)
        crossing_records.append(
            {
                "D55_wall": list(label),
                "target_representative": neighbor.index,
                "target_physical_wall_set_sha256": wall_hash,
                "physical_matches_among_known_states": physical_matches,
                "crossing_denominator": denominator,
                "target_interior": record,
                "verified": verified,
            }
        )

    returned = crossing_states[-1]
    return_is_physical_d86 = (
        _physical_wall_set(atlas, returned) == _physical_wall_set(atlas, d86)
        and returned.row_isometry == d86.row_isometry
    )
    distinct_crossing_wall_sets = len(
        {_physical_wall_set(atlas, state) for state in crossing_states}
    ) == len(crossing_states)
    entry_wall = _restricted_coefficients(physical_orbits[16][1], basis)
    d86_entry = _restricted_coefficients(atlas.state_wall_orbits(d86)[20][0], basis)
    entry_wall_opposite = entry_wall == tuple(-value for value in d86_entry)

    facet_sets = {frozenset(item) for item in D55_FACET_CLASSES}
    computed_sets = {frozenset(item) for item in classes}
    redundant_sets = computed_sets - facet_sets
    redundant_representatives = {
        record.members[0]
        for record in class_records
        if frozenset(record.members) not in facet_sets
    }
    partition_verified = (
        classes == D55_POSITIVE_PROPORTIONAL_CLASSES
        and len(coefficients) == 41
        and len(classes) == 33
        and len(facet_sets) == 11
        and len(redundant_sets) == len(D55_REDUNDANCY_SUPPORTS) == 22
        and {target for target, _ in D55_REDUNDANCY_SUPPORTS}
        == redundant_representatives
    )
    fingerprints_verified = (
        _canonical_sha(d55.row_isometry) == D55_STATE_SHA256
        and _canonical_sha(physical_rows) == D55_PHYSICAL_WALL_ROWS_SHA256
        and _physical_wall_set_sha(atlas, d55) == D55_PHYSICAL_WALL_SET_SHA256
        and _canonical_sha(restricted_rows) == D55_RESTRICTED_ROWS_SHA256
    )
    root_filter_verified = (
        len(root_labels) == 18
        and root_filter["trace_minus_two_verified"] == 18
        and root_filter["field_integral_reflection_roots"] == 17
        and root_filter["field_span_rank"] == 10
        and set(root_labels) - surviving_root_labels == {(3, 1)}
        and root_facet_labels <= surviving_root_labels
    )
    verified = (
        partition_verified
        and fingerprints_verified
        and witnesses_verified
        and dependences_verified
        and generic_audit.passes
        and root_filter_verified
        and crossings_verified
        and distinct_crossing_wall_sets
        and return_is_physical_d86
        and entry_wall_opposite
    )

    return {
        "verified": verified,
        "source_sha256": digest,
        "physical_state": {
            "entered_from": {"representative": 86, "wall": [20, 0]},
            "representative": 55,
            "state_sha256": D55_STATE_SHA256,
            "physical_wall_rows_sha256": D55_PHYSICAL_WALL_ROWS_SHA256,
            "physical_wall_set_sha256": D55_PHYSICAL_WALL_SET_SHA256,
            "restricted_rows_sha256": D55_RESTRICTED_ROWS_SHA256,
            "fingerprints_verified": fingerprints_verified,
        },
        "restricted_section": {
            "ambient_primitive_walls": len(coefficients),
            "positive_proportional_classes": len(classes),
            "irredundant_facets": len(D55_FACET_CLASSES),
            "root_facets": len(D55_ROOT_FACET_CLASSES),
            "internal_facets": len(D55_INTERNAL_FACET_CLASSES),
            "root_facet_classes": [
                [list(label) for label in wall_class] for wall_class in D55_ROOT_FACET_CLASSES
            ],
            "internal_facet_classes": [
                [list(label) for label in wall_class]
                for wall_class in D55_INTERNAL_FACET_CLASSES
            ],
            "positive_proportionality_ratios": {
                f"{record.members[index][0]},{record.members[index][1]}":
                record.ratios_to_first[index].to_json()
                for record in class_records
                for index in range(1, len(record.members))
            },
            "facet_witnesses": witness_records,
            "redundancy_dependences": dependence_records,
            "partition_verified": partition_verified,
            "witnesses_verified": witnesses_verified,
            "dependences_verified": dependences_verified,
            "generic_exact_halfspace_kernel_verified": generic_audit.passes,
        },
        "root_filter": {
            "ambient_minus_two_walls": len(root_labels),
            "positive_field_roots": root_filter["positive_at_both_places"],
            "integral_field_roots": root_filter["field_integral_reflection_roots"],
            "shell_distribution": root_filter["field_integral_shell_distribution"],
            "field_span_rank": root_filter["field_span_rank"],
            "non_field_root_walls": [
                list(label) for label in root_labels if label not in surviving_root_labels
            ],
            "all_root_facet_members_are_field_roots": root_facet_labels <= surviving_root_labels,
            "verified": root_filter_verified,
        },
        "internal_first_star": {
            "crossings": crossing_records,
            "representative_targets": [record["target_representative"] for record in crossing_records],
            "all_crossings_verified": crossings_verified,
            "crossing_physical_wall_sets_pairwise_distinct": distinct_crossing_wall_sets,
            "wall_16_member_1_is_opposite_D86_entry_wall": entry_wall_opposite,
            "wall_16_member_1_returns_to_same_physical_D86": return_is_physical_d86,
        },
        "scope": "complete_local_D55_section_and_internal_first_star_not_global_eigen_nef_completeness",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze_d55(args.disc11), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
