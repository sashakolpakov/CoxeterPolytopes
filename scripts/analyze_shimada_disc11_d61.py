"""Exact reduction of the common physical D_61 golden K3 eigensection.

The state is reached both from physical D_30 wall (22, 0) and physical D_27
wall (23, 0).  Equality is checked on cumulative row isometries and complete
transported wall sets before the local section is reduced and crossed.
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
from scripts.analyze_shimada_disc11_d30 import (
    D30_FACET_WITNESSES,
    D30_INTERIOR,
    _physical_states,
)
from scripts.analyze_shimada_disc11_d55 import (
    _canonical_sha,
    _physical_wall_set,
    _physical_wall_set_sha,
    _solve_fixed_support,
)
from scripts.analyze_shimada_disc11_eigen_nef import (
    _point_record,
    _restricted_coefficients,
)
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    DISC11_SOURCE_SHA256,
    analyze_conjugated_roots,
    expected_shimada_gram,
    transported_identity_eigenbasis,
)


def _g(a: int, b: int = 0, denominator: int = 1) -> Golden:
    return Golden(Fraction(a, denominator), Fraction(b, denominator))


D61_STATE_SHA256 = "a428b0c6df63f0d40aa5708b2773d992e4b39eaa42ccb28958ea05ca3214df13"
D61_PHYSICAL_WALL_ROWS_SHA256 = (
    "8daa0ee3d09e3ef7563277b8a7932e151d322c252389f888fd0baf68083075b8"
)
D61_PHYSICAL_WALL_SET_SHA256 = (
    "a14712879a75637ed4330d2a9b1c9a8629c940195636f231ff952b1d3bdba785"
)
D61_RESTRICTED_ROWS_SHA256 = (
    "9f185c3af6a4c8982346d091e7f137d4320ab5fffb15282b4deb1bb39b32f2d8"
)


D61_POSITIVE_PROPORTIONAL_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0), (2, 0)),
    ((1, 0), (8, 0)),
    ((3, 0), (6, 0)),
    ((4, 0), (15, 0)),
    ((5, 0), (12, 0)),
    ((7, 0),),
    ((9, 0), (17, 0)),
    ((10, 0),),
    ((11, 0),),
    ((13, 0), (14, 0)),
    ((16, 0),),
    ((18, 0),),
    ((19, 0),),
    ((20, 0),),
    ((21, 0),),
    ((22, 0),),
    ((23, 0),),
    ((24, 0),),
    ((25, 0),),
    ((26, 0),),
    ((27, 0),),
    ((28, 0),),
    ((29, 0),),
    ((30, 0),),
    ((31, 0),),
)


D61_ROOT_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0), (2, 0)),
    ((1, 0), (8, 0)),
    ((3, 0), (6, 0)),
    ((4, 0), (15, 0)),
    ((5, 0), (12, 0)),
    ((9, 0), (17, 0)),
    ((13, 0), (14, 0)),
)
D61_INTERNAL_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((18, 0),),
    ((19, 0),),
    ((20, 0),),
    ((23, 0),),
)
D61_FACET_CLASSES = D61_ROOT_FACET_CLASSES + D61_INTERNAL_FACET_CLASSES
D61_FACET_REPRESENTATIVES = tuple(item[0] for item in D61_FACET_CLASSES)


D61_FACET_WITNESSES: tuple[
    tuple[tuple[tuple[int, int], ...], tuple[Golden, ...]], ...
] = (
    (((0, 0), (2, 0)), (
        _g(-14), _g(12, -13, 2), _g(-2, -7, 2), _g(-6, -5, 2),
        _g(-12, -3, 2), _g(-10, -5, 2), _g(-10, -7, 2),
        _g(-10, -7, 2), _g(-6, -11, 2), _g(-8, -11, 2),
    )),
    (((1, 0), (8, 0)), (
        _g(-6, -6), _g(4, -6), _g(-1, -4), _g(-1, -4),
        _g(-2, -4), _g(-3, -4), _g(-3, -5), _g(-5, -4),
        _g(-5, -5), _g(-4, -6),
    )),
    (((3, 0), (6, 0)), (
        _g(-109, 3, 5), _g(51, -52, 5), _g(52, -79, 15),
        _g(41, -47, 15), _g(0), _g(-4, -32, 15),
        _g(-23, -64, 15), _g(-14, -27, 5), _g(-31, -128, 15),
        _g(-14, -32, 5),
    )),
    (((4, 0), (15, 0)), (
        _g(-49, -37, 5), _g(1, -22, 5), _g(-43, -24, 15),
        _g(-14, -12, 15), _g(0), _g(-14, -27, 15),
        _g(-43, -54, 15), _g(-24, -22, 5), _g(-71, -108, 15),
        _g(-9, -37, 5),
    )),
    (((5, 0), (12, 0)), (
        _g(-21, -1), _g(9, -10), _g(9, -16, 3), _g(6, -8, 3),
        _g(0), _g(0, -8, 3), _g(-6, -13, 3), _g(-3, -6),
        _g(-9, -26, 3), _g(-3, -7),
    )),
    (((9, 0), (17, 0)), (
        _g(-109, 3, 5), _g(51, -52, 5), _g(19, -28, 5),
        _g(12, -14, 5), _g(0), _g(-3, -9, 5), _g(-11, -18, 5),
        _g(-14, -27, 5), _g(-12, -41, 5), _g(-14, -32, 5),
    )),
    (((13, 0), (14, 0)), (
        _g(-67, -31, 5), _g(13, -31, 5), _g(-9, -47, 15),
        _g(-12, -16, 15), _g(0), _g(-12, -31, 15),
        _g(-39, -62, 15), _g(-22, -26, 5), _g(-63, -124, 15),
        _g(-17, -36, 5),
    )),
    (((18, 0),), (
        _g(-79, -27, 5), _g(21, -37, 5), _g(17, -64, 15),
        _g(16, -32, 15), _g(0), _g(-14, -32, 15),
        _g(-43, -64, 15), _g(-24, -27, 5), _g(-71, -128, 15),
        _g(-19, -37, 5),
    )),
    (((19, 0),), (
        _g(-18, 2), _g(18, -17, 2), _g(12, -15, 2),
        _g(12, -15, 2), _g(10, -15, 2), _g(8, -15, 2),
        _g(4, -15, 2), _g(0, -13, 2), _g(0, -15, 2),
        _g(2, -17, 2),
    )),
    (((20, 0),), (
        _g(-18, 2), _g(10, -9), _g(5, -7), _g(5, -7), _g(4, -7),
        _g(3, -7), _g(1, -7), _g(-1, -6), _g(-1, -7), _g(0, -8),
    )),
    (((23, 0),), (
        _g(-9, -47, 5), _g(14, -28, 5), _g(4, -28, 5),
        _g(4, -28, 5), _g(-1, -28, 5), _g(-6, -28, 5),
        _g(-16, -28, 5), _g(-26, -23, 5), _g(-26, -28, 5),
        _g(-21, -33, 5),
    )),
)


F0, F1, F3, F4, F5, F9, F13, F18, F19, F20, F23 = (
    D61_FACET_REPRESENTATIVES
)
D61_REDUNDANCY_SUPPORTS: tuple[
    tuple[tuple[int, int], tuple[tuple[int, int], ...]], ...
] = (
    ((7, 0), (F18, F19, F20, F23)),
    ((10, 0), (F19, F23)),
    ((11, 0), (F20, F23)),
    ((16, 0), (F0, F1, F3, F4, F5, F9, F13, F18)),
    ((21, 0), (F0, F3, F4, F9, F13, F18, F19, F20, F23)),
    ((22, 0), (F0, F1, F3, F4, F5, F9, F13, F18, F19, F20)),
    ((24, 0), (F0, F1, F3, F4, F5, F9, F13, F18)),
    ((25, 0), (F0, F1, F3, F4, F5, F9, F13, F18, F19, F20)),
    ((26, 0), (F0, F3, F4, F9, F13, F18, F19, F20, F23)),
    ((27, 0), (F0, F1, F3, F4, F5, F9, F13, F18)),
    ((28, 0), (F0, F1, F3, F4, F5, F9, F13, F18, F19, F20)),
    ((29, 0), (F0, F1, F3, F4, F5, F9, F13, F18)),
    ((30, 0), (F0, F1, F3, F4, F5, F9, F13, F18)),
    ((31, 0), (F19, F20, F23)),
)


D61_INTERNAL_CROSSINGS = (
    ((18, 0), 84, 3, "fbcb407412e6020202a0128c4b36bced36ab46fb49e0b767d057adc8a63f1c43"),
    ((19, 0), 30, 2, "b6b09e5748d52040e223d32b2ffcc52a4435c94b488abfb0b72e9d1b8f6845be"),
    ((20, 0), 27, 3, "8272b57f633c4a6d3a0474cf53a5c5e279b897d7a1d3d9e513d988bb88935735"),
    ((23, 0), 37, 3, "0998d5e9845d9cce86a9b8d8c41a828aa50c42a01efb9b8bc57a324b68fe536c"),
)


def _pair(row: Sequence[Golden], point: Sequence[Golden]) -> Golden:
    return sum((a * b for a, b in zip(row, point)), Golden(0))


def analyze_d61(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")

    atlas = ShimadaDisc11Atlas(path)
    basis = transported_identity_eigenbasis()
    states = _physical_states(atlas)
    d30 = states["d30_via_d181"]
    d86 = states["d86"]
    d27 = atlas.neighboring_state(d86, 24, 0)
    if d27 is None:
        raise ValueError("the comparison D_27 state is unavailable")
    d61_from_d30 = atlas.neighboring_state(d30, 22, 0)
    d61_from_d27 = atlas.neighboring_state(d27, 23, 0)
    if d61_from_d30 is None or d61_from_d27 is None:
        raise ValueError("a D_61 entry path is unavailable")
    path_state_equality = d61_from_d30.row_isometry == d61_from_d27.row_isometry
    path_wall_set_equality = (
        _physical_wall_set(atlas, d61_from_d30)
        == _physical_wall_set(atlas, d61_from_d27)
    )
    d61 = d61_from_d30
    if d61.index != 61:
        raise ValueError("the common physical state is not D_61")

    physical_orbits = atlas.state_wall_orbits(d61)
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
    for target, support in D61_REDUNDANCY_SUPPORTS:
        solution = _solve_fixed_support(
            coefficients[target], tuple(coefficients[label] for label in support)
        )
        verified = all(value.sign() > 0 for value in solution)
        dependences_verified &= verified
        dependencies[target] = tuple(zip(solution, support))
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
    for target_class, witness in D61_FACET_WITNESSES:
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

    d30_face = next(
        witness for wall_class, witness in D30_FACET_WITNESSES
        if wall_class == ((22, 0),)
    )
    d61_interior = tuple(
        d30_face[index] - (D30_INTERIOR[index] - d30_face[index]) / 3
        for index in range(10)
    )
    generic_audit = audit_golden_halfspace_certificate(
        coefficients,
        facet_classes=D61_FACET_CLASSES,
        facet_witnesses={
            wall_class[0]: witness for wall_class, witness in D61_FACET_WITNESSES
        },
        redundancies=dependencies,
        interior=d61_interior,
        quadratic_gram=tuple(
            tuple(K3_TRACE_SCALAR * entry for entry in row) for row in field_gram()
        ),
    )

    chamber = atlas.chamber(61)
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
                raise ValueError("a D_61 root wall recovered a nonintegral root")
            root = tuple(int(entry) for entry in primal)
            if (Matrix(1, 20, root) * gram * Matrix(20, 1, root))[0] != -2:
                raise ValueError("a D_61 root wall did not have square -2")
            root_labels.append((orbit_index, member_index))
            roots.append(root)
    root_filter = analyze_conjugated_roots(tuple(roots))
    surviving_root_labels = {
        root_labels[index] for index in root_filter["surviving_zero_based_indices"]
    }
    root_facet_labels = {label for wall_class in D61_ROOT_FACET_CLASSES for label in wall_class}

    known_states = {
        "D68": states["d68"],
        "D86": d86,
        "D181": states["d181"],
        "D30": d30,
        "D27": d27,
        "D61": d61,
    }
    witness_by_label = {
        wall_class[0]: witness
        for wall_class, witness in D61_FACET_WITNESSES
        if wall_class in D61_INTERNAL_FACET_CLASSES
    }
    crossing_records = []
    crossing_states = []
    crossings_verified = True
    for label, expected_target, denominator, expected_wall_hash in D61_INTERNAL_CROSSINGS:
        neighbor = atlas.neighboring_state(d61, *label)
        if neighbor is None:
            raise ValueError("a D_61 internal facet became a K3 boundary")
        face = witness_by_label[label]
        crossing_point = tuple(
            face[index] - (d61_interior[index] - face[index]) / denominator
            for index in range(10)
        )
        point_record = _point_record(atlas, neighbor, crossing_point, basis)
        wall_hash = _physical_wall_set_sha(atlas, neighbor)
        physical_matches = [
            name
            for name, state in known_states.items()
            if _physical_wall_set(atlas, neighbor) == _physical_wall_set(atlas, state)
        ]
        verified = (
            neighbor.index == expected_target
            and wall_hash == expected_wall_hash
            and point_record["q_square_positive"]
            and point_record["negative_wall_pairings"] == 0
            and point_record["zero_wall_pairings"] == 0
        )
        crossings_verified &= verified
        crossing_states.append(neighbor)
        crossing_records.append(
            {
                "D61_wall": list(label),
                "target_representative": neighbor.index,
                "target_physical_wall_set_sha256": wall_hash,
                "physical_matches_among_known_states": physical_matches,
                "crossing_denominator": denominator,
                "target_interior": point_record,
                "verified": verified,
            }
        )

    return_to_d30 = (
        crossing_states[1].row_isometry == d30.row_isometry
        and _physical_wall_set(atlas, crossing_states[1]) == _physical_wall_set(atlas, d30)
    )
    return_to_d27 = (
        crossing_states[2].row_isometry == d27.row_isometry
        and _physical_wall_set(atlas, crossing_states[2]) == _physical_wall_set(atlas, d27)
    )
    d30_entry = _restricted_coefficients(atlas.state_wall_orbits(d30)[22][0], basis)
    d27_entry = _restricted_coefficients(atlas.state_wall_orbits(d27)[23][0], basis)
    return_walls_opposite = (
        coefficients[(19, 0)] == tuple(-value for value in d30_entry)
        and coefficients[(20, 0)] == tuple(-value for value in d27_entry)
    )

    facet_sets = {frozenset(item) for item in D61_FACET_CLASSES}
    redundant_representatives = {
        record.members[0]
        for record in class_records
        if frozenset(record.members) not in facet_sets
    }
    partition_verified = (
        classes == D61_POSITIVE_PROPORTIONAL_CLASSES
        and len(coefficients) == 32
        and len(classes) == 25
        and len(facet_sets) == 11
        and len(redundant_representatives) == len(D61_REDUNDANCY_SUPPORTS) == 14
        and {target for target, _ in D61_REDUNDANCY_SUPPORTS}
        == redundant_representatives
    )
    fingerprints_verified = (
        _canonical_sha(d61.row_isometry) == D61_STATE_SHA256
        and _canonical_sha(physical_rows) == D61_PHYSICAL_WALL_ROWS_SHA256
        and _physical_wall_set_sha(atlas, d61) == D61_PHYSICAL_WALL_SET_SHA256
        and _canonical_sha(restricted_rows) == D61_RESTRICTED_ROWS_SHA256
    )
    root_filter_verified = (
        len(root_labels) == 18
        and root_filter["trace_minus_two_verified"] == 18
        and root_filter["field_integral_reflection_roots"] == 17
        and root_filter["field_span_rank"] == 10
        and set(root_labels) - surviving_root_labels == {(16, 0)}
        and root_facet_labels <= surviving_root_labels
    )
    verified = (
        path_state_equality
        and path_wall_set_equality
        and partition_verified
        and fingerprints_verified
        and witnesses_verified
        and dependences_verified
        and generic_audit.passes
        and root_filter_verified
        and crossings_verified
        and return_to_d30
        and return_to_d27
        and return_walls_opposite
    )

    return {
        "verified": verified,
        "source_sha256": digest,
        "physical_state": {
            "representative": 61,
            "D30_wall_22_and_D27_wall_23_same_row_isometry": path_state_equality,
            "D30_wall_22_and_D27_wall_23_same_physical_wall_set": path_wall_set_equality,
            "state_sha256": D61_STATE_SHA256,
            "physical_wall_rows_sha256": D61_PHYSICAL_WALL_ROWS_SHA256,
            "physical_wall_set_sha256": D61_PHYSICAL_WALL_SET_SHA256,
            "restricted_rows_sha256": D61_RESTRICTED_ROWS_SHA256,
            "fingerprints_verified": fingerprints_verified,
        },
        "restricted_section": {
            "ambient_primitive_walls": len(coefficients),
            "positive_proportional_classes": len(classes),
            "irredundant_facets": len(D61_FACET_CLASSES),
            "root_facets": len(D61_ROOT_FACET_CLASSES),
            "internal_facets": len(D61_INTERNAL_FACET_CLASSES),
            "root_facet_classes": [
                [list(label) for label in wall_class] for wall_class in D61_ROOT_FACET_CLASSES
            ],
            "internal_facet_classes": [
                [list(label) for label in wall_class]
                for wall_class in D61_INTERNAL_FACET_CLASSES
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
            "returns_to_same_physical_D30": return_to_d30,
            "returns_to_same_physical_D27": return_to_d27,
            "return_walls_are_opposite_entry_walls": return_walls_opposite,
            "new_physical_neighbors": [84, 37],
            "verified": crossings_verified,
        },
        "scope": "complete_local_D61_section_and_internal_first_star_not_global_eigen_nef_completeness",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze_d61(args.disc11), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
