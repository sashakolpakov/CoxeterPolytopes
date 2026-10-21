"""Exact reduction and internal traversal of the physical D_27 eigensection.

The marked state is reached from D_68 through D_86 by crossing D_86 wall
(24, 0).  One restricted internal facet is supported by two ambient walls;
its crossing is certified by two four-step Conway paths with the same exact
physical endpoint.  All chamber comparisons use transported wall sets.
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
from scripts.analyze_shimada_disc11_d55 import (
    _canonical_sha,
    _physical_wall_set,
    _physical_wall_set_sha,
    _solve_fixed_support,
)
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


D27_STATE_SHA256 = "a58b37537a1856b35d640442a3ea69d8dc6fa998190c236a1c066c7263a8e783"
D27_PHYSICAL_WALL_ROWS_SHA256 = (
    "d4ecdae87f5bdf7e8c59de6d68fc88e605c3ea34e20c8bbe75abf8221c4a79e2"
)
D27_PHYSICAL_WALL_SET_SHA256 = (
    "8272b57f633c4a6d3a0474cf53a5c5e279b897d7a1d3d9e513d988bb88935735"
)
D27_RESTRICTED_ROWS_SHA256 = (
    "6445047372ce254d4e1b5cff128ab540cf5c884b036e39ce157eb178860bc568"
)
D27_OPPOSITE_STATE_SHA256 = (
    "5864d28f2ba531a1c60bc7d9e66864c80df7b978fefee5470939197f363e9000"
)
D27_OPPOSITE_PHYSICAL_WALL_SET_SHA256 = (
    "2008055693fd68ed2a0c5ac7778e18edda27c83be5d71c56dbea7b92168b663b"
)


D27_POSITIVE_PROPORTIONAL_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0),),
    ((1, 0), (6, 0)),
    ((2, 0), (9, 0)),
    ((3, 0),),
    ((4, 0), (17, 0)),
    ((5, 0), (12, 0)),
    ((7, 0), (16, 0)),
    ((8, 0), (10, 0)),
    ((11, 0), (15, 0)),
    ((13, 0), (18, 0)),
    ((14, 0),),
    ((19, 0), (31, 0)),
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
    ((32, 0),),
)


D27_ROOT_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((1, 0), (6, 0)),
    ((2, 0), (9, 0)),
    ((4, 0), (17, 0)),
    ((5, 0), (12, 0)),
    ((7, 0), (16, 0)),
    ((8, 0), (10, 0)),
    ((11, 0), (15, 0)),
    ((13, 0), (18, 0)),
)
D27_INTERNAL_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((19, 0), (31, 0)),
    ((20, 0),),
    ((21, 0),),
    ((23, 0),),
)
D27_FACET_CLASSES = D27_ROOT_FACET_CLASSES + D27_INTERNAL_FACET_CLASSES
D27_FACET_REPRESENTATIVES = tuple(item[0] for item in D27_FACET_CLASSES)


D27_FACET_WITNESSES: tuple[
    tuple[tuple[tuple[int, int], ...], tuple[Golden, ...]], ...
] = (
    (((1, 0), (6, 0)), (
        _g(-5, -3), _g(-1), _g(-3), _g(0, -2), _g(2, -4),
        _g(1, -4), _g(1, -5), _g(-1, -4), _g(-2, -4), _g(2, -7),
    )),
    (((2, 0), (9, 0)), (
        _g(15, -15), _g(-1), _g(-4, 1), _g(-1, -1), _g(0, -2),
        _g(-1, -2), _g(-3, -2), _g(-5, -1), _g(-8), _g(-5, -2),
    )),
    (((4, 0), (17, 0)), (
        _g(-2, -5), _g(-1), _g(-1, -1), _g(2, -3), _g(6, -6),
        _g(5, -6), _g(3, -6), _g(1, -5), _g(-2, -4), _g(2, -7),
    )),
    (((5, 0), (12, 0)), (
        _g(19, -17), _g(-1), _g(-3), _g(0, -2), _g(2, -4),
        _g(0, -3), _g(-3, -2), _g(-5, -1), _g(-9, 1), _g(-5, -2),
    )),
    (((7, 0), (16, 0)), (
        _g(11, -13), _g(-1), _g(-3), _g(0, -2), _g(2, -4),
        _g(1, -4), _g(-2, -3), _g(-4, -2), _g(-7, -1), _g(-3, -4),
    )),
    (((8, 0), (10, 0)), (
        _g(17, -16), _g(-1), _g(-4, 1), _g(-2), _g(-1, -1),
        _g(-2, -1), _g(-4, -1), _g(-6), _g(-9, 1), _g(-6, -1),
    )),
    (((11, 0), (15, 0)), (
        _g(-2, -5), _g(-1), _g(-5, 1), _g(-4), _g(-4, -1),
        _g(-3, -2), _g(-3, -3), _g(-3, -3), _g(-4, -3), _g(-2, -5),
    )),
    (((13, 0), (18, 0)), (
        _g(-2, -7), _g(-1), _g(1, -5), _g(4, -7), _g(6, -9),
        _g(5, -9), _g(3, -9), _g(1, -8), _g(-2, -7), _g(2, -10),
    )),
    (((19, 0), (31, 0)), (
        _g(5, -9), _g(0), _g(-2), _g(1, -2), _g(3, -4),
        _g(2, -4), _g(0, -4), _g(-2, -3), _g(-5, -2), _g(-1, -5),
    )),
    (((20, 0),), (
        _g(4, -9), _g(-1), _g(-1, -1), _g(2, -3), _g(4, -5),
        _g(3, -5), _g(1, -5), _g(-1, -4), _g(-4, -3), _g(0, -6),
    )),
    (((21, 0),), (
        _g(-8, -3), _g(-1), _g(-7, 2, 3), _g(-2, 1, 3), _g(0),
        _g(1, -5, 3), _g(-1, -10, 3), _g(-1, -4),
        _g(-8, -14, 3), _g(0, -6),
    )),
    (((23, 0),), (
        _g(-9, -4), _g(5, -9, 2), _g(1, -9, 2), _g(7, -13, 2),
        _g(11, -17, 2), _g(9, -17, 2), _g(5, -17, 2),
        _g(1, -15, 2), _g(-5, -13, 2), _g(3, -19, 2),
    )),
)


F1, F2, F4, F5, F7, F8, F11, F13, F19, F20, F21, F23 = (
    D27_FACET_REPRESENTATIVES
)
D27_REDUNDANCY_SUPPORTS: tuple[
    tuple[tuple[int, int], tuple[tuple[int, int], ...]], ...
] = (
    ((0, 0), (F19, F20, F21)),
    ((3, 0), (F1, F2, F4, F5, F8, F11, F19, F21)),
    ((14, 0), (F1, F2, F4, F5, F7, F8, F11, F21)),
    ((22, 0), (F1, F2, F4, F5, F7, F8, F11, F19, F21)),
    ((24, 0), (F1, F2, F4, F5, F7, F8, F11, F21)),
    ((25, 0), (F2, F4, F5, F8, F11, F13, F19, F20, F21)),
    ((26, 0), (F1, F2, F4, F5, F7, F8, F11, F19, F21)),
    ((27, 0), (F1, F2, F4, F5, F7, F8, F11, F19, F20, F21)),
    ((28, 0), (F1, F2, F4, F5, F7, F8, F11, F21)),
    ((29, 0), (F1, F2, F4, F5, F7, F8, F11, F19, F21)),
    ((30, 0), (F1, F2, F4, F5, F7, F8, F11, F21)),
    ((32, 0), (F1, F2, F4, F5, F7, F8, F11, F21)),
)


D27_SIMPLE_INTERNAL_CROSSINGS = (
    ((20, 0), 86, 2, "3e8dccdc9aedb61f0e32c9b8ffb5f61e893986a545a0d6108814f362494eb3b9"),
    ((21, 0), 14, 2, "24cbab8eba94d9ecb0c569f85e22b3569f4739dff28803cdd8209ecd259058e9"),
    ((23, 0), 61, 2, "a14712879a75637ed4330d2a9b1c9a8629c940195636f231ff952b1d3bdba785"),
)
D27_SIMULTANEOUS_PATHS = (
    ((19, 0), (24, 0), (20, 1), (26, 0)),
    ((31, 0), (25, 0), (24, 0), (16, 1)),
)
D27_SIMULTANEOUS_REPRESENTATIVE_SEQUENCES = (
    (27, 11, 11, 74, 27),
    (27, 74, 11, 11, 27),
)


def _pair(row: Sequence[Golden], point: Sequence[Golden]) -> Golden:
    return sum((a * b for a, b in zip(row, point)), Golden(0))


def analyze_d27(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")

    atlas = ShimadaDisc11Atlas(path)
    basis = transported_identity_eigenbasis()
    d68 = atlas.representative_state(68)
    d86 = atlas.neighboring_state(d68, 11, 1)
    if d86 is None:
        raise ValueError("the physical D_86 state is unavailable")
    d27 = atlas.neighboring_state(d86, 24, 0)
    if d27 is None or d27.index != 27:
        raise ValueError("D_86 wall (24, 0) did not reach D_27")

    physical_orbits = atlas.state_wall_orbits(d27)
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
    for target, support in D27_REDUNDANCY_SUPPORTS:
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
    for target_class, witness in D27_FACET_WITNESSES:
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
        if wall_class == ((24, 0),)
    )
    d86_interior = _crossing_point(D68_INTERNAL_FACET_WITNESSES[1][1], 2)
    d27_interior = tuple(
        d86_face[index] - Golden(1) / 2 * (d86_interior[index] - d86_face[index])
        for index in range(10)
    )
    generic_audit = audit_golden_halfspace_certificate(
        coefficients,
        facet_classes=D27_FACET_CLASSES,
        facet_witnesses={
            wall_class[0]: witness for wall_class, witness in D27_FACET_WITNESSES
        },
        redundancies=dependencies,
        interior=d27_interior,
        quadratic_gram=tuple(
            tuple(K3_TRACE_SCALAR * entry for entry in row) for row in field_gram()
        ),
    )

    chamber = atlas.chamber(27)
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
                raise ValueError("a D_27 root wall recovered a nonintegral root")
            root = tuple(int(entry) for entry in primal)
            if (Matrix(1, 20, root) * gram * Matrix(20, 1, root))[0] != -2:
                raise ValueError("a D_27 root wall did not have square -2")
            root_labels.append((orbit_index, member_index))
            roots.append(root)
    root_filter = analyze_conjugated_roots(tuple(roots))
    surviving_root_labels = {
        root_labels[index] for index in root_filter["surviving_zero_based_indices"]
    }
    root_facet_labels = {label for wall_class in D27_ROOT_FACET_CLASSES for label in wall_class}

    witness_by_label = {
        wall_class[0]: witness
        for wall_class, witness in D27_FACET_WITNESSES
        if wall_class in D27_INTERNAL_FACET_CLASSES
    }
    d55 = atlas.neighboring_state(d86, 20, 0)
    if d55 is None:
        raise ValueError("the comparison D_55 state is unavailable")
    d14_from_d55 = atlas.neighboring_state(d55, 13, 1)
    if d14_from_d55 is None:
        raise ValueError("the comparison D_14 state is unavailable")
    known_states = {
        "D68": d68,
        "D86": d86,
        "D55": d55,
        "D27": d27,
        "D14_from_D55": d14_from_d55,
    }
    simple_records = []
    simple_states = []
    simple_verified = True
    for label, expected_target, denominator, expected_wall_hash in D27_SIMPLE_INTERNAL_CROSSINGS:
        neighbor = atlas.neighboring_state(d27, *label)
        if neighbor is None:
            raise ValueError("a simple D_27 internal facet became a K3 boundary")
        face = witness_by_label[label]
        crossing_point = tuple(
            face[index] - Golden(1) / denominator * (d27_interior[index] - face[index])
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
        simple_verified &= verified
        simple_states.append(neighbor)
        simple_records.append(
            {
                "D27_wall": list(label),
                "target_representative": neighbor.index,
                "target_physical_wall_set_sha256": wall_hash,
                "physical_matches_among_known_states": physical_matches,
                "crossing_denominator": denominator,
                "target_interior": point_record,
                "verified": verified,
            }
        )

    simultaneous_face = witness_by_label[(19, 0)]
    simultaneous_crossing_point = tuple(
        simultaneous_face[index]
        - Golden(1) / 3 * (d27_interior[index] - simultaneous_face[index])
        for index in range(10)
    )
    path_records = []
    path_endpoints = []
    simultaneous_verified = True
    for path_specification, expected_sequence in zip(
        D27_SIMULTANEOUS_PATHS, D27_SIMULTANEOUS_REPRESENTATIVE_SEQUENCES
    ):
        states = [d27]
        current = d27
        for label in path_specification:
            neighbor = atlas.neighboring_state(current, *label)
            if neighbor is None:
                raise ValueError("the simultaneous D_27 path hit a K3 wall")
            states.append(neighbor)
            current = neighbor
        face_records = [_point_record(atlas, state, simultaneous_face, basis) for state in states]
        final_record = _point_record(atlas, current, simultaneous_crossing_point, basis)
        sequence = tuple(state.index for state in states)
        verified = (
            sequence == expected_sequence
            and all(record["q_square_positive"] for record in face_records)
            and all(record["negative_wall_pairings"] == 0 for record in face_records)
            and all(record["zero_wall_pairings"] == 2 for record in face_records)
            and final_record["q_square_positive"]
            and final_record["negative_wall_pairings"] == 0
            and final_record["zero_wall_pairings"] == 0
            and _canonical_sha(current.row_isometry) == D27_OPPOSITE_STATE_SHA256
            and _physical_wall_set_sha(atlas, current)
            == D27_OPPOSITE_PHYSICAL_WALL_SET_SHA256
        )
        simultaneous_verified &= verified
        path_endpoints.append(current)
        path_records.append(
            {
                "path": [list(label) for label in path_specification],
                "representative_sequence": list(sequence),
                "face_closure_records": face_records,
                "opposite_side_interior": final_record,
                "endpoint_state_sha256": _canonical_sha(current.row_isometry),
                "endpoint_physical_wall_set_sha256": _physical_wall_set_sha(atlas, current),
                "verified": verified,
            }
        )
    simultaneous_paths_same_physical_endpoint = (
        path_endpoints[0].row_isometry == path_endpoints[1].row_isometry
        and _physical_wall_set(atlas, path_endpoints[0])
        == _physical_wall_set(atlas, path_endpoints[1])
    )
    simultaneous_endpoint_is_new_physical_D27 = (
        path_endpoints[0].index == d27.index
        and _physical_wall_set(atlas, path_endpoints[0])
        != _physical_wall_set(atlas, d27)
    )

    d86_return = simple_states[0]
    d14_from_d27 = simple_states[1]
    return_is_physical_d86 = (
        d86_return.row_isometry == d86.row_isometry
        and _physical_wall_set(atlas, d86_return) == _physical_wall_set(atlas, d86)
    )
    physical_d14_diamond = (
        d14_from_d27.row_isometry == d14_from_d55.row_isometry
        and _physical_wall_set(atlas, d14_from_d27)
        == _physical_wall_set(atlas, d14_from_d55)
    )
    entry_wall = coefficients[(20, 0)]
    d86_entry = _restricted_coefficients(atlas.state_wall_orbits(d86)[24][0], basis)
    entry_wall_opposite = entry_wall == tuple(-value for value in d86_entry)

    facet_sets = {frozenset(item) for item in D27_FACET_CLASSES}
    redundant_representatives = {
        record.members[0]
        for record in class_records
        if frozenset(record.members) not in facet_sets
    }
    partition_verified = (
        classes == D27_POSITIVE_PROPORTIONAL_CLASSES
        and len(coefficients) == 33
        and len(classes) == 24
        and len(facet_sets) == 12
        and len(redundant_representatives) == len(D27_REDUNDANCY_SUPPORTS) == 12
        and {target for target, _ in D27_REDUNDANCY_SUPPORTS}
        == redundant_representatives
    )
    fingerprints_verified = (
        _canonical_sha(d27.row_isometry) == D27_STATE_SHA256
        and _canonical_sha(physical_rows) == D27_PHYSICAL_WALL_ROWS_SHA256
        and _physical_wall_set_sha(atlas, d27) == D27_PHYSICAL_WALL_SET_SHA256
        and _canonical_sha(restricted_rows) == D27_RESTRICTED_ROWS_SHA256
    )
    root_filter_verified = (
        len(root_labels) == 19
        and root_filter["trace_minus_two_verified"] == 19
        and root_filter["field_integral_reflection_roots"] == 18
        and root_filter["field_span_rank"] == 9
        and set(root_labels) - surviving_root_labels == {(14, 0)}
        and root_facet_labels <= surviving_root_labels
    )
    traversal_verified = (
        simple_verified
        and simultaneous_verified
        and simultaneous_paths_same_physical_endpoint
        and simultaneous_endpoint_is_new_physical_D27
        and return_is_physical_d86
        and physical_d14_diamond
        and entry_wall_opposite
    )
    verified = (
        partition_verified
        and fingerprints_verified
        and witnesses_verified
        and dependences_verified
        and generic_audit.passes
        and root_filter_verified
        and traversal_verified
    )

    return {
        "verified": verified,
        "source_sha256": digest,
        "physical_state": {
            "entered_from": {"representative": 86, "wall": [24, 0]},
            "representative": 27,
            "state_sha256": D27_STATE_SHA256,
            "physical_wall_rows_sha256": D27_PHYSICAL_WALL_ROWS_SHA256,
            "physical_wall_set_sha256": D27_PHYSICAL_WALL_SET_SHA256,
            "restricted_rows_sha256": D27_RESTRICTED_ROWS_SHA256,
            "fingerprints_verified": fingerprints_verified,
        },
        "restricted_section": {
            "ambient_primitive_walls": len(coefficients),
            "positive_proportional_classes": len(classes),
            "irredundant_facets": len(D27_FACET_CLASSES),
            "root_facets": len(D27_ROOT_FACET_CLASSES),
            "internal_facets": len(D27_INTERNAL_FACET_CLASSES),
            "root_facet_classes": [
                [list(label) for label in wall_class] for wall_class in D27_ROOT_FACET_CLASSES
            ],
            "internal_facet_classes": [
                [list(label) for label in wall_class]
                for wall_class in D27_INTERNAL_FACET_CLASSES
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
            "simple_crossings": simple_records,
            "simple_representative_targets": [
                record["target_representative"] for record in simple_records
            ],
            "simultaneous_facet_class": [[19, 0], [31, 0]],
            "simultaneous_paths": path_records,
            "simultaneous_paths_same_physical_endpoint":
            simultaneous_paths_same_physical_endpoint,
            "simultaneous_endpoint_is_distinct_physical_D27":
            simultaneous_endpoint_is_new_physical_D27,
            "wall_20_is_opposite_D86_entry_wall": entry_wall_opposite,
            "wall_20_returns_to_same_physical_D86": return_is_physical_d86,
            "D27_and_D55_reach_same_physical_D14": physical_d14_diamond,
            "verified": traversal_verified,
        },
        "scope": "complete_local_D27_section_and_internal_first_star_not_global_eigen_nef_completeness",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze_d27(args.disc11), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
