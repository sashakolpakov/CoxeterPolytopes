"""Exact D_37 reduction after an exact D_84/D_37 root-facet screen.

Both states are the genuinely new neighbors of physical D_61.  The branch
screen proves exactly that D_84 has six root facets while D_37 has nine, so
D_37 is selected.  Its complete restricted section and both internal
crossings are then certified over Q(phi).
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
from scripts.analyze_shimada_disc11_d61 import D61_FACET_WITNESSES
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


D37_STATE_SHA256 = "05cd8cd6b937051b00917a0ea7eef69ce0108b0a11e18a7cddde065376a9008e"
D37_PHYSICAL_WALL_ROWS_SHA256 = (
    "82cc6a1b480332399a9756a758ac84a427907338fcf6e439de724d14ee7aaf57"
)
D37_PHYSICAL_WALL_SET_SHA256 = (
    "0998d5e9845d9cce86a9b8d8c41a828aa50c42a01efb9b8bc57a324b68fe536c"
)
D37_RESTRICTED_ROWS_SHA256 = (
    "7ee3d170213a5a5880c756607fb5e9f48bd75d4a638201363d0622ba6f9306be"
)


D37_POSITIVE_PROPORTIONAL_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0), (13, 0)),
    ((1, 0), (5, 0)),
    ((2, 0), (19, 0)),
    ((3, 0), (12, 0)),
    ((4, 0), (16, 0)),
    ((6, 0), (17, 0)),
    ((7, 0), (18, 0)),
    ((8, 0), (14, 0)),
    ((9, 0), (15, 0)),
    ((10, 0),),
    ((11, 0),),
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
)
D37_ROOT_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0), (13, 0)),
    ((1, 0), (5, 0)),
    ((2, 0), (19, 0)),
    ((3, 0), (12, 0)),
    ((4, 0), (16, 0)),
    ((6, 0), (17, 0)),
    ((7, 0), (18, 0)),
    ((8, 0), (14, 0)),
    ((9, 0), (15, 0)),
)
D37_INTERNAL_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((20, 0),),
    ((21, 0),),
)
D37_FACET_CLASSES = D37_ROOT_FACET_CLASSES + D37_INTERNAL_FACET_CLASSES
D37_FACET_REPRESENTATIVES = tuple(item[0] for item in D37_FACET_CLASSES)


D37_FACET_WITNESSES: tuple[
    tuple[tuple[tuple[int, int], ...], tuple[Golden, ...]], ...
] = (
    (((0, 0), (13, 0)), (
        _g(59, -46), _g(-52, 29), _g(23, -20), _g(21, -19),
        _g(18, -18), _g(19, -19), _g(16, -18), _g(16, -18),
        _g(15, -18), _g(14, -18),
    )),
    (((1, 0), (5, 0)), (
        _g(-4, -7), _g(4, -6), _g(-6, -1), _g(-9, 1),
        _g(-13, 3), _g(-11, 1), _g(-11), _g(-10, -1),
        _g(-8, -3), _g(-10, -2),
    )),
    (((2, 0), (19, 0)), (
        _g(-18, -39, 5), _g(28, -36, 5), _g(-22, -11, 5),
        _g(-37, -1, 5), _g(-57, 9, 5), _g(-47, -1, 5),
        _g(-47, -6, 5), _g(-42, -11, 5), _g(-42, -16, 5),
        _g(-52, -11, 5),
    )),
    (((3, 0), (12, 0)), (
        _g(-4, -37, 5), _g(54, -48, 5), _g(4, -23, 5),
        _g(-11, -13, 5), _g(-31, -3, 5), _g(-26, -8, 5),
        _g(-41, -3, 5), _g(-36, -8, 5), _g(-41, -8, 5),
        _g(-51, -3, 5),
    )),
    (((4, 0), (16, 0)), (
        _g(-9, -47, 5), _g(29, -38, 5), _g(-11, -18, 5),
        _g(-26, -8, 5), _g(-46, 2, 5), _g(-36, -8, 5),
        _g(-46, -8, 5), _g(-41, -13, 5), _g(-41, -18, 5),
        _g(-51, -13, 5),
    )),
    (((6, 0), (17, 0)), (
        _g(-21, -33, 5), _g(16, -27, 5), _g(-24, -7, 5),
        _g(-29, -2, 5), _g(-39, 3, 5), _g(-29, -7, 5),
        _g(-39, -7, 5), _g(-34, -12, 5), _g(-34, -17, 5),
        _g(-34, -17, 5),
    )),
    (((7, 0), (18, 0)), (
        _g(-19, -37, 5), _g(24, -33, 5), _g(-16, -13, 5),
        _g(-31, -3, 5), _g(-41, 2, 5), _g(-31, -8, 5),
        _g(-41, -8, 5), _g(-36, -13, 5), _g(-36, -18, 5),
        _g(-46, -13, 5),
    )),
    (((8, 0), (14, 0)), (
        _g(-8, -39, 5), _g(48, -46, 5), _g(-7, -16, 5),
        _g(-22, -6, 5), _g(-47, 9, 5), _g(-37, -1, 5),
        _g(-47, -1, 5), _g(-42, -6, 5), _g(-42, -11, 5),
        _g(-57, -1, 5),
    )),
    (((9, 0), (15, 0)), (
        _g(55, -45), _g(-51, 28), _g(19, -18), _g(16, -16),
        _g(12, -14), _g(14, -16), _g(12, -16), _g(13, -17),
        _g(13, -18), _g(11, -17),
    )),
    (((20, 0),), (
        _g(-17, -41, 5), _g(32, -39, 5), _g(-8, -19, 5),
        _g(-23, -9, 5), _g(-43, 1, 5), _g(-33, -9, 5),
        _g(-43, -9, 5), _g(-38, -14, 5), _g(-38, -19, 5),
        _g(-48, -14, 5),
    )),
    (((21, 0),), (
        _g(-18, -39, 5), _g(28, -36, 5), _g(-12, -16, 5),
        _g(-27, -6, 5), _g(-47, 4, 5), _g(-37, -6, 5),
        _g(-47, -6, 5), _g(-42, -11, 5), _g(-42, -16, 5),
        _g(-52, -11, 5),
    )),
)


F0, F1, F2, F3, F4, F6, F7, F8, F9, F20, F21 = D37_FACET_REPRESENTATIVES
D37_REDUNDANCY_SUPPORTS: tuple[
    tuple[tuple[int, int], tuple[tuple[int, int], ...]], ...
] = (
    ((10, 0), (F4, F9, F20, F21)),
    ((11, 0), (F0, F1, F2, F3, F6, F7, F8, F21)),
    ((22, 0), (F0, F1, F2, F3, F4, F6, F7, F8, F20, F21)),
    ((23, 0), (F0, F1, F2, F3, F6, F7, F8, F21)),
    ((24, 0), (F0, F1, F2, F3, F6, F7, F8, F9, F20, F21)),
    ((25, 0), (F0, F3, F4, F6, F7, F8, F9, F20, F21)),
    ((26, 0), (F0, F3, F4, F6, F7, F8, F9, F20, F21)),
    ((27, 0), (F0, F1, F2, F3, F6, F7, F8, F21)),
    ((28, 0), (F0, F1, F2, F3, F6, F7, F8, F21)),
    ((29, 0), (F4, F9, F20)),
    ((30, 0), (F0, F1, F2, F3, F6, F7, F8, F21)),
)
D37_INTERNAL_CROSSINGS = (
    ((20, 0), 61, 2, "a14712879a75637ed4330d2a9b1c9a8629c940195636f231ff952b1d3bdba785"),
    ((21, 0), 22, 2, "33ed32de6e1015ec053042afc6998693081814c85dfd8e3c6f0044214cedba39"),
)


# Exact root-only screen for the competing physical D_84 branch.
D84_ROOT_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0), (6, 0)),
    ((1, 0), (11, 0)),
    ((4, 0), (9, 0)),
    ((7, 0), (13, 0)),
    ((8, 0), (17, 0)),
    ((10, 0), (14, 0)),
)
D84_REDUNDANT_ROOT_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((2, 0),), ((3, 0),), ((5, 0),), ((12, 0),), ((15, 0),), ((16, 0),),
)
D84_ROOT_FACET_WITNESSES: tuple[
    tuple[tuple[tuple[int, int], ...], tuple[Golden, ...]], ...
] = (
    (((0, 0), (6, 0)), (
        _g(-123, 21, 5), _g(67, -59, 5), _g(104, -103, 15),
        _g(97, -74, 15), _g(0), _g(22, -44, 15),
        _g(-16, -58, 15), _g(-3, -29, 5), _g(-17, -116, 15),
        _g(-18, -24, 5),
    )),
    (((1, 0), (11, 0)), (
        _g(-2, -31, 5), _g(77, -59, 5), _g(32, -34, 5),
        _g(17, -24, 5), _g(-28, 1, 5), _g(-23, -4, 5),
        _g(-38, 1, 5), _g(-33, -4, 5), _g(-38, -4, 5),
        _g(-63, 11, 5),
    )),
    (((4, 0), (9, 0)), (
        _g(-27, 5), _g(15, -13), _g(26, -24, 3), _g(19, -15, 3),
        _g(0), _g(4, -9, 3), _g(-4, -12, 3), _g(-1, -6),
        _g(-5, -24, 3), _g(-4, -5),
    )),
    (((7, 0), (13, 0)), (
        _g(-147, 29, 5), _g(83, -71, 5), _g(136, -127, 15),
        _g(113, -86, 15), _g(0), _g(8, -41, 15),
        _g(-14, -67, 15), _g(-7, -31, 5), _g(-13, -134, 15),
        _g(-22, -26, 5),
    )),
    (((8, 0), (17, 0)), (
        _g(-29, 10), _g(37, -28, 2), _g(19, -18, 2),
        _g(13, -14, 2), _g(-5, -4, 2), _g(-1, -8, 2),
        _g(-1, -10, 2), _g(1, -12, 2), _g(1, -14, 2),
        _g(-9, -8, 2),
    )),
    (((10, 0), (14, 0)), (
        _g(-36, 15), _g(49, -35, 2), _g(33, -27, 2),
        _g(29, -25, 2), _g(13, -17, 2), _g(15, -19, 2),
        _g(9, -17, 2), _g(9, -17, 2), _g(7, -17, 2),
        _g(-1, -13, 2),
    )),
)
D84_REDUNDANT_ROOT_SUPPORTS = (
    ((2, 0), ((18, 0), (24, 0))),
    ((3, 0), ((0, 0), (1, 0), (4, 0), (7, 0), (8, 0), (10, 0), (20, 0), (27, 0))),
    ((5, 0), ((0, 0), (1, 0), (4, 0), (7, 0), (8, 0), (10, 0), (20, 0))),
    ((12, 0), ((19, 0), (23, 0))),
    ((15, 0), ((20, 0), (27, 0))),
    ((16, 0), ((18, 0), (23, 0))),
)


def _pair(row: Sequence[Golden], point: Sequence[Golden]) -> Golden:
    return sum((a * b for a, b in zip(row, point)), Golden(0))


def _d84_root_screen(atlas: ShimadaDisc11Atlas, d61, basis) -> dict[str, object]:
    d84 = atlas.neighboring_state(d61, 18, 0)
    if d84 is None or d84.index != 84:
        raise ValueError("D_61 wall (18, 0) did not reach D_84")
    rows = {
        (orbit_index, member_index): _restricted_coefficients(wall, basis)
        for orbit_index, orbit in enumerate(atlas.state_wall_orbits(d84))
        for member_index, wall in enumerate(orbit)
    }
    root_rows = {
        label: row for label, row in rows.items()
        if atlas.chamber(84).adjacencies[label[0]].kind == "_minustwo"
    }
    root_classes = tuple(item.members for item in positive_proportional_classes(root_rows))
    expected_classes = D84_ROOT_FACET_CLASSES + D84_REDUNDANT_ROOT_CLASSES
    witnesses_verified = True
    for target_class, witness in D84_ROOT_FACET_WITNESSES:
        target = set(target_class)
        pairings = {label: _pair(row, witness) for label, row in rows.items()}
        square = K3_TRACE_SCALAR * field_inner_product(field_gram(), witness, witness)
        witnesses_verified &= (
            all(pairings[label] == 0 for label in target)
            and all(value.sign() > 0 for label, value in pairings.items() if label not in target)
            and square.sign() > 0
        )
    dependences_verified = True
    dependence_records = []
    for target, support in D84_REDUNDANT_ROOT_SUPPORTS:
        solution = _solve_fixed_support(rows[target], tuple(rows[label] for label in support))
        verified = all(value.sign() > 0 for value in solution)
        dependences_verified &= verified
        dependence_records.append(
            {
                "target": list(target),
                "support": [list(label) for label in support],
                "coefficients": [value.to_json() for value in solution],
                "verified": verified,
            }
        )
    verified = (
        {frozenset(item) for item in root_classes}
        == {frozenset(item) for item in expected_classes}
        and len(root_classes) == len(expected_classes)
        and len(root_rows) == 18
        and witnesses_verified
        and dependences_verified
    )
    return {
        "verified": verified,
        "representative": 84,
        "ambient_root_walls": len(root_rows),
        "restricted_root_classes": len(root_classes),
        "root_facets": len(D84_ROOT_FACET_CLASSES),
        "root_facet_classes": [
            [list(label) for label in wall_class] for wall_class in D84_ROOT_FACET_CLASSES
        ],
        "redundant_root_classes": [
            [list(label) for label in wall_class] for wall_class in D84_REDUNDANT_ROOT_CLASSES
        ],
        "facet_witnesses_verified": witnesses_verified,
        "redundant_root_dependences": dependence_records,
    }


def analyze_d37(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")

    atlas = ShimadaDisc11Atlas(path)
    basis = transported_identity_eigenbasis()
    states = _physical_states(atlas)
    d30 = states["d30_via_d181"]
    d61 = atlas.neighboring_state(d30, 22, 0)
    if d61 is None:
        raise ValueError("the common physical D_61 state is unavailable")
    d37 = atlas.neighboring_state(d61, 23, 0)
    if d37 is None or d37.index != 37:
        raise ValueError("D_61 wall (23, 0) did not reach D_37")

    d84_screen = _d84_root_screen(atlas, d61, basis)
    branch_selection_verified = (
        d84_screen["verified"]
        and d84_screen["root_facets"] == 6
        and len(D37_ROOT_FACET_CLASSES) == 9
    )

    physical_orbits = atlas.state_wall_orbits(d37)
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
    for target, support in D37_REDUNDANCY_SUPPORTS:
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
    for target_class, witness in D37_FACET_WITNESSES:
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
    d61_face = next(
        witness for wall_class, witness in D61_FACET_WITNESSES
        if wall_class == ((23, 0),)
    )
    d37_interior = tuple(
        d61_face[index] - (d61_interior[index] - d61_face[index]) / 3
        for index in range(10)
    )
    generic_audit = audit_golden_halfspace_certificate(
        coefficients,
        facet_classes=D37_FACET_CLASSES,
        facet_witnesses={
            wall_class[0]: witness for wall_class, witness in D37_FACET_WITNESSES
        },
        redundancies=dependencies,
        interior=d37_interior,
        quadratic_gram=tuple(
            tuple(K3_TRACE_SCALAR * entry for entry in row) for row in field_gram()
        ),
    )

    chamber = atlas.chamber(37)
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
                raise ValueError("a D_37 root wall recovered a nonintegral root")
            root = tuple(int(entry) for entry in primal)
            if (Matrix(1, 20, root) * gram * Matrix(20, 1, root))[0] != -2:
                raise ValueError("a D_37 root wall did not have square -2")
            root_labels.append((orbit_index, member_index))
            roots.append(root)
    root_filter = analyze_conjugated_roots(tuple(roots))
    surviving_root_labels = {
        root_labels[index] for index in root_filter["surviving_zero_based_indices"]
    }
    root_facet_labels = {label for wall_class in D37_ROOT_FACET_CLASSES for label in wall_class}

    witness_by_label = {
        wall_class[0]: witness
        for wall_class, witness in D37_FACET_WITNESSES
        if wall_class in D37_INTERNAL_FACET_CLASSES
    }
    crossing_records = []
    crossing_states = []
    crossings_verified = True
    for label, expected_target, denominator, expected_wall_hash in D37_INTERNAL_CROSSINGS:
        neighbor = atlas.neighboring_state(d37, *label)
        if neighbor is None:
            raise ValueError("a D_37 internal facet became a K3 boundary")
        face = witness_by_label[label]
        crossing_point = tuple(
            face[index] - (d37_interior[index] - face[index]) / denominator
            for index in range(10)
        )
        point_record = _point_record(atlas, neighbor, crossing_point, basis)
        wall_hash = _physical_wall_set_sha(atlas, neighbor)
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
                "D37_wall": list(label),
                "target_representative": neighbor.index,
                "target_physical_wall_set_sha256": wall_hash,
                "crossing_denominator": denominator,
                "target_interior": point_record,
                "verified": verified,
            }
        )
    return_to_d61 = (
        crossing_states[0].row_isometry == d61.row_isometry
        and _physical_wall_set(atlas, crossing_states[0]) == _physical_wall_set(atlas, d61)
    )
    entry_wall = coefficients[(20, 0)]
    d61_entry = _restricted_coefficients(atlas.state_wall_orbits(d61)[23][0], basis)
    entry_wall_opposite = entry_wall == tuple(-value for value in d61_entry)

    facet_sets = {frozenset(item) for item in D37_FACET_CLASSES}
    redundant_representatives = {
        record.members[0]
        for record in class_records
        if frozenset(record.members) not in facet_sets
    }
    partition_verified = (
        classes == D37_POSITIVE_PROPORTIONAL_CLASSES
        and len(coefficients) == 31
        and len(classes) == 22
        and len(facet_sets) == 11
        and len(redundant_representatives) == len(D37_REDUNDANCY_SUPPORTS) == 11
        and {target for target, _ in D37_REDUNDANCY_SUPPORTS}
        == redundant_representatives
    )
    fingerprints_verified = (
        _canonical_sha(d37.row_isometry) == D37_STATE_SHA256
        and _canonical_sha(physical_rows) == D37_PHYSICAL_WALL_ROWS_SHA256
        and _physical_wall_set_sha(atlas, d37) == D37_PHYSICAL_WALL_SET_SHA256
        and _canonical_sha(restricted_rows) == D37_RESTRICTED_ROWS_SHA256
    )
    root_filter_verified = (
        len(root_labels) == 20
        and root_filter["trace_minus_two_verified"] == 20
        and root_filter["field_integral_reflection_roots"] == 19
        and root_filter["field_span_rank"] == 10
        and set(root_labels) - surviving_root_labels == {(11, 0)}
        and root_facet_labels <= surviving_root_labels
    )
    verified = (
        branch_selection_verified
        and partition_verified
        and fingerprints_verified
        and witnesses_verified
        and dependences_verified
        and generic_audit.passes
        and root_filter_verified
        and crossings_verified
        and return_to_d61
        and entry_wall_opposite
    )

    return {
        "verified": verified,
        "source_sha256": digest,
        "branch_selection": {
            "D84_root_screen": d84_screen,
            "D37_root_facets": len(D37_ROOT_FACET_CLASSES),
            "selected": 37,
            "verified": branch_selection_verified,
        },
        "physical_state": {
            "entered_from": {"representative": 61, "wall": [23, 0]},
            "representative": 37,
            "state_sha256": D37_STATE_SHA256,
            "physical_wall_rows_sha256": D37_PHYSICAL_WALL_ROWS_SHA256,
            "physical_wall_set_sha256": D37_PHYSICAL_WALL_SET_SHA256,
            "restricted_rows_sha256": D37_RESTRICTED_ROWS_SHA256,
            "fingerprints_verified": fingerprints_verified,
        },
        "restricted_section": {
            "ambient_primitive_walls": len(coefficients),
            "positive_proportional_classes": len(classes),
            "irredundant_facets": len(D37_FACET_CLASSES),
            "root_facets": len(D37_ROOT_FACET_CLASSES),
            "internal_facets": len(D37_INTERNAL_FACET_CLASSES),
            "root_facet_classes": [
                [list(label) for label in wall_class] for wall_class in D37_ROOT_FACET_CLASSES
            ],
            "internal_facet_classes": [
                [list(label) for label in wall_class]
                for wall_class in D37_INTERNAL_FACET_CLASSES
            ],
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
            "field_root_nonfacet_walls": [[10, 0]],
            "verified": root_filter_verified,
        },
        "internal_first_star": {
            "crossings": crossing_records,
            "representative_targets": [61, 22],
            "returns_to_same_physical_D61": return_to_d61,
            "return_wall_is_opposite_D61_entry_wall": entry_wall_opposite,
            "new_physical_neighbor": 22,
            "verified": crossings_verified,
        },
        "scope": "exact_D84_root_screen_complete_local_D37_section_and_internal_first_star",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze_d37(args.disc11), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
