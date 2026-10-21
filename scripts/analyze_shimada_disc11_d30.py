"""Exact reduction of the physical discriminant-11 chamber D30.

The same physical state is reached by the two certified K3 chamber walks

``D68 -> D181 -> D30`` and ``D68 -> D86 -> D30``.

This verifier proves that equality at the level of Shimada row transports,
restricts all 34 primitive walls to the fixed transported golden eigenspace,
and audits a complete irredundant half-space certificate over ``Q(phi)``.
It then crosses every genuine internal facet once.  The result is local: it
does not assert global eigen-nef completeness, finite volume, or compactness.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from sympy import Matrix

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.geometric.golden_halfspaces import (
    audit_golden_halfspace_certificate,
    golden_dot,
)
from coxeter_polytopes.geometric.k3 import K3_TRACE_SCALAR, field_inner_product
from scripts.analyze_shimada_disc11_eigen_nef import (
    _point_record,
    _restricted_coefficients,
)
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.shimada_disc11_atlas import ShimadaChamberState, ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    DISC11_SOURCE_SHA256,
    analyze_conjugated_roots,
    expected_shimada_gram,
    transported_identity_eigenbasis,
)


WallLabel = tuple[int, int]
WallClass = tuple[WallLabel, ...]
GoldenVector = tuple[Golden, ...]


D30_STATE_SHA256 = "89ee0f1b46508b99b776d83b38d9c70bc7293a142f05c84f7d5256ba21810a6c"
D30_PHYSICAL_WALL_ROWS_SHA256 = (
    "296be6eaa150b5b3876151254022d71a3b7cc89e5b1a6994d0961a8b0c31524c"
)
D30_RESTRICTED_ROWS_SHA256 = (
    "bd49c25c64439562060cb9b32f3aa81071037756b442f6762ecd3127aa36284a"
)


def _g(a: int, b: int = 0, denominator: int = 1) -> Golden:
    return Golden(Fraction(a, denominator), Fraction(b, denominator))


# Complete positive-proportional partition of the 34 primitive physical walls.
D30_POSITIVE_PROPORTIONAL_CLASSES: tuple[WallClass, ...] = (
    ((0, 0), (2, 0)),
    ((1, 0), (15, 0)),
    ((3, 0), (16, 0)),
    ((4, 0), (9, 0)),
    ((5, 0), (6, 0)),
    ((7, 0), (17, 0)),
    ((8, 0), (14, 0)),
    ((10, 0),),
    ((11, 0), (12, 0)),
    ((13, 0),),
    *(((index, 0),) for index in range(18, 34)),
)


D30_ROOT_FACET_CLASSES: tuple[WallClass, ...] = (
    ((0, 0), (2, 0)),
    ((1, 0), (15, 0)),
    ((3, 0), (16, 0)),
    ((4, 0), (9, 0)),
    ((5, 0), (6, 0)),
    ((7, 0), (17, 0)),
    ((8, 0), (14, 0)),
    ((11, 0), (12, 0)),
)
D30_INTERNAL_FACET_CLASSES: tuple[WallClass, ...] = (
    ((19, 0),),
    ((20, 0),),
    ((22, 0),),
    ((26, 0),),
)
D30_FACET_CLASSES = D30_ROOT_FACET_CLASSES + D30_INTERNAL_FACET_CLASSES
D30_FACET_REPRESENTATIVES: tuple[WallLabel, ...] = tuple(
    wall_class[0] for wall_class in D30_FACET_CLASSES
)


# An exact strict point obtained on the D30 side of the D181/D30 face.
D30_INTERIOR: GoldenVector = (
    _g(-3353, 1041, 120),
    _g(1041, -982, 120),
    _g(-2639, 1493, 120),
    _g(-2639, 1493, 120),
    _g(-2399, 1253, 120),
    _g(-2159, 1013, 120),
    _g(-1679, 533, 120),
    _g(-1559, 413, 120),
    _g(-1199, 53, 120),
    _g(-1079, -67, 120),
)


# Every vector is timelike, vanishes on precisely its target class, and is
# strictly positive on all other primitive restricted walls.  The witnesses
# were discovered numerically but are archived and checked only exactly.
D30_FACET_WITNESSES: tuple[tuple[WallClass, GoldenVector], ...] = (
    (((0, 0), (2, 0)), (
        _g(-12), _g(-9, 0, 2), _g(-1, 0, 2), _g(0), _g(0),
        _g(-3, 0, 2), _g(-9, 0, 2), _g(-11, 0, 2),
        _g(-15, 0, 2), _g(-7),
    )),
    (((1, 0), (15, 0)), (
        _g(-10), _g(-11, 0, 3), _g(-2), _g(-7, 0, 3),
        _g(-7, 0, 2), _g(-4), _g(-16, 0, 3), _g(-17, 0, 3),
        _g(-7), _g(-15, 0, 2),
    )),
    (((3, 0), (16, 0)), (
        _g(-25, 0, 3), _g(-8, 0, 3), _g(-3, 0, 2), _g(-1),
        _g(-5, 0, 3), _g(-7, 0, 3), _g(-4), _g(-9, 0, 2),
        _g(-17, 0, 3), _g(-17, 0, 3),
    )),
    (((4, 0), (9, 0)), (
        _g(-25, 0, 2), _g(-9, 0, 2), _g(-1), _g(0), _g(0),
        _g(-2), _g(-9, 0, 2), _g(-6), _g(-8), _g(-15, 0, 2),
    )),
    (((5, 0), (6, 0)), (
        _g(-17, -17, 4), _g(-17, 0, 4), _g(-2), _g(-9, 0, 4),
        _g(-7, 0, 2), _g(-17, 0, 4), _g(-6), _g(-19, 0, 3),
        _g(-15, 0, 2), _g(-33, 0, 4),
    )),
    (((7, 0), (17, 0)), (
        _g(-34, 0, 3), _g(-17, 0, 4), _g(-2, 0, 3), _g(0),
        _g(0), _g(-3, 0, 2), _g(-4), _g(-21, 0, 4), _g(-7),
        _g(-7),
    )),
    (((8, 0), (14, 0)), (
        _g(-10), _g(-11, 0, 3), _g(-7, 0, 3), _g(-3),
        _g(-14, 0, 3), _g(-5), _g(-19, 0, 3), _g(-19, 0, 3),
        _g(-7), _g(-8),
    )),
    (((11, 0), (12, 0)), (
        _g(-28, 0, 3), _g(-10, 0, 3), _g(-2), _g(-7, 0, 3),
        _g(-7, 0, 2), _g(-23, 0, 6), _g(-16, 0, 3),
        _g(-17, 0, 3), _g(-19, 0, 3), _g(-7),
    )),
    (((19, 0),), (
        _g(43, -43, 2), _g(-5), _g(-2, 0, 3), _g(0), _g(0),
        _g(-2), _g(-14, 0, 3), _g(-19, 0, 3), _g(-17, 0, 2),
        _g(-8),
    )),
    (((20, 0),), (
        _g(-21, 0, 2), _g(21, -21, 4), _g(-2), _g(-7, 0, 3),
        _g(-7, 0, 2), _g(-13, 0, 3), _g(-6), _g(-19, 0, 3),
        _g(-15, 0, 2), _g(-8),
    )),
    (((22, 0),), (
        _g(45, -45, 2), _g(-9, 0, 2), _g(-6), _g(-19, 0, 3),
        _g(-15, 0, 2), _g(-25, 0, 3), _g(-10), _g(-21, 0, 2),
        _g(-35, 0, 3), _g(-37, 0, 3),
    )),
    (((26, 0),), (
        _g(35, -35, 2), _g(-7, 0, 2), _g(-2), _g(-2),
        _g(-10, 0, 3), _g(-4), _g(-6), _g(-19, 0, 3),
        _g(-22, 0, 3), _g(-8),
    )),
)


F0, F1, F3, F4, F5, F7, F8, F11, F19, F20, F22, F26 = (
    D30_FACET_REPRESENTATIVES
)


# One positive exact Q(phi)-linear expression for every nonfacet class.
D30_REDUNDANCY_DEPENDENCES: tuple[
    tuple[WallLabel, tuple[tuple[Golden, WallLabel], ...]], ...
] = (
    ((10, 0), (
        (3 + 4 * PHI, F0), (1 + 3 * PHI, F1), (4 + 8 * PHI, F3),
        (4 + 2 * PHI, F4), (4 + 6 * PHI, F7), (2 + 2 * PHI, F8),
        (5 + 4 * PHI, F11), (1 + 4 * PHI, F19),
    )),
    ((13, 0), ((PHI - 1, F19), (PHI - 1, F26))),
    ((18, 0), ((PHI, F5), (PHI, F20), (PHI - 1, F22))),
    ((21, 0), (
        (-3 + 7 * PHI, F0), (3 + PHI, F1), (8 + 4 * PHI, F3),
        (-4 + 6 * PHI, F4), (6 + 4 * PHI, F7), (-2 + 4 * PHI, F8),
        (-5 + 9 * PHI, F11), (-1 + 5 * PHI, F19), (PHI - 1, F26),
    )),
    ((23, 0), (
        (8 - 2 * PHI, F0), (-2 + 2 * PHI, F1), (-3 + 7 * PHI, F3),
        (8 - 3 * PHI, F4), (-1 + 5 * PHI, F7), (6 - 2 * PHI, F8),
        (12 - 4 * PHI, F11), (Golden(4), F19), (2 - PHI, F26),
    )),
    ((24, 0), (
        (1 + 3 * PHI, F0), (1 + 2 * PHI, F1), (3 + 5 * PHI, F3),
        (2 + 2 * PHI, F4), (Golden(1), F5), (3 + 3 * PHI, F7),
        (2 + PHI, F8), (3 + 3 * PHI, F11), (3 * PHI, F19),
        (PHI, F20),
    )),
    ((25, 0), (
        (2 + 4 * PHI, F0), (2 + 2 * PHI, F1), (5 + 7 * PHI, F3),
        (1 + 4 * PHI, F4), (4 + 5 * PHI, F7), (1 + 3 * PHI, F8),
        (2 + 6 * PHI, F11), (4 * PHI, F19),
    )),
    ((27, 0), (
        (1 + PHI, F0), (_g(0, 3, 2), F3), (4 + 5 * PHI, F5),
        (_g(4, 5, 2), F7), (_g(1, 0, 2), F8), (Golden(1), F11),
        (3 * PHI, F19), (_g(7, 13, 2), F20), (_g(1, 2, 2), F22),
    )),
    ((28, 0), (
        (3 * PHI, F0), (1 + PHI, F1), (3 + 3 * PHI, F3),
        (2 * PHI, F4), (1 + 2 * PHI, F5), (3 + 3 * PHI, F7),
        (PHI, F8), (3 * PHI, F11), (3 * PHI, F19),
        (1 + 3 * PHI, F20),
    )),
    ((29, 0), (
        (2 + 2 * PHI, F0), (2 * PHI, F1), (2 + 5 * PHI, F3),
        (3 + PHI, F4), (2 + 3 * PHI, F7), (2 + PHI, F8),
        (4 + 2 * PHI, F11), (2 * PHI, F19),
    )),
    ((30, 0), (
        (3 + 3 * PHI, F0), (1 + 2 * PHI, F1), (3 + 7 * PHI, F3),
        (3 + 2 * PHI, F4), (3 + 5 * PHI, F7), (3 + PHI, F8),
        (5 + 3 * PHI, F11), (1 + 3 * PHI, F19),
    )),
    ((31, 0), ((Golden(1), F5), (PHI, F20))),
    ((32, 0), (
        (2 + 6 * PHI, F0), (2 + 3 * PHI, F1), (7 + 9 * PHI, F3),
        (1 + 5 * PHI, F4), (5 + 7 * PHI, F7), (1 + 4 * PHI, F8),
        (2 + 8 * PHI, F11), (1 + 5 * PHI, F19),
    )),
    ((33, 0), (
        (4 + 5 * PHI, F0), (1 + 4 * PHI, F1), (5 + 10 * PHI, F3),
        (5 + 3 * PHI, F4), (5 + 8 * PHI, F7), (4 + 2 * PHI, F8),
        (7 + 5 * PHI, F11), (2 + 5 * PHI, F19), (Golden(1), F26),
    )),
)


D30_INTERNAL_CROSSINGS = (
    (19, 16, 10),
    (20, 86, 3),
    (22, 61, 3),
    (26, 181, 2),
)


def _canonical_sha(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _physical_states(atlas: ShimadaDisc11Atlas) -> dict[str, ShimadaChamberState]:
    d68 = atlas.representative_state(68)
    d181 = atlas.neighboring_state(d68, 14, 0)
    d86 = atlas.neighboring_state(d68, 11, 1)
    if d181 is None or d181.index != 181:
        raise ValueError("D68 wall (14,0) did not enter physical D181")
    if d86 is None or d86.index != 86:
        raise ValueError("D68 wall (11,1) did not enter physical D86")
    d30_via_d181 = atlas.neighboring_state(d181, 10, 0)
    d30_via_d86 = atlas.neighboring_state(d86, 23, 0)
    if d30_via_d181 is None or d30_via_d181.index != 30:
        raise ValueError("D181 wall (10,0) did not enter physical D30")
    if d30_via_d86 is None or d30_via_d86.index != 30:
        raise ValueError("D86 wall (23,0) did not enter physical D30")
    return {
        "d68": d68,
        "d181": d181,
        "d86": d86,
        "d30_via_d181": d30_via_d181,
        "d30_via_d86": d30_via_d86,
    }


def _restricted_rows(
    atlas: ShimadaDisc11Atlas,
    state: ShimadaChamberState,
) -> dict[WallLabel, GoldenVector]:
    basis = transported_identity_eigenbasis()
    return {
        (orbit_index, member_index): _restricted_coefficients(wall, basis)
        for orbit_index, orbit in enumerate(atlas.state_wall_orbits(state))
        for member_index, wall in enumerate(orbit)
    }


def _proportionality_record(wall_class: WallClass, rows):
    base = rows[wall_class[0]]
    pivot = next(index for index, entry in enumerate(base) if entry)
    return {
        "members": [list(label) for label in wall_class],
        "ratios_to_first": [
            (rows[label][pivot] / base[pivot]).to_json() for label in wall_class
        ],
    }


def _face_record(wall_class: WallClass, witness: GoldenVector, rows):
    target = set(wall_class)
    pairings = {label: golden_dot(row, witness) for label, row in rows.items()}
    square = K3_TRACE_SCALAR * field_inner_product(
        field_gram(), witness, witness
    )
    return {
        "wall_class": [list(label) for label in wall_class],
        "witness": [entry.to_json() for entry in witness],
        "q_square": square.to_json(),
        "timelike": square.sign() > 0,
        "target_zero": all(pairings[label] == 0 for label in target),
        "all_others_positive": all(
            value.sign() > 0
            for label, value in pairings.items()
            if label not in target
        ),
    }


def _dependence_record(target: WallLabel, terms, rows):
    value = tuple(
        sum(
            (coefficient * rows[source][index] for coefficient, source in terms),
            Golden(0),
        )
        for index in range(10)
    )
    positive = all(coefficient.sign() > 0 for coefficient, _ in terms)
    return {
        "target": list(target),
        "terms": [
            {
                "coefficient": coefficient.to_json(),
                "facet": list(source),
            }
            for coefficient, source in terms
        ],
        "all_coefficients_positive": positive,
        "exact_identity": value == rows[target],
        "verified": positive and value == rows[target],
    }


def _root_filter(atlas: ShimadaDisc11Atlas, state: ShimadaChamberState):
    chamber = atlas.chamber(30)
    physical_orbits = atlas.state_wall_orbits(state)
    gram = Matrix(expected_shimada_gram())
    inverse_gram = gram.inv()
    labels: list[WallLabel] = []
    roots: list[tuple[int, ...]] = []
    for orbit_index, adjacency in enumerate(chamber.adjacencies):
        if adjacency.kind != "_minustwo":
            continue
        for member_index, wall in enumerate(physical_orbits[orbit_index]):
            primal = Matrix(1, 20, wall) * inverse_gram
            if any(not entry.is_Integer for entry in primal):
                raise ValueError("a D30 root wall recovered a nonintegral class")
            root = tuple(int(entry) for entry in primal)
            square = (Matrix(1, 20, root) * gram * Matrix(20, 1, root))[0]
            if square != -2:
                raise ValueError("a D30 _minustwo wall did not have square -2")
            labels.append((orbit_index, member_index))
            roots.append(root)
    result = analyze_conjugated_roots(tuple(roots))
    survivors = [labels[index] for index in result["surviving_zero_based_indices"]]
    return {
        "ambient_minus_two_walls": len(labels),
        "positive_field_roots": result["positive_at_both_places"],
        "integral_field_roots": result["field_integral_reflection_roots"],
        "shell_distribution": result["field_integral_shell_distribution"],
        "field_span_rank": result["field_span_rank"],
        "non_field_root_walls": [
            list(label) for label in labels if label not in survivors
        ],
        "field_root_nonfacet_walls": [[13, 0]],
        "root_walls_on_section_facets": 16,
    }


def analyze_d30(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")

    atlas = ShimadaDisc11Atlas(path)
    states = _physical_states(atlas)
    d30 = states["d30_via_d181"]
    path_equality = d30 == states["d30_via_d86"]
    chamber = atlas.chamber(30)
    physical_orbits = atlas.state_wall_orbits(d30)
    rows = _restricted_rows(atlas, d30)
    physical_rows = [
        [[orbit_index, member_index], list(wall)]
        for orbit_index, orbit in enumerate(physical_orbits)
        for member_index, wall in enumerate(orbit)
    ]
    restricted_rows = [
        [[*label], [entry.to_json() for entry in row]]
        for label, row in sorted(rows.items())
    ]

    from coxeter_polytopes.geometric.golden_halfspaces import (
        positive_proportional_classes,
    )

    classes = positive_proportional_classes(rows)
    computed_classes = tuple(item.members for item in classes)
    face_records = [
        _face_record(wall_class, witness, rows)
        for wall_class, witness in D30_FACET_WITNESSES
    ]
    dependence_records = [
        _dependence_record(target, terms, rows)
        for target, terms in D30_REDUNDANCY_DEPENDENCES
    ]
    witnesses = {
        wall_class[0]: witness for wall_class, witness in D30_FACET_WITNESSES
    }
    q_gram = tuple(
        tuple(K3_TRACE_SCALAR * entry for entry in row) for row in field_gram()
    )
    generic_audit = audit_golden_halfspace_certificate(
        rows,
        facet_classes=D30_FACET_CLASSES,
        facet_witnesses=witnesses,
        redundancies=dict(D30_REDUNDANCY_DEPENDENCES),
        interior=D30_INTERIOR,
        quadratic_gram=q_gram,
    )

    crossing_records = []
    crossings_verified = True
    witness_by_orbit = {
        wall_class[0][0]: witness
        for wall_class, witness in D30_FACET_WITNESSES
        if wall_class in D30_INTERNAL_FACET_CLASSES
    }
    for orbit, expected_target, denominator in D30_INTERNAL_CROSSINGS:
        neighbor = atlas.neighboring_state(d30, orbit, 0)
        if neighbor is None:
            raise ValueError("a D30 internal facet became a K3 boundary")
        face = witness_by_orbit[orbit]
        crossing = tuple(
            face[index] - (D30_INTERIOR[index] - face[index]) / denominator
            for index in range(10)
        )
        face_current = _point_record(
            atlas, d30, face, transported_identity_eigenbasis()
        )
        face_neighbor = _point_record(
            atlas, neighbor, face, transported_identity_eigenbasis()
        )
        target_point = _point_record(
            atlas, neighbor, crossing, transported_identity_eigenbasis()
        )
        verified = (
            neighbor.index == expected_target
            and face_current["negative_wall_pairings"] == 0
            and face_current["zero_wall_pairings"] == 1
            and face_current["q_square_positive"]
            and face_neighbor["negative_wall_pairings"] == 0
            and face_neighbor["zero_wall_pairings"] == 1
            and face_neighbor["q_square_positive"]
            and target_point["negative_wall_pairings"] == 0
            and target_point["zero_wall_pairings"] == 0
            and target_point["q_square_positive"]
        )
        crossings_verified &= verified
        crossing_records.append(
            {
                "D30_orbit": orbit,
                "adjacency_kind": chamber.adjacencies[orbit].kind,
                "target_representative": neighbor.index,
                "crossing_denominator": denominator,
                "shared_face_in_D30": face_current,
                "shared_face_in_target": face_neighbor,
                "target_strict_interior": target_point,
                "returns_to_same_physical_D86": neighbor == states["d86"],
                "returns_to_same_physical_D181": neighbor == states["d181"],
                "verified": verified,
            }
        )

    root_orbits = [
        index
        for index, adjacency in enumerate(chamber.adjacencies)
        if adjacency.kind == "_minustwo"
    ]
    fingerprints_verified = (
        _canonical_sha(d30.row_isometry) == D30_STATE_SHA256
        and _canonical_sha(physical_rows) == D30_PHYSICAL_WALL_ROWS_SHA256
        and _canonical_sha(restricted_rows) == D30_RESTRICTED_ROWS_SHA256
    )
    root_filter = _root_filter(atlas, d30)
    verified = (
        path_equality
        and chamber.order_aut == 1
        and len(rows) == 34
        and computed_classes == D30_POSITIVE_PROPORTIONAL_CLASSES
        and len(classes) == 26
        and len(D30_FACET_CLASSES) == 12
        and all(record["timelike"] for record in face_records)
        and all(record["target_zero"] for record in face_records)
        and all(record["all_others_positive"] for record in face_records)
        and all(record["verified"] for record in dependence_records)
        and generic_audit.passes
        and root_orbits == list(range(18))
        and root_filter["integral_field_roots"] == 17
        and root_filter["non_field_root_walls"] == [[10, 0]]
        and crossings_verified
        and fingerprints_verified
    )
    return {
        "verified": verified,
        "source_sha256": digest,
        "physical_state": {
            "representative": 30,
            "D68_D181_D30_equals_D68_D86_D30": path_equality,
            "chamber_stabilizer_order": chamber.order_aut,
            "state_sha256": D30_STATE_SHA256,
            "physical_wall_rows_sha256": D30_PHYSICAL_WALL_ROWS_SHA256,
            "restricted_rows_sha256": D30_RESTRICTED_ROWS_SHA256,
            "fingerprints_verified": fingerprints_verified,
        },
        "strict_interior": _point_record(
            atlas, d30, D30_INTERIOR, transported_identity_eigenbasis()
        ),
        "restricted_section": {
            "ambient_primitive_walls": len(rows),
            "positive_proportional_classes": len(classes),
            "classes": [
                _proportionality_record(item.members, rows) for item in classes
            ],
            "irredundant_facets": len(D30_FACET_CLASSES),
            "root_facets": len(D30_ROOT_FACET_CLASSES),
            "internal_facets": len(D30_INTERNAL_FACET_CLASSES),
            "root_facet_classes": [
                [list(label) for label in wall_class]
                for wall_class in D30_ROOT_FACET_CLASSES
            ],
            "internal_facet_classes": [
                [list(label) for label in wall_class]
                for wall_class in D30_INTERNAL_FACET_CLASSES
            ],
            "facet_witnesses": face_records,
            "redundancy_dependences": dependence_records,
            "generic_exact_halfspace_kernel_verified": generic_audit.passes,
        },
        "root_filter": root_filter,
        "internal_first_star": {
            "crossings": crossing_records,
            "all_crossings_verified": crossings_verified,
            "representative_targets": [
                record["target_representative"] for record in crossing_records
            ],
        },
        "scope": (
            "complete_local_D30_section_and_internal_first_star_"
            "not_global_eigen_nef_completeness_volume_or_compactness"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze_d30(args.disc11), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
