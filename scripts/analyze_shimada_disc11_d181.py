"""Exact restricted-facet audit of the physical discriminant-11 chamber D181.

The physical state is reached from the certified D68 eigensection by crossing
the irredundant internal wall ``(14, 0)``.  This script restricts every one of
its 26 primitive Shimada walls to the *fixed* transported golden eigenspace,
collapses positive-proportional copies, and proves the exact irredundant facet
list by two independent kinds of finite certificate:

* a timelike face point for every facet; and
* a positive Q(phi)-linear dependence for every redundant half-space.

It then crosses both genuine internal facets and verifies strict timelike
interior points in the neighboring physical D30 and D68 chambers.  The result
is local to this two-star patch; it is not a global eigen-nef completeness,
finite-volume, or compactness certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Sequence

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.geometric.k3 import K3_TRACE_SCALAR, field_inner_product
from coxeter_polytopes.geometric.k3_presentation import _inverse_golden
from coxeter_polytopes.geometric.golden_halfspaces import (
    audit_golden_halfspace_certificate,
)
from scripts.analyze_shimada_disc11_eigen_nef import (
    _crossing_point,
    _pair,
    _point_record,
    _restricted_coefficients,
)
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.shimada_disc11_atlas import ShimadaChamberState, ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    D68_INTERNAL_FACET_WITNESSES,
    DISC11_SOURCE_SHA256,
    transported_identity_eigenbasis,
)


WallLabel = tuple[int, int]
WallClass = tuple[WallLabel, ...]
GoldenVector = tuple[Golden, ...]


D181_EXPECTED_PROPORTIONAL_CLASSES: tuple[WallClass, ...] = (
    ((0, 0), (0, 1)),
    ((1, 0), (1, 1)),
    ((2, 0), (6, 0)),
    ((3, 0), (3, 1)),
    ((4, 0), (4, 1)),
    ((5, 0), (5, 1)),
    ((7, 0), (7, 1)),
    ((8, 0), (8, 1)),
    ((9, 0), (9, 1)),
    ((10, 0),),
    ((10, 1),),
    ((11, 0),),
    ((12, 0),),
    ((12, 1),),
    ((13, 0),),
    ((14, 0),),
    ((15, 0),),
)


D181_ROOT_FACET_CLASSES: tuple[WallClass, ...] = (
    ((0, 0), (0, 1)),
    ((1, 0), (1, 1)),
    ((2, 0), (6, 0)),
    ((3, 0), (3, 1)),
    ((4, 0), (4, 1)),
    ((5, 0), (5, 1)),
    ((7, 0), (7, 1)),
    ((8, 0), (8, 1)),
    ((9, 0), (9, 1)),
)


D181_INTERNAL_FACET_CLASSES: tuple[WallClass, ...] = (
    ((10, 0),),
    ((11, 0),),
)


D181_FACET_CLASSES = D181_ROOT_FACET_CLASSES + D181_INTERNAL_FACET_CLASSES
D181_FACET_REPRESENTATIVES: tuple[WallLabel, ...] = tuple(
    wall_class[0] for wall_class in D181_FACET_CLASSES
)


def _g(a: int, b: int = 0, denominator: int = 1) -> Golden:
    return Golden(Fraction(a, denominator), Fraction(b, denominator))


# Each target functional is the displayed positive Q(phi)-linear combination
# of irredundant facet representatives.  These six targets are exactly the
# six nonfacet positive-proportional classes.
D181_REDUNDANCY_DEPENDENCES: tuple[
    tuple[WallLabel, tuple[tuple[Golden, WallLabel], ...]], ...
] = (
    (
        (10, 1),
        (
            (1 + PHI, (2, 0)),
            (PHI - 1, (10, 0)),
            (2 * PHI, (11, 0)),
        ),
    ),
    (
        (12, 0),
        (
            (1 + PHI, (0, 0)),
            (_g(1, 3, 2), (1, 0)),
            (_g(5, 9, 2), (2, 0)),
            (_g(2, -1, 2), (4, 0)),
            (_g(3, 5, 2), (7, 0)),
            (_g(2, 3), (8, 0)),
            (Golden(1), (9, 0)),
            (_g(2, 5, 2), (10, 0)),
            (_g(4, 5), (11, 0)),
        ),
    ),
    (
        (12, 1),
        (
            (_g(1, 4, 3), (0, 0)),
            (_g(2, 4, 3), (1, 0)),
            (_g(10, 16, 3), (2, 0)),
            (_g(4, -2, 3), (3, 0)),
            (_g(-1, 1, 3), (5, 0)),
            (_g(2, 2), (7, 0)),
            (_g(8, 7, 3), (8, 0)),
            (Golden(1), (9, 0)),
            (_g(-2, 10, 3), (10, 0)),
            (_g(12, 20, 3), (11, 0)),
        ),
    ),
    (
        (13, 0),
        (
            (2 * PHI, (0, 0)),
            (2 + 2 * PHI, (1, 0)),
            (1 + PHI, (2, 0)),
            (PHI, (3, 0)),
            (PHI, (4, 0)),
            (2 + 2 * PHI, (7, 0)),
            (2 + 2 * PHI, (8, 0)),
            (2 * PHI, (9, 0)),
            (2 * PHI, (10, 0)),
            (2 * PHI, (11, 0)),
        ),
    ),
    (
        (14, 0),
        (
            (_g(0, 5, 3), (0, 0)),
            (_g(4, 4, 3), (1, 0)),
            (_g(6, 11, 3), (2, 0)),
            (_g(0, 2, 3), (3, 0)),
            (_g(1, 1, 3), (5, 0)),
            (_g(2, 2), (7, 0)),
            (_g(7, 7, 3), (8, 0)),
            (PHI, (9, 0)),
            (_g(0, 8, 3), (10, 0)),
            (_g(7, 15, 3), (11, 0)),
        ),
    ),
    (
        (15, 0),
        (
            (Golden(1), (2, 0)),
            (PHI, (11, 0)),
        ),
    ),
)


def _physical_d181_state(atlas: ShimadaDisc11Atlas) -> ShimadaChamberState:
    state = atlas.neighboring_state(atlas.representative_state(68), 14, 0)
    if state is None or state.index != 181:
        raise ValueError("D68 wall (14,0) did not enter physical D181")
    return state


def _d181_interior_point() -> GoldenVector:
    """The exact point already certified just beyond the D68/D181 wall."""

    return _crossing_point(D68_INTERNAL_FACET_WITNESSES[2][1], 10)


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


def _positive_proportional_classes(
    rows: dict[WallLabel, GoldenVector],
) -> tuple[WallClass, ...]:
    classes: list[list[WallLabel]] = []
    for label, row in rows.items():
        for wall_class in classes:
            base = rows[wall_class[0]]
            pivot = next(index for index, entry in enumerate(base) if entry)
            ratio = row[pivot] / base[pivot]
            if ratio.sign() > 0 and row == tuple(ratio * entry for entry in base):
                wall_class.append(label)
                break
        else:
            classes.append([label])
    return tuple(tuple(wall_class) for wall_class in classes)


def _proportionality_record(
    wall_class: WallClass,
    rows: dict[WallLabel, GoldenVector],
) -> dict[str, object]:
    base = rows[wall_class[0]]
    pivot = next(index for index, entry in enumerate(base) if entry)
    ratios = tuple(rows[label][pivot] / base[pivot] for label in wall_class)
    return {
        "members": [list(label) for label in wall_class],
        "ratios_to_first": [ratio.to_json() for ratio in ratios],
        "all_ratios_positive": all(ratio.sign() > 0 for ratio in ratios),
    }


def _lorentzian_face_witness(interior: GoldenVector, functional: GoldenVector):
    """Project ``interior`` exactly to the hyperplane of ``functional``."""

    inverse = _inverse_golden(field_gram())
    normal = tuple(
        sum((inverse[i][j] * functional[j] for j in range(10)), Golden(0))
        for i in range(10)
    )
    normal_square = _pair(functional, normal)
    if not normal_square:
        raise ValueError("a proposed facet functional has isotropic dual normal")
    scale = _pair(functional, interior) / normal_square
    return tuple(interior[i] - scale * normal[i] for i in range(10))


def _face_record(
    wall_class: WallClass,
    witness: GoldenVector,
    rows: dict[WallLabel, GoldenVector],
) -> dict[str, object]:
    target = set(wall_class)
    pairings = {label: _pair(row, witness) for label, row in rows.items()}
    q_square = K3_TRACE_SCALAR * field_inner_product(
        field_gram(), witness, witness
    )
    return {
        "restricted_wall_class": [list(label) for label in wall_class],
        "face_witness": [entry.to_json() for entry in witness],
        "face_witness_q_square": q_square.to_json(),
        "face_witness_timelike": q_square.sign() > 0,
        "target_pairings_zero": all(pairings[label] == 0 for label in target),
        "all_other_primitive_wall_pairings_positive": all(
            value.sign() > 0
            for label, value in pairings.items()
            if label not in target
        ),
    }


def _dependence_record(
    target: WallLabel,
    terms: Sequence[tuple[Golden, WallLabel]],
    rows: dict[WallLabel, GoldenVector],
) -> dict[str, object]:
    value = tuple(
        sum((coefficient * rows[label][i] for coefficient, label in terms), Golden(0))
        for i in range(10)
    )
    positive = all(coefficient.sign() > 0 for coefficient, _ in terms)
    return {
        "redundant_wall": list(target),
        "terms": [
            {
                "positive_coefficient": coefficient.to_json(),
                "facet_wall": list(label),
            }
            for coefficient, label in terms
        ],
        "all_coefficients_positive": positive,
        "exact_vector_identity": value == rows[target],
        "verified": positive and value == rows[target],
    }


def _cross_face(
    face: GoldenVector,
    interior: GoldenVector,
    denominator: int = 2,
) -> GoldenVector:
    return tuple(
        face[index] - (interior[index] - face[index]) / denominator
        for index in range(10)
    )


def _internal_crossing_record(
    atlas: ShimadaDisc11Atlas,
    state: ShimadaChamberState,
    label: WallLabel,
    face: GoldenVector,
    interior: GoldenVector,
) -> dict[str, object]:
    orbit, member = label
    neighbor = atlas.neighboring_state(state, orbit, member)
    if neighbor is None:
        raise ValueError("an asserted internal facet is a K3 boundary")
    basis = transported_identity_eigenbasis()
    face_current = _point_record(atlas, state, face, basis)
    face_neighbor = _point_record(atlas, neighbor, face, basis)
    opposite = _cross_face(face, interior)
    opposite_record = _point_record(atlas, neighbor, opposite, basis)
    verified = (
        face_current["negative_wall_pairings"] == 0
        and face_current["zero_wall_pairings"] == 1
        and face_current["q_square_positive"]
        and face_neighbor["negative_wall_pairings"] == 0
        and face_neighbor["zero_wall_pairings"] == 1
        and face_neighbor["q_square_positive"]
        and opposite_record["negative_wall_pairings"] == 0
        and opposite_record["zero_wall_pairings"] == 0
        and opposite_record["q_square_positive"]
    )
    return {
        "verified": verified,
        "from_representative": state.index,
        "wall": list(label),
        "adjacency_kind": atlas.chamber(state.index).adjacencies[orbit].kind,
        "to_physical_representative": neighbor.index,
        "shared_face_in_D181": face_current,
        "shared_face_in_neighbor": face_neighbor,
        "opposite_side_strict_interior": opposite_record,
        "returns_to_original_physical_D68": (
            neighbor == atlas.representative_state(68)
        ),
    }


def analyze_d181(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")

    atlas = ShimadaDisc11Atlas(path)
    state = _physical_d181_state(atlas)
    chamber = atlas.chamber(181)
    rows = _restricted_rows(atlas, state)
    interior = _d181_interior_point()
    classes = _positive_proportional_classes(rows)

    expected_class_set = {frozenset(value) for value in D181_EXPECTED_PROPORTIONAL_CLASSES}
    computed_class_set = {frozenset(value) for value in classes}
    classes_verified = computed_class_set == expected_class_set

    root_orbits = {
        index
        for index, adjacency in enumerate(chamber.adjacencies)
        if adjacency.kind == "_minustwo"
    }
    root_classes_verified = all(
        all(label[0] in root_orbits for label in wall_class)
        for wall_class in D181_ROOT_FACET_CLASSES
    )
    internal_classes_verified = all(
        all(label[0] not in root_orbits for label in wall_class)
        for wall_class in D181_INTERNAL_FACET_CLASSES
    )

    face_records = []
    face_witnesses: dict[WallLabel, GoldenVector] = {}
    for wall_class in D181_FACET_CLASSES:
        representative = wall_class[0]
        witness = _lorentzian_face_witness(interior, rows[representative])
        face_witnesses[representative] = witness
        face_records.append(_face_record(wall_class, witness, rows))

    dependence_records = [
        _dependence_record(target, terms, rows)
        for target, terms in D181_REDUNDANCY_DEPENDENCES
    ]
    generic_halfspace_audit = audit_golden_halfspace_certificate(
        rows,
        facet_classes=D181_FACET_CLASSES,
        facet_witnesses=face_witnesses,
        redundancies=dict(D181_REDUNDANCY_DEPENDENCES),
        interior=interior,
        quadratic_gram=tuple(
            tuple(K3_TRACE_SCALAR * entry for entry in row)
            for row in field_gram()
        ),
    )
    redundant_class_set = computed_class_set - {
        frozenset(value) for value in D181_FACET_CLASSES
    }
    expected_redundant_set = {
        frozenset((target,)) for target, _ in D181_REDUNDANCY_DEPENDENCES
    }

    crossing_records = [
        _internal_crossing_record(
            atlas, state, wall_class[0], face_witnesses[wall_class[0]], interior
        )
        for wall_class in D181_INTERNAL_FACET_CLASSES
    ]

    interior_record = _point_record(
        atlas, state, interior, transported_identity_eigenbasis()
    )
    verified = (
        len(rows) == 26
        and classes_verified
        and root_classes_verified
        and internal_classes_verified
        and interior_record["negative_wall_pairings"] == 0
        and interior_record["zero_wall_pairings"] == 0
        and interior_record["q_square_positive"]
        and all(record["face_witness_timelike"] for record in face_records)
        and all(record["target_pairings_zero"] for record in face_records)
        and all(
            record["all_other_primitive_wall_pairings_positive"]
            for record in face_records
        )
        and redundant_class_set == expected_redundant_set
        and all(record["verified"] for record in dependence_records)
        and generic_halfspace_audit.passes
        and all(record["verified"] for record in crossing_records)
    )
    return {
        "verified": verified,
        "analysis_id": "shimada_disc11_physical_D181_restricted_section",
        "source_sha256": digest,
        "atlas_chambers": atlas.chamber_count,
        "physical_state": {
            "reached_from": 68,
            "crossed_D68_wall": [14, 0],
            "representative": state.index,
            "primitive_wall_count": chamber.wall_count,
            "wall_orbit_count": len(chamber.wall_orbits),
            "root_wall_orbits": sorted(root_orbits),
        },
        "strict_interior": interior_record,
        "positive_proportional_restrictions": {
            "class_count": len(classes),
            "classes": [
                _proportionality_record(wall_class, rows) for wall_class in classes
            ],
        },
        "irredundant_section": {
            "facet_count": len(D181_FACET_CLASSES),
            "root_facets": len(D181_ROOT_FACET_CLASSES),
            "internal_facets": len(D181_INTERNAL_FACET_CLASSES),
            "root_facet_classes": [
                [list(label) for label in wall_class]
                for wall_class in D181_ROOT_FACET_CLASSES
            ],
            "internal_facet_classes": [
                [list(label) for label in wall_class]
                for wall_class in D181_INTERNAL_FACET_CLASSES
            ],
            "positive_timelike_face_witnesses": face_records,
            "positive_redundancy_dependences": dependence_records,
            "generic_exact_halfspace_kernel_verified": generic_halfspace_audit.passes,
        },
        "one_star_farther": crossing_records,
        "scope": (
            "complete_local_D181_restriction_and_internal_neighbors_only_"
            "not_global_eigen_nef_completeness_volume_or_compactness"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze_d181(args.disc11), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
