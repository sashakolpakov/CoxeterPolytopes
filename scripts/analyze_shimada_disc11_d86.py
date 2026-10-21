"""Exact reduction and one-star traversal of the physical D_86 eigensection.

The physical state is reached from the certified D_68 state by crossing
orbit 11, member 1.  This chamber-specific verifier consumes Shimada's
external ``disc11.txt`` supplement, restricts every primitive wall to the
fixed transported golden eigenspace, and verifies archived exact
``Q(phi)`` facet and redundancy certificates.

Nothing in this file asserts global eigen-nef completeness.  It certifies
one physical induced chamber and crosses each of its genuine internal
facets once.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from sympy import Matrix

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.geometric.golden_halfspaces import (
    audit_golden_halfspace_certificate,
)
from coxeter_polytopes.geometric.k3 import K3_TRACE_SCALAR, field_inner_product
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


D86_STATE_SHA256 = "9a1a2c3d5a57379afcfced29b458447aec6638188a185e6da722b8ceafeb3d45"
D86_PHYSICAL_WALL_ROWS_SHA256 = (
    "ff6ffd9f03b0f4311e4fa550cb001f4d12ee0bd21853b16808ada305dc90839c"
)
D86_RESTRICTED_ROWS_SHA256 = (
    "2228890afd06e69618ebb95667660e28ebabab53f200c4ffb6cbd6313c7b5a7d"
)


# Complete positive-proportional partition of the 35 physical walls.
D86_POSITIVE_PROPORTIONAL_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0), (13, 0)),
    ((1, 0), (4, 0)),
    ((2, 0), (16, 0)),
    ((3, 0), (14, 0)),
    ((5, 0),),
    ((6, 0), (12, 0)),
    ((7, 0), (10, 0)),
    ((8, 0), (9, 0)),
    ((11, 0),),
    ((15, 0),),
    ((17, 0), (34, 0)),
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
    ((32, 0),),
    ((33, 0),),
)


D86_ROOT_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((0, 0), (13, 0)),
    ((1, 0), (4, 0)),
    ((2, 0), (16, 0)),
    ((3, 0), (14, 0)),
    ((6, 0), (12, 0)),
    ((7, 0), (10, 0)),
    ((8, 0), (9, 0)),
)
D86_INTERNAL_FACET_CLASSES: tuple[tuple[tuple[int, int], ...], ...] = (
    ((20, 0),),
    ((23, 0),),
    ((24, 0),),
    ((26, 0),),
)
D86_FACET_CLASSES = D86_ROOT_FACET_CLASSES + D86_INTERNAL_FACET_CLASSES
D86_FACET_REPRESENTATIVES = tuple(wall_class[0] for wall_class in D86_FACET_CLASSES)


# Every witness is timelike, vanishes on exactly its target class, and is
# strictly positive on the other 25 or 26 restricted classes.
D86_FACET_WITNESSES: tuple[
    tuple[tuple[tuple[int, int], ...], tuple[Golden, ...]], ...
] = (
    (((0, 0), (13, 0)), (
        _g(17, -16), _g(-1), _g(-7, 3), _g(-6, 2), _g(-9, 3),
        _g(-8, 2), _g(-8, 1), _g(-8, 1), _g(-9, 1), _g(-10, 1),
    )),
    (((1, 0), (4, 0)), (
        _g(14, -15), _g(-1), _g(-10, 5, 3), _g(-2, 1, 3), _g(0),
        _g(-2, -2, 3), _g(-7, -4, 3), _g(-4, -1),
        _g(-20, -2, 3), _g(-6, -1),
    )),
    (((2, 0), (16, 0)), (
        _g(38, -28), _g(-1), _g(-5, 2), _g(-2), _g(-3),
        _g(-5, 1), _g(-8, 2), _g(-10, 3), _g(-14, 5), _g(-13, 4),
    )),
    (((3, 0), (14, 0)), (
        _g(30, -24), _g(-1), _g(-5, 2), _g(-2), _g(-3), _g(-4),
        _g(-7, 1), _g(-9, 2), _g(-12, 3), _g(-11, 2),
    )),
    (((6, 0), (12, 0)), (
        _g(14, -14), _g(-1), _g(-5, 2), _g(-2), _g(-3), _g(-4),
        _g(-4, -1), _g(-6), _g(-7), _g(-6, -1),
    )),
    (((7, 0), (10, 0)), (
        _g(45, -33), _g(-1), _g(-4, 2), _g(0), _g(0), _g(-3, 1),
        _g(-7, 2), _g(-11, 4), _g(-16, 6), _g(-14, 5),
    )),
    (((8, 0), (9, 0)), (
        _g(11, -14), _g(-1), _g(-3, 2), _g(0, 1), _g(1, 1),
        _g(0), _g(-2, -1), _g(-4, -1), _g(-7, -1), _g(-4, -2),
    )),
    (((20, 0),), (
        _g(-37, 15), _g(27, -19), _g(-2, 0, 3), _g(5, -3, 3),
        _g(0), _g(-1, -3, 3), _g(-5, -6, 3), _g(-3, -2),
        _g(-16, -6, 3), _g(-5, -2),
    )),
    (((23, 0),), (
        _g(-92, 50), _g(71, -46), _g(-5, 2), _g(-2), _g(-3),
        _g(-4), _g(-6), _g(-8, 1), _g(-11, 2), _g(-10, 1),
    )),
    (((24, 0),), (
        _g(24, -21), _g(-1), _g(-4, 1), _g(-1, -1), _g(-2, -1),
        _g(-3, -1), _g(-5, -1), _g(-7), _g(-10, 1), _g(-9),
    )),
    (((26, 0),), (
        _g(8, -12), _g(-1), _g(-2, 0, 3), _g(5, -3, 3), _g(0),
        _g(-1, -3, 3), _g(-5, -6, 3), _g(-3, -2),
        _g(-16, -6, 3), _g(-5, -2),
    )),
)


F0, F1, F2, F3, F5, F6, F7, F13, F16, F17, F19 = D86_FACET_REPRESENTATIVES


# One positive exact dependence for every nonfacet proportionality class.
D86_REDUNDANCY_DEPENDENCES: tuple[
    tuple[tuple[int, int], tuple[tuple[Golden, tuple[int, int]], ...]], ...
] = (
    ((5, 0), (
        (PHI, F0), (2 * PHI, F1), (2 + 2 * PHI, F2), (PHI, F5),
        (2 + 2 * PHI, F6), (2 * PHI, F7), (2 * PHI, F13),
        (PHI, F17), (PHI, F19),
    )),
    ((11, 0), ((PHI - 1, F13), (PHI - 1, F19))),
    ((15, 0), (
        (2 * PHI, F1), (PHI, F2), ((PHI - 1) / 3, F3),
        ((4 - 2 * PHI) / 3, F5), ((2 + 4 * PHI) / 3, F6),
        ((1 + 4 * PHI) / 3, F7), ((-1 + 8 * PHI) / 3, F13),
        (_g(4, 6, 3), F16), (_g(6, 10, 3), F17),
        (_g(4, 8, 3), F19),
    )),
    ((17, 0), ((2 - PHI, F17), (2 - PHI, F19))),
    ((18, 0), ((Golden(1), F16), (PHI, F17), (Golden(1), F19))),
    ((19, 0), ((PHI - 1, F16), (PHI - 1, F17), (PHI - 1, F19))),
    ((21, 0), (
        (Golden(1), F0), (1 + PHI, F1), (2 * PHI, F2),
        (Golden(1), F5), (2 * PHI, F6), (Golden(2), F7),
        (2 * PHI, F13), (1 + 2 * PHI, F16), (2 + 3 * PHI, F17),
        (2 + 2 * PHI, F19),
    )),
    ((22, 0), (
        (1 - PHI / 2, F0), (_g(2, 3, 2), F1), (PHI, F2),
        (_g(1, 3, 2), F6), (1 + PHI, F7), (1 + 2 * PHI, F13),
        (_g(1, 4, 2), F16), (_g(4, 5, 2), F17),
        (2 + 2 * PHI, F19),
    )),
    ((25, 0), (
        (-2 + 2 * PHI, F1), (Golden(1), F2), (_g(1, 0, 3), F3),
        (_g(-2, 2, 3), F5), (_g(4, 0, 3), F6),
        (_g(-5, 5, 3), F7), (_g(-4, 7, 3), F13),
        (_g(3, 4, 3), F16), (_g(1, 10, 3), F17),
        (_g(-1, 9, 3), F19),
    )),
    ((27, 0), (
        (Golden(2), F1), (PHI, F2), (PHI / 3, F3),
        (_g(2, 0, 3), F5), (4 * PHI / 3, F6), (_g(5, 0, 3), F7),
        (_g(4, 3, 3), F13), (_g(1, 4, 3), F16),
        (_g(1, 8, 3), F17), (_g(3, 5, 3), F19),
    )),
    ((28, 0), (
        (PHI - 1, F0), (2 * PHI - 2, F1), (Golden(2), F2),
        (PHI - 1, F5), (Golden(2), F6), (2 * PHI - 2, F7),
        (2 * PHI - 2, F13), (1 + PHI, F16), (1 + 2 * PHI, F17),
        (2 * PHI, F19),
    )),
    ((29, 0), (
        (_g(1, 0, 2), F0), (_g(1, 4, 2), F1), (PHI, F2),
        (3 * PHI / 2, F6), (1 + PHI, F7), (3 * PHI, F13),
        (_g(3, 5, 2), F16), (_g(3, 10, 2), F17),
        (1 + 4 * PHI, F19),
    )),
    ((30, 0), (
        (2 * PHI, F1), (1 + PHI, F2), ((1 + PHI) / 3, F3),
        (2 * PHI / 3, F5), (_g(4, 4, 3), F6), (5 * PHI / 3, F7),
        (7 * PHI / 3, F13), (_g(4, 2, 3), F16),
        (_g(5, 6, 3), F17), (_g(5, 5, 3), F19),
    )),
    ((31, 0), (
        (_g(3, -1, 2), F0), (_g(1, 3, 2), F1), (2 * PHI, F2),
        (Golden(1), F5), (_g(1, 4, 2), F6), (1 + PHI, F7),
        (2 * PHI, F13), (_g(2, 3, 2), F16), (_g(3, 5, 2), F17),
        (1 + 2 * PHI, F19),
    )),
    ((32, 0), (
        ((PHI - 1) / 2, F0), (_g(-3, 7, 2), F1), (Golden(2), F2),
        (PHI - 1, F5), (_g(5, 0, 2), F6), (-2 + 3 * PHI, F7),
        (-1 + 4 * PHI, F13), (_g(4, 5, 2), F16),
        (_g(5, 9, 2), F17), (2 + 4 * PHI, F19),
    )),
    ((33, 0), (
        (_g(-1, 2, 2), F0), (_g(1, 4, 2), F1), (2 + PHI, F2),
        (PHI - 1, F5), (_g(4, 3, 2), F6), (2 * PHI, F7),
        (3 * PHI, F13), (_g(3, 5, 2), F16), (_g(5, 8, 2), F17),
        (2 + 3 * PHI, F19),
    )),
)


D86_INTERNAL_CROSSINGS = (
    (20, 55, 16),
    (23, 30, 4),
    (24, 27, 2),
    (26, 68, 2),
)


def _pair(row: tuple[Golden, ...], point: tuple[Golden, ...]) -> Golden:
    return sum((row[index] * point[index] for index in range(10)), Golden(0))


def _canonical_sha(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _positive_proportional_classes(coefficients):
    classes: list[list[tuple[int, int]]] = []
    ratios: dict[tuple[int, int], Golden] = {}
    for label, row in coefficients.items():
        for wall_class in classes:
            base = coefficients[wall_class[0]]
            pivot = next(index for index, entry in enumerate(base) if entry)
            ratio = row[pivot] / base[pivot]
            if ratio.sign() > 0 and row == tuple(ratio * entry for entry in base):
                wall_class.append(label)
                ratios[label] = ratio
                break
        else:
            classes.append([label])
            ratios[label] = Golden(1)
    return tuple(tuple(wall_class) for wall_class in classes), ratios


def analyze_d86(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")

    atlas = ShimadaDisc11Atlas(path)
    basis = transported_identity_eigenbasis()
    d68 = atlas.representative_state(68)
    d86 = atlas.neighboring_state(d68, 11, 1)
    if d86 is None or d86.index != 86:
        raise ValueError("D_68 orbit 11 member 1 did not reach D_86")

    physical_orbits = atlas.state_wall_orbits(d86)
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
    classes, ratios = _positive_proportional_classes(coefficients)

    witness_records = []
    witnesses_verified = True
    for target_class, witness in D86_FACET_WITNESSES:
        pairings = {label: _pair(row, witness) for label, row in coefficients.items()}
        target = set(target_class)
        target_zero = all(pairings[label] == 0 for label in target_class)
        others_positive = all(
            value.sign() > 0 for label, value in pairings.items() if label not in target
        )
        q_square = K3_TRACE_SCALAR * field_inner_product(
            field_gram(), witness, witness
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

    dependence_records = []
    dependences_verified = True
    for target, terms in D86_REDUNDANCY_DEPENDENCES:
        total = tuple(
            sum(
                (scalar * coefficients[facet][index] for scalar, facet in terms),
                Golden(0),
            )
            for index in range(10)
        )
        verified = (
            all(scalar.sign() > 0 for scalar, _ in terms)
            and total == coefficients[target]
        )
        dependences_verified &= verified
        dependence_records.append(
            {
                "target": list(target),
                "support": [list(facet) for _, facet in terms],
                "all_coefficients_positive": all(
                    scalar.sign() > 0 for scalar, _ in terms
                ),
                "identity_verified": total == coefficients[target],
            }
        )

    chamber = atlas.chamber(86)
    root_orbits = [
        index for index, adjacency in enumerate(chamber.adjacencies)
        if adjacency.kind == "_minustwo"
    ]
    gram = Matrix(expected_shimada_gram())
    inverse_gram = gram.inv()
    root_labels = []
    roots = []
    for orbit_index in root_orbits:
        for member_index, wall in enumerate(physical_orbits[orbit_index]):
            primal = Matrix(1, 20, wall) * inverse_gram
            if any(not entry.is_Integer for entry in primal):
                raise ValueError("a D_86 root wall recovered a nonintegral root")
            root = tuple(int(entry) for entry in primal)
            if (Matrix(1, 20, root) * gram * Matrix(20, 1, root))[0] != -2:
                raise ValueError("a D_86 root wall did not have square -2")
            root_labels.append((orbit_index, member_index))
            roots.append(root)
    root_filter = analyze_conjugated_roots(tuple(roots))
    surviving_root_labels = [
        root_labels[index] for index in root_filter["surviving_zero_based_indices"]
    ]

    d86_interior = _crossing_point(D68_INTERNAL_FACET_WITNESSES[1][1], 2)
    generic_halfspace_audit = audit_golden_halfspace_certificate(
        coefficients,
        facet_classes=D86_FACET_CLASSES,
        facet_witnesses={
            wall_class[0]: witness
            for wall_class, witness in D86_FACET_WITNESSES
        },
        redundancies=dict(D86_REDUNDANCY_DEPENDENCES),
        interior=d86_interior,
        quadratic_gram=tuple(
            tuple(K3_TRACE_SCALAR * entry for entry in row)
            for row in field_gram()
        ),
    )
    witness_by_orbit = {
        wall_class[0][0]: witness
        for wall_class, witness in D86_FACET_WITNESSES
        if wall_class in D86_INTERNAL_FACET_CLASSES
    }
    crossing_records = []
    crossings_verified = True
    for orbit, expected_target, denominator in D86_INTERNAL_CROSSINGS:
        neighbor = atlas.neighboring_state(d86, orbit, 0)
        if neighbor is None:
            raise ValueError("a D_86 internal facet became a K3 boundary")
        face = witness_by_orbit[orbit]
        epsilon = Golden(1) / denominator
        crossing_point = tuple(
            face[index] - epsilon * (d86_interior[index] - face[index])
            for index in range(10)
        )
        record = _point_record(atlas, neighbor, crossing_point, basis)
        verified = (
            neighbor.index == expected_target
            and record["q_square_positive"]
            and record["negative_wall_pairings"] == 0
            and record["zero_wall_pairings"] == 0
        )
        crossings_verified &= verified
        crossing_records.append(
            {
                "D86_orbit": orbit,
                "target_representative": neighbor.index,
                "crossing_denominator": denominator,
                "target_interior": record,
                "verified": verified,
            }
        )

    # The fourth internal facet is exactly the wall by which the walk entered
    # D_86, with the opposite orientation.
    d68_wall = _restricted_coefficients(
        atlas.state_wall_orbits(d68)[11][1], basis
    )
    return_wall = coefficients[(26, 0)]
    return_wall_verified = return_wall == tuple(-entry for entry in d68_wall)

    facet_class_set = {frozenset(item) for item in D86_FACET_CLASSES}
    computed_class_set = {frozenset(item) for item in classes}
    redundant_class_set = computed_class_set - facet_class_set
    expected_redundant = {
        frozenset(item)
        for item in D86_POSITIVE_PROPORTIONAL_CLASSES
        if frozenset(item) not in facet_class_set
    }
    partition_verified = (
        classes == D86_POSITIVE_PROPORTIONAL_CLASSES
        and len(coefficients) == 35
        and len(classes) == 27
        and len(D86_FACET_CLASSES) == 11
        and redundant_class_set == expected_redundant
        and len(D86_REDUNDANCY_DEPENDENCES) == len(redundant_class_set) == 16
    )
    fingerprints_verified = (
        _canonical_sha(d86.row_isometry) == D86_STATE_SHA256
        and _canonical_sha(physical_rows) == D86_PHYSICAL_WALL_ROWS_SHA256
        and _canonical_sha(restricted_rows) == D86_RESTRICTED_ROWS_SHA256
    )
    verified = (
        partition_verified
        and fingerprints_verified
        and witnesses_verified
        and dependences_verified
        and generic_halfspace_audit.passes
        and len(root_orbits) == 17
        and root_filter["trace_minus_two_verified"] == 17
        and root_filter["field_integral_reflection_roots"] == 16
        and root_filter["field_span_rank"] == 9
        and (15, 0) not in surviving_root_labels
        and crossings_verified
        and return_wall_verified
    )

    return {
        "verified": verified,
        "source_sha256": digest,
        "physical_state": {
            "entered_from": {"representative": 68, "orbit": 11, "member": 1},
            "representative": 86,
            "state_sha256": D86_STATE_SHA256,
            "physical_wall_rows_sha256": D86_PHYSICAL_WALL_ROWS_SHA256,
            "restricted_rows_sha256": D86_RESTRICTED_ROWS_SHA256,
            "fingerprints_verified": fingerprints_verified,
        },
        "restricted_section": {
            "ambient_primitive_walls": len(coefficients),
            "positive_proportional_classes": len(classes),
            "irredundant_facets": len(D86_FACET_CLASSES),
            "root_facets": len(D86_ROOT_FACET_CLASSES),
            "internal_facets": len(D86_INTERNAL_FACET_CLASSES),
            "root_facet_classes": [
                [list(label) for label in wall_class]
                for wall_class in D86_ROOT_FACET_CLASSES
            ],
            "internal_facet_classes": [
                [list(label) for label in wall_class]
                for wall_class in D86_INTERNAL_FACET_CLASSES
            ],
            "positive_proportionality_ratios": {
                f"{label[0]},{label[1]}": ratio.to_json()
                for label, ratio in ratios.items()
                if ratio != 1
            },
            "facet_witnesses": witness_records,
            "redundancy_dependences": dependence_records,
            "partition_verified": partition_verified,
            "witnesses_verified": witnesses_verified,
            "dependences_verified": dependences_verified,
            "generic_exact_halfspace_kernel_verified": generic_halfspace_audit.passes,
        },
        "root_filter": {
            "ambient_minus_two_walls": len(root_labels),
            "positive_field_roots": root_filter["positive_at_both_places"],
            "integral_field_roots": root_filter["field_integral_reflection_roots"],
            "shell_distribution": root_filter["field_integral_shell_distribution"],
            "field_span_rank": root_filter["field_span_rank"],
            "pair_sign_distribution": root_filter["pair_sign_distribution"],
            "non_field_root_walls": [
                list(label) for label in root_labels if label not in surviving_root_labels
            ],
            "all_ample_pairings_positive": root_filter["all_ample_pairings_positive"],
        },
        "internal_first_star": {
            "crossings": crossing_records,
            "all_crossings_verified": crossings_verified,
            "orbit_26_is_opposite_D68_entry_wall": return_wall_verified,
            "representative_targets": [record["target_representative"] for record in crossing_records],
        },
        "scope": "complete_local_D86_section_and_internal_first_star_not_global_eigen_nef_completeness",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze_d86(args.disc11), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
