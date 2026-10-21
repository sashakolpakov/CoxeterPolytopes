"""Exact first-star traversal of the discriminant-11 K3 eigen-nef section.

This script consumes Shimada's external ``disc11.txt`` atlas.  Starting from
the certified ample eigensection in chamber D_68, it crosses the three
internal facets of ``D_68 cap E`` and proves that the neighboring physical
induced chambers have timelike interior points in the fixed eigenspace.

The external source is not vendored.  Its expected SHA-256 is recorded in
``scripts.shimada_disc11_conjugate``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

from sympy import Matrix

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.exact.matrix import rank
from coxeter_polytopes.geometric.k3 import K3_TRACE_SCALAR, field_inner_product
from coxeter_polytopes.geometric.k3_presentation import _inverse_golden
from scripts.coupled_discriminant11_candidate import field_gram
from scripts.shimada_disc11_atlas import ShimadaChamberState, ShimadaDisc11Atlas
from scripts.shimada_disc11_conjugate import (
    AMPLE_FIELD_COEFFICIENTS,
    D68_INTERNAL_FACET_WITNESSES,
    DISC11_SOURCE_SHA256,
    conjugated_phi_action_shimada,
    transported_identity_eigenbasis,
    verify_d68_restricted_facet_certificate,
)


def _restricted_coefficients(wall, basis):
    return tuple(
        sum((basis[k][i] * wall[i] for i in range(20)), Golden(0))
        for k in range(10)
    )


def _pair(coefficient, point):
    return sum(
        (coefficient[index] * point[index] for index in range(10)),
        Golden(0),
    )


def _crossing_point(face, denominator: int):
    scalar = Golden(1) / denominator
    return tuple(
        face[index]
        - scalar * (AMPLE_FIELD_COEFFICIENTS[index] - face[index])
        for index in range(10)
    )


def _state_pairings(atlas, state, point, basis):
    return tuple(
        _pair(_restricted_coefficients(wall, basis), point)
        for orbit in atlas.state_wall_orbits(state)
        for wall in orbit
    )


def _point_record(atlas, state, point, basis):
    pairings = _state_pairings(atlas, state, point, basis)
    square = K3_TRACE_SCALAR * field_inner_product(field_gram(), point, point)
    return {
        "representative_chamber": state.index,
        "primitive_wall_count": len(pairings),
        "q_square": square.to_json(),
        "q_square_positive": square.sign() > 0,
        "negative_wall_pairings": sum(value.sign() < 0 for value in pairings),
        "zero_wall_pairings": sum(not value for value in pairings),
        "positive_wall_pairings": sum(value.sign() > 0 for value in pairings),
        "minimum_wall_pairing": min(pairings).to_json(),
    }


def _induced_eigen_action(row_isometry, basis):
    transported = [
        [
            sum(
                (basis[i][k] * row_isometry[k][j] for k in range(20)),
                Golden(0),
            )
            for j in range(20)
        ]
        for i in range(10)
    ]
    columns = next(
        subset
        for subset in combinations(range(20), 10)
        if rank([[basis[i][j] for j in subset] for i in range(10)]) == 10
    )
    square_basis = [[basis[i][j] for j in columns] for i in range(10)]
    inverse = _inverse_golden(square_basis)
    action = [
        [
            sum(
                (
                    transported[i][columns[k]] * inverse[k][j]
                    for k in range(10)
                ),
                Golden(0),
            )
            for j in range(10)
        ]
        for i in range(10)
    ]
    if any(
        sum((action[i][k] * basis[k][j] for k in range(10)), Golden(0))
        != transported[i][j]
        for i in range(10)
        for j in range(20)
    ):
        raise ValueError("the chamber transport did not preserve the eigenspace")
    return action


def analyze_first_star(path: Path) -> dict[str, object]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != DISC11_SOURCE_SHA256:
        raise ValueError("unexpected disc11.txt checksum")

    atlas = ShimadaDisc11Atlas(path)
    basis = transported_identity_eigenbasis()
    start = atlas.representative_state(68)

    # The two simple crossings enter physical copies of D_86 and D_181.
    face_11 = D68_INTERNAL_FACET_WITNESSES[1][1]
    state_11 = atlas.neighboring_state(start, 11, 1)
    if state_11 is None:
        raise ValueError("D_68 orbit 11 unexpectedly became a K3 boundary")
    point_11 = _crossing_point(face_11, 2)

    face_14 = D68_INTERNAL_FACET_WITNESSES[2][1]
    state_14 = atlas.neighboring_state(start, 14, 0)
    if state_14 is None:
        raise ValueError("D_68 orbit 14 unexpectedly became a K3 boundary")
    point_14 = _crossing_point(face_14, 10)

    # The first internal facet is a codimension-four Conway crossing.  The
    # fixed eigenspace passes directly from D_68 to another physical copy of
    # D_68; the three intermediate chambers meet E only along the face.
    face_9 = D68_INTERNAL_FACET_WITNESSES[0][1]
    path_specification = ((9, 0), (12, 0), (8, 1), (15, 0))
    path_states = [start]
    current = start
    for orbit, member in path_specification:
        neighbor = atlas.neighboring_state(current, orbit, member)
        if neighbor is None:
            raise ValueError("the simultaneous internal path hit a K3 boundary")
        path_states.append(neighbor)
        current = neighbor
    final_state = current
    point_9 = _crossing_point(face_9, 2)

    final_matrix = Matrix(final_state.row_isometry)
    gram = Matrix(atlas.gram)
    tau = Matrix(conjugated_phi_action_shimada())
    eigen_action = _induced_eigen_action(final_state.row_isometry, basis)
    eigen_identity = [
        [eigen_action[i][j] - int(i == j) for j in range(10)]
        for i in range(10)
    ]
    eigen_negative_identity = [
        [eigen_action[i][j] + int(i == j) for j in range(10)]
        for i in range(10)
    ]
    wall_9 = _restricted_coefficients(
        atlas.chamber(68).wall_orbits[9][0], basis
    )
    wall_image = tuple(
        sum((eigen_action[i][j] * wall_9[j] for j in range(10)), Golden(0))
        for i in range(10)
    )

    state_11_record = _point_record(atlas, state_11, point_11, basis)
    state_14_record = _point_record(atlas, state_14, point_14, basis)
    final_record = _point_record(atlas, final_state, point_9, basis)
    face_records = [
        _point_record(atlas, state, face_9, basis) for state in path_states
    ]
    verified = (
        verify_d68_restricted_facet_certificate()["verified"]
        and state_11_record["negative_wall_pairings"] == 0
        and state_11_record["zero_wall_pairings"] == 0
        and state_11_record["q_square_positive"]
        and state_14_record["negative_wall_pairings"] == 0
        and state_14_record["zero_wall_pairings"] == 0
        and state_14_record["q_square_positive"]
        and final_state.index == 68
        and final_record["negative_wall_pairings"] == 0
        and final_record["zero_wall_pairings"] == 0
        and final_record["q_square_positive"]
        and all(record["negative_wall_pairings"] == 0 for record in face_records)
        and all(record["q_square_positive"] for record in face_records)
        and final_matrix.det() == 1
        and final_matrix * gram * final_matrix.T == gram
        and final_matrix * final_matrix == Matrix.eye(20)
        and final_matrix.T * tau == tau * final_matrix.T
        and rank(eigen_identity) == 1
        and rank(eigen_negative_identity) == 9
        and wall_image == tuple(-entry for entry in wall_9)
    )
    return {
        "verified": verified,
        "source_sha256": digest,
        "atlas_chambers": atlas.chamber_count,
        "local_D68_facets": {
            "root": 8,
            "internal": 3,
            "internal_wall_classes": [
                [[9, 0], [16, 0]],
                [[11, 1]],
                [[14, 0]],
            ],
        },
        "simple_internal_crossings": {
            "D68_wall_11_member_1_to_physical_D86": state_11_record,
            "D68_wall_14_to_physical_D181": state_14_record,
        },
        "simultaneous_internal_crossing": {
            "path": [
                {
                    "from_representative": path_states[index].index,
                    "orbit": orbit,
                    "member": member,
                    "to_representative": path_states[index + 1].index,
                }
                for index, (orbit, member) in enumerate(path_specification)
            ],
            "representative_sequence": [state.index for state in path_states],
            "face_closure_records": face_records,
            "opposite_side_interior": final_record,
            "returns_to_physical_copy_of_D68": final_state.index == 68,
            "transport_determinant": str(final_matrix.det()),
            "transport_is_involution": final_matrix * final_matrix == Matrix.eye(20),
            "transport_is_NS_isometry": final_matrix * gram * final_matrix.T == gram,
            "transport_centralizes_tau": final_matrix.T * tau == tau * final_matrix.T,
            "restricted_fixed_rank": 9,
            "restricted_negative_rank": 1,
            "restricted_action_negates_wall_normal": wall_image
            == tuple(-entry for entry in wall_9),
            "conclusion": "enhanced_K3_automorphism_reflects_across_internal_eigen_wall",
        },
        "scope": "first_star_traversal_not_global_eigen_nef_wall_completeness",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("disc11", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze_first_star(args.disc11), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
