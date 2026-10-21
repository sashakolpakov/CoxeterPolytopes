#!/usr/bin/env python3
"""Certify the mod-three obstruction for one ACT real ``H^10`` slice.

The exploratory ACT script fixes the antiunitary involution which reverses
the two ``EE8`` chains and swaps the final hyperbolic pair.  Its fixed module
has Gram matrix ``2*diag(1^10,-1)``.  This certificate turns the apparently
large Vinberg search into a six-dimensional ternary-code calculation.

No call to cddlib occurs here.  The complete obstruction is already visible
    in the residue code: every extendable primitive reflection root has norm one,
    and no two such roots are orthogonal.  Hence an acute simple-root diagram
    has neither an elliptic subsystem above rank one nor a parabolic subsystem
    above ``affine A1``.  Vinberg's first finite-volume gate in dimension ten
    requires an elliptic rank-ten vertex or a parabolic rank-nine cusp.
    Therefore this antiunitary branch cannot produce a finite-volume chamber,
    compact or otherwise.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import itertools
import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np

from scripts.explore_act_real_cubic_h10_reflections import (
    DIMENSION,
    RANK,
    _act_data,
    _extension_mask,
    _stream_gram_summary,
    _vinberg_levels,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "act_real_cubic_h10_congruence_obstruction.json"
)
SCHEMA = "coxeter.act_real_cubic_h10_congruence_obstruction.v1"

EXPECTED_CONSTRAINTS = [
    [0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1, 2, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2],
]


def _rref_mod_three(matrix: np.ndarray) -> np.ndarray:
    reduced = matrix.copy() % 3
    pivot_row = 0
    for column in range(reduced.shape[1]):
        pivot = next(
            (
                row
                for row in range(pivot_row, reduced.shape[0])
                if reduced[row, column]
            ),
            None,
        )
        if pivot is None:
            continue
        reduced[[pivot_row, pivot]] = reduced[[pivot, pivot_row]]
        reduced[pivot_row] = (
            reduced[pivot_row]
            * pow(int(reduced[pivot_row, column]), -1, 3)
        ) % 3
        for row in range(reduced.shape[0]):
            if row != pivot_row and reduced[row, column]:
                reduced[row] = (
                    reduced[row]
                    - reduced[row, column] * reduced[pivot_row]
                ) % 3
        pivot_row += 1
    return reduced[:pivot_row]


def _extension_constraints(data: dict[str, object]) -> np.ndarray:
    """Return independent rows ``a`` for the conditions ``a*r=0 mod 3``."""

    adjoint = data["adjoint_mod3"]
    if not isinstance(adjoint, np.ndarray):
        raise TypeError("the ACT adjoint congruence matrix is not an array")
    # _extension_mask uses r^T Q J.  The constant coefficient of J vanishes;
    # its epsilon coefficient, after applying Q, supplies the five equations.
    coefficient = adjoint[:, :, 1].copy()
    coefficient[-1, :] *= -1
    return _rref_mod_three(coefficient.T)


def _code_basis() -> np.ndarray:
    """Basis of F3 e0 + T4 + T4 + <e9+e10>."""

    rows = np.zeros((6, RANK), dtype=np.int16)
    rows[0, 0] = 1
    rows[1, 1:5] = (2, 2, 1, 0)
    rows[2, 1:5] = (2, 1, 0, 1)
    rows[3, 5:9] = (2, 2, 1, 0)
    rows[4, 5:9] = (2, 1, 0, 1)
    rows[5, 9:11] = (1, 1)
    return rows


def _tetracode() -> set[tuple[int, int, int, int]]:
    return {
        ((-a - b) % 3, (-a + b) % 3, a, b)
        for a, b in itertools.product(range(3), repeat=2)
    }


def _tetracode_monomial_automorphisms() -> list[
    tuple[tuple[int, ...], tuple[int, ...]]
]:
    code = _tetracode()
    automorphisms = []
    for permutation in itertools.permutations(range(4)):
        for signs in itertools.product((1, -1), repeat=4):
            image = {
                tuple(
                    (signs[index] * word[permutation[index]]) % 3
                    for index in range(4)
                )
                for word in code
            }
            if image == code:
                automorphisms.append((permutation, signs))
    return automorphisms


def _block_orbit(
    vector: Sequence[int],
    automorphisms: Sequence[tuple[tuple[int, ...], tuple[int, ...]]],
) -> set[tuple[int, ...]]:
    return {
        tuple(
            signs[index] * int(vector[permutation[index]])
            for index in range(4)
        )
        for permutation, signs in automorphisms
    }


def _wall_orbit_census(walls: np.ndarray) -> dict[str, Any]:
    automorphisms = _tetracode_monomial_automorphisms()

    def block_key(vector: Sequence[int]) -> tuple[int, ...]:
        return min(_block_orbit(vector, automorphisms))

    def wall_key(vector: np.ndarray) -> tuple[Any, ...]:
        first = block_key(vector[1:5])
        second = block_key(vector[5:9])
        return (
            int(vector[0]),
            min(first, second),
            max(first, second),
            int(vector[9]),
            int(vector[10]),
        )

    orbits: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, wall in enumerate(walls):
        orbits[wall_key(wall)].append(index)

    completion_failures = []
    for key, indices in orbits.items():
        wall = walls[indices[0]]
        first_orbit = _block_orbit(wall[1:5], automorphisms)
        second_orbit = _block_orbit(wall[5:9], automorphisms)
        same_block_orbit = block_key(wall[1:5]) == block_key(wall[5:9])
        predicted = len(first_orbit) * len(second_orbit) * (
            1 if same_block_orbit else 2
        )
        if predicted != len(indices):
            completion_failures.append(
                {"key": repr(key), "observed": len(indices), "predicted": predicted}
            )

    height_orbits = Counter(int(key[-1]) for key in orbits)
    size_histogram = Counter(len(indices) for indices in orbits.values())
    representatives = [
        {
            "height": int(key[-1]),
            "representative": list(map(int, walls[indices[0]])),
            "orbit_size": len(indices),
        }
        for key, indices in sorted(orbits.items(), key=lambda item: repr(item[0]))
    ]
    return {
        "tetracode_monomial_automorphism_order": len(automorphisms),
        "basepoint_and_oriented_initial_wall_stabilizer_order": 2 * len(automorphisms) ** 2,
        "stabilizer_structure": "(MonAut(T4) x MonAut(T4)) semidirect C2",
        "wall_orbit_count_through_height_three": len(orbits),
        "wall_orbits_by_height": {
            str(height): count for height, count in sorted(height_orbits.items())
        },
        "wall_orbit_size_histogram": {
            str(size): count for size, count in sorted(size_histogram.items())
        },
        "every_computed_orbit_is_complete": not completion_failures,
        "completion_failures": completion_failures,
        "orbit_representatives": representatives,
    }


def build_certificate() -> dict[str, Any]:
    data = _act_data()
    constraints = _extension_constraints(data)
    if constraints.tolist() != EXPECTED_CONSTRAINTS:
        raise ArithmeticError("the ACT extension constraint matrix changed")

    residues = np.array(
        list(itertools.product(range(3), repeat=RANK)), dtype=np.int16
    )
    exact_extension = _extension_mask(residues, data)
    code_extension = np.all((residues @ constraints.T) % 3 == 0, axis=1)
    if not np.array_equal(exact_extension, code_extension):
        raise ArithmeticError("the five linear equations are not the exact gate")

    codewords = residues[code_extension]
    q_mod_three = (
        np.sum(codewords[:, :DIMENSION] ** 2, axis=1)
        - codewords[:, -1] ** 2
    ) % 3
    if not np.array_equal(q_mod_three, codewords[:, 0] ** 2 % 3):
        raise ArithmeticError("Q restricted to the extension code is not x0^2")

    basis = _code_basis()
    if not np.all((basis @ constraints.T) % 3 == 0):
        raise ArithmeticError("the displayed code basis violates a constraint")
    lorentz = np.diag([1] * DIMENSION + [-1]).astype(np.int16)
    restricted_gram = (basis @ lorentz @ basis.T) % 3
    expected_restricted = np.zeros((6, 6), dtype=np.int16)
    expected_restricted[0, 0] = 1
    if not np.array_equal(restricted_gram, expected_restricted):
        raise ArithmeticError("the code bilinear form did not reduce to x0*y0")

    tetracode = _tetracode()
    tetracode_weights = Counter(sum(coordinate != 0 for coordinate in word) for word in tetracode)
    if tetracode_weights != {0: 1, 3: 8}:
        raise ArithmeticError("the length-four code is not the ternary tetracode")

    walls, raw_ledger = _vinberg_levels(data, max_height=3, chunk_size=512)
    ledger = [
        {key: value for key, value in row.items() if key != "seconds"}
        for row in raw_ledger
    ]
    gram_summary = _stream_gram_summary(walls, block_size=256)
    symmetry = _wall_orbit_census(walls)
    if not (
        len(walls) == 580
        and ledger[-1]["walls"] == 580
        and symmetry["wall_orbit_count_through_height_three"] == 10
        and symmetry["every_computed_orbit_is_complete"]
        and gram_summary["positive_offdiagonal"] == 0
        and gram_summary["zero_offdiagonal_pairs"] == 0
    ):
        raise ArithmeticError("the bounded Vinberg/orbit replay changed")

    accepted_q_histogram = Counter(map(int, q_mod_three))
    certificate = {
        "schema": SCHEMA,
        "status": "PROVED_GLOBAL_INFINITE_COVOLUME_OBSTRUCTION_FOR_THIS_ANTIUNITARY_CLASS",
        "source_lattice_and_convention": {
            "ACT_lattice": "(3) + EE8 + EE8 + [[0,theta],[bar(theta),0]]",
            "Hermitian_convention": "h(x,y)=x^T H conjugate(y), linear in x",
            "antiunitary_involution": "reverse each EE8 block, swap the final pair, then conjugate",
            "fixed_Z_lattice_Gram": "2*diag(1,1,1,1,1,1,1,1,1,1,-1)",
            "fixed_lattice_signature": [10, 1],
            "coordinate_order": [
                "x0",
                "T4a_0",
                "T4a_1",
                "T4a_2",
                "T4a_3",
                "T4b_0",
                "T4b_1",
                "T4b_2",
                "T4b_3",
                "x9",
                "t",
            ],
        },
        "exact_extension_congruence": {
            "ring": "E/(3)=F3[epsilon]/(epsilon^2)",
            "independent_linear_constraints_mod_3": constraints.tolist(),
            "equations": [
                "x1+x3+x4=0",
                "x2+x3-x4=0",
                "x5+x7+x8=0",
                "x6+x7-x8=0",
                "x9-t=0",
            ],
            "accepted_residue_count": int(np.sum(code_extension)),
            "ambient_residue_count": 3 ** RANK,
            "code_dimension": 6,
            "direct_sum_decomposition": "F3*e0 + T4 + T4 + <e9+t>",
            "code_basis_rows_mod_3": basis.tolist(),
            "tetracode_weight_histogram": {
                str(weight): count for weight, count in sorted(tetracode_weights.items())
            },
            "restricted_bilinear_Gram_mod_3": restricted_gram.tolist(),
            "quadratic_identity_on_every_codeword": "Q(r)=x0^2 mod 3",
            "bilinear_identity_on_every_pair": "B(r,s)=x0(r)*x0(s) mod 3",
            "accepted_Q_mod_3_histogram": {
                str(value): count
                for value, count in sorted(accepted_q_histogram.items())
            },
        },
        "primitive_reflection_root_gate": {
            "integrality_argument": "for primitive r, Q(r) divides every 2*r_i, hence Q(r) divides 2",
            "only_possible_positive_norms": [1, 2],
            "norm_two_excluded_by_code_identity": True,
            "every_extendable_root_has_norm": 1,
            "every_extendable_root_has_x0_nonzero_mod_3": True,
            "orthogonal_pair_count_in_the_complete_root_set": 0,
            "orthogonal_pair_proof": "B(r,s)=x0(r)*x0(s) is nonzero mod 3 for any two norm-one roots",
        },
        "safe_symmetry_compression": symmetry,
        "bounded_vinberg_corroboration": {
            "logically_required_for_global_theorem": False,
            "levels_through_height_three": ledger,
            "partial_wall_count": len(walls),
            "streamed_Gram_summary": gram_summary,
            "bounded_exact_cdd_observation_2026_08_09": {
                "wall_count": 580,
                "extreme_rays": 1680,
                "finite_timelike": 0,
                "ideal_lightlike": 0,
                "hyperideal_spacelike": 1680,
                "scope": "partial chamber only; corroboration, not the proof",
            },
        },
        "global_finite_volume_theorem": {
            "dimension": 10,
            "maximum_elliptic_rank": 1,
            "maximum_parabolic_Coxeter_rank": 1,
            "finite_vertex_gate": "requires an elliptic rank-10 subdiagram; none exists",
            "ideal_vertex_gate": "requires a componentwise-affine total rank-9 subdiagram; none exists",
            "Vinberg_first_finite_volume_gate_passes": False,
            "ridge_angle_gate": "an integral norm-one Coxeter ridge requires B(r,s)=0 (angle pi/2)",
            "ridge_existence_gate": "every finite-volume 10-polytope has codimension-two faces",
            "contradiction": "the complete extendable root set has no orthogonal pair",
            "compact_chamber_exists": False,
            "finite_volume_noncompact_chamber_exists": False,
            "new_hyperbolic_example": False,
            "Vinberg_search_status": "TERMINATED_BY_GLOBAL_CONGRUENCE_OBSTRUCTION",
        },
        "scope_and_nonclaims": {
            "proved_for_all_heights": True,
            "proved_for_this_explicit_antiunitary_class_only": True,
            "classification_of_all_ACT_antiunitary_classes": False,
            "real_cubic_threefold_component_realization_proved": False,
            "algebro_geometric_leverage_delta": 0,
            "large_cdd_needed": False,
            "height_four_cdd_permitted": False,
        },
        "all_checks_pass": True,
    }
    return certificate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    certificate = build_certificate()
    if arguments.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if certificate != archived:
            raise ArithmeticError("the frozen ACT congruence archive changed")
        print("archive replay: PASS")
    else:
        print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
