"""Exact resumable screen of the 575 fresh weighted-wall fixed sections.

The complete weighted-wall symmetry inventory found 575 chambers whose full
exact weighted-wall graph group strictly extends Shimada's archived
``Aut(D)``.  Each contributes one proper involution section above H9: one
section in H13 and 574 in H18.  This script performs the missing exact
facet/angle gate, prioritizing the unique H13 case and then increasing wall
count.

The lifted graph actions are rational Lorentz isometries and need not preserve
the integral Neron--Severi lattice.  Accordingly, the arithmetic objects are
kept separate:

* the *saturated fixed divisor lattice* is ``NS(X) intersect Fix_Q(H)``, found
  as an integral Smith kernel after clearing rational denominators;
* the restricted wall mirrors live in its rational Lorentz space; and
* reflection stability of the saturated divisor lattice is tested mirror by
  mirror and is never assumed merely from weighted-diagram symmetry.

For each case every published wall is folded, all conic redundancies are
replayed over Q using only the final facets, every retained hyperbolic facet
gets a strict exact timelike witness, and every normalized pair product is
checked exactly.  Numerical optimization may discover a witness/support, but
only exact rational replay is archived.  A discovery failure is unresolved,
never a negative conclusion.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any

import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_shimada_disc16_d5_full_subgroup_fishing import (
    _dependence_on_retained_facets,
    _pair_gate,
)
from scripts.certify_shimada_disc16_d5_h14_fixed_lattice import (
    _smith_decomposition,
)
from scripts.inventory_shimada_disc16_weighted_wall_symmetries import (
    EXPECTED_SOURCE_SHA256,
    OUTPUT as INVENTORY,
    SOURCE,
)
from scripts.screen_shimada_disc16_atlas_archived_high_sections import (
    _audit_witness,
    _fixed_coordinates_of_ample,
    _primitive,
    _strict_timelike_witness,
)
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _exact_nonnegative_dependence,
    _matrix_payload,
    _orbits,
    _rational_text,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas, _row_product
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CASE_DIRECTORY = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_fresh_weighted_high_sections"
)
MANIFEST = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_fresh_weighted_high_section_census.json"
)
SCHEMA = "coxeter.shimada_disc16_fresh_weighted_high_section_case.v1"
MANIFEST_SCHEMA = "coxeter.shimada_disc16_fresh_weighted_high_section_census.v1"
AMBIENT_RANK = 20
MIN_FIXED_RANK = 11
EXPECTED_CASE_COUNT = 575


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _case_id(chamber_index: int, class_index: int) -> str:
    return f"W{chamber_index:04d}_C{class_index:03d}"


def _case_path(directory: Path, chamber_index: int, class_index: int) -> Path:
    return directory / f"{_case_id(chamber_index, class_index)}.json"


def _matrix(payload: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _saturated_rational_fixed_divisor_lattice(
    primal_actions: tuple[sp.Matrix, ...],
    ambient_gram: sp.Matrix,
    expected_rank: int,
) -> dict[str, Any]:
    """Return a basis of Z^20 intersecting the common rational fixed space."""

    rational_constraints = sp.Matrix.vstack(
        *(matrix.T - sp.eye(AMBIENT_RANK) for matrix in primal_actions)
    )
    denominator = int(
        sp.ilcm(*(sp.Rational(value).q for value in rational_constraints))
    )
    integral_constraints = denominator * rational_constraints
    if any(not value.is_Integer for value in integral_constraints):
        raise ArithmeticError("failed to clear rational fixed-space denominators")
    diagonal, _, right = _smith_decomposition(integral_constraints)
    zero_columns = tuple(
        column
        for column in range(AMBIENT_RANK)
        if all(diagonal[row, column] == 0 for row in range(diagonal.rows))
    )
    if len(zero_columns) != expected_rank:
        raise ArithmeticError("the saturated divisor intersection has wrong rank")
    smith_basis = right[:, list(zero_columns)].T
    basis = smith_basis.lll(delta=sp.Rational(3, 4))
    if not (
        basis.rank() == expected_rank
        and basis * rational_constraints.T
        == sp.zeros(expected_rank, rational_constraints.rows)
        and all(basis * matrix == basis for matrix in primal_actions)
    ):
        raise ArithmeticError("the saturated rational fixed kernel failed replay")
    _, pivot_columns = smith_basis.rref()
    change = basis[:, list(pivot_columns)] * smith_basis[:, list(pivot_columns)].inv()
    if not (change * smith_basis == basis and abs(int(change.det())) == 1):
        raise ArithmeticError("LLL changed the saturated fixed divisor lattice")

    fixed_gram = basis * ambient_gram * basis.T
    if exact_inertia(fixed_gram) != (1, expected_rank - 1, 0):
        raise ArithmeticError("the fixed divisor intersection is not Lorentzian")
    smith = smith_normal_form(fixed_gram, domain=ZZ)
    smith_diagonal = [
        abs(int(smith[index, index])) for index in range(expected_rank)
    ]
    return {
        "rational_constraints": rational_constraints,
        "integral_constraints": integral_constraints,
        "clearing_denominator": denominator,
        "smith_constraint_diagonal": [
            int(diagonal[index, index])
            for index in range(min(diagonal.rows, diagonal.cols))
        ],
        "basis": basis,
        "gram": fixed_gram,
        "determinant": int(fixed_gram.det()),
        "smith_diagonal_absolute": smith_diagonal,
    }


def _weighted_actions(
    chamber_record: dict[str, Any],
) -> dict[Permutation, sp.Matrix]:
    actions: dict[Permutation, sp.Matrix] = {}
    for record in chamber_record["weighted_group_elements"]:
        permutation = tuple(int(value) for value in record["permutation"])
        action = _matrix(record["dual_action"])
        if permutation in actions:
            raise ArithmeticError("a weighted group permutation is duplicated")
        actions[permutation] = action
    if len(actions) != chamber_record["weighted_wall_group_order"]:
        raise ArithmeticError("the weighted group action list is incomplete")
    return actions


def _zero_weyl_wall_orientation(
    atlas: ShimadaDisc11Atlas,
    chamber_index: int,
    orbit_index: int,
    raw_wall: tuple[int, ...],
    dual_gram: sp.Matrix,
) -> tuple[int, str]:
    """Orient a wall on which the projected published Weyl vector vanishes."""

    adjacency = atlas.chamber(chamber_index).adjacencies[orbit_index]
    if adjacency.kind == "_minustwo":
        # Shimada's effective (-2)-wall convention points into the nef side.
        return 1, "effective_minustwo_orientation_convention"
    if adjacency.kind == "_isom":
        adjacent_weyl = atlas.assignment(
            atlas._chamber_name(chamber_index, "adjacentweyl", orbit_index)
        )
    elif adjacency.kind in {"_backto", "_newcham"} and adjacency.target is not None:
        # These adjacency records use the identity local transport.
        adjacent_weyl = atlas.chamber(adjacency.target).weyl
    else:
        raise ArithmeticError("a zero-Weyl wall has no usable adjacent chamber")
    adjacent_dual = sp.Matrix(
        1,
        AMBIENT_RANK,
        _row_product(adjacent_weyl, atlas.projection_to_s_dual),
    )
    adjacent_pairing = sp.Rational(
        (adjacent_dual * dual_gram * sp.Matrix(raw_wall))[0]
    )
    if adjacent_pairing == 0:
        raise ArithmeticError("both chamber Weyl vectors lie on a recorded wall")
    return (
        -1 if adjacent_pairing > 0 else 1,
        f"opposite_adjacent_{adjacency.kind}_Weyl_side",
    )


def _strict_interior_from_oriented_walls(
    published_ample_dual: sp.Matrix,
    ambient_gram: sp.Matrix,
    oriented_wall_matrix: sp.Matrix,
    zero_positions: tuple[int, ...],
) -> tuple[sp.Matrix, dict[str, Any]]:
    """Perturb a timelike published Weyl vector off its zero walls exactly."""

    dual_gram = ambient_gram.inv()
    if not zero_positions:
        return published_ample_dual, {
            "method": "published_projected_Weyl_vector",
            "original_zero_wall_positions": [],
        }
    zero_walls = oriented_wall_matrix[list(zero_positions), :]
    solution = sp.linsolve(
        (zero_walls, sp.ones(len(zero_positions), 1))
    )
    if solution is sp.EmptySet:
        raise ArithmeticError("zero-Weyl walls admit no common inward perturbation")
    raw = tuple(next(iter(solution)))
    parameters = sorted(
        set().union(*(value.free_symbols for value in raw)), key=str
    )
    direction_primal = sp.Matrix(
        [sp.Rational(value.subs({parameter: 0 for parameter in parameters})) for value in raw]
    )
    if zero_walls * direction_primal != sp.ones(len(zero_positions), 1):
        raise ArithmeticError("the exact zero-wall perturbation failed replay")
    published_primal = published_ample_dual * dual_gram
    for exponent in range(1, 121):
        epsilon = sp.Rational(1, 2**exponent)
        candidate_primal = published_primal + epsilon * direction_primal.T
        pairings = candidate_primal * oriented_wall_matrix.T
        candidate_dual = candidate_primal * ambient_gram
        square = sp.Rational((candidate_dual * dual_gram * candidate_dual.T)[0])
        if square > 0 and all(value > 0 for value in pairings):
            return candidate_dual, {
                "method": "published_Weyl_plus_exact_dyadic_zero_wall_perturbation",
                "original_zero_wall_positions": [index + 1 for index in zero_positions],
                "epsilon": _rational_text(epsilon),
                "perturbation_primal_coordinates": [
                    _rational_text(sp.Rational(value)) for value in direction_primal
                ],
            }
    raise ArithmeticError("the exact zero-wall perturbation did not enter the chamber")


def _reflection_integrality(
    primitive_normal: sp.Matrix, fixed_gram: sp.Matrix
) -> tuple[bool, sp.Matrix]:
    square = sp.Rational((primitive_normal.T * fixed_gram * primitive_normal)[0])
    if square >= 0:
        raise ArithmeticError("reflection requested for a nonspacelike mirror")
    reflection = (
        sp.eye(fixed_gram.rows)
        - sp.Rational(2, 1)
        * primitive_normal
        * primitive_normal.T
        * fixed_gram
        / square
    )
    if not (
        reflection * reflection == sp.eye(fixed_gram.rows)
        and reflection.T * fixed_gram * reflection == fixed_gram
    ):
        raise ArithmeticError("a restricted mirror reflection is not Lorentzian")
    return all(value.is_Integer for value in reflection), reflection


def _simplicial_exact_facet_witness(
    target: int,
    fixed_gram: sp.Matrix,
    pairing_matrix: sp.Matrix,
) -> tuple[sp.Matrix, dict[str, Any]] | None:
    """Construct a strict witness from exact slack coordinates when square.

    If the folded cone is simplicial, ``A x`` gives its complete slack
    coordinates.  A timelike coordinate ray on the target face can be
    perturbed by an exact dyadic amount toward every other coordinate ray.
    This avoids treating a numerically ill-conditioned phase-I solve as a
    geometric non-face.
    """

    inequality = pairing_matrix.T
    if inequality.rows != inequality.cols or inequality.det() == 0:
        return None
    inverse = inequality.inv()
    other = tuple(index for index in range(inequality.rows) if index != target)
    for seed in other:
        slack_seed = sp.zeros(inequality.rows, 1)
        slack_seed[seed] = 1
        seed_vector = inverse * slack_seed
        if (seed_vector.T * fixed_gram * seed_vector)[0] <= 0:
            continue
        for exponent in range(1, 81):
            epsilon = sp.Rational(1, 2**exponent)
            slacks = sp.zeros(inequality.rows, 1)
            for index in other:
                slacks[index] = 1 if index == seed else epsilon
            witness = inverse * slacks
            square = sp.Rational((witness.T * fixed_gram * witness)[0])
            if square <= 0:
                continue
            if inequality * witness != slacks:
                raise ArithmeticError("the simplicial exact witness failed replay")
            return witness, {
                "method": "exact_simplicial_slack_ray_plus_dyadic_interior_perturbation",
                "timelike_seed_slack_orbit_index": seed + 1,
                "epsilon": _rational_text(epsilon),
                "witness_square": _rational_text(square),
                "minimum_strict_pairing": _rational_text(epsilon),
            }
    return None


def _exact_face_ray_facet_witness(
    target: int,
    fixed_gram: sp.Matrix,
    pairing_matrix: sp.Matrix,
) -> tuple[sp.Matrix, dict[str, Any]] | None:
    """Enumerate the target face rays exactly and perturb a timelike one.

    An extreme ray in rank ``r`` contains ``r-1`` independent active
    inequalities.  Since the target inequality is fixed active, enumerating
    all ``r-2``-subsets of the other rows is exhaustive.  This fallback is
    reserved for a zero-bad preliminary angle gate, so hard negative cases do
    not pay its combinatorial cost.
    """

    inequality = pairing_matrix.T
    row_count, rank = inequality.shape
    if row_count < rank or rank < 2:
        return None
    other = tuple(index for index in range(row_count) if index != target)
    rays: dict[tuple[int, ...], sp.Matrix] = {}
    for subset in itertools.combinations(other, rank - 2):
        active = (target,) + subset
        system = inequality[list(active), :]
        if system.rank() != rank - 1:
            continue
        kernel = system.nullspace()
        if len(kernel) != 1:
            continue
        vector = kernel[0]
        slacks = inequality * vector
        if all(value >= 0 for value in slacks):
            pass
        elif all(value <= 0 for value in slacks):
            vector = -vector
            slacks = -slacks
        else:
            continue
        denominator = int(sp.ilcm(*(sp.Rational(value).q for value in vector)))
        integers = [int(sp.Rational(value) * denominator) for value in vector]
        divisor = math.gcd(*(abs(value) for value in integers))
        primitive = tuple(value // divisor for value in integers)
        primitive_vector = sp.Matrix(primitive)
        if any(value < 0 for value in inequality * primitive_vector):
            raise ArithmeticError("a primitive exact face ray changed cone orientation")
        rays[primitive] = primitive_vector
    if not rays:
        return None

    interior = sum(rays.values(), sp.zeros(rank, 1))
    interior_slacks = inequality * interior
    if not (
        interior_slacks[target] == 0
        and all(
            value > 0
            for index, value in enumerate(interior_slacks)
            if index != target
        )
    ):
        return None
    interior_square = sp.Rational((interior.T * fixed_gram * interior)[0])
    if interior_square > 0:
        return interior, {
            "method": "exhaustive_exact_face_ray_sum",
            "enumerated_face_extreme_ray_count": len(rays),
            "witness_square": _rational_text(interior_square),
            "minimum_strict_pairing": _rational_text(
                min(
                    value
                    for index, value in enumerate(interior_slacks)
                    if index != target
                )
            ),
        }

    seeds = list(rays.values())
    seed: sp.Matrix | None = next(
        (
            vector
            for vector in seeds
            if (vector.T * fixed_gram * vector)[0] > 0
        ),
        None,
    )
    if seed is None:
        seed = next(
            (
                left + right
                for position, left in enumerate(seeds)
                for right in seeds[position + 1 :]
                if ((left + right).T * fixed_gram * (left + right))[0] > 0
            ),
            None,
        )
    if seed is None:
        return None
    seed_square = sp.Rational((seed.T * fixed_gram * seed)[0])
    for exponent in range(1, 121):
        epsilon = sp.Rational(1, 2**exponent)
        witness = seed + epsilon * interior
        square = sp.Rational((witness.T * fixed_gram * witness)[0])
        slacks = inequality * witness
        if not (
            square > 0
            and slacks[target] == 0
            and all(
                value > 0
                for index, value in enumerate(slacks)
                if index != target
            )
        ):
            continue
        return witness, {
            "method": "exhaustive_exact_face_rays_plus_dyadic_interior_perturbation",
            "enumerated_face_extreme_ray_count": len(rays),
            "timelike_seed_square": _rational_text(seed_square),
            "epsilon": _rational_text(epsilon),
            "witness_square": _rational_text(square),
            "minimum_strict_pairing": _rational_text(
                min(value for index, value in enumerate(slacks) if index != target)
            ),
        }
    return None


def screen_case(
    atlas: ShimadaDisc11Atlas,
    inventory_hash: str,
    chamber_record: dict[str, Any],
    class_record: dict[str, Any],
) -> dict[str, Any]:
    chamber_index = int(chamber_record["chamber_index"])
    class_index = int(class_record["class_index"])
    chamber = atlas.chamber(chamber_index)
    if not (
        chamber_record["class_source"]
        == "fresh_full_weighted_wall_group_enumeration"
        and chamber_record["groups_equal"] is False
        and chamber.wall_count == chamber_record["wall_count"]
        and chamber.order_aut == chamber_record["archived_stabilizer_order"]
        and class_record["order"] > 1
        and class_record["fixed_rank"] >= MIN_FIXED_RANK
    ):
        raise ArithmeticError("the selected fresh weighted case changed")

    ambient_gram = sp.Matrix(atlas.gram)
    dual_gram = ambient_gram.inv()
    raw_walls_tuple = tuple(wall for orbit in chamber.wall_orbits for wall in orbit)
    raw_wall_matrix = sp.Matrix(raw_walls_tuple)
    if len(set(raw_walls_tuple)) != len(raw_walls_tuple):
        raise ArithmeticError("the published chamber has duplicate walls")
    published_ample_tuple = _row_product(chamber.weyl, atlas.projection_to_s_dual)
    published_ample_dual = sp.Matrix(1, AMBIENT_RANK, published_ample_tuple)
    raw_ample_pairings = tuple(
        sp.Rational(value)
        for value in published_ample_dual * dual_gram * raw_wall_matrix.T
    )
    # Shimada's wallorbit rows are mirror vectors, but _isom wall records are
    # not uniformly oriented.  Orient every half-space toward the published
    # Weyl interior point before folding or applying the Coxeter angle gate.
    if chamber.order_aut != 1 or any(len(orbit) != 1 for orbit in chamber.wall_orbits):
        raise ArithmeticError("a fresh strict-extension chamber is not archived-trivial")
    orientation_methods: list[str] = []
    orientation_signs_list: list[int] = []
    for orbit_index, (value, wall) in enumerate(
        zip(raw_ample_pairings, raw_walls_tuple, strict=True)
    ):
        if value > 0:
            sign, method = 1, "positive_published_Weyl_pairing"
        elif value < 0:
            sign, method = -1, "negative_published_Weyl_pairing_reversed"
        else:
            sign, method = _zero_weyl_wall_orientation(
                atlas, chamber_index, orbit_index, wall, dual_gram
            )
        orientation_signs_list.append(sign)
        orientation_methods.append(method)
    orientation_signs = tuple(orientation_signs_list)
    walls_tuple = tuple(
        tuple(sign * value for value in wall)
        for sign, wall in zip(orientation_signs, raw_walls_tuple, strict=True)
    )
    walls = tuple(sp.Matrix(1, AMBIENT_RANK, wall) for wall in walls_tuple)
    wall_matrix = sp.Matrix(walls_tuple)
    source_interior_dual, source_interior_record = _strict_interior_from_oriented_walls(
        published_ample_dual,
        ambient_gram,
        wall_matrix,
        tuple(index for index, value in enumerate(raw_ample_pairings) if value == 0),
    )

    action_lookup = _weighted_actions(chamber_record)
    subgroup_permutations: Subgroup = frozenset(
        tuple(int(value) for value in permutation)
        for permutation in class_record["elements"]
    )
    if not (
        len(subgroup_permutations) == class_record["order"]
        and subgroup_permutations.issubset(action_lookup)
    ):
        raise ArithmeticError("the fresh weighted subgroup is incomplete")
    dual_actions = tuple(action_lookup[element] for element in subgroup_permutations)
    primal_actions = tuple(
        ambient_gram * action * dual_gram for action in dual_actions
    )
    action_integrality = []
    for permutation, action, primal in zip(
        subgroup_permutations, dual_actions, primal_actions, strict=True
    ):
        expected_raw_walls = raw_wall_matrix[
            [permutation[index] for index in range(len(walls))], :
        ]
        expected_walls = wall_matrix[
            [permutation[index] for index in range(len(walls))], :
        ]
        if not (
            raw_wall_matrix * action == expected_raw_walls
            and all(
                orientation_signs[index] == orientation_signs[permutation[index]]
                for index in range(len(walls))
            )
            and wall_matrix * action == expected_walls
            and action * dual_gram * action.T == dual_gram
            and primal * ambient_gram * primal.T == ambient_gram
        ):
            raise ArithmeticError("a rational weighted-wall action failed exact replay")
        action_integrality.append(all(value.is_Integer for value in primal))

    fixed = _saturated_rational_fixed_divisor_lattice(
        primal_actions, ambient_gram, int(class_record["fixed_rank"])
    )
    fixed_basis: sp.Matrix = fixed["basis"]
    fixed_gram: sp.Matrix = fixed["gram"]
    fixed_rank = fixed_basis.rows

    fixed_ample_dual = sum(
        (source_interior_dual * action for action in dual_actions),
        sp.zeros(1, AMBIENT_RANK),
    ) / len(dual_actions)
    if any(fixed_ample_dual * action != fixed_ample_dual for action in dual_actions):
        raise ArithmeticError("the averaged Weyl interior vector is not fixed")
    ample = _fixed_coordinates_of_ample(
        fixed_ample_dual, ambient_gram, fixed_basis
    )
    if (ample.T * fixed_gram * ample)[0] <= 0:
        raise ArithmeticError("the averaged Weyl vector is not timelike")

    wall_orbits = _orbits(subgroup_permutations, len(walls))
    orbit_normals = tuple(
        sum((walls[index] for index in orbit), sp.zeros(1, AMBIENT_RANK))
        for orbit in wall_orbits
    )
    if any(
        normal * action != normal
        for normal in orbit_normals
        for action in dual_actions
    ):
        raise ArithmeticError("a folded wall normal is not fixed")
    pairing_matrix = sp.Matrix.hstack(
        *(fixed_basis * normal.T for normal in orbit_normals)
    )
    if pairing_matrix.rank() != fixed_rank:
        raise ArithmeticError("the restricted walls do not span the fixed divisor dual")
    inverse_fixed_gram = fixed_gram.inv()
    orbit_gram = pairing_matrix.T * inverse_fixed_gram * pairing_matrix
    ambient_orbit_gram = sp.Matrix.vstack(*orbit_normals) * dual_gram * sp.Matrix.vstack(
        *orbit_normals
    ).T
    if not (
        orbit_gram == ambient_orbit_gram
        and exact_inertia(orbit_gram)
        == (1, fixed_rank - 1, len(wall_orbits) - fixed_rank)
    ):
        raise ArithmeticError("the fixed mirror Gram is not the complete Lorentz span")
    coordinates = pairing_matrix.T
    ample_pairings = tuple(sp.Rational(value) for value in ample.T * pairing_matrix)
    if any(value <= 0 for value in ample_pairings):
        raise ArithmeticError("the averaged Weyl vector is not strictly in the folded cone")

    conic_facets_list: list[int] = []
    conic_redundant_list: list[int] = []
    for target in range(len(wall_orbits)):
        dependence = _exact_nonnegative_dependence(coordinates, target)
        (conic_redundant_list if dependence is not None else conic_facets_list).append(
            target
        )
    conic_facets = tuple(conic_facets_list)
    conic_redundant = tuple(conic_redundant_list)
    conic_redundancy_certificates = []
    for target in conic_redundant:
        dependence = _dependence_on_retained_facets(
            coordinates, target, conic_facets
        )
        conic_redundancy_certificates.append(
            {
                "orbit_index": target + 1,
                "terms": [
                    {
                        "conic_facet_orbit_index": support + 1,
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
    for target in conic_facets:
        norm = sp.Rational(orbit_gram[target, target])
        if norm >= 0:
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
            continue
        witness, record = result
        square, pairings = _audit_witness(
            witness, target, fixed_gram, pairing_matrix
        )
        record |= {
            "orbit_index": target + 1,
            "fixed_divisor_lattice_coordinates": [
                _rational_text(sp.Rational(value)) for value in witness
            ],
            "witness_square": _rational_text(square),
            "minimum_strict_pairing": _rational_text(
                min(value for index, value in enumerate(pairings) if index != target)
            ),
        }
        active_facets.append(target)
        witnesses.append(record)

    preliminary_facet_gram = orbit_gram.extract(active_facets, active_facets)
    preliminary_gate = _pair_gate(preliminary_facet_gram)
    if unresolved and preliminary_gate["non_coxeter_pair_count"] == 0:
        still_unresolved = []
        for target in unresolved:
            result = _exact_face_ray_facet_witness(
                target, fixed_gram, pairing_matrix
            )
            if result is None:
                still_unresolved.append(target)
                continue
            witness, record = result
            square, pairings = _audit_witness(
                witness, target, fixed_gram, pairing_matrix
            )
            record |= {
                "orbit_index": target + 1,
                "fixed_divisor_lattice_coordinates": [
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
            active_facets.append(target)
            witnesses.append(record)
        unresolved = still_unresolved

    active_facets.sort()
    witnesses.sort(key=lambda record: int(record["orbit_index"]))
    facet_gram = orbit_gram.extract(active_facets, active_facets)
    pair_gate = _pair_gate(facet_gram)
    if pair_gate["non_coxeter_pair_count"]:
        decision = "not_coxeter"
    elif unresolved:
        decision = "unresolved_facet_causality"
    else:
        decision = "coxeter_candidate_requires_volume_audit"

    orbit_records = []
    reflection_integrality: dict[int, bool] = {}
    for orbit_index, (orbit, normal) in enumerate(
        zip(wall_orbits, orbit_normals, strict=True)
    ):
        fixed_normal = inverse_fixed_gram * pairing_matrix[:, orbit_index]
        primitive, scale = _primitive(fixed_normal)
        norm = sp.Rational(orbit_gram[orbit_index, orbit_index])
        reflection_integral: bool | None = None
        if norm < 0:
            reflection_integral, _ = _reflection_integrality(primitive, fixed_gram)
            reflection_integrality[orbit_index] = reflection_integral
        orbit_records.append(
            {
                "orbit_index": orbit_index + 1,
                "source_wall_indices": [index + 1 for index in orbit],
                "ambient_dual_normal": [int(value) for value in normal],
                "fixed_divisor_pairing_column": [
                    _rational_text(sp.Rational(value))
                    for value in pairing_matrix[:, orbit_index]
                ],
                "rational_fixed_mirror_normal_coordinates": [
                    _rational_text(sp.Rational(value)) for value in fixed_normal
                ],
                "primitive_integral_same_mirror_direction_coordinates": [
                    int(value) for value in primitive
                ],
                "presentation_scale_from_primitive_direction": _rational_text(scale),
                "square": _rational_text(norm),
                "reflection_integral_on_saturated_fixed_divisor_lattice": reflection_integral,
                "conic_status": (
                    "facet" if orbit_index in conic_facets else "redundant"
                ),
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

    active_reflection_flags = [reflection_integrality[index] for index in active_facets]
    divisor_lattice_reflection_stable = bool(active_facets) and all(
        active_reflection_flags
    )
    nonidentity_actions = [
        (permutation, action, primal, integral)
        for permutation, action, primal, integral in zip(
            subgroup_permutations,
            dual_actions,
            primal_actions,
            action_integrality,
            strict=True,
        )
        if action != sp.eye(AMBIENT_RANK)
    ]

    return {
        "schema": SCHEMA,
        "id": _case_id(chamber_index, class_index),
        "source": {
            "primary_disc16_sha256": EXPECTED_SOURCE_SHA256,
            "inventory_path": str(INVENTORY.relative_to(ROOT)),
            "inventory_sha256": inventory_hash,
            "class_source": chamber_record["class_source"],
        },
        "scope": {
            "group": "subgroup of the full exact weighted-wall graph automorphism group",
            "strict_extension_beyond_archived_AutD": True,
            "rational_Lorentz_lifts_allowed": True,
            "surface_automorphism_interpretation_used_as_gate": False,
            "stock_bounded_to_one_complete_published_induced_chamber": True,
            "fixed_divisor_lattice_and_mirror_reflection_lattice_conflated": False,
        },
        "case": {
            "chamber_index": chamber_index,
            "class_index": class_index,
            "published_wall_count": chamber.wall_count,
            "archived_stabilizer_order": chamber.order_aut,
            "weighted_wall_group_order": chamber_record["weighted_wall_group_order"],
            "subgroup_order": class_record["order"],
            "subgroup_conjugacy_class_size": class_record["conjugacy_class_size"],
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "subgroup_wall_permutations": [
                list(permutation) for permutation in sorted(subgroup_permutations)
            ],
            "published_wall_orientation_sign_histogram": {
                str(sign): count
                for sign, count in sorted(Counter(orientation_signs).items())
            },
            "published_wall_orientation_method_histogram": dict(
                sorted(Counter(orientation_methods).items())
            ),
            "wall_halfspaces_oriented_by_published_Weyl_vector": True,
            "strict_source_interior": source_interior_record,
            "nonidentity_rational_actions": [
                {
                    "wall_permutation": list(permutation),
                    "dual_row_action": _matrix_payload(action),
                    "primal_row_action": _matrix_payload(primal),
                    "primal_action_integral": integral,
                }
                for permutation, action, primal, integral in nonidentity_actions
            ],
        },
        "saturated_fixed_divisor_lattice": {
            "definition": "NS(X) intersect Fix_Q(H)",
            "basis_primal_rows": _matrix_payload(fixed_basis),
            "gram": _matrix_payload(fixed_gram),
            "determinant": fixed["determinant"],
            "determinant_absolute": abs(fixed["determinant"]),
            "inertia": list(exact_inertia(fixed_gram)),
            "smith_diagonal_absolute": fixed["smith_diagonal_absolute"],
            "rational_constraint_clearing_denominator": fixed[
                "clearing_denominator"
            ],
            "smith_constraint_diagonal": fixed["smith_constraint_diagonal"],
            "averaged_ample_fixed_coordinates": [
                _rational_text(sp.Rational(value)) for value in ample
            ],
            "averaged_ample_square": _rational_text(
                sp.Rational((ample.T * fixed_gram * ample)[0])
            ),
            "minimum_folded_wall_pairing": _rational_text(min(ample_pairings)),
            "saturated_smith_kernel_after_denominator_clearing": True,
        },
        "restricted_mirror_geometry": {
            "wall_orbit_count": len(wall_orbits),
            "wall_orbit_size_histogram": {
                str(size): count
                for size, count in sorted(Counter(map(len, wall_orbits)).items())
            },
            "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
            "orbit_records": orbit_records,
            "retained_mirror_reflection_integrality_histogram": {
                str(flag).lower(): count
                for flag, count in sorted(Counter(active_reflection_flags).items())
            },
            "saturated_fixed_divisor_lattice_is_stable_under_all_retained_mirror_reflections": divisor_lattice_reflection_stable,
            "reflection_stable_mirror_lattice_status": (
                "saturated_fixed_divisor_lattice itself is reflection-stable"
                if divisor_lattice_reflection_stable
                else "distinct reflection-stable mirror lattice not constructed by this angle screen"
            ),
        },
        "exact_reduction": {
            "conic_facet_count": len(conic_facets),
            "conic_redundant_count": len(conic_redundant),
            "conic_facet_orbit_indices": [index + 1 for index in conic_facets],
            "conic_redundancy_certificates": conic_redundancy_certificates,
            "causally_inactive_conic_facet_count": len(causal_inactive),
            "causally_inactive_orbit_indices": [
                index + 1 for index in causal_inactive
            ],
            "unresolved_conic_facet_count": len(unresolved),
            "unresolved_orbit_indices": [index + 1 for index in unresolved],
            "hyperbolic_facet_count": len(active_facets),
            "hyperbolic_facet_orbit_indices": [index + 1 for index in active_facets],
            "strict_timelike_facet_witnesses": witnesses,
        },
        "coxeter_gate": pair_gate
        | {
            "facet_gram": _matrix_payload(facet_gram),
            "positive_wall_gram": _matrix_payload(-facet_gram),
        },
        "decision": {
            "classification": decision,
            "all_hyperbolic_facets_resolved": not unresolved,
            "all_resolved_facet_pairs_coxeter": pair_gate["all_pairs_coxeter"],
            "finite_volume_claimed": False,
        },
        "checks": {
            "fresh_weighted_group_actions_replayed_exactly": True,
            "rational_primal_fixed_space_intersected_with_NS_integrally": True,
            "saturated_fixed_divisor_lattice_certified": True,
            "all_published_walls_folded": True,
            "orbit_normals_span_complete_fixed_divisor_dual": True,
            "ambient_and_fixed_mirror_grams_agree": True,
            "all_conic_redundancies_use_final_conic_facets": True,
            "every_resolved_hyperbolic_facet_witness_replayed_exactly": True,
            "every_retained_mirror_reflection_checked_on_divisor_lattice": True,
            "every_resolved_facet_pair_checked_exactly": True,
        },
    }


def _inventory_cases(
    inventory: dict[str, Any],
) -> tuple[tuple[dict[str, Any], dict[str, Any]], ...]:
    cases = []
    for chamber in inventory["nontrivial_weighted_group_chambers"]:
        if chamber["class_source"] != "fresh_full_weighted_wall_group_enumeration":
            continue
        for class_record in chamber["classes"]:
            if class_record["order"] > 1 and class_record["fixed_rank"] >= MIN_FIXED_RANK:
                cases.append((chamber, class_record))
    return tuple(
        sorted(
            cases,
            key=lambda pair: (
                0 if pair[1]["hyperbolic_dimension"] == 13 else 1,
                pair[0]["wall_count"],
                pair[0]["chamber_index"],
                pair[1]["class_index"],
            ),
        )
    )


def _valid_case_artifact(path: Path, inventory_hash: str) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not (
        payload.get("schema") == SCHEMA
        and payload.get("source", {}).get("primary_disc16_sha256")
        == EXPECTED_SOURCE_SHA256
        and payload.get("source", {}).get("inventory_sha256") == inventory_hash
        and payload.get("source", {}).get("class_source")
        == "fresh_full_weighted_wall_group_enumeration"
    ):
        return None
    return payload


def _summary_row(payload: dict[str, Any], case_directory: Path) -> dict[str, Any]:
    case = payload["case"]
    reduction = payload["exact_reduction"]
    gate = payload["coxeter_gate"]
    mirror = payload["restricted_mirror_geometry"]
    return {
        "id": payload["id"],
        "chamber_index": case["chamber_index"],
        "class_index": case["class_index"],
        "wall_count": case["published_wall_count"],
        "subgroup_order": case["subgroup_order"],
        "fixed_rank": case["fixed_rank"],
        "hyperbolic_dimension": case["hyperbolic_dimension"],
        "wall_orbit_count": mirror["wall_orbit_count"],
        "conic_facet_count": reduction["conic_facet_count"],
        "hyperbolic_facet_count": reduction["hyperbolic_facet_count"],
        "unresolved_facet_count": reduction["unresolved_conic_facet_count"],
        "non_coxeter_pair_count": gate["non_coxeter_pair_count"],
        "non_coxeter_ratio_histogram": gate["non_coxeter_ratio_histogram"],
        "divisor_lattice_reflection_stable": mirror[
            "saturated_fixed_divisor_lattice_is_stable_under_all_retained_mirror_reflections"
        ],
        "classification": payload["decision"]["classification"],
        "case_certificate": str(
            _case_path(
                case_directory, case["chamber_index"], case["class_index"]
            ).relative_to(ROOT)
        ),
    }


def _audit_case_payload(payload: dict[str, Any]) -> None:
    """Replay one frozen no-go case using only exact archived arithmetic."""

    fixed_gram = _matrix(payload["saturated_fixed_divisor_lattice"]["gram"])
    records = payload["restricted_mirror_geometry"]["orbit_records"]
    normals = {
        int(record["orbit_index"]): sp.Matrix(
            [
                sp.Rational(value)
                for value in record["rational_fixed_mirror_normal_coordinates"]
            ]
        )
        for record in records
    }
    pairings = {
        int(record["orbit_index"]): sp.Matrix(
            [sp.Rational(value) for value in record["fixed_divisor_pairing_column"]]
        )
        for record in records
    }
    for record in records:
        index = int(record["orbit_index"])
        if not (
            fixed_gram * normals[index] == pairings[index]
            and (normals[index].T * fixed_gram * normals[index])[0]
            == sp.Rational(record["square"])
        ):
            raise ArithmeticError(f"{payload['id']} mirror {index} failed replay")

    reduction = payload["exact_reduction"]
    conic = set(map(int, reduction["conic_facet_orbit_indices"]))
    redundant = {
        int(record["orbit_index"]): record
        for record in reduction["conic_redundancy_certificates"]
    }
    if not (
        len(conic) == reduction["conic_facet_count"]
        and len(redundant) == reduction["conic_redundant_count"]
        and conic.isdisjoint(redundant)
        and conic | set(redundant) == set(normals)
    ):
        raise ArithmeticError(f"{payload['id']} conic partition failed replay")
    for target, certificate in redundant.items():
        recovered = sum(
            (
                sp.Rational(term["coefficient"])
                * normals[int(term["conic_facet_orbit_index"])]
                for term in certificate["terms"]
            ),
            sp.zeros(fixed_gram.rows, 1),
        )
        if recovered != normals[target]:
            raise ArithmeticError(f"{payload['id']} redundancy failed replay")

    active = tuple(map(int, reduction["hyperbolic_facet_orbit_indices"]))
    witnesses = {
        int(record["orbit_index"]): record
        for record in reduction["strict_timelike_facet_witnesses"]
    }
    if set(active) != set(witnesses):
        raise ArithmeticError(f"{payload['id']} witness list is incomplete")
    for target, record in witnesses.items():
        witness = sp.Matrix(
            [
                sp.Rational(value)
                for value in record["fixed_divisor_lattice_coordinates"]
            ]
        )
        values = {index: (witness.T * pairing)[0] for index, pairing in pairings.items()}
        if not (
            (witness.T * fixed_gram * witness)[0]
            == sp.Rational(record["witness_square"])
            > 0
            and values[target] == 0
            and all(value > 0 for index, value in values.items() if index != target)
        ):
            raise ArithmeticError(f"{payload['id']} facet witness failed replay")

    facet_gram = sp.Matrix(
        [
            [(normals[left].T * pairings[right])[0] for right in active]
            for left in active
        ]
    )
    gate = payload["coxeter_gate"]
    if not (
        _matrix_payload(facet_gram) == gate["facet_gram"]
        and _pair_gate(facet_gram) == {
            key: gate[key]
            for key in (
                "pair_count",
                "coxeter_edge_counts",
                "non_coxeter_pair_count",
                "non_coxeter_ratio_histogram",
                "non_coxeter_pairs",
                "all_pairs_coxeter",
            )
        }
        and gate["non_coxeter_pair_count"] >= 1
        and payload["decision"]["classification"] == "not_coxeter"
    ):
        raise ArithmeticError(f"{payload['id']} exact no-go pair gate failed replay")


def _audit_no_go_gate_payload(payload: dict[str, Any]) -> None:
    """Replay the decisive exact Gram/pair obstruction with low memory."""

    gate = payload["coxeter_gate"]
    facet_gram = _matrix(gate["facet_gram"])
    replay = _pair_gate(facet_gram)
    archived = {
        key: gate[key]
        for key in (
            "pair_count",
            "coxeter_edge_counts",
            "non_coxeter_pair_count",
            "non_coxeter_ratio_histogram",
            "non_coxeter_pairs",
            "all_pairs_coxeter",
        )
    }
    active = set(payload["exact_reduction"]["hyperbolic_facet_orbit_indices"])
    witnessed = {
        record["orbit_index"]
        for record in payload["exact_reduction"]["strict_timelike_facet_witnesses"]
    }
    if not (
        replay == archived
        and replay["non_coxeter_pair_count"] >= 1
        and payload["decision"]["classification"] == "not_coxeter"
        and active == witnessed
        and facet_gram.rows == len(active)
    ):
        raise ArithmeticError(f"{payload['id']} decisive no-go gate failed replay")


def audit_archive(
    inventory_hash: str,
    cases: tuple[tuple[dict[str, Any], dict[str, Any]], ...],
    case_directory: Path = CASE_DIRECTORY,
    manifest_path: Path = MANIFEST,
) -> dict[str, Any]:
    manifest = _build_manifest(inventory_hash, cases, case_directory)
    if not (
        manifest["summary"]["completed_case_count"] == EXPECTED_CASE_COUNT
        and manifest["summary"]["pending_case_count"] == 0
        and manifest["summary"]["classification_histogram"]
        == {"not_coxeter": EXPECTED_CASE_COUNT}
        and manifest["summary"]["coxeter_candidate_count_requiring_volume_audit"]
        == 0
        and "0" not in manifest["summary"]["bad_pair_count_histogram"]
    ):
        raise ArithmeticError("the complete fresh weighted no-go denominator changed")
    full_replay_ids = {"W0600_C001", "W0584_C001", "W1097_C001", "W2884_C001"}
    full_replay_payloads = []
    for row in manifest["cases"]:
        path = ROOT / row["case_certificate"]
        payload = _valid_case_artifact(path, inventory_hash)
        if payload is None:
            raise ArithmeticError(f"invalid case artifact {path}")
        _audit_no_go_gate_payload(payload)
        if payload["id"] in full_replay_ids:
            full_replay_payloads.append(payload)
    if {payload["id"] for payload in full_replay_payloads} != full_replay_ids:
        raise ArithmeticError("the selected full arithmetic replay cases are incomplete")
    for payload in full_replay_payloads:
        _audit_case_payload(payload)
    archived = json.loads(manifest_path.read_text(encoding="utf-8"))
    if archived != manifest:
        raise ArithmeticError("the frozen fresh weighted manifest differs from replay")
    return manifest


def _build_manifest(
    inventory_hash: str,
    cases: tuple[tuple[dict[str, Any], dict[str, Any]], ...],
    case_directory: Path,
) -> dict[str, Any]:
    rows = []
    for chamber, record in cases:
        path = _case_path(
            case_directory, chamber["chamber_index"], record["class_index"]
        )
        payload = _valid_case_artifact(path, inventory_hash)
        if payload is None:
            rows.append(
                {
                    "id": _case_id(chamber["chamber_index"], record["class_index"]),
                    "chamber_index": chamber["chamber_index"],
                    "class_index": record["class_index"],
                    "wall_count": chamber["wall_count"],
                    "fixed_rank": record["fixed_rank"],
                    "hyperbolic_dimension": record["hyperbolic_dimension"],
                    "status": "pending",
                }
            )
        else:
            rows.append(_summary_row(payload, case_directory) | {"status": "complete"})
    completed = [row for row in rows if row["status"] == "complete"]
    pending = [row for row in rows if row["status"] == "pending"]
    survivors = [
        row
        for row in completed
        if row["classification"] == "coxeter_candidate_requires_volume_audit"
    ]
    bad_ratio_histogram: Counter[str] = Counter()
    bad_pair_count_histogram: Counter[int] = Counter()
    for row in completed:
        bad_pair_count_histogram[row["non_coxeter_pair_count"]] += 1
        bad_ratio_histogram.update(row["non_coxeter_ratio_histogram"])
    return {
        "schema": MANIFEST_SCHEMA,
        "source": {
            "primary_disc16_sha256": EXPECTED_SOURCE_SHA256,
            "inventory_path": str(INVENTORY.relative_to(ROOT)),
            "inventory_sha256": inventory_hash,
            "class_source": "fresh_full_weighted_wall_group_enumeration",
        },
        "scope": {
            "denominator": "575 proper above-H9 involution classes from strict extensions of archived Aut(D) to full exact weighted-wall graph groups",
            "dimension_histogram": {"13": 1, "18": 574},
            "surface_automorphism_interpretation_used_as_gate": False,
            "case_order": "unique H13 first, then increasing wall count/chamber/class for H18",
            "resumable_per_case_artifacts": True,
            "finite_volume_gate_performed": False,
            "saturated_fixed_divisor_and_reflection_mirror_lattices_distinguished": True,
        },
        "summary": {
            "inventory_case_count": len(cases),
            "completed_case_count": len(completed),
            "pending_case_count": len(pending),
            "classification_histogram": dict(
                sorted(Counter(row["classification"] for row in completed).items())
            ),
            "bad_pair_count_histogram": {
                str(count): value
                for count, value in sorted(bad_pair_count_histogram.items())
            },
            "bad_normalized_ratio_histogram": dict(sorted(bad_ratio_histogram.items())),
            "coxeter_candidate_count_requiring_volume_audit": len(survivors),
        },
        "survivors_requiring_volume_audit": survivors,
        "cases": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--inventory", type=Path, default=INVENTORY)
    parser.add_argument("--case-directory", type=Path, default=CASE_DIRECTORY)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--max-cases", type=int)
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        help="restrict to CHAMBER:CLASS; may be repeated",
    )
    parser.add_argument("--recompute", action="store_true")
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument("--audit-archive", action="store_true")
    arguments = parser.parse_args()

    if _sha256(arguments.source) != EXPECTED_SOURCE_SHA256:
        raise SystemExit("the primary disc16 source hash changed")
    inventory_hash = _sha256(arguments.inventory)
    inventory = json.loads(arguments.inventory.read_text(encoding="utf-8"))
    if not (
        inventory["source"]["sha256"] == EXPECTED_SOURCE_SHA256
        and inventory["summary"]["chamber_count"] == 4539
        and inventory["summary"]["weighted_group_strict_extension_count"] == 575
    ):
        raise SystemExit("the frozen weighted-wall inventory changed")
    all_cases = _inventory_cases(inventory)
    if not (
        len(all_cases) == EXPECTED_CASE_COUNT
        and Counter(record["hyperbolic_dimension"] for _, record in all_cases)
        == {13: 1, 18: 574}
    ):
        raise SystemExit(f"expected the 1+574 fresh high-rank queue, got {len(all_cases)}")
    requested = {tuple(map(int, text.split(":"))) for text in arguments.case}
    if any(len(pair) != 2 for pair in requested):
        raise SystemExit("--case must have CHAMBER:CLASS form")

    if arguments.audit_archive:
        manifest = audit_archive(
            inventory_hash, all_cases, arguments.case_directory, arguments.manifest
        )
        print("fresh weighted high-section archive replay passed")
        print(json.dumps(manifest["summary"], indent=2, sort_keys=True))
        return

    arguments.case_directory.mkdir(parents=True, exist_ok=True)
    atlas = ShimadaDisc11Atlas(arguments.source)
    processed = 0
    for chamber, record in all_cases:
        chamber_index = int(chamber["chamber_index"])
        class_index = int(record["class_index"])
        if requested and (chamber_index, class_index) not in requested:
            continue
        path = _case_path(arguments.case_directory, chamber_index, class_index)
        if not arguments.recompute and _valid_case_artifact(path, inventory_hash) is not None:
            continue
        if arguments.max_cases is not None and processed >= arguments.max_cases:
            break
        try:
            payload = screen_case(atlas, inventory_hash, chamber, record)
        except Exception as error:
            print(
                json.dumps(
                    {
                        "id": _case_id(chamber_index, class_index),
                        "status": "ERROR_UNRESOLVED",
                        "error_type": type(error).__name__,
                        "error": str(error),
                    },
                    sort_keys=True,
                ),
                flush=True,
            )
            if arguments.fail_fast:
                raise
            processed += 1
            continue
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        row = _summary_row(payload, arguments.case_directory)
        print(json.dumps(row, sort_keys=True), flush=True)
        if row["classification"] == "coxeter_candidate_requires_volume_audit":
            print("SURVIVOR " + json.dumps(row, sort_keys=True), flush=True)
        processed += 1
        manifest = _build_manifest(inventory_hash, all_cases, arguments.case_directory)
        arguments.manifest.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

    manifest = _build_manifest(inventory_hash, all_cases, arguments.case_directory)
    arguments.manifest.parent.mkdir(parents=True, exist_ok=True)
    arguments.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
