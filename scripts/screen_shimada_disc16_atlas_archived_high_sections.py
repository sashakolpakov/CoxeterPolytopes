"""Resumable exact screen of the archived high-rank disc16 fixed sections.

The frozen inventory contains 107 proper subgroup conjugacy classes with
fixed rank at least 11 among Shimada's 4,539 published ``disc16`` chambers.
Three are the already closed D5 classes.  This script screens the remaining
104 classes, smallest wall stocks first.

Scope is deliberately precise: these are subgroups of Shimada's archived
chamber stabilizers ``Aut(D)``.  Equality with the full abstract exact
weighted-wall automorphism group is a separate atlas-wide audit and is not
claimed here.

For each class the verifier reconstructs the integral subgroup matrices from
the inventory's faithful regular action, computes a saturated Smith-kernel
basis of the fixed lattice, folds every published wall orbit, performs exact
conic reduction, proves every retained hyperbolic facet with an exact strict
timelike witness, and checks all normalized products exactly.  SciPy is used
only to discover supports or rational witness neighbourhoods; every archived
certificate is replayable over Q.

Each case is frozen separately so an interrupted full run resumes without
losing completed work.  A failed numerical discovery is recorded as
``unresolved`` and is never converted into a negative conclusion.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.optimize import linprog, minimize
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form

from scripts.certify_shimada_disc16_d5_full_subgroup_fishing import (
    OUTPUT as D5_CENSUS,
    _dependence_on_retained_facets,
    _pair_gate,
)
from scripts.certify_shimada_disc16_d5_h14_fixed_lattice import (
    _smith_decomposition,
)
from scripts.inventory_shimada_disc16_atlas_fixed_sections import (
    EXPECTED_SOURCE_SHA256,
    OUTPUT as INVENTORY,
    SOURCE,
    _regular_action,
)
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _derive_wall_permutations,
    _exact_nonnegative_dependence,
    _matrix_payload,
    _orbits,
    _rational_text,
)
from scripts.shimada_disc11_atlas import ShimadaDisc11Atlas, _row_product
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CASE_DIRECTORY = (
    ROOT / "data" / "enumerations" / "shimada_disc16_atlas_archived_high_sections"
)
MANIFEST = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_atlas_archived_high_section_census.json"
)
SCHEMA = "coxeter.shimada_disc16_atlas_archived_high_section_case.v1"
MANIFEST_SCHEMA = "coxeter.shimada_disc16_atlas_archived_high_section_census.v1"
AMBIENT_RANK = 20
MIN_FIXED_RANK = 11


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _case_id(chamber_index: int, class_index: int) -> str:
    return f"D{chamber_index:04d}_C{class_index:03d}"


def _case_path(directory: Path, chamber_index: int, class_index: int) -> Path:
    return directory / f"{_case_id(chamber_index, class_index)}.json"


def _matrix(payload: list[list[str]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _primitive(vector: sp.Matrix) -> tuple[sp.Matrix, sp.Rational]:
    denominators = [int(sp.Rational(value).q) for value in vector]
    denominator = int(sp.ilcm(*denominators))
    integers = [int(sp.Rational(value) * denominator) for value in vector]
    common = math.gcd(*(abs(value) for value in integers))
    if common == 0:
        raise ArithmeticError("zero fixed normal")
    primitive = [value // common for value in integers]
    first = next(value for value in primitive if value)
    if first < 0:
        primitive = [-value for value in primitive]
    primitive_matrix = sp.Matrix(primitive)
    scale = next(
        sp.Rational(vector[index], primitive_matrix[index])
        for index in range(vector.rows)
        if primitive_matrix[index]
    )
    if vector != scale * primitive_matrix:
        raise ArithmeticError("primitive normalization failed")
    return primitive_matrix, scale


def _saturated_fixed_lattice(
    matrices: tuple[sp.Matrix, ...], ambient_gram: sp.Matrix, expected_rank: int
) -> dict[str, Any]:
    constraints = sp.Matrix.vstack(
        *(matrix.T - sp.eye(AMBIENT_RANK) for matrix in matrices)
    )
    diagonal, _, right = _smith_decomposition(constraints)
    zero_columns = tuple(
        column
        for column in range(AMBIENT_RANK)
        if all(diagonal[row, column] == 0 for row in range(diagonal.rows))
    )
    if len(zero_columns) != expected_rank:
        raise ArithmeticError("the Smith fixed rank disagrees with the inventory")
    smith_basis = right[:, list(zero_columns)].T
    basis = smith_basis.lll(delta=sp.Rational(3, 4))
    if not (
        basis.rank() == expected_rank
        and basis * constraints.T == sp.zeros(expected_rank, constraints.rows)
        and all(basis * matrix == basis for matrix in matrices)
    ):
        raise ArithmeticError("the saturated fixed basis failed replay")
    _, pivot_columns = smith_basis.rref()
    change = (
        basis[:, list(pivot_columns)]
        * smith_basis[:, list(pivot_columns)].inv()
    )
    if not (change * smith_basis == basis and abs(int(change.det())) == 1):
        raise ArithmeticError("LLL changed the saturated fixed lattice")
    fixed_gram = basis * ambient_gram * basis.T
    smith = smith_normal_form(fixed_gram, domain=ZZ)
    smith_diagonal = [
        abs(int(smith[index, index])) for index in range(expected_rank)
    ]
    if exact_inertia(fixed_gram) != (1, expected_rank - 1, 0):
        raise ArithmeticError("the saturated fixed lattice is not Lorentzian")
    return {
        "constraint_matrix": constraints,
        "smith_constraint_diagonal": [
            int(diagonal[index, index])
            for index in range(min(diagonal.rows, diagonal.cols))
        ],
        "basis": basis,
        "gram": fixed_gram,
        "determinant": int(fixed_gram.det()),
        "smith_diagonal_absolute": smith_diagonal,
    }


def _fixed_coordinates_of_ample(
    ample_dual: sp.Matrix,
    ambient_gram: sp.Matrix,
    fixed_basis: sp.Matrix,
) -> sp.Matrix:
    # A primal row x has dual-basis coordinates xG.  Thus the ample class in
    # primal coordinates is ample_dual G^{-1}.
    primal = ample_dual * ambient_gram.inv()
    solution = sp.linsolve((fixed_basis.T, primal.T))
    if solution is sp.EmptySet:
        raise ArithmeticError("the published Weyl vector is not fixed")
    coordinates = sp.Matrix(tuple(next(iter(solution))))
    if any(not value.is_Rational for value in coordinates):
        raise ArithmeticError("the ample fixed coordinates are not rational")
    if coordinates.T * fixed_basis != primal:
        raise ArithmeticError("the ample fixed coordinates failed replay")
    return coordinates


def _audit_witness(
    coordinates: sp.Matrix,
    target: int,
    gram: sp.Matrix,
    pairing_matrix: sp.Matrix,
) -> tuple[sp.Rational, tuple[sp.Rational, ...]]:
    square = sp.Rational((coordinates.T * gram * coordinates)[0])
    pairings = tuple(
        sp.Rational(value) for value in list(coordinates.T * pairing_matrix)
    )
    if not (
        square > 0
        and pairings[target] == 0
        and all(value > 0 for index, value in enumerate(pairings) if index != target)
    ):
        raise ArithmeticError("a proposed strict timelike witness failed exact replay")
    return square, pairings


def _projection_or_tangent_witness(
    target: int,
    gram: sp.Matrix,
    pairing_matrix: sp.Matrix,
    orbit_gram: sp.Matrix,
    ample: sp.Matrix,
) -> tuple[sp.Matrix, dict[str, Any]] | None:
    norm = sp.Rational(orbit_gram[target, target])
    if norm >= 0:
        return None
    normal = gram.inv() * pairing_matrix[:, target]
    pairing = sp.Rational((ample.T * pairing_matrix[:, target])[0])
    witness = ample - pairing / norm * normal
    square = sp.Rational((witness.T * gram * witness)[0])
    pairings = tuple(
        sp.Rational(value) for value in list(witness.T * pairing_matrix)
    )
    if (
        square > 0
        and pairings[target] == 0
        and all(value > 0 for index, value in enumerate(pairings) if index != target)
    ):
        return witness, {
            "method": "orthogonal_projection_of_published_weyl_vector",
            "witness_square": _rational_text(square),
            "minimum_strict_pairing": _rational_text(
                min(value for index, value in enumerate(pairings) if index != target)
            ),
        }
    negative = tuple(
        index
        for index, value in enumerate(pairings)
        if index != target and value < 0
    )
    zero = tuple(
        index
        for index, value in enumerate(pairings)
        if index != target and value == 0
    )
    if square <= 0 or pairings[target] != 0 or negative or not zero:
        return None
    selected = (target,) + zero
    system = pairing_matrix[:, list(selected)].T
    values = sp.Matrix([0] + [1] * len(zero))
    solution = sp.linsolve((system, values))
    if solution is sp.EmptySet:
        return None
    raw = tuple(next(iter(solution)))
    parameters = sorted(set().union(*(value.free_symbols for value in raw)), key=str)
    direction = sp.Matrix(
        [sp.Rational(value.subs({parameter: 0 for parameter in parameters})) for value in raw]
    )
    for exponent in range(1, 81):
        epsilon = sp.Rational(1, 2**exponent)
        candidate = witness + epsilon * direction
        try:
            candidate_square, candidate_pairings = _audit_witness(
                candidate, target, gram, pairing_matrix
            )
        except ArithmeticError:
            continue
        return candidate, {
            "method": "orthogonal_projection_plus_exact_tangent_perturbation",
            "initial_extra_zero_orbit_indices": [index + 1 for index in zero],
            "epsilon": _rational_text(epsilon),
            "witness_square": _rational_text(candidate_square),
            "minimum_strict_pairing": _rational_text(
                min(
                    value
                    for index, value in enumerate(candidate_pairings)
                    if index != target
                )
            ),
        }
    return None


def _exactify_equalities_near(
    approximate: np.ndarray,
    equality_matrix: sp.Matrix,
    equality_values: sp.Matrix,
    denominators: Iterable[int] = (100, 1_000, 10_000, 1_000_000, 1_000_000_000),
) -> Iterable[sp.Matrix]:
    pivots = tuple(equality_matrix.rref()[1])
    if len(pivots) != equality_matrix.rows:
        raise ArithmeticError("the witness normalization equalities are dependent")
    free = tuple(index for index in range(equality_matrix.cols) if index not in pivots)
    pivot_matrix = equality_matrix[:, list(pivots)]
    for denominator in denominators:
        values: list[sp.Rational | None] = [None] * equality_matrix.cols
        for index in free:
            values[index] = sp.Rational(
                Fraction(float(approximate[index])).limit_denominator(denominator)
            )
        residual = equality_values - equality_matrix[:, list(free)] * sp.Matrix(
            [values[index] for index in free]
        )
        pivot_values = pivot_matrix.inv() * residual
        for index, value in zip(pivots, pivot_values, strict=True):
            values[index] = sp.Rational(value)
        yield sp.Matrix(values)


def _numerically_discovered_exact_witness(
    target: int,
    gram: sp.Matrix,
    pairing_matrix: sp.Matrix,
    ample: sp.Matrix,
) -> tuple[sp.Matrix, dict[str, Any]] | None:
    rank, orbit_count = pairing_matrix.shape
    gram_float = np.array(gram.tolist(), dtype=float)
    pairing_float = np.array(pairing_matrix.tolist(), dtype=float)
    normalization = gram * ample
    normalization_float = np.array(normalization, dtype=float).reshape(-1)
    target_float = pairing_float[:, target]
    others = tuple(index for index in range(orbit_count) if index != target)

    # Phase I: maximize a common strict slack on the target face.
    inequality = np.zeros((len(others), rank + 1))
    inequality[:, :rank] = -pairing_float[:, list(others)].T
    inequality[:, rank] = 1.0
    equalities = np.zeros((2, rank + 1))
    equalities[0, :rank] = target_float
    equalities[1, :rank] = normalization_float
    objective = np.zeros(rank + 1)
    objective[rank] = -1.0
    phase_one = linprog(
        objective,
        A_ub=inequality,
        b_ub=np.zeros(len(others)),
        A_eq=equalities,
        b_eq=np.array([0.0, 1.0]),
        bounds=[(None, None)] * (rank + 1),
        method="highs",
    )
    if not phase_one.success or phase_one.x[rank] <= 1e-10:
        return None
    start = phase_one.x[:rank]
    strict_margin = phase_one.x[rank] / 2.0
    exact_equalities = sp.Matrix.vstack(
        pairing_matrix[:, target].T, normalization.T
    )
    equality_values = sp.Matrix([0, 1])

    # The phase-I point itself is sometimes already timelike.
    for candidate in _exactify_equalities_near(
        start, exact_equalities, equality_values
    ):
        try:
            square, pairings = _audit_witness(
                candidate, target, gram, pairing_matrix
            )
        except ArithmeticError:
            continue
        return candidate, {
            "method": "exact_rationalization_of_strict_phase_I_point",
            "phase_I_numerical_margin": repr(float(phase_one.x[rank])),
            "witness_square": _rational_text(square),
            "minimum_strict_pairing": _rational_text(
                min(value for index, value in enumerate(pairings) if index != target)
            ),
        }

    # On <x,ample>=1 the Lorentz square is strictly concave.  Therefore this
    # is a convex minimization of -q on the linearly constrained target face.
    other_matrix = pairing_float[:, list(others)].T
    constraints = [
        {
            "type": "eq",
            "fun": lambda vector: float(target_float @ vector),
            "jac": lambda vector: target_float,
        },
        {
            "type": "eq",
            "fun": lambda vector: float(normalization_float @ vector - 1.0),
            "jac": lambda vector: normalization_float,
        },
        {
            "type": "ineq",
            "fun": lambda vector: other_matrix @ vector - strict_margin,
            "jac": lambda vector: other_matrix,
        },
    ]
    optimization = minimize(
        lambda vector: float(-(vector @ gram_float @ vector)),
        start,
        jac=lambda vector: -2.0 * gram_float @ vector,
        constraints=constraints,
        method="SLSQP",
        options={"ftol": 1e-12, "maxiter": 2_000},
    )
    if not optimization.success:
        return None
    for candidate in _exactify_equalities_near(
        optimization.x, exact_equalities, equality_values
    ):
        try:
            square, pairings = _audit_witness(
                candidate, target, gram, pairing_matrix
            )
        except ArithmeticError:
            continue
        return candidate, {
            "method": "convex_Lorentz_QP_discovery_then_exact_rational_replay",
            "phase_I_numerical_margin": repr(float(phase_one.x[rank])),
            "QP_numerical_square": repr(
                float(optimization.x @ gram_float @ optimization.x)
            ),
            "witness_square": _rational_text(square),
            "minimum_strict_pairing": _rational_text(
                min(value for index, value in enumerate(pairings) if index != target)
            ),
        }
    return None


def _strict_timelike_witness(
    target: int,
    gram: sp.Matrix,
    pairing_matrix: sp.Matrix,
    orbit_gram: sp.Matrix,
    ample: sp.Matrix,
) -> tuple[sp.Matrix, dict[str, Any]] | None:
    simple = _projection_or_tangent_witness(
        target, gram, pairing_matrix, orbit_gram, ample
    )
    if simple is not None:
        return simple
    return _numerically_discovered_exact_witness(
        target, gram, pairing_matrix, ample
    )


def _reconstruct_inventory_subgroup(
    atlas: ShimadaDisc11Atlas,
    chamber_index: int,
    class_record: dict[str, Any],
) -> tuple[
    tuple[tuple[tuple[int, ...], ...], ...],
    tuple[Permutation, ...],
    dict[Permutation, tuple[tuple[int, ...], ...]],
    tuple[tuple[tuple[int, ...], ...], ...],
]:
    raw_group = atlas.stabilizer_group(chamber_index)
    regular_group, lookup = _regular_action(raw_group)
    requested = tuple(
        tuple(int(value) for value in permutation)
        for permutation in class_record["regular_action_elements"]
    )
    if any(permutation not in lookup for permutation in requested):
        raise ArithmeticError("an inventory regular-action element is absent")
    subgroup_matrices = tuple(lookup[permutation] for permutation in requested)
    if len(subgroup_matrices) != class_record["order"]:
        raise ArithmeticError("the inventory subgroup order changed")
    return raw_group, regular_group, lookup, subgroup_matrices


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
        chamber.wall_count == chamber_record["wall_count"]
        and chamber.order_aut == chamber_record["stabilizer_order"]
        and class_record["order"] > 1
        and class_record["fixed_rank"] >= MIN_FIXED_RANK
    ):
        raise ArithmeticError("the selected inventory case changed")

    ambient_gram = sp.Matrix(atlas.gram)
    dual_gram = ambient_gram.inv()
    walls_tuple = tuple(wall for orbit in chamber.wall_orbits for wall in orbit)
    walls = tuple(sp.Matrix(1, AMBIENT_RANK, wall) for wall in walls_tuple)
    if len(set(walls_tuple)) != len(walls_tuple):
        raise ArithmeticError("the published chamber has duplicate walls")
    ample_tuple = _row_product(chamber.weyl, atlas.projection_to_s_dual)
    ample_dual = sp.Matrix(1, AMBIENT_RANK, ample_tuple)

    raw_group, _, lookup, subgroup_raw = _reconstruct_inventory_subgroup(
        atlas, chamber_index, class_record
    )
    raw_index = {matrix: index for index, matrix in enumerate(raw_group)}
    dual_group = tuple(sp.Matrix(atlas.dual_action(matrix)) for matrix in raw_group)
    source = {
        "gram": dual_gram,
        "ample": ample_dual,
        "stabilizer": dual_group,
        "wall_vectors": walls_tuple,
    }
    archived_wall_group = _derive_wall_permutations(source)
    subgroup_wall: Subgroup = frozenset(
        archived_wall_group[raw_index[matrix]] for matrix in subgroup_raw
    )
    if len(subgroup_wall) != class_record["order"]:
        raise ArithmeticError("the subgroup wall action is not faithful")

    subgroup_matrices = tuple(sp.Matrix(matrix) for matrix in subgroup_raw)
    fixed = _saturated_fixed_lattice(
        subgroup_matrices, ambient_gram, int(class_record["fixed_rank"])
    )
    fixed_basis: sp.Matrix = fixed["basis"]
    fixed_gram: sp.Matrix = fixed["gram"]
    fixed_rank = fixed_basis.rows
    ample = _fixed_coordinates_of_ample(ample_dual, ambient_gram, fixed_basis)
    if (ample.T * fixed_gram * ample)[0] <= 0:
        raise ArithmeticError("the published fixed Weyl vector is not timelike")

    wall_orbits = _orbits(subgroup_wall, len(walls))
    orbit_normals = tuple(
        sum((walls[index] for index in orbit), sp.zeros(1, AMBIENT_RANK))
        for orbit in wall_orbits
    )
    pairing_matrix = sp.Matrix.hstack(
        *(fixed_basis * normal.T for normal in orbit_normals)
    )
    if pairing_matrix.rank() != fixed_rank:
        raise ArithmeticError("the restricted walls do not span the fixed dual")
    inverse_fixed_gram = fixed_gram.inv()
    orbit_gram = pairing_matrix.T * inverse_fixed_gram * pairing_matrix
    if exact_inertia(orbit_gram) != (
        1,
        fixed_rank - 1,
        len(wall_orbits) - fixed_rank,
    ):
        raise ArithmeticError("the orbit-wall Gram is not the complete fixed Lorentz span")
    coordinates = pairing_matrix.T

    conic_facets_list: list[int] = []
    conic_redundant_list: list[int] = []
    for target in range(len(wall_orbits)):
        dependence = _exact_nonnegative_dependence(coordinates, target)
        (conic_redundant_list if dependence is not None else conic_facets_list).append(target)
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
            unresolved.append(target)
            continue
        witness, record = result
        square, pairings = _audit_witness(
            witness, target, fixed_gram, pairing_matrix
        )
        record |= {
            "orbit_index": target + 1,
            "fixed_lattice_coordinates": [
                _rational_text(sp.Rational(value)) for value in witness
            ],
            "witness_square": _rational_text(square),
            "minimum_strict_pairing": _rational_text(
                min(value for index, value in enumerate(pairings) if index != target)
            ),
        }
        active_facets.append(target)
        witnesses.append(record)

    facet_gram = orbit_gram.extract(active_facets, active_facets)
    pair_gate = _pair_gate(facet_gram)
    if pair_gate["non_coxeter_pair_count"]:
        decision = "not_coxeter"
    elif unresolved:
        decision = "unresolved_facet_causality"
    else:
        decision = "coxeter_candidate_requires_volume_audit"

    orbit_records = []
    for orbit_index, (orbit, normal) in enumerate(zip(wall_orbits, orbit_normals, strict=True)):
        fixed_normal = inverse_fixed_gram * pairing_matrix[:, orbit_index]
        primitive, scale = _primitive(fixed_normal)
        orbit_records.append(
            {
                "orbit_index": orbit_index + 1,
                "source_wall_indices": [index + 1 for index in orbit],
                "ambient_dual_normal": [int(value) for value in normal],
                "fixed_pairing_column": [
                    int(value) for value in pairing_matrix[:, orbit_index]
                ],
                "fixed_normal_coordinates": [
                    _rational_text(sp.Rational(value)) for value in fixed_normal
                ],
                "primitive_fixed_normal_coordinates": [int(value) for value in primitive],
                "presentation_scale_from_primitive": _rational_text(scale),
                "square": _rational_text(sp.Rational(orbit_gram[orbit_index, orbit_index])),
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

    return {
        "schema": SCHEMA,
        "id": _case_id(chamber_index, class_index),
        "source": {
            "primary_disc16_sha256": EXPECTED_SOURCE_SHA256,
            "inventory_path": str(INVENTORY.relative_to(ROOT)),
            "inventory_sha256": inventory_hash,
        },
        "scope": {
            "group": "Shimada archived chamber stabilizer Aut(D)",
            "full_abstract_weighted_wall_automorphism_group_audited_here": False,
            "surface_automorphism_interpretation_used_as_gate": False,
            "stock_bounded_to_one_complete_published_induced_chamber": True,
        },
        "case": {
            "chamber_index": chamber_index,
            "class_index": class_index,
            "published_wall_count": chamber.wall_count,
            "archived_stabilizer_order": chamber.order_aut,
            "subgroup_order": class_record["order"],
            "subgroup_conjugacy_class_size": class_record["conjugacy_class_size"],
            "fixed_rank": fixed_rank,
            "hyperbolic_dimension": fixed_rank - 1,
            "regular_action_elements": class_record["regular_action_elements"],
            "raw_primal_row_matrices": [
                _matrix_payload(matrix) for matrix in subgroup_matrices
            ],
            "dual_row_matrices": [
                _matrix_payload(sp.Matrix(atlas.dual_action(tuple(tuple(int(v) for v in row) for row in matrix.tolist()))))
                for matrix in subgroup_matrices
            ],
        },
        "saturated_fixed_lattice": {
            "basis_primal_rows": _matrix_payload(fixed_basis),
            "gram": _matrix_payload(fixed_gram),
            "determinant": fixed["determinant"],
            "determinant_absolute": abs(fixed["determinant"]),
            "inertia": list(exact_inertia(fixed_gram)),
            "smith_diagonal_absolute": fixed["smith_diagonal_absolute"],
            "smith_constraint_diagonal": fixed["smith_constraint_diagonal"],
            "ample_fixed_coordinates": [
                _rational_text(sp.Rational(value)) for value in ample
            ],
            "ample_square": _rational_text(
                sp.Rational((ample.T * fixed_gram * ample)[0])
            ),
            "saturated_smith_kernel": True,
        },
        "folded_walls": {
            "wall_orbit_count": len(wall_orbits),
            "wall_orbit_size_histogram": {
                str(size): count
                for size, count in sorted(Counter(map(len, wall_orbits)).items())
            },
            "orbit_gram_inertia": list(exact_inertia(orbit_gram)),
            "orbit_records": orbit_records,
        },
        "exact_reduction": {
            "conic_facet_count": len(conic_facets),
            "conic_redundant_count": len(conic_redundant),
            "conic_facet_orbit_indices": [index + 1 for index in conic_facets],
            "conic_redundancy_certificates": conic_redundancy_certificates,
            "causally_inactive_conic_facet_count": len(causal_inactive),
            "causally_inactive_orbit_indices": [index + 1 for index in causal_inactive],
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
            "inventory_subgroup_reconstructed_from_regular_action": True,
            "archived_stabilizer_wall_action_faithful": True,
            "saturated_fixed_lattice_certified": True,
            "all_published_walls_folded": True,
            "orbit_normals_span_complete_fixed_dual": True,
            "all_conic_redundancies_use_final_conic_facets": True,
            "every_archived_hyperbolic_facet_witness_replayed_exactly": True,
            "every_resolved_facet_pair_checked_exactly": True,
        },
    }


def _inventory_cases(inventory: dict[str, Any]) -> tuple[tuple[dict[str, Any], dict[str, Any]], ...]:
    cases = []
    for chamber in inventory["chambers"]:
        for class_record in chamber["classes"]:
            if class_record["order"] > 1 and class_record["fixed_rank"] >= MIN_FIXED_RANK:
                cases.append((chamber, class_record))
    return tuple(
        sorted(
            cases,
            key=lambda pair: (
                pair[0]["wall_count"],
                pair[0]["chamber_index"],
                pair[1]["class_index"],
            ),
        )
    )


def _d5_skip_rows(
    cases: tuple[tuple[dict[str, Any], dict[str, Any]], ...]
) -> dict[str, dict[str, Any]]:
    d5 = json.loads(D5_CENSUS.read_text(encoding="utf-8"))
    mapping = {
        11: ("class_03_H10", "finite_volume_noncompact"),
        12: ("class_04_H11", "not_coxeter"),
        15: ("class_02_H14", "finite_volume_noncompact_known_alexeev"),
    }
    result = {}
    d5_cases = [pair for pair in cases if pair[0]["chamber_index"] == 5]
    if len(d5_cases) != 3:
        raise ArithmeticError("the inventory no longer has exactly three proper D5 classes")
    for chamber, record in d5_cases:
        key, decision = mapping[int(record["fixed_rank"])]
        target = d5["above_H9"]["complete_section_screens"][key]
        if not (
            chamber["wall_count"] == 165
            and chamber["stabilizer_order"] == 20
            and record["order"] == 2
            and target.get("hyperbolic_dimension") == record["hyperbolic_dimension"]
        ):
            raise ArithmeticError("a D5 inventory class did not map exactly")
        case_id = _case_id(chamber["chamber_index"], record["class_index"])
        result[case_id] = {
            "id": case_id,
            "chamber_index": 5,
            "class_index": record["class_index"],
            "wall_count": 165,
            "fixed_rank": record["fixed_rank"],
            "hyperbolic_dimension": record["hyperbolic_dimension"],
            "status": "skipped_exactly_mapped_to_closed_D5_class",
            "mapped_D5_key": key,
            "mapped_decision": decision,
            "mapped_certificate": str(D5_CENSUS.relative_to(ROOT)),
        }
    return result


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
    ):
        return None
    return payload


def _summary_row(payload: dict[str, Any]) -> dict[str, Any]:
    case = payload["case"]
    reduction = payload["exact_reduction"]
    gate = payload["coxeter_gate"]
    return {
        "id": payload["id"],
        "chamber_index": case["chamber_index"],
        "class_index": case["class_index"],
        "wall_count": case["published_wall_count"],
        "subgroup_order": case["subgroup_order"],
        "fixed_rank": case["fixed_rank"],
        "hyperbolic_dimension": case["hyperbolic_dimension"],
        "wall_orbit_count": payload["folded_walls"]["wall_orbit_count"],
        "conic_facet_count": reduction["conic_facet_count"],
        "hyperbolic_facet_count": reduction["hyperbolic_facet_count"],
        "unresolved_facet_count": reduction["unresolved_conic_facet_count"],
        "non_coxeter_pair_count": gate["non_coxeter_pair_count"],
        "classification": payload["decision"]["classification"],
        "case_certificate": str(
            _case_path(
                CASE_DIRECTORY, case["chamber_index"], case["class_index"]
            ).relative_to(ROOT)
        ),
    }


def _build_manifest(
    inventory_hash: str,
    cases: tuple[tuple[dict[str, Any], dict[str, Any]], ...],
    case_directory: Path,
) -> dict[str, Any]:
    skips = _d5_skip_rows(cases)
    rows = []
    for chamber, record in cases:
        case_id = _case_id(chamber["chamber_index"], record["class_index"])
        if case_id in skips:
            rows.append(skips[case_id])
            continue
        path = _case_path(
            case_directory, chamber["chamber_index"], record["class_index"]
        )
        payload = _valid_case_artifact(path, inventory_hash)
        if payload is None:
            rows.append(
                {
                    "id": case_id,
                    "chamber_index": chamber["chamber_index"],
                    "class_index": record["class_index"],
                    "wall_count": chamber["wall_count"],
                    "fixed_rank": record["fixed_rank"],
                    "hyperbolic_dimension": record["hyperbolic_dimension"],
                    "status": "pending",
                }
            )
        else:
            rows.append(_summary_row(payload) | {"status": "complete"})
    completed = [row for row in rows if row["status"] == "complete"]
    pending = [row for row in rows if row["status"] == "pending"]
    skipped = [row for row in rows if row["status"].startswith("skipped")]
    survivors = [
        row
        for row in completed
        if row["classification"] == "coxeter_candidate_requires_volume_audit"
    ]
    return {
        "schema": MANIFEST_SCHEMA,
        "source": {
            "primary_disc16_sha256": EXPECTED_SOURCE_SHA256,
            "inventory_path": str(INVENTORY.relative_to(ROOT)),
            "inventory_sha256": inventory_hash,
        },
        "scope": {
            "denominator": "107 proper above-H9 subgroup conjugacy classes of Shimada archived chamber stabilizers Aut(D)",
            "full_weighted_wall_graph_automorphism_groups": "separate audit; not claimed complete here",
            "surface_automorphism_interpretation_used_as_gate": False,
            "case_order": "increasing published wall count, then chamber and class index",
            "resumable_per_case_artifacts": True,
            "finite_volume_gate_performed": False,
        },
        "summary": {
            "inventory_case_count": len(cases),
            "exactly_mapped_closed_D5_case_count": len(skipped),
            "remaining_case_count": len(cases) - len(skipped),
            "completed_remaining_case_count": len(completed),
            "pending_remaining_case_count": len(pending),
            "classification_histogram": dict(
                sorted(Counter(row["classification"] for row in completed).items())
            ),
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
    arguments = parser.parse_args()

    if _sha256(arguments.source) != EXPECTED_SOURCE_SHA256:
        raise SystemExit("the primary disc16 source hash changed")
    inventory_hash = _sha256(arguments.inventory)
    inventory = json.loads(arguments.inventory.read_text(encoding="utf-8"))
    if not (
        inventory["source"]["sha256"] == EXPECTED_SOURCE_SHA256
        and inventory["summary"]["chamber_count"] == 4539
    ):
        raise SystemExit("the frozen high-section inventory changed")
    all_cases = _inventory_cases(inventory)
    if len(all_cases) != 107:
        raise SystemExit(f"expected 107 proper high-rank cases, got {len(all_cases)}")
    requested = {
        tuple(map(int, text.split(":"))) for text in arguments.case
    }
    if any(len(pair) != 2 for pair in requested):
        raise SystemExit("--case must have CHAMBER:CLASS form")

    arguments.case_directory.mkdir(parents=True, exist_ok=True)
    atlas = ShimadaDisc11Atlas(arguments.source)
    processed = 0
    for chamber, record in all_cases:
        chamber_index = int(chamber["chamber_index"])
        class_index = int(record["class_index"])
        if chamber_index == 5:
            continue
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
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        row = _summary_row(payload)
        print(json.dumps(row, sort_keys=True), flush=True)
        if row["classification"] == "coxeter_candidate_requires_volume_audit":
            print("SURVIVOR " + json.dumps(row, sort_keys=True), flush=True)
        processed += 1
        manifest = _build_manifest(inventory_hash, all_cases, arguments.case_directory)
        arguments.manifest.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    manifest = _build_manifest(inventory_hash, all_cases, arguments.case_directory)
    arguments.manifest.parent.mkdir(parents=True, exist_ok=True)
    arguments.manifest.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
