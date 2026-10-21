#!/usr/bin/env python3
"""Raw lattice-to-adjunction Stage 0 for rational AG carriers.

The mathematical input is only an integral unimodular Lorentzian Gram
matrix and a full-span list of primitive oriented negative wall vectors.
No canonical class, physical scale, boundary set, surface label, or carrier
verdict is accepted by :func:`solve_raw_input`.

For a putative smooth rational curve ``C=d*r`` the adjunction equation and
characteristic parity leave ``d`` in ``{1,2}``.  A minimum-ambiguity wall
basis makes the resulting finite enumeration at most ``2**rho``.  Every
integral characteristic solution of square ``10-rho`` is then checked
against every wall.  The output distinguishes the direct NS-integral weak
branch from a standard log-boundary branch on its inferred adapted module.
"""

from __future__ import annotations

import argparse
import copy
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "data/fixtures/raw_adjunction_log_stage0_inputs.json"
ARCHIVE = ROOT / "data/enumerations/raw_adjunction_log_stage0.json"
PB_FIXTURES = ROOT / "data/fixtures/pb_atlas_numerical_stage_a_inputs.json"
TORIC_SOURCE = ROOT / "data/enumerations/toric_full_nef_coxeter_factory.json"
F4_SOURCE = ROOT / "data/enumerations/rational_nef_reflection_integrality_gate.json"

INPUT_SCHEMA = "coxeter.raw_adjunction_log_stage0.input.v1"
FIXTURE_SCHEMA = "coxeter.raw_adjunction_log_stage0.fixtures.v1"
OUTPUT_SCHEMA = "coxeter.raw_adjunction_log_stage0.certificate.v1"
CALIBRATION_SCHEMA = "coxeter.raw_adjunction_log_stage0.calibration.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _fraction(value: Any) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, str) and "/" in value:
        numerator, denominator = value.split("/", 1)
        return Fraction(int(numerator), int(denominator))
    return Fraction(value)


def _serialize(value: Any) -> Any:
    if isinstance(value, Fraction):
        return (
            int(value)
            if value.denominator == 1
            else f"{value.numerator}/{value.denominator}"
        )
    if isinstance(value, sp.Rational):
        return _serialize(Fraction(int(value.p), int(value.q)))
    if isinstance(value, sp.Integer):
        return int(value)
    if isinstance(value, sp.MatrixBase):
        return [[_serialize(value[i, j]) for j in range(value.cols)] for i in range(value.rows)]
    if isinstance(value, tuple):
        return [_serialize(item) for item in value]
    if isinstance(value, list):
        return [_serialize(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _serialize(item) for key, item in value.items()}
    return value


def _vector_list(vector: sp.Matrix) -> list[Any]:
    if vector.cols != 1:
        raise ValueError("expected a column vector")
    return [_serialize(vector[index]) for index in range(vector.rows)]


def _dot(g: sp.Matrix, left: sp.Matrix, right: sp.Matrix) -> sp.Rational:
    return sp.Rational((left.T * g * right)[0])


def _integer_vector(vector: sp.Matrix) -> bool:
    return all(sp.Rational(entry).q == 1 for entry in vector)


def _primitive(vector: Sequence[int]) -> bool:
    return math.gcd(*(abs(int(entry)) for entry in vector)) == 1


def _inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Exact inertia by rational symmetric congruence elimination."""

    a = sp.Matrix(matrix)
    positive = negative = zero = 0
    while a.rows:
        pivot = next((i for i in range(a.rows) if a[i, i] != 0), None)
        if pivot is None:
            pair = next(
                (
                    (i, j)
                    for i in range(a.rows)
                    for j in range(i + 1, a.rows)
                    if a[i, j] != 0
                ),
                None,
            )
            if pair is None:
                zero += a.rows
                break
            i, j = pair
            transform = sp.eye(a.rows)
            transform[:, 0], transform[:, 1] = (
                sp.eye(a.rows)[:, i] + sp.eye(a.rows)[:, j],
                sp.eye(a.rows)[:, i] - sp.eye(a.rows)[:, j],
            )
            used = {i, j}
            remaining = [k for k in range(a.rows) if k not in used]
            for column, old in enumerate(remaining, start=2):
                transform[:, column] = sp.eye(a.rows)[:, old]
            a = transform.T * a * transform
            pivot = 0
        if pivot != 0:
            permutation = [pivot, *[i for i in range(a.rows) if i != pivot]]
            a = a.extract(permutation, permutation)
        value = sp.Rational(a[0, 0])
        if value > 0:
            positive += 1
        else:
            negative += 1
        if a.rows == 1:
            break
        column = a[1:, 0]
        a = a[1:, 1:] - column * column.T / value
    return positive, negative, zero


def _reflection(g: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = _dot(g, root, root)
    return sp.eye(g.rows) - sp.Rational(2, 1) * root * (root.T * g) / square


def _matrix_integral(matrix: sp.Matrix) -> bool:
    return all(sp.Rational(entry).q == 1 for entry in matrix)


def _minimum_ambiguity_basis(walls: Sequence[sp.Matrix], squares: Sequence[int]) -> list[int]:
    """Return a minimum-weight vector-matroid basis.

    A wall has weight one exactly when its primitive square is odd, because
    only then can ``d=2`` survive characteristic parity.  Greedy selection in
    nondecreasing weight is the matroid greedy algorithm, so the number of
    binary choices is minimal among all wall bases.
    """

    rank = walls[0].rows
    chosen: list[int] = []
    span = sp.zeros(rank, 0)
    old_rank = 0
    order = sorted(range(len(walls)), key=lambda index: (squares[index] % 2, index))
    for index in order:
        trial = span.row_join(walls[index])
        new_rank = trial.rank()
        if new_rank > old_rank:
            chosen.append(index)
            span = trial
            old_rank = new_rank
        if old_rank == rank:
            return chosen
    return chosen


def _allowed_scale_pairings(q: int) -> tuple[tuple[int, int], ...]:
    """Return ``(d, K.r)`` pairs surviving adjunction and parity."""

    rows = [(1, q - 2)]
    if q % 2 == 1:
        rows.append((2, 2 * q - 1))
    return tuple(rows)


def _characteristic(g: sp.Matrix, canonical: sp.Matrix) -> bool:
    covector = g * canonical
    return all((int(covector[i]) - int(g[i, i])) % 2 == 0 for i in range(g.rows))


def _lcm(values: Iterable[int]) -> int:
    result = 1
    for value in values:
        result = math.lcm(result, int(value))
    return result


def _adapted_module_certificate(
    g: sp.Matrix,
    walls: Sequence[sp.Matrix],
    physical: Sequence[sp.Matrix],
    boundary_indices: Sequence[int],
    ramification_indices: Sequence[int],
) -> dict[str, Any]:
    """Check all wall reflections on ``L + sum Z*(C_i/m_i)`` exactly."""

    rank = g.rows
    boundary_generators = [
        physical[index] / ramification_indices[position]
        for position, index in enumerate(boundary_indices)
    ]
    generators = [sp.eye(rank)[:, index] for index in range(rank)] + boundary_generators
    denominators = [sp.Rational(entry).q for vector in generators for entry in vector]
    common_denominator = _lcm(denominators)
    integer_columns = [
        sp.Matrix([int(sp.Rational(entry) * common_denominator) for entry in vector])
        for vector in generators
    ]
    hnf = hermite_normal_form(sp.Matrix.hstack(*integer_columns))
    if hnf.rows != rank or hnf.cols != rank or hnf.det() == 0:
        raise ArithmeticError("the inferred adapted module lost full rank")
    hnf_inverse = hnf.inv()

    def member(vector: sp.Matrix) -> bool:
        scaled = vector * common_denominator
        if not _integer_vector(scaled):
            return False
        coordinates = hnf_inverse * scaled
        return _integer_vector(coordinates)

    reflection_rows = []
    for index, wall in enumerate(walls):
        matrix = _reflection(g, wall)
        preserves = all(member(matrix * generator) for generator in generators)
        reflection_rows.append(
            {
                "wall_index_zero_based": index,
                "integral_on_NS": _matrix_integral(matrix),
                "preserves_inferred_adapted_module": preserves,
            }
        )

    ell = _lcm(ramification_indices) if ramification_indices else 1
    scaled_pairings_integral = all(
        sp.Rational(ell) * _dot(g, left, right) in sp.ZZ
        for left in generators
        for right in generators
    )
    return {
        "module_generators_after_NS_basis": [
            _vector_list(generator) for generator in boundary_generators
        ],
        "common_coordinate_denominator": common_denominator,
        "integer_HNF_at_common_denominator": hnf,
        "index_of_scaled_module_in_Zrho": abs(int(hnf.det())),
        "adapted_cover_degree_lcm": ell,
        "ell_scaled_form_is_integral": scaled_pairings_integral,
        "reflection_checks": reflection_rows,
        "all_wall_reflections_preserve_inferred_adapted_module": all(
            row["preserves_inferred_adapted_module"] for row in reflection_rows
        ),
    }


def _candidate_certificate(
    g: sp.Matrix,
    walls: Sequence[sp.Matrix],
    primitive_squares: Sequence[int],
    canonical: sp.Matrix,
    scales: Sequence[int],
) -> dict[str, Any]:
    rank = g.rows
    physical = [scale * wall for scale, wall in zip(scales, walls)]
    physical_squares = [int(_dot(g, curve, curve)) for curve in physical]
    pair_gram = sp.Matrix(
        [[_dot(g, left, right) for right in physical] for left in physical]
    )
    weak = all(square in (-1, -2) for square in physical_squares)
    boundary_indices = [
        index for index, square in enumerate(physical_squares) if square not in (-1, -2)
    ]
    standard_boundary_shapes = bool(boundary_indices) and all(
        physical_squares[index] <= -4 and physical_squares[index] % 2 == 0
        for index in boundary_indices
    )
    ramification_indices = (
        [-physical_squares[index] // 2 for index in boundary_indices]
        if standard_boundary_shapes
        else []
    )
    delta = sp.zeros(rank, 1)
    if standard_boundary_shapes:
        for position, index in enumerate(boundary_indices):
            m = ramification_indices[position]
            delta += sp.Rational(m - 1, m) * physical[index]
    log_residual = canonical + delta
    disjoint = standard_boundary_shapes and all(
        _dot(g, physical[left], physical[right]) == 0
        for offset, left in enumerate(boundary_indices)
        for right in boundary_indices[offset + 1 :]
    )
    module = _adapted_module_certificate(
        g, walls, physical, boundary_indices, ramification_indices
    )
    ns_integral = all(row["integral_on_NS"] for row in module["reflection_checks"])
    log_equation = standard_boundary_shapes and all(entry == 0 for entry in log_residual)
    standard_log = (
        standard_boundary_shapes
        and disjoint
        and log_equation
        and module["ell_scaled_form_is_integral"]
        and module["all_wall_reflections_preserve_inferred_adapted_module"]
    )
    weak_accepted = weak and ns_integral
    accepted = weak_accepted or standard_log
    if weak_accepted:
        branch = "weak_minus_one_minus_two"
    elif standard_log and all(index == 2 for index in ramification_indices):
        branch = "coble_index_two_log_boundary"
    elif standard_log:
        branch = "standard_log_boundary_tame_cover_obstructed_unless_index_two"
    else:
        branch = "adjunction_solution_rejected_by_branch_gate"

    boundary_set = set(boundary_indices)
    incidence_rows = []
    if standard_boundary_shapes:
        for index, curve in enumerate(physical):
            if index in boundary_set:
                continue
            incidence_rows.append(
                {
                    "wall_index_zero_based": index,
                    "physical_square": physical_squares[index],
                    "total_boundary_intersection": sum(
                        _dot(g, physical[boundary], curve) for boundary in boundary_indices
                    ),
                }
            )

    square_histogram = {
        str(square): physical_squares.count(square) for square in sorted(set(physical_squares))
    }
    return {
        "canonical_class": _vector_list(canonical),
        "canonical_square": _dot(g, canonical, canonical),
        "physical_scales": list(scales),
        "physical_squares": physical_squares,
        "physical_square_histogram": square_histogram,
        "physical_wall_gram": pair_gram,
        "weak_branch": weak_accepted,
        "boundary_indices_zero_based_inferred": boundary_indices,
        "ramification_indices_inferred": ramification_indices,
        "standard_boundary_shapes": standard_boundary_shapes,
        "boundary_pairwise_disjoint": disjoint,
        "log_delta": _vector_list(delta),
        "log_residual_K_plus_Delta": _vector_list(log_residual),
        "log_equation": log_equation,
        "other_wall_boundary_incidence": incidence_rows,
        "reflection_module": module,
        "all_reflections_integral_on_NS": ns_integral,
        "standard_log_branch": standard_log,
        "accepted": accepted,
        "branch": branch,
        "tame_adapted_cover_consequence": (
            "K3 double cover and Picard ceiling 20 in characteristic zero / 22 in odd characteristic"
            if standard_log and all(index == 2 for index in ramification_indices)
            else (
                "the adapted-cover theorem excludes an actual tame rational-surface carrier"
                if standard_log
                else "not applicable"
            )
        ),
    }


def _invalid_certificate(reason: str, validation: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema": OUTPUT_SCHEMA,
        "status": "INVALID_INPUT",
        "reason": reason,
        "validation": dict(validation),
        "canonical_candidates": [],
        "accepted_candidates": [],
    }


def solve_raw_input(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Exhaust the raw numerical Stage-0 candidates for one input."""

    allowed_keys = {"schema", "picard_gram", "wall_vectors"}
    extra = sorted(set(payload) - allowed_keys)
    missing = sorted(allowed_keys - set(payload))
    validation: dict[str, Any] = {
        "input_keys_exact": not extra and not missing,
        "extra_keys": extra,
        "missing_keys": missing,
    }
    if extra or missing:
        return _invalid_certificate("the raw input key set is not exact", validation)
    if payload["schema"] != INPUT_SCHEMA:
        return _invalid_certificate("unknown input schema", validation)
    try:
        gram_rows = payload["picard_gram"]
        wall_rows = payload["wall_vectors"]
        if not gram_rows or any(len(row) != len(gram_rows) for row in gram_rows):
            raise ValueError("the Gram matrix is not nonempty and square")
        rank = len(gram_rows)
        if any(not isinstance(entry, int) or isinstance(entry, bool) for row in gram_rows for entry in row):
            raise ValueError("the Gram matrix is not integral")
        if not wall_rows or any(len(row) != rank for row in wall_rows):
            raise ValueError("the wall array has the wrong shape")
        if any(not isinstance(entry, int) or isinstance(entry, bool) for row in wall_rows for entry in row):
            raise ValueError("the wall array is not integral")
        g = sp.Matrix(gram_rows)
        walls = [sp.Matrix(row) for row in wall_rows]
    except (TypeError, ValueError) as error:
        return _invalid_certificate(str(error), validation)

    determinant = int(g.det())
    inertia = _inertia(g) if g == g.T else None
    primitive = [_primitive(row) for row in wall_rows]
    squares = [int(-_dot(g, wall, wall)) for wall in walls]
    wall_rank = sp.Matrix.hstack(*walls).rank()
    oriented_tuples = [tuple(row) for row in wall_rows]
    distinct_vectors = len(set(oriented_tuples)) == len(oriented_tuples)
    unoriented_keys = [min(vector, tuple(-entry for entry in vector)) for vector in oriented_tuples]
    distinct_rays = len(set(unoriented_keys)) == len(unoriented_keys)
    off_diagonal_nonnegative = all(
        _dot(g, walls[i], walls[j]) >= 0
        for i in range(len(walls))
        for j in range(i + 1, len(walls))
    )
    validation.update(
        {
            "rank": rank,
            "wall_count": len(walls),
            "symmetric": g == g.T,
            "determinant": determinant,
            "unimodular": abs(determinant) == 1,
            "inertia": inertia,
            "lorentzian_signature_one_positive": inertia == (1, rank - 1, 0),
            "every_wall_primitive": all(primitive),
            "primitive_flags": primitive,
            "every_wall_negative": all(square > 0 for square in squares),
            "primitive_negative_square_magnitudes": squares,
            "wall_span_rank": wall_rank,
            "full_span": wall_rank == rank,
            "distinct_oriented_vectors": distinct_vectors,
            "distinct_unoriented_wall_rays": distinct_rays,
            "pairwise_nonnegative_oriented_intersections": off_diagonal_nonnegative,
        }
    )
    required = (
        validation["symmetric"]
        and validation["unimodular"]
        and validation["lorentzian_signature_one_positive"]
        and validation["every_wall_primitive"]
        and validation["every_wall_negative"]
        and validation["full_span"]
        and validation["distinct_oriented_vectors"]
        and validation["distinct_unoriented_wall_rays"]
        and validation["pairwise_nonnegative_oriented_intersections"]
    )
    if not required:
        return _invalid_certificate("a raw lattice/wall validation gate failed", validation)

    basis_indices = _minimum_ambiguity_basis(walls, squares)
    if len(basis_indices) != rank:
        return _invalid_certificate("no full-rank wall basis exists", validation)
    covectors = sp.Matrix([list(walls[index].T * g) for index in basis_indices])
    inverse = covectors.inv()
    basis_options = [_allowed_scale_pairings(squares[index]) for index in basis_indices]
    enumeration_bound = math.prod(len(options) for options in basis_options)

    candidates_by_canonical: dict[tuple[int, ...], dict[str, Any]] = {}
    raw_solution_count = integral_solution_count = 0
    for choices in itertools.product(*basis_options):
        pairing_vector = sp.Matrix([pairing for _, pairing in choices])
        canonical = inverse * pairing_vector
        raw_solution_count += 1
        if not _integer_vector(canonical):
            continue
        integral_solution_count += 1
        scales: list[int] = []
        compatible = True
        for q, wall in zip(squares, walls):
            value = int(_dot(g, canonical, wall))
            matching = [scale for scale, pairing in _allowed_scale_pairings(q) if pairing == value]
            if len(matching) != 1:
                compatible = False
                break
            scales.append(matching[0])
        if not compatible:
            continue
        if not _characteristic(g, canonical):
            continue
        if _dot(g, canonical, canonical) != 10 - rank:
            continue
        key = tuple(int(entry) for entry in canonical)
        candidates_by_canonical[key] = _candidate_certificate(
            g, walls, squares, canonical, scales
        )

    candidates = [candidates_by_canonical[key] for key in sorted(candidates_by_canonical)]
    accepted = [candidate for candidate in candidates if candidate["accepted"]]
    branch_histogram: dict[str, int] = {}
    for candidate in candidates:
        branch_histogram[candidate["branch"]] = branch_histogram.get(candidate["branch"], 0) + 1

    result = {
        "schema": OUTPUT_SCHEMA,
        "status": "PASS",
        "input": {
            "schema": INPUT_SCHEMA,
            "picard_gram": [list(row) for row in payload["picard_gram"]],
            "wall_vectors": [list(row) for row in payload["wall_vectors"]],
        },
        "validation": validation,
        "finite_enumeration": {
            "basis_wall_indices_zero_based": basis_indices,
            "ambiguous_basis_wall_count": sum(squares[index] % 2 for index in basis_indices),
            "enumeration_upper_bound": enumeration_bound,
            "basis_pairing_assignments_checked": raw_solution_count,
            "integral_basis_solutions": integral_solution_count,
            "canonical_candidates_after_all_gates": len(candidates),
            "proof": (
                "full span makes K unique from a wall basis; adjunction forces d|2, "
                "and characteristic parity permits d=2 only for odd -r^2"
            ),
        },
        "canonical_candidates": candidates,
        "accepted_candidates": accepted,
        "accepted_candidate_count": len(accepted),
        "branch_histogram": branch_histogram,
        "scope": {
            "proved": [
                "finite exhaustive numerical scaling and canonical-class enumeration",
                "weak NS-integral branch detection",
                "standard disjoint log-boundary equation and inferred-module reflection test",
            ],
            "not_proved": [
                "effectivity or irreducibility of a wall class",
                "existence or nonemptiness of a rational surface",
                "completeness of the Mori cone or absence of hidden walls",
                "same-polytope return or finite volume from the numerical data alone",
            ],
        },
    }
    return _serialize(result)


def _raw_projection(source_input: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema": INPUT_SCHEMA,
        "picard_gram": copy.deepcopy(source_input["picard_gram"]),
        "wall_vectors": copy.deepcopy(source_input["wall_vectors"]),
    }


def build_fixtures_from_sources() -> dict[str, Any]:
    pb = _load(PB_FIXTURES)
    toric = _load(TORIC_SOURCE)["h6_endpoint"]
    f4 = _load(F4_SOURCE)["counterexample"]
    aliases = {case["calibration_alias"]: case for case in pb["cases"]}
    cases = []
    specifications = [
        ("coble_type_ii_h10_30", "H10/30", "coble_index_two_log_boundary", 1),
        ("coble_type_i_h11_26", "H11/26", "coble_index_two_log_boundary", 1),
        ("coble_type_vi_char3_h14_60", "H14/60", "coble_index_two_log_boundary", 1),
        ("coble_type_mii_char3_h17_48", "H17/48", "coble_index_two_log_boundary", 1),
    ]
    for case_id, alias, branch, accepted_count in specifications:
        source = aliases[alias]
        cases.append(
            {
                "case_id": case_id,
                "raw_input": _raw_projection(source["input"]),
                "expected": {
                    "accepted_branch": branch,
                    "accepted_candidate_count": accepted_count,
                    "canonical_class": source["input"]["canonical_class"],
                },
            }
        )
    cases.append(
        {
            "case_id": "toric_h6_9",
            "raw_input": {
                "schema": INPUT_SCHEMA,
                "picard_gram": toric["intersection_matrix"],
                "wall_vectors": [toric["walls"][name] for name in toric["wall_order"]],
            },
            "expected": {
                "accepted_branch": "weak_minus_one_minus_two",
                "accepted_candidate_count": 1,
                "canonical_class": toric["canonical_class"],
            },
        }
    )
    cases.append(
        {
            "case_id": "blowup_F4_maximal_triangle_negative_control",
            "raw_input": {
                "schema": INPUT_SCHEMA,
                "picard_gram": f4["NS_gram"],
                "wall_vectors": [
                    f4["Mori_generators"][name] for name in ("B", "E", "G=F-E")
                ],
            },
            "expected": {
                "accepted_branch": None,
                "accepted_candidate_count": 0,
                "canonical_class": f4["canonical_coordinates"],
                "rejection_gate": "nonzero K+Delta",
            },
        }
    )
    return {
        "schema": FIXTURE_SCHEMA,
        "provenance": {
            str(PB_FIXTURES.relative_to(ROOT)): _sha256(PB_FIXTURES),
            str(TORIC_SOURCE.relative_to(ROOT)): _sha256(TORIC_SOURCE),
            str(F4_SOURCE.relative_to(ROOT)): _sha256(F4_SOURCE),
        },
        "cases": cases,
    }


def _invariant_summary(certificate: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "accepted_candidate_count": certificate["accepted_candidate_count"],
        "branch_histogram": certificate["branch_histogram"],
        "candidates": sorted(
            (
                candidate["branch"],
                tuple(candidate["physical_scales"]),
                tuple(sorted(candidate["physical_squares"])),
                tuple(candidate["ramification_indices_inferred"]),
            )
            for candidate in certificate["canonical_candidates"]
        ),
    }


def _basis_mutation(payload: Mapping[str, Any]) -> tuple[dict[str, Any], sp.Matrix]:
    rank = len(payload["picard_gram"])
    change = sp.eye(rank)
    change[0, 1] = 1
    change[rank - 1, 0] = -1
    if abs(int(change.det())) != 1:
        raise ArithmeticError("the deterministic mutation is not unimodular")
    inverse = change.inv()
    gram = sp.Matrix(payload["picard_gram"])
    walls = [sp.Matrix(row) for row in payload["wall_vectors"]]
    mutated = {
        "schema": INPUT_SCHEMA,
        "picard_gram": _serialize(change.T * gram * change),
        "wall_vectors": [_vector_list(inverse * wall) for wall in walls],
    }
    return mutated, inverse


def _wall_reversal(payload: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema": INPUT_SCHEMA,
        "picard_gram": copy.deepcopy(payload["picard_gram"]),
        "wall_vectors": list(reversed(copy.deepcopy(payload["wall_vectors"]))),
    }


def build_calibration() -> dict[str, Any]:
    fixtures = _load(FIXTURES)
    if fixtures["schema"] != FIXTURE_SCHEMA:
        raise ArithmeticError("the raw Stage-0 fixture schema changed")
    cases = []
    for case in fixtures["cases"]:
        raw = case["raw_input"]
        certificate = solve_raw_input(raw)
        mutated_payload, inverse = _basis_mutation(raw)
        mutated = solve_raw_input(mutated_payload)
        reversed_certificate = solve_raw_input(_wall_reversal(raw))
        expected = case["expected"]
        recovered = [
            candidate
            for candidate in certificate["canonical_candidates"]
            if candidate["canonical_class"] == expected["canonical_class"]
        ]
        accepted_branches = [candidate["branch"] for candidate in certificate["accepted_candidates"]]

        # Transport mutated K-vectors back to the original basis.  If new
        # coordinates are k', then k'=P^{-1}k, so compare P*k'.
        change = inverse.inv()
        transported_mutated_k = sorted(
            tuple(int(entry) for entry in change * sp.Matrix(candidate["canonical_class"]))
            for candidate in mutated["canonical_candidates"]
        )
        original_k = sorted(
            tuple(candidate["canonical_class"])
            for candidate in certificate["canonical_candidates"]
        )

        checks = {
            "raw_input_keys_exact": set(raw) == {"schema", "picard_gram", "wall_vectors"},
            "solver_passes": certificate["status"] == "PASS",
            "expected_canonical_recovered_once": len(recovered) == 1,
            "expected_accepted_count": certificate["accepted_candidate_count"]
            == expected["accepted_candidate_count"],
            "expected_accepted_branch": (
                accepted_branches == [expected["accepted_branch"]]
                if expected["accepted_branch"] is not None
                else accepted_branches == []
            ),
            "basis_mutation_preserves_canonical_set": transported_mutated_k == original_k,
            "basis_mutation_preserves_abstract_summary": _invariant_summary(mutated)
            == _invariant_summary(certificate),
            "wall_permutation_preserves_candidate_count_and_branches": (
                reversed_certificate["accepted_candidate_count"]
                == certificate["accepted_candidate_count"]
                and reversed_certificate["branch_histogram"] == certificate["branch_histogram"]
            ),
            "finite_search_respects_2_power_rank_bound": certificate["finite_enumeration"][
                "enumeration_upper_bound"
            ]
            <= 2 ** certificate["validation"]["rank"],
        }
        if case["case_id"] == "blowup_F4_maximal_triangle_negative_control":
            checks.update(
                {
                    "F4_unique_adjunction_solution": len(certificate["canonical_candidates"]) == 1,
                    "F4_inferred_boundary_is_square_minus_four": recovered[0][
                        "physical_square_histogram"
                    ]
                    == {"-4": 1, "-1": 2},
                    "F4_standard_extension_reflections_integral": recovered[0][
                        "reflection_module"
                    ]["all_wall_reflections_preserve_inferred_adapted_module"],
                    "F4_standard_extension_scaled_form_integral": recovered[0][
                        "reflection_module"
                    ]["ell_scaled_form_is_integral"],
                    "F4_rejected_exactly_at_log_equation_after_shape_and_disjointness": (
                        recovered[0]["standard_boundary_shapes"]
                        and recovered[0]["boundary_pairwise_disjoint"]
                        and not recovered[0]["log_equation"]
                        and recovered[0]["log_residual_K_plus_Delta"] != [0, 0, 0]
                    ),
                }
            )
        cases.append(
            {
                "case_id": case["case_id"],
                "expected": expected,
                "certificate": certificate,
                "mutation_checks": checks,
                "all_checks_pass": all(checks.values()),
            }
        )
    result = {
        "schema": CALIBRATION_SCHEMA,
        "status": "PASS" if all(case["all_checks_pass"] for case in cases) else "FAIL",
        "fixture_sha256": _sha256(FIXTURES),
        "cases": cases,
        "all_checks_pass": all(case["all_checks_pass"] for case in cases),
    }
    if not result["all_checks_pass"]:
        failed = {
            case["case_id"]: [key for key, value in case["mutation_checks"].items() if not value]
            for case in cases
            if not case["all_checks_pass"]
        }
        raise ArithmeticError(f"raw Stage-0 calibration failed: {failed}")
    return _serialize(result)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, help="solve one raw input JSON")
    parser.add_argument("--write-fixtures", action="store_true")
    parser.add_argument("--check-fixtures", action="store_true")
    parser.add_argument("--write-archive", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()

    if args.write_fixtures or args.check_fixtures:
        fixtures = build_fixtures_from_sources()
        encoded = json.dumps(_serialize(fixtures), indent=2, sort_keys=True) + "\n"
        if args.write_fixtures:
            FIXTURES.parent.mkdir(parents=True, exist_ok=True)
            FIXTURES.write_text(encoded, encoding="utf-8")
        if args.check_fixtures:
            assert FIXTURES.read_text(encoding="utf-8") == encoded, "fixture projection differs"

    if args.input:
        payload = _load(args.input)
        certificate = solve_raw_input(payload)
    else:
        certificate = build_calibration()
    encoded = json.dumps(_serialize(certificate), indent=2, sort_keys=True) + "\n"
    if args.write_archive:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        ARCHIVE.write_text(encoded, encoding="utf-8")
    if args.check_archive:
        assert ARCHIVE.read_text(encoding="utf-8") == encoded, "archive differs from rebuild"
    print(encoded, end="")


if __name__ == "__main__":
    main()
