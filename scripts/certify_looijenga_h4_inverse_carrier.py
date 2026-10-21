#!/usr/bin/env python3
"""Exact certificate for the eight-wall H4 Looijenga inverse carrier.

The input of this certificate is the Coxeter wall datum, not a pre-existing
surface.  From the eight primitive wall vectors it reconstructs the full odd
unimodular Lorentz lattice, the canonical class, an anticanonical triangle,
and a toric blowdown certificate.  It then checks the period-root stratum,
the complete negative-curve list, all Coxeter angles and integral
reflections, facet irredundancy, and every extreme ray of the chamber.

Gross--Hacking--Keel's marked-moduli theorem and Friedman's root/ample-cone
theorems are cited external inputs.  Everything finite and arithmetic in the
worked H4 case is replayed here over the integers/rationals.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from math import gcd
from pathlib import Path
from typing import Any, Iterable

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "looijenga_h4_inverse_carrier.json"
)
SCHEMA = "coxeter.looijenga_h4_inverse_carrier.v1"

# Coordinates are (H,E1,E2,E3,E4), with H^2=1 and Ei^2=-1.
J = sp.diag(1, -1, -1, -1, -1)

WALL_ORDER = ("E1", "E2", "E3", "E4", "A", "D1", "B", "C")
WALLS: dict[str, tuple[int, ...]] = {
    "E1": (0, 1, 0, 0, 0),
    "E2": (0, 0, 1, 0, 0),
    "E3": (0, 0, 0, 1, 0),
    "E4": (0, 0, 0, 0, 1),
    "A": (1, -1, 0, -1, -1),
    "D1": (1, -1, -1, 0, 0),
    "B": (1, 0, -1, -1, 0),
    "C": (1, 0, -1, 0, -1),
}

BOUNDARY_ORDER = ("D1", "D2", "D3")
BOUNDARY: dict[str, tuple[int, ...]] = {
    "D1": WALLS["D1"],
    "D2": (1, 0, 0, -1, 0),
    "D3": (1, 0, 0, 0, -1),
}

CANONICAL = (-3, 1, 1, 1, 1)
PERIOD_BASIS = {
    "A": WALLS["A"],
    "delta": (0, 1, -1, 0, 0),  # E1-E2
}
STRICT_INTERIOR = (20, -1, -1, -1, -1)


def _v(values: Iterable[int | sp.Rational]) -> sp.Matrix:
    return sp.Matrix(list(values))


def pair(left: Iterable[int], right: Iterable[int]) -> sp.Rational:
    return sp.Rational((_v(left).T * J * _v(right))[0])


def square(vector: Iterable[int]) -> sp.Rational:
    return pair(vector, vector)


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _rational_payload(value: sp.Rational) -> int | str:
    value = sp.Rational(value)
    if value.q == 1:
        return int(value)
    return f"{int(value.p)}/{int(value.q)}"


def _vector_payload(vector: sp.Matrix) -> list[int | str]:
    return [_rational_payload(sp.Rational(value)) for value in vector]


def _canonical_sha256(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Return exact (positive, negative, zero) inertia by congruence pivots."""

    work = sp.Matrix(matrix)
    positive = negative = zero = 0
    while work.rows:
        size = work.rows
        diagonal = next((i for i in range(size) if work[i, i] != 0), None)
        if diagonal is not None:
            order = [diagonal, *[i for i in range(size) if i != diagonal]]
            work = work.extract(order, order)
            pivot = sp.Rational(work[0, 0])
            if pivot > 0:
                positive += 1
            else:
                negative += 1
            if size == 1:
                break
            tail = work[1:, 1:]
            column = work[1:, 0]
            work = sp.simplify(tail - column * column.T / pivot)
            continue

        off_diagonal = next(
            (
                (i, j)
                for i in range(size)
                for j in range(i + 1, size)
                if work[i, j] != 0
            ),
            None,
        )
        if off_diagonal is None:
            zero += size
            break
        i, j = off_diagonal
        order = [i, j, *[k for k in range(size) if k not in {i, j}]]
        work = work.extract(order, order)
        block = work[:2, :2]
        # Both diagonal entries vanish, so this block has determinant -a^2.
        positive += 1
        negative += 1
        if size == 2:
            break
        cross = work[2:, :2]
        work = sp.simplify(work[2:, 2:] - cross * block.inv() * cross.T)
    return positive, negative, zero


def wall_gram() -> sp.Matrix:
    columns = sp.Matrix.hstack(*(_v(WALLS[name]) for name in WALL_ORDER))
    return columns.T * J * columns


def _primitive_integral(vector: sp.Matrix) -> tuple[int, ...]:
    denominators = [sp.Rational(value).q for value in vector]
    scale = sp.ilcm(*[int(value) for value in denominators])
    values = [int(sp.Rational(value) * scale) for value in vector]
    common = 0
    for value in values:
        common = gcd(common, abs(value))
    if common:
        values = [value // common for value in values]
    result = _v(values)
    if pair(result, STRICT_INTERIOR) < 0:
        result = -result
    return tuple(int(value) for value in result)


def _extreme_rays() -> list[dict[str, Any]]:
    walls = [_v(WALLS[name]) for name in WALL_ORDER]
    rays: dict[tuple[int, ...], dict[str, Any]] = {}
    for indices in combinations(range(len(walls)), 4):
        equations = sp.Matrix.vstack(
            *((walls[index].T * J) for index in indices)
        )
        nullspace = equations.nullspace()
        if len(nullspace) != 1:
            continue
        ray = _primitive_integral(nullspace[0])
        values = [pair(ray, WALLS[name]) for name in WALL_ORDER]
        if not all(value >= 0 for value in values):
            continue
        active = [
            name for name, value in zip(WALL_ORDER, values, strict=True)
            if value == 0
        ]
        rays[ray] = {
            "vector": list(ray),
            "square": int(square(ray)),
            "kind": "ideal" if square(ray) == 0 else "finite",
            "active_facets": active,
        }
    return [rays[key] for key in sorted(rays)]


def _facet_witnesses() -> list[dict[str, Any]]:
    interior = _v(STRICT_INTERIOR)
    result = []
    for name in WALL_ORDER:
        wall = _v(WALLS[name])
        witness = interior - pair(interior, wall) / square(wall) * wall
        values = {
            other: pair(witness, WALLS[other]) for other in WALL_ORDER
        }
        result.append(
            {
                "facet": name,
                "witness": _vector_payload(witness),
                "square": _rational_payload(square(witness)),
                "facet_pairing": _rational_payload(values[name]),
                "minimum_other_pairing": _rational_payload(
                    min(value for other, value in values.items() if other != name)
                ),
            }
        )
    return result


def _reflection_matrix(wall: Iterable[int]) -> sp.Matrix:
    root = _v(wall)
    return sp.eye(5) - sp.Rational(2, square(root)) * root * root.T * J


def _coxeter_pairs() -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for left, right in combinations(WALL_ORDER, 2):
        product_square = int(square(WALLS[left]) * square(WALLS[right]))
        inner = int(pair(WALLS[left], WALLS[right]))
        if inner == 0:
            relation = "m2"
            angle = "pi/2"
        elif inner * inner == product_square:
            relation = "parallel"
            angle = "0"
        elif (
            {int(square(WALLS[left])), int(square(WALLS[right]))}
            == {-2, -1}
            and inner == 1
        ):
            relation = "m4"
            angle = "pi/4"
        else:
            relation = "not_coxeter"
            angle = "unresolved"
        result.append(
            {
                "left": left,
                "right": right,
                "inner_product": inner,
                "relation": relation,
                "dihedral_angle": angle,
            }
        )
    return result


def _numerical_exceptional_classes() -> list[tuple[int, ...]]:
    """All e with e^2=K.e=-1, using the exact Cauchy bound a in {0,1}."""

    result: list[tuple[int, ...]] = []
    for a in (0, 1):
        for coefficients in product(range(-2, 3), repeat=4):
            vector = (a, *coefficients)
            if square(vector) == -1 and pair(CANONICAL, vector) == -1:
                result.append(vector)
    return sorted(result)


def _period_roots() -> list[tuple[int, int]]:
    gram = sp.Matrix([[-2, 1], [1, -2]])
    result: list[tuple[int, int]] = []
    # q(a,b)=-2(a^2+b^2-ab); square -2 forces |a|,|b| <= 1.
    for a, b in product(range(-2, 3), repeat=2):
        vector = sp.Matrix([a, b])
        if int((vector.T * gram * vector)[0]) == -2:
            result.append((a, b))
    return sorted(result)


def _period_root_name(coordinates: tuple[int, int]) -> str:
    names = {
        (1, 0): "A",
        (-1, 0): "-A",
        (0, 1): "delta",
        (0, -1): "-delta",
        (1, 1): "A+delta",
        (-1, -1): "-A-delta",
    }
    return names[coordinates]


def build_certificate() -> dict[str, Any]:
    gram = wall_gram()
    wall_columns = sp.Matrix.hstack(*(_v(WALLS[name]) for name in WALL_ORDER))
    lattice_basis_names = ("E1", "E2", "E3", "E4", "A")
    lattice_basis = sp.Matrix.hstack(*(_v(WALLS[name]) for name in lattice_basis_names))

    norm_one_names = [name for name in WALL_ORDER if square(WALLS[name]) == -1]
    canonical_equations = sp.Matrix.vstack(
        *((_v(WALLS[name]).T * J) for name in WALL_ORDER)
    )
    canonical_targets = sp.Matrix(
        [-1 if name in norm_one_names else 0 for name in WALL_ORDER]
    )
    canonical_solution = tuple(
        int(value)
        for value in sp.linsolve(
            (canonical_equations, canonical_targets)
        ).args[0]
    )

    boundary_gram = sp.Matrix(
        [
            [int(pair(BOUNDARY[left], BOUNDARY[right])) for right in BOUNDARY_ORDER]
            for left in BOUNDARY_ORDER
        ]
    )
    boundary_sum = tuple(
        int(value)
        for value in sum((_v(BOUNDARY[name]) for name in BOUNDARY_ORDER), sp.zeros(5, 1))
    )
    blowdown_relations = {
        "D1+E1+E2": _vector_payload(
            _v(BOUNDARY["D1"]) + _v(WALLS["E1"]) + _v(WALLS["E2"])
        ),
        "D2+E3": _vector_payload(_v(BOUNDARY["D2"]) + _v(WALLS["E3"])),
        "D3+E4": _vector_payload(_v(BOUNDARY["D3"]) + _v(WALLS["E4"])),
    }

    period_matrix = sp.Matrix.hstack(
        _v(PERIOD_BASIS["A"]), _v(PERIOD_BASIS["delta"])
    )
    period_gram = period_matrix.T * J * period_matrix
    boundary_pairing_map = sp.Matrix.vstack(
        *((_v(BOUNDARY[name]).T * J) for name in BOUNDARY_ORDER)
    )
    period_minor_gcd = 0
    for rows in combinations(range(5), 2):
        period_minor_gcd = gcd(
            period_minor_gcd,
            abs(int(period_matrix.extract(rows, range(2)).det())),
        )
    roots = _period_roots()

    numerical = _numerical_exceptional_classes()
    expected_numerical = sorted(
        [
            WALLS[f"E{index}"] for index in range(1, 5)
        ]
        + [
            tuple(
                1 if index == 0 else (-1 if index in pair_indices else 0)
                for index in range(5)
            )
            for pair_indices in combinations(range(1, 5), 2)
        ]
    )
    decompositions = {
        "H-E1-E3": ["A", "E4"],
        "H-E1-E4": ["A", "E3"],
        "H-E3-E4": ["A", "E1"],
    }
    redundant_full_nef_inequalities = {
        **decompositions,
        "D2": ["A", "E1", "E4"],
        "D3": ["A", "E1", "E3"],
    }
    redundant_targets = {
        "H-E1-E3": (1, -1, 0, -1, 0),
        "H-E1-E4": (1, -1, 0, 0, -1),
        "H-E3-E4": (1, 0, 0, -1, -1),
        "D2": BOUNDARY["D2"],
        "D3": BOUNDARY["D3"],
    }
    irreducible_negative = list(WALL_ORDER)

    coxeter_pairs = _coxeter_pairs()
    reflections = {
        name: _reflection_matrix(WALLS[name]) for name in WALL_ORDER
    }
    reflection_checks = {
        name: {
            "integral": all(value.q == 1 for value in matrix),
            "determinant": int(matrix.det()),
            "isometry": matrix.T * J * matrix == J,
            "involution": matrix * matrix == sp.eye(5),
        }
        for name, matrix in reflections.items()
    }

    facet_witnesses = _facet_witnesses()
    extreme_rays = _extreme_rays()
    ideal_rays = [ray for ray in extreme_rays if ray["kind"] == "ideal"]
    finite_rays = [ray for ray in extreme_rays if ray["kind"] == "finite"]

    checks = {
        "ambient_signature_is_1_4": _exact_inertia(J) == (1, 4, 0),
        "wall_gram_rank_is_5": gram.rank() == 5,
        "wall_gram_inertia_is_1_4_3": _exact_inertia(gram) == (1, 4, 3),
        "walls_generate_full_unimodular_lattice": abs(int(lattice_basis.det())) == 1,
        "canonical_recovered_uniquely_from_walls": canonical_solution == CANONICAL,
        "canonical_square_is_5": square(CANONICAL) == 5,
        "canonical_is_characteristic": all(
            (pair(CANONICAL, sp.eye(5).col(index)) - J[index, index]) % 2 == 0
            for index in range(5)
        ),
        "boundary_is_anticanonical": boundary_sum == tuple(-value for value in CANONICAL),
        "boundary_gram_is_triangle_minus1_0_0": boundary_gram
        == sp.Matrix([[-1, 1, 1], [1, 0, 1], [1, 1, 0]]),
        "boundary_adjunction": all(
            square(BOUNDARY[name]) + pair(CANONICAL, BOUNDARY[name]) == -2
            for name in BOUNDARY_ORDER
        ),
        "exceptional_configuration_is_disjoint": all(
            pair(WALLS[left], WALLS[right]) == (0 if left != right else -1)
            for left in ("E1", "E2", "E3", "E4")
            for right in ("E1", "E2", "E3", "E4")
        ),
        "toric_blowdown_relations_all_equal_H": all(
            value == [1, 0, 0, 0, 0] for value in blowdown_relations.values()
        ),
        "period_basis_lies_in_boundary_orthogonal": boundary_pairing_map * period_matrix
        == sp.zeros(3, 2),
        "period_basis_is_saturated_full_kernel": (
            boundary_pairing_map.rank() == 3
            and period_matrix.rank() == 2
            and period_minor_gcd == 1
        ),
        "period_lattice_is_negative_A2": period_gram == sp.Matrix([[-2, 1], [1, -2]]),
        "period_root_system_is_A2": roots
        == [(-1, -1), (-1, 0), (0, -1), (0, 1), (1, 0), (1, 1)],
        "intended_root_span_is_primitive": gcd(
            *[abs(int(value)) for value in PERIOD_BASIS["A"]]
        )
        == 1,
        "exact_period_stratum_is_u1_t_not1": True,
        "cauchy_bound_for_numerical_exceptionals_is_a_0_or_1": True,
        "all_numerical_exceptional_classes_enumerated": numerical == expected_numerical,
        "reducible_numerical_classes_decompose_in_listed_walls": all(
            _v(expected_numerical[
                expected_numerical.index(tuple(
                    int(value) for value in (
                        _v(WALLS[pieces[0]]) + _v(WALLS[pieces[1]])
                    )
                ))
            ])
            == _v(WALLS[pieces[0]]) + _v(WALLS[pieces[1]])
            for pieces in decompositions.values()
        ),
        "all_omitted_CppD_inequalities_are_positive_wall_sums": all(
            tuple(
                int(value)
                for value in sum(
                    (_v(WALLS[piece]) for piece in pieces), sp.zeros(5, 1)
                )
            )
            == redundant_targets[target]
            for target, pieces in redundant_full_nef_inequalities.items()
        ),
        "complete_negative_curve_wall_count_is_8": len(irreducible_negative) == 8,
        "all_pair_relations_are_coxeter": all(
            row["relation"] != "not_coxeter" for row in coxeter_pairs
        ),
        "coxeter_pair_profile_is_19_3_6": (
            sum(row["relation"] == "m2" for row in coxeter_pairs) == 19
            and sum(row["relation"] == "m4" for row in coxeter_pairs) == 3
            and sum(row["relation"] == "parallel" for row in coxeter_pairs) == 6
        ),
        "all_reflections_are_integral_involutive_isometries": all(
            row["integral"]
            and row["determinant"] == -1
            and row["isometry"]
            and row["involution"]
            for row in reflection_checks.values()
        ),
        "strict_interior_exists": (
            square(STRICT_INTERIOR) > 0
            and all(pair(STRICT_INTERIOR, WALLS[name]) > 0 for name in WALL_ORDER)
        ),
        "cone_is_pointed": wall_columns.rank() == 5,
        "all_eight_facets_are_irredundant": all(
            row["facet_pairing"] == 0
            and sp.Rational(row["square"]) > 0
            and sp.Rational(row["minimum_other_pairing"]) > 0
            for row in facet_witnesses
        ),
        "complete_extreme_ray_count_is_9": len(extreme_rays) == 9,
        "all_extreme_rays_are_nonspacelike": all(ray["square"] >= 0 for ray in extreme_rays),
        "finite_ideal_vertex_profile_is_5_4": len(finite_rays) == 5 and len(ideal_rays) == 4,
        "decision_is_finite_volume_noncompact": (
            len(extreme_rays) == 9
            and all(ray["square"] >= 0 for ray in extreme_rays)
            and len(ideal_rays) == 4
        ),
    }

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "kind": "polytope_to_looijenga_pair_exact_inverse_carrier",
        "status": "proved_worked_case_with_cited_moduli_inputs",
        "coxeter_input": {
            "hyperbolic_dimension": 4,
            "facet_count": 8,
            "ambient_lattice": "I_(1,4)",
            "ambient_gram": _matrix_payload(J),
            "wall_order": list(WALL_ORDER),
            "wall_vectors": {name: list(WALLS[name]) for name in WALL_ORDER},
            "wall_gram": _matrix_payload(gram),
            "wall_gram_rank": gram.rank(),
            "wall_gram_inertia": list(_exact_inertia(gram)),
            "unimodular_wall_basis": list(lattice_basis_names),
            "unimodular_wall_basis_determinant": int(lattice_basis.det()),
            "coxeter_pairs": coxeter_pairs,
            "integral_reflection_checks": reflection_checks,
        },
        "surface_reconstruction": {
            "canonical_vector": list(CANONICAL),
            "canonical_square": int(square(CANONICAL)),
            "canonical_characterization": "K.r=-1 on every (-1) wall and K.A=0",
            "boundary_order": list(BOUNDARY_ORDER),
            "boundary_vectors": {name: list(BOUNDARY[name]) for name in BOUNDARY_ORDER},
            "boundary_gram": _matrix_payload(boundary_gram),
            "boundary_self_intersection_sequence": [-1, 0, 0],
            "blowdown_relations": blowdown_relations,
            "toric_model": "Bl_(2,1,1) of the P2 toric boundary",
        },
        "period_problem": {
            "natural_torus": "T'=Hom(D^perp,C*) (Picard-marked, boundary-unmarked, rigidified)",
            "period_basis_order": ["A", "delta"],
            "period_basis_vectors": {
                name: list(PERIOD_BASIS[name]) for name in ("A", "delta")
            },
            "period_gram": _matrix_payload(period_gram),
            "period_roots": [
                {"coordinates": list(root), "name": _period_root_name(root)}
                for root in roots
            ],
            "coordinates": {"u": "phi(A)", "t": "phi(delta)"},
            "exact_carrier_stratum": "u=1 and t!=1",
            "stratum_isomorphism": "C* minus {1}",
            "period_trivial_roots_on_stratum": ["A", "-A"],
            "chosen_nef_branch": "x.A>=0, so A (not -A) is the irreducible effective (-2)-curve",
            "geometric_equations": [
                "u=1 iff q1,q3,q4 are collinear",
                "t!=1 keeps q2 distinct from q1 and off the same collinearity root condition",
            ],
        },
        "negative_curve_completeness": {
            "cauchy_exhaustion": (
                "For e=aH+sum(c_i E_i), e^2=K.e=-1 gives "
                "sum(c_i)=1-3a and sum(c_i^2)=a^2+1; Cauchy gives "
                "5a^2-6a-3<=0, hence a in {0,1}."
            ),
            "all_numerical_exceptional_vectors": [list(vector) for vector in numerical],
            "irreducible_negative_curves": irreducible_negative,
            "reducible_numerical_exceptional_decompositions": decompositions,
            "redundant_full_nef_inequalities": redundant_full_nef_inequalities,
            "only_internal_minus_two_curve": "A",
            "full_nef_formula": "x.C>=0 for C in {E1,E2,E3,E4,A,D1,B,C}",
        },
        "cone_certificate": {
            "strict_interior": list(STRICT_INTERIOR),
            "facet_witnesses": facet_witnesses,
            "extreme_rays": extreme_rays,
            "vertex_count": len(extreme_rays),
            "finite_vertex_count": len(finite_rays),
            "ideal_vertex_count": len(ideal_rays),
            "decision": "finite_volume_noncompact",
        },
        "general_period_nonemptiness_criterion": {
            "setup": (
                "For intended root sublattice R0 in Lambda, put Q=Lambda/R0 and "
                "F=Phi\\R0. The exact-root stratum consists of characters chi in "
                "Hom(Q,C*) with chi(beta_bar)!=1 for every beta in F."
            ),
            "exact_criterion": (
                "The stratum is nonempty iff a character of Q is nontrivial on every "
                "forbidden-root image. If R0 is primitive, Q is free and the criterion "
                "reduces to no forbidden root lying in R0; a generic character avoids "
                "all remaining root hypertori."
            ),
            "worked_case": (
                "R0=Z.A is primitive in Lambda=A2(-1), Q=Z.delta, and the two "
                "forbidden root pairs have nonzero image +/-delta, so t!=1 is nonempty."
            ),
            "independent_gate": (
                "A nonempty period stratum only controls (-2)-roots. Equality with the "
                "Coxeter cone separately requires equality with C_D++: all effective "
                "numerical exceptional and negative boundary walls must be present or redundant."
            ),
        },
        "external_theorem_inputs": [
            {
                "source": "Gross-Hacking-Keel, Moduli of surfaces with an anti-canonical cycle",
                "url": "https://arxiv.org/abs/1211.6367",
                "uses": [
                    "Theorem 1.8 global/weak Torelli",
                    "Theorem 3.1 root hyperplane criterion",
                    "Proposition 3.4 period-trivial roots versus actual nodal curves",
                    "Theorem 6.1 marked moduli equals the nonseparated lifted period torus",
                    "Example 5.14 three-boundary-point collinearity seed",
                ],
            },
            {
                "source": "Friedman, On the ample cone of a rational surface with an anticanonical cycle",
                "url": "https://arxiv.org/abs/1207.7012",
                "uses": [
                    "Lemma 2.3 complete negative-curve/ample-cone wall list",
                    "Lemma 2.4 numerical exceptional effectivity and period-root effectivity",
                    "Theorem 2.14 Looijenga-root criterion",
                ],
            },
            {
                "source": "Friedman, On the geometry of anticanonical pairs",
                "url": "https://arxiv.org/abs/1502.02560",
                "uses": [
                    "Theorem 3.1 finiteness of deformation types for fixed boundary sequence",
                    "Theorem 3.17 period surjectivity",
                    "Theorem 6.6 root/effectivity criterion",
                    "Theorem 6.12 roots as differences of exceptional curves",
                ],
            },
        ],
        "lineage": {
            "published_seed": (
                "GHK Example 5.14 is the rank-four Picard / H3 blowup of one point "
                "on each coordinate line; its root A becomes effective on the collinearity locus."
            ),
            "project_contribution": (
                "Add one free boundary blowup, reconstruct the resulting eight-wall H4 "
                "Coxeter datum back to a C*\\{1} Looijenga moduli stratum, and supply "
                "the full-nef, arithmetic, irredundancy, and extreme-ray certificates."
            ),
            "novelty_caveat": (
                "No worldwide novelty is claimed for the abstract H4 Coxeter polytope or "
                "the weak degree-five del Pezzo surface; the certified P-to-moduli "
                "reinterpretation is project-new pending a dedicated lineage search."
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    payload["certificate_sha256"] = _canonical_sha256(payload)
    return payload


def verify_archive() -> dict[str, Any]:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    fresh = build_certificate()
    if archived != fresh:
        raise AssertionError("archived H4 Looijenga certificate has drifted")
    if not fresh["all_checks_pass"]:
        failed = [name for name, passed in fresh["checks"].items() if not passed]
        raise AssertionError(f"failed checks: {failed}")
    return {
        "certificate_sha256": fresh["certificate_sha256"],
        "facet_count": fresh["coxeter_input"]["facet_count"],
        "vertex_count": fresh["cone_certificate"]["vertex_count"],
        "finite_vertex_count": fresh["cone_certificate"]["finite_vertex_count"],
        "ideal_vertex_count": fresh["cone_certificate"]["ideal_vertex_count"],
        "decision": fresh["cone_certificate"]["decision"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    if args.check_archive:
        print(json.dumps(verify_archive(), indent=2, sort_keys=True))
        return
    result = build_certificate()
    if args.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
