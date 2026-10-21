#!/usr/bin/env python3
"""Certify the full-NS chamber of Kondo's Type-MI Coble surface.

The surface is the blow-up of ``P1 x P1`` in the ten intersection points of
Kondo's two characteristic-three boundary curves.  The thirty published
``(1,1)`` curves are reconstructed from the duad--syntheme--triad design, so
their divisor classes and every incidence are obtained without numerical
coordinates over ``F9``.  The script then constructs the saturated
Coble--Mukai lattice, enumerates all primitive cusps of its forty-root
chamber, and tests every integral boundary lift.

This is a surface-to-polytope certificate.  It does not claim an unmarked
polytope-to-surface existence theorem or worldwide priority for the literal
polytope.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
from typing import Any, Iterable, Sequence

import sympy as sp
from sympy.matrices.normalforms import (
    hermite_normal_form,
    smith_normal_decomp,
    smith_normal_form,
)
from sympy.polys.domains import ZZ

from scripts.certify_coble_one_boundary_full_ns_h10 import _cdd_rays
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/coble_type_mi_full_ns_h11.json"
SCHEMA = "coxeter.coble_type_mi_full_ns_h11.v1"
FERMAT_LINEAGE = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_lineage.json"
)
FERMAT_VINBERG = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_vinberg_certificate.json"
)
FERMAT_FIXED = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_survivor_fixed_lattices.json"
)
FERMAT_COXETER = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_coxeter_screen.json"
)
FERMAT_EXACT_CONE = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_exact_cone_certificate.json"
)
FERMAT_EXACT_RAYS = (
    ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_exact_cone_rays.jsonl.gz"
)

PINNED = {
    "kondo_source_tar": "6989648564a80167f2daa35da64b90a063a0a090b9f23ff843b7fabb9d96a98d",
    "kondo_tex": "8e23df4e77a2666a7a4f9a634c1ec2f36066353bf0bdcd0a72b11a0a70bfb137",
    "fermat_lineage": "e1b20c853d3f2f254900c0dd46cdf187093292fb4356e1d6e9de23ca839faf89",
    "fermat_vinberg": "22f2104680afc9f20f9a2a6fd38f7486037f9890b57cd881c1fda6c96027c9bb",
    "fermat_fixed": "431ddb789fcb30e843efe870a67f9ce276ccd75c4fb65ea4f01f51afdc836124",
    "fermat_coxeter": "6bfeac5a6f7530d641484d920091e11539565ba608353d9eb221eb76010b9921",
    "fermat_exact_cone": "c3be5effcdc72cb9c5c50a2cdfc4280083f113b42bbb0d24508b921fbcbc23c9",
    "fermat_exact_rays": "07f7b280aa339bcc4eaadc7b885b1492434b79a01e7464957f4db55546eb06bc",
}

# ``TYPE_MI_TO_FERMAT`` sends our order
#
#   B1,B2; 30 (1,1) curves; 10 exceptional curves; 164 cusp lifts
#
# to the oriented mirror order in the independently constructed
# ``F3-C2-H11-206`` Fermat-quartic chamber.  The equality of all 42,436 Gram
# entries is checked below; facet count or row fingerprints are not used as a
# substitute for a marked isometry.
TYPE_MI_TO_FERMAT = (
    60, 26, 3, 4, 55, 9, 16, 48, 21, 38, 28, 18, 40, 29, 5, 17, 46,
    32, 35, 37, 30, 36, 39, 10, 2, 1, 34, 49, 50, 43, 45, 47, 44, 11,
    53, 58, 13, 52, 57, 24, 15, 61, 92, 94, 75, 104, 103, 136, 160,
    164, 82, 123, 120, 132, 102, 137, 158, 167, 71, 115, 125, 128, 161,
    165, 159, 169, 152, 157, 117, 131, 144, 168, 111, 114, 86, 100, 72,
    78, 64, 108, 84, 126, 90, 95, 80, 122, 150, 156, 116, 127, 88, 101,
    98, 134, 77, 105, 70, 113, 142, 166, 110, 112, 73, 79, 96, 135, 118,
    133, 124, 129, 66, 109, 85, 130, 162, 196, 151, 193, 93, 188, 99,
    185, 121, 182, 146, 204, 89, 187, 83, 175, 153, 195, 149, 201, 119,
    183, 140, 203, 76, 181, 69, 173, 145, 190, 139, 199, 163, 197, 155,
    194, 147, 205, 141, 202, 87, 184, 81, 171, 74, 180, 65, 172, 68, 170,
    63, 178, 91, 189, 154, 192, 67, 174, 107, 176, 97, 186, 143, 191,
    62, 179, 106, 177, 148, 198, 138, 200, 0, 22, 41, 59, 19, 25, 7,
    12, 42, 56, 20, 23, 8, 14, 33, 54, 31, 51, 6, 27,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("an integral matrix became nonintegral")
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _rational_rows(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(value) if value.q == 1 else str(value) for value in matrix.row(row)]
        for row in range(matrix.rows)
    ]


def _half_integral_basis(
    generators: Sequence[sp.Matrix], rank: int = 12
) -> sp.Matrix:
    raw = sp.Matrix.hstack(*(2 * vector for vector in generators))
    doubled = sp.Matrix(
        [[int(raw[row, column]) for column in range(raw.cols)] for row in range(raw.rows)]
    )
    if doubled.rows != rank or doubled.rank() != rank:
        raise ArithmeticError("the half-integral generators do not have full rank")
    if any(value.q != 1 for value in doubled):
        raise ArithmeticError("a proposed generator has denominator greater than two")
    return hermite_normal_form(doubled) / 2


def _perfect_matchings(values: tuple[int, ...]) -> list[tuple[tuple[int, int], ...]]:
    """Return the fifteen synthemes on six letters."""

    if not values:
        return [tuple()]
    first = values[0]
    result: list[tuple[tuple[int, int], ...]] = []
    for position in range(1, len(values)):
        second = values[position]
        remainder = values[1:position] + values[position + 1 :]
        for tail in _perfect_matchings(remainder):
            result.append(tuple(sorted(((first, second), *tail))))
    return sorted(set(result))


def _surface_model() -> dict[str, Any]:
    """Build Bl_10(P1 x P1) and its two boundaries and forty root curves."""

    gram = sp.zeros(12)
    gram[0, 1] = gram[1, 0] = 1
    for index in range(2, 12):
        gram[index, index] = -1
    identity = sp.eye(12)
    f1, f2 = identity[:, 0], identity[:, 1]
    exceptional = [identity[:, 2 + index] for index in range(10)]
    total_exceptional = sum(exceptional, sp.zeros(12, 1))
    boundaries = [
        f1 + 3 * f2 - total_exceptional,
        3 * f1 + f2 - total_exceptional,
    ]
    canonical = -2 * f1 - 2 * f2 + total_exceptional

    letters = tuple(range(6))
    duads = list(itertools.combinations(letters, 2))
    synthemes = _perfect_matchings(letters)
    triads = [frozenset((0, *pair)) for pair in itertools.combinations(range(1, 6), 2)]
    if not (len(duads) == len(synthemes) == 15 and len(triads) == 10):
        raise ArithmeticError("the duad--syntheme--triad design changed")

    def triad_contains_duad(triad: frozenset[int], duad: tuple[int, int]) -> bool:
        complement = set(letters) - set(triad)
        return set(duad).issubset(triad) or set(duad).issubset(complement)

    duad_blocks = [
        tuple(index for index, triad in enumerate(triads) if triad_contains_duad(triad, duad))
        for duad in duads
    ]
    syntheme_blocks = [
        tuple(
            index
            for index, triad in enumerate(triads)
            if all(len(set(pair) & set(triad)) == 1 for pair in syntheme)
        )
        for syntheme in synthemes
    ]
    blocks = [*duad_blocks, *syntheme_blocks]
    if any(len(block) != 4 for block in blocks) or len(set(blocks)) != 30:
        raise ArithmeticError("the thirty four-point blocks changed")

    root_curves = [
        f1 + f2 - sum((exceptional[index] for index in block), sp.zeros(12, 1))
        for block in blocks
    ]
    root_labels = [
        *(f"D_{left + 1}{right + 1}" for left, right in duads),
        *(
            "S_" + "_".join(f"{left + 1}{right + 1}" for left, right in syntheme)
            for syntheme in synthemes
        ),
    ]
    exceptional_labels = [
        "E_T" + "".join(str(value + 1) for value in sorted(triad))
        for triad in triads
    ]

    boundary_matrix = sp.Matrix.hstack(*boundaries)
    curve_matrix = sp.Matrix.hstack(*root_curves)
    design_intersections = curve_matrix.T * gram * curve_matrix
    boundary_curve_intersections = boundary_matrix.T * gram * curve_matrix
    exceptional_curve_intersections = sp.Matrix.hstack(*exceptional).T * gram * curve_matrix
    checks = {
        "NS_is_odd_unimodular_I_1_11": gram.det() == -1
        and exact_inertia(gram) == (1, 11, 0)
        and any(int(gram[index, index]) % 2 for index in range(12)),
        "anti_bicanonical_boundary": sum(boundaries, sp.zeros(12, 1)) == -2 * canonical,
        "two_disjoint_minus_four_boundaries": boundary_matrix.T * gram * boundary_matrix
        == -4 * sp.eye(2),
        "canonical_square_minus_two": (canonical.T * gram * canonical)[0] == -2,
        "thirty_minus_two_curves": Counter(design_intersections.diagonal())
        == Counter({-2: 30}),
        "root_curves_disjoint_from_boundary": boundary_curve_intersections == sp.zeros(2, 30),
        "four_exceptional_incidences_per_root_curve": all(
            sum(exceptional_curve_intersections[:, column]) == 4
            and all(value in (0, 1) for value in exceptional_curve_intersections[:, column])
            for column in range(30)
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"the Type-MI surface model failed: {checks}")

    pic_tilde_basis = _half_integral_basis(
        [*identity.columnspace(), *(boundary / 2 for boundary in boundaries)]
    )
    return {
        "gram": gram,
        "f1": f1,
        "f2": f2,
        "exceptional": exceptional,
        "boundaries": boundaries,
        "canonical": canonical,
        "root_curves": root_curves,
        "root_labels": root_labels,
        "exceptional_labels": exceptional_labels,
        "blocks": blocks,
        "duads": duads,
        "synthemes": synthemes,
        "triads": triads,
        "pic_tilde_basis": pic_tilde_basis,
        "checks": checks,
    }


def _cm_chamber(surface: dict[str, Any], executable: str) -> dict[str, Any]:
    """Construct saturated CM(S), enumerate its rays, and test cusp lifts."""

    gram: sp.Matrix = surface["gram"]
    boundaries: list[sp.Matrix] = surface["boundaries"]
    exceptional: list[sp.Matrix] = surface["exceptional"]
    root_curves: list[sp.Matrix] = surface["root_curves"]
    pic_basis: sp.Matrix = surface["pic_tilde_basis"]

    constraint = 2 * sp.Matrix.hstack(*boundaries).T * gram * pic_basis
    if any(value.q != 1 for value in constraint):
        raise ArithmeticError("the CM orthogonality constraint became nonintegral")
    diagonal, left_transform, right_transform = smith_normal_decomp(
        constraint, domain=ZZ
    )
    if left_transform * constraint * right_transform != diagonal:
        raise ArithmeticError("the CM Smith decomposition failed")
    rank = constraint.rank()
    cm_basis = pic_basis * right_transform[:, rank:]
    cm_gram = cm_basis.T * gram * cm_basis
    if cm_basis.shape != (12, 10) or exact_inertia(cm_gram) != (1, 9, 0):
        raise ArithmeticError("the saturated Coble--Mukai lattice changed")

    folded_exceptional = [
        boundaries[0] / 2 + boundaries[1] / 2 + 2 * curve
        for curve in exceptional
    ]
    cm_roots = [*root_curves, *folded_exceptional]
    root_labels = [*surface["root_labels"], *surface["exceptional_labels"]]
    root_matrix = sp.Matrix.hstack(*cm_roots)
    root_coordinates = cm_basis.gauss_jordan_solve(root_matrix)[0]
    if any(value.q != 1 for value in root_coordinates):
        raise ArithmeticError("a published Type-MI root left CM integrally")
    root_gram = root_matrix.T * gram * root_matrix
    root_span_smith = [
        abs(int(value))
        for value in smith_normal_form(root_coordinates, domain=ZZ).diagonal()
    ]

    inequalities = root_coordinates.T * cm_gram
    rays = sorted(_cdd_rays(inequalities, executable))
    counts: Counter[str] = Counter()
    cusps: list[tuple[int, ...]] = []
    for raw in rays:
        vector = sp.Matrix(raw)
        square = int((vector.T * cm_gram * vector)[0])
        kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        counts[kind] += 1
        if square == 0:
            cusps.append(tuple(int(value) for value in raw))

    cusp_records = []
    lifted_curves: list[sp.Matrix] = []
    lifted_labels: list[str] = []
    for cusp_index, cusp in enumerate(cusps):
        fiber = cm_basis * sp.Matrix(cusp)
        active = [
            index
            for index, value in enumerate(inequalities * sp.Matrix(cusp))
            if value == 0
        ]
        membership = []
        for boundary_index, boundary in enumerate(boundaries):
            candidate = fiber - boundary / 2
            integral = all(value.q == 1 for value in candidate)
            membership.append(integral)
            if integral:
                lifted_curves.append(candidate)
                lifted_labels.append(f"E_cusp_{cusp_index + 1}_B{boundary_index + 1}")
        cusp_records.append(
            {
                "CM_coordinates": list(cusp),
                "active_root_indices_zero_based": active,
                "active_root_labels": [root_labels[index] for index in active],
                "single_boundary_integral_lift_membership": membership,
            }
        )

    return {
        "basis": cm_basis,
        "gram": cm_gram,
        "root_coordinates": root_coordinates,
        "root_gram": root_gram,
        "root_labels": root_labels,
        "root_span_smith": root_span_smith,
        "inequalities": inequalities,
        "rays": rays,
        "counts": counts,
        "cusps": cusps,
        "cusp_records": cusp_records,
        "lifted_curves": lifted_curves,
        "lifted_labels": lifted_labels,
    }


def _complete_walls(surface: dict[str, Any], cm: dict[str, Any]) -> dict[str, Any]:
    gram: sp.Matrix = surface["gram"]
    walls: list[sp.Matrix] = [
        *surface["boundaries"],
        *surface["root_curves"],
        *surface["exceptional"],
        *cm["lifted_curves"],
    ]
    labels = [
        "B1",
        "B2",
        *surface["root_labels"],
        *surface["exceptional_labels"],
        *cm["lifted_labels"],
    ]
    matrix = sp.Matrix.hstack(*walls)
    negative_gram = matrix.T * gram * matrix
    basis_positions = list(matrix.rref()[1])
    if len(basis_positions) != 12:
        raise ArithmeticError("the complete Type-MI wall stock does not span NS")
    basis_gram = negative_gram.extract(basis_positions, basis_positions)

    relation_histogram: Counter[str] = Counter()
    coxeter = True
    for left in range(len(walls)):
        for right in range(left):
            pairing = negative_gram[left, right]
            if pairing < 0:
                coxeter = False
                continue
            left_norm = -negative_gram[left, left]
            right_norm = -negative_gram[right, right]
            ratio = sp.cancel(4 * pairing**2 / (left_norm * right_norm))
            if ratio == 0:
                relation_histogram["m2"] += 1
            elif ratio == 1:
                relation_histogram["m3"] += 1
            elif ratio == 2:
                relation_histogram["m4"] += 1
            elif ratio == 3:
                relation_histogram["m6"] += 1
            elif ratio == 4:
                relation_histogram["parallel"] += 1
            elif ratio > 4:
                relation_histogram["ultraparallel"] += 1
            else:
                coxeter = False

    # A single exact ample class gives a timelike relative-interior witness
    # for every facet by orthogonal projection.  This avoids mistaking one of
    # the 164 cusp-lift inequalities for a redundant nef inequality.
    strict = (
        3 * surface["f1"]
        + 3 * surface["f2"]
        - sum(surface["exceptional"], sp.zeros(12, 1))
    )
    strict_square = (strict.T * gram * strict)[0]
    strict_pairings = [(strict.T * gram * wall)[0] for wall in walls]
    witness_squares = []
    witness_margins = []
    for index, wall in enumerate(walls):
        square = (wall.T * gram * wall)[0]
        witness = strict - strict_pairings[index] / square * wall
        pairings = [(witness.T * gram * other)[0] for other in walls]
        if not (
            pairings[index] == 0
            and all(value > 0 for other, value in enumerate(pairings) if other != index)
            and (witness.T * gram * witness)[0] > 0
        ):
            raise ArithmeticError(f"wall {index} failed its exact facet witness")
        witness_squares.append((witness.T * gram * witness)[0])
        witness_margins.append(
            min(value for other, value in enumerate(pairings) if other != index)
        )

    physical_rows = _rows(matrix.T)
    return {
        "labels": labels,
        "walls": walls,
        "wall_matrix": matrix,
        "negative_gram": negative_gram,
        "basis_positions": basis_positions,
        "basis_gram": basis_gram,
        "physical_NS_wall_vectors": physical_rows,
        "negative_wall_gram_sha256": hashlib.sha256(
            json.dumps(_rows(negative_gram), separators=(",", ":")).encode()
        ).hexdigest(),
        "rank": matrix.rank(),
        "gram_inertia": exact_inertia(basis_gram),
        "full_gram_nullity": len(walls) - matrix.rank(),
        "norm_histogram": Counter(negative_gram.diagonal()),
        "relation_histogram": relation_histogram,
        "all_pairs_coxeter": coxeter,
        "strict_interior": strict,
        "strict_interior_square": strict_square,
        "strict_interior_minimum_wall_pairing": min(strict_pairings),
        "all_facets_irredundant": True,
        "facet_witness_minimum_square": min(witness_squares),
        "facet_witness_minimum_other_wall_margin": min(witness_margins),
        "all_walls_distinct": len({tuple(row) for row in physical_rows}) == len(walls),
    }


def _arithmetic_lattice(surface: dict[str, Any], full: dict[str, Any]) -> dict[str, Any]:
    gram: sp.Matrix = surface["gram"]
    pic_basis: sp.Matrix = surface["pic_tilde_basis"]
    arithmetic_gram = 2 * pic_basis.T * gram * pic_basis
    arithmetic_normals = [
        surface["boundaries"][0] / 2,
        surface["boundaries"][1] / 2,
        *full["walls"][2:],
    ]
    root_coordinates = pic_basis.inv() * sp.Matrix.hstack(*arithmetic_normals)
    if any(value.q != 1 for value in root_coordinates):
        raise ArithmeticError("an arithmetic Type-MI mirror left doubled Pic~")
    if any(
        math.gcd(*(abs(int(value)) for value in root_coordinates[:, index])) != 1
        for index in range(root_coordinates.cols)
    ):
        raise ArithmeticError("an arithmetic Type-MI mirror is imprimitive")

    smith_matrix = smith_normal_form(root_coordinates, domain=ZZ)
    root_span_smith = [
        abs(int(smith_matrix[index, index])) for index in range(12)
    ]
    lattice_smith_matrix = smith_normal_form(arithmetic_gram, domain=ZZ)
    lattice_smith = [
        abs(int(lattice_smith_matrix[index, index])) for index in range(12)
    ]
    root_types: Counter[tuple[int, int]] = Counter()
    integral_reflections = []
    for index in range(root_coordinates.cols):
        root = root_coordinates[:, index]
        square = int((root.T * arithmetic_gram * root)[0])
        divisibility = math.gcd(
            *(abs(int(value)) for value in arithmetic_gram * root)
        )
        root_types[(square, divisibility)] += 1
        reflection = (
            sp.eye(12)
            - sp.Rational(2, square) * root * (root.T * arithmetic_gram)
        )
        integral_reflections.append(
            all(value.q == 1 for value in reflection)
            and reflection.T * arithmetic_gram * reflection == arithmetic_gram
            and reflection**2 == sp.eye(12)
        )

    ns_reflection_integrality = []
    for wall in full["walls"]:
        square = (wall.T * gram * wall)[0]
        reflection = sp.eye(12) - sp.Rational(2, square) * wall * (wall.T * gram)
        ns_reflection_integrality.append(all(value.q == 1 for value in reflection))

    dual = arithmetic_gram.inv()
    delta = int(any(value.q != 1 for value in dual.diagonal()))
    return {
        "basis": pic_basis,
        "gram": arithmetic_gram,
        "root_coordinates": root_coordinates,
        "determinant": int(arithmetic_gram.det()),
        "smith": lattice_smith,
        "inertia": exact_inertia(arithmetic_gram),
        "even": all(int(arithmetic_gram[index, index]) % 2 == 0 for index in range(12)),
        "nikulin_delta": delta,
        "root_span_smith": root_span_smith,
        "root_types": root_types,
        "all_reflections_integral": all(integral_reflections),
        "NS_nonintegral_reflection_positions": [
            index for index, valid in enumerate(ns_reflection_integrality) if not valid
        ],
    }


def _orient_existing_roots(root_gram: sp.Matrix) -> tuple[list[int], sp.Matrix]:
    signs: list[int | None] = [None] * root_gram.rows
    signs[0] = 1
    queue: deque[int] = deque([0])
    while queue:
        left = queue.popleft()
        for right in range(root_gram.rows):
            if left == right or root_gram[left, right] == 0:
                continue
            proposed = signs[left] * (1 if root_gram[left, right] > 0 else -1)
            if signs[right] is None:
                signs[right] = proposed
                queue.append(right)
            elif signs[right] != proposed:
                raise ArithmeticError("the existing Fermat mirror Gram is not acutely orientable")
    if any(sign is None for sign in signs):
        raise ArithmeticError("the existing Fermat mirror Gram is disconnected")
    integer_signs = [int(sign) for sign in signs]
    diagonal = sp.diag(*integer_signs)
    oriented = diagonal * root_gram * diagonal
    if any(oriented[left, right] < 0 for left in range(oriented.rows) for right in range(left)):
        raise ArithmeticError("the Fermat orientation replay failed")
    return integer_signs, oriented


def _fermat_match(full: dict[str, Any], arithmetic: dict[str, Any]) -> dict[str, Any]:
    upstream = {
        FERMAT_LINEAGE: PINNED["fermat_lineage"],
        FERMAT_VINBERG: PINNED["fermat_vinberg"],
        FERMAT_FIXED: PINNED["fermat_fixed"],
        FERMAT_COXETER: PINNED["fermat_coxeter"],
        FERMAT_EXACT_CONE: PINNED["fermat_exact_cone"],
        FERMAT_EXACT_RAYS: PINNED["fermat_exact_rays"],
    }
    for path, expected in upstream.items():
        if _sha256(path) != expected:
            raise ArithmeticError(f"a pinned Fermat input changed: {path}")
    lineage = json.loads(FERMAT_LINEAGE.read_text(encoding="utf-8"))
    vinberg = json.loads(FERMAT_VINBERG.read_text(encoding="utf-8"))
    fixed = json.loads(FERMAT_FIXED.read_text(encoding="utf-8"))
    coxeter = json.loads(FERMAT_COXETER.read_text(encoding="utf-8"))
    exact = json.loads(FERMAT_EXACT_CONE.read_text(encoding="utf-8"))
    fixed_case = next(case for case in fixed["cases"] if case["id"] == "char3-C2-H11")
    fixed_gram = sp.Matrix(fixed_case["gram"])
    mirror_records = lineage["induced_206_wall_reflection_subgroup"]["mirrors"]
    existing_roots = sp.Matrix.hstack(
        *(sp.Matrix(record["fixed_basis_coordinates"]) for record in mirror_records)
    )
    existing_raw_gram = existing_roots.T * fixed_gram * existing_roots
    orientation_signs, existing_oriented_gram = _orient_existing_roots(existing_raw_gram)
    existing_oriented_roots = existing_roots * sp.diag(*orientation_signs)

    if sorted(TYPE_MI_TO_FERMAT) != list(range(206)):
        raise ArithmeticError("the Type-MI/Fermat facet map is not a permutation")
    matched_roots = existing_oriented_roots[:, list(TYPE_MI_TO_FERMAT)]
    type_mi_gram = arithmetic["root_coordinates"].T * arithmetic["gram"] * arithmetic["root_coordinates"]
    if type_mi_gram != matched_roots.T * fixed_gram * matched_roots:
        raise ArithmeticError("the complete Type-MI and Fermat 206-wall Grams differ")

    basis_positions = list(arithmetic["root_coordinates"].rref()[1])
    source_basis = arithmetic["root_coordinates"][:, basis_positions]
    target_basis = matched_roots[:, basis_positions]
    isometry = target_basis * source_basis.inv()
    if not (
        all(value.q == 1 for value in isometry)
        and abs(int(isometry.det())) == 1
        and isometry * arithmetic["root_coordinates"] == matched_roots
        and isometry.T * fixed_gram * isometry == arithmetic["gram"]
    ):
        raise ArithmeticError("the Type-MI/Fermat marked lattice isometry failed")

    if not (
        coxeter["status"]["all_facet_pairs_are_coxeter"]
        and vinberg["finite_volume"]
        and not vinberg["compact"]
        and vinberg["bad_edge_count"] == 0
        and vinberg["invalid_full_zero_neighborhood_count"] == 0
        and exact["finite_volume_certified"]
        and exact["finite_ray_count"] == vinberg["finite_vertex_count"]
        and exact["ideal_ray_count"] == vinberg["ideal_vertex_count"]
        and exact["hyperideal_ray_count"] == 0
    ):
        raise ArithmeticError("the transported Fermat finite-volume certificate failed")

    return {
        "facet_permutation_type_MI_to_Fermat_zero_based": list(TYPE_MI_TO_FERMAT),
        "facet_permutation_sha256": hashlib.sha256(
            json.dumps(TYPE_MI_TO_FERMAT, separators=(",", ":")).encode()
        ).hexdigest(),
        "existing_orientation_signs": orientation_signs,
        "marked_isometry_from_doubled_Pic_tilde_to_Fermat_fixed_basis": _rows(isometry),
        "marked_isometry_determinant": int(isometry.det()),
        "all_42436_Gram_entries_match": True,
        "same_oriented_polytope": True,
        "exact_cone": {
            "finite_vertex_count": exact["finite_ray_count"],
            "ideal_vertex_count": exact["ideal_ray_count"],
            "hyperideal_vertex_count": exact["hyperideal_ray_count"],
            "extreme_ray_count": exact["ray_count"],
            "active_facet_count_histogram": exact["active_facet_count_histogram"],
            "ray_archive_sha256": _sha256(FERMAT_EXACT_RAYS),
            "certificate_sha256": _sha256(FERMAT_EXACT_CONE),
        },
        "literal_vinberg": {
            "elliptic_edge_rank": vinberg["elliptic_edge_rank"],
            "elliptic_edge_count": vinberg["elliptic_edge_count"],
            "finite_vertex_count": vinberg["finite_vertex_count"],
            "ideal_vertex_count": vinberg["ideal_vertex_count"],
            "endpoint_count_histogram": vinberg["endpoint_count_histogram"],
            "endpoint_type_histogram": vinberg["endpoint_type_histogram"],
            "connected_elliptic_component_count": vinberg["connected_elliptic_component_count"],
            "connected_affine_component_count": vinberg["connected_affine_component_count"],
            "bad_edge_count": vinberg["bad_edge_count"],
            "unrecognized_endpoint_edge_count": vinberg["unrecognized_endpoint_edge_count"],
            "invalid_full_zero_neighborhood_count": vinberg["invalid_full_zero_neighborhood_count"],
            "criterion": vinberg["method"],
        },
        "known_lineage": {
            "existing_polytope_id": lineage["candidate"]["id"],
            "ambient_lattice": lineage["fixed_lattice_identification"]["target_lattice"],
            "nikulin_invariants": lineage["fixed_lattice_identification"]["nikulin_invariants"],
            "known_maximal_chamber": "Alexeev (12,10,1), 14 walls",
            "literal_worldwide_priority": lineage["lineage_decision"]["subgroup_chamber_priority"],
        },
        "upstream_sha256": {
            str(path.relative_to(ROOT)): digest for path, digest in upstream.items()
        },
    }


def build_certificate(cdd_executable: str | None = None) -> dict[str, Any]:
    executable = cdd_executable or shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for the independent CM cusp census")
    surface = _surface_model()
    cm = _cm_chamber(surface, executable)
    full = _complete_walls(surface, cm)
    arithmetic = _arithmetic_lattice(surface, full)
    match = _fermat_match(full, arithmetic)

    membership_histogram = Counter(
        tuple(record["single_boundary_integral_lift_membership"])
        for record in cm["cusp_records"]
    )
    checks = {
        **surface["checks"],
        "CM_is_E10": int(cm["gram"].det()) == -1
        and exact_inertia(cm["gram"]) == (1, 9, 0),
        "forty_complete_CM_roots_span_index_two": cm["root_span_smith"]
        == [1] * 9 + [2],
        "CM_exact_cone_6210_247_0": cm["counts"]
        == Counter({"finite": 6210, "ideal": 247}),
        "exactly_82_boundary_split_CM_cusps": membership_histogram
        == Counter({(False, False): 165, (True, True): 82}),
        "complete_negative_curve_split_2_30_174": full["norm_histogram"]
        == Counter({-1: 174, -2: 30, -4: 2}),
        "full_wall_rank_signature_non_degenerate": full["rank"] == 12
        and full["gram_inertia"] == (1, 11, 0)
        and full["full_gram_nullity"] == 194,
        "all_206_facets_irredundant": full["all_facets_irredundant"],
        "all_pairwise_angles_are_Coxeter": full["all_pairs_coxeter"],
        "doubled_extended_Picard_is_12_10_1": arithmetic["determinant"] == -1024
        and arithmetic["smith"] == [1, 1, *([2] * 10)]
        and arithmetic["inertia"] == (1, 11, 0)
        and arithmetic["even"]
        and arithmetic["nikulin_delta"] == 1,
        "all_mirrors_span_and_reflect_integrally": arithmetic["root_span_smith"]
        == [1] * 12
        and arithmetic["all_reflections_integral"],
        "arithmetic_mirror_types_match": arithmetic["root_types"]
        == Counter({(-2, 2): 164, (-4, 2): 30, (-2, 1): 12}),
        "only_boundaries_fail_to_reflect_NS_integrally": arithmetic[
            "NS_nonintegral_reflection_positions"
        ]
        == [0, 1],
        "marked_same_polytope_as_F3_C2_H11_206": match["same_oriented_polytope"],
        "exact_full_cone_327566_18354_0": match["exact_cone"][
            "finite_vertex_count"
        ]
        == 327566
        and match["exact_cone"]["ideal_vertex_count"] == 18354
        and match["exact_cone"]["hyperideal_vertex_count"] == 0,
        "literal_Vinberg_zero_failures": match["literal_vinberg"]["bad_edge_count"]
        == 0
        and match["literal_vinberg"]["unrecognized_endpoint_edge_count"] == 0
        and match["literal_vinberg"]["invalid_full_zero_neighborhood_count"] == 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Type-MI full-NS H11 checks failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "direction": "SURFACE_TO_POLYTOPE_with_exact_match_to_known_project_polytope",
        "claim": {
            "theorem": (
                "The projectivized full Neron--Severi nef cone of Kondo's unique "
                "characteristic-three Type-MI Coble surface is a 206-facet "
                "finite-volume noncompact arithmetic Coxeter polytope in H11, "
                "and is marked-isometric to the independently constructed "
                "F3-C2-H11-206 Fermat-quartic chamber."
            ),
            "novelty": (
                "project-new complete Type-MI full-nef realization and exact "
                "surface/K3 lineage bridge for an existing project polytope"
            ),
            "not_claimed": [
                "new abstract hyperbolic Coxeter polytope",
                "new arithmetic lattice, maximal reflection group, or commensurability class",
                "worldwide priority for the Type-MI full-nef interpretation",
                "a broad unmarked polytope-to-surface existence theorem",
            ],
        },
        "sources": {
            "kondo": {
                "title": "Coble surfaces with finite automorphism group",
                "arxiv": "2102.00591v2",
                "source_tar_sha256": PINNED["kondo_source_tar"],
                "tex_sha256": PINNED["kondo_tex"],
                "used": [
                    "Propositions 3.1, 3.4, 3.7, and 3.8",
                    "Lemma 3.11 and the primitive-isotropic/fibration correspondence",
                    "Example 4.17, Lemma 4.18, and Proposition 4.19",
                    "Section 6 Type MI and Section 7 uniqueness",
                ],
            },
            "surface_model": (
                "Blow up the ten intersections of the (1,3) and (3,1) boundary "
                "curves on P1 x P1; reconstruct the thirty (1,1) curves from "
                "the duad--syntheme--triad incidence design."
            ),
        },
        "surface": {
            "classification": "Kondo Type MI, characteristic 3, two boundaries",
            "model": "Bl_10(P1 x P1)",
            "NS_identification": "U + <-1>^10 = I_(1,11)",
            "NS_gram": _rows(surface["gram"]),
            "canonical_class": _rows(surface["canonical"].T)[0],
            "boundary_vectors": _rows(sp.Matrix.hstack(*surface["boundaries"]).T),
            "duads": [[left + 1, right + 1] for left, right in surface["duads"]],
            "synthemes": [
                [[left + 1, right + 1] for left, right in syntheme]
                for syntheme in surface["synthemes"]
            ],
            "triads": [[value + 1 for value in sorted(triad)] for triad in surface["triads"]],
            "four_point_blocks_zero_based": [list(block) for block in surface["blocks"]],
        },
        "coble_mukai": {
            "basis_in_NS_coordinates": _rational_rows(cm["basis"]),
            "gram": _rows(cm["gram"]),
            "determinant": int(cm["gram"].det()),
            "root_count": len(cm["root_labels"]),
            "root_labels": cm["root_labels"],
            "root_coordinates": _rows(cm["root_coordinates"].T),
            "root_span_smith_diagonal": cm["root_span_smith"],
            "finite_vertex_count": cm["counts"]["finite"],
            "ideal_vertex_count": cm["counts"]["ideal"],
            "hyperideal_vertex_count": cm["counts"]["hyperideal"],
            "primitive_cusp_count": len(cm["cusps"]),
            "cusp_lift_membership_histogram": {
                "none": membership_histogram[(False, False)],
                "both_boundaries": membership_histogram[(True, True)],
            },
            "cusp_records": cm["cusp_records"],
        },
        "negative_curve_exhaustion": {
            "complete_count": len(full["walls"]),
            "split": {
                "boundary_minus4": 2,
                "minus2_curves": 30,
                "minus1_meeting_both_boundaries_once": 10,
                "minus1_meeting_one_boundary_twice": 164,
            },
            "argument": [
                "Kondo proves that the displayed 40 classes are all effective irreducible CM roots: thirty are actual (-2)-curves and ten fold the ten exceptional (-1)-curves meeting both boundaries.",
                "Every remaining negative curve is a (-1)-curve meeting one boundary twice; E maps to the primitive CM cusp f=E+B_i/2.",
                "Conversely Kondo's primitive-isotropic/fibration and Halphen blow-up results turn every integral f-B_i/2 into that exceptional component.",
                "Exact cdd enumeration gives 247 primitive CM cusps; 82 and only 82 split over both boundaries, producing 164 distinct additional curves.",
            ],
        },
        "full_nef_chamber": {
            "hyperbolic_dimension": 11,
            "facet_count": len(full["walls"]),
            "labels": full["labels"],
            "NS_wall_vectors": full["physical_NS_wall_vectors"],
            "negative_wall_gram_sha256": full["negative_wall_gram_sha256"],
            "wall_rank": full["rank"],
            "basis_facet_positions_zero_based": full["basis_positions"],
            "basis_gram_inertia": list(full["gram_inertia"]),
            "full_gram_nullity": full["full_gram_nullity"],
            "norm_histogram": {str(key): value for key, value in sorted(full["norm_histogram"].items())},
            "Coxeter_relation_histogram": dict(sorted(full["relation_histogram"].items())),
            "strict_interior_vector": _rows(full["strict_interior"].T)[0],
            "strict_interior_square": int(full["strict_interior_square"]),
            "strict_interior_minimum_wall_pairing": int(full["strict_interior_minimum_wall_pairing"]),
            "all_facets_irredundant": full["all_facets_irredundant"],
            "facet_witness_minimum_square": str(full["facet_witness_minimum_square"]),
            "facet_witness_minimum_other_wall_margin": str(full["facet_witness_minimum_other_wall_margin"]),
            "exact_double_description": match["exact_cone"],
            "literal_vinberg": match["literal_vinberg"],
            "classification": "finite_volume_noncompact_H11_206",
        },
        "arithmetic": {
            "lattice": "2*<NS,B1/2,B2/2> = U + E8(-2) + A1(-2)^2",
            "nikulin_invariants": {"r": 12, "a": 10, "delta": 1},
            "basis_in_NS_coordinates": _rational_rows(arithmetic["basis"]),
            "gram": _rows(arithmetic["gram"]),
            "determinant": arithmetic["determinant"],
            "smith_diagonal": arithmetic["smith"],
            "arithmetic_mirror_vectors": _rows(arithmetic["root_coordinates"].T),
            "mirror_span_smith_diagonal": arithmetic["root_span_smith"],
            "primitive_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(arithmetic["root_types"].items())
            },
            "all_206_reflections_integral_involutive_isometries": arithmetic["all_reflections_integral"],
            "NS_nonintegral_reflection_positions_zero_based": arithmetic["NS_nonintegral_reflection_positions"],
        },
        "exact_Fermat_polytope_identity": match,
        "lineage": {
            "surface": "known: Kondo Type MI",
            "cover": "known: supersingular Fermat quartic, Artin invariant 1",
            "abstract_polytope": "known in project as F3-C2-H11-206",
            "arithmetic_lattice_and_maximal_group": "known: Alexeev (12,10,1)",
            "project_gain": "first complete Type-MI full-nef realization and exact quotient/K3 chamber identification",
            "worldwide_literal_priority": "open",
            "inverse_track_score": "N/A/L0: this certificate starts from the surface",
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cdd", default=None)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate(args.cdd)
    if args.check_archive:
        if not args.output.exists() or json.loads(args.output.read_text(encoding="utf-8")) != certificate:
            raise SystemExit("archived Type-MI H11 certificate differs")
        print("archived Type-MI H11 certificate matches")
        return
    args.output.write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(args.output)


if __name__ == "__main__":
    main()
