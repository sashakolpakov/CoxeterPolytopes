#!/usr/bin/env python3
"""Certify the full-NS H11 chamber of the two-boundary Type-I Coble surface.

This is deliberately a SURFACE -> POLYTOPE certificate.  Kondō's complete
Type-I Coble--Mukai root chamber and Halphen-fibration theorem are used to
recover every negative curve in the full Neron--Severi lattice.  The output
does not assert a converse polytope-to-surface theorem.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import tarfile
from typing import Any, Sequence

import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.domains import ZZ

from coxeter_polytopes.exact.matrix import inertia
from coxeter_polytopes.verify.vinberg_volume import (
    _coerce_wall_gram,
    _is_incident_ideal_extension,
    _parabolic_vertex_neighbourhoods,
)
from scripts.certify_coble_one_boundary_full_ns_h10 import (
    _cdd_rays,
    _finite_special_masks,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _euler_from_connected_components,
    _relation_matrix,
)
from scripts.certify_shimada_disc16_d5_h10_lineage import (
    _connected_elliptic_masks,
)


ROOT = Path(__file__).resolve().parents[1]
ENRIQUES = ROOT / "data/known_examples/enriques_finite_automorphism_nef_chambers.json"
ENRIQUES_RAYS = ROOT / "data/enumerations/enriques_finite_automorphism_nef_chamber_rays.jsonl.gz"
BASELINE = ROOT / "data/enumerations/positive_example_lineage_header_controls.json"
KNOWN_CLASS = ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_lineage.json"
OUTPUT = ROOT / "data/enumerations/coble_type_i_two_boundary_full_ns_h11.json"

COXITER = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ZIP = Path("/private/tmp/coxiter-current.zip")
DOUBLING = Path("graphs/testing_index2.tar.gz")

SCHEMA = "coxeter.coble_type_i_two_boundary_full_ns_h11.v1"
PINNED = {
    "enriques": "245715b7352219a7ed1bfc34a43b5b37921e289d39829df51015c347b97072bc",
    "enriques_rays": "b72a2a5479a773c8151ed92cd8b698e71fe4b8cb0b4bf527245847df39cb46c5",
    "baseline": "58e00b546f50a9f3a3bcd4a21e16414aad6bb0d65f471475b6bfc199a57421aa",
    "known_class": "e1b20c853d3f2f254900c0dd46cdf187093292fb4356e1d6e9de23ca839faf89",
    "coxiter_zip": "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3",
    "doubling": "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.is_Integer is not True for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return [[int(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]


def _rational_rows(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(value) if value.q == 1 else str(value) for value in matrix.row(i)]
        for i in range(matrix.rows)
    ]


def _primitive(vector: sp.Matrix) -> sp.Matrix:
    denominator = sp.ilcm(*(sp.denom(value) for value in vector))
    values = [int(value * denominator) for value in vector]
    divisor = math.gcd(*(abs(value) for value in values))
    result = sp.Matrix([value // divisor for value in values])
    return result


def _basis_from_half_integral_generators(generators: Sequence[sp.Matrix]) -> sp.Matrix:
    doubled = sp.Matrix(
        [
            [int(2 * generators[column][row]) for column in range(len(generators))]
            for row in range(12)
        ]
    )
    if doubled.rank() != 12:
        raise ArithmeticError("the proposed divisor generators do not have full rank")
    return hermite_normal_form(doubled) / 2


def _load_type_i() -> tuple[sp.Matrix, list[sp.Matrix], list[sp.Matrix], dict[str, Any]]:
    source = json.loads(ENRIQUES.read_text(encoding="utf-8"))
    record = next(row for row in source["types"] if row["nikulin_kondo_type"] == "I")
    e10 = sp.Matrix(source["lattice"]["gram"])
    roots = [sp.Matrix(row) for row in record["wall_vectors_in_L10_basis"]]
    root_matrix = sp.Matrix(record["wall_vectors_in_L10_basis"])
    pivots = [value - 1 for value in record["exact_cone"]["coordinate_basis_facet_positions"]]
    cusp_basis = root_matrix[pivots, :]
    cusps: list[sp.Matrix] = []
    with gzip.open(ENRIQUES_RAYS, "rt", encoding="utf-8") as handle:
        next(handle)
        for line in handle:
            payload = json.loads(line)
            if payload["type"] != "I" or payload["kind"] != "ideal":
                continue
            cusp = _primitive(cusp_basis.T * sp.Matrix(payload["ray"]))
            if any(value < 0 for value in root_matrix * e10 * cusp):
                cusp = -cusp
            if (cusp.T * e10 * cusp)[0] != 0:
                raise ArithmeticError("a stored Type-I cusp is not isotropic")
            cusps.append(cusp)
    if len(roots) != 12 or len(cusps) != 9:
        raise ArithmeticError("the complete Type-I CM chamber changed")
    # The two circle vertices in Kondō Figure 2 are exactly the last two:
    # they have only even pairings with every remaining root.
    circle = [
        index
        for index, root in enumerate(roots)
        if all(
            int((root.T * e10 * roots[other])[0]) % 2 == 0
            for other in range(12)
            if other != index
        )
    ]
    if circle != [10, 11]:
        raise ArithmeticError("the two Type-I (-1)-root vertices changed")
    return e10, roots, cusps, record


def _surface_lattices(
    e10: sp.Matrix, roots: Sequence[sp.Matrix]
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix, sp.Matrix, dict[str, Any]]:
    ambient = sp.zeros(12)
    ambient[:10, :10] = e10
    ambient[10, 10] = ambient[11, 11] = -1
    h1 = sp.Matrix([0] * 10 + [1, 0])
    h2 = sp.Matrix([0] * 10 + [0, 1])
    boundaries = [2 * h1, 2 * h2]
    actual_roots = [sp.Matrix([*roots[index], 0, 0]) for index in range(10)]
    crossing_exceptionals = [
        sp.Matrix([*(roots[index] / 2), sp.Rational(-1, 2), sp.Rational(-1, 2)])
        for index in (10, 11)
    ]
    canonical = -h1 - h2
    curve_generators = [*boundaries, *actual_roots, *crossing_exceptionals, canonical]
    curve_basis = _basis_from_half_integral_generators(curve_generators)
    curve_gram = curve_basis.T * ambient * curve_basis
    if abs(int(curve_gram.det())) != 4:
        raise ArithmeticError("the visible curve lattice no longer has determinant four")

    # Any full NS of a smooth rational surface is unimodular.  Enumerate all
    # integral index-two overlattices of the visible curve lattice.  There is
    # exactly one, so no geometric glue choice is being guessed.
    glues: list[sp.Matrix] = []
    glue_bits: list[list[int]] = []
    for bits in itertools.product((0, 1), repeat=12):
        if not any(bits):
            continue
        column = sp.Matrix(bits)
        if any(int(value) % 2 for value in curve_gram * column):
            continue
        if sp.denom((column.T * curve_gram * column)[0] / 4) != 1:
            continue
        glues.append(curve_basis * column / 2)
        glue_bits.append(list(bits))
    if len(glues) != 1:
        raise ArithmeticError("the full NS glue is no longer unique")
    ns_basis = _basis_from_half_integral_generators([*curve_basis.columnspace(), glues[0]])
    ns_gram = ns_basis.T * ambient * ns_basis
    if not (
        ns_gram.det() == -1
        and str(inertia(_rows(ns_gram))) == "Inertia(positive=1, negative=11, zero=0)"
        and any(int(ns_gram[i, i]) % 2 for i in range(12))
    ):
        raise ArithmeticError("the reconstructed NS is not I_(1,11)")

    # The boundary reflections need not preserve NS.  The canonical arithmetic
    # module is Pic~ = <NS,B1/2,B2/2>, with the form doubled to make it integral.
    pic_basis = _basis_from_half_integral_generators([*ns_basis.columnspace(), h1])
    pic_gram = 2 * pic_basis.T * ambient * pic_basis
    pic_smith = smith_normal_form(pic_gram, domain=ZZ)
    smith = [abs(int(pic_smith[i, i])) for i in range(12)]
    dual = pic_gram.inv()
    delta = int(any(sp.denom(dual[i, i]) != 1 for i in range(12)))
    if not (
        pic_gram.det() == -1024
        and all(int(pic_gram[i, i]) % 2 == 0 for i in range(12))
        and smith == [1, 1, *([2] * 10)]
        and delta == 1
    ):
        raise ArithmeticError("the doubled extended Picard lattice changed class")
    audit = {
        "visible_curve_lattice_determinant": int(curve_gram.det()),
        "unique_integral_index_two_overlattice": True,
        "unique_glue_bits_in_visible_basis": glue_bits[0],
        "unique_glue_in_CM_h1_h2_coordinates": [
            int(value) if value.q == 1 else str(value) for value in glues[0]
        ],
        "NS_identification": "odd unimodular I_(1,11)",
        "doubled_extended_Picard_nikulin_invariants": {"r": 12, "a": 10, "delta": 1},
        "doubled_extended_Picard_identification": "U + E8(-2) + <-2>^2",
    }
    return ambient, ns_basis, ns_gram, pic_basis, pic_gram, audit


def _complete_walls(
    roots: Sequence[sp.Matrix],
    cusps: Sequence[sp.Matrix],
    ns_basis: sp.Matrix,
) -> tuple[list[str], list[sp.Matrix], dict[str, Any]]:
    inverse = ns_basis.inv()
    h1 = sp.Matrix([0] * 10 + [1, 0])
    h2 = sp.Matrix([0] * 10 + [0, 1])
    labels = ["B1", "B2"]
    walls = [2 * h1, 2 * h2]
    labels.extend(f"R{index + 1}" for index in range(10))
    walls.extend(sp.Matrix([*roots[index], 0, 0]) for index in range(10))
    labels.extend(("E_cross_1", "E_cross_2"))
    walls.extend(
        sp.Matrix([*(roots[index] / 2), sp.Rational(-1, 2), sp.Rational(-1, 2)])
        for index in (10, 11)
    )
    split_cusps: list[int] = []
    membership: list[list[bool]] = []
    for cusp_index, cusp in enumerate(cusps):
        row: list[bool] = []
        candidates = []
        for boundary_index in range(2):
            candidate = sp.Matrix([*cusp, 0, 0])
            candidate[10 + boundary_index] = -1
            integral = all(sp.denom(value) == 1 for value in inverse * candidate)
            row.append(integral)
            candidates.append(candidate)
        membership.append(row)
        if row == [True, True]:
            split_cusps.append(cusp_index)
            for boundary_index, candidate in enumerate(candidates):
                labels.append(f"E_cusp_{cusp_index + 1}_B{boundary_index + 1}")
                walls.append(candidate)
        elif row != [False, False]:
            raise ArithmeticError("a cusp split the two boundary integrality tests")
    if split_cusps != [1, 3, 4, 5, 7, 8] or len(walls) != 26:
        raise ArithmeticError("the complete two-boundary wall count changed")
    coordinates = [inverse * wall for wall in walls]
    if any(any(sp.denom(value) != 1 for value in wall) for wall in coordinates):
        raise ArithmeticError("a full-nef wall left NS")
    return labels, coordinates, {
        "CM_cusp_count": len(cusps),
        "boundary_split_membership_by_CM_cusp": membership,
        "split_CM_cusps_zero_based": split_cusps,
        "together_CM_cusps_zero_based": [index for index in range(9) if index not in split_cusps],
        "negative_curve_count": 26,
        "negative_curve_types": {"boundary_minus4": 2, "minus2": 10, "minus1": 14},
        "minus1_decomposition": {"meeting_both_boundaries_once": 2, "meeting_one_boundary_twice": 12},
    }


def _cone_audit(
    gram: sp.Matrix, walls: Sequence[sp.Matrix], executable: str
) -> dict[str, Any]:
    matrix = sp.Matrix.hstack(*walls)
    wall_gram = matrix.T * gram * matrix
    pivots = list(wall_gram.rref()[1])
    if len(pivots) != 12:
        raise ArithmeticError("the full walls do not span rank twelve")
    inequalities = wall_gram[:, pivots]
    rays = _cdd_rays(inequalities, executable)
    basis_gram = wall_gram.extract(pivots, pivots)
    counts: Counter[str] = Counter()
    face_rows: list[list[list[int]]] = [[] for _ in walls]
    for ray in rays:
        vector = sp.Matrix(ray)
        slacks = inequalities * vector
        square = int((vector.T * basis_gram * vector)[0])
        counts["finite" if square > 0 else "ideal" if square == 0 else "hyperideal"] += 1
        for facet, slack in enumerate(slacks):
            if slack != 0 or len(face_rows[facet]) == 11:
                continue
            trial = [*face_rows[facet], list(ray)]
            if sp.Matrix(trial).rank() > len(face_rows[facet]):
                face_rows[facet].append(list(ray))
    ray_hash = hashlib.sha256(
        json.dumps(sorted(rays), separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "coordinate_basis_facets_zero_based": pivots,
        "extreme_ray_count": len(rays),
        "finite_vertex_count": counts["finite"],
        "ideal_vertex_count": counts["ideal"],
        "hyperideal_vertex_count": counts["hyperideal"],
        "all_facet_face_ranks": [len(rows) for rows in face_rows],
        "all_facets_irredundant": all(len(rows) == 11 for rows in face_rows),
        "primitive_ray_list_sha256": ray_hash,
    }


def _literal_vinberg(negative_gram: sp.Matrix) -> dict[str, Any]:
    positive = [
        [-int(negative_gram[i, j]) for j in range(negative_gram.cols)]
        for i in range(negative_gram.rows)
    ]
    relations = _relation_matrix(positive)
    connected, search = _connected_elliptic_masks(positive)
    special = _finite_special_masks(relations, connected, maximum_rank=11)
    solomon = _euler_from_connected_components(relations, connected)
    matrix = _coerce_wall_gram(positive)
    ideal_vertices, affine_audits = _parabolic_vertex_neighbourhoods(matrix, 11)
    finite_vertices = special[11]
    endpoints: Counter[str] = Counter()
    bad_edges = 0
    for edge_mask in special[10]:
        edge = tuple(index for index in range(negative_gram.rows) if edge_mask & (1 << index))
        finite_count = sum(edge_mask & ~vertex == 0 for vertex in finite_vertices)
        ideal_count = sum(
            _is_incident_ideal_extension(edge, endpoint) for endpoint in ideal_vertices
        )
        endpoints["F" * finite_count + "I" * ideal_count] += 1
        bad_edges += finite_count + ideal_count != 2
    return {
        "elliptic_edge_count": len(special[10]),
        "finite_vertex_count": len(finite_vertices),
        "ideal_vertex_count": len(ideal_vertices),
        "connected_affine_component_audit_count": len(affine_audits),
        "endpoint_type_histogram": dict(sorted(endpoints.items())),
        "bad_edge_count": bad_edges,
        "bad_parabolic_neighborhood_count": sum(
            not audit.is_rank_dimension_minus_one_parabolic for audit in affine_audits
        ),
        "finite_special_subset_counts_by_rank": [len(special[rank]) for rank in range(12)],
        "connected_elliptic_search": search,
        "rational_coxeter_euler_characteristic": solomon["rational_euler_characteristic"],
    }


def _arithmetic(
    ambient: sp.Matrix,
    ns_basis: sp.Matrix,
    ns_gram: sp.Matrix,
    pic_basis: sp.Matrix,
    pic_gram: sp.Matrix,
    ns_walls: Sequence[sp.Matrix],
) -> dict[str, Any]:
    h1 = sp.Matrix([0] * 10 + [1, 0])
    h2 = sp.Matrix([0] * 10 + [0, 1])
    physical = [ns_basis * wall for wall in ns_walls]
    pic_normals = [h1, h2, *physical[2:]]
    pic_walls = [pic_basis.inv() * wall for wall in pic_normals]
    if any(any(sp.denom(value) != 1 for value in wall) for wall in pic_walls):
        raise ArithmeticError("an arithmetic wall normal left doubled Pic~")

    def reflection_integral(gram: sp.Matrix, wall: sp.Matrix) -> bool:
        square = (wall.T * gram * wall)[0]
        reflection = sp.eye(12) - sp.Rational(2, square) * wall * (wall.T * gram)
        return (
            all(value.is_Integer is True for value in reflection)
            and reflection.T * gram * reflection == gram
            and reflection**2 == sp.eye(12)
        )

    ns_integral = [reflection_integral(ns_gram, wall) for wall in ns_walls]
    pic_integral = [reflection_integral(pic_gram, wall) for wall in pic_walls]
    matrix = sp.Matrix.hstack(*pic_walls)
    smith_matrix = smith_normal_form(matrix, domain=ZZ)
    smith = [abs(int(smith_matrix[i, i])) for i in range(12)]
    root_types: Counter[tuple[int, int]] = Counter()
    for wall in pic_walls:
        square = int((wall.T * pic_gram * wall)[0])
        divisibility = math.gcd(*(abs(int(value)) for value in pic_gram * wall))
        root_types[(square, divisibility)] += 1
    return {
        "arithmetic_lattice": "2*Pic~(S) = U + E8(-2) + <-2>^2, Nikulin (12,10,1)",
        "arithmetic_lattice_basis_in_CM_h1_h2_coordinates": _rational_rows(pic_basis),
        "arithmetic_lattice_gram": _rows(pic_gram),
        "arithmetic_wall_vectors": _rows(sp.Matrix.hstack(*pic_walls).T),
        "primitive_square_divisibility_histogram": {
            f"{square},{divisibility}": count
            for (square, divisibility), count in sorted(root_types.items())
        },
        "wall_span_smith_diagonal": smith,
        "walls_span_whole_arithmetic_lattice": smith == [1] * 12,
        "all_26_reflections_integral_on_arithmetic_lattice": all(pic_integral),
        "NS_reflection_integrality": ns_integral,
        "only_boundary_reflections_fail_on_NS": [
            index for index, valid in enumerate(ns_integral) if not valid
        ] == [0, 1],
        "boundary_denominator_explanation": (
            "B_i^2=-4 and E_cross.B_i=1, so s_Bi is not integral on NS; "
            "B_i/2 is integral in Pic~, and doubling the form makes the full "
            "26-reflection system integral"
        ),
    }


def _parse_header(text: str) -> tuple[int, int]:
    line = next(
        content for raw in text.splitlines() if (content := raw.split("#", 1)[0].strip())
    )
    facets, dimension = map(int, line.split()[:2])
    return dimension, facets


def _lineage() -> dict[str, Any]:
    if _sha256(BASELINE) != PINNED["baseline"] or _sha256(KNOWN_CLASS) != PINNED["known_class"]:
        raise ArithmeticError("a frozen project-lineage input changed")
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    prior_h11 = [row for row in baseline["ledger"]["rows"] if row["dimension"] == 11]
    known = json.loads(KNOWN_CLASS.read_text(encoding="utf-8"))
    if not (
        all(row["facets"] != 26 for row in prior_h11)
        and known["fixed_lattice_identification"]["nikulin_invariants"]
        == {"r": 12, "a": 10, "delta": 1}
        and known["known_maximal_reflection_chamber"]["source"]["paper_statements"]["facet_root_count"] == 14
    ):
        raise ArithmeticError("the bounded project lineage changed")
    if _sha256(COXITER_ZIP) != PINNED["coxiter_zip"]:
        raise ArithmeticError("the pinned CoxIter archive changed")
    doubling_path = COXITER / DOUBLING
    if _sha256(doubling_path) != PINNED["doubling"]:
        raise ArithmeticError("the pinned doubling archive changed")
    ordinary: Counter[tuple[int, int]] = Counter()
    ordinary_count = 0
    for path in COXITER.rglob("*.coxiter"):
        ordinary[_parse_header(path.read_text(encoding="utf-8"))] += 1
        ordinary_count += 1
    doubling: Counter[tuple[int, int]] = Counter()
    doubling_count = 0
    with tarfile.open(doubling_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError("a CoxIter doubling member disappeared")
            doubling[_parse_header(handle.read().decode("utf-8"))] += 1
            doubling_count += 1
    if not (
        ordinary_count == 1214
        and doubling_count == 2959
        and ordinary[(11, 26)] == 0
        and doubling[(11, 26)] == 0
    ):
        raise ArithmeticError("the bounded H11/26 CoxIter census changed")
    return {
        "known_arithmetic_lattice": True,
        "known_full_reflection_overgroup": "Alexeev (12,10,1), 14 walls",
        "numeric_subgroup_index": "open: both H11 Coxeter Euler characteristics are zero",
        "frozen_preinsertion_project_H11_rows": len(prior_h11),
        "frozen_preinsertion_project_H11_facet_counts": sorted(row["facets"] for row in prior_h11),
        "no_preinsertion_project_H11_26_row": True,
        "coxiter_commit": "b800be48d8240fbabd2d1c9542e61054584ad7fe",
        "ordinary_record_count": ordinary_count,
        "doubling_record_count": doubling_count,
        "ordinary_H11_26_header_count": ordinary[(11, 26)],
        "doubling_H11_26_header_count": doubling[(11, 26)],
        "bounded_status": "new exact chamber in the project baseline; absent by header from both pinned CoxIter lanes",
        "worldwide_literal_priority": "open",
        "not_new": ["surface", "arithmetic lattice", "maximal reflection group", "commensurability class"],
    }


def build_certificate(cdd_executable: str | None = None) -> dict[str, Any]:
    executable = cdd_executable or shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required")
    if _sha256(ENRIQUES) != PINNED["enriques"] or _sha256(ENRIQUES_RAYS) != PINNED["enriques_rays"]:
        raise ArithmeticError("a Type-I chamber input changed")
    e10, roots, cusps, source_record = _load_type_i()
    ambient, ns_basis, ns_gram, pic_basis, pic_gram, lattice_audit = _surface_lattices(e10, roots)
    labels, walls, curve_audit = _complete_walls(roots, cusps, ns_basis)
    wall_matrix = sp.Matrix.hstack(*walls)
    negative_gram = wall_matrix.T * ns_gram * wall_matrix
    if any(negative_gram[i, j] < 0 for i in range(26) for j in range(i)):
        raise ArithmeticError("the negative curves do not define an acute chamber")
    cone = _cone_audit(ns_gram, walls, executable)
    literal = _literal_vinberg(negative_gram)
    arithmetic = _arithmetic(ambient, ns_basis, ns_gram, pic_basis, pic_gram, walls)
    lineage = _lineage()
    checks = {
        "surface_rank_gives_H11": ns_gram.rows == 12,
        "complete_curve_count_26": curve_audit["negative_curve_count"] == 26,
        "curve_norms_2_10_14": Counter(int(negative_gram[i, i]) for i in range(26))
        == Counter({-4: 2, -2: 10, -1: 14}),
        "all_pairwise_angles_are_coxeter": all(
            relation in {"1", "2", "3", "4", "0"}
            for row in _relation_matrix([[-int(negative_gram[i, j]) for j in range(26)] for i in range(26)])
            for relation in row
        ),
        "cdd_508_42_0": (
            cone["finite_vertex_count"], cone["ideal_vertex_count"], cone["hyperideal_vertex_count"]
        ) == (508, 42, 0),
        "all_26_facets_irredundant": cone["all_facets_irredundant"],
        "literal_counts": (
            literal["finite_vertex_count"],
            literal["ideal_vertex_count"],
            literal["elliptic_edge_count"],
            literal["connected_affine_component_audit_count"],
        ) == (508, 42, 3641, 106),
        "literal_endpoints": literal["endpoint_type_histogram"]
        == {"FF": 2144, "FI": 1300, "II": 197},
        "literal_zero_failures": literal["bad_edge_count"] == 0
        and literal["bad_parabolic_neighborhood_count"] == 0,
        "integral_arithmetic_reflections": arithmetic["all_26_reflections_integral_on_arithmetic_lattice"],
        "boundary_denominator_honest": arithmetic["only_boundary_reflections_fail_on_NS"],
        "bounded_lineage_honest": lineage["worldwide_literal_priority"] == "open",
    }
    if not all(checks.values()):
        raise ArithmeticError(f"two-boundary Coble H11 checks failed: {checks}")
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-08",
        "direction": "SURFACE_TO_POLYTOPE_only",
        "claim": {
            "theorem": (
                "The projectivized full Neron-Severi nef cone of Kondo's complex "
                "Type-I Coble surface with two boundary components is a 26-facet "
                "finite-volume noncompact arithmetic Coxeter polytope in H11."
            ),
            "not_claimed": [
                "polytope-to-surface converse",
                "new surface",
                "new arithmetic lattice or maximal reflection group",
                "worldwide priority for the literal H11/26 chamber",
            ],
        },
        "sources": {
            "kondo": {
                "title": "Coble surfaces with finite automorphism group",
                "arxiv": "2102.00591",
                "used": [
                    "Theorem 1.1/Table 2",
                    "Proposition 3.1",
                    "Lemmas 3.2 and 3.3",
                    "Propositions 3.7 and 3.8",
                    "Lemma 3.11",
                    "Example 4.1 and Proposition 4.2",
                    "Section 6, Type I",
                ],
            },
            "alexeev": {"title": "Reflective hyperbolic 2-elementary lattices", "arxiv": "2209.09110"},
            "input_sha256": {
                str(ENRIQUES.relative_to(ROOT)): PINNED["enriques"],
                str(ENRIQUES_RAYS.relative_to(ROOT)): PINNED["enriques_rays"],
                str(BASELINE.relative_to(ROOT)): PINNED["baseline"],
                str(KNOWN_CLASS.relative_to(ROOT)): PINNED["known_class"],
            },
        },
        "source_CM_chamber": {
            "root_count": 12,
            "cusp_count": 9,
            "wall_sublattice_index": source_record["wall_vector_checks"]["wall_sublattice_index_in_L10"],
            "circle_vertices_zero_based": [10, 11],
        },
        "full_NS_lattice": {
            "basis_in_CM_h1_h2_coordinates": _rational_rows(ns_basis),
            "gram": _rows(ns_gram),
            **lattice_audit,
        },
        "curve_exhaustion": {
            **curve_audit,
            "argument": [
                "Kondo: the only negative curves are boundaries, (-2)-curves, and (-1)-curves",
                "B.E=2 splits every (-1)-curve into incidence (1,1), (2,0), or (0,2)",
                "the (1,1) curves are exactly the two complete effective (-1)-roots",
                "a (2,0)/(0,2) curve is E=f-B_i/2 for a primitive CM cusp f",
                "the unique NS glue makes exactly six of nine cusps boundary-split, yielding twelve curves",
                "Kondo's Halphen blow-up theorem gives the converse cusp-to-exceptional-curve construction",
            ],
        },
        "facets": {
            "labels": labels,
            "NS_wall_vectors": _rows(sp.Matrix.hstack(*walls).T),
            "negative_wall_gram": _rows(negative_gram),
        },
        "exact_double_description": cone,
        "literal_vinberg": literal,
        "arithmetic": arithmetic,
        "lineage": lineage,
        "classification": "finite_volume_noncompact_full_nef_H11_26_arithmetic_Coxeter_polytope",
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate()
    if args.check_archive:
        if not args.output.exists() or json.loads(args.output.read_text(encoding="utf-8")) != certificate:
            raise SystemExit("archived two-boundary Coble H11 certificate differs")
        print("archived two-boundary Coble H11 certificate matches")
        return
    args.output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
