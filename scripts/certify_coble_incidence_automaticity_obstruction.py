#!/usr/bin/env python3
"""Certify the two main automatic-incidence obstructions for Coble carriers.

The cheap/default certificate is self contained apart from two pinned project
archives.  ``--recompute-cones`` reruns the exact cdd calculations for the
four projectivity stabilizers of a four-point set on P1; the full S4 row takes
roughly two minutes on a laptop.  ``--recompute-cm-rays`` reruns the much
smaller cdd calculation showing that both one-boundary 12-root chambers have
three, not one, non-simple ideal vertices.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
import shutil
from typing import Any, Iterable, Sequence

import sympy as sp

from scripts.certify_coble_one_boundary_full_ns_h10 import _cdd_rays
from scripts.certify_coble_polytope_first_blowdown_incidence_bridge import (
    _cone_summary,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/coble_incidence_automaticity_obstruction.json"
MII = ROOT / "data/enumerations/coble_type_mii_full_ns_h17.json"
ONE_BOUNDARY = ROOT / "data/enumerations/coble_one_boundary_full_ns_h10.json"
SCHEMA = "coxeter.coble_incidence_automaticity_obstruction.v1"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _partitions(total: int, maximum: int | None = None) -> Iterable[tuple[int, ...]]:
    if total == 0:
        yield ()
        return
    upper = min(total, maximum if maximum is not None else total)
    for first in range(upper, 0, -1):
        for rest in _partitions(total - first, first):
            yield (first, *rest)


def _sextic_rows() -> list[dict[str, Any]]:
    rows = []
    for degrees in _partitions(6):
        boundary_count = len(degrees)
        internal_nodes = sum((d - 1) * (d - 2) // 2 for d in degrees)
        pairwise_points = sum(
            degrees[i] * degrees[j]
            for i in range(boundary_count)
            for j in range(i)
        )
        blowups = internal_nodes + pairwise_points
        component_squares = []
        for d in degrees:
            other_intersections = d * (6 - d)
            delta = (d - 1) * (d - 2) // 2
            component_squares.append(d * d - 4 * delta - other_intersections)
        rows.append(
            {
                "degrees": list(degrees),
                "boundary_count": boundary_count,
                "internal_node_count": internal_nodes,
                "pairwise_intersection_count": pairwise_points,
                "blowup_count": blowups,
                "expected_blowup_count_9_plus_b": 9 + boundary_count,
                "component_self_intersections": component_squares,
                "K_square": 9 - blowups,
                "hyperbolic_dimension": 9 + boundary_count,
                "finite_Aut_boundary_count_survives_Kondo": boundary_count
                in {1, 2, 5, 8},
                "all_components_can_be_smooth": all(d <= 2 for d in degrees),
            }
        )
    return rows


def _cross_ratio_stabilizers() -> dict[str, list[tuple[int, ...]]]:
    lam = sp.symbols("lambda")
    points = [
        sp.Matrix([1, 0]),  # infinity
        sp.Matrix([0, 1]),  # zero
        sp.Matrix([1, 1]),  # one
        sp.Matrix([lam, 1]),
    ]

    def determinant(left: sp.Matrix, right: sp.Matrix) -> sp.Expr:
        return left[0] * right[1] - left[1] * right[0]

    def cross_ratio(order: Sequence[int]) -> sp.Expr:
        a, b, c, d = (points[index] for index in order)
        return sp.cancel(
            determinant(a, c)
            * determinant(b, d)
            / (determinant(a, d) * determinant(b, c))
        )

    base = cross_ratio((0, 1, 2, 3))
    permutations = list(itertools.permutations(range(4)))
    stabilizers: dict[str, list[tuple[int, ...]]] = {}
    for kind in ("generic", "harmonic", "equianharmonic", "characteristic_3"):
        accepted = []
        for permutation in permutations:
            numerator = sp.together(cross_ratio(permutation) - base).as_numer_denom()[0]
            if kind == "generic":
                keep = sp.Poly(numerator, lam).is_zero
            elif kind == "harmonic":
                keep = sp.expand(numerator.subs(lam, -1)) == 0
            elif kind == "equianharmonic":
                keep = sp.rem(
                    sp.Poly(numerator, lam), sp.Poly(lam**2 - lam + 1, lam)
                ).is_zero
            else:
                keep = sp.Poly(numerator, lam, modulus=3).eval(-1) % 3 == 0
            if keep:
                accepted.append(permutation)
        stabilizers[kind] = accepted
    return stabilizers


def _permutation_label(permutation: Sequence[int]) -> str:
    return "".join(str(value + 1) for value in permutation)


def _archived_grid_cones() -> dict[str, dict[str, Any]]:
    return {
        "generic": {
            "group": "V4",
            "graph_wall_count": 4,
            "wall_count": 28,
            "extreme_ray_count": 5678,
            "finite": 2000,
            "ideal": 1178,
            "hyperideal": 2500,
            "primitive_ray_list_sha256": "06e43b07ff892c26daed2f7f006af6de4594be7e237e2f04256bec4726cde847",
        },
        "harmonic": {
            "group": "D4",
            "graph_wall_count": 8,
            "wall_count": 32,
            "extreme_ray_count": 15748,
            "finite": 10024,
            "ideal": 2380,
            "hyperideal": 3344,
            "primitive_ray_list_sha256": "3014969c410438f8c41476410a2a75d96e4229489afec8121b69bcf6c2be6b09",
        },
        "equianharmonic": {
            "group": "A4",
            "graph_wall_count": 12,
            "wall_count": 36,
            "extreme_ray_count": 13702,
            "finite": 5360,
            "ideal": 3002,
            "hyperideal": 5340,
            "primitive_ray_list_sha256": "85c5a07dac12d6e50488e99b4fba081ef4928cee73b0675cdb462961c8386159",
        },
        "characteristic_3": {
            "group": "S4",
            "graph_wall_count": 24,
            "wall_count": 48,
            "extreme_ray_count": 196052,
            "finite": 190808,
            "ideal": 5244,
            "hyperideal": 0,
            "primitive_ray_list_sha256": "04d6df9b251c6a264060eb559d21ee51b80c294fcc95e264f75856a648a4da57",
        },
    }


def _recompute_grid_cones(
    stabilizers: dict[str, list[tuple[int, ...]]], executable: str
) -> dict[str, dict[str, Any]]:
    payload = json.loads(MII.read_text())
    gram = sp.Matrix(payload["full_NS"]["gram"])
    walls = [sp.Matrix(row) for row in payload["facets"]["physical_NS_wall_vectors"]]
    permutations = list(itertools.permutations(range(4)))
    summaries = {}
    for kind, selected in stabilizers.items():
        chosen = [*walls[:24]]
        chosen.extend(walls[24 + permutations.index(permutation)] for permutation in selected)
        summary, _ = _cone_summary(gram, chosen, executable, check_facets=False)
        summaries[kind] = summary
    return summaries


def _cm_ray_rows(executable: str) -> dict[str, dict[str, Any]]:
    payload = json.loads(ONE_BOUNDARY.read_text())
    output = {}
    for name, case in payload["cases"].items():
        gram = sp.Matrix(case["lattice_gram"])
        walls = [sp.Matrix(row) for row in case["wall_vectors"]]
        boundary = walls[0]
        roots = walls[1:13]
        cm_basis = sp.Matrix.hstack(*(boundary.T * gram).nullspace())
        cm_gram = cm_basis.T * gram * cm_basis
        root_coordinates = [cm_basis.gauss_jordan_solve(root)[0] for root in roots]
        inequalities = sp.Matrix.hstack(*root_coordinates).T * cm_gram
        rays = _cdd_rays(inequalities, executable)
        histogram: Counter[int] = Counter()
        kind_histogram: Counter[str] = Counter()
        for ray in rays:
            vector = sp.Matrix(ray)
            active = sum(value == 0 for value in inequalities * vector)
            square = (vector.T * cm_gram * vector)[0]
            kind = "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
            histogram[active] += 1
            kind_histogram[f"{kind}_active_{active}"] += 1
        output[name] = {
            "ray_count": len(rays),
            "active_facet_histogram": {
                str(key): value for key, value in sorted(histogram.items())
            },
            "kind_and_active_histogram": dict(sorted(kind_histogram.items())),
            "non_simple_ideal_vertex_count": kind_histogram["ideal_active_10"],
        }
    return output


def build_certificate(
    *, recompute_cones: bool = False, recompute_cm_rays: bool = False
) -> dict[str, Any]:
    sextics = _sextic_rows()
    stabilizers = _cross_ratio_stabilizers()
    grid_cones = _archived_grid_cones()
    executable = shutil.which("cddexec_gmp")
    if (recompute_cones or recompute_cm_rays) and executable is None:
        raise RuntimeError("cddexec_gmp is required for exact cone recomputation")
    if recompute_cones:
        fresh = _recompute_grid_cones(stabilizers, executable)
        for kind, summary in fresh.items():
            grid_cones[kind].update(
                {
                    key: summary[key]
                    for key in (
                        "wall_count",
                        "extreme_ray_count",
                        "finite",
                        "ideal",
                        "hyperideal",
                        "primitive_ray_list_sha256",
                    )
                }
            )

    cm_rays = {
        "I": {
            "ray_count": 23,
            "active_facet_histogram": {"9": 20, "10": 3},
            "kind_and_active_histogram": {
                "finite_active_9": 14,
                "ideal_active_9": 6,
                "ideal_active_10": 3,
            },
            "non_simple_ideal_vertex_count": 3,
        },
        "II": {
            "ray_count": 55,
            "active_facet_histogram": {"9": 52, "10": 3},
            "kind_and_active_histogram": {
                "finite_active_9": 38,
                "ideal_active_9": 14,
                "ideal_active_10": 3,
            },
            "non_simple_ideal_vertex_count": 3,
        },
    }
    if recompute_cm_rays:
        cm_rays = _cm_ray_rows(executable)

    for kind, permutations in stabilizers.items():
        grid_cones[kind]["permutations"] = [
            _permutation_label(permutation) for permutation in permutations
        ]
        grid_cones[kind]["finite_volume"] = grid_cones[kind]["hyperideal"] == 0

    kondo_plane_survivors = [
        row["degrees"]
        for row in sextics
        if row["finite_Aut_boundary_count_survives_Kondo"]
    ]
    smooth_partitions = [row["degrees"] for row in sextics if row["all_components_can_be_smooth"]]

    checks = {
        "there_are_11_degree_partitions_of_six": len(sextics) == 11,
        "every_sextic_cluster_has_N_equals_9_plus_b": all(
            row["blowup_count"] == row["expected_blowup_count_9_plus_b"]
            for row in sextics
        ),
        "every_boundary_component_has_square_minus_four": all(
            set(row["component_self_intersections"]) == {-4} for row in sextics
        ),
        "every_sextic_cluster_has_K_square_minus_b": all(
            row["K_square"] == -row["boundary_count"] for row in sextics
        ),
        "Kondo_prunes_plane_partitions_to_five": kondo_plane_survivors
        == [[6], [5, 1], [4, 2], [3, 3], [2, 1, 1, 1, 1]],
        "smooth_component_partitions_are_exactly_four": smooth_partitions
        == [[2, 2, 2], [2, 2, 1, 1], [2, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]],
        "four_point_stabilizer_orders_are_4_8_12_24": [
            len(stabilizers[kind])
            for kind in ("generic", "harmonic", "equianharmonic", "characteristic_3")
        ]
        == [4, 8, 12, 24],
        "only_characteristic_3_grid_cone_has_no_hyperideal_ray": [
            kind for kind, row in grid_cones.items() if row["finite_volume"]
        ]
        == ["characteristic_3"],
        "both_known_H9_12_root_chambers_have_three_non_simple_ideal_vertices": all(
            row["non_simple_ideal_vertex_count"] == 3 for row in cm_rays.values()
        ),
    }

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "direction": "POLYTOPE_TO_SURFACE_ESSENTIAL_SURJECTIVITY_GATE_AUDIT",
        "result": (
            "The natural toric/Looijenga, smooth plane-sextic, and complete "
            "P1xP1-grid factories do not yield a non-singleton automatic finite-volume "
            "Coble denominator. The grid lane has an exact finite-volume iff: among "
            "the four projectivity stabilizers it survives only for S4 in characteristic 3."
        ),
        "pinned_inputs": {
            "type_MII_full_NS": {"path": str(MII.relative_to(ROOT)), "sha256": _sha256(MII)},
            "one_boundary_full_NS": {
                "path": str(ONE_BOUNDARY.relative_to(ROOT)),
                "sha256": _sha256(ONE_BOUNDARY),
            },
        },
        "plane_rational_sextic_factory": {
            "rows": sextics,
            "strong_Coble_argument": (
                "For a reduced union of rational nodal plane curves of total degree 6, "
                "blow up every internal node and every pairwise intersection. Then "
                "N=9+b, every strict component has square -4, and their reduced sum "
                "B equals -2K. Since B.B_i=-4, every member of |-2K| contains every "
                "B_i, hence |-2K|={B}; if -K were effective, twice it would equal the "
                "reduced B, impossible."
            ),
            "Kondo_finite_Aut_surviving_degree_partitions": kondo_plane_survivors,
            "surviving_hyperbolic_dimensions": [10, 11, 14],
            "smooth_component_degree_partitions": smooth_partitions,
            "smooth_component_finite_volume_consequence": (
                "A smooth rational plane component has degree 1 or 2. Kondo's "
                "finite-Aut classification therefore leaves only b=5, dimension 14, "
                "and Type VI in characteristic 3; this lane has at most one exact "
                "finite-volume isometry type, not a strict non-singleton class."
            ),
            "compactness": (
                "Every branch is noncompact: after the first point blowup, pull back "
                "the pencil of lines through that point. Its primitive nef class has square 0."
            ),
        },
        "ordinary_simple_point_interpolation_obstruction": {
            "degree_d_projective_dimension": "d(d+3)/2",
            "minus_two_simple_points_required": "d^2+2",
            "minus_two_excess_conditions": "(d^2-3d+4)/2 > 0 for every integer d",
            "minus_one_simple_points_required": "d^2+1",
            "minus_one_expected_dimension": "-(d-1)(d-2)/2",
            "conclusion": (
                "No additional (-2)-curve is forced on the ordinary general-point "
                "locus. Finite-volume root stocks require special tangency, collinearity, "
                "finite-field, or other determinantal incidence."
            ),
        },
        "P1xP1_complete_grid": {
            "anti_bicanonical_equation": (
                "rF1+sF2-2 sum E_ij = -2K = 4F1+4F2-2 sum E_ij"
            ),
            "forced_grid": "r=s=4; b=8; Picard rank 18; hyperbolic dimension 17",
            "four_point_projectivity_stabilizers_and_exact_cones": grid_cones,
            "theorem": (
                "For the complete-projectivity grid family in characteristic not 2, "
                "finite volume is equivalent to full S4 projectivity, which is equivalent "
                "to characteristic 3 and cross-ratio lambda=-1=2. The unique survivor is "
                "the known MII H17/48 chamber."
            ),
        },
        "Tumarkin_n_plus_3_route": {
            "published_scope": (
                "Tumarkin classifies the compact n+3 case and noncompact pyramids, "
                "and proves the noncompact dimension bound n<=16. The general "
                "noncompact finite-volume n+3 classification is not complete."
            ),
            "Roberts_scope": (
                "Roberts classifies precisely one non-simple vertex; he explicitly "
                "leaves simple polytopes and nonpyramidal polytopes with more than one "
                "non-simple vertex."
            ),
            "known_D1_split_CM_chambers": cm_rays,
            "verdict": (
                "Both known H9/12 root chambers have three non-simple ideal vertices, "
                "so neither lies in Roberts' classified stratum. No published exhaustive "
                "list presently proves that D_1_split contains only Types I and II."
            ),
        },
        "toric_and_Looijenga_no_go": (
            "A toric boundary and, by definition, a Looijenga-pair boundary are effective "
            "anticanonical divisors. Boundary blowups preserve an effective -K divisor. "
            "They therefore cannot directly be a strong Coble carrier, which requires "
            "|-K| empty, without an additional quotient/cover/degeneration operation."
        ),
        "universality_warning": (
            "Rank-three matroid/point-line realization spaces have Mnëv universality. "
            "Consequently local intersection counts, expected dimension, or a fixed "
            "Enriques diagram cannot force nonemptiness for arbitrary plane-incidence "
            "compiler branches; one needs a globally realizable structural matroid class "
            "or an exact incidence certificate."
        ),
        "next_finite_computation": {
            "one_boundary": (
                "Enumerate H9/12 Gale diagrams in the restricted simply-laced integral "
                "E10 stratum with Gram rank 10/corank 2, three non-simple ideal vertices, "
                "and complete cusp lifts; solve the finitely many determinant systems "
                "for dotted-edge integers and run literal Vinberg tests."
            ),
            "plane": (
                "Start from the five Kondo-surviving sextic degree partitions [6], "
                "[5,1], [4,2], [3,3], [2,1,1,1,1]; enumerate multiplicity-two "
                "contact/proximity diagrams and low-degree special curves, then solve "
                "their saturated incidence ideals before running D_CobVin return."
            ),
            "grid": "closed negatively as a source of a non-singleton class",
        },
        "strict_contract_score": {
            "new_non_singleton_class_with_forced_carrier": "not obtained",
            "exact_structural_criterion": (
                "obtained for the complete P1xP1 four-point grid, but it has one "
                "finite-volume survivor"
            ),
            "compact_implication": "none; every rational blowup branch here is noncompact",
            "L4_status": "still open",
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--recompute-cones", action="store_true")
    parser.add_argument("--recompute-cm-rays", action="store_true")
    args = parser.parse_args()
    certificate = build_certificate(
        recompute_cones=args.recompute_cones,
        recompute_cm_rays=args.recompute_cm_rays,
    )
    if args.check_archive:
        archived = json.loads(OUTPUT.read_text())
        if certificate != archived:
            raise SystemExit("archive mismatch")
        print("archive OK")
    else:
        print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
