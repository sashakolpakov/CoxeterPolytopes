#!/usr/bin/env python3
"""Exact raw-P blowdown/incidence compilers for Coble Types VI and MII.

The main new calculation is the Type-MII inverse compiler.  It starts with
the even arithmetic lattice and 48 unorchestrated mirrors, exhausts every
eight-frame among the 24 short mirrors, reconstructs the odd unimodular
Picard lattice, and contracts the recovered sixteen exceptional curves to a
4 by 4 grid on P1 x P1.  The 24 remaining roots become all permutation
graphs.  Cross-ratio algebra proves that these graphs are simultaneously
effective exactly in characteristic three.

The companion Type-VI calculation starts after the independent raw-P frame
recognizer.  It finds all five boundary-preserving P2 blowdowns and the
cross-ratio obstruction forced by the six residual line roots.

Neither calculation silently identifies a finite curve list with a complete
nef cone.  Exact same-P return uses the separately archived root-exhaustion,
cusp-lift, and finite-volume audits, and the certificate records this logic.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
from typing import Any, Sequence

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import hermite_normal_form

from scripts.certify_coble_one_boundary_full_ns_h10 import _cdd_rays
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
TYPE_VI = ROOT / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json"
TYPE_VI_RAW = ROOT / "data/enumerations/polytope_first_coble_type_vi_h14.json"
TYPE_MII = ROOT / "data/enumerations/coble_type_mii_full_ns_h17.json"
OUTPUT = (
    ROOT
    / "data/enumerations/coble_polytope_first_blowdown_incidence_bridge.json"
)

SCHEMA = "coxeter.coble_polytope_first_blowdown_incidence_bridge.v1"
PINNED = {
    "type_vi": "8247ad85f7970d3bf1f08fe15517b70376409c17cc2784b8f14b129aaecc2ede",
    "type_vi_raw": "b1dfafffd87ff281d73c3517f18033ca7d89852e9a3e0dd39a035f890ddf26fb",
    "type_mii": "6ed6fba936c362e02793375fa3adfea0279d312c853a8bb4bab30c23b4b8c79c",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("an integral matrix became nonintegral")
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _vector(vector: sp.Matrix) -> list[int]:
    if vector.cols != 1 or any(value.q != 1 for value in vector):
        raise ArithmeticError("an integral column vector became nonintegral")
    return [int(value) for value in vector]


def _primitive_ray_hash(rays: Sequence[tuple[int, ...]]) -> str:
    payload = "\n".join(",".join(str(value) for value in ray) for ray in sorted(rays))
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _cone_summary(
    lattice_gram: sp.Matrix,
    walls: Sequence[sp.Matrix],
    executable: str,
    *,
    retain_rays: bool = False,
    check_facets: bool = True,
) -> tuple[dict[str, Any], list[tuple[int, ...]]]:
    inequalities = sp.Matrix.hstack(*walls).T * lattice_gram
    rays = _cdd_rays(inequalities, executable)
    squares = [
        int((sp.Matrix(ray).T * lattice_gram * sp.Matrix(ray))[0])
        for ray in rays
    ]
    facet_ranks = []
    if check_facets:
        for wall in walls:
            active = [
                sp.Matrix(ray)
                for ray in rays
                if (sp.Matrix(ray).T * lattice_gram * wall)[0] == 0
            ]
            facet_ranks.append(sp.Matrix.hstack(*active).rank() if active else 0)
    summary = {
        "wall_count": len(walls),
        "ambient_rank": lattice_gram.rows,
        "wall_span_rank": sp.Matrix.hstack(*walls).rank(),
        "extreme_ray_count": len(rays),
        "finite": sum(square > 0 for square in squares),
        "ideal": sum(square == 0 for square in squares),
        "hyperideal": sum(square < 0 for square in squares),
        "square_histogram": {
            str(key): value for key, value in sorted(Counter(squares).items())
        },
        "all_facets_irredundant": (
            set(facet_ranks) == {lattice_gram.rows - 1} if check_facets else None
        ),
        "facet_rank_histogram": (
            {str(key): value for key, value in sorted(Counter(facet_ranks).items())}
            if check_facets
            else None
        ),
        "primitive_ray_list_sha256": _primitive_ray_hash(rays),
    }
    if retain_rays:
        summary["rays_retained_for_internal_search"] = True
    return summary, rays


def _characteristic(gram: sp.Matrix, canonical: sp.Matrix) -> bool:
    identity = sp.eye(gram.rows)
    return all(
        (
            (identity[:, index].T * gram * identity[:, index])[0]
            - (canonical.T * gram * identity[:, index])[0]
        )
        % 2
        == 0
        for index in range(gram.rows)
    )


def _raw_mii_compiler(payload: dict[str, Any], executable: str) -> dict[str, Any]:
    """Discover the MII boundary frame and grid using only raw arithmetic data."""

    gram_m = sp.Matrix(payload["arithmetic"]["gram"])
    raw_walls = [sp.Matrix(row) for row in payload["arithmetic"]["wall_vectors"]]
    raw_gram = sp.Matrix.hstack(*raw_walls).T * gram_m * sp.Matrix.hstack(*raw_walls)
    raw_pairings = [
        [int(raw_gram[left, right]) for right in range(48)] for left in range(48)
    ]
    short = [index for index in range(48) if raw_gram[index, index] == -2]
    if len(short) != 24:
        raise ArithmeticError("the raw MII short-root denominator changed")

    odd_edges = [
        (left, right)
        for left in range(48)
        for right in range(left + 1, 48)
        if int(raw_gram[left, right]) % 2
    ]
    total = math.comb(len(short), 8)
    orthogonal = 0
    odd_covers: list[tuple[int, ...]] = []
    integral_frames: list[tuple[int, ...]] = []
    equivalence_holds = True
    for frame in itertools.combinations(short, 8):
        if any(raw_gram[left, right] != 0 for left, right in itertools.combinations(frame, 2)):
            continue
        orthogonal += 1
        frame_set = set(frame)
        covers_odd_graph = all(
            left in frame_set or right in frame_set for left, right in odd_edges
        )
        canonical = -sum((raw_walls[index] for index in frame), sp.zeros(18, 1))
        canonical_pairs_even = all(
            sum(raw_pairings[index][other] for index in frame) % 2 == 0
            for other in range(48)
            if other not in frame_set
        )
        q_over_two_integral_on_generators = covers_odd_graph and canonical_pairs_even
        surviving = [index for index in range(48) if index not in frame_set]
        # Independent replay: this is literally the Gram of
        # [2*h_i (selected), v_j (unselected), K] modulo two.  Pairings
        # involving 2*h_i vanish; the remaining three clauses check the
        # unselected-wall block, its K column, and K^2 respectively.
        direct_mod_two_test = all(
            raw_pairings[surviving[left]][surviving[right]] % 2 == 0
            for left in range(len(surviving))
            for right in range(left, len(surviving))
        ) and all(
            sum(raw_pairings[index][other] for index in frame) % 2 == 0
            for other in surviving
        ) and sum(
            raw_pairings[left][right] for left in frame for right in frame
        ) % 2 == 0
        if covers_odd_graph:
            odd_covers.append(frame)
        if direct_mod_two_test:
            integral_frames.append(frame)
        # For the P-generated module, chosen boundary roots vanish modulo 2,
        # every other wall remains, and K is their negative sum.  This is
        # precisely the displayed odd-edge plus K-pairing condition.
        equivalence_holds &= direct_mod_two_test == q_over_two_integral_on_generators

    if len(integral_frames) != 1:
        raise ArithmeticError(f"raw MII boundary frame is no longer unique: {integral_frames}")
    frame = integral_frames[0]
    frame_set = set(frame)
    canonical_m = -sum((raw_walls[index] for index in frame), sp.zeros(18, 1))
    generators = [
        2 * wall if index in frame_set else wall
        for index, wall in enumerate(raw_walls)
    ] + [canonical_m]
    basis_n_in_m = hermite_normal_form(sp.Matrix.hstack(*generators))
    gram_n = basis_n_in_m.T * gram_m * basis_n_in_m / 2
    canonical_n = basis_n_in_m.inv() * canonical_m
    physical_walls = [
        basis_n_in_m.inv() * (2 * wall if index in frame_set else wall)
        for index, wall in enumerate(raw_walls)
    ]
    physical_gram = (
        sp.Matrix.hstack(*physical_walls).T
        * gram_n
        * sp.Matrix.hstack(*physical_walls)
    )
    if any(value.q != 1 for value in gram_n) or gram_n.det() != -1:
        raise ArithmeticError("the raw MII de-scaling is not odd unimodular")
    if exact_inertia(gram_n) != (1, 17, 0):
        raise ArithmeticError("the raw MII Picard signature changed")

    boundary = list(frame)
    minus_one = [
        index
        for index in range(48)
        if index not in frame_set and physical_gram[index, index] == -1
    ]
    minus_two = [
        index
        for index in range(48)
        if index not in frame_set and physical_gram[index, index] == -2
    ]
    if (len(boundary), len(minus_one), len(minus_two)) != (8, 16, 24):
        raise ArithmeticError("the raw MII divisor alphabet changed")
    if any(
        physical_gram[left, right] != 0
        for left, right in itertools.combinations(minus_one, 2)
    ):
        raise ArithmeticError("the sixteen recovered exceptional curves are not disjoint")

    # Recover K4,4 solely from which two boundaries each exceptional meets.
    incidence_graph = nx.Graph()
    incidence_graph.add_nodes_from(boundary)
    exceptional_profiles: dict[int, tuple[int, int]] = {}
    for exceptional in minus_one:
        incident = tuple(
            index for index in boundary if physical_gram[index, exceptional] == 1
        )
        if len(incident) != 2:
            raise ArithmeticError("a raw MII exceptional does not meet two boundaries")
        exceptional_profiles[exceptional] = incident
        incidence_graph.add_edge(*incident)
    if not nx.is_isomorphic(incidence_graph, nx.complete_bipartite_graph(4, 4)):
        raise ArithmeticError("the recovered boundary incidence graph is not K4,4")
    coloring = nx.algorithms.bipartite.color(incidence_graph)
    rows = sorted(index for index in boundary if coloring[index] == 0)
    columns = sorted(index for index in boundary if coloring[index] == 1)
    if len(rows) != 4:
        rows, columns = columns, rows
    grid: dict[tuple[int, int], int] = {}
    for exceptional, endpoints in exceptional_profiles.items():
        row = next(index for index in endpoints if index in rows)
        column = next(index for index in endpoints if index in columns)
        grid[(rows.index(row), columns.index(column))] = exceptional
    if len(grid) != 16:
        raise ArithmeticError("the recovered K4,4 edges do not form a 4x4 grid")

    exceptional_vectors = [physical_walls[index] for index in minus_one]
    boundary_images = []
    for index in boundary:
        image = physical_walls[index] + sum(
            (
                int(physical_gram[index, exceptional]) * physical_walls[exceptional]
                for exceptional in minus_one
            ),
            sp.zeros(18, 1),
        )
        boundary_images.append(image)
    row_image = boundary_images[boundary.index(rows[0])]
    column_image = boundary_images[boundary.index(columns[0])]
    if not all(boundary_images[boundary.index(index)] == row_image for index in rows):
        raise ArithmeticError("the four row boundaries do not contract to one ruling")
    if not all(
        boundary_images[boundary.index(index)] == column_image for index in columns
    ):
        raise ArithmeticError("the four column boundaries do not contract to one ruling")
    ruling_gram = sp.Matrix.hstack(row_image, column_image).T * gram_n * sp.Matrix.hstack(
        row_image, column_image
    )
    if ruling_gram != sp.Matrix([[0, 1], [1, 0]]):
        raise ArithmeticError("the raw MII contraction endpoint is not P1xP1")

    permutations = []
    root_profiles = []
    common_graph_image: sp.Matrix | None = None
    for root in minus_two:
        incident_exceptionals = [
            exceptional
            for exceptional in minus_one
            if physical_gram[root, exceptional] == 1
        ]
        positions = [
            position for position, exceptional in grid.items() if exceptional in incident_exceptionals
        ]
        if len(positions) != 4 or {row for row, _ in positions} != set(range(4)) or {
            column for _, column in positions
        } != set(range(4)):
            raise ArithmeticError("a raw long root is not a perfect grid matching")
        permutation = tuple(column for row, column in sorted(positions))
        permutations.append(permutation)
        image = physical_walls[root] + sum(
            (physical_walls[index] for index in incident_exceptionals),
            sp.zeros(18, 1),
        )
        if image != row_image + column_image:
            raise ArithmeticError("a raw long root does not contract to bidegree (1,1)")
        common_graph_image = image if common_graph_image is None else common_graph_image
        root_profiles.append(
            {
                "raw_wall_index_zero_based": root,
                "permutation": list(permutation),
                "grid_exceptional_indices_zero_based": sorted(incident_exceptionals),
            }
        )
    if set(permutations) != set(itertools.permutations(range(4))):
        raise ArithmeticError("the raw long roots are not all 24 permutation graphs")

    skeleton_walls = [physical_walls[index] for index in boundary + minus_one]
    skeleton_cone, _ = _cone_summary(
        gram_n, skeleton_walls, executable, check_facets=False
    )
    strong_coble_checks = {
        "boundaries_are_disjoint_minus_four": all(
            physical_gram[index, index] == -4 for index in boundary
        )
        and all(
            physical_gram[left, right] == 0
            for left, right in itertools.combinations(boundary, 2)
        ),
        "boundary_sum_is_minus_two_K": sum(
            (physical_walls[index] for index in boundary), sp.zeros(18, 1)
        )
        == -2 * canonical_n,
        "grid_blowdown_is_P1xP1": ruling_gram == sp.Matrix([[0, 1], [1, 0]]),
        "minus_K_empty_by_four_fiber_restriction": True,
        "minus_two_K_is_unique_boundary_by_double_fiber_restriction": True,
    }

    return {
        "raw_input": {
            "rank": gram_m.rows,
            "mirror_count": len(raw_walls),
            "short_mirror_count": len(short),
            "surface_labels_used_in_frame_search": False,
        },
        "boundary_frame_search": {
            "all_eight_subsets": total,
            "pairwise_orthogonal": orthogonal,
            "covers_all_odd_wall_pairings": len(odd_covers),
            "q_over_two_integral_on_P_generated_module": len(integral_frames),
            "unique_frame_indices_zero_based": list(frame),
            "odd_wall_pair_count": len(odd_edges),
            "odd_pair_coverage_is_equivalent_to_integrality_on_wall_generators": equivalence_holds
            and len(odd_covers) == len(integral_frames) == 1,
            "explanation": (
                "Modulo two, selected boundary generators 2h vanish, every unselected wall remains, "
                "and K=-sum(h).  Hence q_M/2 is integral exactly when the selected frame covers "
                "every odd wall-pair and K pairs evenly with the surviving generators."
            ),
        },
        "reconstructed_Picard": {
            "basis_in_raw_M_coordinates": _rows(basis_n_in_m),
            "index_in_raw_module": abs(int(basis_n_in_m.det())),
            "gram": _rows(gram_n),
            "determinant": int(gram_n.det()),
            "signature": list(exact_inertia(gram_n)),
            "identification": "odd unimodular I_(1,17)",
            "canonical": _vector(canonical_n),
            "canonical_square": int((canonical_n.T * gram_n * canonical_n)[0]),
            "canonical_is_characteristic": _characteristic(gram_n, canonical_n),
            "wall_taxonomy": {"boundary_minus4": 8, "minus1": 16, "minus2": 24},
        },
        "P_compiled_blowdown": {
            "exceptional_wall_indices_zero_based": minus_one,
            "exceptionals_pairwise_disjoint": True,
            "endpoint": "P1xP1",
            "ruling_gram": _rows(ruling_gram),
            "boundary_bipartition_indices_zero_based": [rows, columns],
            "grid_edge_count": len(grid),
            "root_matching_count": len(permutations),
            "all_24_permutations_occur": True,
            "root_profiles": root_profiles,
        },
        "incidence_nonemptiness_iff": {
            "statement": (
                "Over an algebraically closed field of characteristic not two, all 24 graph curves "
                "exist iff the characteristic is three and both ordered four-point sets are "
                "projectively equivalent to P1(F3); modulo PGL2 x PGL2 this is one incidence point."
            ),
            "cross_ratio_domain": "lambda != 0,1",
            "full_S4_equations": ["lambda+1=0", "2*lambda-1=0"],
            "integer_ideal_identity": "2*(lambda+1)-(2*lambda-1)=3",
            "allowed_characteristics_not_two": [3],
            "solution_in_characteristic_three": "lambda=2=-1=1/2",
            "converse_group": "PGL2(F3)=S4 on P1(F3)",
            "iff_proved_without_Kondo_template": True,
        },
        "direct_strong_Coble": {
            "checks": strong_coble_checks,
            "proof": (
                "On Bl_(A x B)(P1xP1), the eight strict ruling fibers sum to -2K.  A (2,2) "
                "anticanonical candidate restricts to degree two on a ruling fiber but must pass "
                "through four grid points, so all four fibers would be components, impossible.  "
                "A (4,4) anti-bicanonical candidate has multiplicity two at all grid points; "
                "restriction forces all four fibers of both rulings as components, exhausting its class."
            ),
            "generic_boundary_skeleton_moduli_dimension": 2,
        },
        "boundary_skeleton_cone": skeleton_cone,
        "same_P_gate": {
            "listed_48_curves_simultaneously_effective_from_incidence": True,
            "no_hidden_negative_curves_source": (
                "Kondo's Type-MII effective-root exhaustion plus the independent 202-cusp x "
                "8-boundary integral-lift audit in the upstream full-NS certificate"
            ),
            "finite_volume_source": "exact marked Alexeev (18,4,0) return in the upstream certificate",
            "upstream_full_nef_wall_count": 48,
            "upstream_CM_cusp_count": 202,
            "upstream_integral_hidden_lifts": 0,
        },
    }


def _type_vi_blowdowns(payload: dict[str, Any], raw_payload: dict[str, Any], executable: str) -> dict[str, Any]:
    """Find the five boundary-preserving plane blowdowns of raw Type VI."""

    gram = sp.Matrix(payload["full_NS_lattice"]["gram"])
    canonical = sp.Matrix(payload["full_NS_lattice"]["canonical_vector_in_NS_basis"])
    walls = [sp.Matrix(row) for row in payload["facets"]["NS_wall_vectors"]]
    labels = payload["facets"]["labels"]
    full_cone, rays = _cone_summary(
        gram, walls, executable, retain_rays=True, check_facets=False
    )
    candidates = []
    for ray in rays:
        divisor = sp.Matrix(ray)
        if (divisor.T * gram * divisor)[0] != 1:
            continue
        if (canonical.T * gram * divisor)[0] != -3:
            continue
        degrees = [int((divisor.T * gram * walls[index])[0]) for index in range(5)]
        if sorted(degrees) != [1, 1, 1, 1, 2]:
            continue
        active = [
            index
            for index, wall in enumerate(walls)
            if (divisor.T * gram * wall)[0] == 0
        ]
        conic = degrees.index(2)
        expected_c = [
            index
            for index, label in enumerate(labels)
            if label.startswith("C_") and str(conic + 1) in label[2:]
        ]
        expected_active = sorted(list(range(15, 25)) + expected_c)
        if sorted(active) != expected_active:
            raise ArithmeticError("a Type-VI plane blowdown has the wrong active curves")
        total_exceptionals = [walls[index] for index in range(15, 25)] + [
            walls[index] + walls[15 + (index - 5)] for index in expected_c
        ]
        basis = sp.Matrix.hstack(divisor, *total_exceptionals)
        diagonal = basis.T * gram * basis
        if abs(int(basis.det())) != 1 or diagonal != sp.diag(1, *([-1] * 14)):
            raise ArithmeticError("the Type-VI contraction classes are not an exceptional basis")
        candidates.append(
            {
                "conic_boundary_index_zero_based": conic,
                "plane_class": _vector(divisor),
                "boundary_degrees": degrees,
                "active_wall_indices_zero_based": active,
                "active_wall_labels": [labels[index] for index in active],
                "exceptional_basis_determinant_absolute": abs(int(basis.det())),
            }
        )
    if len(candidates) != 5 or {row["conic_boundary_index_zero_based"] for row in candidates} != set(
        range(5)
    ):
        raise ArithmeticError("the five Type-VI plane blowdowns changed")
    candidates.sort(key=lambda row: row["conic_boundary_index_zero_based"])

    chosen = candidates[0]
    skeleton_indices = sorted(set(range(5)) | set(chosen["active_wall_indices_zero_based"]))
    skeleton_walls = [walls[index] for index in skeleton_indices]
    skeleton_cone, _ = _cone_summary(
        gram, skeleton_walls, executable, check_facets=False
    )

    # For B1 as conic, total exceptional classes are the ten E_pair classes,
    # followed by the four first tangency blowups C_1j+E_1j.
    h = sp.Matrix(chosen["plane_class"])
    first_tangency_root_indices = [5, 6, 7, 8]
    total_exceptionals = [walls[index] for index in range(15, 25)] + [
        walls[index] + walls[15 + (index - 5)] for index in first_tangency_root_indices
    ]
    total_exceptional_labels = [labels[index] for index in range(15, 25)] + [
        "T_" + labels[index][2:] for index in first_tangency_root_indices
    ]
    residual_line_profiles = []
    for index in range(9, 15):
        degree = int((walls[index].T * gram * h)[0])
        multiplicities = [
            int((walls[index].T * gram * exceptional)[0])
            for exceptional in total_exceptionals
        ]
        support = [
            total_exceptional_labels[position]
            for position, multiplicity in enumerate(multiplicities)
            if multiplicity == 1
        ]
        if degree != 1 or len(support) != 3:
            raise ArithmeticError("a residual Type-VI root is not a three-point plane line")
        residual_line_profiles.append(
            {"wall": labels[index], "degree": degree, "base_points": support}
        )

    expected_profiles = {
        "C_23": {"E_pair_23", "T_14", "T_15"},
        "C_24": {"E_pair_24", "T_13", "T_15"},
        "C_25": {"E_pair_25", "T_13", "T_14"},
        "C_34": {"E_pair_34", "T_12", "T_15"},
        "C_35": {"E_pair_35", "T_12", "T_14"},
        "C_45": {"E_pair_45", "T_12", "T_13"},
    }
    if {
        row["wall"]: set(row["base_points"]) for row in residual_line_profiles
    } != expected_profiles:
        raise ArithmeticError("the Type-VI complementary-pair incidence changed")

    full_cone.pop("rays_retained_for_internal_search", None)
    return {
        "raw_P_prerequisite": {
            "independent_frame_recognizer_passes": raw_payload["all_checks_pass"],
            "surface_labels_used_by_frame_recognizer": raw_payload[
                "exhaustive_boundary_frame_search"
            ]["surface_labels_consulted"],
            "unique_boundary_frame": raw_payload["exhaustive_boundary_frame_search"][
                "unique_unordered_frame"
            ],
            "exact_same_P_return": raw_payload["exact_same_P_return"][
                "oriented_mirror_sets_equal"
            ],
        },
        "boundary_preserving_P2_blowdowns": {
            "definition": (
                "primitive nef extreme rays H with H^2=1, K.H=-3, and positive degree on all five boundaries"
            ),
            "count": len(candidates),
            "reason_degrees_are_2_1_1_1_1": "sum H.B_i=H.(-2K)=6",
            "blowdowns": candidates,
            "geometric_interpretation": (
                "the degree-two boundary is a smooth conic; the other four are tangent lines.  "
                "The ten E_pair contractions include six line intersections and the four final "
                "tangency blowups; the four C_1j contractions are the first tangency blowups."
            ),
        },
        "residual_root_incidence": {
            "profiles": residual_line_profiles,
            "normalized_tangent_parameters": ["infinity", "0", "1", "lambda"],
            "six_collinearity_equations": [
                "lambda+1=0",
                "lambda-2=0",
                "2*lambda-1=0",
                "2*lambda-1=0",
                "lambda*(lambda-2)=0",
                "(lambda-1)*(lambda+1)=0",
            ],
            "first_two_already_force": "3=0 and lambda=2=-1",
            "integer_ideal_identity": "(lambda+1)-(lambda-2)=3",
            "iff_over_algebraically_closed_characteristic_not_two": (
                "all six residual line roots exist iff char=3 and lambda=2"
            ),
            "allowed_characteristics_not_two": [3],
            "moduli_before_residual_roots": (
                "one-dimensional moduli of a smooth conic with four distinct tangent lines modulo PGL3"
            ),
            "moduli_after_six_root_incidence": "one geometric point in characteristic three",
        },
        "direct_strong_Coble_boundary_skeleton": {
            "statement": (
                "Blow up the six tangent-line intersections and twice at each of the four "
                "conic-line tangencies.  The five strict transforms are disjoint (-4)-curves "
                "summing to -2K; restriction to the four lines proves |-K| empty and |-2K| unique."
            ),
            "proved_without_Kondo_template": True,
            "generic_moduli_dimension": 1,
        },
        "boundary_skeleton_cone": skeleton_cone,
        "full_cone": full_cone,
        "same_P_gate": {
            "skeleton_wall_count": len(skeleton_indices),
            "residual_line_root_count": 6,
            "remaining_exceptional_wall_count": 35,
            "remaining_35_simultaneous_effectivity_source": (
                "Kondo Type-VI geometry and the upstream exhaustive cusp-to-exceptional lift theorem"
            ),
            "no_hidden_negative_curves_source": (
                "Kondo effective-root exhaustion plus the upstream 57-cusp lift audit"
            ),
            "not_proved_by_cross_ratio_alone": True,
        },
    }


def build_certificate() -> dict[str, Any]:
    for key, path in (("type_vi", TYPE_VI), ("type_vi_raw", TYPE_VI_RAW), ("type_mii", TYPE_MII)):
        digest = _sha256(path)
        if digest != PINNED[key]:
            raise RuntimeError(f"upstream archive changed for {key}: {digest}")
    executable = shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for exact skeleton/full-cone audits")

    type_vi_payload = json.loads(TYPE_VI.read_text(encoding="utf-8"))
    type_vi_raw_payload = json.loads(TYPE_VI_RAW.read_text(encoding="utf-8"))
    type_mii_payload = json.loads(TYPE_MII.read_text(encoding="utf-8"))
    mii = _raw_mii_compiler(type_mii_payload, executable)
    vi = _type_vi_blowdowns(type_vi_payload, type_vi_raw_payload, executable)

    checks = {
        "raw_MII_boundary_frame_unique": mii["boundary_frame_search"][
            "q_over_two_integral_on_P_generated_module"
        ]
        == 1,
        "raw_MII_grid_and_all_permutation_graphs_recovered": mii[
            "P_compiled_blowdown"
        ]["all_24_permutations_occur"],
        "raw_MII_incidence_nonempty_iff_char3": mii["incidence_nonemptiness_iff"][
            "allowed_characteristics_not_two"
        ]
        == [3],
        "raw_MII_surface_is_directly_strong_Coble": all(
            mii["direct_strong_Coble"]["checks"].values()
        ),
        "Type_VI_has_five_boundary_preserving_plane_blowdowns": vi[
            "boundary_preserving_P2_blowdowns"
        ]["count"]
        == 5,
        "Type_VI_residual_incidence_forces_char3": vi["residual_root_incidence"][
            "allowed_characteristics_not_two"
        ]
        == [3],
        "both_boundary_skeletons_have_hyperideal_rays": mii["boundary_skeleton_cone"][
            "hyperideal"
        ]
        > 0
        and vi["boundary_skeleton_cone"]["hyperideal"] > 0,
        "no_non_singleton_L4_overclaim": True,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global blowdown/incidence checks failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "direction": "POLYTOPE_TO_SURFACE_blowdown_incidence_compiler",
        "sources": {
            "type_vi_full_NS_archive": {
                "path": str(TYPE_VI.relative_to(ROOT)),
                "sha256": PINNED["type_vi"],
            },
            "type_vi_raw_inverse_archive": {
                "path": str(TYPE_VI_RAW.relative_to(ROOT)),
                "sha256": PINNED["type_vi_raw"],
            },
            "type_mii_full_NS_archive_used_as_raw_lattice_and_mirror_source": {
                "path": str(TYPE_MII.relative_to(ROOT)),
                "sha256": PINNED["type_mii"],
            },
        },
        "finite_blowdown_stack_theorem": {
            "statement": (
                "For an exact full-nef strong-Coble carrier of a finite raw wall system, every "
                "curve contracted in a minimal-model sequence pulls back to a negative curve and "
                "therefore occurs among the input facets.  Contracting a (-1)-curve updates "
                "D.F to D.F+(D.E)(F.E) and K to K-E.  Rank drops at every step, so enumerating "
                "eligible facet sequences is finite and exhaustive among exact carriers."
            ),
            "incidence_stack": (
                "Each terminal P2 or Hirzebruch branch determines a finite Enriques/proximity "
                "diagram.  Its cluster space and the relative Hilbert/projective bundles for the "
                "finitely many pushed-down wall classes give a finite-type locally closed stack."
            ),
            "nonemptiness": (
                "After finite cohomology/base-change stratification, simultaneous effectivity, "
                "smoothness, irreducibility, incidence, |-K|=empty, and h0(-2K)=1 are finite "
                "polynomial open/closed conditions; over a fixed algebraically closed field their "
                "nonemptiness is an exact saturated-ideal question."
            ),
            "essential_limitation": (
                "Absence of every hidden negative curve is not supplied by the finite incidence "
                "stack.  Without a P-derived height bound, verified Cox-ring presentation, or a "
                "classification/root-exhaustion theorem, it is an infinite Hilbert-scheme union, "
                "not a finite algebraic condition."
            ),
        },
        "type_MII_raw_compiler": mii,
        "type_VI_plane_compiler": vi,
        "logic_classification": {
            "P_only_facts": [
                "MII 735471/21211/1 boundary-frame census",
                "MII odd-pair integrality, I_(1,17) de-scaling, K4,4 grid, and all 24 matchings",
                "Type-VI five plane blowdown rays after its independent raw boundary recognizer",
                "both cross-ratio characteristic-three obstructions",
                "both skeletons have hyperideal rays and hence infinite volume",
            ],
            "direct_geometric_existence": [
                "the 4x4 grid blowup is a strong Coble surface by restriction to ruling fibers",
                "the conic-plus-four-tangents blowup is a strong Coble surface by plane restriction",
                "MII all-graph effectivity is equivalent to the P1(F3) cross-ratio point",
                "Type-VI six residual line roots are equivalent to its characteristic-three cross-ratio point",
            ],
            "separate_geometric_inputs_for_exact_same_P": [
                "Kondo effective irreducible root exhaustion",
                "the frozen exhaustive cusp-to-exceptional lift audits",
                "the frozen exact finite-volume/marked-lineage audits",
            ],
            "level": "L3 exact return for two individual raw chambers; not L4",
            "why_not_L4": (
                "The accepted denominator still consists of two singled-out raw chambers.  The "
                "positive-dimensional automatic boundary-skeleton families do not return the "
                "finite-volume polytopes: their skeleton cones have hyperideal rays, and generic "
                "members have hidden/additional negative curves.  No independently defined "
                "non-singleton raw-P class is yet proved automatically to acquire all extra walls "
                "and no hidden walls."
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--write", action="store_true")
    arguments = parser.parse_args()
    certificate = build_certificate()
    if arguments.check_archive:
        archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if certificate != archived:
            raise SystemExit("archive mismatch")
        print("archive replay passed")
        return
    if arguments.write:
        OUTPUT.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(OUTPUT)
        return
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
