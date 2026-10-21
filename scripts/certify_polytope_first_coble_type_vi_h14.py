#!/usr/bin/env python3
"""Recognize the five-boundary Type-VI Coble package from raw Coxeter data.

Only the even arithmetic lattice ``M`` and its sixty oriented primitive wall
vectors are used for discovery.  In particular, this replay does not read the
surface labels, the marked boundary fold, the archived Coble--Mukai roots, or
the archived cusp lifts in the forward Type-VI certificate.

The key inverse operation is an exhaustive Coble de-scaling.  Candidate
boundary halves are enumerated among the divisibility-one short roots.  For a
five-frame H={h_i}, put K=-sum(h_i) and

    N_H = < all other main roots, 2h_1,...,2h_5,K >,
    (x,y)_N = (x,y)_M/2.

The unique frame for which N_H is an odd unimodular Picard lattice and the
walls satisfy adjunction recovers the full divisor package.  The program then
reconstructs the Coble--Mukai chamber, enumerates every primitive cusp, and
checks exact same-wall return.

Kondo's algebraic classification is invoked only after the raw package has
been recovered.  Failure is an obstruction only to this specified Coble
carrier, never to the abstract Coxeter polytope.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any, Sequence

import networkx as nx
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp


ROOT = Path(__file__).resolve().parents[1]
RAW = (
    ROOT
    / "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json"
)
OUTPUT = ROOT / "data/enumerations/polytope_first_coble_type_vi_h14.json"

SCHEMA = "coxeter.polytope_first_coble_type_vi_h14.v1"
RAW_SCHEMA = "coxeter.coble_type_vi_char3_five_boundary_full_ns_h14.v1"
RAW_SHA256 = "8247ad85f7970d3bf1f08fe15517b70376409c17cc2784b8f14b129aaecc2ede"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("an integral matrix was expected")
    return [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def _vector(vector: sp.Matrix) -> list[int]:
    if any(value.q != 1 for value in vector):
        raise ArithmeticError("an integral vector was expected")
    return [int(value) for value in vector]


def _exact_inertia(matrix: sp.Matrix) -> tuple[int, int, int]:
    """Return inertia by exact symmetric elimination."""

    work = sp.Matrix(matrix)
    positive = negative = 0
    while work.rows:
        pivot = next((index for index in range(work.rows) if work[index, index]), None)
        if pivot is None:
            pair = next(
                (
                    (left, right)
                    for left in range(work.rows)
                    for right in range(left)
                    if work[left, right]
                ),
                None,
            )
            if pair is None:
                return positive, negative, work.rows
            left, right = pair
            transform = sp.eye(work.rows)
            transform[:, left] += transform[:, right]
            work = transform.T * work * transform
            pivot = left
        order = [pivot, *range(pivot), *range(pivot + 1, work.rows)]
        work = work.extract(order, order)
        value = work[0, 0]
        positive += int(bool(value > 0))
        negative += int(bool(value < 0))
        if work.rows == 1:
            break
        column = work[1:, 0]
        work = work[1:, 1:] - column * column.T / value
    return positive, negative, 0


def _primitive_integer(vector: sp.Matrix) -> bool:
    if any(value.q != 1 for value in vector):
        return False
    return math.gcd(*(abs(int(value)) for value in vector)) == 1


def _divisibility(gram: sp.Matrix, vector: sp.Matrix) -> int:
    pairings = gram * vector
    return math.gcd(*(abs(int(value)) for value in pairings))


def _characteristic(gram: sp.Matrix, vector: sp.Matrix) -> bool:
    diagonal = sp.Matrix([gram[index, index] for index in range(gram.rows)])
    return all(int(value) % 2 == 0 for value in gram * vector - diagonal)


def _smith_diagonal(matrix: sp.Matrix) -> list[int]:
    smith = smith_normal_form(matrix, domain=ZZ)
    return [
        abs(int(smith[index, index]))
        for index in range(min(smith.rows, smith.cols))
        if smith[index, index]
    ]


def _integer_kernel_basis(matrix: sp.Matrix) -> sp.Matrix:
    """Give a saturated Z-basis of the kernel of an integral matrix."""

    diagonal, left, right = smith_normal_decomp(
        DomainMatrix.from_Matrix(matrix).convert_to(ZZ)
    )
    smith = diagonal.to_Matrix()
    left_matrix = left.to_Matrix()
    right_matrix = right.to_Matrix()
    if left_matrix * matrix * right_matrix != smith:
        raise ArithmeticError("Smith kernel decomposition failed")
    rank = matrix.rank()
    kernel = right_matrix[:, rank:]
    if matrix * kernel != sp.zeros(matrix.rows, matrix.cols - rank):
        raise ArithmeticError("Smith kernel columns are not in the kernel")
    return kernel


def _primitive_positive(vector: Sequence[Fraction]) -> tuple[int, ...]:
    denominator = math.lcm(*(value.denominator for value in vector))
    values = [value.numerator * (denominator // value.denominator) for value in vector]
    divisor = math.gcd(*(abs(value) for value in values))
    if divisor == 0:
        raise ArithmeticError("cddlib returned the zero direction")
    return tuple(value // divisor for value in values)


def _cdd_rays(inequalities: sp.Matrix, executable: str) -> list[tuple[int, ...]]:
    lines = [
        "H-representation",
        "begin",
        f"{inequalities.rows} {inequalities.cols + 1} rational",
    ]
    lines.extend(
        "0 " + " ".join(str(inequalities[row, column]) for column in range(inequalities.cols))
        for row in range(inequalities.rows)
    )
    lines.append("end")
    completed = subprocess.run(
        [executable, "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=False,
        timeout=600,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr[-4000:] or completed.stdout[-4000:])
    output = completed.stdout.splitlines()
    begin = output.index("begin")
    count, width, kind = output[begin + 1].split()
    if kind != "rational" or int(width) != inequalities.cols + 1:
        raise ArithmeticError("cddlib returned an unexpected generator format")
    rays: list[tuple[int, ...]] = []
    for line in output[begin + 2 : begin + 2 + int(count)]:
        tokens = [Fraction(token) for token in line.split()]
        if tokens[0] != 0:
            raise ArithmeticError("the homogeneous cone acquired a vertex row")
        rays.append(_primitive_positive(tokens[1:]))
    if len(rays) != len(set(rays)):
        raise ArithmeticError("cddlib returned duplicate oriented rays")
    return sorted(rays)


def _cone_audit(
    gram: sp.Matrix,
    inequalities: sp.Matrix,
    executable: str,
    expected_face_rank: int,
) -> dict[str, Any]:
    rays = _cdd_rays(inequalities, executable)
    counts: Counter[str] = Counter()
    incidence: Counter[int] = Counter()
    face_spans: list[list[sp.Matrix]] = [[] for _ in range(inequalities.rows)]
    for raw in rays:
        ray = sp.Matrix(raw)
        square = (ray.T * gram * ray)[0]
        counts["finite" if square > 0 else "ideal" if square == 0 else "hyperideal"] += 1
        slacks = inequalities * ray
        incidence[sum(value == 0 for value in slacks)] += 1
        for facet, slack in enumerate(slacks):
            if slack or len(face_spans[facet]) == expected_face_rank:
                continue
            trial = [*face_spans[facet], ray]
            matrix = sp.Matrix.hstack(*trial)
            if matrix.rank() > len(face_spans[facet]):
                face_spans[facet].append(ray)
    timelike = next(
        (
            sp.Matrix(raw)
            for raw in rays
            if (sp.Matrix(raw).T * gram * sp.Matrix(raw))[0] > 0
        ),
        None,
    )
    if timelike is None:
        raise ArithmeticError("the cone has no timelike ray to orient its causal sheet")
    time_covector = timelike.T * gram
    time_pairings = [
        sum(time_covector[column] * raw[column] for column in range(gram.cols))
        for raw in rays
    ]
    same_time_sheet = all(value > 0 for value in time_pairings)
    return {
        "ray_count": len(rays),
        "finite_vertex_count": counts["finite"],
        "ideal_vertex_count": counts["ideal"],
        "hyperideal_vertex_count": counts["hyperideal"],
        "active_facet_count_histogram": {
            str(key): value for key, value in sorted(incidence.items())
        },
        "all_facet_face_ranks": [len(vectors) for vectors in face_spans],
        "all_facets_irredundant": all(
            len(vectors) == expected_face_rank for vectors in face_spans
        ),
        "all_causal_rays_in_one_time_sheet": same_time_sheet,
        "minimum_time_orientation_pairing": int(min(time_pairings)),
        "primitive_ray_list_sha256": _payload_sha256(rays),
        "rays": rays,
    }


def _coxeter_gate(wall_gram: sp.Matrix) -> dict[str, Any]:
    allowed = {
        Fraction(0): "2",
        Fraction(1, 4): "3",
        Fraction(1, 2): "4",
        Fraction(3, 4): "6",
    }
    types: Counter[str] = Counter()
    forbidden: list[dict[str, Any]] = []
    for left in range(wall_gram.rows):
        for right in range(left):
            product = wall_gram[left, right]
            ratio = Fraction(
                int(product**2),
                int(wall_gram[left, left] * wall_gram[right, right]),
            )
            if ratio in allowed:
                types[allowed[ratio]] += 1
            elif ratio == 1:
                types["infinity"] += 1
            elif ratio > 1:
                types["ultraparallel"] += 1
            else:
                forbidden.append(
                    {
                        "facets_zero_based": [right, left],
                        "normalized_square": str(ratio),
                    }
                )
    return {
        "pair_type_histogram": dict(sorted(types.items())),
        "forbidden_pairs": forbidden,
        "all_pairs_are_coxeter": not forbidden,
    }


def _reflection_integral(gram: sp.Matrix, wall: sp.Matrix) -> bool:
    square = (wall.T * gram * wall)[0]
    reflection = sp.eye(gram.rows) - sp.Rational(2, square) * wall * (wall.T * gram)
    return (
        all(value.q == 1 for value in reflection)
        and reflection.T * gram * reflection == gram
        and reflection**2 == sp.eye(gram.rows)
    )


def _type_vi_template() -> tuple[sp.Matrix, list[str]]:
    """The labelled Gamma_VI Gram, used only as the final AG package gate."""

    pairs = list(itertools.combinations(range(5), 2))
    gram = -2 * sp.eye(20)
    for left, pair in enumerate(pairs):
        gram[left, 10 + left] = gram[10 + left, left] = 2
        for right, other in enumerate(pairs[:left]):
            if set(pair).isdisjoint(other):
                gram[left, right] = gram[right, left] = 1
            if len(set(pair) & set(other)) == 1:
                gram[10 + left, 10 + right] = gram[10 + right, 10 + left] = 1
    return gram, ["minus2"] * 10 + ["minus1_root"] * 10


def _complete_coloured_graph(gram: sp.Matrix, kinds: Sequence[str]) -> nx.Graph:
    graph = nx.Graph()
    for index, kind in enumerate(kinds):
        graph.add_node(index, kind=kind, square=int(gram[index, index]))
    for left in range(gram.rows):
        for right in range(left):
            graph.add_edge(right, left, pairing=int(gram[right, left]))
    return graph


def _type_vi_match(root_gram: sp.Matrix, kinds: Sequence[str]) -> dict[str, Any]:
    target_gram, target_kinds = _type_vi_template()
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _complete_coloured_graph(root_gram, kinds),
        _complete_coloured_graph(target_gram, target_kinds),
        node_match=lambda left, right: left == right,
        edge_match=lambda left, right: left == right,
    )
    mappings = [tuple(mapping[index] for index in range(20)) for mapping in matcher.isomorphisms_iter()]
    if not mappings:
        raise ArithmeticError("the recovered CM roots are not the labelled Type-VI package")
    return {
        "labelled_gram_isomorphic_to_Kondo_Gamma_VI": True,
        "labelled_isometry_count": len(mappings),
        "lexicographically_first_source_to_template_zero_based": list(min(mappings)),
        "template_automorphism_group": "S5" if len(mappings) == 120 else "unexpected",
    }


def _frame_search(
    main_gram: sp.Matrix,
    main_indices: Sequence[int],
    wall_coordinates: Sequence[sp.Matrix],
) -> tuple[dict[str, int], dict[str, Any]]:
    rank = main_gram.rows
    identity = sp.eye(rank)
    counts: Counter[str] = Counter()
    survivors: list[dict[str, Any]] = []

    for frame in itertools.combinations(range(rank), 5):
        counts["all_five_subsets"] += 1
        if any(main_gram[left, right] for left, right in itertools.combinations(frame, 2)):
            continue
        counts["pairwise_orthogonal"] += 1
        remaining = [index for index in range(rank) if index not in frame]
        canonical = -sum((identity[:, index] for index in frame), sp.zeros(rank, 1))
        ns_basis = sp.Matrix.hstack(
            *(
                [identity[:, index] for index in remaining]
                + [2 * identity[:, index] for index in frame[:-1]]
                + [canonical]
            )
        )
        ns_gram = ns_basis.T * main_gram * ns_basis / 2
        if any(value.q != 1 for value in ns_gram):
            continue
        counts["half_form_integral"] += 1
        if abs(int(ns_gram.det())) != 1 or all(
            int(ns_gram[index, index]) % 2 == 0 for index in range(rank)
        ):
            continue
        counts["odd_unimodular"] += 1
        inverse = ns_basis.inv()
        canonical_ns = inverse * canonical
        if (canonical_ns.T * ns_gram * canonical_ns)[0] != -5:
            continue
        counts["canonical_square_minus_five"] += 1
        if not _characteristic(ns_gram, canonical_ns):
            continue
        counts["canonical_characteristic"] += 1

        boundary_wall_indices = {main_indices[index] for index in frame}
        physical_walls: list[sp.Matrix] = []
        walls_descend = True
        for wall_index, coordinates in enumerate(wall_coordinates):
            physical = 2 * coordinates if wall_index in boundary_wall_indices else coordinates
            ns_coordinates = inverse * physical
            if any(value.q != 1 for value in ns_coordinates):
                walls_descend = False
                break
            physical_walls.append(ns_coordinates)
        if not walls_descend:
            continue
        counts["all_physical_walls_descend"] += 1

        boundaries = [inverse * (2 * identity[:, index]) for index in frame]
        taxonomy: dict[int, dict[str, Any]] = {}
        kinds: Counter[str] = Counter()
        alphabet_pass = True
        for wall_index, wall in enumerate(physical_walls):
            if wall_index in boundary_wall_indices:
                kind = "boundary_minus4"
                profile = [
                    int((wall.T * ns_gram * boundary)[0]) for boundary in boundaries
                ]
            else:
                square = int((wall.T * ns_gram * wall)[0])
                profile = [
                    int((wall.T * ns_gram * boundary)[0]) for boundary in boundaries
                ]
                ordered = sorted(profile)
                if square == -2 and ordered == [0, 0, 0, 0, 0]:
                    kind = "minus2_curve"
                elif square == -1 and ordered == [0, 0, 0, 1, 1]:
                    kind = "crossing_minus1_curve"
                elif square == -1 and ordered == [0, 0, 0, 0, 2]:
                    kind = "cusp_lift_minus1_curve"
                else:
                    alphabet_pass = False
                    break
            kinds[kind] += 1
            taxonomy[wall_index] = {
                "kind": kind,
                "boundary_intersection_profile": profile,
            }
        if not alphabet_pass:
            continue
        if any(
            (canonical_ns.T * ns_gram * wall)[0]
            != -2 - (wall.T * ns_gram * wall)[0]
            for wall in physical_walls
        ):
            continue
        counts["adjunction_alphabet"] += 1
        survivors.append(
            {
                "frame_main_positions": list(frame),
                "frame_wall_indices": sorted(boundary_wall_indices),
                "canonical_main_coordinates": canonical,
                "ns_basis_in_main_coordinates": ns_basis,
                "ns_gram": ns_gram,
                "canonical_ns_coordinates": canonical_ns,
                "boundary_ns_coordinates": boundaries,
                "physical_wall_ns_coordinates": physical_walls,
                "taxonomy": taxonomy,
                "taxonomy_counts": kinds,
            }
        )

    if len(survivors) != 1:
        raise ArithmeticError(f"expected a unique five-boundary frame, found {len(survivors)}")
    counts["surviving_frames"] = len(survivors)
    return dict(counts), survivors[0]


def _recognize(
    arithmetic_gram: sp.Matrix,
    walls: Sequence[sp.Matrix],
    executable: str,
) -> dict[str, Any]:
    rank = arithmetic_gram.rows
    if not (
        arithmetic_gram == arithmetic_gram.T
        and rank == 15
        and arithmetic_gram.cols == rank
        and _exact_inertia(arithmetic_gram) == (1, 14, 0)
        and int(arithmetic_gram.det()) == 128
    ):
        raise ArithmeticError("the raw arithmetic lattice is not the target H14 lattice")
    if len(walls) != 60 or not all(_primitive_integer(wall) for wall in walls):
        raise ArithmeticError("the raw input is not a 60-wall primitive presentation")

    wall_matrix = sp.Matrix.hstack(*walls)
    wall_gram = wall_matrix.T * arithmetic_gram * wall_matrix
    if any(
        wall_gram[left, right] < 0
        for left in range(60)
        for right in range(left)
    ):
        raise ArithmeticError("the raw walls do not have a common acute orientation")
    arithmetic_types: Counter[tuple[int, int]] = Counter()
    main_indices: list[int] = []
    for index, wall in enumerate(walls):
        square = int(wall_gram[index, index])
        divisibility = _divisibility(arithmetic_gram, wall)
        arithmetic_types[(square, divisibility)] += 1
        if (square, divisibility) == (-2, 1):
            main_indices.append(index)
    expected_types = Counter({(-2, 1): 15, (-2, 2): 35, (-4, 2): 10})
    if arithmetic_types != expected_types:
        raise ArithmeticError(f"the primitive root alphabet changed: {arithmetic_types}")
    if not all(_reflection_integral(arithmetic_gram, wall) for wall in walls):
        raise ArithmeticError("a raw wall reflection is not integral")
    coxeter = _coxeter_gate(wall_gram)
    if not coxeter["all_pairs_are_coxeter"]:
        raise ArithmeticError("the raw wall system contains a forbidden angle")

    main_basis = sp.Matrix.hstack(*(walls[index] for index in main_indices))
    if abs(int(main_basis.det())) != 1:
        raise ArithmeticError("the fifteen short divisibility-one roots are not a basis")
    main_gram = main_basis.T * arithmetic_gram * main_basis
    inverse_main = main_basis.inv()
    wall_coordinates = [inverse_main * wall for wall in walls]
    if any(any(value.q != 1 for value in coordinates) for coordinates in wall_coordinates):
        raise ArithmeticError("a wall left the main-root basis")

    frame_counts, frame = _frame_search(main_gram, main_indices, wall_coordinates)
    expected_frame_counts = {
        "all_five_subsets": 3003,
        "pairwise_orthogonal": 338,
        "half_form_integral": 1,
        "odd_unimodular": 1,
        "canonical_square_minus_five": 1,
        "canonical_characteristic": 1,
        "all_physical_walls_descend": 1,
        "adjunction_alphabet": 1,
        "surviving_frames": 1,
    }
    if frame_counts != expected_frame_counts:
        raise ArithmeticError(f"the exhaustive frame census changed: {frame_counts}")
    expected_taxonomy = Counter(
        {
            "boundary_minus4": 5,
            "minus2_curve": 10,
            "crossing_minus1_curve": 10,
            "cusp_lift_minus1_curve": 35,
        }
    )
    if frame["taxonomy_counts"] != expected_taxonomy:
        raise ArithmeticError(f"the divisor taxonomy changed: {frame['taxonomy_counts']}")

    ns_gram: sp.Matrix = frame["ns_gram"]
    ns_basis: sp.Matrix = frame["ns_basis_in_main_coordinates"]
    inverse_ns = ns_basis.inv()
    boundary_positions: list[int] = frame["frame_main_positions"]
    identity = sp.eye(15)
    half_boundaries = [identity[:, position] for position in boundary_positions]
    boundary_indices = set(frame["frame_wall_indices"])
    canonical = frame["canonical_main_coordinates"]
    boundaries_ns: list[sp.Matrix] = frame["boundary_ns_coordinates"]
    if not (
        _exact_inertia(ns_gram) == (1, 14, 0)
        and int(ns_gram.det()) == 1
        and any(int(ns_gram[index, index]) % 2 for index in range(15))
        and all(_primitive_integer(boundary) for boundary in boundaries_ns)
        and sp.Matrix.hstack(*boundaries_ns).T
        * ns_gram
        * sp.Matrix.hstack(*boundaries_ns)
        == -4 * sp.eye(5)
        and sum(boundaries_ns, sp.zeros(15, 1))
        == -2 * frame["canonical_ns_coordinates"]
    ):
        raise ArithmeticError("the reconstructed Picard lattice failed C0")

    # Adding every B_i/2 to N must recover M exactly.  The form on that
    # extension is the input even form, not the Picard form.
    extension_generators = sp.Matrix.hstack(ns_basis, *half_boundaries)
    extension_basis = hermite_normal_form(extension_generators)
    if extension_basis.shape != (15, 15) or abs(int(extension_basis.det())) != 1:
        raise ArithmeticError("the boundary halves do not recover the raw arithmetic lattice")
    smith = _smith_diagonal(main_gram)
    dual = main_gram.inv()
    delta = int(any(value.q != 1 for value in dual.diagonal()))

    # Fold the physical walls into CM roots, without using the archived CM
    # labels.  The full CM lattice is computed independently as a Smith
    # kernel, so root-span saturation is a genuine test.
    boundary_rows = sp.Matrix.vstack(
        *(half.T * main_gram for half in half_boundaries)
    )
    cm_basis = _integer_kernel_basis(boundary_rows)
    if cm_basis.shape != (15, 10):
        raise ArithmeticError("the recovered CM lattice does not have rank ten")
    cm_gram = cm_basis.T * main_gram * cm_basis / 2
    if not (
        all(value.q == 1 for value in cm_gram)
        and all(int(cm_gram[index, index]) % 2 == 0 for index in range(10))
        and int(cm_gram.det()) == -1
        and _exact_inertia(cm_gram) == (1, 9, 0)
    ):
        raise ArithmeticError("the recovered CM lattice is not E10")

    cm_roots_main: list[sp.Matrix] = []
    cm_root_source_indices: list[int] = []
    cm_root_kinds: list[str] = []
    crossing_boundary_pairs: dict[int, tuple[int, int]] = {}
    for wall_index, record in frame["taxonomy"].items():
        coordinates = wall_coordinates[wall_index]
        if record["kind"] == "minus2_curve":
            root = coordinates
            kind = "minus2"
        elif record["kind"] == "crossing_minus1_curve":
            pair = tuple(
                index
                for index, value in enumerate(record["boundary_intersection_profile"])
                if value == 1
            )
            if len(pair) != 2:
                raise ArithmeticError("a crossing curve lost its two boundaries")
            root = 2 * coordinates + half_boundaries[pair[0]] + half_boundaries[pair[1]]
            crossing_boundary_pairs[wall_index] = (pair[0], pair[1])
            kind = "minus1_root"
        else:
            continue
        if any((half.T * main_gram * root)[0] for half in half_boundaries):
            raise ArithmeticError("a folded root is not boundary-orthogonal")
        cm_roots_main.append(root)
        cm_root_source_indices.append(wall_index)
        cm_root_kinds.append(kind)
    if Counter(cm_root_kinds) != Counter({"minus2": 10, "minus1_root": 10}):
        raise ArithmeticError("the recovered CM root marking changed")
    cm_root_coordinates = cm_basis.gauss_jordan_solve(
        sp.Matrix.hstack(*cm_roots_main)
    )[0]
    if any(value.q != 1 for value in cm_root_coordinates):
        raise ArithmeticError("a folded root left the saturated CM lattice")
    root_span = hermite_normal_form(cm_root_coordinates)
    if root_span.shape != (10, 10) or abs(int(root_span.det())) != 1:
        raise ArithmeticError("the twenty roots no longer span the whole CM lattice")
    cm_root_gram = cm_root_coordinates.T * cm_gram * cm_root_coordinates
    type_vi = _type_vi_match(cm_root_gram, cm_root_kinds)
    if type_vi["labelled_isometry_count"] != 120:
        raise ArithmeticError("the labelled Type-VI symmetry order changed")

    cm_inequalities = (cm_gram * cm_root_coordinates).T
    cm_cone = _cone_audit(cm_gram, cm_inequalities, executable, 9)
    if not (
        cm_cone["ray_count"] == 657
        and cm_cone["finite_vertex_count"] == 600
        and cm_cone["ideal_vertex_count"] == 57
        and cm_cone["hyperideal_vertex_count"] == 0
        and cm_cone["all_facets_irredundant"]
        and cm_cone["all_causal_rays_in_one_time_sheet"]
    ):
        raise ArithmeticError("the complete CM cone changed")

    raw_cusp_lifts = {
        tuple(_vector(frame["physical_wall_ns_coordinates"][index]))
        for index, record in frame["taxonomy"].items()
        if record["kind"] == "cusp_lift_minus1_curve"
    }
    reconstructed_cusp_lifts: dict[tuple[int, ...], tuple[int, int]] = {}
    lift_masks: Counter[str] = Counter()
    boundary_lift_counts: Counter[int] = Counter()
    ideal_rays = [
        raw
        for raw in cm_cone["rays"]
        if (sp.Matrix(raw).T * cm_gram * sp.Matrix(raw))[0] == 0
    ]
    for cusp_index, raw in enumerate(ideal_rays):
        cusp = cm_basis * sp.Matrix(raw)
        mask: list[str] = []
        for boundary_index, half in enumerate(half_boundaries):
            exceptional_main = cusp - half
            exceptional_ns = inverse_ns * exceptional_main
            integral = all(value.q == 1 for value in exceptional_ns)
            mask.append("1" if integral else "0")
            if not integral:
                continue
            if not (
                (exceptional_ns.T * ns_gram * exceptional_ns)[0] == -1
                and sorted(
                    [
                    int((exceptional_ns.T * ns_gram * boundary)[0])
                    for boundary in boundaries_ns
                    ]
                )
                == [0, 0, 0, 0, 2]
                and _primitive_integer(exceptional_ns)
            ):
                raise ArithmeticError("an integral cusp lift failed adjunction")
            key = tuple(_vector(exceptional_ns))
            if key in reconstructed_cusp_lifts:
                raise ArithmeticError("two cusp-boundary pairs produced one exceptional class")
            reconstructed_cusp_lifts[key] = (cusp_index, boundary_index)
            boundary_lift_counts[boundary_index] += 1
        lift_masks["".join(mask)] += 1
    if not (
        len(reconstructed_cusp_lifts) == 35
        and set(reconstructed_cusp_lifts) == raw_cusp_lifts
        and sorted(boundary_lift_counts.values()) == [7, 7, 7, 7, 7]
        and lift_masks["00000"] == 22
        and sum(count for mask, count in lift_masks.items() if mask.count("1") == 1) == 35
        and all(mask.count("1") <= 1 for mask in lift_masks)
    ):
        raise ArithmeticError("the complete cusp-to-exceptional lift test changed")

    # Reconstruct every arithmetic mirror from the decoded divisor package.
    recovered_main: list[sp.Matrix] = [*half_boundaries]
    for wall_index, record in frame["taxonomy"].items():
        if record["kind"] == "minus2_curve":
            recovered_main.append(wall_coordinates[wall_index])
        elif record["kind"] == "crossing_minus1_curve":
            root_position = cm_root_source_indices.index(wall_index)
            root = cm_roots_main[root_position]
            left, right = crossing_boundary_pairs[wall_index]
            recovered_main.append(
                (root - half_boundaries[left] - half_boundaries[right]) / 2
            )
    recovered_main.extend(
        ns_basis * sp.Matrix(key) for key in sorted(reconstructed_cusp_lifts)
    )
    if any(any(value.q != 1 for value in vector) for vector in recovered_main):
        raise ArithmeticError("same-P reconstruction produced a nonintegral mirror")
    raw_set = {tuple(_vector(vector)) for vector in wall_coordinates}
    recovered_set = {tuple(_vector(vector)) for vector in recovered_main}
    if len(recovered_main) != 60 or recovered_set != raw_set:
        raise ArithmeticError("the decoded package did not return the same sixty mirrors")

    full_inequalities = (main_gram * sp.Matrix.hstack(*wall_coordinates)).T
    full_cone = _cone_audit(main_gram, full_inequalities, executable, 14)
    if not (
        full_cone["ray_count"] == 28_872
        and full_cone["finite_vertex_count"] == 27_845
        and full_cone["ideal_vertex_count"] == 1_027
        and full_cone["hyperideal_vertex_count"] == 0
        and full_cone["all_facets_irredundant"]
        and full_cone["all_causal_rays_in_one_time_sheet"]
    ):
        raise ArithmeticError("the full returned H14 cone changed")

    physical_matrix = sp.Matrix.hstack(*frame["physical_wall_ns_coordinates"])
    physical_gram = physical_matrix.T * ns_gram * physical_matrix
    taxonomy_payload = {
        kind: sorted(index + 1 for index, record in frame["taxonomy"].items() if record["kind"] == kind)
        for kind in expected_taxonomy
    }
    result = {
        "raw_arithmetic_lattice": {
            "rank": 15,
            "signature": [1, 14],
            "determinant": int(main_gram.det()),
            "smith_diagonal": smith,
            "two_elementary_invariants": {"r": 15, "a": 7, "delta": delta},
            "primitive_root_type_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(arithmetic_types.items())
            },
            "fifteen_divisibility_one_short_roots_form_basis": True,
            "all_wall_reflections_integral": True,
            "coxeter_gate": coxeter,
        },
        "exhaustive_boundary_frame_search": {
            "search_denominator": "all five-subsets of the 15 short divisibility-one basis roots",
            "prefix_survivor_counts": frame_counts,
            "unique_unordered_frame": True,
            "recovered_boundary_wall_indices_one_based": [
                index + 1 for index in frame["frame_wall_indices"]
            ],
            "recovered_boundary_half_vectors_in_main_root_basis": [
                _vector(half) for half in half_boundaries
            ],
            "surface_labels_consulted": False,
        },
        "reconstructed_Picard_lattice": {
            "identification": "odd unimodular I_(1,14)",
            "basis_in_main_root_coordinates": _rows(ns_basis),
            "gram": _rows(ns_gram),
            "determinant": int(ns_gram.det()),
            "signature": list(_exact_inertia(ns_gram)),
            "canonical_vector": _vector(frame["canonical_ns_coordinates"]),
            "canonical_square": int(
                (frame["canonical_ns_coordinates"].T * ns_gram * frame["canonical_ns_coordinates"])[0]
            ),
            "canonical_is_characteristic": True,
            "boundary_vectors": [_vector(boundary) for boundary in boundaries_ns],
            "minus_two_K_equals_boundary_sum": True,
            "boundary_frame_pairwise_orthogonal_square_minus_four": True,
            "Picard_index_in_extended_module": abs(int(ns_basis.det())),
        },
        "wall_taxonomy": {
            "counts": dict(sorted(frame["taxonomy_counts"].items())),
            "indices_one_based": taxonomy_payload,
            "adjunction_identity_K_dot_C_equals_minus_two_minus_C_square": True,
            "physical_negative_wall_gram_sha256": _payload_sha256(_rows(physical_gram)),
        },
        "extended_Picard_return": {
            "module": "M=<N,B_1/2,...,B_5/2>",
            "index_M_over_N": abs(int(ns_basis.det())),
            "doubled_form_is_exactly_the_raw_arithmetic_form": True,
            "boundary_halves_generate_M_over_N": True,
            "two_elementary_identification": "(15,7,1)",
        },
        "Coble_Mukai_reconstruction": {
            "rank": 10,
            "signature": list(_exact_inertia(cm_gram)),
            "even_unimodular": True,
            "identification": "E10=U+E8(-1)",
            "basis_in_main_root_coordinates": _rows(cm_basis),
            "gram": _rows(cm_gram),
            "roots_span_whole_CM": True,
            "root_count": 20,
            "root_marking_counts": dict(sorted(Counter(cm_root_kinds).items())),
            "root_gram": _rows(cm_root_gram),
            "Type_VI_gate": type_vi,
            "exact_cone": {key: value for key, value in cm_cone.items() if key != "rays"},
        },
        "cusp_lift_exhaustion": {
            "primitive_CM_cusps": len(ideal_rays),
            "lift_mask_histogram": dict(sorted(lift_masks.items())),
            "integral_cusp_boundary_pairs": len(reconstructed_cusp_lifts),
            "integral_lifts_per_boundary": {
                str(index + 1): boundary_lift_counts[index] for index in range(5)
            },
            "cusps_with_no_integral_lift": lift_masks["00000"],
            "all_and_only_integral_lifts_are_raw_minus_one_walls": True,
        },
        "exact_same_P_return": {
            "reconstructed_mirror_count": len(recovered_main),
            "raw_mirror_count": len(raw_set),
            "oriented_mirror_sets_equal": recovered_set == raw_set,
            "canonical_oriented_mirror_set_sha256": _payload_sha256(sorted(raw_set)),
            "full_H14_cone": {key: value for key, value in full_cone.items() if key != "rays"},
            "finite_volume_noncompact": True,
        },
    }
    return result


def build_certificate(
    raw_path: Path = RAW,
    cdd_executable: str | None = None,
) -> dict[str, Any]:
    executable = cdd_executable or shutil.which("cddexec_gmp")
    if executable is None:
        raise RuntimeError("cddexec_gmp is required for the exact cone replays")
    raw_hash = _sha256(raw_path)
    if raw_path == RAW and raw_hash != RAW_SHA256:
        raise ArithmeticError(
            f"the frozen Type-VI input changed: expected {RAW_SHA256}, got {raw_hash}"
        )
    archived = json.loads(raw_path.read_text(encoding="utf-8"))
    if archived.get("schema") != RAW_SCHEMA:
        raise ArithmeticError("the Type-VI input schema changed")

    # These are deliberately the only mathematical discovery fields read.
    arithmetic_gram = sp.Matrix(archived["arithmetic"]["gram"])
    walls = [sp.Matrix(row) for row in archived["arithmetic"]["wall_vectors"]]
    recognition = _recognize(arithmetic_gram, walls, executable)
    all_checks = {
        "raw_input_is_only_lattice_and_oriented_mirrors": True,
        "boundary_frame_unique_by_exhaustive_search": recognition[
            "exhaustive_boundary_frame_search"
        ]["unique_unordered_frame"],
        "odd_unimodular_NS_reconstructed": recognition["reconstructed_Picard_lattice"][
            "determinant"
        ]
        == 1,
        "twenty_roots_are_labelled_Type_VI": recognition["Coble_Mukai_reconstruction"][
            "Type_VI_gate"
        ]["labelled_gram_isomorphic_to_Kondo_Gamma_VI"],
        "all_cusp_lifts_exhausted": recognition["cusp_lift_exhaustion"][
            "all_and_only_integral_lifts_are_raw_minus_one_walls"
        ],
        "exact_same_P_return": recognition["exact_same_P_return"][
            "oriented_mirror_sets_equal"
        ],
        "finite_volume_noncompact": recognition["exact_same_P_return"][
            "finite_volume_noncompact"
        ],
    }
    if not all(all_checks.values()):
        raise ArithmeticError(f"the inverse recognition checks failed: {all_checks}")
    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "direction": "POLYTOPE_TO_SURFACE_inverse_recognizer",
        "input": {
            "path": str(raw_path.relative_to(ROOT)) if raw_path.is_relative_to(ROOT) else str(raw_path),
            "sha256": raw_hash,
            "schema": archived["schema"],
            "discovery_fields_used": ["arithmetic.gram", "arithmetic.wall_vectors"],
            "surface_fields_not_used": [
                "facets.labels",
                "arithmetic.boundary_frame_fold",
                "Coble_Mukai_chamber",
                "curve_exhaustion",
                "sources.kondo",
            ],
        },
        **recognition,
        "algebro_geometric_gate": {
            "first_point_where_classification_is_used": (
                "after unique NS/boundary/root/cusp reconstruction and exact same-P return"
            ),
            "Kondo_row": "Type VI, characteristic 3, five boundaries, twenty roots",
            "characteristic": 3,
            "surface_isomorphism_classes": 1,
            "automorphism_group": "S5",
            "R_invariant": "(E6 + D5, Z/2Z)",
            "existence_and_uniqueness_source": "Kondo, arXiv:2102.00591, Main Theorem and Type-VI example/uniqueness",
            "complete_nef_reason": (
                "Kondo-Vinberg exhausts the twenty effective irreducible CM roots; "
                "adjunction leaves only boundary, (-2), and (-1) curves; the Halphen/cusp "
                "dictionary and the exhaustive integral lift test recover every (-1) wall"
            ),
            "full_nef_carrier_exists": True,
        },
        "scope": {
            "proved": (
                "necessary-and-sufficient recognition of the unique Type-VI five-boundary "
                "package inside the stated raw (15,7,1), 60-wall arithmetic datum"
            ),
            "not_proved": [
                "automatic Coble realization for every numerical package passing C0-C2",
                "a theorem about all arithmetic Coxeter polytopes",
                "nonexistence of the abstract chamber when another carrier species is allowed",
                "a new abstract H14 polytope (the chamber is Alexeev's known (15,7,1) chamber)",
            ],
            "nine_class_status": (
                "this is one calibration for the finite nine-template decision compiler; "
                "template-free essential surjectivity remains a separate missing lemma"
            ),
        },
        "checks": all_checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=RAW)
    parser.add_argument("--cdd", default=None)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate(arguments.input, arguments.cdd)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(
        json.dumps(
            {
                "boundary_frame_prefix_counts": result[
                    "exhaustive_boundary_frame_search"
                ]["prefix_survivor_counts"],
                "wall_taxonomy": result["wall_taxonomy"]["counts"],
                "same_P": result["exact_same_P_return"]["oriented_mirror_sets_equal"],
                "carrier": result["algebro_geometric_gate"]["full_nef_carrier_exists"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
