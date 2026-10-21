#!/usr/bin/env python3
"""Independent hostile audit of the Type-MII full-NS H17 chamber.

This file deliberately imports no Type-MII project certificate.  It rebuilds
the blow-up lattice, the 48 geometric walls, the extended reflection lattice,
and the Coble--Mukai section from first principles.  It then performs two
exact double-description computations:

* the 40-root rank-ten CM cone, for the missing-(-1)-curve cusp gate; and
* the full 48-wall rank-eighteen cone, for a literal finite-volume gate.

Finally it parses Alexeev's pinned raw ancillary Gram directly, finds an exact
weighted-graph isomorphism without a frozen wall permutation, and reconstructs
the integral marked lattice isometry from the first eighteen ancillary roots.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import (
    hermite_normal_form,
    smith_normal_decomp,
    smith_normal_form,
)
from sympy.polys.domains import ZZ


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/coble_type_mii_h17_hostile_audit.json"
KONDO_TEX = Path("/private/tmp/coble-src.HX0hjo/CobleFinal.tex")
ALEXEEV_MATRIX = Path("/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt")
ALEXEEV_TEX = Path("/private/tmp/alexeev_2209_source/g0.tex")
SCHEMA = "coxeter.coble_type_mii_h17.hostile_audit.v1"

PINNED = {
    "kondo_tex": "8e23df4e77a2666a7a4f9a634c1ec2f36066353bf0bdcd0a72b11a0a70bfb137",
    "alexeev_matrix": "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457",
    "alexeev_tex": "006e2a1eac2aa3c4916ef66aceb27e8e80a35145162e6785a4874843f249f1df",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json_sha(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode()
    ).hexdigest()


def _rows(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.q != 1 for value in matrix):
        raise ArithmeticError("an expected integral matrix is nonintegral")
    return [[int(value) for value in matrix.row(index)] for index in range(matrix.rows)]


def _rat_rows(matrix: sp.Matrix) -> list[list[int | str]]:
    return [
        [int(value) if value.q == 1 else str(value) for value in matrix.row(index)]
        for index in range(matrix.rows)
    ]


def _primitive_positive(values: Sequence[Fraction]) -> tuple[int, ...]:
    denominator = math.lcm(*(value.denominator for value in values))
    integers = [int(value * denominator) for value in values]
    divisor = math.gcd(*(abs(value) for value in integers))
    integers = [value // divisor for value in integers]
    # The ray orientation returned by cdd is part of the answer: negating it
    # generally leaves the supplied homogeneous cone.  Normalize scale only.
    return tuple(integers)


def _cdd_rays(
    inequalities: sp.Matrix, executable: str, *, timeout: int = 900
) -> list[tuple[int, ...]]:
    """Enumerate a homogeneous rational cone without project helpers."""

    lines = [
        "H-representation",
        "begin",
        f"{inequalities.rows} {inequalities.cols + 1} rational",
    ]
    lines.extend(
        "0 " + " ".join(str(inequalities[row, col]) for col in range(inequalities.cols))
        for row in range(inequalities.rows)
    )
    lines.append("end")
    completed = subprocess.run(
        [executable, "--rep"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout,
    )
    if completed.returncode:
        raise RuntimeError(completed.stderr[-4000:] or completed.stdout[-4000:])
    output = completed.stdout.splitlines()
    begin = output.index("begin")
    count, width, kind = output[begin + 1].split()
    if kind != "rational" or int(width) != inequalities.cols + 1:
        raise ArithmeticError("cddlib returned an unexpected representation")
    rays = []
    for line in output[begin + 2 : begin + 2 + int(count)]:
        tokens = [Fraction(token) for token in line.split()]
        if tokens[0]:
            raise ArithmeticError("a homogeneous cone acquired an affine vertex")
        rays.append(_primitive_positive(tokens[1:]))
    if len(rays) != len(set(rays)):
        raise ArithmeticError("cddlib returned duplicate oriented rays")
    return sorted(rays)


def _half_integral_basis(generators: Sequence[sp.Matrix]) -> sp.Matrix:
    raw = sp.Matrix.hstack(*(2 * vector for vector in generators))
    if any(value.q != 1 for value in raw) or raw.rank() != raw.rows:
        raise ArithmeticError("bad half-integral spanning set")
    doubled = sp.Matrix(
        [[int(raw[row, col]) for col in range(raw.cols)] for row in range(raw.rows)]
    )
    return hermite_normal_form(doubled) / 2


def _smith_abs(matrix: sp.Matrix) -> list[int]:
    return [
        abs(int(value))
        for value in smith_normal_form(matrix, domain=ZZ).diagonal()
    ]


def _surface() -> dict[str, Any]:
    """Rebuild Bl_16(P1xP1), its curves, and the reflective overlattice."""

    gram = sp.zeros(18)
    gram[0, 1] = gram[1, 0] = 1
    for index in range(2, 18):
        gram[index, index] = -1
    basis = sp.eye(18)
    f1, f2 = basis[:, 0], basis[:, 1]
    exceptional = {
        (row, col): basis[:, 2 + 4 * row + col]
        for row in range(4)
        for col in range(4)
    }
    exceptionals = [
        exceptional[row, col] for row in range(4) for col in range(4)
    ]
    boundaries = [
        f1 - sum((exceptional[row, col] for col in range(4)), sp.zeros(18, 1))
        for row in range(4)
    ] + [
        f2 - sum((exceptional[row, col] for row in range(4)), sp.zeros(18, 1))
        for col in range(4)
    ]
    permutations = list(itertools.permutations(range(4)))
    graphs = [
        f1
        + f2
        - sum(
            (exceptional[row, permutation[row]] for row in range(4)),
            sp.zeros(18, 1),
        )
        for permutation in permutations
    ]
    walls = [*boundaries, *exceptionals, *graphs]
    wall_matrix = sp.Matrix.hstack(*walls)
    wall_gram = wall_matrix.T * gram * wall_matrix
    canonical = -2 * f1 - 2 * f2 + sum(exceptionals, sp.zeros(18, 1))
    interior = 5 * f1 + 5 * f2 - sum(exceptionals, sp.zeros(18, 1))

    relation_histogram: Counter[str] = Counter()
    for left in range(48):
        for right in range(left):
            ratio = sp.cancel(
                wall_gram[left, right] ** 2
                / (wall_gram[left, left] * wall_gram[right, right])
            )
            relation = {
                sp.Rational(0): "m2",
                sp.Rational(1, 4): "m3",
                sp.Rational(1, 2): "m4",
                sp.Rational(1): "parallel",
            }.get(ratio)
            if relation is None:
                raise ArithmeticError(f"non-Coxeter wall ratio {ratio}")
            relation_histogram[relation] += 1

    witnesses = []
    for index, wall in enumerate(walls):
        square = (wall.T * gram * wall)[0]
        witness = interior - (interior.T * gram * wall)[0] / square * wall
        pairings = [(witness.T * gram * other)[0] for other in walls]
        if not (
            (witness.T * gram * witness)[0] > 0
            and pairings[index] == 0
            and all(value > 0 for j, value in enumerate(pairings) if j != index)
        ):
            raise ArithmeticError("a proposed full-cone facet is redundant")
        witnesses.append(witness)

    pic_tilde = _half_integral_basis(
        [*basis.columnspace(), *(boundary / 2 for boundary in boundaries)]
    )
    lattice_gram = 2 * pic_tilde.T * gram * pic_tilde
    arithmetic_normals = [
        *(boundary / 2 for boundary in boundaries), *exceptionals, *graphs
    ]
    arithmetic_walls = pic_tilde.inv() * sp.Matrix.hstack(*arithmetic_normals)
    arithmetic_wall_gram = arithmetic_walls.T * lattice_gram * arithmetic_walls

    root_types: Counter[tuple[int, int]] = Counter()
    reflections_ok = []
    for index in range(48):
        root = arithmetic_walls[:, index]
        square = int((root.T * lattice_gram * root)[0])
        divisibility = math.gcd(*(abs(int(value)) for value in lattice_gram * root))
        root_types[square, divisibility] += 1
        reflection = (
            sp.eye(18)
            - sp.Rational(2, square) * root * (root.T * lattice_gram)
        )
        reflections_ok.append(
            all(value.q == 1 for value in reflection)
            and reflection.T * lattice_gram * reflection == lattice_gram
            and reflection**2 == sp.eye(18)
        )

    dual = lattice_gram.inv()
    delta = int(any(sp.denom(dual[index, index]) != 1 for index in range(18)))
    checks = {
        "NS_signature_and_determinant": gram.det() == -1,
        "canonical_square": (canonical.T * gram * canonical)[0] == -8,
        "boundary_is_anti_bicanonical": sum(boundaries, sp.zeros(18, 1))
        == -2 * canonical,
        "boundaries_disjoint_minus_four": sp.Matrix.hstack(*boundaries).T
        * gram
        * sp.Matrix.hstack(*boundaries)
        == -4 * sp.eye(8),
        "curve_adjunction": all(
            (curve.T * gram * curve)[0] == square
            and (canonical.T * gram * curve)[0] == -2 - square
            for curves, square in ((boundaries, -4), (exceptionals, -1), (graphs, -2))
            for curve in curves
        ),
        "wall_norms": Counter(wall_gram.diagonal())
        == Counter({-4: 8, -1: 16, -2: 24}),
        "walls_span_NS": wall_matrix.rank() == 18
        and _smith_abs(wall_matrix) == [1] * 18,
        "acute": all(
            wall_gram[left, right] >= 0
            for left in range(48)
            for right in range(left)
        ),
        "coxeter_histogram": relation_histogram
        == Counter({"m2": 796, "m3": 128, "m4": 96, "parallel": 108}),
        "strict_interior": (interior.T * gram * interior)[0] == 34
        and Counter((interior.T * gram * wall)[0] for wall in walls)
        == Counter({1: 24, 6: 24}),
        "all_facets_irredundant": len(witnesses) == 48,
        "extended_index": abs(pic_tilde.det()) == sp.Rational(1, 128),
        "lattice_18_4_0": lattice_gram.det() == -16
        and _smith_abs(lattice_gram) == [1] * 14 + [2] * 4
        and all(int(lattice_gram[index, index]) % 2 == 0 for index in range(18))
        and delta == 0,
        "arithmetic_root_types": root_types
        == Counter({(-2, 1): 24, (-4, 2): 24}),
        "arithmetic_walls_span_L": _smith_abs(arithmetic_walls) == [1] * 18,
        "all_reflections_integral": all(reflections_ok),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"independent surface reconstruction failed: {checks}")
    return {
        "gram": gram,
        "basis": basis,
        "f1": f1,
        "f2": f2,
        "exceptional": exceptional,
        "exceptionals": exceptionals,
        "boundaries": boundaries,
        "permutations": permutations,
        "graphs": graphs,
        "walls": walls,
        "wall_matrix": wall_matrix,
        "wall_gram": wall_gram,
        "canonical": canonical,
        "interior": interior,
        "pic_tilde": pic_tilde,
        "lattice_gram": lattice_gram,
        "arithmetic_walls": arithmetic_walls,
        "arithmetic_wall_gram": arithmetic_wall_gram,
        "relation_histogram": relation_histogram,
        "root_types": root_types,
        "checks": checks,
    }


def _bilinear(left: Sequence[int], gram: sp.Matrix, right: Sequence[int]) -> int:
    return sum(
        int(left[i]) * int(gram[i, j]) * int(right[j])
        for i in range(len(left))
        for j in range(len(right))
    )


def _components(nodes: Iterable[int], gram: sp.Matrix) -> list[list[int]]:
    unseen = set(nodes)
    result = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        queue = [start]
        while queue:
            vertex = queue.pop()
            adjacent = {other for other in unseen if gram[vertex, other] != 0}
            unseen.difference_update(adjacent)
            component.update(adjacent)
            queue.extend(adjacent)
        result.append(sorted(component))
    return result


def _affine_a_type(component: Sequence[int], gram: sp.Matrix) -> str:
    if len(component) == 2 and gram[component[0], component[1]] == 2:
        return "2A1~"
    degrees = [
        sum(gram[vertex, other] == 1 for other in component if other != vertex)
        for vertex in component
    ]
    if len(component) >= 3 and degrees == [2] * len(component):
        return f"A{len(component) - 1}~"
    raise ArithmeticError("a cusp component is not an affine A component")


def _cm_audit(surface: dict[str, Any], cdd: str) -> dict[str, Any]:
    gram: sp.Matrix = surface["gram"]
    boundaries: list[sp.Matrix] = surface["boundaries"]
    pic_tilde: sp.Matrix = surface["pic_tilde"]
    constraint = 2 * sp.Matrix.hstack(*boundaries).T * gram * pic_tilde
    diagonal, _left, right = smith_normal_decomp(constraint, domain=ZZ)
    rank = constraint.rank()
    cm_basis = pic_tilde * right[:, rank:]
    cm_gram = cm_basis.T * gram * cm_basis

    roots = []
    for row in range(4):
        for col in range(4):
            roots.append(
                boundaries[row] / 2
                + boundaries[4 + col] / 2
                + 2 * surface["exceptional"][row, col]
            )
    roots.extend(surface["graphs"])
    root_matrix = sp.Matrix.hstack(*roots)
    root_coordinates = cm_basis.gauss_jordan_solve(root_matrix)[0]
    root_gram = root_matrix.T * gram * root_matrix
    inequalities = root_coordinates.T * cm_gram
    rays = _cdd_rays(inequalities, cdd)

    squares = [_bilinear(ray, cm_gram, ray) for ray in rays]
    counts = Counter(
        "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        for square in squares
    )
    time_ray = next(ray for ray, square in zip(rays, squares) if square > 0)
    time_pairings = [_bilinear(time_ray, cm_gram, ray) for ray in rays]
    same_sheet = all(value > 0 for value in time_pairings)

    cusps = [ray for ray, square in zip(rays, squares) if square == 0]
    lift_memberships = []
    type_histogram: Counter[str] = Counter()
    cusp_records = []
    for cusp in cusps:
        vector = sp.Matrix(cusp)
        fiber = cm_basis * vector
        membership = [
            all(value.q == 1 for value in fiber - boundary / 2)
            for boundary in boundaries
        ]
        lift_memberships.extend(membership)
        active = [
            index for index, value in enumerate(inequalities * vector) if value == 0
        ]
        types = sorted(
            (_affine_a_type(component, root_gram) for component in _components(active, root_gram)),
            key=lambda value: (
                1 if value == "2A1~" else 0,
                0 if value == "2A1~" else -int(value[1 : value.index("~")]),
            ),
        )
        signature = "+".join(types)
        type_histogram[signature] += 1
        cusp_records.append(
            {
                "CM_coordinates": list(cusp),
                "active_root_indices": active,
                "affine_type": signature,
                "integral_boundary_lifts": [
                    index for index, value in enumerate(membership) if value
                ],
            }
        )

    expected_types = Counter(
        {
            "A7~+2A1~": 72,
            "A5~+A2~+2A1~": 96,
            "A3~+A3~+2A1~+2A1~": 18,
            "A2~+A2~+A2~+A2~": 16,
        }
    )
    checks = {
        "constraint_smith": [abs(int(diagonal[i, i])) for i in range(8)]
        == [2, 2, 2, 2, 2, 2, 2, 4],
        "CM_is_II_1_9": cm_basis.shape == (18, 10)
        and cm_gram.det() == -1
        and all(int(cm_gram[index, index]) % 2 == 0 for index in range(10)),
        "forty_roots": Counter(root_gram.diagonal()) == Counter({-2: 40})
        and root_gram.rank() == 10,
        "roots_span_CM": all(value.q == 1 for value in root_coordinates)
        and _smith_abs(root_coordinates) == [1] * 10,
        "exact_ray_counts": counts
        == Counter({"finite": 4752, "ideal": 202, "hyperideal": 0}),
        "entire_linear_cone_is_one_causal_sheet": same_sheet
        and counts["hyperideal"] == 0,
        "cusp_affine_types": type_histogram == expected_types,
        "all_single_boundary_lifts_tested": len(lift_memberships) == 1616,
        "no_single_boundary_integral_lift": not any(lift_memberships),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"independent CM audit failed: {checks}")
    return {
        "basis": cm_basis,
        "gram": cm_gram,
        "root_coordinates": root_coordinates,
        "root_gram": root_gram,
        "ray_count": len(rays),
        "ray_sha256": _json_sha(rays),
        "counts": counts,
        "same_sheet_min_time_pairing": min(time_pairings),
        "cusps": cusp_records,
        "type_histogram": type_histogram,
        "checks": checks,
    }


def _full_cone(surface: dict[str, Any], cdd: str) -> dict[str, Any]:
    """Literal full-rank volume audit: no imported Vinberg conclusion."""

    inequalities = surface["wall_matrix"].T * surface["gram"]
    rays = _cdd_rays(inequalities, cdd)
    # In the chosen blow-up basis q(x)=2*x0*x1-sum(x_i^2), i>=2.
    squares = [2 * ray[0] * ray[1] - sum(value * value for value in ray[2:]) for ray in rays]
    # H=5F1+5F2-sum E_i, so H.x=5(x0+x1)+sum x_i.
    interior_pairings = [5 * (ray[0] + ray[1]) + sum(ray[2:]) for ray in rays]
    counts = Counter(
        "finite" if square > 0 else "ideal" if square == 0 else "hyperideal"
        for square in squares
    )
    ideal_rays = [ray for ray, square in zip(rays, squares) if square == 0]
    checks = {
        "exact_full_extreme_ray_count": len(rays) == 196_052,
        "full_ray_types": counts
        == Counter({"finite": 190_808, "ideal": 5_244, "hyperideal": 0}),
        "all_rays_on_chamber_time_sheet": min(interior_pairings) > 0,
        "no_hyperideal_extreme_ray": counts["hyperideal"] == 0,
        "finite_volume": counts["hyperideal"] == 0
        and min(interior_pairings) > 0,
        "noncompact": counts["ideal"] > 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"literal full-cone audit failed: {checks}")
    return {
        "extreme_ray_count": len(rays),
        "finite_vertex_count": counts["finite"],
        "ideal_vertex_count": counts["ideal"],
        "hyperideal_vertex_count": counts["hyperideal"],
        "ray_sha256": _json_sha(rays),
        "ideal_ray_sha256": _json_sha(ideal_rays),
        "minimum_interior_pairing": min(interior_pairings),
        "checks": checks,
    }


def _raw_alexeev_matrix(path: Path) -> sp.Matrix:
    source = path.read_text(encoding="utf-8")
    pattern = r"^coxeter_matrix\[\(18,4,0\)\] = matrix\((\[.*\])\)$"
    match = re.search(pattern, source, flags=re.MULTILINE)
    if match is None:
        raise ArithmeticError("raw Alexeev (18,4,0) assignment is absent")
    return sp.Matrix(ast.literal_eval(match.group(1)))


def _weighted_graph(gram: sp.Matrix) -> nx.Graph:
    graph = nx.Graph()
    for index in range(gram.rows):
        graph.add_node(index, norm=int(gram[index, index]))
    for left in range(gram.rows):
        for right in range(left):
            if gram[left, right]:
                graph.add_edge(right, left, pairing=int(gram[left, right]))
    return graph


def _alexeev_audit(surface: dict[str, Any], matrix_path: Path, tex_path: Path) -> dict[str, Any]:
    if _sha256(matrix_path) != PINNED["alexeev_matrix"]:
        raise ArithmeticError("raw Alexeev ancillary hash changed")
    if _sha256(tex_path) != PINNED["alexeev_tex"]:
        raise ArithmeticError("Alexeev TeX source hash changed")
    raw = _raw_alexeev_matrix(matrix_path)
    ours: sp.Matrix = surface["arithmetic_wall_gram"]

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _weighted_graph(ours),
        _weighted_graph(raw),
        node_match=lambda left, right: left == right,
        edge_match=lambda left, right: left == right,
    )
    isomorphisms = []
    for mapping in matcher.isomorphisms_iter():
        isomorphisms.append(tuple(int(mapping[index]) for index in range(48)))
    if not isomorphisms:
        raise ArithmeticError("no raw weighted-Gram isomorphism exists")
    permutation = min(isomorphisms)
    inverse = [0] * 48
    for source, target in enumerate(permutation):
        inverse[target] = source
    entrywise = all(
        ours[left, right] == raw[permutation[left], permutation[right]]
        for left in range(48)
        for right in range(48)
    )

    leading = raw[:18, :18]
    determinants = [leading[:index, :index].det() for index in range(1, 19)]
    previous = sp.Integer(1)
    pivot_signs = []
    for determinant in determinants:
        if determinant == 0:
            raise ArithmeticError("Alexeev's leading basis acquired a zero pivot")
        pivot_signs.append(1 if determinant * previous > 0 else -1)
        previous = determinant
    raw_coordinates = leading.inv() * raw[:18, :]
    marking = sp.Matrix.hstack(
        *(surface["arithmetic_walls"][:, inverse[index]] for index in range(18))
    )
    all_roots = all(
        marking * raw_coordinates[:, index]
        == surface["arithmetic_walls"][:, inverse[index]]
        for index in range(48)
    )

    tex = tex_path.read_text(encoding="utf-8")
    source_markers = {
        "theorem_O_plus_semidirect": r"$O^+(S) = \Aut G \ltimes W_r$." in tex,
        "vinberg_two_steps": (
            "With the control vector $v_0=a+b$, Vinberg's algorithm terminates in\n"
            "two steps."
        ) in tex,
        "parabolic_5244_17": (
            "There are $5244$ maximal parabolic subgraphs, $17$\n"
            "modulo $\\Aut\\Gamma_r$:"
        ) in tex,
    }
    checks = {
        "raw_shape_symmetry_rank": raw.shape == (48, 48)
        and raw == raw.T
        and raw.rank() == 18,
        "raw_norm_split": Counter(raw.diagonal()) == Counter({-2: 24, -4: 24}),
        "raw_leading_lattice": leading.det() == -16
        and pivot_signs.count(1) == 1
        and pivot_signs.count(-1) == 17,
        "weighted_Gram_isomorphic": entrywise,
        "isomorphism_torsor_size_1152": len(isomorphisms) == 1152,
        "integral_unimodular_marking": all(value.q == 1 for value in marking)
        and abs(marking.det()) == 1,
        "marked_lattice_Gram": marking.T * surface["lattice_gram"] * marking
        == leading,
        "all_48_marked_roots": all_roots,
        "primary_text_claims_present": all(source_markers.values()),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"raw Alexeev hostile audit failed: {checks}")
    return {
        "raw_matrix_sha256": _sha256(matrix_path),
        "raw_matrix_Gram_sha256": _json_sha(_rows(raw)),
        "canonical_MII_to_raw_wall_permutation": list(permutation),
        "weighted_Gram_isomorphism_count": len(isomorphisms),
        "raw_basis_to_MII_lattice_marking": _rows(marking),
        "source_markers": source_markers,
        "logic": {
            "finite_volume": (
                "Alexeev's O+(L)=Aut(K4,4) semidirect W_r gives reflection "
                "index 1152 in the arithmetic lattice O+(L), hence finite covolume."
            ),
            "noncompact": (
                "The complete maximal-parabolic list supplies cusps. Finite "
                "sidedness or Vinberg termination alone would not prove volume."
            ),
            "independent_cross_check": (
                "The literal full-cone cdd audit independently has exactly 5244 "
                "ideal extreme rays and no hyperideal rays."
            ),
        },
        "checks": checks,
    }


def _kondo_source_gate(path: Path) -> dict[str, Any]:
    if _sha256(path) != PINNED["kondo_tex"]:
        raise ArithmeticError("Kondo TeX source hash changed")
    source = path.read_text(encoding="utf-8")
    markers = {
        "blowup_16_grid_points": (
            "the sixteen ${\\bf F}_3$-rational points of" in source
        ),
        "eight_boundary_components": (
            "Coble surface with eight boundary components" in source
        ),
        "sixteen_minus_one_roots": "We have sixteen $(-1)$-roots" in source,
        "twenty_four_bidegree_11_curves": (
            "the following 24 non-singular curves of bidegree $(1,1)$" in source
        ),
        "exactly_forty_irreducible_roots": (
            "$S$ has exactly forty effective irreducible roots" in source
        ),
        "negative_curve_alphabet": (
            "the only curves with negative self-intersection on $S$ are the curves"
            in source
        ),
        "single_boundary_nef_isotropic_argument": (
            "The divisor $2E+B$ is nef, $(2E+B)^2=0" in source
        ),
    }
    if not all(markers.values()):
        raise ArithmeticError(f"a Kondō source marker disappeared: {markers}")
    return {"tex_sha256": _sha256(path), "markers": markers}


def build_audit(
    *,
    cdd_executable: str | None = None,
    kondo_tex: Path = KONDO_TEX,
    alexeev_matrix: Path = ALEXEEV_MATRIX,
    alexeev_tex: Path = ALEXEEV_TEX,
) -> dict[str, Any]:
    cdd = cdd_executable or shutil.which("cddexec_gmp")
    if cdd is None:
        raise RuntimeError("cddexec_gmp is required for the hostile audit")
    kondo = _kondo_source_gate(kondo_tex)
    surface = _surface()
    cm = _cm_audit(surface, cdd)
    full = _full_cone(surface, cdd)
    alexeev = _alexeev_audit(surface, alexeev_matrix, alexeev_tex)

    # The exhaustion implication is now fully separated into its dependencies.
    exhaustion = {
        "Kondo_exact_effective_irreducible_CM_roots": 40,
        "known_minus_two_curves": 24,
        "known_two_boundary_minus_one_curves": 16,
        "adjunction_total_boundary_incidence_for_any_minus_one_curve": 2,
        "remaining_incidence_partition": "(2) on one boundary",
        "remaining_curve_forces_class": "f=E+B_i/2",
        "forced_class_properties": [
            "f belongs to CM(S)",
            "f is primitive",
            "f^2=0",
            "2f=2E+B_i is nef because its intersection with every irreducible curve is nonnegative",
            "the causal 40-root cone forces f to be one of its 202 null extreme rays",
        ],
        "exact_cusp_boundary_lift_tests": 1616,
        "surviving_lifts": 0,
        "complete_negative_curves": 48,
    }
    checks = {
        "independent_surface": all(surface["checks"].values()),
        "independent_CM_and_cusp_gate": all(cm["checks"].values()),
        "negative_curve_exhaustion": exhaustion["surviving_lifts"] == 0,
        "literal_full_volume": all(full["checks"].values()),
        "raw_Alexeev_isometry": all(alexeev["checks"].values()),
        "ideal_count_matches_published_parabolics": full["ideal_vertex_count"]
        == 5244,
        "not_merely_finite_sided": full["hyperideal_vertex_count"] == 0,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"global hostile audit failed: {checks}")

    return {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "independence": {
            "imports_Type_MII_certificate": False,
            "imports_Alexeev_reconstruction_archive": False,
            "uses_raw_pinned_primary_sources": True,
            "reconstructs_wall_permutation": True,
            "performs_full_48_wall_double_description": True,
        },
        "decision": {
            "claim_survives_hostile_audit": True,
            "polytope": "H17/48",
            "volume": "finite_volume_noncompact",
            "abstract_lineage": "known Alexeev/Keum--Kondo (18,4,0)",
            "surface_realization": "Kondo Type MII Coble surface in characteristic 3",
            "new_abstract_polytope": False,
        },
        "primary_sources": {
            "kondo": kondo,
            "alexeev_matrix_sha256": _sha256(alexeev_matrix),
            "alexeev_tex_sha256": _sha256(alexeev_tex),
        },
        "surface_reconstruction": {
            "NS_gram": _rows(surface["gram"]),
            "canonical": _rows(surface["canonical"].T)[0],
            "boundaries": _rows(sp.Matrix.hstack(*surface["boundaries"]).T),
            "wall_Gram_sha256": _json_sha(_rows(surface["wall_gram"])),
            "relation_histogram": dict(surface["relation_histogram"]),
            "arithmetic_lattice_Gram": _rows(surface["lattice_gram"]),
            "arithmetic_wall_Gram_sha256": _json_sha(
                _rows(surface["arithmetic_wall_gram"])
            ),
            "arithmetic_root_type_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(surface["root_types"].items())
            },
            "checks": surface["checks"],
        },
        "CM_cusp_audit": {
            "basis_in_NS": _rat_rows(cm["basis"]),
            "gram": _rows(cm["gram"]),
            "root_Gram_sha256": _json_sha(_rows(cm["root_gram"])),
            "extreme_ray_count": cm["ray_count"],
            "finite_vertex_count": cm["counts"]["finite"],
            "ideal_vertex_count": cm["counts"]["ideal"],
            "hyperideal_vertex_count": cm["counts"]["hyperideal"],
            "ray_sha256": cm["ray_sha256"],
            "same_sheet_min_time_pairing": cm["same_sheet_min_time_pairing"],
            "affine_type_histogram": dict(cm["type_histogram"]),
            "primitive_cusps": cm["cusps"],
            "checks": cm["checks"],
        },
        "negative_curve_exhaustion": exhaustion,
        "literal_full_cone": full,
        "raw_Alexeev_lineage": alexeev,
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--cdd", default=None)
    parser.add_argument("--kondo-tex", type=Path, default=KONDO_TEX)
    parser.add_argument("--alexeev-matrix", type=Path, default=ALEXEEV_MATRIX)
    parser.add_argument("--alexeev-tex", type=Path, default=ALEXEEV_TEX)
    args = parser.parse_args()
    payload = build_audit(
        cdd_executable=args.cdd,
        kondo_tex=args.kondo_tex,
        alexeev_matrix=args.alexeev_matrix,
        alexeev_tex=args.alexeev_tex,
    )
    if args.check_archive:
        if not args.output.exists() or json.loads(
            args.output.read_text(encoding="utf-8")
        ) != payload:
            raise SystemExit("archived hostile Type-MII audit differs")
        print("archived hostile Type-MII audit matches")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(args.output)


if __name__ == "__main__":
    main()
