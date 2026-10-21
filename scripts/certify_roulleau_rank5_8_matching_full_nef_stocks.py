#!/usr/bin/env python3
"""Freeze the rank-5--8 Roulleau stocks that meet open ADEG headers.

This is a deliberately header-scoped source compiler, not a classification of
all low-rank K3 nef chambers.  Roulleau's final table is filtered by the six
rank/root-count headers occurring among the previously header-empty stable
ADEG triality packages:

    (5,5), (6,6), (7,7), (8,8), (8,9), (8,10).

For every resulting stock, the published complete irreducible (-2)-curve
graph and a strict ample class are replayed exactly.  A small exhaustive
implementation of Vinberg's edge/cusp criterion independently certifies
finite volume (and records compact versus noncompact).  The output is a
machine denominator for the separate polytope-first promotion gate.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
from pathlib import Path
from typing import Any, Iterable

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/roulleau_rank5_8_matching_full_nef_stocks.json"
SCHEMA = "coxeter.roulleau_rank5_8_matching_full_nef_stocks.v1"
ROULLEAU_ARXIV = "https://arxiv.org/abs/2003.08985"
ROULLEAU_SOURCE_TAR_SHA256 = (
    "3e445ed01446b47ded6c007b395308de34e566b736f0ba982470bc497adeb1e9"
)
ROULLEAU_SOURCE_TEX_SHA256 = (
    "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3"
)
TARGET_HEADERS = {(5, 5), (6, 6), (7, 7), (8, 8), (8, 9), (8, 10)}


def _canonical_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _payload(matrix: sp.Matrix) -> list[list[int]]:
    if not all(value.is_Integer for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return [
        [int(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def _gram(
    count: int,
    ordinary_one_based: Iterable[tuple[int, int]],
    parallel_one_based: Iterable[tuple[int, int]] = (),
) -> sp.Matrix:
    result = -2 * sp.eye(count)
    for left, right in ordinary_one_based:
        result[left - 1, right - 1] = result[right - 1, left - 1] = 1
    for left, right in parallel_one_based:
        result[left - 1, right - 1] = result[right - 1, left - 1] = 2
    return result


def _case(
    stock_id: str,
    lattice: str,
    rank: int,
    count: int,
    ordinary: Iterable[tuple[int, int]],
    ample: tuple[int, ...],
    ample_square: int,
    ample_pairings: tuple[int, ...],
    tex_lines: str,
    branch_model: str,
    parallel: Iterable[tuple[int, int]] = (),
) -> dict[str, Any]:
    return {
        "id": stock_id,
        "lattice": lattice,
        "rank": rank,
        "gram": _gram(count, ordinary, parallel),
        "basis_indices": tuple(range(rank)),
        "ample": ample,
        "ample_square": ample_square,
        "ample_pairings": ample_pairings,
        "tex_lines": tex_lines,
        "branch_model": branch_model,
    }


def _cases() -> list[dict[str, Any]]:
    return [
        _case(
            "R5-UA3-NEF-H4-5", "U+A3(-1)", 5, 5,
            [(1, 2), (2, 3), (2, 4), (3, 5), (4, 5)],
            (5, 11, 9, 9, 8), 50, (1, 1, 1, 1, 2), "3257-3320",
            "double cover of F2 with an A3 branch singularity",
        ),
        _case(
            "R5-P4D4-NEF-H4-5", "[4]+D4(-1)", 5, 5,
            [(1, 5), (5, 3), (1, 2), (1, 4), (3, 2), (3, 4)],
            (4, 3, 3, 3, 3), 22, (1, 1, 3, 1, 1), "3321-3392",
            "double plane branched over a sextic with three collinear nodes",
        ),
        _case(
            "R6-UA4-NEF-H5-6", "U+A4(-1)", 6, 6,
            [(3, 5), (3, 2), (2, 1), (1, 6), (6, 5), (3, 4)],
            (8, 9, 11, 5, 9, 8), 50, (1,) * 6, "4030-4100",
            "double cover of F4 with an A4 branch singularity",
        ),
        _case(
            "R6-U2D4-NEF-H5-6", "U(2)+D4(-1)", 6, 6,
            [(6, 1), (6, 2), (6, 3), (6, 4), (6, 5)],
            (3, 3, 3, 3, 3, 7), 22, (1,) * 6, "4507-4580",
            "double plane branched over a transverse line and smooth quintic",
        ),
        _case(
            "R6-UD4-NEF-H5-6", "U+D4(-1)", 6, 6,
            [(6, 5), (5, 1), (1, 3), (1, 2), (1, 4)],
            (21, 10, 10, 10, 13, 6), 70, (1,) * 6, "4581-4660",
            "double cover of F4 with a D4 branch singularity",
        ),
        _case(
            "R7-UA5-NEF-H6-7", "U+A5(-1)", 7, 7,
            [(5, 4), (4, 3), (3, 2), (2, 7), (7, 6), (6, 5), (2, 1)],
            (7, 15, 12, 10, 9, 10, 12), 84,
            (1, 1, 1, 1, 2, 1, 1), "5216-5292",
            "double cover of F4 with an A5 branch singularity",
        ),
        _case(
            "R7-UD5-NEF-H6-7", "U+D5(-1)", 7, 7,
            [(1, 2), (2, 3), (3, 5), (5, 6), (3, 4), (5, 7)],
            (8, 17, 27, 13, 25, 12, 12), 114, (1,) * 7, "5293-5353",
            "double cover of F4 with a D5 branch singularity",
        ),
        _case(
            "R8-UD6-NEF-H7-8", "U+D6(-1)", 8, 8,
            [(1, 2), (2, 4), (4, 5), (5, 6), (6, 7), (4, 3), (6, 8)],
            (11, 24, 18, 38, 35, 33, 16, 16), 220,
            (2, 1, 2, 1, 1, 1, 1, 1), "5354-5413",
            "double cover of F4 with a D6 branch singularity",
        ),
        _case(
            "R8-UD4A1x2-NEF-H7-10", "U+D4(-1)+A1(-1)^2", 8, 10,
            [(1, 2), (2, 5), (5, 6), (2, 3), (2, 4), (6, 7), (6, 8)],
            (6, 13, 6, 6, 9, 6, 2, 2), 54,
            (1, 1, 1, 1, 1, 1, 2, 2, 4, 4), "5414-5555",
            "double plane branched over a line and a nodal quintic",
            parallel=[(7, 9), (8, 10)],
        ),
        _case(
            "R8-UA2x3-NEF-H7-10", "U+A2(-1)^3", 8, 10,
            [
                (1, 2), (1, 3), (1, 4),
                (2, 5), (2, 6), (5, 6),
                (3, 7), (3, 9), (7, 9),
                (4, 10), (4, 8), (10, 8),
            ],
            (3, 5, 1, 1, 4, 4, 0, 0), 18, (1,) * 10, "5795-5898",
            "double plane with a D4 branch singularity and splitting tangents",
        ),
        _case(
            "R8-UA3x2-NEF-H7-9", "U+A3(-1)^2", 8, 9,
            [
                (4, 5), (5, 6),
                (4, 2), (4, 3), (2, 1), (3, 1),
                (6, 7), (6, 8), (7, 9), (8, 9),
            ],
            (5, 6, 6, 8, 5, 3, 1, 1), 40,
            (2, 1, 1, 1, 1, 1, 1, 1, 2), "5899-5969",
            "double cover of F4 with two A3 branch singularities",
        ),
        _case(
            "R8-UA2A4-NEF-H7-9", "U+A2(-1)+A4(-1)", 8, 9,
            [
                (5, 1), (1, 2), (2, 3), (3, 4), (4, 5),
                (5, 6), (6, 7), (7, 8), (7, 9), (8, 9),
            ],
            (7, 6, 6, 7, 9, 5, 2, 0), 42,
            (1, 1, 1, 1, 1, 1, 1, 2, 2), "5970-6037",
            "double cover of F4 with A4 and A2 branch singularities",
        ),
        _case(
            "R8-UA1A5-NEF-H7-9", "U+A1(-1)+A5(-1)", 8, 9,
            [
                (7, 6), (6, 5), (5, 4), (4, 9), (9, 8), (8, 7),
                (4, 3), (3, 2),
            ],
            (10, 12, 7, 3, 0, -2, -3, -2), 76,
            (4, 3, 1, 1, 1, 1, 2, 1, 1), "6038-6105",
            "double cover of F4 with A5 and nodal branch singularities",
            parallel=[(2, 1)],
        ),
        _case(
            "R8-UA6-NEF-H7-8", "U+A6(-1)", 8, 8,
            [(2, 8), (8, 7), (7, 6), (6, 5), (5, 4), (4, 3), (3, 2), (2, 1)],
            (7, 15, 12, 10, 9, 9, 10, 12), 84, (1,) * 8, "6106-6182",
            "double cover of F4 with an A6 branch singularity",
        ),
        _case(
            "R8-UA2D4-NEF-H7-9", "U+A2(-1)+D4(-1)", 8, 9,
            [
                (2, 4), (4, 5), (5, 6), (6, 7), (4, 1), (4, 3),
                (7, 8), (7, 9), (8, 9),
            ],
            (7, 7, 7, 15, 10, 6, 3, 1), 56,
            (1, 1, 1, 1, 1, 1, 1, 1, 4), "6183-6247",
            "double cover of F4 with D4 and A2 branch singularities",
        ),
        _case(
            "R8-UA1D5-NEF-H7-9", "U+A1(-1)+D5(-1)", 8, 9,
            [
                (2, 3), (3, 4), (4, 5), (5, 7), (5, 6),
                (7, 8), (7, 9),
            ],
            (9, 12, 8, 5, 3, 1, 1, 0), 96,
            (6, 2, 1, 1, 1, 1, 1, 1, 1), "6248-6314",
            "double cover of F4 with D5 and nodal branch singularities",
            parallel=[(1, 2)],
        ),
        _case(
            "R8-UE6-NEF-H7-8", "U+E6(-1)", 8, 8,
            [(1, 2), (2, 3), (3, 4), (4, 5), (5, 7), (4, 6), (6, 8)],
            (12, 25, 39, 54, 35, 35, 17, 17), 234, (1,) * 8, "6315-6378",
            "double cover of F4 with an E6 branch singularity",
        ),
    ]


def _smith(matrix: sp.Matrix) -> tuple[int, ...]:
    diagonal = smith_normal_form(matrix, domain=ZZ).diagonal()
    return tuple(abs(int(value)) for value in diagonal if value)


def _wall_coordinates(gram: sp.Matrix, basis_indices: tuple[int, ...]) -> sp.Matrix:
    basis = gram.extract(basis_indices, basis_indices)
    coordinates = basis.inv() * gram.extract(basis_indices, range(gram.cols))
    if not all(value.is_Integer for value in coordinates):
        raise ArithmeticError("published wall list does not lie integrally in its basis")
    return coordinates


def _components(positive: sp.Matrix, vertices: tuple[int, ...]) -> list[tuple[int, ...]]:
    graph = nx.Graph()
    graph.add_nodes_from(vertices)
    for offset, left in enumerate(vertices):
        for right in vertices[:offset]:
            if positive[left, right] != 0:
                graph.add_edge(left, right)
    return [tuple(sorted(component)) for component in nx.connected_components(graph)]


def _is_elliptic(positive: sp.Matrix, vertices: tuple[int, ...]) -> bool:
    if not vertices:
        return True
    submatrix = positive.extract(vertices, vertices)
    return exact_inertia(submatrix) == (len(vertices), 0, 0)


def _is_affine_component(positive: sp.Matrix, vertices: tuple[int, ...]) -> bool:
    if len(vertices) < 2 or len(_components(positive, vertices)) != 1:
        return False
    submatrix = positive.extract(vertices, vertices)
    if exact_inertia(submatrix) != (len(vertices) - 1, 0, 1):
        return False
    return all(
        _is_elliptic(positive, vertices[:index] + vertices[index + 1 :])
        for index in range(len(vertices))
    )


def _vinberg(negative: sp.Matrix, dimension: int) -> dict[str, Any]:
    """Brute-force Vinberg edge/cusp criterion for at most ten facets."""
    positive = -negative
    facets = positive.rows
    subsets = [
        tuple(combination)
        for size in range(1, facets + 1)
        for combination in itertools.combinations(range(facets), size)
    ]
    elliptic = {vertices for vertices in subsets if _is_elliptic(positive, vertices)}
    affine = {
        vertices for vertices in subsets if _is_affine_component(positive, vertices)
    }
    edges = sorted(vertices for vertices in elliptic if len(vertices) == dimension - 1)
    finite_vertices = {vertices for vertices in elliptic if len(vertices) == dimension}

    bad_neighborhoods: list[dict[str, Any]] = []
    cusp_components: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    for component in sorted(affine):
        neighborhood = tuple(
            vertex
            for vertex in range(facets)
            if vertex in component
            or all(positive[vertex, member] == 0 for member in component)
        )
        parts = _components(positive, neighborhood)
        if all(part in affine for part in parts) and sum(len(part) - 1 for part in parts) == dimension - 1:
            cusp_components[neighborhood] = parts
        else:
            bad_neighborhoods.append(
                {
                    "affine_component_one_based": [value + 1 for value in component],
                    "zero_neighborhood_one_based": [value + 1 for value in neighborhood],
                }
            )

    endpoint_histogram: Counter[tuple[int, int]] = Counter()
    bad_edges: list[dict[str, Any]] = []
    for edge in edges:
        edge_set = set(edge)
        finite_count = sum(edge_set.issubset(vertex) for vertex in finite_vertices)
        ideal_count = 0
        for neighborhood, parts in cusp_components.items():
            for omitted in itertools.product(*parts):
                cusp_edge = set(neighborhood) - set(omitted)
                if cusp_edge == edge_set:
                    ideal_count += 1
        endpoint_histogram[(finite_count, ideal_count)] += 1
        if finite_count + ideal_count != 2:
            bad_edges.append(
                {
                    "edge_one_based": [value + 1 for value in edge],
                    "finite_endpoints": finite_count,
                    "ideal_endpoints": ideal_count,
                }
            )

    finite_volume = not bad_neighborhoods and not bad_edges
    compact = finite_volume and not cusp_components
    return {
        "criterion": (
            "every elliptic rank-(n-1) edge has exactly two finite/ideal "
            "endpoints, and every connected affine zero-neighborhood is a "
            "full rank-(n-1) parabolic diagram"
        ),
        "dimension": dimension,
        "facet_count": facets,
        "elliptic_edge_count": len(edges),
        "finite_vertex_count": len(finite_vertices),
        "connected_affine_component_count": len(affine),
        "ideal_vertex_count": len(cusp_components),
        "endpoint_type_histogram": {
            f"{finite},{ideal}": count
            for (finite, ideal), count in sorted(endpoint_histogram.items())
        },
        "bad_parabolic_neighborhoods": bad_neighborhoods,
        "bad_edges": bad_edges,
        "decision": (
            "compact" if compact else "finite_volume_noncompact" if finite_volume else "infinite_volume"
        ),
        "all_checks_pass": finite_volume,
    }


def _matrix_certificate(case: dict[str, Any]) -> dict[str, Any]:
    gram: sp.Matrix = case["gram"]
    rank = int(case["rank"])
    indices = tuple(case["basis_indices"])
    basis = gram.extract(indices, indices)
    coordinates = _wall_coordinates(gram, indices)
    ample = sp.Matrix(1, rank, case["ample"])
    pairings = tuple(int(value) for value in ample * basis * coordinates)
    square = int((ample * basis * ample.T)[0])
    checks = {
        "header_is_targeted": (rank, gram.rows) in TARGET_HEADERS,
        "symmetric_minus_two_walls": gram == gram.T
        and all(gram[index, index] == -2 for index in range(gram.rows)),
        "nonnegative_curve_intersections": all(
            gram[left, right] >= 0
            for left in range(gram.rows)
            for right in range(left)
        ),
        "rank_and_inertia": gram.rank() == rank
        and exact_inertia(gram) == (1, rank - 1, gram.rows - rank),
        "published_basis_is_unimodular_wall_span": abs(int(coordinates[:, indices].det())) == 1,
        "ample_square": square == case["ample_square"],
        "ample_pairings": pairings == case["ample_pairings"],
        "ample_is_strict": min(pairings) > 0,
    }
    if not all(checks.values()):
        failed = [name for name, passed in checks.items() if not passed]
        raise ArithmeticError(f"{case['id']} matrix checks failed: {failed}; pairings={pairings}; square={square}")
    return {
        "wall_gram": _payload(gram),
        "wall_gram_sha256": _canonical_sha256(_payload(gram)),
        "basis_indices_one_based": [value + 1 for value in indices],
        "basis_determinant": int(basis.det()),
        "basis_smith_invariants": list(_smith(basis)),
        "basis_inertia": list(exact_inertia(basis)),
        "wall_rank": gram.rank(),
        "wall_nullity": gram.rows - gram.rank(),
        "ample_coordinates_in_basis": list(case["ample"]),
        "ample_square": square,
        "ample_wall_pairings": list(pairings),
        "checks": checks,
        "all_checks_pass": True,
    }


def analyze() -> dict[str, Any]:
    cases = []
    seen_headers: Counter[tuple[int, int]] = Counter()
    for source in _cases():
        matrix = _matrix_certificate(source)
        vinberg = _vinberg(source["gram"], source["rank"] - 1)
        if not vinberg["all_checks_pass"]:
            raise ArithmeticError(f"{source['id']} failed the exact Vinberg criterion")
        seen_headers[(source["rank"], source["gram"].rows)] += 1
        cases.append(
            {
                "id": source["id"],
                "lattice": source["lattice"],
                "picard_rank": source["rank"],
                "hyperbolic_dimension": source["rank"] - 1,
                "facet_count": source["gram"].rows,
                "tex_lines": source["tex_lines"],
                "branch_model": source["branch_model"],
                "complete_curve_source": (
                    "Roulleau's stated curve configuration together with the final "
                    "table's complete irreducible (-2)-curve count"
                ),
                "complete_full_nef_reason": (
                    "the published walls are the complete irreducible (-2)-curve "
                    "set; the K3 ample-cone theorem makes their positive-cone "
                    "chamber the full nef cone"
                ),
                "matrix_certificate": matrix,
                "literal_vinberg": vinberg,
                "lineage": {
                    "published_status": (
                        "published Roulleau complete K3 root chamber; this compiler "
                        "makes no abstract-polytope priority claim"
                    ),
                    "primary_source": ROULLEAU_ARXIV,
                    "primary_source_tex_lines": source["tex_lines"],
                },
                "classification": (
                    ("compact" if vinberg["decision"] == "compact" else "finite_volume_noncompact")
                    + "_full_nef_K3_Coxeter_polytope"
                ),
                "all_checks_pass": True,
            }
        )

    expected_headers = {
        (5, 5): 2,
        (6, 6): 3,
        (7, 7): 2,
        (8, 8): 3,
        (8, 9): 5,
        (8, 10): 2,
    }
    if dict(seen_headers) != expected_headers or len(cases) != 17:
        raise ArithmeticError(f"Roulleau targeted header census changed: {seen_headers}")
    decisions = Counter(case["literal_vinberg"]["decision"] for case in cases)
    payload = {
        "schema": SCHEMA,
        "scope": {
            "selection": (
                "all Roulleau final-table rank-5--8 complete root systems whose "
                "(Picard rank, root count) header lies in the six previously "
                "header-empty stable ADEG target headers"
            ),
            "target_headers_rank_root_count": [list(header) for header in sorted(TARGET_HEADERS)],
            "exclusion_boundary": (
                "not all rank-5--8 finite-automorphism K3 lattices, not an all-rank "
                "K3 chamber census, and not an abstract-polytope novelty claim"
            ),
        },
        "primary_source": {
            "citation": (
                "Xavier Roulleau, An atlas of K3 surfaces with finite automorphism group, "
                "EPIGA 6 (2022), arXiv:2003.08985"
            ),
            "url": ROULLEAU_ARXIV,
            "source_tar_sha256": ROULLEAU_SOURCE_TAR_SHA256,
            "source_tex_sha256": ROULLEAU_SOURCE_TEX_SHA256,
            "final_table_curve_count_statement_tex_line": 10977,
        },
        "cases": cases,
        "summary": {
            "stock_count": len(cases),
            "header_stock_counts": {
                f"{rank},{roots}": count
                for (rank, roots), count in sorted(seen_headers.items())
            },
            "decision_counts": dict(sorted(decisions.items())),
            "all_complete_full_nef": True,
            "all_literal_vinberg_finite_volume": True,
        },
        "all_checks_pass": True,
    }
    payload["canonical_sha256"] = _canonical_sha256(payload)
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    arguments = parser.parse_args()
    payload = analyze()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.write:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(encoded, encoding="utf-8")
    if arguments.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != encoded:
            raise SystemExit("archive does not match exact replay")
    if not arguments.quiet:
        print(encoded, end="")


if __name__ == "__main__":
    main()
