#!/usr/bin/env python3
"""Certify the integral and literature lineage of the D5 H10 30-wall chamber.

The complete D5 subgroup census supplies the involution, its 30 irredundant
facet orbits, and the exact facet Gram matrix.  This certificate computes the
saturated integral fixed lattice, primitive-normalizes every mirror, and
places all reflections in a standard odd-unimodular Lorentzian lattice.

The resulting lattice is ``2*I_{1,10}``.  Its full root-reflection chamber is
Vinberg's classical ``P10``.  We replay this identification in three pinned
data sources: Alexeev's ``(11,11,1)`` matrix, Scharlau--Kirschmer List10 entry
1, and CoxIter's unit-form record.  Solomon's formula then gives the exact
index of the 30-wall reflection subgroup in the full P10 group.

The corpus exclusion is deliberately bounded.  It is not a worldwide
priority claim for the literal 30-wall chamber.
"""

from __future__ import annotations

import argparse
import ast
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import tarfile
from typing import Any, Sequence

import networkx as nx
import sympy as sp
from sympy import ZZ
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp

from scripts.audit_jacobian_kummer_h12_parallel_cusps import _positive_definite
from scripts.certify_jacobian_kummer_h10_maximal_reflection_lineage import (
    _parse_coxiter,
    _permutation,
)
from scripts.certify_jacobian_kummer_h10_reflection_index import (
    _connected_masks_direct,
    _euler_from_connected_components,
    _factorization,
    _relation_matrix,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
CENSUS = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_full_subgroup_fishing.json"
)
CONE = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h10_30_exact_cone_certificate.json"
)
VINBERG = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h10_30_vinberg_face_certificate.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h10_30_lineage.json"
)

COXITER_CORPUS = Path("/private/tmp/coxiter-current/CoxIter-master")
COXITER_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
COXITER_ARCHIVE_SHA256 = (
    "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
)
DOUBLING_ARCHIVE_RELATIVE = Path("graphs/testing_index2.tar.gz")
DOUBLING_ARCHIVE_SHA256 = (
    "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"
)
P10_COXITER_RELATIVE = Path(
    "graphs/group units qf/10-1,1,1,1,1,1,1,1,1,1,1.coxiter"
)
P10_COXITER_SHA256 = (
    "a8e8618529200271e4be249442ec283c65f222372ece8d1910476433f57ce43e"
)

SK_DIRECTORY = Path("/private/tmp/Scharlau-Kirschmer-Lorentz")
SK_ARCHIVE = Path("/private/tmp/Scharlau-Kirschmer-Lorentz.zip")
SK_ARCHIVE_SHA256 = (
    "74a4dd3e534f1e793cba45765ff33575dbd331d5cdce6f34cdd64340b3f0d4a5"
)
SK_LIST10_SHA256 = (
    "9e627c4e1303015fd96acdcf6ccbeb264ebadba49abfcd761705f54d211b2e21"
)
SK_GRAM10_SHA256 = (
    "b8186e1fa2771c517a6040998f3b234568c3edde3ec8b7bb8412865d7710af18"
)

ALEXEEV_DIRECTORY = Path("/private/tmp/alexeev_2209_source")
ALEXEEV_ARCHIVE = Path("/private/tmp/alexeev-2209.09110v4.tar")
ALEXEEV_ARCHIVE_SHA256 = (
    "10e4002f6518f8bf6d669946775bbe42731f69277a52847c5602fa60ea33c2a6"
)
ALEXEEV_MATRIX_FILE_SHA256 = (
    "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457"
)
ALEXEEV_TEX_SHA256 = (
    "006e2a1eac2aa3c4916ef66aceb27e8e80a35145162e6785a4874843f249f1df"
)
ALEXEEV_11_11_1_MATRIX_SHA256 = (
    "c41902c3ed657ee99cd018a8883b9c8bcb5e313ddb14a39f0bdf448497014717"
)


# Columns are the standard ``I_{1,10}`` basis in the reduced saturated fixed
# basis used below.  Thus Q^T * (G_fixed/2) * Q = diag(1,-1,...,-1).
FIXED_TO_I_1_10 = sp.Matrix(
    [
        [24, 1, 0, 0, 1, 1, 0, -3, 0, -6, -6],
        [28, 0, 1, -1, 0, 0, 0, 0, 5, -7, -7],
        [-21, 0, 0, 1, 0, 0, 0, 0, -4, 5, 5],
        [24, 0, 0, 0, 1, 1, 0, -3, 0, -6, -6],
        [6, 0, 0, 0, 0, 0, 0, -1, 0, -2, -1],
        [24, 0, 0, 0, 0, 1, 0, -3, 0, -6, -6],
        [4, 0, 0, 0, 0, 0, 1, 0, 0, -1, -1],
        [-5, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1],
        [-8, 0, 0, 0, 0, 0, 0, 0, -1, 2, 2],
        [-4, 0, 0, 0, 0, 0, 0, 0, -1, 1, 1],
        [-8, 0, 0, 0, 0, 0, 0, 1, 0, 2, 2],
    ]
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _matrix(payload: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(entry) for entry in row] for row in payload])


def _matrix_payload(matrix: sp.Matrix) -> list[list[int]]:
    if any(value.is_Integer is not True for value in matrix):
        raise ArithmeticError("expected an integral matrix")
    return [[int(value) for value in row] for row in matrix.tolist()]


def _matrix_sha256(matrix: sp.Matrix) -> str:
    payload = ";".join(
        ",".join(str(int(matrix[row, column])) for column in range(matrix.cols))
        for row in range(matrix.rows)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _smith_decomposition(matrix: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    diagonal, left, right = smith_normal_decomp(
        DomainMatrix.from_Matrix(matrix).convert_to(ZZ)
    )
    return diagonal.to_Matrix(), left.to_Matrix(), right.to_Matrix()


def _primitive_integer_column(vector: sp.Matrix) -> tuple[sp.Matrix, sp.Rational]:
    denominator = int(sp.ilcm(*[sp.Rational(value).q for value in vector]))
    entries = [int(sp.Rational(value) * denominator) for value in vector]
    content = math.gcd(*(abs(entry) for entry in entries))
    if content == 0:
        raise ArithmeticError("zero mirror vector")
    primitive = sp.Matrix([entry // content for entry in entries])
    if next(entry for entry in primitive if entry) < 0:
        primitive = -primitive
    scale = next(
        sp.Rational(vector[index], primitive[index])
        for index in range(vector.rows)
        if primitive[index]
    )
    if vector != scale * primitive:
        raise ArithmeticError("primitive normalization did not replay")
    return primitive, scale


def _root_type(gram: sp.Matrix, root: sp.Matrix) -> tuple[int, int]:
    square = int((root.T * gram * root)[0])
    divisibility = math.gcd(*(abs(int(value)) for value in gram * root))
    return square, divisibility


def _reflection_is_integral(gram: sp.Matrix, root: sp.Matrix) -> bool:
    square = int((root.T * gram * root)[0])
    reflection = sp.eye(gram.rows) - sp.Rational(2, square) * root * (
        root.T * gram
    )
    return (
        all(value.is_Integer is True for value in reflection)
        and reflection.T * gram * reflection == gram
        and reflection * reflection == sp.eye(gram.rows)
        and reflection * root == -root
    )


def _connected_elliptic_masks(
    positive_gram_payload: Sequence[Sequence[str]],
) -> tuple[set[int], list[dict[str, int]]]:
    rational = [
        [sp.Rational(entry) for entry in row] for row in positive_gram_payload
    ]
    scale = math.lcm(*(int(value.q) for row in rational for value in row))
    gram = [[int(scale * value) for value in row] for row in rational]
    size = len(gram)
    compatible = [
        {
            other
            for other in range(size)
            if other != vertex
            and (
                gram[vertex][other] == 0
                or gram[vertex][vertex] * gram[other][other]
                > gram[vertex][other] * gram[vertex][other]
            )
        }
        for vertex in range(size)
    ]
    finite_neighbors = [
        {other for other in compatible[vertex] if gram[vertex][other] != 0}
        for vertex in range(size)
    ]

    elliptic: set[tuple[int, ...]] = {(vertex,) for vertex in range(size)}
    masks = {1 << vertex for vertex in range(size)}
    order_counts: list[dict[str, int]] = []
    for order in range(2, 13):
        candidates: set[tuple[int, ...]] = set()
        for facets in elliptic:
            facet_set = set(facets)
            possible = set.intersection(*(compatible[node] for node in facets))
            connected_extensions = set.union(
                *(finite_neighbors[node] for node in facets)
            )
            for vertex in possible & connected_extensions - facet_set:
                candidates.add(tuple(sorted((*facets, vertex))))
        next_elliptic: set[tuple[int, ...]] = set()
        for facets in candidates:
            subgram = [
                [gram[left][right] for right in facets] for left in facets
            ]
            if _positive_definite(subgram):
                next_elliptic.add(facets)
                masks.add(sum(1 << vertex for vertex in facets))
        order_counts.append(
            {
                "order": order,
                "candidate_count": len(candidates),
                "elliptic_count": len(next_elliptic),
            }
        )
        elliptic = next_elliptic
        if not elliptic:
            break
    return masks, order_counts


def _parse_coxiter_text(text: str) -> tuple[int, int, list[list[str]]]:
    lines = [
        content
        for line in text.splitlines()
        if (content := line.split("#", 1)[0].strip())
    ]
    facets, dimension = map(int, lines[0].split()[:2])
    matrix = [
        ["1" if left == right else "2" for right in range(facets)]
        for left in range(facets)
    ]
    vertex_line = next(
        (line for line in lines[1:] if line.startswith("vertices labels:")), None
    )
    labels = (
        vertex_line.split(":", 1)[1].split()
        if vertex_line is not None
        else [str(index + 1) for index in range(facets)]
    )
    if len(labels) != facets or len(set(labels)) != facets:
        raise ArithmeticError("invalid CoxIter vertex labels")
    index_of = {label: index for index, label in enumerate(labels)}
    for line in lines[1:]:
        if line.startswith("vertices labels:"):
            continue
        tokens = line.split()
        if len(tokens) != 3:
            raise ArithmeticError(f"unsupported CoxIter line: {line}")
        left, right = index_of[tokens[0]], index_of[tokens[1]]
        matrix[left][right] = matrix[right][left] = tokens[2]
    return facets, dimension, matrix


def _corpus_audit(
    corpus: Path, candidate_relations: Sequence[Sequence[str]]
) -> dict[str, Any]:
    ordinary_count = 0
    ordinary_h10: Counter[int] = Counter()
    ordinary_candidates: list[dict[str, Any]] = []
    ordinary_matches: list[dict[str, Any]] = []
    for path in sorted(corpus.rglob("*.coxiter")):
        facets, dimension, relations = _parse_coxiter(path)
        ordinary_count += 1
        if dimension == 10:
            ordinary_h10[facets] += 1
        if (facets, dimension) != (30, 10):
            continue
        relative = path.relative_to(corpus).as_posix()
        ordinary_candidates.append({"path": relative})
        permutation = _permutation(candidate_relations, relations)
        if permutation is not None:
            ordinary_matches.append(
                {
                    "path": relative,
                    "candidate_to_corpus_zero_based_permutation": permutation,
                    "entrywise_relation_matrix_match": True,
                }
            )

    doubling_path = corpus / DOUBLING_ARCHIVE_RELATIVE
    doubling_count = 0
    doubling_h10: Counter[int] = Counter()
    doubling_candidates: list[dict[str, Any]] = []
    doubling_matches: list[dict[str, Any]] = []
    with tarfile.open(doubling_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"failed to read {member.name}")
            text = handle.read().decode("utf-8")
            facets, dimension, relations = _parse_coxiter_text(text)
            doubling_count += 1
            if dimension == 10:
                doubling_h10[facets] += 1
            if (facets, dimension) != (30, 10):
                continue
            doubling_candidates.append({"path": member.name})
            permutation = _permutation(candidate_relations, relations)
            if permutation is not None:
                doubling_matches.append(
                    {
                        "path": member.name,
                        "candidate_to_corpus_zero_based_permutation": permutation,
                        "entrywise_relation_matrix_match": True,
                    }
                )
    return {
        "repository": "https://github.com/rgugl/CoxIter",
        "commit": COXITER_COMMIT,
        "ordinary_file_count": ordinary_count,
        "ordinary_H10_facet_histogram": {
            str(facets): count for facets, count in sorted(ordinary_h10.items())
        },
        "ordinary_H10_30_candidates": ordinary_candidates,
        "ordinary_exact_weighted_diagram_matches": ordinary_matches,
        "doubling_archive_relative_path": DOUBLING_ARCHIVE_RELATIVE.as_posix(),
        "doubling_archive_sha256": _sha256(doubling_path),
        "doubling_file_count": doubling_count,
        "doubling_H10_facet_histogram": {
            str(facets): count for facets, count in sorted(doubling_h10.items())
        },
        "doubling_H10_30_candidates": doubling_candidates,
        "doubling_exact_weighted_diagram_matches": doubling_matches,
    }


def _alexeev_matrix(path: Path) -> sp.Matrix:
    prefix = "coxeter_matrix[(11,11,1)] = matrix("
    matching = [
        line for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(prefix)
    ]
    if len(matching) != 1 or not matching[0].endswith(")"):
        raise ArithmeticError("Alexeev's (11,11,1) matrix line changed")
    matrix = sp.Matrix(ast.literal_eval(matching[0][len(prefix) : -1]))
    if matrix.shape != (12, 12) or matrix != matrix.T:
        raise ArithmeticError("Alexeev's (11,11,1) matrix is malformed")
    return matrix


def _weighted_gram_permutation(source: sp.Matrix, target: sp.Matrix) -> list[int] | None:
    def graph(matrix: sp.Matrix) -> nx.Graph:
        result = nx.Graph()
        for index in range(matrix.rows):
            result.add_node(index, square=int(matrix[index, index]))
        for left in range(matrix.rows):
            for right in range(left + 1, matrix.cols):
                result.add_edge(left, right, product=int(matrix[left, right]))
        return result

    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph(source),
        graph(target),
        node_match=nx.algorithms.isomorphism.categorical_node_match(
            "square", None
        ),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match(
            "product", None
        ),
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    if mapping is None:
        return None
    permutation = [mapping[index] for index in range(source.rows)]
    if any(
        source[left, right] != target[permutation[left], permutation[right]]
        for left in range(source.rows)
        for right in range(source.cols)
    ):
        raise ArithmeticError("weighted Gram permutation did not replay")
    return permutation


def _p10_roots() -> sp.Matrix:
    basis = [sp.eye(11)[:, index] for index in range(11)]
    roots = [basis[index] - basis[index + 1] for index in range(1, 10)]
    roots.append(basis[10])
    roots.append(basis[0] - basis[1] - basis[2] - basis[3])
    roots.append(
        3 * basis[0]
        - sum((basis[index] for index in range(1, 11)), sp.zeros(11, 1))
    )
    return sp.Matrix.hstack(*roots)


def _parse_sk_list10_first(path: Path) -> dict[str, Any]:
    line = next(
        content.strip()
        for content in path.read_text(encoding="utf-8").splitlines()
        if '"I_{10,1}()"' in content
    )
    match = re.fullmatch(
        r'<(\d+),\s*"([^"]+)",\s*\[\s*([^]]+)\s*\],\s*(\d+)>,?',
        line,
    )
    if match is None:
        raise ArithmeticError("Scharlau--Kirschmer List10 entry 1 changed format")
    return {
        "determinant": int(match.group(1)),
        "genus_symbol": match.group(2),
        "finite_face_counts_facets_to_vertices": [
            int(value.strip()) for value in match.group(3).split(",")
        ],
        "cusp_count": int(match.group(4)),
    }


def build_certificate(
    coxiter_corpus: Path = COXITER_CORPUS,
    coxiter_archive: Path = COXITER_ARCHIVE,
    sk_directory: Path = SK_DIRECTORY,
    sk_archive: Path = SK_ARCHIVE,
    alexeev_directory: Path = ALEXEEV_DIRECTORY,
    alexeev_archive: Path = ALEXEEV_ARCHIVE,
) -> dict[str, Any]:
    census = json.loads(CENSUS.read_text(encoding="utf-8"))
    cone = json.loads(CONE.read_text(encoding="utf-8"))
    vinberg = json.loads(VINBERG.read_text(encoding="utf-8"))
    section = census["above_H9"]["complete_section_screens"]["class_03_H10"]
    if not (
        section["fixed_rank"] == 11
        and section["hyperbolic_dimension"] == 10
        and section["facet_count"] == 30
        and section["pair_gate"]["all_pairs_coxeter"]
        and cone["decision"]["finite_volume"]
        and not cone["decision"]["compact"]
        and cone["exact_double_description"]["hyperideal_vertex_count"] == 0
        and vinberg["decision"]["finite_volume"]
        and vinberg["literal_endpoint_audit"]["bad_edge_count"] == 0
        and vinberg["literal_endpoint_audit"]["bad_parabolic_neighborhood_count"]
        == 0
        and _sha256(CENSUS) == cone["source"]["census_sha256"]
    ):
        raise ArithmeticError("an upstream H10-30 gate failed")

    ambient_gram = _matrix(census["ambient_lattice"]["gram"])
    raw_matrices = [
        _matrix(payload)
        for payload in section["representative_subgroup"][
            "raw_primal_row_matrices"
        ]
    ]
    involution = next(matrix for matrix in raw_matrices if matrix != sp.eye(20))
    if not (
        involution**2 == sp.eye(20)
        and involution * ambient_gram * involution.T == ambient_gram
    ):
        raise ArithmeticError("the archived involution is not an integral isometry")

    constraint = involution.T - sp.eye(20)
    smith_constraint, _, smith_right = _smith_decomposition(constraint)
    equation_rank = constraint.rank()
    smith_basis = smith_right[:, equation_rank:].T
    reduced_basis = smith_basis.lll(delta=sp.Rational(3, 4))
    _, pivot_columns = smith_basis.rref()
    basis_change = (
        reduced_basis[:, list(pivot_columns)]
        * smith_basis[:, list(pivot_columns)].inv()
    )
    fixed_gram = reduced_basis * ambient_gram * reduced_basis.T
    fixed_snf = smith_normal_form(fixed_gram, domain=ZZ)
    fixed_snf_diagonal = [
        abs(int(fixed_snf[index, index])) for index in range(11)
    ]
    half_gram = fixed_gram / 2
    target_gram = sp.diag(1, *([-1] * 10))
    isometry = FIXED_TO_I_1_10
    if not (
        equation_rank == 9
        and smith_basis.rank() == 11
        and smith_basis * involution == smith_basis
        and reduced_basis * involution == reduced_basis
        and basis_change * smith_basis == reduced_basis
        and abs(int(basis_change.det())) == 1
        and fixed_gram.det() == 2**11
        and exact_inertia(fixed_gram) == (1, 10, 0)
        and fixed_snf_diagonal == [2] * 11
        and all(value.is_Integer is True for value in half_gram)
        and half_gram.det() == 1
        and any(int(half_gram[index, index]) % 2 for index in range(11))
        and abs(int(isometry.det())) == 1
        and isometry.T * half_gram * isometry == target_gram
    ):
        raise ArithmeticError("the saturated fixed-lattice identification failed")

    walls = tuple(
        sp.Matrix(1, 20, wall)
        for wall in census["chamber"]["published_primitive_dual_walls"]
    )
    wall_orbits = tuple(
        tuple(index - 1 for index in orbit) for orbit in section["wall_orbits"]
    )
    facets = tuple(index - 1 for index in section["facet_orbit_indices"])
    mirror_records: list[dict[str, Any]] = []
    fixed_types: Counter[tuple[int, int]] = Counter()
    half_types: Counter[tuple[int, int]] = Counter()
    presentation_types: Counter[tuple[str, int, str]] = Counter()
    target_roots: list[sp.Matrix] = []
    for facet_position, orbit_index in enumerate(facets):
        presented_dual_normal = sum(
            (walls[index] for index in wall_orbits[orbit_index]),
            sp.zeros(1, 20),
        )
        fixed_coordinates_rational = fixed_gram.inv() * (
            reduced_basis * presented_dual_normal.T
        )
        primitive, scale = _primitive_integer_column(fixed_coordinates_rational)
        if (
            fixed_coordinates_rational.T * reduced_basis
            != presented_dual_normal * ambient_gram.inv()
        ):
            raise ArithmeticError("a mirror was not recovered in the fixed span")
        target_coordinates = isometry.inv() * primitive
        if any(value.is_Integer is not True for value in target_coordinates):
            raise ArithmeticError("a mirror left I_{1,10}")
        presented_square = sp.Rational(
            (
                presented_dual_normal
                * ambient_gram.inv()
                * presented_dual_normal.T
            )[0]
        )
        fixed_square, fixed_divisibility = _root_type(fixed_gram, primitive)
        half_square, half_divisibility = _root_type(half_gram, primitive)
        reflection_valid = _reflection_is_integral(fixed_gram, primitive)
        if not (
            presented_square == scale**2 * fixed_square
            and fixed_square == 2 * half_square
            and fixed_divisibility == 2 * half_divisibility
            and reflection_valid
        ):
            raise ArithmeticError("a primitive mirror integrality check failed")
        fixed_types[(fixed_square, fixed_divisibility)] += 1
        half_types[(half_square, half_divisibility)] += 1
        presentation_types[(str(presented_square), fixed_square, str(scale))] += 1
        target_roots.append(target_coordinates)
        mirror_records.append(
            {
                "facet_position": facet_position + 1,
                "wall_orbit_index": orbit_index + 1,
                "source_wall_indices": [
                    index + 1 for index in wall_orbits[orbit_index]
                ],
                "presented_dual_normal": [
                    int(value) for value in presented_dual_normal
                ],
                "presented_square": str(presented_square),
                "presented_equals_scale_times_primitive": str(scale),
                "primitive_fixed_coordinates": [int(value) for value in primitive],
                "primitive_I_1_10_coordinates": [
                    int(value) for value in target_coordinates
                ],
                "fixed_lattice_square": fixed_square,
                "fixed_lattice_divisibility": fixed_divisibility,
                "half_scaled_square": half_square,
                "half_scaled_divisibility": half_divisibility,
                "reflection_integral_involutive_isometry": reflection_valid,
            }
        )
    target_root_matrix = sp.Matrix.hstack(*target_roots)
    root_span_snf = smith_normal_form(target_root_matrix, domain=ZZ)
    root_span_diagonal = [
        abs(int(root_span_snf[index, index])) for index in range(11)
    ]
    if not (
        fixed_types == Counter({(-2, 2): 18, (-4, 2): 12})
        and half_types == Counter({(-1, 1): 18, (-2, 1): 12})
        and target_root_matrix.rank() == 11
        and root_span_diagonal == [1] * 11
    ):
        raise ArithmeticError("the 30-root type or integral-span audit failed")

    p10_roots = _p10_roots()
    p10_gram = p10_roots.T * (2 * target_gram) * p10_roots
    p10_types = Counter(
        _root_type(2 * target_gram, p10_roots[:, index])
        for index in range(p10_roots.cols)
    )
    p10_reflections_integral = all(
        _reflection_is_integral(2 * target_gram, p10_roots[:, index])
        for index in range(p10_roots.cols)
    )
    p10_relations = _relation_matrix(
        [[str(-p10_gram[row, column]) for column in range(12)] for row in range(12)]
    )
    p10_path = coxiter_corpus / P10_COXITER_RELATIVE
    p10_facets, p10_dimension, coxiter_p10_relations = _parse_coxiter(p10_path)
    coxiter_p10_permutation = _permutation(p10_relations, coxiter_p10_relations)

    alexeev_matrix_path = alexeev_directory / "all_cases_in_2el_paper.txt"
    alexeev_tex_path = alexeev_directory / "g0.tex"
    alexeev_matrix = _alexeev_matrix(alexeev_matrix_path)
    alexeev_permutation = _weighted_gram_permutation(p10_gram, alexeev_matrix)
    alexeev_tex = " ".join(alexeev_tex_path.read_text(encoding="utf-8").split())

    sk_list_path = sk_directory / "List10"
    sk_gram_path = sk_directory / "Gram10"
    sk_entry = _parse_sk_list10_first(sk_list_path)
    sk_grams = ast.literal_eval(sk_gram_path.read_text(encoding="utf-8"))
    sk_gram_entry1 = sp.Matrix(11, 11, sk_grams[0])
    sk_target_permutation = list(range(1, 11)) + [0]
    sk_signed_target = (-target_gram).extract(
        sk_target_permutation, sk_target_permutation
    )

    candidate_relations = _relation_matrix(section["positive_wall_gram"])
    connected_masks, component_search = _connected_elliptic_masks(
        section["positive_wall_gram"]
    )
    candidate_euler = _euler_from_connected_components(
        candidate_relations, connected_masks
    )
    p10_connected_masks = _connected_masks_direct(p10_relations)
    p10_euler = _euler_from_connected_components(p10_relations, p10_connected_masks)
    candidate_chi = Fraction(candidate_euler["rational_euler_characteristic"])
    p10_chi = Fraction(p10_euler["rational_euler_characteristic"])
    index = candidate_chi / p10_chi
    volume_coefficient = -candidate_chi * Fraction(32, 945)
    full_face_vector = [
        cone["exact_double_description"]["finite_vertex_count"]
        + cone["exact_double_description"]["ideal_vertex_count"],
        *reversed(candidate_euler["finite_special_subset_count_by_rank"][1:10]),
        1,
    ]

    corpus_audit = _corpus_audit(coxiter_corpus, candidate_relations)
    expected_h10_ordinary = {"11": 14, "12": 8, "13": 4, "14": 1}
    expected_h10_doubling = {
        "12": 3,
        "13": 4,
        "14": 3,
        "15": 1,
        "16": 3,
        "19": 1,
        "20": 3,
        "28": 2,
        "44": 2,
    }
    expected_component_counts = {
        1: 30,
        2: 45,
        3: 90,
        4: 126,
        5: 114,
        6: 120,
        7: 108,
        8: 132,
        9: 96,
    }
    actual_component_counts = Counter(mask.bit_count() for mask in connected_masks)
    expected_candidate_special = [
        1,
        30,
        321,
        1_890,
        7_072,
        17_784,
        30_452,
        34_857,
        25_206,
        10_030,
        1_333,
    ]
    expected_p10_special = [1, 12, 65, 210, 450, 672, 714, 540, 285, 99, 17]
    expected_face_vector = [
        1_540,
        10_030,
        25_206,
        34_857,
        30_452,
        17_784,
        7_072,
        1_890,
        321,
        30,
        1,
    ]
    expected_alexeev_permutation = [2, 3, 4, 5, 6, 7, 8, 10, 0, 11, 1, 9]
    checks = {
        "census_and_both_volume_certificates_replay": True,
        "fixed_kernel_is_saturated_rank_11": (
            smith_basis.rank() == 11 and abs(int(basis_change.det())) == 1
        ),
        "fixed_lattice_is_2_times_odd_unimodular_I_1_10": (
            isometry.T * half_gram * isometry == target_gram
        ),
        "explicit_isometry_has_determinant_minus_one": int(isometry.det()) == -1,
        "all_30_mirrors_are_integral_reflecting_roots": all(
            record["reflection_integral_involutive_isometry"]
            for record in mirror_records
        ),
        "thirty_roots_generate_the_whole_fixed_lattice": root_span_diagonal
        == [1] * 11,
        "p10_explicit_root_types_are_expected": p10_types
        == Counter({(-4, 2): 10, (-2, 2): 2}),
        "all_p10_reflections_are_integral": p10_reflections_integral,
        "p10_coxiter_record_is_12_wall_H10": (
            p10_facets == 12 and p10_dimension == 10
        ),
        "p10_coxiter_relation_matrix_matches_entrywise": coxiter_p10_permutation
        == list(range(12)),
        "alexeev_sources_are_pinned": (
            _sha256(alexeev_archive) == ALEXEEV_ARCHIVE_SHA256
            and _sha256(alexeev_matrix_path) == ALEXEEV_MATRIX_FILE_SHA256
            and _sha256(alexeev_tex_path) == ALEXEEV_TEX_SHA256
            and _matrix_sha256(alexeev_matrix)
            == ALEXEEV_11_11_1_MATRIX_SHA256
        ),
        "alexeev_calls_lattice_I_1_10_2": (
            "I_{1,10}(2)" in alexeev_tex and "has $12$ roots" in alexeev_tex
        ),
        "p10_gram_matches_alexeev_11_11_1_entrywise": alexeev_permutation
        == expected_alexeev_permutation,
        "scharlau_kirschmer_sources_are_pinned": (
            _sha256(sk_archive) == SK_ARCHIVE_SHA256
            and _sha256(sk_list_path) == SK_LIST10_SHA256
            and _sha256(sk_gram_path) == SK_GRAM10_SHA256
        ),
        "scharlau_kirschmer_gram10_entry1_is_signed_I_1_10": (
            sk_gram_entry1 == sk_signed_target
        ),
        "scharlau_kirschmer_list10_entry1_matches_p10_faces": (
            sk_entry
            == {
                "determinant": 1,
                "genus_symbol": "I_{10,1}()",
                "finite_face_counts_facets_to_vertices": expected_p10_special[1:],
                "cusp_count": 2,
            }
        ),
        "connected_elliptic_census_is_expected": actual_component_counts
        == expected_component_counts,
        "candidate_finite_special_counts_match_vinberg_faces": (
            candidate_euler["finite_special_subset_count_by_rank"]
            == expected_candidate_special
            and expected_candidate_special[9]
            == vinberg["literal_endpoint_audit"]["edge_count"]
            and expected_candidate_special[10]
            == vinberg["literal_endpoint_audit"]["finite_vertex_count"]
        ),
        "p10_finite_special_counts_match_scharlau_kirschmer": (
            p10_euler["finite_special_subset_count_by_rank"]
            == expected_p10_special
        ),
        "candidate_euler_is_expected": candidate_chi
        == Fraction(-527, 483_840),
        "p10_euler_is_expected": p10_chi == Fraction(-1, 11_147_673_600),
        "euler_ratio_is_positive_integer": index.denominator == 1 and index > 0,
        "index_is_expected": index == 12_142_080,
        "volume_is_expected": volume_coefficient == Fraction(527, 14_288_400),
        "full_face_vector_reconstructs": full_face_vector == expected_face_vector,
        "coxiter_sources_are_pinned": (
            _sha256(coxiter_archive) == COXITER_ARCHIVE_SHA256
            and _sha256(p10_path) == P10_COXITER_SHA256
            and corpus_audit["doubling_archive_sha256"]
            == DOUBLING_ARCHIVE_SHA256
        ),
        "ordinary_corpus_scope_is_expected": (
            corpus_audit["ordinary_file_count"] == 1_214
            and corpus_audit["ordinary_H10_facet_histogram"]
            == expected_h10_ordinary
        ),
        "doubling_corpus_scope_is_expected": (
            corpus_audit["doubling_file_count"] == 2_959
            and corpus_audit["doubling_H10_facet_histogram"]
            == expected_h10_doubling
        ),
        "no_bounded_corpus_H10_30_candidate_or_match": (
            corpus_audit["ordinary_H10_30_candidates"] == []
            and corpus_audit["ordinary_exact_weighted_diagram_matches"] == []
            and corpus_audit["doubling_H10_30_candidates"] == []
            and corpus_audit["doubling_exact_weighted_diagram_matches"] == []
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"D5 H10-30 lineage checks failed: {checks}")

    return {
        "schema": "shimada_disc16_d5_h10_30_lineage_v1",
        "candidate": {
            "id": "D5-C2-H10-30",
            "source_chamber": "Shimada disc16 induced chamber D5",
            "hyperbolic_dimension": 10,
            "facet_count": 30,
            "finite_vertex_count": 1_333,
            "ideal_vertex_count": 207,
            "compact": False,
            "finite_volume": True,
        },
        "sources": {
            "subgroup_census": str(CENSUS.relative_to(ROOT)),
            "subgroup_census_sha256": _sha256(CENSUS),
            "exact_cone_certificate": str(CONE.relative_to(ROOT)),
            "exact_cone_certificate_sha256": _sha256(CONE),
            "literal_vinberg_certificate": str(VINBERG.relative_to(ROOT)),
            "literal_vinberg_certificate_sha256": _sha256(VINBERG),
            "vinberg_primary_source": (
                "https://doi.org/10.1070/SM1972v016n01ABEH001346"
            ),
            "alexeev_primary_source": "https://arxiv.org/abs/2209.09110",
            "alexeev_archive_sha256": _sha256(alexeev_archive),
            "alexeev_matrix_file_sha256": _sha256(alexeev_matrix_path),
            "scharlau_kirschmer_source": (
                "https://www.math.rwth-aachen.de/~Markus.Kirschmer/lorentz/index.html"
            ),
            "scharlau_kirschmer_archive_sha256": _sha256(sk_archive),
            "coxiter_repository_commit": COXITER_COMMIT,
            "coxiter_archive_sha256": _sha256(coxiter_archive),
        },
        "saturated_fixed_lattice": {
            "definition": "Fix(M) = ker_Z(M^T-I) for Shimada row action M",
            "constraint_smith_diagonal": [
                int(smith_constraint[index, index]) for index in range(20)
            ],
            "rank": 11,
            "saturated": True,
            "smith_kernel_basis_rows_in_NS": _matrix_payload(smith_basis),
            "reduced_basis_rows_in_NS": _matrix_payload(reduced_basis),
            "reduced_to_smith_basis_change": _matrix_payload(basis_change),
            "basis_change_determinant": int(basis_change.det()),
            "gram": _matrix_payload(fixed_gram),
            "inertia": list(exact_inertia(fixed_gram)),
            "determinant": int(fixed_gram.det()),
            "smith_diagonal": fixed_snf_diagonal,
            "even": True,
            "half_scaled_gram": _matrix_payload(half_gram),
            "half_scaled_determinant": int(half_gram.det()),
            "half_scaled_odd_unimodular": True,
            "target_lattice": "2*I_{1,10}",
            "target_half_scaled_gram": _matrix_payload(target_gram),
            "explicit_isometry_convention": (
                "Q^T * (G_fixed/2) * Q = diag(1,-1,...,-1)"
            ),
            "explicit_isometry_matrix": _matrix_payload(isometry),
            "explicit_isometry_determinant": int(isometry.det()),
            "integrally_isometric": True,
            "orthogonal_group_unchanged_by_scalar_2": True,
        },
        "primitive_mirrors": {
            "facet_count": len(mirror_records),
            "fixed_lattice_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(fixed_types.items())
            },
            "half_scaled_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(half_types.items())
            },
            "presentation_histogram": {
                f"presented_{presented}__primitive_{square}__scale_{scale}": count
                for (presented, square, scale), count in sorted(
                    presentation_types.items()
                )
            },
            "all_reflections_integral_involutive_isometries": True,
            "target_root_matrix_rank": target_root_matrix.rank(),
            "target_root_span_smith_diagonal": root_span_diagonal,
            "roots_generate_entire_fixed_lattice": True,
            "records": mirror_records,
        },
        "known_maximal_reflection_chamber": {
            "name": "Vinberg P10",
            "lattice": "I_{1,10}, equivalently 2*I_{1,10}",
            "dimension": 10,
            "facet_count": 12,
            "finite_vertex_count": 17,
            "ideal_vertex_count": 2,
            "explicit_roots_in_I_1_10_coordinates": _matrix_payload(p10_roots.T),
            "root_gram_in_2_I_1_10": _matrix_payload(p10_gram),
            "primitive_square_divisibility_histogram": {
                f"{square},{divisibility}": count
                for (square, divisibility), count in sorted(p10_types.items())
            },
            "all_reflections_integral": p10_reflections_integral,
            "normalized_relation_matrix": p10_relations,
            "coxiter": {
                "relative_path": P10_COXITER_RELATIVE.as_posix(),
                "sha256": _sha256(p10_path),
                "p10_root_to_coxiter_zero_based_permutation": (
                    coxiter_p10_permutation
                ),
                "entrywise_relation_matrix_match": True,
            },
            "alexeev_11_11_1": {
                "lattice_name": "I_{1,10}(2)",
                "root_count": 12,
                "isotropic_orbit_count": 2,
                "matrix_sha256": _matrix_sha256(alexeev_matrix),
                "p10_root_to_alexeev_zero_based_permutation": (
                    alexeev_permutation
                ),
                "entrywise_weighted_gram_match": True,
            },
            "scharlau_kirschmer_list10_entry1": {
                **sk_entry,
                "gram10_entry1": _matrix_payload(sk_gram_entry1),
                "signed_target_to_sk_zero_based_permutation": (
                    sk_target_permutation
                ),
                "entrywise_lattice_gram_match": True,
            },
            "full_root_reflection_group": True,
        },
        "solomon_index": {
            "method": "complete Solomon sum over every finite special subgroup",
            "connected_elliptic_components": [
                [index + 1 for index in range(30) if mask & (1 << index)]
                for mask in sorted(connected_masks)
            ],
            "connected_component_search_by_order": component_search,
            "candidate_30_wall_group": candidate_euler,
            "vinberg_p10_group": p10_euler,
            "inclusion": "W_30 < W(P10) by the explicit integral root coordinates",
            "index_formula": "chi(W_30) / chi(W(P10))",
            "index": index.numerator,
            "index_prime_factorization": _factorization(index.numerator),
            "hyperbolic_covolume_formula": (
                "vol = -(32*pi^5/945)*chi in H10"
            ),
            "hyperbolic_covolume": "pi^5 * 527/14288400",
            "full_face_vector": full_face_vector,
        },
        "bounded_weighted_diagram_audit": {
            "scope": (
                "pinned CoxIter ordinary corpus and its separate reflection-"
                "doubling archive only"
            ),
            "candidate_relation_matrix": candidate_relations,
            "coxiter_snapshot": corpus_audit,
            "literal_30_facet_record_present": False,
            "does_not_prove_worldwide_absence": True,
        },
        "lineage_decision": {
            "new_lattice": False,
            "new_maximal_reflection_group": False,
            "new_arithmetic_commensurability_class": False,
            "positive_algebraic_construction": True,
            "literal_30_wall_subgroup_chamber_priority": (
                "open after bounded corpus exclusion"
            ),
            "worldwide_priority_claimed": False,
            "novelty_boundary": (
                "at most the divisor-involution realization and possibly the "
                "literal index-12142080 subgroup chamber"
            ),
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coxiter-corpus", type=Path, default=COXITER_CORPUS)
    parser.add_argument("--coxiter-archive", type=Path, default=COXITER_ARCHIVE)
    parser.add_argument("--sk-directory", type=Path, default=SK_DIRECTORY)
    parser.add_argument("--sk-archive", type=Path, default=SK_ARCHIVE)
    parser.add_argument("--alexeev-directory", type=Path, default=ALEXEEV_DIRECTORY)
    parser.add_argument("--alexeev-archive", type=Path, default=ALEXEEV_ARCHIVE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate(
        arguments.coxiter_corpus,
        arguments.coxiter_archive,
        arguments.sk_directory,
        arguments.sk_archive,
        arguments.alexeev_directory,
        arguments.alexeev_archive,
    )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if (
            not arguments.output.is_file()
            or arguments.output.read_text(encoding="utf-8") != rendered
        ):
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(arguments.output)


if __name__ == "__main__":
    main()
