#!/usr/bin/env python3
"""Independent hostile audit of two K3-carried H10 redoubling rays.

The construction files are deliberately not imported.  This verifier starts
from the two frozen 13-wall lattice-root archives, reconstructs finite path
unions exactly, proves the closed root recurrences used for the all-N formulas,
and compares the first descendants with the pinned CoxIter snapshots.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import tarfile
from typing import Any, Iterable, Sequence

import networkx as nx
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

from coxeter_polytopes.lattices import (
    block_diagonal,
    e8_cartan_gram,
    hyperbolic_plane,
)


ROOT = Path(__file__).resolve().parents[1]
F_SOURCE = ROOT / "data/enumerations/k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.json"
F_COMPILER = ROOT / "data/enumerations/ue8m2a1_common_k3_moduli_carrier.json"
S_SOURCE = ROOT / "data/enumerations/jacobian_kummer_h10_maximal_reflection_lineage.json"
S_COMPILER = ROOT / "data/enumerations/ue8m2a4_common_k3_moduli_carrier.json"
OUTPUT = ROOT / "data/enumerations/h10_two_k3_redoubling_rays_hostile_audit.json"

CORPUS = Path("/private/tmp/coxiter-current/CoxIter-master")
CORPUS_ARCHIVE = Path("/private/tmp/coxiter-current.zip")
DOUBLING_RELATIVE = Path("graphs/testing_index2.tar.gz")

SCHEMA = "coxeter.h10.two_k3_redoubling_rays.hostile_audit.v1"
EXPECTED_SOURCE_SHA256 = {
    "F_seed": "93cec35a3c395482f4f99be827b3f5a6ba44764fc7defa22bf23eb6ee951acf0",
    "F_compiler": "1a2b7f8e5a2869d6a0e23b5079c8551b526453ac4764f2a52573909722ffe28e",
    "S_seed": "3cb8b0e9120d13479f9b9e3a32cfb63922e03f41f2f0fc1dc501b0a9f6654c29",
    "S_compiler": "99c02da329d34974bc3e6315d23f6d6a7cfee3854f874b4bdfa86159db8e5d64",
}
CORPUS_ARCHIVE_SHA256 = "ae04648b291c1e5f3d0e47d666c846dc805b244469a20486bca59674b4d171c3"
DOUBLING_ARCHIVE_SHA256 = "3330ca41fb3e047123ace718f08542bc1dd0019b35922d4fa836a217647e7036"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_sha256(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _pair(left: sp.Matrix, gram: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.expand((left.T * gram * right)[0])


def _reflect_raw(vector: sp.Matrix, root: sp.Matrix, gram: sp.Matrix) -> sp.Matrix:
    return (
        vector
        - 2 * _pair(vector, gram, root) / _pair(root, gram, root) * root
    ).applyfunc(sp.expand)


def _reflect(vector: sp.Matrix, root: sp.Matrix, gram: sp.Matrix) -> sp.Matrix:
    result = _reflect_raw(vector, root, gram)
    if any(entry.is_Integer is not True for entry in result):
        raise ArithmeticError("selected reflection left the lattice")
    return result


def _transform(
    word: tuple[int, ...],
    vector: sp.Matrix,
    selected: tuple[sp.Matrix, sp.Matrix],
    gram: sp.Matrix,
) -> sp.Matrix:
    for letter in reversed(word):
        vector = _reflect(vector, selected[letter], gram)
    return vector


def _oriented_key(vector: sp.Matrix) -> tuple[int, ...]:
    values = tuple(int(entry) for entry in vector)
    for value in values:
        if value:
            return values if value > 0 else tuple(-entry for entry in values)
    raise ArithmeticError("zero root")


def _alternating_path(copy_count: int, start: int) -> tuple[tuple[int, ...], ...]:
    words: list[tuple[int, ...]] = [()]
    word: tuple[int, ...] = ()
    for position in range(copy_count - 1):
        word += ((start + position) % 2,)
        words.append(word)
    return tuple(words)


def _path_boundary(
    roots: tuple[sp.Matrix, ...],
    gram: sp.Matrix,
    selected_indices: tuple[int, int],
    copy_count: int,
    start: int,
) -> tuple[sp.Matrix, ...]:
    selected = tuple(roots[index] for index in selected_indices)
    subtree = _alternating_path(copy_count, start)
    internal: set[tuple[int, ...]] = set()
    for word in subtree[1:]:
        prefix, letter = word[:-1], word[-1]
        internal.add(
            _oriented_key(
                _transform(prefix, selected[letter], selected, gram)
            )
        )
    boundary: dict[tuple[int, ...], sp.Matrix] = {}
    for word in subtree:
        for root in roots:
            image = _transform(word, root, selected, gram)
            key = _oriented_key(image)
            if key not in internal:
                boundary.setdefault(key, image)
    return tuple(boundary.values())


def _label(left: sp.Matrix, right: sp.Matrix, gram: sp.Matrix) -> str:
    cross = _pair(left, gram, right)
    if cross < 0:
        raise ArithmeticError(f"non-acute oriented wall pair: {cross}")
    if cross == 0:
        return "2"
    ratio = sp.factor(
        cross**2
        / (_pair(left, gram, left) * _pair(right, gram, right))
    )
    labels = {
        sp.Rational(1, 4): "3",
        sp.Rational(1, 2): "4",
        sp.Rational(3, 4): "6",
        sp.Rational(1): "0",
    }
    if ratio in labels:
        return labels[ratio]
    if ratio > 1:
        return "1"
    raise ArithmeticError(f"unsupported Coxeter ratio {ratio}")


def _relation_matrix(
    roots: Sequence[sp.Matrix], gram: sp.Matrix
) -> list[list[str]]:
    return [
        [
            "1" if left == right else _label(alpha, beta, gram)
            for right, beta in enumerate(roots)
        ]
        for left, alpha in enumerate(roots)
    ]


def _relation_histogram(matrix: Sequence[Sequence[str]]) -> dict[str, int]:
    counts = Counter(
        matrix[left][right]
        for right in range(len(matrix))
        for left in range(right)
    )
    return {label: counts[label] for label in sorted(counts)}


def _node_profiles(matrix: Sequence[Sequence[str]]) -> list[dict[str, Any]]:
    profiles = Counter(
        tuple(
            sorted(
                Counter(
                    matrix[index][other]
                    for other in range(len(matrix))
                    if other != index
                ).items()
            )
        )
        for index in range(len(matrix))
    )
    return [
        {"multiplicity": multiplicity, "labels": dict(profile)}
        for profile, multiplicity in sorted(profiles.items())
    ]


def _graph(matrix: Sequence[Sequence[str]]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(matrix)))
    for right in range(len(matrix)):
        for left in range(right):
            graph.add_edge(left, right, label=matrix[left][right])
    return graph


EDGE_MATCH = nx.algorithms.isomorphism.categorical_edge_match("label", None)


def _isomorphic(
    left: Sequence[Sequence[str]], right: Sequence[Sequence[str]]
) -> bool:
    return nx.is_isomorphic(_graph(left), _graph(right), edge_match=EDGE_MATCH)


def _content(vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(entry)) for entry in vector if entry))


def _divisibility(gram: sp.Matrix, vector: sp.Matrix) -> int:
    return math.gcd(*(abs(int(entry)) for entry in gram * vector))


def _mirror_histogram(
    roots: Sequence[sp.Matrix], gram: sp.Matrix
) -> dict[str, int]:
    counts: Counter[tuple[int, int]] = Counter()
    for root in roots:
        square = int(_pair(root, gram, root))
        divisibility = _divisibility(gram, root)
        counts[(square, divisibility)] += 1
        reflection = (
            sp.eye(gram.rows)
            - sp.Rational(2, square) * root * (root.T * gram)
        )
        if not (
            _content(root) == 1
            and all(entry.is_Integer for entry in reflection)
            and reflection.det() == -1
            and reflection.T * gram * reflection == gram
        ):
            raise ArithmeticError("a boundary wall failed its lattice gate")
    return {
        f"square_{square}_div_{divisibility}": count
        for (square, divisibility), count in sorted(counts.items())
    }


def _doubling_walls(matrix: Sequence[Sequence[str]]) -> list[int]:
    result = []
    for index, row in enumerate(matrix):
        if all(
            value in {"0", "1", "2"}
            or (value.isdigit() and int(value) % 2 == 0)
            for value in row
        ):
            result.append(index)
    return result


def _gluing_cost(matrix: Sequence[Sequence[str]], wall: int) -> int:
    return 2 + sum(
        matrix[wall][other] == "2"
        for other in range(len(matrix))
        if other != wall
    )


def _span_data(roots: Sequence[sp.Matrix]) -> tuple[int, list[int], int | None]:
    matrix = sp.Matrix.hstack(*roots).T
    smith = smith_normal_form(matrix, domain=ZZ)
    diagonal = [
        abs(int(smith[index, index]))
        for index in range(min(smith.shape))
        if smith[index, index]
    ]
    rank = matrix.rank()
    index = math.prod(diagonal) if rank == matrix.cols else None
    return rank, diagonal, index


def _s_gram() -> sp.Matrix:
    e8 = [[-2 * entry for entry in row] for row in e8_cartan_gram()]
    return sp.Matrix(block_diagonal(hyperbolic_plane(), e8, [[-4]]))


def _expected_relations(ray: str, copies: int, start: int) -> dict[str, int]:
    ultraparallel = (copies - 1) * (copies + 2) // 2
    if ray == "F":
        if copies % 2:
            return {
                "0": (copies + 3) // 2,
                "1": ultraparallel,
                "2": (19 * copies + 111) // 2,
                "3": copies + 8,
                "4": 2,
            }
        if start == 0:
            return {
                "0": (copies + 4) // 2,
                "1": ultraparallel,
                "2": (19 * copies + 112) // 2,
                "3": copies + 8,
                "4": 1,
            }
        return {
            "0": (copies + 2) // 2,
            "1": ultraparallel,
            "2": (19 * copies + 110) // 2,
            "3": copies + 8,
            "4": 3,
        }
    if copies % 2:
        return {
            "0": (copies + 3) // 2,
            "1": ultraparallel,
            "2": (19 * copies + 111) // 2,
            "3": (copies + 17) // 2,
            "4": (copies + 3) // 2,
        }
    if start == 0:
        return {
            "0": (copies + 2) // 2,
            "1": ultraparallel,
            "2": (19 * copies + 110) // 2,
            "3": (copies + 18) // 2,
            "4": (copies + 4) // 2,
        }
    return {
        "0": (copies + 4) // 2,
        "1": ultraparallel,
        "2": (19 * copies + 112) // 2,
        "3": (copies + 16) // 2,
        "4": (copies + 2) // 2,
    }


def _expected_mirrors(ray: str, copies: int, start: int) -> dict[str, int]:
    if ray == "F":
        if copies % 2:
            values = {
                (-4, 2): (copies + 17) // 2,
                (-2, 1): (copies + 3) // 2,
                (-2, 2): 2,
            }
        elif start == 0:
            values = {
                (-4, 2): (copies + 18) // 2,
                (-2, 1): (copies + 2) // 2,
                (-2, 2): 2,
            }
        else:
            values = {
                (-4, 2): (copies + 16) // 2,
                (-2, 1): (copies + 4) // 2,
                (-2, 2): 2,
            }
    elif copies % 2:
        values = {
            (-4, 2): (copies + 21) // 2,
            (-4, 4): 1,
            (-2, 1): (copies + 1) // 2,
        }
    elif start == 0:
        values = {
            (-4, 2): copies // 2 + 12,
            (-2, 1): copies // 2,
        }
    else:
        values = {
            (-4, 2): copies // 2 + 9,
            (-4, 4): 2,
            (-2, 1): copies // 2 + 1,
        }
    return {
        f"square_{square}_div_{divisibility}": count
        for (square, divisibility), count in sorted(values.items())
    }


def _closed_boundary(
    ray: str,
    roots: tuple[sp.Matrix, ...],
    selected_indices: tuple[int, int],
    copies: int,
    start: int,
) -> tuple[sp.Matrix, ...]:
    a, b = (roots[index] for index in selected_indices)
    delta = a + b
    k = copies // 2
    e = lambda index: a + index * delta
    if ray == "F":
        fixed = [roots[index] for index in (*range(8), 11)]
        c, d = roots[8], roots[10]
        c_at = lambda index: c + 2 * index * a + 2 * index * (index - 1) * delta
        d_at = lambda index: d + 2 * index * a + 2 * index * index * delta
        if start == 0:
            c_range = range(k + 1)
            d_range = range(k + 1 if copies % 2 else k)
            endpoints = [-e(-1), e(2 * k if copies % 2 else 2 * k - 1)]
        else:
            c_range = range(-k, 1) if copies % 2 else range(-(k - 1), 1)
            d_range = range(-k, 1)
            endpoints = [e(0), -e(-(2 * k + 1) if copies % 2 else -2 * k)]
    else:
        fixed = [roots[index] for index in (*range(1, 9), 11)]
        c, d = roots[0], roots[12]
        c_at = lambda index: c + index * a + index * index * delta
        d_at = lambda index: d + 2 * index * a + 2 * index * (index - 1) * delta
        if start == 0:
            c_range = range(k + 1 if copies % 2 else k)
            d_range = range(k + 1)
            endpoints = [-e(-1), e(2 * k if copies % 2 else 2 * k - 1)]
        else:
            c_range = range(-k, 1)
            d_range = range(-k, 1) if copies % 2 else range(-(k - 1), 1)
            endpoints = [e(0), -e(-(2 * k + 1) if copies % 2 else -2 * k)]
    return tuple(
        [*fixed, *(c_at(index) for index in c_range), *(d_at(index) for index in d_range), *endpoints]
    )


def _symbolic_recurrence_green(
    ray: str,
    roots: tuple[sp.Matrix, ...],
    gram: sp.Matrix,
    selected_indices: tuple[int, int],
) -> bool:
    a, b = (roots[index] for index in selected_indices)
    delta = a + b
    variable = sp.symbols("j", integer=True)
    if ray == "F":
        c, d = roots[8], roots[10]
        c_at = lambda j: c + 2 * j * a + 2 * j * (j - 1) * delta
        d_at = lambda j: d + 2 * j * a + 2 * j * j * delta
    else:
        c, d = roots[0], roots[12]
        c_at = lambda j: c + j * a + j * j * delta
        d_at = lambda j: d + 2 * j * a + 2 * j * (j - 1) * delta
    transform = lambda vector: _reflect_raw(
        _reflect_raw(vector, b, gram), a, gram
    )
    identities = [
        transform(delta) - delta,
        transform(a) - (a + 2 * delta),
        transform(c_at(variable)) - c_at(variable + 1),
        transform(d_at(variable)) - d_at(variable + 1),
    ]
    return all(
        sp.simplify(entry) == 0
        for identity in identities
        for entry in identity
    )


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
        (line for line in lines[1:] if line.startswith("vertices labels:")),
        None,
    )
    labels = (
        vertex_line.split(":", 1)[1].split()
        if vertex_line is not None
        else [str(index + 1) for index in range(facets)]
    )
    index_of = {label: index for index, label in enumerate(labels)}
    if len(index_of) != facets:
        raise ArithmeticError("invalid CoxIter vertex labels")
    for line in lines[1:]:
        if line.startswith("vertices labels:"):
            continue
        left_label, right_label, value = line.split()
        left, right = index_of[left_label], index_of[right_label]
        matrix[left][right] = matrix[right][left] = value
    return facets, dimension, matrix


def _corpus_records(corpus: Path, corpus_archive: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if _sha256(corpus_archive) != CORPUS_ARCHIVE_SHA256:
        raise ArithmeticError("ordinary CoxIter archive changed")
    doubling_path = corpus / DOUBLING_RELATIVE
    if _sha256(doubling_path) != DOUBLING_ARCHIVE_SHA256:
        raise ArithmeticError("CoxIter doubling archive changed")
    records: list[dict[str, Any]] = []
    ordinary_count = 0
    for path in sorted(corpus.rglob("*.coxiter")):
        facets, dimension, matrix = _parse_coxiter_text(
            path.read_text(encoding="utf-8")
        )
        ordinary_count += 1
        if facets in {14, 15, 16}:
            records.append(
                {
                    "lane": "ordinary",
                    "path": path.relative_to(corpus).as_posix(),
                    "facets": facets,
                    "dimension": dimension,
                    "matrix": matrix,
                }
            )
    doubling_count = 0
    with tarfile.open(doubling_path, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"failed to read {member.name}")
            facets, dimension, matrix = _parse_coxiter_text(
                handle.read().decode("utf-8")
            )
            doubling_count += 1
            if facets in {14, 15, 16}:
                records.append(
                    {
                        "lane": "doubling",
                        "path": member.name,
                        "facets": facets,
                        "dimension": dimension,
                        "matrix": matrix,
                    }
                )
    metadata = {
        "repository": "https://github.com/rgugl/CoxIter",
        "commit": "b800be48d8240fbabd2d1c9542e61054584ad7fe",
        "ordinary_archive_sha256": CORPUS_ARCHIVE_SHA256,
        "doubling_archive_sha256": DOUBLING_ARCHIVE_SHA256,
        "ordinary_record_count": ordinary_count,
        "doubling_record_count": doubling_count,
        "same_node_records_by_lane_and_facets": {
            f"{lane}_{facets}": count
            for (lane, facets), count in sorted(
                Counter((record["lane"], record["facets"]) for record in records).items()
            )
        },
    }
    return records, metadata


def _lineage_record(
    identifier: str,
    matrix: list[list[str]],
    records: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    histogram = _relation_histogram(matrix)
    same_nodes = [record for record in records if record["facets"] == len(matrix)]
    histogram_matches = [
        record
        for record in same_nodes
        if _relation_histogram(record["matrix"]) == histogram
    ]
    exact_matches = [
        record
        for record in histogram_matches
        if _isomorphic(matrix, record["matrix"])
    ]
    def compact(record: dict[str, Any]) -> dict[str, Any]:
        return {
            "lane": record["lane"],
            "path": record["path"],
            "declared_dimension": record["dimension"],
        }
    return {
        "id": identifier,
        "facets": len(matrix),
        "relation_histogram": histogram,
        "same_node_record_count_all_dimensions": len(same_nodes),
        "histogram_matches": [compact(record) for record in histogram_matches],
        "exact_0_1_distinguished_labelled_graph_matches": [
            compact(record) for record in exact_matches
        ],
    }


def build_certificate(
    corpus: Path = CORPUS,
    corpus_archive: Path = CORPUS_ARCHIVE,
) -> dict[str, Any]:
    sources = {
        "F_seed": F_SOURCE,
        "F_compiler": F_COMPILER,
        "S_seed": S_SOURCE,
        "S_compiler": S_COMPILER,
    }
    actual_hashes = {key: _sha256(path) for key, path in sources.items()}
    if actual_hashes != EXPECTED_SOURCE_SHA256:
        raise ArithmeticError(f"a pinned source changed: {actual_hashes}")
    f_source = json.loads(F_SOURCE.read_text(encoding="utf-8"))
    f_compiler = json.loads(F_COMPILER.read_text(encoding="utf-8"))
    s_source = json.loads(S_SOURCE.read_text(encoding="utf-8"))
    s_compiler = json.loads(S_COMPILER.read_text(encoding="utf-8"))
    if not all(
        document["all_checks_pass"]
        for document in (f_source, f_compiler, s_source, s_compiler)
    ):
        raise ArithmeticError("an upstream carrier theorem is not green")

    f_gram = sp.Matrix(f_source["full_Picard_lattice"]["standard_geometric_gram"])
    f_roots = tuple(
        sp.Matrix(record["primitive_geometric_standard_F_coordinates"])
        for record in f_source["transported_H10_13_chamber"]["wall_records"]
    )
    s_gram = _s_gram()
    s_roots = tuple(
        sp.Matrix(root)
        for root in s_source["vinberg_maximal_integral_reflection_chamber"]["roots_in_vinberg_order"]
    )
    specifications = {
        "F": {
            "gram": f_gram,
            "roots": f_roots,
            "selected": (9, 12),
            "expected_doubling": [9, 12],
            "span_index": 1,
            "carrier_level": "explicit nine-parameter equation family",
        },
        "S": {
            "gram": s_gram,
            "roots": s_roots,
            "selected": (9, 10),
            "expected_doubling": [0, 9, 10],
            "span_index": 2,
            "carrier_level": "nine-dimensional period locus; no common equation family",
        },
    }

    ray_records: dict[str, Any] = {}
    candidate_matrices: dict[str, list[list[str]]] = {}
    ray_checks: dict[str, bool] = {}
    for ray, specification in specifications.items():
        gram = specification["gram"]
        roots = specification["roots"]
        selected_indices = specification["selected"]
        seed_matrix = _relation_matrix(roots, gram)
        selected = tuple(roots[index] for index in selected_indices)
        delta = selected[0] + selected[1]
        omitted_span = [
            _span_data([root for index, root in enumerate(roots) if index != omit])
            for omit in selected_indices
        ]
        sample_records = []
        sample_green = True
        for copies in range(1, 13):
            starts = (0,) if copies == 1 else (0, 1)
            for start in starts:
                boundary = _path_boundary(
                    roots, gram, selected_indices, copies, start
                )
                closed = _closed_boundary(
                    ray, roots, selected_indices, copies, start
                )
                matrix = _relation_matrix(boundary, gram)
                span_rank, span_smith, span_index = _span_data(boundary)
                record_green = (
                    len(boundary) == copies + 12
                    and {_oriented_key(root) for root in boundary}
                    == {_oriented_key(root) for root in closed}
                    and _relation_histogram(matrix)
                    == {
                        label: count
                        for label, count in _expected_relations(
                            ray, copies, start
                        ).items()
                        if count
                    }
                    and _mirror_histogram(boundary, gram)
                    == _expected_mirrors(ray, copies, start)
                    and span_rank == 11
                    and span_index == specification["span_index"]
                )
                sample_green &= record_green
                if copies in {1, 2, 3, 4, 11, 12}:
                    sample_records.append(
                        {
                            "copies": copies,
                            "start": start,
                            "facets": len(boundary),
                            "relations": _relation_histogram(matrix),
                            "mirrors": _mirror_histogram(boundary, gram),
                            "span_smith": span_smith,
                            "span_index": span_index,
                            "green": record_green,
                        }
                    )
                if copies in {2, 3, 4} and (copies % 2 == 0 or start == 0):
                    suffix = "A" if start == 0 else "B"
                    candidate_matrices[f"{ray}_N{copies}_{suffix}"] = matrix
        odd_phase_green = all(
            _isomorphic(
                _relation_matrix(
                    _path_boundary(roots, gram, selected_indices, copies, 0),
                    gram,
                ),
                _relation_matrix(
                    _path_boundary(roots, gram, selected_indices, copies, 1),
                    gram,
                ),
            )
            for copies in range(3, 12, 2)
        )
        even_phase_distinct = all(
            _expected_relations(ray, copies, 0)
            != _expected_relations(ray, copies, 1)
            for copies in range(2, 13, 2)
        )
        checks = {
            "seed_has_13_facets": len(roots) == 13,
            "doubling_wall_inventory_exact": _doubling_walls(seed_matrix)
            == specification["expected_doubling"],
            "selected_walls_are_distinct_parallel": (
                selected[0] != selected[1]
                and _pair(selected[0], gram, selected[1]) ** 2
                == _pair(selected[0], gram, selected[0])
                * _pair(selected[1], gram, selected[1])
                and delta != sp.zeros(11, 1)
                and _pair(delta, gram, delta) == 0
                and all(_pair(delta, gram, root) == 0 for root in selected)
            ),
            "both_gluing_costs_are_12": [
                _gluing_cost(seed_matrix, index) for index in selected_indices
            ]
            == [12, 12],
            "symbolic_parabolic_orbit_recurrence": _symbolic_recurrence_green(
                ray, roots, gram, selected_indices
            ),
            "selected_wall_omissions_keep_expected_span": all(
                rank == 11 and index == specification["span_index"]
                for rank, _, index in omitted_span
            ),
            "closed_form_boundary_relations_mirrors_and_spans_N1_to_N12": sample_green,
            "odd_phases_are_congruent": odd_phase_green,
            "even_phases_have_distinct_relation_histograms": even_phase_distinct,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"{ray} ray failed: {checks}")
        ray_checks.update({f"{ray}_{key}": value for key, value in checks.items()})
        ray_records[ray] = {
            "lattice": (
                "U+E8(-2)+<-2>" if ray == "F" else "U+E8(-2)+<-4>"
            ),
            "carrier_level": specification["carrier_level"],
            "selected_doubling_walls_zero_based": list(selected_indices),
            "all_seed_doubling_walls_zero_based": _doubling_walls(seed_matrix),
            "selected_wall_types": [
                {
                    "square": int(_pair(root, gram, root)),
                    "divisibility": _divisibility(gram, root),
                }
                for root in selected
            ],
            "selected_raw_inner_product": int(
                _pair(selected[0], gram, selected[1])
            ),
            "selected_raw_inner_product_square": int(
                _pair(selected[0], gram, selected[1]) ** 2
            ),
            "selected_normalized_square_ratio": "1",
            "common_isotropic_vector": [int(entry) for entry in delta],
            "gluing_costs": [
                _gluing_cost(seed_matrix, index) for index in selected_indices
            ],
            "seed_relations": _relation_histogram(seed_matrix),
            "seed_mirrors": _mirror_histogram(roots, gram),
            "wall_span_index_for_every_path": specification["span_index"],
            "all_N_conclusion": (
                "For every N>=1, an N-vertex path in the infinite-dihedral "
                "Cayley line gives a genuine finite-volume noncompact H10 "
                "Coxeter polyhedron with index N, volume N times the seed, "
                "and N+12 irredundant facets."
            ),
            "sample_formula_replays": sample_records,
            "checks": checks,
        }

    f_seed_matrix = _relation_matrix(f_roots, f_gram)
    s_seed_matrix = _relation_matrix(s_roots, s_gram)
    seed_isomorphic = _isomorphic(f_seed_matrix, s_seed_matrix)
    seed_profiles_distinct = _node_profiles(f_seed_matrix) != _node_profiles(s_seed_matrix)

    corpus_records, corpus_metadata = _corpus_records(corpus, corpus_archive)
    lineage = {
        identifier: _lineage_record(identifier, matrix, corpus_records)
        for identifier, matrix in sorted(candidate_matrices.items())
    }
    expected_exact = {
        "F_N2_A": [],
        "F_N2_B": [],
        "F_N3_A": [],
        "F_N4_A": [],
        "F_N4_B": [],
        "S_N2_A": [
            "testing_index2/10-2,1,1,1,1,1,1,1,1,1,1-12-13_1.coxiter"
        ],
        "S_N2_B": [
            "testing_index2/10-2,1,1,1,1,1,1,1,1,1,1-13-12_1.coxiter"
        ],
        "S_N3_A": [],
        "S_N4_A": [
            "testing_index2/10-2,1,1,1,1,1,1,1,1,1,1-12-13_2.coxiter"
        ],
        "S_N4_B": [
            "testing_index2/10-2,1,1,1,1,1,1,1,1,1,1-13-12_2.coxiter"
        ],
    }
    actual_exact = {
        identifier: [
            record["path"]
            for record in value["exact_0_1_distinguished_labelled_graph_matches"]
        ]
        for identifier, value in lineage.items()
    }
    f_histogram_survivor = lineage["F_N2_B"]["histogram_matches"]
    collision_profiles_distinct = (
        _node_profiles(candidate_matrices["F_N2_B"])
        != _node_profiles(candidate_matrices["S_N2_A"])
    )
    lineage_checks = {
        "ordinary_and_doubling_corpus_counts_1214_2959": (
            corpus_metadata["ordinary_record_count"],
            corpus_metadata["doubling_record_count"],
        )
        == (1214, 2959),
        "exact_match_inventory_frozen": actual_exact == expected_exact,
        "F_N2_B_has_one_histogram_survivor_but_no_exact_match": (
            [record["path"] for record in f_histogram_survivor]
            == expected_exact["S_N2_A"]
            and not lineage["F_N2_B"][
                "exact_0_1_distinguished_labelled_graph_matches"
            ]
        ),
        "node_profiles_separate_the_F_S_N2_collision": collision_profiles_distinct,
        "F_and_S_seeds_are_not_relation_diagram_isomorphic": not seed_isomorphic,
        "seed_node_profiles_already_separate_them": seed_profiles_distinct,
        "S_seed_is_classical_qf2": (
            s_source["lineage_decision"]["maximal_chamber_is_known_qf2_chamber"]
            is True
        ),
    }
    if not all(lineage_checks.values()):
        raise ArithmeticError(f"lineage audit failed: {lineage_checks}")

    checks = {**ray_checks, **lineage_checks}
    result: dict[str, Any] = {
        "schema": SCHEMA,
        "status_date": "2026-08-08",
        "theorem_transfer": {
            "source": (
                "Daniel Allcock, Infinitely many hyperbolic Coxeter groups "
                "through dimension 19, Geom. Topol. 10 (2006), Lemma 1.3, "
                "Theorem 4.1, and Corollary 4.2"
            ),
            "doi": "https://doi.org/10.2140/gt.2006.10.737",
            "hypothesis_audited": (
                "A doubling wall has only even finite Coxeter labels.  Two "
                "selected doubling walls need only be disjoint in H^10; "
                "parallel walls qualify even though their closures meet ideally."
            ),
            "application": (
                "The selected reflections generate the infinite dihedral "
                "group. Every finite N-vertex subtree is a path; its union "
                "is a Coxeter fundamental polyhedron of a reflection subgroup "
                "of index N. Interiors are disjoint, so volume scales by N."
            ),
        },
        "sources": {
            key: {
                "path": path.relative_to(ROOT).as_posix(),
                "sha256": actual_hashes[key],
            }
            for key, path in sources.items()
        },
        "rays": ray_records,
        "closed_relation_formulas": {
            "semantics": {
                "0": "parallel",
                "1": "generic ultraparallel",
                "2": "m=2",
                "3": "m=3",
                "4": "m=4",
            },
            "F": {
                "odd_N": "0=(N+3)/2, 1=(N-1)(N+2)/2, 2=(19N+111)/2, 3=N+8, 4=2",
                "even_A": "0=(N+4)/2, 1=(N-1)(N+2)/2, 2=(19N+112)/2, 3=N+8, 4=1",
                "even_B": "0=(N+2)/2, 1=(N-1)(N+2)/2, 2=(19N+110)/2, 3=N+8, 4=3",
            },
            "S": {
                "odd_N": "0=(N+3)/2, 1=(N-1)(N+2)/2, 2=(19N+111)/2, 3=(N+17)/2, 4=(N+3)/2",
                "even_A": "0=(N+2)/2, 1=(N-1)(N+2)/2, 2=(19N+110)/2, 3=(N+18)/2, 4=(N+4)/2",
                "even_B": "0=(N+4)/2, 1=(N-1)(N+2)/2, 2=(19N+112)/2, 3=(N+16)/2, 4=(N+2)/2",
            },
        },
        "ray_distinction": {
            "seed_global_histograms_equal": _relation_histogram(f_seed_matrix)
            == _relation_histogram(s_seed_matrix),
            "seed_relation_diagrams_isomorphic": seed_isomorphic,
            "F_seed_node_profiles": _node_profiles(f_seed_matrix),
            "S_seed_node_profiles": _node_profiles(s_seed_matrix),
            "only_descendant_histogram_collision_through_the_closed_formulas": (
                "F even-B at N=2 equals S even-A at N=2; exact node profiles differ"
            ),
        },
        "bounded_lineage": {
            "scope": (
                "All same-node records in the pinned ordinary and index-two-"
                "doubling CoxIter snapshots were scanned, regardless of declared "
                "dimension. This is not a worldwide publication-priority theorem."
            ),
            "corpus": corpus_metadata,
            "candidates": lineage,
            "conclusion": (
                "The S ray is the classical qf2/Allcock-CoxIter redoubling "
                "ray and its N=2,N=4 branches match the pinned descendants. "
                "The F N=2,N=3,N=4 descendants have no exact labelled-diagram "
                "match. F_N2_B requires exact graph comparison because one S "
                "record has the same global label histogram."
            ),
        },
        "carrier_scope": {
            "F": (
                "Every boundary root is primitive integral in F and the walls "
                "span F with index one. The same explicit nine-parameter K3 "
                "equation family is a weak full-Picard divisor carrier."
            ),
            "S": (
                "Every boundary root is primitive integral in S and the wall "
                "span has index two. The same nine-dimensional S-polarized K3 "
                "period locus is a weak full-Picard divisor carrier; no common "
                "equation family is supplied."
            ),
            "excluded_claims": [
                "full nef cone",
                "effective or irreducible representative for every wall",
                "surface automorphisms realizing the reflections",
                "worldwide novelty of the F descendants",
                "abstract novelty of the classical S descendants",
            ],
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    result["payload_sha256"] = _canonical_sha256(result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--corpus-archive", type=Path, default=CORPUS_ARCHIVE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate(arguments.corpus, arguments.corpus_archive)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(arguments.output)
    print(result["payload_sha256"])


if __name__ == "__main__":
    main()
