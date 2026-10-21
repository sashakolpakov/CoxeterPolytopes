#!/usr/bin/env python3
"""Certify the maximal reflection chamber containing JK-H10-0007.

The script independently runs Vinberg's algorithm on the intrinsic lattice
``U+E8(-2)+<-4>``, applies the exact finite-volume stopping test, and compares
the resulting normalized Coxeter matrix against the pinned official CoxIter
corpus.  The unique corpus match is the chamber for
``-2*x0^2+x1^2+...+x10^2``.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import tarfile
from typing import Any, Sequence

import networkx as nx

from coxeter_polytopes.search.jk_h10_maximal_vinberg import (
    divisibility,
    ordered_level_keys,
    pair,
    run_search,
)


ROOT = Path(__file__).resolve().parents[1]
CORPUS = Path("/private/tmp/coxiter-current/CoxIter-master")
ARCHIVE = Path("/private/tmp/coxiter-current.zip")
OUTPUT = ROOT / "data" / "enumerations" / "jacobian_kummer_h10_maximal_reflection_lineage.json"
COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"
EXPECTED_MATCH = "graphs/group units qf/10-2,1,1,1,1,1,1,1,1,1,1.coxiter"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relation(left: Sequence[int], right: Sequence[int]) -> str:
    product = pair(left, right)
    if product == 0:
        return "2"
    ratio = Fraction(product * product, pair(left, left) * pair(right, right))
    labels = {
        Fraction(1, 4): "3",
        Fraction(1, 2): "4",
        Fraction(3, 4): "6",
        Fraction(1): "0",
    }
    if ratio not in labels:
        raise ArithmeticError(f"unsupported normalized product {ratio}")
    return labels[ratio]


def _root_relation_matrix(roots: Sequence[Sequence[int]]) -> list[list[str]]:
    matrix = [["1" if i == j else "2" for j in range(len(roots))] for i in range(len(roots))]
    for i in range(len(roots)):
        for j in range(i + 1, len(roots)):
            matrix[i][j] = matrix[j][i] = _relation(roots[i], roots[j])
    return matrix


def _parse_coxiter(path: Path) -> tuple[int, int, list[list[str]]]:
    lines = [
        content
        for line in path.read_text(encoding="utf-8").splitlines()
        if (content := line.split("#", 1)[0].strip())
    ]
    facets, dimension = map(int, lines[0].split()[:2])
    matrix = [["1" if i == j else "2" for j in range(facets)] for i in range(facets)]
    vertex_line = next((line for line in lines[1:] if line.startswith("vertices labels:")), None)
    if vertex_line is None:
        labels = [str(index + 1) for index in range(facets)]
    else:
        labels = vertex_line.split(":", 1)[1].split()
    if len(labels) != facets or len(set(labels)) != facets:
        raise ValueError(f"invalid CoxIter vertex labels in {path}")
    index_of = {label: index for index, label in enumerate(labels)}
    for line in lines[1:]:
        if line.startswith("vertices labels:"):
            continue
        tokens = line.split()
        if len(tokens) != 3:
            raise ValueError(f"unsupported CoxIter line in {path}: {line}")
        left, right = index_of[tokens[0]], index_of[tokens[1]]
        matrix[left][right] = matrix[right][left] = tokens[2]
    return facets, dimension, matrix


def _graph(matrix: Sequence[Sequence[str]]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(matrix)))
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            graph.add_edge(i, j, relation=matrix[i][j])
    return graph


def _permutation(source: Sequence[Sequence[str]], target: Sequence[Sequence[str]]) -> list[int] | None:
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        _graph(source),
        _graph(target),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("relation", None),
    )
    mapping = next(matcher.isomorphisms_iter(), None)
    if mapping is None:
        return None
    permutation = [mapping[index] for index in range(len(source))]
    if any(source[i][j] != target[permutation[i]][permutation[j]] for i in range(len(source)) for j in range(len(source))):
        raise ArithmeticError("reported graph permutation is not entrywise exact")
    return permutation


def _scan_doubling_archive(path: Path, corpus: Path) -> dict[str, Any]:
    total = 0
    header_histogram: Counter[tuple[int, int]] = Counter()
    h10_84_records: list[str] = []
    with tarfile.open(path, "r:gz") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".coxiter"):
                continue
            handle = archive.extractfile(member)
            if handle is None:
                raise ArithmeticError(f"failed to read {member.name}")
            header = ""
            for raw_line in handle:
                content = raw_line.decode("utf-8").split("#", 1)[0].strip()
                if content:
                    header = content
                    break
            facets, dimension = map(int, header.split()[:2])
            total += 1
            header_histogram[(dimension, facets)] += 1
            if (dimension, facets) == (10, 84):
                h10_84_records.append(member.name)
    return {
        "path": str(path.relative_to(corpus)),
        "sha256": _sha256(path),
        "coxiter_file_count": total,
        "h10_facet_count_histogram": {
            str(facets): count
            for (dimension, facets), count in sorted(header_histogram.items())
            if dimension == 10
        },
        "h10_84_facet_records": h10_84_records,
    }


def build_certificate(corpus: Path, archive: Path | None = None) -> dict[str, Any]:
    search = run_search(3, attempt_stopping=True)
    roots = search.roots
    relation_matrix = _root_relation_matrix(roots)
    candidates: list[dict[str, Any]] = []
    matches: list[dict[str, Any]] = []
    for path in sorted(corpus.rglob("*.coxiter")):
        header = next(
            content
            for line in path.read_text(encoding="utf-8").splitlines()
            if (content := line.split("#", 1)[0].strip())
        )
        header_facets, header_dimension = map(int, header.split()[:2])
        if (header_facets, header_dimension) != (13, 10):
            continue
        facets, dimension, matrix = _parse_coxiter(path)
        if (facets, dimension) != (13, 10):
            raise ArithmeticError("CoxIter header changed during parsing")
        relative = path.relative_to(corpus).as_posix()
        candidates.append({"path": relative})
        permutation = _permutation(relation_matrix, matrix)
        if permutation is not None:
            matches.append(
                {
                    "path": relative,
                    "source_root_index_to_corpus_zero_based_vertex": permutation,
                    "entrywise_relation_matrix_match": True,
                }
            )

    type_counts = Counter((pair(root, root), divisibility(root)) for root in roots)
    doubling_archive = _scan_doubling_archive(
        corpus / "graphs" / "testing_index2.tar.gz", corpus
    )
    accepted_levels = [
        {
            "control": level.control,
            "root_square": -level.norm_abs,
            "candidate_count": level.candidate_count,
            "accepted_count": len(level.accepted_roots),
            "accepted_roots": [list(root) for root in level.accepted_roots],
        }
        for level in search.levels
        if level.accepted_roots
    ]
    checks = {
        "level_order_exhausted_through_height_9_over_2": ordered_level_keys(3)[-1] == (6, 8),
        "wall_count_is_13": len(roots) == 13,
        "finite_volume_stopping_fired": search.stopping.fired,
        "finite_vertex_count_is_29": search.stopping.finite_vertex_count == 29,
        "cusp_count_is_3": search.stopping.cusp_count == 3,
        "no_hyperideal_or_opposite_sheet_rays": search.stopping.hyperideal_ray_count == 0
        and search.stopping.opposite_sheet_ray_count == 0,
        "root_type_histogram": type_counts == Counter({(-4, 2): 11, (-2, 1): 1, (-4, 4): 1}),
        "unique_corpus_match": [match["path"] for match in matches] == [EXPECTED_MATCH],
        "doubling_archive_has_2959_records": doubling_archive["coxiter_file_count"]
        == 2_959,
        "doubling_archive_has_no_h10_84_facet_record": doubling_archive[
            "h10_84_facet_records"
        ]
        == [],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"maximal-reflection lineage checks failed: {checks}")

    return {
        "schema": "jacobian_kummer_h10_maximal_reflection_lineage_v1",
        "candidate": {
            "id": "JK-H10-0007",
            "intrinsic_lattice": "U + E8(-2) + <-4>",
            "induced_chamber_facets": 84,
        },
        "vinberg_maximal_integral_reflection_chamber": {
            "dimension": 10,
            "control_vector": [1, 1, *([0] * 9)],
            "height_zero_type": "A1 + E8 + A1",
            "reflection_norms_exhausted": [-2, -4, -8],
            "height_equation": "d^2 + 4*v_E8^2 + 8*z^2 = n^2 + 2*q",
            "height_cutoff": "9/2",
            "accepted_positive_levels": accepted_levels,
            "roots_in_vinberg_order": [list(root) for root in roots],
            "primitive_root_type_histogram": {
                f"square_{square}_div_{root_div}": count
                for (square, root_div), count in sorted(type_counts.items())
            },
            "normalized_coxeter_relation_matrix": relation_matrix,
            "facet_count": len(roots),
            "finite_vertices": search.stopping.finite_vertex_count,
            "cusps": search.stopping.cusp_count,
            "hyperideal_rays": search.stopping.hyperideal_ray_count,
            "vinberg_stopping_fired": search.stopping.fired,
        },
        "coxiter_snapshot": {
            "repository": "https://github.com/rgugliel/CoxIter",
            "commit": COXITER_COMMIT,
            "archive_sha256": _sha256(archive) if archive is not None else None,
            "h10_13_facet_candidates": candidates,
            "exact_matches": matches,
            "reflection_doubling_archive": doubling_archive,
        },
        "lineage_decision": {
            "maximal_chamber_is_known_qf2_chamber": True,
            "known_quadratic_form": "-2*x0^2 + x1^2 + ... + x10^2",
            "induced_84_wall_group_is_a_finite_index_reflection_subgroup": True,
            "reason_for_finite_index": "all 84 facet reflections are integral reflections of L and both chambers have finite covolume",
            "does_not_settle_whether_the_84_wall_subgroup_chamber_was_printed_before": True,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--archive", type=Path, default=ARCHIVE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate(arguments.corpus, arguments.archive if arguments.archive.is_file() else None)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")
    print(arguments.output)


if __name__ == "__main__":
    main()
