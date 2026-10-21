#!/usr/bin/env python3
"""Certify the Coxeter rows in D'Urso's Cremona-nef fundamental domains.

The imported algebraic-geometric theorem is D'Urso's description of a
fundamental cone for the Cremona action on the K-negative nef-effective cone
of the blow-up of P2 at n very general points.  This script independently
checks the finite hyperbolic part of the statement:

* the physical divisor normals in I_(1,n) and their exact Gram matrices;
* the complete finite/ideal/hyperideal ray census for n=10,...,14;
* irredundancy and finite volume;
* the exact criterion that the terminal angle is Coxeter iff n=10,11,13;
* integrality (or nonintegrality) of every divisor-lattice reflection; and
* exact labelled-graph lineage against frozen CoxIter records.

The H13 row is Tumarkin's first 13-dimensional n+2 pyramid, not Vinberg's
odd-unimodular P13 chamber.  Their different ray splits give an independent
nonisometry witness.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Iterable

import networkx as nx
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.certify_other_algebraic_reservoirs import (
    _extreme_ray_census,
    _inertia_symmetric,
    _matrix_payload,
)


DEFAULT_OUTPUT = (
    ROOT / "data" / "enumerations" / "durso_cremona_nef_coxeter_transfer.json"
)
DEFAULT_COXITER_ROOT = Path("/private/tmp/coxiter-current/CoxIter-master")

SOURCE = {
    "author": "Luize D'Urso",
    "title": "On the Nef cones of blowups of the projective plane",
    "arxiv": "https://arxiv.org/abs/2412.15460",
    "version": "v2",
    "version_date": "2025-07-09",
    "pdf_sha256": "f5b48d32f67c789a843e98aed1f05664cb35ff2526b8ef2606d41946c9399b39",
    "used_results": ["Theorem 1.6", "Proposition 6.1", "Theorem 6.2", "Proposition 6.3"],
}

COXITER_COMMIT = "b800be48d8240fbabd2d1c9542e61054584ad7fe"


def _edge_key(left: int, right: int, relation: int) -> tuple[int, int, int]:
    return min(left, right), max(left, right), relation


FROZEN_LINEAGE: dict[int, dict[str, Any]] = {
    10: {
        "name": "Vinberg P10 / odd-unimodular full chamber",
        "relative_path": "graphs/group units qf/10-1,1,1,1,1,1,1,1,1,1,1.coxiter",
        "sha256": "a8e8618529200271e4be249442ec283c65f222372ece8d1910476433f57ce43e",
        "edges": [
            *[_edge_key(i, i + 1, 3) for i in range(1, 9)],
            _edge_key(9, 10, 4),
            _edge_key(3, 11, 3),
            _edge_key(10, 12, 0),
        ],
    },
    11: {
        "name": "Vinberg P11 / odd-unimodular full chamber",
        "relative_path": "graphs/group units qf/11-1,1,1,1,1,1,1,1,1,1,1,1.coxiter",
        "sha256": "fbec0cc4d22509b42f8aa2b692b3cd568745c8f869b05f564f30b21a9372534b",
        "edges": [
            *[_edge_key(i, i + 1, 3) for i in range(1, 10)],
            _edge_key(10, 11, 4),
            _edge_key(3, 12, 3),
            _edge_key(11, 13, 4),
        ],
    },
    13: {
        "name": "Tumarkin 13-tum04_13_01",
        "relative_path": "graphs/Tumarkin pyramids n+2/13-tum04_13_01.coxiter",
        "sha256": "d1de070e5aa7cc90292f21424625ee8ef7dfbe718d004fb8434cc8beee0b35c8",
        "edges": [
            *[_edge_key(i, i + 1, 3) for i in range(1, 12)],
            _edge_key(12, 13, 4),
            _edge_key(13, 14, 3),
            _edge_key(3, 15, 3),
        ],
    },
}

FROZEN_PROJECT_P13 = {
    "name": "Vinberg P13 / Tumarkin 13-tum04_13_02",
    "relative_path": "graphs/Tumarkin pyramids n+2/13-tum04_13_02.coxiter",
    "sha256": "871d05dd8911ee8f2ef3ca305753e3dcc1e7c7d1489f9a2fa005f45e6ce5b8d1",
    "edges": [
        *[_edge_key(i, i + 1, 3) for i in range(1, 12)],
        _edge_key(12, 13, 4),
        _edge_key(11, 14, 3),
        _edge_key(3, 15, 3),
    ],
    "ray_split": [43, 3, 0],
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _graph(node_count: int, edges: Iterable[tuple[int, int, int]]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(1, node_count + 1))
    for left, right, relation in edges:
        graph.add_edge(left, right, relation=int(relation))
    return graph


def _isomorphic(
    left_nodes: int,
    left_edges: Iterable[tuple[int, int, int]],
    right_nodes: int,
    right_edges: Iterable[tuple[int, int, int]],
) -> bool:
    if left_nodes != right_nodes:
        return False
    matcher = nx.algorithms.isomorphism.categorical_edge_match("relation", None)
    return nx.is_isomorphic(
        _graph(left_nodes, left_edges),
        _graph(right_nodes, right_edges),
        edge_match=matcher,
    )


def _parse_coxiter(path: Path) -> tuple[int, int, list[tuple[int, int, int]]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    lines = [line for line in lines if line and not line.startswith("vertices labels:")]
    facets, dimension = (int(value) for value in lines[0].split()[:2])
    edges = [_edge_key(*(int(value) for value in line.split()[:3])) for line in lines[1:]]
    return facets, dimension, sorted(edges)


def _verify_external_pins(coxiter_root: Path | None, source_pdf: Path | None) -> None:
    if source_pdf is not None:
        if _sha256(source_pdf) != SOURCE["pdf_sha256"]:
            raise AssertionError("D'Urso source PDF SHA-256 changed")
    if coxiter_root is None:
        return
    records = [*FROZEN_LINEAGE.values(), FROZEN_PROJECT_P13]
    for record in records:
        path = coxiter_root / record["relative_path"]
        if _sha256(path) != record["sha256"]:
            raise AssertionError(f"CoxIter source SHA-256 changed: {path}")
        facets, dimension, edges = _parse_coxiter(path)
        expected_dimension = 13 if "13-tum" in record["relative_path"] else int(
            record["relative_path"].split("/")[-1].split("-")[0]
        )
        if facets != expected_dimension + 2 or dimension != expected_dimension:
            raise AssertionError(f"CoxIter header changed: {path}")
        if edges != sorted(record["edges"]):
            raise AssertionError(f"CoxIter graph changed: {path}")


def _primitive_rays(ray_census: dict[str, Any]) -> list[sp.Matrix]:
    return [sp.Matrix(row["vector"]) for row in ray_census["primitive_rays"]]


def _reflection(intersection: sp.Matrix, root: sp.Matrix) -> sp.Matrix:
    square = (root.T * intersection * root)[0]
    return sp.simplify(sp.eye(intersection.rows) - sp.Rational(2, square) * root * root.T * intersection)


def _is_integral_matrix(matrix: sp.Matrix) -> bool:
    return all(sp.denom(value) == 1 for value in matrix)


def _candidate_edges(n: int) -> list[tuple[int, int, int]]:
    if n not in (10, 11, 13):
        raise ValueError("only the three Coxeter rows have a relation graph")
    terminal_relation = {10: 0, 11: 4, 13: 3}[n]
    return sorted(
        [
            *[_edge_key(i, i + 1, 3) for i in range(1, n - 1)],
            _edge_key(n - 1, n, 4),
            _edge_key(3, n + 1, 3),
            _edge_key(n, n + 2, terminal_relation),
        ]
    )


def _package(n: int) -> dict[str, Any]:
    if n < 10:
        raise ValueError("D'Urso's K-negative package starts at n=10")
    intersection = sp.diag(1, *([-1] * n))
    labels: list[str] = []
    roots: list[list[int]] = []

    for index in range(1, n):
        root = [0] * (n + 1)
        root[index] = 1
        root[index + 1] = -1
        roots.append(root)
        labels.append(f"E{index}-E{index + 1}")

    last = [0] * (n + 1)
    last[n] = 1
    roots.append(last)
    labels.append(f"E{n}")

    roots.append([1, -1, -1, -1] + [0] * (n - 3))
    labels.append("H-E1-E2-E3")

    roots.append([3] + [-1] * n)
    labels.append("-K=3H-E1-...-En")

    columns = sp.Matrix(roots).T
    gram = columns.T * intersection * columns
    ray_census = _extreme_ray_census(intersection, columns)
    rays = _primitive_rays(ray_census)
    inequality_matrix = columns.T * intersection

    facet_ranks = []
    for wall in range(columns.cols):
        on_wall = [ray for ray in rays if (inequality_matrix[wall, :] * ray)[0] == 0]
        facet_ranks.append(sp.Matrix.hstack(*on_wall).rank() if on_wall else 0)

    reflection_rows = []
    for index, label in enumerate(labels):
        reflection = _reflection(intersection, columns[:, index])
        reflection_rows.append(
            {
                "wall": label,
                "square": int(gram[index, index]),
                "integral_on_NS": _is_integral_matrix(reflection),
                "determinant": int(reflection.det()),
                "isometry": reflection.T * intersection * reflection == intersection,
            }
        )

    terminal_cosine_squared = sp.Rational(1, n - 9)
    terminal_relation: int | str | None
    if n == 10:
        terminal_relation = "infinity"
    elif n == 11:
        terminal_relation = 4
    elif n == 13:
        terminal_relation = 3
    else:
        terminal_relation = None

    expected = {
        "total": 9 * n - 71,
        "finite": 9 * n - 73,
        "ideal": 2,
        "hyperideal": 0,
    }
    span_basis = columns[:, : n + 1]
    checks = {
        "ambient_signature": _inertia_symmetric(intersection) == (1, n, 0),
        "wall_gram_inertia": _inertia_symmetric(gram) == (1, n, 1),
        "wall_rank": columns.rank() == n + 1,
        "physical_walls_span_NS_index_one": abs(int(span_basis.det())) == 1,
        "all_facets_irredundant": facet_ranks == [n] * (n + 2),
        "exact_ray_formula": all(ray_census[key] == value for key, value in expected.items()),
        "finite_volume_from_no_hyperideal_rays": ray_census["hyperideal"] == 0,
        "terminal_pairing_exact": (
            gram[n - 1, n + 1] == 1
            and gram[n - 1, n - 1] == -1
            and gram[n + 1, n + 1] == -(n - 9)
        ),
        "canonical_reflection_integral_iff_n_minus_9_divides_2": (
            reflection_rows[-1]["integral_on_NS"] == ((n - 9) in (1, 2))
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"D'Urso package H{n} changed: {checks}")

    result: dict[str, Any] = {
        "dimension": n,
        "facet_count": n + 2,
        "surface": f"blow-up S_{n} of P2 at {n} very general points",
        "algebraic_cone": "Nef^e(S_n) intersect {D.K_S <= 0}, modulo the Cremona action",
        "labels": labels,
        "physical_divisor_columns": _matrix_payload(columns),
        "physical_gram": _matrix_payload(gram),
        "wall_norm_histogram": {
            str(key): value for key, value in sorted(Counter(int(gram[i, i]) for i in range(gram.rows)).items())
        },
        "terminal_cosine_squared": str(terminal_cosine_squared),
        "terminal_relation": terminal_relation,
        "is_coxeter": n in (10, 11, 13),
        "facet_support_ranks": facet_ranks,
        "ray_census": ray_census,
        "reflection_rows": reflection_rows,
        "checks": checks,
        "all_checks_pass": True,
    }
    if n in FROZEN_LINEAGE:
        edges = _candidate_edges(n)
        lineage = FROZEN_LINEAGE[n]
        matches = _isomorphic(n + 2, edges, n + 2, lineage["edges"])
        if not matches:
            raise ArithmeticError(f"H{n} no longer matches its frozen lineage graph")
        result["coxeter_edges"] = [list(edge) for edge in edges]
        result["lineage"] = {
            "name": lineage["name"],
            "coxiter_relative_path": lineage["relative_path"],
            "coxiter_sha256": lineage["sha256"],
            "exact_weighted_graph_isomorphism": True,
        }
    return result


def build_payload() -> dict[str, Any]:
    packages = {str(n): _package(n) for n in range(10, 15)}
    durso_h13 = packages["13"]
    durso_edges = [tuple(edge) for edge in durso_h13["coxeter_edges"]]
    distinct_from_project_p13 = not _isomorphic(
        15, durso_edges, 15, FROZEN_PROJECT_P13["edges"]
    )
    durso_split = [
        durso_h13["ray_census"]["finite"],
        durso_h13["ray_census"]["ideal"],
        durso_h13["ray_census"]["hyperideal"],
    ]
    if not (
        distinct_from_project_p13
        and durso_split == [44, 2, 0]
        and FROZEN_PROJECT_P13["ray_split"] == [43, 3, 0]
    ):
        raise ArithmeticError("the two H13/15 lineages are no longer separated")

    return {
        "schema": "durso_cremona_nef_coxeter_transfer_v1",
        "source": SOURCE,
        "source_theorem_import": {
            "statement": (
                "For every n>=10 the displayed cone is a finite-volume fundamental "
                "domain for the Cremona action on the K-negative nef-effective cone; "
                "it is Coxeter exactly for n=10,11,13."
            ),
            "independently_replayed_here": False,
        },
        "exact_arithmetic_replay": {
            "tested_dimensions": [10, 11, 12, 13, 14],
            "packages": packages,
            "general_terminal_angle": "cos(theta_n)^2 = 1/(n-9)",
            "all_n_ge_10_coxeter_criterion": [10, 11, 13],
            "criterion_proof_partition": {
                "n=10": "theta=0 (parallel)",
                "n=11": "theta=pi/4",
                "n=12": "pi/4 < theta < pi/3, hence theta is not pi/m",
                "n=13": "theta=pi/3",
                "n>=14": "pi/3 < theta < pi/2, hence theta is not pi/m",
            },
        },
        "h13_lineage_separation": {
            "durso": {
                "name": FROZEN_LINEAGE[13]["name"],
                "ray_split": durso_split,
            },
            # Keep the in-memory payload JSON-native as well as the archived
            # serialization.  Otherwise tuple-valued frozen edges compare
            # unequal to their list-valued JSON replay despite encoding the
            # same labelled graph.
            "project_existing": {
                **{
                    key: value
                    for key, value in FROZEN_PROJECT_P13.items()
                    if key != "edges"
                },
                "edges": [list(edge) for edge in FROZEN_PROJECT_P13["edges"]],
            },
            "weighted_graphs_nonisomorphic": distinct_from_project_p13,
            "ray_splits_different": True,
        },
        "coxiter_snapshot": {
            "repository": "https://github.com/rgugliel/CoxIter",
            "commit": COXITER_COMMIT,
            "records": {
                str(n): {
                    "relative_path": row["relative_path"],
                    "sha256": row["sha256"],
                }
                for n, row in FROZEN_LINEAGE.items()
            },
        },
        "claim_boundary": {
            "new_abstract_polytopes": False,
            "new_project_lineage_interpretation": True,
            "h13_abstract_source": "Tumarkin 2004, already in the pinned CoxIter corpus",
            "full_nef_cone_claim": False,
            "precise_cone_claim": "Cremona fundamental sector of the K-negative nef-effective cone",
            "reflection_integral_on_NS_in_H13": False,
            "strict_inverse_L4": False,
            "reason_not_L4": (
                "the surface theorem supplies the cone, while the Coxeter angle test "
                "itself selects 10,11,13; no independent AG classification pushes "
                "back to a new theorem on an a-priori polytope class"
            ),
        },
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--verify-corpus", action="store_true")
    parser.add_argument("--coxiter-root", type=Path, default=DEFAULT_COXITER_ROOT)
    parser.add_argument("--source-pdf", type=Path)
    arguments = parser.parse_args()

    _verify_external_pins(
        arguments.coxiter_root if arguments.verify_corpus else None,
        arguments.source_pdf,
    )
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")

    summary = {
        "coxeter_dimensions": payload["exact_arithmetic_replay"]["all_n_ge_10_coxeter_criterion"],
        "ray_splits": {
            n: [
                row["ray_census"]["finite"],
                row["ray_census"]["ideal"],
                row["ray_census"]["hyperideal"],
            ]
            for n, row in payload["exact_arithmetic_replay"]["packages"].items()
            if row["is_coxeter"]
        },
        "h13_lineage": payload["h13_lineage_separation"]["durso"]["name"],
        "all_checks_pass": payload["all_checks_pass"],
    }
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
