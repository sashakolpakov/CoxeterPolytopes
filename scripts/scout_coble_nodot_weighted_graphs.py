#!/usr/bin/env python3
"""Scout small-edge ``{0,1,2}`` Coble--Mukai Gram diagrams.

This is deliberately a scout rather than a certificate.  It streams
unlabelled connected graphs from nauty ``geng``, promotes up to three of
their edges from root pairing 1 to root pairing 2, and retains matrices
whose negative Cartan matrix has inertia ``(9,1,2)`` and exact rank 10.

The Type-I Coble--Mukai diagram occurs in the degree-at-most-three,
three-promoted-edge lane.  The point of the scout is to measure how much
larger that finite Lorentzian denominator is before the face/volume and
``E10`` saturation gates are imposed.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
import subprocess
from typing import Iterator

import networkx as nx
import numpy as np
import sympy as sp
from flint import fmpz_mat


def _graphs(maximum_degree: int, edge_range: str) -> Iterator[nx.Graph]:
    command = ["geng", "-q", "-c", f"-D{maximum_degree}", "12", edge_range]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=False)
    assert process.stdout is not None
    for line in process.stdout:
        yield nx.from_graph6_bytes(line.strip())
    return_code = process.wait()
    if return_code:
        raise RuntimeError(f"geng exited with status {return_code}")


def _matrix(graph: nx.Graph, promoted: tuple[tuple[int, int], ...]) -> np.ndarray:
    matrix = 2 * np.eye(12, dtype=np.int64)
    for left, right in graph.edges():
        matrix[left, right] = matrix[right, left] = -1
    for left, right in promoted:
        matrix[left, right] = matrix[right, left] = -2
    return matrix


def _numerical_signature(matrix: np.ndarray) -> tuple[int, int, int]:
    eigenvalues = np.linalg.eigvalsh(matrix.astype(np.float64))
    tolerance = 1e-7
    return (
        int(np.count_nonzero(eigenvalues > tolerance)),
        int(np.count_nonzero(eigenvalues < -tolerance)),
        int(np.count_nonzero(np.abs(eigenvalues) <= tolerance)),
    )


def _radical_quotient_determinant(matrix: np.ndarray) -> int:
    """Return the determinant on ``Z^12 / rad`` using a unimodular basis."""

    gram = fmpz_mat(matrix.tolist())
    _, transform = gram.transpose().hnf(transform=True)
    basis = transform.transpose()
    quotient = basis.transpose() * gram * basis
    if any(quotient[i, j] for i in range(12) for j in range(12) if i >= 10 or j >= 10):
        raise ArithmeticError("column HNF failed to isolate the saturated radical")
    block = fmpz_mat(
        10,
        10,
        [quotient[i, j] for i in range(10) for j in range(10)],
    )
    return int(block.det())


def build(
    maximum_degree: int,
    edge_range: str,
    maximum_parallel: int,
    sample_limit: int,
    minimum_parallel: int = 1,
    sample_all_exact: bool = False,
) -> dict:
    graph_count = 0
    trials: Counter[tuple[int, int]] = Counter()
    numerical: Counter[tuple[int, int]] = Counter()
    exact_counts: Counter[tuple[int, int]] = Counter()
    discriminants: Counter[tuple[int, int, int]] = Counter()
    discriminant_four_counts: Counter[tuple[int, int]] = Counter()
    exact_samples: list[dict] = []
    for graph_index, graph in enumerate(_graphs(maximum_degree, edge_range)):
        graph_count += 1
        edges = tuple(tuple(sorted(edge)) for edge in graph.edges())
        for parallel_count in range(minimum_parallel, maximum_parallel + 1):
            for promoted in combinations(edges, parallel_count):
                key = (len(edges), parallel_count)
                trials[key] += 1
                matrix = _matrix(graph, promoted)
                if _numerical_signature(matrix) != (9, 1, 2):
                    continue
                numerical[key] += 1
                exact_matrix = sp.Matrix(matrix.tolist())
                if exact_matrix.rank() != 10:
                    continue
                exact_counts[key] += 1
                determinant = _radical_quotient_determinant(matrix)
                discriminants[(len(edges), parallel_count, determinant)] += 1
                if abs(determinant) == 4:
                    discriminant_four_counts[key] += 1
                if (sample_all_exact or abs(determinant) == 4) and len(exact_samples) < sample_limit:
                    exact_samples.append(
                        {
                            "source_graph_index_zero_based": graph_index,
                            "underlying_graph6": nx.to_graph6_bytes(
                                graph, header=False
                            ).decode().strip(),
                            "underlying_edge_count": len(edges),
                            "promoted_parallel_edges_zero_based": promoted,
                            "radical_quotient_determinant": determinant,
                            "negative_cartan": matrix.tolist(),
                        }
                    )
    return {
        "scope": {
            "vertices": 12,
            "connected": True,
            "maximum_degree": maximum_degree,
            "underlying_edge_range": edge_range,
            "maximum_parallel_edges": maximum_parallel,
            "minimum_parallel_edges": minimum_parallel,
            "target_positive_negative_zero_inertia": [9, 1, 2],
        },
        "underlying_graph_count": graph_count,
        "trial_count_by_edges_and_parallel": {
            f"{edges}:{parallel}": count
            for (edges, parallel), count in sorted(trials.items())
        },
        "numerical_survivor_count_by_edges_and_parallel": {
            f"{edges}:{parallel}": count
            for (edges, parallel), count in sorted(numerical.items())
        },
        "exact_labeled_survivor_count_by_edges_and_parallel": {
            f"{edges}:{parallel}": count
            for (edges, parallel), count in sorted(exact_counts.items())
        },
        "exact_labeled_survivor_count": sum(exact_counts.values()),
        "radical_quotient_determinant_histogram": {
            f"{edges}:{parallel}:{determinant}": count
            for (edges, parallel, determinant), count in sorted(discriminants.items())
        },
        "absolute_discriminant_four_count_by_edges_and_parallel": {
            f"{edges}:{parallel}": count
            for (edges, parallel), count in sorted(discriminant_four_counts.items())
        },
        "exact_sample_limit": sample_limit,
        "exact_sample_filter": (
            "all_exact_lorentzian_corank_two"
            if sample_all_exact
            else "absolute_radical_quotient_determinant_four"
        ),
        "exact_labeled_survivor_samples": exact_samples,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--maximum-degree", type=int, default=3)
    parser.add_argument("--edge-range", default="10:18")
    parser.add_argument("--maximum-parallel", type=int, default=3)
    parser.add_argument("--minimum-parallel", type=int, default=1)
    parser.add_argument("--sample-limit", type=int, default=20)
    parser.add_argument("--sample-all-exact", action="store_true")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()
    payload = build(
        args.maximum_degree,
        args.edge_range,
        args.maximum_parallel,
        args.sample_limit,
        args.minimum_parallel,
        args.sample_all_exact,
    )
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        # The scout normally writes through the caller's shell redirection; the
        # option is retained for local exploratory use only.
        from pathlib import Path

        Path(args.output).write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
