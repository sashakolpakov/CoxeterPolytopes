"""Literal exact Vinberg volume audit for Kummer H10 class_0007.

The verifier constructs every elliptic rank-9 edge and rank-10 finite
vertex as an orthogonal union of connected elliptic components.  Ideal
endpoints come from the complete connected-affine/full-zero-neighborhood
census.  Vinberg's criterion is satisfied exactly when every edge has two
finite/ideal endpoints and every affine neighborhood has parabolic rank 9.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import product
import json
from pathlib import Path
import time
from typing import Any

from scripts.enumerate_jacobian_kummer_h10_class0007_components import (
    OUTPUT as COMPONENT_AUDIT,
    _components,
    _integral_positive_gram,
)
from scripts.screen_jacobian_kummer_h10_class0007 import OUTPUT as SCREEN
from scripts.verify_jacobian_kummer_h12_vinberg import _elliptic_unions


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h10_class0007_vinberg_certificate.json"
)
DIMENSION = 10


def _ideal_neighborhoods(
    gram: list[list[int]], affine_components: list[tuple[int, ...]]
) -> dict[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    neighborhoods: dict[tuple[int, ...], tuple[tuple[int, ...], ...]] = {}
    for component in affine_components:
        component_set = set(component)
        neighborhood = tuple(
            vertex
            for vertex in range(len(gram))
            if vertex in component_set
            or all(gram[vertex][node] == 0 for node in component)
        )
        neighborhoods.setdefault(neighborhood, _components(gram, neighborhood))
    return neighborhoods


def verify(
    screen_path: Path = SCREEN,
    component_path: Path = COMPONENT_AUDIT,
) -> dict[str, Any]:
    metric = json.loads(screen_path.read_text(encoding="utf-8"))
    component_audit = json.loads(component_path.read_text(encoding="utf-8"))
    gram, scale = _integral_positive_gram(metric)
    dimension = int(metric["group_class"]["hyperbolic_dimension"])
    if not (
        metric["status"]["all_facet_pairs_are_coxeter"]
        and component_audit["dimension"] == dimension
        and component_audit["facet_count"]
        == len(gram)
        == metric["facet_count"]
        and component_audit["global_integral_gram_scale"] == scale
    ):
        raise ArithmeticError("the class_0007 Vinberg inputs are inconsistent")
    connected_elliptic = [
        tuple(vertex - 1 for vertex in component)
        for component in component_audit["connected_elliptic_components"]
    ]
    connected_affine = [
        tuple(vertex - 1 for vertex in component)
        for component in component_audit["connected_affine_components"]
    ]

    endpoint_counts: defaultdict[tuple[int, ...], list[int]] = defaultdict(
        lambda: [0, 0]
    )
    started = time.monotonic()
    finite_vertex_count = 0
    for vertex in _elliptic_unions(gram, connected_elliptic, dimension):
        finite_vertex_count += 1
        for deleted in range(dimension):
            edge = vertex[:deleted] + vertex[deleted + 1 :]
            endpoint_counts[edge][0] += 1
    finite_seconds = time.monotonic() - started
    print(
        f"finite vertices {finite_vertex_count}; "
        f"incidences {finite_vertex_count * dimension}; "
        f"{finite_seconds:.3f}s",
        flush=True,
    )

    ideal_vertices = _ideal_neighborhoods(gram, connected_affine)
    for components in ideal_vertices.values():
        if not (
            sum(len(component) - 1 for component in components)
            == dimension - 1
        ):
            raise ArithmeticError("an archived ideal neighborhood has wrong rank")
        for deleted_nodes in product(*components):
            deleted = set(deleted_nodes)
            edge = tuple(
                sorted(
                    vertex
                    for component in components
                    for vertex in component
                    if vertex not in deleted
                )
            )
            if len(edge) != dimension - 1:
                raise ArithmeticError("a cusp deletion did not produce a codimension-one edge")
            endpoint_counts[edge][1] += 1

    endpoint_type_histogram: Counter[str] = Counter()
    endpoint_count_histogram: Counter[int] = Counter()
    bad_edges: list[dict[str, Any]] = []
    edge_count = 0
    started = time.monotonic()
    for edge in _elliptic_unions(
        gram, connected_elliptic, dimension - 1
    ):
        edge_count += 1
        finite_count, ideal_count = endpoint_counts.pop(edge, [0, 0])
        endpoint_count = finite_count + ideal_count
        endpoint_count_histogram[endpoint_count] += 1
        endpoint_type_histogram["F" * finite_count + "I" * ideal_count] += 1
        if endpoint_count != 2 and len(bad_edges) < 100:
            bad_edges.append(
                {
                    "facets": [vertex + 1 for vertex in edge],
                    "finite_endpoint_count": finite_count,
                    "ideal_endpoint_count": ideal_count,
                }
            )
    edge_seconds = time.monotonic() - started
    print(f"elliptic edges {edge_count}; {edge_seconds:.3f}s", flush=True)

    unrecognized = {
        edge: counts for edge, counts in endpoint_counts.items() if any(counts)
    }
    finite_incidence_count = finite_vertex_count * dimension
    ideal_incidence_count = sum(
        endpoint_type.count("I") * count
        for endpoint_type, count in endpoint_type_histogram.items()
    )
    invalid_neighborhood_count = component_audit[
        "invalid_full_zero_neighborhood_count"
    ]
    valid = (
        invalid_neighborhood_count == 0
        and endpoint_count_histogram == Counter({2: edge_count})
        and not bad_edges
        and not unrecognized
        and edge_count > 0
        and (finite_vertex_count > 0 or ideal_vertices)
    )
    return {
        "schema": "jacobian_kummer_h10_class0007_literal_vinberg_volume_certificate_v1",
        "method": (
            "complete connected elliptic/affine component enumeration and "
            "Vinberg edge-endpoint/full-zero-neighborhood criterion"
        ),
        "logically_independent_of_exact_double_description": True,
        "dimension": dimension,
        "facet_count": len(gram),
        "global_integral_gram_scale": scale,
        "connected_elliptic_component_count": len(connected_elliptic),
        "connected_affine_component_count": len(connected_affine),
        "finite_vertex_count": finite_vertex_count,
        "ideal_vertex_count": len(ideal_vertices),
        "elliptic_rank_9_edge_count": edge_count,
        "finite_endpoint_incidence_count": finite_incidence_count,
        "ideal_endpoint_incidence_count": ideal_incidence_count,
        "endpoint_count_histogram": {
            str(key): value for key, value in sorted(endpoint_count_histogram.items())
        },
        "endpoint_type_histogram": dict(sorted(endpoint_type_histogram.items())),
        "bad_edge_count": sum(
            count
            for endpoints, count in endpoint_count_histogram.items()
            if endpoints != 2
        ),
        "first_bad_edges": bad_edges,
        "unrecognized_endpoint_edge_count": len(unrecognized),
        "invalid_full_zero_neighborhood_count": invalid_neighborhood_count,
        "finite_volume": valid,
        "compact": False if valid and ideal_vertices else (True if valid else None),
        "vinberg_decision": (
            "finite_volume_noncompact"
            if valid and ideal_vertices
            else "compact"
            if valid
            else "infinite_volume"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen", type=Path, default=SCREEN)
    parser.add_argument("--components", type=Path, default=COMPONENT_AUDIT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = verify(arguments.screen, arguments.components)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not arguments.output.is_file() or arguments.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {arguments.output}")
    else:
        arguments.output.write_text(rendered, encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
