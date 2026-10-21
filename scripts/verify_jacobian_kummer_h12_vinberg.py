"""Literal exact Vinberg audit for the 172-facet H12 Kummer chamber.

This verifier is diagram-internal.  It constructs every elliptic rank-11
edge and rank-12 finite vertex as an orthogonal union of connected elliptic
components, constructs every ideal endpoint from the complete affine
zero-neighbourhood census, and checks that each edge has exactly two
endpoints.  It does not use the independent cdd extreme-ray certificate.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from itertools import product
import json
from pathlib import Path
import time
from typing import Iterator

from scripts.enumerate_jacobian_kummer_h12_affine_components import (
    INPUT,
    OUTPUT as AFFINE_AUDIT,
    _components,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_vinberg_volume_certificate.json"
)


def _elliptic_unions(
    gram: list[list[int]],
    connected_components: list[tuple[int, ...]],
    target_order: int,
) -> Iterator[tuple[int, ...]]:
    """Yield every elliptic diagram of ``target_order`` exactly once.

    Every elliptic Coxeter diagram is the orthogonal union of its uniquely
    determined connected elliptic components.  Components are selected in
    increasing order of their least vertex, which makes the enumeration
    duplication-free.
    """

    size = len(gram)
    all_vertices = (1 << size) - 1
    zero_masks = [
        sum(1 << other for other in range(size) if other != vertex and gram[vertex][other] == 0)
        for vertex in range(size)
    ]
    records_by_minimum: list[list[tuple[int, int, int]]] = [
        [] for _ in range(size)
    ]
    for component in connected_components:
        mask = sum(1 << vertex for vertex in component)
        orthogonal = all_vertices
        for vertex in component:
            orthogonal &= zero_masks[vertex]
        records_by_minimum[min(component)].append(
            (mask, orthogonal, len(component))
        )

    def recurse(
        remaining: int,
        minimum: int,
        allowed: int,
        selected: int,
    ) -> Iterator[tuple[int, ...]]:
        for least in range(minimum, size):
            if not (allowed >> least) & 1:
                continue
            for mask, orthogonal, order in records_by_minimum[least]:
                if order > remaining or mask & ~allowed:
                    continue
                combined = selected | mask
                if order == remaining:
                    yield tuple(
                        vertex for vertex in range(size) if (combined >> vertex) & 1
                    )
                    continue
                next_allowed = allowed & orthogonal
                next_minimum = least + 1
                available = next_allowed & ~((1 << next_minimum) - 1)
                if available.bit_count() >= remaining - order:
                    yield from recurse(
                        remaining - order,
                        next_minimum,
                        next_allowed,
                        combined,
                    )

    yield from recurse(target_order, 0, all_vertices, 0)


def _ideal_neighborhoods(
    gram: list[list[int]], affine_components: list[tuple[int, ...]]
) -> dict[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    size = len(gram)
    neighborhoods: dict[tuple[int, ...], tuple[tuple[int, ...], ...]] = {}
    for component in affine_components:
        component_set = set(component)
        neighborhood = tuple(
            vertex
            for vertex in range(size)
            if vertex in component_set
            or all(gram[vertex][node] == 0 for node in component)
        )
        neighborhoods.setdefault(neighborhood, _components(gram, neighborhood))
    return neighborhoods


def verify() -> dict[str, object]:
    metric = json.loads(INPUT.read_text())
    affine_audit = json.loads(AFFINE_AUDIT.read_text())
    gram = [[int(entry) for entry in row] for row in metric["positive_wall_gram"]]
    connected_elliptic = [
        tuple(vertex - 1 for vertex in component)
        for component in affine_audit["connected_elliptic_components"]
    ]
    connected_affine = [
        tuple(vertex - 1 for vertex in component)
        for component in affine_audit["connected_affine_components"]
    ]

    endpoint_counts: defaultdict[tuple[int, ...], list[int]] = defaultdict(
        lambda: [0, 0]
    )
    started = time.monotonic()
    finite_vertex_count = 0
    for vertex in _elliptic_unions(gram, connected_elliptic, 12):
        finite_vertex_count += 1
        for deleted in range(12):
            edge = vertex[:deleted] + vertex[deleted + 1 :]
            endpoint_counts[edge][0] += 1
    finite_seconds = time.monotonic() - started
    print(
        f"finite vertices {finite_vertex_count}; "
        f"incidences {sum(counts[0] for counts in endpoint_counts.values())}; "
        f"{finite_seconds:.3f}s",
        flush=True,
    )

    ideal_vertices = _ideal_neighborhoods(gram, connected_affine)
    for components in ideal_vertices.values():
        for deleted_nodes in product(*components):
            deleted = set(deleted_nodes)
            edge = tuple(
                vertex
                for component in components
                for vertex in component
                if vertex not in deleted
            )
            endpoint_counts[tuple(sorted(edge))][1] += 1

    endpoint_type_histogram: Counter[str] = Counter()
    endpoint_count_histogram: Counter[int] = Counter()
    missing_edges: list[tuple[int, ...]] = []
    edge_count = 0
    started = time.monotonic()
    for edge in _elliptic_unions(gram, connected_elliptic, 11):
        edge_count += 1
        finite_count, ideal_count = endpoint_counts.pop(edge, [0, 0])
        endpoint_count_histogram[finite_count + ideal_count] += 1
        endpoint_type_histogram["F" * finite_count + "I" * ideal_count] += 1
        if finite_count + ideal_count != 2 and len(missing_edges) < 100:
            missing_edges.append(edge)
    edge_seconds = time.monotonic() - started
    print(f"elliptic edges {edge_count}; {edge_seconds:.3f}s", flush=True)

    unrecognized_endpoint_edges = {
        edge: counts for edge, counts in endpoint_counts.items() if any(counts)
    }
    finite_incidence_count = finite_vertex_count * 12
    ideal_incidence_count = sum(
        endpoint_type.count("I") * edge_type_count
        for endpoint_type, edge_type_count in endpoint_type_histogram.items()
    )
    valid = (
        affine_audit["invalid_full_zero_neighborhood_count"] == 0
        and endpoint_count_histogram == Counter({2: edge_count})
        and not missing_edges
        and not unrecognized_endpoint_edges
    )
    return {
        "schema": "jacobian_kummer_h12_literal_vinberg_volume_certificate_v1",
        "method": (
            "complete connected elliptic/affine component enumeration and "
            "Vinberg edge-endpoint/full-zero-neighborhood criterion"
        ),
        "logically_independent_of_cdd_extreme_ray_certificate": True,
        "dimension": 12,
        "facet_count": len(gram),
        "connected_elliptic_component_count": len(connected_elliptic),
        "connected_affine_component_count": len(connected_affine),
        "finite_vertex_count": finite_vertex_count,
        "ideal_vertex_count": len(ideal_vertices),
        "elliptic_rank_11_edge_count": edge_count,
        "finite_endpoint_incidence_count": finite_incidence_count,
        "ideal_endpoint_incidence_count": ideal_incidence_count,
        "endpoint_count_histogram": dict(sorted(endpoint_count_histogram.items())),
        "endpoint_type_histogram": dict(sorted(endpoint_type_histogram.items())),
        "bad_edge_count": sum(
            count for endpoints, count in endpoint_count_histogram.items() if endpoints != 2
        ),
        "unrecognized_endpoint_edge_count": len(unrecognized_endpoint_edges),
        "invalid_full_zero_neighborhood_count": affine_audit[
            "invalid_full_zero_neighborhood_count"
        ],
        "finite_volume": valid,
        "compact": False if valid and ideal_vertices else None,
        "vinberg_decision": "finite_volume_noncompact" if valid else "infinite_volume",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = verify()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check_archive:
        if not OUTPUT.is_file() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"archived artifact differs: {OUTPUT}")
    else:
        OUTPUT.write_text(rendered, encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"checked {OUTPUT}" if arguments.check_archive else f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
