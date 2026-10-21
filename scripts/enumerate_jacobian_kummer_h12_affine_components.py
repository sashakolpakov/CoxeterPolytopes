"""Enumerate all connected affine components in the H12 Kummer diagram.

The search is exact.  Apart from affine A1 (parallel pairs), every proper
two-node subdiagram of an irreducible affine diagram is elliptic.  Therefore
the non-A1 search grows only connected elliptic subsets and permits a new
node only when every new two-node subdiagram is elliptic or orthogonal.  This
removes ultraparallel and parallel branches before matrix arithmetic.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import json
from pathlib import Path
import time

from scripts.audit_jacobian_kummer_h12_parallel_cusps import (
    INPUT,
    _affine_component,
    _positive_definite,
)


OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_affine_neighborhood_audit.json"
)


def _components(matrix: list[list[int]], vertices: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Connected components of an induced nonzero-edge graph."""

    unseen = set(vertices)
    result: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        component = {start}
        frontier = [start]
        while frontier:
            vertex = frontier.pop()
            neighbors = {other for other in unseen if matrix[vertex][other] != 0}
            unseen.difference_update(neighbors)
            component.update(neighbors)
            frontier.extend(neighbors)
        result.append(tuple(sorted(component)))
    return tuple(result)


def enumerate_components(path: Path = INPUT) -> dict[str, object]:
    payload = json.loads(path.read_text())
    gram = [[int(entry) for entry in row] for row in payload["positive_wall_gram"]]
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
        {
            other
            for other in compatible[vertex]
            if gram[vertex][other] != 0
        }
        for vertex in range(size)
    ]

    elliptic: set[tuple[int, ...]] = {(vertex,) for vertex in range(size)}
    all_elliptic: list[tuple[int, ...]] = sorted(elliptic)
    affine: list[tuple[int, ...]] = []
    order_counts: list[dict[str, int | float]] = []
    for order in range(2, 15):
        started = time.monotonic()
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
        next_affine: list[tuple[int, ...]] = []
        for facets in sorted(candidates):
            subgram = [[gram[left][right] for right in facets] for left in facets]
            if _positive_definite(subgram):
                next_elliptic.add(facets)
            elif _affine_component(subgram):
                next_affine.append(facets)
        affine.extend(next_affine)
        elapsed = time.monotonic() - started
        order_counts.append(
            {
                "order": order,
                "candidate_count": len(candidates),
                "elliptic_count": len(next_elliptic),
                "affine_count": len(next_affine),
            }
        )
        print({**order_counts[-1], "seconds": round(elapsed, 6)}, flush=True)
        all_elliptic.extend(sorted(next_elliptic))
        elliptic = next_elliptic
        if not elliptic:
            break

    parallel_pairs = tuple(
        (right, left)
        for left in range(size)
        for right in range(left)
        if gram[left][left] * gram[right][right]
        == gram[left][right] * gram[left][right]
    )
    all_affine = (*parallel_pairs, *affine)

    components_by_neighborhood: dict[tuple[int, ...], list[tuple[int, ...]]] = defaultdict(list)
    for component in all_affine:
        component_set = set(component)
        neighborhood = tuple(
            vertex
            for vertex in range(size)
            if vertex in component_set
            or all(gram[vertex][node] == 0 for node in component)
        )
        components_by_neighborhood[neighborhood].append(component)

    affine_cache: dict[tuple[int, ...], bool] = {}
    invalid_neighborhoods: list[dict[str, object]] = []
    neighborhood_signatures: Counter[tuple[int, ...]] = Counter()
    for neighborhood, source_components in sorted(components_by_neighborhood.items()):
        local_components = _components(gram, neighborhood)
        signature = tuple(sorted(len(component) - 1 for component in local_components))
        neighborhood_signatures[signature] += 1
        for component in local_components:
            if component not in affine_cache:
                component_gram = [
                    [gram[left][right] for right in component] for left in component
                ]
                affine_cache[component] = _affine_component(component_gram)
        if not (
            len(neighborhood) == 11 + len(local_components)
            and all(affine_cache[component] for component in local_components)
        ):
            invalid_neighborhoods.append(
                {
                    "neighborhood": [vertex + 1 for vertex in neighborhood],
                    "component_orders": [len(component) for component in local_components],
                    "source_affine_component_count": len(source_components),
                }
            )

    return {
        "schema": "jacobian_kummer_h12_connected_affine_components_v1",
        "scope": "complete exact connected-affine and full zero-neighborhood Vinberg audit",
        "search_order_counts": order_counts,
        "parallel_a1_count": len(parallel_pairs),
        "non_a1_affine_count": len(affine),
        "affine_order_histogram": dict(sorted(Counter(map(len, all_affine)).items())),
        "connected_affine_components": [
            [vertex + 1 for vertex in component] for component in all_affine
        ],
        "connected_elliptic_components": [
            [vertex + 1 for vertex in component] for component in all_elliptic
        ],
        "distinct_full_zero_neighborhood_count": len(components_by_neighborhood),
        "full_zero_neighborhood_size_histogram": dict(
            sorted(Counter(map(len, components_by_neighborhood)).items())
        ),
        "affine_components_per_neighborhood_histogram": dict(
            sorted(Counter(map(len, components_by_neighborhood.values())).items())
        ),
        "neighborhood_component_rank_signature_histogram": {
            "+".join(map(str, signature)): count
            for signature, count in sorted(neighborhood_signatures.items())
        },
        "distinct_neighborhood_components_checked": len(affine_cache),
        "invalid_full_zero_neighborhood_count": len(invalid_neighborhoods),
        "invalid_full_zero_neighborhoods": invalid_neighborhoods,
    }


def main() -> None:
    result = enumerate_components()
    OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    summary = {
        key: result[key]
        for key in (
            "parallel_a1_count",
            "non_a1_affine_count",
            "affine_order_histogram",
            "distinct_full_zero_neighborhood_count",
            "full_zero_neighborhood_size_histogram",
            "invalid_full_zero_neighborhood_count",
        )
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
