"""Audit every parallel-pair cusp neighbourhood of the H12 Kummer chamber.

This is a deliberately partial Vinberg audit: every affine A1 component
(a parallel pair of facets) is treated exhaustively, and its *full* zero
neighbourhood is checked for parabolic rank 11.  Affine components without a
parallel pair still have to be enumerated before finite volume can be claimed.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import json
from pathlib import Path

from coxeter_polytopes.exact.diagram import connected_components


ROOT = Path(__file__).resolve().parents[1]
INPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h12_interim_facet_gram.json"
)


def _determinant(matrix: list[list[int]]) -> int:
    """Exact fraction-free determinant for a small integral matrix."""

    size = len(matrix)
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (
                    row
                    for row in range(pivot_index + 1, size)
                    if work[row][pivot_index] != 0
                ),
                None,
            )
            if swap is None:
                return 0
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                if numerator % previous:
                    raise ArithmeticError("Bareiss division was not exact")
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def _positive_definite(matrix: list[list[int]]) -> bool:
    """Exact Sylvester test using integral leading principal minors."""

    return all(
        _determinant([row[:order] for row in matrix[:order]]) > 0
        for order in range(1, len(matrix) + 1)
    )


def _affine_component(matrix: list[list[int]]) -> bool:
    """Intrinsic exact affine test: det zero and every deletion elliptic."""

    size = len(matrix)
    return (
        size >= 2
        and _determinant(matrix) == 0
        and all(
            _positive_definite(
                [
                    [matrix[left][right] for right in range(size) if right != deleted]
                    for left in range(size)
                    if left != deleted
                ]
            )
            for deleted in range(size)
        )
    )


def audit(path: Path = INPUT) -> dict[str, object]:
    payload = json.loads(path.read_text())
    gram = [[int(entry) for entry in row] for row in payload["positive_wall_gram"]]
    size = len(gram)

    pairs_by_neighborhood: dict[tuple[int, ...], list[tuple[int, int]]] = defaultdict(list)
    for left in range(size):
        for right in range(left):
            if (
                gram[left][left] * gram[right][right]
                != gram[left][right] * gram[left][right]
            ):
                continue
            neighborhood = tuple(
                vertex
                for vertex in range(size)
                if vertex in (left, right)
                or (gram[vertex][left] == 0 and gram[vertex][right] == 0)
            )
            pairs_by_neighborhood[neighborhood].append((right, left))

    invalid: list[dict[str, object]] = []
    component_signatures: Counter[tuple[int, ...]] = Counter()
    multiplicities: Counter[int] = Counter()
    affine_cache: dict[tuple[int, ...], bool] = {}
    for neighborhood, pairs in sorted(pairs_by_neighborhood.items()):
        local_gram = [[gram[left][right] for right in neighborhood] for left in neighborhood]
        local_components = connected_components(local_gram)
        signature = tuple(sorted(len(component) - 1 for component in local_components))
        component_signatures[signature] += 1
        multiplicities[len(pairs)] += 1
        global_components = tuple(
            tuple(neighborhood[local] for local in component)
            for component in local_components
        )
        for component in global_components:
            if component not in affine_cache:
                component_gram = [
                    [gram[left][right] for right in component] for left in component
                ]
                affine_cache[component] = _affine_component(component_gram)
        valid = (
            len(neighborhood) == 11 + len(local_components)
            and all(affine_cache[component] for component in global_components)
        )
        if not valid:
            invalid.append(
                {
                    "neighborhood": [vertex + 1 for vertex in neighborhood],
                    "parallel_pairs": [
                        [left + 1, right + 1] for left, right in pairs
                    ],
                    "component_orders": [len(component) for component in local_components],
                }
            )

    return {
        "schema": "jacobian_kummer_h12_parallel_a1_cusp_audit_v1",
        "scope": (
            "exhaustive over affine A1 (parallel-pair) components only; "
            "not by itself a complete Vinberg finite-volume certificate"
        ),
        "parallel_pair_count": sum(map(len, pairs_by_neighborhood.values())),
        "distinct_full_zero_neighborhood_count": len(pairs_by_neighborhood),
        "neighborhood_size_histogram": dict(
            sorted(Counter(map(len, pairs_by_neighborhood)).items())
        ),
        "parallel_pair_multiplicity_histogram": dict(sorted(multiplicities.items())),
        "affine_component_rank_signature_histogram": {
            "+".join(map(str, signature)): count
            for signature, count in sorted(component_signatures.items())
        },
        "invalid_neighborhood_count": len(invalid),
        "invalid_neighborhoods": invalid,
        "distinct_connected_components_checked": len(affine_cache),
        "non_affine_connected_component_count": sum(
            not valid for valid in affine_cache.values()
        ),
    }


def main() -> None:
    result = audit()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
