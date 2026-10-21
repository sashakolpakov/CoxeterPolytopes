#!/usr/bin/env python3
"""Certify the exact index of the 84-wall JK-H10 reflection subgroup.

For a Coxeter group ``W`` the rational Euler characteristic is

    chi(W) = sum_{T finite special} (-1)^|T| / |W_T|.

The 84-wall computation consumes the already exhaustive connected-elliptic
component census and enumerates every mutually orthogonal product of those
components.  The 13-wall maximal chamber is small enough to enumerate
directly from its exact relation matrix.  The quotient of the two Euler
characteristics is the reflection-subgroup index.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
SCREEN = ROOT / "data" / "enumerations" / "jacobian_kummer_h10_class0007_screen.json"
COMPONENTS = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h10_class0007_component_audit.json"
)
MAXIMAL = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h10_maximal_reflection_lineage.json"
)
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "jacobian_kummer_h10_reflection_index.json"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fraction_payload(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def _relation_matrix(gram_payload: Sequence[Sequence[str]]) -> list[list[str]]:
    gram = [[Fraction(entry) for entry in row] for row in gram_payload]
    size = len(gram)
    result = [["1" if left == right else "2" for right in range(size)] for left in range(size)]
    labels = {
        Fraction(1, 4): "3",
        Fraction(1, 2): "4",
        Fraction(3, 4): "6",
        Fraction(1): "0",
    }
    for left in range(size):
        if len(gram[left]) != size or gram[left][left] <= 0:
            raise ArithmeticError("invalid positive wall Gram")
        for right in range(left + 1, size):
            if gram[left][right] != gram[right][left] or gram[left][right] > 0:
                raise ArithmeticError("wall Gram is not symmetric acute")
            if gram[left][right] == 0:
                relation = "2"
            else:
                ratio = gram[left][right] ** 2 / (gram[left][left] * gram[right][right])
                relation = labels.get(ratio, "1" if ratio > 1 else "")
                if not relation:
                    raise ArithmeticError(f"unsupported Coxeter ratio {ratio}")
            result[left][right] = result[right][left] = relation
    return result


def _vertices(mask: int) -> list[int]:
    result: list[int] = []
    while mask:
        bit = mask & -mask
        result.append(bit.bit_length() - 1)
        mask -= bit
    return result


def _finite_connected_order(mask: int, relations: Sequence[Sequence[str]]) -> int:
    """Return the order of a connected crystallographic finite Coxeter group."""

    vertices = _vertices(mask)
    rank = len(vertices)
    if rank == 1:
        return 2

    edges: list[tuple[int, int, int]] = []
    degree = {vertex: 0 for vertex in vertices}
    for offset, left in enumerate(vertices):
        for right in vertices[offset + 1 :]:
            relation = relations[left][right]
            if relation in {"0", "1"}:
                raise ValueError("parallel or ultraparallel pair")
            if relation != "2":
                label = int(relation)
                edges.append((left, right, label))
                degree[left] += 1
                degree[right] += 1
    if len(edges) != rank - 1:
        raise ValueError("finite connected Coxeter diagram must be a tree")

    labels = sorted(label for _, _, label in edges)
    if labels == [3] * (rank - 1):
        if max(degree.values()) <= 2:
            return factorial(rank + 1)  # A_rank
        degrees = sorted(degree.values())
        if degrees[-1] != 3 or degrees.count(1) != 3:
            raise ValueError("unclassified simply-laced tree")
        branch = next(vertex for vertex, value in degree.items() if value == 3)
        adjacency = {vertex: [] for vertex in vertices}
        for left, right, _ in edges:
            adjacency[left].append(right)
            adjacency[right].append(left)
        arm_lengths: list[int] = []
        for neighbor in adjacency[branch]:
            previous, current, length = branch, neighbor, 1
            while degree[current] == 2:
                following = next(vertex for vertex in adjacency[current] if vertex != previous)
                previous, current = current, following
                length += 1
            arm_lengths.append(length)
        arm_lengths.sort()
        if arm_lengths == [1, 1, rank - 3]:
            return (2 ** (rank - 1)) * factorial(rank)  # D_rank
        exceptional = {
            (6, (1, 2, 2)): 51_840,
            (7, (1, 2, 3)): 2_903_040,
            (8, (1, 2, 4)): 696_729_600,
        }
        try:
            return exceptional[(rank, tuple(arm_lengths))]
        except KeyError as error:
            raise ValueError("unclassified simply-laced finite diagram") from error

    if labels.count(4) == 1 and all(label in {3, 4} for label in labels):
        if max(degree.values()) > 2:
            raise ValueError("non-path diagram with a 4-edge")
        four_edge = next((left, right) for left, right, label in edges if label == 4)
        if degree[four_edge[0]] == 1 or degree[four_edge[1]] == 1:
            return (2**rank) * factorial(rank)  # B_rank = C_rank
        if rank == 4:
            return 1_152  # F4
        raise ValueError("internal 4-edge outside F4")

    if rank == 2 and len(labels) == 1 and labels[0] >= 3:
        return 2 * labels[0]  # I_2(m), including G2
    raise ValueError("unclassified connected finite Coxeter diagram")


def _connected_masks_direct(relations: Sequence[Sequence[str]]) -> set[int]:
    """Enumerate all finite connected masks; intended for the 13-wall base."""

    size = len(relations)
    result: set[int] = set()
    for mask in range(1, 1 << size):
        vertices = _vertices(mask)
        if any(
            relations[left][right] in {"0", "1"}
            for offset, left in enumerate(vertices)
            for right in vertices[offset + 1 :]
        ):
            continue
        adjacency = {vertex: set() for vertex in vertices}
        for offset, left in enumerate(vertices):
            for right in vertices[offset + 1 :]:
                if relations[left][right] != "2":
                    adjacency[left].add(right)
                    adjacency[right].add(left)
        seen = {vertices[0]}
        frontier = [vertices[0]]
        while frontier:
            vertex = frontier.pop()
            new = adjacency[vertex] - seen
            seen.update(new)
            frontier.extend(new)
        if len(seen) != len(vertices):
            continue
        try:
            _finite_connected_order(mask, relations)
        except ValueError:
            continue
        result.add(mask)
    return result


def _connected_masks_by_nonright_components(
    relations: Sequence[Sequence[str]],
) -> set[int]:
    """Enumerate finite connected masks componentwise in the nonright graph.

    A connected Coxeter diagram is connected by its labels greater than two;
    consequently it lies inside one connected component of that sparse graph.
    This is exactly equivalent to :func:`_connected_masks_direct`, but avoids
    scanning all ``2^facets`` masks when a many-wall chamber has a small
    nonright core plus numerous isolated nodes.
    """

    size = len(relations)
    adjacency = [set() for _ in range(size)]
    for left in range(size):
        if len(relations[left]) != size:
            raise ValueError("relation matrix must be square")
        for right in range(left):
            if relations[left][right] in {"3", "4", "6"}:
                adjacency[left].add(right)
                adjacency[right].add(left)

    unseen = set(range(size))
    graph_components: list[list[int]] = []
    while unseen:
        start = min(unseen)
        component = {start}
        frontier = [start]
        unseen.remove(start)
        while frontier:
            vertex = frontier.pop()
            new = adjacency[vertex] & unseen
            unseen.difference_update(new)
            component.update(new)
            frontier.extend(new)
        graph_components.append(sorted(component))

    result: set[int] = set()
    for component in graph_components:
        if len(component) > 24:
            raise ArithmeticError(
                "nonright component too large for componentwise subset enumeration"
            )
        for local_mask in range(1, 1 << len(component)):
            vertices = [
                component[offset]
                for offset in range(len(component))
                if local_mask & (1 << offset)
            ]
            if any(
                relations[left][right] in {"0", "1"}
                for offset, left in enumerate(vertices)
                for right in vertices[offset + 1 :]
            ):
                continue
            if len(vertices) > 1:
                seen = {vertices[0]}
                frontier = [vertices[0]]
                vertex_set = set(vertices)
                while frontier:
                    vertex = frontier.pop()
                    new = (adjacency[vertex] & vertex_set) - seen
                    seen.update(new)
                    frontier.extend(new)
                if seen != vertex_set:
                    continue
            mask = sum(1 << vertex for vertex in vertices)
            try:
                _finite_connected_order(mask, relations)
            except ValueError:
                continue
            result.add(mask)
    return result


def _euler_from_connected_components(
    relations: Sequence[Sequence[str]], connected_masks: set[int]
) -> dict[str, Any]:
    size = len(relations)
    orders = {
        mask: _finite_connected_order(mask, relations) for mask in connected_masks
    }
    compatible = [0] * size
    nonorthogonal = [0] * size
    for left in range(size):
        for right in range(size):
            if left == right:
                continue
            relation = relations[left][right]
            if relation in {"2", "3", "4", "6"}:
                compatible[left] |= 1 << right
            if relation in {"3", "4", "6"}:
                nonorthogonal[left] |= 1 << right

    counts = [0] * (size + 1)
    contributions = [Fraction(0) for _ in range(size + 1)]
    counts[0] = 1
    contributions[0] = Fraction(1)

    def recurse(
        components: tuple[int, ...],
        candidates: int,
        rank: int,
        reciprocal_order: Fraction,
    ) -> None:
        while candidates:
            bit = candidates & -candidates
            candidates -= bit
            vertex = bit.bit_length() - 1
            touching: list[int] = []
            untouched: list[int] = []
            for component in components:
                if nonorthogonal[vertex] & component:
                    touching.append(component)
                else:
                    untouched.append(component)
            merged = bit
            for component in touching:
                merged |= component
            if merged not in connected_masks:
                continue
            new_reciprocal = reciprocal_order
            for component in touching:
                new_reciprocal *= orders[component]
            new_reciprocal /= orders[merged]
            new_rank = rank + 1
            counts[new_rank] += 1
            contributions[new_rank] += (
                -new_reciprocal if new_rank % 2 else new_reciprocal
            )
            recurse(
                tuple((*untouched, merged)),
                candidates & compatible[vertex],
                new_rank,
                new_reciprocal,
            )

    recurse(tuple(), (1 << size) - 1, 0, Fraction(1))
    last_rank = max(index for index, count in enumerate(counts) if count)
    euler = sum(contributions, Fraction(0))
    connected_rank_histogram = Counter(mask.bit_count() for mask in connected_masks)
    connected_order_histogram = Counter(orders.values())
    return {
        "facet_count": size,
        "connected_elliptic_count": len(connected_masks),
        "connected_elliptic_rank_histogram": {
            str(rank): count for rank, count in sorted(connected_rank_histogram.items())
        },
        "connected_finite_group_order_histogram": {
            str(order): count for order, count in sorted(connected_order_histogram.items())
        },
        "finite_special_subset_count_by_rank": counts[: last_rank + 1],
        "euler_contribution_by_rank": [
            _fraction_payload(value) for value in contributions[: last_rank + 1]
        ],
        "rational_euler_characteristic": _fraction_payload(euler),
    }


def _factorization(value: int) -> dict[str, int]:
    result: dict[str, int] = {}
    divisor = 2
    while divisor * divisor <= value:
        while value % divisor == 0:
            result[str(divisor)] = result.get(str(divisor), 0) + 1
            value //= divisor
        divisor += 1
    if value > 1:
        result[str(value)] = result.get(str(value), 0) + 1
    return result


def build_certificate() -> dict[str, Any]:
    screen = json.loads(SCREEN.read_text(encoding="utf-8"))
    components = json.loads(COMPONENTS.read_text(encoding="utf-8"))
    maximal = json.loads(MAXIMAL.read_text(encoding="utf-8"))
    if not (
        screen["status"]["exact_irredundant_facets"]
        and screen["status"]["exact_full_coxeter_pair_gate"]
        and components["invalid_full_zero_neighborhood_count"] == 0
        and maximal["all_checks_pass"]
    ):
        raise ArithmeticError("an upstream exact gate failed")

    subgroup_relations = _relation_matrix(screen["positive_wall_gram"])
    subgroup_connected = {
        sum(1 << (vertex - 1) for vertex in component)
        for component in components["connected_elliptic_components"]
    }
    subgroup = _euler_from_connected_components(
        subgroup_relations, subgroup_connected
    )

    base_relations = maximal["vinberg_maximal_integral_reflection_chamber"][
        "normalized_coxeter_relation_matrix"
    ]
    base_connected = _connected_masks_direct(base_relations)
    base = _euler_from_connected_components(base_relations, base_connected)

    subgroup_euler = Fraction(subgroup["rational_euler_characteristic"])
    base_euler = Fraction(base["rational_euler_characteristic"])
    index = subgroup_euler / base_euler
    volume_coefficient = -subgroup_euler * Fraction(32, 945)
    coxiter_f_vector = [
        16_560,
        107_233,
        263_403,
        347_634,
        281_244,
        146_599,
        49_643,
        10_650,
        1_350,
        84,
        1,
    ]
    reconstructed_f_vector = [
        subgroup["finite_special_subset_count_by_rank"][10] + 862,
        *reversed(subgroup["finite_special_subset_count_by_rank"][1:10]),
        1,
    ]
    checks = {
        "subgroup_connected_count_matches_archive": len(subgroup_connected)
        == components["connected_elliptic_components"].__len__(),
        "subgroup_rank_10_count_matches_finite_vertices": subgroup[
            "finite_special_subset_count_by_rank"
        ][10]
        == 15_698,
        "subgroup_rank_9_count_matches_vinberg_edges": subgroup[
            "finite_special_subset_count_by_rank"
        ][9]
        == 107_233,
        "subgroup_euler_is_expected": subgroup_euler == Fraction(-527, 11_520),
        "base_has_13_facets": base["facet_count"] == 13,
        "base_euler_matches_coxiter": base_euler == Fraction(-31, 11_147_673_600),
        "euler_ratio_is_positive_integer": index.denominator == 1 and index > 0,
        "index_is_expected": index == 16_450_560,
        "volume_coefficient_is_expected": volume_coefficient == Fraction(527, 340_200),
        "independent_coxiter_f_vector_agrees": reconstructed_f_vector
        == coxiter_f_vector,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"reflection-index checks failed: {checks}")

    return {
        "schema": "jacobian_kummer_h10_reflection_index_v1",
        "sources": {
            "subgroup_screen": str(SCREEN.relative_to(ROOT)),
            "subgroup_screen_sha256": _sha256(SCREEN),
            "connected_component_audit": str(COMPONENTS.relative_to(ROOT)),
            "connected_component_audit_sha256": _sha256(COMPONENTS),
            "maximal_chamber_audit": str(MAXIMAL.relative_to(ROOT)),
            "maximal_chamber_audit_sha256": _sha256(MAXIMAL),
        },
        "method": "Solomon rational Euler sum over every finite special subgroup",
        "subgroup_84_wall_chamber": subgroup,
        "vinberg_f2_h10_13_wall_chamber": base,
        "index": {
            "formula": "chi(W_84) / chi(W_f2_H10)",
            "value": index.numerator,
            "prime_factorization": _factorization(index.numerator),
        },
        "hyperbolic_covolume": {
            "gauss_bonnet_formula": "vol = -(32*pi^5/945)*chi in H10",
            "value": f"pi^5 * {_fraction_payload(volume_coefficient)}",
        },
        "independent_coxiter_1_1_8_replay": {
            "input_exporter": "scripts/export_jacobian_kummer_h12_coxiter.py --input data/enumerations/jacobian_kummer_h10_class0007_screen.json",
            "finite_covolume": True,
            "cocompact": False,
            "f_vector": coxiter_f_vector,
            "vertices_at_infinity": 862,
            "rational_euler_characteristic": "-527/11520",
            "covolume": "pi^5 * 527/340200",
            "agreement_with_solomon_enumeration": True,
        },
        "lineage_decision": {
            "maximal_reflection_group": "Vinberg f_2^10",
            "subgroup_is_proper": True,
            "subgroup_index_is_exact": True,
            "literal_84_wall_subgroup_priority_remains_open": True,
        },
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    result = build_certificate()
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
