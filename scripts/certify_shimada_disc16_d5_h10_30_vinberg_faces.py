"""Literal face-driven Vinberg audit of the disc16/D5 H10 chamber.

The blind Vinberg scan of all 30-choose-9 and 30-choose-10 subsets is wasteful.
This verifier is nevertheless exhaustive.  It uses exact pair compatibility
to enumerate every elliptic rank-9 and rank-10 Coxeter subdiagram, and exact
clique/deletion/determinant tests to enumerate every connected affine
subdiagram.  Independently, the GMP double-description archive supplies the
complete vertex and face incidence lattice.  The two enumerations must agree:

* every elliptic rank-9 diagram is one geometric edge with two endpoints;
* every elliptic rank-10 diagram is one finite endpoint; and
* the full zero-neighbourhood of every affine component is one rank-9
  parabolic ideal endpoint.

Thus this is Vinberg's edge/endpoint and full-parabolic-neighbourhood test,
not a count inferred only from pairwise Coxeter compatibility.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from functools import lru_cache
import gzip
import json
from pathlib import Path
from typing import Any, Iterator

import sympy as sp

from scripts.certify_shimada_disc16_d5_full_subgroup_fishing import (
    OUTPUT as CENSUS_OUTPUT,
    ROOT,
)
from scripts.certify_shimada_disc16_d5_h10_30_exact_cone import (
    RAY_OUTPUT,
    SUMMARY_OUTPUT as CONE_CERTIFICATE,
    load_ray_archive,
)


OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h10_30_vinberg_face_certificate.json"
)
EDGE_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d5_h10_30_vinberg_edges.jsonl.gz"
)
DIMENSION = 10
FACET_COUNT = 30
CANDIDATE_KEY = "class_03_H10"


def _mask(indices: tuple[int, ...]) -> int:
    return sum(1 << index for index in indices)


def _indices(mask: int) -> tuple[int, ...]:
    return tuple(index for index in range(FACET_COUNT) if mask & (1 << index))


class DiagramEnumerator:
    """Exact finite/affine classifier for the rational Coxeter diagram."""

    def __init__(self, positive_gram: sp.Matrix, dimension: int = DIMENSION):
        if positive_gram.rows != positive_gram.cols:
            raise ValueError("the Vinberg Gram must be square")
        self.gram = positive_gram
        self.facet_count = positive_gram.rows
        self.dimension = dimension
        self.labels = [[2] * self.facet_count for _ in range(self.facet_count)]
        self.compatible_adjacency: list[int] = []
        self.nonorthogonal_adjacency: list[int] = []
        for left in range(self.facet_count):
            compatible = 0
            nonorthogonal = 0
            for right in range(self.facet_count):
                if left == right:
                    continue
                product = sp.Rational(positive_gram[left, right])
                if product:
                    nonorthogonal |= 1 << right
                ratio = sp.factor(
                    product**2
                    / (
                        positive_gram[left, left]
                        * positive_gram[right, right]
                    )
                )
                label = (
                    2
                    if product == 0
                    else 3
                    if ratio == sp.Rational(1, 4)
                    else 4
                    if ratio == sp.Rational(1, 2)
                    else 6
                    if ratio == sp.Rational(3, 4)
                    else "infinity"
                    if ratio == 1
                    else "ultraparallel"
                    if ratio > 1
                    else "forbidden"
                )
                self.labels[left][right] = label
                if label in {2, 3, 4, 6}:
                    compatible |= 1 << right
            self.compatible_adjacency.append(compatible)
            self.nonorthogonal_adjacency.append(nonorthogonal)
        self.finite_nonorthogonal_adjacency = tuple(
            sum(
                1 << right
                for right in range(self.facet_count)
                if self.labels[left][right] in {3, 4, 6}
            )
            for left in range(self.facet_count)
        )

    def _components(self, mask: int, *, all_nonorthogonal: bool) -> tuple[int, ...]:
        adjacency = (
            self.nonorthogonal_adjacency
            if all_nonorthogonal
            else self.finite_nonorthogonal_adjacency
        )
        remaining = mask
        result = []
        while remaining:
            start = remaining & -remaining
            remaining ^= start
            component = 0
            frontier = start
            while frontier:
                node_bit = frontier & -frontier
                frontier ^= node_bit
                node = node_bit.bit_length() - 1
                component |= node_bit
                new = adjacency[node] & remaining
                remaining ^= new
                frontier |= new
            result.append(component)
        return tuple(result)

    def _finite_connected_type(self, mask: int) -> str | None:
        nodes = tuple(
            index for index in range(self.facet_count) if mask & (1 << index)
        )
        order = len(nodes)
        if order == 1:
            return "A1"
        adjacency = {
            node: {
                other
                for other in nodes
                if other != node and self.labels[node][other] in {3, 4, 6}
            }
            for node in nodes
        }
        edges = [
            (left, right, self.labels[left][right])
            for position, left in enumerate(nodes)
            for right in nodes[position + 1 :]
            if self.labels[left][right] in {3, 4, 6}
        ]
        if len(edges) != order - 1:
            return None
        degrees = {node: len(adjacency[node]) for node in nodes}
        six_edges = [(left, right) for left, right, label in edges if label == 6]
        four_edges = [(left, right) for left, right, label in edges if label == 4]
        if six_edges:
            return "G2" if order == 2 and len(six_edges) == 1 else None
        if max(degrees.values()) <= 2:
            if not four_edges:
                return f"A{order}"
            if len(four_edges) != 1:
                return None
            left, right = four_edges[0]
            if degrees[left] == 1 or degrees[right] == 1:
                return f"B{order}"
            return "F4" if order == 4 else None
        if four_edges or max(degrees.values()) > 3:
            return None
        branch = [node for node in nodes if degrees[node] == 3]
        if len(branch) != 1:
            return None
        branch_node = branch[0]
        arm_lengths = []
        for neighbour in adjacency[branch_node]:
            previous, current, length = branch_node, neighbour, 1
            while degrees[current] == 2:
                following = next(
                    node for node in adjacency[current] if node != previous
                )
                previous, current = current, following
                length += 1
            arm_lengths.append(length)
        arms = tuple(sorted(arm_lengths))
        if arms[:2] == (1, 1):
            return f"D{order}"
        return {
            (1, 2, 2): "E6",
            (1, 2, 3): "E7",
            (1, 2, 4): "E8",
        }.get(arms)

    @lru_cache(maxsize=None)
    def finite_signature(self, mask: int) -> tuple[str, ...] | None:
        components = self._components(mask, all_nonorthogonal=False)
        types = tuple(self._finite_connected_type(component) for component in components)
        if any(kind is None for kind in types):
            return None
        return tuple(sorted(kind for kind in types if kind is not None))

    def _cliques(self, order: int) -> Iterator[int]:
        def recurse(chosen: int, candidates: int, needed: int) -> Iterator[int]:
            if needed == 0:
                yield chosen
                return
            while candidates.bit_count() >= needed:
                bit = candidates & -candidates
                candidates ^= bit
                vertex = bit.bit_length() - 1
                yield from recurse(
                    chosen | bit,
                    candidates & self.compatible_adjacency[vertex],
                    needed - 1,
                )

        yield from recurse(0, (1 << self.facet_count) - 1, order)

    def finite_subdiagrams(self, order: int) -> dict[int, tuple[str, ...]]:
        return {
            mask: signature
            for mask in self._cliques(order)
            if (signature := self.finite_signature(mask)) is not None
        }

    def connected_affine_subdiagrams(self) -> dict[int, dict[str, Any]]:
        """Enumerate all irreducible affine subsets of Coxeter rank <= 11.

        For order at least three, every proper two-node subset of an affine
        Coxeter diagram is elliptic, hence the nodes form a compatibility
        clique.  Every vertex deletion is a finite Coxeter diagram.  These
        necessary filters leave a tiny list on which determinant and rank are
        checked exactly.  Parallel two-node affine A1 diagrams are added
        separately.  Orders through n+2 are scanned, matching Vinberg's bad-
        neighbourhood guard rather than assuming the desired rank bound.
        """

        affine: dict[int, dict[str, Any]] = {}
        for right in range(self.facet_count):
            for left in range(right):
                if self.labels[left][right] == "infinity":
                    mask = (1 << left) | (1 << right)
                    subgram = self.gram.extract((left, right), (left, right))
                    if subgram.det() != 0 or subgram.rank() != 1:
                        raise ArithmeticError("a parallel pair is not affine A1")
                    affine[mask] = {
                        "order": 2,
                        "coxeter_rank": 1,
                        "deleted_finite_signatures": [["A1"], ["A1"]],
                    }

        for order in range(3, self.dimension + 3):
            for mask in self._cliques(order):
                if len(self._components(mask, all_nonorthogonal=True)) != 1:
                    continue
                nodes = tuple(
                    index
                    for index in range(self.facet_count)
                    if mask & (1 << index)
                )
                deleted = tuple(
                    self.finite_signature(mask ^ (1 << node)) for node in nodes
                )
                if any(signature is None for signature in deleted):
                    continue
                subgram = self.gram.extract(nodes, nodes)
                if subgram.det() != 0 or subgram.rank() != order - 1:
                    continue
                affine[mask] = {
                    "order": order,
                    "coxeter_rank": order - 1,
                    "deleted_finite_signatures": [list(signature) for signature in deleted],
                }
        return affine


def build_certificate(
    census_path: Path = CENSUS_OUTPUT,
    cone_path: Path = CONE_CERTIFICATE,
    ray_path: Path = RAY_OUTPUT,
) -> tuple[dict[str, Any], bytes]:
    census = json.loads(census_path.read_text(encoding="utf-8"))
    cone = json.loads(cone_path.read_text(encoding="utf-8"))
    ray_header, ray_records = load_ray_archive(ray_path)
    screen = census["above_H9"]["complete_section_screens"][CANDIDATE_KEY]
    if not (
        screen["facet_count"] == FACET_COUNT
        and screen["pair_gate"]["all_pairs_coxeter"] is True
        and cone["decision"]["finite_volume"] is True
        and cone["exact_double_description"]["extreme_ray_count"] == 1540
        and ray_header["ray_count"] == 1540
    ):
        raise ArithmeticError("the pinned H10 cone inputs changed")
    raw_gram = sp.Matrix(
        [[sp.Rational(value) for value in row] for row in screen["facet_gram"]]
    )
    positive_gram = -raw_gram
    enumerator = DiagramEnumerator(positive_gram)

    active_sets = tuple(
        tuple(index - 1 for index in record[1]) for record in ray_records
    )
    kinds = tuple(str(record[3]) for record in ray_records)
    active_masks = tuple(_mask(active) for active in active_sets)
    finite_endpoint_masks = {
        mask for mask, kind in zip(active_masks, kinds, strict=True) if kind == "finite"
    }
    ideal_endpoint_masks = {
        mask for mask, kind in zip(active_masks, kinds, strict=True) if kind == "ideal"
    }

    finite_rank_n_minus_one = enumerator.finite_subdiagrams(DIMENSION - 1)
    finite_rank_n = enumerator.finite_subdiagrams(DIMENSION)
    if set(finite_rank_n) != finite_endpoint_masks:
        raise ArithmeticError("elliptic rank-10 diagrams do not equal finite endpoints")

    edge_endpoint_pairs: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for right in range(len(active_masks)):
        for left in range(right):
            intersection = active_masks[left] & active_masks[right]
            if intersection.bit_count() >= DIMENSION - 1:
                edge_endpoint_pairs[intersection].append((left, right))
    if not (
        set(edge_endpoint_pairs) == set(finite_rank_n_minus_one)
        and all(mask.bit_count() == DIMENSION - 1 for mask in edge_endpoint_pairs)
        and all(len(pairs) == 1 for pairs in edge_endpoint_pairs.values())
    ):
        raise ArithmeticError("elliptic edge diagrams do not have exactly two endpoints")

    affine = enumerator.connected_affine_subdiagrams()
    affine_full_neighborhoods: dict[int, int] = {}
    for component_mask, record in affine.items():
        component = _indices(component_mask)
        zero_neighborhood = tuple(
            facet
            for facet in range(FACET_COUNT)
            if facet not in component
            and all(positive_gram[facet, node] == 0 for node in component)
        )
        neighborhood_mask = component_mask | _mask(zero_neighborhood)
        if neighborhood_mask not in ideal_endpoint_masks:
            raise ArithmeticError(
                "an affine component has a non-parabolic full zero-neighborhood"
            )
        affine_full_neighborhoods[component_mask] = neighborhood_mask
        record["full_zero_neighborhood_facet_positions"] = [
            index + 1 for index in _indices(neighborhood_mask)
        ]

    cusp_component_order_histogram: Counter[tuple[int, ...]] = Counter()
    cusp_component_count = Counter()
    for cusp_mask in ideal_endpoint_masks:
        components = enumerator._components(cusp_mask, all_nonorthogonal=True)
        if not (
            sum(component.bit_count() - 1 for component in components)
            == DIMENSION - 1
            and all(component in affine for component in components)
            and all(
                affine_full_neighborhoods[component] == cusp_mask
                for component in components
            )
        ):
            raise ArithmeticError("an ideal endpoint is not a full rank-9 parabolic neighborhood")
        orders = tuple(sorted(component.bit_count() for component in components))
        cusp_component_order_histogram[orders] += 1
        cusp_component_count[len(components)] += 1

    endpoint_kind_histogram: Counter[str] = Counter()
    edge_type_histogram: Counter[str] = Counter()
    edge_lines = []
    for edge_mask in sorted(edge_endpoint_pairs):
        left, right = edge_endpoint_pairs[edge_mask][0]
        endpoint_code = "".join(sorted((kinds[left][0], kinds[right][0])))
        endpoint_kind_histogram[endpoint_code] += 1
        signature = finite_rank_n_minus_one[edge_mask]
        edge_type_histogram[" + ".join(signature)] += 1
        edge_lines.append(
            json.dumps(
                [
                    [index + 1 for index in _indices(edge_mask)],
                    [left + 1, right + 1],
                    [kinds[left], kinds[right]],
                    list(signature),
                ],
                separators=(",", ":"),
            )
        )
    edge_header = {
        "schema": "shimada_disc16_d5_h10_30_vinberg_edges_v1",
        "record_format": "[edge_facet_positions, endpoint_ray_indices, endpoint_types, finite_diagram_components]",
        "edge_count": len(edge_lines),
        "ray_indexing": "one_based archive order",
    }
    edge_archive = gzip.compress(
        (
            json.dumps(edge_header, sort_keys=True, separators=(",", ":"))
            + "\n"
            + "\n".join(edge_lines)
            + "\n"
        ).encode("utf-8"),
        compresslevel=9,
        mtime=0,
    )

    compatible_clique_counts = {
        str(order): sum(1 for _ in enumerator._cliques(order))
        for order in range(1, DIMENSION + 3)
    }
    summary = {
        "schema": "shimada_disc16_d5_h10_30_vinberg_face_certificate_v1",
        "source": {
            "census": str(census_path.relative_to(ROOT)),
            "cone_certificate": str(cone_path.relative_to(ROOT)),
            "ray_archive": str(ray_path.relative_to(ROOT)),
            "candidate_key": CANDIDATE_KEY,
        },
        "method": {
            "criterion": "Vinberg edge-endpoint and full parabolic zero-neighborhood criterion",
            "enumeration": "all exact pair-compatible cliques through order 12, finite Coxeter classification, and exact affine determinant/rank",
            "face_completeness": "independent GMP rational double-description ray/incidence archive",
            "blind_binomial_subset_scan_used": False,
            "all_relevant_subdiagrams_enumerated": True,
        },
        "diagram_enumeration": {
            "pair_compatible_clique_count_by_order": compatible_clique_counts,
            "elliptic_rank_9_edge_diagram_count": len(finite_rank_n_minus_one),
            "elliptic_rank_10_finite_endpoint_diagram_count": len(finite_rank_n),
            "connected_affine_component_count": len(affine),
            "connected_affine_order_histogram": {
                str(order): count
                for order, count in sorted(
                    Counter(record["order"] for record in affine.values()).items()
                )
            },
            "edge_finite_type_histogram": dict(sorted(edge_type_histogram.items())),
        },
        "literal_endpoint_audit": {
            "edge_count": len(edge_endpoint_pairs),
            "edge_endpoint_count_histogram": {"2": len(edge_endpoint_pairs)},
            "edge_endpoint_type_histogram": dict(sorted(endpoint_kind_histogram.items())),
            "finite_vertex_count": len(finite_endpoint_masks),
            "ideal_vertex_count": len(ideal_endpoint_masks),
            "cusp_affine_component_count_histogram": {
                str(count): value for count, value in sorted(cusp_component_count.items())
            },
            "cusp_component_order_histogram": {
                "+".join(map(str, orders)): count
                for orders, count in sorted(cusp_component_order_histogram.items())
            },
            "bad_edge_count": 0,
            "bad_parabolic_neighborhood_count": 0,
        },
        "decision": {
            "finite_volume": True,
            "compact": False,
            "classification": "finite_volume_noncompact",
            "reason": "all elliptic rank-9 edges have exactly two finite/ideal endpoints and every affine component has a full rank-9 parabolic zero-neighborhood",
        },
        "checks": {
            "all_133463_pair_compatible_rank_9_cliques_screened": compatible_clique_counts["9"] == 133463,
            "all_10030_elliptic_rank_9_diagrams_equal_exact_geometric_edges": len(finite_rank_n_minus_one) == 10030,
            "all_1333_elliptic_rank_10_diagrams_equal_exact_finite_vertices": len(finite_rank_n) == 1333,
            "every_edge_has_exactly_two_endpoints": True,
            "all_436_connected_affine_components_enumerated": len(affine) == 436,
            "every_affine_full_zero_neighborhood_is_a_certified_cusp": True,
            "every_ideal_vertex_is_componentwise_affine_of_total_rank_9": True,
            "ideal_vertices_exist": len(ideal_endpoint_masks) == 207,
        },
    }
    expected_endpoint_types = Counter({"ff": 4423, "fi": 4484, "ii": 1123})
    if not (
        len(edge_endpoint_pairs) == 10030
        and len(finite_endpoint_masks) == 1333
        and len(ideal_endpoint_masks) == 207
        and endpoint_kind_histogram == expected_endpoint_types
        and len(affine) == 436
    ):
        raise ArithmeticError("the literal Vinberg census changed")
    return summary, edge_archive


def _write_or_check_text(path: Path, text: str, check: bool) -> None:
    if check:
        if not path.is_file() or path.read_text(encoding="utf-8") != text:
            raise SystemExit(f"archived artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _write_or_check_bytes(path: Path, payload: bytes, check: bool) -> None:
    if check:
        if not path.is_file() or path.read_bytes() != payload:
            raise SystemExit(f"archived artifact differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path, default=CENSUS_OUTPUT)
    parser.add_argument("--cone-certificate", type=Path, default=CONE_CERTIFICATE)
    parser.add_argument("--ray-archive", type=Path, default=RAY_OUTPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--edge-output", type=Path, default=EDGE_OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    summary, edge_archive = build_certificate(
        arguments.census, arguments.cone_certificate, arguments.ray_archive
    )
    _write_or_check_text(
        arguments.output,
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        arguments.check_archive,
    )
    _write_or_check_bytes(arguments.edge_output, edge_archive, arguments.check_archive)
    print(json.dumps(summary["decision"] | summary["literal_endpoint_audit"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
