"""Literal Vinberg face/cusp audit of the D0004/C007 H10 chamber.

The verifier independently combines two complete exact enumerations:

* the Coxeter diagram side enumerates every elliptic rank-9 edge diagram,
  every elliptic rank-10 finite endpoint, and every connected affine
  subdiagram (including the full zero-neighbourhood test); and
* the GMP double-description archive supplies every projective vertex ray
  and its complete wall-incidence set.

The two sides must agree.  In particular, every elliptic edge has exactly two
finite/ideal endpoints, and every affine component has a full parabolic
rank-9 cusp neighbourhood.  Thus finite volume is not inferred from a
pairwise Coxeter gate or from a vertex count alone.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_shimada_disc16_d4_c7_h10_27_exact_cone import (
    CASE_INPUT,
    FACET_COUNT,
    RAY_OUTPUT,
    ROOT,
    SUMMARY_OUTPUT as CONE_CERTIFICATE,
    load_ray_archive,
)
from scripts.certify_shimada_disc16_d5_h10_30_vinberg_faces import (
    DiagramEnumerator,
)


OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d4_c7_h10_27_vinberg_face_certificate.json"
)
EDGE_OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "shimada_disc16_d4_c7_h10_27_vinberg_edges.jsonl.gz"
)
CASE_ID = "D0004_C007"
DIMENSION = 10


def _mask(indices: tuple[int, ...]) -> int:
    return sum(1 << index for index in indices)


def _indices(mask: int) -> tuple[int, ...]:
    return tuple(index for index in range(FACET_COUNT) if mask & (1 << index))


def build_certificate(
    case_path: Path = CASE_INPUT,
    cone_path: Path = CONE_CERTIFICATE,
    ray_path: Path = RAY_OUTPUT,
) -> tuple[dict[str, Any], bytes]:
    case = json.loads(case_path.read_text(encoding="utf-8"))
    cone = json.loads(cone_path.read_text(encoding="utf-8"))
    ray_header, ray_records = load_ray_archive(ray_path)
    if not (
        case["id"] == CASE_ID
        and case["exact_reduction"]["hyperbolic_facet_count"] == FACET_COUNT
        and case["coxeter_gate"]["all_pairs_coxeter"] is True
        and cone["decision"]["finite_volume"] is True
        and cone["exact_double_description"]["extreme_ray_count"] == 962
        and cone["exact_double_description"]["finite_vertex_count"] == 895
        and cone["exact_double_description"]["ideal_vertex_count"] == 67
        and cone["exact_double_description"]["hyperideal_vertex_count"] == 0
        and ray_header["ray_count"] == 962
    ):
        raise ArithmeticError("the pinned D0004/C007 H10 inputs changed")

    positive_gram = sp.Matrix(
        [
            [sp.Rational(value) for value in row]
            for row in case["coxeter_gate"]["positive_wall_gram"]
        ]
    )
    enumerator = DiagramEnumerator(positive_gram, dimension=DIMENSION)

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

    # A pair of complete endpoint incidence sets determines an edge exactly
    # when their common active mask has Coxeter rank n-1.  We retain every
    # intersection of size at least n-1, so a hidden nonsimple edge would be
    # detected rather than silently discarded.
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
            raise ArithmeticError(
                "an ideal endpoint is not a full rank-9 parabolic neighborhood"
            )
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
        "schema": "shimada_disc16_d4_c7_h10_27_vinberg_edges_v1",
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
        "schema": "shimada_disc16_d4_c7_h10_27_vinberg_face_certificate_v1",
        "source": {
            "case": str(case_path.relative_to(ROOT)),
            "cone_certificate": str(cone_path.relative_to(ROOT)),
            "ray_archive": str(ray_path.relative_to(ROOT)),
            "case_id": CASE_ID,
            "scope": "subgroup of Shimada archived chamber stabilizer Aut(D)",
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
                str(count): value
                for count, value in sorted(cusp_component_count.items())
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
            "all_55451_pair_compatible_rank_9_cliques_screened": compatible_clique_counts["9"]
            == 55451,
            "all_6221_elliptic_rank_9_diagrams_equal_exact_geometric_edges": len(
                finite_rank_n_minus_one
            )
            == 6221,
            "all_895_elliptic_rank_10_diagrams_equal_exact_finite_vertices": len(
                finite_rank_n
            )
            == 895,
            "every_edge_has_exactly_two_endpoints": True,
            "all_198_connected_affine_components_enumerated": len(affine) == 198,
            "every_affine_full_zero_neighborhood_is_a_certified_cusp": True,
            "every_ideal_vertex_is_componentwise_affine_of_total_rank_9": True,
            "ideal_vertices_exist": len(ideal_endpoint_masks) == 67,
        },
    }
    expected_endpoint_types = Counter({"ff": 3161, "fi": 2628, "ii": 432})
    expected_cliques = {
        "1": 27,
        "2": 264,
        "3": 1446,
        "4": 5238,
        "5": 13683,
        "6": 27177,
        "7": 42516,
        "8": 53658,
        "9": 55451,
        "10": 47223,
        "11": 33078,
        "12": 18870,
    }
    if not (
        len(edge_endpoint_pairs) == 6221
        and len(finite_endpoint_masks) == 895
        and len(ideal_endpoint_masks) == 67
        and endpoint_kind_histogram == expected_endpoint_types
        and len(affine) == 198
        and compatible_clique_counts == expected_cliques
    ):
        raise ArithmeticError("the literal D0004/C007 Vinberg census changed")
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
    parser.add_argument("--case", type=Path, default=CASE_INPUT)
    parser.add_argument("--cone-certificate", type=Path, default=CONE_CERTIFICATE)
    parser.add_argument("--ray-archive", type=Path, default=RAY_OUTPUT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--edge-output", type=Path, default=EDGE_OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    arguments = parser.parse_args()
    summary, edge_archive = build_certificate(
        arguments.case, arguments.cone_certificate, arguments.ray_archive
    )
    _write_or_check_text(
        arguments.output,
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        arguments.check_archive,
    )
    _write_or_check_bytes(arguments.edge_output, edge_archive, arguments.check_archive)
    print(
        json.dumps(
            summary["decision"] | summary["literal_endpoint_audit"],
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
