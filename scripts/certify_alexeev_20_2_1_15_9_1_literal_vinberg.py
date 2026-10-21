"""Literal Vinberg endpoint/parabolic audit of three Alexeev fixed sections.

The independent GMP cone census supplies the complete exact endpoint and
incidence list.  This verifier separately enumerates every elliptic
rank-(n-1) Coxeter edge diagram, every elliptic rank-n finite endpoint, and
every connected affine component through order n+2.  It then applies
Vinberg's full edge/endpoint and affine zero-neighbourhood criterion.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
from typing import Any

import sympy as sp

from scripts.certify_alexeev_20_2_1_15_9_1_exact_cones import (
    OUTPUT as CONE_CENSUS,
    ROOT,
)
from scripts.certify_shimada_disc16_d5_h10_30_vinberg_faces import (
    DiagramEnumerator,
)
from scripts.screen_alexeev_20_2_1_15_9_1_fixed_sections import (
    OUTPUT as SECTION_CENSUS,
    _matrix,
)


OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_literal_vinberg_census.json"
)
EDGE_DIRECTORY = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_20_2_1_15_9_1_vinberg_edges"
)
SCHEMA = "coxeter.alexeev_20_2_1_15_9_1_literal_vinberg_census.v1"


def _mask(indices: tuple[int, ...]) -> int:
    return sum(1 << index for index in indices)


def _indices(mask: int, facet_count: int) -> tuple[int, ...]:
    return tuple(index for index in range(facet_count) if mask & (1 << index))


def _load_rays(path: Path) -> tuple[dict[str, Any], tuple[tuple[Any, ...], ...]]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        header = json.loads(next(stream))
        records = tuple(tuple(json.loads(line)) for line in stream)
    if len(records) != int(header["ray_count"]):
        raise ArithmeticError("exact ray archive is truncated")
    return header, records


def _audit_case(
    section_case: dict[str, Any], cone_case: dict[str, Any]
) -> tuple[dict[str, Any], bytes]:
    case_id = section_case["id"]
    dimension = int(section_case["case"]["hyperbolic_dimension"])
    facet_gram = _matrix(section_case["restricted_mirrors"]["facet_gram"])
    positive_gram = -facet_gram
    facet_count = facet_gram.rows
    if not (
        cone_case["id"] == case_id
        and cone_case["decision"]["finite_volume"] is True
        and cone_case["decision"]["compact"] is False
        and cone_case["hyperbolic_dimension"] == dimension
        and cone_case["facet_count"] == facet_count
        and section_case["coxeter_gate"]["all_pairs_coxeter"] is True
    ):
        raise ArithmeticError(f"pinned input for {case_id} changed")
    ray_path = ROOT / cone_case["artifacts"]["ray_archive_path"]
    ray_header, ray_records = _load_rays(ray_path)
    if ray_header["case_id"] != case_id:
        raise ArithmeticError("ray archive belongs to a different case")
    active_sets = tuple(
        tuple(int(index) - 1 for index in record[1]) for record in ray_records
    )
    kinds = tuple(str(record[3]) for record in ray_records)
    active_masks = tuple(_mask(active) for active in active_sets)
    finite_endpoint_masks = {
        mask for mask, kind in zip(active_masks, kinds, strict=True) if kind == "finite"
    }
    ideal_endpoint_masks = {
        mask for mask, kind in zip(active_masks, kinds, strict=True) if kind == "ideal"
    }
    if any(kind == "hyperideal" for kind in kinds):
        raise ArithmeticError("literal Vinberg input unexpectedly has a hyperideal endpoint")

    enumerator = DiagramEnumerator(positive_gram, dimension=dimension)
    finite_edges = enumerator.finite_subdiagrams(dimension - 1)
    finite_vertices = enumerator.finite_subdiagrams(dimension)
    if set(finite_vertices) != finite_endpoint_masks:
        raise ArithmeticError("elliptic rank-n diagrams do not equal finite endpoints")

    edge_endpoint_pairs: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for right in range(len(active_masks)):
        for left in range(right):
            intersection = active_masks[left] & active_masks[right]
            if intersection.bit_count() >= dimension - 1:
                edge_endpoint_pairs[intersection].append((left, right))
    if not (
        set(edge_endpoint_pairs) == set(finite_edges)
        and all(mask.bit_count() == dimension - 1 for mask in edge_endpoint_pairs)
        and all(len(pairs) == 1 for pairs in edge_endpoint_pairs.values())
    ):
        raise ArithmeticError("an elliptic edge lacks exactly two geometric endpoints")

    affine = enumerator.connected_affine_subdiagrams()
    affine_full_neighborhoods = {}
    for component_mask, record in affine.items():
        component = _indices(component_mask, facet_count)
        zero_neighborhood = tuple(
            facet
            for facet in range(facet_count)
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
            index + 1 for index in _indices(neighborhood_mask, facet_count)
        ]

    cusp_component_order_histogram: Counter[tuple[int, ...]] = Counter()
    cusp_component_count = Counter()
    for cusp_mask in ideal_endpoint_masks:
        components = enumerator._components(cusp_mask, all_nonorthogonal=True)
        if not (
            sum(component.bit_count() - 1 for component in components)
            == dimension - 1
            and all(component in affine for component in components)
            and all(
                affine_full_neighborhoods[component] == cusp_mask
                for component in components
            )
        ):
            raise ArithmeticError("an ideal endpoint is not a full rank-(n-1) cusp")
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
        signature = finite_edges[edge_mask]
        edge_type_histogram[" + ".join(signature)] += 1
        edge_lines.append(
            json.dumps(
                [
                    [index + 1 for index in _indices(edge_mask, facet_count)],
                    [left + 1, right + 1],
                    [kinds[left], kinds[right]],
                    list(signature),
                ],
                separators=(",", ":"),
            )
        )
    edge_header = {
        "schema": "coxeter.alexeev_fixed_section_vinberg_edges.v1",
        "case_id": case_id,
        "record_format": "[edge_facet_positions,endpoint_ray_indices,endpoint_types,finite_diagram_components]",
        "edge_count": len(edge_lines),
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
        for order in range(1, dimension + 3)
    }
    return (
        {
            "id": case_id,
            "ambient_key": section_case["ambient_key"],
            "hyperbolic_dimension": dimension,
            "facet_count": facet_count,
            "method": {
                "criterion": "Vinberg edge-endpoint and full affine zero-neighborhood criterion",
                "diagram_enumeration": "all exact compatible cliques through order n+2, finite component recognition, affine determinant/rank",
                "geometric_face_completeness": "independent GMP exact cone ray/incidence archive",
            },
            "diagram_enumeration": {
                "pair_compatible_clique_count_by_order": compatible_clique_counts,
                "elliptic_rank_n_minus_1_edge_diagram_count": len(finite_edges),
                "elliptic_rank_n_finite_endpoint_diagram_count": len(finite_vertices),
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
                "edge_endpoint_type_histogram": dict(
                    sorted(endpoint_kind_histogram.items())
                ),
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
                "reason": "every elliptic rank-(n-1) edge has exactly two finite/ideal endpoints and every affine component has a full rank-(n-1) parabolic zero-neighborhood",
            },
            "checks": {
                "all_elliptic_edges_equal_exact_geometric_edges": True,
                "all_elliptic_rank_n_diagrams_equal_exact_finite_vertices": True,
                "every_edge_has_exactly_two_endpoints": True,
                "every_affine_full_zero_neighborhood_is_a_certified_cusp": True,
                "every_ideal_vertex_is_componentwise_affine_of_total_rank_n_minus_1": True,
                "ideal_vertices_exist": bool(ideal_endpoint_masks),
            },
        },
        edge_archive,
    )


def build() -> dict[str, Any]:
    sections = json.loads(SECTION_CENSUS.read_text())
    cones = json.loads(CONE_CENSUS.read_text())
    section_by_id = {row["id"]: row for row in sections["cases"]}
    cone_by_id = {row["id"]: row for row in cones["cases"]}
    if set(section_by_id) != set(cone_by_id) or len(section_by_id) != 3:
        raise ArithmeticError("the section and cone censuses do not cover the same cases")
    EDGE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    records = []
    for case_id in sorted(section_by_id):
        record, archive = _audit_case(section_by_id[case_id], cone_by_id[case_id])
        edge_path = EDGE_DIRECTORY / f"{case_id}.jsonl.gz"
        edge_path.write_bytes(archive)
        record["edge_archive_path"] = str(edge_path.relative_to(ROOT))
        records.append(record)
        print(
            json.dumps(
                {
                    "id": case_id,
                    "edges": record["literal_endpoint_audit"]["edge_count"],
                    "finite_vertices": record["literal_endpoint_audit"][
                        "finite_vertex_count"
                    ],
                    "ideal_vertices": record["literal_endpoint_audit"][
                        "ideal_vertex_count"
                    ],
                    "affine_components": record["diagram_enumeration"][
                        "connected_affine_component_count"
                    ],
                    "decision": record["decision"]["classification"],
                }
            ),
            flush=True,
        )
    return {
        "schema": SCHEMA,
        "source": {
            "fixed_section_census": str(SECTION_CENSUS.relative_to(ROOT)),
            "exact_cone_census": str(CONE_CENSUS.relative_to(ROOT)),
        },
        "summary": {
            "case_count": len(records),
            "finite_volume_noncompact_count": sum(
                row["decision"]["classification"] == "finite_volume_noncompact"
                for row in records
            ),
            "bad_edge_count": sum(
                row["literal_endpoint_audit"]["bad_edge_count"] for row in records
            ),
            "bad_parabolic_neighborhood_count": sum(
                row["literal_endpoint_audit"]["bad_parabolic_neighborhood_count"]
                for row in records
            ),
        },
        "cases": records,
    }


def main() -> None:
    payload = build()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
