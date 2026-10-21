#!/usr/bin/env python3
"""Freeze the five finite-volume H9 sections parked below the old H10 cutoff."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import shutil
from typing import Any, Sequence

import networkx as nx
import sympy as sp

from scripts.certify_alexeev_20_2_1_15_9_1_exact_cones import (
    _analyze_case,
    _parse_cdd_rays,
    _replay_case,
    _run_cdd,
)
from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
from scripts.explore_alexeev_14_10_1_rank10_c1 import (
    DENOMINATOR,
    _literal_vinberg,
    _time_orientation,
    build as build_section,
)
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    SOURCE,
)


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "scripts" / "explore_alexeev_14_10_1_rank10_c1.py"
FACTORY = ROOT / "data" / "enumerations" / "alexeev_k3_h9_14_factory.json"
DOUBLING = ROOT / "data" / "enumerations" / "alexeev_h9_doubling_cubes.json"
OUTPUT = (
    ROOT / "data" / "enumerations" / "alexeev_parked_rank10_h9_examples.json"
)
INEQUALITY_DIRECTORY = (
    ROOT / "data" / "enumerations" / "alexeev_parked_rank10_h9_inequalities"
)
RAY_DIRECTORY = (
    ROOT / "data" / "enumerations" / "alexeev_parked_rank10_h9_rays"
)

EXPECTED_SHA256 = {
    "denominator": "9c518c02083588cffaa15f38f5887b6f9c2bc109eed7d93611b065a0e0c6fdff",
    "helper": "bb0a921276fddbb4137712c6be78df9abfad2a1f7af96a47a327b8f2fa9bda6a",
    "source": "4cdf39f61f01238a883129551c0181427e5821b7fed2fabfea436cbf3627f457",
    "factory": "b1ee6fd769527ae6dab01a2bb58099e7c7d1776d0417b6a8d710316697e8a775",
    "doubling": "4919728a5c5ff8adebd80b250cfedfb562a040c8bd1e72e551b398b1f13e7e4c",
}

SPECS: tuple[dict[str, Any], ...] = (
    {
        "key": (14, 8, 0),
        "class_index": 1,
        "id": "A14_8_0_C01",
        "raw": 12,
        "redundant": 0,
        "facets": 12,
        "determinant": -1024,
        "smith": [1, 1, 2, 2, 2, 2, 2, 2, 4, 4],
        "relations": {"2": 54, "3": 8, "4": 4},
        "finite": 44,
        "ideal": 5,
        "edges": 260,
        "affine": 10,
        "time_min": "9",
        "time_max": "70",
        "source_norms": {"-4": 7, "-2": 5},
        "folded_norms": {"-8": 3, "-4": 6, "-2": 3},
        "doubling_walls": [],
    },
    {
        "key": (14, 8, 1),
        "class_index": 2,
        "id": "A14_8_1_C02",
        "raw": 17,
        "redundant": 1,
        "facets": 16,
        "determinant": -1024,
        "smith": [1, 1, 2, 2, 2, 2, 2, 2, 4, 4],
        "relations": {"2": 96, "3": 9, "4": 10, "parallel": 5},
        "finite": 165,
        "ideal": 23,
        "edges": 1112,
        "affine": 55,
        "time_min": "21",
        "time_max": "294",
        "source_norms": {"-4": 7, "-2": 9},
        "folded_norms": {"-8": 3, "-4": 7, "-2": 6},
        "doubling_walls": [7, 15, 16],
    },
    {
        "key": (14, 10, 1),
        "class_index": 1,
        "id": "A14_10_1_C01",
        "raw": 14,
        "redundant": 1,
        "facets": 13,
        "determinant": -4096,
        "smith": [2, 2, 2, 2, 2, 2, 2, 2, 4, 4],
        "relations": {"2": 63, "3": 9, "4": 6},
        "finite": 70,
        "ideal": 7,
        "edges": 393,
        "affine": 14,
        "time_min": "12",
        "time_max": "100",
        "source_norms": {"-4": 10, "-2": 3},
        "folded_norms": {"-8": 3, "-4": 10},
        "doubling_walls": [3, 5, 11],
    },
    {
        "key": (14, 10, 1),
        "class_index": 2,
        "id": "A14_10_1_C02",
        "raw": 15,
        "redundant": 0,
        "facets": 15,
        "determinant": -4096,
        "smith": [2, 2, 2, 2, 2, 2, 2, 2, 4, 4],
        "relations": {"2": 85, "3": 8, "4": 9, "parallel": 3},
        "finite": 123,
        "ideal": 15,
        "edges": 778,
        "affine": 35,
        "time_min": "11",
        "time_max": "111",
        "source_norms": {"-4": 9, "-2": 6},
        "folded_norms": {"-8": 4, "-4": 7, "-2": 4},
        "doubling_walls": [1, 3, 6, 15],
    },
    {
        "key": (16, 14, 1),
        "class_index": 1,
        "id": "A16_14_1_C01",
        "raw": 18,
        "redundant": 3,
        "facets": 15,
        "determinant": -16384,
        "smith": [2, 2, 2, 2, 2, 2, 4, 4, 4, 4],
        "relations": {"2": 85, "3": 8, "4": 6, "parallel": 6},
        "finite": 100,
        "ideal": 9,
        "edges": 588,
        "affine": 23,
        "time_min": "36",
        "time_max": "304",
        "source_norms": {"-4": 11, "-2": 4},
        "folded_norms": {"-8": 7, "-4": 4, "-2": 4},
        "doubling_walls": [1, 8, 14, 15],
    },
)


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _matrix(payload: Sequence[Sequence[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _relation_token(gram: sp.Matrix, left: int, right: int) -> str:
    ratio = sp.cancel(
        gram[left, right] ** 2 / (gram[left, left] * gram[right, right])
    )
    lookup = {
        sp.Rational(0): "2",
        sp.Rational(1, 4): "3",
        sp.Rational(1, 2): "4",
        sp.Rational(3, 4): "6",
        sp.Rational(1): "parallel",
    }
    if ratio not in lookup:
        raise ArithmeticError(f"unexpected normalized Coxeter product {ratio}")
    return lookup[ratio]


def _relations(gram: sp.Matrix) -> tuple[tuple[str, ...], ...]:
    return tuple(
        tuple(
            "1" if left == right else _relation_token(gram, left, right)
            for right in range(gram.rows)
        )
        for left in range(gram.rows)
    )


def _relation_histogram(relations: Sequence[Sequence[str]]) -> dict[str, int]:
    return dict(
        sorted(
            Counter(
                relations[left][right]
                for left in range(len(relations))
                for right in range(left)
            ).items()
        )
    )


def _graph(relations: Sequence[Sequence[str]]) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(range(len(relations)))
    for left in range(len(relations)):
        for right in range(left):
            if relations[left][right] != "2":
                graph.add_edge(left, right, label=relations[left][right])
    return graph


def _isomorphic(
    left: Sequence[Sequence[str]], right: Sequence[Sequence[str]]
) -> bool:
    return nx.is_isomorphic(
        _graph(left),
        _graph(right),
        edge_match=nx.algorithms.isomorphism.categorical_edge_match("label", None),
    )


def _prior_relations() -> list[tuple[str, tuple[tuple[str, ...], ...]]]:
    factory = json.loads(FACTORY.read_text(encoding="utf-8"))
    doubling = json.loads(DOUBLING.read_text(encoding="utf-8"))
    rows: list[tuple[str, tuple[tuple[str, ...], ...]]] = []
    for record in factory["constructions"]:
        gram = _matrix(record["section"]["restricted_mirrors"]["facet_gram"])
        rows.append((record["id"], _relations(gram)))
    for record in doubling["certified_types"]:
        rows.append((record["id"], _relations(_matrix(record["facet_gram"]))))
    return rows


def build(executable: str | None = None) -> tuple[dict[str, Any], dict[Path, bytes]]:
    actual_hashes = {
        "denominator": _sha(DENOMINATOR),
        "helper": _sha(HELPER),
        "source": _sha(SOURCE),
        "factory": _sha(FACTORY),
        "doubling": _sha(DOUBLING),
    }
    if actual_hashes != EXPECTED_SHA256:
        raise ArithmeticError("a pinned parked-rank-ten source changed")
    denominator = json.loads(DENOMINATOR.read_text(encoding="utf-8"))
    rank_ten_rows = sorted(
        (tuple(example["key"]), int(row["class_index"]))
        for example in denominator["examples"]
        for row in example["subgroup_classes"]
        if int(row["fixed_rank"]) == 10
    )
    expected_denominator = sorted(
        [
            ((14, 8, 0), 1),
            ((14, 8, 1), 1),
            ((14, 8, 1), 2),
            ((14, 10, 1), 1),
            ((14, 10, 1), 2),
            ((16, 14, 1), 1),
        ]
    )
    if rank_ten_rows != expected_denominator:
        raise ArithmeticError("the parked rank-ten subgroup denominator changed")

    cdd = executable or shutil.which("cddexec_gmp") or "cddexec_gmp"
    artifacts: dict[Path, bytes] = {}
    records: list[dict[str, Any]] = []
    new_relations: list[tuple[tuple[str, ...], ...]] = []
    source_text = SOURCE.read_text(encoding="utf-8")
    for spec in SPECS:
        section = build_section(spec["key"], spec["class_index"])
        reduction = section["exact_reduction"]
        fixed = section["saturated_fixed_lattice"]
        facet_gram = _matrix(section["restricted_mirrors"]["facet_gram"])
        relations = _relations(facet_gram)
        relation_histogram = _relation_histogram(relations)
        doubling_walls = [
            index + 1
            for index in range(facet_gram.rows)
            if all(
                relations[index][other] not in {"3", "5", "7", "9"}
                for other in range(facet_gram.rows)
                if other != index
            )
        ]

        source_wall_gram = _extract_matrix(source_text, spec["key"])
        active = [
            row
            for row in section["restricted_mirrors"]["records"]
            if row["hyperbolic_status"] == "active_facet"
        ]
        source_norms: Counter[str] = Counter()
        for row in active:
            orbit_norms = {
                str(source_wall_gram[index - 1, index - 1])
                for index in row["source_wall_indices"]
            }
            if len(orbit_norms) != 1:
                raise ArithmeticError("a wall orbit mixes source norms")
            source_norms[orbit_norms.pop()] += 1
        folded_norms = Counter(str(row["square"]) for row in active)

        if not (
            section["id"] == spec["id"]
            and all(section["checks"].values())
            and section["case"]["fixed_rank"] == 10
            and section["case"]["hyperbolic_dimension"] == 9
            and section["case"]["wall_orbit_count"] == spec["raw"]
            and reduction["conic_redundant_count"] == spec["redundant"]
            and reduction["active_hyperbolic_facet_count"] == spec["facets"]
            and reduction["causally_inactive_orbit_indices"] == []
            and reduction["unresolved_orbit_indices"] == []
            and fixed["determinant"] == spec["determinant"]
            and fixed["smith_diagonal_absolute"] == spec["smith"]
            and relation_histogram == spec["relations"]
            and dict(source_norms) == spec["source_norms"]
            and dict(folded_norms) == spec["folded_norms"]
            and doubling_walls == spec["doubling_walls"]
            and section["restricted_mirrors"][
                "all_active_reflections_integral_on_saturated_fixed_lattice"
            ]
            and section["coxeter_gate"]["non_coxeter_pair_count"] == 0
            and len(reduction["strict_facet_witnesses"]) == spec["facets"]
        ):
            raise ArithmeticError(f"the exact section gate failed for {spec['id']}")

        cone = _replay_case(section)
        cdd_output = _run_cdd(cone["inequality_text"], cdd)
        rays = _parse_cdd_rays(cdd_output, int(cone["rank"]))
        volume, ray_archive = _analyze_case(section, cone, rays)
        time_orientation = _time_orientation(section, cone, rays)
        vinberg = _literal_vinberg(section, cone, rays)
        exact = volume["exact_double_description"]
        if not (
            volume["decision"]["classification"] == "finite_volume_noncompact"
            and exact["finite_vertex_count"] == spec["finite"]
            and exact["ideal_vertex_count"] == spec["ideal"]
            and exact["hyperideal_vertex_count"] == 0
            and time_orientation["minimum_strict_interior_pairing"]
            == spec["time_min"]
            and time_orientation["maximum_strict_interior_pairing"]
            == spec["time_max"]
            and vinberg["elliptic_edge_count"] == spec["edges"]
            and vinberg["connected_affine_component_count"] == spec["affine"]
            and vinberg["bad_edge_count"] == 0
            and vinberg["bad_parabolic_neighborhood_count"] == 0
        ):
            raise ArithmeticError(f"the exact volume gate failed for {spec['id']}")

        inequality_path = INEQUALITY_DIRECTORY / f"{spec['id']}.ine"
        ray_path = RAY_DIRECTORY / f"{spec['id']}.jsonl.gz"
        inequality_bytes = cone["inequality_text"].encode()
        artifacts[inequality_path] = inequality_bytes
        artifacts[ray_path] = ray_archive
        record = {
            "id": spec["id"],
            "ambient_key": list(spec["key"]),
            "class_index": spec["class_index"],
            "section": section,
            "relation_histogram": relation_histogram,
            "doubling_walls_one_based": doubling_walls,
            "source_wall_audit": {
                "active_facet_source_norm_histogram": dict(source_norms),
                "folded_primitive_facet_norm_histogram": dict(folded_norms),
            },
            "exact_cone": volume,
            "time_orientation": time_orientation,
            "literal_vinberg": vinberg,
            "artifacts": {
                "inequality_path": str(inequality_path.relative_to(ROOT)),
                "inequality_sha256": sha256(inequality_bytes).hexdigest(),
                "ray_archive_path": str(ray_path.relative_to(ROOT)),
                "ray_archive_sha256": sha256(ray_archive).hexdigest(),
            },
        }
        records.append(record)
        new_relations.append(relations)

    pairwise = {
        f"{records[left]['id']}__{records[right]['id']}": _isomorphic(
            new_relations[left], new_relations[right]
        )
        for left in range(len(records))
        for right in range(left + 1, len(records))
    }
    prior = _prior_relations()
    prior_matches = {
        records[index]["id"]: [
            prior_id
            for prior_id, prior_relation in prior
            if len(prior_relation) == len(new_relations[index])
            and _isomorphic(new_relations[index], prior_relation)
        ]
        for index in range(len(records))
    }
    if any(pairwise.values()) or any(prior_matches.values()):
        raise ArithmeticError("a parked H9 diagram duplicates this project batch")

    payload: dict[str, Any] = {
        "schema": "coxeter.alexeev_parked_rank10_h9_examples.v1",
        "source": {
            "file_sha256": actual_hashes,
            "complete_parked_rank10_subgroup_rows": [
                [list(key), class_index] for key, class_index in rank_ten_rows
            ],
            "previously_frozen_row": [[14, 8, 1], 1],
        },
        "scope": {
            "five_abstract_hyperbolic_coxeter_polytopes_certified": True,
            "finite_volume_and_noncompactness_certified": True,
            "integral_fixed_lattice_reflections_certified": True,
            "weak_full_picard_k3_divisor_carriers_certified": True,
            "surface_automorphism_or_nef_claimed": False,
            "worldwide_novelty_claimed": False,
        },
        "cases": records,
        "pairwise_relation_diagram_isomorphisms": pairwise,
        "matches_with_prior_project_H9_factory_or_cube_diagrams": prior_matches,
        "summary": {
            "new_project_ledger_diagram_count": 5,
            "finite_volume_noncompact_count": 5,
            "compact_count": 0,
            "facet_count_histogram": {"12": 1, "13": 1, "15": 2, "16": 1},
            "total_bad_edge_count": 0,
            "total_bad_parabolic_neighborhood_count": 0,
        },
        "decision": {
            "classification": "five_pairwise_nonisomorphic_finite_volume_noncompact_H9_fixed_sections",
            "lineage_status": "worldwide_priority_pending_independent_audit",
        },
    }
    digest_payload = json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode()
    payload["certificate_payload_sha256"] = sha256(digest_payload).hexdigest()
    return payload, artifacts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-archive", action="store_true")
    parser.add_argument("--cdd", default=None)
    args = parser.parse_args()
    payload, artifacts = build(args.cdd)
    rendered = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    if args.check_archive:
        if OUTPUT.read_bytes() != rendered:
            raise ArithmeticError("parked-rank-ten H9 JSON archive changed")
        for path, content in artifacts.items():
            if path.read_bytes() != content:
                raise ArithmeticError(f"parked-rank-ten artifact changed: {path}")
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(rendered)
        for path, content in artifacts.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
    print(json.dumps(payload["summary"], sort_keys=True))
    print(OUTPUT)


if __name__ == "__main__":
    main()
