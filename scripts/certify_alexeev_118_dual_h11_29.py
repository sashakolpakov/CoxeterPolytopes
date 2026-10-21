#!/usr/bin/env python3
"""Certify the dual Alexeev 118-wall constructions of an H11/29 chamber.

The Alexeev ancillary file contains two complete rank-16 mixed-root chambers,
``(16,6,1)`` and ``(16,12,1)``, each with 118 walls.  Their edge-coloured
automorphism groups are too large for NetworkX's generic iterator.  Here a
colour-preserving incidence gadget is passed to Traces/dreadnaut instead.
The resulting complete order-720 groups are lifted entrywise to integral
Lorentz isometries, and GAP enumerates every subgroup conjugacy class.

Each stock has exactly one proper subgroup class with fixed rank at least 11:
an involution class of size 15 and fixed rank 12.  This script pins a
canonical representative, constructs its saturated fixed lattice, and runs
the exact strict-facet/Coxeter gate.  Finite volume is deliberately handled
by the companion literal Vinberg certificate.

Only weighted wall/lattice symmetries are asserted.  No surface automorphism
membership is inferred.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any

import sympy as sp

from scripts.certify_alexeev_high_dim_reflection_examples import _extract_matrix
from scripts.inventory_alexeev_20_2_1_15_9_1_wall_symmetries import (
    SOURCE,
    SOURCE_SHA256,
    _ambient_reconstruction,
    _gap_subgroup_classes,
    _greedy_generators,
    _lift_actions,
)
from scripts.screen_alexeev_20_2_1_15_9_1_fixed_sections import _screen_case
from scripts.search_apery_fermi_chamber_folds import (
    Permutation,
    Subgroup,
    _fixed_rank,
    _generated_subgroup,
    _identity,
    _permutation_order,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "alexeev_118_dual_h11_29_construction.json"
)
SCHEMA = "coxeter.alexeev_118_dual_h11_29.construction.v1"
MIN_FIXED_RANK = 11


CONFIGS: dict[tuple[int, int, int], dict[str, Any]] = {
    (16, 6, 1): {
        "id": "A16_6_1_C01",
        "rank": 16,
        "wall_count": 118,
        "ambient_determinant": -64,
        "dreadnaut_generator_count": 4,
        "fixed_determinant": -256,
        "fixed_smith": [1] * 6 + [2] * 4 + [4] * 2,
        "interior_square": 37418,
        "interior_minimum_pairing": 27,
        "sorted_group_elements_sha256": "02978aec3460d53c6bab1d747d429aa69445c1f77949f38ac3d0a03f1938d9ad",
        "subgroup_fingerprint_sha256": "0fce2a69b6feb1850e41247e0a0196c613ee61f92563d131b9e92463d64256fa",
        "positive_facet_gram_sha256": "ce514e067c7ceed35f0c9d39518054885eca771df23f3bcf8c863522878ac17e",
    },
    (16, 12, 1): {
        "id": "A16_12_1_C01",
        "rank": 16,
        "wall_count": 118,
        "ambient_determinant": -4096,
        "dreadnaut_generator_count": 3,
        "fixed_determinant": -16384,
        "fixed_smith": [2] * 10 + [4] * 2,
        "interior_square": 62866,
        "interior_minimum_pairing": 54,
        "sorted_group_elements_sha256": "3c71760042d12bba61c876aa541d864c5a6f1a8e1af67cfb2318cf33e630c989",
        "subgroup_fingerprint_sha256": "0fce2a69b6feb1850e41247e0a0196c613ee61f92563d131b9e92463d64256fa",
        "positive_facet_gram_sha256": "59b2c4c11937173d8dc9852c4777a4629c0c8d02a15014205fa995bfc72c1eab",
    },
}


EXPECTED_PAIR_COUNTS = {
    "2": 293,
    "3": 28,
    "4": 54,
    "infinity": 19,
    "ultraparallel": 12,
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _dreadnaut_weighted_gram_group(
    gram: sp.Matrix, executable: str
) -> tuple[tuple[Permutation, ...], dict[str, Any]]:
    """Return the complete weighted-Gram group via a coloured pair gadget.

    Original vertices are coloured by norm.  For every unordered pair of
    original vertices, a new vertex is joined to its two endpoints and
    coloured by their Gram entry.  Colour-preserving gadget automorphisms are
    therefore exactly automorphisms of the complete edge-coloured Gram.
    """

    wall_count = gram.rows
    cells: defaultdict[tuple[str, int], list[int]] = defaultdict(list)
    for index in range(wall_count):
        cells[("original", int(gram[index, index]))].append(index)

    pair_rows: list[tuple[int, int, int]] = []
    vertex_count = wall_count
    for right in range(wall_count):
        for left in range(right):
            pair_rows.append((vertex_count, left, right))
            cells[("pair", int(gram[left, right]))].append(vertex_count)
            vertex_count += 1

    adjacency: list[list[int]] = [[] for _ in range(vertex_count)]
    for gadget, left, right in pair_rows:
        adjacency[gadget] = [left, right]
        adjacency[left].append(gadget)
        adjacency[right].append(gadget)

    partition = "|".join(
        ",".join(map(str, cells[colour]))
        for colour in sorted(cells, key=lambda value: (value[0], value[1]))
    )
    rows = ["At", f"n={vertex_count} g"]
    rows.extend(
        " ".join(map(str, neighbours))
        + (";" if index + 1 < vertex_count else ".")
        for index, neighbours in enumerate(adjacency)
    )
    rows.extend([f"f=[{partition}]", "l=100000", "+a +p", "x", "q"])
    completed = subprocess.run(
        [executable],
        input="\n".join(rows) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )

    gadget_generators = []
    for line in completed.stdout.splitlines():
        match = re.fullmatch(r"Gen(?:\(A\))? #\d+:\s*((?:\d+\s*)+)", line.strip())
        if match:
            permutation = tuple(map(int, match.group(1).split()))
            if len(permutation) != vertex_count:
                raise ArithmeticError("dreadnaut emitted a truncated generator")
            gadget_generators.append(permutation)
    summary_match = re.search(
        r"(\d+) orbits?; grpsize=([0-9.e+]+); (\d+) gens", completed.stdout
    )
    if not gadget_generators or summary_match is None:
        raise ArithmeticError("could not parse the dreadnaut certificate")
    if "e" in summary_match.group(2) or "+" in summary_match.group(2):
        raise ArithmeticError("unexpected nonintegral group-order rendering")

    generators = tuple(
        permutation[:wall_count] for permutation in gadget_generators
    )
    if any(
        gram.extract(generator, generator) != gram for generator in generators
    ):
        raise ArithmeticError("a restricted gadget generator changed the wall Gram")
    group = tuple(sorted(_generated_subgroup(generators)))
    reported_order = int(summary_match.group(2))
    if len(group) != reported_order:
        raise ArithmeticError("restricted generators do not recover the reported group")
    canonical_generators = _greedy_generators(group)
    if any(
        gram.extract(generator, generator) != gram
        for generator in canonical_generators
    ):
        raise ArithmeticError("a canonical group generator changed the wall Gram")
    return group, {
        "engine": "Traces through dreadnaut",
        "gadget_vertex_count": vertex_count,
        "gadget_incidence_edge_count": 2 * len(pair_rows),
        "gadget_colour_cell_count": len(cells),
        "gadget_orbit_count": int(summary_match.group(1)),
        "reported_group_order": reported_order,
        "reported_generator_count": int(summary_match.group(3)),
        "canonical_greedy_generators_zero_based": [
            list(row) for row in canonical_generators
        ],
        "raw_dreadnaut_generators_replay_every_gram_entry": True,
        "canonical_generators_replay_every_gram_entry": True,
        "completeness_reason": (
            "original vertices are norm-coloured and every unordered wall "
            "pair has an incidence vertex coloured by its exact Gram entry"
        ),
    }


def _subgroup_fingerprint(
    subgroup: Subgroup,
    class_size: int,
    fixed_rank: int,
) -> tuple[Any, ...]:
    cyclic = any(
        _permutation_order(element) == len(subgroup) for element in subgroup
    )
    histogram = tuple(
        sorted(Counter(_permutation_order(element) for element in subgroup).items())
    )
    return len(subgroup), class_size, fixed_rank, cyclic, histogram


def _gap_structure_description(generators: tuple[Permutation, ...]) -> str:
    gap_generators = ",".join(
        "PermList(" + str([value + 1 for value in generator]) + ")"
        for generator in generators
    )
    program = (
        f"G:=Group([{gap_generators}]);;"
        'Print(Size(G),"|",StructureDescription(G),"\\n");;QUIT;\n'
    )
    completed = subprocess.run(
        ["gap", "-q"], input=program, text=True, capture_output=True, check=True
    )
    if completed.stdout.strip() != "720|S6":
        raise ArithmeticError("GAP no longer identifies the wall group as S6")
    return "S6"


def _fingerprint_payload(fingerprint: tuple[Any, ...], multiplicity: int) -> dict[str, Any]:
    order, class_size, fixed_rank, cyclic, histogram = fingerprint
    return {
        "subgroup_order": order,
        "subgroup_conjugacy_class_size": class_size,
        "fixed_rank": fixed_rank,
        "hyperbolic_dimension": fixed_rank - 1,
        "cyclic": cyclic,
        "element_order_histogram": {str(key): value for key, value in histogram},
        "classes_with_this_fingerprint": multiplicity,
    }


def _inventory(
    group: tuple[Permutation, ...],
    actions: dict[Permutation, sp.Matrix],
    rank: int,
) -> tuple[dict[str, Any], Permutation]:
    generators = _greedy_generators(group)
    structure_description = _gap_structure_description(generators)
    rows = _gap_subgroup_classes(group, generators)
    fingerprints: Counter[tuple[Any, ...]] = Counter()
    high_rows = []
    for row in rows:
        subgroup: Subgroup = row["representative_subgroup"]
        representative_generators = row["representative_generators"]
        matrices = tuple(actions[element] for element in representative_generators)
        fixed_rank = _fixed_rank(matrices, rank) if matrices else rank
        fingerprint = _subgroup_fingerprint(
            subgroup, int(row["conjugacy_class_size"]), fixed_rank
        )
        fingerprints[fingerprint] += 1
        if len(subgroup) > 1 and fixed_rank >= MIN_FIXED_RANK:
            high_rows.append(fingerprint)

    high_involutions = tuple(
        element
        for element in group
        if _permutation_order(element) == 2
        and _fixed_rank((actions[element],), rank) == 12
    )
    if not (
        len(rows) == 56
        and len(high_rows) == 1
        and high_rows[0][:4] == (2, 15, 12, True)
        and len(high_involutions) == 15
    ):
        raise ArithmeticError("the complete high-subgroup denominator changed")
    canonical = min(high_involutions)
    records = [
        _fingerprint_payload(fingerprint, multiplicity)
        for fingerprint, multiplicity in sorted(fingerprints.items(), key=str)
    ]
    return {
        "gap_subgroup_conjugacy_class_count": len(rows),
        "gap_structure_description": structure_description,
        "stable_subgroup_class_fingerprint_histogram": records,
        "stable_fingerprint_histogram_sha256": _payload_sha256(records),
        "proper_fixed_rank_at_least_11_class_count": 1,
        "unique_high_class": {
            "subgroup_order": 2,
            "conjugacy_class_size": 15,
            "fixed_rank": 12,
            "hyperbolic_dimension": 11,
            "cyclic": True,
            "element_order_histogram": {"1": 1, "2": 1},
            "all_fixed_rank_12_involutions_in_group": 15,
        },
        "canonical_representative_rule": (
            "lexicographically least zero-based wall permutation among all "
            "15 order-two elements of fixed rank 12"
        ),
        "canonical_involution_zero_based": list(canonical),
    }, canonical


def _all_wall_sum_interior(
    wall_rows: sp.Matrix,
    ambient_gram: sp.Matrix,
    expected_square: int,
    expected_minimum: int,
) -> sp.Matrix:
    interior = sp.ones(1, wall_rows.rows) * wall_rows
    square = int((interior * ambient_gram * interior.T)[0])
    minimum = int(min(interior * ambient_gram * wall_rows.T))
    if (square, minimum) != (expected_square, expected_minimum):
        raise ArithmeticError("the strict all-wall-sum interior changed")
    return interior


def _audit_case(case: dict[str, Any], config: dict[str, Any]) -> dict[str, bool]:
    checks = {
        "identifier": case["id"] == config["id"],
        "unique_high_class": (
            case["case"]["subgroup_order"] == 2
            and case["case"]["subgroup_conjugacy_class_size"] == 15
            and case["case"]["fixed_rank"] == 12
            and case["case"]["hyperbolic_dimension"] == 11
        ),
        "wall_reduction": (
            case["case"]["source_wall_count"] == 118
            and case["case"]["wall_orbit_count"] == 68
            and case["case"]["wall_orbit_size_histogram"] == {"1": 18, "2": 50}
            and case["exact_reduction"]["conic_facet_count"] == 29
            and case["exact_reduction"]["conic_redundant_count"] == 39
            and case["exact_reduction"]["active_hyperbolic_facet_count"] == 29
            and case["exact_reduction"]["causally_inactive_orbit_indices"] == []
            and case["exact_reduction"]["unresolved_orbit_indices"] == []
        ),
        "fixed_lattice": (
            case["saturated_fixed_lattice"]["determinant"]
            == config["fixed_determinant"]
            and case["saturated_fixed_lattice"]["smith_diagonal_absolute"]
            == config["fixed_smith"]
            and case["saturated_fixed_lattice"]["inertia"] == [1, 11, 0]
        ),
        "strict_integral_coxeter": (
            case["decision"]["classification"]
            == "coxeter_candidate_requires_volume_audit"
            and case["coxeter_gate"]["all_pairs_coxeter"] is True
            and case["coxeter_gate"]["coxeter_edge_counts"]
            == EXPECTED_PAIR_COUNTS
            and case["restricted_mirrors"]
            ["all_active_reflections_integral_on_saturated_fixed_lattice"]
            is True
        ),
        "positive_facet_gram_hash": (
            _payload_sha256(case["restricted_mirrors"]["positive_wall_gram"])
            == config["positive_facet_gram_sha256"]
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"construction replay changed for {case['id']}: {checks}")
    return checks


def _one_case(
    source_text: str,
    key: tuple[int, int, int],
    config: dict[str, Any],
    dreadnaut: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    wall_gram = _extract_matrix(source_text, key)
    rank = int(config["rank"])
    if not (
        wall_gram.shape == (118, 118)
        and wall_gram.rank() == rank
        and wall_gram == wall_gram.T
    ):
        raise ArithmeticError(f"Alexeev matrix {key} changed")
    ambient_gram, wall_rows = _ambient_reconstruction(wall_gram, rank)
    if not (
        int(ambient_gram.det()) == config["ambient_determinant"]
        and exact_inertia(ambient_gram) == (1, 15, 0)
    ):
        raise ArithmeticError("ambient lattice invariants changed")

    group, dreadnaut_record = _dreadnaut_weighted_gram_group(wall_gram, dreadnaut)
    if not (
        len(group) == 720
        and dreadnaut_record["gadget_vertex_count"] == 7021
        and dreadnaut_record["gadget_incidence_edge_count"] == 13806
        and dreadnaut_record["gadget_orbit_count"] == 41
        and dreadnaut_record["reported_generator_count"]
        == config["dreadnaut_generator_count"]
    ):
        raise ArithmeticError("the complete weighted wall group changed")
    actions = _lift_actions(group, wall_rows, ambient_gram)
    inventory, involution = _inventory(group, actions, rank)
    group_hash = _payload_sha256([list(row) for row in group])
    if not (
        group_hash == config["sorted_group_elements_sha256"]
        and inventory["stable_fingerprint_histogram_sha256"]
        == config["subgroup_fingerprint_sha256"]
    ):
        raise ArithmeticError("the frozen group/subgroup inventory changed")
    subgroup = tuple(sorted((_identity(118), involution)))
    example = {
        "ambient": {"rank": rank},
        "weighted_diagram_group": {"elements": [list(row) for row in group]},
    }
    class_record = {
        "class_index": 1,
        "order": 2,
        "conjugacy_class_size": 15,
        "cyclic": True,
        "fixed_rank": 12,
        "representative_elements": [list(row) for row in subgroup],
    }
    interior = _all_wall_sum_interior(
        wall_rows,
        ambient_gram,
        int(config["interior_square"]),
        int(config["interior_minimum_pairing"]),
    )
    case = _screen_case(
        key=key,
        example=example,
        class_record=class_record,
        wall_gram=wall_gram,
        source_interior=interior,
    )
    case["strict_fixed_interior"]["source_method"] = (
        "sum of all 118 primitive wall rows in the leading-root basis"
    )
    case["construction_checks"] = _audit_case(case, config)
    stock = {
        "key": list(key),
        "source_wall_count": 118,
        "ambient_rank": rank,
        "ambient_hyperbolic_dimension": 15,
        "ambient_determinant": int(ambient_gram.det()),
        "ambient_inertia": list(exact_inertia(ambient_gram)),
        "weighted_wall_group": dreadnaut_record
        | {
            "order": len(group),
            "all_720_actions_integral_unimodular_Lorentz": True,
            "sorted_group_elements_sha256": group_hash,
        },
        "subgroup_inventory": inventory,
        "geometric_distinction": {
            "full_weighted_wall_group_used": True,
            "integral_lattice_actions_certified": True,
            "surface_automorphism_membership_claimed": False,
        },
    }
    return stock, case


def build(source: Path = SOURCE, dreadnaut: str | None = None) -> dict[str, Any]:
    if _sha256(source) != SOURCE_SHA256:
        raise ArithmeticError("Alexeev ancillary source hash changed")
    selected = dreadnaut or shutil.which("dreadnaut")
    if selected is None:
        raise RuntimeError("Traces/dreadnaut is required for the 118-wall groups")
    source_text = source.read_text(encoding="utf-8")
    stocks = []
    cases = []
    for key, config in CONFIGS.items():
        print(json.dumps({"starting": list(key)}), flush=True)
        stock, case = _one_case(source_text, key, config, selected)
        stocks.append(stock)
        cases.append(case)
        print(
            json.dumps(
                {
                    "finished": list(key),
                    "dimension": case["case"]["hyperbolic_dimension"],
                    "facets": case["exact_reduction"]["active_hyperbolic_facet_count"],
                }
            ),
            flush=True,
        )
    return {
        "schema": SCHEMA,
        "source": {
            "path": str(source),
            "sha256": SOURCE_SHA256,
            "arxiv": "2209.09110v4",
            "attachment": "all_cases_in_2el_paper.txt",
        },
        "scope": {
            "ambient_keys": [list(key) for key in CONFIGS],
            "complete_weighted_wall_groups": True,
            "all_subgroup_conjugacy_classes": True,
            "minimum_fixed_rank": MIN_FIXED_RANK,
            "finite_volume_claimed": False,
            "surface_automorphism_membership_claimed": False,
        },
        "stocks": stocks,
        "cases": cases,
        "summary": {
            "stock_count": 2,
            "weighted_group_orders": [720, 720],
            "subgroup_conjugacy_class_counts": [56, 56],
            "proper_high_class_counts": [1, 1],
            "candidate_dimensions": [11, 11],
            "candidate_facet_counts": [29, 29],
            "both_strict_integral_coxeter_candidates": True,
        },
        "all_checks_pass": all(
            all(case["construction_checks"].values()) for case in cases
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--dreadnaut")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build(args.source, args.dreadnaut)
    if args.check_archive:
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt H11/29 construction archive changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
