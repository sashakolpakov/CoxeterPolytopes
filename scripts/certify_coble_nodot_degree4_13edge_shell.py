#!/usr/bin/env python3
"""Certify the exact degree-four, thirteen-edge nodot shell.

This enlarges the frozen maximum-degree-three denominator in one intrinsic
polytope direction: connected twelve-node supports have exactly thirteen
ordinary edges and maximum degree exactly four.  From zero through three
ordinary edges may be promoted from root pairing one to pairing two.  Every
presentation with exact root-Gram inertia ``(1,9,2)`` is retained, without a
discriminant or surface filter, deduplicated as a weighted diagram, and sent
through the exact quotient-cone and literal Vinberg tests.

The result is deliberately a surface-free negative shell theorem.  Every
Lorentzian class has a hyperideal extreme ray, so none has finite hyperbolic
volume and none reaches an algebraic carrier gate.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import shutil
from typing import Any

import networkx as nx

from scripts.certify_coble_nodot_13edge_lane import (
    _analyze_cone,
    _deduplicate_weighted,
    _literal_vinberg,
    _matrix,
    _quotient_lattice,
)
from scripts.certify_coble_nodot_maxdegree3_closure import (
    _perron_strict_interior_certificate,
)
from scripts.scout_coble_nodot_weighted_graphs import build as build_scout


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "coble_nodot_degree4_13edge_shell.json"
)
SCHEMA = "coxeter.coble_nodot_degree4_13edge_shell.v1"

PINNED = {
    ROOT / "data/enumerations/coble_nodot_13edge_lane.json":
        "103f4867ef89bca47f96d0ee806ae4c6fa157a62058294ed5852b112bfad6ffa",
    ROOT / "data/enumerations/coble_nodot_maxdegree3_closure.json":
        "56a0646809fff76514729a58887172eb916d28cbda699128ed2f96a23556e82d",
}

EXPECTED_MARKED = {0: 1, 1: 26, 2: 189, 3: 422}
EXPECTED_CLASSES = {0: 1, 1: 24, 2: 145, 3: 330}
EXPECTED_MAX_DEGREE_AT_MOST_FOUR_SUPPORTS = 16_977
EXPECTED_MAX_DEGREE_AT_MOST_THREE_SUPPORTS = 2_888
EXPECTED_DEGREE_FOUR_SUPPORTS = 14_089


def _digest(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _load_pinned() -> dict[Path, dict[str, Any]]:
    loaded: dict[Path, dict[str, Any]] = {}
    for path, expected in PINNED.items():
        if sha256(path.read_bytes()).hexdigest() != expected:
            raise ArithmeticError(f"pinned archive changed: {path}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not payload.get("all_checks_pass"):
            raise ArithmeticError(f"pinned archive is not green: {path}")
        loaded[path] = payload
    return loaded


def _support_maximum_degree(graph6: str) -> int:
    graph = nx.from_graph6_bytes(graph6.encode())
    return max(dict(graph.degree()).values())


def _class_record(
    pairing_two_count: int,
    class_index: int,
    diagram_class: dict[str, Any],
    cdd: str,
) -> dict[str, Any]:
    representative = diagram_class["representative"]
    root_gram = -_matrix(representative["negative_cartan"])
    quotient = _quotient_lattice(root_gram)
    cone = _analyze_cone(
        quotient["gram"], quotient["root_coordinates"], cdd
    )
    interior = _perron_strict_interior_certificate(
        quotient["gram"], quotient["root_coordinates"]
    )
    vinberg = _literal_vinberg(
        quotient["gram"], quotient["root_coordinates"], cone
    )
    hyperideal_records = [
        record for record in cone["rays"] if record["type"] == "hyperideal"
    ]
    if not (
        cone["all_input_walls_irredundant"]
        and cone["all_causal_rays_same_time_cone"]
        and cone["hyperideal_vertex_count"] == len(hyperideal_records) > 0
        and interior["all_wall_slacks_strictly_positive"]
        and interior["square"] > 0
        and not vinberg["literal_vinberg_passes"]
        and vinberg["stopped_at_hyperideal_endpoint"]
    ):
        raise ArithmeticError(
            f"degree-four class k={pairing_two_count}, i={class_index} "
            "failed an exact chamber or hyperideal gate"
        )
    return {
        "pairing_two_count": pairing_two_count,
        "class_index_within_pairing_count_zero_based": class_index,
        "labelled_member_count": len(diagram_class["members"]),
        "support_graph6": representative["underlying_graph6"],
        "promoted_pairing_two_edges_zero_based": representative[
            "promoted_parallel_edges_zero_based"
        ],
        "negative_cartan": representative["negative_cartan"],
        "radical_quotient_determinant": quotient["determinant"],
        "finite_ideal_hyperideal": [
            cone["finite_vertex_count"],
            cone["ideal_vertex_count"],
            cone["hyperideal_vertex_count"],
        ],
        "classification": "infinite_volume_hyperideal",
        "all_input_walls_irredundant": cone[
            "all_input_walls_irredundant"
        ],
        "positive_strict_interior": interior,
        "literal_vinberg": vinberg,
        "one_exact_hyperideal_ray_witness": hyperideal_records[0],
        "complete_ray_record_sha256": _digest(cone["rays"]),
        "hyperideal_ray_records_sha256": _digest(hyperideal_records),
    }


def build(cdd: str | None = None) -> dict[str, Any]:
    cdd = cdd or shutil.which("cddexec_gmp") or "cddexec_gmp"
    pinned = _load_pinned()
    prior_13 = pinned[ROOT / "data/enumerations/coble_nodot_13edge_lane.json"]
    subcubic_support_count = prior_13["scout_denominator"][
        "underlying_unlabelled_graph_count"
    ]
    if subcubic_support_count != EXPECTED_MAX_DEGREE_AT_MOST_THREE_SUPPORTS:
        raise ArithmeticError("pinned subcubic thirteen-edge support count changed")

    counts: dict[str, dict[str, Any]] = {}
    records: list[dict[str, Any]] = []
    determinant_histogram: Counter[int] = Counter()
    finite_ideal_hyperideal_totals = [0, 0, 0]

    for pairing_two_count in range(4):
        scout = build_scout(
            maximum_degree=4,
            edge_range="13:13",
            maximum_parallel=pairing_two_count,
            minimum_parallel=pairing_two_count,
            sample_limit=100_000,
            sample_all_exact=True,
        )
        if scout["underlying_graph_count"] != (
            EXPECTED_MAX_DEGREE_AT_MOST_FOUR_SUPPORTS
        ):
            raise ArithmeticError("maximum-degree-four support count changed")
        if len(scout["exact_labeled_survivor_samples"]) != scout[
            "exact_labeled_survivor_count"
        ]:
            raise ArithmeticError("the exact scout sample was truncated")
        marked = [
            sample
            for sample in scout["exact_labeled_survivor_samples"]
            if _support_maximum_degree(sample["underlying_graph6"]) == 4
        ]
        classes = _deduplicate_weighted(marked)
        if len(marked) != EXPECTED_MARKED[pairing_two_count]:
            raise ArithmeticError(
                f"degree-four marked count changed at k={pairing_two_count}"
            )
        if len(classes) != EXPECTED_CLASSES[pairing_two_count]:
            raise ArithmeticError(
                f"degree-four weighted class count changed at "
                f"k={pairing_two_count}"
            )
        trials = EXPECTED_DEGREE_FOUR_SUPPORTS * comb(13, pairing_two_count)
        counts[str(pairing_two_count)] = {
            "weighted_trials": trials,
            "exact_Lorentzian_marked_presentations": len(marked),
            "weighted_diagram_isometry_classes": len(classes),
        }
        for class_index, diagram_class in enumerate(classes):
            record = _class_record(
                pairing_two_count, class_index, diagram_class, cdd
            )
            records.append(record)
            determinant_histogram[record["radical_quotient_determinant"]] += 1
            for index, value in enumerate(record["finite_ideal_hyperideal"]):
                finite_ideal_hyperideal_totals[index] += value

    if len(records) != sum(EXPECTED_CLASSES.values()):
        raise ArithmeticError("combined degree-four class count changed")
    if any(record["classification"] != "infinite_volume_hyperideal" for record in records):
        raise ArithmeticError("a degree-four finite-volume survivor appeared")

    exact_census = {
        "maximum_degree_at_most_four_support_count": (
            EXPECTED_MAX_DEGREE_AT_MOST_FOUR_SUPPORTS
        ),
        "maximum_degree_at_most_three_support_count": (
            EXPECTED_MAX_DEGREE_AT_MOST_THREE_SUPPORTS
        ),
        "maximum_degree_exactly_four_support_count": (
            EXPECTED_DEGREE_FOUR_SUPPORTS
        ),
        "weighted_trials": sum(
            row["weighted_trials"] for row in counts.values()
        ),
        "exact_Lorentzian_marked_presentations": sum(EXPECTED_MARKED.values()),
        "weighted_diagram_isometry_classes": len(records),
        "finite_volume_classes": 0,
        "infinite_volume_hyperideal_classes": len(records),
    }
    expected_census = {
        "maximum_degree_at_most_four_support_count": 16_977,
        "maximum_degree_at_most_three_support_count": 2_888,
        "maximum_degree_exactly_four_support_count": 14_089,
        "weighted_trials": 5_325_642,
        "exact_Lorentzian_marked_presentations": 638,
        "weighted_diagram_isometry_classes": 500,
        "finite_volume_classes": 0,
        "infinite_volume_hyperideal_classes": 500,
    }
    if exact_census != expected_census:
        raise ArithmeticError(f"degree-four exact census changed: {exact_census}")

    result = {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "intrinsic_surface_free_shell": {
            "name": "S_Delta=4,e=13,kappa<=3",
            "definition": (
                "twelve norm-minus-two wall symbols; connected ordinary "
                "support with exactly thirteen edges and maximum degree "
                "exactly four; off-diagonal root pairings in {0,1,2}; at "
                "most three pairing-two edges; exact root-Gram inertia "
                "(1,9,2); no determinant, E10, surface, or carrier filter"
            ),
            "membership_is_intrinsic_to_weighted_Coxeter_Gram_data": True,
            "surface_or_carrier_success_in_membership": False,
            "relation_to_previous_denominator": (
                "the maximum-degree bound is enlarged from at most three to "
                "an exact degree-four shell while the Type-I thirteen-edge "
                "support count and pairing bound are held fixed"
            ),
        },
        "counts_by_pairing_two_count": counts,
        "exact_census": exact_census,
        "quotient_determinant_histogram": {
            str(key): value for key, value in sorted(determinant_histogram.items())
        },
        "extreme_ray_type_totals_finite_ideal_hyperideal": (
            finite_ideal_hyperideal_totals
        ),
        "class_records": records,
        "theorem": {
            "statement": (
                "No weighted diagram in S_Delta=4,e=13,kappa<=3 defines a "
                "finite-volume hyperbolic root chamber: each of its 500 exact "
                "Lorentzian weighted-diagram classes has an exact hyperideal "
                "extreme ray in the saturated radical quotient."
            ),
            "vinberg_volume_test_applied_to_every_class": True,
            "compact_classes": 0,
            "finite_volume_noncompact_classes": 0,
            "infinite_volume_hyperideal_classes": 500,
            "new_hyperbolic_finite_volume_examples": [],
        },
        "polytope_to_surface_transition_audit": {
            "chain_attempted": (
                "intrinsic weighted Gram shell -> exact Lorentzian quotient "
                "cone -> Vinberg volume gate -> carrier enumeration -> "
                "independent algebraic geometry -> pushback to P"
            ),
            "last_completed_gate": "exact Vinberg hyperideal-ray exclusion",
            "carrier_enumeration_reached": False,
            "algebraic_geometry_invoked": False,
            "exact_same_P_return_reached": False,
            "independent_AG_pushback_reached": False,
            "strict_L4_flagship_gate_crossed": False,
            "rating": "SCOPED_P_SIDE_NEGATIVE_L2_L3_NOT_L4",
            "honest_value": (
                "a certified search boundary: this natural degree-four shell "
                "cannot contain the desired new finite-volume polytope, but "
                "the exclusion is already visible on the P side and therefore "
                "does not constitute polytope-to-surface theoretical leverage"
            ),
        },
        "next_completion_gate": (
            "Enlarge a different intrinsic shell until at least one exact "
            "finite-volume class survives Vinberg; then exhaust all relevant "
            "surface carriers and require an independent AG theorem to decide "
            "a nonthin subset before returning exactly to the original P."
        ),
        "pinned_sources": {
            str(path.relative_to(ROOT)): {"sha256": digest}
            for path, digest in PINNED.items()
        },
        "all_checks_pass": True,
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cdd", default=None)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build(args.cdd)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check_archive:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit("archive mismatch")
        print(sha256(rendered.encode()).hexdigest())
    elif args.write:
        OUTPUT.write_text(rendered, encoding="utf-8")
        print(OUTPUT)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
