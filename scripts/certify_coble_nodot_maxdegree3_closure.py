#!/usr/bin/env python3
"""Close the maximum-degree-three one-boundary nodot denominator.

The new computation exhausts the two tails (11--13 and 16--18 support
edges) for twelve roots, at most three pairing-two edges, and exact inertia
(9,1,2).  The already frozen 14- and 15-edge certificates are pinned and
combined with the tail result.  The same artifact also records the classical
P10 counterexample to automatic Coble incidence on all of D_CobVin.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil
from typing import Any

import sympy as sp

from scripts.certify_coble_nodot_13edge_lane import (
    _analyze_cone,
    _deduplicate_weighted,
    _known_root_graphs,
    _literal_vinberg,
    _matrix,
    _quotient_lattice,
    _weighted_isomorphic,
)
from scripts.scout_coble_nodot_weighted_graphs import build as build_scout


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "coble_nodot_maxdegree3_closure.json"
)
SCHEMA = "coxeter.coble_nodot_maxdegree3_closure.v1"

PINNED = {
    ROOT / "data/enumerations/coble_nodot_13edge_lane.json":
        "103f4867ef89bca47f96d0ee806ae4c6fa157a62058294ed5852b112bfad6ffa",
    ROOT / "data/enumerations/coble_nodot_14edge_lane.json":
        "e1a14bfed652de62979baf61afaabf262bd96482e0ec7e3376994267884a03aa",
    ROOT / "data/enumerations/coble_nodot_15edge_lane.json":
        "c9e2da6bcc08736fbc9848d64ca3fd273c23d535394ef727b671f369e7b20dae",
    ROOT / "data/enumerations/coble_nodot_15edge_one_pair2_lane.json":
        "1e3e20026c418c6857ceb7a834693fd9217479e05506ddcab93d16304c38b084",
    ROOT / "data/enumerations/coble_nodot_15edge_two_pair2_lane.json":
        "85a2a347b84954b20e3f796958a244755fd994a76508fbfa3bad84f077c15167",
    ROOT / "data/enumerations/coble_nodot_15edge_three_pair2_lane.json":
        "581480510c340cb77ad21f017aa48ace1e072d0f07ba35c6f54eb895a5b6be2c",
    ROOT / "data/enumerations/coble_one_boundary_raw_k3_promotion.json":
        "1c8b8699da8f8747fdad22d520743bb4640b41b8e988c16844ecfe5e27c1958c",
    ROOT / "data/enumerations/polytope_first_coble_h10_hostile_audit.json":
        "12c23f9e1818f962175232a55d7497c9d18f63554db7a571539e8e8d6ce9d0f7",
    ROOT / "data/enumerations/coble_vinberg_incidence_inverse_transfer.json":
        "ee4280e2c9cbee2a730f2333c63566f8cc36bf7dab125f9b168ce437260028f8",
}


def _digest(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _load_pinned() -> dict[Path, dict[str, Any]]:
    loaded = {}
    for path, expected in PINNED.items():
        if sha256(path.read_bytes()).hexdigest() != expected:
            raise ArithmeticError(f"pinned archive changed: {path}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not payload.get("all_checks_pass"):
            raise ArithmeticError(f"pinned archive is not green: {path}")
        loaded[path] = payload
    return loaded


def _perron_strict_interior_certificate(
    gram: sp.Matrix,
    normal_coordinates: sp.Matrix,
) -> dict[str, Any]:
    """Give an exact timelike interior using Perron iteration.

    If W is the connected acute wall Gram, A=W+2I is its nonnegative
    weighted adjacency matrix.  The Perron eigenvalue of A is greater than
    two because W has one positive direction.  Iteration by A+I is primitive
    and converges projectively to its positive Perron vector.  Hence, after
    finitely many exact integer steps, W*a is entrywise positive and
    a^t*W*a is positive.  Then x=N*a is a strict timelike chamber point.
    """

    wall_gram = normal_coordinates.T * gram * normal_coordinates
    size = wall_gram.rows
    adjacency_with_loops = wall_gram + 3 * sp.eye(size)
    coefficients = sp.ones(size, 1)
    for iteration in range(512):
        coordinates = normal_coordinates * coefficients
        slacks = wall_gram * coefficients
        square = int((coefficients.T * wall_gram * coefficients)[0])
        if square > 0 and all(value > 0 for value in slacks):
            return {
                "construction": "exact_primitive_Perron_iteration_on_W_plus_3I",
                "iteration_count": iteration,
                "positive_wall_coefficients": [int(value) for value in coefficients],
                "coordinates": [int(value) for value in coordinates],
                "square": square,
                "all_wall_slacks_strictly_positive": True,
            }
        coefficients = adjacency_with_loops * coefficients
        divisor = 0
        for value in coefficients:
            divisor = math.gcd(divisor, abs(int(value)))
        if divisor > 1:
            coefficients /= divisor
    raise ArithmeticError("exact Perron iteration did not reach a timelike interior")


def _tail(
    edge_range: str,
    expected_marked: dict[int, int],
    expected_classes: dict[int, int],
    cdd: str,
) -> dict[str, Any]:
    known = _known_root_graphs()
    records = []
    counts = {}
    volume = Counter()
    finite = []
    trials = 0
    marked_total = 0
    classes_total = 0
    for pairing_two_count in range(4):
        scout = build_scout(
            maximum_degree=3,
            edge_range=edge_range,
            maximum_parallel=pairing_two_count,
            minimum_parallel=pairing_two_count,
            sample_limit=10_000,
            sample_all_exact=True,
        )
        marked = scout["exact_labeled_survivor_samples"]
        classes = _deduplicate_weighted(marked)
        if len(marked) != expected_marked[pairing_two_count]:
            raise ArithmeticError(
                f"{edge_range}, k={pairing_two_count}: marked count changed"
            )
        if len(classes) != expected_classes[pairing_two_count]:
            raise ArithmeticError(
                f"{edge_range}, k={pairing_two_count}: class count changed"
            )
        trials_here = sum(scout["trial_count_by_edges_and_parallel"].values())
        trials += trials_here
        marked_total += len(marked)
        classes_total += len(classes)
        counts[str(pairing_two_count)] = {
            "weighted_trials": trials_here,
            "exact_Lorentzian_marked_presentations": len(marked),
            "weighted_isometry_classes": len(classes),
            "marked_by_edge_count": scout[
                "exact_labeled_survivor_count_by_edges_and_parallel"
            ],
        }
        for class_index, diagram_class in enumerate(classes):
            representative = diagram_class["representative"]
            root_gram = -_matrix(representative["negative_cartan"])
            quotient = _quotient_lattice(root_gram)
            cone = _analyze_cone(
                quotient["gram"], quotient["root_coordinates"], cdd
            )
            interior = _perron_strict_interior_certificate(
                quotient["gram"], quotient["root_coordinates"]
            )
            if not (
                cone["all_input_walls_irredundant"]
                and cone["all_causal_rays_same_time_cone"]
                and interior["all_wall_slacks_strictly_positive"]
                and interior["square"] > 0
            ):
                raise ArithmeticError("an exact root cone failed its chamber gate")
            hyperideal = cone["hyperideal_vertex_count"]
            if hyperideal:
                classification = "infinite_volume_hyperideal"
                vinberg_pass = False
            else:
                vinberg = _literal_vinberg(
                    quotient["gram"], quotient["root_coordinates"], cone
                )
                vinberg_pass = bool(vinberg["literal_vinberg_passes"])
                classification = (
                    "finite_volume_noncompact"
                    if vinberg_pass and cone["ideal_vertex_count"]
                    else "compact"
                    if vinberg_pass
                    else "failed_literal_vinberg"
                )
            volume[classification] += 1
            known_matches = [
                label
                for label, graph in known.items()
                if _weighted_isomorphic(diagram_class["graph"], graph)
            ]
            record = {
                "pairing_two_count": pairing_two_count,
                "class_index_within_pairing_count_zero_based": class_index,
                "labelled_member_count": len(diagram_class["members"]),
                "support_edge_count": representative["underlying_edge_count"],
                "support_graph6": representative["underlying_graph6"],
                "promoted_pairing_two_edges_zero_based": representative[
                    "promoted_parallel_edges_zero_based"
                ],
                "radical_quotient_determinant": quotient["determinant"],
                "finite_ideal_hyperideal": [
                    cone["finite_vertex_count"],
                    cone["ideal_vertex_count"],
                    hyperideal,
                ],
                "classification": classification,
                "literal_vinberg_passes_when_no_hyperideal": vinberg_pass,
                "known_root_diagram_matches": known_matches,
                "positive_strict_interior": interior,
                "complete_ray_record_sha256": _digest(cone["rays"]),
            }
            records.append(record)
            if not hyperideal:
                finite.append(record)
    return {
        "edge_range": edge_range,
        "counts_by_pairing_two_count": counts,
        "weighted_trials": trials,
        "exact_Lorentzian_marked_presentations": marked_total,
        "weighted_isometry_classes": classes_total,
        "volume_ledger": dict(sorted(volume.items())),
        "no_hyperideal_records": finite,
        "class_records": records,
    }


def _p10_counterexample(
    hostile: dict[str, Any],
) -> dict[str, Any]:
    rows = [
        row
        for row in hostile["presentations"]
        if row["chamber_id"] == "R12-E8A1x2-C2-H10-12"
    ]
    if len(rows) != 1 or len(rows[0]["candidate_audits"]) != 1:
        raise ArithmeticError("the pinned P10 presentation changed")
    audit = rows[0]["candidate_audits"][0]
    if not (
        audit["coble_incidence_pass"]
        and audit["complete_package_regenerated_from_P"]
        and audit["b_perp"]["is_E10"]
        and audit["cusp_regeneration"]["sets_equal"]
        and audit["surface_classification_gate"]["status"]
        == "Coble_carrier_obstructed"
        and len(audit["root_facet_indices_one_based"]) == 10
    ):
        raise ArithmeticError("the P10 automaticity counterexample changed")
    return {
        "raw_chamber": "classical P10 = R12-E8A1x2-C2-H10-12",
        "unique_characteristic_norm_minus_one_boundary": True,
        "boundary_orthogonal_lattice": "even unimodular E10",
        "root_count": 10,
        "root_cone_vertices_finite_ideal": [9, 1],
        "complete_cusp_lift_count": 1,
        "complete_raw_package_regenerated": True,
        "D_CobVin_arithmetic_gates_pass": True,
        "base_field_scope": "complex numbers (characteristic zero)",
        "exact_complex_one_boundary_strong_Coble_full_nef_carrier_exists": False,
        "AG_obstruction": (
            "finite root-reflection index would force finite Aut, but Kondo's "
            "complex one-boundary finite-Aut Types I and II both have twelve "
            "effective irreducible roots, whereas P10 has ten"
        ),
        "consequence": "R_CobVin(C) is a proper subset of D_CobVin",
    }


def build(cdd: str | None = None) -> dict[str, Any]:
    cdd = cdd or shutil.which("cddexec_gmp") or "cddexec_gmp"
    pinned = _load_pinned()
    low = _tail(
        "11:13",
        expected_marked={0: 0, 1: 4, 2: 53, 3: 99},
        expected_classes={0: 0, 1: 4, 2: 47, 3: 79},
        cdd=cdd,
    )
    high = _tail(
        "16:18",
        expected_marked={0: 7, 1: 34, 2: 158, 3: 209},
        expected_classes={0: 7, 1: 32, 2: 121, 3: 142},
        cdd=cdd,
    )
    if not (
        len(low["no_hyperideal_records"]) == 1
        and low["no_hyperideal_records"][0]["known_root_diagram_matches"]
        == ["Type I"]
        and high["no_hyperideal_records"] == []
    ):
        raise ArithmeticError("tail finite-volume classification changed")

    p13 = pinned[ROOT / "data/enumerations/coble_nodot_13edge_lane.json"]
    p14 = pinned[ROOT / "data/enumerations/coble_nodot_14edge_lane.json"]
    p15 = pinned[ROOT / "data/enumerations/coble_nodot_15edge_lane.json"]
    p151 = pinned[
        ROOT / "data/enumerations/coble_nodot_15edge_one_pair2_lane.json"
    ]
    p152 = pinned[
        ROOT / "data/enumerations/coble_nodot_15edge_two_pair2_lane.json"
    ]
    p153 = pinned[
        ROOT / "data/enumerations/coble_nodot_15edge_three_pair2_lane.json"
    ]
    promotion = pinned[
        ROOT / "data/enumerations/coble_one_boundary_raw_k3_promotion.json"
    ]
    hostile = pinned[
        ROOT / "data/enumerations/polytope_first_coble_h10_hostile_audit.json"
    ]
    if not (
        p13["actual_full_D_1_split_isomorphism_class_count"] == 1
        and p14["denominator"]["weighted_diagram_isometry_class_count"] == 194
        and p15["full_D_1_split_isomorphism_class_count"] == 1
        and p151["full_D_1_split_isomorphism_class_count"] == 0
        and p152["full_D_1_split_isomorphism_class_count"] == 0
        and p153["full_D_1_split_isomorphism_class_count"] == 0
        and promotion["classification_corollary"][
            "accepted_underlying_polytope_isometry_classes"
        ]
        == 2
    ):
        raise ArithmeticError("a pinned middle-lane conclusion changed")

    combined = {
        "weighted_trials": (
            low["weighted_trials"]
            + p14["denominator"]["weighted_trial_count"]
            + p15["denominator"]["weighted_trial_count"]
            + p151["denominator"]["weighted_trial_count"]
            + p152["denominator"]["weighted_trial_count"]
            + p153["denominator"]["weighted_trial_count"]
            + high["weighted_trials"]
        ),
        "exact_Lorentzian_marked_presentations": (
            low["exact_Lorentzian_marked_presentations"]
            + p14["denominator"]["exact_Lorentzian_corank_two_labelled_count"]
            + p15["denominator"]["exact_Lorentzian_corank_two_labelled_count"]
            + p151["denominator"]["exact_Lorentzian_corank_two_marked_trial_count"]
            + p152["denominator"]["exact_Lorentzian_corank_two_marked_trial_count"]
            + p153["denominator"]["exact_Lorentzian_corank_two_marked_trial_count"]
            + high["exact_Lorentzian_marked_presentations"]
        ),
        "weighted_root_isometry_classes": (
            low["weighted_isometry_classes"]
            + p14["denominator"]["weighted_diagram_isometry_class_count"]
            + p15["denominator"]["weighted_diagram_isometry_class_count"]
            + p151["denominator"]["weighted_diagram_isometry_class_count"]
            + p152["denominator"]["weighted_diagram_isometry_class_count"]
            + p153["denominator"]["weighted_diagram_isometry_class_count"]
            + high["weighted_isometry_classes"]
        ),
        "finite_volume_root_isometry_classes": 2,
        "finite_volume_root_types": ["Type I", "Type II"],
        "complete_cusp_lifted_polytope_isometry_classes": 2,
        "exact_Coble_same_P_returns": 2,
    }
    if combined != {
        "weighted_trials": 10_456_309,
        "exact_Lorentzian_marked_presentations": 1_038,
        "weighted_root_isometry_classes": 837,
        "finite_volume_root_isometry_classes": 2,
        "finite_volume_root_types": ["Type I", "Type II"],
        "complete_cusp_lifted_polytope_isometry_classes": 2,
        "exact_Coble_same_P_returns": 2,
    }:
        raise ArithmeticError(f"combined denominator changed: {combined}")

    p10 = _p10_counterexample(hostile)
    result = {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "surface_free_denominator": {
            "name": "D_1_split^nodot,Delta<=3,kappa<=3",
            "definition": (
                "twelve primitive norm-minus-two E10 roots; connected support "
                "of maximum degree at most three; pairings in {0,1,2} with at "
                "most three pairing-two pairs; exact root Gram inertia (1,9,2); "
                "finite-volume root cone; characteristic boundary and complete "
                "primitive cusp lift as in D_1_split"
            ),
            "surface_or_carrier_success_in_membership": False,
            "support_edge_range_is_exhaustive": (
                "connected twelve-vertex maximum-degree-three supports have "
                "between 11 and 18 edges"
            ),
        },
        "new_tail_censuses": {"low_11_13": low, "high_16_18": high},
        "combined_exact_census": combined,
        "same_P_carrier_theorem": {
            "statement": (
                "Every finite-volume member of the declared raw denominator is "
                "integrally isometric to Type I or Type II, and the exhaustive "
                "split-K3/Coble compiler constructs a complex one-boundary strong "
                "Coble surface whose full marked nef chamber is the same P."
            ),
            "positive_survivors": [
                {"root": "Type I H9/12", "full": "H10/22"},
                {"root": "Type II H9/12", "full": "H10/30"},
            ],
            "all_carriers_exhausted_by_existing_compiler": True,
            "hidden_walls_exhausted": True,
        },
        "strict_L4_audit": {
            "class_wide_carrier_nonemptiness": "PASS on the declared denominator",
            "exact_same_P_return": "PASS",
            "independent_AG_input": (
                "Kondo classifies complex finite-Aut one-boundary Coble surfaces "
                "as Types I and II and supplies their existence"
            ),
            "nonthin_new_AG_pushback": "FAIL",
            "reason": (
                "the exact P-side cone census has already isolated the same two "
                "root/full-wall classes before Kondo is invoked; AG supplies "
                "existence, carrier interpretation and moduli lineage, but no "
                "new dimension bound or additional polytope exclusion"
            ),
            "rating": "L3_PLUS_CLASSWIDE_SCOPED_ESSENTIAL_SURJECTIVITY_NOT_STRICT_L4",
        },
        "maximal_automaticity_counterexample": p10,
        "sharp_remaining_L4_gate": (
            "Find a natural raw condition broader than the answer-isolating "
            "max-degree-three census, excluding P10 for a structural reason, "
            "which forces incidence-stack nonemptiness while leaving multiple "
            "P-side possibilities for an independent AG theorem to eliminate."
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
