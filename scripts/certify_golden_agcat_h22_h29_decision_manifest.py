#!/usr/bin/env python3
"""Freeze the exact Q(sqrt(5)), H22--H29 AGCat decision ledger."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coxeter_polytopes.classification.golden_agcat_manifest import (  # noqa: E402
    FRAME_BOUND,
    carrier_logic_records,
    cycle_box_record,
    golden_closed_cauchy_edge_box,
    golden_coxeter_frame_pair_alphabet,
    golden_finite_angle_triangle_packets,
    golden_finite_pair_values,
    golden_frame_edge_box,
    rank_manifest_rows,
)
from coxeter_polytopes.exact.golden import Golden  # noqa: E402


OUTPUT = ROOT / "data/enumerations/golden_agcat_h22_h29_decision_manifest.json"


def _element(value) -> dict[str, int]:
    if not value.is_integral:
        raise AssertionError("manifest element left Z[phi]")
    return {"a": int(value.a), "b": int(value.b)}


def _digest_payload(value: object) -> str:
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("utf-8")
    return sha256(encoded).hexdigest()


def _artifact(relative: str) -> dict[str, object]:
    path = ROOT / relative
    data = path.read_bytes()
    return {
        "path": relative,
        "sha256": sha256(data).hexdigest(),
        "byte_count": len(data),
    }


def build_manifest() -> dict[str, object]:
    closed_edge_box = golden_closed_cauchy_edge_box()
    edge_box = golden_frame_edge_box()
    pair_alphabet = golden_coxeter_frame_pair_alphabet()
    finite_pairs = golden_finite_pair_values()
    triangle_packets = golden_finite_angle_triangle_packets()
    admissible_triangles = tuple(row for row in triangle_packets if row.admissible_lorentzian)
    edge_payload = [_element(value) for value in edge_box]
    alphabet_payload = [_element(value) for value in pair_alphabet]
    alphabet_set = set(pair_alphabet)
    invalid_subparallel = tuple(
        value
        for value in edge_box
        if (value - 4).sign(0) < 0 and value not in alphabet_set
    )
    dotted = tuple(value for value in pair_alphabet if (value - 4).sign(0) > 0)
    cauchy_equality = tuple(value for value in closed_edge_box if (value - 4).sign(1) == 0)

    if len(closed_edge_box) != 6_878 or len(edge_box) != 6_877:
        raise AssertionError("golden edge-box count changed")
    if cauchy_equality != (Golden(4),) or Golden(4) in edge_box:
        raise AssertionError("independent-frame Cauchy equality was not removed")
    if tuple(map(str, invalid_subparallel)) != ("2 - 1*phi", "3 - 1*phi"):
        raise AssertionError("unexpected invalid subparallel values")
    if len(pair_alphabet) != 6_875 or len(dotted) != 6_869:
        raise AssertionError("golden Coxeter pair alphabet census changed")
    if len(triangle_packets) != 10 or len(admissible_triangles) != 5:
        raise AssertionError("golden finite-angle triangle census changed")
    if tuple(row.labels for row in admissible_triangles) != (
        (3, 3, 5),
        (3, 5, 5),
        (4, 4, 5),
        (5, 5, 5),
        (5, 10, 10),
    ):
        raise AssertionError("golden admissible triangle list changed")

    ranks = rank_manifest_rows()
    if len(ranks) != 8 or sum(row["required_cycle_box_count"] for row in ranks) != 204:
        raise AssertionError("rank/cycle application count changed")

    source_artifacts = {
        "universal_arrangement_AGCat": _artifact(
            "results/exact_agcat_functor_and_ag_finiteness_pushback_synthesis.md"
        ),
        "global_maximal_finiteness": _artifact(
            "results/agcat_global_maximal_arithmetic_coxeter_finiteness.md"
        ),
        "computational_architecture": _artifact(
            "results/computational_agcat_coxeter_classification_architecture.md"
        ),
        "coherent_full_H2_bound": _artifact(
            "results/coherent_full_h2_real_multiplication_bound.md"
        ),
        "golden_rank20_Hodge_gate": _artifact(
            "results/golden_rank20_full_h2_extension_gate.md"
        ),
        "maximal_overgroup_existence_argument": _artifact(
            "results/higher_dimensional_coxeter_orbifold_ag_quotient_no_go.md"
        ),
        "finite_index_field_invariance_argument": _artifact(
            "results/unconditional_maximal_arithmetic_reflection_chamber_finiteness.md"
        ),
        "golden_Dn_prefix": _artifact(
            "data/enumerations/golden_dn_h22_h29_through_phi2.json"
        ),
        "golden_A1_Golay_cutter": _artifact(
            "data/enumerations/golden_a1_uniform_ray_cutter.json"
        ),
        "golden_Leech_norm128": _artifact(
            "data/enumerations/golden_leech_low_norm_manifest.json"
        ),
        "golden_Leech_cleanroom": _artifact(
            "data/enumerations/golden_leech_low_norm_manifest_cleanroom.json"
        ),
        "golden_odd_Leech_section_H23": _artifact(
            "data/enumerations/golden_odd_leech_section_h23_first_packets.json"
        ),
        "golden_full_triangle_packets": _artifact(
            "data/enumerations/golden_agcat_triangle_packets.json"
        ),
        "golden_four_node_at_most_one_dotted_packets": _artifact(
            "data/enumerations/golden_agcat_four_node_at_most_one_dotted.json"
        ),
        "golden_four_node_at_most_one_dotted_independent_hostile_audit": _artifact(
            "data/enumerations/golden_agcat_four_node_at_most_one_dotted_hostile_audit.json"
        ),
        "golden_two_dotted_four_node_trees": _artifact(
            "data/enumerations/golden_agcat_two_dotted_four_node_trees.json"
        ),
        "golden_two_dotted_four_node_trees_independent_hostile_audit": _artifact(
            "data/enumerations/golden_agcat_two_dotted_four_node_trees_independent_audit.json"
        ),
    }

    return {
        "schema": "golden-agcat-h22-h29-decision-manifest-0.1.0",
        "field": {
            "name": "Q(sqrt(5))",
            "generator": "phi",
            "minimal_polynomial": "t^2-t-1",
            "ring_of_integers": "Z[phi]",
            "discriminant": 5,
            "real_embeddings": ["phi=(1+sqrt(5))/2", "phi'=(1-sqrt(5))/2"],
            "class_number": 1,
            "narrow_class_number": 1,
        },
        "decision_target": {
            "compact": True,
            "arithmetic": True,
            "hyperbolic_dimensions": list(range(22, 30)),
            "native_ranks": list(range(23, 31)),
            "primary_enumeration_scope": {
                "candidate_queue": (
                    "all full native hosts R(L) reconstructed from retained golden frame/order/lattice packages"
                ),
                "every_candidate_full_host_must_be_completely_classified": True,
                "direct_terminal_empty_set": (
                    "compact reflective full-host survivors among all classified candidate R(L)"
                ),
                "direct_terminal_empty_set_must_be_empty": True,
                "all_full_native_hosts_themselves_must_be_absent": False,
                "direct_route_reason": (
                    "for every hypothetical compact W, the invariant-lattice construction "
                    "already gives W subset R(L) at finite index; therefore R(L) is a compact "
                    "reflective full-host survivor, without first maximalizing W"
                ),
                "equivalent_maximal_completion_cross_check": (
                    "one may also maximalize W; for Wmax the same construction gives "
                    "Wmax subset R(L), and reflection maximality forces equality"
                ),
            },
            "reduction_to_all_compact_chambers": {
                "status": "proved_by_the_declared_standard_lattice_arguments",
                "fixed_dimension": True,
                "fixed_ambient": (
                    "the single Lie group Isom(H^n) containing the fixed input W; no common "
                    "integral lattice or preselected commensurator is assumed"
                ),
                "maximal_overgroup_exists": (
                    "for every discrete reflection R containing the fixed finite-covolume W, "
                    "R is a lattice and [R:W]=vol(W\\H^n)/vol(R\\H^n); the fixed-dimensional "
                    "Kazhdan-Margulis lower covolume makes this set of indices bounded. The "
                    "set is nonempty (R=W), so it has a maximum; if Rmax had any strict "
                    "overgroup R', then [R':W]=[R':Rmax][Rmax:W] would be larger. No nesting "
                    "among the other candidate overgroups is required"
                ),
                "finite_index": True,
                "compactness_preserved": True,
                "arithmeticity_preserved": True,
                "minimal_Vinberg_field_preserved": (
                    "finite-index commensurability invariance via Vinberg's least adjoint "
                    "ring/field and cyclic Cartan-product criterion"
                ),
                "dimension_preserved": True,
                "conclusion": (
                    "empty maximal-host list at this field and rank implies no nonmaximal "
                    "compact arithmetic chamber at that field and rank"
                ),
            },
        },
        "carrier_logic": list(carrier_logic_records()),
        "universal_length_two_box": {
            "frame_bound": FRAME_BOUND,
            "weak_closed_Cauchy_superset_definition": (
                "b in Z[phi], 0<=sigma0(b)<62^2 and 0<=sigma1(b)<=4"
            ),
            "weak_closed_Cauchy_superset_count": len(closed_edge_box),
            "Cauchy_equality_elements": [_element(value) for value in cauchy_equality],
            "Cauchy_equality_excluded_reason": (
                "positive-definite equality makes the two nonselected embedded roots "
                "proportional, contradicting independence of a frame pair"
            ),
            "exact_independent_frame_box_definition": (
                "b in Z[phi], 0<=sigma0(b)<62^2 and 0<=sigma1(b)<4"
            ),
            "complete_box_count": len(edge_box),
            "complete_box_sha256": _digest_payload(edge_payload),
            "complete_box_elements": edge_payload,
            "finite_Coxeter_angle_products": [
                {"label": row.label, "product": _element(row.product)}
                for row in finite_pairs
            ],
            "parallel_product_present": False,
            "parallel_product_excluded": _element(Golden(4)),
            "compactness_also_excludes_parallel_facets": True,
            "invalid_subparallel_non_Coxeter_values": [
                _element(value) for value in invalid_subparallel
            ],
            "dotted_product_count": len(dotted),
            "exact_Coxeter_pair_alphabet_count": len(pair_alphabet),
            "exact_Coxeter_pair_alphabet_sha256": _digest_payload(alphabet_payload),
            "exact_Coxeter_pair_alphabet_elements": alphabet_payload,
            "status": "enumerated_exactly",
        },
        "finite_angle_triangle_cycle_tranche": {
            "scope": (
                "complete three-edge frames whose three pair relations are finite and "
                "nonorthogonal; modulo all triangle vertex permutations"
            ),
            "raw_unordered_label_multisets": 35,
            "integral_golden_cycle_packet_count": len(triangle_packets),
            "admissible_lorentzian_packet_count": len(admissible_triangles),
            "packets": [
                {
                    "labels": list(row.labels),
                    "edge_products": [_element(value) for value in row.edge_products],
                    "directed_three_cycle": _element(row.directed_three_cycle),
                    "gram_determinant": _element(row.gram_determinant),
                    "selected_determinant_sign": row.selected_determinant_sign,
                    "conjugate_determinant_sign": row.conjugate_determinant_sign,
                    "admissible_lorentzian": row.admissible_lorentzian,
                }
                for row in triangle_packets
            ],
            "admissible_label_triples": [list(row.labels) for row in admissible_triangles],
            "status": "enumerated_exactly",
            "universal_scope_warning": (
                "triangles containing a dotted edge and non-complete graph triples are not "
                "covered by this finite-angle tranche"
            ),
        },
        "full_triangle_cycle_tranche": {
            "scope": (
                "every complete nonorthogonal three-node frame over the full "
                "6875-element edge alphabet, including all dotted values, modulo S3"
            ),
            "raw_unordered_nonzero_edge_multisets": 54_158_527_500,
            "integral_golden_square_product_packet_count": 4_768_141,
            "conjugate_positive_packet_count": 2_633_885,
            "retained_selected_lorentzian_packet_count": 2_633_885,
            "raw_selected_degenerate_packet_count": 1,
            "retained_selected_degenerate_packet_count": 0,
            "certificate": source_artifacts["golden_full_triangle_packets"],
            "status": "enumerated_exactly",
            "high_rank_scope_warning": (
                "complete nonorthogonal triangles are not forced in a connected "
                "rank-23--30 graph; orthogonal/tree triples, overlapping compatibility, "
                "longer cycles, gauges, orders, conductors, and hosts remain"
            ),
        },
        "four_node_at_most_one_dotted_tranche": {
            "scope": (
                "every connected four-node packet with only finite/orthogonal edges, "
                "and every connected packet with exactly one dotted edge ranging over "
                "all 6869 dotted values; simultaneous compatibility on every present "
                "triangle and Hamilton four-cycle"
            ),
            "finite_only_connected_orbits_before_cycle_gate": 2_380,
            "one_dotted_connected_orbits_before_cycle_gate": 14_287_520,
            "finite_only_cycle_compatible_orbits": 327,
            "one_dotted_cycle_compatible_orbits": 532_734,
            "finite_only_two_place_retained_orbits": 119,
            "one_dotted_two_place_retained_orbits": 180_402,
            "overlapping_cycle_topologies_retained": [
                "diamond_K4_minus_edge",
                "complete_K4",
            ],
            "global_one_dotted_label3_path_family": {
                "native_ranks": list(range(23, 31)),
                "admissible_dotted_value_counts": [
                    3_590, 3_584, 3_577, 3_571, 3_566, 3_561, 3_556, 3_552
                ],
                "uniform_sample_dotted_product": _element(Golden(2, 2)),
                "status": "exact_full_rank_frame_packets",
                "ambient_arithmetic_lattices_are_cocompact": True,
                "cocompactness_reason": (
                    "conjugate positive definiteness forbids a nonzero isotropic "
                    "Q(sqrt(5))-vector"
                ),
                "compactness_warning": (
                    "the associated n+1-wall simplicial cones are not finite-volume; "
                    "full Vinberg completion is required before any reflective or "
                    "compact example can be claimed"
                ),
            },
            "certificate": source_artifacts[
                "golden_four_node_at_most_one_dotted_packets"
            ],
            "independent_hostile_audit": source_artifacts[
                "golden_four_node_at_most_one_dotted_independent_hostile_audit"
            ],
            "independent_hostile_audit_status": "PASS",
            "status": "enumerated_exactly",
            "high_rank_scope_warning": (
                "four-node packets with at least two dotted edges and all compatible "
                "connected rank-23--30 graphs remain unenumerated; this local tranche "
                "does not imply high-rank nonexistence"
            ),
        },
        "four_node_exactly_two_dotted_tree_tranche": {
            "scope": (
                "every connected four-node tree with exactly two dotted edges over "
                "the complete 6869-value dotted alphabet and one of the five "
                "nonorthogonal finite golden Coxeter products; modulo graph automorphisms"
            ),
            "canonical_tree_packets_before_gram_gate": 471_865_955,
            "two_place_retained_packets": 97_172_716,
            "retained_by_topology": {
                "star": 25_023_281,
                "disjoint_dotted_path": 81,
                "adjacent_dotted_path": 72_149_354,
            },
            "integral_tree_gauge": (
                "orient each product c as (-1,-c); tree propagation supplies a "
                "totally positive integral diagonal symmetrizer"
            ),
            "production_certificate": source_artifacts[
                "golden_two_dotted_four_node_trees"
            ],
            "independent_hostile_audit": source_artifacts[
                "golden_two_dotted_four_node_trees_independent_hostile_audit"
            ],
            "independent_hostile_audit_status": "PASS_AFTER_CORRECTING_DISJOINT_PATH_LOCUS",
            "correction_history": (
                "the first production draft incorrectly rejected every disjoint-dotted "
                "path; the hostile audit found 81 survivors. The current production "
                "artifact, theorem count, and exact Fenwick dominance join include them"
            ),
            "status": "enumerated_exactly",
            "high_rank_scope_warning": (
                "cyclic two-dotted four-node packets, packets with at least three "
                "dotted edges, and their compatible joins into ranks 23--30 remain open"
            ),
        },
        "cycle_boxes": [cycle_box_record(length) for length in range(2, 31)],
        "rank_rows": list(ranks),
        "existing_normalized_diagonal_and_root_block_screens": [
            {
                "family": "<-2phi> plus D_n",
                "dimensions": list(range(22, 30)),
                "coverage": "norm-two controls 1, phi, phi^2 and exact current cones",
                "outcome": "all prefixes have hyperideal rays",
                "global_host_decision": False,
                "remaining_debt": "complete mirror-norm audit and every later Vinberg level",
            },
            {
                "family": "<-2phi> plus N(A1^24)",
                "dimensions": [24],
                "coverage": (
                    "complete norm class; exact controls through phi^3; accepted 6072-root "
                    "cutter family against all earlier complete shells"
                ),
                "outcome": "replacement hyperideal extreme ray survives the audited prefix",
                "global_host_decision": False,
                "remaining_debt": "full phi^3 accepted packet and all later Vinberg levels",
            },
            {
                "family": "<-a> plus Leech",
                "dimensions": [24],
                "coverage": (
                    "all 47 squarefree coefficient-ideal directions of norm <=128, complete "
                    "local norm classes, first necessary packets"
                ),
                "outcome": "strict first-prefix witnesses for all 47 directions",
                "global_host_decision": False,
                "remaining_debt": (
                    "actual Leech realization/acceptance, later levels, nonunit square "
                    "multiples, and coefficient directions beyond norm 128"
                ),
            },
            {
                "family": "<-a> plus S23 odd-Leech section",
                "dimensions": [23],
                "coverage": (
                    "complete local norms, height-zero exclusion, and first positive "
                    "packets for a=phi and a=2phi"
                ),
                "outcome": (
                    "553104 and 53682 walls respectively; strict hyperideal points "
                    "survive both exact prefixes"
                ),
                "global_host_decision": False,
                "remaining_debt": "every later Vinberg level and other coefficient ideals",
            },
        ],
        "finite_enumeration_debt": {
            "rank_field_pairs": 8,
            "cycle_box_definitions": 29,
            "rank_cycle_box_applications": 204,
            "connected_graph_manifests_completed": 0,
            "complete_golden_cycle_packet_manifests_completed": 0,
            "local_finite_angle_triangle_cycle_tranches_completed": 1,
            "local_full_triangle_cycle_tranches_completed": 1,
            "local_four_node_at_most_one_dotted_tranches_completed": 1,
            "local_four_node_exactly_two_dotted_tree_tranches_completed": 1,
            "golden_ideal_gauge_backend_completed": False,
            "reflection_order_conductor_manifests_completed": 0,
            "complete_full_host_Vinberg_transcripts_completed": 0,
            "compact_reflective_full_host_survivor_empty_certificates_completed": 0,
            "terminal_maximal_completion_empty_certificates_completed": 0,
            "canonical_host_comparison_and_empty_certificates_completed": 0,
        },
        "verdict": {
            "existential_finiteness_of_reflection_maximal_slice_proved": True,
            "effective_finite_queue_architecture_proved": True,
            "universal_pair_box_materialized": True,
            "enumerated_emptiness": False,
            "nonexistence_in_dimensions_22_29_over_Qsqrt5": "OPEN",
            "compact_or_finite_volume_survivor_found": False,
            "prohibited_inference": (
                "finite decision space does not imply that the decision space is empty"
            ),
        },
        "source_artifacts": source_artifacts,
    }


def verify_manifest(payload: dict[str, object]) -> None:
    rebuilt = build_manifest()
    if payload != rebuilt:
        raise AssertionError("manifest differs from an exact rebuild")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        verify_manifest(json.loads(args.output.read_text(encoding="utf-8")))
        print(json.dumps({"manifest": str(args.output), "status": "PASS"}, sort_keys=True))
        return
    payload = build_manifest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "edge_box": payload["universal_length_two_box"]["complete_box_count"],
                "pair_alphabet": payload["universal_length_two_box"][
                    "exact_Coxeter_pair_alphabet_count"
                ],
                "rank_rows": len(payload["rank_rows"]),
                "nonexistence": payload["verdict"][
                    "nonexistence_in_dimensions_22_29_over_Qsqrt5"
                ],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
