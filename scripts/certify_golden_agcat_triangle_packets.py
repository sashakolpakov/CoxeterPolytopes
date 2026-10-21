#!/usr/bin/env python3
"""Freeze the full exact golden three-node AGCat packet census."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.classification.golden_triangle_packets import (
    golden_triangle_packet_census,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data/enumerations/golden_agcat_triangle_packets.json"


def build_payload() -> dict[str, object]:
    census = golden_triangle_packet_census()
    if (
        census.alphabet_count != 6_875
        or census.zero_edge_value_count != 1
        or census.nonzero_edge_value_count != 6_874
        or census.squareclass_bucket_count != 3_661
        or census.maximum_squareclass_bucket_size != 109
        or census.compatible_squareclass_bucket_triples != 27_088
        or census.raw_unordered_nonzero_edge_multisets != 54_158_527_500
        or census.square_product_edge_multisets != 4_768_141
        or census.conjugate_positive_packet_count != 2_633_885
        or census.selected_positive_packet_count != 0
        or census.selected_lorentzian_packet_count != 2_633_885
        or census.raw_selected_zero_determinant_count != 1
        or census.retained_selected_zero_determinant_count != 0
        or census.conjugate_nonpositive_packet_count != 2_134_256
    ):
        raise ArithmeticError("the frozen golden triangle census changed")
    return {
        "schema": "golden-agcat-full-triangle-packets-0.1.0",
        "scope": {
            "field": "Q(sqrt(5))",
            "frame_bound": 62,
            "edge_alphabet": (
                "the complete 6875-element Coxeter-filtered independent-frame "
                "length-two alphabet from the golden H22--H29 AGCat manifest"
            ),
            "graph": "complete nonorthogonality graph on three nodes",
            "canonicalization": (
                "unordered edge-product multiset, equivalently modulo all S3 "
                "vertex permutations"
            ),
            "strict_local_stage_warning": (
                "This is the complete local three-node tranche only. It does not "
                "enumerate compatible rank-23--30 connected graphs, longer cycles, "
                "ideal gauges, orders, conductor lattices, or Vinberg hosts. A "
                "connected high-rank frame need not contain a complete nonorthogonal "
                "triangle; three-node induced subgraphs can be trees or contain zero "
                "edges, and no universal triangular three-face is asserted."
            ),
        },
        "algorithm": {
            "zero_edge_removed": True,
            "reason": "a complete three-node graph has three nonorthogonal edges",
            "principal_ideal_squareclass_factorization": (
                "exact rational-norm factorization, inert/ramified valuations, "
                "split-prime Hensel valuations, and the conjugate unit-sign class"
            ),
            "squareclass_join": (
                "group 6874 nonzero edge values by element squareclass; join bucket "
                "triples by XOR instead of materializing the cubic edge product"
            ),
            "golden_square_root_solver": (
                "for A+B*phi=(m+n*phi)^2, t=n^2 solves "
                "5*t^2-(4*A+2*B)*t+B^2=0; every root is multiplied back exactly"
            ),
            "directed_cycle_orientation": (
                "choose the identity-positive square root and negate it for the "
                "acute directed three-cycle"
            ),
            "gram_determinant": "8-2*(c12+c23+c31)+2*b123",
            "stream_materialized": False,
            "stream_digests_freeze_every_canonical_record": True,
        },
        "counts": {
            "edge_alphabet": census.alphabet_count,
            "zero_edge_values": census.zero_edge_value_count,
            "nonzero_edge_values": census.nonzero_edge_value_count,
            "element_squareclass_buckets": census.squareclass_bucket_count,
            "maximum_bucket_size": census.maximum_squareclass_bucket_size,
            "compatible_bucket_triples": census.compatible_squareclass_bucket_triples,
            "raw_unordered_nonzero_edge_multisets": census.raw_unordered_nonzero_edge_multisets,
            "integral_golden_square_product_packets": census.square_product_edge_multisets,
            "conjugate_positive_definite_packets": census.conjugate_positive_packet_count,
            "conjugate_nonpositive_rejections": census.conjugate_nonpositive_packet_count,
            "retained_selected_positive_definite": census.selected_positive_packet_count,
            "retained_selected_lorentzian": census.selected_lorentzian_packet_count,
            "raw_selected_degenerate": census.raw_selected_zero_determinant_count,
            "retained_selected_degenerate": census.retained_selected_zero_determinant_count,
        },
        "determinant_logic": {
            "conjugate_requirement": (
                "strictly positive determinant; the conjugate form is positive "
                "definite and the three frame roots remain independent"
            ),
            "selected_universal_possibilities": [
                "positive definite", "degenerate", "Lorentzian"
            ],
            "selected_degenerate_compact_rule": (
                "count universally, but reject for a compact reflection target: an "
                "affine/parabolic Coxeter subdiagram produces a parabolic element, "
                "whereas a cocompact lattice contains no parabolic element"
            ),
            "observed_retained_result": (
                "all 2633885 packets passing conjugate positivity are already "
                "selected-Lorentzian; no selected-degenerate packet survives"
            ),
        },
        "strata": [
            {
                "edge_kinds": list(row.edge_kinds),
                "identity_positive_root_conjugate_sign": row.root_conjugate_sign,
                "selected_determinant_sign": row.selected_determinant_sign,
                "conjugate_determinant_sign": row.conjugate_determinant_sign,
                "count": row.count,
            }
            for row in census.strata
        ],
        "stream_digests": {
            "all_integral_square_product_packets_sha256": census.canonical_stream_sha256,
            "conjugate_positive_packets_sha256": census.retained_stream_sha256,
            "record_format": (
                "sorted edges c1;c2;c3;negative directed cycle;determinant, each "
                "golden integer encoded as decimal a,b and terminated by newline"
            ),
        },
        "outcome": {
            "complete_local_triangle_tranche": True,
            "dotted_edges_included": True,
            "new_high_rank_nonexistence_proof": False,
            "selected_positive_zero_is_not_a_no_go": (
                "complete nonorthogonal triangles are not forced in a connected "
                "high-rank frame, while a selected-Lorentzian three-node restriction "
                "is locally allowed inside a global index-one form"
            ),
            "Qsqrt5_H22_H29_status": "OPEN",
        },
    }


def encoded_payload() -> bytes:
    return (json.dumps(build_payload(), indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    encoded = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()
    if args.verify:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != encoded:
            raise SystemExit("FAIL: frozen golden triangle JSON differs from recomputation")
        status = "PASS"
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(encoded)
        status = "WROTE"
    print(
        json.dumps(
            {
                "status": status,
                "output": str(OUTPUT),
                "sha256": sha256(encoded).hexdigest(),
                "square_product_packets": payload["counts"][
                    "integral_golden_square_product_packets"
                ],
                "retained_packets": payload["counts"][
                    "conjugate_positive_definite_packets"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
