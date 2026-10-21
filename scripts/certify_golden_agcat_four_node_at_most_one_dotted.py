#!/usr/bin/env python3
"""Freeze the exact golden four-node at-most-one-dotted packet census."""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path

from coxeter_polytopes.classification.golden_four_node_packets import (
    golden_four_node_packet_census,
    golden_one_dotted_path_rank_packets,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data/enumerations/golden_agcat_four_node_at_most_one_dotted.json"
)


def _topology_rows(census) -> list[dict[str, object]]:
    counts: Counter[tuple[str, str]] = Counter()
    retained: Counter[tuple[str, str]] = Counter()
    for row in census.strata:
        key = row.tranche, row.topology
        counts[key] += row.count
        conjugate_positive = (
            row.conjugate_leading_triple_sign > 0
            and row.conjugate_full_determinant_sign > 0
        )
        if row.tranche == "exactly_one_dotted":
            selected_allowed = (
                row.selected_leading_triple_sign < 0
                and row.selected_full_determinant_sign < 0
            )
        else:
            selected_allowed = row.selected_full_determinant_sign < 0 or (
                row.selected_leading_triple_sign > 0
                and row.selected_full_determinant_sign > 0
            )
        if conjugate_positive and selected_allowed:
            retained[key] += row.count
    return [
        {
            "tranche": tranche,
            "topology": topology,
            "cycle_compatible": count,
            "two_place_retained": retained[(tranche, topology)],
        }
        for (tranche, topology), count in sorted(counts.items())
    ]


def build_payload() -> dict[str, object]:
    census = golden_four_node_packet_census()
    path_rows = golden_one_dotted_path_rank_packets()
    frozen = {
        "finite_value_count": 6,
        "dotted_value_count": 6_869,
        "finite_only_connected_orbit_count": 2_380,
        "one_dotted_complement_orbit_count": 2_080,
        "one_dotted_raw_orbit_count": 14_287_520,
        "finite_only_cycle_compatible_count": 327,
        "one_dotted_cycle_compatible_count": 532_734,
        "finite_only_conjugate_positive_count": 119,
        "one_dotted_conjugate_positive_count": 180_402,
        "finite_only_retained_count": 119,
        "finite_only_retained_selected_positive_count": 5,
        "finite_only_retained_selected_lorentzian_count": 114,
        "one_dotted_retained_count": 180_402,
        "one_dotted_retained_selected_lorentzian_count": 180_402,
        "one_dotted_selected_nonlorentzian_count": 0,
        "finite_only_selected_leading_triple_degenerate_count": 21,
        "finite_only_selected_full_degenerate_count": 3,
        "one_dotted_selected_leading_triple_degenerate_count": 0,
        "one_dotted_selected_full_degenerate_count": 0,
        "finite_only_conjugate_leading_triple_degenerate_count": 21,
        "finite_only_conjugate_full_degenerate_count": 3,
        "one_dotted_conjugate_leading_triple_degenerate_count": 0,
        "one_dotted_conjugate_full_degenerate_count": 0,
        "finite_only_compatible_stream_sha256": (
            "687762835339ac1f56fed84bcc41975e51da471abba4d97a7ab4cafc64590e7a"
        ),
        "finite_only_retained_stream_sha256": (
            "a8a04a852ce054ebfa1668156200998263939d7417d7c437662130b2320c50e8"
        ),
        "one_dotted_compatible_stream_sha256": (
            "0c988e10b93d804457dcfb9f238a832c546fd749b0a66751e67c875462153cc5"
        ),
        "one_dotted_retained_stream_sha256": (
            "7de03d51631e744ba34c248a23a933da855e841256f7f27a33c1cf09b8a0b0a1"
        ),
    }
    for field, expected in frozen.items():
        if getattr(census, field) != expected:
            raise ArithmeticError(
                f"golden four-node frozen invariant {field} changed: "
                f"{getattr(census, field)!r} != {expected!r}"
            )
    expected_path_counts = (3_590, 3_584, 3_577, 3_571, 3_566, 3_561, 3_556, 3_552)
    expected_path_hashes = (
        "41abc5d3199c0d6c9c9a8e3d43062e71b9eef438c902105f7d013d540877d9fd",
        "213ab61fe4f5553a11bbf57a28a95c8ac77625ca90ed13051b79cc72212bb1b5",
        "b9565103be25f316cf66ffe861c291de0dc3636b537d8f0d48c70bd8f8db440a",
        "ec319cc01937d0f71e862c094c26b4c9aa89e40117dacce8e1600d76ca82ef5b",
        "82fa892077eca7a86e7306d9974d19a46e74234996ee32f48ba80e107d289f05",
        "50c7a2c02c6ac629e916a8600ccf56e3a9e7bf9b2dd8e424baf64dafdbb8e37f",
        "b67576b086ff11c0124207ead54a3382e54c12fdf0ea7aa42faa8dbf09cb8e6c",
        "aa9ddb53c6aad025f3adc785cc7fdc0c8f784711a974b5e146d749011fbb7f96",
    )
    if tuple(row.admissible_dotted_value_count for row in path_rows) != expected_path_counts:
        raise ArithmeticError("golden full-rank path counts changed")
    if tuple(row.admissible_dotted_stream_sha256 for row in path_rows) != expected_path_hashes:
        raise ArithmeticError("golden full-rank path hashes changed")

    return {
        "schema": "golden-agcat-four-node-at-most-one-dotted-0.1.0",
        "scope": {
            "field": "Q(sqrt(5))",
            "native_target_ranks": list(range(23, 31)),
            "frame_bound": 62,
            "four_node_graph": "connected",
            "tranches": [
                {
                    "name": "finite_only",
                    "pair_slots": (
                        "all six slots independently range over the complete six-value "
                        "finite/orthogonal golden Coxeter alphabet"
                    ),
                },
                {
                    "name": "exactly_one_dotted",
                    "pair_slots": (
                        "the unique dotted slot ranges over every one of the 6869 "
                        "dotted values; the other five slots independently range over "
                        "all six finite/orthogonal values"
                    ),
                },
            ],
            "not_covered": (
                "four-node packets with two or more dotted edges; compatible connected "
                "rank-23--30 graphs; ideal gauges; reflection orders; conductor lattices; "
                "complete Vinberg hosts"
            ),
            "high_rank_warning": (
                "This exact local four-node tranche is not a proof that every high-rank "
                "golden frame lies in it, and it is not a high-dimensional no-go."
            ),
        },
        "canonicalization": {
            "finite_only": "quotient all connected assignments by S4",
            "exactly_one_dotted": (
                "send the unique dotted edge to 01, then quotient its finite complement "
                "by the order-four setwise stabilizer S2 x S2"
            ),
            "finite_only_connected_orbits_before_cycles": (
                census.finite_only_connected_orbit_count
            ),
            "one_dotted_finite_complement_orbits_before_cycles": (
                census.one_dotted_complement_orbit_count
            ),
            "one_dotted_orbits_after_expanding_all_dotted_values_before_cycles": (
                census.one_dotted_raw_orbit_count
            ),
        },
        "cycle_compatibility": {
            "cycles_checked_simultaneously": (
                "all four triangles and all three Hamilton four-cycles that are present"
            ),
            "condition": (
                "the XOR of golden element squareclasses on every graph cycle is zero"
            ),
            "sufficiency": (
                "for a connected graph, propagate a totally positive symmetrizer along "
                "a spanning tree; the zero cycle classes make every chord quotient a "
                "square in Q(sqrt(5)). This is sufficiency for the field-level "
                "symmetrizable cycle packet, not for an integral reflection gauge. "
                "Integrality of cycle roots follows because Z[phi] is integrally closed"
            ),
            "one_dotted_reduction": (
                "cycles avoiding the unique dotted edge must already have trivial finite "
                "squareclass; all cycles containing it must prescribe one common dotted "
                "squareclass. If it is a bridge, all 6869 dotted values survive the cycle gate"
            ),
            "finite_only_compatible_orbits": census.finite_only_cycle_compatible_count,
            "one_dotted_compatible_orbits": census.one_dotted_cycle_compatible_count,
        },
        "two_place_gram_gate": {
            "normalization": "diagonal entries 2 and pair products c_ij",
            "four_node_determinant": (
                "16 - 4*sum(c_e) + sum_opposite(c_e*c_f) "
                "- 4*sum_triangle(sqrt(product c_e)) "
                "- 2*sum_four_cycle(sqrt(product c_e))"
            ),
            "root_orientation": (
                "identity-positive roots; odd acute triangle cycles are negative and "
                "even acute four-cycles are positive"
            ),
            "conjugate_requirement": (
                "strict positive definiteness, checked by exact leading 3x3 and 4x4 "
                "Sylvester determinants; every 2x2 determinant is already positive"
            ),
            "selected_finite_only_requirement": (
                "positive definite or Lorentzian; degeneracy and two-negative inertia rejected"
            ),
            "selected_one_dotted_requirement": (
                "Lorentzian; the dotted 2x2 minor already supplies one negative direction"
            ),
            "finite_only_conjugate_positive": census.finite_only_conjugate_positive_count,
            "one_dotted_conjugate_positive": census.one_dotted_conjugate_positive_count,
            "finite_only_retained": census.finite_only_retained_count,
            "finite_only_retained_selected_positive": (
                census.finite_only_retained_selected_positive_count
            ),
            "finite_only_retained_selected_lorentzian": (
                census.finite_only_retained_selected_lorentzian_count
            ),
            "one_dotted_retained": census.one_dotted_retained_count,
            "one_dotted_retained_selected_lorentzian": (
                census.one_dotted_retained_selected_lorentzian_count
            ),
            "observed_one_dotted_selected_nonlorentzian": (
                census.one_dotted_selected_nonlorentzian_count
            ),
            "zero_determinant_counts": {
                "finite_only_selected_leading_triple": (
                    census.finite_only_selected_leading_triple_degenerate_count
                ),
                "finite_only_selected_full_four_node": (
                    census.finite_only_selected_full_degenerate_count
                ),
                "finite_only_conjugate_leading_triple": (
                    census.finite_only_conjugate_leading_triple_degenerate_count
                ),
                "finite_only_conjugate_full_four_node": (
                    census.finite_only_conjugate_full_degenerate_count
                ),
                "one_dotted_selected_leading_triple": (
                    census.one_dotted_selected_leading_triple_degenerate_count
                ),
                "one_dotted_selected_full_four_node": (
                    census.one_dotted_selected_full_degenerate_count
                ),
                "one_dotted_conjugate_leading_triple": (
                    census.one_dotted_conjugate_leading_triple_degenerate_count
                ),
                "one_dotted_conjugate_full_four_node": (
                    census.one_dotted_conjugate_full_degenerate_count
                ),
            },
        },
        "topology_counts": _topology_rows(census),
        "global_one_dotted_label3_path_family": {
            "scope": (
                "for each native rank 23--30, the connected path has one terminal "
                "dotted edge product c and every other nonzero edge has Coxeter label 3"
            ),
            "cycle_compatibility": (
                "vacuous but exact: a tree has no graph cycles, so every edge-product "
                "assignment is symmetrizable over Q(sqrt(5))"
            ),
            "cartan_realization": (
                "A_01=-1, A_10=-c, every later adjacent pair is (-1,-1), and all "
                "other off-diagonal entries vanish; symmetrizer diag(c,1,...,1)"
            ),
            "leading_determinants": "D_j=2*j-(j-1)*c",
            "selected_inertia": (
                "c>4 gives D_1>0 and D_j<0 for every j>=2, hence exactly one "
                "negative direction by exact LDL pivot signs"
            ),
            "conjugate_positive_definiteness": (
                "D_j'>0 for every j<=r iff c'<2*r/(r-1), since these thresholds "
                "strictly decrease with j"
            ),
            "rows": [
                {
                    "native_rank": row.native_rank,
                    "hyperbolic_dimension": row.hyperbolic_dimension,
                    "admissible_dotted_value_count": row.admissible_dotted_value_count,
                    "admissible_dotted_stream_sha256": (
                        row.admissible_dotted_stream_sha256
                    ),
                    "sample_dotted_product": list(row.sample_dotted_product),
                    "sample_full_determinant": list(row.sample_full_determinant),
                }
                for row in path_rows
            ],
            "uniform_sample": (
                "c=2+2*phi works in every rank 23--30; it gives integral Cartan "
                "reflections in an arithmetic Lorentzian lattice"
            ),
            "arithmetic_host_cocompactness": (
                "the symmetrized form is positive definite at the conjugate real place, "
                "so it has no nonzero isotropic vector over Q(sqrt(5)); hence its full "
                "arithmetic orthogonal lattice is uniform/cocompact"
            ),
            "minimal_vinberg_field": (
                "exactly Q(sqrt(5)): every dotted alphabet value is nonrational (a "
                "rational value cannot be both >4 selected and <4 conjugate), and the "
                "length-two cyclic Cartan product is c"
            ),
            "compactness_warning": (
                "these are exact full-rank seeds inside cocompact arithmetic orthogonal "
                "lattices, not compact Coxeter chambers. The displayed n+1-wall "
                "simplicial cone has a dotted proper subdiagram and is not a finite-volume "
                "simplex; a complete Vinberg run must decide whether the full reflection "
                "subgroup has finite index"
            ),
        },
        "strata": [
            {
                "tranche": row.tranche,
                "topology": row.topology,
                "nonzero_edge_count": row.nonzero_edge_count,
                "degree_sequence": list(row.degree_sequence),
                "cycle_rank": row.cycle_rank,
                "selected_leading_triple_sign": row.selected_leading_triple_sign,
                "selected_full_determinant_sign": row.selected_full_determinant_sign,
                "conjugate_leading_triple_sign": row.conjugate_leading_triple_sign,
                "conjugate_full_determinant_sign": row.conjugate_full_determinant_sign,
                "count": row.count,
            }
            for row in census.strata
        ],
        "stream_digests": {
            "finite_only_cycle_compatible_sha256": (
                census.finite_only_compatible_stream_sha256
            ),
            "finite_only_two_place_retained_sha256": (
                census.finite_only_retained_stream_sha256
            ),
            "one_dotted_cycle_compatible_sha256": (
                census.one_dotted_compatible_stream_sha256
            ),
            "one_dotted_two_place_retained_sha256": (
                census.one_dotted_retained_stream_sha256
            ),
            "record_format": (
                "six canonical edge products; leading 012 determinant; full determinant; "
                "each golden integer as decimal a,b, semicolon-separated and newline-terminated"
            ),
        },
        "outcome": {
            "complete_four_node_finite_only_tranche": True,
            "complete_four_node_exactly_one_dotted_tranche": True,
            "all_6869_dotted_values_in_input": True,
            "overlapping_triangle_and_four_cycle_compatibility_checked": True,
            "new_cocompact_arithmetic_host_lattice_family_with_integral_reflection_frames": True,
            "host_family_reflectivity_status": "OPEN",
            "new_compact_example": False,
            "new_high_rank_nonexistence_proof": False,
            "Qsqrt5_H22_H29_status": "OPEN",
        },
    }


def encoded_payload() -> bytes:
    return (json.dumps(build_payload(), indent=2, sort_keys=True) + "\n").encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    encoded = encoded_payload()
    if args.verify:
        if not OUTPUT.exists() or OUTPUT.read_bytes() != encoded:
            raise SystemExit("FAIL: frozen golden four-node JSON differs from recomputation")
        status = "PASS"
    else:
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_bytes(encoded)
        status = "WROTE"
    payload = json.loads(encoded)
    print(
        json.dumps(
            {
                "status": status,
                "output": str(OUTPUT),
                "sha256": sha256(encoded).hexdigest(),
                "finite_only_retained": payload["two_place_gram_gate"][
                    "finite_only_retained"
                ],
                "one_dotted_retained": payload["two_place_gram_gate"][
                    "one_dotted_retained"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
