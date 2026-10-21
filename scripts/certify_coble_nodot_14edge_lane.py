#!/usr/bin/env python3
"""Exhaust the 14-edge, degree-three, small-pairing Coble--Mukai lane.

The lane consists of connected ordinary support graphs on twelve vertices
with exactly fourteen edges and maximum degree three.  Zero through three
ordinary edges may be promoted from pairing one to pairing two.  Every exact
Lorentzian corank-two presentation is retained, every possible even
unimodular ``E10`` saturation is enumerated through the full discriminant
quadratic form, and every integral branch receives exact GMP double
description.  Literal Vinberg tests are run on every branch without a
hyperideal endpoint; the certified lane turns out to have none.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import shutil
from typing import Any

from scripts.certify_coble_nodot_13edge_lane import (
    _all_even_unimodular_extensions,
    _analyze_cone,
    _deduplicate_weighted,
    _integer_rows,
    _literal_vinberg,
    _matrix,
    _quotient_lattice,
)
from scripts.scout_coble_nodot_weighted_graphs import build as build_scout


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "enumerations" / "coble_nodot_14edge_lane.json"
SCHEMA = "coxeter.coble_nodot_14edge_lane.v1"


def _digest(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def _public_cone(cone: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in cone.items()
        if key not in {"rays", "cusp_audits", "facet_active_span_ranks"}
    }


def build(cdd: str | None = None) -> dict[str, Any]:
    cdd = cdd or shutil.which("cddexec_gmp") or "cddexec_gmp"
    scout = build_scout(
        maximum_degree=3,
        edge_range="14:14",
        maximum_parallel=3,
        sample_limit=1000,
        minimum_parallel=0,
        sample_all_exact=True,
    )
    raw = scout["exact_labeled_survivor_samples"]
    expected_trials = {
        "14:0": 5_291,
        "14:1": 74_074,
        "14:2": 481_481,
        "14:3": 1_925_924,
    }
    expected_exact = {"14:0": 2, "14:1": 12, "14:2": 90, "14:3": 132}
    if not (
        scout["underlying_graph_count"] == 5_291
        and scout["trial_count_by_edges_and_parallel"] == expected_trials
        and scout["exact_labeled_survivor_count_by_edges_and_parallel"]
        == expected_exact
        and scout["exact_labeled_survivor_count"] == 236
        and scout["exact_sample_filter"] == "all_exact_lorentzian_corank_two"
        and len(raw) == 236
    ):
        raise ArithmeticError("the exhaustive 14-edge scout denominator changed")

    classes = _deduplicate_weighted(raw)
    if len(classes) != 194:
        raise ArithmeticError("the 14-edge weighted-isomorphism quotient changed")
    class_records = []
    branch_records_flat = []
    for class_index, diagram_class in enumerate(classes):
        representative = diagram_class["representative"]
        root_gram = -_matrix(representative["negative_cartan"])
        quotient = _quotient_lattice(root_gram)
        extensions = _all_even_unimodular_extensions(
            quotient["gram"], quotient["root_coordinates"]
        )
        extension_records = []
        for extension_index, extension in enumerate(extensions):
            cone = _analyze_cone(
                extension["gram"], extension["root_coordinates"], cdd
            )
            vinberg = _literal_vinberg(
                extension["gram"], extension["root_coordinates"], cone
            )
            no_hyperideal = cone[
                "complete_cdd_V_representation_has_no_hyperideal_ray"
            ]
            valid_polyhedral_gate = (
                cone["all_input_walls_irredundant"]
                and cone["strict_interior_from_ray_sum"]
                and cone["strict_interior_square"] > 0
                and cone["all_causal_rays_same_time_cone"]
            )
            if no_hyperideal:
                finite_volume_gate = (
                    valid_polyhedral_gate
                    and cone[
                        "all_ideal_supports_parabolic_rank_n_minus_one"
                    ]
                    and vinberg["literal_vinberg_passes"]
                )
                classification = (
                    "compact"
                    if finite_volume_gate
                    and cone["compact_by_cdd_vertex_types"]
                    else "finite_volume_noncompact"
                    if finite_volume_gate
                    else "unresolved_no_hyperideal_but_failed_exact_gate"
                )
            else:
                finite_volume_gate = False
                classification = "infinite_volume_hyperideal"
            hyperideal_witness = next(
                (
                    record
                    for record in cone["rays"]
                    if record["type"] == "hyperideal"
                ),
                None,
            )
            branch = {
                "extension_index_zero_based": extension_index,
                "index_over_root_quotient": extension[
                    "index_over_root_quotient"
                ],
                "isotropic_subgroup": extension["isotropic_subgroup"],
                "e10_gram_sha256": _digest(_integer_rows(extension["gram"])),
                "root_coordinates_sha256": _digest(
                    _integer_rows(extension["root_coordinates"])
                ),
                "cone": _public_cone(cone),
                "literal_vinberg": vinberg,
                "finite_volume_gate_pass": finite_volume_gate,
                "classification": classification,
                "hyperideal_witness": hyperideal_witness,
                "full_cusp_lift_attempted": False,
                "full_cusp_lift_stop_reason": (
                    "the Coble--Mukai root chamber already has a hyperideal "
                    "extreme ray"
                    if not no_hyperideal
                    else None
                ),
                "split_K3_compiler_attempted": False,
                "split_K3_stop_reason": (
                    "D_1_split requires a finite-volume Coble--Mukai chamber"
                    if not finite_volume_gate
                    else None
                ),
            }
            extension_records.append(branch)
            branch_records_flat.append((class_index, branch))
        class_records.append(
            {
                "class_index_zero_based": class_index,
                "labelled_member_count": len(diagram_class["members"]),
                "parallel_pair_count": sum(
                    root_gram[left, right] == 2
                    for left in range(12)
                    for right in range(left + 1, 12)
                ),
                "negative_cartan": _integer_rows(-root_gram),
                "radical_quotient": {
                    "determinant": quotient["determinant"],
                    "smith_diagonal_absolute": quotient[
                        "smith_diagonal_absolute"
                    ],
                    "even_unimodular_E10_extension_count": len(extensions),
                },
                "extension_branches": extension_records,
            }
        )

    branch_count = len(branch_records_flat)
    classifications = Counter(
        branch["classification"] for _, branch in branch_records_flat
    )
    determinant_class_histogram = Counter(
        row["radical_quotient"]["determinant"] for row in class_records
    )
    determinant_branch_histogram = Counter(
        class_records[class_index]["radical_quotient"]["determinant"]
        for class_index, _ in branch_records_flat
    )
    hyperideal_histogram = Counter(
        branch["cone"]["hyperideal_vertex_count"]
        for _, branch in branch_records_flat
    )
    if not (
        branch_count == 169
        and classifications == {"infinite_volume_hyperideal": 169}
        and min(hyperideal_histogram) == 5
        and max(hyperideal_histogram) == 56
    ):
        raise ArithmeticError("the complete 14-edge volume conclusion changed")

    payload = {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "scope": {
            "vertices": 12,
            "connected_ordinary_support": True,
            "ordinary_support_edge_count": 14,
            "ordinary_support_maximum_degree": 3,
            "pairing_alphabet": [0, 1, 2],
            "parallel_pair_counts": [0, 1, 2, 3],
            "all_exact_Lorentzian_corank_two_determinants_retained": True,
            "all_even_unimodular_E10_extensions_retained": True,
            "not_all_14_edge_nodot": (
                "ordinary degree above three and four or more pairing-2 pairs "
                "are outside the declared finite lane"
            ),
            "not_the_full_D_1_split_nodot_denominator": True,
        },
        "denominator": {
            "underlying_unlabelled_graph_count": scout[
                "underlying_graph_count"
            ],
            "weighted_trial_count": sum(
                scout["trial_count_by_edges_and_parallel"].values()
            ),
            "trial_count_by_parallel_pairs": scout[
                "trial_count_by_edges_and_parallel"
            ],
            "exact_Lorentzian_corank_two_labelled_count": len(raw),
            "exact_count_by_parallel_pairs": scout[
                "exact_labeled_survivor_count_by_edges_and_parallel"
            ],
            "weighted_diagram_isometry_class_count": len(classes),
            "E10_extendable_weighted_class_count": sum(
                bool(row["extension_branches"]) for row in class_records
            ),
            "integral_E10_extension_branch_count": branch_count,
            "determinant_histogram_on_weighted_classes": {
                str(key): value
                for key, value in sorted(determinant_class_histogram.items())
            },
            "determinant_histogram_on_E10_branches": {
                str(key): value
                for key, value in sorted(determinant_branch_histogram.items())
            },
        },
        "class_records": class_records,
        "volume_ledger": {
            "compact_branch_count": classifications["compact"],
            "finite_volume_noncompact_branch_count": classifications[
                "finite_volume_noncompact"
            ],
            "infinite_volume_hyperideal_branch_count": classifications[
                "infinite_volume_hyperideal"
            ],
            "hyperideal_vertex_count_histogram": {
                str(key): value
                for key, value in sorted(hyperideal_histogram.items())
            },
            "minimum_hyperideal_vertex_count": min(hyperideal_histogram),
            "maximum_hyperideal_vertex_count": max(hyperideal_histogram),
            "finite_volume_survivors": [],
        },
        "theorem": {
            "statement": (
                "Every one of the 169 integral E10 branches in the declared "
                "14-edge lane has a hyperideal extreme ray and therefore "
                "infinite volume; there are no compact or finite-volume "
                "noncompact survivors."
            ),
            "literal_Vinberg_policy": (
                "the literal edge/endpoint and affine-zero-neighborhood audit "
                "is run on every no-hyperideal branch; there are zero such "
                "branches, while each of the 169 branches is already rejected "
                "by an explicit cdd hyperideal endpoint"
            ),
            "AG_or_split_K3_pushback": (
                "none is logically reached because no raw root chamber passes "
                "the finite-volume gate"
            ),
        },
        "claims": {
            "new_hyperbolic_example": False,
            "strict_L4": False,
            "all_nodot": False,
            "scoped_negative_classification": True,
        },
        "all_checks_pass": True,
    }
    return payload


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
