#!/usr/bin/env python3
"""Certify the 15-edge, no-parallel Coble--Mukai root lane.

The surface-free denominator consists of connected ordinary support graphs on
twelve roots, with fifteen pairing-one edges, no pairing-two edge, and maximum
ordinary degree three.  Every Lorentzian corank-two Gram is retained, followed
by every even unimodular E10 extension, exact GMP cone enumeration, literal
Vinberg tests, the complete one-boundary cusp lift, and comparison with the
two archived strong-Coble return packages.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import shutil
from typing import Any

import sympy as sp

from scripts.certify_coble_nodot_13edge_lane import (
    _all_even_unimodular_extensions,
    _analyze_cone,
    _deduplicate_weighted,
    _full_cusp_lift,
    _full_isomorphic,
    _integer_rows,
    _known_root_graphs,
    _literal_vinberg,
    _matrix,
    _quotient_lattice,
    _weighted_isomorphic,
)
from scripts.scout_coble_nodot_weighted_graphs import build as build_scout


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "data" / "enumerations" / "coble_nodot_15edge_lane.json"
SOURCE_FULL = (
    ROOT / "data" / "enumerations" / "coble_one_boundary_full_ns_h10.json"
)
SOURCE_PROMOTION = (
    ROOT
    / "data"
    / "enumerations"
    / "coble_one_boundary_raw_k3_promotion.json"
)
SOURCE_SHA256 = {
    SOURCE_FULL: "1e90fcd5f9b1ba8538937b91d2f26ea401f6952e8484a597624939edaf2a873e",
    SOURCE_PROMOTION: "1c8b8699da8f8747fdad22d520743bb4640b41b8e988c16844ecfe5e27c1958c",
}
SCHEMA = "coxeter.coble_nodot_15edge_lane.v1"


def _public_cone(cone: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in cone.items()
        if key not in {"rays", "cusp_audits", "facet_active_span_ranks"}
    }


def _digest(value: Any) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def build(cdd: str | None = None) -> dict[str, Any]:
    cdd = cdd or shutil.which("cddexec_gmp") or "cddexec_gmp"
    if any(
        sha256(path.read_bytes()).hexdigest() != expected
        for path, expected in SOURCE_SHA256.items()
    ):
        raise ArithmeticError("a pinned one-boundary source archive changed")
    scout = build_scout(
        maximum_degree=3,
        edge_range="15:15",
        maximum_parallel=0,
        sample_limit=100,
        minimum_parallel=0,
        sample_all_exact=True,
    )
    raw = scout["exact_labeled_survivor_samples"]
    if not (
        scout["underlying_graph_count"] == 5_741
        and scout["trial_count_by_edges_and_parallel"] == {"15:0": 5_741}
        and scout["exact_labeled_survivor_count_by_edges_and_parallel"]
        == {"15:0": 2}
        and scout["exact_labeled_survivor_count"] == 2
        and scout["radical_quotient_determinant_histogram"]
        == {"15:0:-4": 1, "15:0:-1": 1}
        and len(raw) == 2
    ):
        raise ArithmeticError("the exhaustive 15-edge scout denominator changed")

    classes = _deduplicate_weighted(raw)
    if len(classes) != 2:
        raise ArithmeticError("the two 15-edge Grams became isomorphic")
    known_root_graphs = _known_root_graphs()
    known_full = json.loads(SOURCE_FULL.read_text(encoding="utf-8"))

    class_records: list[dict[str, Any]] = []
    live_full: list[dict[str, Any]] = []
    branch_classifications: Counter[str] = Counter()
    for class_index, diagram_class in enumerate(classes):
        representative = diagram_class["representative"]
        root_gram = -_matrix(representative["negative_cartan"])
        quotient = _quotient_lattice(root_gram)
        extensions = _all_even_unimodular_extensions(
            quotient["gram"], quotient["root_coordinates"]
        )
        extension_records = []
        for extension_index, extension in enumerate(extensions):
            cm = _analyze_cone(
                extension["gram"], extension["root_coordinates"], cdd
            )
            cm_vinberg = _literal_vinberg(
                extension["gram"], extension["root_coordinates"], cm
            )
            polyhedral_gate = (
                cm["all_input_walls_irredundant"]
                and cm["strict_interior_from_ray_sum"]
                and cm["strict_interior_square"] > 0
                and cm["all_causal_rays_same_time_cone"]
            )
            finite_volume = (
                polyhedral_gate
                and cm["complete_cdd_V_representation_has_no_hyperideal_ray"]
                and cm["all_ideal_supports_parabolic_rank_n_minus_one"]
                and cm_vinberg["literal_vinberg_passes"]
            )
            if finite_volume:
                classification = (
                    "compact"
                    if cm["compact_by_cdd_vertex_types"]
                    else "finite_volume_noncompact"
                )
            elif cm["hyperideal_vertex_count"]:
                classification = "infinite_volume_hyperideal"
            else:
                classification = "failed_exact_volume_gate"
            branch_classifications[classification] += 1

            full_record = None
            if finite_volume:
                full = _full_cusp_lift(
                    extension["gram"], extension["root_coordinates"], cm, cdd
                )
                full_cone = full["cone"]
                full_vinberg = full["literal_vinberg"]
                full_pass = (
                    full_cone["all_input_walls_irredundant"]
                    and full_cone[
                        "complete_cdd_V_representation_has_no_hyperideal_ray"
                    ]
                    and full_cone["strict_interior_from_ray_sum"]
                    and full_cone["strict_interior_square"] > 0
                    and full_cone["all_causal_rays_same_time_cone"]
                    and full_cone[
                        "all_ideal_supports_parabolic_rank_n_minus_one"
                    ]
                    and full_vinberg["literal_vinberg_passes"]
                    and full["characteristic_wall_positions_zero_based"] == [0]
                    and full["wall_coordinate_span_index_in_I_1_10"] == 1
                )
                full_matches = [
                    f"Type {case_key}"
                    for case_key, case in known_full["cases"].items()
                    if _full_isomorphic(
                        full["wall_gram"], _matrix(case["negative_wall_gram"])
                    )
                ]
                full_record = {
                    "gate_pass": full_pass,
                    "wall_count": full["normal_coordinates"].cols,
                    "known_full_chamber_matches": full_matches,
                    "wall_coordinate_span_index_in_I_1_10": full[
                        "wall_coordinate_span_index_in_I_1_10"
                    ],
                    "characteristic_wall_positions_zero_based": full[
                        "characteristic_wall_positions_zero_based"
                    ],
                    "wall_gram_sha256": _digest(_integer_rows(full["wall_gram"])),
                    "cone": _public_cone(full_cone),
                    "literal_vinberg": full_vinberg,
                }
                if full_pass:
                    live_full.append(
                        {
                            "class_index_zero_based": class_index,
                            "extension_index_zero_based": extension_index,
                            "matches": full_matches,
                            "wall_gram": full["wall_gram"],
                        }
                    )

            extension_records.append(
                {
                    "extension_index_zero_based": extension_index,
                    "index_over_root_quotient": extension[
                        "index_over_root_quotient"
                    ],
                    "isotropic_subgroup": extension["isotropic_subgroup"],
                    "e10_gram_sha256": _digest(_integer_rows(extension["gram"])),
                    "root_coordinates_sha256": _digest(
                        _integer_rows(extension["root_coordinates"])
                    ),
                    "cm_cone": _public_cone(cm),
                    "cm_literal_vinberg": cm_vinberg,
                    "classification": classification,
                    "full_cusp_lift": full_record,
                }
            )

        known_root_matches = [
            label
            for label, graph in known_root_graphs.items()
            if _weighted_isomorphic(diagram_class["graph"], graph)
        ]
        class_records.append(
            {
                "class_index_zero_based": class_index,
                "source_graph6": representative["underlying_graph6"],
                "labelled_member_count": len(diagram_class["members"]),
                "known_root_diagram_matches": known_root_matches,
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

    promotion = json.loads(SOURCE_PROMOTION.read_text(encoding="utf-8"))
    promotion_pass = {
        str(case["raw_id"]): bool(case["all_checks_pass"])
        for case in promotion["cases"]
    }
    survivor_ledger = []
    for survivor in live_full:
        matches = survivor["matches"]
        if matches == ["Type I"]:
            raw_id = "raw_22"
        elif matches == ["Type II"]:
            raw_id = "raw_30"
        else:
            raw_id = None
        survivor_ledger.append(
            {
                key: survivor[key]
                for key in ("class_index_zero_based", "extension_index_zero_based")
            }
            | {
                "known_full_chamber_matches": matches,
                "full_split_K3_Coble_return_pass": (
                    promotion_pass.get(raw_id) if raw_id is not None else None
                ),
                "return_scope": (
                    f"inherited exact {raw_id} compiler after literal full-Gram match"
                    if raw_id is not None
                    else "no known-return inference"
                ),
            }
        )

    full_isomorphism_classes: list[sp.Matrix] = []
    for survivor in live_full:
        wall_gram = survivor["wall_gram"]
        if not any(
            _full_isomorphic(wall_gram, representative)
            for representative in full_isomorphism_classes
        ):
            full_isomorphism_classes.append(wall_gram)

    payload = {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "scope": {
            "vertices": 12,
            "connected_ordinary_support": True,
            "ordinary_support_edge_count": 15,
            "ordinary_support_maximum_degree": 3,
            "pairing_alphabet": [0, 1],
            "parallel_pair_count": 0,
            "all_exact_Lorentzian_corank_two_determinants_retained": True,
            "all_even_unimodular_E10_extensions_retained": True,
            "not_the_full_D_1_split_nodot_denominator": True,
        },
        "denominator": {
            "underlying_unlabelled_graph_count": scout[
                "underlying_graph_count"
            ],
            "weighted_trial_count": sum(
                scout["trial_count_by_edges_and_parallel"].values()
            ),
            "exact_Lorentzian_corank_two_labelled_count": len(raw),
            "weighted_diagram_isometry_class_count": len(classes),
            "integral_E10_extension_branch_count": sum(
                len(record["extension_branches"]) for record in class_records
            ),
            "radical_quotient_determinant_histogram": scout[
                "radical_quotient_determinant_histogram"
            ],
        },
        "class_records": class_records,
        "volume_ledger": dict(sorted(branch_classifications.items())),
        "full_D_1_split_branch_count": len(live_full),
        "full_D_1_split_isomorphism_class_count": len(full_isomorphism_classes),
        "full_survivor_ledger": survivor_ledger,
        "claims": {
            "strict_L4_for_full_D_1_split_nodot": False,
            "reason": (
                "this theorem exhausts only the 15-edge, no-parallel, "
                "maximum-degree-three lane"
            ),
            "full_split_compiler_only_inferred_after_literal_known_gram_match": True,
        },
    }
    payload["all_checks_pass"] = (
        len(class_records) == 2
        and sum(row["labelled_member_count"] for row in class_records) == 2
        and sorted(
            (
                row["radical_quotient"]["determinant"],
                row["radical_quotient"]["even_unimodular_E10_extension_count"],
            )
            for row in class_records
        )
        == [(-4, 2), (-1, 1)]
        and branch_classifications
        == {"finite_volume_noncompact": 2, "infinite_volume_hyperideal": 1}
        and len(live_full) == 2
        and len(full_isomorphism_classes) == 1
        and all(row["known_full_chamber_matches"] == ["Type II"] for row in survivor_ledger)
        and all(
            row["full_split_K3_Coble_return_pass"] is True
            for row in survivor_ledger
        )
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
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
