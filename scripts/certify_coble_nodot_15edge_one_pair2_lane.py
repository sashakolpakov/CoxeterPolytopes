#!/usr/bin/env python3
"""Certify the 15-edge Coble--Mukai lane with one pairing-two edge.

The P-only denominator consists of connected ordinary support graphs on
twelve roots, with exactly fifteen support edges, maximum support degree
three, exactly one root pairing equal to two, and all fourteen other
pairings equal to one.  Every marked edge of every unlabelled support is
tested by exact rank and inertia arithmetic.  Every even unimodular E10
extension of every Lorentzian corank-two quotient is then sent through exact
GMP cone enumeration and literal Vinberg tests.  Finite-volume branches
receive the complete one-boundary cusp lift and same-P split-K3/Coble return
test when their full Gram has an archived compiler.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any

from flint import fmpz_mat
import networkx as nx
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
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "data"
    / "enumerations"
    / "coble_nodot_15edge_one_pair2_lane.json"
)
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
SCHEMA = "coxeter.coble_nodot_15edge_one_pair2_lane.v1"


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


def _ray_hashes(cone: dict[str, Any]) -> dict[str, str]:
    return {
        "extreme_ray_records_sha256": _digest(cone["rays"]),
        "hyperideal_ray_records_sha256": _digest(
            [record for record in cone["rays"] if record["type"] == "hyperideal"]
        ),
    }


def _volume_classification(
    cone: dict[str, Any], vinberg: dict[str, Any]
) -> tuple[str, bool]:
    polyhedral_gate = (
        cone["all_input_walls_irredundant"]
        and cone["strict_interior_from_ray_sum"]
        and cone["strict_interior_square"] > 0
        and cone["all_causal_rays_same_time_cone"]
    )
    finite_volume = (
        polyhedral_gate
        and cone["complete_cdd_V_representation_has_no_hyperideal_ray"]
        and cone["all_ideal_supports_parabolic_rank_n_minus_one"]
        and vinberg["literal_vinberg_passes"]
    )
    if finite_volume:
        return (
            "compact"
            if cone["compact_by_cdd_vertex_types"]
            else "finite_volume_noncompact"
        ), True
    if cone["hyperideal_vertex_count"]:
        return "infinite_volume_hyperideal", False
    return "failed_exact_volume_gate", False


def _marked_cartan(graph: nx.Graph, promoted: tuple[int, int]) -> list[list[int]]:
    matrix = [[0] * 12 for _ in range(12)]
    for index in range(12):
        matrix[index][index] = 2
    for left, right in graph.edges():
        matrix[left][right] = matrix[right][left] = -1
    left, right = promoted
    matrix[left][right] = matrix[right][left] = -2
    return matrix


def _exact_denominator(geng: str) -> dict[str, Any]:
    """Enumerate every support and every marked edge without a float filter."""

    process = subprocess.Popen(
        [geng, "-q", "-c", "-D3", "12", "15:15"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdout is not None
    graph_count = 0
    trial_count = 0
    rank_histogram: Counter[int] = Counter()
    rank_ten_inertia_histogram: Counter[tuple[int, int, int]] = Counter()
    determinant_histogram: Counter[int] = Counter()
    survivors: list[dict[str, Any]] = []
    for graph_index, raw in enumerate(process.stdout):
        graph6 = raw.strip().decode()
        graph = nx.from_graph6_bytes(raw.strip())
        graph_count += 1
        edges = tuple(tuple(sorted(edge)) for edge in graph.edges())
        if len(edges) != 15 or max(dict(graph.degree()).values()) > 3:
            raise ArithmeticError("geng emitted a support outside the denominator")
        for promoted in edges:
            trial_count += 1
            matrix_rows = _marked_cartan(graph, promoted)
            rank = int(fmpz_mat(matrix_rows).rank())
            rank_histogram[rank] += 1
            if rank != 10:
                continue
            matrix = sp.Matrix(matrix_rows)
            inertia = exact_inertia(matrix)
            rank_ten_inertia_histogram[inertia] += 1
            if inertia != (9, 1, 2):
                continue
            quotient = _quotient_lattice(-matrix)
            determinant = int(quotient["determinant"])
            determinant_histogram[determinant] += 1
            survivors.append(
                {
                    "source_graph_index_zero_based": graph_index,
                    "underlying_graph6": graph6,
                    "underlying_edge_count": 15,
                    "promoted_parallel_edges_zero_based": [list(promoted)],
                    "radical_quotient_determinant": determinant,
                    "negative_cartan": matrix_rows,
                }
            )
    stderr = process.stderr.read().decode() if process.stderr is not None else ""
    return_code = process.wait()
    if return_code:
        raise RuntimeError(stderr or f"geng exited with status {return_code}")
    if trial_count != sum(rank_histogram.values()):
        raise ArithmeticError("rank histogram does not cover every marked trial")
    return {
        "underlying_graph_count": graph_count,
        "weighted_trial_count": trial_count,
        "exact_rank_histogram": {
            str(rank): count for rank, count in sorted(rank_histogram.items())
        },
        "rank_ten_exact_inertia_histogram": {
            ",".join(map(str, inertia)): count
            for inertia, count in sorted(rank_ten_inertia_histogram.items())
        },
        "radical_quotient_determinant_histogram": {
            str(determinant): count
            for determinant, count in sorted(determinant_histogram.items())
        },
        "survivors": survivors,
    }


def build(cdd: str | None = None, geng: str | None = None) -> dict[str, Any]:
    cdd = cdd or shutil.which("cddexec_gmp") or "cddexec_gmp"
    geng = geng or shutil.which("geng") or "geng"
    if any(
        sha256(path.read_bytes()).hexdigest() != expected
        for path, expected in SOURCE_SHA256.items()
    ):
        raise ArithmeticError("a pinned one-boundary source archive changed")

    denominator = _exact_denominator(geng)
    raw = denominator.pop("survivors")
    classes = _deduplicate_weighted(raw)
    known_root_graphs = _known_root_graphs()
    known_full = json.loads(SOURCE_FULL.read_text(encoding="utf-8"))

    class_records: list[dict[str, Any]] = []
    live_full: list[dict[str, Any]] = []
    root_cone_classifications: Counter[str] = Counter()
    branch_classifications: Counter[str] = Counter()
    for class_index, diagram_class in enumerate(classes):
        representative = diagram_class["representative"]
        root_gram = -_matrix(representative["negative_cartan"])
        quotient = _quotient_lattice(root_gram)
        root_cone = _analyze_cone(
            quotient["gram"], quotient["root_coordinates"], cdd
        )
        root_vinberg = _literal_vinberg(
            quotient["gram"], quotient["root_coordinates"], root_cone
        )
        root_classification, _ = _volume_classification(root_cone, root_vinberg)
        root_cone_classifications[root_classification] += 1
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
            classification, finite_volume = _volume_classification(cone, vinberg)
            branch_classifications[classification] += 1

            full_record = None
            if finite_volume:
                full = _full_cusp_lift(
                    extension["gram"], extension["root_coordinates"], cone, cdd
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
                    "cone": _public_cone(cone),
                    "cone_ray_hashes": _ray_hashes(cone),
                    "literal_vinberg": vinberg,
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
                "promoted_pair_zero_based": representative[
                    "promoted_parallel_edges_zero_based"
                ][0],
                "marked_trial_member_count": len(diagram_class["members"]),
                "known_root_diagram_matches": known_root_matches,
                "negative_cartan": _integer_rows(-root_gram),
                "radical_quotient": {
                    "determinant": quotient["determinant"],
                    "smith_diagonal_absolute": quotient[
                        "smith_diagonal_absolute"
                    ],
                    "even_unimodular_E10_extension_count": len(extensions),
                },
                "root_cone": _public_cone(root_cone),
                "root_cone_ray_hashes": _ray_hashes(root_cone),
                "root_literal_vinberg": root_vinberg,
                "root_cone_classification": root_classification,
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
                "class_index_zero_based": survivor["class_index_zero_based"],
                "extension_index_zero_based": survivor[
                    "extension_index_zero_based"
                ],
                "known_full_chamber_matches": matches,
                "full_split_K3_Coble_return_pass": (
                    promotion_pass.get(raw_id) if raw_id is not None else None
                ),
                "return_scope": (
                    f"inherited exact {raw_id} compiler after literal full-Gram match"
                    if raw_id is not None
                    else "no archived compiler; same-P surface return not established"
                ),
            }
        )

    full_isomorphism_classes: list[sp.Matrix] = []
    for survivor in live_full:
        wall_gram = survivor["wall_gram"]
        if not any(
            _full_isomorphic(wall_gram, other)
            for other in full_isomorphism_classes
        ):
            full_isomorphism_classes.append(wall_gram)

    nonsquare_determinant_classes = sum(
        math.isqrt(abs(record["radical_quotient"]["determinant"])) ** 2
        != abs(record["radical_quotient"]["determinant"])
        for record in class_records
    )
    square_determinant_classes_without_extension = sum(
        math.isqrt(abs(record["radical_quotient"]["determinant"])) ** 2
        == abs(record["radical_quotient"]["determinant"])
        and not record["extension_branches"]
        for record in class_records
    )

    payload = {
        "schema": SCHEMA,
        "status_date": "2026-08-09",
        "scope": {
            "vertices": 12,
            "connected_ordinary_support": True,
            "ordinary_support_edge_count": 15,
            "ordinary_support_maximum_degree": 3,
            "pairing_alphabet": [0, 1, 2],
            "pairing_two_edge_count": 1,
            "all_other_nonzero_pairings_equal_one": True,
            "all_marked_edges_of_every_unlabelled_support_tested": True,
            "weighted_trials_are_edge_markings_not_vertex_labellings": True,
            "exact_rank_on_every_marked_trial_without_float_prefilter": True,
            "all_Lorentzian_corank_two_determinants_retained": True,
            "all_even_unimodular_E10_extensions_retained": True,
            "not_the_full_D_1_split_nodot_denominator": True,
        },
        "denominator": denominator
        | {
            "exact_Lorentzian_corank_two_marked_trial_count": len(raw),
            "weighted_diagram_isometry_class_count": len(classes),
            "integral_E10_extension_branch_count": sum(
                len(record["extension_branches"]) for record in class_records
            ),
            "classes_without_even_unimodular_E10_extension": sum(
                not record["extension_branches"] for record in class_records
            ),
            "classes_rejected_by_nonsquare_quotient_determinant": (
                nonsquare_determinant_classes
            ),
            "square_determinant_classes_without_E10_extension": (
                square_determinant_classes_without_extension
            ),
        },
        "class_records": class_records,
        "root_cone_volume_ledger": dict(sorted(root_cone_classifications.items())),
        "E10_extension_volume_ledger": dict(
            sorted(branch_classifications.items())
        ),
        "full_D_1_split_branch_count": len(live_full),
        "full_D_1_split_isomorphism_class_count": len(full_isomorphism_classes),
        "full_survivor_ledger": survivor_ledger,
        "claims": {
            "strict_L4_for_full_D_1_split_nodot": False,
            "new_abstract_hyperbolic_polytope_claimed": False,
            "reason": (
                "this theorem exhausts only the 15-edge, exactly-one-pairing-two, "
                "maximum-degree-three lane"
            ),
            "full_split_compiler_only_inferred_after_literal_known_gram_match": True,
            "classes_without_E10_saturation_are_lattice_obstructions": True,
            "hyperideal_branches_are_infinite_volume": True,
            "every_Lorentzian_root_cone_has_a_hyperideal_ray": True,
            "no_finite_volume_survivor_in_the_scoped_lane": True,
        },
    }
    payload["all_checks_pass"] = (
        denominator["underlying_graph_count"] == 5_741
        and denominator["weighted_trial_count"] == 86_115
        and len(raw) == 29
        and len(classes) == 28
        and sum(record["marked_trial_member_count"] for record in class_records)
        == 29
        and sum(
            len(record["extension_branches"]) for record in class_records
        )
        == 21
        and sum(not record["extension_branches"] for record in class_records)
        == 10
        and nonsquare_determinant_classes == 10
        and square_determinant_classes_without_extension == 0
        and root_cone_classifications == {"infinite_volume_hyperideal": 28}
        and branch_classifications == {"infinite_volume_hyperideal": 21}
        and all(
            record["root_cone"]["all_input_walls_irredundant"]
            and record["root_cone"]["strict_interior_from_ray_sum"]
            and record["root_cone"]["strict_interior_square"] > 0
            and record["root_cone"]["all_causal_rays_same_time_cone"]
            and record["root_cone"]["hyperideal_vertex_count"] > 0
            for record in class_records
        )
        and all(
            (
                branch["cone"]["finite_vertex_count"],
                branch["cone"]["ideal_vertex_count"],
                branch["cone"]["hyperideal_vertex_count"],
            )
            == (
                record["root_cone"]["finite_vertex_count"],
                record["root_cone"]["ideal_vertex_count"],
                record["root_cone"]["hyperideal_vertex_count"],
            )
            for record in class_records
            for branch in record["extension_branches"]
        )
        and not live_full
        and not full_isomorphism_classes
        and not survivor_ledger
    )
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cdd", default=None)
    parser.add_argument("--geng", default=None)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build(args.cdd, args.geng)
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
