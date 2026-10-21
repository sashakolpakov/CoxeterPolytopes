#!/usr/bin/env python3
"""Freeze the rank >= 11 Roulleau K3-nef source denominator.

The certificate is deliberately a source inventory, not yet a transcription of
the fourteen rank 11--14 wall matrices.  With ``--check-source`` it checks the
pinned EPIGA TeX and ancillary Magma files in /private/tmp against their hashes
and exact line-local witnesses.  The default mode is reproducible without the
copyrighted source archive and writes the frozen manifest used by the tests.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "data/enumerations/roulleau_rank11_14_source_stock_inventory.json"
)
DEFAULT_TEX = Path("/private/tmp/roulleau_source/epiga_vol6_Roulleau.tex")
DEFAULT_EXAMPLES = Path("/private/tmp/roulleau_source/anc/ExamplesMagma.txt")
DEFAULT_ALGO = Path("/private/tmp/roulleau_source/anc/FiniteAutAlgoMagma.txt")

SOURCE_HASHES = {
    "pdf": "a8514b074cd955d487c1aa4e8798ce8d29c3d9274e004ba2b0a1bafdf922ad87",
    "source_tar": "3e445ed01446b47ded6c007b395308de34e566b736f0ba982470bc497adeb1e9",
    "tex": "e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3",
    "examples_magma": "ba2bddee3f5946d4603ab76e2a1db7039738f8bcd64dfe60b3b2c5481676fe8b",
    "algorithm_magma": "79c60e186a603508558e610141e3120d8822ce3d87158b11f2f8e390653fe02e",
}


def _case(
    key: str,
    rank: int,
    lattice: str,
    tex_token: str,
    walls: int,
    discriminant: int,
    aut: str,
    action_tag: str,
    section_lines: tuple[int, int],
    table_lines: tuple[int, int],
    extraction: str,
    completeness: str,
    ledger_status: str,
    ledger_evidence: str,
    ancillary_lines: tuple[int, int] | None = None,
    tex_table_token: str | None = None,
) -> dict:
    return {
        "key": key,
        "picard_rank": rank,
        "hyperbolic_dimension": rank - 1,
        "lattice": lattice,
        "tex_lattice_token": tex_token,
        "tex_table_lattice_token": tex_table_token or tex_token,
        "absolute_discriminant": discriminant,
        "wall_count": walls,
        "general_surface_automorphism_group": aut,
        "roulleau_ns_action_tag": action_tag,
        "source": {
            "section_lines": list(section_lines),
            "table_lines": list(table_lines),
            "ancillary_examples_magma_lines": (
                list(ancillary_lines) if ancillary_lines else None
            ),
        },
        "wall_data_extraction": extraction,
        "exact_wall_gram_extractable": True,
        "source_completeness_method": completeness,
        "reconstruction_status": (
            "already_exactly_archived"
            if ledger_status == "represented_exact_full_nef_stock"
            else "source_complete_transcription_and_exact_lattice_checks_pending"
        ),
        "current_machine_ledger_status": ledger_status,
        "ledger_evidence": ledger_evidence,
    }


CASES = [
    _case(
        "R11_UE8A1", 11, "U+E8(-1)+A1(-1)", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{A}_{1}",
        12, 2, "C2", "t", (8053, 8119), (10908, 10908),
        "complete_weighted_dual_graph", "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock", "exact full stock, trivial weighted group, and intrinsic full-chamber lineage in data/enumerations/roulleau_rank11_direct_graph_denominator.json",
    ),
    _case(
        "R11_UD8A1", 11, "U+D8(-1)+A1(-1)", r"U\oplus\mathbf{D}_{8}\oplus\mathbf{A}_{1}",
        14, 8, "C2", "t", (8120, 8205), (10910, 10910),
        "complete_weighted_dual_graph", "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock", "exact full stock and complete C2 denominator; the proper fixed rank is 10, so no fold above H9, in data/enumerations/roulleau_rank11_direct_graph_denominator.json",
    ),
    _case(
        "R11_UD4D4A1", 11, "U+D4(-1)^2+A1(-1)", r"U\oplus\mathbf{D}_{4}\oplus\mathbf{D}_{4}\oplus\mathbf{A}_{1}",
        22, 32, "C2", "t", (8206, 8318), (10912, 10912),
        "basis_weighted_graph_plus_linear_equivalences_for_all_extra_curves",
        "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock",
        "ledger records the prior 22-wall K3 nef chamber on U+D4(-1)^2+A1(-1)",
        tex_table_token=r"U\oplus\mathbf{D}_{4}^{\oplus2}\oplus\mathbf{A}_{1}",
    ),
    _case(
        "R11_UD4A1x5", 11, "U+D4(-1)+A1(-1)^5", r"U\oplus\mathbf{D}_{4}\oplus\mathbf{A}_{1}^{\oplus5}",
        90, 128, "C2xC2", "tn", (8319, 8542), (10914, 10914),
        "basis_weighted_graph_plus_closed_subset_formula_for_all_70_extra_curves",
        "induced_chamber_completion_inside_II_1_17",
        "represented_exact_full_nef_stock", "exact 90-wall matrix, S8 x C2 weighted group, explicit cusp, and complete no-fold-above-H9 proof in data/enumerations/roulleau_r11_ud4a1x5_90_denominator.json",
        (76, 142),
    ),
    _case(
        "R12_UE8A1x2", 12, "U+E8(-1)+A1(-1)^2", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{A}_{1}^{\oplus2}",
        14, 4, "C2", "t", (8545, 8647), (10922, 10922),
        "complete_weighted_dual_graph", "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock", "exact full stock, complete C2 subgroup denominator, and one cofinite fold in data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json",
    ),
    _case(
        "R12_UD8A1x2", 12, "U+D8(-1)+A1(-1)^2", r"U\oplus\mathbf{D}_{8}\oplus\mathbf{A}_{1}^{\oplus2}",
        19, 16, "C2", "t", (8648, 8803), (10924, 10924),
        "complete_weighted_dual_graph_split_across_two_tikz_panels",
        "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock", "exact full stock, all 10 D12 subgroup classes, and two cofinite folds in data/enumerations/roulleau_rank12_graph_all_subgroup_folds.json",
    ),
    _case(
        "R12_UD4x2A1x2", 12, "U+D4(-1)^2+A1(-1)^2", r"U\oplus\mathbf{D}_{4}^{\oplus2}\oplus\mathbf{A}_{1}^{\oplus2}",
        59, 64, "C2xC2", "tn", (8804, 9033), (10926, 10926),
        "basis_weighted_graph_plus_closed_Bij_Cij_formulas_for_all_40_extra_curves",
        "induced_chamber_completion_inside_II_1_17",
        "represented_exact_full_nef_stock",
        "exact 59-wall matrix, C2^2 x S6 weighted group, all 914 subgroup classes, and two cofinite H10 folds in data/enumerations/roulleau_r12_ud4x2a1x2_59_all_subgroup_folds.json",
    ),
    _case(
        "R12_UA2E8", 12, "U+A2(-1)+E8(-1)", r"U\oplus\mathbf{A}_{2}\oplus\mathbf{E}_{8}",
        13, 3, "C2", "n", (9034, 9125), (10928, 10928),
        "complete_weighted_dual_graph", "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock", "exact full stock and geometric involution fold in data/enumerations/roulleau_non2elementary_involution_folds.json",
    ),
    _case(
        "R13_UE8A1x3", 13, "U+E8(-1)+A1(-1)^3", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{A}_{1}^{\oplus3}",
        17, 8, "C2", "t", (9129, 9249), (10937, 10937),
        "complete_weighted_dual_graph", "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock", "exact full stock, all four S3 subgroup classes, and two known full reflection chambers in data/enumerations/roulleau_ue8a1x3_all_subgroup_folds.json",
    ),
    _case(
        "R13_UD8A1x3", 13, "U+D8(-1)+A1(-1)^3", r"U\oplus\mathbf{D}_{8}\oplus\mathbf{A}_{1}^{\oplus3}",
        39, 32, "C2xC2", "tn", (9250, 9490), (10941, 10941),
        "basis_weighted_graph_plus_closed_AJ_BJ_formula_for_all_20_extra_curves",
        "induced_chamber_completion_inside_II_1_17",
        "represented_exact_full_nef_stock",
        "exact coordinate/formula stock, all 57 diagram subgroup classes, and two cofinite fold types in data/enumerations/roulleau_r13_ud8a1x3_all_subgroup_folds.json",
        (145, 214),
    ),
    _case(
        "R13_UE8A3", 13, "U+E8(-1)+A3(-1)", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{A}_{3}",
        14, 4, "C2", "n", (9491, 9583), (10945, 10945),
        "complete_weighted_dual_graph", "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock", "exact full stock and geometric involution fold in data/enumerations/roulleau_non2elementary_involution_folds.json",
    ),
    _case(
        "R14_UE8D4", 14, "U+E8(-1)+D4(-1)", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{D}_{4}",
        15, 4, "C2", "t", (9586, 9668), (10957, 10957),
        "complete_weighted_dual_graph", "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock",
        "exact full stock and all four diagram subgroup classes in data/enumerations/roulleau_ue8d4_all_subgroup_folds.json; both proper fold diagrams are known full chambers",
    ),
    _case(
        "R14_UD8D4", 14, "U+D8(-1)+D4(-1)", r"U\oplus\mathbf{D}_{8}\oplus\mathbf{D}_{4}",
        20, 16, "C2xC2", "tn", (9669, 9827), (10962, 10962),
        "complete_weighted_dual_graph", "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock", "exact full stock, all 57 diagram subgroup classes, and five cofinite fold types in data/enumerations/roulleau_ud8d4_all_subgroup_folds.json",
    ),
    _case(
        "R14_UE8A1x4", 14, "U+E8(-1)+A1(-1)^4", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{A}_{1}^{\oplus4}",
        27, 16, "C2xC2", "t", (9828, 10185), (10966, 10966),
        "basis_weighted_graph_plus_explicit_coordinates_and_pair_relations_for_all_extra_curves",
        "direct_curve_enumeration_and_completeness_test",
        "represented_exact_full_nef_stock",
        "exact coordinate/equivalence stock, all 98 diagram subgroups, and four cofinite fold diagrams in data/enumerations/roulleau_r14_ue8a1x4_all_subgroup_folds.json",
    ),
    _case(
        "R15_UE8D4A1", 15, "U+E8(-1)+D4(-1)+A1(-1)", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{D}_{4}\oplus\mathbf{A}_{1}",
        21, 8, "C2xC2", "tn", (10186, 10324), (10938, 10938),
        "complete_weighted_wall_matrix_already_archived", "published_complete_curve_configuration",
        "represented_exact_full_nef_stock", "data/known_examples/k3_high_rank_nef_chambers.json rank-15 stock",
    ),
    _case(
        "R16_UE8D6", 16, "U+E8(-1)+D6(-1)", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{D}_{6}",
        19, 4, "C2xC2", "tn", (10325, 10425), (10942, 10942),
        "complete_weighted_wall_matrix_already_archived", "published_complete_curve_configuration",
        "represented_exact_full_nef_stock", "data/known_examples/k3_high_rank_nef_chambers.json rank-16 stock",
    ),
    _case(
        "R17_UE8E7", 17, "U+E8(-1)+E7(-1)", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{E}_{7}",
        19, 2, "C2xC2", "tn", (10426, 10520), (10946, 10946),
        "complete_weighted_wall_matrix_already_archived", "published_complete_curve_configuration",
        "represented_exact_full_nef_stock", "data/known_examples/k3_high_rank_nef_chambers.json rank-17 stock",
    ),
    _case(
        "R18_UE8E8", 18, "U+E8(-1)^2", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{E}_{8}",
        19, 1, "C2xC2", "tn", (10521, 10611), (10950, 10950),
        "complete_weighted_wall_matrix_already_archived", "published_complete_curve_configuration",
        "represented_exact_full_nef_stock", "data/known_examples/k3_high_rank_nef_chambers.json rank-18 stock",
    ),
    _case(
        "R19_UE8E8A1", 19, "U+E8(-1)^2+A1(-1)", r"U\oplus\mathbf{E}_{8}\oplus\mathbf{E}_{8}\oplus\mathbf{A}_{1}",
        24, 2, "S3xC2", "tn", (10612, 10707), (10954, 10954),
        "complete_weighted_wall_matrix_already_archived", "published_complete_curve_configuration",
        "represented_exact_full_nef_stock", "data/known_examples/k3_high_rank_nef_chambers.json rank-19 stock",
    ),
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _check_line_range(lines: list[str], case: dict) -> None:
    start, end = case["source"]["section_lines"]
    excerpt = "\n".join(lines[start - 1 : end])
    assert case["tex_lattice_token"] in excerpt, case["key"]
    wall_tokens = [
        f"contains ${case['wall_count']}$",
        f"contains ${case['wall_count']}$ $(-2)$-curves",
        f"contains ${case['wall_count']}$ $\\cu$-curves",
    ]
    numbered_list_witness = (
        "$(-2)$-curves" in excerpt
        and f"A_{{{case['wall_count']}}}" in excerpt
    )
    assert any(token in excerpt for token in wall_tokens) or numbered_list_witness, case["key"]
    tstart, tend = case["source"]["table_lines"]
    table_excerpt = "\n".join(lines[tstart - 1 : tend])
    assert case["tex_table_lattice_token"] in table_excerpt, case["key"]
    assert f"${case['wall_count']}$" in table_excerpt, case["key"]


def check_sources(tex: Path, examples: Path, algo: Path) -> dict:
    assert tex.exists() and examples.exists() and algo.exists()
    assert sha256(tex) == SOURCE_HASHES["tex"]
    assert sha256(examples) == SOURCE_HASHES["examples_magma"]
    assert sha256(algo) == SOURCE_HASHES["algorithm_magma"]
    tex_lines = tex.read_text().splitlines()
    for case in CASES:
        _check_line_range(tex_lines, case)
    examples_text = examples.read_text()
    assert 'Lattice ["U","D4","A1","A1","A1","A1","A1"]' in examples_text
    assert "We obtain 90 (-2)-curves" in examples_text
    assert 'Lattice "U","D8","A1","A1","A1"' in examples_text
    assert "SE:=[17,16,15,6,7,8" in examples_text
    algo_text = algo.read_text()
    for function in (
        "SmoothRationalCurves",
        "SmoothRationalCurvesWithKnowledge",
        "IsAmplePol",
        "IntersectionForCartan",
        "IntersectionPolarizationForCartan",
        "MatrixOfIntersections",
    ):
        assert f"function {function}" in algo_text
    return {"checked": True, "all_checks_pass": True}


def build_artifact(source_check: dict | None = None) -> dict:
    rank11_14 = [c for c in CASES if 11 <= c["picard_rank"] <= 14]
    represented = [
        c for c in CASES
        if c["current_machine_ledger_status"] == "represented_exact_full_nef_stock"
    ]
    queue = [
        c for c in rank11_14
        if c["current_machine_ledger_status"] == "scouted_not_yet_machine_represented"
    ]
    formula_cases = [
        c for c in rank11_14 if "plus_" in c["wall_data_extraction"]
    ]
    assert len(CASES) == 19
    assert len(rank11_14) == 14
    assert len(represented) == 19
    assert len(queue) == 0
    assert len(formula_cases) == 5
    assert sum(c["wall_count"] for c in rank11_14) == 375
    assert all(c["hyperbolic_dimension"] == c["picard_rank"] - 1 for c in CASES)
    return {
        "schema": "coxeter.roulleau_rank11_14_source_stock_inventory.v1",
        "artifact_id": "roulleau_rank11_14_source_stock_inventory",
        "all_checks_pass": True,
        "source": {
            "title": "An atlas of K3 surfaces with finite automorphism group",
            "author": "Xavier Roulleau",
            "publication": "EPIGA 6 (2022), article 19",
            "classification_scope": "118 finite-automorphism complex K3 Neron-Severi lattices of rank at least 3",
            "local_paths": {
                "pdf": "/private/tmp/roulleau-k3-atlas.pdf",
                "source_tar": "/private/tmp/roulleau_k3_atlas_source.tar",
                "tex": str(DEFAULT_TEX),
                "examples_magma": str(DEFAULT_EXAMPLES),
                "algorithm_magma": str(DEFAULT_ALGO),
            },
            "sha256": SOURCE_HASHES,
            "source_conventions_lines": [527, 565],
            "completeness_discussion_lines": [729, 820],
            "rank_at_least_11_table_lines": [10906, 10966],
            "algorithm_function_lines": {
                "SmoothRationalCurves": 26,
                "SmoothRationalCurvesWithKnowledge": 102,
                "IsAmplePol": 283,
                "IntersectionForCartan": 300,
                "IntersectionPolarizationForCartan": 347,
                "MatrixOfIntersections": 378,
            },
        },
        "denominator": {
            "rank_at_least_11_stock_count": len(CASES),
            "rank11_14_stock_count": len(rank11_14),
            "rank15_19_already_archived_stock_count": 5,
            "rank11_14_total_wall_count": sum(c["wall_count"] for c in rank11_14),
            "rank11_14_direct_complete_graph_count": 9,
            "rank11_14_basis_graph_plus_formula_count": len(formula_cases),
            "rank11_14_induced_chamber_completion_count": sum(
                c["source_completeness_method"] == "induced_chamber_completion_inside_II_1_17"
                for c in rank11_14
            ),
            "already_exactly_machine_represented_count": len(represented),
            "new_transcription_queue_count": len(queue),
            "new_transcription_queue_keys": [c["key"] for c in queue],
        },
        "mathematical_status": {
            "full_nef_chambers_are_positive_examples": True,
            "volume_type": "finite-volume_noncompact",
            "reason": (
                "finite automorphism group makes the K3 Weyl group finite index in the "
                "arithmetic orthogonal group; every listed lattice contains U and hence has cusps"
            ),
            "priority_warning": (
                "the full nef chambers are published known polytopes, not new abstract examples; "
                "the transcription queue is closed and the two new 59-wall folds have only "
                "bounded literal-priority claims"
            ),
        },
        "source_check": source_check or {"checked": False},
        "stocks": CASES,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check-source", action="store_true")
    parser.add_argument("--tex", type=Path, default=DEFAULT_TEX)
    parser.add_argument("--examples", type=Path, default=DEFAULT_EXAMPLES)
    parser.add_argument("--algo", type=Path, default=DEFAULT_ALGO)
    args = parser.parse_args()
    source_check = (
        check_sources(args.tex, args.examples, args.algo) if args.check_source else None
    )
    artifact = build_artifact(source_check)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact["denominator"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
