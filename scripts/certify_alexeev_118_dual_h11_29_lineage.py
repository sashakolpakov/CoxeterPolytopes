#!/usr/bin/env python3
"""Exact duality and bounded lineage audit for the Alexeev H11/29 catch.

This certificate proves that the two 29-facet sections are one abstract
hyperbolic Coxeter polytope by matching every normalized squared Gram entry.
It also performs four deliberately bounded comparisons: the current project
H11 positive ledger, every full Alexeev H11 matrix in the pinned ancillary
file, all 27 Scharlau--Kirschmer List11 entries, and the pinned ordinary and
index-two-doubling CoxIter snapshots.  Absence from these corpora is not a
worldwide priority theorem.
"""

from __future__ import annotations

from collections import Counter
import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Any

import sympy as sp

from scripts.certify_alexeev_20_2_1_15_9_1_lineage_bounds import (
    _rational_comparisons,
    _scan_coxiter,
)


ROOT = Path(__file__).resolve().parents[1]
CONSTRUCTION = ROOT / "data/enumerations/alexeev_118_dual_h11_29_construction.json"
VINBERG = ROOT / "data/enumerations/alexeev_118_dual_h11_29_vinberg.json"
STABLE = ROOT / "data/enumerations/alexeev_118_dual_h11_29_stable_overlattices.json"
ALEXEEV_MATRIX_SOURCE = Path(
    "/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt"
)
ALEXEEV_TEX_SOURCE = Path("/private/tmp/alexeev_2209_source/g0.tex")
SK_LIST11 = Path("/private/tmp/Scharlau-Kirschmer-Lorentz/List11")
OUTPUT = ROOT / "data/enumerations/alexeev_118_dual_h11_29_lineage.json"
SCHEMA = "coxeter.alexeev_118_dual_h11_29.lineage.v1"


SOURCE_TO_DUAL_FACET_PERMUTATION_ONE_BASED = [
    1,
    19,
    2,
    18,
    10,
    11,
    15,
    3,
    4,
    12,
    27,
    7,
    6,
    22,
    8,
    24,
    14,
    13,
    17,
    5,
    9,
    25,
    20,
    23,
    29,
    16,
    26,
    28,
    21,
]
NORMALIZED_SQUARED_GRAM_SHA256 = (
    "348ef97ede4d2cf17719e6f434dc96539e446db3a6f2ec943e64e3bb7de27e53"
)


PROJECT_H11_SOURCES = {
    "F3-C2-H11-206": (
        ROOT / "data/enumerations/ss_fermat_quartic_char3_c2_h11_vinberg_certificate.json",
        206,
    ),
    "F3-V4-H11-328": (
        ROOT / "data/enumerations/ss_fermat_quartic_char3_v4_h11_exact_cone_certificate.json",
        328,
    ),
    "A20-C2b-H11-20": (
        ROOT / "data/enumerations/alexeev_20_2_1_15_9_1_fixed_section_census.json",
        20,
    ),
    "SS2-3B-H11-22": (
        ROOT / "data/enumerations/borcherds_char2_3b_h11_22_certificate.json",
        22,
    ),
    "A18-C2a-H11-28": (
        ROOT / "data/enumerations/alexeev_remaining_high_fixed_section_construction.json",
        28,
    ),
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _payload_sha256(payload: Any) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(raw).hexdigest()


def _matrix(payload: list[list[str | int]]) -> sp.Matrix:
    return sp.Matrix([[sp.Rational(value) for value in row] for row in payload])


def _normalized_squared_gram(gram: sp.Matrix) -> list[list[str]]:
    if any(
        gram[left, right] > 0
        for left in range(gram.rows)
        for right in range(left)
    ):
        raise ArithmeticError("positive-wall Gram has a positive off-diagonal")
    return [
        [
            str(
                sp.factor(
                    sp.Rational(gram[left, right]) ** 2
                    / (sp.Rational(gram[left, left]) * sp.Rational(gram[right, right]))
                )
            )
            for right in range(gram.rows)
        ]
        for left in range(gram.rows)
    ]


def _duality(construction: dict[str, Any]) -> dict[str, Any]:
    cases = {case["id"]: case for case in construction["cases"]}
    source = _matrix(cases["A16_6_1_C01"]["restricted_mirrors"]["positive_wall_gram"])
    dual = _matrix(cases["A16_12_1_C01"]["restricted_mirrors"]["positive_wall_gram"])
    source_normalized = _normalized_squared_gram(source)
    dual_normalized = _normalized_squared_gram(dual)
    mapping = [value - 1 for value in SOURCE_TO_DUAL_FACET_PERMUTATION_ONE_BASED]
    reordered = [
        [dual_normalized[mapping[left]][mapping[right]] for right in range(29)]
        for left in range(29)
    ]
    if not (
        sorted(mapping) == list(range(29))
        and source_normalized == reordered
        and _payload_sha256(source_normalized) == NORMALIZED_SQUARED_GRAM_SHA256
    ):
        raise ArithmeticError("the exact dual 29-wall Gram identification changed")
    lower_histogram = Counter(
        source_normalized[left][right]
        for left in range(29)
        for right in range(left)
    )
    lattices = {
        identifier: {
            "determinant": cases[identifier]["saturated_fixed_lattice"]["determinant"],
            "smith_diagonal_absolute": cases[identifier]["saturated_fixed_lattice"]
            ["smith_diagonal_absolute"],
        }
        for identifier in cases
    }
    return {
        "canonical_realization": "A16_6_1_C01 active-facet order",
        "canonical_normalized_squared_gram": source_normalized,
        "canonical_normalized_squared_gram_sha256": NORMALIZED_SQUARED_GRAM_SHA256,
        "source_to_dual_facet_permutation_one_based": (
            SOURCE_TO_DUAL_FACET_PERMUTATION_ONE_BASED
        ),
        "all_841_normalized_squared_entries_match": True,
        "normalized_lower_triangle_histogram": dict(sorted(lower_histogram.items())),
        "same_abstract_hyperbolic_Coxeter_polytope": True,
        "fixed_lattice_realizations": lattices,
        "fixed_lattices_not_integrally_isometric": (
            lattices["A16_6_1_C01"]["determinant"]
            != lattices["A16_12_1_C01"]["determinant"]
        ),
    }


def _project_positive_shapes() -> dict[str, Any]:
    records = []
    for identifier, (path, expected_facets) in PROJECT_H11_SOURCES.items():
        payload = json.loads(path.read_text(encoding="utf-8"))
        if identifier == "F3-C2-H11-206":
            dimension, facets = payload["dimension"], payload["facet_count"]
        elif identifier == "F3-V4-H11-328":
            dimension = payload["cone"]["hyperbolic_dimension"]
            facets = payload["cone"]["facet_count"]
        elif identifier == "SS2-3B-H11-22":
            dimension = payload["decision"]["hyperbolic_dimension"]
            facets = payload["decision"]["facet_count"]
        else:
            cases = payload["cases"]
            target = next(
                case
                for case in cases
                if case["case"]["hyperbolic_dimension"] == 11
                and case["exact_reduction"]["active_hyperbolic_facet_count"]
                == expected_facets
            )
            dimension = target["case"]["hyperbolic_dimension"]
            facets = target["exact_reduction"]["active_hyperbolic_facet_count"]
        if (dimension, facets) != (11, expected_facets):
            raise ArithmeticError(f"project comparison target changed: {identifier}")
        records.append(
            {
                "id": identifier,
                "dimension": dimension,
                "facet_count": facets,
                "archive_path": str(path.relative_to(ROOT)),
                "archive_sha256": _sha256(path),
            }
        )
    return {
        "records": records,
        "facet_count_histogram": {
            str(key): value
            for key, value in sorted(Counter(row["facet_count"] for row in records).items())
        },
        "no_prior_project_positive_has_29_facets": all(
            row["facet_count"] != 29 for row in records
        ),
    }


def _alexeev_full_h11() -> dict[str, Any]:
    source = ALEXEEV_MATRIX_SOURCE.read_text(encoding="utf-8")
    keys = [
        tuple(map(int, groups))
        for groups in re.findall(
            r"coxeter_matrix\[\((\d+),(\d+),(\d+)\)\]\s*=", source
        )
    ]
    rank12 = [key for key in keys if key[0] == 12]
    marker = "coxeter_matrix[(12,10,1)] = matrix(["
    start = source.index(marker) + len(marker) - 1
    depth = 0
    end = start
    while end < len(source):
        if source[end] == "[":
            depth += 1
        elif source[end] == "]":
            depth -= 1
            if depth == 0:
                end += 1
                break
        end += 1
    import ast

    rows = ast.literal_eval(source[start:end])
    if rank12 != [(12, 10, 1)] or len(rows) != 14:
        raise ArithmeticError("the pinned Alexeev full H11 inventory changed")
    return {
        "source_path": str(ALEXEEV_MATRIX_SOURCE),
        "source_sha256": _sha256(ALEXEEV_MATRIX_SOURCE),
        "full_rank_12_keys": [list(key) for key in rank12],
        "full_H11_chambers": [{"key": [12, 10, 1], "facet_count": 14}],
        "no_full_Alexeev_H11_chamber_has_29_facets": True,
    }


def _scharlau_kirschmer() -> dict[str, Any]:
    source = SK_LIST11.read_text(encoding="utf-8")
    rows = re.findall(
        r'<\s*(\d+)\s*,\s*"([^"]+)"\s*,\s*\[\s*(\d+)', source
    )
    facets = [int(row[2]) for row in rows]
    if len(rows) != 27 or 29 in facets:
        raise ArithmeticError("Scharlau--Kirschmer List11 comparison changed")
    return {
        "source_path": str(SK_LIST11),
        "source_sha256": _sha256(SK_LIST11),
        "entry_count": len(rows),
        "full_chamber_facet_count_histogram": {
            str(key): value for key, value in sorted(Counter(facets).items())
        },
        "no_List11_full_chamber_has_29_facets": True,
    }


def _alexeev_published_duality() -> dict[str, Any]:
    source = ALEXEEV_TEX_SOURCE.read_text(encoding="utf-8")
    checks = {
        "both_parent_subsections_present": (
            "$(16,6,1)$" in source and "$(16,12,1)$" in source
        ),
        "both_parent_chambers_have_118_roots": source.count("$118$ roots") >= 2,
        "parent_group_S6_stated": source.count(r"\Aut\Gamma_r = S_6") >= 2,
        "doubled_dual_relation_stated": (
            "If $L=(16,6,1)$ then $(16,12,1)$ is" in source
            and "the even sublattice of $L^*(2)$" in source
        ),
        "derived_29_wall_section_not_stated": (
            "29-facet" not in source and "29 facet" not in source
        ),
    }


def _normalized_integral_lineage(stable: dict[str, Any]) -> dict[str, Any]:
    """Freeze the stronger literal lineage proved by the companion audit."""

    if not stable["all_checks_pass"]:
        raise ArithmeticError("stable-overlattice companion certificate failed")
    common = stable["common_normalized_reflection_system"]
    maximal = stable["maximal_reflection_lineage"]
    collision = stable["same_header_collision_with_A18_4_0_C03"]
    checks = {
        "same_normalized_integral_lattice": (
            common["identity"] == "A^T (G_12/2) A = M_6"
            and common["checks"]["A_integral"]
            and common["checks"]["A_unimodular"]
            and common["checks"]["A_transports_normalized12_to_maximal6"]
        ),
        "same_29_wall_subgroup": (
            common["checks"]["all_29_oriented_walls_map_exactly"]
            and common["checks"]["all_29_reflections_conjugate_exactly"]
            and len(common["source_root_factors_before_A"]) == 29
        ),
        "known_containing_full_reflection_chamber": (
            maximal["scharlau_kirschmer_list"] == "List11"
            and maximal["entry_one_based"] == 5
            and maximal["genus_symbol"] == "I_{11,1}(2^{2})"
            and maximal["full_chamber_facets"] == 14
            and all(maximal["checks"].values())
        ),
        "not_a_new_arithmetic_commensurability_class": (
            maximal["new_arithmetic_commensurability_class"] is False
            and stable["summary"]["new_arithmetic_commensurability_class_count"]
            == 0
        ),
        "A18_C03_integral_collision_excluded": (
            all(collision["checks"].values())
            and collision["checks"]["different_discriminant_quadratic_histograms"]
            and collision["checks"]["different_maximal_stable_determinants"]
        ),
        "only_numerical_full_group_index_open": (
            stable["summary"]["numerical_full_group_index_open"]
            and maximal["section_to_full_index_status"].startswith("finite;")
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"normalized integral lineage changed: {checks}")
    return {
        "normalizations": {
            "A16_6_1_C01": common["A16_6_model"],
            "A16_12_1_C01": common["A16_12_model"],
        },
        "integral_isometry_identity": common["identity"],
        "source_to_dual_facet_permutation_one_based": common[
            "source_to_dual_facet_permutation_one_based"
        ],
        "source_root_factors_before_isometry": common[
            "source_root_factors_before_A"
        ],
        "conclusion": common["conclusion"],
        "containing_full_reflection_chamber": {
            "census": maximal["scharlau_kirschmer_list"],
            "entry_one_based": maximal["entry_one_based"],
            "genus_symbol": maximal["genus_symbol"],
            "facet_count": maximal["full_chamber_facets"],
            "finite_vertex_count": maximal["full_chamber_finite_vertices"],
            "cusp_count": maximal["full_chamber_cusps"],
            "section_to_full_index_status": maximal[
                "section_to_full_index_status"
            ],
        },
        "A18_4_0_C03_noncollision": {
            "conclusion": collision["conclusion"],
            "A16_6_discriminant_q_mod_2_histogram": collision[
                "A16_6_discriminant_q_mod_2_histogram"
            ],
            "A18_C03_discriminant_q_mod_2_histogram": collision[
                "A18_C03_discriminant_q_mod_2_histogram"
            ],
        },
        "checks": checks,
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Alexeev source-text lineage changed: {checks}")
    return {
        "source_path": str(ALEXEEV_TEX_SOURCE),
        "source_sha256": _sha256(ALEXEEV_TEX_SOURCE),
        "published_parent_relation": (
            "(16,12,1) is the even sublattice of (16,6,1)^*(2); "
            "their 118-wall parabolic diagrams correspond with -2/-4 roots switched"
        ),
        "checks": checks,
        "interpretation": (
            "the parent duality is published; the common transposition-fixed "
            "H11/29 chamber is the new derived object tested here"
        ),
    }


def build(
    construction_path: Path = CONSTRUCTION,
    vinberg_path: Path = VINBERG,
    stable_path: Path = STABLE,
) -> dict[str, Any]:
    construction = json.loads(construction_path.read_text(encoding="utf-8"))
    vinberg = json.loads(vinberg_path.read_text(encoding="utf-8"))
    stable = json.loads(stable_path.read_text(encoding="utf-8"))
    if not construction["all_checks_pass"] or not vinberg["all_checks_pass"]:
        raise ArithmeticError("upstream construction/volume certificate failed")
    duality = _duality(construction)
    project = _project_positive_shapes()
    alexeev = _alexeev_full_h11()
    sk = _scharlau_kirschmer()
    coxiter = _scan_coxiter({(11, 29)})
    target = coxiter["targets"]["H11_facets_29"]
    parent = _alexeev_published_duality()
    integral_lineage = _normalized_integral_lineage(stable)

    rational_rows = {}
    for case in construction["cases"]:
        lattice = _matrix(case["saturated_fixed_lattice"]["gram"])
        comparisons = _rational_comparisons(lattice)
        rational_rows[case["id"]] = {
            key: value["rationally_similar"] for key, value in comparisons.items()
        }
    checks = {
        "one_abstract_polytope_two_integral_realizations": (
            duality["same_abstract_hyperbolic_Coxeter_polytope"]
            and duality["fixed_lattices_not_integrally_isometric"]
        ),
        "both_volume_audits_pass": vinberg["summary"]["finite_volume_count"] == 2,
        "not_a_full_parent_chamber": all(
            case["case"]["source_wall_count"] == 118
            and case["exact_reduction"]["active_hyperbolic_facet_count"] == 29
            for case in construction["cases"]
        ),
        "absent_from_prior_project_H11_positives_by_facet_count": project[
            "no_prior_project_positive_has_29_facets"
        ],
        "absent_from_full_Alexeev_H11_by_facet_count": alexeev[
            "no_full_Alexeev_H11_chamber_has_29_facets"
        ],
        "absent_from_SK_List11_by_facet_count": sk[
            "no_List11_full_chamber_has_29_facets"
        ],
        "absent_from_pinned_CoxIter_snapshots_by_facet_count": target[
            "absent_from_both_snapshots_by_dimension_and_facet_count"
        ],
        "both_in_odd_unimodular_rational_class": all(
            row["f1_odd_unimodular"]
            and not row["f2"]
            and not row["f3_McLeod"]
            for row in rational_rows.values()
        ),
        "same_normalized_integral_lattice_and_29_reflections": (
            integral_lineage["checks"]["same_normalized_integral_lattice"]
            and integral_lineage["checks"]["same_29_wall_subgroup"]
        ),
        "known_full_reflection_lineage_not_new_arithmetic_class": (
            integral_lineage["checks"][
                "known_containing_full_reflection_chamber"
            ]
            and integral_lineage["checks"][
                "not_a_new_arithmetic_commensurability_class"
            ]
        ),
        "A18_C03_same_header_collision_excluded": integral_lineage["checks"][
            "A18_C03_integral_collision_excluded"
        ],
    }
    if not all(checks.values()):
        raise ArithmeticError(f"bounded H11/29 lineage checks failed: {checks}")
    return {
        "schema": SCHEMA,
        "source": {
            "construction_path": str(construction_path.relative_to(ROOT)),
            "construction_sha256": _sha256(construction_path),
            "vinberg_path": str(vinberg_path.relative_to(ROOT)),
            "vinberg_sha256": _sha256(vinberg_path),
            "stable_overlattice_path": str(stable_path.relative_to(ROOT)),
            "stable_overlattice_sha256": _sha256(stable_path),
        },
        "duality": duality,
        "normalized_integral_lineage": integral_lineage,
        "published_parent_duality": parent,
        "rational_similarity": rational_rows,
        "bounded_comparisons": {
            "prior_project_H11_positives": project,
            "full_Alexeev_H11_chambers": alexeev,
            "Scharlau_Kirschmer_List11": sk,
            "CoxIter": coxiter,
        },
        "priority_statement": {
            "proved": (
                "a finite-volume noncompact H11/29 Coxeter polytope with two "
                "exact dual Alexeev fixed-section realizations which normalize "
                "to the same integral I_(11,1)(2^2) lattice and the same "
                "29-reflection subgroup; its containing full reflection chamber "
                "is Scharlau--Kirschmer List11 entry 5; the literal 29-wall "
                "chamber is absent from the four pinned comparison denominators; "
                "the apparent A18-C03 determinant/Smith collision is excluded by "
                "the discriminant quadratic form"
            ),
            "not_proved": (
                "worldwide literal novelty, the numerical index in the known "
                "full reflection group, or realization of the diagram "
                "involution by an automorphism of an algebraic variety"
            ),
        },
        "checks": checks,
        "all_checks_pass": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--construction", type=Path, default=CONSTRUCTION)
    parser.add_argument("--vinberg", type=Path, default=VINBERG)
    parser.add_argument("--stable", type=Path, default=STABLE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check-archive", action="store_true")
    args = parser.parse_args()
    payload = build(args.construction, args.vinberg, args.stable)
    if args.check_archive:
        if payload != json.loads(args.output.read_text(encoding="utf-8")):
            raise ArithmeticError("rebuilt H11/29 lineage archive changed")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["checks"], sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
